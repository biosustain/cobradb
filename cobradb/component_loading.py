# -*- coding: utf-8 -*-

from pathlib import Path
from typing import List, Set
from cobradb.models import AlreadyLoadedError, Gene, Genome, Chromosome, Synonym
from cobradb.util import scrub_gene_id, get_or_create_data_source, get_or_create, timing
from cobradb import ncbi_data

from sqlalchemy import select, func
import logging
import gzip


class BadGenomeError(Exception):
    pass


def _load_gb_file(genbank_file_handle):
    """Load the Genbank file.

    Arguments
    ---------

    genbank_file_handle: The handle to the genbank file.

    """
    # imports
    from Bio import SeqIO

    # load the genbank file
    logging.debug("Loading file: %s" % genbank_file_handle.name)
    try:
        gb_file = SeqIO.parse(genbank_file_handle, "gb")
    except IOError:
        raise BadGenomeError("File '%s' not found" % genbank_file_handle.name)
    except Exception as e:
        raise BadGenomeError(
            'BioPython failed to parse %s with error "%s"'
            % (genbank_file_handle.name, getattr(e, "message", "?"))
        )
    return gb_file


def collect_gene_synonym(synonym_collection, gene_db_id, synonym, data_source_id):
    if not data_source_id in synonym_collection:
        synonym_collection[data_source_id] = set()
    synonym_collection[data_source_id].add((synonym, gene_db_id))


def insert_synonyms(session, synonym_collection):
    # syns = []
    for data_source_id, synonyms in synonym_collection.items():
        data_source_db = get_or_create_data_source(session, data_source_id)
        for syn, gene_db_id in synonyms:
            synonym = Synonym(
                type="gene",
                ome_id=gene_db_id,
                synonym=syn,
            )
            data_source_db.synonyms.append(synonym)


def _nonempty_str(s):
    s = s.strip()
    return None if s == "" else s


def _get_qual(feat, name):
    """Get a non-null attribute from the feature."""
    try:
        qual = feat.qualifiers[name]
    except KeyError:
        return []

    return [y for y in (_nonempty_str(x) for x in qual) if y is not None]


def _get_first_qual(feat, name):
    """Get a non-null attribute from the feature."""
    try:
        qual = feat.qualifiers[name]
    except KeyError:
        return None

    return _nonempty_str(qual[0])


def _get_geneid(feature):
    """Get the value of GeneID from db_xref, or else return None."""
    db_xref = _get_qual(feature, "db_xref")
    if not db_xref:
        return None
    for ref in db_xref:
        splitrefs = [x.strip() for x in ref.split(":")]
        if len(splitrefs) == 2 and splitrefs[0].lower() == "geneid":
            return splitrefs[1]
    return None


@timing
def load_assembly(assembly_id, assembly_path, chromosome_accessions, session):
    """Load the genome and chromosomes."""

    accession_type, accession_value = assembly_id

    # check that the genome doesn't already exist
    if (
        session.scalar(
            select(func.count(Genome.id))
            .filter(Genome.accession_type == accession_type)
            .filter(Genome.accession_value == accession_value)
        )
        > 0
    ):
        raise AlreadyLoadedError(f"Assembly {assembly_id} already loaded")

    logging.debug(f"Adding new genome: {assembly_id}")
    genome_db = Genome(accession_type=accession_type, accession_value=accession_value)
    session.add(genome_db)
    session.commit()
    genome_db_id = genome_db.id

    if assembly_path is None:
        logging.warning(f"No assembly file found for {assembly_id}")
    else:
        open_f = open
        if Path(assembly_path).suffix == ".gz":
            open_f = gzip.open
        with open_f(assembly_path, "rt") as f:
            gb_file = _load_gb_file(f)
            for i, record in enumerate(gb_file):
                if (
                    chromosome_accessions is not None
                    and record.id not in chromosome_accessions
                ):
                    logging.warning(f"Skipping chromosome [{i+1} of ?] {record.id}")
                    continue
                logging.info(f"Loading chromosome [{i+1} of ?] {record.id}")
                load_chromosome(session, record, genome_db_id)

    genome_db = session.get(Genome, genome_db_id)
    if not genome_db.organism and accession_type == "ncbi_assembly":
        organism_info = ncbi_data.get_organism_for_ncbi_assembly_accession(
            genome_db.accession_value
        )
        if organism_info is not None:
            new_organism, new_tax_id, new_strain = organism_info
            genome_db.organism = new_organism
            genome_db.taxon_id = new_tax_id
            genome_db.strain = new_strain
    session.commit()
    session.close()


def first_uppercase(val):
    if val is None:
        return val
    return str(val[0]).upper()


@timing
def load_chromosome(session, record, genome_db_id: int):
    genome_db = session.get(Genome, genome_db_id)
    chromosome_db = session.scalars(
        select(Chromosome)
        .filter(Chromosome.ncbi_accession == record.id)
        .filter(Chromosome.genome_id == genome_db.id)
        .limit(1)
    ).first()
    if not chromosome_db:
        logging.debug("Loading new chromosome: {}".format(record.id))
        chromosome_db = Chromosome(ncbi_accession=record.id)
        genome_db.chromosomes.append(chromosome_db)
        session.commit()
    else:
        logging.debug("Chromosome already loaded: %s" % record.id)

    # update genome
    if genome_db.organism is None:
        logging.warning(f"Organism: {record.annotations['organism']}")
        genome_db.organism = record.annotations["organism"]
        session.commit()

    bigg_id_warnings = 0
    duplicate_genes_warnings = 0
    warning_num = 5
    synonym_collection = {}
    added_bigg_ids = set(
        session.scalars(
            select(Gene.bigg_id).filter(Gene.chromosome_id == chromosome_db.id)
        ).all()
    )
    for i, feature in enumerate(record.features):
        # update genome with the source information
        if feature.type == "source":
            if genome_db.strain is None:
                strain = _get_first_qual(feature, "strain")
                if strain is not None:
                    genome_db.strain = strain
            if genome_db.taxon_id is None:
                db_xref = _get_qual(feature, "db_xref")
                if db_xref is not None:
                    for ref in db_xref:
                        if "taxon" == ref.split(":")[0]:
                            genome_db.taxon_id = int(ref.split(":")[1])
                            break
            session.commit()

        # only read in CDSs
        if feature.type != "CDS":
            continue

        # bigg_id required
        bigg_id = None
        gene_name = None
        refseq_name = None
        locus_tag = None

        # get bigg_id if possible from locus_tag or and GeneID
        found_tag = _get_first_qual(feature, "locus_tag")
        found_gene_id = _get_geneid(feature)
        if found_tag is not None:
            locus_tag = found_tag
            bigg_id = scrub_gene_id(found_tag)
        elif found_gene_id is not None:
            bigg_id = scrub_gene_id(found_gene_id)

        # get name
        found_name = _get_first_qual(feature, "gene")
        if found_name is not None:
            gene_name = found_name
            refseq_name = found_name

        # warn about no locus_tag / bigg_id
        if gene_name is not None and bigg_id is None:
            if bigg_id_warnings <= warning_num:
                msg = (
                    "No locus_tag for gene. Using Gene name as bigg_id: %s" % gene_name
                )
                if bigg_id_warnings == warning_num:
                    msg += " (Warnings limited to %d)" % warning_num
                logging.warning(msg)
                bigg_id_warnings += 1
            bigg_id = scrub_gene_id(gene_name)
            gene_name = bigg_id
        elif bigg_id is None:
            logging.warning(
                (
                    "No locus_tag or gene name for gene %d in chromosome "
                    "%s" % (i, chromosome_db.ncbi_accession)
                )
            )
            continue

        gene_db = None
        if bigg_id in added_bigg_ids:
            if duplicate_genes_warnings <= warning_num:
                msg = "Duplicate genes %s on chromosome %s" % (
                    bigg_id,
                    chromosome_db.ncbi_accession,
                )
                if duplicate_genes_warnings == warning_num:
                    msg += " (Warnings limited to %d)" % warning_num
                logging.warning(msg)
                duplicate_genes_warnings += 1

            gene_db = session.scalars(
                select(Gene)
                .filter(Gene.bigg_id == bigg_id)
                .filter(Gene.chromosome_id == chromosome_db.id)
                .limit(1)
            ).first()
        else:
            # get the strand and positions
            strand = None
            if feature.location.strand == 1:
                strand = "+"
            elif feature.location.strand == -1:
                strand = "-"
            leftpos = int(feature.location.start) + 1
            rightpos = int(feature.location.end)

            dna_sequence = str(feature.extract(record.seq)).upper()
            protein_sequence = first_uppercase(
                feature.qualifiers.get("translation", None)
            )

            # finally, create the gene
            gene_db = Gene(
                bigg_id=bigg_id,
                locus_tag=locus_tag,
                name=gene_name,
                leftpos=leftpos,
                rightpos=rightpos,
                strand=strand,
                dna_sequence=dna_sequence,
                protein_sequence=protein_sequence,
                mapped_to_genbank=True,
                chromosome_id=chromosome_db.id,
            )
            session.add(gene_db)
            session.commit()
            added_bigg_ids.add(bigg_id)

        gene_db_id = gene_db.id
        # session.expunge(gene_db)
        # load the synonyms for the gene
        if locus_tag is not None:
            collect_gene_synonym(
                synonym_collection, gene_db_id, locus_tag, "refseq_locus_tag"
            )

        if refseq_name is not None:
            collect_gene_synonym(
                synonym_collection, gene_db_id, refseq_name, "refseq_name"
            )

        for ref in _get_qual(feature, "gene_synonym"):
            synonyms = [x.strip() for x in ref.split(";")]
            for syn in synonyms:
                collect_gene_synonym(
                    synonym_collection, gene_db_id, syn, "refseq_synonym"
                )

        for ref in _get_qual(feature, "db_xref"):
            splitrefs = [x.strip() for x in ref.split(":")]
            if len(splitrefs) == 2:
                collect_gene_synonym(
                    synonym_collection, gene_db_id, splitrefs[1], splitrefs[0]
                )

        for ref in _get_qual(feature, "old_locus_tag"):
            for syn in [x.strip() for x in ref.split(";")]:
                collect_gene_synonym(
                    synonym_collection, gene_db_id, syn, "refseq_old_locus_tag"
                )

        for ref in _get_qual(feature, "note"):
            for value in [x.strip() for x in ref.split(";")]:
                sp = value.split(":")
                if len(sp) == 2 and sp[0] == "ORF_ID":
                    collect_gene_synonym(
                        synonym_collection, gene_db_id, sp[1], "refseq_orf_id"
                    )
    session.commit()
    session.close()
    insert_synonyms(session, synonym_collection)
    session.commit()
    session.close()
