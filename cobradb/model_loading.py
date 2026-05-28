import json
from pathlib import Path
from typing import Any, Dict, Optional, Union
from sqlalchemy.orm import Session, contains_eager, joinedload, subqueryload, aliased
from cobradb.api import metabolites
from cobradb.api.bigg_ids import create_component_bigg_id
from cobradb.api.metabolites import (
    create_collection_specific_metabolite,
    get_or_create_compartmentalized_component,
    get_universal_component_by_bigg_id,
)
from cobradb.api.reactions import (
    ReactionParticipantInfo,
    get_or_create_reaction_for_universal_reaction,
    get_or_create_universal_reaction,
    get_universal_reaction_participants_list,
)

from cobradb.models import (
    AlreadyLoadedError,
    Chromosome,
    Compartment,
    CompartmentalizedComponent,
    Component,
    ComponentReferenceMapping,
    Gene,
    GeneReactionMatrix,
    Genome,
    Model,
    ModelCollection,
    ModelCompartmentalizedComponent,
    ModelCount,
    ModelGene,
    ModelReaction,
    Publication,
    PublicationModel,
    Reaction,
    ReactionMatrix,
    Synonym,
    Taxon,
    UniversalCompartmentalizedComponent,
    UniversalComponent,
    UniversalReaction,
)
from cobradb import settings
from cobradb import parse
from cobradb.api import utils
from cobradb.util import (
    format_formula,
    load_tsv,
    ref_str_to_tuple,
    scrub_name,
    check_none,
    timing,
)

from sqlalchemy import func, select, update
import re
import logging
from collections import defaultdict
import os
from difflib import SequenceMatcher


class GenbankNotFound(Exception):
    pass


def load_model_assemblies(
    model_genome_path: Union[str, os.PathLike], refseq_dir: Union[str, os.PathLike]
):
    model_genome_path = Path(model_genome_path)
    if model_genome_path.suffix == ".tsv":
        return load_model_assemblies_tsv(model_genome_path, refseq_dir)
    elif model_genome_path.suffix == ".json":
        return load_model_assemblies_json(model_genome_path, refseq_dir)
    raise ValueError("model_genome_path should be a .tsv or .json file")


def load_model_assemblies_json(
    model_genome_path: Union[str, os.PathLike], refseq_dir: Union[str, os.PathLike]
):
    with open(model_genome_path, "r") as f:
        data = json.load(f)
    models_list = data.get("models", [])
    collection_list = data.get("collections", [])
    assemblies = {}
    chromosomes = {}

    for model_entry in models_list:
        if "filename" not in model_entry:
            logging.error(f"Skipping, no 'filename' defined in entry '{model_entry}'.")
            continue

        pub_ref = model_entry.get("pub_ref")
        if pub_ref is not None:
            if isinstance(pub_ref, str):
                pub_ref = pub_ref.split(":", maxsplit=1)
        model_entry["pub_ref"] = pub_ref

        model_assembly = model_entry.get("assembly")
        if model_assembly is not None:
            if ":" in model_assembly:
                model_assembly = tuple(model_assembly.split(":", maxsplit=1))
            else:
                model_assembly = ("ncbi_assembly", model_assembly)
        if model_assembly is not None and model_assembly not in assemblies:
            if model_assembly[0] == "ncbi_assembly":
                assembly_file = (
                    Path(refseq_dir)
                    / "ncbi_dataset"
                    / "data"
                    / model_assembly[1]
                    / "genomic.gbff"
                )
                if (gz_file := assembly_file.with_suffix(".gbff.gz")).is_file():
                    assemblies[model_assembly] = str(gz_file)
                elif assembly_file.is_file():
                    assemblies[model_assembly] = str(assembly_file)
                else:
                    print(f"Assembly file not found: {assembly_file}")
                    assemblies[model_assembly] = None
            elif model_assembly[0] == "pankb_assembly":
                logging.warning(f"PanKB assembly: {model_assembly[1]}")
                assemblies[model_assembly] = None
            elif model_assembly[0] == "patric_assembly":
                logging.warning(f"PATRIC assembly: {model_assembly[1]}")
                assemblies[model_assembly] = None

        model_entry["assembly"] = model_assembly

        model_chromosomes = model_entry.get("chromosomes", [])
        for model_chromosome in model_chromosomes:
            if model_assembly is not None:
                if not model_assembly in chromosomes:
                    chromosomes[model_assembly] = set()
                chromosomes[model_assembly].add(model_chromosome)
        model_entry["chromosomes"] = model_chromosomes

    return models_list, assemblies, chromosomes, collection_list


def load_model_assemblies_tsv(
    model_genome_path: Union[str, os.PathLike], refseq_dir: Union[str, os.PathLike]
):
    lines = load_tsv(model_genome_path, required_column_num=None)
    models_list = []

    assemblies = {}

    def get_assembly_filenames(assembly_string):
        assembly_ids = [x.strip() for x in assembly_string.split(",")]
        r = []
        for assembly in assembly_ids:
            r.append(assembly)
            if assembly in assemblies:
                continue
            assembly_file = (
                Path(refseq_dir)
                / "ncbi_dataset"
                / "data"
                / assembly
                / "genomic.gbff.gz"
            )
            if assembly_file.is_file():
                assemblies[assembly] = str(assembly_file)
            else:
                print(f"Assembly file not found: {assembly_file}")
                assemblies[assembly] = None
        return r

    chromosomes = {}
    for line in lines:
        if len(line) > 4:
            line = line[:4]
        model_filename, pub_ref_string, assembly_ref_string, genome_ref_string_plus = (
            line
        )
        if assembly_ref_string is None:
            model_assemblies = None
        else:
            model_assemblies = get_assembly_filenames(assembly_ref_string)
        if genome_ref_string_plus is None:
            model_chromosomes = None
        else:
            model_chromosomes = set(
                [x.strip() for x in genome_ref_string_plus.split(",")]
            )
            if not model_assemblies is None:
                for assembly in model_assemblies:
                    if not assembly in chromosomes:
                        chromosomes[assembly] = set()
                    for chr in model_chromosomes:
                        chromosomes[assembly].add(chr)
        if pub_ref_string is None or pub_ref_string.strip() == "None":
            pub_ref = None
        else:
            pub_ref = ref_str_to_tuple(pub_ref_string)
        models_list.append(
            {
                "filename": model_filename,
                "pub_ref": pub_ref,
                "assemblies": model_assemblies,
                "assembly": None if model_assemblies is None else model_assemblies[0],
                "chromosomes": model_chromosomes,
                "genome_ref": (
                    None
                    if model_assemblies is None
                    else ("ncbi_assembly", model_assemblies[0])
                ),
            }
        )
    return models_list, assemblies, chromosomes, []


def get_model_list():
    """Get the models that are available, as SBML, in ome_data/models"""
    return [
        x.replace(".xml", "").replace(".mat", "")
        for x in os.listdir(settings.model_directory)
        if ".xml" in x or ".mat" in x
    ]


def check_for_model(name):
    """Check for model, case insensitive, and ignore periods and underscores"""

    def min_name(n):
        return n.lower().replace(".", "").replace(" ", "").replace("_", "")

    for x in get_model_list():
        if min_name(name) == min_name(x):
            return x
    return None


def _sim(name1, name2):
    """Return true if the names are similar"""
    clean1, clean2 = [n.lower().replace(" ", "") for n in [name1, name2]]
    return SequenceMatcher(None, clean1, clean2).ratio() > 0.7


def improve_name(session, db, new_name):
    """If the new_name is a better descriptive name for the reaction or metabolite,
    then update.

    """
    cur_name = db.name
    cur_id = db.id
    # New name is not None and not similar to bigg_id
    if (
        new_name is not None
        and not _sim(new_name, cur_id)
        and (cur_name is None or _sim(cur_name, cur_id))
    ):
        logging.debug("Replacing name %s with %s" % (cur_name, new_name))
        db.name = new_name
    # session.commit()


@timing
def load_model(
    session: Session,
    model_data: Dict[str, Any],
    model_dir: Optional[Union[str, os.PathLike]] = None,
):
    """Load a model into the database. Returns the bigg_id for the new model.

    Arguments
    ---------

    model_filepath: the path to the file where model is stored.

    pub_ref: a publication PMID or doi for the model, as a tuple like this:

        ('doi', '10.1128/ecosalplus.10.2.1')

        ('pmid', '21988831')

        Can be None

    genome_ref: A tuple specifying the genome accession type and value. The
    first element can be ncbi_accession, ncbi_assembly, or organism.

    session: An instance of Session.

    """
    # apply id normalization
    logging.debug("Parsing SBML")
    if model_dir is None:
        full_path = model_data["filename"]
    else:
        full_path = Path(model_dir) / model_data["filename"]
    model, old_parsed_ids = parse.load_and_normalize(full_path)

    if model_data.get("prefix") is not None:
        model.id = f"{model_data['prefix']}{model.id}"

    # Fallback: if model.id is empty in the model file, derive bigg_id from the filename
    if not model.id:
        from pathlib import Path as _P
        derived_id = _P(model_data["filename"]).stem.replace(".", "_")
        while derived_id.endswith(("_json", "_xml", "_mat", "_sbml")):
            import re as _re
            derived_id = _re.sub(r"_(json|xml|mat|sbml)$", "", derived_id)
        logging.warning(
            f"Empty model.id for {model_data['filename']}, deriving bigg_id={derived_id}"
        )
        model.id = derived_id

    model_bigg_id = model.id

    # check that the model doesn't already exist
    if (
        model_db := utils.get_object_by_bigg_id(session, model_bigg_id, Model)
    ) is not None:
        raise AlreadyLoadedError(
            f"Model {model_bigg_id} already loaded",
            bigg_id=model_db.bigg_id,
            db_id=model_db.id,
        )

    organism = None
    genome_id = None
    tax_id = None
    collection_id = None

    # check for a genome annotation for this model
    if (assembly := model_data.get("assembly")) is not None:
        genome_db = session.scalars(
            select(Genome)
            .filter(Genome.accession_type == assembly[0])
            .filter(Genome.accession_value == assembly[1])
            .limit(1)
        ).first()
        if genome_db is None:
            raise GenbankNotFound(
                "Genome for model {} not found with assembly {}".format(
                    model_bigg_id, assembly
                )
            )
        genome_id = genome_db.id
        organism = genome_db.organism
        tax_id = genome_db.taxon_id

    if (model_organism := model_data.get("organism")) is not None:
        organism = model_organism

    if (model_tax_id := model_data.get("tax_id")) is not None:
        tax_id = model_tax_id

    if organism is not None and tax_id is None:
        tax_ids = session.scalars(select(Taxon.id).filter(Taxon.name == organism)).all()
        if len(tax_ids) > 1:
            logging.error(
                f"Could not determine tax_id for {model_bigg_id}, because name '{organism}' is associated with multiple tax ids."
            )
        elif len(tax_ids) == 1:
            tax_id = tax_ids[0]

    if not organism and not genome_id:
        logging.info(
            "No Genome reference or organism provided for model {}".format(
                model_bigg_id
            )
        )

    collection_db = None
    if (model_collection_bigg_id := model_data.get("collection")) is not None:
        collection_db = utils.get_object_by_bigg_id(
            session, model_collection_bigg_id, ModelCollection
        )
    else:
        # Create single-model collection
        model_collection_bigg_id = model_bigg_id
    if collection_db is None:
        collection_db = ModelCollection(bigg_id=model_collection_bigg_id)
        session.add(collection_db)
        session.commit()
    collection_id = collection_db.id

    # Load the model objects. Remember: ORDER MATTERS! So don't mess around.
    logging.debug("Loading objects for model {}".format(model.id))
    published_filename = os.path.basename(model_data["filename"])
    model_db_id = load_new_model(
        session,
        model,
        genome_id,
        model_data.get("pub_ref"),
        published_filename,
        organism,
        tax_id=tax_id,
        collection_id=collection_id,
    )

    load_metabolites(
        session,
        model_db_id,
        model,
    )

    # # reactions
    load_reactions(
        session,
        model_db_id,
        model,
    )
    #
    # genes
    load_genes(session, model_db_id, model)

    # count model objects for the model summary web page
    load_model_count(session, model_db_id)

    session.commit()

    return model_bigg_id, model_db_id


@timing
def load_new_model(
    session,
    model,
    genome_db_id,
    pub_ref,
    published_filename,
    organism,
    collection_id,
    tax_id=None,
):
    """Load the model.

    Arguments:
    ---------

    session: A SQLAlchemy session.

    model: A COBRApy model.

    genome_db_id: The database ID of the genome. Can be None.

    pub_ref: a publication PMID or doi for the model, as a string like this:

        doi:10.1128/ecosalplus.10.2.1

        pmid:21988831

        Can be None

    organism: The organism. Can be None.

    Returns:
    -------

    The database ID of the new model row.

    """
    if tax_id is not None and session.get(Taxon, tax_id) is None:
        logging.error(f"Could not find tax id {tax_id} in database.")
        tax_id = None
    model_db = Model(
        bigg_id=model.id,
        genome_id=genome_db_id,
        published_filename=published_filename,
        organism=organism,
        taxon_id=tax_id,
        collection_id=collection_id,
    )
    session.add(model_db)
    if pub_ref is not None:
        session.commit()
        ref_type, ref_id = pub_ref
        publication_db = session.scalars(
            select(Publication)
            .filter(Publication.reference_type == ref_type)
            .filter(Publication.reference_id == ref_id)
            .limit(1)
        ).first()
        if publication_db is None:
            publication_db = Publication(reference_type=ref_type, reference_id=ref_id)
        publication_model_db = None
        if publication_db.id is not None and model_db.id is not None:
            publication_model_db = session.scalars(
                select(PublicationModel)
                .filter(PublicationModel.publication_id == publication_db.id)
                .filter(PublicationModel.model_id == model_db.id)
                .limit(1)
            ).first()
        if publication_model_db is None:
            publication_model_db = PublicationModel(publication=publication_db)
            model_db.publication_models.append(publication_model_db)
    session.commit()
    return model_db.id


@timing
def load_metabolites(session, model_db_id, model):
    """Load the metabolites as components and model components.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_id: The database ID for the model.

    model: The COBRApy model.

    old_metabolite_ids: A dictionary where keys are new IDs and values are old
    IDs for compartmentalized metabolites.
    Returns
    -------

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    final_metabolite_ids: A new dictionary where keys are original
    compartmentalized metabolite IDs from the model and values are the new
    compartmentalized metabolite IDs.

    """

    model_db = session.get(Model, model_db_id)
    collection_db = session.get(ModelCollection, model_db.collection_id)

    # for each metabolite in the model
    for metabolite in model.metabolites:
        # model_specific_id = metabolite.id
        # metabolite_id = model_specific_id.removeprefix(f"__{model_db.bigg_id}__")
        # metabolite_id = parse.remove_duplicate_tag(metabolite_id)
        # print(f"Metabolite: {model_specific_id} -> {metabolite_id}")
        metabolite_id = metabolite.id
        clean_metabolite_id = parse.remove_duplicate_tag(metabolite_id)
        print(f"Metabolite: {metabolite_id}")

        try:
            component_bigg_id, compartment_bigg_id = parse.split_compartment(
                clean_metabolite_id
            )
        except Exception:
            logging.error(
                (
                    "Could not find compartment for metabolite %s in"
                    "model %s, assuming compartment 'c'."
                    % (clean_metabolite_id, model_db.bigg_id)
                )
            )
            component_bigg_id = clean_metabolite_id
            compartment_bigg_id = "c"

        if (
            universal_component_db := get_universal_component_by_bigg_id(
                session, component_bigg_id, model_collection_id=collection_db.id
            )
        ) is not None:
            new_universal_bigg_id = universal_component_db.bigg_id
        else:
            new_universal_bigg_id = component_bigg_id

        print(f"# {new_universal_bigg_id}: {component_bigg_id}")
        print(f"## <> {universal_component_db}")

        # look for the formula in these places
        formula_fns = [
            lambda m: getattr(m, "formula", None),  # support cobra v0.3 and 0.4
            lambda m: m.notes.get("FORMULA", None),
            lambda m: m.notes.get("FORMULA1", None),
        ]
        # Cast to string, but not for None
        strip_str_or_none = lambda v: str(v).strip() if v is not None else None
        # Ignore the empty string
        ignore_empty_str = lambda s: s if s != "" else None
        # Use a generator for lazy evaluation
        values = (
            ignore_empty_str(strip_str_or_none(formula_fn(metabolite)))
            for formula_fn in formula_fns
        )
        # Get the first non-null result. Otherwise formula = None.
        formula = format_formula(next(filter(None, values), None))
        # Check for non-valid formulas
        if parse.invalid_formula(formula):
            logging.warning(
                "Invalid formula %s for metabolite %s in model %s"
                % (formula, metabolite_id, model.id)
            )
            formula = None
        _is_orig, formula = utils.fix_explicit_formula(formula, allow_R=True)

        # get charge
        try:
            charge = int(metabolite.charge)
            # check for float charge
            if charge != metabolite.charge:
                logging.warning(
                    "Could not load charge {} for {} in model {}".format(
                        metabolite.charge, metabolite_id, model.id
                    )
                )
                charge = None
        except Exception:
            if hasattr(metabolite, "charge") and metabolite.charge is not None:
                logging.debug(
                    "Could not convert charge to integer for metabolite {} in model {}: {}".format(
                        metabolite_id, model.id, metabolite.charge
                    )
                )
            charge = None

        if charge is None:
            charge = 0

        new_biggr_id = create_component_bigg_id(new_universal_bigg_id, charge=charge)

        print(f"% {new_biggr_id}")
        if (
            metabolite_db := session.scalars(
                select(Component)
                .filter(Component.bigg_id == new_biggr_id)
                .filter(
                    (Component.collection_id == None)
                    | (Component.collection_id == collection_db.id)
                )
                .limit(1)
            ).first()
        ) is not None:
            if (
                metabolite_db.charge != charge
                or not utils.are_explicit_formulae_equivalent(
                    metabolite_db.formula, formula
                )
            ):
                logging.warning(
                    f"Found component, but charge or formula did not match: {metabolite_db} (charge: {metabolite_db.charge}, formula: {metabolite_db.formula}) != (charge: {charge}, formula: {formula})"
                )
                metabolite_db = None

        if metabolite_db is None and formula is not None:
            _universal_component_db = session.scalars(
                select(UniversalComponent)
                .options(
                    subqueryload(
                        UniversalComponent.components.and_(
                            Component.collection_id == None
                        )
                    )
                    .subqueryload(Component.reference_mappings)
                    .joinedload(ComponentReferenceMapping.reference_compound)
                )
                .filter(UniversalComponent.bigg_id == new_universal_bigg_id)
                .limit(1)
            ).first()
            comp_formula = utils.Formula(formula)
            if (
                _universal_component_db is not None
                and _universal_component_db.allow_flexible_variants
            ):
                formula_delta_variant_mapping = {}
                for _comp_db in _universal_component_db.components:
                    if len(_comp_db.reference_mappings) != 1:
                        continue
                    # delta_formula = _comp_db.reference_mappings[
                    #     0
                    # ].reference_formula_delta
                    delta_formula = utils.Formula(_comp_db.formula) - utils.Formula(
                        _comp_db.reference_mappings[0].reference_compound.formula
                    )
                    normalized_delta = delta_formula - utils.Formula(
                        {
                            "H": (
                                _comp_db.charge
                                - int(
                                    _comp_db.reference_mappings[
                                        0
                                    ].reference_compound.charge
                                )
                            )
                        }
                    )
                    normalized_delta_str = str(normalized_delta)
                    if normalized_delta_str not in formula_delta_variant_mapping:
                        formula_delta_variant_mapping[normalized_delta_str] = (
                            _comp_db.variant
                        )
                    if _comp_db.charge != charge:
                        continue
                    if comp_formula == utils.Formula(_comp_db.formula):
                        metabolite_db = _comp_db
                        break
                if metabolite_db is None:
                    for _comp_db in _universal_component_db.components:
                        if len(_comp_db.reference_mappings) != 1:
                            continue
                        reference_db = _comp_db.reference_mappings[0].reference_compound
                        ref_formula = utils.Formula(reference_db.formula)
                        delta_formula = comp_formula - ref_formula

                        normalized_delta = delta_formula - utils.Formula(
                            {"H": (charge - int(reference_db.charge))}
                        )

                        if (
                            normalized_delta_str := str(normalized_delta)
                        ) in formula_delta_variant_mapping:
                            variant = formula_delta_variant_mapping[
                                normalized_delta_str
                            ]
                        else:
                            variant = max(formula_delta_variant_mapping.values()) + 1

                        full_bid = create_component_bigg_id(
                            new_universal_bigg_id, variant=variant, charge=charge
                        )
                        metabolite_db = metabolites._create_component(
                            session,
                            full_bid,
                            _universal_component_db,
                            formula=formula,
                            variant=variant,
                            charge=charge,
                            name=reference_db.name,
                        )
                        if metabolite_db is None:
                            raise ValueError()
                        component_ref_mapping_db = ComponentReferenceMapping(
                            universal_component=_universal_component_db,
                            reference_compound=reference_db,
                            reference_formula_delta=str(delta_formula),
                        )
                        metabolite_db.reference_mappings.append(
                            component_ref_mapping_db
                        )
                        session.add(metabolite_db)
                        session.flush()
                        break

        # if necessary, add the new metabolite, and keep track of the ID
        if metabolite_db is None:
            print("Creating a model-specific metabolite entry.")
            # make the new metabolite
            name = scrub_name(getattr(metabolite, "name", None))
            if name is None:
                name = ""
            universal_component_db, metabolite_db = (
                create_collection_specific_metabolite(
                    bigg_id=new_universal_bigg_id,
                    collection_db=collection_db,
                    charge=charge,
                    formula=formula,
                    name=name,
                    session=session,
                )
            )

        print(f"%% <> {metabolite_db} ({universal_component_db})")
        if metabolite_db is None or universal_component_db is None:
            raise ValueError("Error parsing metabolites")

        compartment_db = utils.get_object_by_bigg_id(
            session, compartment_bigg_id, Compartment
        )
        if compartment_db is None:
            raise ValueError(f"Could not find compartment {compartment_bigg_id}")

        universal_comp_comp_id = create_component_bigg_id(
            universal_component_db.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
        )
        comp_comp_id = create_component_bigg_id(
            universal_component_db.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
            charge=metabolite_db.charge,
            variant=metabolite_db.variant,
        )

        # if there is no universal compartmentalized component, add a new one
        if (
            universal_comp_component_db := utils.get_object_by_bigg_id(
                session, universal_comp_comp_id, UniversalCompartmentalizedComponent
            )
        ) is None:
            universal_comp_component_db = UniversalCompartmentalizedComponent(
                bigg_id=universal_comp_comp_id,
                universal_component=universal_component_db,
                compartment=compartment_db,
            )
            session.add(universal_comp_component_db)

        # if there is no compartmentalized component, add a new one
        if (
            comp_component_db := utils.get_object_by_bigg_id(
                session, comp_comp_id, CompartmentalizedComponent
            )
        ) is None:
            comp_component_db = CompartmentalizedComponent(
                bigg_id=comp_comp_id,
                component=metabolite_db,
                universal_compartmentalized_component=universal_comp_component_db,
                compartment=compartment_db,
            )
            session.add(comp_component_db)

        # Check if there is an existing entry that exactly matches
        model_comp_comp_db = session.scalars(
            select(ModelCompartmentalizedComponent)
            .filter(ModelCompartmentalizedComponent.model == model_db)
            .filter(
                ModelCompartmentalizedComponent.id_in_original_model == metabolite_id
            )
            .limit(1)
        ).first()
        if model_comp_comp_db is None:
            model_comp_comp_db_all = session.scalars(
                select(ModelCompartmentalizedComponent)
                .options(
                    contains_eager(
                        ModelCompartmentalizedComponent.compartmentalized_component
                    ).options(
                        contains_eager(
                            CompartmentalizedComponent.universal_compartmentalized_component
                        ),
                        joinedload(CompartmentalizedComponent.component).joinedload(
                            Component.universal_component
                        ),
                    )
                )
                .join(ModelCompartmentalizedComponent.compartmentalized_component)
                .join(CompartmentalizedComponent.universal_compartmentalized_component)
                .filter(
                    CompartmentalizedComponent.universal_compartmentalized_component
                    == universal_comp_component_db
                )
                .filter(ModelCompartmentalizedComponent.model == model_db)
            ).all()
            mcc_duplicates = 0
            for mcc_db in model_comp_comp_db_all:
                if not ":" in mcc_db.bigg_id:
                    mcc_db.bigg_id = create_component_bigg_id(
                        mcc_db.compartmentalized_component.component.universal_component.bigg_id,
                        compartment_bigg_id=compartment_db.bigg_id,
                        charge=mcc_db.compartmentalized_component.component.charge,
                        variant=mcc_db.compartmentalized_component.component.variant,
                    )
                if mcc_db.compartmentalized_component == comp_component_db:
                    mcc_duplicates += 1
            if mcc_duplicates == 0:
                if len(model_comp_comp_db_all) > 0:
                    model_comp_comp_id = create_component_bigg_id(
                        universal_component_db.bigg_id,
                        compartment_bigg_id=compartment_db.bigg_id,
                        charge=metabolite_db.charge,
                        variant=metabolite_db.variant,
                    )
                else:
                    model_comp_comp_id = universal_comp_component_db.bigg_id
                model_comp_comp_db = ModelCompartmentalizedComponent(
                    bigg_id=model_comp_comp_id,
                    model=model_db,
                    compartmentalized_component=comp_component_db,
                    id_in_original_model=metabolite_id,
                )
                session.add(model_comp_comp_db)
            else:
                print(
                    f"Duplicate metabolite in model, creating a model-specific metabolite entry. (metabolite_id: {metabolite_id}, {new_universal_bigg_id})"
                )
                metabolite_db = None
                base_universal_bigg_id = new_universal_bigg_id

                # Try to use the same duplicate component as in the same model, but with different compartment
                component_bigg_id_pattern = create_component_bigg_id(
                    f"{base_universal_bigg_id}_DUP%",
                    charge=charge,
                    collection_bigg_id=collection_db.bigg_id,
                )
                model_compartmentalized_components_db = session.scalars(
                    select(ModelCompartmentalizedComponent)
                    .options(
                        contains_eager(
                            ModelCompartmentalizedComponent.compartmentalized_component
                        ).contains_eager(CompartmentalizedComponent.component)
                    )
                    .join(ModelCompartmentalizedComponent.compartmentalized_component)
                    .join(CompartmentalizedComponent.component)
                    .filter(Component.bigg_id.like(component_bigg_id_pattern))
                    .filter(ModelCompartmentalizedComponent.model_id == model_db.id)
                ).all()
                for model_comp_comp_db in model_compartmentalized_components_db:
                    if (
                        model_comp_comp_db.compartmentalized_component.compartment_id
                        == compartment_db.id
                    ):
                        continue
                    other_bigg_id, other_compartment_id = parse.split_compartment(
                        model_comp_comp_db.id_in_original_model
                    )
                    if other_bigg_id == base_universal_bigg_id:
                        metabolite_db = (
                            model_comp_comp_db.compartmentalized_component.component
                        )
                        break

                if metabolite_db is None:
                    # Try to use the same duplicate component as in other models in the collection
                    component_bigg_id_pattern = create_component_bigg_id(
                        f"{base_universal_bigg_id}_DUP%",
                        charge=charge,
                        collection_bigg_id=collection_db.bigg_id,
                    )

                    met_same_model_subq = (
                        select(Component.id.label("comp_id"))
                        .filter(Component.bigg_id.like(component_bigg_id_pattern))
                        .join(Component.compartmentalized_components)
                        .join(
                            CompartmentalizedComponent.model_compartmentalized_components
                        )
                        .filter(
                            ModelCompartmentalizedComponent.id_in_original_model
                            == metabolite_id
                        )
                        .filter(ModelCompartmentalizedComponent.model_id == model_db.id)
                        .group_by(Component.id)
                        .having(func.count(ModelCompartmentalizedComponent.id) > 0)
                    ).subquery()
                    met_query = (
                        select(Component.id)
                        .filter(Component.bigg_id.like(component_bigg_id_pattern))
                        .join(Component.compartmentalized_components)
                        .join(
                            CompartmentalizedComponent.model_compartmentalized_components
                        )
                        .filter(
                            ModelCompartmentalizedComponent.id_in_original_model
                            == metabolite_id
                        )
                        .filter(ModelCompartmentalizedComponent.model_id != model_db.id)
                        .group_by(Component.id)
                        .having(func.count(ModelCompartmentalizedComponent.id) > 0)
                        .join(
                            met_same_model_subq,
                            met_same_model_subq.c.comp_id == Component.id,
                            isouter=True,
                        )
                        .filter(met_same_model_subq.c.comp_id == None)
                        .limit(1)
                    )
                    metabolite_db_id = session.scalars(met_query).first()
                    if metabolite_db_id is not None:
                        metabolite_db = session.get(Component, metabolite_db_id)

                if metabolite_db is None:
                    print("Creating new duplicate")
                    new_universal_bigg_id = f"{base_universal_bigg_id}_DUP1"
                    for i in range(1, 100):
                        new_universal_bigg_id = f"{base_universal_bigg_id}_DUP{i}"
                        new_bigg_id = create_component_bigg_id(
                            new_universal_bigg_id,
                            collection_bigg_id=collection_db.bigg_id,
                        )
                        temp_uni_comp_db = utils.get_object_by_bigg_id(
                            session, new_bigg_id, UniversalComponent
                        )
                        if temp_uni_comp_db is None:
                            break

                    # make the new metabolite
                    universal_component_db, metabolite_db = (
                        create_collection_specific_metabolite(
                            bigg_id=new_universal_bigg_id,
                            collection_db=collection_db,
                            charge=charge,
                            formula=formula,
                            name="",
                            session=session,
                        )
                    )
                    if universal_component_db is None or metabolite_db is None:
                        raise ValueError()
                    comp_comp_id = create_component_bigg_id(
                        universal_component_db.bigg_id,
                        compartment_bigg_id=compartment_db.bigg_id,
                        charge=metabolite_db.charge,
                        variant=metabolite_db.variant,
                    )
                    universal_comp_comp_id = create_component_bigg_id(
                        universal_component_db.bigg_id,
                        compartment_bigg_id=compartment_db.bigg_id,
                    )
                    universal_comp_component_db = UniversalCompartmentalizedComponent(
                        bigg_id=universal_comp_comp_id,
                        universal_component=universal_component_db,
                        compartment=compartment_db,
                    )
                    session.add(universal_comp_component_db)
                    comp_component_db = CompartmentalizedComponent(
                        bigg_id=comp_comp_id,
                        component=metabolite_db,
                        universal_compartmentalized_component=universal_comp_component_db,
                        compartment=compartment_db,
                    )
                    session.add(comp_component_db)
                else:
                    universal_component_db = session.get(
                        UniversalComponent, metabolite_db.universal_component_id
                    )
                    comp_component_db = get_or_create_compartmentalized_component(
                        session, metabolite_db, compartment_db
                    )

                model_comp_comp_id = create_component_bigg_id(
                    universal_component_db.bigg_id,
                    compartment_bigg_id=compartment_db.bigg_id,
                    # charge=metabolite_db.charge,
                )
                model_comp_comp_db = ModelCompartmentalizedComponent(
                    bigg_id=model_comp_comp_id,
                    model=model_db,
                    compartmentalized_component=comp_component_db,
                    id_in_original_model=metabolite_id,
                )
                session.add(model_comp_comp_db)
    session.commit()


def compartmentalized_id_to_universal_compartmentalized_id(comp_comp_id):
    try:
        universal_comp_comp_id, _charge = comp_comp_id.rsplit(":", maxsplit=1)
        return universal_comp_comp_id
    except:
        return comp_comp_id


@timing
def load_reactions(
    session,
    model_db_id,
    model,
):
    """Load the reactions and stoichiometries into the model.

    TODO if the reaction is already loaded, we need to check the stoichometry
    has. If that doesn't match, then add a new reaction with an incremented ID
    (e.g. ACALD_1)

    Arguments
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    old_reaction_ids: A dictionary where keys are new IDs and values are old IDs
    for reactions.

    comp_comp_db_ids: A dictionary where keys are the original compartmentalized
    metabolite ids and the values are the database IDs for the compartmentalized
    components.

    final_metabolite_ids: A new dictionary where keys are original
    compartmentalized metabolite IDs from the model and values are the new
    compartmentalized metabolite IDs.

    Returns
    -------

    A dictionary with keys for reaction BiGG IDs in the model and values for the
    associated ModelReaction.id in the database.

    """

    model_db = session.get(Model, model_db_id)
    collection_db = session.get(ModelCollection, model_db.collection_id)
    model_db_rxn_ids = {}
    for reaction in model.reactions:
        # Drop duplicates label
        # model_specific_id = reaction.id
        # reaction_id = model_specific_id.removeprefix(f"__{model_db.bigg_id}__")
        # reaction_id = parse.remove_duplicate_tag(reaction_id)
        reaction_id = reaction.id
        clean_reaction_id = parse.remove_duplicate_tag(reaction_id)
        clean_reaction_id, model_reaction_copy_number = parse.split_id_and_copy_tag(
            clean_reaction_id
        )

        reaction_name = rn if (rn := reaction.name.strip()) != "" else None

        # TODO: Make sure to keep copy number equal model and biggr

        met_mappings = {
            m.id: (
                session.scalars(
                    select(ModelCompartmentalizedComponent)
                    .options(
                        joinedload(
                            ModelCompartmentalizedComponent.compartmentalized_component
                        ).options(
                            joinedload(CompartmentalizedComponent.component),
                            joinedload(
                                CompartmentalizedComponent.universal_compartmentalized_component
                            ).joinedload(
                                UniversalCompartmentalizedComponent.universal_component
                            ),
                        )
                    )
                    .filter(
                        (ModelCompartmentalizedComponent.model == model_db)
                        & (
                            ModelCompartmentalizedComponent.id_in_original_model == m.id
                        ),
                    )
                    .limit(1)
                ).first()
            )
            for m in reaction.metabolites.keys()
        }
        if any(x is None for x in met_mappings.values()):
            raise ValueError()

        participants = [
            ReactionParticipantInfo(
                coefficient=coeff,
                compartmentalized_component=met_mappings[
                    m.id
                ].compartmentalized_component,
            )
            for m, coeff in reaction.metabolites.items()
        ]
        # print(f"{reaction_id}: {participants}")
        reaction_hash = Reaction.generate_hash(participants)
        # print(f"reaction hash 0: {reaction_hash}")

        # Get the reaction
        reaction_db = session.scalars(
            select(Reaction)
            .filter(
                (Reaction.hash == reaction_hash)
                & (
                    (Reaction.collection_id == None)
                    | (Reaction.collection_id == collection_db.id)
                )
            )
            .join(Reaction.universal_reaction)
            .limit(1)
        ).first()

        if reaction_db is None:
            universal_participants = get_universal_reaction_participants_list(
                session, participants
            )
            universal_reaction_hash = UniversalReaction.generate_hash(
                universal_participants
            )
            print(f"universal hash 0: {universal_reaction_hash}")
            universal_reaction_db = session.scalars(
                select(UniversalReaction)
                .filter(UniversalReaction.hash == universal_reaction_hash)
                .filter(
                    (UniversalReaction.collection_id == None)
                    | (UniversalReaction.collection_id == collection_db.id)
                )
                .limit(1)
            ).first()
            if universal_reaction_db is not None:
                logging.warning(
                    f"\t UniversalReaction {reaction_id}: {universal_reaction_db}"
                )
                print("! Creating new reaction variant.")
                get_or_create_reaction_for_universal_reaction(
                    session, universal_reaction_db, participants
                )
                session.flush()
                print(f"reaction hash 3: {reaction_hash}")
                reaction_db = session.scalars(
                    select(Reaction)
                    .filter(
                        (Reaction.hash == reaction_hash)
                        & (
                            (Reaction.collection_id == None)
                            | (Reaction.collection_id == collection_db.id)
                        )
                    )
                    .join(Reaction.universal_reaction)
                    .limit(1)
                ).first()
            else:
                # Check for exchange reactions
                reaction_collection_id = None
                is_exchange = False
                if (
                    len(participants) == 1
                    and abs(float(participants[0].coefficient)) == 1
                ):
                    universal_compartmentalized_component_bigg_id = (
                        universal_participants[
                            0
                        ].universal_compartmentalized_component.bigg_id
                    )
                    if (
                        clean_reaction_id
                        != f"EX_{universal_compartmentalized_component_bigg_id}"
                    ):
                        print(f"Wrong name for exchange reaction: {clean_reaction_id}")
                        clean_reaction_id = (
                            f"EX_{universal_compartmentalized_component_bigg_id}"
                        )
                    is_exchange = True
                    if (
                        collection_id := universal_participants[
                            0
                        ].universal_compartmentalized_component.universal_component.collection_id
                    ) is not None:
                        reaction_collection_id = collection_id

                if not is_exchange:
                    print("! Creating collection-specific reaction.")
                    clean_reaction_id = (
                        f"__{collection_db.bigg_id}__{clean_reaction_id}"
                    )
                    reaction_collection_id = collection_db.id
                else:
                    print("! Creating exchange reaction.")

                universal_reaction_db = get_or_create_universal_reaction(
                    session,
                    clean_reaction_id,
                    universal_participants,
                    reaction_collection_id=reaction_collection_id,
                    exchange_reaction=is_exchange,
                    reaction_name=reaction_name,
                )
                session.flush()
                reaction_db = get_or_create_reaction_for_universal_reaction(
                    session, universal_reaction_db, participants
                )

                session.flush()
                print(f"Collection id: {reaction_collection_id} ({collection_db.id})")
                if universal_reaction_db is not None:
                    print(f"Created universal reaction")
                    print(
                        f"Hash of created universal_reaction: '{universal_reaction_db.hash}'"
                    )
                    print(
                        f"Expected hash of universal_reaction: '{universal_reaction_hash}'"
                    )
                else:
                    print("Did not create universal_reaction")

                if reaction_db is not None:
                    print(f"Created reaction")
                    print(f"Hash of created reaction: '{reaction_db.hash}'")
                    print(f"Expected hash of reaction: '{reaction_hash}'")
                else:
                    print("Did not create reaction")

                reaction_db = session.scalars(
                    select(Reaction)
                    .filter(
                        (Reaction.hash == reaction_hash)
                        & (
                            (Reaction.collection_id == None)
                            | (Reaction.collection_id == collection_db.id)
                        )
                    )
                    .join(Reaction.universal_reaction)
                    .limit(1)
                ).first()
                print(f"Retreived reaction: {reaction_db}")

        if reaction_db is None:
            print("ERROR: Reaction was not correctly created.")

        if reaction_db.universal_reaction.name is None and reaction_name is not None:
            reaction_db.universal_reaction.name = reaction_name

        logging.warning(f"Reaction {clean_reaction_id}: {reaction_db}")

        reaction_matrix_db = session.scalars(
            select(ReactionMatrix)
            .join(ReactionMatrix.universal_reaction_matrix)
            .join(ReactionMatrix.compartmentalized_component)
            .filter(ReactionMatrix.reaction == reaction_db)
        ).all()
        model_reaction_is_reversed = False
        for rm_db in reaction_matrix_db:
            # TODO: May fail in some cases where metabolites occur at both sides
            # for m, coeff in reaction.metabolites.items():
            for p in participants:
                comp_comp_id = p.compartmentalized_component.bigg_id
                coeff = p.coefficient
                if comp_comp_id == rm_db.compartmentalized_component.bigg_id:
                    if rm_db.universal_reaction_matrix.coefficient == coeff:
                        pass
                    elif rm_db.universal_reaction_matrix.coefficient == -1 * coeff:
                        model_reaction_is_reversed = True
                        break
                    else:
                        logging.error(
                            f"Coefficients do not match for {clean_reaction_id}, {comp_comp_id}"
                        )
            break

        # If the reaction is reversed, then switch upper and lower bound
        lower_bound = (
            -reaction.upper_bound
            if model_reaction_is_reversed
            else reaction.lower_bound
        )
        upper_bound = (
            -reaction.lower_bound
            if model_reaction_is_reversed
            else reaction.upper_bound
        )

        # TODO: Flip reaction in model if necessary

        # subsystem
        subsystem = check_none(reaction.subsystem.strip())

        copy_number = (
            session.scalars(
                select(func.count(ModelReaction.id))
                .join(ModelReaction.reaction)
                .filter(
                    Reaction.universal_reaction_id == reaction_db.universal_reaction_id
                )
                .filter(ModelReaction.model == model_db)
            ).first()
            + 1
        )

        # make a new reaction
        model_reaction_id = ModelReaction.create_id(
            model_db.bigg_id, reaction_db.universal_reaction.bigg_id, copy_number
        )
        model_reaction_db = ModelReaction(
            bigg_id=model_reaction_id,
            model=model_db,
            reaction=reaction_db,
            reversed=model_reaction_is_reversed,
            gene_reaction_rule=reaction.gene_reaction_rule,
            original_gene_reaction_rule=reaction.gene_reaction_rule,
            upper_bound=upper_bound,
            lower_bound=lower_bound,
            objective_coefficient=reaction.objective_coefficient,
            copy_number=copy_number,
            subsystem=subsystem,
            id_in_original_model=reaction_id,
        )
        session.add(model_reaction_db)
        session.flush()
    session.commit()
    return model_db_rxn_ids


# find gene functions
def _match_gene_by_fns(fn_list, session, gene_id, chromosome_ids):
    """Go through each funciton and look for a match."""
    for fn in fn_list:
        match, is_alternative_transcript = fn(session, gene_id, chromosome_ids)
        if len(match) > 0:
            if len(match) > 1:
                logging.warning(
                    "Multiple matches for gene {} with function {}. Using the first match.".format(
                        gene_id, fn.__name__
                    )
                )
            return match[0], is_alternative_transcript
    return None, False


def _by_bigg_id(session, gene_id, chromosome_ids):
    # look for a matching model gene
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.bigg_id) == func.lower(gene_id))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _by_name(session, gene_id, chromosome_ids):
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.name) == func.lower(gene_id))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _by_synonym(session, gene_id, chromosome_ids):
    gene_db = session.scalars(
        select(Gene)
        .join(Synonym, Synonym.ome_id == Gene.id)
        .filter(Gene.chromosome_id.in_(chromosome_ids))
        .filter(func.lower(Synonym.synonym) == func.lower(gene_id))
    ).all()
    return gene_db, False


def _by_alternative_transcript(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Gene.bigg_id) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_alternative_transcript_name(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Gene.name) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_alternative_transcript_synonym(session, gene_id, chromosome_ids):
    """Function to check for the alternative transcript match."""
    check = re.match(r"(.*)_AT[0-9]{1,2}$", gene_id)
    if not check:
        gene_db = []
    else:
        # find the old gene
        gene_db = session.scalars(
            select(Gene)
            .join(Synonym, Synonym.ome_id == Gene.id)
            .filter(Gene.chromosome_id.in_(chromosome_ids))
            .filter(func.lower(Synonym.synonym) == func.lower(check.group(1)))
            .filter(Gene.alternative_transcript_of.is_(None))
        ).all()
    return gene_db, True


def _by_bigg_id_no_underscore(session, gene_id, chromosome_ids):
    """Matches for T maritima genes"""
    # look for a matching model gene
    gene_db = session.scalars(
        select(Gene)
        .filter(func.lower(Gene.bigg_id) == func.lower(gene_id.replace("_", "")))
        .filter(Gene.chromosome_id.in_(chromosome_ids))
    ).all()
    return gene_db, False


def _replace_gene_str(rule, old_gene, new_gene):
    return re.sub(r"\b" + old_gene + r"\b", new_gene, rule)


@timing
def load_genes(session, model_db_id, model):
    """Load the genes for this model.

    Arguments:
    ---------

    session: An SQLAlchemy session.

    model_db_id: The database ID for the model.

    model: The COBRApy model.

    model_db_rxn_ids: A dictionary with keys for reactions in the model and
    values for the associated ModelReaction.id in the database.

    old_gene_ids: A dictionary where keys are new IDs and values are old IDs for
    genes.

    """

    # find the model in the db
    model_db = session.get(Model, model_db_id)

    # find the chromosomes in the db
    chromosome_ids = session.scalars(
        select(Chromosome.id).filter(Chromosome.genome_id == model_db.genome_id)
    ).all()
    chromosome_ids = [c for c in chromosome_ids]
    if len(chromosome_ids) == 0:
        logging.warning("No chromosomes for model %s" % model_db.bigg_id)

    context = {}

    # keep track of the gene-reaction associations
    gene_bigg_id_to_model_reaction_db_ids = defaultdict(set)
    for reaction in model.reactions:
        # find the ModelReaction that corresponds to this particular reaction in
        # the model
        # if reaction.id not in model_db_rxn_ids:
        #     continue
        # universal_reaction_id, copy_number = ModelReaction.interpret_id(reaction.id)
        # model_reaction_db = session.scalars(
        #     select(ModelReaction)
        #     .join(ModelReaction.reaction)
        #     .join(Reaction.universal_reaction)
        #     .filter(
        #         (UniversalReaction.bigg_id == universal_reaction_id)
        #         & (ModelReaction.copy_number == copy_number)
        #         & (ModelReaction.model == model_db)
        #     )
        #     .limit(1)
        # ).first()
        model_reaction_db = session.scalars(
            select(ModelReaction)
            .filter(ModelReaction.model == model_db)
            .filter(ModelReaction.id_in_original_model == reaction.id)
            .limit(1)
        ).first()

        if model_reaction_db is None:
            logging.error(
                "Could not find ModelReaction {} in model {}. Cannot load GeneReactionMatrix entries".format(
                    reaction.id, model_db.bigg_id
                )
            )
            continue
        for gene in reaction.genes:
            gene_bigg_id_to_model_reaction_db_ids[gene.id].add(model_reaction_db.id)

    # load the genes
    for gene in model.genes:
        if len(chromosome_ids) == 0:
            gene_db = None
            is_alternative_transcript = False
        else:
            # find a matching gene
            fns = [
                _by_bigg_id,
                _by_name,
                # _by_synonym,
                # _by_alternative_transcript,
                # _by_alternative_transcript_name,
                # _by_alternative_transcript_synonym,
                _by_bigg_id_no_underscore,
            ]
            gene_db, is_alternative_transcript = _match_gene_by_fns(
                fns, session, gene.id, chromosome_ids
            )

        dup_gene_alternative_transcript = (not gene_db) and is_alternative_transcript
        if not gene_db:
            # add
            if len(chromosome_ids) > 0:
                logging.warning(
                    "Gene not in genbank file: {} from model {}".format(
                        gene.id, model.id
                    )
                )
            gene_db = session.scalars(
                select(Gene).filter(Gene.bigg_id == gene.id).limit(1)
            ).first()
            if gene_db is None:
                gene_db = Gene(
                    bigg_id=gene.id,
                    # name is optional in cobra 0.4b2. This will probably change back.
                    name=scrub_name(getattr(gene, "name", None)),
                    mapped_to_genbank=False,
                )
                session.add(gene_db)
            # session.commit()

        context[gene.id] = {
            "gene": gene,
            "gene_db": gene_db,
            "dup_gene_alternative_transcript": dup_gene_alternative_transcript,
        }
        if dup_gene_alternative_transcript:
            # duplicate gene for the alternative transcript
            old_gene_db = gene_db
            ome_gene = {}
            ome_gene["bigg_id"] = gene.bigg_id
            ome_gene["name"] = old_gene_db.name
            ome_gene["leftpos"] = old_gene_db.leftpos
            ome_gene["rightpos"] = old_gene_db.rightpos
            ome_gene["chromosome_id"] = old_gene_db.chromosome_id
            ome_gene["strand"] = old_gene_db.strand
            ome_gene["mapped_to_genbank"] = True
            ome_gene["alternative_transcript_of"] = old_gene_db.id
            gene_db = Gene(**ome_gene)
            session.add(gene_db)
            # session.commit()
            context[gene.id]["old_gene_db"] = old_gene_db
    session.commit()

    for gene_id, ctx in context.items():
        # add the deprecated id if necessary
        gene_db, old_gene_db, dup_gene_alternative_transcript = (
            ctx["gene_db"],
            ctx.get("old_gene_db"),
            ctx["dup_gene_alternative_transcript"],
        )

        # if dup_gene_alternative_transcript:
        #     # duplicate all the synonyms
        #     synonyms_db = (
        #         session.query(Synonym).filter(Synonym.ome_id == old_gene_db.id).all()
        #     )
        #     for syn_db in synonyms_db:
        #         # add a new synonym
        #         ome_synonym = {}
        #         ome_synonym["type"] = syn_db.type
        #         ome_synonym["ome_id"] = gene_db.id
        #         ome_synonym["synonym"] = syn_db.synonym
        #         ome_synonym["data_source_id"] = syn_db.data_source_id
        #         synonym_object = Synonym(**ome_synonym)
        #         session.add(synonym_object)
        #
        # add model gene
        model_gene_db = session.scalars(
            select(ModelGene)
            .filter(ModelGene.gene == gene_db)
            .filter(ModelGene.model == model_db)
            .limit(1)
        ).first()
        if model_gene_db is None:
            model_gene_db = ModelGene(gene_id=gene_db.id, model=model_db)
            session.add(model_gene_db)
            # session.commit()
        ctx["model_gene_db"] = model_gene_db

    session.commit()
    for gene_id, ctx in context.items():
        gene_db, model_gene_db = ctx["gene_db"], ctx["model_gene_db"]
        # add old gene synonym
        old_bigg_synonyms = {}
        # for old_bigg_id in old_gene_ids[gene_id]:
        #     synonym_db = (
        #         session.query(Synonym)
        #         .filter(Synonym.type == "gene")
        #         .filter(Synonym.ome_id == gene_db.id)
        #         .filter(Synonym.synonym == old_bigg_id)
        #         .filter(Synonym.data_source_id == data_source_id)
        #         .first()
        #     )
        #     if synonym_db is None:
        #         synonym_db = Synonym(
        #             type="gene",
        #             ome_id=gene_db.id,
        #             synonym=old_bigg_id,
        #             data_source_id=data_source_id,
        #         )
        #         session.add(synonym_db)
        #         # session.commit()
        #     old_bigg_synonyms[old_bigg_id] = synonym_db
        ctx["old_bigg_synonyms"] = old_bigg_synonyms
    session.commit()
    for gene_id, ctx in context.items():
        gene_db, model_gene_db, old_bigg_synonyms = (
            ctx["gene_db"],
            ctx["model_gene_db"],
            ctx["old_bigg_synonyms"],
        )
        # for old_bigg_id, synonym_db in old_bigg_synonyms.items():
        #     # add OldIDSynonym
        #     old_id_db = (
        #         session.query(OldIDSynonym)
        #         .filter(OldIDSynonym.type == "model_gene")
        #         .filter(OldIDSynonym.ome_id == model_gene_db.id)
        #         .filter(OldIDSynonym.synonym_id == synonym_db.id)
        #         .first()
        #     )
        #     if old_id_db is None:
        #         old_id_db = OldIDSynonym(
        #             type="model_gene", ome_id=model_gene_db.id, synonym_id=synonym_db.id
        #         )
        #         session.add(old_id_db)
        #         # session.commit()

        # find model reaction
        try:
            model_reaction_db_ids = gene_bigg_id_to_model_reaction_db_ids[gene_id]
        except KeyError:
            # error message above
            continue

        for mr_db_id in model_reaction_db_ids:
            # add to the GeneReactionMatrix, if not already present
            found_gene_reaction_row = (
                session.scalars(
                    select(func.count(GeneReactionMatrix.id))
                    .filter(GeneReactionMatrix.model_gene == model_gene_db)
                    .filter(GeneReactionMatrix.model_reaction_id == mr_db_id)
                ).first()
                > 0
            )
            if not found_gene_reaction_row:
                new_object = GeneReactionMatrix(
                    model_gene=model_gene_db, model_reaction_id=mr_db_id
                )
                session.add(new_object)

            # update the gene_reaction_rule if the gene id has changed
            if gene_id != gene_db.bigg_id:
                mr = session.get(ModelReaction, mr_db_id)
                new_rule = _replace_gene_str(
                    mr.gene_reaction_rule, gene_id, gene_db.bigg_id
                )
                session.execute(
                    update(ModelReaction),
                    [{"id": mr_db_id, "gene_reaction_rule": new_rule}],
                )
    session.commit()


@timing
def load_model_count(session: Session, model_db_id: int):
    metabolite_count = session.scalars(
        select(func.count(ModelCompartmentalizedComponent.id)).filter(
            ModelCompartmentalizedComponent.model_id == model_db_id
        )
    ).first()
    reaction_count = session.scalars(
        select(func.count(ModelReaction.id)).filter(
            ModelReaction.model_id == model_db_id
        )
    ).first()
    gene_count = session.scalars(
        select(func.count(ModelGene.id)).filter(ModelGene.model_id == model_db_id)
    ).first()
    mc = ModelCount(
        model_id=model_db_id,
        gene_count=gene_count,
        metabolite_count=metabolite_count,
        reaction_count=reaction_count,
    )
    session.add(mc)
