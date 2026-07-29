import logging
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import joinedload
from cobradb.models import (
    Annotation,
    AnnotationLink,
    AnnotationProperty,
    CompartmentalizedComponent,
    Component,
    ComponentAnnotationMapping,
    ComponentIDMapping,
    ComponentReferenceMapping,
    InChI,
    Reaction,
    ReactionAnnotationMapping,
    ReferenceCompound,
    UniversalComponent,
    UniversalReaction,
)
from cobradb.data_sources import get_data_source_id
from modelseedpy.biochem import from_local
import subprocess
import pandas as pd

SEED_METABOLITE_PROPERTIES = [
    "name",
    "mass",
    "is_core",
    "is_obsolete",
    "is_cofactor",
    "delta_g",
    "delta_g_err",
    "pka",
    "pkb",
    "is_abstract",
    "smiles",
]
SEED_METABOLITE_ALIASES = {
    "seed.compound": "seed.compound",
    "KEGG": "kegg.compound",
    "MetaCyc": "metacyc.compound",
    "metanetx.chemical": "metanetx.chemical",
}

SEED_REACTION_PROPERTIES = [
    "name",
    "is_obsolete",
    "deltag",
    "deltagerr",
    "is_abstract",
]
SEED_REACTION_ALIASES = {
    "seed.reaction": "seed.reaction",
    "KEGG": "kegg.reaction",
    "MetaCyc": "metacyc.reaction",
    "metanetx.reaction": "metanetx.reaction",
}


def download_or_update(directory):
    directory = Path(directory)
    git_directory = directory / "ModelSEEDDatabase"

    if git_directory.exists():
        process = subprocess.run(
            ["git", "pull"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=git_directory,
        )
    else:
        process = subprocess.run(
            [
                "git",
                "clone",
                # "-b",
                # "dev",
                "https://github.com/ModelSEED/ModelSEEDDatabase.git",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            cwd=directory,
        )
    if process.returncode != 0:
        print("Could not download or update the ModelSEED Database.")


def load_modelseed_database(directory):
    directory = Path(directory)
    db = from_local(str(directory / "ModelSEEDDatabase"))
    return db


def push_model_seed_metabolites(modelseed_db, session):
    for cpd_id in modelseed_db.compounds:
        cpd = modelseed_db.get_seed_compound(cpd_id)

        aliases = cpd.aliases
        cpd_old_bigg_ids = aliases.get("BiGG", [])
        if isinstance(cpd_old_bigg_ids, str):
            cpd_old_bigg_ids = [cpd_old_bigg_ids]

        inchi_key = cpd.inchikey
        inchi_key = None if pd.isna(inchi_key) or len(inchi_key) == 0 else inchi_key

        charge = cpd.data.get("charge")
        charge = None if pd.isna(charge) else int(charge)

        if charge is None:
            charge = 0

        cpd_bigg_ids = set()
        bigg_component_db = set()
        for cpd_bigg_id in cpd_old_bigg_ids:
            bc_db = session.scalars(
                select(Component)
                .join(Component.universal_component)
                .join(UniversalComponent.old_bigg_ids)
                .filter(ComponentIDMapping.old_bigg_id == cpd_bigg_id)
                .filter(Component.charge == charge)
                .limit(1)
            ).first()
            if bc_db is None:
                bc_db = session.scalars(
                    select(Component)
                    .join(Component.universal_component)
                    .filter(UniversalComponent.bigg_id == cpd_bigg_id)
                    .filter(Component.charge == charge)
                    .limit(1)
                ).first()
            if bc_db is not None:
                bigg_component_db.add(bc_db)
                cpd_bigg_ids.add(bc_db.bigg_id)
            else:
                cpd_bigg_ids.add(cpd_bigg_id)
        bigg_component_db = {x.bigg_id: x for x in bigg_component_db}

        inchi_key_component_db = None
        if inchi_key is not None:
            inchi_key_component_db = session.execute(
                select(Component, InChI)
                .join(Component.reference_mappings)
                .join(ComponentReferenceMapping.reference_compound)
                .join(ReferenceCompound.inchi)
                .filter(Component.charge == charge)
                .filter(InChI.key == inchi_key)
            ).all()

            inchi = cpd.inchi
            inchi = None if pd.isna(inchi) or len(inchi) == 0 else inchi

            # Select proper match based on full InChI, if possible
            if inchi is None:
                if len(inchi_key_component_db) == 1:
                    inchi_key_component_db = inchi_key_component_db[0][0]
                else:
                    inchi_key_component_db = None
            else:
                for cdb, ikdb in inchi_key_component_db:
                    if ikdb.to_string() == inchi:
                        inchi_key_component_db = cdb
                        break
                else:
                    inchi_key_component_db = None
        inchi_key_component_db = (
            {}
            if inchi_key_component_db is None
            else {inchi_key_component_db.bigg_id: inchi_key_component_db}
        )

        if not bigg_component_db and not inchi_key_component_db:
            continue

        cpd_names = cpd.data.get("aliases", "")
        cpd_names = cpd_names.split("|")
        cpd_names = [x[5:] for x in cpd_names if x.startswith("Name:")]
        cpd_names = cpd_names[0] if len(cpd_names) > 0 else ""
        cpd_names = [x.strip() for x in cpd_names.split(";")]

        annotation_bigg_id = f"seed.metabolite:{cpd_id}"
        annotation_db = session.scalars(
            select(Annotation).filter(Annotation.bigg_id == annotation_bigg_id).limit(1)
        ).first()
        if annotation_db is None:
            annotation_db = Annotation(
                bigg_id=annotation_bigg_id,
                type="seed",
                default_data_source_id=get_data_source_id("seed.compound", session),
            )
            session.add(annotation_db)

            for property in SEED_METABOLITE_PROPERTIES:
                prop_val = getattr(cpd, property, None)
                if prop_val is None:
                    prop_val = cpd.data.get(property)
                if isinstance(prop_val, str):
                    prop_val = prop_val.strip()
                    if len(prop_val) == 0:
                        prop_val = None
                if prop_val is None:
                    continue
                if pd.isna(prop_val):
                    continue
                if isinstance(prop_val, int) and property.startswith("is_"):
                    prop_val = bool(prop_val)
                if property == "is_obsolete":
                    annotation_db.is_obsolete = bool(prop_val)
                    continue
                prop_db = AnnotationProperty(key=property)
                prop_db.value = prop_val
                annotation_db.properties.append(prop_db)

            for name in cpd_names:
                if name == getattr(cpd, "name", None):
                    continue
                prop_db = AnnotationProperty(key="name")
                prop_db.value = name
                annotation_db.properties.append(prop_db)

            for property, namespace in SEED_METABOLITE_ALIASES.items():
                if namespace == "seed.compound":
                    prop_val = cpd_id
                else:
                    prop_val = aliases.get(property)
                    if isinstance(prop_val, str):
                        prop_val = prop_val.strip()
                        if len(prop_val) == 0:
                            prop_val = None
                    if prop_val is None:
                        continue
                if pd.isna(prop_val):
                    continue
                if isinstance(prop_val, str):
                    prop_val = {prop_val}

                for val in prop_val:
                    alias_db = AnnotationLink(
                        identifier=val,
                        data_source_id=get_data_source_id(namespace, session),
                    )
                    annotation_db.links.append(alias_db)

        for k in set(bigg_component_db.keys()) | set(inchi_key_component_db):
            bc_db = bigg_component_db.get(k)
            ikc_db = inchi_key_component_db.get(k)

            bigg_id_match = True if bc_db is not None else None
            if bigg_id_match is None and k in cpd_bigg_ids:
                bigg_id_match = False
            inchi_match = True if ikc_db is not None else None
            if inchi_match is None and inchi_key:
                inchi_match = False

            annotation_mapping = ComponentAnnotationMapping(
                component=bc_db if bigg_id_match else ikc_db,
                bigg_id_match=bigg_id_match,
                inchi_match=inchi_match,
            )
            found_match = False
            for cmap in annotation_db.component_mappings:
                if (
                    (cmap.component_id == annotation_mapping.component_id)
                    and (cmap.bigg_id_match == annotation_mapping.bigg_id_match)
                    and (cmap.inchi_match == annotation_mapping.inchi_match)
                ):
                    found_match = True
                    break
            if not found_match:
                annotation_db.component_mappings.append(annotation_mapping)
    session.commit()


def push_model_seed_reactions(modelseed_db, session):
    seed_compound_data_source_id = get_data_source_id("seed.compound", session)
    for rxn_id in modelseed_db.reactions:
        rxn = modelseed_db.get_seed_reaction(rxn_id)

        aliases = rxn.aliases
        rxn_bigg_ids = aliases.get("BiGG", [])
        if isinstance(rxn_bigg_ids, str):
            rxn_bigg_ids = [rxn_bigg_ids]

        bigg_reaction_db = set()
        for rxn_bigg_id in rxn_bigg_ids:
            # bigg_reaction_db = session.scalars(
            #     select(Reaction)
            #     .join(Reaction.universal_reaction)
            #     .join(UniversalComponent.old_bigg_ids)
            #     .filter(ComponentIDMapping.old_bigg_id == cpd_bigg_id)
            #     .filter(Component.charge == charge)
            #     .limit(1)
            # ).first()
            # if bigg_component_db is None:
            br_db = session.scalars(
                select(Reaction)
                .join(Reaction.universal_reaction)
                .filter(UniversalReaction.bigg_id == rxn_bigg_id)
                .limit(1)
            ).first()
            if br_db is not None:
                bigg_reaction_db.add(br_db)
        bigg_reaction_db = {x.bigg_id: x for x in bigg_reaction_db}

        pattern_matched_reaction_db = []
        front = [([], {})]
        try:  # modelseedpy will throw an error in some cases
            cstoichiometry = rxn.cstoichiometry
        except:
            cstoichiometry = {}
        if cstoichiometry:
            for (cpd, seed_comp), coefficient in cstoichiometry.items():
                possible_metabolites = session.scalars(
                    select(CompartmentalizedComponent)
                    .join(CompartmentalizedComponent.component)
                    .join(Component.annotation_mappings)
                    .join(ComponentAnnotationMapping.annotation)
                    .join(Annotation.links)
                    .filter(
                        AnnotationLink.data_source_id == seed_compound_data_source_id
                    )
                    .filter(AnnotationLink.identifier == cpd)
                    .distinct()
                ).all()
                if not possible_metabolites:
                    front = []
                    break
                old_front = front
                front = []
                for m in possible_metabolites:
                    for f_reactants, f_comp_map in old_front:
                        if (
                            seed_comp in f_comp_map
                            and f_comp_map[seed_comp] != m.compartment_id
                        ):
                            continue
                        new_comp_map = f_comp_map | {seed_comp: m.compartment_id}
                        new_reactants = f_reactants + [
                            {
                                "compartmentalized_component": m,
                                "coefficient": coefficient,
                            }
                        ]
                        if len(new_comp_map.values()) == len(
                            set(new_comp_map.values())
                        ):
                            front.append((new_reactants, new_comp_map))
            if front:
                hashes = [Reaction.generate_hash(f[0]) for f in front]
                pattern_matched_reaction_db = session.scalars(
                    select(Reaction).filter(Reaction.hash.in_(hashes))
                ).all()
            pattern_matched_reaction_db = {
                x.bigg_id: x for x in set(pattern_matched_reaction_db)
            }

        if not bigg_reaction_db and not pattern_matched_reaction_db:
            continue

        rxn_names = str(rxn.data.get("aliases", ""))
        rxn_names = rxn_names.split("|")
        rxn_names = [x[5:] for x in rxn_names if x.startswith("Name:")]
        rxn_names = rxn_names[0] if len(rxn_names) > 0 else ""
        rxn_names = [x.strip() for x in rxn_names.split(";")]

        annotation_bigg_id = f"seed.reaction:{rxn_id}"
        annotation_db = session.scalars(
            select(Annotation).filter(Annotation.bigg_id == annotation_bigg_id).limit(1)
        ).first()
        if annotation_db is None:
            annotation_db = Annotation(
                bigg_id=annotation_bigg_id,
                type="seed",
                default_data_source_id=get_data_source_id("seed.reaction", session),
            )
            session.add(annotation_db)

            for property in SEED_REACTION_PROPERTIES:
                prop_val = getattr(rxn, property, None)
                if prop_val is None:
                    prop_val = rxn.data.get(property)
                if isinstance(prop_val, str):
                    prop_val = prop_val.strip()
                    if len(prop_val) == 0:
                        prop_val = None
                if prop_val is None:
                    continue
                if pd.isna(prop_val):
                    continue
                if isinstance(prop_val, int) and property.startswith("is_"):
                    prop_val = bool(prop_val)
                if property == "is_obsolete":
                    annotation_db.is_obsolete = bool(prop_val)
                prop_db = AnnotationProperty(key=property)
                prop_db.value = prop_val
                annotation_db.properties.append(prop_db)

            for name in rxn_names:
                if name == getattr(rxn, "name", None):
                    continue
                prop_db = AnnotationProperty(key="name")
                prop_db.value = name
                annotation_db.properties.append(prop_db)

            for property, namespace in SEED_REACTION_ALIASES.items():
                if namespace == "seed.reaction":
                    prop_val = rxn_id
                else:
                    prop_val = aliases.get(property)
                    if isinstance(prop_val, str):
                        prop_val = prop_val.strip()
                        if len(prop_val) == 0:
                            prop_val = None
                    if prop_val is None:
                        continue
                    if pd.isna(prop_val):
                        continue
                if isinstance(prop_val, str):
                    prop_val = {prop_val}

                for val in prop_val:
                    alias_db = AnnotationLink(
                        identifier=val,
                        data_source_id=get_data_source_id(namespace, session),
                    )
                    annotation_db.links.append(alias_db)

        for k in set(bigg_reaction_db.keys()) | set(pattern_matched_reaction_db.keys()):
            br_db = bigg_reaction_db.get(k)
            pmr_db = pattern_matched_reaction_db.get(k)

            bigg_id_match = True if br_db is not None else None
            if bigg_id_match is None and k in rxn_bigg_ids:
                bigg_id_match = False

            pattern_match = pmr_db is not None
            annotation_mapping = ReactionAnnotationMapping(
                reaction=br_db if bigg_id_match else pmr_db,
                bigg_id_match=bigg_id_match,
                pattern_match=pattern_match,
            )

            found_match = False
            for rmap in annotation_db.reaction_mappings:
                if (
                    (rmap.reaction_id == annotation_mapping.reaction_id)
                    and (rmap.bigg_id_match == annotation_mapping.bigg_id_match)
                    and (rmap.pattern_match == annotation_mapping.pattern_match)
                ):
                    found_match = True
                    break
            if not found_match:
                annotation_db.reaction_mappings.append(annotation_mapping)
    session.commit()

    fill_universal_reaction_names_from_annotations(session)


# Values that appear in annotation name properties but carry no information.
PLACEHOLDER_NAMES = {"", "-", "null", "na", "n/a", "none", "unknown"}


def _is_usable_name(value):
    if value is None:
        return False
    stripped = value.strip()
    return len(stripped) > 2 and stripped.lower() not in PLACEHOLDER_NAMES


def fill_universal_reaction_names_from_annotations(session):
    """Give nameless universal reactions a name from their annotations.

    A universal reaction only gets a name at creation time, from either its
    reference reaction or the model file it came from
    (api.reactions.get_or_create_universal_reaction). Reactions that match
    neither are left with name=None, even when ModelSEED annotations carrying
    perfectly good names are attached later.

    This fills those in with the first usable annotation name, skipping
    placeholders like "-" and "" that ModelSEED uses for missing values.
    Annotation order is preserved, so the primary name wins over the aliases
    that follow it.

    Only touches rows where name IS NULL, so it never overwrites a name that
    came from a reference reaction or a model.

    Returns the number of reactions updated.
    """
    rows = session.execute(
        select(UniversalReaction.id, AnnotationProperty.value_str)
        .join(Reaction, Reaction.universal_reaction_id == UniversalReaction.id)
        .join(
            ReactionAnnotationMapping,
            ReactionAnnotationMapping.reaction_id == Reaction.id,
        )
        .join(
            AnnotationProperty,
            AnnotationProperty.annotation_id == ReactionAnnotationMapping.annotation_id,
        )
        .filter(UniversalReaction.name.is_(None))
        .filter(AnnotationProperty.key == "name")
        .order_by(UniversalReaction.id, AnnotationProperty.id)
    ).all()

    names = {}
    for universal_reaction_id, value in rows:
        if universal_reaction_id in names:
            continue
        if _is_usable_name(value):
            names[universal_reaction_id] = value.strip()

    for universal_reaction_id, name in names.items():
        session.get(UniversalReaction, universal_reaction_id).name = name
    session.commit()

    logging.info(f"Named {len(names)} universal reactions from annotations")
    return len(names)
