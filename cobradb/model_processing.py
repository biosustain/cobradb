# -*- coding: utf-8 -*-
import gzip
from os import PathLike
from pathlib import Path
from typing import Union
from cobra.io import save_json_model, save_yaml_model, write_sbml_model
from sqlalchemy.orm import Session, contains_eager, joinedload

from cobradb.models import (
    CompartmentalizedComponent,
    Component,
    ComponentReferenceMapping,
    Model,
    ModelCompartmentalizedComponent,
    ModelReaction,
    Reaction,
    ReferenceReaction,
    UniversalReaction,
)
from cobradb import parse
from cobradb.api import utils
from cobradb.util import (
    timing,
)

from sqlalchemy import select
import logging


def create_model_base_filename(session: Session, bigg_id: str) -> str:
    base_filename = "".join(
        [x for x in bigg_id if x.isalpha() or x.isdigit() or x in "_"]
    )
    i = 0
    while True:
        if i == 0:
            proposed_filename = base_filename
        else:
            import hashlib

            filename_hash = hashlib.sha256(f"{base_filename}{i}".encode()).hexdigest()
            proposed_filename = f"{base_filename}_{filename_hash[:8]}"
        model_db = session.scalars(
            select(Model)
            .filter(Model.base_filename == proposed_filename)
            .filter(Model.bigg_id != bigg_id)
            .limit(1)
        ).first()
        if not model_db:
            return proposed_filename
        i += 1
    raise Exception("Should not occur.")


@timing
def process_model(session: Session, model_data, model_filepath: Union[str, PathLike]):
    # apply id normalization
    logging.debug("Parsing SBML")
    model, old_parsed_ids = parse.load_and_normalize(model_filepath)
    if model_data.get("prefix") is not None:
        model.id = f"{model_data['prefix']}{model.id}"
    model_bigg_id = model.id

    # check that the model exists
    if (model_db := utils.get_object_by_bigg_id(session, model_bigg_id, Model)) is None:
        return None, None
    model_db_id = model_db.id

    process_metabolites(
        session,
        model,
        model_db_id,
    )

    process_reactions(
        session,
        model,
        model_db_id,
    )

    base_filename = create_model_base_filename(session, model_db.bigg_id)
    model_output_path = Path("/models/models/") / base_filename
    logging.warning(f"Writing corrected model to: {model_output_path}")
    write_sbml_model(model, model_output_path.with_suffix(".biggr.sbml"))
    save_json_model(model, model_output_path.with_suffix(".biggr.json"))
    save_yaml_model(model, str(model_output_path.with_suffix(".biggr.yaml")))
    for ext in [".biggr.sbml", ".biggr.json", ".biggr.yaml"]:
        with open(model_output_path.with_suffix(ext), "rb") as f_in, gzip.open(
            model_output_path.with_suffix(f"{ext}.gz"), "wb"
        ) as f_out:
            f_out.writelines(f_in)
    model_db.base_filename = base_filename
    session.commit()

    return model_bigg_id, model_db_id


@timing
def process_metabolites(session, model, model_db_id):
    model_db = session.get(Model, model_db_id)

    # Make sure there are no clashes when renaming
    for metabolite in model.metabolites:
        metabolite.id = f"_____{metabolite.id}"

    # for each metabolite in the model
    for metabolite in model.metabolites:
        metabolite_id = metabolite.id[5:]
        print(f"Metabolite: {metabolite_id}")

        # Check if there is an existing entry that exactly matches
        model_comp_comp_db = session.scalars(
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
            .filter(ModelCompartmentalizedComponent.model == model_db)
            .filter(
                ModelCompartmentalizedComponent.id_in_original_model == metabolite_id
            )
            .limit(1)
        ).first()
        if model_comp_comp_db is None:
            raise ValueError()

        comp_component_db = model_comp_comp_db.compartmentalized_component
        universal_comp_component_db = (
            comp_component_db.universal_compartmentalized_component
        )
        metabolite_db = comp_component_db.component

        metabolite.id = model_comp_comp_db.bigg_id
        metabolite.name = metabolite_db.name
        metabolite.charge = float(metabolite_db.charge)
        metabolite.formula = metabolite_db.formula

        metabolite.annotation["biggr"] = comp_component_db.bigg_id
        metabolite.annotation["sbo"] = "SBO:0000247"

        comp_ref_map = session.scalars(
            select(ComponentReferenceMapping)
            .join(ComponentReferenceMapping.reference_compound)
            .filter(ComponentReferenceMapping.component == metabolite_db)
        ).all()
        if comp_ref_map:
            chebis = []
            for crm in comp_ref_map:
                if crm.reference_compound.bigg_id.startswith("CHEBI:"):
                    chebis.append(crm.reference_compound.bigg_id)
                metabolite.annotation["sbo"] = crm.reference_compound.get_sbo()

            if chebis:
                metabolite.annotation["chebi"] = chebis


@timing
def process_reactions(
    session,
    model,
    model_db_id,
):
    model_db = session.get(Model, model_db_id)

    # Make sure there are no clashes when renaming
    for reaction in model.reactions:
        reaction.id = f"_____{reaction.id}"

    for reaction in model.reactions:
        reaction_id = reaction.id[5:]
        print(f"Reaction: {reaction_id}")

        model_reaction_db = session.scalars(
            select(ModelReaction)
            .options(
                joinedload(ModelReaction.reaction).joinedload(
                    Reaction.universal_reaction
                )
            )
            .filter(
                (ModelReaction.model == model_db)
                & (ModelReaction.id_in_original_model == reaction_id)
            )
            .limit(1)
        ).first()
        if model_reaction_db is None:
            raise ValueError()

        reaction.id = model_reaction_db.bigg_id
        reaction_name = model_reaction_db.reaction.universal_reaction.name
        reaction.name = "" if reaction_name is None else reaction_name
        # reaction.bounds = (model_reaction_db.lower_bound, model_reaction_db.upper_bound)
        # TODO: Handle flipping reactions

        reaction.annotation.clear()
        reaction.annotation["biggr"] = model_reaction_db.reaction.bigg_id

        reaction.annotation["sbo"] = "SBO:0000176"

        uni_ref_db = session.execute(
            select(UniversalReaction, ReferenceReaction)
            .outerjoin(
                ReferenceReaction,
                UniversalReaction.reference_id == ReferenceReaction.id,
            )
            .filter(
                UniversalReaction.id == model_reaction_db.reaction.universal_reaction_id
            )
            .limit(1)
        ).first()

        if uni_ref_db is not None:
            universal_db, reference_db = uni_ref_db
            reaction.annotation["sbo"] = universal_db.get_sbo(reference_db)
            if reference_db is not None:
                if reference_db.bigg_id.startswith("RHEA:"):
                    reaction.annotation["rhea"] = reference_db.bigg_id
