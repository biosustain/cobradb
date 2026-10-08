from typing import List, Optional
from sqlalchemy.orm import Session
from cobradb.api import metabolites
from cobradb.api import reactions
from cobradb.api.bigg_ids import create_component_bigg_id
from cobradb.api.reactions import (
    ReactionParticipantInfo,
    get_default_reaction_participants_list,
    get_universal_reaction_participants_list,
    is_reaction_charge_balanced,
    is_reaction_mass_balanced,
)
from cobradb.api.utils import get_object_by_bigg_id
from cobradb.models import (
    Compartment,
    CompartmentalizedComponent,
    Component,
    Model,
    ModelCollection,
    Reaction,
    ReferenceReaction,
    UniversalReaction,
)
from cobradb.util import timing

from sqlalchemy import select
import logging
import json


def load_bigg_id_data(filename):
    with open(filename, "r") as f:
        return json.load(f)


def parse_reaction_participants(
    session: Session, participants
) -> Optional[List[ReactionParticipantInfo]]:
    parsed_participants = []
    for i, side in enumerate([-1, 1]):
        for p in participants[i]:
            full_id = p[1]
            if full_id == "":
                continue
            if ":" in full_id:
                universal_comp_comp_id, charge = full_id.rsplit(":", maxsplit=1)
            else:
                universal_comp_comp_id, charge = full_id, None
            universal_id, compartment = universal_comp_comp_id.rsplit("_", maxsplit=1)

            universal_component_db = metabolites.get_universal_component_by_bigg_id(
                session, universal_id
            )
            if universal_component_db is None:
                # logging.warning(f"Participant {full_id}: No universal component found.")
                return None
            compartment_db = get_object_by_bigg_id(session, compartment, Compartment)
            if compartment_db is None:
                # logging.warning(f"Participant {full_id}: No compartment found.")
                return None

            if charge is None:
                component_db = universal_component_db.default_component
            else:
                component_db = session.scalars(
                    select(Component)
                    .filter(Component.charge == float(charge))
                    .filter(Component.variant == 0)
                    .filter(
                        Component.universal_component_id == universal_component_db.id
                    )
                    .limit(1)
                ).first()
            if component_db is None:
                # logging.warning(f"Participant {full_id}: No component found.")
                return None

            compartmentalized_component_db = session.scalars(
                select(CompartmentalizedComponent)
                .filter(CompartmentalizedComponent.component_id == component_db.id)
                .filter(CompartmentalizedComponent.compartment_id == compartment_db.id)
                .limit(1)
            ).first()
            if compartmentalized_component_db is None:
                universal_compartmentalized_component_db = (
                    metabolites.get_or_create_universal_compartmentalized_component(
                        session, universal_component_db, compartment_db
                    )
                )
                cc_bigg_id = create_component_bigg_id(
                    base_bigg_id=universal_component_db.bigg_id,
                    charge=component_db.charge,
                    variant=component_db.variant,
                    compartment_bigg_id=compartment_db.bigg_id,
                )
                compartmentalized_component_db = CompartmentalizedComponent(
                    bigg_id=cc_bigg_id,
                    component=component_db,
                    compartment=compartment_db,
                    universal_compartmentalized_component=universal_compartmentalized_component_db,
                )
                session.add(compartmentalized_component_db)
            participant = ReactionParticipantInfo(
                coefficient=side * p[0],
                compartmentalized_component=compartmentalized_component_db,
            )
            parsed_participants.append(participant)
    return parsed_participants


def match_curated_reaction(session: Session, bigg_id, reaction_data):
    """Parse and balance a curated reaction and look it up in the database.

    Returns None when the reaction cannot be loaded. Otherwise returns a tuple
    (reaction_participants, reaction_hash, collection_cond, reaction_db,
    universal_reaction_db), where reaction_db is an existing reaction with the
    same participants and universal_reaction_db an existing universal reaction
    with the same universal participants (looked up only when reaction_db is
    None); either is None when no match exists.
    """
    reaction_participants = parse_reaction_participants(
        session, reaction_data["participants"]
    )

    # logging.warning(f"reaction participants: {reaction_participants}")
    if reaction_participants is None:
        # logging.error(
        #     f"Could not parse reaction participants for '{bigg_id}' ({reaction_data['participants']})"
        # )
        return None

    reaction_collection_bigg_id = reaction_data.get("collection_bigg_id")
    if reaction_collection_bigg_id is None:
        reaction_collection = None
    else:
        reaction_collection = session.scalars(
            select(ModelCollection)
            .filter(ModelCollection.bigg_id == reaction_collection_bigg_id)
            .limit(1)
        ).first()

    # logging.warning(f"reaction collection: {reaction_collection}")

    if reaction_collection is None:
        collection_cond = lambda x: x.collection_id == None
    else:
        reaction_collection_id = reaction_collection.id
        collection_cond = lambda x: (x.collection_id == None) | (
            x.collection_id == reaction_collection_id
        )

    if not (charge_balanced := is_reaction_charge_balanced(reaction_participants)):
        logging.error(
            f"Reaction template is unbalanced for reaction {bigg_id}, fixing."
        )
        universal_reaction_participants = get_universal_reaction_participants_list(
            session, reaction_participants, balance_charge=True
        )
        reaction_participants = get_default_reaction_participants_list(
            session, universal_reaction_participants
        )
        if reaction_participants is None:
            logging.error(
                f"Failed to update reaction participants to correct charge of reaction {bigg_id}."
            )
            return None

    if not (charge_balanced := is_reaction_charge_balanced(reaction_participants)):
        logging.error(f"Failed to correct charge of reaction {bigg_id}.")
        return None

    if not (mass_balanced := is_reaction_mass_balanced(reaction_participants)):
        logging.error(f"Reaction {bigg_id} is not mass balanced.")
        return None

    reaction_hash = Reaction.generate_hash(reaction_participants)
    # logging.warning(f"Reaction hash: {reaction_hash}")

    # Get the reaction
    reaction_db = session.scalars(
        select(Reaction)
        .filter(Reaction.hash == reaction_hash)
        .filter(collection_cond(Reaction))
        .join(Reaction.universal_reaction)
        .limit(1)
    ).first()

    # logging.warning(f"reaction_db: {reaction_db}")
    universal_reaction_db = None
    if reaction_db is None:
        universal_participants = reactions.get_universal_reaction_participants_list(
            session, reaction_participants
        )
        # logging.warning(f"universal participants: {universal_participants}")
        universal_reaction_hash = UniversalReaction.generate_hash(
            universal_participants
        )
        # logging.warning(f"universal hash: {universal_reaction_hash}")
        universal_reaction_db = session.scalars(
            select(UniversalReaction)
            .filter(UniversalReaction.hash == universal_reaction_hash)
            .filter(collection_cond(UniversalReaction))
            .limit(1)
        ).first()

    return (
        reaction_participants,
        reaction_hash,
        collection_cond,
        reaction_db,
        universal_reaction_db,
    )


@timing
def push_reactions(session: Session, data):
    logging.warning("function push reaction")
    for n, (bigg_id, reaction_data) in enumerate(data.items()):
        match = match_curated_reaction(session, bigg_id, reaction_data)
        if match is None:
            continue
        (
            reaction_participants,
            reaction_hash,
            collection_cond,
            reaction_db,
            universal_reaction_db,
        ) = match

        if reaction_db is not None:
            logging.warning(
                f"Reaction '{bigg_id}' already in database as '{reaction_db.bigg_id}', skipping."
            )
            # Keep the redundant BiGG ID as a superseded identifier.
            reactions.create_reaction_id_mapping(
                session, bigg_id, reaction_db.universal_reaction
            )
            session.commit()
            continue

        # logging.warning(f"universal reaction_db: {universal_reaction_db}")
        if universal_reaction_db is None:
            # logging.warning("get_or_create_universal_reaction")
            universal_participants = (
                reactions.get_universal_reaction_participants_list(
                    session, reaction_participants
                )
            )
            universal_reaction_db = reactions.get_or_create_universal_reaction(
                session,
                bigg_id,
                universal_participants,
            )
            session.commit()
        else:
            # A proton-stoichiometry variant of an existing universal reaction.
            reactions.create_reaction_id_mapping(
                session, bigg_id, universal_reaction_db
            )
        # logging.warning(f"universal reaction_db: {universal_reaction_db}")
        reaction_db = reactions.get_or_create_reaction_for_universal_reaction(
            session, universal_reaction_db, reaction_participants
        )
        session.commit()

        reaction_db = session.scalars(
            select(Reaction)
            .filter(Reaction.hash == reaction_hash)
            .filter(collection_cond(Reaction))
            .limit(1)
        ).first()
        if reaction_db is None:
            logging.error(f"Reaction {bigg_id} was not correctly created.")
            continue

        session.commit()


@timing
def load_reactions(session: Session, curated_reactions_filepath):
    logging.debug("Loading Curated Reactions reference data")

    data = load_bigg_id_data(curated_reactions_filepath)

    push_reactions(session, data)
