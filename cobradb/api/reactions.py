import logging
from operator import itemgetter
from typing import Dict, List, Optional, Tuple, Union
from sqlalchemy import func, select
from sqlalchemy.orm import Session

from cobradb.api import metabolites
from cobradb.data_sources import get_data_source_id
from cobradb.models import (
    Annotation,
    AnnotationLink,
    Compartment,
    CompartmentalizedComponent,
    ComponentReferenceMapping,
    ModelCollection,
    Reaction,
    ReactionIDMapping,
    ReactionMatrix,
    ReferenceCompound,
    ReferenceReaction,
    ReferenceReactionAnnotationMapping,
    ReferenceReactionParticipant,
    UniversalCompartmentalizedComponent,
    UniversalComponent,
    UniversalReaction,
    UniversalReactionMatrix,
)
from cobradb.api import utils

COMP_IN_OUT = [
    ("c", "p"),
    ("c", "e"),
    ("p", "e"),
    ("n", "c"),
    ("n", "p"),
    ("n", "e"),
]


class ReactionParticipantInfo:
    def __init__(
        self,
        coefficient: Union[float, int],
        compartmentalized_component: CompartmentalizedComponent,
    ):
        self.coefficient = coefficient
        self.compartmentalized_component = compartmentalized_component

    @property
    def compartment_id(self) -> int:
        return self.compartmentalized_component.compartment_id

    @property
    def charge(self):
        return self.compartmentalized_component.component.charge

    def __repr__(self):
        return f"<ReactionParticipantInfo cc='{self.compartmentalized_component.bigg_id} coeff='{self.coefficient}' charge='{self.charge}'>"


class UniversalReactionParticipantInfo:
    def __init__(
        self,
        coefficient: Union[float, int],
        universal_compartmentalized_component: UniversalCompartmentalizedComponent,
    ):
        self.coefficient = coefficient
        self.universal_compartmentalized_component = (
            universal_compartmentalized_component
        )

    @property
    def compartment_id(self) -> int:
        return self.universal_compartmentalized_component.compartment_id

    @property
    def charge(self):
        if (
            self.universal_compartmentalized_component.universal_component.default_component
            is None
        ):
            return 0
        return (
            self.universal_compartmentalized_component.universal_component.default_component.charge
        )

    def __repr__(self):
        return f"<UniversalReactionParticipantInfo ucc='{self.universal_compartmentalized_component.bigg_id} coeff='{self.coefficient}' charge='{self.charge}'>"


class ReferenceReactionParticipantInfo:
    def __init__(
        self,
        reference_compound: ReferenceCompound,
        universal_reaction_participant_info: UniversalReactionParticipantInfo,
        reference_reaction_participant: Optional[ReferenceReactionParticipant] = None,
    ):
        self.reference_compound = reference_compound
        self.universal_reaction_participant_info = universal_reaction_participant_info
        self.reference_reaction_participant = reference_reaction_participant

    @property
    def side(self):
        return "L" if self.universal_reaction_participant_info.coefficient < 0 else "R"

    @property
    def compartment_id(self) -> int:
        return self.universal_reaction_participant_info.compartment_id


def is_reaction_charge_balanced(
    reaction_participants: List[ReactionParticipantInfo],
) -> bool:
    overall_charge = 0
    for reaction_participant in reaction_participants:
        overall_charge += reaction_participant.coefficient * reaction_participant.charge
    return overall_charge == 0


def is_reaction_mass_balanced(
    reaction_participants: List[ReactionParticipantInfo],
) -> bool:
    overall_counts = {}
    for reaction_participant in reaction_participants:
        formula = reaction_participant.compartmentalized_component.component.formula
        d_formula = utils.formula_to_dict(formula)
        for k, v in d_formula.items():
            n = v * reaction_participant.coefficient
            if k in overall_counts:
                overall_counts[k] += n
            else:
                overall_counts[k] = n
    return all(x == 0 for x in overall_counts.values())


def _get_compartment_charge_dict(
    participants: Union[
        List[ReactionParticipantInfo], List[UniversalReactionParticipantInfo]
    ],
) -> Dict[int, int]:
    participants_charge = {}
    for p in participants:
        compartment_id = p.compartment_id
        if compartment_id not in participants_charge:
            participants_charge[compartment_id] = 0
        participants_charge[compartment_id] += p.charge * p.coefficient
    return participants_charge


def _corrected_charge_diffs(charge_dict_1, charge_dict_2, offset=0):
    charge_diffs = {
        c: int(charge_dict_1[c] - charge_dict_2[c])
        for c in (set(charge_dict_1.keys()) | set(charge_dict_2.keys()))
    }

    if (total_charge_diff := (sum(charge_diffs.values()) + offset)) == 0:
        return None

    # Add hydrons in compartments with the biggest difference first.
    corr_charge_diffs = {}
    remaining_charge_diff = total_charge_diff
    for c, charge_diff in reversed(sorted(charge_diffs.items(), key=itemgetter(1))):
        if charge_diff == 0:
            continue
        if (charge_diff > 0) != (total_charge_diff > 0):
            continue
        if abs(charge_diff) < abs(remaining_charge_diff):
            corr_charge_diffs[c] = charge_diff
            remaining_charge_diff -= charge_diff
        else:
            corr_charge_diffs[c] = remaining_charge_diff
            remaining_charge_diff = 0
        if remaining_charge_diff == 0:
            break

    if remaining_charge_diff != 0:
        if len(charge_diffs) == 1:
            k = next(iter(charge_diffs.keys()))
            if k not in corr_charge_diffs:
                corr_charge_diffs[k] = 0
            corr_charge_diffs[k] += remaining_charge_diff
            remaining_charge_diff = 0

    return corr_charge_diffs


def get_default_reaction_participants_list(
    session: Session,
    universal_reaction_participants: List[UniversalReactionParticipantInfo],
) -> Optional[List[ReactionParticipantInfo]]:

    reaction_participants = []
    for participant in universal_reaction_participants:
        dc_db = (
            participant.universal_compartmentalized_component.universal_component.default_component
        )
        if dc_db is None:
            return None
        cc_db = metabolites.get_or_create_compartmentalized_component(
            session,
            dc_db,
            participant.universal_compartmentalized_component.compartment,
        )
        rp = ReactionParticipantInfo(
            compartmentalized_component=cc_db,
            coefficient=participant.coefficient,
        )
        reaction_participants.append(rp)
    return reaction_participants


def get_universal_reaction_participants_list(
    session: Session,
    reaction_participants: List[ReactionParticipantInfo],
    balance_charge: bool = False,
) -> List[UniversalReactionParticipantInfo]:
    universal_participants = [
        UniversalReactionParticipantInfo(
            universal_compartmentalized_component=participant.compartmentalized_component.universal_compartmentalized_component,
            coefficient=participant.coefficient,
        )
        for participant in reaction_participants
    ]

    # Find the difference in charge and correct this using H+.
    participants_charge = _get_compartment_charge_dict(reaction_participants)
    universal_participants_charge = _get_compartment_charge_dict(universal_participants)

    offset = 0
    if balance_charge:
        offset = -sum(participants_charge.values())

    corr_charge_diffs = _corrected_charge_diffs(
        participants_charge, universal_participants_charge, offset
    )
    if corr_charge_diffs is None:
        return universal_participants

    universal_participants_copy = universal_participants[:]
    for compartment_id, charge_diff in corr_charge_diffs.items():
        if charge_diff == 0:
            continue
        for p in universal_participants:
            # Update existing entry if possible
            ucp = p.universal_compartmentalized_component
            if (
                ucp.compartment_id == compartment_id
                and ucp.universal_component.bigg_id == "h"
            ):
                p.coefficient += charge_diff
                break
        else:
            # Create a new H+ entry
            hydron_db = session.scalars(
                select(UniversalCompartmentalizedComponent)
                .join(UniversalCompartmentalizedComponent.universal_component)
                .filter(
                    UniversalCompartmentalizedComponent.compartment_id == compartment_id
                )
                .filter(UniversalComponent.bigg_id == "h")
                .limit(1)
            ).first()
            if not hydron_db:
                universal_participants = universal_participants_copy
                break
            universal_participants.append(
                UniversalReactionParticipantInfo(
                    universal_compartmentalized_component=hydron_db,
                    coefficient=charge_diff,
                )
            )
    # Remove any entries with coefficient 0
    universal_participants = [p for p in universal_participants if p.coefficient != 0]
    return universal_participants


def get_reference_reaction_participants_list(
    session: Session,
    universal_reaction_participants: List[UniversalReactionParticipantInfo],
) -> Optional[List[ReferenceReactionParticipantInfo]]:

    reference_reaction_participants = []
    for universal_participant in universal_reaction_participants:
        universal_component_db = (
            universal_participant.universal_compartmentalized_component.universal_component
        )
        if universal_component_db.reference_mapping is None:
            return None

        reference_reaction_participants.append(
            ReferenceReactionParticipantInfo(
                reference_compound=universal_component_db.reference_mapping.mapping.reference_compound,
                universal_reaction_participant_info=universal_participant,
            )
        )
    return reference_reaction_participants


def match_reference_reaction_participants_with_reference_reaction(
    session: Session,
    reference_reaction_participants: List[ReferenceReactionParticipantInfo],
    reference_reaction_db: ReferenceReaction,
):
    reaction_mappings = [({}, {}, reference_reaction_participants[:], [])]
    participants_query = select(
        ReferenceReactionParticipant,
    ).filter(ReferenceReactionParticipant.reaction_id == reference_reaction_db.id)

    if (
        len(reference_reaction_participants)
        != session.execute(
            select(func.count()).select_from(participants_query.subquery())
        ).scalar_one()
    ):
        return None
    participants_db = session.scalars(participants_query)

    for participant in participants_db:
        requires_n = "n" in str(participant.coefficient).lower()

        # print(f"Participant: {participant} -> {compound}")
        crmappings_db = session.scalars(
            select(ComponentReferenceMapping)
            .join(
                ComponentReferenceMapping.universal_component,
            )
            .join(ComponentReferenceMapping.component)
            .filter(
                ComponentReferenceMapping.reference_compound_id
                == participant.compound_id,
            )
        ).all()
        participant_mappings = []
        for crmapping in crmappings_db:
            coefficient = participant.coefficient
            if requires_n:
                if crmapping.reference_n is None:
                    continue
                coefficient = coefficient.lower().replace(
                    "n", str(crmapping.reference_n)
                )
            coefficient = float(coefficient)

            ms = [
                p
                for p in reference_reaction_participants
                if p.universal_reaction_participant_info.universal_compartmentalized_component.universal_component_id
                == crmapping.universal_component_id
                and float(abs(p.universal_reaction_participant_info.coefficient))
                == coefficient
                and (
                    float(p.universal_reaction_participant_info.charge)
                    == crmapping.component.charge
                )
            ]
            if not ms:
                continue
            for m in ms:
                participant_mappings.append((participant, m, crmapping))
        if not participant_mappings:
            reaction_mappings = []
            break
        old_front = reaction_mappings
        reaction_mappings = []
        for front in old_front:
            for (
                participant_db,
                reference_reaction_participant,
                crmapping,
            ) in participant_mappings:
                if reference_reaction_participant not in front[2]:
                    continue
                new_item = (dict(front[0]), dict(front[1]), front[2][:], front[3][:])
                for i, att in enumerate(
                    [("side", "side"), ("compartment", "compartment_id")]
                ):
                    if new_item is None:
                        break
                    att_m = getattr(participant_db, att[0])
                    att_r = getattr(reference_reaction_participant, att[1])
                    if att_m in front[i]:
                        if att_r != front[i][att_m]:
                            print(
                                f"Inconsistent {att} mapping: prev: {att_m}:{front[i][att_m]}, new: {att_m}:{att_r}"
                            )
                            new_item = None
                    elif att_r in front[i].values():
                        print(
                            f"Inconsistent {att} mapping: value:{att_r} already in {front[i]}"
                        )
                        new_item = None
                    else:
                        new_item[i][att_m] = att_r
                if new_item is not None:
                    new_item[3].append(
                        (
                            participant_db,
                            reference_reaction_participant,
                            crmapping,
                        )
                    )
                    new_item[2].remove(reference_reaction_participant)
                    reaction_mappings.append(new_item)
    successful_mappings = []
    for m in reaction_mappings:
        # TODO: reference_n should also be consistent
        if m[2]:
            continue
        if "in" in m[1] and "out" in m[1]:
            if not (m[1]["in"], m[1]["out"]) in COMP_IN_OUT:
                print(
                    f"mapping in: {m[1]['in']}, out: {m[1]['out']} not a known in/out pair."
                )
                continue
        successful_mappings.append(m)
    if not successful_mappings:
        print(f"No successful mappings: {reaction_mappings}")
        return None
    elif len(successful_mappings) > 1:
        print(f"Multiple successful mappings: {reaction_mappings}")
        return None
    else:
        result = []
        for (
            participant_db,
            reference_reaction_participant,
            crmapping,
        ) in successful_mappings[0][3]:
            reference_reaction_participant.reference_reaction_participant = (
                participant_db
            )
            result.append(reference_reaction_participant)
        return result


def find_reference_reaction(
    session: Session,
    reference_reaction_participants: List[ReferenceReactionParticipantInfo],
):
    reference_db, mapping = None, None

    hash_pattern = ReferenceReaction.generate_hash(
        reference_reaction_participants, pattern=True
    )
    logging.warning(f"hash pattern: {hash_pattern}")

    db_entries = session.scalars(
        select(ReferenceReaction).where(
            ReferenceReaction.hash.regexp_match(hash_pattern)
        )
    ).all()

    logging.warning(f"db entries: {db_entries}")
    if not db_entries:
        return None, None
    db_entries = list(db_entries)
    # print(db_entries)

    for db_entry in db_entries:
        # print(f"Matching {db_entry.id}: {db_entry.equation}")
        m = match_reference_reaction_participants_with_reference_reaction(
            session, reference_reaction_participants, db_entry
        )
        print(f"match: {m}")
        # print(m)
        if m is not None:
            reference_db, mapping = db_entry, m
            break

    return reference_db, mapping


def get_or_create_universal_reaction(
    session: Session,
    bigg_id: str,
    universal_reaction_participants: List[UniversalReactionParticipantInfo],
    reaction_collection_id: Optional[int] = None,
    exchange_reaction: bool = False,
    reaction_name: Optional[str] = None,
) -> UniversalReaction:
    universal_reaction_hash = UniversalReaction.generate_hash(
        universal_reaction_participants
    )

    universal_reaction_db = session.scalars(
        select(UniversalReaction)
        .filter(UniversalReaction.hash == universal_reaction_hash)
        .filter(
            (UniversalReaction.collection_id == None)
            | (UniversalReaction.collection_id == reaction_collection_id)
        )
        .limit(1)
    ).first()
    if universal_reaction_db is not None:
        return universal_reaction_db

    reference_db = None
    reference_reaction_participants = get_reference_reaction_participants_list(
        session, universal_reaction_participants
    )
    logging.warning(
        f"reference reaction participants: {reference_reaction_participants}"
    )
    if reference_reaction_participants is not None:
        logging.warning("find_reference_reaction")
        reference_db, reference_reaction_participants = find_reference_reaction(
            session, reference_reaction_participants
        )
    logging.warning(f"reference_db: {reference_db}")

    if reaction_collection_id is None and exchange_reaction:
        reaction_collection_id = universal_reaction_participants[
            0
        ].universal_compartmentalized_component.universal_component.collection_id

    if reaction_collection_id is None:
        reaction_collection = None
    else:
        reaction_collection = session.get(ModelCollection, reaction_collection_id)

    is_biomass_reaction = reference_db is None and "BIOMASS" in bigg_id.upper()
    if is_biomass_reaction and reaction_collection is not None:
        bare_bigg_id = bigg_id.removeprefix(f"__{reaction_collection.bigg_id}__")
        other_universal_reaction_db = session.scalars(
            select(UniversalReaction)
            .filter(UniversalReaction.bigg_id == bare_bigg_id)
            .filter(
                (UniversalReaction.collection_id != None)
                & (UniversalReaction.collection_id != reaction_collection.id)
            )
            .limit(1)
        ).first()
        if not other_universal_reaction_db:
            bigg_id = bare_bigg_id
        biomass_reference_db = session.scalars(
            select(ReferenceReaction)
            .filter(ReferenceReaction.bigg_id == "BiGGr:BIOMASS")
            .limit(1)
        ).first()
        if biomass_reference_db is None:
            biomass_reference_db = ReferenceReaction(
                bigg_id="BiGGr:BIOMASS",
                name="Biomass reaction",
                hash="#BIOMASS",
            )
            session.add(biomass_reference_db)
        reference_db = biomass_reference_db

    for i in range(100):
        proposed_bigg_id = bigg_id if i == 0 else f"{bigg_id}_{i}"
        universal_reaction_db_by_bigg_id = session.scalars(
            select(UniversalReaction)
            .filter(UniversalReaction.bigg_id == proposed_bigg_id)
            .limit(1)
        ).first()
        if universal_reaction_db_by_bigg_id is None:
            bigg_id = proposed_bigg_id
            break

    print("Creating new universal reaction")
    if reference_db is not None:
        if reference_db.name is not None and reference_db.name.strip() != "":
            reaction_name = reference_db.name
    universal_reaction_db = UniversalReaction(
        bigg_id=bigg_id,
        name=reaction_name,
        reference=reference_db,
        hash=universal_reaction_hash,
        collection=reaction_collection,
        is_exchange=exchange_reaction,
        is_pseudo=is_biomass_reaction,
    )

    compartment_ids = set()
    if reference_reaction_participants is None:
        for participant in universal_reaction_participants:
            ucc = participant.universal_compartmentalized_component
            compartment_ids.add(ucc.compartment_id)
            urm_db = UniversalReactionMatrix(
                universal_compartmentalized_component=ucc,
                reference_reaction_participant=None,
                coefficient=participant.coefficient,
            )
            universal_reaction_db.matrix.append(urm_db)
    else:
        for participant in reference_reaction_participants:
            ucc = (
                participant.universal_reaction_participant_info.universal_compartmentalized_component
            )
            compartment_ids.add(ucc.compartment_id)
            urm_db = UniversalReactionMatrix(
                universal_compartmentalized_component=ucc,
                reference_reaction_participant=participant.reference_reaction_participant,
                coefficient=participant.universal_reaction_participant_info.coefficient,
            )
            universal_reaction_db.matrix.append(urm_db)
    if len(compartment_ids) > 1:
        universal_reaction_db.is_transport = True
    session.add(universal_reaction_db)
    return universal_reaction_db


def get_or_create_reaction_for_universal_reaction(
    session: Session,
    universal_reaction_db: UniversalReaction,
    reaction_participants: List[ReactionParticipantInfo],
) -> Optional[Reaction]:
    reaction_hash = Reaction.generate_hash(reaction_participants)

    # Get the reaction
    reaction_db = session.scalars(
        select(Reaction)
        .filter(
            (Reaction.hash == reaction_hash)
            & (
                (Reaction.collection_id == None)
                | (Reaction.collection == universal_reaction_db.collection)
            )
        )
        .limit(1)
    ).first()
    if reaction_db is not None:
        logging.warning(f"Reaction '{reaction_hash}' already exists.")
        return reaction_db

    universal_matrix_orig = {}
    for x in universal_reaction_db.matrix:
        ucc_id = x.universal_compartmentalized_component_id
        universal_matrix_orig[ucc_id] = universal_matrix_orig.get(ucc_id, []) + [x]

    reaction_matrix = []
    for is_reversed in [False, True]:
        universal_matrix = universal_matrix_orig.copy()
        for participant in reaction_participants:
            compartmentalized_component_db = participant.compartmentalized_component
            original_coefficient = participant.coefficient
            coefficient = -original_coefficient if is_reversed else original_coefficient

            matched_matrices = universal_matrix.get(
                compartmentalized_component_db.universal_compartmentalized_component_id,
                [],
            )

            matched_matrix = None
            for mm in matched_matrices:
                if mm.coefficient == coefficient:
                    matched_matrix = mm
                    break
            else:
                logging.warning(
                    f"Coefficients do not match: {[matched_matrix.coefficient for matched_matrix in matched_matrices]} vs. {coefficient}"
                )
                matched_matrix = None

            if matched_matrix is None:
                # No match is only allowed for H+
                if (
                    c_bigg_id := compartmentalized_component_db.component.bigg_id
                ) != "h:1":
                    # print(
                    #     f"No correct match for {compartmentalized_component_db.bigg_id} when creating a reaction for univeral reaction {universal_reaction_db.bigg_id}."
                    # )
                    logging.error(f"No correct match for {c_bigg_id}.")
                    reaction_matrix = []
                    break
                matched_matrix = UniversalReactionMatrix(
                    universal_compartmentalized_component=compartmentalized_component_db.universal_compartmentalized_component,
                    coefficient=coefficient,
                )
            else:
                ucc_id = (
                    compartmentalized_component_db.universal_compartmentalized_component_id
                )
                universal_matrix[ucc_id].remove(matched_matrix)
                if len(universal_matrix[ucc_id]) == 0:
                    del universal_matrix[ucc_id]
            reaction_matrix.append(
                ReactionMatrix(
                    universal_reaction_matrix=matched_matrix,
                    compartmentalized_component=compartmentalized_component_db,
                )
            )
        else:
            for unmatched_matrix in universal_matrix.values():
                if (
                    uc_bigg_id := unmatched_matrix[
                        0
                    ].universal_compartmentalized_component.universal_component.bigg_id
                ) != "h":
                    # print(
                    #     f"No match for {unmatched_matrix.universal_compartmentalized_component} when creating a reaction for univeral reaction {universal_reaction_db.bigg_id}."
                    # )
                    logging.warning(f"Unmatched entries that are not 'h': {uc_bigg_id}")
                    reaction_matrix = []
                    break
        if reaction_matrix:
            break

    if not reaction_matrix:
        logging.error("No reaction matrix could be constructed.")
        return None

    copy_number = (
        session.scalar(
            select(func.count(Reaction.id)).filter(
                Reaction.universal_reaction == universal_reaction_db
            )
        )
        + 1
    )
    reaction_bigg_id = Reaction.create_id(universal_reaction_db.bigg_id, copy_number)

    logging.info(f"Creating new reaction '{reaction_bigg_id}'")
    reaction_db = Reaction(
        bigg_id=reaction_bigg_id,
        hash=reaction_hash,
        collection_id=universal_reaction_db.collection_id,
        copy_number=copy_number,
        universal_reaction=universal_reaction_db,
    )
    reaction_db.matrix.extend(reaction_matrix)
    session.add(reaction_db)
    return reaction_db


def _create_reference_reaction(
    session: Session,
    bigg_id: str,
    equation: str,
    participants: List[Tuple[int, float, str, Compartment]],
    annotations: Optional[Dict[str, List[str]]],
) -> Optional[ReferenceReaction]:
    reaction_db = ReferenceReaction(
        bigg_id=bigg_id,
        equation=equation,
    )
    for side_n, coefficient, compound_id, compartment in participants:
        compound_db = utils.get_object_by_bigg_id(
            session, compound_id, ReferenceCompound
        )
        participant_db = ReferenceReactionParticipant(
            compound=compound_db,
            side="L" if side_n < 0 else "R",
            coefficient=coefficient,
            compartment=compartment,
        )
        reaction_db.reaction_participants.append(participant_db)
    if len(reaction_db.reaction_participants) > 0:
        reaction_db.update_hash()
        session.add(reaction_db)

        add_reaction_annotations(session, reaction_db, annotations)
        return reaction_db
    return None


def add_reference_exchange_reactions(
    session: Session, compound_db: ReferenceCompound
) -> Optional[ReferenceReaction]:
    reaction_bigg_id = f"EX:{compound_db.bigg_id}"
    reference_reaction_db = session.scalars(
        select(ReferenceReaction)
        .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
        .limit(1)
    ).first()
    if reference_reaction_db is not None:
        return
    compound_name = compound_db.name if compound_db.name else compound_db.id
    reference_reaction_db = ReferenceReaction(
        bigg_id=reaction_bigg_id,
        name=f"Exchange of {compound_name}.",
        equation=f"{compound_name} = ∅",
    )
    lhs_part = ReferenceReactionParticipant(
        compound=compound_db,
        side="L",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(lhs_part)
    reference_reaction_db.update_hash()
    session.add(reference_reaction_db)
    session.commit()
    return reference_reaction_db


def add_reaction_annotations(
    session: Session,
    reaction_db: ReferenceReaction,
    annotations: Optional[Dict[str, List[str]]],
):
    if annotations is None:
        return
    default_data_source_id = get_data_source_id("rhea", session)
    if default_data_source_id is None:
        print("Could not find RHEA data source.")
        return
    annotation_db = Annotation(
        bigg_id=reaction_db.bigg_id,
        default_data_source_id=default_data_source_id,
        type="rhea",
    )
    mapping = ReferenceReactionAnnotationMapping(
        reference_reaction=reaction_db,
    )
    annotation_db.reference_reaction_mappings.append(mapping)

    for namespace, identifiers in annotations.items():
        data_source_id = get_data_source_id(namespace, session)
        if not data_source_id:
            print(f"Unknown data source: {namespace}")
            continue
        for identifier in identifiers:
            link = AnnotationLink(
                data_source_id=data_source_id,
                identifier=identifier,
            )
            annotation_db.links.append(link)
    session.add(annotation_db)


def create_reaction_id_mapping(
    session: Session, old_bigg_id: str, universal_reaction_db: UniversalReaction
) -> Optional[ReactionIDMapping]:
    """Record old_bigg_id as a superseded identifier of universal_reaction_db."""
    if old_bigg_id == universal_reaction_db.bigg_id:
        return None
    if session.get(ReactionIDMapping, old_bigg_id) is not None:
        return None
    # Check that we're not making an entry inaccessible.
    if utils.get_object_by_bigg_id(session, old_bigg_id, UniversalReaction) is not None:
        return None
    id_mapping_db = ReactionIDMapping(
        old_bigg_id=old_bigg_id,
        new_universal_reaction=universal_reaction_db,
    )
    session.add(id_mapping_db)
    return id_mapping_db
