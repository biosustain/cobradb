from pathlib import Path
from typing import Dict, Optional, Tuple
from libchebipy import ChebiEntity as lcpChebiEntity

from libchebipy._chebi_entity import parsers
import libchebipy._parsers as chebi_parsers

import pandas as pd
from sqlalchemy import select
from sqlalchemy.orm import Session

from cobradb import settings
from cobradb.models import (
    ReferenceCompound,
    ReferenceReaction,
    ReferenceReactionParticipant,
)
from cobradb.api import utils

parsers.set_auto_update(False)
parsers.set_download_cache_path(Path(settings.chebi_cache_dir) / "libChEBI")


def load_default_chebi_mapping() -> Dict[str, str]:
    try:
        df = pd.read_csv(settings.default_chebis, sep="\t", index_col=0, header=0)
        df.index = "CHEBI:" + df.index.astype(str)
        df["CHEBI_PH7_3"] = "CHEBI:" + df["CHEBI_PH7_3"].astype(str)
        return df["CHEBI_PH7_3"].to_dict()
    except:
        return {}


DEFAULT_CHEBI_MAPPING = load_default_chebi_mapping()

MAIN_RELATIONS = [
    "is_conjugate_acid_of",
    "is_conjugate_base_of",
    "is_tautomer_of",
]

CHEBI_PROTONATION_RELATIONS = ["is_conjugate_acid_of", "is_conjugate_base_of"]


def parse_structures():
    """Fixes compatibility issue of libchebipy with information on ftp server."""
    filename = chebi_parsers.get_file("structures.tsv.gz")

    df = pd.read_csv(filename, sep="\t", index_col=0, header=0)
    for _id, row in df.iterrows():
        cpd_id = int(row["compound_id"])
        inchi_key = row["standard_inchi_key"]
        if not pd.isna(inchi_key):
            chebi_parsers.__INCHI_KEYS[cpd_id] = chebi_parsers.Structure(
                inchi_key, chebi_parsers.Structure.InChIKey, 1
            )
        smiles = row["smiles"]
        if not pd.isna(smiles):
            chebi_parsers.__SMILES[cpd_id] = chebi_parsers.Structure(
                smiles, chebi_parsers.Structure.SMILES, 1
            )
    print("Successfully loaded all CHEBI InChI Keys and SMILES structures.")


def get_chebi_for_inchikey(inchikey: str):
    for chebi_id, structure in chebi_parsers.__INCHI_KEYS.items():
        s_obj = structure.get_structure()
        if not isinstance(s_obj, str):
            continue
        if s_obj.startswith(inchikey):
            return chebi_id
    return None


class ChebiEntity(lcpChebiEntity):
    """Class changing some default behavior of ChebiEntity, since the original often fails."""

    # def __init__(self, chebi_id):
    #     super().__init__(
    #         chebi_id,
    #         parser="filesystem",
    #         auto_update=False,
    #         download_dir="/chebi/libChEBI",
    #     )

    def __init__(self, chebi_id):
        super().__init__(chebi_id)


def get_related_chebis(
    chebi: str, related_chebis: Optional[Dict[str, Optional[Tuple[str, str]]]] = None
) -> Dict[str, Optional[Tuple[str, str]]]:
    chebi_entity = ChebiEntity(chebi)

    if related_chebis is None:
        related_chebis = {chebi: None}

    new_chebis = set()

    outgoing_rel = chebi_entity.get_outgoings()
    for rel in outgoing_rel:
        if rel._Relation__typ in MAIN_RELATIONS:
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            if rel_chebi in related_chebis:
                continue
            related_chebis[rel_chebi] = (rel._Relation__typ, chebi)
            new_chebis.add(rel_chebi)

    for new_chebi in new_chebis:
        get_related_chebis(new_chebi, related_chebis)

    return related_chebis


def create_hierarchical_conversion_reaction(
    session: Session, lhs_chebi: ChebiEntity, rhs_chebi: ChebiEntity
) -> Optional[ReferenceReaction]:
    if lhs_chebi.get_charge() != rhs_chebi.get_charge():
        return
    reaction_bigg_id = f"HIER:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
    reference_reaction_db = session.scalars(
        select(ReferenceReaction)
        .filter(ReferenceReaction.bigg_id == reaction_bigg_id)
        .limit(1)
    ).first()
    if reference_reaction_db is not None:
        return
    reference_reaction_db = ReferenceReaction(
        bigg_id=reaction_bigg_id,
        name=f"Conversion of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}, because one is an instance of the other.",
        equation=f"{lhs_chebi.get_name()} = {rhs_chebi.get_name()}",
    )
    lhs_compound = session.scalars(
        select(ReferenceCompound)
        .filter(ReferenceCompound.bigg_id == lhs_chebi.get_id())
        .limit(1)
    ).first()
    lhs_part = ReferenceReactionParticipant(
        compound=lhs_compound,
        side="L",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(lhs_part)
    rhs_compound = session.scalars(
        select(ReferenceCompound)
        .filter(ReferenceCompound.bigg_id == rhs_chebi.get_id())
        .limit(1)
    ).first()
    rhs_part = ReferenceReactionParticipant(
        compound=rhs_compound,
        side="R",
        coefficient="1",
        compartment="0",
    )
    reference_reaction_db.reaction_participants.append(rhs_part)
    reference_reaction_db.update_hash()
    session.add(reference_reaction_db)
    return reference_reaction_db


def add_reference_conversion_reactions(
    session: Session, compound_db: ReferenceCompound
) -> None:
    if not compound_db.bigg_id.startswith("CHEBI:"):
        return
    chebi_entity = ChebiEntity(compound_db.bigg_id)
    main_charge = chebi_entity.get_charge()

    outgoing_rel = chebi_entity.get_outgoings()
    for rel in outgoing_rel:
        if rel._Relation__typ in CHEBI_PROTONATION_RELATIONS:
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = utils.get_object_by_bigg_id(
                session, rel_chebi, ReferenceCompound
            )
            if rel_chebi_db is None:
                continue
            rel_chebi_entity = ChebiEntity(rel_chebi)
            rel_charge = rel_chebi_entity.get_charge()
            if main_charge == rel_charge:
                print("Same charge")
                continue
            elif main_charge < rel_charge:
                lhs_chebi = chebi_entity
                rhs_chebi = rel_chebi_entity
            else:
                lhs_chebi = rel_chebi_entity
                rhs_chebi = chebi_entity
            n_h_plus = rhs_chebi.get_charge() - lhs_chebi.get_charge()
            reaction_bigg_id = f"PROT:{lhs_chebi.get_id()}_{rhs_chebi.get_id()}"
            reference_reaction_db = utils.get_object_by_bigg_id(
                session, reaction_bigg_id, ReferenceReaction
            )
            if reference_reaction_db is not None:
                continue
            reference_reaction_db = ReferenceReaction(
                bigg_id=reaction_bigg_id,
                name=f"Protonation of {lhs_chebi.get_id()} to {rhs_chebi.get_id()}.",
                equation=f"{lhs_chebi.get_name()} + {n_h_plus} H+ = {rhs_chebi.get_name()}",
            )
            lhs_compound = utils.get_object_by_bigg_id(
                session, lhs_chebi.get_id(), ReferenceCompound
            )
            lhs_part = ReferenceReactionParticipant(
                compound=lhs_compound,
                side="L",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(lhs_part)
            rhs_compound = utils.get_object_by_bigg_id(
                session, rhs_chebi.get_id(), ReferenceCompound
            )
            rhs_part = ReferenceReactionParticipant(
                compound=rhs_compound,
                side="R",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(rhs_part)
            proton_compound = utils.get_object_by_bigg_id(
                session, "CHEBI:15378", ReferenceCompound
            )
            proton_part = ReferenceReactionParticipant(
                compound=proton_compound,
                side="L",
                coefficient="1",
                compartment="0",
            )
            reference_reaction_db.reaction_participants.append(proton_part)
            reference_reaction_db.update_hash()
            session.add(reference_reaction_db)
            session.commit()
        elif rel._Relation__typ == "is_a":
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = utils.get_object_by_bigg_id(
                session, rel_chebi, ReferenceCompound
            )
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_outgoings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = utils.get_object_by_bigg_id(
                            session, rel_2_chebi, ReferenceCompound
                        )
                        if rel_chebi_db is None:
                            continue
                        rel_chebi_entity = ChebiEntity(rel_2_chebi)
                        create_hierarchical_conversion_reaction(
                            session, chebi_entity, rel_chebi_entity
                        )
                session.commit()
                continue
            create_hierarchical_conversion_reaction(
                session, chebi_entity, rel_chebi_entity
            )
            session.commit()

    incoming_rel = chebi_entity.get_incomings()
    for rel in incoming_rel:
        if rel._Relation__typ == "is_a":
            rel_chebi = f"CHEBI:{rel._Relation__target_chebi_id}"
            rel_chebi_db = utils.get_object_by_bigg_id(
                session, rel_chebi, ReferenceCompound
            )
            rel_chebi_entity = ChebiEntity(rel_chebi)
            if rel_chebi_db is None:
                # Look for relations where one step is skipped in our DB.
                for rel_2 in rel_chebi_entity.get_incomings():
                    if rel_2._Relation__typ == "is_a":
                        rel_2_chebi = f"CHEBI:{rel_2._Relation__target_chebi_id}"
                        rel_chebi_db = utils.get_object_by_bigg_id(
                            session, rel_2_chebi, ReferenceCompound
                        )
                        if rel_chebi_db is None:
                            continue
                        rel_chebi_entity = ChebiEntity(rel_2_chebi)
                        create_hierarchical_conversion_reaction(
                            session, rel_chebi_entity, chebi_entity
                        )
                session.commit()
                continue
            create_hierarchical_conversion_reaction(
                session, rel_chebi_entity, chebi_entity
            )
            session.commit()
