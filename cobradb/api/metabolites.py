import logging
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar, Union

from sqlalchemy import select
from sqlalchemy.orm import Bundle, Session, contains_eager, subqueryload

from cobradb.data_sources import get_data_source_id
from cobradb.models import (
    Annotation,
    AnnotationLink,
    AnnotationProperty,
    BiGGBase,
    Compartment,
    CompartmentalizedComponent,
    Component,
    ComponentIDMapping,
    ComponentReferenceMapping,
    InChI,
    Model,
    ModelCollection,
    ReferenceCompound,
    ReferenceCompoundAnnotationMapping,
    ReferenceReactivePart,
    ReferenceReactivePartMatrix,
    UniversalCompartmentalizedComponent,
    UniversalComponent,
    UniversalComponentReferenceMapping,
)
from cobradb import chebi
from cobradb.chebi import DEFAULT_CHEBI_MAPPING, ChebiEntity, get_related_chebis
from cobradb.api import utils, reactions
import cobradb.api.bigg_ids as bigg_ids_api

CHEBI_METABOLITE_PROPERTIES = {
    ("CHEBI", "CHEBI"): "CHEBI",
    ("KEGG COMPOUND accession", "KEGG COMPOUND"): "kegg.compound",
    ("DrugBank accession", "DrugBank"): "drugbank",
    ("KEGG DRUG accession", "KEGG DRUG"): "kegg.drug",
    ("Wikipedia accession", "Wikipedia"): "wikipedia.en",
    ("MetaCyc accession", "MetaCyc"): "metacyc.compound",
    ("HMDB accession", "HMDB"): "hmdb",
}


def get_universal_component_by_bigg_id(
    session: Session,
    universal_bigg_id: str,
    model_id: Optional[int] = None,
    model_collection_id: Optional[int] = None,
) -> Optional[UniversalComponent]:
    """Get the univeral component object from a universal BiGG ID. This
    function maps deprecated IDs to the correct new BiGG IDs.

    Arguments
    ---------

    session: SQLAlchemy session

    universal_bigg_id: The BiGG ID to look up
    """

    universal_component_db = session.scalars(
        select(UniversalComponent)
        .join(UniversalComponent.old_bigg_ids)
        .filter(ComponentIDMapping.old_bigg_id == universal_bigg_id)
        .limit(1)
    ).first()
    if universal_component_db is None:
        if model_id is not None:
            model_collection_id = session.scalars(
                select(ModelCollection.id)
                .join(ModelCollection.models)
                .filter(Model.id == model_id)
                .limit(1)
            ).first()
        if model_collection_id is None:
            model_collection_sel = UniversalComponent.collection_id == None
        else:
            model_collection_sel = (UniversalComponent.collection_id == None) | (
                UniversalComponent.collection_id == model_collection_id
            )
        universal_component_db = session.scalars(
            select(UniversalComponent)
            .filter(UniversalComponent.bigg_id == universal_bigg_id)
            .filter(model_collection_sel)
            .limit(1)
        ).first()
    return universal_component_db


def _create_component(
    session: Session,
    bigg_id: str,
    universal_component_db: UniversalComponent,
    reference_compound_db: Optional[ReferenceCompound] = None,
    charge: Optional[Union[int, str]] = None,
    formula: Optional[str] = None,
    variant: Optional[int] = None,
    name: Optional[str] = None,
    collection_db: Optional[ModelCollection] = None,
    require_formula: bool = True,
) -> Optional[Component]:
    if reference_compound_db is not None:
        charge = reference_compound_db.charge
        formula = reference_compound_db.formula
        if name is None:
            name = reference_compound_db.name
    try:
        int_charge = int(str(charge))
    except:
        return None

    if not formula or formula == "nan":
        if require_formula:
            return None
        else:
            formula = ""

    if not ":" in bigg_id:
        bigg_id = bigg_ids_api.create_component_bigg_id(bigg_id, charge=int_charge)

    component_db = Component(
        bigg_id=bigg_id,
        name=name,
        formula=formula,
        charge=int_charge,
        collection=collection_db,
        variant=variant,
    )
    universal_component_db.components.append(component_db)
    return component_db


def _create_universal_component(
    session: Session,
    bigg_id: str,
    name: Optional[str] = None,
    collection_db: Optional[ModelCollection] = None,
    default_component: Optional[Component] = None,
    allow_flexible_variants: bool = False,
):
    universal_component_db = UniversalComponent(
        bigg_id=bigg_id,
        name=name,
        collection=collection_db,
        default_component=default_component,
        allow_flexible_variants=allow_flexible_variants,
    )
    session.add(universal_component_db)
    return universal_component_db


def _create_inchi_object(
    session: Session,
    inchi: str,
) -> Optional[InChI]:
    """Create an InChI object from an InChI string if possible, or retrieve the corresponding object from the database.

    Arguments
    ---------

    session: SQLAlchemy session

    inchi: InChI string
    """
    if not inchi:
        return None

    inchi_db = None
    inchi_obj = InChI.from_string(inchi)
    if inchi_obj is not None:
        inchi_db = session.scalars(
            select(InChI)
            .filter(
                (InChI.key_major == inchi_obj.key_major)
                & (InChI.key_minor == inchi_obj.key_minor)
                & (InChI.key_proton == inchi_obj.key_proton)
            )
            .limit(1)
        ).first()
        if inchi_obj != inchi_db:
            inchi_db = None
        if inchi_db is None:
            inchi_db = inchi_obj
    return inchi_db


SMRT = TypeVar("SMRT", bound=BiGGBase)


def get_or_create_small_molecule_reference(
    session: Session,
    chebi: str,
    cpd_cls: Type[SMRT] = ReferenceCompound,
) -> Tuple[bool, Optional[SMRT]]:
    """Create a small molecule reference based on a ChEBI ID. Can either create a ReferenceCompound or ReferenceReactivePart.

    Arguments
    ---------
    session: SQLAlchemy session

    chebi: ChEBI ID as a string (CHEBI:00000 format)

    cpd_cls: Class of reference object to create. Default: ReferenceCompound."""
    chebi_db = utils.get_object_by_bigg_id(session, chebi, cpd_cls)
    if chebi_db:
        return True, chebi_db

    chebi_entity = ChebiEntity(chebi)
    if not chebi_entity:
        return False, None

    inchi_db = _create_inchi_object(session, chebi_entity.get_inchi())

    cpd_dict = dict(
        bigg_id=chebi,
        name=chebi_entity.get_name(),
        html_name=chebi_entity.get_name(),
        charge=str(chebi_entity.get_charge()),
        formula=chebi_entity.get_formula(),
        inchi=inchi_db,
    )
    if cpd_cls is ReferenceCompound:
        cpd_dict["compound_type"] = "small_molecule"

    chebi_db = cpd_cls(**cpd_dict)
    session.add(chebi_db)

    if isinstance(
        chebi_db, ReferenceCompound
    ):  # Typechecker is confused when using cpd_cls.
        add_chebi_annotations(session, chebi_db, chebi_entity)

    return False, chebi_db


def create_small_molecule_reference_for_inchi(
    session: Session,
    ref_id: str,
    inchi: InChI,
    name: Optional[str],
    cpd_cls: Type[SMRT] = ReferenceCompound,
) -> SMRT:
    """Create a small molecule reference based on a ChEBI ID. Can either create a ReferenceCompound or ReferenceReactivePart.

    Arguments
    ---------
    session: SQLAlchemy session

    chebi: ChEBI ID as a string (CHEBI:00000 format)

    cpd_cls: Class of reference object to create. Default: ReferenceCompound."""

    n_protons = inchi.n_protons()
    if n_protons == 0:
        formula = inchi.formula
    else:
        delta_h = utils.Formula({"H": n_protons})
        formula = str(utils.Formula(inchi.formula) + delta_h)

    cpd_dict = dict(
        bigg_id=ref_id,
        name=name,
        html_name=name,
        charge=str(inchi.charge()),
        formula=formula,
        inchi=inchi,
    )
    if cpd_cls is ReferenceCompound:
        cpd_dict["compound_type"] = "small_molecule"

    ref_db = cpd_cls(**cpd_dict)
    session.add(ref_db)

    return ref_db


def add_chebi_annotations(
    session: Session,
    reference_compound_db: ReferenceCompound,
    chebi_entity: ChebiEntity,
) -> None:
    """Add annotations from ChEBI to the reference compound.

    Arguments
    ---------
    session: SQLAlchemy session

    reference_compound_db: ReferenceCompound object to annotate.

    chebi_entity: ChebiEntity object to retrieve annotations from."""
    chebi = reference_compound_db.bigg_id
    annotation_db = utils.get_object_by_bigg_id(
        session,
        chebi,
        Annotation,
        opts=(subqueryload(Annotation.links), subqueryload(Annotation.properties)),
    )
    if not annotation_db:
        annotation_db = Annotation(
            bigg_id=chebi,
            type="chebi",
            default_data_source_id=get_data_source_id("CHEBI", session),
        )
        session.add(annotation_db)

    for property, func in {
        "definition": chebi_entity.get_definition,
        "star": chebi_entity.get_star,
        "mass": chebi_entity.get_mass,
        "smiles": chebi_entity.get_smiles,
    }.items():
        prop_val = func()
        if isinstance(prop_val, str):
            prop_val = prop_val.strip()
            if len(prop_val) == 0:
                prop_val = None
        if prop_val is None:
            continue
        if isinstance(prop_val, int) and property.startswith("is_"):
            prop_val = bool(prop_val)
        prop_db = AnnotationProperty(key=property)
        prop_db.value = prop_val
        annotation_db.properties.append(prop_db)

    for name in chebi_entity.get_names():
        prop_db = AnnotationProperty(key="name")
        prop_db.value = name.get_name()
        annotation_db.properties.append(prop_db)

    alias_db = AnnotationLink(
        identifier=chebi.removeprefix("CHEBI:"),
        data_source_id=get_data_source_id("CHEBI", session),
    )
    annotation_db.links.append(alias_db)

    for database_accession in chebi_entity.get_database_accessions():
        data_source_bigg_id = CHEBI_METABOLITE_PROPERTIES.get(
            (database_accession.get_type(), database_accession.get_source())
        )
        if data_source_bigg_id is None:
            continue
        prop_val = database_accession.get_accession_number()
        if isinstance(prop_val, str):
            prop_val = prop_val.strip()
            if len(prop_val) == 0:
                prop_val = None
        if prop_val is None:
            continue
        if isinstance(prop_val, str):
            prop_val = {prop_val}

        for val in prop_val:
            alias_db = AnnotationLink(
                identifier=val,
                data_source_id=get_data_source_id(data_source_bigg_id, session),
            )
            annotation_db.links.append(alias_db)

    annotation_mapping = ReferenceCompoundAnnotationMapping(
        reference_compound=reference_compound_db,
    )
    annotation_db.reference_compound_mappings.append(annotation_mapping)


def _create_reactive_part(
    session: Session,
    bigg_id: str,
    name: str,
    formula: str,
    charge: str,
    html_name: Optional[str] = None,
) -> Optional[ReferenceReactivePart]:
    if bigg_id.startswith("CHEBI:"):
        _existed, reactive_part_db = get_or_create_small_molecule_reference(
            session, bigg_id, cpd_cls=ReferenceReactivePart
        )
    else:
        reactive_part_db = ReferenceReactivePart(
            bigg_id=bigg_id,
            name=name,
            html_name=html_name,
            formula=formula,
            charge=charge,
        )
        session.add(reactive_part_db)
    return reactive_part_db


def _create_reference_compound(
    session: Session,
    bigg_id: str,
    name: str,
    compound_type: str,
    formula: Optional[str] = None,
    charge: Optional[str] = None,
    html_name: Optional[str] = None,
    reactive_parts: Optional[List[str]] = None,
):
    if bigg_id.startswith("CHEBI:"):
        _existed, compound_db = get_or_create_small_molecule_reference(session, bigg_id)
    else:
        compound_db = ReferenceCompound(
            bigg_id=bigg_id,
            name=name,
            html_name=html_name,
            formula=formula,
            charge=charge,
            compound_type=compound_type,
        )
        reactive_parts_formula_sum = utils.Formula({})
        reactive_parts_charge_sum = 0
        if reactive_parts is not None:
            for reactive_part in reactive_parts:
                reactive_part_db = utils.get_object_by_bigg_id(
                    session, reactive_part, ReferenceReactivePart
                )
                if reactive_part_db is None:
                    return
                reactive_parts_formula_sum += reactive_part_db.formula
                reactive_parts_charge_sum += (
                    0
                    if reactive_part_db.charge is None
                    else int(reactive_part_db.charge)
                )
                reactive_part_matrix_db = ReferenceReactivePartMatrix(
                    reactive_part=reactive_part_db,
                )
                compound_db.reactive_part_matrix.append(reactive_part_matrix_db)
        if compound_db.formula is None or compound_db.charge is None:
            compound_db.formula = str(reactive_parts_formula_sum)
            compound_db.charge = str(reactive_parts_charge_sum)
        session.add(compound_db)
    session.commit()
    if compound_db is not None:
        chebi.add_reference_conversion_reactions(session, compound_db)
        reactions.add_reference_exchange_reactions(session, compound_db)
    return compound_db


def create_metabolite(
    session: Session,
    proposed_bigg_id: str,
    input_chebi: str,
    assure_present: Optional[Tuple[Any, Any]] = None,
    reference_n: Optional[int] = None,
):
    if not bigg_ids_api.is_bigg_id_valid(
        proposed_bigg_id, bigg_ids_api.BiGGIDType.UNIVERSAL_COMPONENT
    ):
        return {"status": "error", "message": "invalid bigg id"}

    universal_component_db = get_universal_component_by_bigg_id(
        session, proposed_bigg_id
    )
    if universal_component_db is not None:
        if universal_component_db.bigg_id == proposed_bigg_id:
            return {"status": "error", "message": "bigg id already exists"}
        else:
            return {
                "status": "error",
                "message": "bigg id is a deprecated id",
                "new_id": universal_component_db.bigg_id,
            }

    if reference_n is None:
        all_chebis = get_related_chebis(input_chebi)
        all_chebis = list({input_chebi} | set(all_chebis.keys()))
    else:
        all_chebis = [input_chebi]

    component_ref_mapping_db = session.execute(
        select(
            ComponentReferenceMapping,
            Bundle("ref_comp", ReferenceCompound.bigg_id),
            Bundle("uni_comp", UniversalComponent.bigg_id),
        )
        .join(ComponentReferenceMapping.reference_compound)
        .join(ComponentReferenceMapping.component)
        .join(Component.universal_component)
        .filter(ReferenceCompound.bigg_id.in_(all_chebis))
        .filter(ComponentReferenceMapping.reference_n == reference_n)
        .limit(1)
    ).one_or_none()
    if component_ref_mapping_db:
        return {
            "status": "error",
            "message": "chebi already associated",
            "chebi": component_ref_mapping_db.ref_comp.bigg_id,
            "bigg_id": component_ref_mapping_db.uni_comp.bigg_id,
        }

    references_db = {
        x: ref
        for x in all_chebis
        if (ref := get_or_create_small_molecule_reference(session, x)[1]) is not None
    }

    if assure_present:
        assure_charge, assure_formula = assure_present
        assure_f = utils.Formula(assure_formula)

        if reference_n is None:
            comp = any(
                (float(str(ref.charge)) == assure_charge)
                and (utils.Formula(ref.formula) == assure_f)
                for ref in references_db.values()
            )
        else:
            comp = any(
                (utils.NCharge(ref.charge).fill(reference_n) == assure_charge)
                and (utils.NFormula(ref.formula).fill(reference_n) == assure_f)
                for ref in references_db.values()
            )
        if not comp:
            return {
                "status": "error",
                "message": "charge + formula combination not present",
            }

    default_chebi = None
    for ch in all_chebis:
        default_chebi = DEFAULT_CHEBI_MAPPING.get(ch)
        if default_chebi is not None:
            break
    if default_chebi not in references_db:
        default_chebi = None
    if default_chebi is None:
        default_chebi = input_chebi
    default_chebi_entity = ChebiEntity(default_chebi)

    universal_component_db = _create_universal_component(
        session,
        bigg_id=proposed_bigg_id,
        name=default_chebi_entity.get_name(),
    )
    session.commit()

    successfully_added = []

    full_bids_created = {}
    for chebi, reference_db in references_db.items():
        if reference_n is None:
            explicit_charge = reference_db.charge
            explicit_formula = reference_db.formula
        else:
            explicit_charge = utils.NCharge(reference_db.charge).fill(reference_n)
            explicit_formula = str(
                utils.NFormula(reference_db.formula).fill(reference_n)
            )
        full_bid = bigg_ids_api.create_component_bigg_id(
            proposed_bigg_id, charge=explicit_charge
        )
        if full_bid not in full_bids_created:
            if (
                component_db := _create_component(
                    session,
                    full_bid,
                    universal_component_db,
                    name=reference_db.name,
                    charge=explicit_charge,
                    formula=explicit_formula,
                )
            ) is None:
                continue
            full_bids_created[full_bid] = component_db
        component_db = full_bids_created[full_bid]

        component_ref_mapping_db = ComponentReferenceMapping(
            universal_component=universal_component_db,
            reference_compound=reference_db,
            reference_n=reference_n,
        )
        component_db.reference_mappings.append(component_ref_mapping_db)

        if chebi == default_chebi or default_chebi is None:
            universal_component_ref_mapping_db = session.scalars(
                select(UniversalComponentReferenceMapping)
                .filter(
                    UniversalComponentReferenceMapping.id == universal_component_db.id
                )
                .limit(1)
            ).first()
            if not universal_component_ref_mapping_db:
                universal_component_ref_mapping_db = UniversalComponentReferenceMapping(
                    id=universal_component_db.id, mapping=component_ref_mapping_db
                )
                session.add(universal_component_ref_mapping_db)

            universal_component_db.default_component = component_db

        successfully_added.append((full_bid, reference_db.bigg_id))
        session.commit()

    if not full_bids_created:
        session.commit()
        session.delete(universal_component_db)
    session.commit()

    return {"status": "success", "components_added": successfully_added}


def create_complex_metabolite(
    session: Session,
    proposed_bigg_id: str,
    ref_identifier: str,
):
    if not bigg_ids_api.is_bigg_id_valid(
        proposed_bigg_id, bigg_ids_api.BiGGIDType.UNIVERSAL_COMPONENT
    ):
        return {"status": "error", "message": "invalid bigg id"}

    universal_component_db = get_universal_component_by_bigg_id(
        session, proposed_bigg_id
    )
    if universal_component_db is not None:
        if universal_component_db.bigg_id == proposed_bigg_id:
            return {"status": "error", "message": "bigg id already exists"}
        else:
            return {
                "status": "error",
                "message": "bigg id is a deprecated id",
                "new_id": universal_component_db.bigg_id,
            }

    component_ref_mapping_db = session.execute(
        select(
            ComponentReferenceMapping,
            Bundle("ref_comp", ReferenceCompound.bigg_id),
            Bundle("uni_comp", UniversalComponent.bigg_id),
        )
        .join(ComponentReferenceMapping.reference_compound)
        .join(ComponentReferenceMapping.component)
        .join(Component.universal_component)
        .filter(ReferenceCompound.bigg_id == ref_identifier)
        .limit(1)
    ).one_or_none()
    if component_ref_mapping_db:
        return {
            "status": "error",
            "message": "reference identifier already associated",
            "chebi": component_ref_mapping_db.ref_comp.bigg_id,
            "bigg_id": component_ref_mapping_db.uni_comp.bigg_id,
        }

    reference_db = session.scalars(
        select(ReferenceCompound)
        .filter(ReferenceCompound.bigg_id == ref_identifier)
        .limit(1)
    ).first()

    if reference_db is None:
        return {"status": "error", "message": "Reference could not be found"}

    universal_component_db = _create_universal_component(
        session,
        bigg_id=proposed_bigg_id,
        name=reference_db.name,
        allow_flexible_variants=True,
    )
    session.commit()

    charge = 0 if reference_db.charge is None else int(reference_db.charge)
    full_bid = bigg_ids_api.create_component_bigg_id(proposed_bigg_id, charge=charge)
    component_db = _create_component(
        session, full_bid, universal_component_db, reference_db
    )
    if component_db is None:
        return {
            "status": "error",
            "message": "could not create component",
        }
    component_ref_mapping_db = ComponentReferenceMapping(
        universal_component=universal_component_db,
        reference_compound=reference_db,
        reference_formula_delta="",
    )
    component_db.reference_mappings.append(component_ref_mapping_db)
    universal_component_ref_mapping_db = UniversalComponentReferenceMapping(
        id=universal_component_db.id, mapping=component_ref_mapping_db
    )
    session.add(universal_component_ref_mapping_db)

    universal_component_db.default_component = component_db
    session.commit()

    return {"status": "success", "components_added": [full_bid]}


def create_metabolites_for_inchis(
    session: Session, bid: str, inchis: List[Dict[str, str]]
):
    inchi_objects = []
    for inchi_d in inchis:
        inchi = inchi_d["string"]
        inchi_obj = InChI.from_string(inchi)
        if inchi_obj is None:
            logging.error(f"Could not parse InChI string: '{inchi}'")
            continue
        possible_matches = session.scalars(
            select(InChI)
            .filter(InChI.key_major == inchi_obj.key_major)
            .filter(InChI.key_minor == inchi_obj.key_minor)
            .filter(InChI.key_proton == inchi_obj.key_proton)
        ).all()
        if any(x == inchi_obj for x in possible_matches):
            logging.error(f"InChI is already in database: '{inchi}'")
            continue
        inchi_objects.append((inchi_obj, inchi_d.get("name")))
    existing_components = session.scalars(
        select(Component)
        .options(contains_eager(Component.universal_component))
        .join(Component.universal_component)
        .filter(UniversalComponent.bigg_id == bid)
    ).all()
    common_key_major = None
    common_key_minor = None
    existing_charges = set()
    for existing_component in existing_components:
        if len(existing_component.reference_mappings) == 0:
            logging.error(
                f"Existing component '{existing_component.bigg_id}' does not have a reference"
            )
            return
        for ref_map in existing_component.reference_mappings:
            if ref_map.reference_compound.inchi is None:
                logging.error(
                    f"Existing component '{existing_component.bigg_id}' with reference {ref_map.reference_compound.bigg_id} does not have an InChI associated."
                )
                return
            if common_key_major is None:
                common_key_major = ref_map.reference_compound.inchi.key_major
                common_key_minor = ref_map.reference_compound.inchi.key_minor
            else:
                if (
                    common_key_major != ref_map.reference_compound.inchi.key_major
                    or common_key_minor != ref_map.reference_compound.inchi.key_minor
                ):
                    logging.error(
                        f"Existing component '{existing_component.bigg_id}' with reference {ref_map.reference_compound.bigg_id} has conflicting InChI associated."
                    )
                    # return
            new_charge = ref_map.reference_compound.inchi.charge()
            existing_charges.add(new_charge)

    universal_component_db = None
    if len(existing_components) > 0:
        universal_component_db = existing_components[0].universal_component
    for inchi_obj, name in inchi_objects:
        if common_key_major is None:
            common_key_major = inchi_obj.key_major
            common_key_minor = inchi_obj.key_minor
        else:
            if (
                common_key_major != inchi_obj.key_major
                or common_key_minor != inchi_obj.key_minor
            ):
                logging.error(
                    f"Provided InChI {inchi_obj.to_string()} conflicts with existing InChI references. Creating reference anyway."
                )
                # return
        new_charge = inchi_obj.charge()
        if new_charge in existing_charges:
            logging.error(
                f"An InChI reference with charge {new_charge} already exists for this molecule."
            )
            return

        if universal_component_db is None:
            universal_component_db = _create_universal_component(
                session,
                bigg_id=bid,
                name=name,
            )
            session.commit()

        full_bid = bigg_ids_api.create_component_bigg_id(bid, charge=new_charge)
        reference_bigg_id = f"BIGGR:{full_bid}"
        reference_db = create_small_molecule_reference_for_inchi(
            session, reference_bigg_id, inchi_obj, name
        )
        component_db = _create_component(
            session, full_bid, universal_component_db, reference_db
        )
        session.commit()
        if component_db is None:
            continue

        component_ref_mapping_db = ComponentReferenceMapping(
            universal_component=universal_component_db,
            reference_compound=reference_db,
        )
        component_db.reference_mappings.append(component_ref_mapping_db)

        if len(universal_component_db.components) == 1:
            universal_component_ref_mapping_db = UniversalComponentReferenceMapping(
                id=universal_component_db.id, mapping=component_ref_mapping_db
            )
            session.add(universal_component_ref_mapping_db)
            universal_component_db.default_component = component_db
        session.commit()

        existing_charges.add(new_charge)


def create_collection_specific_metabolite(
    session: Session,
    bigg_id: str,
    collection_db: ModelCollection,
    charge: Optional[Any],
    name: str,
    formula: str,
) -> Tuple[Optional[UniversalComponent], Optional[Component]]:
    if charge is None:
        charge = 0
    new_universal_id = bigg_ids_api.create_component_bigg_id(
        bigg_id,
        collection_bigg_id=collection_db.bigg_id,
    )
    new_bigg_id = bigg_ids_api.create_component_bigg_id(
        bigg_id,
        charge=charge,
        collection_bigg_id=collection_db.bigg_id,
    )

    if (
        universal_component_db := utils.get_object_by_bigg_id(
            session, new_universal_id, UniversalComponent
        )
    ) is None:
        universal_component_db = _create_universal_component(
            session,
            bigg_id=new_universal_id,
            name=name,
            collection_db=collection_db,
        )

    if (
        component_db := utils.get_object_by_bigg_id(session, new_bigg_id, Component)
    ) is None:
        component_db = _create_component(
            session,
            bigg_id=new_bigg_id,
            universal_component_db=universal_component_db,
            name=name,
            formula=formula,
            charge=charge,
            collection_db=collection_db,
            require_formula=False,
        )
        if universal_component_db.default_component is None:
            universal_component_db.default_component = component_db
    return universal_component_db, component_db


def create_component_id_mapping(
    session: Session, old_bigg_id: str, new_bigg_id: str
) -> Optional[ComponentIDMapping]:
    id_mapping_db = session.scalars(
        select(ComponentIDMapping)
        .filter(ComponentIDMapping.old_bigg_id == old_bigg_id)
        .limit(1)
    ).first()
    if id_mapping_db is not None:
        return None

    # Check that we're not making an entry inaccessible.
    if (
        utils.get_object_by_bigg_id(session, old_bigg_id, UniversalComponent)
        is not None
    ):
        return None

    if (
        new_component_db := utils.get_object_by_bigg_id(
            session, new_bigg_id, UniversalComponent
        )
    ) is None:
        return None
    id_mapping_db = ComponentIDMapping(
        old_bigg_id=old_bigg_id,
        new_universal_component=new_component_db,
    )
    session.add(id_mapping_db)
    return id_mapping_db


def get_or_create_universal_compartmentalized_component(
    session: Session,
    universal_component_db: UniversalComponent,
    compartment_db: Compartment,
) -> UniversalCompartmentalizedComponent:
    universal_compartmentalized_component_db = session.scalars(
        select(UniversalCompartmentalizedComponent)
        .filter(UniversalCompartmentalizedComponent.compartment_id == compartment_db.id)
        .filter(
            UniversalCompartmentalizedComponent.universal_component
            == universal_component_db
        )
    ).first()
    if universal_compartmentalized_component_db is None:
        ucc_bigg_id = bigg_ids_api.create_component_bigg_id(
            base_bigg_id=universal_component_db.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
        )
        universal_compartmentalized_component_db = UniversalCompartmentalizedComponent(
            bigg_id=ucc_bigg_id,
            universal_component=universal_component_db,
            compartment=compartment_db,
        )
        session.add(universal_compartmentalized_component_db)
    return universal_compartmentalized_component_db


def get_or_create_compartmentalized_component(
    session: Session,
    component_db: Component,
    compartment_db: Compartment,
) -> CompartmentalizedComponent:
    compartmentalized_component_db = session.scalars(
        select(CompartmentalizedComponent)
        .filter(CompartmentalizedComponent.compartment_id == compartment_db.id)
        .filter(CompartmentalizedComponent.component == component_db)
    ).first()
    if compartmentalized_component_db is None:
        cc_bigg_id = bigg_ids_api.create_component_bigg_id(
            base_bigg_id=component_db.universal_component.bigg_id,
            compartment_bigg_id=compartment_db.bigg_id,
            charge=component_db.charge,
            variant=component_db.variant,
        )
        compartmentalized_component_db = CompartmentalizedComponent(
            bigg_id=cc_bigg_id,
            component=component_db,
            compartment=compartment_db,
        )
        session.add(compartmentalized_component_db)
    return compartmentalized_component_db
