from collections.abc import Mapping
import gzip
import json
import re
import cobra
from sqlalchemy import select
from cobradb.api import utils
import memote
import memote.suite.api as api
from memote.suite.results import ResultManager

memote_version = getattr(memote, "__version__", None)

from cobradb.models import (
    AlreadyLoadedError,
    CompartmentalizedComponent,
    MemoteResult,
    MemoteRun,
    MemoteTest,
    Model,
    ModelCompartmentalizedComponent,
    ModelReaction,
    Reaction,
    UniversalCompartmentalizedComponent,
    UniversalReaction,
)

TEST_TYPES = {
    "test_absolute_extreme_coefficient_ratio": "general",
    "test_biomass_consistency": "per_reaction",
    "test_biomass_default_production": "per_reaction",
    "test_biomass_open_production": "per_reaction",
    # "test_biomass_precursors_default_production": "per_reaction_count_metabolite",
    # "test_biomass_precursors_open_production": "per_reaction_count_metabolite",
    "test_biomass_presence": "count_reaction",
    "test_biomass_specific_sbo_presence": {"type": "general", "invert": True},
    "test_blocked_reactions": "count_reaction",
    # "test_compartments_presence": "compartment_count",
    "test_degrees_of_freedom": "general",
    "test_demand_specific_sbo_presence": {"type": "count_reaction", "invert": True},
    # "test_detect_energy_generating_cycles": "TODO",
    # "test_direct_metabolites_in_biomass": "per_reaction_count_metabolite",
    # "test_essential_precursors_not_in_biomass": "per_reaction_count_metabolite_MNX?",
    "test_exchange_specific_sbo_presence": {
        "type": "count_reaction",
        "invert": True,
        "invert_result": True,
    },
    "test_fast_growth_default": "per_reaction",
    "test_fbc_presence": "general",
    # "test_find_constrained_pure_metabolic_reactions": "TODO",
    "test_find_constrained_transport_reactions": {
        "type": "count_reaction",
        "field": "data_count",
    },
    "test_find_deadends": "count_metabolite",
    "test_find_disconnected": {"type": "count_metabolite", "invert": True},
    "test_find_duplicate_metabolites_in_compartments": {
        "type": "count_metabolite",
        "field": "data_count",
    },
    "test_find_duplicate_reactions": {"type": "count_reaction"},
    "test_find_medium_metabolites": {
        "type": "count_metabolite",
        "field": "data_count",
        "invert_result": True,
    },
    "test_find_metabolites_not_consumed_with_open_bounds": "count_metabolite",
    "test_find_metabolites_not_produced_with_open_bounds": "count_metabolite",
    "test_find_orphans": "count_metabolite",
    "test_find_pure_metabolic_reactions": {
        "type": "count_reaction",
        "field": "data_count",
    },
    "test_find_reactions_unbounded_flux_default_condition": {
        "type": "count_reaction",
        "invert": True,
    },
    # "test_find_reactions_with_identical_genes": "TODO",
    # "test_find_reactions_with_partially_identical_annotations": "TODO",
    "test_find_reversible_oxygen_reactions": "count_reaction",
    "test_find_stoichiometrically_balanced_cycles": "count_reaction",
    "test_find_transport_reactions": {"type": "count_reaction", "field": "data_count"},
    # "test_find_unique_metabolites": "count_metabolite", (NO compartment)
    "test_gam_in_biomass": "per_reaction",
    # "test_gene_essentiality_from_data_qualitative": "TODO",
    # "test_gene_product_annotation_overview": "TODO",
    # "test_gene_product_annotation_presence": "TODO",
    # "test_gene_product_annotation_wrong_ids": "TODO",
    "test_gene_protein_reaction_rule_presence": "count_reaction",
    # "test_gene_sbo_presence": "count_gene",
    # "test_gene_specific_sbo_presence": "count_gene",
    # "test_genes_presence": "count_gene",
    # "test_growth_from_data_qualitative": "TODO",
    "test_inconsistent_min_stoichiometry": "count_metabolite",
    "test_matrix_rank": "general",
    # "test_metabolic_coverage": "TODO",
    "test_metabolic_reaction_specific_sbo_presence": {
        "type": "count_reaction",
        "invert": True,
        "invert_result": True,
    },
    # "test_metabolite_annotation_overview": "TODO",
    "test_metabolite_annotation_presence": "count_metabolite",
    # "test_metabolite_annotation_wrong_ids": "TODO",
    "test_metabolite_id_namespace_consistency": "count_metabolite",
    "test_metabolite_sbo_presence": {
        "type": "count_metabolite",
        "invert": True,
        "invert_result": True,
    },
    "test_metabolite_specific_sbo_presence": {
        "type": "count_metabolite",
        "invert": True,
        "invert_result": True,
    },
    "test_metabolites_charge_presence": {
        "type": "count_metabolite",
        "field": "data_count",
    },
    "test_metabolites_formula_presence": {
        "type": "count_metabolite",
        "field": "data_count",
    },
    "test_metabolites_presence": {"type": "count_metabolite", "field": "data_count"},
    # "test_model_id_presence": "general",
    "test_ngam_presence": "count_reaction",
    "test_number_independent_conservation_relations": "general",
    "test_protein_complex_presence": "count_reaction",
    # "test_reaction_annotation_overview": "TODO",
    "test_reaction_annotation_presence": "count_reaction",
    # "test_reaction_annotation_wrong_ids": "TODO",
    "test_reaction_charge_balance": {"type": "count_reaction", "invert": True},
    "test_reaction_id_namespace_consistency": "count_reaction",
    "test_reaction_mass_balance": {"type": "count_reaction", "invert": True},
    "test_reaction_sbo_presence": {
        "type": "count_reaction",
        "invert": True,
        "invert_result": True,
    },
    "test_reactions_presence": {"type": "count_reaction", "field": "data_count"},
    # "test_sbml_level": "general (STRING)",
    "test_sink_specific_sbo_presence": {"type": "count_reaction", "invert": True},
    "test_stoichiometric_consistency": {"type": "general", "invert": True},
    "test_transport_reaction_gpr_presence": "count_reaction",
    "test_transport_reaction_specific_sbo_presence": {
        "type": "count_reaction",
        "invert": True,
        "invert_result": True,
    },
    "test_unconserved_metabolites": {"type": "count_metabolite", "field": "data_count"},
}

MEMOTE_REACTION_REGEX = re.compile(
    r"^([0-9\.]+ \<\= )?(?P<rname>[_a-zA-Z0-9\:]+)( \<\= [0-9\.]+)?$"
)


SOLVER = "glpk"


def record_memote_run(session, model_db, results, skip=None):
    """Record which MEMOTE version produced this model's results.

    MEMOTE does not write its own version into the result file -- meta.packages
    lists cobrapy and the rest but leaves "memote" as None -- so read it off the
    installed package instead.

    Updates the row in place when one exists, so re-running a model does not
    leave a stale version behind.
    """
    cobrapy_version = results.get("meta", {}).get("packages", {}).get("cobra")
    skipped_tests = ", ".join(sorted(skip)) if skip else None

    run_db = session.scalars(
        select(MemoteRun).filter(MemoteRun.model_id == model_db.id).limit(1)
    ).first()
    if run_db is None:
        run_db = MemoteRun(model=model_db)
        session.add(run_db)

    run_db.memote_version = memote_version
    run_db.cobrapy_version = cobrapy_version or cobra.__version__
    run_db.solver = SOLVER
    run_db.skipped_tests = skipped_tests


def run_memote(model_filename, result_filename, skip=None):
    """Run the MEMOTE suite and store the result file.

    skip names tests to leave out. GLPK aborts the process at the C level on
    the largest models (RECON1 dies in test_inconsistent_min_stoichiometry),
    and since that is not a Python exception it cannot be caught -- the only
    way to get a report for such a model is to not run the test. Whatever is
    skipped is recorded by load_memote_results so a missing test is not
    mistaken for a failing one.
    """
    config = cobra.Configuration()
    config.solver = SOLVER

    # Check if the model can be loaded at all.
    model, sbml_ver, notifications = api.validate_model(model_filename)
    if model is None:
        raise Exception(
            "The model could not be loaded due to the following SBML errors."
        )

    code, result = api.test_model(
        model=model,
        sbml_version=sbml_ver,
        results=True,
        skip=skip,
        solver_timeout=10,
    )

    manager = ResultManager()
    manager.store(result, filename=result_filename)


def load_memote_results(model_bigg_id, filename, session, skip=None):
    print(f"Loading memote results for {model_bigg_id}")

    with gzip.open(filename, "rb") as f:
        results = json.load(f)

    model_db = utils.get_object_by_bigg_id(session, model_bigg_id, Model)
    if model_db is None:
        raise ValueError(f"Model {model_bigg_id} could not be found.")
    model_db_id = model_db.id

    existing_result_db = session.scalars(
        select(MemoteResult).filter(MemoteResult.model == model_db).limit(1)
    ).first()
    if existing_result_db is not None:
        raise AlreadyLoadedError(
            f"There are already Memote results for {model_bigg_id} in the database",
            bigg_id=model_bigg_id,
        )

    record_memote_run(session, model_db, results, skip=skip)

    for test_func, test_result in results["tests"].items():
        test_bigg_id = test_func
        test_type = TEST_TYPES.get(test_bigg_id)
        if test_type is None:
            continue

        invert = False
        invert_result = False
        field = "metric"
        if isinstance(test_type, dict):
            field = test_type.get("field", field)
            invert = test_type.get("invert", invert)
            invert_result = test_type.get("invert_result", invert_result)
            test_type = test_type["type"]

        test_db = session.scalars(
            select(MemoteTest).filter(MemoteTest.bigg_id == test_bigg_id).limit(1)
        ).first()
        if test_db is None:
            test_db = MemoteTest(
                bigg_id=test_bigg_id,
                name=test_result["title"],
                summary=test_result.get("summary"),
                format_type=test_result["format_type"],
                invert=invert,
                invert_result=invert_result,
                field=field,
            )
            session.add(test_db)
            session.commit()

        general_result = {
            "model_id": model_db_id,
            "test_id": test_db.id,
        }
        specific_results = {}

        for prop in ["data", "duration", "metric", "result", "message"]:
            prop_val = test_result.get(prop)
            if prop_val is None:
                continue
            if isinstance(prop_val, (int, float)):
                general_result[prop] = float(prop_val)
            elif isinstance(prop_val, str):
                general_result[prop] = prop_val
            elif isinstance(prop_val, bool):
                general_result[prop] = 1.0 if prop_val else 0.0
            elif (is_dict := isinstance(prop_val, dict)) or isinstance(prop_val, list):
                if prop == "data":
                    general_result["data_count"] = len(prop_val)
                for k in prop_val:
                    if isinstance(k, list) and len(k) == 1:
                        k = k[0]
                    if not isinstance(k, str):
                        print(f"Key {k} not a string.")
                        continue
                    if k not in specific_results:
                        specific_results[k] = {
                            "model_id": model_db_id,
                            "test_id": test_db.id,
                        }
                    if test_type in ["per_reaction", "count_reaction"]:
                        reaction_match = MEMOTE_REACTION_REGEX.match(k)
                        if not reaction_match:
                            universal_reaction_id = k
                        else:
                            universal_reaction_id = reaction_match.group("rname")
                        # copy_number = 1
                        # if ":" in universal_reaction_id:
                        #     universal_reaction_id, copy_number = (
                        #         universal_reaction_id.rsplit(":", maxsplit=1)
                        #     )
                        #     copy_number = int(copy_number)
                        model_reaction_db = session.scalars(
                            select(ModelReaction)
                            .join(Reaction.universal_reaction)
                            .filter(
                                (ModelReaction.bigg_id == universal_reaction_id)
                                & (ModelReaction.model_id == model_db_id)
                            )
                            .limit(1)
                        ).first()
                        if model_reaction_db is None:
                            print(f"Reaction not found: {k} ({universal_reaction_id})")
                            continue
                        specific_results[k]["model_reaction_id"] = model_reaction_db.id
                    elif test_type in ["per_metabolite", "count_metabolite"]:
                        metabolite_id = k
                        if ":" in metabolite_id:
                            model_comp_comp_db = session.scalars(
                                select(ModelCompartmentalizedComponent)
                                .join(
                                    ModelCompartmentalizedComponent.compartmentalized_component
                                )
                                .filter(
                                    (
                                        CompartmentalizedComponent.bigg_id
                                        == metabolite_id
                                    )
                                    & (
                                        ModelCompartmentalizedComponent.model_id
                                        == model_db_id
                                    )
                                )
                                .limit(1)
                            ).first()
                        else:
                            model_comp_comp_db = session.scalars(
                                select(ModelCompartmentalizedComponent)
                                .join(
                                    ModelCompartmentalizedComponent.compartmentalized_component
                                )
                                .join(
                                    CompartmentalizedComponent.universal_compartmentalized_component
                                )
                                .filter(
                                    (
                                        UniversalCompartmentalizedComponent.bigg_id
                                        == metabolite_id
                                    )
                                    & (
                                        ModelCompartmentalizedComponent.model_id
                                        == model_db_id
                                    )
                                )
                                .limit(1)
                            ).first()
                        if model_comp_comp_db is None:
                            print(f"Metabolite not found: {k}")
                            continue
                        specific_results[k][
                            "model_compartmentalized_component_id"
                        ] = model_comp_comp_db.id
                    else:
                        print(f"Type of test currently not implemented: {test_type}")
                        print(test_result)

                    if is_dict:
                        v = prop_val[k]
                        if v is None:
                            continue
                        if prop in ["data", "metric", "duration"]:
                            specific_results[k][prop] = float(v)
                        else:
                            specific_results[k][prop] = v

        general_result = MemoteResult(**general_result)
        session.add(general_result)

        for res in specific_results.values():
            res_db = MemoteResult(**res)
            session.add(res_db)
        session.commit()
        session.close()
