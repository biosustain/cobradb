import logging
from pathlib import Path
from typing import List, Optional

from cobradb.component_loading import load_assembly
from cobradb.memote import load_memote_results, run_memote
from cobradb.model_processing import process_model
from cobradb.models import AlreadyLoadedError, Session


def process_model_workflow(
    model_dir: Path, model_filename: str, model_data, skip_memote: bool = False
):
    with Session() as session:
        model_bigg_id = None
        try:
            model_bigg_id, _model_db_id = process_model(
                session,
                model_data,
                model_dir / model_filename,
            )
        except Exception as e:
            print(
                f"model_dir: '{model_dir}', model_path: '{model_dir / model_filename}'"
            )
            logging.error("Could not load model %s." % model_filename)
            logging.exception(e)
        logging.warning(f"Model bigg_id: {model_bigg_id}")

        if not skip_memote and model_bigg_id is not None:
            model_biggr_filename = f"/models/models/{model_bigg_id}.biggr.sbml"
            memote_result_filename = f"/models/memote/{model_bigg_id}.biggr.json.gz"
            Path(memote_result_filename).parent.mkdir(parents=True, exist_ok=True)

            successful_memote_run = False
            try:
                # if not Path(memote_result_filename).exists():
                run_memote(model_biggr_filename, memote_result_filename)
                successful_memote_run = True
            except Exception as e:
                logging.error("Memote failed to run.")
                logging.exception(e)

            if successful_memote_run:
                try:
                    load_memote_results(
                        model_bigg_id,
                        memote_result_filename,
                        session,
                    )
                except Exception as e:
                    logging.error("Could not load memote results %s." % model_filename)
                    logging.exception(e)
        else:
            logging.warning("Skipping MEMOTE")


def load_genomes_workflow(
    assembly_id: str, assembly_path: Path, chromosome_accessions: Optional[List[str]]
):
    logging.info(f"Loading assembly {assembly_id}")
    with Session() as session:
        try:
            load_assembly(
                assembly_id=assembly_id,
                assembly_path=assembly_path,
                chromosome_accessions=chromosome_accessions,
                session=session,
            )
            session.commit()
        except AlreadyLoadedError as e:
            logging.info(str(e))
        except Exception as e:
            logging.exception(e)
