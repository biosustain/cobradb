from typing import Optional, Tuple, Union

from urllib.parse import quote
import requests

import time


def rate_limit(calls_per_second: float):
    seconds_per_call = 1 / calls_per_second

    def decorator(func):
        last_call = time.time()

        def wrapper(*args, **kwargs):
            nonlocal last_call
            elapsed = time.time() - last_call
            if (wait_time := (seconds_per_call - elapsed)) > 0:
                time.sleep(wait_time)

            last_call = time.time()
            return func(*args, **kwargs)

        return wrapper

    return decorator


@rate_limit(calls_per_second=5)
def get_organism_for_ncbi_assembly_accession(
    accession: str,
) -> Optional[Tuple[str, int, Optional[str]]]:
    r = requests.get(
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{quote(accession)}/dataset_report",
        headers={"accept": "application/json"},
    )
    if r.status_code != 200:
        return None

    response = r.json()
    if response.get("total_count") != 1:
        return None

    report = response.get("reports", [None])[0]
    if report is None:
        return None

    organism = report.get("organism")
    if organism is None:
        return None

    organism_name = organism.get("organism_name")
    tax_id = organism.get("tax_id")

    if organism_name is None or tax_id is None:
        return None

    strain = organism.get("strain")
    if (
        strain is None
        and (infraspecific := organism.get("infraspecific_names")) is not None
    ):
        strain = infraspecific.get("strain")

    return organism_name, tax_id, strain


def resolve_organism(accession: str) -> Optional[Tuple[str, int, Optional[str]]]:
    """Look up an assembly accession, retrying without the version suffix.

    RefSeq versions get superseded (GCF_x.1 -> GCF_x.2) and the old versioned
    accession then returns nothing, while the unversioned one resolves to the
    current version. Accessions that are not GenBank/RefSeq assemblies (PATRIC
    genome ids, for instance) are skipped rather than sent to the API.
    """
    if accession is None or not accession.startswith(("GCA_", "GCF_")):
        return None

    info = get_organism_for_ncbi_assembly_accession(accession)
    if info is None and "." in accession:
        info = get_organism_for_ncbi_assembly_accession(accession.rsplit(".", 1)[0])
    return info
