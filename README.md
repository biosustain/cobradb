# Cobradb

Cobradb loads genome-scale metabolic models (GEMs) and genome annotations into a
relational PostgreSQL database. It works both as an independent ETL pipeline
that pipes data into [BiGGr Models](https://biggr.org) database, and as an
installed package providing ORM models that the
[BiGGr Models](https://biggr.org) web app uses to query the database.

## What cobradb does as the ETL pipeline

1. Reads SBML/JSON metabolic models, NCBI assemblies and model-genome-set.json
   as input, and refers to Rhea, ChEBI, curated annotation files, etc. in the
   mounted directories.
2. Normalizes IDs, deduplicates entities, applies curation rules.
3. Writes data including models, genomes, GPRs, cross references to a PostgreSQL
   database.
4. Produces standardized SBML files and runs MEMOTE quality checks on each
   model, writing both results to the output folder.

## Running the pipeline

### 1. Set up configuration and insert input files

- adjust variables in `settings.ini` when necessary
- create a new model-genome-set.json file with unique traceable name (e.g.
  `/app/model-genome-set-8.json`), follow practices and add the new models
- insert newly added models and corresponding NCBI assemblies (when applicable)
  into `/data/biggr-input-models` and `/data/biggr-assemblies/ncbi_dataset/data`
  directories

### 2. Create new output folders and run docker compose

The recommended way to run the pipeline is to use full update, which uses the
default CMD ["bin/load_db", "--drop-all"] defined in the Dockerfile. The docker
compose would set up an ephemeral pipeline service and a new local PostgreSQL
instance. Note that in the same VM, the BiGGr web app and the production
PostgreSQL instance are also running. The BiGGr web app reads model files from
`/data/biggr-bigg-models`, and the production PostgreSQL instance stores its
data in `/data/biggr-postgres`. Therefore, it is important to use separate
output directories to avoid writing into the production database. You can create
new folders with

```bash
mkdir -p /data/biggr-bigg-models-new /data/biggr-postgres-new
```

and then run docker compose with

```bash
docker compose up --build --force-recreate -d
```

### 3. Switch to production and archive ETL artifacts

After the pipeline is done, stop the ETL stack so the cobradb PostgreSQL
container releases its lock on `/data/biggr-postgres-new` (two postgres
instances cannot share the same data directory):

```bash
docker compose down
```

Then test the new database with the web app before discarding the old data.
Edit `../biggr_models/docker-compose.yml` to point the `biggr-web` and
`biggr-db` volumes at the `-new` paths:

```yaml
biggr-web:
  volumes:
    - /data/biggr-bigg-models-new/:/models
biggr-db:
  volumes:
    - /data/biggr-postgres-new:/var/lib/postgresql/data
```

Bring the web app down and back up against the new data:

```bash
cd ../biggr_models && docker compose --profile prod down
docker compose --profile prod up -d
```

Once the results are verified, promote the new directories and remove the
old ones. Keep the old data on disk until you're sure — this step is
irreversible.

```bash
sudo rm -rf /data/biggr-bigg-models /data/biggr-postgres
sudo mv /data/biggr-bigg-models-new /data/biggr-bigg-models
sudo mv /data/biggr-postgres-new   /data/biggr-postgres
```

Revert the volume paths in `biggr_models/docker-compose.yml` back to the
non-`-new` paths and restart:

```bash
docker compose --profile prod down
docker compose --profile prod up -d
```

### 4. Archive the run to Azure Blob Storage

It is good practice to ensure reproducibility and upload all ETL related input,
output, tools, config files to a backup storage every time you do a full update.
Archived runs are stored in the `backups` container of the
[biggr storage account](https://portal.azure.com/#@dtudk.onmicrosoft.com/resource/subscriptions/aee8556f-d2fd-4efd-a6bd-f341a90fa76e/resourceGroups/rg-recon/providers/Microsoft.Storage/storageAccounts/biggr/overview),
one folder per run (`backups/<YYYYMMDD>/`).

This archive is the data-level backup. For the VM-level backup (disk
snapshots), see the "Backup and restore" section in the
[biggr_models README](https://github.com/biosustain/bigg_models).

`bin/upload_archive` uploads:

| What | How | Source |
|---|---|---|
| Raw input models | `azcopy`, copied as-is | `/data/biggr-input-models` |
| Standardized models | `azcopy`, copied as-is | `/data/biggr-bigg-models` |
| Genome assemblies, data files, ChEBI cache | `azcopy`, copied as-is | `/data/biggr-assemblies`, `/data/biggr-datafiles`, `/data/biggr-chebi-cache` |
| Model list used for the run | `azcopy` | `$USED_JSON_FILE` |
| PostgreSQL database | `pg_dump -Fc` → `biggr_<YYYYMMDD>.dump` | `biggr` database in the `biggr-postgres` container |
| PostgreSQL roles (without passwords) | `pg_dumpall --globals-only --no-role-passwords` → `biggr_globals_<YYYYMMDD>.sql` | same |


To run it, update `SAS` (a SAS token for the `backups` container) and
`USED_JSON_FILE` in `.env`, then run the script from the `cobradb/` directory:

```bash
cd ../cobradb && bin/upload_archive
```

A full run takes a while (about 30 min for the database dump, plus about 70 GB
of uploads), so run it outside peak hours. `azcopy` and `pg_dump` run at the
lowest CPU/IO priority and uploads are capped at 400 Mbps, so that the website
stays responsive.

The archive contains no passwords: the roles file is dumped with
`--no-role-passwords`. The `biggr` storage account does not allow public
access, so the archive can only be read with a SAS token or Azure login. 

If you don't have the azcopy installed already, you can install the non-snap
version with the following command

```bash
wget -O /tmp/azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux
tar xf /tmp/azcopy.tar.gz -C /tmp
sudo mv /tmp/azcopy_linux_amd64_*/azcopy /usr/local/bin/azcopy
```

This standalone binary can run unsandboxed and read any path the calling user
has access to.

### Restoring the database from an archive

Download `biggr_<YYYYMMDD>.dump` and `biggr_globals_<YYYYMMDD>.sql` from the
archive folder, then on an empty PostgreSQL instance restore the roles first
and the data second:

```bash
psql -f biggr_globals_<YYYYMMDD>.sql              # roles
createdb biggr                                     # skip if it already exists
pg_restore -d biggr -j 4 biggr_<YYYYMMDD>.dump     # data, 4 parallel jobs
```

The roles file has no passwords. With the Docker setup, the `biggr-postgres`
container sets the password from `POSTGRES_PASSWORD` in `.env` on first start.

Run them against the PostgreSQL instance that will serve the website, e.g. the
`biggr-db` container started on an empty data directory.

## model-genome-set.json reference

The model list json file has the following shape:
`{"models": [entry1, entry2, ...]}`. Each entry supports:

| Field       | Required | Description                                                                                                                                                                                      |
| ----------- | -------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| filename    | yes      | Model file (.xml or .json), resolved under model_directory in settings.ini                                                                                                                       |
| assembly    | no       | source:accession. Supported sources: ncbi_assembly (default), pankb_assembly, patric_assembly. ncbi_assembly needs a genomic.gbff[.gz] under /data/biggr-assemblies/ncbi_dataset/data/accession/ |
| chromosomes | no       | Replicon accessions, e.g. ["NC_000913.3"]                                                                                                                                                        |
| organism    | no       | Overrides organism inferred from assembly                                                                                                                                                        |
| tax_id      | no       | Overrides tax id inferred from assembly or organism                                                                                                                                              |
| collection  | no       | BiGG collection id; defaults to a single-model collection                                                                                                                                        |
| prefix      | no       | Prepended to the model's BiGG id at load time                                                                                                                                                    |
| pub_ref     | no       | type:value, e.g. pmid:21988831 or doi:10.5281/zenodo.17581962                                                                                                                                    |
