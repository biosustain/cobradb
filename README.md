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

It is good practice to ensure reproducibility and upload all ETL related input,
output, tools, config files to a backup storage every time you do a full update.
Archived runs are stored at
[biggr storage account](https://portal.azure.com/#@dtudk.onmicrosoft.com/resource/subscriptions/aee8556f-d2fd-4efd-a6bd-f341a90fa76e/resourceGroups/rg-recon/providers/Microsoft.Storage/storageAccounts/biggr/overview)
You can update the SAS token and the model list json file in the .env file. Then
run the upload script with

```bash
cd ../cobradb && bin/upload_archive
```

If you don't have the azcopy installed already, you can install the non-snap
version with the following command

```bash
wget -O /tmp/azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux
tar xf /tmp/azcopy.tar.gz -C /tmp
sudo mv /tmp/azcopy_linux_amd64_*/azcopy /usr/local/bin/azcopy
```

This standalone binary can run unsandboxed and read any path the calling user
has access to.

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
