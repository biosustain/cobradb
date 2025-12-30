# COBRAdb

COBRAdb loads genome-scale metabolic models and genome annotations into a
relational database. It powers [BiGGr Models](https://biggr.org), and
it is available under the MIT license. The bin/load_db script is the entrypoint
for building a database. The cobradb python library is also used by the
[biggr_models](https://github.org/biosustain/biggr_models) python package.

## Installation

It is recommended to run the program using Docker compose.

1. Install and set up Docker.

2. Set up your configuration. Copy the file `settings.ini.docker` and rename it
   `settings.ini`. Change the `postgres_password` property to a password of choice.
   Change other options (like the number of processes to use) as desired.

3. Similarly, copy the `env.docker` file to `.env` and set the same postgres
   password here. (Make sure to neither put the `settings.ini` or `.env` file under
   source control or share them in any other way.)

4. Create the following directories:
   - `/data/biggr-assemblies`
   - `/data/biggr-input-models`
   - `/data/biggr-chebi-cache`
   - `/data/biggr-datafiles`

5. Populate these directories with the required datafiles:
   `/data/biggr-chebi-cache` can be left empty. For `/data/biggr-assemblies` and
   `/data/biggr-input-models`, see the [biggr_models_data](https://github.org/biosustain/biggr_models_data) documentation. `/data/biggr-datafiles` files are
   available internally on the `biggr` storage account.

6. Build the docker container using `docker compose build` and run using `docker compose up`.
