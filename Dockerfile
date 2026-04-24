FROM python:3.10
ENV PYTHONUNBUFFERED=1
RUN mkdir /app
WORKDIR /app

# Install dependencies

RUN apt-get update && apt-get install -y \
  libxml2-dev
COPY requirements.txt /app
RUN pip install -r requirements.txt

RUN git clone https://github.com/biosustain/bigg_models_data.git bigg_models_data

RUN git clone https://github.com/biosustain/biggr_maps.git biggr_maps
WORKDIR /app/biggr_maps
RUN python setup.py install

WORKDIR /app
COPY ./ /app
RUN python setup.py install
# Install reqs  

#COPY settings.ini /app/settings.ini
#COPY settings.small.ini /app/settings.ini

# CMD ["bin/ecmdb_helper"]
# CMD ["bin/metanetx_helper"]
# CMD ["bin/acp_helper"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-taxonomy", "--skip-curated-reactions"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-model-processing"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-taxonomy", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-genomes", "--skip-compartments"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-taxonomy", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-models", "--skip-model-processing", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps"]
# CMD ["bin/load_db", "--drop-models", "--skip-rhea", "--skip-taxonomy", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-model-processing", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps"]
# CMD ["bin/load_db", "--drop-all", "--skip-model-processing", "--skip-collections", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps"]
# CMD ["bin/load_db", "--drop-models", "--skip-rhea", "--skip-curated-reactions", "--skip-curated-metabolites", "--skip-taxonomy", "--skip-model-processing", "--skip-collections", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps"]
# CMD ["bin/load_db", "--skip-models", "--skip-rhea", "--skip-curated-reactions", "--skip-curated-metabolites", "--skip-collections", "--skip-taxonomy", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps"]
# CMD ["bin/load_db", "--drop-models", "--drop-genomes", "--skip-memote", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-taxonomy", "--skip-maps"]
# CMD ["bin/load_db", "--drop-genomes", "--skip-compartments", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-taxonomy", "--skip-models", "--skip-collections", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-maps", "--skip-model-processing", "--skip-memote"]
# CMD ["bin/load_db", "--drop-all", "--skip-genomes", "--skip-curated-reactions", "--skip-models"]
# CMD ["bin/load_db", "--skip-genomes"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-curated-metabolites"]
# CMD ["bin/load_db", "--skip-rhea", "--skip-curated-metabolites"]
# CMD ["bin/load_db", "--drop-models", "--skip-compartments", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-taxonomy", "--skip-genomes"]
# CMD ["bin/load_db", "--drop-maps", "--skip-compartments", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-taxonomy", "--skip-genomes", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-model-processing", "--skip-collections"]
# CMD ["bin/load_db", "--skip-compartments", "--skip-rhea", "--skip-curated-metabolites", "--skip-taxonomy", "--skip-genomes", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-model-processing", "--skip-collections", "--skip-maps"]
# CMD ["bin/load_db", "--drop-all", "--skip-compartments", "--skip-curated-metabolites", "--skip-taxonomy", "--skip-genomes", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-model-processing", "--skip-collections", "--skip-maps"]
# CMD ["bin/load_db", "--skip-compartments", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-taxonomy", "--skip-genomes", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-collections", "--skip-maps"]
# CMD ["bin/load_db", "--drop-models", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-genomes", "--skip-seed-metabolites", "--skip-seed-reactions"]
# CMD ["bin/load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-models", "--skip-seed-metabolites", "--skip-seed-reactions"]
# CMD ["bin/load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-models"]
# CMD ["bin/load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions"]
# CMD ["bin/load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-models", "--skip-memote", "--skip-compartments", "--skip-seed-metabolites"]
CMD ["bin/load_db", "--drop-all"]
# CMD ["bin/load_db", "--drop-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-compartments", "--skip-taxonomy"]
# CMD ["bin/load_db", "--drop-all", "--skip-genomes", "--skip-models", "--skip-model-processing", "--skip-collections", "--skip-maps", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-taxonomy"]
# CMD ["bin/load_db", "--skip-rhea"]

# CMD ["bin/load_db", "--skip-genomes", "--skip-rhea", "--skip-curated-metabolites", "--skip-curated-reactions", "--skip-models", "--skip-memote", "--skip-compartments", "--skip-seed-metabolites", "--skip-seed-reactions", "--skip-model-processing"]

# CMD ["bin/test_parts"]
