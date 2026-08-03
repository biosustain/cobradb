# -*- coding: utf-8 -*-

"""Load whole-database entity counts into the database_summary_count table.

These counts back the BiGGr front page cards. Each count MUST match the number
of rows the user sees when clicking through to the corresponding list page, so
the filters here mirror the DataHandler.pre_filter of each list view in
biggr_models.
"""

import logging
from datetime import datetime

from sqlalchemy import func, select

from cobradb.models import (
    Compartment,
    DatabaseSummaryCount,
    Genome,
    Model,
    ModelCollection,
    UniversalComponent,
    UniversalReaction,
)

# Canonical entity-type keys. The BiGGr front page reads these exact strings,
# so do not rename one without updating biggr_models/queries/summary_queries.py.
COLLECTIONS = "collections"
MODELS = "models"
METABOLITES = "metabolites"
REACTIONS = "reactions"
GENOMES = "genomes"
COMPARTMENTS = "compartments"

ENTITY_TYPES = (COLLECTIONS, MODELS, METABOLITES, REACTIONS, GENOMES, COMPARTMENTS)


def _count_statements():
    """Return {entity_type: select statement} for every summarized entity.

    Metabolites and reactions are restricted to the universal namespace
    (collection_id IS NULL) to match /universal/metabolites/ and
    /universal/reactions/, which filter the same way.
    """
    return {
        COLLECTIONS: select(func.count(ModelCollection.id)),
        MODELS: select(func.count(Model.id)),
        METABOLITES: select(func.count(UniversalComponent.id)).filter(
            UniversalComponent.collection_id.is_(None)
        ),
        REACTIONS: select(func.count(UniversalReaction.id)).filter(
            UniversalReaction.collection_id.is_(None)
        ),
        GENOMES: select(func.count(Genome.id)),
        COMPARTMENTS: select(func.count(Compartment.id)),
    }


def load_summary_counts(session):
    """Recompute every database summary count and upsert it.

    Idempotent: running this twice leaves exactly one row per entity type.
    Safe to run on a partially loaded database; entity types whose count query
    fails (e.g. a table that was dropped) are logged and skipped rather than
    aborting the whole step.
    """
    now = datetime.now()
    results = {}

    for entity_type, statement in _count_statements().items():
        try:
            count = session.scalar(statement)
        except Exception as e:
            # A failed statement aborts the postgres transaction, so roll back
            # before counting the next entity.
            session.rollback()
            logging.warning("Could not count %s: %s", entity_type, e)
            continue

        if count is None:
            count = 0

        row = (
            session.query(DatabaseSummaryCount)
            .filter(DatabaseSummaryCount.entity_type == entity_type)
            .first()
        )
        if row is None:
            row = DatabaseSummaryCount(
                entity_type=entity_type, count=count, date_time=now
            )
            session.add(row)
        else:
            row.count = count
            row.date_time = now

        results[entity_type] = count

    session.commit()
    logging.info("Loaded database summary counts: %s", results)
    return results