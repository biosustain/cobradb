#!/usr/bin/env python
"""One-off: add the model.method column and set it to "CarveMe" for the models
in the iPKB_Clostridiaceae collection.

These models are already in the database, so this bypasses the ETL entirely.
The ETL will set method at load time in a later change; until then this script
is the only record of where these values came from.

There is no migration framework in this repo (schema comes from
Base.metadata.create_all, which creates missing tables but not missing
columns), so the ALTER TABLE lives here alongside the UPDATE to keep the two
from drifting apart. Both statements are idempotent — rerunning is a no-op.
"""
import logging

from sqlalchemy import text

from cobradb.models import Model, ModelCollection, Session

logging.basicConfig(level=logging.INFO)

COLLECTION = "iPKB_Clostridiaceae"
METHOD = "CarveMe"

with Session() as session:
    session.execute(text("ALTER TABLE model ADD COLUMN IF NOT EXISTS method varchar"))
    session.commit()
    print("model.method column present")

    collection = (
        session.query(ModelCollection)
        .filter(ModelCollection.bigg_id == COLLECTION)
        .first()
    )
    if collection is None:
        raise SystemExit(f"No collection found with BiGG ID {COLLECTION}")

    models = session.query(Model).filter(Model.collection_id == collection.id).all()
    for m in models:
        m.method = METHOD
    session.commit()
    print(f"Done. Set method={METHOD!r} on {len(models)} models in {COLLECTION}")
