# Embedded license notices in model files

Audit date: 2026-10-22

## Result

**No model served on bigg.bio carries its own license notice.**

One model file in the bigg_models_data repository carries an embedded notice, but it is not
loaded into BiGG:

| Model id | File | Notice text | License | Served on bigg.bio | Notice in SBML / JSON downloads |
|---|---|---|---|---|---|
| yeast_7.00 | `bigg_models_data/original_models/yeast_7.00.xml` | "This SBML representation of the yeast metabolic network is made available under the Creative Commons Attribution-Share Alike 3.0 Unported Licence (see www.creativecommons.org)." (model-level SBML `<notes>`) | CC BY-SA 3.0 | No | Not applicable (no download exists) |

## What was searched

Case-insensitive search for `license`, `licence`, `creativecommons`,
`copyright`, `CC BY` / `CC-BY`, `dc:rights`, `dcterms:rights` and `"rights"` in:

| Location | Files | Hits |
|---|---|---|
| Source models loaded by cobradb (`biggr-input-models`: 3892 JSON, 1006 SBML, 2 MAT) | 4900 | 0 |
| `bigg_models_data/original_models/` and `bigg_models_data/iAF987_conv.xml` | 213 | 1 (`yeast_7.00.xml`) |
| `cobradb/test_data/` | 10 | 0 |
| Exported downloads (`*.biggr.sbml`, `*.biggr.json`) for all 4881 served models | 9762 | 0 |

The two MAT files (`iHN637.mat`, `iLJ478.mat`) store strings in compressed
blocks, so they were additionally loaded with `scipy.io.loadmat` and every
string field was searched; no hits.

Model-level metadata in the served source files was also read directly, to
rule out notices phrased without any of the search terms:

- 4801 of the 4900 files have no model-level `<notes>` or `<annotation>` at
  all (JSON files only contain `id`, `version` and `compartments` at the top
  level).
- 97 CarveMe-built SBML files (`GCF_*.xml`) have a model-level annotation that
  only says "This model was built with CarveMe version 1.6.4".
- `ecoli_core_model.xml` has the note "COBRA_Toolbox".
- `iJN1462.json` has the note "Added Detect Reactions and corrected metabolite
  names / Added Levulinic Acid Catabolism".

Served models were taken from the `model` table of the production database
(4881 models). `yeast_7.00` is not among them (the served yeast models are
`iMM904` and `iND750`).

## Does the cobradb export keep model notes and annotations?

**No, not reliably.** `cobradb/model_processing.py` reads each source model
with cobrapy (`read_sbml_model` / `load_json_model` / `load_matlab_model`),
modifies it and writes it with `write_sbml_model` and `save_json_model`.

Round-trip test with `yeast_7.00.xml` and cobrapy 0.29.1 (the version in the
production containers; `requirements.txt` allows `>=0.14.2,<0.30.0`):

- `model.notes` after reading: `{}`. cobrapy only keeps notes of the form
  `<p>key: value</p>`; free-text notes such as this license sentence are dropped.
- `model.annotation` after reading: only CV terms (here `taxonomy`). The RDF
  `dc:creator` block is dropped.
- The written SBML and JSON files do not contain the notice.

So if a model with an embedded license notice is added to BiGG in the future,
its notice will **not** appear in the downloads unless the export is changed
(for example by copying the source model's notes into `model.notes` as a
`License: ...` entry before writing). Re-run this audit whenever new source
models are added.
