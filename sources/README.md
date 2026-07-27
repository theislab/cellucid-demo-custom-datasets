# Synthetic source design

Every value in `exports/` is created by `scripts/build_datasets.py`; no source
data is downloaded.

The generator uses fixed arithmetic sequences rather than a random-number
generator. It records a fixed UTC creation timestamp, uses exact dependency
versions, invokes Cellucid’s current `prepare()` writer, generates the current
catalog with `generate_datasets_manifest()`, and writes a sorted SHA-256
inventory. `--check` constructs a clean temporary tree and compares every byte.

The names beginning with `SYN_` are teaching labels, not measured genes.
Population, lineage, batch, condition, quality, and response labels are
invented. The examples must not be used for biological inference.
