# Contributing

Thank you for helping make custom Cellucid datasets easier to understand.

This repository contains only deterministic synthetic examples. Do not add
patient-derived, private, embargoed, identifying, credential-bearing, or
licensed third-party data. A proposed dataset must be generated entirely by a
reviewable script, stay small, demonstrate a distinct current Cellucid
capability, and make no biological claim.

Before opening a pull request:

1. Install the exact environment in `requirements-build.txt`.
2. Change `scripts/build_datasets.py` and the expected contract together.
3. Run `python scripts/build_datasets.py --force`.
4. Run `python scripts/build_datasets.py --check`.
5. Run `python scripts/validate_exports.py`.
6. Run `python -m unittest discover -s tests -p "test_*.py" -v`.
7. Run `node --test tests/catalog-contract.test.mjs`.
8. Run `python -m ruff check scripts tests`.
9. Run `python -m ruff format --check scripts tests`.
10. Inspect the complete diff and explain the capability and byte-size effect.

Keep dataset IDs stable after publication. A changed scientific row order or
meaning requires a new explicit dataset ID, because sessions, links, and
annotations depend on identity.

Report security concerns through `SECURITY.md`. General Cellucid product bugs
belong in the repository that owns the affected web, Python, R, dataset, or
annotation behavior.
