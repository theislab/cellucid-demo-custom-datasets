# Contributing to the Cellucid custom-dataset example

Contributions are welcome — corrections to the publishing guide, fixes to the
generator, and issues about an example dataset that will not open.

This file focuses on `cellucid-demo-custom-datasets`, the worked example of
publishing your own prepared datasets on GitHub. It describes contributing **to
this example**. If you are publishing your own data, you do not contribute
here: you build your own repository, and `README.md` is your guide.

By participating, you agree to follow the project’s Code of Conduct:
- `CODE_OF_CONDUCT.md`

If you’re reporting a security issue, please follow:
- `SECURITY.md`

---

## Which repo should I contribute to?

Cellucid is split by responsibility:

| Repo | What it is | Contribute here when you… |
|---|---|---|
| `cellucid` | Web app (UI + state + WebGL rendering) | are fixing how the viewer connects to GitHub, or how it draws what it loaded |
| `cellucid-python` | Python package + CLI + Sphinx docs | are fixing `prepare`/`serve`, or any documentation page for any repository |
| `cellucid-r` | R package exporter | are changing `cellucid_prepare()` |
| `cellucid-demo-custom-datasets` (this repo) | Worked example of publishing your own datasets | are changing this guide or its synthetic examples |
| `cellucid-datasets` | The published demo catalog | are correcting a published generation or its catalog entry |
| `cellucid-annotation` | Reference layout for annotation repositories | are changing the annotation schema, validation, or workflows |

A dataset that renders wrongly is usually a viewer issue, not an example issue.
Load the same folder with the **Prepared** local input first: if it is wrong
there too, it is the export; if it is only wrong over GitHub, it is the
connection.

The publishing workflow itself is documented once, in
[Publish a Custom Dataset Repository](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/11_custom_dataset_repository.html).
`README.md` here is the short path through it. A rule that changes in one has to
change in the other, in the same pass.

---

## What is authored here, and what is generated

Nothing under `exports/` is written by hand. The whole tree is produced by
`generate_datasets.py`, which holds the synthetic source data, the export
parameters, and the catalog call in one file, and replaces `exports/` whole.

So a fix to an example dataset is a fix to `generate_datasets.py` plus the
regenerated tree — never an edit to a payload file or a manifest. A hand-edited
export is invisible to review, and vanishes at the next generation.

---

## Testing & validation

There is deliberately **no CI in this repository**, so this check is yours to
run before you open a pull request. Nothing else will run it for you:

```bash
python generate_datasets.py --check
```

It rebuilds the complete catalog into a temporary directory and compares every
generated path and byte against the checked-in `exports/`, without modifying the
repository. It must print `PASS` and exit 0.

If you changed the synthetic inputs or the export parameters on purpose,
regenerate with `--force` and then re-run `--check`. README.md's *Regenerate
these example datasets* section is the single description of what those two
commands guarantee, which pins they enforce, and how to recover from an
interrupted swap — follow it rather than a second copy here.

Two consequences are yours as a contributor rather than the reader's:

- Regenerating with a different `cellucid` release changes the bytes and is a
  version bump, not a fix. Update the inline pins, `CELLUCID_RELEASE`, and the
  install command in `README.md` in the same change.
- The pins are what make the bytes reproducible, so a PR that regenerates
  `exports/` without saying which environment produced it cannot be reviewed.
  Say so, and paste the `--check` output.

---

## PR guidelines

- Keep PRs small and focused (one dataset, or one section of the guide).
- Include what changed, why, and the `--check` output.
- Regenerate rather than patch: attach the `generate_datasets.py` diff alongside
  the changed `exports/` tree.
- Keep the root minimal. `README.md` enumerates this repository's committed
  top-level contents exactly, and the documentation page above describes that
  root to readers who copy it, so a new root file means updating both in the
  same change. There are intentionally no tests, no validation harness, no
  package scaffolding, no workflows, no `.github/` directory, no requirements
  file, and no site assets.
- There is intentionally no `CITATION.cff` either. The datasets here are
  synthetic teaching material, not a result anyone should cite; the citable
  companions are the Python (`theislab/cellucid-python`) and R
  (`theislab/cellucid-r`) packages, whose `CITATION.cff` files carry `version:`.
  Adding one here would put a "Cite this repository" button on a tutorial, and
  would name this repository's author in every copy of the example.
