# Support

Thanks for using the Cellucid custom-dataset example.

## Getting help

- Documentation: https://cellucid.readthedocs.io
- The complete publishing workflow, including the repository contract and the
  fast-diagnosis table:
  https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/11_custom_dataset_repository.html
- One of the three synthetic examples here that will not open, or a step in
  `README.md` that does not work as written, is an issue here:
  https://github.com/theislab/cellucid-demo-custom-datasets/issues
- Anything about how the viewer connects, draws, filters, or exports what it
  loaded belongs in the web app: https://github.com/theislab/cellucid/issues
- Anything about `prepare()` — the export it wrote, its arguments, or an error
  it raised — belongs in the Python package:
  https://github.com/theislab/cellucid-python/issues, or in
  https://github.com/theislab/cellucid-r/issues for `cellucid_prepare()`.
- For security issues, see `SECURITY.md`.

## Your own repository

Most questions here are about a reader's own dataset repository rather than
this one. Two checks settle nearly all of them, and both are worth running
before opening an issue anywhere:

- Load the prepared folder locally, with **Prepared** under **Local data:**. If
  it is already wrong there, the problem is in the export, not in GitHub.
- Open your raw catalog URL in a browser tab:
  `https://raw.githubusercontent.com/<owner>/<repo>/<branch>/<path>/datasets.json`.
  It must return the JSON itself — not an HTML sign-in page, a Git LFS pointer,
  or a 404. The most common mistake is giving Cellucid the repository instead of
  the exports root inside it.

Quote the exact value you typed into the **GitHub data:** field and the exact
message Cellucid showed; the notification names both the plain cause and the
URL it actually requested.

## Data & privacy

The datasets in this repository are synthetic teaching material and contain no
patient, donor, or clinical information. Please don’t attach private data
(including patient-derived data) to issues or logs here. Prefer minimal
synthetic examples, anonymized screenshots, or small reproducible exports.

Publishing a prepared dataset publishes the coordinates, metadata, gene values,
graphs, vector fields, and provenance you included in it. `README.md` has the
privacy checklist to work through before a repository goes public.

## Code of Conduct

All project spaces are covered by `CODE_OF_CONDUCT.md`.
