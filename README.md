# Cellucid custom dataset examples

Welcome! This repository is a small, inspectable model for publishing your own
[Cellucid](https://www.cellucid.com/) prepared datasets.

It contains three fully synthetic examples, their complete browser-ready
artifacts, the one script that generated them, strict validators, and a
cross-platform CI workflow. You can clone it to learn the format, replace the
synthetic generator with your own export script, or copy its repository
structure into a new public data repository.

No example is real biology. No patient, donor, clinical, or private source data
is present, and the `SYN_` gene names are teaching labels.

## Use the public examples

This complete tree is published as
`theislab/cellucid-demo-custom-datasets` on the `main` branch. The GitHub loader
reads the public repository directly from `raw.githubusercontent.com`. Because
this catalog contains three datasets, the shareable URL always names the exact
dataset:

- [2D cell-type islands](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-cell-types-2d)
- [3D branching development](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-development-3d)
- [1D trajectory](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-trajectory-1d)

The exact repository value to paste into Cellucid’s **GitHub** connection is:

```text
theislab/cellucid-demo-custom-datasets/exports
```

Then choose one of the three IDs shown above. The compact form resolves only
the `main` branch; Cellucid does not probe other branches. The explicit
`theislab/cellucid-demo-custom-datasets@main/exports` form is equivalent.

The same files can be served from GitHub Pages. The exact roots are:

```text
https://theislab.github.io/cellucid-demo-custom-datasets/
https://theislab.github.io/cellucid-demo-custom-datasets/exports/
https://theislab.github.io/cellucid-demo-custom-datasets/exports/datasets.json
```

This direct Cellucid URL points its catalog transport at the Pages exports
root:

```text
https://www.cellucid.com/?exportsBaseUrl=https%3A%2F%2Ftheislab.github.io%2Fcellucid-demo-custom-datasets%2Fexports%2F&dataset=synthetic-cell-types-2d
```

## What each tiny dataset teaches

| Dataset ID | Cells | Genes | Obs metadata | View | Graph | Vector field | Exact dataset bytes |
|---|---:|---:|---|---|---:|---|---:|
| `synthetic-cell-types-2d` | 72 | 8 | cell type, batch, quality, library size | 2D Planar | 144 weighted edges | — | 5,006 |
| `synthetic-development-3d` | 96 | 10 | lineage, timepoint, condition, differentiation, response | 3D Orbit | 187 weighted edges | 3D `velocity_umap` | 7,675 |
| `synthetic-trajectory-1d` | 48 | 6 | stage, replicate, pseudotime, activity | 1D Planar | 93 weighted edges | 1D `velocity_umap` | 4,898 |

The complete `exports/` tree is 27,085 bytes in 72 files, including the catalog
and checksum inventory. Payloads use gzip level 6, 8-bit gene and continuous
metadata quantization, and `uint8` categorical codes. A browser requests only
the selected dataset and fields, so adding more catalog entries does not make
every visitor download every gene.

## Repository map

```text
.
├── .github/workflows/test.yml
├── .gitattributes
├── .gitignore
├── .nojekyll
├── exports/
│   ├── datasets.json
│   ├── SHA256SUMS
│   ├── synthetic-cell-types-2d/
│   ├── synthetic-development-3d/
│   └── synthetic-trajectory-1d/
├── scripts/
│   ├── __init__.py
│   ├── build_datasets.py
│   └── validate_exports.py
├── sources/
│   └── README.md
├── tests/
│   ├── catalog-contract.test.mjs
│   └── test_contract.py
├── CODE_OF_CONDUCT.md
├── CONTRIBUTING.md
├── LICENSE
├── README.md
├── SECURITY.md
├── SUPPORT.md
├── favicon.svg
├── index.html
└── requirements-build.txt
```

Within each dataset:

- `dataset_identity.json` is the canonical identity and feature inventory.
- `points_1d.bin.gz`, `points_2d.bin.gz`, or `points_3d.bin.gz` stores aligned
  Float32 coordinates.
- `obs_manifest.json` and `obs/` describe and store cell metadata.
- `var_manifest.json` and `var/` describe and store lazy gene values.
- `connectivity_manifest.json` and `connectivity/` store one canonical weighted
  undirected graph.
- `vectors/` is present only when `dataset_identity.json` declares an aligned
  vector field.

Every array uses the same cell order. Changing that order changes the dataset’s
meaning even if its ID and filenames stay the same.

## Rebuild these examples exactly

Cellucid 0.9.1 requires Python 3.11–3.14. Use a clean environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install --requirement requirements-build.txt
```

On Windows PowerShell, activate with:

```powershell
.venv\Scripts\Activate.ps1
```

Validate the checked-in bytes without changing them:

```bash
python scripts/validate_exports.py
python -m unittest discover -s tests -p "test_*.py" -v
node --test tests/catalog-contract.test.mjs
python -m ruff check scripts tests
python -m ruff format --check scripts tests
python scripts/build_datasets.py --check
```

`--check` creates a clean temporary generation, compares every path and byte
with `exports/`, and removes the temporary tree. To intentionally publish a
complete new local generation after editing the source:

```bash
python scripts/build_datasets.py --force
python scripts/build_datasets.py --check
```

The requirements file pins Cellucid 0.9.1 to the exact audited
`cellucid-python` source commit that generated these bytes, along with exact
NumPy, pandas, SciPy, and Ruff versions. The builder refuses an inexact data
dependency version before touching `exports/`. It stages and validates the
whole catalog before an atomic directory replacement. `SHA256SUMS` owns every
catalog, manifest, embedding, metadata, gene, graph, and vector file.

## Create your own dataset from Python

Start from aligned arrays. The number and order of cells must agree in:

- every embedding;
- every `obs` column;
- every row of `gene_expression`;
- both axes of `connectivities`;
- every vector field.

The number and order of genes must agree between `var` and the columns of
`gene_expression`. Use finite, real numeric values. Make categorical columns
explicit pandas categoricals so their category order is intentional.

This is a complete current 2D example:

```python
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import sparse

from cellucid import prepare

n_cells = 40
n_genes = 3

embedding = np.column_stack(
    [
        np.linspace(-1, 1, n_cells, dtype=np.float32),
        np.repeat([-0.5, 0.5], n_cells // 2).astype(np.float32),
    ]
)
obs = pd.DataFrame(
    {
        "cell_type": pd.Categorical(
            ["type-a"] * 20 + ["type-b"] * 20,
            categories=["type-a", "type-b"],
            ordered=True,
        ),
        "quality_score": np.linspace(0.7, 0.95, n_cells, dtype=np.float32),
    }
)
var = pd.DataFrame(index=pd.Index(["GENE_A", "GENE_B", "GENE_C"]))
expression = np.column_stack(
    [
        np.linspace(0, 5, n_cells, dtype=np.float32),
        np.linspace(5, 0, n_cells, dtype=np.float32),
        np.tile(np.asarray([0, 1], dtype=np.float32), n_cells // 2),
    ]
)

graph = np.zeros((n_cells, n_cells), dtype=np.float64)
for index in range(n_cells - 1):
    graph[index, index + 1] = 1.0
    graph[index + 1, index] = 1.0

velocity = np.column_stack(
    [
        np.full(n_cells, 0.04, dtype=np.float32),
        np.zeros(n_cells, dtype=np.float32),
    ]
)

prepare(
    latent_space=embedding,
    obs=obs,
    var=var,
    gene_expression=expression,
    connectivities=sparse.csr_matrix(graph),
    out_dir="exports/my-study-v1",
    obs_keys=["cell_type", "quality_score"],
    centroid_outlier_quantile=0.95,
    centroid_min_points=6,
    force=False,
    var_quantization=8,
    obs_continuous_quantization=8,
    compression=6,
    obs_categorical_dtype="uint8",
    dataset_name="My study",
    dataset_id="my-study-v1",
    created_at="2026-07-27T00:00:00Z",
    dataset_description="A public, documented example dataset.",
    X_umap_2d=embedding,
    vector_fields={"velocity_umap_2d": velocity},
)
```

Important current contracts:

- Provide at least one exact `X_umap_1d`, `X_umap_2d`, or `X_umap_3d` array.
  Each must have exactly 1, 2, or 3 columns.
- Vector keys use `<field>_umap_<dimension>d`, such as
  `velocity_umap_2d`. The suffix must match the vector’s column count and an
  available embedding dimension.
- `connectivities` is a square, exactly symmetric, finite, non-negative graph
  with an exact zero diagonal. Sparse input must not store zeros or duplicate
  coordinates.
- `dataset_id`, observation keys, and gene identifiers are portable file
  components. Keep them short, stable, case-distinct, and limited to letters,
  digits, `.`, `_`, and `-`.
- `uint8` supports at most 255 categories per field; choose `uint16` explicitly
  when a field needs more.
- 8-bit quantization is compact and well suited to interactive color mapping,
  but it is lossy. Use 16-bit or unquantized Float32 when the visual precision
  requirement justifies the extra bytes.
- The prepared export is a visualization and exploration artifact. Preserve
  the original analysis object and provenance separately.

After exporting one or more dataset directories, create the required catalog:

```python
from cellucid.prepare_data import generate_datasets_manifest

generate_datasets_manifest(
    "./exports",
    default_dataset="my-study-v1",
)
```

The default must be an exact ID present in the catalog. A catalog with several
datasets still requires an exact `dataset=` value in a shareable GitHub or
remote-server URL; Cellucid does not choose an arbitrary entry.

## Create an export from R

The R package writes the same current prepared contract. A minimal aligned
example is:

```r
library(cellucid)

n_cells <- 40L
embedding <- cbind(
  seq(-1, 1, length.out = n_cells),
  rep(c(-0.5, 0.5), each = n_cells / 2L)
)
obs <- data.frame(
  cell_type = factor(
    c(rep("type-a", 20L), rep("type-b", 20L)),
    levels = c("type-a", "type-b")
  ),
  quality_score = seq(0.7, 0.95, length.out = n_cells)
)
expression <- cbind(
  GENE_A = seq(0, 5, length.out = n_cells),
  GENE_B = seq(5, 0, length.out = n_cells),
  GENE_C = rep(c(0, 1), n_cells / 2L)
)
var <- data.frame(row.names = colnames(expression))
graph <- matrix(0, nrow = n_cells, ncol = n_cells)
for (index in seq_len(n_cells - 1L)) {
  graph[index, index + 1L] <- 1
  graph[index + 1L, index] <- 1
}

cellucid_prepare(
  latent_space = embedding,
  obs = obs,
  var = var,
  gene_expression = expression,
  connectivities = graph,
  out_dir = "exports/my-r-study-v1",
  obs_keys = c("cell_type", "quality_score"),
  centroid_min_points = 6L,
  var_quantization = 8L,
  obs_continuous_quantization = 8L,
  obs_categorical_dtype = "uint8",
  compression = 6L,
  dataset_name = "My R study",
  dataset_description = "A public, documented R example dataset.",
  dataset_id = "my-r-study-v1",
  X_umap_2d = embedding
)
```

Create `exports/datasets.json` after all R datasets are present. You can use
the Python `generate_datasets_manifest()` call above or create the exact
version-1 object yourself and run this repository’s validator. Do not edit a
generated dataset manifest or binary by hand.

## Validate your own repository

Adjust the `EXPECTED` inventory in `scripts/validate_exports.py` and
`tests/catalog-contract.test.mjs` to describe your deliberate datasets. The
validators then check:

- exact catalog, identity, obs, var, connectivity, and vector schemas;
- catalog-to-identity names, descriptions, cell counts, and gene counts;
- canonical root-relative paths with no traversal;
- exact embedding, observation, gene, graph, and vector byte lengths;
- finite Float32/Float64 payloads;
- categorical code bounds;
- strictly ordered unique connectivity pairs with `source < destination`;
- truthful maximum graph degree;
- deterministic gzip timestamps;
- a complete sorted SHA-256 inventory;
- no undeclared artifact in a dataset or at the exports root.

Run both Python and Node validations because they challenge the same artifacts
through independent implementations. The included workflow executes static
contracts on Ubuntu, macOS, and Windows with Python 3.11/3.14 and Node 20/26,
then performs one byte-identical clean rebuild.

For a real project, also load every dataset through the visible Cellucid
workflow before publication. Check field replacement, filtering, highlighting,
connectivity, vectors, 1D/2D/3D navigation, browser console output, and network
requests.

## Test locally in Cellucid

For the prepared-data server:

```bash
cellucid serve ./exports --port 8765
```

Open the exact **Viewer URL** printed by the command. A prepared server URL
ends in:

```text
/?source=remote
```

The server’s dataset endpoint is:

```text
http://127.0.0.1:8765/_cellucid/datasets
```

If you intentionally change the port, use the exact Viewer URL and endpoint
printed by the server. This repository has three datasets, so choose the exact
ID in the UI rather than expecting one to be selected by position.

You can also open <https://www.cellucid.com/>, choose the prepared-directory
control, and select one individual directory below `exports/`. The folder must
start at the directory containing `dataset_identity.json`, not the parent
catalog directory.

## Publish with GitHub Pages

GitHub’s current branch-publishing flow supports a repository root as the Pages
source:

1. Create a public repository named `cellucid-demo-custom-datasets`.
2. Put this complete tree on the `main` branch.
3. Open **Settings → Pages**.
4. Under **Build and deployment**, select **Deploy from a branch**.
5. Select `main`, select `/(root)`, and save.
6. Wait for the Pages deployment workflow to finish.

`.nojekyll` ensures the prepared paths are published unchanged. GitHub warns
that Pages sites are public even for many private-repository plans, so perform
the privacy review in `SECURITY.md` before enabling it.

Verify the exact response contract after deployment:

```bash
curl -sS -D - -o /dev/null \
  https://theislab.github.io/cellucid-demo-custom-datasets/exports/datasets.json

curl -sS -D - -o /dev/null -H "Range: bytes=0-31" \
  https://theislab.github.io/cellucid-demo-custom-datasets/exports/synthetic-cell-types-2d/points_2d.bin.gz
```

The catalog request must return `200`. The range request must return `206` with
`Content-Range`. Both responses must include:

```text
Access-Control-Allow-Origin: *
Accept-Ranges: bytes
```

These headers let an HTTPS Cellucid page fetch JSON and stream binary bytes
from the Pages origin. A custom host must provide the same cross-origin and
range behavior for `GET` and `HEAD`; a browser cannot repair a server that
rejects the request.

GitHub’s official setup and privacy details are in
[Configuring a publishing source](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

## Host somewhere other than GitHub

Any static HTTPS host is suitable when it:

- preserves every case-sensitive path under `exports/`;
- serves `datasets.json` and JSON manifests as JSON;
- serves gzip payload bytes without applying a second content encoding;
- permits public `GET` and `HEAD`;
- supports byte ranges;
- sends `Access-Control-Allow-Origin: *`, or names the exact Cellucid origin;
- does not require a cookie, bearer token, redirect login, or expiring query
  parameter for public data.

Point Cellucid at the canonical absolute directory URL ending in `/`:

```text
https://www.cellucid.com/?exportsBaseUrl=https%3A%2F%2Fdata.example.org%2Fexports%2F&dataset=my-study-v1
```

If the data must remain private, do not put it on a public static host. Use the
Cellucid server over a controlled network or SSH tunnel and follow your
institution’s data-governance rules.

## Size and performance guidance

This repository is intentionally only 27,085 bytes. Real exports can be much
larger because each selected gene is a separate lazy payload.

- Export only the gene identifiers users need to color and query.
- Keep metadata fields purposeful; category labels are downloaded with the
  selected dataset.
- Use gzip and deliberate quantization.
- Keep embedding and vector dimensions aligned.
- Avoid committing intermediate AnnData, notebooks with embedded data, caches,
  or duplicate export generations.
- Validate actual transfer sizes and request counts in browser developer tools.

GitHub warns above 50 MiB and blocks regular Git files above 100 MiB.
[Git LFS is not supported by GitHub Pages](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-git-large-file-storage),
so use Cellucid server mode or a static/object host designed for large files
when any payload approaches GitHub’s limit. The
[GitHub large-file guide](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github)
also recommends keeping repositories small.

## Privacy checklist before publishing

- Confirm every source permits public redistribution.
- Remove direct and indirect identifiers from obs values and category labels.
- Remove secrets and private URLs from scripts, notebooks, logs, manifests, and
  Git history.
- Confirm that rare combinations of metadata cannot identify a person.
- Remember that compression, quantization, and synthetic display names do not
  anonymize real measurements.
- Open the exact public catalog and every dataset from a signed-out browser.
- Review the host’s access logging, retention, and deletion policies.

See `SECURITY.md` for the repository policy.

## Troubleshooting

### “datasets.json not found”

The GitHub path must point to the directory containing `datasets.json`:

```text
theislab/cellucid-demo-custom-datasets/exports
```

Do not point at the repository root or one dataset folder. Confirm the raw URL:

```text
https://raw.githubusercontent.com/theislab/cellucid-demo-custom-datasets/main/exports/datasets.json
```

### “Multiple datasets require a selection”

Add the exact ID:

```text
&dataset=synthetic-development-3d
```

The `default` catalog field is for the configured sample catalog. A GitHub or
remote catalog with several entries requires an explicit choice.

### CORS or network error

Inspect the failed request in browser developer tools. Open the same URL with
`curl -D -` and confirm status, `Access-Control-Allow-Origin`, and range
headers. Fix the host configuration; do not add a second browser transport.

### Field fails to load

Run `python scripts/validate_exports.py`. A path, manifest tuple, dtype,
quantization setting, or payload length may disagree. Re-export the complete
dataset from its source script rather than editing generated bytes.

### Vector control is absent

Open `dataset_identity.json` and confirm `vector_fields` exists. The selected
view dimension must be listed under that field’s `available_dimensions`.
`velocity_umap_2d` does not describe a 3D vector.

### Category or legend error

Confirm every categorical code is either a valid category index or the
declared missing code. Keep category order stable and replace the entire
prepared generation atomically.

### Repository is too large

Use fewer exported genes, intentional quantization, and gzip. For large
scientific data, use `cellucid serve` or an HTTPS object host with CORS and
range support instead of GitHub Pages.

## Cellucid ecosystem

- [Cellucid web viewer](https://www.cellucid.com/)
- [Web viewer source](https://github.com/theislab/cellucid)
- [Python package and complete documentation](https://cellucid.readthedocs.io/)
- [Python package source](https://github.com/theislab/cellucid-python)
- [R package source](https://github.com/theislab/cellucid-r)
- [Canonical public sample catalog](https://github.com/theislab/cellucid-datasets)
- [Community annotation template](https://github.com/theislab/cellucid-annotation)
- [Data loading overview](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/01_loading_options_overview.html)
- [Dataset identity and reproducibility](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/06_dataset_identity_why_it_matters.html)
- [Python server reference](https://cellucid.readthedocs.io/en/latest/user_guide/python_package/g_api_reference_coverage/api/server.html)

## License

Code and generated synthetic examples are available under the BSD 3-Clause
License in `LICENSE`. The datasets make no biological or clinical claim.
