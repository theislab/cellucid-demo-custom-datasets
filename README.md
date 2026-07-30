# Publish your own Cellucid datasets on GitHub

This is a deliberately small example repository. Its public data surface
contains only:

- this guide; and
- an `exports/` folder with three tiny datasets that Cellucid can open.

This repository also has a small, dependency-free validator, adversarial
tests, and cross-platform maintenance CI that protect the checked-in examples.
Those maintenance files are not needed by Cellucid and do not need to be
copied. Your own repository can still be just as simple as the two items above.

## Try this repository

Open [Cellucid](https://www.cellucid.com), choose the **GitHub** data source,
and enter:

```text
theislab/cellucid-demo-custom-datasets/exports
```

Then choose one of these small synthetic examples:

- `synthetic-cell-types-2d`
- `synthetic-development-3d`
- `synthetic-trajectory-1d`

You can also open them directly:

- [2-D cell types](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-cell-types-2d)
- [3-D development](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-development-3d)
- [1-D trajectory](https://www.cellucid.com/?github=theislab/cellucid-demo-custom-datasets/exports&dataset=synthetic-trajectory-1d)

These are teaching datasets, not biological results. They contain no patient,
donor, or clinical information.

## What your repository needs

The finished repository needs only this shape:

```text
your-repository/
├── README.md
└── exports/
    ├── datasets.json
    └── my-dataset/
        ├── dataset_identity.json
        ├── points_2d.bin.gz
        ├── obs_manifest.json
        ├── obs/
        ├── var_manifest.json
        └── var/
```

Do not create those data files by hand. The Cellucid Python package writes
them from your AnnData object.

## Step 1: prepare your dataset in Python

Install Cellucid in the Python environment used by your notebook:

```bash
python -m pip install "cellucid @ https://github.com/theislab/cellucid-python/archive/eedd3fca1dbb0f57a3ec9468c4a460003bda570a.zip"
```

This immutable revision is the verified Cellucid 0.9.1 source, including the
current weighted-connectivity writer. The archive pin keeps this tutorial
reproducible and does not require Git.

The example below starts from an H5AD file produced by a typical Scanpy
workflow. It assumes the file contains:

- `X_pca` in `adata.obsm`;
- a 2-D UMAP in `adata.obsm["X_umap"]`; and
- a neighbor graph in `adata.obsp["connectivities"]`.

In a Jupyter notebook, run:

```python
from pathlib import Path

import anndata as ad

from cellucid import prepare
from cellucid.prepare_data import generate_datasets_manifest

adata = ad.read_h5ad("my_data.h5ad")

dataset_id = "my-dataset"
exports_dir = Path("exports")

prepare(
    latent_space=adata.obsm["X_pca"],
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    connectivities=adata.obsp["connectivities"],
    X_umap_2d=adata.obsm["X_umap"],
    out_dir=exports_dir / dataset_id,
    dataset_id=dataset_id,
    dataset_name="My dataset",
    dataset_description="A short description of this dataset.",
    obs_categorical_dtype="uint16",
    compression=6,
    var_quantization=8,
    obs_continuous_quantization=8,
)

generate_datasets_manifest(
    exports_dir,
    default_dataset=dataset_id,
)
```

This creates:

```text
exports/
├── datasets.json
└── my-dataset/
    └── ...Cellucid data files...
```

The `dataset_id` is the stable name used in links. Use lowercase letters,
numbers, and hyphens, and do not change it after sharing the dataset.

If your AnnData uses different keys, inspect them first:

```python
print(list(adata.obsm.keys()))
print(list(adata.obsp.keys()))
print(list(adata.obs.columns))
```

Then replace the keys in the `prepare(...)` call with the exact keys in your
file. Do not invent a missing embedding or graph.

For a dataset without a connectivity graph, omit the `connectivities=...`
line. For 1-D or 3-D data, pass the real array as `X_umap_1d` or
`X_umap_3d`. The
[complete preparation guide](https://cellucid.readthedocs.io/en/latest/user_guide/python_package/c_data_preparation_api/index.html)
explains velocity fields, multiple dimensions, metadata selection, and large
datasets.

## Step 2: check the folder locally

Before uploading anything:

1. Open [Cellucid](https://www.cellucid.com).
2. Choose the prepared-folder input.
3. Select `exports/my-dataset`.
4. Confirm the cell count, embedding, metadata, and gene values.

This keeps GitHub setup separate from data-preparation problems.

## Step 3: create a public GitHub repository

The GitHub data source reads public files directly. A private repository will
not work with this public, token-free input.

1. Sign in to [GitHub](https://github.com).
2. Choose **New repository**.
3. Give it a clear name, such as `my-lab-cellucid-datasets`.
4. Select **Public**.
5. Create the repository.
6. Choose **Add file → Upload files**.
7. Upload your `README.md` and the complete `exports/` folder without changing
   its internal paths.
8. Commit the upload to the `main` branch.

After the upload, verify that GitHub shows:

```text
README.md
exports/datasets.json
exports/my-dataset/dataset_identity.json
```

## Step 4: open the GitHub data in Cellucid

In Cellucid:

1. Choose the **GitHub** data source.
2. Enter `OWNER/REPOSITORY/exports`.
3. Connect.
4. Choose `my-dataset`.

For example:

```text
my-lab/my-lab-cellucid-datasets/exports
```

The shareable link has this form:

```text
https://www.cellucid.com/?github=OWNER/REPOSITORY/exports&dataset=my-dataset
```

The short repository form uses the `main` branch. If you intentionally publish
from another branch, include it explicitly:

```text
OWNER/REPOSITORY@BRANCH/exports
```

## Step 5: turn on GitHub Pages

The GitHub input above works without Pages. Pages is useful when you also want
a stable static URL for the prepared files.

1. Open the repository on GitHub.
2. Choose **Settings**.
3. In the left sidebar, choose **Pages**.
4. Under **Build and deployment**, set **Source** to
   **Deploy from a branch**.
5. Choose the `main` branch and `/(root)`.
6. Choose **Save**.
7. Wait until GitHub reports that the site is live.

Your catalog should then be available at:

```text
https://OWNER.github.io/REPOSITORY/exports/datasets.json
```

Open that URL in a browser. You should see JSON, not a 404 page.

The corresponding direct Cellucid URL is:

```text
https://www.cellucid.com/?exportsBaseUrl=https%3A%2F%2FOWNER.github.io%2FREPOSITORY%2Fexports%2F&dataset=my-dataset
```

Replace `OWNER`, `REPOSITORY`, and `my-dataset` with your values.

## Validate this example repository

The maintenance workflow runs on Linux with Python 3.11 and 3.14, and on
macOS and Windows with Python 3.14. It verifies the catalog, identities,
manifest-declared artifacts, bounded gzip payloads, semantic binary contents,
and an exact checksum inventory without downloading data or installing Python
packages.

To run the same checks locally with Python 3.11 or newer:

```bash
python scripts/validate_exports.py
python -m unittest discover -s tests -p "test_*.py" -v
```

## Add more datasets

Repeat the complete Step 1 cell for the second AnnData object. Give it a new
`dataset_id`, a new `dataset_name`, and a new output such as
`Path("exports") / "second-dataset"`. After every dataset folder exists, call
`generate_datasets_manifest("exports", default_dataset="my-dataset")` once
more. Upload the updated `exports/` folder to GitHub. Cellucid will read the
new `datasets.json` and show both datasets.

## Update an existing dataset

Keep the same `dataset_id` when the dataset is still the same shared object.
To replace an existing local export intentionally, rerun the complete
`prepare(...)` call with `force=True`, regenerate `datasets.json`, inspect the
result locally, and upload the changed `exports/` folder.

Use a new `dataset_id` when cell order, gene identity, or biological meaning
changes enough that old saved sessions should not be applied to the new data.

## Privacy checklist

Before making the repository public:

- inspect every published `adata.obs` column;
- remove patient names, clinical identifiers, private sample IDs, and internal
  notes;
- confirm that gene and cell labels are safe to share;
- load the prepared folder locally and inspect what appears in Cellucid; and
- ask the person responsible for the data if public release is permitted.

If the data cannot be public, do not upload it to this type of repository. Use
the [Cellucid server workflow](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/04_server_tutorial.html)
on an appropriately protected system instead.

## Common problems

**Cellucid cannot find `datasets.json`**

Enter the folder that directly contains it. If the file is
`exports/datasets.json`, use `OWNER/REPOSITORY/exports`.

**The repository connects but no dataset appears**

Open `exports/datasets.json` on GitHub and confirm that each listed folder
contains `dataset_identity.json`.

**GitHub Pages returns 404**

Confirm Pages uses `main` and `/(root)`, wait for the deployment to finish, and
open the exact `/exports/datasets.json` URL.

**A file is too large for GitHub**

GitHub rejects ordinary Git files larger than 100 MiB. Use the
[Cellucid server workflow](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/04_server_tutorial.html)
instead of Git LFS pointers.

**The wrong dataset opens**

Include the exact dataset ID in the link:

```text
&dataset=my-dataset
```

For the full data-source map, see the
[Cellucid loading documentation](https://cellucid.readthedocs.io/en/latest/user_guide/web_app/b_data_loading/index.html).
