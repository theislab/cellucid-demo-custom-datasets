# Publish your own Cellucid datasets on GitHub

This is a deliberately artifact-only example repository. Cellucid reads AnnData
directly, so publishing a dataset on GitHub means committing a **Cellucid-ready
AnnData** — an ordinary AnnData that `cellucid.prepare()` (Python) or
`cellucid_prepare()` (R) laid out for the browser — and sharing one link. Its
committed top-level contents are exactly:

- `README.md` — this guide;
- `generate_datasets.py` — the single, self-contained reproducibility file;
- `catalog/` — a **dataset catalog**: `datasets.json` beside three tiny
  synthetic datasets, each a Zarr store folder (`.zarr/`), the form to publish
  on GitHub;
- `files/` — the 2-D dataset again as the two single files other hosts serve:
  an `.h5ad` and a Zarr ZIP (`.zarr.zip`);
- `.nojekyll` — an empty file that makes GitHub Pages publish the Zarr stores
  whole (see [Step 6](#step-6-optional-turn-on-github-pages));
- `.gitignore` and `.gitattributes` — the git metadata that keeps a checkout on
  any platform byte-identical to the tree `--check` compares against: the
  files inside a Zarr folder have no extension, so `.gitattributes` keeps them
  byte-exact (no line-ending conversion), and pins `datasets.json` to LF; and
- `LICENSE`, `CONTRIBUTING.md`, `SUPPORT.md`, `SECURITY.md` and
  `CODE_OF_CONDUCT.md` — how *this* repository is licensed, changed, asked
  about, and reported to. They govern this example; they are not part of the
  shape you copy. Your own repository needs its own licence and its own
  contact routes, not these.

Nothing else belongs in this repository: there are no tests, validation
harness, package scaffolding, extra scripts, workflows, `.github/` directory,
requirements file, or site assets. `generate_datasets.py` is the sole
maintenance script and regenerates the complete checked-in `catalog/` and
`files/` trees.
There is deliberately no `CITATION.cff`: these are synthetic teaching datasets,
not a result to cite. To cite Cellucid itself, see
[Citing Cellucid](https://docs.cellucid.com/en/latest/citing.html).

## Try this repository

Open [Cellucid](https://www.cellucid.com), expand **Session**, paste this into
the **Link or server** address bar, and press **Open**:

```text
theislab/cellucid-demo-custom-datasets/catalog
```

The catalog opens on its default dataset, and **Sample dataset** lists all
three small synthetic examples:

| Dataset id | Contents | Open |
| --- | --- | --- |
| `synthetic-cell-types-2d` | 72 cells, 8 genes; two 2-D embeddings of the same cells (`X_islands`, the default, and `X_marker_scores`), so the embedding picker has a choice | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog&dataset=synthetic-cell-types-2d) |
| `synthetic-development-3d` | 96 cells, 10 genes; a 3-D embedding `X_branches` with a velocity field `velocity_branches` drawn over it | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog&dataset=synthetic-development-3d) |
| `synthetic-trajectory-1d` | 48 cells, 6 genes; a 1-D embedding `X_trajectory` with a velocity field `velocity_trajectory` | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog&dataset=synthetic-trajectory-1d) |

Every dataset also has a weighted neighbour graph (`obsp['connectivities']`).
Its genes are indexed by a synthetic accession (`SYNC000001`, as a real
dataset is indexed by Ensembl ids) and named by their `SYN_…` symbols, which
the viewer shows and searches by default; **Session › Data mapping › Gene
names** switches to the accessions.

The velocity fields follow scVelo's naming: `obsm['velocity_<basis>']` is drawn
over `obsm['X_<basis>']`, for any basis — here `X_branches` and `X_trajectory`,
not a UMAP.

### One dataset, without the catalog

Each Zarr folder in `catalog/` opens on its own, and `files/` holds the 2-D
dataset as one file in each single-file form. Every form is read in place with
range requests, so only what a view needs is downloaded:

| Form | `?data=` value | Open |
| --- | --- | --- |
| Zarr folder (for GitHub) | `gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr` | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/catalog/synthetic-cell-types-2d.zarr) |
| `.h5ad` (for single-file hosts) | `gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad` | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.h5ad) |
| Zarr ZIP (one file) | `gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.zarr.zip` | [![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:theislab/cellucid-demo-custom-datasets/files/synthetic-cell-types-2d.zarr.zip) |

The same files from GitHub Pages, which this repository has turned on:

| What | `?data=` value |
| --- | --- |
| The catalog | [`https://theislab.github.io/cellucid-demo-custom-datasets/catalog/`](https://www.cellucid.com/?data=https%3A%2F%2Ftheislab.github.io%2Fcellucid-demo-custom-datasets%2Fcatalog%2F&dataset=synthetic-development-3d) |
| One Zarr folder | [`https://theislab.github.io/cellucid-demo-custom-datasets/catalog/synthetic-trajectory-1d.zarr/`](https://www.cellucid.com/?data=https%3A%2F%2Ftheislab.github.io%2Fcellucid-demo-custom-datasets%2Fcatalog%2Fsynthetic-trajectory-1d.zarr%2F) |

A `gh:` value without `@branch` reads the repository's default branch (GitHub
resolves `HEAD`, so `main`, `master` or any other works), and a link opens what
is on that branch now. A file's `raw.githubusercontent.com` address, or its
page on github.com, works as well as the `gh:` form.

These are teaching datasets, not biological results. They contain no patient,
donor, or clinical information.

## Cellucid documentation

The Cellucid documentation uses this repository as its public reference:

- [Publish a Custom Dataset Repository](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/11_custom_dataset_repository.html)
  is the complete catalog contract and publishing workflow.
- [Host a Dataset on GitHub](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/15_host_on_github.html)
  explains why a Zarr folder is the way to publish on GitHub, with measured
  numbers, the size limits, `.nojekyll`, and Git LFS as the alternative.
- [Host a Store for Sharing](https://docs.cellucid.com/en/latest/user_guide/r_package/d_viewing_loading/02_host_stores_for_sharing.html)
  is the same path from R.
- [Cellucid-ready AnnData specification](https://docs.cellucid.com/en/latest/user_guide/python_package/c_data_preparation_api/09_cellucid_ready_anndata_specification.html)
  is the exact contract the files here satisfy.
- [Troubleshoot Data Loading](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/08_troubleshooting_data_loading.html)
  uses this repository as the known-good connection when diagnosing GitHub
  paths.

## What your repository needs

One dataset needs only its Zarr folder. Several datasets, with a list the
viewer shows and one it opens on, need a catalog folder:

```text
your-repository/
├── README.md
├── .nojekyll                 # only if you turn on GitHub Pages
└── catalog/
    ├── datasets.json
    ├── my-dataset.zarr/
    └── second-dataset.zarr/
```

Do not write the Zarr folders by hand: `cellucid.prepare()` writes them from
your data.

## Step 1: write a Cellucid-ready AnnData

Any AnnData opens in Cellucid as it is. Writing it with `prepare()` makes it
read much faster over the internet (the expression matrix gene-major, arrays
chunked for range requests) and lets it carry its own name, description,
citation and launch defaults in `uns['cellucid']`.

Install Cellucid in the Python environment used by your notebook:

```bash
python -m pip install "cellucid==0.9.5"  # CELLUCID_VERSION
```

The version is pinned so this tutorial stays reproducible:
`0.9.5` is the version these example files were built with. <!-- CELLUCID_VERSION -->

The example below starts from an `.h5ad` produced by a typical Scanpy
workflow. It assumes the file contains:

- a 2-D UMAP in `adata.obsm["X_umap"]` and `X_pca` in `adata.obsm`;
- a neighbour graph in `adata.obsp["connectivities"]`; and
- gene ids in `adata.var_names` and symbols in `adata.var["gene_symbols"]`.

In a Jupyter notebook, run:

```python
import anndata as ad
import cellucid

adata = ad.read_h5ad("my_data.h5ad")

cellucid.prepare(
    "catalog/my-dataset.zarr",
    format="zarr",
    embeddings={"X_umap": adata.obsm["X_umap"], "X_pca": adata.obsm["X_pca"]},
    dataset_id="my-dataset",
    dataset_name="My dataset",
    dataset_description="A short description of this dataset.",
    obs=adata.obs,
    var=adata.var,
    gene_expression=adata.X,
    var_gene_id_column="gene_symbols",
    connectivities=adata.obsp["connectivities"],
    latent_key="X_pca",
)
```

This writes the folder `catalog/my-dataset.zarr/`: a standard AnnData Zarr
store that `anndata.read_zarr()` reads too. The `dataset_id` is the stable name
links, saved sessions and community annotations use. Use letters, digits,
`.`, `_` and `-`, and do not change it after sharing the dataset.

### Which `format`

`format` is required, because the right one depends on where the data goes:

| `format` | Writes | Use it for |
| --- | --- | --- |
| `"zarr"` | a `.zarr` folder of small files | **GitHub** (no file over 100 MB, no Git LFS; many small files at once over HTTP/2) |
| `"h5ad"` | one `.h5ad` file | everything else: single-file hosts (Zenodo, Figshare, Hugging Face, a bucket, a lab server), local files, `cellucid serve`, Jupyter |
| `"zarr.zip"` | the Zarr store as one `.zarr.zip` file | wherever an `.h5ad` goes; the single file cellucid-r writes |

From R, `cellucid_prepare("catalog/my-dataset.zarr", format = "zarr", ...)`
writes the same Zarr store (R writes `"zarr"` and `"zarr.zip"`).

### Embeddings and vector fields

`embeddings=` maps each `obsm` key to its array; every column is stored, and
the viewer plots any 1 to 3 columns of any of them.
`embedding_defaults={"2d": "X_umap"}` names the embedding a dimension opens on,
and `{"3d": ("X_pca", [0, 1, 2])}` also picks its columns. `latent_key` names the latent space
(scVI, PCA; never the 2D/3D layout) in which the viewer's outlier filter
measures each cell's distance to its category's centre; `prepare()` computes
those distances for every categorical obs column and stores them in the store
(`uns['cellucid']['outlierQuantiles']`, a few MB per column), so the filter
works at once instead of reading the whole latent space in the browser.

A velocity is stored under scVelo's name, `obsm['velocity_<basis>']`, and
drawn over `obsm['X_<basis>']`:

```python
cellucid.prepare(
    ...,
    embeddings={"X_umap": adata.obsm["X_umap"]},
    vector_fields={"velocity": {"X_umap": adata.obsm["velocity_umap"]}},
)
```

### Gene names your readers can search

`var_gene_id_column=` names the `var` column whose values the viewer shows and
searches. Point it at the symbols your readers type (`"gene_symbols"` above):
in an accession-indexed dataset left on the index, a reader searching `CD8A`
finds nothing. Every `var` column is stored, so a reader can switch to another
column under **Session › Data mapping › Gene names**.

Cellucid performs no symbol lookup and ships no mapping of any kind. If your
source is indexed by accession and you want symbols, resolve them yourself,
before calling `prepare()`, and pass the resulting column.

If your AnnData uses different keys, inspect them first:

```python
print(list(adata.obsm.keys()))
print(list(adata.obsp.keys()))
print(list(adata.obs.columns))
print(adata.var.head())
```

Then use the exact keys in your file. Do not invent a missing embedding or
graph; for a dataset without a neighbour graph, omit `connectivities=`. The
[complete preparation guide](https://docs.cellucid.com/en/latest/user_guide/python_package/c_data_preparation_api/index.html)
explains every argument.

## Step 2: write `datasets.json` (several datasets)

A catalog folder lists its datasets in `datasets.json`, a few lines written by
hand:

```json
{
  "schema": 2,
  "default": "my-dataset",
  "datasets": [
    { "id": "my-dataset", "path": "my-dataset.zarr/" },
    { "id": "second-dataset", "path": "second-dataset.zarr/" }
  ]
}
```

`schema` is exactly `2`; `default` is the id the catalog opens on; each `id` is
the dataset's `dataset_id`; each `path` is relative to the catalog folder — a
Zarr folder ends in `/`, and an `.h5ad` or `.zarr.zip` file can be listed too.
Names, descriptions and counts are not repeated here: Cellucid reads them from
each store. [`catalog/datasets.json`](catalog/datasets.json) is this
repository's.

## Step 3: check it locally

Before uploading anything:

1. Open [Cellucid](https://www.cellucid.com).
2. Expand **Session**, click **Open…** under **File or folder**, and choose
   **Zarr folder…** (in Safari, drag the folder onto the page instead: its
   folder dialog leaves out the store's hidden metadata files).
3. Select `catalog/my-dataset.zarr`.
4. Confirm the cell and gene counts, the embeddings, the metadata, a gene, and
   each vector field.

To check the whole catalog, serve the catalog folder with a static server that
sends CORS headers and answers range requests, and open
`https://www.cellucid.com/?data=http://127.0.0.1:8000/`. This keeps GitHub
setup separate from data problems.

## Step 4: publish on GitHub

Cellucid reads public GitHub repositories directly; a private repository will
not work with this public, token-free input.

1. Create a **public** repository on [GitHub](https://github.com), for example
   `my-lab-cellucid-datasets`.
2. Commit the catalog folder with Git from the command line (a Zarr folder is
   many small files, and the github.com upload page takes at most 100 files at
   a time):

   ```bash
   git add catalog .nojekyll
   git commit -m "Add datasets"
   git push
   ```

3. Check that the catalog is published: open
   `https://raw.githubusercontent.com/OWNER/REPOSITORY/HEAD/catalog/datasets.json`
   in a browser tab. It must show the JSON itself, not a 404 page.

## Step 5: open it in Cellucid

In **Session › Link or server**, enter `OWNER/REPOSITORY/catalog` and press
**Open**, then choose a dataset in the dataset list. For one Zarr folder, enter its path,
`OWNER/REPOSITORY/catalog/my-dataset.zarr`.

The shareable links have this form:

```text
https://www.cellucid.com/?data=gh:OWNER/REPOSITORY/catalog&dataset=my-dataset
https://www.cellucid.com/?data=gh:OWNER/REPOSITORY/catalog/my-dataset.zarr
```

To give a link a one-click badge in your README, wrap it in the
"Open in Cellucid" badge (the **Share** button beside **Dataset** in the
Session panel writes this snippet for you):

```markdown
[![Open in Cellucid](https://www.cellucid.com/assets/img/open-in-cellucid.svg)](https://www.cellucid.com/?data=gh:OWNER/REPOSITORY/catalog&dataset=my-dataset)
```

The short form reads the repository's default branch, whatever it is called.
To publish from another branch, tag or commit, put it after the repository
name:

```text
OWNER/REPOSITORY@BRANCH/catalog
```

## Step 6 (optional): turn on GitHub Pages

The links above work without Pages. Pages answers a first request faster than
`raw.githubusercontent.com` and gives the files a stable static address.

GitHub Pages runs Jekyll by default, which publishes no file whose name starts
with `.` and no folder whose name starts with `_` — so a Zarr store would be
served without its metadata (`.zgroup`, `.zattrs`, `.zmetadata`) and its cell
names (`obs/_index`). An empty `.nojekyll` file at the root of the published
branch turns Jekyll off; this repository has one.

1. Open the repository's **Settings › Pages**.
2. Under **Build and deployment**, set **Source** to **Deploy from a branch**.
3. Choose the default branch and `/ (root)`, and **Save**.
4. When GitHub reports that the site is live, open
   `https://OWNER.github.io/REPOSITORY/catalog/datasets.json` and
   `https://OWNER.github.io/REPOSITORY/catalog/my-dataset.zarr/.zmetadata`.
   Both must show JSON, not a 404 page.

A Pages address is linked URL-encoded, ending in `/`:

```text
https://www.cellucid.com/?data=https%3A%2F%2FOWNER.github.io%2FREPOSITORY%2Fcatalog%2F&dataset=my-dataset
```

## Or publish one file elsewhere

For a host that serves single files, write one file instead of a folder:

```python
cellucid.prepare("my-dataset.h5ad", format="h5ad", ...)        # Zenodo, Figshare, Hugging Face, a bucket
cellucid.prepare("my-dataset.zarr.zip", format="zarr.zip", ...)  # the Zarr store as one file
```

Upload the file and open its address the same way (on GitHub,
`OWNER/REPOSITORY/path/to/my-dataset.h5ad`). GitHub refuses files over 100 MB
unless they are stored with [Git LFS](https://git-lfs.com), whose transfers
count against the account's LFS quota; Cellucid reads an LFS file from GitHub's
media host, but on GitHub a Zarr folder needs no LFS at all.

## Regenerate these example datasets

[`generate_datasets.py`](generate_datasets.py) contains the complete synthetic
source data, exact dependency pins, generation logic, and safe replacement
workflow. It is the only Python file in this repository. Every dataset is
written by `cellucid.prepare()`: the three Zarr folders in `catalog/` with
`format="zarr"`, and the 2-D dataset's single files in `files/` with
`format="h5ad"` and `format="zarr.zip"`; it writes `catalog/datasets.json`
itself. Its `_gene_frame()` helper is the worked version of the gene-name
section above, and the 3-D and 1-D datasets pass
`vector_fields={"velocity": {"X_branches": velocity}}` and
`{"velocity": {"X_trajectory": velocity}}`. With
[uv](https://docs.astral.sh/uv/), its inline dependency metadata makes the
checked-in files reproducible without a requirements file. `prepare()` writes
through `h5py` and `zarr`, so `anndata`, `h5py` and `zarr` are pinned as well.

Rebuild into a temporary directory and compare every generated path and byte
with `catalog/` and `files/`, without changing the repository:

```bash
uv run generate_datasets.py --check
```

To intentionally renew the complete generation after editing the synthetic
inputs, use:

```bash
uv run generate_datasets.py --force
uv run generate_datasets.py --check
```

If the exact dependencies embedded at the top of the file are already
installed from their recorded sources, the equivalent commands are
`python generate_datasets.py --force` and
`python generate_datasets.py --check`.

Before touching `catalog/` or `files/`, the script refuses to run unless the
environment is the exact one the pins name: the Python version is in range,
`cellucid`, `numpy`, `pandas`, `scipy`, `anndata`, `h5py` and `zarr` are
installed at exactly the pinned versions, the `cellucid` metadata and runtime
versions agree, and the imported `cellucid` is the pinned distribution rather
than another copy shadowing it. If `cellucid` is an editable checkout instead
of a released install, it additionally requires that checkout to have no
uncommitted changes under `src/cellucid` or `pyproject.toml` — it does not pin
a commit, so a clean checkout of a different commit is accepted, and re-running
`--check` after changing `cellucid` is what catches that. `--check` also
requires `.nojekyll` to be present and empty, and reads every written dataset
back with `anndata` to confirm its id, counts and `obsm` keys. `--force` stages
and verifies a complete generation before replacing the old directories, and
restores them when an ordinary Python error interrupts that swap. Portable
directory replacement is not one filesystem-atomic operation: do not run a
local viewer concurrently, and if the process or machine is forcibly
terminated during the swap, restore the `.catalog-backup-*` and
`.files-backup-*` directories before retrying. There is intentionally no CI or
separate test/validation framework in this artifact-only repository; inspect
the renewed files in Cellucid before publishing them.

## Add more datasets

Run `prepare()` again for the second AnnData object, with a new `dataset_id`,
a new `dataset_name`, and a new output such as `catalog/second-dataset.zarr`.
Add its entry to `datasets.json`, commit both, and push. Cellucid reads the
new `datasets.json` and lists both datasets.

## Update an existing dataset

Keep the same `dataset_id` when the dataset is still the same shared object.
To replace a store you already wrote, rerun `prepare(...)` with `force=True`,
inspect the result locally, and commit the changed folder; Git removes the
chunk files that are gone when you `git add` the folder.

Use a new `dataset_id` when cell order, gene identity, or biological meaning
changes enough that old saved sessions should not be applied to the new data.

## Privacy checklist

Before making the repository public:

- inspect every published `adata.obs` column;
- remove patient names, clinical identifiers, private sample IDs, and internal
  notes;
- confirm that gene and cell labels are safe to share;
- open the written store locally and inspect what appears in Cellucid; and
- ask the person responsible for the data if public release is permitted.

If the data cannot be public, do not upload it to this type of repository. Use
the [Cellucid server workflow](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/04_server_tutorial.html)
on an appropriately protected system instead.

## Common problems

**Cellucid says the folder holds no `datasets.json`**

Enter the folder that directly contains it. If the file is
`catalog/datasets.json`, use `OWNER/REPOSITORY/catalog` — or link one Zarr
folder directly.

**`Invalid datasets.json: …`**

The message names the key or value: an unknown key, a missing `schema`,
`default` or `datasets` (a catalog must hold all three), a `schema` other than
`2`, a path that is not a Zarr folder (ending in `/`), `.h5ad` or `.zarr.zip`,
or a repeated id or path.

**A dataset is listed but does not open**

The message names the dataset and the reason: a missing file (404), an id that
differs from the store's `uns['cellucid']['id']`, or a folder that is not an
AnnData store.

**It opens from `raw.githubusercontent.com` but not from GitHub Pages**

Add the empty `.nojekyll` file at the root of the published branch, so Pages
publishes the store's dot-files and its `_index` folders.

**A file is too large for GitHub**

GitHub rejects ordinary Git files larger than 100 MiB. Write `format="zarr"`:
each chunk file of a Zarr folder is small. For a single file, use Git LFS or
another host.

**The wrong dataset opens**

Include the exact dataset ID in the link:

```text
&dataset=my-dataset
```

For the full data-source map, see the
[Cellucid loading documentation](https://docs.cellucid.com/en/latest/user_guide/web_app/b_data_loading/index.html).
