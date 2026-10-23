#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.15"
# dependencies = [
#   "cellucid==0.9.5",  # CELLUCID_VERSION
#   "numpy==2.5.1",
#   "pandas==2.3.3",
#   "scipy==1.18.0",
#   "anndata==0.12.19",
#   "h5py==3.16.0",
#   "zarr==3.3.0",
# ]
# ///
"""Generate the deterministic synthetic Cellucid datasets of this repository.

Every dataset is a Cellucid-ready AnnData written by ``cellucid.prepare()``:

- ``catalog/`` is a dataset catalog: ``datasets.json`` beside three Zarr store
  folders, the form to publish on GitHub (small files, no LFS);
- ``files/`` holds the 2-D dataset again as the two single files other hosts
  serve: an ``.h5ad`` and a ``.zarr.zip``.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import os
import shutil
import subprocess
import sys
import tempfile
import urllib.parse
import urllib.request
import uuid
from pathlib import Path
from typing import Any

import anndata
import cellucid
import numpy as np
import pandas as pd
from cellucid import prepare
from scipy import sparse

REPOSITORY_ROOT = Path(__file__).resolve().parent
SOURCE_NAME = "Cellucid deterministic synthetic examples"
SOURCE_URL = "https://github.com/theislab/cellucid-demo-custom-datasets"
# The Cellucid these files were built with. Keep in step with the inline
# dependency block above and with the install command in README.md.
CELLUCID_RELEASE = "0.9.5"  # CELLUCID_VERSION
# The paths within cellucid-python that determine what `prepare` writes. Only
# consulted for a maintainer's editable checkout; a released install has no
# working tree to be dirty.
PACKAGED_SOURCE_PATHS = ("src/cellucid", "pyproject.toml")
EXPECTED_VERSIONS = {
    "cellucid": CELLUCID_RELEASE,
    "numpy": "2.5.1",
    "pandas": "2.3.3",
    "scipy": "1.18.0",
    # prepare() writes through h5py and zarr-python, whose bytes are only
    # reproducible at fixed versions of them.
    "anndata": "0.12.19",
    "h5py": "3.16.0",
    "zarr": "3.3.0",
}

#: The catalog folder and the single-file folder, each replaced as a whole and
#: compared byte by byte.
CATALOG = "catalog"
FILES = "files"
GENERATED_DIRECTORIES = (CATALOG, FILES)
#: GitHub Pages runs Jekyll unless this empty file is at the site root, and
#: Jekyll publishes no file starting with `.` (.zgroup, .zattrs, .zmetadata)
#: and no folder starting with `_` (obs/_index): a Zarr store would be served
#: without its metadata and its cell names.
NOJEKYLL = ".nojekyll"

#: What each dataset must contain once written: its counts, its obsm keys with
#: their widths, and which of them are vector fields.
EXPECTED_DATASETS: dict[str, dict[str, Any]] = {
    "synthetic-cell-types-2d": {
        "n_cells": 72,
        "n_genes": 8,
        "obsm": {"X_islands": 2, "X_marker_scores": 2},
    },
    "synthetic-development-3d": {
        "n_cells": 96,
        "n_genes": 10,
        "obsm": {"X_branches": 3, "velocity_branches": 3},
    },
    "synthetic-trajectory-1d": {
        "n_cells": 48,
        "n_genes": 6,
        "obsm": {"X_trajectory": 1, "velocity_trajectory": 1},
    },
}
DEFAULT_DATASET = "synthetic-cell-types-2d"
# The dataset written again as single files, and the files it becomes.
SINGLE_FILE_DATASET = "synthetic-cell-types-2d"
SINGLE_FILES = {
    f"{SINGLE_FILE_DATASET}.h5ad": "h5ad",
    f"{SINGLE_FILE_DATASET}.zarr.zip": "zarr.zip",
}


def _require_cellucid_provenance(
    distribution: importlib.metadata.Distribution,
) -> None:
    """Require Cellucid to be the released distribution, or a clean checkout of it.

    The published recipe installs `cellucid` from PyPI at the exact version in
    EXPECTED_VERSIONS, which the caller has already checked. A release installed
    from an index records no ``direct_url.json``, so its absence is the expected
    case and is accepted.

    A maintainer working inside a clone of cellucid-python installs the package
    editable instead. That is accepted too, but only while the packaged sources
    are clean: uncommitted edits there change what `prepare` writes, and the
    files would no longer match the version this file names.
    """
    direct_url_text = distribution.read_text("direct_url.json")
    if direct_url_text is None:
        # Installed from an index. This is the published path.
        _require_distribution_import_ownership(distribution)
        return

    try:
        direct_url = json.loads(direct_url_text)
    except json.JSONDecodeError as error:
        raise RuntimeError(
            "The Cellucid installation records unreadable provenance; reinstall "
            f'it with `python -m pip install "cellucid=={CELLUCID_RELEASE}"`.'
        ) from error

    directory_info = (
        direct_url.get("dir_info") if isinstance(direct_url, dict) else None
    )
    url = direct_url.get("url") if isinstance(direct_url, dict) else None
    parsed_url = urllib.parse.urlparse(url) if isinstance(url, str) else None
    if (
        not isinstance(directory_info, dict)
        or directory_info.get("editable") is not True
        or parsed_url is None
        or parsed_url.scheme != "file"
    ):
        raise RuntimeError(
            "Cellucid must be the released distribution, installed with "
            f'`python -m pip install "cellucid=={CELLUCID_RELEASE}"`, or an '
            "editable checkout of it."
        )

    source_root = Path(
        urllib.request.url2pathname(urllib.parse.unquote(parsed_url.path))
    ).resolve()
    imported_package = Path(cellucid.__file__ or "").resolve()
    try:
        imported_package.relative_to(source_root / "src" / "cellucid")
    except ValueError as error:
        raise RuntimeError(
            "The editable Cellucid metadata does not own the imported package."
        ) from error

    try:
        modified_sources = subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "status",
                "--porcelain",
                "--untracked-files=all",
                "--",
                *PACKAGED_SOURCE_PATHS,
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise RuntimeError(
            "Cannot verify the editable Cellucid source checkout."
        ) from error
    if modified_sources:
        differing = ", ".join(
            sorted({line[3:] for line in modified_sources.splitlines() if line[3:]})
        )
        raise RuntimeError(
            "The editable Cellucid packaged sources have uncommitted changes: "
            f"{differing}. Commit or restore them, or install the released "
            f"cellucid=={CELLUCID_RELEASE}, before generating."
        )


def _require_distribution_import_ownership(
    distribution: importlib.metadata.Distribution,
) -> None:
    expected_package = Path(distribution.locate_file("cellucid")).resolve()
    imported_package = Path(cellucid.__file__ or "").resolve().parent
    if imported_package != expected_package:
        raise RuntimeError(
            "The imported Cellucid package is shadowing the pinned distribution: "
            f"imported={imported_package}, expected={expected_package}."
        )


def _require_exact_environment() -> None:
    if not (sys.version_info >= (3, 11) and sys.version_info < (3, 15)):
        raise RuntimeError(
            "Dataset generation requires Python >=3.11,<3.15; "
            f"found {sys.version.split()[0]}."
        )

    cellucid_distribution = importlib.metadata.distribution("cellucid")
    installed = {
        "cellucid": cellucid_distribution.version,
        **{
            name: importlib.metadata.version(name)
            for name in EXPECTED_VERSIONS
            if name != "cellucid"
        },
    }
    if installed != EXPECTED_VERSIONS:
        details = ", ".join(
            f"{name}={installed[name]} (required {version})"
            for name, version in EXPECTED_VERSIONS.items()
        )
        raise RuntimeError(f"Build environment is not exact: {details}")

    if cellucid.__version__ != cellucid_distribution.version:
        raise RuntimeError(
            "Cellucid package metadata and runtime version disagree: "
            f"metadata={cellucid_distribution.version}, "
            f"runtime={cellucid.__version__}."
        )
    _require_cellucid_provenance(cellucid_distribution)


def _categorical(values: list[str], categories: list[str]) -> pd.Categorical:
    return pd.Categorical(values, categories=categories, ordered=True)


def _graph(
    n_cells: int,
    edges: list[tuple[int, int, float]],
) -> sparse.csr_matrix:
    matrix = np.zeros((n_cells, n_cells), dtype=np.float32)
    for source, destination, weight in edges:
        if not 0 <= source < destination < n_cells:
            raise ValueError(
                "Synthetic connectivity edges require source < destination."
            )
        if matrix[source, destination] != 0:
            raise ValueError("Synthetic connectivity edges must be unique.")
        matrix[source, destination] = weight
        matrix[destination, source] = weight
    return sparse.csr_matrix(matrix)


def _gene_frame(symbols: list[str], *, prefix: str) -> pd.DataFrame:
    """The ``var`` frame of one dataset, laid out as a Scanpy workflow leaves it.

    The index holds one stable accession per gene (``SYNC000001``, the way a
    real dataset is indexed by Ensembl ids) and the ``symbol`` column holds the
    names a reader types. ``prepare(var_gene_id_column="symbol")`` makes the
    symbols the names the viewer shows and searches by default; every ``var``
    column is stored, so a reader can switch to the accessions under
    Session > Data mapping > Gene names.

    ``prepare()`` never invents a name: it performs no symbol lookup and ships
    no mapping of any kind, so resolving accessions to symbols -- if that is
    what your source needs -- belongs here, in the repository that owns the
    data, before ``prepare()`` is called.
    """
    if len(set(symbols)) != len(symbols):
        raise ValueError("Every gene needs one distinct name.")
    return pd.DataFrame(
        {"symbol": list(symbols)},
        index=pd.Index(
            [f"{prefix}{position + 1:06d}" for position in range(len(symbols))]
        ),
    )


def _prepare(
    out_path: Path,
    *,
    output_format: str,
    dataset_id: str,
    dataset_name: str,
    description: str,
    embeddings: dict[str, np.ndarray],
    embedding_defaults: dict[str, str],
    latent_key: str,
    obs: pd.DataFrame,
    genes: list[str],
    gene_id_prefix: str,
    expression: np.ndarray,
    connectivities: sparse.csr_matrix,
    vector_fields: dict[str, dict[str, np.ndarray]] | None = None,
) -> None:
    prepare(
        out_path,
        format=output_format,
        embeddings=embeddings,
        dataset_id=dataset_id,
        dataset_name=dataset_name,
        dataset_description=description,
        source_name=SOURCE_NAME,
        source_url=SOURCE_URL,
        obs=obs,
        var=_gene_frame(genes, prefix=gene_id_prefix),
        gene_expression=expression,
        # The var column whose values name the genes in the viewer: the
        # symbols a reader types, not the accessions in var.index.
        var_gene_id_column="symbol",
        connectivities=connectivities,
        vector_fields=vector_fields,
        embedding_defaults=embedding_defaults,
        latent_key=latent_key,
        obs_keys=list(obs.columns),
        centroid_min_points=6,
    )


def _trajectory_1d(out_path: Path, output_format: str) -> None:
    n_cells = 48
    index = np.arange(n_cells, dtype=np.float32)
    pseudotime = index / np.float32(n_cells - 1)
    wobble = ((np.arange(n_cells) % 4) - 1.5).astype(np.float32) / np.float32(80)
    embedding = (np.float32(2) * pseudotime - np.float32(1) + wobble)[:, None]

    stage_categories = ["early", "transition", "committed", "mature"]
    stages = [stage_categories[min(value // 12, 3)] for value in range(n_cells)]
    replicate_categories = ["replicate-a", "replicate-b"]
    replicates = [replicate_categories[value % 2] for value in range(n_cells)]
    transcriptional_activity = (
        np.float32(0.25)
        + np.float32(0.6) * pseudotime
        + (np.arange(n_cells) % 3).astype(np.float32) / np.float32(30)
    )
    obs = pd.DataFrame(
        {
            "stage": _categorical(stages, stage_categories),
            "replicate": _categorical(replicates, replicate_categories),
            "pseudotime": pseudotime,
            "transcriptional_activity": transcriptional_activity,
        },
        index=[f"trajectory-cell-{value:03d}" for value in range(n_cells)],
    )

    center_peak = np.float32(1) - np.abs(np.float32(2) * pseudotime - np.float32(1))
    saw = (np.arange(n_cells) % 8).astype(np.float32) / np.float32(7)
    expression = np.column_stack(
        [
            np.float32(6) * (np.float32(1) - pseudotime),
            np.float32(1) + np.float32(7) * center_peak,
            np.float32(6) * pseudotime,
            np.float32(2) + np.float32(2) * saw,
            np.float32(1) + transcriptional_activity,
            np.float32(0.5) + np.float32(3) * pseudotime * pseudotime,
        ]
    ).astype(np.float32)
    # The default name of each gene, in matrix-column order -- see _gene_frame().
    genes = [
        "SYN_EARLY",
        "SYN_TRANSITION",
        "SYN_MATURE",
        "SYN_SAWTOOTH",
        "SYN_ACTIVITY",
        "SYN_LATE",
    ]

    edges: list[tuple[int, int, float]] = []
    for source in range(n_cells - 1):
        edges.append((source, source + 1, 1.0))
    for source in range(n_cells - 2):
        edges.append((source, source + 2, 0.5))
    velocity = (
        np.float32(0.07) + (np.arange(n_cells) % 5).astype(np.float32) / np.float32(500)
    )[:, None]

    _prepare(
        out_path,
        output_format=output_format,
        dataset_id="synthetic-trajectory-1d",
        dataset_name="Synthetic trajectory — 1D",
        description=(
            "A synthetic 48-cell progression for learning planar 1D navigation, "
            "pseudotime coloring, gene trends, connectivity, and a velocity-style overlay."
        ),
        embeddings={"X_trajectory": embedding},
        embedding_defaults={"1d": "X_trajectory"},
        latent_key="X_trajectory",
        obs=obs,
        genes=genes,
        gene_id_prefix="SYNT",
        expression=expression,
        connectivities=_graph(n_cells, edges),
        # Stored as obsm['velocity_trajectory'], scVelo's name for a velocity
        # on the basis X_trajectory: the viewer draws it over that embedding.
        vector_fields={"velocity": {"X_trajectory": velocity}},
    )


def _cell_types_2d(out_path: Path, output_format: str) -> None:
    n_groups = 3
    cells_per_group = 24
    n_cells = n_groups * cells_per_group
    centers = np.asarray(
        [
            [-0.70, -0.35],
            [0.65, -0.25],
            [0.00, 0.68],
        ],
        dtype=np.float32,
    )
    offsets = []
    for local_index in range(cells_per_group):
        row = local_index // 6
        column = local_index % 6
        offsets.append(
            [
                np.float32(column - 2.5) / np.float32(18),
                np.float32(row - 1.5) / np.float32(12),
            ]
        )
    offsets_array = np.asarray(offsets, dtype=np.float32)
    embedding = np.vstack(
        [centers[group] + offsets_array for group in range(n_groups)]
    ).astype(np.float32)

    cell_type_categories = ["type-a", "type-b", "type-c"]
    cell_types = [
        cell_type_categories[index // cells_per_group] for index in range(n_cells)
    ]
    batch_categories = ["batch-1", "batch-2", "batch-3"]
    batches = [batch_categories[index % 3] for index in range(n_cells)]
    quality_score = np.float32(0.72) + (np.arange(n_cells) % 9).astype(
        np.float32
    ) / np.float32(50)
    library_size = np.float32(900) + (np.arange(n_cells) % 12).astype(
        np.float32
    ) * np.float32(35)
    obs = pd.DataFrame(
        {
            "cell_type": _categorical(cell_types, cell_type_categories),
            "batch": _categorical(batches, batch_categories),
            "quality_score": quality_score,
            "library_size": library_size,
        },
        index=[f"island-cell-{value:03d}" for value in range(n_cells)],
    )

    group_codes = np.repeat(np.arange(n_groups), cells_per_group)
    local = np.tile(np.arange(cells_per_group), n_groups).astype(np.float32)
    expression = np.column_stack(
        [
            np.where(group_codes == 0, np.float32(8), np.float32(1))
            + local / np.float32(48),
            np.where(group_codes == 1, np.float32(8), np.float32(1))
            + local / np.float32(48),
            np.where(group_codes == 2, np.float32(8), np.float32(1))
            + local / np.float32(48),
            np.float32(2) + (local % 6) / np.float32(3),
            np.float32(1) + quality_score,
            np.float32(1) + library_size / np.float32(500),
            np.float32(0.5) + embedding[:, 0] * embedding[:, 0],
            np.float32(0.5) + embedding[:, 1] * embedding[:, 1],
        ]
    ).astype(np.float32)
    # A second 2D layout of the same cells, read off the marker genes: the
    # horizontal axis contrasts marker A with marker B, the vertical one sets
    # marker C against both, and the shared and quality genes spread each
    # island out. With two 2D embeddings, `embedding_defaults` names the one
    # the dataset opens on; the viewer's embedding picker offers both.
    marker_scores = np.column_stack(
        [
            expression[:, 0] - expression[:, 1] + expression[:, 3] / np.float32(4),
            expression[:, 2]
            - (expression[:, 0] + expression[:, 1]) / np.float32(2)
            + expression[:, 4],
        ]
    ).astype(np.float32)
    # The default name of each gene, in matrix-column order -- see _gene_frame().
    genes = [
        "SYN_MARKER_A",
        "SYN_MARKER_B",
        "SYN_MARKER_C",
        "SYN_SHARED",
        "SYN_QUALITY",
        "SYN_LIBRARY",
        "SYN_AXIS_X",
        "SYN_AXIS_Y",
    ]

    edges: list[tuple[int, int, float]] = []
    for group in range(n_groups):
        start = group * cells_per_group
        for local_index in range(cells_per_group):
            current = start + local_index
            right = start + (local_index // 6) * 6 + ((local_index + 1) % 6)
            down = start + ((local_index + 6) % cells_per_group)
            for neighbor, weight in ((right, 1.0), (down, 0.5)):
                source, destination = sorted((current, neighbor))
                edge = (source, destination, weight)
                if edge not in edges:
                    edges.append(edge)

    _prepare(
        out_path,
        output_format=output_format,
        dataset_id="synthetic-cell-types-2d",
        dataset_name="Synthetic cell-type islands — 2D",
        description=(
            "Three synthetic 2D populations for learning categorical and continuous "
            "coloring, marker genes, filtering, highlighting, and weighted connectivity."
        ),
        embeddings={"X_islands": embedding, "X_marker_scores": marker_scores},
        embedding_defaults={"2d": "X_islands"},
        latent_key="X_islands",
        obs=obs,
        genes=genes,
        gene_id_prefix="SYNC",
        expression=expression,
        connectivities=_graph(n_cells, edges),
    )


def _development_3d(out_path: Path, output_format: str) -> None:
    cells_per_lineage = 48
    n_cells = cells_per_lineage * 2
    local_index = np.tile(np.arange(cells_per_lineage), 2).astype(np.float32)
    progress = local_index / np.float32(cells_per_lineage - 1)
    lineage_code = np.repeat(np.asarray([-1, 1], dtype=np.float32), cells_per_lineage)
    centered = np.float32(2) * progress - np.float32(1)
    embedding = np.column_stack(
        [
            centered,
            lineage_code * (np.float32(0.18) + np.float32(0.72) * progress),
            lineage_code * centered * centered + np.float32(0.12) * centered,
        ]
    ).astype(np.float32)

    lineage_categories = ["lineage-left", "lineage-right"]
    lineages = [
        lineage_categories[index // cells_per_lineage] for index in range(n_cells)
    ]
    timepoint_categories = ["t0", "t1", "t2", "t3"]
    timepoints = [
        timepoint_categories[min((index % cells_per_lineage) // 12, 3)]
        for index in range(n_cells)
    ]
    condition_categories = ["control", "stimulated"]
    conditions = [condition_categories[(index // 8) % 2] for index in range(n_cells)]
    response_score = (
        np.float32(0.2)
        + np.float32(0.65) * progress
        + np.where(lineage_code > 0, np.float32(0.08), np.float32(0))
    ).astype(np.float32)
    obs = pd.DataFrame(
        {
            "lineage": _categorical(lineages, lineage_categories),
            "timepoint": _categorical(timepoints, timepoint_categories),
            "condition": _categorical(conditions, condition_categories),
            "differentiation_score": progress,
            "response_score": response_score,
        },
        index=[f"development-cell-{value:03d}" for value in range(n_cells)],
    )

    left = (lineage_code < 0).astype(np.float32)
    right = (lineage_code > 0).astype(np.float32)
    triangular = np.float32(1) - np.abs(np.float32(2) * progress - np.float32(1))
    expression = np.column_stack(
        [
            np.float32(7) * (np.float32(1) - progress),
            np.float32(1) + np.float32(6) * triangular,
            np.float32(1) + np.float32(7) * progress,
            np.float32(1) + np.float32(6) * left * progress,
            np.float32(1) + np.float32(6) * right * progress,
            np.float32(1) + response_score,
            np.float32(1) + embedding[:, 0] * embedding[:, 0],
            np.float32(1) + np.abs(embedding[:, 1]),
            np.float32(1) + np.abs(embedding[:, 2]),
            np.float32(1) + (local_index % 7) / np.float32(3),
        ]
    ).astype(np.float32)
    # The default name of each gene, in matrix-column order -- see _gene_frame().
    genes = [
        "SYN_PROGENITOR",
        "SYN_COMMITMENT",
        "SYN_DIFFERENTIATED",
        "SYN_LEFT_FATE",
        "SYN_RIGHT_FATE",
        "SYN_RESPONSE",
        "SYN_CURVATURE_X",
        "SYN_BRANCH_Y",
        "SYN_BRANCH_Z",
        "SYN_LOCAL_PHASE",
    ]

    edges: list[tuple[int, int, float]] = []
    for lineage in range(2):
        start = lineage * cells_per_lineage
        for source in range(start, start + cells_per_lineage - 1):
            edges.append((source, source + 1, 1.0))
        for source in range(start, start + cells_per_lineage - 2):
            edges.append((source, source + 2, 0.5))
    edges.append((0, cells_per_lineage, 0.25))

    velocity = np.column_stack(
        [
            np.full(n_cells, np.float32(0.05), dtype=np.float32),
            lineage_code * np.float32(0.018),
            lineage_code * centered * np.float32(0.05) + np.float32(0.003),
        ]
    ).astype(np.float32)

    _prepare(
        out_path,
        output_format=output_format,
        dataset_id="synthetic-development-3d",
        dataset_name="Synthetic branching development — 3D",
        description=(
            "A synthetic two-lineage 3D progression for orbit navigation, temporal "
            "metadata, fate markers, graph context, and a velocity-style vector field."
        ),
        embeddings={"X_branches": embedding},
        embedding_defaults={"3d": "X_branches"},
        latent_key="X_branches",
        obs=obs,
        genes=genes,
        gene_id_prefix="SYND",
        expression=expression,
        connectivities=_graph(n_cells, edges),
        # Stored as obsm['velocity_branches'], scVelo's name for a velocity on
        # the basis X_branches: the viewer draws it over that embedding.
        vector_fields={"velocity": {"X_branches": velocity}},
    )


#: Each dataset's writer, in catalog order.
DATASET_WRITERS = {
    "synthetic-cell-types-2d": _cell_types_2d,
    "synthetic-development-3d": _development_3d,
    "synthetic-trajectory-1d": _trajectory_1d,
}


def _catalog_document() -> dict[str, Any]:
    """``datasets.json``: which datasets the folder holds, and where.

    Names, descriptions and counts are not repeated here: Cellucid reads them
    from each store's ``uns['cellucid']`` when it lists the catalog.
    """
    return {
        "schema": 2,
        "default": DEFAULT_DATASET,
        "datasets": [
            {"id": dataset_id, "path": f"{dataset_id}.zarr/"}
            for dataset_id in DATASET_WRITERS
        ],
    }


def _read_store(path: Path) -> anndata.AnnData:
    if path.suffix == ".h5ad":
        return anndata.read_h5ad(path)
    return anndata.read_zarr(path)


def _require_dataset(adata: anndata.AnnData, dataset_id: str, where: Path) -> None:
    """Read a written dataset back with anndata and check what it must contain."""
    expected = EXPECTED_DATASETS[dataset_id]
    meta = adata.uns.get("cellucid", {})
    widths = {key: int(value.shape[1]) for key, value in adata.obsm.items()}
    problems = []
    if meta.get("id") != dataset_id:
        problems.append(f"uns['cellucid']['id'] is {meta.get('id')!r}")
    if adata.n_obs != expected["n_cells"] or adata.n_vars != expected["n_genes"]:
        problems.append(f"shape is {adata.n_obs} x {adata.n_vars}")
    if widths != expected["obsm"]:
        problems.append(f"obsm is {widths}, not {expected['obsm']}")
    if list(adata.obsp) != ["connectivities"]:
        problems.append(f"obsp is {list(adata.obsp)}")
    if "symbol" not in adata.var:
        problems.append("var has no 'symbol' column")
    if problems:
        raise RuntimeError(f"{where} is not the expected dataset: {'; '.join(problems)}.")


def _check_generation(root: Path) -> None:
    catalog_root = root / CATALOG
    expected_entries = {"datasets.json", *(f"{key}.zarr" for key in DATASET_WRITERS)}
    entries = {path.name for path in catalog_root.iterdir()}
    if entries != expected_entries:
        raise RuntimeError(
            f"{CATALOG}/ holds {sorted(entries)}, not {sorted(expected_entries)}."
        )
    catalog = json.loads((catalog_root / "datasets.json").read_text(encoding="utf-8"))
    if catalog != _catalog_document():
        raise RuntimeError(f"{CATALOG}/datasets.json differs from the catalog it lists.")
    for dataset_id in DATASET_WRITERS:
        store = catalog_root / f"{dataset_id}.zarr"
        _require_dataset(_read_store(store), dataset_id, store)

    files_root = root / FILES
    names = sorted(path.name for path in files_root.iterdir())
    if names != sorted(SINGLE_FILES):
        raise RuntimeError(f"{FILES}/ holds {names}, not {sorted(SINGLE_FILES)}.")
    single = files_root / f"{SINGLE_FILE_DATASET}.h5ad"
    _require_dataset(_read_store(single), SINGLE_FILE_DATASET, single)
    with tempfile.TemporaryDirectory(prefix="cellucid-unzip-") as temporary:
        shutil.unpack_archive(
            files_root / f"{SINGLE_FILE_DATASET}.zarr.zip", temporary, format="zip"
        )
        _require_dataset(
            anndata.read_zarr(temporary), SINGLE_FILE_DATASET, files_root / "*.zarr.zip"
        )


def build_generation(root: Path) -> None:
    """Write ``catalog/`` and ``files/`` under ``root`` with ``cellucid.prepare()``."""
    catalog_root = root / CATALOG
    catalog_root.mkdir(parents=True, exist_ok=False)
    for dataset_id, write in DATASET_WRITERS.items():
        write(catalog_root / f"{dataset_id}.zarr", "zarr")
    (catalog_root / "datasets.json").write_text(
        json.dumps(_catalog_document(), indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    files_root = root / FILES
    files_root.mkdir(parents=True, exist_ok=False)
    for name, output_format in SINGLE_FILES.items():
        DATASET_WRITERS[SINGLE_FILE_DATASET](files_root / name, output_format)
    _check_generation(root)


def _tree_entries(root: Path) -> dict[str, tuple[str, bytes]]:
    entries: dict[str, tuple[str, bytes]] = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        if path.is_symlink():
            raise RuntimeError(f"Generated entries must not be symbolic links: {relative}")
        if path.is_dir():
            entries[relative] = ("directory", b"")
        elif path.is_file():
            entries[relative] = ("file", path.read_bytes())
        else:
            raise RuntimeError(
                f"Generated entries must be regular files or directories: {relative}"
            )
    return entries


def _require_nojekyll() -> None:
    marker = REPOSITORY_ROOT / NOJEKYLL
    if marker.is_symlink() or not marker.is_file() or marker.stat().st_size != 0:
        raise RuntimeError(
            f"{marker} must be an empty file: without it GitHub Pages publishes "
            "no Zarr metadata (.zgroup, .zattrs, .zmetadata) and no obs/_index."
        )


def _check_reproducibility() -> None:
    _require_nojekyll()
    for name in GENERATED_DIRECTORIES:
        expected_root = REPOSITORY_ROOT / name
        if expected_root.is_symlink() or not expected_root.is_dir():
            raise FileNotFoundError(
                f"Checked-in {name}/ must be a real directory: {expected_root}"
            )
    _check_generation(REPOSITORY_ROOT)
    with tempfile.TemporaryDirectory(prefix="cellucid-rebuild-") as temporary:
        rebuilt = Path(temporary)
        build_generation(rebuilt)
        for name in GENERATED_DIRECTORIES:
            expected = _tree_entries(REPOSITORY_ROOT / name)
            actual = _tree_entries(rebuilt / name)
            if expected != actual:
                missing = sorted(expected.keys() - actual.keys())
                extra = sorted(actual.keys() - expected.keys())
                changed = sorted(
                    path
                    for path in expected.keys() & actual.keys()
                    if expected[path] != actual[path]
                )
                raise RuntimeError(
                    f"Checked-in {name}/ is not byte-reproducible: "
                    f"missing={missing}, extra={extra}, changed={changed}"
                )
    print(
        f"PASS: checked-in {CATALOG}/ and {FILES}/ are byte-identical to a clean rebuild"
    )


def _publish_generation() -> None:
    for name in GENERATED_DIRECTORIES:
        target = REPOSITORY_ROOT / name
        if target.is_symlink():
            raise ValueError(f"{name}/ must not be a symbolic link: {target}")
        if target.exists() and not target.is_dir():
            raise NotADirectoryError(f"{name}/ is not a directory: {target}")

    stage = Path(tempfile.mkdtemp(prefix=".generation-stage-", dir=REPOSITORY_ROOT))
    backups: dict[str, Path] = {}
    placed: list[str] = []
    try:
        build_generation(stage)
        # Move every old directory aside first, so a failure restores them all.
        try:
            for name in GENERATED_DIRECTORIES:
                target = REPOSITORY_ROOT / name
                if target.exists():
                    backups[name] = REPOSITORY_ROOT / f".{name}-backup-{uuid.uuid4().hex}"
                    os.replace(target, backups[name])
            for name in GENERATED_DIRECTORIES:
                os.replace(stage / name, REPOSITORY_ROOT / name)
                placed.append(name)
        except BaseException:
            for name in placed:
                shutil.rmtree(REPOSITORY_ROOT / name)
            for name, backup in backups.items():
                os.replace(backup, REPOSITORY_ROOT / name)
            backups.clear()
            raise
        for backup in backups.values():
            shutil.rmtree(backup)
        backups.clear()
    finally:
        shutil.rmtree(stage, ignore_errors=True)
        leftover = [str(backup) for backup in backups.values() if backup.exists()]
        if leftover:
            raise RuntimeError(
                f"Published generation is complete, but old backup cleanup failed: {leftover}"
            )


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build or verify the deterministic synthetic Cellucid datasets."
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--force",
        action="store_true",
        help=f"Replace this repository's {CATALOG}/ and {FILES}/ from a clean staged generation.",
    )
    action.add_argument(
        "--check",
        action="store_true",
        help="Rebuild in a temporary directory and compare every checked-in byte.",
    )
    return parser.parse_args()


def main() -> int:
    arguments = _parse_arguments()
    _require_exact_environment()
    if arguments.check:
        _check_reproducibility()
    else:
        _publish_generation()
        print(
            f"PASS: published deterministic {CATALOG}/ and {FILES}/ to {REPOSITORY_ROOT}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
