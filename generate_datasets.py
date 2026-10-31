#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.15"
# dependencies = [
#   "cellucid==0.9.6",  # CELLUCID_VERSION
#   "numpy==2.5.1",
#   "pandas==2.3.3",
#   "scipy==1.18.0",
#   "anndata==0.12.19",
#   "h5py==3.16.0",
#   "zarr==3.3.0",
#   "mudata==0.4.1",
# ]
# ///
"""Generate the deterministic synthetic Cellucid datasets of this repository.

Every AnnData and MuData here is Cellucid-ready, written by
``cellucid.prepare()``:

- ``catalog/`` is a dataset catalog: ``datasets.json`` beside three AnnData
  Zarr store folders, a multimodal (CITE-seq) MuData Zarr store folder, and
  the same CITE-seq cells as a Seurat object saved as R's ``saveRDS()`` does
  (``.rds``), the forms to publish on GitHub (small files, no LFS);
- ``files/`` holds the 2-D dataset again as the two single files other hosts
  serve: an ``.h5ad`` and a ``.zarr.zip``.
"""

from __future__ import annotations

import argparse
import gzip
import importlib.metadata
import json
import os
import shutil
import struct
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
CELLUCID_RELEASE = "0.9.6"  # CELLUCID_VERSION
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
    # The MuData container prepare(data=mdata) is given and the store is read
    # back with.
    "mudata": "0.4.1",
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
#: The multimodal (CITE-seq) dataset: the same cells as a Cellucid-ready
#: MuData Zarr store and as a Seurat object in an .rds file.
CITESEQ_MUDATA = "synthetic-citeseq"
CITESEQ_SEURAT = "synthetic-citeseq-seurat"
EXPECTED_MODALITIES = {"rna": 8, "prot": 4}
CITESEQ_CELLS = 60
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


def _citeseq_cells() -> dict[str, Any]:
    """The one synthetic CITE-seq experiment both multimodal files hold.

    Sixty cells of three populations, each measured twice: eight genes (RNA)
    and four surface proteins (antibody-derived tags, ADT), as integer counts.
    The MuData Zarr store in ``catalog/`` and the Seurat object written as an
    ``.rds`` are two encodings of these same counts, embeddings and graph.
    """
    n_groups = 3
    cells_per_group = 20
    n_cells = n_groups * cells_per_group
    group = np.repeat(np.arange(n_groups), cells_per_group)
    local = np.tile(np.arange(cells_per_group), n_groups)
    centers = np.asarray(
        [[-0.70, -0.40], [0.70, -0.40], [0.00, 0.70]], dtype=np.float32
    )
    offsets = np.column_stack(
        [
            ((local % 5) - 2).astype(np.float32) / np.float32(14),
            ((local // 5) - 1.5).astype(np.float32) / np.float32(10),
        ]
    ).astype(np.float32)
    umap = (centers[group] + offsets).astype(np.float32)

    in_group = [(group == index).astype(np.int32) for index in range(n_groups)]
    rna_counts = np.column_stack(
        [
            1 + 11 * in_group[0] + local % 3,
            1 + 11 * in_group[1] + local % 4,
            12 * in_group[2] + local % 2,
            4 + local % 5,
            np.where(local % 6 == 0, 6, 0),
            8 + local % 3,
            np.where(local % 4 == 1, 3, 0),
            2 + (local + group) % 4,
        ]
    ).astype(np.int32)
    # Seurat refuses "_" in feature names (it turns them into "-"), so these
    # names use "-" in both files.
    genes = [
        "SYN-B-MARKER",
        "SYN-T-MARKER",
        "SYN-MONO-MARKER",
        "SYN-SHARED",
        "SYN-CYCLING",
        "SYN-HOUSEKEEPING",
        "SYN-STRESS",
        "SYN-RIBOSOMAL",
    ]
    adt_counts = np.column_stack(
        [
            5 + 55 * in_group[0] + local % 4,
            4 + 56 * in_group[1] + local % 3,
            3 + 47 * in_group[2] + local % 5,
            local % 3,
        ]
    ).astype(np.int32)
    proteins = ["CD19-ADT", "CD3-ADT", "CD14-ADT", "Isotype-ADT"]

    # Seurat's two defaults: LogNormalize (counts per 10,000, log1p) for RNA,
    # and CLR across a cell's proteins for ADT.
    rna_totals = rna_counts.sum(axis=1, keepdims=True).astype(np.float64)
    rna_data = np.log1p(rna_counts / rna_totals * 1e4)
    adt_log = np.log1p(adt_counts.astype(np.float64))
    adt_data = np.log1p(
        adt_counts / np.exp(adt_log.sum(axis=1, keepdims=True) / adt_counts.shape[1])
    )
    pca = np.column_stack(
        [
            rna_data[:, 0] - rna_data[:, 1],
            rna_data[:, 2] - (rna_data[:, 0] + rna_data[:, 1]) / 2,
            rna_data[:, 3] - rna_data[:, 5],
            rna_data[:, 4] + rna_data[:, 6] / 2,
        ]
    ).astype(np.float32)

    cell_types = ["B cell", "T cell", "monocyte"]
    samples = ["donor-1", "donor-2"]
    obs = pd.DataFrame(
        {
            "cell_type": _categorical(
                [cell_types[value] for value in group], cell_types
            ),
            "sample": _categorical([samples[value % 2] for value in local], samples),
        },
        index=[f"citeseq-cell-{value:03d}" for value in range(n_cells)],
    )

    edges: list[tuple[int, int, float]] = []
    for index in range(n_groups):
        start = index * cells_per_group
        for position in range(cells_per_group):
            current = start + position
            right = start + (position // 5) * 5 + ((position + 1) % 5)
            down = start + ((position + 5) % cells_per_group)
            for neighbor, weight in ((right, 1.0), (down, 0.5)):
                source, destination = sorted((current, neighbor))
                if (source, destination, weight) not in edges:
                    edges.append((source, destination, weight))
    return {
        "obs": obs,
        "umap": umap,
        "pca": pca,
        "genes": genes,
        "rna_counts": rna_counts,
        "rna_data": rna_data,
        "proteins": proteins,
        "adt_counts": adt_counts,
        "adt_data": adt_data,
        "graph": _graph(n_cells, edges),
    }


def _citeseq_mudata(out_path: Path, output_format: str) -> None:
    """The CITE-seq experiment as a Cellucid-ready MuData: ``prepare(data=mdata)``."""
    import mudata

    cells = _citeseq_cells()
    obs_names = cells["obs"].index
    rna = anndata.AnnData(
        X=sparse.csr_matrix(cells["rna_data"].astype(np.float32)),
        obs=pd.DataFrame(index=obs_names),
        var=_gene_frame(cells["genes"], prefix="SYNM"),
    )
    prot = anndata.AnnData(
        X=sparse.csr_matrix(cells["adt_data"].astype(np.float32)),
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=pd.Index(cells["proteins"])),
    )
    mdata = mudata.MuData({"rna": rna, "prot": prot})
    for column in cells["obs"].columns:
        mdata.obs[column] = cells["obs"][column].to_numpy()
    mdata.obsm["X_umap"] = cells["umap"]
    mdata.obsm["X_pca"] = cells["pca"]
    mdata.obsp["connectivities"] = cells["graph"]
    prepare(
        out_path,
        format=output_format,
        data=mdata,
        dataset_id=CITESEQ_MUDATA,
        dataset_name="Synthetic CITE-seq — MuData",
        dataset_description=(
            "Sixty synthetic cells of three populations measured twice, as RNA (8 genes) "
            "and as surface proteins (4 antibody tags): two modalities to switch between "
            "under Data mapping."
        ),
        source_name=SOURCE_NAME,
        source_url=SOURCE_URL,
        modality="rna",
        var_gene_id_column="symbol",
        embedding_defaults={"2d": "X_umap"},
        latent_key="X_pca",
        obs_keys=list(cells["obs"].columns),
        centroid_min_points=6,
    )


# --- A Seurat object as R's saveRDS() writes it ------------------------------
#
# The same experiment as a Seurat v5 object (SeuratObject 5.4.0: two Assay5
# assays, a PCA and a UMAP, an SNN graph), serialized here as R's own
# saveRDS() would: R serialization format 3 (XDR, big-endian), gzip. Writing it
# from Python keeps this file the one generator, with no R to install; R reads
# it back with readRDS() as an ordinary Seurat object.

_NILVALUE_SXP = 254
_REFSXP = 255
_SYMSXP = 1
_LISTSXP = 2
_CHARSXP = 9
_SXP_TYPES = {"lgl": 10, "int": 13, "dbl": 14, "str": 16, "list": 19, "s4": 25}
#: The R and SeuratObject releases whose objects this file reproduces.
_R_VERSION = (4, 5, 3)
_SEURATOBJECT_VERSION = (5, 4, 0)


class _RValue:
    """One R value: a vector (``lgl``, ``int``, ``dbl``, ``str``, ``list``) or an S4 object."""

    def __init__(
        self,
        kind: str,
        data: Any = (),
        attributes: list[tuple[str, _RValue]] | None = None,
        *,
        s4: bool = False,
    ) -> None:
        self.kind = kind
        self.data = data
        self.attributes = attributes or []
        self.s4 = s4


def _r_str(*values: str) -> _RValue:
    return _RValue("str", list(values))


def _r_list(
    items: dict[str, _RValue] | None = None,
    extra: list[tuple[str, _RValue]] | None = None,
) -> _RValue:
    items = items or {}
    attributes = [("names", _r_str(*items))] if items else []
    return _RValue("list", list(items.values()), attributes + (extra or []))


def _r_class(name: str, package: str) -> _RValue:
    return _RValue("str", [name], [("package", _r_str(package))])


def _r_s4(name: str, package: str, slots: list[tuple[str, _RValue]]) -> _RValue:
    return _RValue("s4", (), [*slots, ("class", _r_class(name, package))], s4=True)


def _r_factor(values: pd.Categorical, names: list[str] | None = None) -> _RValue:
    attributes = [
        ("levels", _r_str(*map(str, values.categories))),
        ("class", _r_str("factor")),
    ]
    if names is not None:
        attributes.insert(0, ("names", _r_str(*names)))
    return _RValue("int", [int(code) + 1 for code in values.codes], attributes)


def _r_matrix(values: np.ndarray, dimnames: _RValue | None = None) -> _RValue:
    """A dense double matrix (R stores it column by column)."""
    rows, columns = values.shape
    attributes = [("dim", _RValue("int", [rows, columns]))]
    if dimnames is not None:
        attributes.append(("dimnames", dimnames))
    return _RValue("dbl", [float(value) for value in values.T.reshape(-1)], attributes)


def _r_empty_matrix() -> _RValue:
    return _RValue("dbl", [], [("dim", _RValue("int", [0, 0]))])


def _r_null() -> _RValue:
    return _RValue("null")


def _r_dgcmatrix(
    features_by_cells: np.ndarray, *, dimnames: _RValue
) -> list[tuple[str, _RValue]]:
    """The slots of a Matrix ``dgCMatrix``: column-compressed, one column per cell."""
    compressed = sparse.csc_matrix(features_by_cells)
    compressed.sort_indices()
    return [
        ("i", _RValue("int", compressed.indices.tolist())),
        ("p", _RValue("int", compressed.indptr.tolist())),
        ("Dim", _RValue("int", list(features_by_cells.shape))),
        ("Dimnames", dimnames),
        ("x", _RValue("dbl", compressed.data.astype(np.float64).tolist())),
        ("factors", _r_list()),
    ]


def _r_logmap(row_names: list[str], layers: list[str]) -> _RValue:
    """SeuratObject's ``LogMap``: which rows (cells or features) each layer holds."""
    return _RValue(
        "lgl",
        [1] * (len(row_names) * len(layers)),
        [
            ("dim", _RValue("int", [len(row_names), len(layers)])),
            ("class", _r_class("LogMap", "SeuratObject")),
            ("dimnames", _RValue("list", [_r_str(*row_names), _r_str(*layers)])),
        ],
        s4=True,
    )


def _r_data_frame(columns: dict[str, _RValue], row_names: _RValue) -> _RValue:
    return _RValue(
        "list",
        list(columns.values()),
        [
            ("row.names", row_names),
            ("names", _r_str(*columns)),
            ("class", _r_str("data.frame")),
        ],
    )


def _r_assay5(
    counts: np.ndarray,
    data: np.ndarray,
    *,
    cells: list[str],
    features: list[str],
    key: str,
) -> _RValue:
    no_names = _RValue("list", [_r_null(), _r_null()])
    layers = _r_list(
        {
            "counts": _r_s4(
                "dgCMatrix", "Matrix", _r_dgcmatrix(counts.T, dimnames=no_names)
            ),
            "data": _r_s4(
                "dgCMatrix", "Matrix", _r_dgcmatrix(data.T, dimnames=no_names)
            ),
        }
    )
    return _r_s4(
        "Assay5",
        "SeuratObject",
        [
            ("layers", layers),
            ("cells", _r_logmap(cells, ["counts", "data"])),
            ("features", _r_logmap(features, ["counts", "data"])),
            ("default", _RValue("int", [1])),
            ("assay.orig", _r_str()),
            (
                "meta.data",
                _RValue(
                    "list",
                    [],
                    [
                        ("names", _r_str()),
                        ("class", _r_str("data.frame")),
                        ("row.names", _RValue("int", [-(2**31), -len(features)])),
                    ],
                ),
            ),
            # Setting a layer leaves an Assay5's misc a named empty list.
            ("misc", _RValue("list", [], [("names", _r_str())])),
            ("key", _r_str(key)),
        ],
    )


def _r_dimreduc(embedding: np.ndarray, *, cells: list[str], key: str) -> _RValue:
    columns = [f"{key}{index + 1}" for index in range(embedding.shape[1])]
    empty = _r_empty_matrix
    return _r_s4(
        "DimReduc",
        "SeuratObject",
        [
            (
                "cell.embeddings",
                _r_matrix(
                    embedding.astype(np.float64),
                    _RValue("list", [_r_str(*cells), _r_str(*columns)]),
                ),
            ),
            ("feature.loadings", empty()),
            ("feature.loadings.projected", empty()),
            ("assay.used", _r_str("RNA")),
            ("global", _RValue("lgl", [0])),
            ("stdev", _RValue("dbl", [])),
            (
                "jackstraw",
                _r_s4(
                    "JackStrawData",
                    "SeuratObject",
                    [
                        ("empirical.p.values", empty()),
                        ("fake.reduction.scores", empty()),
                        ("empirical.p.values.full", empty()),
                        ("overall.p.values", empty()),
                    ],
                ),
            ),
            ("misc", _r_list()),
            ("key", _r_str(key)),
        ],
    )


def _citeseq_seurat_object() -> _RValue:
    """The Seurat object, slot for slot as SeuratObject 5.4.0 builds it in R 4.5.3.

    ``serialize(seu, NULL, version = 3)`` of the object these R lines build
    from the same values is byte for byte what this writes (before gzip)::

        seu <- CreateSeuratObject(counts = rna_counts, assay = "RNA",
                                  project = "synthetic-citeseq-seurat")
        LayerData(seu, assay = "RNA", layer = "data") <- rna_data
        seu[["ADT"]] <- CreateAssay5Object(counts = adt_counts)
        LayerData(seu, assay = "ADT", layer = "data") <- adt_data
        seu$cell_type <- cell_type; seu$sample <- sample
        Idents(seu) <- "cell_type"
        seu[["pca"]] <- CreateDimReducObject(pca, key = "PC_", assay = "RNA")
        seu[["umap"]] <- CreateDimReducObject(umap, key = "umap_", assay = "RNA")
        seu[["RNA_snn"]] <- graph  # a Graph, DefaultAssay "RNA"
        seu@misc$cellucid <- list(id = ..., name = ..., description = ...,
                                  source = list(...), defaults = list(...))
    """
    cells = _citeseq_cells()
    obs = cells["obs"]
    names = list(obs.index)
    project = CITESEQ_SEURAT
    rna_counts, adt_counts = cells["rna_counts"], cells["adt_counts"]
    meta = _r_data_frame(
        {
            "orig.ident": _r_factor(pd.Categorical([project] * len(names))),
            "nCount_RNA": _RValue(
                "dbl", rna_counts.sum(axis=1).astype(np.float64).tolist()
            ),
            "nFeature_RNA": _RValue("int", (rna_counts > 0).sum(axis=1).tolist()),
            "nCount_ADT": _RValue(
                "dbl", adt_counts.sum(axis=1).astype(np.float64).tolist()
            ),
            "nFeature_ADT": _RValue("int", (adt_counts > 0).sum(axis=1).tolist()),
            "cell_type": _r_factor(obs["cell_type"].array),
            "sample": _r_factor(obs["sample"].array),
        },
        _r_str(*names),
    )
    graph_slots = _r_dgcmatrix(
        cells["graph"].toarray().astype(np.float64),
        dimnames=_RValue("list", [_r_str(*names), _r_str(*names)]),
    )
    cellucid_record = _r_list(
        {
            "id": _r_str(CITESEQ_SEURAT),
            "name": _r_str("Synthetic CITE-seq — Seurat object"),
            "description": _r_str(
                "The synthetic CITE-seq cells as a Seurat object saved with saveRDS(): "
                "RNA and ADT assays, a PCA, a UMAP and an SNN graph."
            ),
            "source": _r_list({"name": _r_str(SOURCE_NAME), "url": _r_str(SOURCE_URL)}),
            "defaults": _r_list(
                {
                    "embedding2d": _r_str("X_umap"),
                    "obsKeys": _r_str("cell_type", "sample"),
                }
            ),
        }
    )
    return _r_s4(
        "Seurat",
        "SeuratObject",
        [
            (
                "assays",
                _r_list(
                    {
                        "RNA": _r_assay5(
                            rna_counts,
                            cells["rna_data"],
                            cells=names,
                            features=cells["genes"],
                            key="rna_",
                        ),
                        "ADT": _r_assay5(
                            adt_counts,
                            cells["adt_data"],
                            cells=names,
                            features=cells["proteins"],
                            key="adt_",
                        ),
                    }
                ),
            ),
            ("meta.data", meta),
            ("active.assay", _r_str("RNA")),
            ("active.ident", _r_factor(obs["cell_type"].array, names)),
            (
                "graphs",
                _r_list(
                    {
                        "RNA_snn": _r_s4(
                            "Graph",
                            "SeuratObject",
                            [("assay.used", _r_str("RNA")), *graph_slots],
                        )
                    }
                ),
            ),
            ("neighbors", _r_list()),
            (
                "reductions",
                _r_list(
                    {
                        "pca": _r_dimreduc(cells["pca"], cells=names, key="PC_"),
                        "umap": _r_dimreduc(cells["umap"], cells=names, key="umap_"),
                    }
                ),
            ),
            ("images", _r_list()),
            ("project.name", _r_str(project)),
            ("misc", _r_list({"cellucid": cellucid_record})),
            (
                "version",
                _RValue(
                    "list",
                    [_RValue("int", list(_SEURATOBJECT_VERSION))],
                    [("class", _r_str("package_version", "numeric_version"))],
                ),
            ),
            ("commands", _r_list()),
            ("tools", _r_list()),
        ],
    )


class _RdsWriter:
    """R serialization format 3, XDR: what ``serialize(x, NULL, version = 3)`` writes."""

    def __init__(self) -> None:
        self.out = bytearray()
        self.symbols: dict[str, int] = {}

    def integer(self, value: int) -> None:
        self.out += struct.pack(">i", value)

    def header(self) -> None:
        major, minor, patch = _R_VERSION
        self.out += b"X\n"
        self.integer(3)
        self.integer(major * 65536 + minor * 256 + patch)
        self.integer(3 * 65536 + 5 * 256)  # R 3.5.0 reads format 3
        self.integer(len(b"UTF-8"))
        self.out += b"UTF-8"

    def flags(
        self,
        sxp: int,
        *,
        is_object: bool = False,
        has_attributes: bool = False,
        has_tag: bool = False,
        levels: int = 0,
    ) -> None:
        self.integer(
            sxp | is_object << 8 | has_attributes << 9 | has_tag << 10 | levels << 12
        )

    def string(self, value: str | None) -> None:
        if value is None:
            self.integer(_CHARSXP)
            self.integer(-1)
            return
        encoded = value.encode("utf-8")
        # R marks a string ASCII (64) or UTF-8 (8).
        self.flags(_CHARSXP, levels=64 if encoded.isascii() else 8)
        self.integer(len(encoded))
        self.out += encoded

    def symbol(self, name: str) -> None:
        # A symbol is written once; later uses refer to its place in the table.
        if name in self.symbols:
            self.integer(self.symbols[name] << 8 | _REFSXP)
            return
        self.symbols[name] = len(self.symbols) + 1
        self.integer(_SYMSXP)
        self.string(name)

    def item(self, value: _RValue) -> None:
        if value.kind == "null":
            self.integer(_NILVALUE_SXP)
            return
        is_object = any(name == "class" for name, _ in value.attributes)
        self.flags(
            _SXP_TYPES[value.kind],
            is_object=is_object,
            has_attributes=bool(value.attributes),
            levels=16 if value.s4 else 0,  # the S4 bit
        )
        if value.kind != "s4":
            self.integer(len(value.data))
            if value.kind in ("lgl", "int"):
                self.out += struct.pack(f">{len(value.data)}i", *value.data)
            elif value.kind == "dbl":
                self.out += struct.pack(f">{len(value.data)}d", *value.data)
            elif value.kind == "str":
                for text in value.data:
                    self.string(text)
            else:
                for element in value.data:
                    self.item(element)
        if value.attributes:
            for name, attribute in value.attributes:
                self.flags(_LISTSXP, has_tag=True)
                self.symbol(name)
                self.item(attribute)
            self.integer(_NILVALUE_SXP)


def _citeseq_rds(out_path: Path) -> None:
    """The Seurat object as ``saveRDS(seu, file)`` writes it: format 3, gzip level 6."""
    writer = _RdsWriter()
    writer.header()
    writer.item(_citeseq_seurat_object())
    # mtime 0 keeps the gzip header free of the time of writing.
    out_path.write_bytes(gzip.compress(bytes(writer.out), compresslevel=6, mtime=0))


#: Each dataset's writer, in catalog order.
DATASET_WRITERS = {
    "synthetic-cell-types-2d": _cell_types_2d,
    "synthetic-development-3d": _development_3d,
    "synthetic-trajectory-1d": _trajectory_1d,
}


def _catalog_document() -> dict[str, Any]:
    """``datasets.json``: which datasets the folder holds, and where.

    Names, descriptions and counts are not repeated here: Cellucid reads them
    from each file's ``uns['cellucid']`` (the Seurat object's
    ``misc$cellucid``) when it lists the catalog.
    """
    return {
        "schema": 2,
        "default": DEFAULT_DATASET,
        "datasets": [
            *(
                {"id": dataset_id, "path": f"{dataset_id}.zarr/"}
                for dataset_id in DATASET_WRITERS
            ),
            {"id": CITESEQ_MUDATA, "path": f"{CITESEQ_MUDATA}.zarr/"},
            {"id": CITESEQ_SEURAT, "path": f"{CITESEQ_SEURAT}.rds"},
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
        raise RuntimeError(
            f"{where} is not the expected dataset: {'; '.join(problems)}."
        )


def _require_citeseq(root: Path) -> None:
    """Read the CITE-seq MuData back with mudata, and check the .rds stream's header."""
    import mudata

    store = root / f"{CITESEQ_MUDATA}.zarr"
    mdata = mudata.read_zarr(store)
    meta = mdata.uns.get("cellucid", {})
    modalities = {name: modality.n_vars for name, modality in mdata.mod.items()}
    problems = []
    if meta.get("id") != CITESEQ_MUDATA:
        problems.append(f"uns['cellucid']['id'] is {meta.get('id')!r}")
    if mdata.n_obs != CITESEQ_CELLS or modalities != EXPECTED_MODALITIES:
        problems.append(f"{mdata.n_obs} cells, modalities {modalities}")
    if sorted(mdata.obsm) != ["X_pca", "X_umap", "prot", "rna"]:
        problems.append(f"obsm is {sorted(mdata.obsm)}")
    if problems:
        raise RuntimeError(
            f"{store} is not the expected MuData: {'; '.join(problems)}."
        )
    rds = root / f"{CITESEQ_SEURAT}.rds"
    if not gzip.decompress(rds.read_bytes()).startswith(b"X\n\x00\x00\x00\x03"):
        raise RuntimeError(f"{rds} is not R serialization format 3 in gzip.")


def _check_generation(root: Path) -> None:
    catalog_root = root / CATALOG
    expected_entries = {
        "datasets.json",
        *(f"{key}.zarr" for key in DATASET_WRITERS),
        f"{CITESEQ_MUDATA}.zarr",
        f"{CITESEQ_SEURAT}.rds",
    }
    entries = {path.name for path in catalog_root.iterdir()}
    if entries != expected_entries:
        raise RuntimeError(
            f"{CATALOG}/ holds {sorted(entries)}, not {sorted(expected_entries)}."
        )
    catalog = json.loads((catalog_root / "datasets.json").read_text(encoding="utf-8"))
    if catalog != _catalog_document():
        raise RuntimeError(
            f"{CATALOG}/datasets.json differs from the catalog it lists."
        )
    for dataset_id in DATASET_WRITERS:
        store = catalog_root / f"{dataset_id}.zarr"
        _require_dataset(_read_store(store), dataset_id, store)
    _require_citeseq(catalog_root)

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
    _citeseq_mudata(catalog_root / f"{CITESEQ_MUDATA}.zarr", "zarr")
    _citeseq_rds(catalog_root / f"{CITESEQ_SEURAT}.rds")
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
            raise RuntimeError(
                f"Generated entries must not be symbolic links: {relative}"
            )
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
                    backups[name] = (
                        REPOSITORY_ROOT / f".{name}-backup-{uuid.uuid4().hex}"
                    )
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
