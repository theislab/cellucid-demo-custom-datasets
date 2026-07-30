#!/usr/bin/env python3
# /// script
# requires-python = ">=3.11,<3.15"
# dependencies = [
#   "cellucid @ https://github.com/theislab/cellucid-python/archive/ce5fb6dab9afc6f2d4580da063773ffe3e641bc7.zip",
#   "numpy==2.5.1",
#   "pandas==2.3.3",
#   "scipy==1.18.0",
# ]
# ///
"""Generate the complete deterministic synthetic Cellucid catalog."""

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

import cellucid
import numpy as np
import pandas as pd
from cellucid import prepare
from cellucid.prepare_data import generate_datasets_manifest
from scipy import sparse

REPOSITORY_ROOT = Path(__file__).resolve().parent
DEFAULT_EXPORTS = REPOSITORY_ROOT / "exports"
CREATED_AT = "2026-07-27T00:00:00Z"
SOURCE_NAME = "Cellucid deterministic synthetic examples"
SOURCE_URL = "https://github.com/theislab/cellucid-demo-custom-datasets"
CELLUCID_COMMIT = "ce5fb6dab9afc6f2d4580da063773ffe3e641bc7"
CELLUCID_ARCHIVE_URL = (
    f"https://github.com/theislab/cellucid-python/archive/{CELLUCID_COMMIT}.zip"
)
EXPECTED_VERSIONS = {
    "cellucid": "0.9.1",
    "numpy": "2.5.1",
    "pandas": "2.3.3",
    "scipy": "1.18.0",
}
EXPECTED_DATASETS = {
    "synthetic-cell-types-2d": {
        "path": "synthetic-cell-types-2d/",
        "n_cells": 72,
        "n_genes": 8,
    },
    "synthetic-development-3d": {
        "path": "synthetic-development-3d/",
        "n_cells": 96,
        "n_genes": 10,
    },
    "synthetic-trajectory-1d": {
        "path": "synthetic-trajectory-1d/",
        "n_cells": 48,
        "n_genes": 6,
    },
}
DEFAULT_DATASET = "synthetic-cell-types-2d"


def _require_cellucid_provenance(
    distribution: importlib.metadata.Distribution,
) -> None:
    try:
        direct_url = json.loads(distribution.read_text("direct_url.json") or "")
    except (TypeError, json.JSONDecodeError) as error:
        raise RuntimeError(
            "Cellucid installation provenance is unavailable; "
            "run this file with `uv run generate_datasets.py ...`."
        ) from error

    url = direct_url.get("url") if isinstance(direct_url, dict) else None
    if url == CELLUCID_ARCHIVE_URL:
        _require_distribution_import_ownership(distribution)
        return

    vcs_info = direct_url.get("vcs_info") if isinstance(direct_url, dict) else None
    if (
        isinstance(url, str)
        and "github.com/theislab/cellucid-python" in url
        and isinstance(vcs_info, dict)
        and vcs_info.get("vcs") == "git"
        and vcs_info.get("commit_id") == CELLUCID_COMMIT
    ):
        _require_distribution_import_ownership(distribution)
        return

    directory_info = (
        direct_url.get("dir_info") if isinstance(direct_url, dict) else None
    )
    parsed_url = urllib.parse.urlparse(url) if isinstance(url, str) else None
    if (
        not isinstance(directory_info, dict)
        or directory_info.get("editable") is not True
        or parsed_url is None
        or parsed_url.scheme != "file"
    ):
        raise RuntimeError(
            f"Cellucid must come from pinned commit {CELLUCID_COMMIT}; "
            "run this file with `uv run generate_datasets.py ...`."
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
        head = subprocess.run(
            ["git", "-C", str(source_root), "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
        source_status = subprocess.run(
            [
                "git",
                "-C",
                str(source_root),
                "status",
                "--porcelain=v1",
                "--untracked-files=all",
                "--",
                "src/cellucid",
                "pyproject.toml",
            ],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as error:
        raise RuntimeError(
            "Cannot verify the editable Cellucid source checkout."
        ) from error
    if head != CELLUCID_COMMIT or source_status:
        raise RuntimeError(
            "Editable Cellucid must be an unchanged checkout of pinned commit "
            f"{CELLUCID_COMMIT}."
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
        "numpy": importlib.metadata.version("numpy"),
        "pandas": importlib.metadata.version("pandas"),
        "scipy": importlib.metadata.version("scipy"),
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
    matrix = np.zeros((n_cells, n_cells), dtype=np.float64)
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


def _common_prepare(
    *,
    out_dir: Path,
    dataset_id: str,
    dataset_name: str,
    description: str,
    embedding_key: str,
    embedding: np.ndarray,
    latent_space: np.ndarray,
    obs: pd.DataFrame,
    genes: list[str],
    expression: np.ndarray,
    connectivities: sparse.csr_matrix,
    vector_fields: dict[str, np.ndarray] | None = None,
) -> None:
    embedding_arguments = {
        "X_umap_1d": None,
        "X_umap_2d": None,
        "X_umap_3d": None,
    }
    embedding_arguments[embedding_key] = embedding
    var = pd.DataFrame(index=pd.Index(genes, name="gene_id"))
    prepare(
        latent_space=latent_space,
        obs=obs,
        var=var,
        gene_expression=expression,
        connectivities=connectivities,
        out_dir=out_dir,
        obs_keys=list(obs.columns),
        centroid_outlier_quantile=0.95,
        centroid_min_points=6,
        force=False,
        var_quantization=8,
        obs_continuous_quantization=8,
        compression=6,
        obs_categorical_dtype="uint8",
        dataset_name=dataset_name,
        dataset_id=dataset_id,
        created_at=CREATED_AT,
        dataset_description=description,
        source_name=SOURCE_NAME,
        source_url=SOURCE_URL,
        vector_fields=vector_fields,
        **embedding_arguments,
    )


def _build_trajectory_1d(out_dir: Path) -> None:
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
    genes = [
        "SYN_EARLY",
        "SYN_TRANSITION",
        "SYN_MATURE",
        "SYN_CYCLE",
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

    _common_prepare(
        out_dir=out_dir,
        dataset_id="synthetic-trajectory-1d",
        dataset_name="Synthetic trajectory — 1D",
        description=(
            "A synthetic 48-cell progression for learning planar 1D navigation, "
            "pseudotime coloring, gene trends, connectivity, and a velocity-style overlay."
        ),
        embedding_key="X_umap_1d",
        embedding=embedding,
        latent_space=np.column_stack([pseudotime, transcriptional_activity]),
        obs=obs,
        genes=genes,
        expression=expression,
        connectivities=_graph(n_cells, edges),
        vector_fields={"velocity_umap_1d": velocity},
    )


def _build_cell_types_2d(out_dir: Path) -> None:
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

    _common_prepare(
        out_dir=out_dir,
        dataset_id="synthetic-cell-types-2d",
        dataset_name="Synthetic cell-type islands — 2D",
        description=(
            "Three synthetic 2D populations for learning categorical and continuous "
            "coloring, marker genes, filtering, highlighting, and weighted connectivity."
        ),
        embedding_key="X_umap_2d",
        embedding=embedding,
        latent_space=embedding,
        obs=obs,
        genes=genes,
        expression=expression,
        connectivities=_graph(n_cells, edges),
    )


def _build_development_3d(out_dir: Path) -> None:
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
        "SYN_PULSE",
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

    _common_prepare(
        out_dir=out_dir,
        dataset_id="synthetic-development-3d",
        dataset_name="Synthetic branching development — 3D",
        description=(
            "A synthetic two-lineage 3D progression for orbit navigation, temporal "
            "metadata, fate markers, graph context, and a velocity-style vector field."
        ),
        embedding_key="X_umap_3d",
        embedding=embedding,
        latent_space=embedding,
        obs=obs,
        genes=genes,
        expression=expression,
        connectivities=_graph(n_cells, edges),
        vector_fields={"velocity_umap_3d": velocity},
    )


def _load_json_object(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise RuntimeError(f"Cannot read generated JSON object {path}.") from error
    if not isinstance(value, dict):
        raise TypeError(f"Generated JSON root must be an object: {path}")
    return value


def _remove_export_locks(exports_root: Path) -> None:
    """Remove the writer coordination files created beside each fresh export."""
    for dataset_id in EXPECTED_DATASETS:
        lock_path = exports_root / f".{dataset_id}.cellucid.lock"
        try:
            lock_path.unlink()
        except FileNotFoundError:
            pass


def _sanity_check_generation(exports_root: Path) -> None:
    expected_root_entries = {"datasets.json", *EXPECTED_DATASETS}
    actual_root_entries = {path.name for path in exports_root.iterdir()}
    if actual_root_entries != expected_root_entries:
        raise RuntimeError(
            "Generated exports root entries differ: "
            f"expected={sorted(expected_root_entries)}, "
            f"actual={sorted(actual_root_entries)}"
        )

    catalog_path = exports_root / "datasets.json"
    if catalog_path.is_symlink() or not catalog_path.is_file():
        raise RuntimeError(f"Generated catalog is not a regular file: {catalog_path}")
    catalog = _load_json_object(catalog_path)
    if catalog.get("version") != 1 or catalog.get("default") != DEFAULT_DATASET:
        raise RuntimeError("Generated catalog version or default dataset differs.")
    entries = catalog.get("datasets")
    if not isinstance(entries, list) or len(entries) != len(EXPECTED_DATASETS):
        raise RuntimeError("Generated catalog does not contain exactly three datasets.")

    entries_by_id: dict[str, dict[str, Any]] = {}
    for entry in entries:
        if not isinstance(entry, dict) or not isinstance(entry.get("id"), str):
            raise TypeError("Generated catalog contains an invalid dataset entry.")
        dataset_id = entry["id"]
        if dataset_id in entries_by_id:
            raise RuntimeError(f"Generated catalog repeats dataset id {dataset_id!r}.")
        entries_by_id[dataset_id] = entry
    if set(entries_by_id) != set(EXPECTED_DATASETS):
        raise RuntimeError(
            f"Generated catalog dataset ids differ: actual={sorted(entries_by_id)}"
        )

    for dataset_id, expected in EXPECTED_DATASETS.items():
        dataset_root = exports_root / dataset_id
        if dataset_root.is_symlink() or not dataset_root.is_dir():
            raise RuntimeError(
                f"Generated dataset path is not a directory: {dataset_root}"
            )
        entry = entries_by_id[dataset_id]
        if (
            entry.get("path") != expected["path"]
            or entry.get("n_cells") != expected["n_cells"]
            or entry.get("n_genes") != expected["n_genes"]
        ):
            raise RuntimeError(
                f"Generated catalog counts or path differ for {dataset_id!r}."
            )

        identity_path = dataset_root / "dataset_identity.json"
        if identity_path.is_symlink() or not identity_path.is_file():
            raise RuntimeError(
                f"Generated dataset identity is not a regular file: {identity_path}"
            )
        identity = _load_json_object(identity_path)
        stats = identity.get("stats")
        if (
            identity.get("id") != dataset_id
            or not isinstance(stats, dict)
            or stats.get("n_cells") != expected["n_cells"]
            or stats.get("n_genes") != expected["n_genes"]
        ):
            raise RuntimeError(
                f"Generated identity counts or id differ for {dataset_id!r}."
            )


def build_generation(exports_root: Path) -> None:
    exports_root.mkdir(parents=True, exist_ok=False)
    _build_cell_types_2d(exports_root / "synthetic-cell-types-2d")
    _build_development_3d(exports_root / "synthetic-development-3d")
    _build_trajectory_1d(exports_root / "synthetic-trajectory-1d")
    generate_datasets_manifest(
        exports_root,
        default_dataset=DEFAULT_DATASET,
    )
    _remove_export_locks(exports_root)
    _sanity_check_generation(exports_root)


def _tree_entries(root: Path) -> dict[str, tuple[str, bytes]]:
    entries: dict[str, tuple[str, bytes]] = {}
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        if path.is_symlink():
            raise RuntimeError(f"Export entries must not be symbolic links: {relative}")
        if path.is_dir():
            entries[relative] = ("directory", b"")
        elif path.is_file():
            entries[relative] = ("file", path.read_bytes())
        else:
            raise RuntimeError(
                f"Export entries must be regular files or directories: {relative}"
            )
    return entries


def _check_reproducibility(expected_root: Path) -> None:
    if expected_root.is_symlink() or not expected_root.is_dir():
        raise FileNotFoundError(
            f"Checked-in exports must be a real directory: {expected_root}"
        )
    _sanity_check_generation(expected_root)
    with tempfile.TemporaryDirectory(prefix="cellucid-rebuild-") as temporary:
        rebuilt = Path(temporary) / "exports"
        build_generation(rebuilt)
        expected = _tree_entries(expected_root)
        actual = _tree_entries(rebuilt)
        if expected != actual:
            missing = sorted(expected.keys() - actual.keys())
            extra = sorted(actual.keys() - expected.keys())
            changed = sorted(
                path
                for path in expected.keys() & actual.keys()
                if expected[path] != actual[path]
            )
            raise RuntimeError(
                "Checked-in exports are not byte-reproducible: "
                f"missing={missing}, extra={extra}, changed={changed}"
            )
    print("PASS: checked-in exports are byte-identical to a clean rebuild")


def _publish_generation() -> None:
    target = DEFAULT_EXPORTS
    if target.is_symlink():
        raise ValueError(f"Export root must not be a symbolic link: {target}")
    if target.exists() and not target.is_dir():
        raise NotADirectoryError(f"Export root is not a directory: {target}")

    target.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(
        tempfile.mkdtemp(
            prefix=f".{target.name}-stage-",
            dir=target.parent,
        )
    )
    stage.rmdir()
    backup: Path | None = None
    try:
        build_generation(stage)
        if target.exists():
            backup = target.parent / f".{target.name}-backup-{uuid.uuid4().hex}"
            os.replace(target, backup)
        try:
            os.replace(stage, target)
        except BaseException:
            if backup is not None:
                os.replace(backup, target)
                backup = None
            raise
        if backup is not None:
            shutil.rmtree(backup)
            backup = None
    finally:
        if stage.exists():
            shutil.rmtree(stage)
        if backup is not None and backup.exists():
            raise RuntimeError(
                f"Published exports are complete, but old backup cleanup failed: {backup}"
            )


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build or verify the deterministic synthetic Cellucid catalog."
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "--force",
        action="store_true",
        help="Replace this repository's exports/ from a clean staged generation.",
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
        _check_reproducibility(DEFAULT_EXPORTS)
    else:
        _publish_generation()
        print(f"PASS: published deterministic exports to {DEFAULT_EXPORTS}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, RuntimeError, TypeError, ValueError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
