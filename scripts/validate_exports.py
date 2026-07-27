#!/usr/bin/env python3
"""Validate the committed Cellucid catalog and every declared payload."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import re
import struct
import sys
import zlib
from pathlib import Path, PurePosixPath
from typing import Any

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPORTS = REPOSITORY_ROOT / "exports"
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")
SAFE_SEGMENT = re.compile(r"^[A-Za-z0-9_.-]+$")
EXPECTED = {
    "synthetic-cell-types-2d": {
        "name": "Synthetic cell-type islands — 2D",
        "n_cells": 72,
        "n_genes": 8,
        "dimensions": [2],
        "obs": {
            "cell_type": ("category", 3),
            "batch": ("category", 3),
            "quality_score": ("continuous", None),
            "library_size": ("continuous", None),
        },
        "vector_dimensions": None,
    },
    "synthetic-development-3d": {
        "name": "Synthetic branching development — 3D",
        "n_cells": 96,
        "n_genes": 10,
        "dimensions": [3],
        "obs": {
            "lineage": ("category", 2),
            "timepoint": ("category", 4),
            "condition": ("category", 2),
            "differentiation_score": ("continuous", None),
            "response_score": ("continuous", None),
        },
        "vector_dimensions": [3],
    },
    "synthetic-trajectory-1d": {
        "name": "Synthetic trajectory — 1D",
        "n_cells": 48,
        "n_genes": 6,
        "dimensions": [1],
        "obs": {
            "stage": ("category", 4),
            "replicate": ("category", 2),
            "pseudotime": ("continuous", None),
            "transcriptional_activity": ("continuous", None),
        },
        "vector_dimensions": [1],
    },
}


class ContractError(ValueError):
    """The export tree violates the one current repository contract."""


def _object_pairs(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"JSON object contains duplicate key {key!r}.")
        result[key] = value
    return result


def _load_json(path: Path) -> Any:
    try:
        return json.loads(
            path.read_text(encoding="utf-8"),
            object_pairs_hook=_object_pairs,
            parse_constant=lambda value: (_ for _ in ()).throw(
                ContractError(f"{path} contains non-finite JSON value {value}.")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise ContractError(f"{path} must contain valid UTF-8 JSON.") from error


def _exact_keys(
    value: Any,
    required: set[str],
    optional: set[str],
    label: str,
) -> dict[str, Any]:
    if type(value) is not dict:
        raise ContractError(f"{label} must be an object.")
    missing = sorted(required - value.keys())
    extra = sorted(value.keys() - required - optional)
    if missing or extra:
        raise ContractError(f"{label} keys differ: missing={missing}, extra={extra}.")
    return value


def _safe_relative(value: Any, *, directory: bool = False) -> str:
    if type(value) is not str or value == "" or value != value.strip():
        raise ContractError("Artifact path must be exact non-empty text.")
    if directory != value.endswith("/"):
        raise ContractError(f"Artifact path has incorrect directory suffix: {value!r}.")
    trimmed = value[:-1] if directory else value
    pure = PurePosixPath(trimmed)
    if (
        pure.is_absolute()
        or pure.as_posix() != trimmed
        or "\\" in value
        or "?" in value
        or "#" in value
        or any(part in {"", ".", ".."} for part in pure.parts)
        or any(not SAFE_SEGMENT.fullmatch(part) for part in pure.parts)
    ):
        raise ContractError(f"Artifact path is not canonical and safe: {value!r}.")
    return value


def _payload(path: Path) -> bytes:
    data = path.read_bytes()
    if path.suffix == ".gz":
        if len(data) < 10 or data[:2] != b"\x1f\x8b":
            raise ContractError(f"{path} is not gzip data.")
        if data[4:8] != b"\x00\x00\x00\x00":
            raise ContractError(f"{path} gzip timestamp is not deterministic zero.")
        try:
            return gzip.decompress(data)
        except (gzip.BadGzipFile, EOFError, zlib.error) as error:
            raise ContractError(f"{path} is corrupt gzip data.") from error
    return data


def _finite_float32(path: Path, count: int) -> None:
    data = _payload(path)
    if len(data) != count * 4:
        raise ContractError(f"{path} has {len(data)} bytes; expected {count * 4}.")
    values = struct.unpack(f"<{count}f", data)
    if not all(math.isfinite(value) for value in values):
        raise ContractError(f"{path} contains non-finite Float32 values.")


def _validate_connectivity(
    dataset_root: Path,
    n_cells: int,
    expected_files: set[str],
) -> int:
    manifest_path = dataset_root / "connectivity_manifest.json"
    manifest = _exact_keys(
        _load_json(manifest_path),
        {
            "format",
            "n_cells",
            "n_edges",
            "max_neighbors",
            "index_bytes",
            "index_dtype",
            "sourcesPath",
            "destinationsPath",
            "weightsPath",
            "weight_dtype",
            "weight_bytes",
            "compression",
        },
        set(),
        f"{manifest_path} connectivity manifest",
    )
    if (
        manifest["format"] != "edge_pairs"
        or manifest["n_cells"] != n_cells
        or manifest["compression"] != 6
        or manifest["index_dtype"] != "uint16"
        or manifest["index_bytes"] != 2
        or manifest["weight_dtype"] != "float64"
        or manifest["weight_bytes"] != 8
    ):
        raise ContractError(
            f"{manifest_path} has incorrect fixed connectivity metadata."
        )
    n_edges = manifest["n_edges"]
    if type(n_edges) is not int or n_edges <= 0:
        raise ContractError(f"{manifest_path} requires a positive edge count.")

    payloads = []
    for key in ("sourcesPath", "destinationsPath", "weightsPath"):
        relative = _safe_relative(manifest[key])
        expected_files.add(relative)
        payloads.append(_payload(dataset_root / relative))
    sources_data, destinations_data, weights_data = payloads
    if len(sources_data) != n_edges * 2 or len(destinations_data) != n_edges * 2:
        raise ContractError(f"{manifest_path} connectivity index sizes disagree.")
    if len(weights_data) != n_edges * 8:
        raise ContractError(f"{manifest_path} connectivity weight size disagrees.")
    sources = struct.unpack(f"<{n_edges}H", sources_data)
    destinations = struct.unpack(f"<{n_edges}H", destinations_data)
    weights = struct.unpack(f"<{n_edges}d", weights_data)
    pairs = list(zip(sources, destinations, strict=True))
    if pairs != sorted(set(pairs)):
        raise ContractError(
            f"{manifest_path} edge pairs are not unique lexicographic data."
        )
    if not all(0 <= source < destination < n_cells for source, destination in pairs):
        raise ContractError(
            f"{manifest_path} requires source < destination within bounds."
        )
    if not all(math.isfinite(weight) and weight > 0 for weight in weights):
        raise ContractError(f"{manifest_path} requires finite positive weights.")
    degrees = [0] * n_cells
    for source, destination in pairs:
        degrees[source] += 1
        degrees[destination] += 1
    if max(degrees) != manifest["max_neighbors"]:
        raise ContractError(f"{manifest_path} max_neighbors is not truthful.")
    expected_files.add("connectivity_manifest.json")
    return n_edges


def _validate_dataset(dataset_root: Path, expected: dict[str, Any]) -> set[str]:
    if not dataset_root.is_dir() or dataset_root.is_symlink():
        raise ContractError(f"Dataset root must be a real directory: {dataset_root}.")
    dataset_id = dataset_root.name
    expected_files = {
        "dataset_identity.json",
        "obs_manifest.json",
        "var_manifest.json",
    }
    identity_path = dataset_root / "dataset_identity.json"
    identity_required = {
        "version",
        "id",
        "name",
        "description",
        "created_at",
        "cellucid_data_version",
        "stats",
        "embeddings",
        "obs_fields",
        "export_settings",
        "source",
    }
    if expected["vector_dimensions"] is not None:
        identity_required.add("vector_fields")
    identity = _exact_keys(
        _load_json(identity_path),
        identity_required,
        set(),
        f"{identity_path} identity",
    )
    if (
        identity["version"] != 2
        or identity["id"] != dataset_id
        or identity["name"] != expected["name"]
        or identity["created_at"] != "2026-07-27T00:00:00Z"
        or identity["cellucid_data_version"] != "0.9.1"
        or type(identity["description"]) is not str
        or identity["description"] == ""
    ):
        raise ContractError(f"{identity_path} identity values differ.")
    settings = _exact_keys(
        identity["export_settings"],
        {
            "compression",
            "var_quantization",
            "obs_continuous_quantization",
            "obs_categorical_dtype",
        },
        set(),
        f"{identity_path} export settings",
    )
    if settings != {
        "compression": 6,
        "var_quantization": 8,
        "obs_continuous_quantization": 8,
        "obs_categorical_dtype": "uint8",
    }:
        raise ContractError(f"{identity_path} export settings differ.")
    source = _exact_keys(
        identity["source"],
        {"name", "url"},
        set(),
        f"{identity_path} source",
    )
    if source != {
        "name": "Cellucid deterministic synthetic examples",
        "url": "https://github.com/theislab/cellucid-demo-custom-datasets",
    }:
        raise ContractError(f"{identity_path} source differs.")

    stats = _exact_keys(
        identity["stats"],
        {
            "n_cells",
            "n_genes",
            "n_obs_fields",
            "n_categorical_fields",
            "n_continuous_fields",
            "has_connectivity",
            "n_edges",
        },
        set(),
        f"{identity_path} stats",
    )
    obs_expected = expected["obs"]
    categorical_count = sum(kind == "category" for kind, _ in obs_expected.values())
    continuous_count = len(obs_expected) - categorical_count
    if (
        stats["n_cells"] != expected["n_cells"]
        or stats["n_genes"] != expected["n_genes"]
        or stats["n_obs_fields"] != len(obs_expected)
        or stats["n_categorical_fields"] != categorical_count
        or stats["n_continuous_fields"] != continuous_count
        or stats["has_connectivity"] is not True
    ):
        raise ContractError(f"{identity_path} statistics differ.")

    embeddings = _exact_keys(
        identity["embeddings"],
        {"available_dimensions", "default_dimension", "files"},
        set(),
        f"{identity_path} embeddings",
    )
    dimensions = expected["dimensions"]
    if embeddings["available_dimensions"] != dimensions or embeddings[
        "default_dimension"
    ] != max(dimensions):
        raise ContractError(f"{identity_path} embedding dimensions differ.")
    expected_embedding_files = {
        f"{dimension}d": f"points_{dimension}d.bin.gz" for dimension in dimensions
    }
    if embeddings["files"] != expected_embedding_files:
        raise ContractError(f"{identity_path} embedding files differ.")
    for dimension, relative in (
        (dimension, expected_embedding_files[f"{dimension}d"])
        for dimension in dimensions
    ):
        _finite_float32(dataset_root / relative, expected["n_cells"] * dimension)
        expected_files.add(relative)

    identity_obs: dict[str, tuple[str, int | None]] = {}
    if type(identity["obs_fields"]) is not list or len(identity["obs_fields"]) != len(
        obs_expected
    ):
        raise ContractError(
            f"{identity_path} obs_fields must be the exact field array."
        )
    for field in identity["obs_fields"]:
        if type(field) is not dict:
            raise ContractError(f"{identity_path} obs field must be an object.")
        required = {"key", "kind"}
        if field.get("kind") == "category":
            required.add("n_categories")
        record = _exact_keys(field, required, set(), "identity obs field")
        key = record["key"]
        if (
            type(key) is not str
            or not SAFE_SEGMENT.fullmatch(key)
            or key in identity_obs
        ):
            raise ContractError(f"{identity_path} has an invalid obs field key.")
        identity_obs[key] = (
            record["kind"],
            record.get("n_categories"),
        )
    if identity_obs != obs_expected:
        raise ContractError(f"{identity_path} observation inventory differs.")

    obs_path = dataset_root / "obs_manifest.json"
    obs_manifest = _exact_keys(
        _load_json(obs_path),
        {
            "_format",
            "n_points",
            "centroid_outlier_quantile",
            "latent_key",
            "compression",
            "_obsSchemas",
            "_continuousFields",
            "_categoricalFields",
        },
        set(),
        f"{obs_path} observation manifest",
    )
    if (
        obs_manifest["_format"] != "compact_v1"
        or obs_manifest["n_points"] != expected["n_cells"]
        or obs_manifest["centroid_outlier_quantile"] != 0.95
        or obs_manifest["latent_key"] != "latent_space"
        or obs_manifest["compression"] != 6
    ):
        raise ContractError(f"{obs_path} fixed metadata differs.")
    schemas = _exact_keys(
        obs_manifest["_obsSchemas"],
        {"continuous", "categorical"},
        set(),
        f"{obs_path} schemas",
    )
    if schemas["continuous"] != {
        "pathPattern": "obs/{key}.values.u8.gz",
        "ext": "u8",
        "dtype": "uint8",
        "quantized": True,
        "quantizationBits": 8,
    }:
        raise ContractError(f"{obs_path} continuous schema differs.")
    categorical_schema = schemas["categorical"]
    if categorical_schema != {
        "codesPathPattern": "obs/{key}.codes.{ext}.gz",
        "outlierPathPattern": "obs/{key}.outliers.u8.gz",
        "outlierExt": "u8",
        "outlierDtype": "uint8",
        "outlierQuantized": True,
    }:
        raise ContractError(f"{obs_path} categorical schema differs.")
    continuous_tuples: dict[str, list[Any]] = {}
    for record in obs_manifest["_continuousFields"]:
        if (
            type(record) is not list
            or len(record) == 0
            or type(record[0]) is not str
            or record[0] in continuous_tuples
        ):
            raise ContractError(f"{obs_path} has an invalid continuous field tuple.")
        continuous_tuples[record[0]] = record
    categorical_tuples: dict[str, list[Any]] = {}
    for record in obs_manifest["_categoricalFields"]:
        if (
            type(record) is not list
            or len(record) == 0
            or type(record[0]) is not str
            or record[0] in categorical_tuples
            or record[0] in continuous_tuples
        ):
            raise ContractError(f"{obs_path} has an invalid categorical field tuple.")
        categorical_tuples[record[0]] = record
    for key, (kind, count) in obs_expected.items():
        if kind == "continuous":
            record = continuous_tuples.get(key)
            if (
                type(record) is not list
                or len(record) != 3
                or not all(type(value) in {int, float} for value in record[1:])
                or not record[1] < record[2]
            ):
                raise ContractError(f"{obs_path} continuous field {key!r} differs.")
            relative = f"obs/{key}.values.u8.gz"
            if len(_payload(dataset_root / relative)) != expected["n_cells"]:
                raise ContractError(f"{relative} has an incorrect cell axis.")
            expected_files.add(relative)
        else:
            record = categorical_tuples.get(key)
            if type(record) is not list or len(record) != 7:
                raise ContractError(f"{obs_path} categorical field {key!r} differs.")
            categories, dtype, missing_code, centroids, outlier_min, outlier_max = (
                record[1:]
            )
            if (
                type(categories) is not list
                or len(categories) != count
                or len(set(categories)) != count
                or not all(
                    type(category) is str and category != "" for category in categories
                )
                or dtype != "uint8"
                or missing_code != 255
                or type(centroids) is not dict
                or set(centroids) != {str(value) for value in dimensions}
                or type(outlier_min) not in {int, float}
                or type(outlier_max) not in {int, float}
                or not math.isfinite(outlier_min)
                or not math.isfinite(outlier_max)
                or not outlier_min < outlier_max
            ):
                raise ContractError(f"{obs_path} categorical metadata {key!r} differs.")
            codes_relative = f"obs/{key}.codes.u8.gz"
            outliers_relative = f"obs/{key}.outliers.u8.gz"
            codes = _payload(dataset_root / codes_relative)
            outliers = _payload(dataset_root / outliers_relative)
            if len(codes) != expected["n_cells"] or any(
                value >= count for value in codes
            ):
                raise ContractError(f"{codes_relative} has invalid category codes.")
            if len(outliers) != expected["n_cells"]:
                raise ContractError(f"{outliers_relative} has an incorrect cell axis.")
            expected_files.update({codes_relative, outliers_relative})
    if set(continuous_tuples) | set(categorical_tuples) != set(obs_expected):
        raise ContractError(f"{obs_path} includes unexpected fields.")

    var_path = dataset_root / "var_manifest.json"
    var_manifest = _exact_keys(
        _load_json(var_path),
        {
            "_format",
            "n_points",
            "var_gene_id_column",
            "compression",
            "quantization",
            "_varSchema",
            "fields",
        },
        set(),
        f"{var_path} var manifest",
    )
    if (
        var_manifest["_format"] != "compact_v1"
        or var_manifest["n_points"] != expected["n_cells"]
        or var_manifest["var_gene_id_column"] is not None
        or var_manifest["compression"] != 6
        or var_manifest["quantization"] != 8
        or var_manifest["_varSchema"]
        != {
            "kind": "continuous",
            "pathPattern": "var/{key}.values.u8.gz",
            "ext": "u8",
            "dtype": "uint8",
            "quantized": True,
            "quantizationBits": 8,
        }
    ):
        raise ContractError(f"{var_path} schema differs.")
    fields = var_manifest["fields"]
    if type(fields) is not list or len(fields) != expected["n_genes"]:
        raise ContractError(f"{var_path} gene count differs.")
    gene_ids: set[str] = set()
    for record in fields:
        if (
            type(record) is not list
            or len(record) != 3
            or type(record[0]) is not str
            or record[0] == ""
            or not SAFE_SEGMENT.fullmatch(record[0])
            or record[0] in gene_ids
            or not all(
                type(value) in {int, float} and math.isfinite(value)
                for value in record[1:]
            )
            or not record[1] < record[2]
        ):
            raise ContractError(f"{var_path} contains an invalid gene tuple.")
        gene_ids.add(record[0])
        relative = f"var/{record[0]}.values.u8.gz"
        if len(_payload(dataset_root / relative)) != expected["n_cells"]:
            raise ContractError(f"{relative} has an incorrect cell axis.")
        expected_files.add(relative)

    edge_count = _validate_connectivity(
        dataset_root,
        expected["n_cells"],
        expected_files,
    )
    if edge_count != stats["n_edges"]:
        raise ContractError(f"{identity_path} n_edges disagrees with its manifest.")

    vector_dimensions = expected["vector_dimensions"]
    if vector_dimensions is not None:
        vectors = _exact_keys(
            identity["vector_fields"],
            {"default_field", "fields"},
            set(),
            f"{identity_path} vectors",
        )
        if vectors["default_field"] != "velocity_umap":
            raise ContractError(f"{identity_path} vector default differs.")
        vector_fields = _exact_keys(
            vectors["fields"],
            {"velocity_umap"},
            set(),
            f"{identity_path} vector fields",
        )
        field = _exact_keys(
            vector_fields["velocity_umap"],
            {"label", "basis", "available_dimensions", "default_dimension", "files"},
            set(),
            f"{identity_path} velocity field",
        )
        expected_vector_files = {
            f"{dimension}d": f"vectors/velocity_umap_{dimension}d.bin.gz"
            for dimension in vector_dimensions
        }
        if field != {
            "label": "velocity_umap",
            "basis": "umap",
            "available_dimensions": vector_dimensions,
            "default_dimension": max(vector_dimensions),
            "files": expected_vector_files,
        }:
            raise ContractError(f"{identity_path} velocity metadata differs.")
        for dimension in vector_dimensions:
            relative = expected_vector_files[f"{dimension}d"]
            _finite_float32(dataset_root / relative, expected["n_cells"] * dimension)
            expected_files.add(relative)

    actual_files = {
        path.relative_to(dataset_root).as_posix()
        for path in dataset_root.rglob("*")
        if path.is_file()
    }
    if actual_files != expected_files:
        raise ContractError(
            f"{dataset_root} files differ: "
            f"missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}."
        )
    return expected_files


def _validate_checksums(exports_root: Path, expected_paths: set[str]) -> None:
    checksum_path = exports_root / "SHA256SUMS"
    checksum_text = checksum_path.read_text(encoding="utf-8")
    if not checksum_text.endswith("\n"):
        raise ContractError("SHA256SUMS must end with one LF record terminator.")
    lines = checksum_text.splitlines()
    records: dict[str, str] = {}
    for line in lines:
        if len(line) < 67 or line[64:66] != "  ":
            raise ContractError("SHA256SUMS contains a malformed record.")
        digest, relative = line[:64], line[66:]
        if not SHA256_PATTERN.fullmatch(digest):
            raise ContractError("SHA256SUMS contains an invalid digest.")
        _safe_relative(relative)
        if relative in records:
            raise ContractError(f"SHA256SUMS duplicates {relative}.")
        records[relative] = digest
    if list(records) != sorted(records):
        raise ContractError("SHA256SUMS paths must be sorted.")
    if set(records) != expected_paths:
        raise ContractError(
            "SHA256SUMS inventory differs: "
            f"missing={sorted(expected_paths - records.keys())}, "
            f"extra={sorted(records.keys() - expected_paths)}."
        )
    for relative, digest in records.items():
        actual = hashlib.sha256((exports_root / relative).read_bytes()).hexdigest()
        if actual != digest:
            raise ContractError(f"SHA256 mismatch for {relative}.")


def validate_exports(exports_root: Path = DEFAULT_EXPORTS) -> dict[str, int]:
    candidate = Path(exports_root)
    if candidate.is_symlink():
        raise ContractError(f"Exports root must not be a symbolic link: {candidate}.")
    exports_root = candidate.resolve()
    if not exports_root.is_dir():
        raise ContractError(f"Exports root must be a real directory: {exports_root}.")
    for path in exports_root.rglob("*"):
        if path.is_symlink() or not (path.is_file() or path.is_dir()):
            raise ContractError(f"Export tree entry must be a real file: {path}.")
    catalog_path = exports_root / "datasets.json"
    catalog = _exact_keys(
        _load_json(catalog_path),
        {"version", "default", "datasets"},
        set(),
        f"{catalog_path} catalog",
    )
    if catalog["version"] != 1 or catalog["default"] != "synthetic-cell-types-2d":
        raise ContractError(f"{catalog_path} version/default differs.")
    entries = catalog["datasets"]
    if type(entries) is not list or len(entries) != len(EXPECTED):
        raise ContractError(f"{catalog_path} dataset count differs.")
    expected_order = sorted(EXPECTED)
    ordered_ids = [
        entry.get("id") if type(entry) is dict else None for entry in entries
    ]
    if ordered_ids != expected_order:
        raise ContractError(f"{catalog_path} dataset ordering differs.")

    expected_paths = {"datasets.json"}
    total_cells = 0
    total_genes = 0
    for entry in entries:
        record = _exact_keys(
            entry,
            {"id", "path", "name", "description", "n_cells", "n_genes"},
            set(),
            "catalog entry",
        )
        dataset_id = record["id"]
        if type(dataset_id) is not str or dataset_id not in EXPECTED:
            raise ContractError(
                f"Catalog dataset ID is not recognized: {dataset_id!r}."
            )
        expected = EXPECTED[dataset_id]
        expected_path = f"{dataset_id}/"
        if (
            _safe_relative(record["path"], directory=True) != expected_path
            or record["name"] != expected["name"]
            or type(record["description"]) is not str
            or record["description"] == ""
            or record["n_cells"] != expected["n_cells"]
            or record["n_genes"] != expected["n_genes"]
        ):
            raise ContractError(f"Catalog entry {dataset_id!r} differs.")
        dataset_root = exports_root / dataset_id
        files = _validate_dataset(dataset_root, expected)
        expected_paths.update(f"{dataset_id}/{relative}" for relative in files)
        total_cells += expected["n_cells"]
        total_genes += expected["n_genes"]

    actual_root_entries = {path.name for path in exports_root.iterdir()}
    if actual_root_entries != set(EXPECTED) | {"datasets.json", "SHA256SUMS"}:
        raise ContractError(
            f"Exports root entries differ: {sorted(actual_root_entries)}."
        )
    _validate_checksums(exports_root, expected_paths)
    return {
        "datasets": len(EXPECTED),
        "cells": total_cells,
        "genes": total_genes,
        "files": len(expected_paths) + 1,
    }


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate every catalog, manifest, and binary payload."
    )
    parser.add_argument(
        "exports_root",
        nargs="?",
        type=Path,
        default=DEFAULT_EXPORTS,
    )
    return parser.parse_args()


def main() -> int:
    result = validate_exports(_parse_arguments().exports_root)
    print(
        "PASS: "
        f"{result['datasets']} datasets, {result['cells']} cells, "
        f"{result['genes']} genes, {result['files']} files"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ContractError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
