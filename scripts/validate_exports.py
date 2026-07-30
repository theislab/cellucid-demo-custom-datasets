#!/usr/bin/env python3
"""Validate the complete, immutable Cellucid example export contract."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import re
import stat
import struct
import sys
import zlib
from pathlib import Path, PurePosixPath
from typing import Any, Iterable

REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EXPORTS = REPOSITORY_ROOT / "exports"
JSON_SIZE_LIMIT = 1_048_576
CHECKSUM_RECORDS = 71
SAFE_SEGMENT = re.compile(r"^[A-Za-z0-9_.-]+$")
SHA256_RECORD = re.compile(r"^([0-9a-f]{64})  ([A-Za-z0-9_./-]+)$")

EXPECTED: dict[str, dict[str, Any]] = {
    "synthetic-cell-types-2d": {
        "name": "Synthetic cell-type islands — 2D",
        "description": (
            "Three synthetic 2D populations for learning categorical and "
            "continuous coloring, marker genes, filtering, highlighting, and "
            "weighted connectivity."
        ),
        "n_cells": 72,
        "n_edges": 144,
        "dimension": 2,
        "genes": (
            "SYN_MARKER_A",
            "SYN_MARKER_B",
            "SYN_MARKER_C",
            "SYN_SHARED",
            "SYN_QUALITY",
            "SYN_LIBRARY",
            "SYN_AXIS_X",
            "SYN_AXIS_Y",
        ),
        "obs": (
            ("quality_score", "continuous", None),
            ("library_size", "continuous", None),
            ("cell_type", "category", ("type-a", "type-b", "type-c")),
            ("batch", "category", ("batch-1", "batch-2", "batch-3")),
        ),
        "vector": False,
    },
    "synthetic-development-3d": {
        "name": "Synthetic branching development — 3D",
        "description": (
            "A synthetic two-lineage 3D progression for orbit navigation, "
            "temporal metadata, fate markers, graph context, and a "
            "velocity-style vector field."
        ),
        "n_cells": 96,
        "n_edges": 187,
        "dimension": 3,
        "genes": (
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
        ),
        "obs": (
            ("differentiation_score", "continuous", None),
            ("response_score", "continuous", None),
            ("lineage", "category", ("lineage-left", "lineage-right")),
            ("timepoint", "category", ("t0", "t1", "t2", "t3")),
            ("condition", "category", ("control", "stimulated")),
        ),
        "vector": True,
    },
    "synthetic-trajectory-1d": {
        "name": "Synthetic trajectory — 1D",
        "description": (
            "A synthetic 48-cell progression for learning planar 1D "
            "navigation, pseudotime coloring, gene trends, connectivity, and "
            "a velocity-style overlay."
        ),
        "n_cells": 48,
        "n_edges": 93,
        "dimension": 1,
        "genes": (
            "SYN_EARLY",
            "SYN_TRANSITION",
            "SYN_MATURE",
            "SYN_CYCLE",
            "SYN_ACTIVITY",
            "SYN_LATE",
        ),
        "obs": (
            ("pseudotime", "continuous", None),
            ("transcriptional_activity", "continuous", None),
            ("stage", "category", ("early", "transition", "committed", "mature")),
            ("replicate", "category", ("replicate-a", "replicate-b")),
        ),
        "vector": True,
    },
}


class ContractError(ValueError):
    """The export tree violates this repository's exact current contract."""


def _duplicate_safe_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise ContractError(f"JSON object contains duplicate key {key!r}.")
        result[key] = value
    return result


def _finite_json_float(path: Path, token: str) -> float:
    value = float(token)
    if not math.isfinite(value):
        raise ContractError(f"{path} contains non-finite JSON number {token!r}.")
    return value


def _reject_nonfinite_constant(path: Path, token: str) -> None:
    raise ContractError(f"{path} contains non-finite JSON constant {token!r}.")


def _reject_unpaired_surrogates(value: Any, path: Path) -> None:
    if isinstance(value, str):
        if any(0xD800 <= ord(character) <= 0xDFFF for character in value):
            raise ContractError(f"{path} contains an unpaired Unicode surrogate.")
    elif type(value) is list:
        for item in value:
            _reject_unpaired_surrogates(item, path)
    elif type(value) is dict:
        for key, item in value.items():
            _reject_unpaired_surrogates(key, path)
            _reject_unpaired_surrogates(item, path)


def _load_json(path: Path) -> Any:
    try:
        size = path.stat().st_size
        if size <= 0 or size > JSON_SIZE_LIMIT:
            raise ContractError(
                f"{path} JSON size {size} is outside 1..{JSON_SIZE_LIMIT} bytes."
            )
        text = path.read_text(encoding="utf-8")
        value = json.loads(
            text,
            object_pairs_hook=_duplicate_safe_object,
            parse_float=lambda token: _finite_json_float(path, token),
            parse_constant=lambda token: _reject_nonfinite_constant(path, token),
        )
    except ContractError:
        raise
    except (
        OSError,
        UnicodeError,
        ValueError,
        json.JSONDecodeError,
        RecursionError,
    ) as error:
        raise ContractError(f"{path} must contain valid finite UTF-8 JSON.") from error
    _reject_unpaired_surrogates(value, path)
    return value


def _exact_keys(value: Any, keys: set[str], label: str) -> dict[str, Any]:
    if type(value) is not dict:
        raise ContractError(f"{label} must be an object.")
    missing = sorted(keys - value.keys())
    extra = sorted(value.keys() - keys)
    if missing or extra:
        raise ContractError(f"{label} keys differ: missing={missing}, extra={extra}.")
    return value


def _exact_int(value: Any, expected: int, label: str) -> None:
    if type(value) is not int or value != expected:
        raise ContractError(f"{label} must be the integer {expected}.")


def _finite_number(value: Any, label: str) -> float:
    if type(value) not in {int, float}:
        raise ContractError(f"{label} must be a finite number.")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ContractError(f"{label} must be a finite number.") from error
    if not math.isfinite(result):
        raise ContractError(f"{label} must be a finite number.")
    return result


def _safe_relative(value: Any, *, directory: bool = False) -> str:
    if type(value) is not str or not value or value != value.strip():
        raise ContractError("Artifact path must be exact non-empty text.")
    if value.endswith("/") is not directory:
        raise ContractError(f"Artifact path has an incorrect suffix: {value!r}.")
    trimmed = value[:-1] if directory else value
    pure = PurePosixPath(trimmed)
    if (
        len(value) > 240
        or pure.is_absolute()
        or pure.as_posix() != trimmed
        or "\\" in value
        or "?" in value
        or "#" in value
        or any(
            part in {"", ".", ".."}
            or len(part) > 100
            or not SAFE_SEGMENT.fullmatch(part)
            for part in pure.parts
        )
    ):
        raise ContractError(f"Artifact path is not canonical and safe: {value!r}.")
    return value


def _assert_casefold_unique(paths: Iterable[str], label: str) -> None:
    seen: dict[str, str] = {}
    for relative in paths:
        folded = relative.casefold()
        previous = seen.get(folded)
        if previous is not None and previous != relative:
            raise ContractError(
                f"{label} has a case-fold collision: {previous!r} and {relative!r}."
            )
        seen[folded] = relative


def _scan_real_tree(root: Path) -> dict[str, str]:
    if root.is_symlink():
        raise ContractError(f"Tree root must not be a symbolic link: {root}.")
    if not root.is_dir():
        raise ContractError(f"Tree root must be a real directory: {root}.")

    entries: dict[str, str] = {}
    pending: list[tuple[Path, str]] = [(root, "")]
    while pending:
        directory, prefix = pending.pop()
        try:
            children = sorted(os.scandir(directory), key=lambda item: item.name)
        except OSError as error:
            raise ContractError(f"Cannot inspect export directory {directory}.") from error
        for entry in children:
            relative = f"{prefix}{entry.name}"
            try:
                mode = entry.stat(follow_symlinks=False).st_mode
            except OSError as error:
                raise ContractError(f"Cannot inspect export entry {relative!r}.") from error
            if entry.is_symlink():
                raise ContractError(f"Export entry must not be a symlink: {relative}.")
            if stat.S_ISDIR(mode):
                entries[relative] = "directory"
                pending.append((Path(entry.path), f"{relative}/"))
            elif stat.S_ISREG(mode):
                entries[relative] = "file"
            else:
                raise ContractError(
                    f"Export entry must be a regular file or directory: {relative}."
                )
    _assert_casefold_unique(entries, "Export tree")
    return entries


def _artifact_path(dataset_root: Path, relative: str) -> Path:
    pure = PurePosixPath(_safe_relative(relative))
    path = dataset_root.joinpath(*pure.parts)
    if path.is_symlink() or not path.is_file():
        raise ContractError(f"Declared artifact must be a real file: {relative}.")
    return path


def _claim_payload(
    dataset_root: Path,
    relative: str,
    expected_size: int,
    expected_files: set[str],
) -> bytes:
    if relative in expected_files:
        raise ContractError(f"Artifact is declared more than once: {relative}.")
    _assert_casefold_unique((*expected_files, relative), "Artifact declarations")
    expected_files.add(relative)
    path = _artifact_path(dataset_root, relative)
    if path.suffix != ".gz":
        raise ContractError(f"Binary artifact must use .gz: {relative}.")
    if type(expected_size) is not int or expected_size <= 0:
        raise ContractError(f"Internal expected size is invalid for {relative}.")

    try:
        compressed_size = path.stat().st_size
        if compressed_size < 18 or compressed_size > expected_size + 1024:
            raise ContractError(
                f"{relative} compressed size {compressed_size} is outside the "
                f"bounded range 18..{expected_size + 1024}."
            )
        compressed = path.read_bytes()
    except ContractError:
        raise
    except OSError as error:
        raise ContractError(f"Cannot read declared artifact {relative}.") from error

    if (
        compressed[:3] != b"\x1f\x8b\x08"
        or compressed[3] != 0
        or compressed[4:8] != b"\x00\x00\x00\x00"
    ):
        raise ContractError(
            f"{relative} must have a deterministic, flag-free gzip header."
        )
    if int.from_bytes(compressed[-4:], "little") != expected_size % (1 << 32):
        raise ContractError(f"{relative} gzip ISIZE disagrees with its manifest.")

    decoder = zlib.decompressobj(16 + zlib.MAX_WBITS)
    try:
        payload = decoder.decompress(compressed, expected_size + 1)
    except zlib.error as error:
        raise ContractError(f"{relative} is corrupt gzip data.") from error
    if len(payload) > expected_size:
        raise ContractError(f"{relative} expands beyond its declared size.")
    if (
        not decoder.eof
        or decoder.unconsumed_tail
        or decoder.unused_data
        or len(payload) != expected_size
    ):
        raise ContractError(
            f"{relative} must be one exact gzip member expanding to "
            f"{expected_size} bytes, got {len(payload)}."
        )
    return payload


def _finite_float32(
    payload: bytes,
    count: int,
    label: str,
    *,
    dimension: int,
    require_each_axis_variation: bool,
) -> tuple[float, ...]:
    if len(payload) != count * 4:
        raise ContractError(f"{label} has an incorrect Float32 byte length.")
    try:
        values = struct.unpack(f"<{count}f", payload)
    except struct.error as error:
        raise ContractError(f"{label} is not aligned Float32 data.") from error
    if not all(math.isfinite(value) for value in values):
        raise ContractError(f"{label} contains non-finite Float32 data.")
    if require_each_axis_variation:
        for axis in range(dimension):
            axis_values = values[axis::dimension]
            if min(axis_values) >= max(axis_values):
                raise ContractError(f"{label} axis {axis} is degenerate.")
    return values


def _validate_centroids(
    centroids: Any,
    categories: tuple[str, ...],
    codes: bytes,
    points: tuple[float, ...],
    dimension: int,
    label: str,
) -> None:
    mapping = _exact_keys(centroids, {str(dimension)}, f"{label} centroids")
    records = mapping[str(dimension)]
    if type(records) is not list or len(records) != len(categories):
        raise ContractError(f"{label} requires one centroid per category.")
    coordinate_ranges = [
        (min(points[axis::dimension]), max(points[axis::dimension]))
        for axis in range(dimension)
    ]
    for index, (record, category) in enumerate(zip(records, categories, strict=True)):
        item = _exact_keys(
            record,
            {"category", "position", "n_points"},
            f"{label} centroid {index}",
        )
        if item["category"] != category:
            raise ContractError(f"{label} centroid category ordering differs.")
        position = item["position"]
        if type(position) is not list or len(position) != dimension:
            raise ContractError(f"{label} centroid position dimension differs.")
        for axis, coordinate in enumerate(position):
            value = _finite_number(coordinate, f"{label} centroid coordinate")
            minimum, maximum = coordinate_ranges[axis]
            if not minimum <= value <= maximum:
                raise ContractError(f"{label} centroid is outside the embedding.")
        n_points = item["n_points"]
        category_total = codes.count(index)
        if (
            type(n_points) is not int
            or n_points <= 0
            or n_points > category_total
        ):
            raise ContractError(f"{label} centroid support is not truthful.")


def _validate_obs(
    dataset_root: Path,
    expected: dict[str, Any],
    points: tuple[float, ...],
    expected_files: set[str],
) -> None:
    path = dataset_root / "obs_manifest.json"
    manifest = _exact_keys(
        _load_json(path),
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
        f"{path} observation manifest",
    )
    _exact_int(manifest["n_points"], expected["n_cells"], f"{path} n_points")
    _exact_int(manifest["compression"], 6, f"{path} compression")
    if (
        manifest["_format"] != "compact_v1"
        or manifest["latent_key"] != "latent_space"
        or manifest["centroid_outlier_quantile"] != 0.95
    ):
        raise ContractError(f"{path} fixed metadata differs.")
    schemas = _exact_keys(
        manifest["_obsSchemas"],
        {"continuous", "categorical"},
        f"{path} schemas",
    )
    if schemas["continuous"] != {
        "pathPattern": "obs/{key}.values.u8.gz",
        "ext": "u8",
        "dtype": "uint8",
        "quantized": True,
        "quantizationBits": 8,
    }:
        raise ContractError(f"{path} continuous schema differs.")
    if schemas["categorical"] != {
        "codesPathPattern": "obs/{key}.codes.{ext}.gz",
        "outlierPathPattern": "obs/{key}.outliers.u8.gz",
        "outlierExt": "u8",
        "outlierDtype": "uint8",
        "outlierQuantized": True,
    }:
        raise ContractError(f"{path} categorical schema differs.")

    obs = expected["obs"]
    continuous_expected = [item for item in obs if item[1] == "continuous"]
    categorical_expected = [item for item in obs if item[1] == "category"]
    continuous = manifest["_continuousFields"]
    categorical = manifest["_categoricalFields"]
    if (
        type(continuous) is not list
        or len(continuous) != len(continuous_expected)
        or type(categorical) is not list
        or len(categorical) != len(categorical_expected)
    ):
        raise ContractError(f"{path} field inventory differs.")

    for record, (key, _kind, _categories) in zip(
        continuous, continuous_expected, strict=True
    ):
        if type(record) is not list or len(record) != 3 or record[0] != key:
            raise ContractError(f"{path} continuous field ordering differs.")
        minimum = _finite_number(record[1], f"{path} {key} minimum")
        maximum = _finite_number(record[2], f"{path} {key} maximum")
        if minimum >= maximum:
            raise ContractError(f"{path} {key} range is degenerate.")
        relative = f"obs/{key}.values.u8.gz"
        payload = _claim_payload(
            dataset_root, relative, expected["n_cells"], expected_files
        )
        if 255 in payload:
            raise ContractError(f"{relative} contains an undeclared missing value.")

    dimension = expected["dimension"]
    for record, (key, _kind, category_values) in zip(
        categorical, categorical_expected, strict=True
    ):
        if type(record) is not list or len(record) != 7 or record[0] != key:
            raise ContractError(f"{path} categorical field ordering differs.")
        categories, dtype, missing, centroids, outlier_min, outlier_max = record[1:]
        expected_categories = tuple(category_values)
        if (
            type(categories) is not list
            or tuple(categories) != expected_categories
            or dtype != "uint8"
        ):
            raise ContractError(f"{path} categorical metadata for {key!r} differs.")
        _exact_int(missing, 255, f"{path} {key} missing code")
        minimum = _finite_number(outlier_min, f"{path} {key} outlier minimum")
        maximum = _finite_number(outlier_max, f"{path} {key} outlier maximum")
        if minimum >= maximum:
            raise ContractError(f"{path} {key} outlier range is degenerate.")
        codes_relative = f"obs/{key}.codes.u8.gz"
        outliers_relative = f"obs/{key}.outliers.u8.gz"
        codes = _claim_payload(
            dataset_root, codes_relative, expected["n_cells"], expected_files
        )
        outliers = _claim_payload(
            dataset_root, outliers_relative, expected["n_cells"], expected_files
        )
        if any(code >= len(expected_categories) for code in codes):
            raise ContractError(f"{codes_relative} contains an invalid category code.")
        if 255 in outliers:
            raise ContractError(
                f"{outliers_relative} contains an undeclared missing value."
            )
        _validate_centroids(
            centroids,
            expected_categories,
            codes,
            points,
            dimension,
            f"{path} {key}",
        )


def _validate_var(
    dataset_root: Path,
    expected: dict[str, Any],
    expected_files: set[str],
) -> None:
    path = dataset_root / "var_manifest.json"
    manifest = _exact_keys(
        _load_json(path),
        {
            "_format",
            "n_points",
            "var_gene_id_column",
            "compression",
            "quantization",
            "_varSchema",
            "fields",
        },
        f"{path} variable manifest",
    )
    _exact_int(manifest["n_points"], expected["n_cells"], f"{path} n_points")
    _exact_int(manifest["compression"], 6, f"{path} compression")
    _exact_int(manifest["quantization"], 8, f"{path} quantization")
    if (
        manifest["_format"] != "compact_v1"
        or manifest["var_gene_id_column"] is not None
        or manifest["_varSchema"]
        != {
            "kind": "continuous",
            "pathPattern": "var/{key}.values.u8.gz",
            "ext": "u8",
            "dtype": "uint8",
            "quantized": True,
            "quantizationBits": 8,
        }
    ):
        raise ContractError(f"{path} schema differs.")
    fields = manifest["fields"]
    genes = expected["genes"]
    if type(fields) is not list or len(fields) != len(genes):
        raise ContractError(f"{path} gene count differs.")
    for record, gene in zip(fields, genes, strict=True):
        if type(record) is not list or len(record) != 3 or record[0] != gene:
            raise ContractError(f"{path} gene ordering differs.")
        minimum = _finite_number(record[1], f"{path} {gene} minimum")
        maximum = _finite_number(record[2], f"{path} {gene} maximum")
        if minimum >= maximum:
            raise ContractError(f"{path} {gene} range is degenerate.")
        relative = f"var/{gene}.values.u8.gz"
        payload = _claim_payload(
            dataset_root, relative, expected["n_cells"], expected_files
        )
        if 255 in payload:
            raise ContractError(f"{relative} contains an undeclared missing value.")


def _validate_connectivity(
    dataset_root: Path,
    expected: dict[str, Any],
    expected_files: set[str],
) -> None:
    path = dataset_root / "connectivity_manifest.json"
    manifest = _exact_keys(
        _load_json(path),
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
        f"{path} connectivity manifest",
    )
    _exact_int(manifest["n_cells"], expected["n_cells"], f"{path} n_cells")
    _exact_int(manifest["n_edges"], expected["n_edges"], f"{path} n_edges")
    _exact_int(manifest["max_neighbors"], 4, f"{path} max_neighbors")
    _exact_int(manifest["index_bytes"], 2, f"{path} index_bytes")
    _exact_int(manifest["weight_bytes"], 8, f"{path} weight_bytes")
    _exact_int(manifest["compression"], 6, f"{path} compression")
    fixed = {
        "format": "edge_pairs",
        "index_dtype": "uint16",
        "sourcesPath": "connectivity/edges.src.bin.gz",
        "destinationsPath": "connectivity/edges.dst.bin.gz",
        "weightsPath": "connectivity/edges.weights.f64.bin.gz",
        "weight_dtype": "float64",
    }
    if any(manifest[key] != value for key, value in fixed.items()):
        raise ContractError(f"{path} fixed connectivity metadata differs.")

    n_edges = expected["n_edges"]
    sources_data = _claim_payload(
        dataset_root, manifest["sourcesPath"], n_edges * 2, expected_files
    )
    destinations_data = _claim_payload(
        dataset_root, manifest["destinationsPath"], n_edges * 2, expected_files
    )
    weights_data = _claim_payload(
        dataset_root, manifest["weightsPath"], n_edges * 8, expected_files
    )
    sources = struct.unpack(f"<{n_edges}H", sources_data)
    destinations = struct.unpack(f"<{n_edges}H", destinations_data)
    weights = struct.unpack(f"<{n_edges}d", weights_data)
    previous = (-1, -1)
    degrees = [0] * expected["n_cells"]
    for source, destination, weight in zip(
        sources, destinations, weights, strict=True
    ):
        pair = (source, destination)
        if not 0 <= source < destination < expected["n_cells"]:
            raise ContractError(f"{path} has an out-of-bounds or self edge.")
        if pair <= previous:
            raise ContractError(f"{path} edges are not unique lexicographic data.")
        if not math.isfinite(weight) or weight <= 0:
            raise ContractError(f"{path} requires finite positive edge weights.")
        degrees[source] += 1
        degrees[destination] += 1
        previous = pair
    if max(degrees) != manifest["max_neighbors"]:
        raise ContractError(f"{path} max_neighbors is not truthful.")


def _validate_identity(
    dataset_root: Path,
    expected: dict[str, Any],
    expected_files: set[str],
) -> tuple[float, ...]:
    dataset_id = dataset_root.name
    path = dataset_root / "dataset_identity.json"
    identity_keys = {
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
    if expected["vector"]:
        identity_keys.add("vector_fields")
    identity = _exact_keys(_load_json(path), identity_keys, f"{path} identity")
    _exact_int(identity["version"], 2, f"{path} version")
    if (
        identity["id"] != dataset_id
        or identity["name"] != expected["name"]
        or identity["description"] != expected["description"]
        or identity["created_at"] != "2026-07-27T00:00:00Z"
        or identity["cellucid_data_version"] != "0.9.1"
    ):
        raise ContractError(f"{path} fixed identity metadata differs.")

    continuous_count = sum(item[1] == "continuous" for item in expected["obs"])
    categorical_count = len(expected["obs"]) - continuous_count
    stats_value = _exact_keys(
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
        f"{path} statistics",
    )
    expected_stats = {
        "n_cells": expected["n_cells"],
        "n_genes": len(expected["genes"]),
        "n_obs_fields": len(expected["obs"]),
        "n_categorical_fields": categorical_count,
        "n_continuous_fields": continuous_count,
        "n_edges": expected["n_edges"],
    }
    for key, value in expected_stats.items():
        _exact_int(stats_value[key], value, f"{path} statistics {key}")
    if stats_value["has_connectivity"] is not True:
        raise ContractError(f"{path} must declare connectivity.")

    if identity["export_settings"] != {
        "compression": 6,
        "var_quantization": 8,
        "obs_continuous_quantization": 8,
        "obs_categorical_dtype": "uint8",
    }:
        raise ContractError(f"{path} export settings differ.")
    if identity["source"] != {
        "name": "Cellucid deterministic synthetic examples",
        "url": "https://github.com/theislab/cellucid-demo-custom-datasets",
    }:
        raise ContractError(f"{path} source differs.")

    expected_obs_fields = []
    for key, kind, categories in expected["obs"]:
        record = {"key": key, "kind": kind}
        if categories is not None:
            record["n_categories"] = len(categories)
        expected_obs_fields.append(record)
    if identity["obs_fields"] != expected_obs_fields:
        raise ContractError(f"{path} observation identity inventory differs.")

    dimension = expected["dimension"]
    point_relative = f"points_{dimension}d.bin.gz"
    embeddings = _exact_keys(
        identity["embeddings"],
        {"available_dimensions", "default_dimension", "files"},
        f"{path} embeddings",
    )
    if (
        type(embeddings["available_dimensions"]) is not list
        or embeddings["available_dimensions"] != [dimension]
        or type(embeddings["available_dimensions"][0]) is not int
    ):
        raise ContractError(f"{path} available dimensions differ.")
    _exact_int(
        embeddings["default_dimension"], dimension, f"{path} default dimension"
    )
    if embeddings["files"] != {f"{dimension}d": point_relative}:
        raise ContractError(f"{path} embedding files differ.")
    point_payload = _claim_payload(
        dataset_root,
        point_relative,
        expected["n_cells"] * dimension * 4,
        expected_files,
    )
    points = _finite_float32(
        point_payload,
        expected["n_cells"] * dimension,
        point_relative,
        dimension=dimension,
        require_each_axis_variation=True,
    )

    if expected["vector"]:
        vector_relative = f"vectors/velocity_umap_{dimension}d.bin.gz"
        vector_fields = _exact_keys(
            identity["vector_fields"],
            {"default_field", "fields"},
            f"{path} vector fields",
        )
        field_map = _exact_keys(
            vector_fields["fields"],
            {"velocity_umap"},
            f"{path} vector field map",
        )
        field = _exact_keys(
            field_map["velocity_umap"],
            {
                "label",
                "basis",
                "available_dimensions",
                "default_dimension",
                "files",
            },
            f"{path} velocity field",
        )
        if (
            vector_fields["default_field"] != "velocity_umap"
            or field["label"] != "velocity_umap"
            or field["basis"] != "umap"
            or field["available_dimensions"] != [dimension]
            or type(field["available_dimensions"][0]) is not int
            or field["files"] != {f"{dimension}d": vector_relative}
        ):
            raise ContractError(f"{path} vector field metadata differs.")
        _exact_int(
            field["default_dimension"], dimension, f"{path} vector dimension"
        )
        vector_payload = _claim_payload(
            dataset_root,
            vector_relative,
            expected["n_cells"] * dimension * 4,
            expected_files,
        )
        vectors = _finite_float32(
            vector_payload,
            expected["n_cells"] * dimension,
            vector_relative,
            dimension=dimension,
            require_each_axis_variation=False,
        )
        if not any(value != 0.0 for value in vectors):
            raise ContractError(f"{vector_relative} is a zero vector field.")
    return points


def _validate_dataset(dataset_root: Path, expected: dict[str, Any]) -> set[str]:
    if dataset_root.is_symlink() or not dataset_root.is_dir():
        raise ContractError(f"Dataset root must be a real directory: {dataset_root}.")
    expected_files = {
        "connectivity_manifest.json",
        "dataset_identity.json",
        "obs_manifest.json",
        "var_manifest.json",
    }
    points = _validate_identity(dataset_root, expected, expected_files)
    _validate_obs(dataset_root, expected, points, expected_files)
    _validate_var(dataset_root, expected, expected_files)
    _validate_connectivity(dataset_root, expected, expected_files)

    actual_entries = _scan_real_tree(dataset_root)
    actual_files = {
        relative for relative, kind in actual_entries.items() if kind == "file"
    }
    if actual_files != expected_files:
        raise ContractError(
            f"{dataset_root} artifact inventory differs: "
            f"missing={sorted(expected_files - actual_files)}, "
            f"extra={sorted(actual_files - expected_files)}."
        )
    return expected_files


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    try:
        with path.open("rb") as handle:
            while chunk := handle.read(131_072):
                digest.update(chunk)
    except OSError as error:
        raise ContractError(f"Cannot hash export artifact {path}.") from error
    return digest.hexdigest()


def _validate_checksums(exports_root: Path, expected_paths: set[str]) -> None:
    checksum_path = exports_root / "SHA256SUMS"
    try:
        checksum_bytes = checksum_path.read_bytes()
        checksum_text = checksum_bytes.decode("ascii")
    except (OSError, UnicodeError) as error:
        raise ContractError("SHA256SUMS must be readable ASCII.") from error
    if (
        not checksum_text.endswith("\n")
        or checksum_text.endswith("\n\n")
        or "\r" in checksum_text
    ):
        raise ContractError("SHA256SUMS must use one LF per record.")
    lines = checksum_text[:-1].split("\n")
    if len(lines) != CHECKSUM_RECORDS:
        raise ContractError(
            f"SHA256SUMS must contain exactly {CHECKSUM_RECORDS} records."
        )

    records: dict[str, str] = {}
    for line in lines:
        match = SHA256_RECORD.fullmatch(line)
        if match is None:
            raise ContractError("SHA256SUMS contains a malformed record.")
        digest, relative = match.groups()
        _safe_relative(relative)
        if relative == "SHA256SUMS":
            raise ContractError("SHA256SUMS must not checksum itself.")
        if relative in records:
            raise ContractError(f"SHA256SUMS duplicates {relative}.")
        records[relative] = digest
    if list(records) != sorted(records):
        raise ContractError("SHA256SUMS paths must be sorted.")
    _assert_casefold_unique(records, "SHA256SUMS")
    if set(records) != expected_paths:
        raise ContractError(
            "SHA256SUMS inventory differs: "
            f"missing={sorted(expected_paths - records.keys())}, "
            f"extra={sorted(records.keys() - expected_paths)}."
        )
    for relative, digest in records.items():
        path = exports_root.joinpath(*PurePosixPath(relative).parts)
        if _sha256(path) != digest:
            raise ContractError(f"SHA256 mismatch for {relative}.")


def validate_exports(exports_root: Path = DEFAULT_EXPORTS) -> dict[str, int]:
    candidate = Path(exports_root)
    if candidate.is_symlink():
        raise ContractError(f"Exports root must not be a symlink: {candidate}.")
    try:
        root = candidate.resolve(strict=True)
    except OSError as error:
        raise ContractError(f"Exports root does not exist: {candidate}.") from error
    tree_entries = _scan_real_tree(root)

    catalog_path = root / "datasets.json"
    catalog = _exact_keys(
        _load_json(catalog_path),
        {"version", "default", "datasets"},
        f"{catalog_path} catalog",
    )
    _exact_int(catalog["version"], 1, f"{catalog_path} version")
    if catalog["default"] != "synthetic-cell-types-2d":
        raise ContractError(f"{catalog_path} default dataset differs.")
    entries = catalog["datasets"]
    expected_ids = sorted(EXPECTED)
    if type(entries) is not list or len(entries) != len(expected_ids):
        raise ContractError(f"{catalog_path} dataset count differs.")

    actual_ids = [
        entry.get("id") if type(entry) is dict else None for entry in entries
    ]
    if actual_ids != expected_ids:
        raise ContractError(f"{catalog_path} dataset ordering differs.")
    _assert_casefold_unique(expected_ids, "Catalog dataset IDs")

    expected_paths = {"datasets.json"}
    total_cells = 0
    total_genes = 0
    for entry, dataset_id in zip(entries, expected_ids, strict=True):
        record = _exact_keys(
            entry,
            {"id", "path", "name", "description", "n_cells", "n_genes"},
            f"{catalog_path} entry {dataset_id}",
        )
        expected = EXPECTED[dataset_id]
        if (
            record["id"] != dataset_id
            or _safe_relative(record["path"], directory=True) != f"{dataset_id}/"
            or record["name"] != expected["name"]
            or record["description"] != expected["description"]
        ):
            raise ContractError(f"Catalog entry {dataset_id!r} metadata differs.")
        _exact_int(
            record["n_cells"], expected["n_cells"], f"Catalog {dataset_id} n_cells"
        )
        _exact_int(
            record["n_genes"],
            len(expected["genes"]),
            f"Catalog {dataset_id} n_genes",
        )
        files = _validate_dataset(root / dataset_id, expected)
        expected_paths.update(f"{dataset_id}/{relative}" for relative in files)
        total_cells += expected["n_cells"]
        total_genes += len(expected["genes"])

    expected_root_entries = set(expected_ids) | {"datasets.json", "SHA256SUMS"}
    actual_root_entries = {path.name for path in root.iterdir()}
    if actual_root_entries != expected_root_entries:
        raise ContractError(
            f"Exports root entries differ: "
            f"missing={sorted(expected_root_entries - actual_root_entries)}, "
            f"extra={sorted(actual_root_entries - expected_root_entries)}."
        )

    expected_tree_files = expected_paths | {"SHA256SUMS"}
    actual_tree_files = {
        relative for relative, kind in tree_entries.items() if kind == "file"
    }
    if actual_tree_files != expected_tree_files:
        raise ContractError(
            "Complete export inventory differs: "
            f"missing={sorted(expected_tree_files - actual_tree_files)}, "
            f"extra={sorted(actual_tree_files - expected_tree_files)}."
        )
    if len(expected_paths) != CHECKSUM_RECORDS:
        raise ContractError("Internal checksum inventory count differs.")
    _validate_checksums(root, expected_paths)
    return {
        "datasets": len(EXPECTED),
        "cells": total_cells,
        "genes": total_genes,
        "files": len(expected_tree_files),
        "checksums": len(expected_paths),
    }


def _parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate every committed Cellucid example artifact."
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
        f"{result['genes']} genes, {result['files']} files, "
        f"{result['checksums']} checksums"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (ContractError, OSError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1) from error
