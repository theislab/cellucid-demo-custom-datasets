"""Adversarial tests for the immutable custom-dataset export contract."""

from __future__ import annotations

import ast
import gzip
import json
import math
import shutil
import struct
import tempfile
import unittest
from pathlib import Path
from unittest import mock

from scripts.validate_exports import (
    ContractError,
    _assert_casefold_unique,
    validate_exports,
)

ROOT = Path(__file__).resolve().parents[1]
EXPORTS = ROOT / "exports"
FINAL_PYTHON_SHA = "eedd3fca1dbb0f57a3ec9468c4a460003bda570a"


class ExportContractTests(unittest.TestCase):
    def _copy_exports(self, temporary: str) -> Path:
        candidate = Path(temporary) / "exports"
        shutil.copytree(EXPORTS, candidate)
        return candidate

    def _load(self, candidate: Path, relative: str) -> dict:
        path = candidate / relative
        return json.loads(path.read_text(encoding="utf-8"))

    def _store(self, candidate: Path, relative: str, value: object) -> None:
        (candidate / relative).write_text(
            json.dumps(value, ensure_ascii=False),
            encoding="utf-8",
        )

    def _rewrite_gzip(
        self,
        candidate: Path,
        relative: str,
        transform,
    ) -> None:
        path = candidate / relative
        payload = bytearray(gzip.decompress(path.read_bytes()))
        transform(payload)
        path.write_bytes(gzip.compress(bytes(payload), compresslevel=6, mtime=0))

    def test_committed_tree_is_valid(self) -> None:
        self.assertEqual(
            validate_exports(EXPORTS),
            {
                "datasets": 3,
                "cells": 216,
                "genes": 24,
                "files": 72,
                "checksums": 71,
            },
        )

    def test_duplicate_json_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / "datasets.json").write_text(
                '{"version":1,"version":1,'
                '"default":"synthetic-cell-types-2d","datasets":[]}',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ContractError, "duplicate key"):
                validate_exports(candidate)

    def test_nonfinite_json_constant_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            text = (candidate / "datasets.json").read_text(encoding="utf-8")
            (candidate / "datasets.json").write_text(
                text.replace('"version": 1', '"version": NaN', 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ContractError, "non-finite JSON constant"):
                validate_exports(candidate)

    def test_overflowing_json_number_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            text = (candidate / "datasets.json").read_text(encoding="utf-8")
            (candidate / "datasets.json").write_text(
                text.replace('"version": 1', '"version": 1e9999', 1),
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ContractError, "non-finite JSON number"):
                validate_exports(candidate)

    def test_enormous_json_integer_is_rejected_without_overflow(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = "synthetic-cell-types-2d/obs_manifest.json"
            manifest = self._load(candidate, relative)
            manifest["_continuousFields"][0][1] = 10**999
            self._store(candidate, relative, manifest)
            with self.assertRaisesRegex(ContractError, "must be a finite number"):
                validate_exports(candidate)

    def test_catalog_path_attacks_are_rejected(self) -> None:
        attacks = (
            "../outside/",
            "/absolute/",
            "synthetic-cell-types-2d//",
            "synthetic-cell-types-2d\\",
            "synthetic-cell-types-2d/?raw=1",
            "synthetic-cell-types-2d/#fragment",
        )
        for attack in attacks:
            with self.subTest(attack=attack), tempfile.TemporaryDirectory() as temporary:
                candidate = self._copy_exports(temporary)
                catalog = self._load(candidate, "datasets.json")
                catalog["datasets"][0]["path"] = attack
                self._store(candidate, "datasets.json", catalog)
                with self.assertRaisesRegex(ContractError, "path"):
                    validate_exports(candidate)

    def test_catalog_order_and_identity_are_exact(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            catalog = self._load(candidate, "datasets.json")
            catalog["datasets"].reverse()
            self._store(candidate, "datasets.json", catalog)
            with self.assertRaisesRegex(ContractError, "dataset ordering"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = "synthetic-cell-types-2d/dataset_identity.json"
            identity = self._load(candidate, relative)
            identity["id"] = "other"
            self._store(candidate, relative, identity)
            with self.assertRaisesRegex(ContractError, "identity metadata"):
                validate_exports(candidate)

    def test_boolean_cannot_impersonate_an_integer(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            catalog = self._load(candidate, "datasets.json")
            catalog["version"] = True
            self._store(candidate, "datasets.json", catalog)
            with self.assertRaisesRegex(ContractError, "integer 1"):
                validate_exports(candidate)

    def test_missing_and_unexpected_artifacts_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (
                candidate
                / "synthetic-trajectory-1d"
                / "obs"
                / "stage.codes.u8.gz"
            ).unlink()
            with self.assertRaisesRegex(ContractError, "real file"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / "undeclared.bin").write_bytes(b"undeclared")
            with self.assertRaisesRegex(ContractError, "root entries differ"):
                validate_exports(candidate)

    def test_casefold_collisions_are_rejected(self) -> None:
        with self.assertRaisesRegex(ContractError, "case-fold collision"):
            _assert_casefold_unique(
                ("dataset/points.bin.gz", "DATASET/POINTS.bin.gz"),
                "test paths",
            )

    def test_symlink_root_is_rejected_when_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            link = Path(temporary) / "exports-link"
            try:
                link.symlink_to(EXPORTS, target_is_directory=True)
            except OSError as error:
                self.skipTest(f"symlinks are unavailable: {error}")
            with self.assertRaisesRegex(ContractError, "must not be a symlink"):
                validate_exports(link)

    def test_symlink_artifact_is_rejected_when_supported(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            artifact = (
                candidate
                / "synthetic-trajectory-1d"
                / "obs"
                / "stage.codes.u8.gz"
            )
            target = Path(temporary) / "target.gz"
            target.write_bytes(artifact.read_bytes())
            artifact.unlink()
            try:
                artifact.symlink_to(target)
            except OSError as error:
                self.skipTest(f"symlinks are unavailable: {error}")
            with self.assertRaisesRegex(ContractError, "must not be a symlink"):
                validate_exports(candidate)

    def test_gzip_timestamp_and_trailing_member_are_rejected(self) -> None:
        relative = "synthetic-trajectory-1d/points_1d.bin.gz"
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            path = candidate / relative
            data = bytearray(path.read_bytes())
            data[4] = 1
            path.write_bytes(data)
            with self.assertRaisesRegex(ContractError, "gzip header"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            path = candidate / relative
            member = path.read_bytes()
            path.write_bytes(member + member)
            with self.assertRaisesRegex(ContractError, "one exact gzip member"):
                validate_exports(candidate)

    def test_gzip_expansion_and_compressed_size_are_bounded(self) -> None:
        relative = "synthetic-trajectory-1d/obs/stage.codes.u8.gz"
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / relative).write_bytes(
                gzip.compress(b"\0" * 49, compresslevel=6, mtime=0)
            )
            with self.assertRaisesRegex(ContractError, "ISIZE"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            high_ratio = bytearray(
                gzip.compress(b"\0" * 500_000, compresslevel=9, mtime=0)
            )
            high_ratio[-4:] = (48).to_bytes(4, "little")
            self.assertLessEqual(len(high_ratio), 48 + 1024)
            (candidate / relative).write_bytes(high_ratio)
            with self.assertRaisesRegex(
                ContractError, "expands beyond its declared size"
            ):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / relative).write_bytes(b"x" * 2048)
            with self.assertRaisesRegex(ContractError, "compressed size"):
                validate_exports(candidate)

    def test_gzip_overflow_never_uses_an_unbounded_flush(self) -> None:
        class OverflowingDecoder:
            eof = False
            unconsumed_tail = b"remaining compressed input"
            unused_data = b""

            @staticmethod
            def decompress(_compressed: bytes, max_length: int) -> bytes:
                return b"\0" * max_length

            @staticmethod
            def flush() -> bytes:
                raise AssertionError("bounded validation must never call flush()")

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            with mock.patch(
                "scripts.validate_exports.zlib.decompressobj",
                return_value=OverflowingDecoder(),
            ):
                with self.assertRaisesRegex(
                    ContractError, "expands beyond its declared size"
                ):
                    validate_exports(candidate)

    def test_nonfinite_embedding_is_rejected_semantically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = "synthetic-trajectory-1d/points_1d.bin.gz"

            def inject_nan(payload: bytearray) -> None:
                struct.pack_into("<f", payload, 0, math.nan)

            self._rewrite_gzip(candidate, relative, inject_nan)
            with self.assertRaisesRegex(ContractError, "non-finite Float32"):
                validate_exports(candidate)

    def test_invalid_category_code_is_rejected_semantically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = "synthetic-trajectory-1d/obs/stage.codes.u8.gz"

            def inject_missing(payload: bytearray) -> None:
                payload[0] = 255

            self._rewrite_gzip(candidate, relative, inject_missing)
            with self.assertRaisesRegex(ContractError, "invalid category code"):
                validate_exports(candidate)

    def test_invalid_connectivity_is_rejected_semantically(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = (
                "synthetic-trajectory-1d/connectivity/edges.src.bin.gz"
            )

            def inject_out_of_bounds(payload: bytearray) -> None:
                struct.pack_into("<H", payload, 0, 65535)

            self._rewrite_gzip(candidate, relative, inject_out_of_bounds)
            with self.assertRaisesRegex(ContractError, "out-of-bounds"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            relative = (
                "synthetic-trajectory-1d/connectivity/"
                "edges.weights.f64.bin.gz"
            )

            def inject_nan(payload: bytearray) -> None:
                struct.pack_into("<d", payload, 0, math.nan)

            self._rewrite_gzip(candidate, relative, inject_nan)
            with self.assertRaisesRegex(ContractError, "finite positive"):
                validate_exports(candidate)

    def test_checksum_mismatch_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            checksum = candidate / "SHA256SUMS"
            text = checksum.read_text(encoding="ascii")
            checksum.write_text(("0" * 64) + text[64:], encoding="ascii")
            with self.assertRaisesRegex(ContractError, "SHA256 mismatch"):
                validate_exports(candidate)

    def test_checksum_format_order_duplicates_and_count_are_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            checksum = candidate / "SHA256SUMS"
            lines = checksum.read_text(encoding="ascii").splitlines()
            lines[0], lines[1] = lines[1], lines[0]
            checksum.write_text("\n".join(lines) + "\n", encoding="ascii")
            with self.assertRaisesRegex(ContractError, "paths must be sorted"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            checksum = candidate / "SHA256SUMS"
            lines = checksum.read_text(encoding="ascii").splitlines()
            lines[-1] = lines[0]
            checksum.write_text("\n".join(lines) + "\n", encoding="ascii")
            with self.assertRaisesRegex(ContractError, "duplicates"):
                validate_exports(candidate)

        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            checksum = candidate / "SHA256SUMS"
            lines = checksum.read_text(encoding="ascii").splitlines()
            checksum.write_text("\n".join(lines[:-1]) + "\n", encoding="ascii")
            with self.assertRaisesRegex(ContractError, "exactly 71 records"):
                validate_exports(candidate)


class RepositoryMaintenanceTests(unittest.TestCase):
    def test_readme_pins_only_the_final_python_revision(self) -> None:
        readme = (ROOT / "README.md").read_text(encoding="utf-8")
        self.assertEqual(readme.count(FINAL_PYTHON_SHA), 1)
        self.assertNotIn("a60eeb4432a7822685fd20bc55da04685e53d5ed", readme)
        self.assertIn("cross-platform maintenance CI", readme)

    def test_workflow_is_the_exact_fast_cross_platform_matrix(self) -> None:
        workflow = (ROOT / ".github" / "workflows" / "test.yml").read_text(
            encoding="utf-8"
        )
        self.assertEqual(workflow.count("uses: actions/checkout@v7"), 1)
        self.assertEqual(workflow.count("uses: actions/setup-python@v7"), 1)
        self.assertEqual(workflow.count("os: ubuntu-latest"), 2)
        self.assertEqual(workflow.count("os: macos-latest"), 1)
        self.assertEqual(workflow.count("os: windows-latest"), 1)
        self.assertEqual(workflow.count('python: "3.11"'), 1)
        self.assertEqual(workflow.count('python: "3.14"'), 3)
        self.assertNotIn("pip install", workflow)

    def test_validator_has_no_remote_or_process_dependencies(self) -> None:
        source = (ROOT / "scripts" / "validate_exports.py").read_text(
            encoding="utf-8"
        )
        tree = ast.parse(source)
        imported_roots = {
            alias.name.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        } | {
            (node.module or "").split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        }
        self.assertTrue(
            imported_roots.isdisjoint(
                {
                    "http",
                    "requests",
                    "socket",
                    "subprocess",
                    "urllib",
                }
            )
        )

    def test_portable_attributes_and_python_ignores_are_exact(self) -> None:
        self.assertEqual(
            (ROOT / ".gitattributes").read_text(encoding="utf-8"),
            "* text=auto eol=lf\n*.bin binary\n*.gz binary\n",
        )
        ignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
        for record in (
            ".pytest_cache/",
            ".ruff_cache/",
            ".venv/",
            "__pycache__/",
            "*.py[cod]",
        ):
            self.assertIn(f"{record}\n", ignore)


if __name__ == "__main__":
    unittest.main()
