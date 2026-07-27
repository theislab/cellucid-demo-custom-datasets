"""Adversarial tests for the committed custom-dataset repository."""

from __future__ import annotations

import json
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.validate_exports import ContractError, validate_exports

ROOT = Path(__file__).resolve().parents[1]
EXPORTS = ROOT / "exports"


class ExportContractTests(unittest.TestCase):
    def test_committed_tree_is_valid(self) -> None:
        self.assertEqual(
            validate_exports(EXPORTS),
            {"datasets": 3, "cells": 216, "genes": 24, "files": 72},
        )

    def _copy_exports(self, temporary: str) -> Path:
        candidate = Path(temporary) / "exports"
        shutil.copytree(EXPORTS, candidate)
        return candidate

    def test_unexpected_artifact_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / "undeclared.bin").write_bytes(b"not declared")
            with self.assertRaisesRegex(ContractError, "root entries differ"):
                validate_exports(candidate)

    def test_catalog_path_traversal_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            catalog_path = candidate / "datasets.json"
            catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            catalog["datasets"][0]["path"] = "../outside/"
            catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaisesRegex(ContractError, "canonical and safe"):
                validate_exports(candidate)

    def test_corrupt_payload_is_rejected_before_checksum_acceptance(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            payload = candidate / "synthetic-trajectory-1d" / "points_1d.bin.gz"
            bytes_value = bytearray(payload.read_bytes())
            bytes_value[-1] ^= 0x01
            payload.write_bytes(bytes_value)
            with self.assertRaises(ContractError):
                validate_exports(candidate)

    def test_noncanonical_catalog_path_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            catalog_path = candidate / "datasets.json"
            catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            catalog["datasets"][0]["path"] = "synthetic-cell-types-2d//"
            catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaisesRegex(ContractError, "canonical and safe"):
                validate_exports(candidate)

    def test_duplicate_identity_field_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            identity_path = (
                candidate / "synthetic-cell-types-2d" / "dataset_identity.json"
            )
            identity = json.loads(identity_path.read_text(encoding="utf-8"))
            identity["obs_fields"].append(identity["obs_fields"][0])
            identity_path.write_text(json.dumps(identity), encoding="utf-8")
            with self.assertRaisesRegex(ContractError, "exact field array"):
                validate_exports(candidate)

    def test_malformed_catalog_entry_is_a_contract_error(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            catalog_path = candidate / "datasets.json"
            catalog = json.loads(catalog_path.read_text(encoding="utf-8"))
            catalog["datasets"][0] = None
            catalog_path.write_text(json.dumps(catalog), encoding="utf-8")
            with self.assertRaisesRegex(ContractError, "dataset ordering differs"):
                validate_exports(candidate)

    def test_duplicate_json_key_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            candidate = self._copy_exports(temporary)
            (candidate / "datasets.json").write_text(
                '{"version":1,"version":1,"default":"synthetic-cell-types-2d",'
                '"datasets":[]}',
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ContractError, "duplicate key"):
                validate_exports(candidate)


if __name__ == "__main__":
    unittest.main()
