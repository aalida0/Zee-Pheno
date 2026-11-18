"""Restricted analytic proposal catalog: no accidental native pass."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Model.audit_neutrino_native import sha
from Scan.proposal_catalog import build

ROOT=Path(__file__).resolve().parents[1]

def template_card():
    parts=[]
    for block in ("FZEEIN","IMFZEEIN","YEE2IN","IMYEE2IN"):
        parts.append(f"BLOCK {block}\n")
        for i in range(1,4):
            for j in range(1,4):
                parts.append(f"  {i} {j} 0.0 # fixture\n")
    return "".join(parts)

def catalog(source):
    return {
        "schema_version":1,
        "snapshot":str(ROOT/"Model/reference/native_neutrino_snapshot.json"),
        "reference_card":str(source),
        "targets":[
            {"point_id":"NO_a","ordering":"NO","mlight_eV":.0093,
             "alpha21_rad":0.,"alpha31_rad":0.},
            {"point_id":"IO_b","ordering":"IO","mlight_eV":.001,
             "alpha21_rad":1.,"alpha31_rad":-1.},
        ]}

class ProposalTests(unittest.TestCase):
    def test_no_io_candidates_are_analytic_only_and_immutable(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source=root/"template.in"
            source.write_text(template_card())
            cfg=root/"catalog.json"
            cfg.write_text(json.dumps(catalog(source)))
            with patch("Scan.proposal_catalog.REFERENCE_CARD_SHA",sha(source)):
                result=build(cfg,root/"out")
            self.assertEqual(result["n_proposals"],2)
            self.assertEqual(result["status"],"ANALYTIC_PROPOSALS_NOT_NATIVE_POINTS")
            for item in result["all_proposals"]:
                self.assertEqual(item["status"],"NATIVE_REVALIDATION_REQUIRED")
                self.assertEqual(item["scalar_and_LFV_gates"],"NOT_EVALUATED")
                self.assertTrue((root/"out"/item["card"]).exists())
                self.assertEqual(item["card_sha256"],sha(root/"out"/item["card"]))
                self.assertLess(item["analytic_forward_residual"],1e-8)
            self.assertFalse((root/"out"/"NO_a"/"result.json").exists())
            self.assertEqual(source.read_text(),template_card())

    def test_requires_attested_template_and_unique_safe_targets(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            source=root/"template.in"
            source.write_text(template_card())
            cfg=root/"catalog.json"
            payload=catalog(source)
            cfg.write_text(json.dumps(payload))
            with self.assertRaisesRegex(ValueError,"Unrecognized"):
                build(cfg,root/"out")
            payload["targets"][1]["point_id"]="NO_a"
            cfg.write_text(json.dumps(payload))
            with patch("Scan.proposal_catalog.REFERENCE_CARD_SHA",sha(source)):
                with self.assertRaisesRegex(ValueError,"Duplicate"):
                    build(cfg,root/"out")
            payload["targets"][1]["point_id"]="../../unsafe"
            cfg.write_text(json.dumps(payload))
            with patch("Scan.proposal_catalog.REFERENCE_CARD_SHA",sha(source)):
                with self.assertRaisesRegex(ValueError,"safe"):
                    build(cfg,root/"out")

if __name__=="__main__":
    unittest.main()
