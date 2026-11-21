"""Manifest-based Higgs screen source gates using real frozen SPheno files."""
import json
import tempfile
import unittest
from pathlib import Path

from ZeeHiggs.spheno_audit import NATIVE
from ZeeHiggs.native_higgs_scan import audit_manifest_source,evaluate_manifest


def entry(order):
    root=NATIVE/f"s00017_{order}_benchmark"
    return {
        "point_id":f"s00017_{order}","ordering":order,
        "card":str(root/"input.in"),
        "native":str(root/"native.json"),
        "spectrum":str(root/"spectrum.slha"),
    }


class HiggsNativeManifestTests(unittest.TestCase):
    def test_source_matched_native_NO_IO_have_per_point_neutral_and_charged_inputs(self):
        for order in ("NO","IO"):
            with self.subTest(order=order):
                p=audit_manifest_source(entry(order),Path("/"))
                self.assertTrue(p["native_theory_gate_pass"] if
                                "native_theory_gate_pass" in p else p["theory_scope_pass"])
                self.assertIn("37",p["masses_GeV"])
                self.assertIn("9000037",p["masses_GeV"])
                self.assertEqual(set(p["neutral_ggF_LHC13_pb"]),{"25","35","36"})
                self.assertGreater(p["neutral_ggF_LHC13_pb"]["25"],0)

    def test_tampered_card_never_becomes_higgs_input(self):
        with tempfile.TemporaryDirectory() as d:
            src=entry("NO")
            card=Path(d)/"tamper.in"
            card.write_text(Path(src["card"]).read_text()+"\n# wrong\n")
            src["card"]=str(card)
            with self.assertRaisesRegex(ValueError,"card|SHA|source"):
                audit_manifest_source(src,Path("/"))

    def test_manifest_records_backend_failure_without_allowed_points(self):
        with tempfile.TemporaryDirectory() as d:
            manifest=Path(d)/"manifest.json"
            manifest.write_text(json.dumps({"schema_version":1,
                                             "points":[entry("NO"),entry("IO")]}))
            # No backend installed in unit test. The source-binding stage
            # must complete; backend failures must remain unevaluated.
            from unittest.mock import patch
            with patch("ZeeHiggs.native_higgs_scan.evaluate_one",
                       side_effect=ValueError("backend unavailable")):
                r=evaluate_manifest(manifest,None,None,None)
            self.assertEqual(r["n_points"],2)
            self.assertEqual(r["n_scoped_evaluated"],0)
            self.assertEqual(r["globally_allowed"],None)
            self.assertTrue(all(p["status"]==
                 "NOT_EVALUATED_SOURCE_OR_BACKEND_ERROR" for p in r["points"]))


if __name__=="__main__":
    unittest.main()
