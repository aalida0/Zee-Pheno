"""Real native and adversarial manifest tests of final pre-scan theory gate."""
import copy
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from ZeeHiggs.spheno_audit import NATIVE,read_slha
from Scan.theory_gate_native import (
    screen_blocks,screen_one,evaluate_manifest,ORIGINAL_BINARY_SHA256
)


class PreScanTheoryGateTests(unittest.TestCase):
    def test_real_native_NO_IO_scoped_complete_theory(self):
        for order in ("NO","IO"):
            with self.subTest(order=order):
                dir=NATIVE/f"s00017_{order}_benchmark"
                entry={"point_id":"s00017_"+order,"ordering":order,
                       "card":str(dir/"input.in"),
                       "native":str(dir/"native.json"),
                       "spectrum":str(dir/"spectrum.slha")}
                result=screen_one(entry,Path("/"))
                self.assertEqual(result["status"],
                                 "DECLARED_LOW_SCALE_THEORY_SCOPE_PASS")
                self.assertTrue(result["passed_in_declared_tree_and_one_loop_scalar_scope"])
                self.assertTrue(all(result["requirements"].values()))
                self.assertEqual(
                    result["BFB"]["full_Zee_quartic_BFB"],
                    "CERTIFIED_SUFFICIENT_PASS")
                self.assertEqual(
                    result["scalar_unitarity"]["tree_high_energy_full_Zee_unitarity"],
                    "PASS")
                self.assertTrue(result["scalar_one_loop_STU"][
                    "source_match_within_declared_absolute_tolerances"])
                self.assertTrue(result["global_tree_vacuum"][
                    "defined_tree_potential_certificate"]["certified_global_tree_EW"])
                self.assertEqual(result["global_tree_vacuum"][
                    "original_native_running_potential_certificate"]["status"],
                    "UNKNOWN_SOURCE_TADPOLE_SCHEME_MISMATCH")
                self.assertEqual(result["physics_scope"]["Yukawa_RGE_above_source_Q"],
                                 "NOT_REQUESTED")
                self.assertEqual(result["physics_scope"][
                    "native_loop_corrected_scalar_potential_globality"],"UNKNOWN")
                self.assertEqual(result["full_physical_model_viability"],"UNKNOWN")
                self.assertIsNone(result["globally_allowed"])

    def test_high_UV_cutoff_fails_closed_instead_of_assuming_RGE(self):
        d=NATIVE/"s00017_NO_benchmark"
        blocks,_=read_slha((d/"spectrum.slha").read_text())
        result=screen_blocks(blocks,point_id="NO-highscale",cutoff_Q_GeV=1e4)
        self.assertFalse(result["passed_in_declared_tree_and_one_loop_scalar_scope"])
        self.assertFalse(result["requirements"]["optional_UV_RGE_cutoff_if_requested"])
        self.assertEqual(result["physics_scope"]["Yukawa_RGE_above_source_Q"],
                         "UNKNOWN_REQUIRES_MODEL_SOURCE_RGE_TRAJECTORY")

    def test_manifest_preserves_real_card_hash_and_invalid_source(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            card=NATIVE/"s00017_NO_benchmark/input.in"
            native=NATIVE/"s00017_NO_benchmark/native.json"
            spectrum=NATIVE/"s00017_NO_benchmark/spectrum.slha"
            rows=[
                {"point_id":"n","ordering":"NO","card":str(card),
                 "native":str(native),"spectrum":str(spectrum)},
                {"point_id":"tampered","ordering":"NO",
                 "card":str(root/"bad.in"),"native":str(native),"spectrum":str(spectrum)},
            ]
            (root/"bad.in").write_text(card.read_text()+"\n# tamper\n")
            (root/"manifest.json").write_text(json.dumps({
                "schema_version":1,"points":rows
            }))
            output=evaluate_manifest(root/"manifest.json")
            self.assertEqual(output["n_points"],2)
            self.assertEqual(output["n_theory_scope_pass"],1)
            self.assertEqual(output["n_source_or_convention_invalid"],1)
            self.assertTrue(output["points"][0]["provenance"]["source_card_matched_native"])
            self.assertEqual(output["points"][1]["status"],
                             "NOT_EVALUATED_SOURCE_OR_CONVENTION_ERROR")
            self.assertEqual(output["n_global_allowed"],None)

    def test_insufficient_native_data_is_not_silently_a_pass(self):
        file=NATIVE/"s00017_NO_benchmark/spectrum.slha"
        blocks,_=read_slha(file.read_text())
        b=copy.deepcopy(blocks)
        b.pop("YEE2")
        with self.assertRaisesRegex(ValueError,"Missing complex Yukawa"):
            screen_blocks(b)
        b=copy.deepcopy(blocks)
        # Scalar source has moved, but MINPAR was not updated.
        b["ZEEQUARTIC"][0][0][1]="40.0"
        with self.assertRaisesRegex(ValueError,"Source MINPAR"):
            screen_blocks(b)
        b=copy.deepcopy(blocks)
        b.pop("SPHENOLOWENERGY")
        with self.assertRaises(Exception):
            screen_blocks(b)


if __name__=="__main__":unittest.main()
