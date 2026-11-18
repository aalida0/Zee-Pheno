"""End-to-end native evidence and fail-closed full-scan regression."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from LFV.test_constraints import sample
from Model.audit_neutrino_native import sha
from Model.model import scalar_parameters
from Neutrino.fit import target
from Scan.native_gates import evaluate_entry
from Scan.native_pipeline import evaluate_manifest
from Scan.prepare_reference_manifest import prepare

ROOT=Path(__file__).resolve().parents[1]
REFERENCE=json.loads((ROOT/"Model/reference/native_neutrino_snapshot.json").read_text())["snapshot"]


def fixture(root,ordering,name,higgs_mass=125.,neutrino_bad=False):
    """Synthetic Takagi card for integration tests, never a real SPheno run."""
    base=copy.deepcopy(REFERENCE)
    mnu,m,u=target(ordering,.0093 if ordering=="NO" else .001)
    if neutrino_bad:
        mnu=mnu*1.5
        m=m*1.5
    indices=[0,1,2] if ordering=="NO" else [2,0,1]
    usorted=u[:,indices]
    base["mat1"]={"real":(mnu.real/1e9).tolist(),
                  "imag":(mnu.imag/1e9).tolist()}
    base["UV_1L"]={"real":usorted.conj().T.real.tolist(),
                    "imag":usorted.conj().T.imag.tolist()}
    base["MFv_1L"]=(m[indices]/1e9).tolist()
    base["ZEL"]={"real":np.eye(3).tolist(),
                 "imag":np.zeros((3,3)).tolist()}
    base["ZER"]=base["ZEL"]
    card=root/f"{name}.in"
    card.write_text(f"synthetic native card {name}\n")
    native=root/f"{name}.json"
    native.write_text(json.dumps({
        "input_sha256":sha(card),"snapshots":[base],
        "spectrum_snapshot_index":0,
    }))
    spectrum=root/f"{name}.spc"
    spectrum.write_text(
        "BLOCK MASS\n"
        " 11 0.00051099893\n 13 0.105658372\n 15 1.77686\n"
        f" 25 {higgs_mass:.8E}\n 35 7.67531354E+02\n 36 7.18038367E+02\n"
        "BLOCK SCALARMIX\n"
        " 1 1 7.25651649E-01\n 1 2 6.88062268E-01\n"
        " 2 1 -6.88062268E-01\n 2 2 7.25651649E-01\n")
    mass=sample()["quark_yukawa_masses_GeV"]
    q=root/"quarks.json"
    q.write_text(json.dumps(mass))
    return {
        "point_id":name,"ordering":ordering,
        "card":card.name,"native":native.name,"spectrum":spectrum.name,
        "quark_loop":q.name,"quark_yukawa":q.name,
        "mass_yukawa_scheme":"SYNTHETIC_TEST_ONLY_NOT_MATCHED_QCD",
    }


class NativeGateTests(unittest.TestCase):
    def test_real_native_sources_can_screen_without_arbitrary_quark_masses(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","source_only",125.)
            for key in ("quark_loop","quark_yukawa","mass_yukawa_scheme"):
                entry.pop(key)
            result=evaluate_entry(entry,root)
            self.assertEqual(result["gates"]["neutrino_1D_3sigma"]["status"],"PASS")
            self.assertEqual(result["verdict"],"UNKNOWN_PENDING_GATES")
            self.assertEqual(result["LFV_status"],
                             "SCALAR_ONE_LOOP_ONLY_QUARK_SCHEME_NOT_SUPPLIED")
            self.assertIsNone(result["lfv_point_sha256"])
            self.assertGreaterEqual(result["LFV_retained_BR"]["mu_to_e_gamma"],0.)
            self.assertEqual(result["gates"]["full_Zee_LFV"]["status"],"UNKNOWN")
            bad=dict(entry,quark_loop="quarks.json")
            with self.assertRaisesRegex(ValueError,"Quark"):
                evaluate_entry(bad,root)

    def test_no_io_reference_manifest_writer_without_quark_guessing(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            for ordering,higgs in (("NO",125.),("IO",484.2844)):
                entry=fixture(root,ordering,ordering.lower(),higgs)
                (root/f"step1_{ordering}_proposal.in").write_bytes(
                    (root/entry["card"]).read_bytes())
                dest=root/f"step1_{ordering}_native"
                dest.mkdir()
                (dest/"result.json").write_bytes((root/entry["native"]).read_bytes())
                (dest/"output.spc").write_bytes((root/entry["spectrum"]).read_bytes())
            manifest=root/"reference_manifest.json"
            status=prepare(root,manifest)
            self.assertEqual(status["LFV_scope"],"SCALAR_ONE_LOOP_ONLY")
            self.assertEqual(status["n_points"],2)
            report=evaluate_manifest(manifest)
            self.assertEqual(report["n_native_1D_pass"],2)
            self.assertEqual(report["n_fail_necessary"],1)
            self.assertEqual(report["n_pending"],1)
            self.assertIsNone(report["n_globally_allowed"])

    def test_native_pass_is_not_globally_allowed(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","no_control",125.)
            result=evaluate_entry(entry,root)
            self.assertEqual(result["gates"]["card_native_binding"]["status"],"PASS")
            self.assertEqual(result["gates"]["neutrino_1D_3sigma"]["status"],"PASS")
            self.assertEqual(result["gates"]["h125_mass_screen"]["status"],"PASS")
            self.assertEqual(result["gates"]["cascade_on_shell"]["status"],"PASS")
            self.assertEqual(result["verdict"],"UNKNOWN_PENDING_GATES")
            self.assertIsNone(result["globally_allowed"])
            self.assertEqual(result["gates"]["native_FZEEIN_propagation"]["status"],"UNKNOWN")
            self.assertEqual(result["gates"]["B_flavor"]["status"],"UNKNOWN")
            self.assertEqual(result["gates"]["full_Zee_LFV"]["status"],"UNKNOWN")
            self.assertEqual(result["charged_photon_topology"],"UNKNOWN_POTENTIAL_NOT_SUPPLIED")
            self.assertGreaterEqual(result["LFV_retained_BR"]["mu_to_e_gamma"],0.)

    def test_higgs_mass_fails_even_if_neutrino_fits_and_lfv_small(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"IO","io_heavy_h",484.2844)
            result=evaluate_entry(entry,root)
            self.assertEqual(result["gates"]["neutrino_1D_3sigma"]["status"],"PASS")
            self.assertEqual(result["gates"]["h125_mass_screen"]["status"],"FAIL")
            self.assertEqual(result["verdict"],"FAIL_NECESSARY_SCREEN")
            self.assertFalse(result["globally_allowed"])
            self.assertIn("h125_mass_screen",result["failed_necessary_screens"])

    def test_bad_neutrino_data_fail_screen(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","no_bad",125.,neutrino_bad=True)
            result=evaluate_entry(entry,root)
            self.assertEqual(result["gates"]["neutrino_1D_3sigma"]["status"],"FAIL")
            self.assertEqual(result["verdict"],"FAIL_NECESSARY_SCREEN")

    def test_manifest_isolates_tampered_native_card_without_physics(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            valid=fixture(root,"NO","good",125.)
            invalid=fixture(root,"IO","tampered",125.)
            (root/invalid["card"]).write_text("modified card\n")
            manifest=root/"manifest.json"
            manifest.write_text(json.dumps({
                "schema_version":1,"points":[valid,invalid]}))
            result=evaluate_manifest(manifest)
            self.assertEqual(result["n_entries"],2)
            self.assertEqual(result["n_native_1D_pass"],1)
            self.assertEqual(result["n_pending"],1)
            self.assertEqual(result["n_invalid_source"],1)
            self.assertIsNone(result["n_globally_allowed"])
            self.assertEqual(result["reports"][1]["verdict"],
                             "NOT_EVALUATED_SOURCE_ERROR")
            self.assertNotIn("LFV_retained_BR",result["reports"][1])

    def test_verified_F_input_response_requires_same_card_and_native_checks(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","probed",125.)
            native_report=root/"input_propagation.json"
            evidence={
                "status":"NATIVE_INPUT_GATE_PASS",
                "input_sha256":sha(root/entry["card"]),
                "native_executed":True,
                "binary_sha256":"a"*64,
                "delta_f":.000001,"delta_imye2":.000001,
                "output_sha256":{k:"b"*64 for k in ("baseline","f_probe","ye2_probe")},
                "checks":{"baseline_F_relative":0.,
                          "baseline_Ye2_relative":0.,
                          "F_antisymmetric_max_abs":0.,
                          "F_probe_output_delta":.000001,
                          "Ye2_probe_output_delta":.000001},
            }
            native_report.write_text(json.dumps(evidence))
            entry["native_input_report"]=native_report.name
            report=evaluate_entry(entry,root)
            self.assertEqual(report["gates"]["native_FZEEIN_propagation"]["status"],"PASS")
            evidence["status"]="NATIVE_INPUT_GATE_FAIL"
            native_report.write_text(json.dumps(evidence))
            with self.assertRaisesRegex(ValueError,"INPUT_PROBE_NOT_NATIVE_PASS"):
                evaluate_entry(entry,root)
            evidence["status"]="NATIVE_INPUT_GATE_PASS"
            evidence["input_sha256"]="0"*64
            native_report.write_text(json.dumps(evidence))
            with self.assertRaisesRegex(ValueError,"INPUT_PROBE_REPORT_CARD_SHA_MISMATCH"):
                evaluate_entry(entry,root)

    def test_scalar_inversion_only_necessary_and_source_consistent(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","with_potential",125.)
            base=evaluate_entry(entry,root)
            p=base["native_scalar_parameters"]
            inv={"m12sq_GeV2": 250000.*np.sin(p["beta"])*np.cos(p["beta"]),
                 "lambda_c":.25,"K1":.3,"K2":.3,"K12":.01}
            x=scalar_parameters(**p,m12sq=inv["m12sq_GeV2"],
                                lc=inv["lambda_c"],k1=inv["K1"],
                                k2=inv["K2"],k12=inv["K12"])
            potential={k:x[k] for k in ("Lambda1","Lambda2","Lambda3",
                "Lambda4","Lambda5","K1","K2","K12","Mu4")}
            (root/"inv.json").write_text(json.dumps(inv))
            (root/"pot.json").write_text(json.dumps(potential))
            entry.update(potential="pot.json",potential_source="synthetic source",
                         scalar_inversion="inv.json")
            result=evaluate_entry(entry,root)
            self.assertEqual(result["gates"]["scalar_tree_inversion"]["status"],"PASS")
            self.assertEqual(result["gates"]["potential_provenance"]["status"],"UNKNOWN")
            self.assertEqual(result["gates"]["full_Zee_LFV"]["status"],"UNKNOWN")
            self.assertEqual(result["charged_photon_topology"],
                             "INCLUDED_FROM_RESTRICTED_POTENTIAL")
            potential["Lambda1"]+=1.
            (root/"pot.json").write_text(json.dumps(potential))
            with self.assertRaisesRegex(ValueError,"inconsistent"):
                evaluate_entry(entry,root)

    def test_no_potential_with_or_without_provenance_and_duplicate_ids(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            entry=fixture(root,"NO","no_alter",125.)
            bad=dict(entry,potential="arbitrary.json")
            with self.assertRaises(ValueError):
                evaluate_entry(bad,root)
            manifest=root/"manifest.json"
            manifest.write_text(json.dumps({
                "schema_version":1,"points":[entry,entry]}))
            with self.assertRaisesRegex(ValueError,"Duplicate"):
                evaluate_manifest(manifest)


if __name__=="__main__":
    unittest.main()
