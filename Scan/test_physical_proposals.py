"""Integration tests for scalar-first, multiple-neutrino-family proposals."""
import copy
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from Model.audit_neutrino_native import sha
from Model.check_native import blocks,mat
from Scan.physical_proposals import generate, cheap_scalar

ROOT=Path(__file__).resolve().parents[1]


def template():
    names=[("MINPAR",3),("ZEEPOT",5),("ZEEQUARTIC",9)]
    out=[]
    for name,n in names:
        out.append(f"Block {name}\n")
        for k in range(1,n+1):
            if name=="MINPAR" and k!=3:
                continue
            out.append(f" {k} 0.0 # source\n")
    for name in ("FZEEIN","IMFZEEIN","YEE2IN","IMYEE2IN"):
        out.append(f"BLOCK {name}\n")
        for i in range(1,4):
            for j in range(1,4):
                out.append(f" {i} {j} 0.0 # source\n")
    return "".join(out)


def arm(name="baseline",frac=.7):
    return {"name":name,"weight":1,"ranges":{
        "tan_beta":[2.,2.],"gamma_rad":[.35,.35],
        "mH2_GeV":[900.,900.],"mH1_fraction":[frac,frac],
        "mH_offset_GeV":[0.,0.],"mA_offset_GeV":[0.,0.],
        "alignment_delta_rad":[0.,0.],"soft_mass_ratio":[1.,1.],
        "lambda_c":[.3,.3],"K1":[.2,.2],"K2":[.2,.2],
        "K12":[.01,.01]}}


def config():
    return {"schema_version":2,"seed":117,
            "scalar_count":2,"families_per_scalar":4,
            "arms":[arm("baseline"),arm("invalid_gap",.95)],
            "neutrino":{"orderings":["NO","IO"],"ranges":{
                "lightest_eV":[.0093,.0093],
                "Majorana_phase_rad":[.2,.2],
                "log10_F_scale":[-4.,-4.],
                "log10_homogeneous_scale":[-4.,-4.],
                "null_t":[.25,.25]}}}


class PhysicalProposalsTests(unittest.TestCase):
    def test_scalar_denominator_full_yukawa_family_and_native_pending(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            cfg=root/"config.json"
            original=config()
            cfg.write_text(json.dumps(original))
            card=root/"reference.in"
            card.write_text(template())
            snap=ROOT/"Model/reference/native_neutrino_snapshot.json"
            with patch("Scan.physical_proposals.REFERENCE_CARD_SHA",sha(card)):
                report=generate(cfg,card,snap,root/"out")
            self.assertEqual(report["n_scalar_attempted"],2)
            self.assertEqual(report["n_scalar_preselected"],1)
            self.assertEqual(report["n_native_cards_written"],4)
            self.assertEqual(report["n_native_points_validated"],0)
            self.assertIsNone(report["n_globally_allowed"])
            self.assertEqual([r["arm"] for r in report["outputs"]],
                             ["baseline","invalid_gap"])
            self.assertEqual(report["outputs"][1]["gates"]["status"],
                             "FAIL_PREDECLARED_KINEMATIC_OR_MASS_PRIOR")
            families=report["outputs"][0]["families"]
            self.assertEqual([(f["ordering"],f["root"]) for f in families],
                             [("NO",0),("IO",0),("NO",1),("IO",1)])
            self.assertEqual(families[0]["homogeneous_coefficients"]["real"],
                             [[0.,0.,0.,0.]])
            self.assertEqual(families[1]["homogeneous_coefficients"]["real"],
                             [[0.,0.,0.,0.]])
            for point in report["outputs"][0]["families"]:
                self.assertLess(point["analytical_closure"],1e-8)
                self.assertEqual(point["status"],"NATIVE_REVALIDATION_REQUIRED")
                self.assertIn(point["ordering"],("NO","IO"))
                path=root/"out"/"cards"/point["card"]
                self.assertEqual(point["card_sha256"],sha(path))
                data=blocks(path.read_text())
                self.assertEqual(len(data["ZEEPOT"]),5)
                self.assertEqual(len(data["ZEEQUARTIC"]),9)
                self.assertAlmostEqual(data["MINPAR"][(3,)],2.)
                f=mat(data,"FZEEIN")
                self.assertAlmostEqual(f[1,2],-f[2,1])
            self.assertGreater(
                abs(report["outputs"][0]["families"][1]["D_mass_basis"]["real"][1][1]-
                    report["outputs"][0]["families"][0]["D_mass_basis"]["real"][1][1]),0.)
            self.assertEqual(json.loads(cfg.read_text()),original)
            self.assertEqual(card.read_text(),template())
            with patch("Scan.physical_proposals.REFERENCE_CARD_SHA",sha(card)):
                with self.assertRaises(FileExistsError):
                    generate(cfg,card,snap,root/"out")

    def test_reject_modified_template_and_undocumented_priors(self):
        with tempfile.TemporaryDirectory() as td:
            root=Path(td)
            card=root/"card.in"
            card.write_text(template())
            snap=ROOT/"Model/reference/native_neutrino_snapshot.json"
            cfg=root/"config.json"
            cfg.write_text(json.dumps(config()))
            with self.assertRaisesRegex(ValueError,"source-attested"):
                generate(cfg,card,snap,root/"bad")
            invalid=config()
            invalid["arms"][0]["ranges"]["unknown"]=[1.,2.]
            cfg.write_text(json.dumps(invalid))
            with patch("Scan.physical_proposals.REFERENCE_CARD_SHA",sha(card)):
                with self.assertRaisesRegex(ValueError,"prior"):
                    generate(cfg,card,snap,root/"bad")
            self.assertFalse((root/"bad").exists())

    def test_sufficient_quartic_bfb_not_called_global_vacuum(self):
        good={"Lambda1":.25,"Lambda2":.3,"Lambda3":.2,
              "Lambda4":-.1,"Lambda5":.01,"LambdaC":.3,
              "K1":.2,"K2":.2,"K12":.01,"Mu4":20.,
              "M112":10.,"M222":20.,"M12Sq":30.,"MCS2":40.}
        valid=cheap_scalar(good)
        self.assertEqual(valid["zee_quartic_BFB_sufficient"],"PASS")
        self.assertEqual(valid["global_vacuum"],"UNKNOWN")
        bad=dict(good,Lambda1=-.1)
        self.assertEqual(cheap_scalar(bad)["doublet_BFB_necessary"],"FAIL")


if __name__=="__main__":
    unittest.main()
