"""Strict synthetic tests + two archived real native spectra.

CI never claims to run HEPfit. A real run requires the local pinned driver.
"""
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from Flavor.bsgamma import (
    native_scalar, ordinary, parse_driver_output, injected_coeff_check,
    run_driver, calculate,
)


def fake_values(br=3.49e-4, mt=164.29, np=None):
    np = np or [[0.]*8 for _ in range(2)]
    a = 0.00868
    p = 0.1
    vals = dict(
        BR_NLO=br, BR_stock_NNLO=br-0.00001,
        cut_GeV=1.6, mt_MSbar_mut=mt, muw_GeV=160.,
        mu_b_GeV=2., mu_c_GeV=1.99999999,
        a_s_muw=a, overall=br / p, P_NLO=p, N=0.,
        high_1_0_7=np[0][6], high_1_0_8=np[0][7],
        high_1_1_4=a*np[1][3],
        high_1_1_7=a*(np[1][6]-4*np[1][3]/9),
        high_1_1_8=a*(np[1][7]-np[1][3]/6),
    )
    return vals


class NativeBridgeTests(unittest.TestCase):
    def test_real_native_pole_masses_and_doublet_projectors(self):
        for order in ("NO", "IO"):
            with self.subTest(order=order):
                p = native_scalar(order)
                self.assertAlmostEqual(p["masses_GeV"]["H1"], 666.305, delta=.002)
                self.assertAlmostEqual(p["masses_GeV"]["H2"], 1184.554, delta=.002)
                self.assertAlmostEqual(sum(p["doublet_weights"]), 1, delta=1e-5)
                self.assertGreater(p["doublet_weights"][0], .97)
                self.assertLess(p["doublet_weights"][1], .03)
                self.assertAlmostEqual(p["legacy_flavorkit_BR"], 3.63471789e-4,
                                       places=10)

    def test_coherent_np_matching_and_raw_order(self):
        point = native_scalar("NO")
        def matching(mass, tanb, mt, mu=160.):
            self.assertAlmostEqual(mu, 160)
            return dict(C7_LO=-300/mass, C8_LO=-400/mass,
                        C4_NLO_raw=100/mass, C7_NLO_raw=200/mass,
                        C8_NLO_raw=250/mass)
        vector, c = ordinary(point, 164., matching)
        direct = sum(w * (-300 / point["masses_GeV"][name])
                     for w,name in zip(point["doublet_weights"], ("H1","H2")))
        self.assertAlmostEqual(c["C7_LO"], direct)
        self.assertEqual(vector[0][6], c["C7_LO"])
        self.assertEqual(vector[0][7], c["C8_LO"])
        self.assertEqual(vector[1][3], c["C4_NLO_raw"])
        self.assertTrue(all(vector[1][i] == 0 for i in (0,1,2,4,5)))
        injected_coeff_check(fake_values(np=vector), vector)
        wrong = fake_values(np=vector)
        wrong["high_1_1_7"] += .1
        with self.assertRaisesRegex(ValueError,"bridge mismatch"):
            injected_coeff_check(wrong, vector)

    def test_output_fails_closed(self):
        vals = fake_values()
        s = "\n".join("VALUE "+k+" "+str(v) for k,v in vals.items())
        self.assertAlmostEqual(parse_driver_output(s)["BR_NLO"], vals["BR_NLO"])
        with self.assertRaisesRegex(ValueError, "silently ignores"):
            parse_driver_output(s+"\nUNKNOWN_INPUT BAD\n")
        with self.assertRaisesRegex(ValueError, "Incomplete"):
            parse_driver_output("VALUE BR_NLO 0.00035\n")
        with self.assertRaisesRegex(ValueError, "Wrong HEPfit photon cut"):
            parse_driver_output(s.replace("VALUE cut_GeV 1.6",
                                          "VALUE cut_GeV 1.9"))

    def test_runner_preserves_raw_receipt_with_zero_exit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fake=SimpleNamespace(returncode=0, stdout="\n".join(
                "VALUE "+k+" "+str(v) for k,v in fake_values().items()))
            with patch("Flavor.bsgamma.subprocess.run",return_value=fake):
                result, receipts=run_driver("sm", "sm", {"muw":160.},
                                            root, root/"fake.exe", {})
            self.assertEqual(receipts["exit_code"],0)
            self.assertTrue((root/"sm.input").exists())
            self.assertEqual(result["BR_NLO"],fake_values()["BR_NLO"])
            with patch("Flavor.bsgamma.subprocess.run",return_value=fake):
                with self.assertRaises(FileExistsError):
                    run_driver("sm", "sm", {"muw":160.},root,root/"fake.exe",{})

    def test_full_orchestration_stays_unknown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            # No native HEPfit is executed in this synthetic test.
            sm_input={"mut":160., "muw":160., "mtop":173.2}
            host=dict(sm_input, mHh2=1., mA2=1.,mHp2=1.,m12_2=1.,
                      bma=1.,logtb=.3,mHl=125.09)
            record={"runs":{"SM_reference":{"input":sm_input},
                            "F13_heavy_0001330_central":{"input":host}}}
            def match(mass,tanb,mt,mu=160):
                return dict(C7_LO=-.04 * (600/mass)**2, C8_LO=-.05*(600/mass)**2,
                            C4_NLO_raw=.001,C7_NLO_raw=.9,C8_NLO_raw=.5)
            def run(label,mode,params,directory,driver,env):
                np = [[params.get(f"NP{o}_{i+1}",0.) for i in range(8)]
                      for o in range(2)]
                br=3.49e-4+(np[0][6]*-1e-4 if label.endswith("_Zee_NLO") else 0)
                return fake_values(br=br,np=np),{"label":label}
            with patch("Flavor.bsgamma.backend",return_value=(record,match,{"synthetic":True})), \
                 patch("Flavor.bsgamma.run_driver",side_effect=run):
                result=calculate(root,root/"driver",root/"source",root/"compat",
                                 root/"rootlib",root/"out")
            self.assertEqual(result["status"],"EXECUTED_CENTRAL_NLO_ONLY_NO_PHYSICAL_EXCLUSION")
            self.assertIsNone(result["globally_allowed"])
            self.assertEqual([p["ordering"] for p in result["points"]],["NO","IO"])
            self.assertTrue(all(p["physical_bsgamma_gate"]=="UNKNOWN"
                                for p in result["points"]))
            self.assertTrue((root/"out/bsgamma_nlo.json").exists())


if __name__ == "__main__":
    unittest.main()
