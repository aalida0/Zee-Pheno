"""Step 2 scalar-only LFV physics tests; all exclusion claims fail closed."""
import json
import tempfile
import unittest
from pathlib import Path

import numpy as np

from LFV.core import (kernels, neutral_dipole, radiative_width, combine,
                      dipole_ledger)
from LFV.pilot import evaluate
from Neutrino.fit import complex_matrix


REFERENCE=Path(__file__).resolve().parents[1]/"Model/reference/native_neutrino_snapshot.json"


class DipoleTests(unittest.TestCase):
    def test_Lavoura_kernel_limits_and_vanishing_flavor(self):
        self.assertEqual(kernels(1.),(1/24,1/3,1/24))
        c=np.diag([.01,.02,.03]).astype(complex)
        amp=neutral_dipole(c,250.,1,0,np.array([.0005,.104,1.75]))
        np.testing.assert_allclose(amp,0.,atol=1e-28)
        self.assertEqual(radiative_width(np.zeros(2,complex),.104,.0005,.0073),0.)

    def test_external_pole_override_recovers_old_equal_mass_scheme(self):
        ml=np.array([.0005,.104,1.75])
        c=np.array([[.001,.0003j,.0001],[-.0001j,.003,.0002j],
                    [.0004,-.0008j,.009]],complex)
        np.testing.assert_array_equal(
            neutral_dipole(c,450.,1,0,ml),
            neutral_dipole(c,450.,1,0,ml,external=ml))
        p=dict(v=246.,alpha=1.39,beta=.8,gamma=.6,
               mh1=125.,mh2=600.,ma=650.,mH1=675.,mH2=1167.)
        ff=np.array([[0,1e-5,2e-5],[-1e-5,0,3e-5],[-2e-5,-3e-5,0]],complex)
        x=dipole_ledger(ff,c,p,ml,1,0)
        y=dipole_ledger(ff,c,p,ml,1,0,external=ml)
        np.testing.assert_array_equal(combine(x),combine(y))

    def test_physical_status_unknown_and_input_binding(self):
        native=json.loads(REFERENCE.read_text())
        snapshot=native["snapshot"]
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            card=root/"card.in"
            card.write_text("test reference only\n")
            from Model.audit_neutrino_native import sha
            packed={"input_sha256":sha(card),"spectrum_snapshot_index":0,
                    "snapshots":[snapshot]}
            payload=root/"native.json"
            payload.write_text(json.dumps(packed))
            spectrum=root/"spectrum.spc"
            spectrum.write_text("BLOCK MASS\n 11 0.00051099893\n 13 0.105658372\n 15 1.77686\n 25 125.0\n 35 600.\n 36 640.\n")
            result=evaluate(payload,card,spectrum,alpha_rad=1.398)
            self.assertEqual(result["status"],"SCALAR_ONE_LOOP_SUBSET_ONLY")
            self.assertEqual(result["physical_lfv"],"UNKNOWN")
            self.assertIsNone(result["global_zee_acceptance"])
            self.assertGreater(len(result["missing_contributions"]),5)
            for v in result["processes"].values():
                self.assertTrue(np.isfinite(v["BR_scalar_one_loop_subset"]))
                self.assertGreaterEqual(v["BR_scalar_one_loop_subset"],0.)
                self.assertEqual(len(v["component_CL_CR"]),7)
                pairs=v["component_CL_CR"].values()
                for chirality in ("C_L","C_R"):
                    total=complex(
                        v[chirality+"_GeV_minus2"]["real"],
                        v[chirality+"_GeV_minus2"]["imag"])
                    summed=sum(complex(x[chirality]["real"],x[chirality]["imag"])
                               for x in pairs)
                    self.assertLess(abs(total-summed),1e-25)
            card.write_text("changed\n")
            with self.assertRaises(ValueError):
                evaluate(payload,card,spectrum,alpha_rad=1.398)
            with self.assertRaises(ValueError):
                evaluate(payload,card,spectrum,alpha_rad=float("nan"))


if __name__=="__main__":
    unittest.main()
