"""Native coupling mapping and effective tbH tabulation tests (no fake exclusion)."""
import math
import unittest
from unittest.mock import patch

from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.charged_effective_production import (
    native_charged_couplings, production_reference,
)


class ChargedEffectiveProductionTests(unittest.TestCase):
    def test_native_spheno_mix_gives_two_orthogonal_doublet_fractions(self):
        for order in ("NO","IO"):
            p=native_charged_couplings(audit(order))
            q=p["mass_eigenstates"]
            self.assertEqual(set(q),{37,9000037})
            self.assertAlmostEqual(
                sum(row["native_doublet_amplitude"]**2 for row in q.values()),
                1.,delta=3e-3)
            self.assertAlmostEqual(p["tan_beta_native"],4.69279,delta=.001)
            for row in q.values():
                self.assertAlmostEqual(
                    row["kappa_t"]/row["kappa_b"],
                    1/p["tan_beta_native"]**2,delta=1e-9)
                self.assertGreater(row["mass_GeV"],200)
                self.assertEqual(row["br_top_to_Hpb"],0)

    def test_backend_quadratic_doublet_scaling_and_both_states(self):
        src=audit("NO")
        def table(mass,kt,kb,top_br):
            self.assertEqual(top_br,0)
            return (kt*kt+0.1*kt*kb+0.3*kb*kb)/mass
        a=production_reference(src,table)
        self.assertFalse(a["full_Zee_charged_production_complete"])
        self.assertEqual(set(a["charged_mass_eigenstates"]),{37,9000037})
        self.assertTrue(all(v>0 for v in a["higgstools_cross_sections"].values()))
        self.assertIsNone(a["full_HiggsBounds_95CL"])

    def test_backend_invalid_rate_is_not_pass(self):
        with self.assertRaisesRegex(ValueError,"nonfinite/invalid"):
            production_reference(audit("NO"),lambda *args:float("nan"))
        with self.assertRaisesRegex(ValueError,"quadratic"):
            production_reference(audit("NO"),lambda mass,kt,kb,br:1.)

    def test_13TeV_reference_does_not_report_full_constraint(self):
        a=production_reference(audit("IO"),
             lambda mass,kt,kb,br: .01*(kt*kt+kb*kb))
        self.assertEqual(a["globally_allowed"],None)
        self.assertIsNone(a["full_HiggsBounds_95CL"])
        self.assertFalse(a["full_Zee_charged_production_complete"])


if __name__=="__main__":
    unittest.main()
