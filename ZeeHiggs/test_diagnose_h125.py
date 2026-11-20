"""Pure-contract tests for official HiggsSignals per-measurement BR diagnosis."""
import unittest
from ZeeHiggs.diagnose_h125 import (
    compare_measurements, particle_report, number
)


class AuditTests(unittest.TestCase):
    def test_correlated_measurement_matching_and_sort_order(self):
        source=[{"id":"100","chi2":3.,"reference":"A","source":"foo","n_subobservables":2},
                {"id":"200","chi2":6.,"reference":"B","source":"bar","n_subobservables":1}]
        eff=[{**source[0],"chi2":5.},{**source[1],"chi2":4.}]
        exp=[{**source[0],"chi2":7.},{**source[1],"chi2":20.}]
        rows=compare_measurements(source,eff,exp)
        self.assertEqual(rows[0]["id"],"200")
        self.assertEqual(rows[0]["delta_explicit_vs_effective"],16.)
        self.assertEqual(rows[0]["delta_explicit_vs_SM"],14.)
        self.assertEqual(rows[1]["delta_effective_vs_SM"],2.)
        self.assertEqual(rows[1]["has_correlations_within_measurement"],
                         "included_by_official_measurement")

    def test_reject_mismatched_measurement_ids_and_duplicates(self):
        x=[{"id":"A","chi2":1.,"reference":"","source":"","n_subobservables":1}]
        y=[{**x[0],"id":"B"}]
        with self.assertRaisesRegex(ValueError,"Inconsistent"):
            compare_measurements(x,x,y)
        with self.assertRaisesRegex(ValueError,"Inconsistent"):
            compare_measurements(x+x,x,x)

    def test_higgs_predicted_width_and_br(self):
        class FakeH:
            def totalWidth(self):return 0.003
            def br(self,d):return {"bb":.5,"WW":.2}.get(d,0.)
        class FakePred:
            def particle(self,id):
                self.assert_id=id
                return FakeH()
        result=particle_report(FakePred(),"25")
        self.assertAlmostEqual(result["partial_widths_GeV"]["bb"],.0015)
        self.assertAlmostEqual(result["branching_ratios"]["WW"],.2)
        self.assertAlmostEqual(result["listed_channel_sum"],.7)

    def test_reject_bad_physical_inputs(self):
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            number(float("nan"),"test")
        class Bad:
            def totalWidth(self):return -1.
            def br(self,d):return 0.
        class Pred:
            def particle(self,id):return Bad()
        with self.assertRaisesRegex(ValueError,"Invalid Higgs width"):
            particle_report(Pred(),"25")


if __name__=="__main__":
    unittest.main()
