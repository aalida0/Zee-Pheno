"""Regression for legacy declared-potential theory flags and source warning."""
import unittest
import math

from Model.model import scalar_parameters
from Scan.theory_gates import check_potential


class DeclaredTheoryTests(unittest.TestCase):
    def setUp(self):
        self.scalar=dict(mh1=125.25,mh2=690.,ma=735.,mH1=666.,mH2=1184.,
                         alpha=-.18,beta=math.atan(4.69),gamma=.22,v=246.21965)
        self.inv=dict(m12sq_GeV2=98976.1653,lambda_c=.18029,
                      K1=.087475,K2=.5648,K12=-.02893)

    def test_declared_only_full_tree_unitarity_and_BFB(self):
        p=scalar_parameters(**self.scalar,m12sq=self.inv["m12sq_GeV2"],
                            lc=self.inv["lambda_c"],k1=self.inv["K1"],
                            k2=self.inv["K2"],k12=self.inv["K12"])
        declared={k:p[k] for k in ("Lambda1","Lambda2","Lambda3",
                  "Lambda4","Lambda5","K1","K2","K12","Mu4")}
        result=check_potential(self.scalar,declared,self.inv)
        self.assertEqual(result["internal_tree_scalar_inversion"],"PASS")
        self.assertIn(result["full_Zee_BFB"],
                      ("CERTIFIED_SUFFICIENT_PASS","FAIL_NECESSARY","UNKNOWN"))
        self.assertIn(result["full_Zee_unitarity"],("PASS","FAIL"))
        self.assertGreater(result["full_Zee_tree_unitarity_eigenvalue"],0)
        self.assertEqual(result["native_scalar_potential_provenance"],"UNKNOWN")
        self.assertEqual(result["full_Zee_vacuum_globality"],"UNKNOWN")
        declared["Lambda1"]+=.03
        with self.assertRaisesRegex(ValueError,"inconsistent"):
            check_potential(self.scalar,declared,self.inv)

    def test_missing_singlet_inversion_cannot_be_synthesized(self):
        p=scalar_parameters(**self.scalar,m12sq=self.inv["m12sq_GeV2"],
                            lc=self.inv["lambda_c"],k1=self.inv["K1"],
                            k2=self.inv["K2"],k12=self.inv["K12"])
        declared={k:p[k] for k in ("Lambda1","Lambda2","Lambda3",
                  "Lambda4","Lambda5","K1","K2","K12","Mu4")}
        inv=dict(self.inv);inv.pop("lambda_c")
        with self.assertRaises(ValueError):
            check_potential(self.scalar,declared,inv)


if __name__=="__main__":
    unittest.main()
