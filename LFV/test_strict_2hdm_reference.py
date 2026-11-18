"""Strict 2HDM decoupling mapping and source-bound reference checks."""
import copy
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np

from LFV.test_constraints import sample, cmat
from LFV.strict_2hdm_reference import (SM_KEYS, cp_even_rotation,
    decoupled_input, fill_type_ii_quarks, evaluate)

SM={
    "alphaQED":1/129,"alphasMZ":.1180,"sw2":.23129,
    "MZ":91.1880,"MW":80.3692,"MTpole":172.4,
    "mbmb":4.183,"mcmc":1.273,
    "ms_at_2GeV":.0935,"md_at_2GeV":.0047,"mu_at_2GeV":.00216
}


def decoupled_sample():
    p=sample()
    p["F"]=cmat(np.zeros((3,3),complex))
    p["scalar"]["gamma"]=0.
    p["internal_lepton_masses_GeV"]=p["pole_lepton_masses_GeV"][:]
    p["potential"]={
        "Lambda1":.3,"Lambda2":.45,"Lambda3":.1,
        "Lambda4":-.1,"Lambda5":-.05,
        "K1":0.,"K2":0.,"K12":0.,"Mu4":0.}
    return p


class StrictReferenceTests(unittest.TestCase):
    def test_mixing_quark_yukawa_and_potential_dictionary(self):
        point=decoupled_sample()
        p=point["scalar"]
        inp=decoupled_input(point,SM)
        self.assertEqual(inp["theta12"],p["beta"]-p["alpha"])
        self.assertEqual(inp["theta13"],0.)
        self.assertEqual(inp["theta23"],0.)
        self.assertEqual(inp["MHp"],p["mH1"])
        self.assertEqual(inp["Mh3"],p["ma"])
        self.assertEqual(len([k for k in inp if k.startswith("absrho")]),27)
        D=np.array(point["D"]["real"])+1j*np.array(point["D"]["imag"])
        for i in range(3):
            for j in range(3):
                val=inp[f"absrhol{i+1}{j+1}"]*np.exp(
                    1j*inp[f"argrhol{i+1}{j+1}"])
                self.assertAlmostEqual(abs(val-D[j,i].conjugate()/math.sqrt(2)),0.,
                                       places=14)
        source=dict(beta=p["beta"],mumt=.002,mc=.5,mt=160.,
                    mdmt=.003,msmt=.05,mb=2.7)
        new=fill_type_ii_quarks(inp,source)
        self.assertAlmostEqual(new["absrhou33"],160/p["v"]/math.tan(p["beta"]))
        self.assertAlmostEqual(new["absrhod33"],2.7/p["v"]*math.tan(p["beta"]))
        self.assertAlmostEqual(new["argrhod33"],math.pi)

    def test_invalid_decoupling_and_scheme_fail_closed(self):
        base=decoupled_sample()
        for key,value in (("gamma",.3),("gamma",-1e-3)):
            point=copy.deepcopy(base)
            point["scalar"][key]=value
            with self.assertRaisesRegex(ValueError,"Strict"):
                decoupled_input(point,SM)
        point=copy.deepcopy(base)
        point["F"]["real"][0][1]=1e-4
        point["F"]["real"][1][0]=-1e-4
        with self.assertRaisesRegex(ValueError,"Strict"):
            decoupled_input(point,SM)
        point=copy.deepcopy(base)
        point["potential"]["K1"]=.1
        with self.assertRaisesRegex(ValueError,"Strict"):
            decoupled_input(point,SM)
        point=copy.deepcopy(base)
        point["potential"]["Mu4"]=3.
        with self.assertRaisesRegex(ValueError,"Strict"):
            decoupled_input(point,SM)
        point=copy.deepcopy(base)
        point["internal_lepton_masses_GeV"][1]=.104
        with self.assertRaisesRegex(ValueError,"mass scheme"):
            decoupled_input(point,SM)
        bad=dict(SM)
        bad.pop("MZ")
        with self.assertRaises(ValueError):
            decoupled_input(base,bad)

    def test_fake_public_reference_contract_and_neutral_signs(self):
        # Synthetic API check only; independent numerical validation
        # requires the unmodified author-supplied thdm.py file.
        fake='''\
import math
class thdm:
    def __init__(self, ip):
        self.ip=ip
        self.mumt=.002;self.mc=.5;self.mt=160.
        self.mdmt=.003;self.msmt=.05;self.mb=2.7
    def q1(self,k):
        t=self.ip["theta12"]
        return (math.cos(t),math.sin(t),0)[k-1]
    def q2(self,k):
        t=self.ip["theta12"]
        return (-math.sin(t),math.cos(t),1j)[k-1]
    def cL_llpg(self,i,j):
        return complex(i,j)
    def cR_llpg(self,i,j):
        return complex(-i,j)
'''
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/"fake_author.py"
            path.write_text(fake)
            result=evaluate(decoupled_sample(),SM,path)
            self.assertEqual(result["status"],"STRICT_2HDM_REFERENCE_EVALUATED_NOT_ZEE")
            self.assertIsNone(result["Zee_full_dipoles"])
            self.assertEqual(result["Zee_gauge_extension"],"UNVALIDATED")
            self.assertLess(result["neutral_vertex_max_residual"],1e-12)
            self.assertEqual(len(result["reference_dipoles"]),3)


if __name__=="__main__":
    unittest.main()
