"""Audited SPheno source controls, malformed input falsifiers, HiggsTools API mock."""
import copy
import json
import math
import sys
import tempfile
import types
import unittest
from pathlib import Path
from unittest.mock import patch

from ZeeHiggs.spheno_audit import NATIVE, audit, read_slha
from ZeeHiggs.higgstools_pilot import infer_ids, clean_missing_zero_couplings, run_one


class SourceAuditTests(unittest.TestCase):
    def test_real_native_NO_IO_audit(self):
        for order in ("NO","IO"):
            with self.subTest(order=order):
                d=audit(order)
                self.assertEqual(d["status"],"SOURCE_ATTESTED_NATIVE_SPHENO_HIGGS_INPUTS_AUDITED")
                self.assertEqual(d["physical_higgs"],"UNKNOWN")
                self.assertFalse(d["charged_production_accepted"])
                self.assertFalse(d["widths_precision_accepted"])
                self.assertEqual(d["neutral_pdgs"],[25,35,36])
                self.assertEqual(d["charged_pdgs"],[37,9000037])
                self.assertAlmostEqual(d["masses_GeV"]["25"],125.25,delta=.002)
                self.assertAlmostEqual(d["masses_GeV"]["37"],666.305,delta=.01)
                self.assertAlmostEqual(d["masses_GeV"]["9000037"],1184.55,delta=.02)
                self.assertAlmostEqual(d["total_widths_GeV"]["25"],.0034521311,delta=5e-7)
                self.assertAlmostEqual(d["charged_cascade_br"],.2781925,delta=1e-5)
                self.assertAlmostEqual(d["neutral_couplings"]["25"]["boson"]["WW"],
                                       .999735666,delta=1e-8)
                self.assertAlmostEqual(d["neutral_couplings"]["25"]["fermion_re_im"]["5"][0],
                                       .891842454,delta=1e-7)
                self.assertAlmostEqual(d["neutral_ggF_LHC13_pb"]["25"],48.3793449,delta=.01)
                for brsum in d["branching_sum"].values():
                    self.assertAlmostEqual(brsum,1.,delta=.001)
                self.assertIn("NOT CALCULATED","; ".join(d["input_limitations"]))

    def test_falsify_source_binary_or_card(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            src=NATIVE/"s00017_NO_benchmark"
            out=root/"s00017_NO_benchmark"
            out.mkdir()
            for name in ("input.in","spectrum.slha","native.json"):
                (out/name).write_bytes((src/name).read_bytes())
            v=json.loads((out/"native.json").read_text())
            v["production_binary_sha256"]="invalid"
            (out/"native.json").write_text(json.dumps(v))
            with self.assertRaisesRegex(ValueError,"Unattested"):
                audit("NO",root)

    def test_falsify_missing_and_nonclosing_decays(self):
        original=(NATIVE/"s00017_NO_benchmark/spectrum.slha").read_text()
        with self.assertRaisesRegex(ValueError,"Duplicate DECAY"):
            read_slha(original+"\nDECAY 25 1.0\n")
        with self.assertRaisesRegex(ValueError,"Nonfinite"):
            read_slha(original.replace("DECAY        25     3.45213110E-03",
                                       "DECAY        25 nan"))
        with self.assertRaises(ValueError):
            read_slha(original.replace("DECAY        25     3.45213110E-03",
                                       "DECAY        25 -1.0"))

    def test_no_extra_mass_ids_inferred_by_approx(self):
        report=audit("NO")
        data={f"m_h{i}":report["masses_GeV"][str(p)]
              for i,p in enumerate([25,35,36],1)}
        data.update({f"m_Hp{i}":report["masses_GeV"][str(p)]
                     for i,p in enumerate([37,9000037],1)})
        ids=infer_ids(data,report)
        self.assertEqual(ids[25],"h1")
        self.assertEqual(ids[9000037],"Hp2")
        with self.assertRaisesRegex(ValueError,"ambiguous"):
            infer_ids({**data,"m_other":data["m_h1"]},report)
        with self.assertRaisesRegex(ValueError,"unmapped"):
            infer_ids({k:v for k,v in data.items() if k!="m_h1"},report)

    def test_clean_not_calculated_zgamma_but_not_gammagamma(self):
        data={"effc_h1_Zgam":0.,"effc_h2_gamZ":0.,"effc_h1_gamgam":1.039,
              "m_h1":125.25}
        c,missing=clean_missing_zero_couplings(data)
        self.assertEqual(len(missing),2)
        self.assertNotIn("effc_h1_Zgam",c)
        self.assertIn("effc_h1_gamgam",c)

    def test_mocked_actual_higgstools_entrypoint_contract(self):
        report=audit("NO")
        data={}
        for pid,i in [(25,"h1"),(35,"h2"),(36,"h3"),(37,"Hp1"),(9000037,"Hp2")]:
            data["m_"+i]=report["masses_GeV"][str(pid)]
        data["effc_h1_Zgam"]=0.
        calls=[]
        class FakeInput:
            def readHB5SLHA(self,path,*,neutralPDGs,chargedPDGs):
                selfn=list(neutralPDGs);selfc=list(chargedPDGs)
                if selfn != [25,35,36] or selfc != [37,9000037]:
                    raise AssertionError("wrong PDGs")
                return dict(data)
            def predictionsFromDict(self,d,neutrals,charged,doubly,*,useExplicitBr):
                self_outer={"dict":d,"neutral":neutrals,"charged":charged,
                            "doubly":doubly,"explicit":useExplicitBr}
                calls.append(self_outer)
                return self_outer
        fake=FakeInput()
        m1=types.ModuleType("Higgs")
        m2=types.ModuleType("Higgs.tools")
        m3=types.ModuleType("Higgs.tools.Input")
        m3.readHB5SLHA=fake.readHB5SLHA
        m3.predictionsFromDict=fake.predictionsFromDict
        m1.tools=m2
        m2.Input=m3
        class HB:
            allowed=True
            selectedLimits={}
        bounds=lambda pred:HB()
        def signals(pred): return 123.4 if pred["explicit"] is True else 120.3
        with patch.dict(sys.modules,{"Higgs":m1,"Higgs.tools":m2,"Higgs.tools.Input":m3}), \
             patch("ZeeHiggs.higgstools_pilot.sm_higgs_reference",return_value=100.0):
            result=run_one("NO",bounds,signals)
        self.assertEqual(result["status"],"OFFICIAL_HIGGSTOOLS_DIAGNOSTIC_ONLY")
        self.assertEqual(result["physical_higgs"],"UNKNOWN")
        self.assertIsNone(result["globally_allowed"])
        self.assertEqual(len(calls),2)
        self.assertEqual(calls[0]["neutral"],["h1","h2","h3"])
        self.assertEqual(calls[0]["charged"],["Hp1","Hp2"])
        self.assertEqual(calls[1]["explicit"]["h1"],False)
        self.assertEqual(calls[1]["explicit"]["Hp2"],True)
        self.assertEqual(calls[0]["dict"]["cxn_h2_ggH_LHC13"],
                         report["neutral_ggF_LHC13_pb"]["35"])
        self.assertNotIn("effc_h1_Zgam",calls[0]["dict"])
        self.assertAlmostEqual(
            result["branches"]["effective_coupling_h125_BR"]["chi2_HiggsSignals"],120.3)
        self.assertEqual(result["higgsSignals_SM_chi2_same_h125_mass"],100.0)
        self.assertAlmostEqual(
            result["branches"]["effective_coupling_h125_BR"]["delta_chi2_vs_same_mass_SM"],20.3)


if __name__=="__main__":
    unittest.main()
