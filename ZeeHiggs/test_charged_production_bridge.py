"""Source-bound HiggsTools Hpmtb rate interface; no invented MadGraph run."""
import copy
import unittest

from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.charged_production_bridge import (
    validate_charged_rates,insert_charged_rates
)


def fixture(order="NO"):
    source=audit(order)
    rows=[]
    for pdg in (37,9000037):
        rows.append({
            "pdg":pdg,
            "mass_GeV":source["masses_GeV"][str(pdg)],
            "mode":"Hpmtb","collider":"LHC13","sqrt_s_TeV":13,
            "flavor_scheme":"5FS","perturbative_order":"LO",
            "charge_inclusive":True,
            "sigma_pb":0.012 if pdg==37 else 0.003,
            "integration_error_pb":0.00001,
            "generator":"MG5_aMC",
            "param_card_sha256":"a"*64,
            "run_card_sha256":"b"*64,
            "process_card_sha256":"c"*64,
            "generator_result_sha256":"d"*64
        })
    return {
        "schema_version":1,"point_id":source["point_id"],
        "source":dict(source["source"]),"rates":rows
    },source


class ChargedProductionReceiptTests(unittest.TestCase):
    def test_exact_point_two_eigenstates_and_preserve_existing_modes(self):
        item,source=fixture()
        validated=validate_charged_rates(item,source)
        self.assertTrue(validated["charged_13TeV_Hpmtb_input_present"])
        self.assertFalse(validated["all_energy_and_topology_production_complete"])
        self.assertIsNone(validated["full_HiggsBounds_model_exclusion"])
        base={"cxn_37_brtHpb_LHC13":0.0, "cxn_25_ggH_LHC13":48.0}
        ids={25:"25",35:"35",36:"36",37:"37",9000037:"9000037"}
        result=insert_charged_rates(base,ids,validated)
        self.assertEqual(result["cxn_37_Hpmtb_LHC13"],0.012)
        self.assertEqual(result["cxn_9000037_Hpmtb_LHC13"],0.003)
        self.assertEqual(result["cxn_37_brtHpb_LHC13"],0)
        self.assertNotIn("cxn_37_Hpmtb_LHC13",base)
        with self.assertRaisesRegex(ValueError,"overwrite"):
            insert_charged_rates(result,ids,validated)

    def test_refuse_cross_order_and_cross_spectrum(self):
        item,source=fixture("NO")
        _,io=fixture("IO")
        with self.assertRaisesRegex(ValueError,"point ID"):
            validate_charged_rates(item,io)
        item["source"]["spectrum_sha256"]="0"*64
        with self.assertRaisesRegex(ValueError,"source mismatch"):
            validate_charged_rates(item,source)

    def test_refuse_13p6_TeV_or_nonchargeinclusive_or_inclusive_4FS(self):
        for key,value in (
            ("sqrt_s_TeV",13.6),("charge_inclusive",False),
            ("flavor_scheme","4FS"),("mode","brtHpb"),
            ("perturbative_order","NNLO")
        ):
            item,source=fixture()
            item["rates"][0][key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,"Incompatible"):
                validate_charged_rates(item,source)

    def test_refuse_masses_duplicates_unverified_zero_and_missing(self):
        for change in ("mass","duplicate","zero","missing","nan","hash"):
            item,source=fixture()
            if change=="mass":item["rates"][0]["mass_GeV"]+=2
            if change=="duplicate":item["rates"][1]["pdg"]=37
            if change=="zero":item["rates"][0]["sigma_pb"]=0
            if change=="missing":item["rates"].pop()
            if change=="nan":item["rates"][0]["sigma_pb"]=float("nan")
            if change=="hash":item["rates"][0]["generator_result_sha256"]="bogus"
            with self.subTest(change=change),self.assertRaises(ValueError):
                validate_charged_rates(item,source)

    def test_explicit_zero_requires_provenance(self):
        item,source=fixture()
        item["rates"][0]["sigma_pb"]=0
        item["rates"][0]["integration_error_pb"]=0
        item["rates"][0]["explicit_zero_reason"]="Audited vanishing doublet coupling"
        result=validate_charged_rates(item,source)
        self.assertEqual(result["higgstools_cross_sections"]["cxn_37_Hpmtb_LHC13"],0)
        self.assertFalse(result["generator_amplitudes_and_raw_rates_independently_verified"])


if __name__=="__main__":
    unittest.main()
