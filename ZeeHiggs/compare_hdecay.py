"""Compare actual HDECAY SM, HiggsTools HXSWG SM, and SPheno Zee h125 widths.

Physics scope: HDECAY is *SM only*. Zee is approximated by channelwise
source-attested effective couplings to test discrepancies; no new HDECAY
Zee-model amplitude or 2HDM/ScannerS has been introduced.

HDECAY h->WW*/ZZ* and the HXSWG Prophecy4f-based SM reference do not
necessarily have identical electroweak/off-shell schemes; differences are
diagnostics, not theoretical nuisance parameters.
"""
import argparse
import json
import math
from pathlib import Path

CHANNELS=("bb","cc","tautau","mumu","WW","ZZ","gg","gamgam")


def combine(hdecay,hxswg):
    if hdecay.get("status")!="INDEPENDENT_HDECAY_SM_REFERENCE_ONLY":
        raise ValueError("Unvalidated HDECAY source")
    if hxswg.get("status")!="ACTUAL_HXSWG_SPHENO_RATE_VALIDATION_ONLY":
        raise ValueError("Unvalidated HXSWG validation inputs")
    rows={}
    xpts={p["point_id"].rsplit("_",1)[-1]:p for p in hxswg["points"]}
    if set(xpts)!={"NO","IO"} or set(hdecay["results"])!={"NO","IO"}:
        raise ValueError("Unexpected NO/IO reference control set")
    for order in ("NO","IO"):
        p=xpts[order]; h=hdecay["results"][order]
        if abs(p["mh_GeV"]-h["native_Zee_mass_GeV"])>1e-6:
            raise ValueError("Incompatible scalar mass")
        lines={}
        for ch in CHANNELS:
            hd_sm=h["partial_width_GeV"][ch]
            k=p["widths"]["channel_table"][ch].get("kappa_squared")
            if not isinstance(k,(int,float)) or not math.isfinite(k):
                raise ValueError("Missing real source-matched effective coupling for "+ch)
            he_sm=p["widths"]["channel_table"][ch]["SM_partial_width_GeV"]
            he_zee=p["widths"]["channel_table"][ch]["HXSWG_normalized_Zee_partial_GeV"]
            sp=p["widths"]["channel_table"][ch]["native_SPheno_partial_GeV"]
            if min(hd_sm,he_sm,he_zee,sp)<=0:
                raise ValueError("Zero/negative source width: "+ch)
            hdecay_zee_proxy=hd_sm*k
            lines[ch]={
                "HDECAY_SM_partial_width_GeV":hd_sm,
                "HXSWG_SM_partial_width_GeV":he_sm,
                "HDECAY_SM_over_HXSWG_SM":hd_sm/he_sm,
                "native_Zee_partial_GeV":sp,
                "HXSWG_kappa_rescaled_Zee_partial_GeV":he_zee,
                "HDECAY_SM_kappa_rescaled_Zee_proxy_GeV":hdecay_zee_proxy,
                "native_Zee_div_HDECAY_proxy":sp/hdecay_zee_proxy,
                "native_Zee_div_HXSWG_reference":sp/he_zee,
                "HDECAY_proxy_is_an_independent_full_Zee_calculation":False,
            }
        rows[order]={
            "mass_GeV":p["mh_GeV"],
            "HDECAY_SM_total_width_GeV":h["width_GeV"],
            "HXSWG_SM_total_width_GeV":p["SM_reference"]["total_width_GeV"],
            "native_Zee_total_width_GeV":p["widths"]["SPheno_native_total_GeV"],
            "HXSWG_kappa_rescaled_Zee_total_width_GeV":p["widths"]["HXSWG_ref_rescaled_total_GeV"],
            "channels":lines,
        }
    return {
        "status":"INDEPENDENT_SM_HDECAY_VS_HXSWG_SPHENO_VALIDATION_DIAGNOSTIC_ONLY",
        "HDECAY_source_commit":hdecay["HDECAY_git_commit"],
        "HDECAY_source_sha256":hdecay["HDECAY_source_sha256"],
        "dataset":hxswg["reference_data"],
        "results":rows,
        "theory_sigma":None,
        "physical_higgs":"UNKNOWN",
        "globally_allowed":None,
        "warnings":[
            "HDECAY calculates only SM; kappa-rescaled HDECAY rates are Zee proxies, not source-matched full BSM results",
            "HDECAY stand-alone WW/ZZ off-shell implementation can differ from HXSWG Prophecy4f 4f reference",
            "SM parameter inputs (mb, mc, alpha_s, MW, electroweak options) differ and were not retuned",
            "A channel difference is NOT a calibrated theory uncertainty",
            "BSM contributions to Zgamma are not calculated in native SPheno"
        ],
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hdecay",type=Path,required=True)
    p.add_argument("--hxswg",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    data=combine(json.loads(a.hdecay.read_text()),
                 json.loads(a.hxswg.read_text()))
    with a.output.open("x") as f:
        json.dump(data,f,indent=2,allow_nan=False);f.write("\n")
    for o,row in data["results"].items():
        print("HDECAY_VS_HXSWG",o,"SM_total_width",row["HDECAY_SM_total_width_GeV"],
              row["HXSWG_SM_total_width_GeV"],
              "native_ratio_bb",row["channels"]["bb"]["native_Zee_div_HDECAY_proxy"],
              "native_ratio_ZZ",row["channels"]["ZZ"]["native_Zee_div_HDECAY_proxy"],
              "phys=UNKNOWN")


if __name__=="__main__":
    main()
