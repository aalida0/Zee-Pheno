"""Source-locked h125 production/partial-width reference validation.

Uses SPheno 4.0.5 native NO/IO inputs for Zee loop couplings, and the
official HiggsTools v1.3 HiggsPredictions SM reference (HXSWG
HDECAY/Prophecy4f-based).  No 2HDMC/ScannerS, no invented BSM HDECAY
amplitudes.  Varies ONLY the h125 13TeV ggF reference normalization;
H/A SPheno rates are held unchanged.  h125 BR method is crossed with
ggF policy to quantify scheme dependence in official HiggsSignals 3.1.

This is NOT a frequentist confidence-level calculation or theory sigma.
"""
import argparse
import json
import math
from pathlib import Path

from ZeeHiggs.spheno_audit import audit
from ZeeHiggs.diagnose_h125 import (
    CHANNELS, make_inputs, predictions_for, particle_report, number
)
from ZeeHiggs.higgstools_pilot import dataset_receipt

COUPLINGS={
    "bb":("fermion",5),
    "cc":("fermion",4),
    "tautau":("fermion",15),
    "mumu":("fermion",13),
    "WW":("boson","WW"),
    "ZZ":("boson","ZZ"),
    "gg":("boson","gg"),
    "gamgam":("boson","gamgam"),
}
# 'ss' tracks the native scalar Type-II down-type coupling but independent
# source-verified kappa_s has not been populated; exclude from coupling closure.


def expected_width_scaling(source, channel):
    """Leading kappa rescaling, NOT a new perturbative calculation."""
    group,index=COUPLINGS[channel]
    h=source["neutral_couplings"]["25"]
    if group=="boson":
        c=number(h["boson"][index],channel)
        return c*c
    components=h["fermion_re_im"][str(index)]
    c_s,c_p=(number(x,channel) for x in components)
    return c_s*c_s+c_p*c_p


def validate_reference(order, signals):
    import Higgs.predictions as HP
    import Higgs.tools.Input as hinput

    source=audit(order)
    data,ids,removed,prod=make_inputs(source,hinput)
    hid=ids[25]
    mh=source["masses_GeV"]["25"]
    sm=HP.Predictions()
    h=sm.addParticle(HP.BsmParticle("SM_h",charge="neutral",cp="even"))
    h.setMass(mh)
    HP.effectiveCouplingInput(h,HP.smLikeEffCouplings)
    smchi=number(signals(sm),"sm_chi2")
    smwidths=particle_report(sm,"SM_h")

    # "historical" keeps SPheno approximate ggF; "HXSWG" lets the
    # official HiggsPredictions SM reference rescale by Zee kappa_gg.
    # Other neutral and charged states are IDENTICAL in the two policies.
    bsm_xs_reference=dict(data)
    bsm_xs_reference.pop(f"cxn_{hid}_ggH_LHC13",None)
    central={
        "SPheno_ggF_with_effective_BR":(data,False),
        "HXSWG_ggF_with_effective_BR":(bsm_xs_reference,False),
        "SPheno_ggF_with_explicit_BR":(data,True),
        "HXSWG_ggF_with_explicit_BR":(bsm_xs_reference,True),
    }
    results={}
    for name,(values,use_explicit) in central.items():
        predictions=predictions_for(values,ids,hinput,explicit=use_explicit)
        hh=predictions.particle(hid)
        xs=number(hh.cxn("LHC13","ggH"),"ggF")
        width=particle_report(predictions,hid)
        chi2=number(signals(predictions),"chi2_"+name)
        results[name]={
            "chi2":chi2,
            "delta_chi2_vs_same_mass_SM":chi2-smchi,
            "higgs_width_GeV":width["width_GeV"],
            "ggF_LHC13_pb":xs,
            "gamma_gamma_partial_GeV":width["partial_widths_GeV"]["gamgam"],
            "model_input":"SPheno source-locked native "+ ("explicit BRs" if use_explicit else "effective BRs") +
                          " with "+("old SPheno ggF absolute" if "SPheno_ggF" in name else
                          "official SM ggF rescaled by native Zee kappa_gg"),
        }
    r_effective=predictions_for(bsm_xs_reference,ids,hinput,explicit=False)
    r_explicit=predictions_for(bsm_xs_reference,ids,hinput,explicit=True)
    eff=particle_report(r_effective,hid)
    explicit=particle_report(r_explicit,hid)
    ref_sm_gg=number(sm.particle("SM_h").cxn("LHC13","ggH"),"SM ggF")
    kap_gg=expected_width_scaling(source,"gg")
    expected_reference_gg=ref_sm_gg*kap_gg
    actual_reference_gg=results["HXSWG_ggF_with_effective_BR"]["ggF_LHC13_pb"]
    if not math.isclose(expected_reference_gg,actual_reference_gg,rel_tol=3e-8):
        raise ValueError("HXSWG ggF does not reproduce SM*Zee kappa_gg^2")
    actual_old=results["SPheno_ggF_with_effective_BR"]["ggF_LHC13_pb"]
    if not math.isclose(actual_old,prod[f"cxn_{hid}_ggH_LHC13"],rel_tol=1e-10):
        raise ValueError("Native SPheno ggF override is not preserved")

    table={}
    for channel in CHANNELS:
        gamma_sm=smwidths["partial_widths_GeV"][channel]
        gamma_eff=eff["partial_widths_GeV"][channel]
        gamma_native=explicit["partial_widths_GeV"][channel]
        native_ratio=(gamma_native/gamma_eff) if gamma_eff>0 else None
        entry={
            "SM_partial_width_GeV":gamma_sm,
            "HXSWG_normalized_Zee_partial_GeV":gamma_eff,
            "native_SPheno_partial_GeV":gamma_native,
            "native_div_reference":native_ratio,
            "native_minus_reference_percent":100*(native_ratio-1) if native_ratio is not None else None,
            "ratio_is_calibrated_uncertainty":False,
        }
        if channel in COUPLINGS:
            ksquared=expected_width_scaling(source,channel)
            predicted=gamma_sm*ksquared
            entry["kappa_squared"]=ksquared
            entry["SM_times_kappa2_GeV"]=predicted
            entry["HXSWG_closure"]=(gamma_eff/predicted) if predicted else None
            # For weak-vector final states h->VV* this is within SM-like
            # kappa_W/V limit; for complex Yukawa phases the scalar+pseudo
            # norm is only a tiny-phase diagnostic.
            if not math.isclose(gamma_eff,predicted,rel_tol=1.5e-5):
                raise ValueError("Reference partial width fails kappa closure: "+channel)
        table[channel]=entry

    return {
        "point_id":source["point_id"],
        "status":"HXSWG_REF_AND_SPHENO_WIDTH_COMPARISON_EXECUTED",
        "source":source["source"],"mh_GeV":mh,
        "SM_reference":{"chi2":smchi,"total_width_GeV":smwidths["width_GeV"],
                        "ggF_LHC13_pb":ref_sm_gg},
        "production":{
            "SPheno_approx_ggF_LHC13_pb":actual_old,
            "HXSWG_normalized_Zee_ggF_LHC13_pb":actual_reference_gg,
            "SPheno_effective_SM_ggF_LHC13_pb":actual_old/kap_gg,
            "HXSWG_SM_ggF_LHC13_pb":ref_sm_gg,
            "kappa_g_squared":kap_gg,
            "normalization_ratio_new_div_old":actual_reference_gg/actual_old,
            "normalization_mismatch_percent":100*(actual_reference_gg/actual_old-1),
            "SPheno_SusHi_version":"1.5.0 approximate 13 TeV rescaling",
            "policy":"HXSWG HiggsPredictions SM normalization times source-proven Zee loop coupling",
            "new_policy_is_NNLO_Zee_production":False,
        },
        "widths":{
            "SPheno_native_total_GeV":explicit["width_GeV"],
            "HXSWG_ref_rescaled_total_GeV":eff["width_GeV"],
            "ratio_native_div_reference":explicit["width_GeV"]/eff["width_GeV"],
            "channel_table":table,
            "reference":"Official HiggsTools v1.3 SMHiggsInterp based on LHCHXSWG HDECAY/Prophecy4f",
            "SM_and_BSM_EW_corrections_independently_2loop_validated":False,
        },
        "higgsSignals_chi2_four_corners":results,
        "primary_exploratory_policy":"HXSWG_ggF_with_effective_BR",
        "crosscheck_policy":"SPheno_ggF_with_explicit_BR",
        "sensitivity": {
            "ggF_normalization_effect_with_effective_BR":
                results["HXSWG_ggF_with_effective_BR"]["chi2"]-
                results["SPheno_ggF_with_effective_BR"]["chi2"],
            "BR_method_effect_with_HXSWG_ggF":
                results["HXSWG_ggF_with_explicit_BR"]["chi2"]-
                results["HXSWG_ggF_with_effective_BR"]["chi2"],
            "are_not_additive_or_gaussian_theory_errors":True
        },
        "uncalculated_Zgamma_keys_removed":removed,
        "physical_higgs":"UNKNOWN","physical_HiggsSignals_95CL":None,
        "globally_allowed":None,
    }


def evaluate(hsdataset):
    import Higgs.signals as HS
    data=dataset_receipt(hsdataset)
    signals=HS.Signals(str(Path(hsdataset).resolve()))
    return {
      "status":"ACTUAL_HXSWG_SPHENO_RATE_VALIDATION_ONLY",
      "reference_data":data,"backend":"HiggsTools v1.3 HiggsSignals v3.1",
      "points":[validate_reference(order,signals) for order in ("NO","IO")],
      "physical_higgs":"UNKNOWN","globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hsdataset",required=True,type=Path)
    p.add_argument("--output",required=True,type=Path)
    a=p.parse_args()
    result=evaluate(a.hsdataset)
    with a.output.open("x") as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write("\n")
    for p in result["points"]:
        print(p["point_id"], "SMggF",p["production"]["HXSWG_SM_ggF_LHC13_pb"],
              "nativeggF",p["production"]["SPheno_approx_ggF_LHC13_pb"],
              "HXSWGggF",p["production"]["HXSWG_normalized_Zee_ggF_LHC13_pb"],
              "chi2ref",p["higgsSignals_chi2_four_corners"]["HXSWG_ggF_with_effective_BR"]["chi2"],
              "phys=UNKNOWN")


if __name__ == "__main__":
    main()
