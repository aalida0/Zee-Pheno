"""Source-bound SPheno Zee scalar theory and correlated PDG 2025 oblique fit.

PDG Review of Particle Physics 2025, "Electroweak Model and Constraints on
New Physics", Table 10.8, independent S,T,U column:
  (S,T,U)=(0.021 +/-0.096, 0.040 +/-0.120, 0.008 +/-0.092)
  correlations rho_ST=+0.91, rho_SU=-0.62, rho_TU=-0.83.
Use the THREE-PARAMETER covariance with U FLOATING, not the separate
PDG U=0 fit or independent marginal 2 sigma cuts.
95% conditional 3D ellipsoid: chi2_3=7.814727903251179.
The actual SPheno output's S,T,U are from its one-loop BSM calculation,
and the audit checks label/field bindings but does NOT independently
rederive the Zee two-charged-scalar self-energies. The full-model EWPO
confidence gate is therefore NOT yet claimed as source-matched validated.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from Scan.zee_scalar_theory import NAMES, evaluate as scalar_evaluate
from ZeeHiggs.spheno_audit import audit, read_slha, sha256

PDG2025={
 "source":"PDG 2025, Review 'Electroweak Model and Constraints on New Physics', Table 10.8",
 "url":"https://pdg.lbl.gov/2025/reviews/rpp2025-rev-standard-model.pdf",
 "reference":"U FLOATING (right half) three-parameter STU column",
 "mean_STU":[.021,.040,.008],
 "sigma_STU":[.096,.120,.092],
 "correlations_STU":[[1.,.91,-.62],[.91,1.,-.83],[-.62,-.83,1.]],
 "chi2_3_95":7.814727903251179,
 "fit_reference_SM_spectrum":"PDG ref SM as defined in review; source-matched one-loop SPheno subtraction still to audit",
}
MINPAR_MAP={i:k for i,k in enumerate(NAMES,1)}
SPHENO_LOWENERGY_KEYS={1:"T",2:"S",3:"U"}


def numeric(v):
    if type(v) not in (int,float) or not math.isfinite(v):
        raise ValueError("Nonfinite numeric observable")
    return float(v)


def block_1d(rows,name):
    out={}
    for fields,raw in rows:
        if len(fields)<2:
            raise ValueError("Bad SLHA "+name+" row")
        try:
            key=int(fields[0])
            value=float(fields[1].replace("D","E").replace("d","e"))
        except (TypeError,ValueError,IndexError) as exc:
            raise ValueError("Malformed SLHA "+name) from exc
        if not math.isfinite(value) or key in out:
            raise ValueError("Nonfinite or duplicate "+name+" row")
        out[key]=(value,raw.lower())
    return out


def read_attested_theory(order):
    # This verifies native spectrum and input SHA against archived
    # original SPheno binary before using ANY scalar/EWPO value.
    report=audit(order)
    text=Path(report["spectrum_path"]).read_text()
    blocks,_=read_slha(text)
    for name in ("MINPAR","ZEEQUARTIC","ZEEPOT","SPHENOLOWENERGY","SPHENO"):
        if name not in blocks:
            raise ValueError("Native SPheno missing required block "+name)
    native_min=block_1d(blocks["MINPAR"],"MINPAR")
    output=block_1d(blocks["ZEEQUARTIC"],"ZEEQUARTIC")
    potential=block_1d(blocks["ZEEPOT"],"ZEEPOT")
    low=block_1d(blocks["SPHENOLOWENERGY"],"SPhenoLowEnergy")
    settings=block_1d(blocks["SPHENO"],"SPheno")

    quartic={}
    for index,name in MINPAR_MAP.items():
        if index not in native_min or index not in output:
            raise ValueError("Missing native/low-scale Zee quartic "+name)
        value=native_min[index][0]
        produced=output[index][0]
        if not math.isclose(value,produced,rel_tol=2e-7,abs_tol=3e-9):
            raise ValueError("Native scalar quartic MINPAR/ZEEQUARTIC mismatch "+name)
        quartic[name]=produced

    # Check the native model's running potential components, including
    # cubic mu4.  These are not new independent inputs.
    for inkey,outkey in ((11,3),(12,4),(13,5)):
        if not math.isclose(native_min[inkey][0],potential[outkey][0],
                            rel_tol=2e-7,abs_tol=3e-5):
            raise ValueError("Native SPheno MINPAR vs ZEEPOT mismatch")
    if settings[57][0]!=1:
        raise ValueError("SPheno did not enable its low-energy observables")

    stu={}
    for index,key in SPHENO_LOWENERGY_KEYS.items():
        value,label=low[index]
        if key.lower()+"-parameter" not in label or "1-loop bsm" not in label:
            raise ValueError("SPheno STU observable label/key mismatch")
        stu[key]=numeric(value)
    return {
        "point_id":"s00017_"+order,
        "source":report["source"],
        "quartics":quartic,
        "quartics_at_Q_GeV":160.,
        "MINPAR_match_running_ZEEQUARTIC":"PASS",
        "MINPAR_match_ZEEPOT":"PASS",
        "spheno_low_energy_enabled":True,
        "STU":{k:stu[k] for k in ("S","T","U")},
        "SPheno_delta_rho_separate":low.get(39,(None,))[0],
        "spheno_output_comments_attest_one_loop_BSM":True,
        "STU_independently_recomputed_with_two_mixed_charged_scalars":False,
    }


def correlated_pdg2025(stu):
    if not isinstance(stu,dict) or set(stu)!={"S","T","U"}:
        raise ValueError("Exactly all three source STU values required")
    v=np.array([numeric(stu[k]) for k in ("S","T","U")])
    center=np.array(PDG2025["mean_STU"])
    sigma=np.array(PDG2025["sigma_STU"])
    rho=np.array(PDG2025["correlations_STU"])
    cov=np.outer(sigma,sigma)*rho
    eigenvalues=np.linalg.eigvalsh(cov)
    if np.min(eigenvalues)<=0 or np.min(np.linalg.eigvalsh(rho))<=0:
        raise ValueError("PDG 2025 STU covariance is not positive definite")
    d=v-center
    chi=float(d@np.linalg.solve(cov,d))
    if chi<0 or not math.isfinite(chi):
        raise ValueError("Bad correlated oblique statistic")
    # This is a fixed three-dimensional Gaussian ellipse, not a
    # model best-fit or joint-Higgs/flavor statistical analysis.
    flag=("PDG2025_STU_CORRELATED_3D_SCREEN_PASS" if
          chi<=PDG2025["chi2_3_95"] else
          "PDG2025_STU_CORRELATED_3D_SCREEN_FAIL")
    return {
        "status":flag,
        "point_STU":dict(stu),
        "PDG2025_3D_correlated_chi2":chi,
        "PDG2025_95pct_threshold_3D":PDG2025["chi2_3_95"],
        "ellipsoid_margin":PDG2025["chi2_3_95"]-chi,
        "PDG2025_U_floats":True,
        "PDG2025_fit_reference":PDG2025,
        "physics_assumptions_validated":{
            "source_SPheno_one_loop_output_binding":True,
            "full_Zee_specific_oblique_selfenergy_recalculation":False,
            "SM_reference_point_and_two_charged_sector_convention_reconciled":False,
            "additional_nonoblique_EW_corrections_included":False,
        },
        "physical_EWPO":"UNKNOWN",
        "globally_allowed":None,
    }


def evaluate(order):
    point=read_attested_theory(order)
    scalar=scalar_evaluate(point["quartics"])
    ew=correlated_pdg2025(point["STU"])
    # Delayed import avoids a cycle: the general Grimus reader reuses
    # this module's strictly checked native SLHA low-energy parser.
    from Scan.zee_oblique_STU_independent import evaluate as independent_STU
    source_comparison=independent_STU(order)
    return {
        "status":"SOURCE_ATTESTED_ZEE_TREE_AND_EWPO_2025_PILOT",
        "point_id":point["point_id"],
        "source":point["source"],
        "quartic_source":{k:point[k] for k in ("quartics","quartics_at_Q_GeV",
          "MINPAR_match_running_ZEEQUARTIC","MINPAR_match_ZEEPOT")},
        "scalar_theory":scalar,
        "EWPO":ew,
        "SPheno_low_energy":{k:point[k] for k in ("STU",
          "SPheno_delta_rho_separate",
          "spheno_output_comments_attest_one_loop_BSM",
          "STU_independently_recomputed_with_two_mixed_charged_scalars")},
        "BFB_quartic_asymptotic":scalar["BFB"]["full_Zee_quartic_BFB"],
        "tree_high_energy_Zee_unitarity":scalar["unitarity"]["tree_high_energy_full_Zee_unitarity"],
        "EWPO_correlated_3D_screen":ew["status"],
        "EWPO_scalar_full_STU_independent_Grimus2008":source_comparison,
        "EWPO_scalar_one_loop_pragmatically_reproduced":
            source_comparison["independent_one_loop_scalar_STU_pragmatic_consistency"],
        "EWPO_full_theory_validated":False,
        "global_vacuum":"UNKNOWN",
        "globally_allowed":None,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    rows={k:evaluate(k) for k in ("NO","IO")}
    result={
        "status":"SOURCE_ATTESTED_SPHENO_ZEE_BFB_UNITARITY_EWPO_PILOT",
        "points":rows,
        "PDG_fit":PDG2025,
        "physical_full_Zee_EWPO":"UNKNOWN",
        "global_vacuum":"UNKNOWN",
        "globally_allowed":None,
    }
    with a.output.open("x") as file:
        json.dump(result,file,indent=2,allow_nan=False);file.write("\n")
    for key,row in rows.items():
        print(row["point_id"],"BFB",row["BFB_quartic_asymptotic"],
              "unitarity",row["tree_high_energy_Zee_unitarity"],
              "STU_chi2",row["EWPO"]["PDG2025_3D_correlated_chi2"],
              "EWPO",row["EWPO_correlated_3D_screen"],
              "global_vacuum=UNKNOWN")


if __name__=="__main__":
    main()
