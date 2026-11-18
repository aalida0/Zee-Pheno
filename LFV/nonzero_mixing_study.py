"""Nonzero-mixing Zee LFV truncation sensitivity; NOT an exclusion scan.

No inverse neutrino fit, phenomenological acceptance, or LFV-driven tuning.
A caller explicitly provides all lepton Yukawas and a scalar inversion scheme.
The reconstructed restricted Zee scalar potential changes deterministically
with masses and angles; it is not independently tuned to rescue any point.

The retained photon subset is coherent. Missing Zee W/Z gauge amplitudes
remain unbounded at every point, even if photon/one-loop ratios are small.
"""
import argparse
import copy
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

import numpy as np

from LFV.constraints import evaluate, matrix
from LFV.gauge_vertices import charged_currents
from Model.model import scalar_parameters


GRID_NAMES=("tan_beta","gamma_rad","mH2_GeV","F_multiplier")
INVERSION_BASE=("lambda_c","K1","K2","K12")
SOFT_MODES=("m12sq_GeV2","soft_mass_sq_GeV2")
POTENTIAL_NAMES=("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
                 "K1","K2","K12","Mu4")
DEFAULT_GRID={
    "tan_beta":[2.0,8.0],
    "gamma_rad":[0.15,0.35,0.55],
    "mH2_GeV":[900.0,1200.0],
    "F_multiplier":[1.0,100.0],
}
Z_MASS_GEV=91.1876


def _validate_scheme(config):
    if not isinstance(config,dict):
        raise ValueError("Explicit scalar inversion controls required")
    modes=set(config)&set(SOFT_MODES)
    optional={"alignment_delta_rad"} if "alignment_delta_rad" in config else set()
    if len(modes)!=1 or set(config)!=set(INVERSION_BASE)|modes|optional:
        raise ValueError("Provide exactly one of m12sq_GeV2 or soft_mass_sq_GeV2")
    if any(type(value) not in (float,int) or not math.isfinite(value)
           for value in config.values()):
        raise ValueError("Scalar inversion controls must be finite")
    if "soft_mass_sq_GeV2" in config and config["soft_mass_sq_GeV2"]<=0:
        raise ValueError("Soft scalar squared mass must be positive")
    if "alignment_delta_rad" in config and abs(config["alignment_delta_rad"])>.35:
        raise ValueError("Alignment control outside declared domain")
    return config


def grid_points(point, inversion, grid):
    """Declare and exhaust all nodes before evaluating any LFV observable."""
    if not isinstance(grid,dict) or set(grid)!=set(GRID_NAMES):
        raise ValueError("Exactly four grid axes required")
    for name,values in grid.items():
        if (not isinstance(values,list) or not 1<=len(values)<=12 or
            any(type(v) not in (float,int) or not math.isfinite(v)
                for v in values) or len(set(values))!=len(values)):
            raise ValueError(f"Invalid grid axis {name}")
    inversion=_validate_scheme(inversion)
    if point.get("potential") is not None:
        raise ValueError("Template potential must be reconstructed, not reused")
    F=matrix(point["F"],"F")
    if not np.max(abs(F))>0:
        raise ValueError("Need a nonzero F texture to audit charged-singlet LFV")
    n=0
    for tanbeta,gamma,mass,scale in itertools.product(
        grid["tan_beta"],grid["gamma_rad"],
        grid["mH2_GeV"],grid["F_multiplier"]
    ):
        if not (0.1<tanbeta<60 and
                0.0<abs(gamma)<math.pi/2-0.01 and
                scale>0 and mass>0):
            raise ValueError("Grid includes invalid tanbeta, mixing, F or mass")
        p=copy.deepcopy(point)
        p["scalar"]["beta"]=math.atan(tanbeta)
        p["scalar"]["gamma"]=gamma
        p["scalar"]["mH2"]=mass
        scalar=p["scalar"]
        if "alignment_delta_rad" in inversion:
            scalar["alpha"]=scalar["beta"]-inversion["alignment_delta_rad"]
        if not scalar["mH2"]>scalar["mH1"]+Z_MASS_GEV:
            raise ValueError("Cascade H2->H1 Z closed in requested grid")
        m12sq=(inversion["m12sq_GeV2"] if "m12sq_GeV2" in inversion
               else inversion["soft_mass_sq_GeV2"]*
                    math.sin(scalar["beta"])*math.cos(scalar["beta"]))
        sp=scalar_parameters(
            **scalar,
            m12sq=m12sq,
            lc=inversion["lambda_c"],
            k1=inversion["K1"],k2=inversion["K2"],k12=inversion["K12"],
        )
        p["potential"]={name:sp[name] for name in POTENTIAL_NAMES}
        p["F"]={
            "real":(F.real*scale).tolist(),
            "imag":(F.imag*scale).tolist(),
        }
        n+=1
        p["point_id"]=f"mix_{n:03d}_tb{tanbeta:g}_g{gamma:+.3f}_M{mass:g}_F{scale:g}"
        yield p,{
            "tan_beta":tanbeta,
            "gamma_rad":gamma,
            "mH2_GeV":mass,
            "F_multiplier":scale,
            "mH1_GeV":scalar["mH1"],
            "charged_cascade_open":True,
            "Mu4_GeV":sp["Mu4"],
            "beta_minus_alpha_rad":scalar["beta"]-scalar["alpha"],
            "m12sq_GeV2":m12sq,
            "max_abs_scalar_coupling":max(abs(sp[k]) for k in (
                "Lambda1","Lambda2","Lambda3","Lambda4","Lambda5","LambdaC",
                "K1","K2","K12")),
            "simple_quartic_magnitude_below_4pi":all(
                abs(sp[k])<4*math.pi for k in (
                "Lambda1","Lambda2","Lambda3","Lambda4","Lambda5","LambdaC",
                "K1","K2","K12")),
        }


def _coefficient(row,name):
    return complex(row[name]["real"],row[name]["imag"])


def evaluate_grid(template,inversion,grid,sw2):
    """Selected-topology computation for genuine gamma!=0 and F!=0 points.

    Returns diagnostics only; do NOT use any subset BR for acceptance.
    """
    if not (0<sw2<1):
        raise ValueError("Finite weak mixing parameter in (0,1) required")
    # Reject at input level rather than skipping failed points. This is
    # important for an unbiased grid and reproducible failure accounting.
    points=list(grid_points(template,inversion,grid))
    results=[]
    for p,axes in points:
        source=evaluate(p)
        if source["full_lfv_gate"]!="UNKNOWN":
            raise AssertionError("Incomplete Zee LFV became a physical gate")
        current=charged_currents(p["scalar"]["gamma"],sw2)
        offdiag=abs(float(current["Z_over_e"][0,1]))
        for name,process in source["processes"].items():
            rows=process["included_components"]
            cf=np.array([
                sum(_coefficient(row,c) for key,row in rows.items()
                    if key in ("1loop_H1_F","1loop_H2_F"))
                for c in ("C_L","C_R")
            ],complex)
            total=np.array([
                _coefficient(process,"C_L_GeV_minus2"),
                _coefficient(process,"C_R_GeV_minus2")
            ])
            topology=process["retained_topology_materiality"]
            if not math.isclose(float(np.linalg.norm(total)),
                               topology["retained_total_norm_GeV_minus2"],
                               rel_tol=1e-11,abs_tol=1e-24):
                raise AssertionError("Included topology ledger mismatch")
            nf=float(np.linalg.norm(cf))
            nt=float(np.linalg.norm(total))
            results.append({
                "point_id":p["point_id"],"channel":name,**axes,
                "Z12_over_e_abs":offdiag,
                "known_F_one_loop_norm_GeV_minus2":nf,
                "known_F_one_loop_over_retained_norm":nf/nt if nt else None,
                "retained_amplitude_norm_GeV_minus2":nt,
                "known_two_loop_over_one_loop":
                    topology["coherent_retained_two_loop_over_one_loop"],
                "retained_partial_BR":process["branching_ratio_included_topologies_only"],
                "retained_cancellation_ratio":process["cancellation_ratio"],
                "full_Zee_gauge_amplitude_upper_bound_GeV_minus2":None,
                "full_Zee_LFV_gate":"UNKNOWN",
                "neutrino_fit_gate":"NOT_EVALUATED",
                "scalar_Higgs_vacuum_gate":"NOT_EVALUATED",
                "2HDM_gauge_control":"OMITTED_NONPHOTON_TERMS_MATERIAL_IN_18PT_CONTROL",
                "method":"SOURCE_AUDITED_ONE_LOOP_PLUS_PHOTON_BZ_ONLY",
            })
    if len(results)!=3*len(points):
        raise AssertionError("Missing LFV channels in mixing sweep")
    summaries={}
    for channel in ("mu_to_e_gamma","tau_to_e_gamma","tau_to_mu_gamma"):
        channel_rows=[r for r in results if r["channel"]==channel]
        known=[r["known_two_loop_over_one_loop"] for r in channel_rows
               if r["known_two_loop_over_one_loop"] is not None]
        fracs=[r["known_F_one_loop_over_retained_norm"] for r in channel_rows
               if r["known_F_one_loop_over_retained_norm"] is not None]
        summaries[channel]={
            "n_points":len(channel_rows),
            "retained_two_loop_over_one_loop_median":
                float(np.median(known)) if known else None,
            "retained_two_loop_over_one_loop_max":
                max(known) if known else None,
            "fraction_retained_photon_two_loop_gt_one_loop":
                sum(value>1 for value in known)/len(known) if known else None,
            "fraction_known_F_one_loop_gt_10pct_retained":
                sum(value>.1 for value in fracs)/len(fracs) if fracs else None,
            "full_Zee_LFV_gate":"UNKNOWN",
        }
    output={
        "status":"NONZERO_MIXING_ZEE_TRUNCATION_DIAGNOSTIC_ONLY",
        "n_points":len(points),"n_channel_rows":len(results),
        "n_simple_quartic_magnitude_below_4pi":sum(
            axes["simple_quartic_magnitude_below_4pi"] for _,axes in points),
        "grid":grid,
        "scalar_inversion":inversion,
        "template_hash_sha256":hashlib.sha256(
            json.dumps(template,sort_keys=True,allow_nan=False).encode()
        ).hexdigest(),
        "fixed_D_texture":True,
        "fixed_F_texture_up_to_declared_multiplier":True,
        "physical_full_Zee_LFV_allowed_points":None,
        "type_ii_2HDM_omitted_gauge_bound":"NOT_TRANSFERABLE_TO_MIXED_ZEE",
        "source_only_LFV_threshold":None,
        "model_constraints_applied":[],
        "scalar_quartic_diagnostic":"ABS_LAMBDA_LT_4PI_ONLY_NOT_UNITARITY_BFB_OR_VACUUM",
        "known_source_controls":[
            "2017 and 2023 Zee partial Barr-Zee conventions; assumption-dependent",
            "2025/26 corrected general 2HDM 18-point control: all tested omitted norms >10% retained",
            "Exact Zee offdiagonal Z charged-scalar current; no Ward-complete Zee 2-loop gauge amplitude",
        ],
        "summaries":summaries,
        "rows":results,
    }
    return output


def main():
    a=argparse.ArgumentParser(description=__doc__)
    for key in ("template","inversion","output-dir"):
        a.add_argument("--"+key,type=Path,required=True)
    a.add_argument("--grid",type=Path)
    a.add_argument("--sin2thetaW",type=float,default=0.23129)
    x=a.parse_args()
    t=json.loads(x.template.read_text())
    inv=json.loads(x.inversion.read_text())
    grid=json.loads(x.grid.read_text()) if x.grid else DEFAULT_GRID
    report=evaluate_grid(t,inv,grid,x.sin2thetaW)
    x.output_dir.mkdir(parents=True,exist_ok=False)
    (x.output_dir/"nonzero_mixing_diagnostics.json").write_text(
        json.dumps(report,indent=2,allow_nan=False)+"\n")
    with (x.output_dir/"nonzero_mixing_diagnostics.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(report["rows"][0]))
        writer.writeheader()
        writer.writerows(report["rows"])
    print(json.dumps({
        "status":report["status"],"n_points":report["n_points"],
        "n_channel_rows":report["n_channel_rows"],
        "n_simple_quartic_magnitude_below_4pi":report["n_simple_quartic_magnitude_below_4pi"],
        "summaries":report["summaries"],
        "sample":report["rows"][:1],
        "full_Zee_LFV_status":"UNKNOWN",
    },indent=2))


if __name__=="__main__":
    main()
