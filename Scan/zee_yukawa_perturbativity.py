"""Native Zee running Yukawa and quartic perturbativity at specified Q.

All entries are *source-output* matrices, not NuFIT fitted template
Yukawas; maximal singular value is basis independent. Complex
antisymmetric FZEE is checked, as is source RG order and Q.

RGE flow is a separate physical condition: 2-loop RGE level in
SPhenoInput 38 DOES NOT mean one spectrum at Q=160 GeV proves
perturbativity at a requested 10 TeV cutoff. To make this final,
supply actual native running solutions with source-matched Q labels
through the scan's RGE trajectory manifest; do not use arbitrary
hand-built 2HDM beta functions that omit the singlet lepton Yukawas.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ZeeHiggs.spheno_audit import audit,read_slha
from Scan.zee_ewpo import block_1d
from Scan.zee_scalar_theory import NAMES,unitarity,bfb

YUKAWAS=("Yu","Yd","Ye","YEE2","FZEE")
MAX_2BODY_YUKAWA=math.sqrt(4*math.pi)
LIMIT_SCALAR=4*math.pi


def y_matrix(blocks,name):
    imag="IM"+name
    if name not in blocks or imag not in blocks:
        raise ValueError("Missing complex Yukawa block "+name)
    def parse(rows,label):
        mat=np.full((3,3),np.nan)
        for f,_ in rows:
            if len(f)!=3:raise ValueError("Non-3-index Yukawa row "+label)
            i,j=int(f[0])-1,int(f[1])-1
            if not(0<=i<3 and 0<=j<3) or np.isfinite(mat[i,j]):
                raise ValueError("Invalid/duplicate Yukawa entry "+label)
            mat[i,j]=float(f[2].replace("D","E"))
        if not np.isfinite(mat).all():
            raise ValueError("Incomplete/nonfinite Yukawa "+label)
        return mat
    return parse(blocks[name],name)+1j*parse(blocks[imag],imag)


def analyze_matrices(blocks,*,cutoff_Q_GeV=None, native_Q_GeV=None):
    if native_Q_GeV is None or not 0<native_Q_GeV<1e20:
        raise ValueError("Missing actual source native Q")
    results={}
    max_sv=0
    for name in YUKAWAS:
        y=y_matrix(blocks,name)
        vals=np.linalg.svd(y,compute_uv=False)
        largest=float(max(vals))
        max_sv=max(max_sv,largest)
        result={"max_singular_value":largest,
                "frobenius_norm":float(np.linalg.norm(y)),
                "spectral_threshold_sqrt_4pi":MAX_2BODY_YUKAWA,
                "local_Q_status":"PASS" if largest<MAX_2BODY_YUKAWA else "FAIL"}
        if name=="FZEE":
            anti=float(np.linalg.norm(y+y.T))
            result["antisymmetry_residual"]=anti
            if anti>max(1e-8,1e-6*np.linalg.norm(y)):
                raise ValueError("Native FZEE violates antisymmetry")
        results[name]=result
    cutoff_Q_GeV=cutoff_Q_GeV or native_Q_GeV
    if cutoff_Q_GeV<native_Q_GeV:
        raise ValueError("Perturbativity cutoff cannot be below running input Q")
    return {
        "status":"LOCAL_RUNNING_YUKAWA_SINGULAR_VALUES_VALIDATED",
        "native_Q_GeV":native_Q_GeV,
        "Yukawa_scheme":"source_SPheno_running_at_native_Q",
        "singular_values":results,
        "max_Yukawa_singular_value":max_sv,
        "local_Yukawa_perturbativity":"PASS" if max_sv<MAX_2BODY_YUKAWA else "FAIL",
        "high_scale_cutoff_GeV":cutoff_Q_GeV,
        "RGE_to_cutoff_proven":cutoff_Q_GeV<=native_Q_GeV*(1+1e-10),
        "RGE_to_cutoff_status":
            "NOT_REQUIRED_LOCAL_SCALE_ONLY" if cutoff_Q_GeV<=native_Q_GeV*(1+1e-10)
            else "UNKNOWN_REQUIRES_SOURCE_RGE_TRAJECTORY",
    }


def read_native(order,cutoff_Q_GeV=10000):
    x=audit(order)
    blocks,_=read_slha(Path(x["spectrum_path"]).read_text())
    sp=block_1d(blocks["SPHENO"],"SPHENO")
    q=sp[33][0]
    if sp[38][0]!=2:
        raise ValueError("Unexpected SPheno RGE loop order")
    source=analyze_matrices(blocks,cutoff_Q_GeV=cutoff_Q_GeV,native_Q_GeV=q)
    p=block_1d(blocks["ZEEQUARTIC"],"ZEEQUARTIC")
    quartics={k:p[i][0] for i,k in enumerate(NAMES,1)}
    maxlambda=max(abs(v) for v in quartics.values())
    return {
        "status":"SOURCE_ATTESTED_RUNNING_PERTURBATIVITY_LOW_SCALE",
        "point_id":"s00017_"+order,"source":x["source"],
        "SPheno_2_loop_RGE_mode_enabled":True,
        "local_quartic_4pi_prior":"PASS" if maxlambda<LIMIT_SCALAR else "FAIL",
        "local_high_energy_scalar_unitarity":unitarity(quartics)["tree_high_energy_full_Zee_unitarity"],
        "local_full_quartic_BFB":bfb(quartics)["full_Zee_quartic_BFB"],
        "max_abs_native_quartic":maxlambda,
        "yukawas":source,
        "full_UV_perturbativity_to_10TeV":
            ("PASS_LOCAL_ONLY_NOT_UV" if source["RGE_to_cutoff_proven"] else
             "UNKNOWN_REQUIRES_SOURCE_RGE_TRAJECTORY"),
        "globally_allowed":None,
    }


def evaluate_trajectory(rows,cutoff_Q_GeV):
    """Validate source-evaluated multislice RGE trajectory (NOT interpolation).

    Each row from a separately authenticated spectrum with a real
    Q label and physical matrices; fail closed on missing intermediate
    values and ANY failed local check. Does NOT certify unobserved
    intervals between Q sample points.
    """
    if not rows or len(rows)<2:
        raise ValueError("Need >=2 native RGE scales")
    qs=[r["native_Q_GeV"] for r in rows]
    if qs!=sorted(set(qs)) or max(qs)<cutoff_Q_GeV*(1-1e-8):
        raise ValueError("Unsorted/insufficient native Q coverage")
    failure=any(r["local_Yukawa_perturbativity"]=="FAIL" for r in rows)
    return {
        "status":"NATIVE_RGE_SAMPLED_YUKAWA_PERTURBATIVITY",
        "native_Q_grid_GeV":qs,
        "highest_Q_GeV":max(qs),
        "requested_cutoff_GeV":cutoff_Q_GeV,
        "Yukawa_samples_all_perturbative":not failure,
        "continuous_RGE_trajectory_and_scalar_couplings_validated":False,
        "RGE_final_status":"FAIL_SAMPLED" if failure else "UNKNOWN_CONTINUOUS_FLOW_NOT_CERTIFIED",
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--cutoff-Gev",type=float,default=10000.)
    a=p.parse_args()
    res={"status":"NATIVE_ZEE_LOW_SCALE_YUKAWA_AND_SCALAR_CHECK",
         "points":{o:read_native(o,a.cutoff_Gev) for o in ("NO","IO")},
         "full_RGE_perturbativity":"UNKNOWN",
         "globally_allowed":None}
    with a.output.open("x") as f:
        json.dump(res,f,indent=2,allow_nan=False);f.write("\n")
    for k,v in res["points"].items():
        print(k,"max_Yukawa_singular",v["yukawas"]["max_Yukawa_singular_value"],
              "local_status",v["yukawas"]["local_Yukawa_perturbativity"],
              "RGE_10TeV",v["full_UV_perturbativity_to_10TeV"])


if __name__=="__main__":main()
