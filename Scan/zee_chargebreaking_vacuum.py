"""Zee 2021 Eq.17 exact gauge-orbit reduced tree vacuum and witness search.

All complex SU(2) doublet and charged-singlet directions are captured:
x=|Phi1|^2=r1^2, y=|Phi2|^2=r2^2,
Phi1†Phi2=r1*r2*cos(theta)*exp(i phi),
|det(Phi1,Phi2)|=r1*r2*sin(theta),
|chi+|=rC. The cubic phase can always be chosen to minimize it:
-2*abs(mu4)*r1*r2*rC*sin(theta).

This is an EXACT minimization of phases, NOT a neutral-only slice.
A converged lower-potential field configuration is a constructive
charge-breaking/cp-breaking witness. Failure to find one after a
finite multistart is NOT a mathematical global-minimum proof.

MOST IMPORTANT: native SPheno's loop-corrected tadpole-derived
m1^2 and m2^2 need not satisfy the tree minimum equations.
We DO NOT compare to a claimed tree normal stationary point if they
do not. That requires a tree-level SPheno tadpole output or an
independent matching calculation on the same source card.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from ZeeHiggs.spheno_audit import audit,read_slha
from Scan.zee_ewpo import block_1d
from Scan.zee_scalar_theory import bfb,validate

VDIM=246.0


def normal_tree_mass_squared(q,vevs,m12sq):
    q=validate(q)
    v1,v2=(float(v) for v in vevs)
    if not(0<v1<1e4 and 0<v2<1e4) or not math.isfinite(m12sq):
        raise ValueError("Invalid electroweak vevs/m12sq")
    ls=q["Lambda3"]+q["Lambda4"]+q["Lambda5"]
    return (
        m12sq*v2/v1 - .5*q["Lambda1"]*v1*v1-.5*ls*v2*v2,
        m12sq*v1/v2 - .5*q["Lambda2"]*v2*v2-.5*ls*v1*v1,
    )


def source_inputs(order):
    r=audit(order)
    blocks,_=read_slha(Path(r["spectrum_path"]).read_text())
    data={n:block_1d(blocks[n],n) for n in
          ("HMIX","ZEEQUARTIC","ZEEPOT")}
    q={key:data["ZEEQUARTIC"][i][0] for i,key in enumerate(
       ("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
        "LambdaC","K1","K2","K12"),1)}
    pot=data["ZEEPOT"]; hmix=data["HMIX"]
    return {"quartics":q,"m1sq":pot[1][0],"m2sq":pot[2][0],
            "m12sq":pot[3][0],"mCsq":pot[4][0],"mu4":pot[5][0],
            "vevs":[hmix[102][0],hmix[103][0]],
            "original_spectrum_sha256":r["source"]["spectrum_sha256"],
            "original_card_sha256":r["source"]["card_sha256"]}


def potential(x,p):
    """Exact minimized cubic-phase tree potential (GeV^4 / VDIM^4).

    x=(r1,r2,rC)/VDIM, theta in [0,pi/2], phi in [0,pi].
    With theta=0, phi=0, rC=0, this is the neutral EW ray.
    """
    a,b,c,theta,phase=np.asarray(x,dtype=float)
    if not np.all(np.isfinite([a,b,c,theta,phase])):
        raise ValueError("Nonfinite fields")
    q=validate(p["quartics"])
    r1,r2,rc=a*VDIM,b*VDIM,c*VDIM
    t=math.cos(theta)
    f=math.cos(phase)
    z1,z2,zc=r1*r1,r2*r2,rc*rc
    overlap=r1*r2*t
    mu=p["mu4"]
    v=(p["m1sq"]*z1+p["m2sq"]*z2+p["mCsq"]*zc
       -2*p["m12sq"]*overlap*f
       +.5*q["Lambda1"]*z1*z1+.5*q["Lambda2"]*z2*z2
       +q["Lambda3"]*z1*z2
       +(q["Lambda4"]+q["Lambda5"]*math.cos(2*phase))*overlap**2
       +q["LambdaC"]*zc*zc
       +(q["K1"]*z1+q["K2"]*z2
          -2*q["K12"]*overlap*f)*zc
       -2*abs(mu)*r1*r2*rc*math.sin(theta))
    return v/VDIM**4


def ew_normal(p):
    v1,v2=p["vevs"]
    x=np.array([v1/(math.sqrt(2)*VDIM),
                v2/(math.sqrt(2)*VDIM),0.,0.,0.])
    return x,float(potential(x,p))


def check_tadpole(p,rel_tol=1e-4):
    expected=normal_tree_mass_squared(
        p["quartics"],p["vevs"],p["m12sq"])
    native=np.array([p["m1sq"],p["m2sq"]],dtype=float)
    diff=native-np.asarray(expected)
    norm=max(float(np.linalg.norm(native)),float(np.linalg.norm(expected)),1.)
    rel=float(np.linalg.norm(diff)/norm)
    return {"tree_expected_m1sq_m2sq_GeV2":list(expected),
            "native_m1sq_m2sq_GeV2":list(native),
            "native_minus_tree_GeV2":list(diff),
            "relative_residual":rel,
            "native_matches_tree_stationarity":rel<=rel_tol,
            "likely_loop_tadpole_effect":
             "POSSIBLE_BUT_NOT_SOURCE_VALIDATED" if rel>rel_tol else "NOT_REQUIRED"}


def search_witness(p, *, nstarts=48, seed=13, max_field_GeV=3000.,
                   numerical_rel_tol=1e-7):
    """Find deeper field configurations, only a one-sided certificate.

    All starts include known neutral and CB directions; L-BFGS-B
    bounded at declared finite radius. No global-minimum PASS label.
    """
    if nstarts<8 or max_field_GeV<VDIM*2:
        raise ValueError("Vacuum search too small for declared protocol")
    if bfb(p["quartics"])["full_Zee_quartic_BFB"]!="CERTIFIED_SUFFICIENT_PASS":
        return {"status":"BFB_CERTIFICATE_REQUIRED_FOR_THIS_FINITE_RADIUS_SEARCH",
                "globally_stable":None}
    x0,vn=ew_normal(p)
    radii=float(max_field_GeV/VDIM)
    bounds=[(0,radii)]*3+[(0,math.pi/2),(0,math.pi)]
    rng=np.random.default_rng(seed)
    starts=[
        x0,
        np.array([x0[0],x0[1],.02, .35, .0]),
        np.array([x0[0],x0[1],1., math.pi/3, .0]),
        np.array([x0[0],x0[1],1., math.pi/4, math.pi/2]),
        np.array([x0[0],x0[1],0., 0., math.pi]),
        np.array([.0,.0,.3, 0., .0]),
        np.array([.5,.0,2.,math.pi/3, .2]),
        np.array([.0,.5,2.,math.pi/3, .2]),
    ]
    for _ in range(nstarts):
        starts.append(np.array([
            float(rng.uniform(.002, min(4,radii))),
            float(rng.uniform(.002, min(4,radii))),
            float(rng.uniform(0,min(4,radii))),
            float(rng.uniform(0,math.pi/2)),
            float(rng.uniform(0,math.pi))]))
    best=(vn,x0,True)
    converged=0
    for start in starts:
        result=minimize(potential,start,args=(p,),bounds=bounds,
                        method="L-BFGS-B",
                        options={"ftol":1e-12,"gtol":1e-8,"maxiter":900})
        if result.success:converged+=1
        value=float(potential(result.x,p))
        if value<best[0]:best=(value,result.x,result.success)
    fv,x,successful=best
    delta=fv-vn
    threshold=numerical_rel_tol*max(1,abs(vn),abs(fv))
    deeper=delta<-threshold
    r1,r2,rc=x[:3]*VDIM
    qbreaking=bool(rc>1e-4 or (r1*r2*math.sin(x[3]))>1e-4)
    return {
       "status":("DEEPER_WITNESS_FOUND" if deeper else
                 "NO_DEEPER_WITNESS_FOUND_NONCERTIFYING"),
       "normal_tree_potential_div_v4":vn,
       "lowest_found_potential_div_v4":fv,
       "delta_potential_div_v4":delta,
       "deeper_than_EW_normal":deeper,
       "witness_fields_GeV":[float(r1),float(r2),float(rc)],
       "witness_angles_radians":[float(x[3]),float(x[4])],
       "charge_breaking_witness":qbreaking,
       "n_local_starts":len(starts),
       "n_optimizer_converged":converged,
       "best_optimization_report_success":bool(successful),
       "max_fields_GeV":max_field_GeV,
       "all_complex_field_directions_included_by_orbit_parameterization":True,
       "certified_exhaustive_global_search":False,
       "global_EW_vacuum":"FAIL_DEEPER_TREE_WITNESS" if deeper else "UNKNOWN",
    }


def evaluate(order,nstarts=48):
    p=source_inputs(order)
    t=check_tadpole(p)
    # Always compare *two separate* questions: does the native running
    # quadratic potential have tree-stationary EW input? What would a
    # tree-consistent reconstructed potential suggest?
    surrogate=dict(p)
    surrogate["m1sq"],surrogate["m2sq"]=t["tree_expected_m1sq_m2sq_GeV2"]
    trial=search_witness(surrogate,nstarts=nstarts)
    return {
       "status":"NATIVE_SPHENO_TREE_TADPOLE_AUDIT_WITH_EXACT_ORBIT_WITNESSES",
       "point_id":"s00017_"+order,
       "source_card_sha256":p["original_card_sha256"],
       "source_spectrum_sha256":p["original_spectrum_sha256"],
       "native_tadpole_consistency":t,
       "tree_consistent_SURROGATE_only":trial,
       "surrogate_masses_not_native_running_values":not t["native_matches_tree_stationarity"],
       "original_native_tree_global_EW_vacuum":
            (trial["global_EW_vacuum"] if t["native_matches_tree_stationarity"]
             else "UNKNOWN_NEEDS_MATCHED_TREE_TADPOLES"),
       "original_native_global_EW_vacuum_fully_certified":False,
       "globally_allowed":None,
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--nstarts",type=int,default=48)
    a=parser.parse_args()
    out={"status":"EXACT_ORBIT_ZEE_VACUUM_WITNESS_AUDIT",
         "points":{o:evaluate(o,a.nstarts) for o in ("NO","IO")},
         "physical_global_vacuum":"UNKNOWN","globally_allowed":None}
    with a.output.open("x") as f:
        json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
    for k,v in out["points"].items():
        print(k,"tadpole_residual",
              v["native_tadpole_consistency"]["relative_residual"],
              "surrogate",v["tree_consistent_SURROGATE_only"]["status"],
              "native_vacuum",v["original_native_tree_global_EW_vacuum"])


if __name__=="__main__":main()
