"""Independent one-loop scalar T in Zee, Grimus et al. NPB801 (2008) Eq.(28).

Uses actual native SPheno charged and neutral mixing matrices; NOT a
reimplementation or recasting of SPheno's printed T.  All physical
charged/neutral scalars, both charged mixing states and Goldstones
are represented. The source uses the finite SM Higgs subtraction.

This validates *T* only. S and U use different loop functions G and
G-hat and CANNOT be inferred from T. A full S,T,U theory gate requires
all three independently computed quantities and consistent EW schemes.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from ZeeHiggs.spheno_audit import audit, read_slha
from Scan.zee_ewpo import block_1d

CITE="Grimus, Lavoura, Ogreid, Osland, Nucl.Phys.B801 (2008) 81-96, Eq.28"


def f_kernel(x,y):
    x=float(x);y=float(y)
    if not math.isfinite(x) or not math.isfinite(y) or min(x,y)<=0:
        raise ValueError("Scalar kernel needs two finite positive squared masses")
    if math.isclose(x,y,rel_tol=1e-7):
        d=(x-y)/(x+y)
        return (x+y)*(d*d/3+d**4/15+d**6/35)
    return .5*(x+y)-x*y/(x-y)*math.log(x/y)


def read_matrix(blocks,name,size):
    if name not in blocks: raise ValueError("Missing native "+name)
    a=np.full((size,size),np.nan)
    for row,_ in blocks[name]:
        if len(row)!=3: raise ValueError("Malformed native mixing matrix")
        i,j=(int(z)-1 for z in row[:2])
        if not (0<=i<size and 0<=j<size) or math.isfinite(a[i,j]):
            raise ValueError("Duplicate or invalid native "+name+" indices")
        a[i,j]=float(row[2].replace("D","E"))
    if not np.isfinite(a).all() or not np.allclose(a@a.T,np.eye(size),rtol=0,atol=2e-4):
        raise ValueError("Nonorthogonal or incomplete native "+name+" matrix")
    return a


def source_point(order):
    base=audit(order)
    blocks,_=read_slha(Path(base["spectrum_path"]).read_text())
    h=read_matrix(blocks,"SCALARMIX",2)
    a=read_matrix(blocks,"PSEUDOSCALARMIX",2)
    c=read_matrix(blocks,"CHARGEMIX",3)
    mass=block_1d(blocks["MASS"],"MASS")
    gauge=block_1d(blocks["GAUGE"],"GAUGE")
    low=block_1d(blocks["SPHENOLOWENERGY"],"SPHENOLOWENERGY")
    for pid in (25,35,36,37,9000037,23,24):
        if pid not in mass:raise ValueError("Missing particle mass "+str(pid))
    # In native SPheno convention each *row* denotes mass state
    # (G0,h,H,A and G+,H1+,H2+ respectively).  This is a complex V
    # with the neutral CP-odd states multiplied by i, Eq.22/24.
    U=c[:,:2].T
    V=np.column_stack((1j*a[0,:],h[0,:],h[1,:],1j*a[1,:]))
    if not np.allclose(V.conj().T@V,(V.conj().T@V).conj().T,atol=1e-9):
        raise ValueError("Non-Hermitian neutral mixing metric")
    masses_c=[None,mass[37][0],mass[9000037][0]]
    masses_n=[None,mass[25][0],mass[35][0],mass[36][0]]
    g1,g2=gauge[1][0],gauge[2][0]
    if not (0<g1<1 and 0<g2<1):raise ValueError("Bad native EW gauge couplings")
    return {"source":base["source"],"U":U,"V":V,
            "masses_charged":masses_c,"masses_neutral":masses_n,
            "mZ":mass[23][0],"mW":mass[24][0],
            "g1":g1,"g2":g2,"native_T":low[1][0],
            "native_S":low[2][0],"native_U":low[3][0],"order":order}


def independent_t(p,weak_angle="running_gauge"):
    """Eq.28 full multi-doublet + charged-singlet scalar T.

    weak_angle='running_gauge' uses source GAUGE g1,g2 at 160 GeV.
    'on_shell' uses 1-mW^2/mZ^2 as a separate scheme sensitivity,
    not a statistical nuisance.
    """
    U,V=p["U"],p["V"]
    if U.shape!=(2,3) or V.shape!=(2,4):
        raise ValueError("Zee charged/neutral state matrices mismatch")
    UDU=U.conj().T@U
    UDV=U.conj().T@V
    VDV=V.conj().T@V
    mc,mn=p["masses_charged"],p["masses_neutral"]
    z,w=p["mZ"],p["mW"]
    if weak_angle=="running_gauge":
        sw2=p["g1"]**2/(p["g1"]**2+p["g2"]**2)
    elif weak_angle=="on_shell":
        sw2=1-(w/z)**2
    else:
        raise ValueError("Unsupported sin^2 thetaW convention")
    if not 0<sw2<1:raise ValueError("Invalid weak mixing angle")
    f=lambda a,b:f_kernel(a*a,b*b)
    one=sum(abs(UDV[a,b])**2*f(mc[a],mn[b])
            for a in (1,2) for b in (1,2,3))
    two=sum(np.imag(VDV[b,bp])**2*f(mn[b],mn[bp])
            for b in (1,2,3) for bp in range(b+1,4))
    three=2*abs(UDU[1,2])**2*f(mc[1],mc[2])
    four=3*sum(np.imag(VDV[0,b])**2*
               (f(z,mn[b])-f(w,mn[b])) for b in (1,2,3))
    five=3*(f(z,mn[1])-f(w,mn[1]))
    result=(one-two-three+four-five)/(16*math.pi*sw2*w*w)
    if not math.isfinite(result):raise ValueError("Nonfinite T independent result")
    return {
      "independent_scalar_T":float(result),
      "SPheno_T":float(p["native_T"]),
      "difference":float(result-p["native_T"]),
      "relative_difference":float((result-p["native_T"])/p["native_T"])
        if p["native_T"] else None,
      "weak_angle_definition":weak_angle,
      "weak_angle_sin2":sw2,
      "parts_GeV2":{"positive_charged_neutral":float(one),
        "neutral_pairs":float(two),"charged_mixing_pairs":float(three),
        "gauge_neutral":float(four),"SM_higgs_subtraction":float(five)},
      "charged_states_included":[37,9000037],
      "neutral_states_included":[25,35,36],
      "independently_validated_STU_complete":False,
      "theory_T_same_mass_and_scheme_certified":False,
    }


def evaluate(order):
    p=source_point(order)
    variants={s:independent_t(p,s) for s in ("running_gauge","on_shell")}
    return {"status":"INDEPENDENT_2008_ZEE_SCALAR_T_NUMERICAL_CROSSCHECK",
       "point_id":"s00017_"+order,
       "source":p["source"],"reference":CITE,
       "results":variants,
       "full_STU_independently_checked":False,
       "physical_full_EWPO":"UNKNOWN"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    a=parser.parse_args()
    out={"status":"INDEPENDENT_ZEE_T_ONLY_VALIDATION",
         "points":{o:evaluate(o) for o in ("NO","IO")},
         "all_STU_independently_validated":False,
         "physical_full_EWPO":"UNKNOWN","globally_allowed":None}
    with a.output.open("x") as f:
        json.dump(out,f,indent=2,allow_nan=False);f.write("\n")
    for k,v in out["points"].items():
        p=v["results"]["running_gauge"]
        print(k,"independent_T",p["independent_scalar_T"],
              "SPheno_T",p["SPheno_T"],
              "delta",p["difference"],"STU_full=UNKNOWN")


if __name__=="__main__":main()
