"""Independent full one-loop scalar S,T,U for Zee, Grimus et al. (2008).

NPB 801 (2008) 81-96, Eq.(28),(30),(31), App.B,C. The general formula
includes arbitrary charged singlets and DOUBLETS and *off-diagonal*
charged mixing via Udagger U, plus neutral CP-even/odd mixing via
Vdagger V. It is exactly applicable to restricted Zee scalar sectors.

G and Ghat are evaluated through INDEPENDENT, stable Feynman-parameter
quadratures (original Eqs B3,B11,C1-C5) rather than hardcoded formulas
with catastrophic cancellation near degenerate masses.

No independent gauge/fermion new-physics self-energy contribution is
silently included. Crosscheck the independent one-loop scalar
predictions against source SPheno's BSM S,T,U from the same physical
spectrum, and treat finite scheme/mass input differences as declared.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from Scan.zee_oblique_T_independent import source_point,independent_t

PAPER="Grimus et al., Nucl.Phys.B801 (2008) 81-96, Eqs.(28),(30),(31),(B3),(B11),(C1)-(C5)"
NUMERICAL_ORDER=64
MAX_DELTA_ABSOLUTE={"S":5e-4,"T":5e-4,"U":1e-4}


def quadrature(order=NUMERICAL_ORDER):
    if order<16 or order>512:raise ValueError("Unsupported quadrature order")
    node,weight=leggauss(order)
    x=(node+1)*.5
    return x,weight*.5


def G(I,J,Q,order=NUMERICAL_ORDER):
    """C.2 as convergent stable Gauss-Legendre integral, dimensionless.

    Eq.B3: A(Q)= (1/(8 pi2))*integral_0^1 Delta*(div-ln Delta) dx.
    Eq.C1: [A(Q)-A(0)]/Q =
        [2-2div+ln I+ln J+G(I,J,Q)]/(96pi2).
    Hence G = -2-ln(I/Q)-ln(J/Q)
      -12 integral [(D-u)ln(D-u)-D ln D] dx,
    with D=(I/Q)(1-x)+(J/Q)x and u=x(1-x). Exact.
    """
    if not all(math.isfinite(z) and z>0 for z in (I,J,Q)):
        raise ValueError("G needs positive finite squared masses")
    i,j=I/Q,J/Q
    x,wt=quadrature(order)
    u=x*(1-x)
    d=i*(1-x)+j*x
    if np.any(d<=u):
        raise ValueError("G kinematics cross physical branch cut")
    integral=np.sum(wt*((d-u)*np.log1p(-u/d)-u*np.log(d)))
    return float(-2-math.log(i)-math.log(j)-12*integral)


def G_tilde(I,J,Q,order=NUMERICAL_ORDER):
    """Eq.C3/C4: Gtilde = 2 integral ln[Delta(Q)/Delta(0)]dx."""
    if not all(math.isfinite(z) and z>0 for z in (I,J,Q)):
        raise ValueError("Gtilde needs positive masses")
    x,wt=quadrature(order)
    d=(I/Q)*(1-x)+(J/Q)*x
    u=x*(1-x)
    if np.any(d<=u):raise ValueError("Gtilde hits physical branch cut")
    return float(2*np.sum(wt*np.log1p(-u/d)))


def G_hat(I,Q,order=NUMERICAL_ORDER):
    """Eq.C5 = G(I,Q,Q)+12*Gtilde(I,Q,Q)."""
    return G(I,Q,Q,order)+12*G_tilde(I,Q,Q,order)


def independent_S_U(p,weak_angle="running_gauge",order=NUMERICAL_ORDER):
    U,V=p["U"],p["V"]
    if U.shape!=(2,3) or V.shape!=(2,4):
        raise ValueError("Unexpected doublet/charged singlet mixing dimension")
    udu=U.conj().T@U
    udv=U.conj().T@V
    vdv=V.conj().T@V
    ma,mn=p["masses_charged"],p["masses_neutral"]
    mW,mZ=p["mW"],p["mZ"]
    if weak_angle=="running_gauge":
        sw2=p["g1"]**2/(p["g1"]**2+p["g2"]**2)
    elif weak_angle=="on_shell":
        sw2=1-(mW/mZ)**2
    else:raise ValueError("Invalid EW weak-angle scheme")
    if not 0<sw2<1:raise ValueError("Invalid weak mixing angle")
    g=lambda x,y,mv:G(x*x,y*y,mv*mv,order)
    h=lambda x,mv:G_hat(x*x,mv*mv,order)
    # Eq.30, bracket in bar S. S=bar S *4sW^2 cW^2/alpha
    # and g^2/alpha=4pi/sW^2 -> bracket/(24pi).
    s1=sum((2*sw2-udu[a,a].real)**2*g(ma[a],ma[a],mZ)
           for a in (1,2))
    s2=2*abs(udu[1,2])**2*g(ma[1],ma[2],mZ)
    s3=sum(np.imag(vdv[b,bp])**2*g(mn[b],mn[bp],mZ)
           for b in (1,2,3) for bp in range(b+1,4))
    s4=-2*sum(udu[a,a].real*math.log(ma[a]**2) for a in (1,2))
    s4+=sum(vdv[b,b].real*math.log(mn[b]**2) for b in (1,2,3))
    s4-=math.log(mn[1]**2)
    s5=sum(np.imag(vdv[0,b])**2*h(mn[b],mZ) for b in (1,2,3))
    s5-=h(mn[1],mZ)
    s=(s1+s2+s3+s4+s5)/(24*math.pi)
    # Eq.31, full bar U bracket. All states and SM subtraction.
    u1=sum(abs(udv[a,b])**2*g(ma[a],mn[b],mW)
           for a in (1,2) for b in (1,2,3))
    u2=-s1-s2-s3
    u3=sum(np.imag(vdv[0,b])**2*(h(mn[b],mW)-h(mn[b],mZ))
           for b in (1,2,3))
    u3+=-h(mn[1],mW)+h(mn[1],mZ)
    uu=(u1+u2+u3)/(24*math.pi)
    if not all(math.isfinite(float(z)) for z in (s,uu)):
        raise ValueError("Nonfinite independent oblique result")
    return {
      "S":float(s),"U":float(uu),
      "weak_angle_scheme":weak_angle,
      "weak_angle_sin2":sw2,
      "S_bracket":{"charged_diagonal":float(s1),"charged_off_diagonal":float(s2),
         "neutral_pairs":float(s3),"SM_subtracted_logs":float(s4),
         "SM_subtracted_neutral_gauge":float(s5)},
      "U_bracket":{"charged_neutral":float(u1),
                   "minus_Z_couplings":float(u2),
                   "SM_subtracted_neutral_gauge":float(u3)},
      "numerical_G_and_Ghat_quadrature_order":order,
      "two_charged_mass_eigenstates_included":True,
      "SM_subtraction_explicit":True,
    }


def evaluate(order):
    p=source_point(order)
    numeric={}
    for scheme in ("running_gauge","on_shell"):
        t=independent_t(p,scheme)
        su=independent_S_U(p,scheme)
        independent={"S":su["S"],"T":t["independent_scalar_T"],"U":su["U"]}
        native={"S":p["native_S"] if "native_S" in p else None,
                "T":p["native_T"],
                "U":p["native_U"] if "native_U" in p else None}
        if any(z is None for z in native.values()):
            raise ValueError("Native S and U source fields unavailable; must not infer")
        residual={key:independent[key]-native[key] for key in ("S","T","U")}
        consistent=all(abs(residual[key])<MAX_DELTA_ABSOLUTE[key]
                       for key in residual)
        numeric[scheme]={
            "source_SPheno_STU":native,
            "independent_Grimus_scalar_STU":independent,
            "absolute_differences":residual,
            "comparison_tolerance":MAX_DELTA_ABSOLUTE,
            "source_and_independent_agree_at_pragmatic_tolerance":consistent,
            "T_parts":t["parts_GeV2"],
            "S_U_parts":su,
            "same_exact_renormalization_scheme_validated":False,
        }
    return {
        "point_id":"s00017_"+order,
        "status":"FULL_ONE_LOOP_SCALAR_STU_INDEPENDENT_CROSSCHECK",
        "reference":PAPER,
        "source":p["source"],
        "two_charged_states":True,
        "results":numeric,
        "independent_one_loop_scalar_STU_pragmatic_consistency":
           numeric["running_gauge"]["source_and_independent_agree_at_pragmatic_tolerance"],
        "on_shell_weak_angle_scheme_is_a_sensitivity_variation":True,
        "oblique_experimental_screen_further_requires_PDG_covariance":True,
        "complete_SPheno_finite_matching_scheme_certified":False,
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    a=p.parse_args()
    rows={o:evaluate(o) for o in ("NO","IO")}
    output={"status":"INDEPENDENT_2008_GRIMUS_ZEE_SCALAR_STU_COMPLETE",
            "points":rows,
            "all_points_pragmatically_match_SPheno":
             all(z["independent_one_loop_scalar_STU_pragmatic_consistency"]
                 for z in rows.values()),
            "model_scalar_one_loop_oblique_screen_ready":True,
            "full_Zee_nonoblique_EW_corrections_included":False,
            "global_vacuum":"UNKNOWN",
            "globally_allowed":None}
    with a.output.open("x") as f:
        json.dump(output,f,indent=2,allow_nan=False);f.write("\n")
    for k,v in rows.items():
        p=v["results"]["running_gauge"]
        print(k,"S/T/U Grimus",p["independent_Grimus_scalar_STU"],
              "SPheno",p["source_SPheno_STU"],
              "crosscheck",p["source_and_independent_agree_at_pragmatic_tolerance"])


if __name__=="__main__":main()
