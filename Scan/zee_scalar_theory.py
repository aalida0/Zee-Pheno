"""Restricted CP-conserving Zee scalar theory at one renormalization scale.

Potential conventions are EXACTLY Eq. (17) in EPJC 81 (2021) 1087,
doi:10.1140/epjc/s10052-021-09956-2 (two SU2 doublets + charged
complex singlet, no quartic lambda6/lambda7).

Analytic FULL-ZEE-quartic BFB SUFFICIENT proof:
  V4 = V4_2HDM(Phi1,Phi2) + lambda_C*(chi*chi)^2
       +(chi*chi) Phi^dagger [[k1,-k12],[-k12,k2]] Phi.
If the 2HDM quartic is strictly BFB, lambda_C>0 and the Hermitian
portal matrix is positive SEMIDEFINITE, then V4 >= 0 on ALL fields.
This is sufficient, not necessary: failure of portal PSD must not
be called physical BFB FAIL without an explicit negative direction.
We also check necessary pure singlet/doublet rays, Eq. (41), (45c).

Full high-energy TREE perturbative unitarity:
Eigenvalues of all published scattering matrices E.22-E.29 with
the SAME source-potential conventions, including 13x13 (Q=0,Y=0).
Condition max(abs(eigenvalues)) < 8pi (not max positive eigenvalue
and NOT the 4pi coupling prior).  Only the asymptotic high-energy
tree constraint is claimed; no finite-energy cubic exchanges, loop
unitarity, perturbative running or global EW vacuum certification.
"""
import math

import numpy as np

NAMES=("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
       "LambdaC","K1","K2","K12")
SOURCE={
 "paper":"Off diagonal charged scalar couplings with the Z boson: Zee-type models as an example",
 "doi":"10.1140/epjc/s10052-021-09956-2",
 "url":"https://link.springer.com/article/10.1140/epjc/s10052-021-09956-2",
 "potential":"Eq.17","BFB":"Eq.41,45; PSD portal sufficient by sum of nonnegative quartic forms",
 "tree_unitarity":"Eq.57 and App.E E.22--E.29",
}


def validate(p):
    if not isinstance(p,dict) or not set(NAMES)<=set(p):
        raise ValueError("All nine source-matched restricted Zee quartics are required")
    out={}
    for k in NAMES:
        v=p[k]
        if type(v) not in (float,int) or not math.isfinite(v):
            raise ValueError("Invalid or nonfinite Zee quartic "+k)
        out[k]=float(v)
    return out


def bfb(p):
    p=validate(p)
    l1,l2,l3,l4,l5,lc,k1,k2,k12=(p[k] for k in NAMES)
    doublet=(l1>0 and l2>0 and
             l3+math.sqrt(max(0,l1*l2))>0 and
             l3+l4-abs(l5)+math.sqrt(max(0,l1*l2))>0)
    pure_singlet=lc>0
    portal_ray_1=(lc>=0 and l1>=0 and
                  k1+math.sqrt(max(0,2*lc*l1))>=0)
    portal_ray_2=(lc>=0 and l2>=0 and
                  k2+math.sqrt(max(0,2*lc*l2))>=0)
    necessary=(doublet and pure_singlet and portal_ray_1 and portal_ray_2)
    determinant=k1*k2-k12*k12
    portal_psd=k1>=0 and k2>=0 and determinant>=0
    # Important: PSD portal + strict doublet BFB is sufficient for
    # ALL complex field configurations, not only a neutral-field ray.
    sufficient=doublet and pure_singlet and portal_psd
    return {
        "doublet_BFB_necessary":"PASS" if doublet else "FAIL",
        "pure_singlet_BFB_necessary":"PASS" if pure_singlet else "FAIL",
        "singlet_doublet_axis_BFB_necessary":
            "PASS" if portal_ray_1 and portal_ray_2 else "FAIL",
        "portal_hermitian_PSD_sufficient":"PASS" if portal_psd else "NOT_PROVEN",
        "portal_matrix_determinant":determinant,
        "full_Zee_quartic_BFB":
            ("CERTIFIED_SUFFICIENT_PASS" if sufficient else
             "FAIL_NECESSARY" if not necessary else "UNKNOWN"),
        "BFB_is_global_EW_vacuum_proof":False,
        "BFB_scheme":"TREE_QUARTIC_AT_GIVEN_SCALE",
    }


def scattering_matrices(p):
    """Published E.22--E.29 in the Eq.17 scalar-potential convention.

    Matrices are REAL symmetric for the restricted CP-conserving model.
    """
    p=validate(p)
    a,b,c,d,e,lc,k1,k2,k12=(p[k] for k in NAMES)
    e22=np.array([
        [a,0,0,e,0,0],
        [0,c+d,0,0,0,0],
        [0,0,k1,0,-k12,0],
        [e,0,0,b,0,0],
        [0,0,-k12,0,k2,0],
        [0,0,0,0,0,2*lc],
    ],dtype=float)
    e23=np.array([
        [a,0,0,e,0,0],
        [0,c,d,0,0,0],
        [0,d,c,0,0,0],
        [e,0,0,b,0,0],
        [0,0,0,0,k1,-k12],
        [0,0,0,0,-k12,k2],
    ],dtype=float)
    e24=np.array([
        [a,0,0,d,0,0],
        [0,c,e,0,0,0],
        [0,e,c,0,0,0],
        [d,0,0,b,0,0],
        [0,0,0,0,k1,-k12],
        [0,0,0,0,-k12,k2],
    ],dtype=float)
    lam34=c+d
    # E.27: 7x7 A
    A=np.array([
        [2*a,0,0,0,lam34,0,0],
        [0,2*e,0,lam34,0,0,0],
        [0,0,0,0,0,0,k1],
        [0,lam34,0,2*e,0,0,0],
        [lam34,0,0,0,2*b,0,0],
        [0,0,0,0,0,0,-k12],
        [0,0,k1,0,0,-k12,0],
    ],dtype=float)
    # E.28: 7x6 B
    B=np.array([
        [0,k1,a,0,0,c],
        [0,-k12,0,e,d,0],
        [-k12,0,0,0,0,0],
        [0,-k12,0,d,e,0],
        [0,k2,c,0,0,b],
        [k2,0,0,0,0,0],
        [0,0,0,0,0,0],
    ],dtype=float)
    # E.29: 6x6 C
    C=np.array([
        [0,0,0,0,0,0],
        [0,4*lc,k1,-k12,-k12,k2],
        [0,k1,2*a,0,0,lam34],
        [0,-k12,0,2*e,lam34,0],
        [0,-k12,0,lam34,2*e,0],
        [0,k2,lam34,0,0,2*b],
    ],dtype=float)
    e26=np.block([[A,B],[B.T,C]])
    result={"Q2_Y1_E22":e22,"Q1_Y1_E23":e23,
            "Q1_Y0_E24":e24,"Q0_Y1_E25":np.zeros((3,3)),
            "Q0_Y0_E26_to_E29":e26}
    for name,m in result.items():
        if not np.allclose(m,m.T,rtol=0,atol=1e-12):
            raise ValueError("Unitarity scattering matrix is not symmetric: "+name)
    return result


def independently_known_eigenvalues(p):
    """15 published analytic independent eigenvalues, Eq. E.31.

    The four nontrivial 13x13-sector roots are tested by full numerical
    diagonalization of the ACTUAL E.26/E.29 matrix, avoiding the
    typographical eta^2 duplicate in the printed quartic polynomial.
    """
    p=validate(p)
    a,b,c,d,e,lc,k1,k2,k12=(p[k] for k in NAMES)
    mixk=math.hypot(k1-k2,2*k12)
    mix4=math.hypot(a-b,2*d)
    mix5=math.hypot(a-b,2*e)
    return [
       .5*(k1+k2-mixk),.5*(k1+k2+mixk),c+d,
       .5*(a+b-mix5),.5*(a+b+mix5),2*lc,c-d,
       .5*(a+b-mix4),.5*(a+b+mix4),c-e,c+e,
       -.5*(k1+k2+mixk),.5*(-k1-k2+mixk),
       e-c,-c-2*d+3*e,
    ]


def unitarity(p):
    p=validate(p)
    matrices=scattering_matrices(p)
    eig={name:sorted(float(v) for v in np.linalg.eigvalsh(m))
         for name,m in matrices.items()}
    full=[v for vals in eig.values() for v in vals]
    highest=max(abs(v) for v in full)
    analytic=independently_known_eigenvalues(p)
    # Published analytic eigenvalues must all be in *some* block spectrum.
    # Some appear with multiplicity; no normalization adjustments allowed.
    unmatched=[v for v in analytic if not any(
        math.isclose(v,w,rel_tol=2e-10,abs_tol=2e-9) for w in full)]
    if unmatched:
        raise ValueError("Independent published unitarity eigenvalues missing: "+str(unmatched))
    return {
        "tree_high_energy_full_Zee_unitarity":
            "PASS" if highest<8*math.pi else "FAIL",
        "all_2body_scattering_blocks_evaluated":list(matrices),
        "dimensions":{key:int(m.shape[0]) for key,m in matrices.items()},
        "max_abs_scattering_eigenvalue":highest,
        "limit_abs_eigenvalue":8*math.pi,
        "margin_to_8pi":8*math.pi-highest,
        "analytic_15_eigenvalues_crosschecked":"PASS",
        "eigenvalues_by_QY_block":eig,
        "eigenvalue_test_is_for_quartics_at_single_scale":True,
        "finite_energy_scattering_and_cubic_exchange_included":False,
        "unitarity_is_absolute_quartic_4pi_test":False,
    }


def evaluate(p):
    p=validate(p)
    z=bfb(p);u=unitarity(p)
    mag=max(abs(v) for v in p.values())
    return {
        "status":"PUBLISHED_FULL_HIGH_ENERGY_TREE_UNITARITY_AND_BFB_SUFFICIENT_AUDIT",
        "source":SOURCE,
        "BFB":z,"unitarity":u,
        "max_abs_quartic":mag,
        "max_abs_quartic_lt_4pi_conventional_prior":"PASS" if mag<4*math.pi else "FAIL",
        "full_tree_global_vacuum":"UNKNOWN",
        "RG_running_perturbativity":"UNKNOWN",
        "full_model_globally_allowed":None,
    }
