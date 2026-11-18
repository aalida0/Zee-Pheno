"""Audited photon Barr-Zee subsets of the restricted Zee model.

Adapted from aalida0/zee module4/{kernels,amplitudes,charged_scalar,couplings}.py,
branch research/zee-full-constraints. See module4/MATCHING_AUDIT.md and
module4/CHARGED_SCALAR_AUDIT.md there. Not complete electroweak LFV.
"""
from functools import lru_cache
import math

import numpy as np
from scipy.integrate import quad


def _ratio(u, z):
    t = (u-z)/z
    if abs(t) < 1e-4:
        return 1-t/2+t*t/3-t**3/4+t**4/5-t**5/6
    return math.log(u/z)/t


@lru_cache(maxsize=2048)
def photon_kernels(z):
    """Positive fermionic f,g and negative scalar F2; numerical errors."""
    z=float(z)
    if not math.isfinite(z) or not 1e-24<=z<=1e10:
        raise ValueError("Barr-Zee mass ratio outside validated integration range")
    def fn(x, k):
        u=x*(1-x)
        return .5*_ratio(u,z)*(1-2*u if k=="f" else 1) if k!="F2" else -u*_ratio(u,z)/(2*z)
    points=[]
    if z<.25:
        points=[2*z/(1+math.sqrt(1-4*z))]
    rows=[tuple(2*t for t in quad(fn,0,.5,args=(k,),epsabs=min(1e-12,5e-10*z),
                  epsrel=2e-12,points=points,limit=500))
          for k in ("f","g","F2")]
    if any(not math.isfinite(x) or err>max(2e-10,abs(x)*2e-9)
           for x,err in rows):
        raise ArithmeticError("Barr-Zee quadrature did not converge")
    return tuple(x for x,_ in rows),tuple(err for _,err in rows)


def quark_coupling(species, mf, scalar, p):
    """Type-II neutral quark vertex -i(C PL+C* PR), no CKM factor."""
    alpha,beta,v=p["alpha"],p["beta"],p["v"]
    if species in ("up","charm","top"):
        return dict(h=mf/v*math.sin(alpha)/math.sin(beta),
                    H=mf/v*math.cos(alpha)/math.sin(beta),
                    A=1j*mf/v/math.tan(beta))[scalar]
    if species in ("down","strange","bottom"):
        return dict(h=mf/v*math.cos(alpha)/math.cos(beta),
                    H=-mf/v*math.sin(alpha)/math.cos(beta),
                    A=1j*mf/v*math.tan(beta))[scalar]
    raise ValueError("Unknown quark species")


QUARKS={"up":(3,2/3),"charm":(3,2/3),"top":(3,2/3),
        "down":(3,-1/3),"strange":(3,-1/3),"bottom":(3,-1/3)}
SPECIES=tuple(QUARKS)+("electron","muon","tau")


def fermion_photon(external, internal, mass, scalar_mass, mi, alpha_em, nc, charge):
    """Signed CL,CR from module4/MATCHING_AUDIT.md, eqs for photon subset."""
    vals=(mass,scalar_mass,mi,alpha_em,nc)
    if any(not math.isfinite(v) or v<=0 for v in vals):
        raise ValueError("Positive finite Barr-Zee masses/constants required")
    ext=np.asarray(external,complex)
    if ext.shape!=(2,) or not np.isfinite(ext).all() or not np.isfinite(internal):
        raise ValueError("Invalid chiral/inner scalar Yukawa")
    (f,g,_),errors=photon_kernels((mass/scalar_mass)**2)
    effective=internal.real*f+1j*internal.imag*g
    pref=-alpha_em*nc*charge**2/(8*math.pi**3*mi*mass)
    return (pref*np.array([ext[0]*effective,ext[1]*effective.conjugate()]),
            abs(pref)*max(map(abs,ext))*(abs(internal.real)*errors[0]+abs(internal.imag)*errors[1]))


def charged_scalar_photon(external, trilinear, mass, scalar_mass, mi, alpha_em):
    """Signed physical charged-scalar photon subset; G from L=-G*s*H-H+."""
    vals=(mass,scalar_mass,mi,alpha_em)
    if any(not math.isfinite(v) or v<=0 for v in vals):
        raise ValueError("Positive finite scalar Barr-Zee input required")
    if abs(complex(trilinear).imag)>1e-12*max(1.,abs(trilinear)):
        raise ValueError("Physical diagonal trilinear is real for CP-even scalar")
    ext=np.asarray(external,complex)
    if ext.shape!=(2,) or not np.isfinite(ext).all():
        raise ValueError("Invalid chiral external coupling")
    (_,_,f2),(_,_,error)=photon_kernels((mass/scalar_mass)**2)
    pref=alpha_em*complex(trilinear).real/(16*math.pi**3*mi*scalar_mass**2)
    return pref*f2*ext,abs(pref)*max(map(abs,ext))*error


def physical_trilinears(p, potential):
    """Diagonal G_h/H of physical H1,H2 from the restricted scalar potential."""
    names=("Lambda1","Lambda2","Lambda3","Lambda4","Lambda5",
           "K1","K2","K12","Mu4")
    if set(potential)!=set(names) or not all(math.isfinite(potential[k]) for k in names):
        raise ValueError("Complete finite restricted scalar potential is required")
    v1,v2=p["v"]*math.cos(p["beta"]),p["v"]*math.sin(p["beta"])
    l1,l2,l3,l45=(potential["Lambda1"],potential["Lambda2"],
                   potential["Lambda3"],potential["Lambda4"]+potential["Lambda5"])
    mu=potential["Mu4"]/math.sqrt(2)
    g1=np.array([[l1*v1,l45*v2/2,0],[l45*v2/2,l3*v1,-mu],
                 [0,-mu,potential["K1"]*v1-potential["K12"]*v2]])
    g2=np.array([[l3*v2,l45*v1/2,mu],[l45*v1/2,l2*v2,0],
                 [mu,0,potential["K2"]*v2-potential["K12"]*v1]])
    ca,sa=math.cos(p["alpha"]),math.sin(p["alpha"])
    cb,sb=math.cos(p["beta"]),math.sin(p["beta"])
    cg,sg=math.cos(p["gamma"]),math.sin(p["gamma"])
    zp=np.array([[cb,sb,0],[-sb*cg,cb*cg,sg],[sb*sg,-cb*sg,cg]])
    return {name:np.diag((zp@matrix@zp.T)[1:,1:]).copy()
            for name,matrix in (("h",ca*g1+sa*g2),("H",-sa*g1+ca*g2))}
