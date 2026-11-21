"""Restricted Zee H/A -> charged lepton pair: scoped LO partial widths.

Required inputs: a tree-level physical scalar dictionary, reconstructed scalar
quartics, and dimensionless 3x3 complex Higgs-basis D matrix in the charged
lepton mass basis. No LFV acceptance, SPheno match, or pp cross section claim.

Convention: h = cos(alpha)*rho1+sin(alpha)*rho2, H = -sin(alpha)*rho1
+cos(alpha)*rho2, S=cos(beta)*Ye+sin(beta)*Ye2=sqrt(2)*M_l/v,
D=-sin(beta)*Ye+cos(beta)*Ye2.
"""
import math
import numpy as np

ML = np.array([0.00051099895, 0.1056583755, 1.77686])


def neutral_yukawas(D, alpha, beta, v=246.21965):
    D = np.asarray(D, dtype=complex)
    if D.shape != (3, 3) or not np.isfinite(D).all():
        raise ValueError("D must be a finite complex 3x3 matrix")
    S = np.diag(math.sqrt(2) * ML / v)
    c, s = math.cos(beta-alpha), math.sin(beta-alpha)
    return {"h": (c*S-s*D)/math.sqrt(2),
            "H": (s*S+c*D)/math.sqrt(2),
            "A": D/math.sqrt(2)}


def inclusive_pair_width(mass, Y, i, j):
    """Charge-summed H/A->l_i l_j, i!=j, massless final leptons."""
    if not (0 <= i < 3 and 0 <= j < 3 and i != j):
        raise ValueError("Distinct physical lepton flavours required")
    return mass/(8*math.pi)*(abs(Y[i, j])**2+abs(Y[j, i])**2)


def fermion_width(mass, yukawa, kinematic_mass, colors, cp):
    if cp not in ("scalar", "pseudo") or colors not in (1, 3):
        raise ValueError("Invalid CP or color")
    if mass <= 2*kinematic_mass:
        return 0.
    beta = math.sqrt(1-4*(kinematic_mass/mass)**2)
    return colors*mass/(8*math.pi)*abs(yukawa)**2*beta**(3 if cp == "scalar" else 1)


def H_hh_trilinear(physical, quartics):
    a, b = physical["alpha"], physical["beta"]
    ca, sa, cb, sb = math.cos(a), math.sin(a), math.cos(b), math.sin(b)
    h, H = (ca, sa), (-sa, ca)
    v1, v2 = physical.get("v", 246.21965)*cb, physical.get("v", 246.21965)*sb
    l1, l2 = quartics["Lambda1"], quartics["Lambda2"]
    l345 = sum(quartics[k] for k in ("Lambda3", "Lambda4", "Lambda5"))
    return (3*l1*v1*H[0]*h[0]**2+3*l2*v2*H[1]*h[1]**2+
            l345*(v1*(H[0]*h[1]**2+2*h[0]*H[1]*h[1])+
                  v2*(H[1]*h[0]**2+2*h[0]*H[0]*h[1])))


def _vv_width(mass, mv, v, identical=False):
    if mass <= 2*mv: return 0.
    x = (mv/mass)**2
    return (mass**3/((32 if identical else 16)*math.pi*v*v)*
            math.sqrt(1-4*x)*(1-4*x+12*x*x))


def _a_zh_width(mass, mh, mZ, v, sin_beta_minus_alpha):
    if mass <= mh+mZ: return 0.
    lam = ((mass-mh-mZ)*(mass+mh+mZ)*(mass-mh+mZ)*(mass+mh-mZ))
    return (sin_beta_minus_alpha**2*lam**1.5/(16*math.pi*v*v*mass**3))


def included_widths(physical, quartics, D, *, mt_yuk=149.4, mb_yuk=2.717,
                    mt_kin=172.4, mb_kin=4.78, mW=80.379, mZ=91.1876):
    """Partial tree rates and OPTIMISTIC mu/tau branching ratio bounds.

    No gg/gammagamma/Zgamma or full off-shell/loop widths, no hadronic sigma.
    Running Yukawa masses and pole kinematics must be given consistently.
    """
    v = physical.get("v", 246.21965)
    b, a = physical["beta"], physical["alpha"]
    mh, mH, mA = (physical[k] for k in ("mh1", "mh2", "ma"))
    sin_ba = math.sin(b-a)
    Y = neutral_yukawas(D, a, b, v)
    tops = {"H": mt_yuk/v*math.cos(a)/math.sin(b),
            "A": mt_yuk/v/math.tan(b)}
    bottoms = {"H": -mb_yuk/v*math.sin(a)/math.cos(b),
               "A": mb_yuk/v*math.tan(b)}
    result = {}
    for name, mass, cp in (("H", mH, "scalar"), ("A", mA, "pseudo")):
        y = Y[name]
        parts = {"tt": fermion_width(mass, tops[name], mt_kin, 3, cp),
                 "bb": fermion_width(mass, bottoms[name], mb_kin, 3, cp)}
        for i, flavour in enumerate(("ee", "mumu", "tautau")):
            parts[flavour] = fermion_width(mass, y[i,i], ML[i], 1, cp)
        for (i,j),flavour in [((0,1),"emu"), ((0,2),"etau"), ((1,2),"mutau")]:
            parts[flavour] = inclusive_pair_width(mass,y,i,j)
        if name == "H":
            g = H_hh_trilinear(physical, quartics)
            parts["hh"] = 0. if mass<=2*mh else g*g/(32*math.pi*mass)*math.sqrt(1-4*(mh/mass)**2)
            parts["WW"] = sin_ba**2*_vv_width(mass,mW,v)
            parts["ZZ"] = sin_ba**2*_vv_width(mass,mZ,v,True)
        else:
            parts["Zh"] = _a_zh_width(mass,mh,mZ,v,sin_ba)
        included = sum(parts.values())
        result[name] = {"mass_GeV":mass,"partial_widths_GeV":parts,
                        "total_included_GeV":included,
                        "BR_mu_tau_optimistic_upper":parts["mutau"]/included,
                        "acceptance":"UNKNOWN"}
    return result
