"""Restricted CP-conserving Zee scalar model, tree level only.

Conventions match the audited ZeeTypeII scalar reconstruction in
https://github.com/aalida0/zee (scripts/gate6_benchmark.py).
No SPheno, LFV, neutrino-mass or vacuum/global-acceptance claim.
"""
import math


def zee_f(f12: complex, f13: complex, f23: complex):
    """Antisymmetric complex lepton-singlet Yukawa matrix (f = -f.T)."""
    entries = (complex(f12), complex(f13), complex(f23))
    if not all(math.isfinite(z.real) and math.isfinite(z.imag) for z in entries):
        raise ValueError("Fzee entries must be finite")
    a, b, c = entries
    return ((0j, a, b), (-a, 0j, c), (-b, -c, 0j))


def scalar_parameters(*, mh1, mh2, ma, mH1, mH2, alpha, beta, gamma,
                      m12sq, lc, k1, k2, k12, v=246.21965):
    """Invert the tree-level restricted Zee scalar spectrum to Lagrangian parameters.

    Inputs: masses and v in GeV; m12sq in GeV^2; angles in radians.
    gamma may be signed. Lambda6/7 are absent by model definition.
    """
    values = (mh1, mh2, ma, mH1, mH2, alpha, beta, gamma,
              m12sq, lc, k1, k2, k12, v)
    if not all(math.isfinite(x) for x in values):
        raise ValueError("All scalar inputs must be finite")
    if not (0 < mh1 <= mh2 and ma > 0 and 0 < mH1 <= mH2 and v > 0
            and 0.01 < beta < math.pi/2 - 0.01
            and -math.pi/2 < gamma < math.pi/2):
        raise ValueError("Invalid masses or mixing angles")
    sa, ca = math.sin(alpha), math.cos(alpha)
    sb, cb = math.sin(beta), math.cos(beta)
    sg, cg = math.sin(gamma), math.cos(gamma)
    vv = v * v
    soft = m12sq / (sb * cb)
    md = mH1*mH1*cg*cg + mH2*mH2*sg*sg
    mn11 = mh1*mh1*ca*ca + mh2*mh2*sa*sa
    mn22 = mh1*mh1*sa*sa + mh2*mh2*ca*ca
    mn12 = (mh1*mh1 - mh2*mh2)*sa*ca
    l1 = (mn11 - m12sq*sb/cb) / (vv*cb*cb)
    l2 = (mn22 - m12sq*cb/sb) / (vv*sb*sb)
    l5 = (soft - ma*ma) / vv
    l4 = 2*(soft - md)/vv - l5
    l3 = (mn12 + m12sq)/(vv*sb*cb) - l4 - l5
    mu4 = -math.sqrt(2)/v * (mH1*mH1 - mH2*mH2)*cg*sg
    mc2 = mH1*mH1*sg*sg + mH2*mH2*cg*cg - vv/2*(k1*cb*cb + k2*sb*sb - 2*k12*sb*cb)
    l345 = l3 + l4 + l5
    m112 = m12sq*sb/cb - vv/2*(l1*cb*cb + l345*sb*sb)
    m222 = m12sq*cb/sb - vv/2*(l2*sb*sb + l345*cb*cb)
    out = dict(Lambda1=l1, Lambda2=l2, Lambda3=l3, Lambda4=l4,
               Lambda5=l5, LambdaC=lc, K1=k1, K2=k2, K12=k12,
               Mu4=mu4, MCS2=mc2, M12Sq=m12sq, M112=m112, M222=m222)
    if not all(math.isfinite(x) for x in out.values()):
        raise ValueError("Nonfinite reconstructed parameter")
    return out
