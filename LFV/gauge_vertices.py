"""Exact tree-level electroweak current matrices for Zee charged scalars.

Basis: (H_D+, S+); physical H1+=c_gamma H_D+ + s_gamma S+,
H2+=-s_gamma H_D+ + c_gamma S+. Q=+1, T3=(1/2,0).
These matrices are NOT a validated LFV two-loop gauge amplitude.
"""
import math
import numpy as np


def charged_currents(gamma, sw2):
    if not math.isfinite(gamma) or not math.isfinite(sw2) or not 0<sw2<1:
        raise ValueError("Require finite gamma, 0<sin²thetaW<1")
    c,s=math.cos(gamma),math.sin(gamma)
    rotation=np.array([[c,s],[-s,c]])
    doublet=rotation[:,0]
    singlet=rotation[:,1]
    # Covariant derivative: g_Z/e = (T3-sW² Q)/(sW cW)
    t3=np.outer(doublet,doublet)/2
    photon=np.eye(2)
    z=(t3-sw2*photon)/math.sqrt(sw2*(1-sw2))
    return {
        "charged_rotation":rotation,
        "doublet_fractions":doublet,
        "singlet_fractions":singlet,
        "T3_mass_basis":t3,
        "photon_over_e":photon,
        "Z_over_e":z,
        "W_doublet_projector":np.outer(doublet,doublet),
        "status":"TREE_EW_VERTEX_IDENTITIES_ONLY",
    }
