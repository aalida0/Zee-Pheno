"""Step 1: one NO/IO Zee neutrino target and inverse-fit card proposal.

NuFIT 6.1 IC24 (including SK) best fits, November 2025 data. This is an
exact central-value target, NOT a NuFIT likelihood or a native SPheno PASS.
"""
import argparse
import json
import math
from pathlib import Path

import numpy as np

from Model.check_native import fix_card
from Model.audit_neutrino_native import sha

# NuFIT 6.1 (2025), IC24 + SK atmospheric, lower half of v61.tbl-parameters.pdf
# dm3l is dm31 for NO, dm32 for IO. Delta CP is in degrees.
NUFIT = {
    "NO": dict(s12=0.3088, s13=0.02248, s23=0.470,
               delta_deg=212., dm21=7.537e-5, dm3l=2.511e-3),
    "IO": dict(s12=0.3088, s13=0.02262, s23=0.550,
               delta_deg=274., dm21=7.537e-5, dm3l=-2.483e-3),
}
SOURCE = "NuFIT 6.1 (2025), IC24 with SK, https://www.nu-fit.org/sites/default/files/v61.tbl-parameters.pdf"
REFERENCE_CARD_SHA = "e83f315928a3906bc14e48aec39b2bd842abf66949f53b6cd78824d4c45a8e07"


def finite_kernel(m1, m2, ml):
    if min(m1, m2) <= 0 or ml < 0:
        raise ValueError("Nonphysical mass")
    if m1 == m2:
        return 0.
    if ml == 0:
        return math.log((m1/m2)**2)
    a, b, z = m1*m1, m2*m2, ml*ml
    if z == a or z == b:
        raise ValueError("Degenerate mass requires finite limiting expression")
    return (math.log(a/b) + z/(a-z)*math.log(a/z)
            - z/(b-z)*math.log(b/z))


def zee_matrix(f, d, leptons, charged, sine):
    leptons = np.asarray(leptons, float)
    weights = 1e9*sine/(32*math.pi**2)*leptons*np.array(
        [finite_kernel(*charged, m) for m in leptons])
    a = (np.asarray(f, complex)*weights[None, :]) @ np.asarray(d, complex)
    return a+a.T


def complex_matrix(value):
    return np.asarray(value["real"])+1j*np.asarray(value["imag"])


def encode(a):
    return {"real": a.real.tolist(), "imag": a.imag.tolist()}


def target(ordering, mlight, alpha21=0., alpha31=0.):
    p = NUFIT[ordering]
    if not math.isfinite(mlight) or not 0 <= mlight <= 0.2:
        raise ValueError("Lightest mass must be finite and between 0 and 0.2 eV")
    if not all(math.isfinite(x) for x in (alpha21, alpha31)):
        raise ValueError("Invalid Majorana phase")
    if ordering == "NO":
        masses = np.array([mlight, math.sqrt(mlight**2+p["dm21"]),
                           math.sqrt(mlight**2+p["dm3l"])])
    else:
        m2 = math.sqrt(mlight**2-p["dm3l"])
        masses = np.array([math.sqrt(m2**2-p["dm21"]), m2, mlight])
    s12,s13,s23 = (math.sqrt(p[k]) for k in ("s12","s13","s23"))
    c12,c13,c23 = (math.sqrt(1-p[k]) for k in ("s12","s13","s23"))
    delta = math.radians(p["delta_deg"])
    e = np.exp(1j*delta)
    u = np.array([
        [c12*c13, s12*c13, s13/e],
        [-s12*c23-c12*s23*s13*e, c12*c23-s12*s23*s13*e, s23*c13],
        [s12*s23-c12*c23*s13*e, -c12*s23-s12*c23*s13*e, c23*c13],
    ], complex) @ np.diag([1., np.exp(.5j*alpha21), np.exp(.5j*alpha31)])
    mnu = u.conj() @ np.diag(masses) @ u.conj().T
    if np.linalg.norm(u.conj().T@u-np.eye(3)) > 1e-12:
        raise ArithmeticError("PMNS input not unitary")
    return mnu, masses, u


def inverse(mnu, leptons, charged, sine, fscale, t=1.+0j):
    """Zee null cone n.T M n=0; minimal-norm D for fixed antisymmetric F."""
    if not np.isfinite(fscale) or fscale <= 0 or abs(sine) < 1e-12:
        raise ValueError("Invalid F scale or scalar mixing")
    # n=(1,x,t) solves n.T @ M @ n=0.
    a = mnu[1, 1]
    b = 2*(mnu[0, 1]+t*mnu[1, 2])
    c = mnu[0, 0]+2*t*mnu[0, 2]+t*t*mnu[2, 2]
    roots = np.roots([a, b, c])
    x = min(roots, key=abs)
    n = np.array([1, x, t], complex)
    f = np.array([[0, n[2], -n[1]], [-n[2], 0, n[0]],
                  [n[1], -n[0], 0]], complex)
    f *= fscale/np.max(abs(f))
    upper = np.triu_indices(3)
    cols = []
    for basis in np.eye(9).reshape(9,3,3):
        cols.append(zee_matrix(f,basis,leptons,charged,sine)[upper])
    matrix = np.array(cols).T
    singular = np.linalg.svd(matrix, compute_uv=False)
    rank = int(np.sum(singular > singular[0]*1e-12))
    if rank != 5:
        raise ValueError(f"Inverse Zee rank is {rank}, expected 5")
    solution = np.linalg.lstsq(matrix,mnu[upper],rcond=1e-12)[0]
    d = solution.reshape(3,3)
    residual = float(np.linalg.norm(zee_matrix(f,d,leptons,charged,sine)-mnu)
                     /np.linalg.norm(mnu))
    if residual > 1e-8 or abs(n.T@mnu@n)/np.linalg.norm(mnu)>1e-9:
        raise ValueError(f"Inverse Zee solution failed ({residual})")
    return f,d,residual


def fit(snapshot_path, ordering, mlight, alpha21=0., alpha31=0.):
    data = json.loads(snapshot_path.read_text())
    snap = data["snapshot"] if "snapshot" in data else data
    leptons = np.asarray(snap["MFe"], float)
    charged = np.asarray(snap["MHm"],float)[1:]
    zp = np.asarray(snap["ZP"],float)
    v = math.hypot(snap["vd"],snap["vu"])
    cb,sb = snap["vd"]/v,snap["vu"]/v
    sine = 2*zp[1,2]*(-sb*zp[1,0]+cb*zp[1,1])
    el,er = complex_matrix(snap["ZEL"]),complex_matrix(snap["ZER"])
    fref = complex_matrix(snap["fzee"])
    sy = cb*complex_matrix(snap["Ye"])+sb*complex_matrix(snap["Ye2"])
    fscale = float(np.max(abs(fref)))
    mnu,masses,u = target(ordering,mlight,alpha21,alpha31)
    fm,dm,residual = inverse(mnu,leptons,charged,sine,fscale)
    # Invert the exact charged-lepton mass-basis conventions in the native audit.
    fg = el.T @ fm @ el
    dg = er.T @ dm @ el
    ye2 = sb*sy+cb*dg
    # This is a *proposal* at the reference running point. Native matching
    # re-solves Ye and running parameters after Ye2,F have been changed.
    report = dict(status="ANALYTIC_TARGET_CLOSED_NATIVE_UNKNOWN",
                  ordering=ordering, nufit_source=SOURCE, nufit_parameters=NUFIT[ordering],
                  reference_snapshot_sha256=sha(snapshot_path),
                  mlight_eV=mlight, alpha21_rad=alpha21,alpha31_rad=alpha31,
                  target_masses_eV=masses.tolist(),target_Mnu_eV=encode(mnu),
                  target_pmns=encode(u), forward_residual=residual,
                  native_revalidation_required=True,
                  f_scale=fscale,max_abs_D=float(np.max(abs(dm))),
                  F_gauge=encode(fg), Ye2_gauge=encode(ye2))
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--snapshot",type=Path,default=Path(__file__).resolve().parents[1]/
                    "Model/reference/native_neutrino_snapshot.json")
    ap.add_argument("--ordering",choices=("NO","IO"),required=True)
    ap.add_argument("--mlight-ev",type=float,required=True)
    ap.add_argument("--alpha21",type=float,default=0.)
    ap.add_argument("--alpha31",type=float,default=0.)
    ap.add_argument("--reference-card",type=Path,help="Attested loop-enabled native card")
    ap.add_argument("--write-card",type=Path,help="Write new proposed input card (must not exist)")
    ap.add_argument("--output",type=Path,help="Write JSON report (must not exist)")
    args=ap.parse_args()
    if bool(args.reference_card)!=bool(args.write_card):
        ap.error("--reference-card and --write-card must be used together")
    report=fit(args.snapshot,args.ordering,args.mlight_ev,args.alpha21,args.alpha31)
    if args.reference_card:
        if sha(args.reference_card)!=REFERENCE_CARD_SHA:
            raise ValueError("Unrecognized loop-enabled source card")
        names={"FZEEIN":complex_matrix(report["F_gauge"]),
               "YEE2IN":complex_matrix(report["Ye2_gauge"])}
        updates={(("IM" if imag else "")+name,(i+1,j+1)):float(
                 (matrix[i,j].imag if imag else matrix[i,j].real))
                 for name,matrix in names.items() for imag in (False,True)
                 for i in range(3) for j in range(3)}
        proposal=fix_card(args.reference_card.read_text(),updates)
        with args.write_card.open("x") as out:
            out.write(proposal)
        report["proposed_card_sha256"]=sha(args.write_card)
    if args.output:
        with args.output.open("x") as out:
            json.dump(report,out,indent=2)
            out.write("\n")
    print(json.dumps({k:v for k,v in report.items() if k not in
                     ("target_Mnu_eV","target_pmns","F_gauge","Ye2_gauge")},indent=2))


if __name__ == "__main__":
    main()
