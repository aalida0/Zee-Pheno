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


def inverse(mnu, leptons, charged, sine, fscale, t=1.+0j, *,
            root=0, homogeneous=None, f_phase=0., chart=(0,1,2)):
    """Two null-cone roots and four-complex-dimensional homogeneous D freedom."""
    if (not np.isfinite(fscale) or fscale <= 0 or abs(sine) < 1e-12 or
        not np.isfinite(t) or type(root) is not int or root not in (0,1) or
        not math.isfinite(f_phase) or
        not isinstance(chart,(tuple,list)) or len(chart)!=3 or
        any(type(i) is not int for i in chart) or sorted(chart)!=[0,1,2]):
        raise ValueError("Invalid F scale, null-cone root, or scalar mixing")
    coeff = np.zeros(4,complex) if homogeneous is None else np.asarray(homogeneous,complex)
    if coeff.shape!=(4,) or not np.isfinite(coeff).all():
        raise ValueError("Homogeneous D requires exactly four finite complex coefficients")
    # Each chart covers a distinct projective F null-vector patch.
    q=np.asarray(mnu)[np.ix_(chart,chart)]
    a = q[1, 1]
    b = 2*(q[0, 1]+t*q[1, 2])
    c = q[0, 0]+2*t*q[0, 2]+t*t*q[2, 2]
    if abs(a) < 1e-15*max(abs(b),abs(c),1e-300):
        if root or abs(b)<1e-15*max(abs(c),1e-300):
            raise ValueError("Unavailable null-cone root; choose another chart or t")
        x=-c/b
    else:
        roots=sorted(np.roots([a,b,c]),key=lambda z:(abs(z),z.real,z.imag))
        x=roots[root]
    n = np.empty(3,complex)
    n[list(chart)]=[1,x,t]
    f = np.array([[0, n[2], -n[1]], [-n[2], 0, n[0]],
                  [n[1], -n[0], 0]], complex)
    f *= (fscale/np.max(abs(f)))*np.exp(1j*f_phase)
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
    if np.any(coeff):
        # Four complex directions in the kernel of the symmetric
        # neutrino-mass map. This explores Yukawa families, not LFV tuning.
        vh=np.linalg.svd(matrix,full_matrices=True)[2]
        directions=vh[rank:].conj().reshape(4,3,3)
        if directions.shape!=(4,3,3) or np.linalg.norm(
            matrix@directions.reshape(4,9).T
        )/np.linalg.norm(matrix)>1e-11:
            raise ValueError("Unresolved four-dimensional neutrino kernel")
        d=d+np.einsum("k,kij->ij",coeff,directions)
    residual = float(np.linalg.norm(zee_matrix(f,d,leptons,charged,sine)-mnu)
                     /np.linalg.norm(mnu))
    if residual > 1e-8 or abs(n.T@mnu@n)/np.linalg.norm(mnu)>1e-9:
        raise ValueError(f"Inverse Zee solution failed ({residual})")
    return f,d,residual


def fit(snapshot_path, ordering, mlight, alpha21=0., alpha31=0.):
    data = json.loads(snapshot_path.read_text())
    # A native diagnostic run may contain multiple OneLoopFv snapshots;
    # use the invocation whose masses actually appear in the MASS block.
    snap = (data["snapshots"][data["spectrum_snapshot_index"]]
            if "snapshots" in data else data["snapshot"] if "snapshot" in data else data)
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
                  sum_mnu_eV=float(sum(masses)),
                  m_beta_eV=float(np.sqrt(np.sum(abs(u[0])**2 * masses**2))),
                  m_bb_eV=float(abs(mnu[0,0])),
                  target_pmns=encode(u), forward_residual=residual,
                  native_revalidation_required=True,
                  f_scale=fscale,max_abs_D=float(np.max(abs(dm))),
                  F_gauge=encode(fg), Ye2_gauge=encode(ye2))
    return report



# NuFIT 6.1 IC24+SK 3-sigma ranges (one-dimensional projections only).
THREE_SIGMA = {
    "NO": dict(s12=(.2893,.3295), s13=(.02064,.02418),
               s23=(.435,.584), dm21=(7.236e-5,7.823e-5),
               dm3l=(2.450e-3,2.576e-3), delta_deg=(125.,365.)),
    "IO": dict(s12=(.2893,.3295), s13=(.02093,.02441),
               s23=(.439,.584), dm21=(7.236e-5,7.822e-5),
               dm3l=(-2.547e-3,-2.421e-3), delta_deg=(203.,335.)),
}


def inspect_native(path, ordering):
    """Native loop Takagi + charged lepton rotation; one-dimensional screen."""
    source = json.loads(path.read_text())
    snap = source["snapshots"][source["spectrum_snapshot_index"]]
    native_mass = complex_matrix(snap["mat1"])*1e9
    u_native = complex_matrix(snap["UV_1L"])
    el = complex_matrix(snap["ZEL"])
    sorted_mass = np.asarray(snap["MFv_1L"], float)*1e9
    if np.max(abs(native_mass-native_mass.T)) > 1e-11 * max(np.max(abs(native_mass)),1e-20):
        raise ValueError("Native Majorana matrix is not symmetric")
    diagonal = u_native.conj()@native_mass@u_native.conj().T
    if np.linalg.norm(diagonal-np.diag(sorted_mass)) > 1e-9*np.linalg.norm(native_mass):
        raise ValueError("Native Takagi rotation does not diagonalize Majorana mass")
    if any(sorted_mass[i] > sorted_mass[i+1] for i in range(2)):
        raise ValueError("Native neutrino masses are not ascending")
    # U_mass = ZEL U_gauge, U_gauge = UV_1L*^T.
    u = el @ u_native.conj().T
    indices = [0,1,2] if ordering=="NO" else [1,2,0]
    m = sorted_mass[indices]
    u = u[:,indices]
    s13 = float(abs(u[0,2])**2)
    if not 0 < s13 < 1:
        raise ValueError("PMNS reactor angle undefined")
    s12 = float(abs(u[0,1])**2/(1-s13))
    s23 = float(abs(u[1,2])**2/(1-s13))
    a,b,c = math.sqrt(s12),math.sqrt(s13),math.sqrt(s23)
    ca,cb,cc = math.sqrt(1-s12),math.sqrt(1-s13),math.sqrt(1-s23)
    j=float(np.imag(u[0,0]*u[1,1]*u[0,1].conjugate()*u[1,0].conjugate()))
    denom=ca*cc*cb**2*a*c*b
    if denom < 1e-12:
        raise ValueError("Unresolved native Dirac CP phase")
    sind=j/denom
    cosd=(abs(u[1,0])**2-a*a*cc*cc-ca*ca*c*c*b*b)/(2*a*ca*c*cc*b)
    delta = math.degrees(math.atan2(sind,cosd))%360
    values=dict(s12=s12,s13=s13,s23=s23,
                dm21=float(m[1]**2-m[0]**2),
                dm3l=float(m[2]**2-m[0 if ordering=="NO" else 1]**2),
                delta_deg=delta)
    allowed=THREE_SIGMA[ordering]
    checks={k: (lo<=value<=hi if k!="delta_deg" else
                 any(lo<=value+shift<=hi for shift in (0.,360.,-360.)))
            for k,value in values.items() for lo,hi in [allowed[k]]}
    return dict(status="NATIVE_3SIGMA_SCREEN_PASS" if all(checks.values()) else
                "NATIVE_3SIGMA_SCREEN_FAIL",ordering=ordering,
                parameters=values,three_sigma_ranges=allowed,checks=checks,
                masses_eV=m.tolist(),source=SOURCE,
                caveat="One-dimensional projected 3-sigma ranges, not a NuFIT joint likelihood or full Zee acceptance")


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
    ap.add_argument("--native-result",type=Path,help="Existing instrumented native result.json")
    ap.add_argument("--checked-card",type=Path,help="Actual native input card to verify against native result hash")
    args=ap.parse_args()
    if bool(args.reference_card)!=bool(args.write_card):
        ap.error("--reference-card and --write-card must be used together")
    report=fit(args.snapshot,args.ordering,args.mlight_ev,args.alpha21,args.alpha31)
    if args.checked_card and not args.native_result:
        ap.error("--checked-card requires --native-result")
    if args.native_result:
        screen=inspect_native(args.native_result,args.ordering)
        report["native_screen"]=screen
        verified=False
        if args.checked_card:
            native_header=json.loads(args.native_result.read_text())
            verified=(sha(args.checked_card)==native_header["input_sha256"])
            if not verified:
                raise ValueError("Native result input SHA-256 does not match checked card")
            report["native_card_sha256"]=sha(args.checked_card)
        if screen["status"]=="NATIVE_3SIGMA_SCREEN_PASS":
            report["status"]=("NATIVE_3SIGMA_PASS_CARD_BOUND" if verified
                              else "NATIVE_3SIGMA_PASS_CARD_UNVERIFIED")
        else:
            report["status"]="NATIVE_3SIGMA_FAIL"
        report["native_revalidation_required"]=not (verified and screen["status"]=="NATIVE_3SIGMA_SCREEN_PASS")
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
