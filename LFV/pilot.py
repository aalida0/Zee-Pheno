"""Step 2 pilot: complete *scalar one-loop subset*, NOT complete Zee LFV.

Uses matching-time native Yukawas and charged mixing; neutral h/H/A masses
and SCALARMIX are read from the same SLHA spectrum. No physical LFV pass/exclusion is claimed.
"""
import argparse
import json
import math
import re
from pathlib import Path

import numpy as np

from LFV.core import (BASIS, MISSING, OPERATOR, combine, dipole_ledger,
                      radiative_width)
from Model.audit_neutrino_native import parse_entries, sha
from Neutrino.fit import complex_matrix

# Historic Module3/PDG-2024 illustrative lifetime normalizations, not LFV limits.
ALPHA_ZERO = 0.0072973525692838015
MU_WIDTH_GEV = 2.995983701908041e-19
TAU_WIDTH_GEV = 2.267350867723045e-12


def as_json(z):
    return {"real": float(z.real), "imag": float(z.imag)}



def neutral_alpha(spectrum):
    """Read alpha, allowing the arbitrary sign of the second real eigenstate."""
    rows, active = {}, False
    for line in spectrum.splitlines():
        header = re.match(r"\s*BLOCK\s+(\w+)\b", line, re.I)
        if header:
            active = header[1].upper() == "SCALARMIX"
            continue
        if not active:
            continue
        fields = line.split("#", 1)[0].split()
        if len(fields) < 3:
            continue
        try:
            i, j = map(int, fields[:2])
            value = float(fields[2].replace("D", "E").replace("d", "e"))
        except ValueError as exc:
            raise ValueError("Malformed native SCALARMIX row") from exc
        if (i, j) in rows or not math.isfinite(value):
            raise ValueError("Duplicate or nonfinite SCALARMIX entry")
        rows[i, j] = value
    if set(rows) != {(1,1),(1,2),(2,1),(2,2)}:
        raise ValueError("Native spectrum requires complete 2x2 SCALARMIX")
    c, s = rows[1,1], rows[1,2]
    second_row = min(max(abs(rows[2,1]+s),abs(rows[2,2]-c)),
                     max(abs(rows[2,1]-s),abs(rows[2,2]+c)))
    if max(abs(c*c+s*s-1),second_row)>1e-5:
        raise ValueError("SCALARMIX does not match expected CP-even rotation")
    return math.atan2(s, c)


def evaluate(native_path, card_path, spectrum_path):
    spectrum=spectrum_path.read_text()
    alpha_rad=neutral_alpha(spectrum)
    record = json.loads(native_path.read_text())
    if "snapshots" not in record or "spectrum_snapshot_index" not in record:
        raise ValueError("Expected the complete instrumented native result.json")
    if sha(card_path) != record["input_sha256"]:
        raise ValueError("Native result is not associated with this input card")
    s = record["snapshots"][record["spectrum_snapshot_index"]]
    vd, vu = float(s["vd"]), float(s["vu"])
    v = math.hypot(vd, vu)
    if min(vd,vu,v) <= 0:
        raise ValueError("Invalid EWSB")
    cb, sb = vd/v, vu/v
    beta = math.atan2(vu,vd)
    el, er = complex_matrix(s["ZEL"]), complex_matrix(s["ZER"])
    f = complex_matrix(s["fzee"])
    d = -sb*complex_matrix(s["Ye"]) + cb*complex_matrix(s["Ye2"])
    fm = el.conj()@f@el.conj().T
    dm = er.conj()@d@el.conj().T
    if np.max(abs(fm+fm.T)) > 1e-10:
        raise ValueError("Non-antisymmetric native F")
    ml = np.asarray(s["MFe"],float)
    charged = np.asarray(s["MHm"],float)[1:]
    z = np.asarray(s["ZP"],float)
    # q_2 is the singlet fraction; h_2 the physical doublet fraction.
    sin_gamma=z[1,2]
    cos_gamma=-sb*z[1,0]+cb*z[1,1]
    if abs(sin_gamma**2+cos_gamma**2-1) > 1e-8:
        raise ValueError("Unexpected native charged mixing convention")
    gamma=math.atan2(sin_gamma,cos_gamma)
    mass = dict(parse_entries(spectrum,"MASS"))
    try:
        external = np.array([mass[k] for k in (11,13,15)], float)
    except KeyError as exc:
        raise ValueError(f"Missing charged-lepton pole mass {exc}") from exc
    if (not np.isfinite(external).all() or not np.isfinite(ml).all() or
            np.any(external<=0) or np.any(ml<=0) or
            not np.isfinite(fm).all() or not np.isfinite(dm).all() or
            not np.isfinite(z).all()):
        raise ValueError("Nonfinite or unphysical native LFV input")
    try:
        mh,mH,ma=(mass[k] for k in (25,35,36))
    except KeyError as exc:
        raise ValueError(f"Missing native neutral scalar mass {exc}") from exc
    p = dict(v=v,beta=beta,alpha=alpha_rad,gamma=gamma,
             mh1=mh,mh2=mH,ma=ma,mH1=charged[0],mH2=charged[1])
    if (not np.isfinite([mh,mH,ma,*charged]).all() or
            min(mh,mH,ma,*charged) < 20*ml[-1]):
        raise ValueError("Scalar mass too light for audited one-loop expansion")
    results={}
    for label,i,j,width in (("mu_to_e_gamma",1,0,MU_WIDTH_GEV),
                             ("tau_to_e_gamma",2,0,TAU_WIDTH_GEV),
                             ("tau_to_mu_gamma",2,1,TAU_WIDTH_GEV)):
        parts=dipole_ledger(fm,dm,p,ml,i,j,external=external)
        total=combine(parts)
        br=radiative_width(total,external[i],external[j],ALPHA_ZERO)/width
        results[label]=dict(
            BR_scalar_one_loop_subset=br, C_L_GeV_minus2=as_json(total[0]),
            C_R_GeV_minus2=as_json(total[1]),
            component_CL_CR={k:{"C_L":as_json(a["C_L"]),"C_R":as_json(a["C_R"])}
                             for k,a in parts.items()},
            missing_sector_not_bounded=True)
    return dict(status="SCALAR_ONE_LOOP_SUBSET_ONLY",
                physical_lfv="UNKNOWN", global_zee_acceptance=None,
                native_input_sha256=sha(card_path),native_diagnostic_sha256=sha(native_path),
                neutral_spectrum_sha256=sha(spectrum_path),
                scalar_inputs=p, internal_running_lepton_masses_GeV=ml.tolist(),
                external_pole_lepton_masses_GeV=external.tolist(),
                alpha_em=ALPHA_ZERO,
                total_lepton_widths_GeV={"mu":MU_WIDTH_GEV,"tau":TAU_WIDTH_GEV},
                dipole_operator=OPERATOR,basis=BASIS,processes=results,
                missing_contributions=MISSING,
                note="A subset BR is not an exclusion/pass. Internal running masses/Yukawas and external pole leptons are distinguished, with charged tree masses and neutral pole masses documented as a mixed-order scalar scheme.")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--native",type=Path,required=True,
                    help="Instrumented SPheno result.json for one fitted NO/IO card")
    ap.add_argument("--card",type=Path,required=True,help="Original NO/IO input card")
    ap.add_argument("--spectrum",type=Path,required=True,
                    help="Neutral Higgs pole masses PDG 25,35,36")
    ap.add_argument("--output",type=Path,required=True,help="Must not already exist")
    args=ap.parse_args()
    report=evaluate(args.native,args.card,args.spectrum)
    with args.output.open("x") as out:
        json.dump(report,out,indent=2)
        out.write("\n")
    print(json.dumps({k:v for k,v in report.items() if k not in
                      ("processes","missing_contributions")},indent=2))
    for k,v in report["processes"].items():
        print(k,"BR_scalar_one_loop_subset",f'{v["BR_scalar_one_loop_subset"]:.8e}')


if __name__=="__main__":
    main()
