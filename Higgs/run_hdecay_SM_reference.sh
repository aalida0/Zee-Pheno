#!/usr/bin/env bash
# Compile and execute public HDECAY 6.62 in SM mode (not Zee physics).
# Git commit pinned, without modifying its physics source.
set -euo pipefail
WORK_ROOT="$1"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$WORK_ROOT"
cd "$WORK_ROOT"
HDECAY_COMMIT=4884213286e1edbdf3b2bb9fe30095b29e336dff
if [[ ! -d hdecay-source/.git ]]; then
  for attempt in 1 2 3; do
    rm -rf hdecay-source
    if git -c http.lowSpeedTime=15 clone https://gitea.psi.ch/ltpth/HDECAY.git hdecay-source; then
      break
    fi
    echo "PSI HDECAY clone attempt $attempt failed (upstream availability); retrying" >&2
    sleep 4
  done
fi
if [[ ! -d hdecay-source/.git ]]; then
  echo "BLOCKED: official PSI HDECAY source unavailable; no HDECAY numerical validation performed" >&2
  exit 69
fi
git -C hdecay-source checkout --detach "$HDECAY_COMMIT"
test "$(git -C hdecay-source rev-parse HEAD)" = "$HDECAY_COMMIT"
gfortran --version | head -n 1
make -C hdecay-source/zzz -B hdecay
test -x hdecay-source/zzz/run
echo "HDECAY_SOURCE_SHA256 $(sha256sum hdecay-source/zzz/hdecay.f | cut -d' ' -f1)"

for row in NO:125.250000 IO:125.250135; do
  ORDER="$(echo "$row" | cut -d: -f1)"
  MASS="$(echo "$row" | cut -d: -f2)"
  RUN_DIR="$WORK_ROOT/hdecay-sm-$ORDER"
  mkdir -p "$RUN_DIR"
  cp hdecay-source/zzz/run "$RUN_DIR/run"
  cp hdecay-source/zzz/hdecay.in "$RUN_DIR/hdecay.in"
  HDECAY_H_MASS="$MASS" HDECAY_INPUT="$RUN_DIR/hdecay.in" python3 - <<'PY'
import os,re
p=os.environ["HDECAY_INPUT"]
s=open(p).read()
for flag,value in (("2HDM","0"),("COUPVAR","0"),("SM4","0"),("HIGGS","0"),("SLHAOUT","0"),("NMA","1")):
    pattern=rf"(?m)^(\s*{flag}\s*=\s*)[^\n]*"
    s,n=re.subn(pattern,lambda m:m.group(1)+value,s,count=1)
    if n!=1:
        raise RuntimeError("Missing HDECAY input flag "+flag)
# HDECAY manual: MABEG/MAEND set the *SM* mass range (not MHL,
# which is used by the hMSSM/2HDM input).
for masskey in ("MABEG","MAEND"):
    patt=rf"(?m)^(\\s*{masskey}\\s*=\\s*)[^\\n]*"
    s,n=re.subn(patt,lambda m:m.group(1)+os.environ["HDECAY_H_MASS"]+"D0",s,count=1)
    if n!=1:
        raise RuntimeError("Missing HDECAY SM mass key "+masskey)
patt=r"(?m)^(\\s*MHL\\s*=\\s*)[^\\n]*"
s,n=re.subn(patt,lambda m:m.group(1)+os.environ["HDECAY_H_MASS"]+"D0",s)
if n<1:
    raise RuntimeError("Missing HDECAY SM mass input")
open(p,"w").write(s)
PY
  if ! (
     cd "$RUN_DIR"
     ./run > run.stdout.log 2>&1
  ); then
     echo "HDECAY run returned nonzero for $ORDER" >&2
     tail -n 45 "$RUN_DIR/run.stdout.log" >&2
     ls -lah "$RUN_DIR" >&2
     exit 1
  fi
  if [[ ! -s "$RUN_DIR/br.sm1" || ! -s "$RUN_DIR/br.sm2" ]]; then
     echo "HDECAY SM result files missing for $ORDER" >&2
     tail -n 45 "$RUN_DIR/run.stdout.log" >&2
     ls -lah "$RUN_DIR" >&2
     exit 1
  fi
  echo "HDECAY_SM_COMPLETED $ORDER $MASS"
  head -n 6 "$RUN_DIR/br.sm2"
done
python3 "$SCRIPT_DIR/../ZeeHiggs/hdecay_sm_reference.py" \
  --source "$WORK_ROOT/hdecay-source" \
  --NO "$WORK_ROOT/hdecay-sm-NO" \
  --IO "$WORK_ROOT/hdecay-sm-IO" \
  --output "$WORK_ROOT/hdecay-sm-reference.json"
