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
  git clone https://gitea.psi.ch/ltpth/HDECAY.git hdecay-source
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
for flag,value in (("2HDM","0"),("COUPVAR","0"),("SM4","0"),("HIGGS","5")):
    pattern=rf"(?m)^(\s*{flag}\s*=\s*)[^\n]*"
    s,n=re.subn(pattern,lambda m:m.group(1)+value,s,count=1)
    if n!=1:
        raise RuntimeError("Missing HDECAY input flag "+flag)
patt=r"(?m)^(\s*MHL\s*=\s*)[^\n]*"
s,n=re.subn(patt,lambda m:m.group(1)+os.environ["HDECAY_H_MASS"]+"D0",s)
if n<1:
    raise RuntimeError("Missing HDECAY SM mass input")
open(p,"w").write(s)
PY
  (
     cd "$RUN_DIR"
     ./run > run.stdout.log 2>&1
  )
  test -s "$RUN_DIR/br.sm1"
  test -s "$RUN_DIR/br.sm2"
  echo "HDECAY_SM_COMPLETED $ORDER $MASS"
  head -n 6 "$RUN_DIR/br.sm2"
done
python3 "$SCRIPT_DIR/../ZeeHiggs/hdecay_sm_reference.py" \
  --source "$WORK_ROOT/hdecay-source" \
  --NO "$WORK_ROOT/hdecay-sm-NO" \
  --IO "$WORK_ROOT/hdecay-sm-IO" \
  --output "$WORK_ROOT/hdecay-sm-reference.json"
