#!/usr/bin/env bash
# Native Zee pilot: official HiggsTools only, no ScannerS/2HDMC.
set -euo pipefail
PREFIX="$HOME/zee-higgstools"
if [[ "$#" -gt 0 ]]; then PREFIX="$1"; fi
mkdir -p "$PREFIX"

clone_pinned() {
    local url="$1" tag="$2" name="$3"
    local dir="$PREFIX/$name"
    if [[ ! -d "$dir/.git" ]]; then
        git clone --depth 1 --branch "$tag" "$url" "$dir"
    fi
    local current wanted
    current="$(git -C "$dir" rev-parse HEAD)"
    wanted="$(git -C "$dir" rev-parse "$tag^{commit}")"
    if [[ "$current" != "$wanted" ]]; then
        echo "ERROR: $name HEAD differs from pinned tag $tag" >&2
        return 1
    fi
    echo "PINNED $name $tag $current"
}

clone_pinned https://gitlab.com/higgsbounds/higgstools.git v1.3 higgstools
clone_pinned https://gitlab.com/higgsbounds/hbdataset.git v1.7 hbdataset
clone_pinned https://gitlab.com/higgsbounds/hsdataset.git v1.1 hsdataset

if [[ ! -x "$PREFIX/venv/bin/python" ]]; then
    python3 -m venv "$PREFIX/venv"
fi
VENV_PY="$PREFIX/venv/bin/python"
"$VENV_PY" -m pip install --upgrade pip setuptools wheel
"$VENV_PY" -m pip install "pylha==0.4" numpy scipy
# Avoid overloading the 24-GiB host.
CMAKE_BUILD_PARALLEL_LEVEL=2 MAKEFLAGS=-j2 "$VENV_PY" -m pip install "$PREFIX/higgstools"
"$VENV_PY" - <<'PY'
import Higgs
import Higgs.predictions
import Higgs.bounds
import Higgs.signals
import Higgs.tools.Input
print("PASS: official HiggsTools Python interface")
PY
echo "READY: $PREFIX/venv/bin/python"
echo "Datasets: $PREFIX/hbdataset and $PREFIX/hsdataset"
