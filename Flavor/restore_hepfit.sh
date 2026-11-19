#!/usr/bin/env bash
# Restore the historical source-qualified runtime without modifying ~/zee.
set -euo pipefail
build_root=${1:-/home/amine/zee-hepfit-build}
upstream=${2:-/home/amine/zee}
task_dir="$upstream/research/full_constraints/bsgamma_hepfit"
root_prefix=/home/amine/root
hepfit_src="$build_root/HEPfit-322bd73a2dadd94ce94a6d71cb79d616e569120a"
gslpp_src="$build_root/gslpp-f5b059e058a47f6797823a3846fd2cf9d17ff12c"
rge_src="$build_root/RGESolver-d1eb75559e03ae74a9ed0532760560757ef3beec"
compat_dir="$build_root/compat/usr/lib/x86_64-linux-gnu"
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$build_root"
# Existing archives are verified, never downloaded again or silently replaced.
python3 - "$build_root" "$task_dir/SOURCE_LOCK.json" <<'PY'
import hashlib,json,pathlib,subprocess,sys
root=pathlib.Path(sys.argv[1]); lock=json.load(open(sys.argv[2]))
urls=dict(zip(('source.tar.gz','gslpp.tar.gz','rgesolver.tar.gz','libtbb2.deb'),
              (lock['source_urls'][k] for k in ('HEPfit','gslpp','RGESolver','libtbb2'))))
for name,url in urls.items():
    p=root/name
    if not p.exists():
        subprocess.run(['curl','--fail','--location','--silent','--show-error',
                        '--retry','2','--connect-timeout','30',url,'-o',str(p)],check=True)
    if hashlib.sha256(p.read_bytes()).hexdigest()!=lock['archives'][name]:
        raise SystemExit('Archive SHA256 mismatch: '+str(p))
    print('Verified archive',name,flush=True)
PY
if [ ! -d "$hepfit_src" ]; then
    tar -xzf "$build_root/source.tar.gz" -C "$build_root"
    tar -xzf "$build_root/gslpp.tar.gz" -C "$build_root"
    tar -xzf "$build_root/rgesolver.tar.gz" -C "$build_root"
    cp -a "$gslpp_src/src" "$hepfit_src/gslpp/"
    cp "$gslpp_src/CMakeLists.txt" "$hepfit_src/gslpp/"
    cp -a "$gslpp_src/src" "$rge_src/gslpp/"
    patch -d "$hepfit_src" -p1 --batch -i "$task_dir/initialize_matching.patch"
fi
python3 "$task_dir/provenance.py" --source "$hepfit_src"
# Verify the reused ROOT/GSL artifacts against the historical build receipt.
python3 - "$task_dir/BUILD_RECEIPT.json" "$task_dir" <<'PY'
import hashlib,json,pathlib,sys
receipt=json.load(open(sys.argv[1]))
task=pathlib.Path(sys.argv[2])
execution=json.loads((task/'execution/HEPFIT_EXECUTION.json').read_text())
lock=json.loads((task/'SOURCE_LOCK.json').read_text())
for p,expected in ((task/'driver.cpp',execution['driver_sha256']),
                   (task/'initialize_matching.patch',lock['initialization_patch_sha256'])):
    if hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
        raise SystemExit('Historical driver/patch SHA256 mismatch: '+str(p))
for name in ('/home/amine/root/bin/root-config','/home/amine/root/lib/libMathCore.so',
             '/home/amine/root/lib/libMathMore.so',
             '/home/amine/codes/Herwig/herwigpp/lib/libgsl.so.25.1.0'):
    key=(name if '/root/' in name else
         '/tmp/zee-hepfit/compat/usr/lib/x86_64-linux-gnu/libgsl.so.25.1.0')
    if hashlib.sha256(pathlib.Path(name).read_bytes()).hexdigest()!=receipt['artifacts'][key]:
        raise SystemExit('Historical dependency SHA256 mismatch: '+name)
PY
cmake -S "$rge_src" -B "$build_root/rge-build" \
    -DCMAKE_INSTALL_PREFIX="$build_root/deps" -DROOT_CONFIG_DIR="$root_prefix/bin"
cmake --build "$build_root/rge-build" --target RGESolver -j2
cmake --install "$build_root/rge-build"
cmake -S "$hepfit_src" -B "$build_root/build" -DNOMCMC=ON -DBAT_INSTALL=OFF \
    -DROOT_CONFIG_DIR="$root_prefix/bin" -DRGESOLVER_CONFIG_DIR="$build_root/deps/bin" \
    -DCMAKE_INSTALL_PREFIX="$build_root/install"
cmake --build "$build_root/build" --target StandardModel THDM Flavour LeptonFlavour LoopFunctions -j2
mkdir -p "$build_root/compat"
dpkg-deb -x "$build_root/libtbb2.deb" "$build_root/compat"
cp /home/amine/codes/Herwig/herwigpp/lib/libgsl.so.25.1.0 "$compat_dir/"
ln -sfn libgsl.so.25.1.0 "$compat_dir/libgsl.so.25"
python3 "$task_dir/build_driver.py" --source "$hepfit_src" --build "$build_root/build" \
    --root-config "$root_prefix/bin/root-config" \
    --rgesolver-config "$build_root/deps/bin/rgesolver-config" \
    --compat-lib-dir "$compat_dir" --output "$build_root/driver"
# The original execution suite uses these small independent comparison clients.
python3 - <<'PY'
import hashlib,pathlib
p=pathlib.Path('/home/amine/Downloads/superiso_v5.0.tgz')
if hashlib.sha256(p.read_bytes()).hexdigest()!='7f8400e98ae80deac50af64486b980891e1630d729a0a4660ca309a4d7706b5f':
    raise SystemExit('SuperIso archive SHA256 mismatch')
PY
if [ ! -d "$build_root/superiso/SuperIso_v5.0" ]; then
    mkdir -p "$build_root/superiso"
    tar -xzf /home/amine/Downloads/superiso_v5.0.tgz -C "$build_root/superiso"
    patch -d "$build_root/superiso/SuperIso_v5.0" -p1 --batch \
        -i "$upstream/module5/patches/d9b_superiso_v5.0.patch"
fi
make -C "$build_root/superiso/SuperIso_v5.0" -j1
# Upstream SuperIso ar-member rules require serial make.
gcc -O2 -I "$build_root/superiso/SuperIso_v5.0/src" \
    "$upstream/scripts/zee_d9b_rate_probe.c" "$build_root/superiso/SuperIso_v5.0/src/libsuperiso.a" \
    -lm -o "$build_root/superiso_probe"
gcc -O2 -I "$build_root/superiso/SuperIso_v5.0/src" \
    "$task_dir/superiso_inputs.c" "$build_root/superiso/SuperIso_v5.0/src/libsuperiso.a" \
    -lm -o "$build_root/superiso_inputs"
python3 "$task_dir/capture_build.py" --scratch "$build_root" --root-prefix "$root_prefix" \
    --output "$build_root/BUILD_RECEIPT.json"
printf 'Build complete: %s\n' "$build_root/driver"
