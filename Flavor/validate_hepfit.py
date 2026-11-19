"""Attest a rebuilt HEPfit binary using the unchanged historical source/run suite.

A changed binary hash requires source verification, live build artifacts and
all 37 historical input/output comparisons. A claimed PASS alone is not proof.
"""
import argparse
import json
import math
from pathlib import Path
import runpy

from Flavor.bsgamma import HEPFIT_COMMIT, MATCHER_SHA, parse_driver_output, sha

# Same absolute operator-map and implementation-identity bounds as provenance.py.
OPERATOR_TOLERANCE = 2e-14
IMPLEMENTATION_TOLERANCE = 2e-13


def compare_runs(historical, rebuilt):
    if len(historical['runs'])!=37 or set(rebuilt['runs'])!=set(historical['runs']):
        raise ValueError('Restored runtime must reproduce all 37 historical runs')
    residuals = {}
    for name,old in historical['runs'].items():
        new=rebuilt['runs'][name]
        if new['input']!=old['input'] or new['mode']!=old['mode']:
            raise ValueError('Historical input or mode changed: '+name)
        if set(new['values'])!=set(old['values']):
            raise ValueError('Historical observable vector changed: '+name)
        for key,expected in old['values'].items():
            if not math.isfinite(new['values'][key]) or not math.isfinite(expected):
                raise ValueError('Nonfinite historical comparison: '+name+'/'+key)
            delta=abs(new['values'][key]-expected)
            tolerance=(OPERATOR_TOLERANCE if key.startswith(('high_','low_'))
                       else IMPLEMENTATION_TOLERANCE)
            if delta>tolerance:
                raise ValueError(f'Historical numerical mismatch: {name}/{key}: {delta}')
        residuals[name]=max(abs(new['values'][k]-v) for k,v in old['values'].items())
    return residuals


def validate(upstream,driver,source,build_path,execution_path):
    base=upstream/'research/full_constraints/bsgamma_hepfit'
    historical_path=base/'execution/HEPFIT_EXECUTION.json'
    lock_path=base/'SOURCE_LOCK.json'
    lock=json.loads(lock_path.read_text())
    build=json.loads(build_path.read_text())
    if lock['HEPfit_commit']!=HEPFIT_COMMIT:
        raise ValueError('Wrong pinned HEPfit commit')
    if sha(upstream/'scripts/zee_module5_bmu_nlo.py')!=MATCHER_SHA:
        raise ValueError('BMU matching source changed')
    provenance=runpy.run_path(str(base/'provenance.py'))
    count=provenance['verify_source'](source)
    if build.get('verified_source_files')!=count:
        raise ValueError('Build source verification count mismatch')
    if sha(base/'initialize_matching.patch')!=lock['initialization_patch_sha256']:
        raise ValueError('Initialization patch changed')
    if build['artifacts'].get(str(driver.resolve()))!=sha(driver):
        raise ValueError('Driver does not match restored build receipt')
    for name,digest in build['artifacts'].items():
        if sha(Path(name))!=digest:
            raise ValueError('Build artifact changed: '+name)
    caches=[p for p in build['artifacts'] if p.endswith('/build/CMakeCache.txt')]
    if len(caches)!=1:
        raise ValueError('Missing unique HEPfit CMake build cache')
    cache=Path(caches[0])
    homes=[line.split('=',1)[1] for line in cache.read_text().splitlines()
           if line.startswith('CMAKE_HOME_DIRECTORY:INTERNAL=')]
    if homes!=[str(source.resolve())]:
        raise ValueError('CMake build is not tied to supplied pinned source')
    historical=provenance['verify_execution'](historical_path)
    rebuilt=provenance['verify_execution'](execution_path)
    if (rebuilt['HEPfit_commit']!=HEPFIT_COMMIT or rebuilt['binary_sha256']!=sha(driver)
        or rebuilt['driver_sha256']!=sha(base/'driver.cpp')):
        raise ValueError('Restored execution does not bind the actual source and binary')
    for name,r in rebuilt['runs'].items():
        raw=parse_driver_output((execution_path.parent/r['stdout']).read_text(),
                                matching_scale=r['input']['muw'])
        if raw!=r['values']:
            raise ValueError('Execution values differ from actual stdout: '+name)
        actual=dict((k,float(v)) for k,v in (line.split() for line in
            (execution_path.parent/(name+'.txt')).read_text().splitlines()))
        if actual!=r['input']:
            raise ValueError('Execution input differs from actual card: '+name)
    residuals=compare_runs(historical,rebuilt)
    return dict(
        schema_version=1,status='RESTORED_HEPFIT_NUMERICALLY_EQUIVALENT',
        HEPfit_commit=HEPFIT_COMMIT,binary_sha256=sha(driver),
        driver_sha256=sha(base/'driver.cpp'),matching_sha256=MATCHER_SHA,
        source_lock_sha256=sha(lock_path),
        historical_execution_sha256=sha(historical_path),
        build_receipt=str(build_path.resolve()),build_receipt_sha256=sha(build_path),
        execution_receipt=str(execution_path.resolve()),execution_receipt_sha256=sha(execution_path),
        verified_source_files=count,validated_runs=len(residuals),
        maximum_output_absolute_residual=max(residuals.values()),
        per_run_maximum_residual=residuals,
        comparison_tolerances=dict(operator_absolute=OPERATOR_TOLERANCE,
                                   other_absolute=IMPLEMENTATION_TOLERANCE),
        original_central_BR={name:rebuilt['runs'][name]['values']['BR_NLO'] for name in
            ('SM_reference','d7_NO_neutral_near_heavy_07_central','F13_heavy_0001330_central')},
        physical_bsgamma='UNKNOWN',theory_sigma=None,experimental_exclusion=None)


def verify_attestation(path,upstream,driver,source):
    record=json.loads(path.read_text())
    if (record.get('schema_version')!=1 or
        record.get('status')!='RESTORED_HEPFIT_NUMERICALLY_EQUIVALENT'):
        raise ValueError('Unrecognized restored HEPfit attestation')
    build_path=Path(record['build_receipt'])
    execution_path=Path(record['execution_receipt'])
    if sha(build_path)!=record['build_receipt_sha256'] or sha(execution_path)!=record['execution_receipt_sha256']:
        raise ValueError('Restored backend receipt changed')
    verified=validate(upstream,driver,source,build_path,execution_path)
    if verified!=record:
        raise ValueError('Restored backend attestation does not match actual evidence')
    return verified


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('upstream','driver','source','build-receipt','execution-receipt','output'):
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    result=validate(a.upstream.resolve(),a.driver.resolve(),a.source.resolve(),
                    a.build_receipt.resolve(),a.execution_receipt.resolve())
    with a.output.open('x') as f:
        json.dump(result,f,indent=2,allow_nan=False)
        f.write('\n')
    print(json.dumps(result['original_central_BR'],indent=2))
    print(result['status'],result['validated_runs'])


if __name__=='__main__':
    main()
