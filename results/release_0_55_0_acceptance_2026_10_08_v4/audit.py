"""Verify fresh distributions and an installed wheel outside the source tree."""
from pathlib import Path
import hashlib,json,os,shutil,subprocess,sys,zipfile
ROOT=Path.cwd();sys.path.insert(0,str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
OUT=ROOT/'results/release_0_55_0_acceptance_2026_10_08_v4'
SCRATCH=Path('/tmp/starplast-release-055-95v7eqjm')
ENV=dict(os.environ,QT_QPA_PLATFORM='offscreen',OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2')
def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def run(name,args,cwd=ROOT):
    result=subprocess.run(args,cwd=cwd,env=ENV,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    (OUT/(name+'.log')).write_text(result.stdout)
    if result.returncode:raise RuntimeError(f'{name} exited {result.returncode}')
    return {'step':name,'status':'passed'}
nb=ExecutedNotebook('Starplast 0.55.0 clean wheel acceptance')
nb.ns.update(ROOT=ROOT,OUT=OUT,SCRATCH=SCRATCH,run=run,sha=sha,python=sys.executable,json=json)
nb.md('613 release tests and documentation build passed in local packet V2. Its final reporting cell lacked json in the notebook namespace; preserve that failure and complete the receipt here without repeating passed checks. Initial packet retains the two missing-docstring failures, corrected without algorithm changes. Export only tracked package/build inputs; restricted untracked files are not included. The initial builder cell contained an interpreted null separator; corrected to chr(0), with its original failure preserved. The local build/ directory was mistaken for the build distribution; neither interpreter had its CLI. Install the real build tool only in the temporary tools environment; retain V2/V3 failures. Explicitly reuse parent runtime dependencies in the isolated wheel environment. Verify installed resources and new functional cards; no fit or biological admission.')
try:
    nb.code('from pathlib import Path','import json,shutil,subprocess,zipfile',
        "prior=ROOT/'results/release_0_55_0_2026_10_08_v2'",
        "assert '613 passed' in (prior/'tests.log').read_text()",
        "assert 'Documentation:' in (prior/'docs.log').read_text()",
        "assert json.loads((prior/'failure.json').read_text())['error']=='NameError'",
        "source=SCRATCH/'source-v4';source.mkdir(exist_ok=False)",
        "files=subprocess.check_output(['git','ls-files','-z','pyproject.toml','README.md','LICENSE','starplast'],cwd=ROOT).decode().split(chr(0))",
        'inputs={}',
        "for name in filter(None,files):\n path=ROOT/name;target=source/name;target.parent.mkdir(parents=True,exist_ok=True)\n shutil.copy2(path,target);inputs[name]=sha(path)",
        "(OUT/'input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\\n')",
        "dist=SCRATCH/'dist-v4';toolpython=str(SCRATCH/'tools/bin/python')",
        "run('build-tools',[toolpython,'-m','pip','install','build>=1.2','setuptools>=77','wheel'])",
        "run('build',[toolpython,'-m','build','--no-isolation','--outdir',str(dist)],cwd=source)",
        "run('wheel-check',[python,'scripts/check_wheel.py',str(dist)])",
        "run('metadata',[toolpython,'-m','twine','check','--strict',*[str(p) for p in sorted(dist.iterdir())]])")
    nb.code("wheel=next(dist.glob('*.whl'))",
        "required=['starplast/functional_pair_class.py','starplast/group_ratios.py','starplast/data/functional_results.json','starplast/data/functional_agreement.json','starplast/data/functional_agreement_uncertainty.json']",
        "with zipfile.ZipFile(wheel) as archive:\n for name in required:assert archive.read(name)==(ROOT/name).read_bytes(),name",
        "installed=SCRATCH/'installed'",
        "run('installed-env',[python,'-m','venv','--system-site-packages',str(installed)])",
        "import sys,sysconfig",
        "site=installed/'lib'/f'python{sys.version_info.major}.{sys.version_info.minor}'/'site-packages'",
        "(site/'release_dependencies.pth').write_text(sysconfig.get_path('purelib')+'\\n')",
        "installedpython=str(installed/'bin/python')",
        "run('install',[installedpython,'-m','pip','install','--no-deps','--force-reinstall',str(wheel)])",
        "code=\"from pathlib import Path; import starplast; from starplast import paths,functional_results as F,functional_agreement as A,functional_pair_class as C; assert starplast.__version__=='0.55.0'; assert '/installed/lib/' in starplast.__file__; assert paths.check()[0]; b,reason=F.shipped('Pf'); assert not reason and len(b)==2; pair=A.find_pair(b[0],b); r=A.compare(*pair); result,view=C.build(pair,r,'major:2'); assert len(result['rows'])==152; assert len(C.addresses(pair))==19; print(starplast.__version__,starplast.__file__,len(result['rows']))\"",
        "run('installed-smoke',[installedpython,'-I','-c',code],cwd=SCRATCH)",
        "run('installed-cli',[str(installed/'bin/starplast-discover'),'--help'],cwd=SCRATCH)",
        'for name,digest in inputs.items():assert sha(ROOT/name)==digest,name',
        "result={'version':'0.55.0','release_checks':613,'documentation':'built','tracked_build_inputs':len(inputs),'wheel_and_sdist':{p.name:{'sha256':sha(p),'size':p.stat().st_size} for p in dist.iterdir()},'installed_wheel':'verified outside source','new_functional_resources':'exact original bytes','broader_goal':'paused','publication':'pending CI and upload'}",
        "(OUT/'acceptance.json').write_text(json.dumps(result,indent=2)+'\\n')",'result')
except Exception as exc:
    (OUT/'failure.json').write_text(json.dumps({'error':type(exc).__name__,'detail':str(exc)})+'\n')
    nb.md('Failed: '+type(exc).__name__+': '+str(exc));raise
finally:
    nb.write(str(OUT/'executed.ipynb'))
    (OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n' for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(json.dumps(nb.ns['result']))
