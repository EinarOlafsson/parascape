"""Bounded local release checks; upload is handled by the main-branch workflow."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys
import tempfile

ROOT = Path.cwd()
sys.path.insert(0,str(ROOT))
from scripts.notebook_runner import ExecutedNotebook

OUT = ROOT/'results/release_0_55_0_2026_10_08_v2'
SCRATCH = Path(tempfile.mkdtemp(prefix='starplast-release-055-'))
ENV = dict(os.environ,QT_QPA_PLATFORM='offscreen',PYQTGRAPH_QT_LIB='PyQt6',
           OPENBLAS_NUM_THREADS='2',OMP_NUM_THREADS='2',NUMBA_NUM_THREADS='2',
           STARPLAST_STATE=str(SCRATCH/'state'))
def run(name,args,cwd=ROOT):
    result=subprocess.run(args,cwd=cwd,env=ENV,stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT,text=True)
    (OUT/(name+'.log')).write_text(result.stdout)
    if result.returncode:
        raise RuntimeError(f'{name} exited {result.returncode}; see preserved log')
    return {'step':name,'status':'passed'}

nb=ExecutedNotebook('Starplast 0.55.0 local release verification')
nb.ns.update(ROOT=ROOT,OUT=OUT,SCRATCH=SCRATCH,run=run,python=sys.executable)
nb.md('User-authorized release despite the unfinished action list. The broader goal remains paused. Execute only explicit release tests; restricted untracked activity-truth files are not collected. Initial release check failed on two missing public API docstrings; the original packet is preserved. Added both docstrings without changing algorithms. Build tools use an isolated temporary environment. CI and a clean installed-wheel check remain separate publication gates.')
try:
    nb.code("version=run('version',[python,'scripts/release.py','check'])",'version')
    nb.code("tests=['test_release.py','test_wheel.py','test_packaging.py','test_readme.py','test_docstrings.py','test_paths.py','test_app_smoke.py','test_app_polish.py','test_discover_cli.py','test_importer.py']",
            "checked=run('tests',[python,'-m','pytest',*[str(ROOT/'tests'/name) for name in tests],'-q'])",'checked')
    nb.code("tools=SCRATCH/'tools'", "run('tools-env',[python,'-m','venv','--system-site-packages',str(tools)])",
            "toolpython=str(tools/'bin/python')", "run('tools-install',[toolpython,'-m','pip','install','pdoc>=14','markdown2>=2.4','twine>=6'])")
    nb.code("docs=run('docs',[toolpython,'scripts/build_docs.py'])",'docs')
    nb.code("result={'version':'0.55.0','scratch':str(SCRATCH),'release_tests':'passed','documentation':'built','broader_goal':'paused','publication':'not performed by this script'}",
            "(OUT/'local_checks.json').write_text(json.dumps(result,indent=2)+'\\n')",'result')
except Exception as exc:
    (OUT/'failure.json').write_text(json.dumps({'error':type(exc).__name__,'detail':str(exc)})+'\n')
    nb.md('Failed: '+type(exc).__name__+': '+str(exc))
    raise
finally:
    nb.write(str(OUT/'executed.ipynb'))
    def sha(path):
        with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
    (OUT/'SHA256SUMS.txt').write_text(''.join(f'{sha(p)}  {p.relative_to(OUT)}\n'
        for p in sorted(OUT.rglob('*')) if p.is_file() and p.name!='SHA256SUMS.txt'))
print(json.dumps(nb.ns['result']))
