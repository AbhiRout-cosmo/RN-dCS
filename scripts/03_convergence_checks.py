"""App. E: convergence of beta*d(omega) at high charge and high l under changes of contour angle, box size and resolution
(Q/M=0.95, 0.98 at l=2; Q/M=0.4, 0.98 at l=6).
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from coef import *
import json
print('--- convergence spot checks (Richardson spread also shown)')
for (l,Q,wg) in [(2,0.95,0.4216928-0.0866587j),(2,0.98,0.4273073-0.0849811j),(6,0.4,1.2310040-0.0956899j),(6,0.98,1.4555508-0.0900726j)]:
    for kw in [dict(th=0.6,S=90,N=400),dict(th=0.5,S=130,N=560),dict(th=0.75,S=110,N=500)]:
        wGR,c,e=coef(1.,Q,l,wg,**kw)
        print(f'l={l} Q={Q} {kw}: wGR={wGR:.8f} beta*dw={c:.5f} spread={e:.1e}',flush=True)
