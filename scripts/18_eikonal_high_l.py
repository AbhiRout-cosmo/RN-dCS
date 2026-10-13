"""Sec. IX A, Fig. 2: grav-led coefficient for l=12..30 at Q=0,0.3,0.6 (reduced ECS, beta=1e7) for the ratio to the
l^3 and l^2/Q laws.  Writes data/eikonal_high_l.json.  ~30 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
import branch as BR
out={}
for Q in (0.0,0.3,0.6):
    rows=[]
    c=json.load(open(DATA/'c3_sweeps.json'))['lscal'][str(Q)]
    w12=complex(c[-1][1],c[-1][2])
    for l in (12,16,20,24,30):
        g=w12.real*(l+0.5)/12.5+1j*w12.imag
        w0=BR.gr(1.,max(Q,1e-12),l,1,g,N=360)
        w0=BR.gr(1.,max(Q,1e-12),l,1,w0,N=360)
        _,c=BR.coef_b(1.,max(Q,1e-12),l,1,w0,B=1e7,N=360)
        rows.append([l,w0,c]); print(f'Q={Q} l={l} w0={w0:.6f} beta*dw={c:.4f}',flush=True)
        out[Q]=rows; save('eikonal_high_l',out)

