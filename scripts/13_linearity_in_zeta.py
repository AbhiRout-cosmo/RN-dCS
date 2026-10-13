"""App. E: beta*dw at seven beta in [5e3,3.2e5]; fit c1+c2/beta+c3/beta^2 -> c1 equals Richardson value to 1e-6.
Writes data/linearity_in_zeta.json.  ~15 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
import branch as BR
C={'A':(2,0.0,0.37367168-0.08896231j,dict()),
   'B':(2,0.6,0.38621747-0.08981368j,dict()),
   'C':(2,0.95,0.42169276-0.08665867j,dict(S=90.,N=400)),
   'D':(6,0.6,1.26604130-0.09617508j,dict()),
   'E':(2,0.3,0.46986497-0.09582639j,dict())}     # E: EM-led branch, off resonance
out={'beta':{}}
t0=time.time()
for k in ('A','B','C','D','E'):
    l,Q,g,base=C[k]; kw=dict(th=0.6,S=60.,N=300); kw.update(base)
    w0=near(evals(op_full(1.,Q,l,1,coupled=False,**kw)),g); rows=[]
    for B in (5e3,1e4,2e4,4e4,8e4,1.6e5,3.2e5):
        w=near(evals(op_full(1.,Q,l,B,**kw)),w0); rows.append([B,B*(w-w0)])
        print(k,'beta',B,f'{B*(w-w0):.7f}',flush=True)
    b=np.array([x[0] for x in rows]); y=np.array([x[1] for x in rows])
    A=np.vstack([np.ones_like(b),1/b,1/b**2]).T.astype(complex); sol,*_=np.linalg.lstsq(A,y,rcond=None)
    out['beta'][k]=dict(l=l,Q=Q,w0=w0,rows=rows,fit=sol,maxres=float(np.max(abs(A@sol-y))))
    print(k,'fit c1=',sol[0],' c2=',sol[1],' maxres',np.max(abs(A@sol-y)),flush=True); save('linearity_in_zeta',out)
save('linearity_in_zeta',out)
