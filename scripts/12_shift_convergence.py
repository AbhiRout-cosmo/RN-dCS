"""App. E: convergence of the shift beta*dw itself under contour angle, box size, resolution and rotation pivot
(l=2 Q=0,0.6,0.95; l=6 Q=0.6; EM-led l=2 Q=0.3).  Writes data/shift_convergence.json.  ~40 min.
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
out={'conv':{}}
t0=time.time()
for k,(l,Q,g,base) in C.items():
    b=dict(th=0.6,S=60.,N=300,rc=0.0); b.update(base)
    vars_=[('base',{})]+[(f'th={x}',dict(th=x)) for x in (0.5,0.7,0.8)]+[('th=0.4,S=90,N=400',dict(th=0.4,S=90.,N=400))]
    vars_+=[(f'S={x}',dict(S=x)) for x in ((45.,80.) if b['S']==60 else (70.,120.))]
    vars_+=[(f'N={x}',dict(N=x)) for x in ((220,260,360) if b['N']==300 else (340,460))]
    vars_+=[(f'pivot={x}',dict(rc=x)) for x in (-2.,2.,5.)]
    res=[]
    for name,ch in vars_:
        kw=dict(b); kw.update(ch)
        w0,c,spr=coef_full(1.,Q,l,g,**kw); res.append([name,w0,c,spr])
        print(k,name,f'{w0:.10f}',f'{c:.7f}',f'{spr:.1e}','%.0fs'%(time.time()-t0),flush=True)
    out['conv'][k]=dict(l=l,Q=Q,rows=res); save('shift_convergence',out)
save('shift_convergence',out)
