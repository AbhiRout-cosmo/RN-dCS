"""Sec. IX C, Fig. 3(c): splitting of the coupled EM-led/pseudoscalar pair at the exact degeneracy
Q_c = 2M sqrt(L+1)/(L+2) (symbolic/07_exact_Qc.py), l=2,3, beta = 1e3..1e10.
Result: |w+ - w-| sqrt(beta) = 0.47995 (l=2), 0.57617 (l=3), constant to 5 digits => O(sqrt(zeta)).
Writes data/avoided.json.  (v1.0.0 used a mis-tracked minimum near Q=0.655; superseded.)  ~3 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
out={}
for l,g in ((2,0.52618687-0.09850443j),(3,0.70858954-0.09775797j)):
    L=l*(l+1); Qc=2*np.sqrt(L+1)/(L+2); rows=[]
    for B in (1e3,1e4,1e5,1e6,1e7,1e8,1e10):
        ev=evals(op_full(1.,Qc,l,B,N=300)); idx=np.argsort(abs(ev-g))[:2]; a,b=ev[idx]
        rows.append([B,a,b,abs(a-b)*np.sqrt(B),(a+b)/2-g])
        print(f'l={l} Qc={Qc:.7f} beta={B:.0e}: |split|*sqrt(beta)={abs(a-b)*np.sqrt(B):.5f}  mean shift={(a+b)/2-g:.2e}',flush=True)
    out[l]=dict(Qc=Qc,rows=rows)
save('avoided',out)
