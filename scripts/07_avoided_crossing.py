"""Sec. IX C: avoided crossing of the EM-led and pseudoscalar modes near Q_c (l=2,3): the pair splitting scales as beta^-1/2
(split*sqrt(beta) ~ 0.47 for l=2, 0.56 for l=3), i.e. O(sqrt(zeta)).  Needs data/spectra.json (script 06).
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from coef import *
import json
from scipy.optimize import minimize_scalar
sp=json.load(open(DATA/'spectra.json'))
out={}
for l in [2,3]:
    q=np.array(sp[str(l)]['Q']); wE=np.array([complex(*x) for x in sp[str(l)]['em']]); wS=np.array([complex(*x) for x in sp[str(l)]['sc']])
    i=np.argmin(abs(wE-wS)); 
    # refine Q_c by minimizing |w_EM - w_Theta|
    def gap(Q):
        H=operator(1.,Q,l,1e300,coupled=False,th=0.6,S=60,N=300); ev=np.sqrt(sla.eigvals(H)); ev=np.where(ev.real<0,-ev,ev)
        idx=np.argsort(abs(ev-(wE[i]+wS[i])/2))[:2]; a,b=ev[idx]; return abs(a-b),a,b
    r=minimize_scalar(lambda Q: gap(Q)[0],bounds=(q[i]-0.015,q[i]+0.015),method='bounded',options=dict(xatol=1e-7))
    Qc=r.x; g,a,b=gap(Qc)
    print(f'l={l}: Q_c={Qc:.6f} min|w_EM-w_Theta|={g:.2e}  w_EM={a:.6f} w_Theta={b:.6f}',flush=True)
    rows=[]
    for beta in [1e2,1e3,1e4,1e5,1e6,1e7,1e8]:
        H=operator(1.,Qc,l,beta,th=0.6,S=60,N=300); ev=np.sqrt(sla.eigvals(H)); ev=np.where(ev.real<0,-ev,ev)
        # the two modes emerging from the pair (a,b)
        cand=ev[np.argsort(abs(ev-(a+b)/2))[:2]]
        d=cand-(a+b)/2
        rows.append([beta,[complex(x) for x in cand]])
        print(f'   beta={beta:.0e}: pair - mean = {d[0]:.3e}, {d[1]:.3e};  |split|*sqrt(beta)={abs(cand[0]-cand[1])*np.sqrt(beta):.4f}',flush=True)
    out[l]=dict(Qc=Qc,rows=[[r_[0],[[z.real,z.imag] for z in r_[1]]] for r_ in rows])
json.dump(out,open(DATA/'avoided.json','w'))
