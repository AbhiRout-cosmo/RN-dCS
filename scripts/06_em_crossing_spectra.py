"""Decoupled RN branch and pseudoscalar spectra vs Q (input to script 07).  NOTE: nearest-eigenvalue tracking can merge the
EM-led and scalar modes beyond the crossing; the paper's Q_c(l) values come from script 05, which uses exact branch identification.
Writes data/spectra.json.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from coef import *
import json
# decoupled spectra: track grav-led, EM-led, scalar with fine Q steps
def spectra(Q,l,guesses,N=300,S=60):
    H=operator(1.,Q,l,1e300,coupled=False,th=0.6,S=S,N=N)
    ev=np.sqrt(sla.eigvals(H)); ev=np.where(ev.real<0,-ev,ev)
    return [ev[np.argmin(abs(ev-g))] for g in guesses]
sch={2:[0.3736717-0.0889623j,0.4575956-0.0950044j,0.4836439-0.0967588j],
     3:[0.5994433-0.0927030j,0.6568989-0.0956162j,0.6753662-0.0964996j],
     4:[0.8091784-0.0941640j,0.8530954-0.0958599j,0.8673656-0.0967549j],
     5:[1.0122953-0.0948705j,1.0479131-0.0959817j,1.0596506-0.0968682j],
     6:[1.2120098-0.0952658j,1.2419962-0.0960512j,1.2519780-0.0969280j]}
out={}
for l in [2,3,4,5,6]:
    g=sch[l]; tab=[]
    # refine scalar Schwarzschild value first
    g=spectra(0.0,l,g)
    for Q in np.round(np.arange(0,0.985,0.01),3):
        N=300 if Q<0.9 else 400; S=60 if Q<0.9 else 90
        g=spectra(Q,l,g,N=N,S=S); tab.append([Q]+[complex(x) for x in g])
    q=np.array([t[0] for t in tab]); wE=np.array([t[2] for t in tab]); wS=np.array([t[3] for t in tab])
    d=wE.real-wS.real; k=np.where(np.sign(d[:-1])!=np.sign(d[1:]))[0]
    Qc=[q[i]-d[i]*(q[i+1]-q[i])/(d[i+1]-d[i]) for i in k]
    print(f'l={l}: scalar(Q=0)={tab[0][3]:.6f}; Re(w_EM)=Re(w_Theta) at Q_c={np.round(Qc,4)}; |dIm| there={[abs(wE[i].imag-wS[i].imag) for i in k]}',flush=True)
    out[l]=dict(Q=q.tolist(),grav=[[x.real,x.imag] for x in [t[1] for t in tab]],em=[[x.real,x.imag] for x in wE],sc=[[x.real,x.imag] for x in wS],Qc=Qc)
json.dump(out,open(DATA/'spectra.json','w'))
