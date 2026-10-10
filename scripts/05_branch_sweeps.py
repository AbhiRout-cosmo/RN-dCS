"""Sec. IX, Fig. 1: (i) EM-led branch shift vs Q (sign reversal at Q_c), l=2..4; (ii) crossing charge Q_c(l), l=2..12;
(iii) large-l scaling of the gravitational-led shift at Q/M=0, 0.3, 0.6 (l^3 vs l^2/Q).  Writes data/c3_sweeps.json.  Runtime ~10 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from branch import *
import json
from scipy.optimize import brentq
out={}
# Schwarzschild seeds: grav (RW), EM, scalar via single-field ECS from WKB-ish guesses refined by continuation in l
def seed(l):
    x=(l+0.5)/np.sqrt(27)
    return x-0.0962j*0.92, x+0.04-0.096j, x+0.06-0.0965j
res={'em':{},'Qc':{},'lscal':{}}
# 1) EM-led linear coefficient, l=2..4, fine Q grid, exact branch identity
for l in [2,3,4]:
    g2=gr(1.,1e-9,l,2,{2:0.4576-0.0950j,3:0.6569-0.0956j,4:0.8531-0.0959j}[l]); rows=[]
    for Q in np.round(np.r_[np.arange(0.02,0.96,0.02),0.98],3):
        N=300 if Q<0.9 else 400; Sb=60 if Q<0.9 else 90
        w0,c=coef_b(1.,Q,l,2,g2,N=N,S=Sb); g2=w0; rows.append([Q,w0.real,w0.imag,c.real,c.imag])
    res['em'][l]=rows; print('EM-led done l=',l,flush=True)
# 2) Q_c(l): Re w_2 = Re w_Theta
for l in range(2,13):
    gb=gr(1.,1e-9,l,2,(l+0.5)/np.sqrt(27)+0.03-0.096j); gs=sc(1.,0.,l,(l+0.5)/np.sqrt(27)+0.05-0.0965j)
    Qs=np.arange(0.01,0.9,0.01); prev=None; Qc=None
    for Q in Qs:
        gb=gr(1.,Q,l,2,gb); gs=sc(1.,Q,l,gs); d=gb.real-gs.real
        if prev is not None and np.sign(d)!=np.sign(prev[1]):
            Qc=prev[0]-prev[1]*(Q-prev[0])/(d-prev[1]); break
        prev=(Q,d)
    res['Qc'][l]=Qc; print(f'l={l}: Q_c={Qc:.4f}',flush=True)
# 3) l-scaling of grav-led coefficient at fixed Q
for Q in [0.0,0.3,0.6]:
    rows=[]; g1=None
    for l in range(2,13):
        g=(l+0.5)/np.sqrt(27)-0.0962j*0.95 if g1 is None else g1*(l+0.5)/(l-0.5)
        w0,c=coef_b(1.,Q,l,1,g,B=2e5 if l>6 else 2e4,N=320); g1=w0; rows.append([l,w0.real,w0.imag,c.real,c.imag])
        print(f'Q={Q} l={l}: w1={w0:.6f}  beta*dw={c:.4f}',flush=True)
    res['lscal'][Q]=rows
json.dump(res,open(DATA/'c3_sweeps.json','w'))
