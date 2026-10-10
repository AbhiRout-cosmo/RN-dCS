"""Sec. XI C / Eq. (75): dCS shift of the fundamental gravitational-led ZDM (l=2) vs sigma=(r+ - r-)/r+.
Prints beta*dw/(kappa_+ b^2) -> C(sigma); quadratic fit gives C0 = 481.7  =>  dw/w = 4*pi*C0*zeta*(1-Q^2/M^2).
Writes data/zdm_shift.json.  Units M = alpha = 1.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from zdm_tools import *
import json
l=2; out=[]
for sig in [0.006,0.012,0.024,0.048]:
    b=sig/(2-sig); Q=np.sqrt(1-b*b)
    s1=RN(Q,Ubranch(Q,l,'G'),1); w0=zdm(s1,2,0)
    cs=[]
    for B in [1e6,2e6]:
        s2=RN(Q,Ucoupled(Q,l,B),2); w=secant(s2.F,w0,w0*(1+1e-6)); cs.append(B*(w-w0))
    c=2*cs[1]-cs[0]
    out.append([sig,b,s1.kap,w0.real,w0.imag,c.real,c.imag])
    print(f'sigma={sig}: b={b:.5f} kappa={s1.kap:.6f}  w0={w0.imag:+.9f}i  beta*dw={c:.4e}  [Richardson spread {abs(cs[1]-cs[0]):.1e}]  beta*dw/(kappa b^2)={c/(s1.kap*b*b):.4f}',flush=True)
json.dump(out,open(DATA/'zdm_shift.json','w'))
