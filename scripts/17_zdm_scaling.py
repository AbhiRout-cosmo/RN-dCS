"""Sec. XI C, Fig. 4: ZDM shift vs b=sqrt(1-Q^2) (l=2 fundamental grav-led), sigma=0.004..0.1; free-exponent fit;
mode tracking: only purely damped roots continuous in sigma are kept (sigma>=0.07 rejected).  Writes data/zdm_scaling.json.  ~8 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
from zdm_tools import RN, Ubranch, secant, zdm, Ucoupled
out={'zdm':[]}
t0=time.time()
l=2
for sig in (0.004,0.006,0.009,0.012,0.018,0.024,0.036,0.048,0.07,0.1):
    b=sig/(2-sig); Q=np.sqrt(1-b*b)
    s1=RN(Q,Ubranch(Q,l,'G'),1); w0=zdm(s1,2,0)
    cs=[]
    for B in (1e6,2e6):
        s2=RN(Q,Ucoupled(Q,l,B),2); w=secant(s2.F,w0,w0*(1+1e-6)); cs.append(B*(w-w0))
    c=2*cs[1]-cs[0]
    out['zdm'].append([sig,b,s1.kap,w0,c,abs(cs[1]-cs[0])])
    print(f'sigma={sig} b={b:.5f} w0/kap={w0/s1.kap:.6f} beta*dw={c:.5e} C={(-c/(1j*s1.kap*b*b)).real:.3f} spread {abs(cs[1]-cs[0]):.1e} %.0fs'%(time.time()-t0),flush=True)
    save('zdm_scaling',out)


# free-exponent fit (only purely damped roots, Re w0 = 0)
rows=[r for r in out['zdm'] if abs(complex(r[3]).real)<1e-6]
b=np.array([r[1] for r in rows]); w0=np.array([r[3] for r in rows]); c=np.array([r[4] for r in rows]); rel=(c/w0).real
print('free-exponent fit dw/w0 ~ b^p: p =',np.polyfit(np.log(b),np.log(rel),1)[0])
print('local exponents:',np.round(np.diff(np.log(rel))/np.diff(np.log(b)),4))
