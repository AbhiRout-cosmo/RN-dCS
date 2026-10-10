"""App. E: validation of the Frobenius--shooting ZDM solver against Zimmerman & Mark, PRD 93, 044033 (2016), Fig. 10.
Prints the scaled residual M (n+1/2)^-1 sigma^-2 Delta(omega_I) for (j=1,l=1), (j=2,l=2), (j=1,l=2).
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from zdm_tools import RN, Ubranch, secant, np
from scipy.optimize import minimize_scalar
def zdm(sv,h,n,width=0.06):
    p=sv.kap*(h+n); g=np.linspace(p*(1-width),p*(1+width),61)
    v=np.array([abs(sv.F(-1j*x)) for x in g]); i=np.argmin(v)
    lo,hi=g[max(i-1,0)],g[min(i+1,len(g)-1)]
    res=minimize_scalar(lambda x: abs(sv.F(-1j*x)),bounds=(lo,hi),method='bounded',options=dict(xatol=1e-13))
    z=secant(sv.F,-1j*res.x,-1j*res.x*(1+1e-7))
    return z, v[i]/np.median(v)
for (l,br,h,lab,zm) in [(1,'EM',3,'ZM j=1,l=1',(3.0,1.35)),(2,'G',2,'ZM j=2,l=2',(2.0,1.0)),(2,'EM',4,'ZM j=1,l=2',(4.0,1.6))]:
    print(lab,' ZM Fig.10 |resid| ~',zm)
    for sig in [0.006,0.012]:
        b=sig/(2-sig); Q=np.sqrt(1-b*b); sv=RN(Q,Ubranch(Q,l,br),1)
        out=[]
        for n in [0,1]:
            z,dip=zdm(sv,h,n); pred=-sv.kap*(h+n)
            out.append(f'n={n}: Mw={z.real:+.1e}{z.imag:+.9f}i resid={(z.imag-pred)/((n+0.5)*sig**2):+.3f} (dip {dip:.0e})')
        print(f'  sigma={sig}: '+' | '.join(out),flush=True)
