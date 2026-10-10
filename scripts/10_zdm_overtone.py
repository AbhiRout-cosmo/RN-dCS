"""Sec. XI C and App. E: (1) robustness of the ZDM shift to the matching radius / Frobenius point;
(2) near-degenerate gravitational-led n=1 ZDM vs scalar n=0 ZDM: gap ~ sigma^2 and first-order shift ~ kappa_+.
Writes data/zdm_n1.json.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from zdm_tools import *
import json
l=2
# (1) solver-parameter robustness at sigma=0.012
sig=0.012; b=sig/(2-sig); Q=np.sqrt(1-b*b)
s1=RN(Q,Ubranch(Q,l,'G'),1); w0=zdm(s1,2,0)
for xm,x0f in [(2.0,0.3),(3.0,0.2),(1.5,0.35)]:
    s2=RN(Q,Ucoupled(Q,l,1e6),2)
    s2.F_=s2.F; F=lambda w: s2.F_(w,xm=xm)
    s2.frob_=s2.frob; s2.frob=lambda w,K=60,x0f_=x0f: s2.frob_(w,K=K,x0f=x0f_)
    w=secant(F,w0,w0*(1+1e-6)); print(f'  robustness xm={xm} x0f={x0f}: beta*dw/(kappa b^2) = {1e6*(w-w0)/(s1.kap*b*b):.3f}',flush=True)
# (2) first overtone n=1 (near-degenerate with scalar ZDM h=3)
print('grav-led n=1 and scalar n=0 (h=3) ZDMs, and the n=1 shift:')
res=[]
for sig in [0.012,0.024,0.048]:
    b=sig/(2-sig); Q=np.sqrt(1-b*b)
    s1=RN(Q,Ubranch(Q,l,'G'),1); w1=zdm(s1,2,1)
    sT=RN(Q,lambda r: np.array([[l*(l+1)*r+2-2*Q*Q/r]]),1); wT=zdm(sT,3,0)
    cs=[]
    for B in [1e6,2e6]:
        s2=RN(Q,Ucoupled(Q,l,B),2); w=secant(s2.F,w1,w1*(1+1e-6)); cs.append(B*(w-w1))
    c=2*cs[1]-cs[0]; res.append([sig,w1.imag,wT.imag,c.real,c.imag])
    print(f'  sigma={sig}: w_G(n=1)={w1.imag:+.9f}i  w_Theta(n=0)={wT.imag:+.9f}i  gap/kappa={(w1.imag-wT.imag)/s1.kap:+.2e}   beta*dw_G(n=1)/kappa={c/s1.kap:.3e}  [spread {abs(cs[1]-cs[0])/s1.kap:.1e}]',flush=True)
json.dump(res,open(DATA/'zdm_n1.json','w'))
