"""Helpers for zero-damping-mode (ZDM) root finding and the branch-reduced {Z_g, Theta} coupling matrix."""
from shoot import RN, Ubranch, secant, np
from scipy.optimize import minimize_scalar
import json
def zdm(sv,h,n,width=0.06):
    p=sv.kap*(h+n); g=np.linspace(p*(1-width),p*(1+width),61)
    v=np.array([abs(sv.F(-1j*x)) for x in g]); i=np.argmin(v)
    res=minimize_scalar(lambda x: abs(sv.F(-1j*x)),bounds=(g[max(i-1,0)],g[min(i+1,60)]),method='bounded',options=dict(xatol=1e-14))
    return secant(sv.F,-1j*res.x,-1j*res.x*(1+1e-7))
def Ucoupled(Q,l,beta,branch='G'):
    mu=np.sqrt((l-1)*(l+2)); s_=np.sqrt(9+4*Q*Q*mu*mu); L=l*(l+1)
    q={'G':3+s_,'EM':3-s_}[branch]; phi=np.arctan((s_-3)/(2*Q*mu))
    def U(r):
        chi=2*Q/(r*mu); w={'G':np.cos(phi)-chi*np.sin(phi),'EM':np.sin(phi)+chi*np.cos(phi)}[branch]
        k=(r-Q*Q)/r**3
        return np.array([[L*r+4*Q*Q/r-q, 96*np.pi*k*w],[6*L*(L-2)*k*w/beta, L*r+2-2*Q*Q/r]])
    return U
