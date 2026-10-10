"""C2: coupled RN axial {Psi, Psi_e, Theta} QNMs in dCS by exterior complex scaling.
Tortoise contour r* = e^{i th} s, s in [-S,S], Chebyshev collocation, Dirichlet ends.
H v = w^2 v,  H = -d2/dr*2 + V + Coff  (linear eigenproblem in w^2)."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import scipy.linalg as sla

def cheb(N):
    x=np.cos(np.pi*np.arange(N+1)/N); c=np.r_[2,np.ones(N-1),2]*(-1)**np.arange(N+1)
    X=np.tile(x,(N+1,1)).T; dX=X-X.T
    D=np.outer(c,1/c)/(dX+np.eye(N+1)); D-=np.diag(D.sum(1)); return D,x

def rgrid(M,Q,th,s):
    rp=M+np.sqrt(M*M-Q*Q); rm=M-np.sqrt(M*M-Q*Q)
    def rstar(r):
        if rm>0: return r+rp**2/(rp-rm)*np.log(r-rp)-rm**2/(rp-rm)*np.log(r-rm)
        return r+2*M*np.log(r-rp)
    r0=brentq(rstar,rp+1e-12,rp+50)
    # u = log(r - r+),  du/ds = e^{i th} (r - r-)/r^2
    def rhs(t,y):
        u=y[0]+1j*y[1]; r=rp+np.exp(u); du=np.exp(1j*th)*(r-rm)/r**2
        return [du.real,du.imag]
    u0=np.log(r0-rp); out=np.empty(len(s),complex)
    for sel in (s>=0, s<0):
        ss=s[sel]; order=np.argsort(np.abs(ss)); tgt=ss[order]
        if len(tgt)==0: continue
        sol=solve_ivp(rhs,(0,tgt[-1]),[u0.real,u0.imag],t_eval=tgt,rtol=1e-13,atol=1e-14,method='DOP853')
        uu=sol.y[0]+1j*sol.y[1]; tmp=np.empty(len(ss),complex); tmp[order]=rp+np.exp(uu); out[sel]=tmp
    return out

def operator(M,Q,l,beta,th=0.6,S=60.,N=240,chi_on=True,coupled=True):
    D,x=cheb(N); s=S*x; D=D/S; D2=D@D
    r=rgrid(M,Q,th,s); Dl=r*r-2*M*r+Q*Q; L=l*(l+1)
    VRW=Dl/r**5*(L*r-6*M+4*Q*Q/r); VEM=Dl/r**5*(L*r+4*Q*Q/r); VT=Dl/r**5*(L*r+2*M-2*Q*Q/r)
    CQ=2*Q*np.sqrt((l-1)*(l+2))*Dl/r**5
    K=(M*r-Q*Q)*Dl/r**8
    CdCS=96*np.pi*K; CTg=6*(l-1)*l*(l+1)*(l+2)*K/beta if coupled else 0*K
    chi=2*Q/(r*np.sqrt((l-1)*(l+2))) if chi_on else 0*r
    n=N+1; I=np.arange(1,N)          # interior (Dirichlet)
    lap=-np.exp(-2j*th)*D2[np.ix_(I,I)]
    dg=lambda v: np.diag(v[I])
    Z=np.zeros((N-1,N-1))
    cd=CdCS if coupled else 0*K
    H=np.block([[lap+dg(VRW), dg(CQ),        dg(cd)],
                [dg(CQ),      lap+dg(VEM),   dg(chi*cd)],
                [dg(CTg),     dg(chi*CTg),   lap+dg(VT)]])
    return H

def qnm(M,Q,l,beta,w_guess,**kw):
    H=operator(M,Q,l,beta,**kw)
    ev=sla.eigvals(H); w=np.sqrt(ev); w=np.where(w.real<0,-w,w)
    return w[np.argmin(abs(w-w_guess))]
