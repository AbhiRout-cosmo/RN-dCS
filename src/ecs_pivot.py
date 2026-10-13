"""ECS helpers with movable rotation pivot r*_c (rc=0 reproduces ecs.rgrid), full 3-field operator, Richardson coefficient, JSON saver."""
import sys,json,time,pathlib,numpy as np,scipy.linalg as sla
DATA=pathlib.Path(__file__).resolve().parents[1]/'data'
from ecs import cheb
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
def rgrid_c(M,Q,th,s,rc=0.0):
    """ECS map with adjustable pivot: r*(s)=rc+e^{i th}(s-rc)  (rc=0 reproduces src/ecs.rgrid)"""
    rp=M+np.sqrt(M*M-Q*Q); rm=M-np.sqrt(M*M-Q*Q)
    def rstar(r):
        if rm>0: return r+rp**2/(rp-rm)*np.log(r-rp)-rm**2/(rp-rm)*np.log(r-rm)
        return r+2*M*np.log(r-rp)
    r0=brentq(lambda r: rstar(r)-rc,rp+1e-12,rp+80)
    def rhs(t,y):
        u=y[0]+1j*y[1]; r=rp+np.exp(u); du=np.exp(1j*th)*(r-rm)/r**2
        return [du.real,du.imag]
    u0=np.log(r0-rp); out=np.empty(len(s),complex); sc=s-rc
    for sel in (sc>=0, sc<0):
        ss=sc[sel]; order=np.argsort(np.abs(ss)); tgt=ss[order]
        if len(tgt)==0: continue
        sol=solve_ivp(rhs,(0,tgt[-1]),[u0.real,u0.imag],t_eval=tgt,rtol=1e-13,atol=1e-14,method='DOP853')
        uu=sol.y[0]+1j*sol.y[1]; tmp=np.empty(len(ss),complex); tmp[order]=rp+np.exp(uu); out[sel]=tmp
    return out
def op_full(M,Q,l,beta,th=0.6,S=60.,N=300,rc=0.0,chi_on=True,coupled=True):
    D,x=cheb(N); s=S*x; D=D/S; D2=D@D
    r=rgrid_c(M,Q,th,s,rc); Dl=r*r-2*M*r+Q*Q; L=l*(l+1)
    VRW=Dl/r**5*(L*r-6*M+4*Q*Q/r); VEM=Dl/r**5*(L*r+4*Q*Q/r); VT=Dl/r**5*(L*r+2*M-2*Q*Q/r)
    CQ=2*Q*np.sqrt((l-1)*(l+2))*Dl/r**5; K=(M*r-Q*Q)*Dl/r**8
    cd=96*np.pi*K if coupled else 0*K; CTg=6*(l-1)*l*(l+1)*(l+2)*K/beta if coupled else 0*K
    chi=2*Q/(r*np.sqrt((l-1)*(l+2))) if chi_on else 0*r
    I=np.arange(1,N); lap=-np.exp(-2j*th)*D2[np.ix_(I,I)]; dg=lambda v: np.diag(v[I])
    return np.block([[lap+dg(VRW),dg(CQ),dg(cd)],[dg(CQ),lap+dg(VEM),dg(chi*cd)],[dg(CTg),dg(chi*CTg),lap+dg(VT)]])
def evals(H):
    w=np.sqrt(sla.eigvals(H)); return np.where(w.real<0,-w,w)
def near(ev,g): return ev[np.argmin(abs(ev-g))]
def coef_full(M,Q,l,g,B=2e4,**kw):
    w0=near(evals(op_full(M,Q,l,1,coupled=False,**kw)),g)
    c1=B*(near(evals(op_full(M,Q,l,B,**kw)),w0)-w0); c2=2*B*(near(evals(op_full(M,Q,l,2*B,**kw)),w0)-w0)
    return w0,2*c2-c1,abs(c2-c1)
def save(name,obj):
    def conv(o):
        if isinstance(o,complex) or isinstance(o,np.complexfloating): return [float(o.real),float(o.imag)]
        if isinstance(o,np.floating): return float(o)
        if isinstance(o,np.ndarray): return [conv(x) for x in o]
        if isinstance(o,dict): return {str(k):conv(v) for k,v in o.items()}
        if isinstance(o,(list,tuple)): return [conv(x) for x in o]
        return o
    json.dump(conv(obj),open(DATA/f'{name}.json','w'),indent=1)
