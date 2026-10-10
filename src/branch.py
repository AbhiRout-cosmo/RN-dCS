"""C3: branch-resolved system. Exact GR decoupling with constant angle phi:
   Z1 = cos(phi) Psi - sin(phi) Psi_e  (grav-led, V1 = D/r^5[Lr+4Q^2/r-q1]),  Z2 = sin(phi) Psi + cos(phi) Psi_e (EM-led, V2 with q2)
   q_{1,2} = 3M +- s,  s = sqrt(9M^2+4Q^2 mu^2),  tan(phi) = (s-3M)/(2 Q mu)
   pseudoscalar couples to Phi = Psi + chi Psi_e = w1 Z1 + w2 Z2,  w1 = cos(phi) - chi sin(phi), w2 = sin(phi) + chi cos(phi)"""
import numpy as np, scipy.linalg as sla
from ecs import cheb, rgrid
def branch_op(M,Q,l,beta,branch,th=0.6,S=60.,N=260,coupled=True,chi_on=True,field=None):
    D,x=cheb(N); s_=S*x; D=D/S; D2=D@D
    r=rgrid(M,Q,th,s_); Dl=r*r-2*M*r+Q*Q; L=l*(l+1); mu=np.sqrt((l-1)*(l+2))
    s=np.sqrt(9*M*M+4*Q*Q*mu*mu); q={1:3*M+s,2:3*M-s}[branch]
    phi=np.arctan((s-3*M)/(2*Q*mu)) if Q>0 else 0.0
    chi=2*Q/(r*mu) if chi_on else 0*r
    w={1:np.cos(phi)-chi*np.sin(phi),2:np.sin(phi)+chi*np.cos(phi)}[branch]
    Vb=Dl/r**5*(L*r+4*Q*Q/r-q); VT=Dl/r**5*(L*r+2*M-2*Q*Q/r)
    K=(M*r-Q*Q)*Dl/r**8; Cg=96*np.pi*K*w; Cs=6*(l-1)*l*(l+1)*(l+2)*K*w/beta
    if not coupled: Cg=Cs=0*K
    I=np.arange(1,N); lap=-np.exp(-2j*th)*D2[np.ix_(I,I)]; dg=lambda v: np.diag(v[I])
    if field=='b': return lap+dg(Vb)
    if field=='s': return lap+dg(VT)
    return np.block([[lap+dg(Vb),dg(Cg)],[dg(Cs),lap+dg(VT)]])
def eig_near(H,g):
    ev=np.sqrt(sla.eigvals(H)); ev=np.where(ev.real<0,-ev,ev); return ev[np.argmin(abs(ev-g))]
def gr(M,Q,l,branch,g,**kw): return eig_near(branch_op(M,Q,l,1,branch,field='b',**kw),g)
def sc(M,Q,l,g,**kw): return eig_near(branch_op(M,Q,l,1,1,field='s',**kw),g)
def coef_b(M,Q,l,branch,g,chi_on=True,B=2e4,**kw):
    w0=gr(M,Q,l,branch,g,**kw)
    c1=B*(eig_near(branch_op(M,Q,l,B,branch,chi_on=chi_on,**kw),w0)-w0)
    c2=2*B*(eig_near(branch_op(M,Q,l,2*B,branch,chi_on=chi_on,**kw),w0)-w0)
    return w0,2*c2-c1
