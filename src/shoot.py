"""Frobenius(horizon, matrix, FFT-Taylor) + complex-ray shooting QNM solver for RN systems
   psi'' + (f'/f) psi' + (w^2 - V(r))/f^2 psi = 0,  V = (Delta/r^5) U(r)  (U: n x n matrix incl. couplings).  M=1 units."""
import numpy as np
from scipy.integrate import solve_ivp
class RN:
    def __init__(s,Q,Ufun,n):
        s.Q=Q; b=np.sqrt(1-Q*Q); s.rp=1+b; s.rm=1-b; s.d=s.rp-s.rm; s.U=Ufun; s.n=n
        s.kap=s.d/(2*s.rp**2)
    def rhs(s,r,y,w):
        n=s.n; p=y[:n]; dp=y[n:]
        D=(r-s.rp)*(r-s.rm); f=D/r**2; fp=(2*r-s.rp-s.rm)/r**2-2*D/r**3
        V=D/r**5*s.U(r)
        ddp=-(fp/f)*dp-((w*w*np.eye(n)-V)@p)/f**2
        return np.concatenate([dp,ddp])
    def frob(s,w,K=60,x0f=0.3):
        """ingoing horizon solutions (n of them), value & derivative at x0 = x0f*d"""
        n=s.n; d=s.d; rp=s.rp; m=256; rc=0.5           # circle |xi|=0.5 in xi=x/d
        th=2*np.pi*np.arange(m)/m; xi=rc*np.exp(1j*th); x=d*xi
        r=rp+x
        pv=1+x/(x+d)-2*x/r
        qv=np.array([r[i]**4*w*w/(x[i]+d)**2*np.eye(n)-x[i]*s.U(r[i])/(r[i]*(x[i]+d)) for i in range(m)])
        sc=rc**np.arange(K)
        pc=(np.fft.fft(pv)/m)[:K]/sc; qc=(np.fft.fft(qv,axis=0)/m)[:K]/sc[:,None,None]
        rho=-1j*w*rp**2/d
        sols=[]
        for j in range(n):
            a=[np.eye(n)[j].astype(complex)]
            for k in range(1,K):
                lhs=((k+rho)*(k+rho-1)+pc[0]*(k+rho))*np.eye(n)+qc[0]
                rhs=-sum((pc[jj]*(k-jj+rho))*a[k-jj]+qc[jj]@a[k-jj] for jj in range(1,k+1))
                a.append(np.linalg.solve(lhs,rhs))
            val=sum(a[k]*x0f**(k+rho) for k in range(K)); der=sum(a[k]*(k+rho)*x0f**(k+rho-1) for k in range(K))/d
            sols.append((val,der))
        return sols, rp+x0f*d
    def F(s,w,xm=2.0,ra=8.,R=None,vt=0.75*np.pi):
        """QNM function: det of [horizon sols | outgoing sols] at r_m = r+ + xm*d"""
        n=s.n; rm_=s.rp+xm*s.d
        sols,r0=s.frob(w)
        H=[]
        for (v,dv) in sols:
            sol=solve_ivp(lambda r,y: s.rhs(r,y,w),(r0,rm_),np.concatenate([v,dv]),rtol=1e-12,atol=1e-14,method='DOP853')
            H.append(sol.y[:,-1])
        # outgoing: start far on ray r = ra + t e^{i vt}, t=R ; integrate inward along ray to ra, then real axis to rm_
        if R is None: R=40./max(abs(w),0.02)
        O=[]
        for j in range(n):
            re=ra+R*np.exp(1j*vt)
            D=(re-s.rp)*(re-s.rm)
            lg=1j*w*re**2/D              # psi'/psi ~ i w dr*/dr  (leading)
            y0=np.concatenate([np.eye(n)[j],lg*np.eye(n)[j]]).astype(complex)
            e=np.exp(1j*vt)
            f1=lambda t,y: e*s.rhs(ra+t*e,y,w)
            sol=solve_ivp(f1,(R,0),y0,rtol=1e-12,atol=1e-300,method='DOP853')
            y=sol.y[:,-1]; y=y/np.max(abs(y))
            sol=solve_ivp(lambda r,yy: s.rhs(r,yy,w),(ra,rm_),y,rtol=1e-12,atol=1e-300,method='DOP853')
            O.append(sol.y[:,-1])
        Mx=np.array(H+O).T
        Mx=Mx/np.max(abs(Mx),axis=0)
        return np.linalg.det(Mx)
def secant(Ffun,z0,z1=None,tol=1e-11,it=60):
    z1=z0*(1+1e-3)+1e-5 if z1 is None else z1
    f0,f1=Ffun(z0),Ffun(z1)
    for _ in range(it):
        z2=z1-f1*(z1-z0)/(f1-f0); z0,f0=z1,f1; z1=z2; f1=Ffun(z1)
        if abs(z1-z0)<tol*max(1,abs(z1)): break
    return z1
def Ubranch(Q,l,branch):
    mu=np.sqrt((l-1)*(l+2)); s_=np.sqrt(9+4*Q*Q*mu*mu); L=l*(l+1)
    q={'G':3+s_,'EM':3-s_}[branch]
    return lambda r: np.array([[L*r+4*Q*Q/r-q]])
if __name__=='__main__':
    sv=RN(0.5,Ubranch(0.5,2,'G'),1)
    print('RN Q=0.5 grav-led l=2:',secant(sv.F,0.38-0.09j),' ref 0.3816771514-0.0896123793j')
    sv=RN(0.,Ubranch(0.,2,'G'),1)
    print('Schw l=2 n=1:',secant(sv.F,0.35-0.27j),' ref 0.3467110-0.2739149j')
