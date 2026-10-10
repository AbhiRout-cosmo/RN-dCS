"""Covariant C-tensor components C_{t phi}, C_{r phi}, C_{theta phi} for d(theta) = (Theta/r) Y e^{-iwt} on RN
(weight-1/2 symmetrization, orientation eps^{0123}=+1/sqrt(-g); flip the sign for the paper's convention). C_{theta phi}=0.
"""
import sympy as sp, itertools, time
t,r,th,ph=sp.symbols('t r theta phi'); M,Q,w,L=sp.symbols('M Q omega L')
X=[t,r,th,ph]; f=1-2*M/r+Q**2/r**2
Th=sp.Function('Th')(r); Yf=sp.Function('Y')(th); E=sp.exp(-sp.I*w*t)
g=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2); gi=g.inv()
G=[[[sp.cancel(sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
R={}
for a,b,c,d in itertools.product(range(4),repeat=4):
    R[a,b,c,d]=sp.cancel(sp.diff(G[a][d][b],X[c])-sp.diff(G[a][c][b],X[d])+sum(G[a][c][e]*G[e][d][b]-G[a][d][e]*G[e][c][b] for e in range(4)))
Ric=sp.Matrix(4,4,lambda b,d: sp.cancel(sum(R[a,b,a,d] for a in range(4))))
Rmix=gi*Ric   # R^nu_sigma
# nabla_tau R^nu_sigma
def dRmix(tau,nu,sig):
    return sp.cancel(sp.diff(Rmix[nu,sig],X[tau])+sum(G[nu][tau][l]*Rmix[l,sig]-G[l][tau][sig]*Rmix[nu,l] for l in range(4)))
sg=r**2*sp.sin(th)
def eps(*i): return sp.LeviCivita(*i)/sg
# R^{ab}_{cd}
Rud={}
for a,b,c,d in itertools.product(range(4),repeat=4):
    Rud[a,b,c,d]=sp.cancel(sum(gi[b,e]*R[a,e,c,d] for e in range(4)))
# dual  *R^{s m n p} = 1/2 eps^{n p t e} R^{s m}_{t e}
def dual(s,m,n,p):
    return sum(eps(n,p,a,b)*Rud[s,m,a,b] for a in range(4) for b in range(4) if a!=b)/2
th_=Th/r*Yf*E
v=[sp.diff(th_,x) for x in X]
vv=[[sp.diff(v[b],X[a])-sum(G[l][a][b]*v[l] for l in range(4)) for b in range(4)] for a in range(4)]
T0=time.time()
def Cup(m,n):
    s=0
    for (mu,nu) in [(m,n),(n,m)]:
        s+=sum(v[rho]*eps(rho,sig,tau,mu)*dRmix(tau,nu,sig) for rho in range(4) for sig in range(4) for tau in range(4) if len({rho,sig,tau,mu})==4)
        s+=sum(vv[rho][sig]*dual(sig,mu,nu,rho) for rho in range(4) for sig in range(4))
    return s/2
Cl={}
for (m,n) in [(0,3),(1,3),(2,3)]:
    cu={}
    Cl[(m,n)]=sum(g[m,a]*g[n,b]*Cup(a,b) for a in range(4) for b in range(4) if g[m,a]!=0 and g[n,b]!=0)
    print((m,n),'built',time.time()-T0,flush=True)
Ys,Yp,Y2,sn,cs=sp.symbols('Ys Yp Y2 sn cs')
T2,T1,T0s=sp.symbols('T2 T1 T0')
out={}
for k,ex in Cl.items():
    ex=ex.subs(t,0).subs(Yf.diff(th,2),Y2).subs(Yf.diff(th),Yp).subs(Yf,Ys)
    ex=ex.subs({Th.diff(r,2):T2}).subs({Th.diff(r):T1}).subs({Th:T0s})
    ex=ex.subs(Y2,-sp.cos(th)/sp.sin(th)*Yp-L*Ys)
    ex=sp.expand(ex.subs({sp.sin(th):sn,sp.cos(th):cs}))
    out[k]={s:sp.cancel(ex.coeff(s)) for s in (Ys,Yp)}
    print(k,out[k],flush=True)
