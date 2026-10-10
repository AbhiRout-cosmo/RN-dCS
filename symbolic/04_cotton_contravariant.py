"""Sec. V, App. A: contravariant C-tensor pieces C_(1) (Ricci) and C_(2) (dual Riemann) in the paper's conventions
(eps_{t r th ph}=+sqrt(-g), weight-1/2 symmetrization, S^phi = d_theta Y / sin theta):
C^{t phi} = [-3(Mr-Q^2) R'/r^6 + 3(Mr-2Q^2) R/r^7] S^phi,  C^{r phi} = -3 i w (Mr-Q^2) R S^phi / r^6.
"""
import sympy as sp, itertools
t,r,th,ph=sp.symbols('t r theta phi'); M,Q,w,L=sp.symbols('M Q omega L')
X=[t,r,th,ph]; f=1-2*M/r+Q**2/r**2
R=sp.Function('R')(r); Yf=sp.Function('Y')(th); E=sp.exp(-sp.I*w*t)
g=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2); gi=g.inv()
G=[[[sp.cancel(sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a in range(4)]
Rm={}
for a,b,c,d in itertools.product(range(4),repeat=4):
    Rm[a,b,c,d]=sp.cancel(sp.diff(G[a][d][b],X[c])-sp.diff(G[a][c][b],X[d])+sum(G[a][c][e]*G[e][d][b]-G[a][d][e]*G[e][c][b] for e in range(4)))
Ric=sp.Matrix(4,4,lambda b,d: sp.cancel(sum(Rm[a,b,a,d] for a in range(4)))); Rmix=gi*Ric
def dRmix(tau,nu,sig): return sp.cancel(sp.diff(Rmix[nu,sig],X[tau])+sum(G[nu][tau][l]*Rmix[l,sig]-G[l][tau][sig]*Rmix[nu,l] for l in range(4)))
sg=r**2*sp.sin(th)
def eps(*i): return -sp.LeviCivita(*i)/sg          # eps^{0123} = -1/sqrt(-g)
Rud={k:sp.cancel(sum(gi[k[1],e]*Rm[k[0],e,k[2],k[3]] for e in range(4))) for k in itertools.product(range(4),repeat=4)}
def dual(s,m,n,p): return sum(eps(n,p,a,b)*Rud[s,m,a,b] for a in range(4) for b in range(4) if a!=b)/2
th_=R*Yf*E; v=[sp.diff(th_,x) for x in X]
vv=[[sp.diff(v[b],X[a])-sum(G[l][a][b]*v[l] for l in range(4)) for b in range(4)] for a in range(4)]
def C1(m,n): return sum(sum(v[rho]*eps(rho,sig,tau,mu)*dRmix(tau,nu,sig) for rho in range(4) for sig in range(4) for tau in range(4) if len({rho,sig,tau,mu})==4) for (mu,nu) in [(m,n),(n,m)])/2
def C2(m,n): return sum(sum(vv[rho][sig]*dual(sig,mu,nu,rho) for rho in range(4) for sig in range(4)) for (mu,nu) in [(m,n),(n,m)])/2
Sphi=sp.Symbol('Sphi')
for name,fn in [('C1',C1),('C2',C2)]:
    for (m,n) in [(0,3),(1,3)]:
        ex=sp.expand(fn(m,n).subs(t,0))
        ex=ex.subs(sp.diff(Yf,th),Sphi*sp.sin(th))
        print(f'{name}^{X[m]}{X[n]} =',sp.factor(sp.simplify(ex)),flush=True)
