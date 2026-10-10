"""Sec. VII A / App. C, Eq. (dstarRR): linearized Pontryagin density for Regge-Wheeler-gauge axial perturbations of RN.
Result: d(*RR) = 24 L (Mr-Q^2)/r^7 [r h0' - 2 h0 + i w r h1] Y_l  in the orientation eps^{0123}=+1/sqrt(-g) used here;
the paper's convention eps_{t r th ph}=+sqrt(-g) flips the overall sign (Eq. dstarRR carries a minus sign). Runtime ~10 s.
"""
import sympy as sp
t,r,th,ph=sp.symbols('t r theta phi')
M,Q,w,L=sp.symbols('M Q omega L')
X=[t,r,th,ph]
f=1-2*M/r+Q**2/r**2
Y=sp.Function('Y')(th)
h0=sp.Function('h0')(r); h1=sp.Function('h1')(r)
E=sp.exp(-sp.I*w*t)
W=sp.sin(th)*sp.diff(Y,th)
g0=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2)
g0i=g0.inv()
h=sp.zeros(4)
h[0,3]=h[3,0]=h0*W*E
h[1,3]=h[3,1]=h1*W*E
def Gam(g,gi):
    return [[[sp.simplify(sum(gi[a,d]*(sp.diff(g[d,b],X[c])+sp.diff(g[d,c],X[b])-sp.diff(g[b,c],X[d])) for d in range(4))/2)
              for c in range(4)] for b in range(4)] for a in range(4)]
G0=Gam(g0,g0i)
# delta Gamma
hu=g0i*h  # h^a_d
dG=[[[sp.expand(sum(g0i[a,d]*(sp.diff(h[d,b],X[c])+sp.diff(h[d,c],X[b])-sp.diff(h[b,c],X[d])) for d in range(4))/2
       - sum(hu[a,d]*G0[d][b][c] for d in range(4))) for c in range(4)] for b in range(4)] for a in range(4)]
def Riem(G):  # R^a_{bcd}
    R={}
    for a in range(4):
     for b in range(4):
      for c in range(4):
       for d in range(c+1,4):
        R[a,b,c,d]=sp.diff(G[a][d][b],X[c])-sp.diff(G[a][c][b],X[d])+sum(G[a][c][e]*G[e][d][b]-G[a][d][e]*G[e][c][b] for e in range(4))
    return R
R0=Riem(G0)
R0={k:sp.simplify(v) for k,v in R0.items()}
dR={}
for a in range(4):
 for b in range(4):
  for c in range(4):
   for d in range(c+1,4):
    dR[a,b,c,d]=sp.expand(sp.diff(dG[a][d][b],X[c])-sp.diff(dG[a][c][b],X[d])
       +sum(dG[a][c][e]*G0[e][d][b]+G0[a][c][e]*dG[e][d][b]-dG[a][d][e]*G0[e][c][b]-G0[a][d][e]*dG[e][c][b] for e in range(4)))
def get(R,a,b,c,d):
    if c==d: return 0
    return R[a,b,c,d] if c<d else -R[a,b,d,c]
sqrtg=r**2*sp.sin(th)
eps=lambda *i: sp.LeviCivita(*i)/sqrtg   # eps^{0123}=+1/sqrt(-g)  (orientation O+)
# delta(*RR) = eps^{rho sig tau eta} R0^nu_{mu rho sig} dR^mu_{nu tau eta}
import itertools
dP=0
for p in itertools.permutations(range(4)):
    rho,sig,tau,eta=p
    s=eps(*p)
    for mu in range(4):
        for nu in range(4):
            a=get(R0,nu,mu,rho,sig)
            if a==0: continue
            b=get(dR,mu,nu,tau,eta)
            if b==0: continue
            dP+=s*a*b
dP=sp.expand(dP/E)
# reduce Y derivatives via Legendre: Y'' = -cot Y' - L Y
def red(ex):
    for k in (4,3,2):
        for _ in range(3):
            ex=ex.subs(sp.diff(Y,th,k), sp.diff(-sp.cos(th)/sp.sin(th)*sp.diff(Y,th)-L*Y,th,k-2))
    return sp.expand(ex)
for _ in range(3): dP=red(dP)
Yp=sp.Symbol('Yp'); Ys=sp.Symbol('Ys')
dP=sp.expand(dP.subs(sp.diff(Y,th),Yp).subs(Y,Ys))
cY=sp.simplify(dP.coeff(Ys)); cYp=sp.simplify(dP.coeff(Yp))
print('coef Y  :',sp.factor(cY))
print('coef Yp :',sp.simplify(cYp))
