"""Linearized axial Einstein-Maxwell equations on RN (t-phi, r-phi, theta-phi and Maxwell phi components), Regge-Wheeler gauge,
background A_t = +Q/r (Gaussian-geometrized: G_mn = 2(F F - g F^2/4)). Coefficients of the axial harmonic W are printed.
"""
import sympy as sp
t,r,th,ph=sp.symbols('t r theta phi')
M,Q,w,L=sp.symbols('M Q omega L')
X=[t,r,th,ph]
f=1-2*M/r+Q**2/r**2
h0=sp.Function('h0')(r); h1=sp.Function('h1')(r); a=sp.Function('a')(r)
E=sp.exp(-sp.I*w*t)
Wf=sp.Function('W')(th)
g0=sp.diag(-f,1/f,r**2,r**2*sp.sin(th)**2); g0i=g0.inv()
h=sp.zeros(4); h[0,3]=h[3,0]=h0*Wf*E; h[1,3]=h[3,1]=h1*Wf*E
G0=[[[sp.simplify(sum(g0i[a_,d]*(sp.diff(g0[d,b],X[c])+sp.diff(g0[d,c],X[b])-sp.diff(g0[b,c],X[d])) for d in range(4))/2) for c in range(4)] for b in range(4)] for a_ in range(4)]
hu=g0i*h
dG=[[[sp.expand(sum(g0i[a_,d]*(sp.diff(h[d,b],X[c])+sp.diff(h[d,c],X[b])-sp.diff(h[b,c],X[d])) for d in range(4))/2 - sum(hu[a_,d]*G0[d][b][c] for d in range(4))) for c in range(4)] for b in range(4)] for a_ in range(4)]
def dRic(b,d):
    return sp.expand(sum(sp.diff(dG[a_][d][b],X[a_])-sp.diff(dG[a_][a_][b],X[d]) +sum(dG[a_][a_][e]*G0[e][d][b]+G0[a_][a_][e]*dG[e][d][b]-dG[a_][d][e]*G0[e][a_][b]-G0[a_][d][e]*dG[e][a_][b] for e in range(4)) for a_ in range(4)))
F0=sp.zeros(4); F0[0,1]=Q/r**2; F0[1,0]=-Q/r**2
dA=[0,0,0,a*Wf*E]
dF=sp.Matrix(4,4,lambda m,n: sp.diff(dA[n],X[m])-sp.diff(dA[m],X[n]))
dgi=-g0i*h*g0i
def dT(m,n):
    s=sum(dF[m,al]*F0[n,be]*g0i[al,be]+F0[m,al]*dF[n,be]*g0i[al,be]+F0[m,al]*F0[n,be]*dgi[al,be] for al in range(4) for be in range(4))
    return sp.expand(2*(s - h[m,n]*(-2*Q**2/r**4)/4))
sg=r**2*sp.sin(th)
def dFu(m,n):
    return sum(g0i[m,al]*g0i[n,be]*dF[al,be]+dgi[m,al]*g0i[n,be]*F0[al,be]+g0i[m,al]*dgi[n,be]*F0[al,be] for al in range(4) for be in range(4))

import time
Ws,Wp,W2s,cs,sn=sp.symbols('Ws Wp W2s cs sn')
H0,H0p,H0pp,H1,H1p,A0,A1,A2=sp.symbols('H0 H0p H0pp H1 H1p A0 A1 A2')
def sym(ex):
    ex=ex.subs(t,0)
    ex=ex.subs(Wf.diff(th,2),W2s).subs(Wf.diff(th),Wp).subs(Wf,Ws)
    ex=ex.subs({h0.diff(r,2):H0pp,h1.diff(r,2):sp.Symbol('H1pp'),a.diff(r,2):A2}).subs({h0.diff(r):H0p,h1.diff(r):H1p,a.diff(r):A1}).subs({h0:H0,h1:H1,a:A0})
    ex=ex.subs(W2s,sp.cos(th)/sp.sin(th)*Wp-L*Ws)
    ex=ex.subs(sp.sin(2*th),2*sp.sin(th)*sp.cos(th)).subs({sp.sin(th):sn,sp.cos(th):cs})
    return sp.expand(ex)
T=time.time()
for k,ex in [('tp',dRic(0,3)-dT(0,3)),('rp',dRic(1,3)-dT(1,3)),('thp',dRic(2,3)-dT(2,3)),('mx',sum(sp.diff(sg*dFu(m,3),X[m]) for m in range(4)))]:
    e=sym(ex)
    cW=sp.collect(sp.expand(sp.cancel(sp.together(e.coeff(Ws)))*r**6),[H0pp,H0p,H0,H1p,H1,A2,A1,A0],sp.factor); cWp=sp.cancel(sp.together(e.coeff(Wp)))
    print(k,time.time()-T,'\n  W*r^6 :',cW,'\n  Wp:',cWp,flush=True)
