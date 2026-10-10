"""Secs. VI-VII, App. B: reduction of the sourced axial Einstein-Maxwell system to master form.
Prints, for P = Psi_RW and the Maxwell function a, the coefficients of d^2/dr*^2 + ... : recovers V_RW, V_EM, C_Q and the
pseudoscalar sources (C_dCS in the Psi equation and chi*C_dCS in the Maxwell equation; alpha=1 without the 16*pi factor).
"""
import sympy as sp
r,M,Q,w,L,beta=sp.symbols('r M Q omega L beta')
p0,p1,p2,p3,a0,a1,a2,t0,t1,t2=sp.symbols('p0 p1 p2 p3 a0 a1 a2 t0 t1 t2')
D=r**2-2*M*r+Q**2
def d(ex):
    return sp.diff(ex,r)+sp.diff(ex,p0)*p1+sp.diff(ex,p1)*p2+sp.diff(ex,p2)*p3+sp.diff(ex,a0)*a1+sp.diff(ex,a1)*a2+sp.diff(ex,t0)*t1+sp.diff(ex,t1)*t2
h1=r**3*p0/D
h0=sp.I/w*(D/r**5)*(r*D*d(h1)+2*(M*r-Q**2)*h1)        # theta-phi (unsourced: C_thph=0)
CrW=3*sp.I*M*t0*w*r/(r**3*D)-3*sp.I*Q**2*t0*w/(r**3*D)  # C_{r phi}|_W  (alpha=1)
rp=-2*sp.I*a0*Q*w*r**6/D - sp.I*h0*w*r**7/D + sp.I*d(h0)*w*r**8/(2*D) - h1*r**4*(2*L*M*r-L*Q**2-L*r**2-4*M*r+2*Q**2+w**2*r**4+2*r**2)/(2*D) + r**6*CrW
mx=a0*r**4*(2*L*M*r-L*Q**2-L*r**2+w**2*r**4)/D + 2*a1*r**3*(M*r-Q**2) + a2*r**4*D + 2*h0*Q*r**3 - d(h0)*Q*r**4 - sp.I*h1*Q*w*r**4
ds=lambda F: D/r**2*d(F)
# P equation: solve for p2
c2=sp.cancel(sp.diff(rp,p2))
P_eq=rp*(D/r**2)**2/c2          # = d2P/dr*2 + (w^2-V)P + ... 
p2sol=sp.solve(rp,p2)[0]
# Maxwell with p2 eliminated (p3 doesn't appear)
mx2=sp.cancel(mx.subs(p2,p2sol))
c2a=sp.cancel(sp.diff(mx2,a2))
A_eq=mx2*(D/r**2)**2/c2a
def show(eq,F0,name):
    rem=eq-ds(ds(F0))
    print(name)
    for s in (p2,p1,p0,a2,a1,a0,t2,t1,t0):
        c=sp.cancel(sp.diff(rem,s))
        if c!=0: print('   ',s,':',c)
show(P_eq,p0,'Psi_RW eq  (d2/dr*2 P + sum coef*var = 0):')
show(A_eq,a0,'a eq       (d2/dr*2 a + sum coef*var = 0):')
