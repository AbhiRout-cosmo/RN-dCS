"""App. B/C consistency check: the unused t-phi Einstein equation, including its C-tensor source, vanishes identically on shell
(exact rational arithmetic at random (M, Q, omega, l)).
"""
import sympy as sp, random
r=sp.symbols('r')
p0,p1,p2,p3,a0,a1,a2,t0,t1,t2,T2s=sp.symbols('p0 p1 p2 p3 a0 a1 a2 t0 t1 t2 T2s')
def run(M,Q,w,L):
    D=r**2-2*M*r+Q**2
    def d(ex):
        return sp.diff(ex,r)+sp.diff(ex,p0)*p1+sp.diff(ex,p1)*p2+sp.diff(ex,p2)*p3+sp.diff(ex,a0)*a1+sp.diff(ex,a1)*a2+sp.diff(ex,t0)*t1+sp.diff(ex,t1)*t2
    ds=lambda F: D/r**2*d(F)
    h1=r**3*p0/D
    h0=sp.I/w*(D/r**5)*(r*D*d(h1)+2*(M*r-Q**2)*h1)
    CrW=3*sp.I*w*(M*r-Q**2)*t0/(r**3*D)
    rp=-2*sp.I*a0*Q*w*r**6/D - sp.I*h0*w*r**7/D + sp.I*d(h0)*w*r**8/(2*D) - h1*r**4*(2*L*M*r-L*Q**2-L*r**2-4*M*r+2*Q**2+w**2*r**4+2*r**2)/(2*D) + r**6*CrW
    mx=a0*r**4*(2*L*M*r-L*Q**2-L*r**2+w**2*r**4)/D + 2*a1*r**3*(M*r-Q**2) + a2*r**4*D + 2*h0*Q*r**3 - d(h0)*Q*r**4 - sp.I*h1*Q*w*r**4
    Ct=(-12*M**2*t0*r**2 + 6*M**2*t1*r**3 + 24*M*Q**2*t0*r - 9*M*Q**2*t1*r**2 + 6*M*t0*r**3 - 3*M*t1*r**4 - 9*Q**4*t0 + 3*Q**4*t1*r - 9*Q**2*t0*r**2 + 3*Q**2*t1*r**3)/r**8
    tp=2*a1*Q*r**2*D + h0*r**2*(L*r**2-4*M*r+2*Q**2)/2 - d(d(h0))*r**4*D/2 - sp.I*h1*w*r**3*D - sp.I*d(h1)*w*r**4*D/2 + r**6*Ct
    VT=D/r**5*(L*r+2*M-2*Q**2/r)
    t2sol=sp.solve(sp.expand(ds(ds(t0))).subs(t2,T2s)-(VT-w**2)*t0,T2s)[0]
    p2sol=sp.solve(rp,p2)[0]
    a2sol=sp.solve(mx,a2)[0]
    p3sol=sp.solve(d(rp),p3)[0]
    ex=tp.subs(p3,p3sol)
    for _ in range(3):
        ex=ex.subs({p2:p2sol,a2:a2sol,t2:t2sol})
    return sp.cancel(sp.together(ex))
random.seed(1)
for k in range(3):
    M=sp.Rational(random.randint(3,9),4); Q=sp.Rational(random.randint(1,5),7)
    w=sp.Rational(random.randint(2,9),5); l=random.randint(2,5); L=l*(l+1)
    print(f'M={M} Q={Q} w={w} l={l}: t-phi residual =',run(M,Q,w,L),flush=True)
