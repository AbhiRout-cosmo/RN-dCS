"""App. D, Table V: local-WKB proxy for the l=2, n=0 shift (no chi channel, local pseudoscalar Green's function):
d(w^2) = C_dCS C_Thg / (w0^2 - V_Theta) at the peak of V_RW,  d(w) = d(w^2)/(2 w0).   Q=0: beta*dw = -0.4466+0.1470i.
Qualitative proxy only: it gets the sign of Im(dw) wrong (exact: -0.7733-0.3946i, script 01)."""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import numpy as np
from scipy.optimize import minimize_scalar
from ecs import operator
import scipy.linalg as sla

M, l = 1.0, 2
L = l*(l+1)
def baseline(Q, g):
    ev = np.sqrt(sla.eigvals(operator(M, Q, l, 1e300, coupled=False, N=260)))
    ev = np.where(ev.real < 0, -ev, ev)
    return ev[np.argmin(abs(ev-g))]
g = 0.3736716844-0.0889623157j
print(" Q/M   r_peak   beta*dw (local WKB)   |beta*dw|")
for Q in [0.0, 0.2, 0.4, 0.6, 0.8]:
    D = lambda r: r*r-2*M*r+Q*Q
    VRW = lambda r: D(r)/r**5*(L*r-6*M+4*Q*Q/r)
    VT = lambda r: D(r)/r**5*(L*r+2*M-2*Q*Q/r)
    rp = minimize_scalar(lambda r: -VRW(r), bounds=(1.8, 5), method='bounded', options=dict(xatol=1e-12)).x
    w0 = baseline(Q, g); g = w0
    K = (M*rp-Q*Q)*D(rp)/rp**8
    CdCS, bCTg = 96*np.pi*K, 6*(l-1)*l*(l+1)*(l+2)*K
    dw2 = CdCS*bCTg/(w0**2-VT(rp))
    dw = dw2/(2*w0)
    print(f" {Q:.1f}   {rp:.3f}   {dw.real:+.4f}{dw.imag:+.4f}i     {abs(dw):.3f}")
