"""Sec. XIII-XIV and App. E: ECS gates (Schwarzschild RW, RN vs Hatsuda-Kimura series) and the neutral-limit coefficient
beta*d(omega) = -0.7733-0.3946i (l=2, M=alpha=1), with its independence of (theta, S, N).  Runtime ~1 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from coef import *
M=1.; wk=[0.3736716844-0.0889623157j,0.02581767285-0.00282403214j,0.02518778870+0.00020532453j,-0.004748170246+0.002508402108j,0.01557265014+0.00041287974j]
w0=0.3736716844-0.0889623157j
print('Gate 1: Schwarzschild l=2 RW:', qnm(M,0.,2,1e300,w0,coupled=False,N=260), ' ref', w0)
for Q in [0.1,0.3]:
    w=qnm(M,Q,2,1e300,w0,th=0.6,S=60,N=260,coupled=False); ser=sum(c*Q**(2*k) for k,c in enumerate(wk))
    print(f'Gate 2: RN Q={Q}: ECS {w:.9f}  Hatsuda-Kimura series {ser:.9f}  |diff| {abs(w-ser):.1e}')
print('Gate 3 / Eq. (82): Q=0 coefficient beta*dw (target -0.7733-0.3946i) and grid independence')

w0=0.3736716844-0.0889623157j
for kw in [dict(th=0.6,S=60,N=200),dict(th=0.6,S=60,N=300),dict(th=0.5,S=80,N=340),dict(th=0.7,S=50,N=300)]:
    wGR,c,e=coef(1.,0.,2,w0,**kw); print('Q=0',kw,f'beta*dw={c:.6f}  (Richardson spread {e:.1e})',flush=True)
for kw in [dict(th=0.6,S=60,N=300),dict(th=0.5,S=80,N=340)]:
    wGR,c,e=coef(1.,0.5,2,0.3817-0.0896j,**kw); print('Q=0.5',kw,f'wGR={wGR:.8f} beta*dw={c:.6f}',flush=True)
