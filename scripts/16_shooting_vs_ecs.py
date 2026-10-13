"""App. E: independent Frobenius/complex-ray shooting of the FULL {Psi,Psi_e,Theta} system (and the reduced system)
vs ECS: l=2, Q=0,0.6,0.9.  Agreement 1.4e-6 in beta*dw, 8e-9 in omega_0.  Writes data/shooting_vs_ecs.json.  ~4 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
from zdm_tools import RN, Ubranch, secant, zdm, Ucoupled
def Ufull(Q,l,beta):
    mu=np.sqrt((l-1)*(l+2)); L=l*(l+1)
    def U(r):
        chi=2*Q/(r*mu); k=(r-Q*Q)/r**3; ct=6*L*(L-2)*k/beta
        return np.array([[L*r-6+4*Q*Q/r, 2*Q*mu, 96*np.pi*k],[2*Q*mu, L*r+4*Q*Q/r, chi*96*np.pi*k],[ct, chi*ct, L*r+2-2*Q*Q/r]])
    return U
out={'shoot':[]}
t0=time.time()
for (l,Q,g) in ((2,0.0,0.37367168-0.08896231j),(2,0.6,0.38621747-0.08981368j),(2,0.9,0.41357050-0.08833301j)):
    Qs=max(Q,1e-9)
    sv0=RN(Qs,Ufull(Qs,l,1e300),3); w0=secant(sv0.F,g)
    cs=[]
    for B in (2e4,4e4):
        sv=RN(Qs,Ufull(Qs,l,B),3); w=secant(sv.F,w0,w0*(1+1e-6)); cs.append(B*(w-w0))
    cf=2*cs[1]-cs[0]
    s1=RN(Qs,Ubranch(Qs,l,'G'),1); w0r=secant(s1.F,g); cr=[]
    for B in (2e4,4e4):
        s2=RN(Qs,Ucoupled(Qs,l,B,'G') if Q>0 else Ufull(Qs,l,B),2 if Q>0 else 3); w=secant(s2.F,w0r,w0r*(1+1e-6)); cr.append(B*(w-w0r))
    crr=2*cr[1]-cr[0]
    kw=dict(N=300,S=60.) if Q<0.9 else dict(N=400,S=90.)
    w0e,ce,_=coef_full(1.,Q,l,g,**kw)
    out['shoot'].append(dict(l=l,Q=Q,w0_shoot=w0,w0_ecs=w0e,c_shoot_full=cf,c_shoot_red=crr,c_ecs_full=ce))
    print(f'l={l} Q={Q}: w0 shoot {w0:.10f} ecs {w0e:.10f} | c shoot(full) {cf:.6f} shoot(red) {crr:.6f} ecs(full) {ce:.6f}  %.0fs'%(time.time()-t0),flush=True)
    save('shooting_vs_ecs',out)
