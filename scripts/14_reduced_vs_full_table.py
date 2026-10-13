"""App. E: reduced {Z_g,Theta} vs full {Psi,Psi_e,Theta} for l=2,4,6, Q=0.2-0.95, plus the beta^-2 test at l=6, Q=0.2.
Writes data/reduced_vs_full_table.json.  ~20 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
import branch as BR
C={'A':(2,0.0,0.37367168-0.08896231j,dict()),
   'B':(2,0.6,0.38621747-0.08981368j,dict()),
   'C':(2,0.95,0.42169276-0.08665867j,dict(S=90.,N=400)),
   'D':(6,0.6,1.26604130-0.09617508j,dict()),
   'E':(2,0.3,0.46986497-0.09582639j,dict())}     # E: EM-led branch, off resonance
out={'red':[]}
t0=time.time()
S=json.load(open(DATA/'sweep.json'))
for l in (2,4,6):
    for x in S[str(l)]:
        if x['Q'] in (0.2,0.4,0.6,0.8,0.9,0.95):
            kw=dict(N=300,S=60.) if x['Q']<0.9 else dict(N=400,S=90.)
            g=complex(*x['wGR'])
            w0f,cf,_=coef_full(1.,x['Q'],l,g,**kw)
            w0r,cr=BR.coef_b(1.,x['Q'],l,1,g,**kw)
            out['red'].append([l,x['Q'],cf,cr,abs(cf-cr)/abs(cf)])
            print('red',l,x['Q'],f'{cf:.6f}',f'{cr:.6f}',f'{abs(cf-cr)/abs(cf):.1e}',flush=True); save('reduced_vs_full_table',out)

g=complex(*[x for x in S['6'] if x['Q']==0.2][0]['wGR']); out['beta_test']=[]
for B in (2e4,2e5,2e6):
    _,cf,_=coef_full(1.,0.2,6,g,B=B); _,cr=BR.coef_b(1.,0.2,6,1,g,B=B,N=300,S=60.)
    out['beta_test'].append([B,cf,cr,abs(cf-cr)/abs(cf)]); print(f'beta-test B={B:.0e}: rel diff {abs(cf-cr)/abs(cf):.1e}',flush=True)
save('reduced_vs_full_table',out)
