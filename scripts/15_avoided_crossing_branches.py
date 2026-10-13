"""Sec. IX C, Fig. 3(a,b): both coupled EM-led/pseudoscalar modes through Q_c (l=2), beta=1e4,1e5,1e6, with exact branch
identification (single-field operators) and the symmetrized Theta weight of each eigenvector.  Writes data/avoided_branches.json.  ~60 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from ecs_pivot import *
l=2; N=260; S=60.; th=0.6
Qs=np.round(np.arange(0.60,0.7201,0.004),4)
out={'Q':Qs,'dec':[], 'cpl':{}}
gE,gS=0.5100-0.0985j,0.5240-0.0985j   # rough guesses near Q~0.6
prevE=prevS=None
for Q in Qs:
    ev=evals(op_full(1.,Q,l,1,coupled=False,th=th,S=S,N=N))
    if prevE is None:
        # identify EM-led and scalar at start by branch-specific single-field operators
        import branch as BR
        prevE=BR.eig_near(BR.branch_op(1.,Q,l,1,2,field='b',N=N,S=S),0.50-0.098j)
        prevS=BR.eig_near(BR.branch_op(1.,Q,l,1,2,field='s',N=N,S=S),0.53-0.098j)
    import branch as BR
    wE=BR.eig_near(BR.branch_op(1.,Q,l,1,2,field='b',N=N,S=S),prevE)   # exact branch identification
    wS=BR.eig_near(BR.branch_op(1.,Q,l,1,2,field='s',N=N,S=S),prevS)
    prevE,prevS=wE,wS; out['dec'].append([wE,wS])
    print('Q',Q,'dec',f'{wE:.6f}',f'{wS:.6f}',flush=True)
for B in (1e4,1e5,1e6):
    rows=[]
    for i,Q in enumerate(Qs):
        H=op_full(1.,Q,l,B,th=th,S=S,N=N)
        lam,V=sla.eig(H); w=np.sqrt(lam); sgn=np.where(w.real<0,-1,1); w=w*sgn
        wE,wS=out['dec'][i]
        pick=[]
        for tgt in (wE,wS):
            idx=np.argsort(abs(w-tgt))
            j=[k for k in idx if k not in pick][0]; pick.append(j)
        n1=N-1
        res=[]
        for j in pick:
            v=V[:,j]; sc_=np.sqrt(16*np.pi*B/((l-1)*l*(l+1)*(l+2))); nt=np.linalg.norm(v[2*n1:])*sc_; ng=np.linalg.norm(v[:n1]); ne=np.linalg.norm(v[n1:2*n1])
            res.append([w[j],nt**2/(ng**2+ne**2+nt**2)])
        rows.append(res)
        print('beta',B,'Q',Q,[f'{x[0]:.6f} th-frac {x[1]:.2f}' for x in res],flush=True)
    out['cpl'][B]=rows; save('avoided_branches',out)
print('DONE')
