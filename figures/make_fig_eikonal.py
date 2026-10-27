"""Fig. 2 (approach to the eikonal laws): data from scripts 05 and 18."""
import json, pathlib, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; DATA = ROOT / "data"; FIG = ROOT / "figures"
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm'})
c=json.load(open(DATA/'c3_sweeps.json'))['lscal']; D=json.load(open(DATA/'eikonal_high_l.json'))
def asym(Q,l):
    if Q==0: return -4*np.sqrt(3)*np.pi/243*l**3
    rph=(3+np.sqrt(9-8*Q*Q))/2; Dl=rph*rph-2*rph+Q*Q; f=Dl/rph**2; Kp=(rph-Q*Q)*Dl/rph**8
    return -72*np.pi*l**2/Q*Kp**2*rph**6/(np.sqrt(f)*Dl)
fig,ax=plt.subplots(figsize=(3.4,2.5)); out={}
for Qk,mk in (('0.0','o'),('0.3','s'),('0.6','^')):
    Q=float(Qk)
    L=[r[0] for r in c[Qk] if r[0]<12]+[r[0] for r in D[Qk]]
    R=[r[3]/asym(Q,r[0]) for r in c[Qk] if r[0]<12]+[complex(*r[2]).real/asym(Q,r[0]) for r in D[Qk]]
    x=np.array(L,float); y=np.array(R); out[Qk]=list(zip(L,R))
    lab='$Q=0$: exact / $(l^3)$' if Q==0 else f'$Q={Q}M$: exact / ($l^2/Q$ law)'
    ax.plot(1/x,y,mk,ms=3.5,label=lab)
    sel=x>=12; p=np.polyfit(1/x[sel],y[sel],1); xx=np.linspace(0,0.1,20); ax.plot(xx,np.polyval(p,xx),'k-',lw=0.5)
    print(Qk,'ratio -> %.4f as 1/l->0 ; slope %.3f'%(p[1],p[0]))
ax.axhline(1,color='gray',lw=0.6,ls='--'); ax.set_xlabel('$1/l$'); ax.set_ylabel('$\\mathrm{Re}\\,\\delta\\omega$ / asymptotic law')
ax.set_xlim(0,0.52); ax.legend(fontsize=6.5,frameon=False); fig.tight_layout(); fig.savefig(FIG/'fig_R3_eikonal.pdf')

