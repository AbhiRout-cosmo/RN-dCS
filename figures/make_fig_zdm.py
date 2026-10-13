"""Fig. 4 (ZDM scaling): data from script 17."""
import json, pathlib, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; DATA = ROOT / "data"; FIG = ROOT / "figures"
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm'})
C=json.load(open(DATA/'zdm_scaling.json'))['zdm']; rows=[r for r in C if abs(r[3][0])<1e-6]
sig=np.array([r[0] for r in rows]); b=np.array([r[1] for r in rows]); kap=np.array([r[2] for r in rows])
w0=np.array([complex(*r[3]) for r in rows]); c=np.array([complex(*r[4]) for r in rows]); rel=(c/w0).real
fig,ax=plt.subplots(1,2,figsize=(5.6,2.4))
ax[0].loglog(b**2,rel,'o',ms=3.5,label='Frobenius--shooting')
p=np.polyfit(np.log(b),np.log(rel),1); ax[0].loglog(b**2,np.exp(np.polyval(p,np.log(b))),'k-',lw=0.7,label=rf'fit: $\propto b^{{{p[0]:.2f}}}$ (local exp. $\to2$)')
ax[0].set_xlabel('$b^2=(M^2-Q^2)/M^2$'); ax[0].set_ylabel('$\\beta\\,\\delta\\omega/\\omega_0$'); ax[0].legend(fontsize=7,frameon=False)
Cs=-c.imag/(kap*b*b); q=np.polyfit(sig,Cs,2); ss=np.linspace(0,0.05,50)
ax[1].plot(sig,Cs,'o',ms=3.5); ax[1].plot(ss,np.polyval(q,ss),'k-',lw=0.7); ax[1].plot([0],[np.polyval(q,0)],'r*',ms=7)
ax[1].set_xlabel('$\\sigma=(r_+-r_-)/r_+$'); ax[1].set_ylabel('$C(\\sigma)=-\\beta\\,\\mathrm{Im}\\,\\delta\\omega/(\\kappa_+b^2)$'); ax[1].set_title(f'$C_0={np.polyval(q,0):.1f}$',fontsize=9)
fig.tight_layout(); fig.savefig(FIG/'fig_R2_zdm.pdf')

