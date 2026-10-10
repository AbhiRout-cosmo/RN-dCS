"""Fig. 1: branch-resolved results, (a) l-scaling, (b) Q_c(l), (c) EM-led shift. Reads data/c3_sweeps.json; writes figures/fig_branches.pdf.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.optimize import brentq
plt.rcParams.update({'font.size':9,'axes.linewidth':0.6,'axes.edgecolor':'#555','axes.labelcolor':'#222','xtick.color':'#444','ytick.color':'#444','font.family':'serif'})
C=['#2a78d6','#eb6834','#1baf7a']; MK=['o','s','^']
R=json.load(open(DATA/'c3_sweeps.json'))
fig,ax=plt.subplots(1,3,figsize=(7.2,2.55),constrained_layout=True)
a=ax[0]
for i,Q in enumerate(['0.0','0.3','0.6']):
    d=np.array(R['lscal'][Q]); a.loglog(d[:,0],-d[:,3],'-',color=C[i],marker=MK[i],ms=3.5,mfc='white',mew=1,lw=1.5,label=f'$Q/M={float(Q):.1f}$')
ll=np.array([6,12]); a.loglog(ll,0.09*ll**3,':',color='#666',lw=1); a.text(5.0,0.09*5.0**3*2.2,'$\\propto l^3$',fontsize=7.5,color='#444')
a.loglog(ll,0.42*ll**2,':',color='#666',lw=1); a.text(8.3,0.42*8.3**2*0.45,'$\\propto l^2$',fontsize=7.5,color='#444')
a.set_xlabel('$l$'); a.set_ylabel(r'$-\mathrm{Re}\,\beta\,\delta\omega_1$'); a.set_title('(a) gravitational-led, large $l$',fontsize=9,loc='left')
a.legend(fontsize=7,frameon=False,loc='upper left'); a.grid(alpha=0.25,lw=0.4,which='both'); a.set_xticks([2,3,4,6,8,12]); a.set_xticklabels(['2','3','4','6','8','12'])
b=ax[1]; ls=np.array([int(k) for k in R['Qc']]); qc=np.array([R['Qc'][str(k)] for k in ls])
lf=np.linspace(2,14,200)
def crit(l):
    mu=np.sqrt((l-1)*(l+2))
    return brentq(lambda Q:2-6*Q*Q/((3+np.sqrt(9-8*Q*Q))/2)+3-np.sqrt(9+4*Q*Q*mu*mu),1e-3,0.999)
b.plot(lf,[crit(x) for x in lf],'-',color=C[0],lw=1.5,label='criterion (eikonal)')
b.plot(lf,2/np.sqrt((lf-1)*(lf+2)),'--',color=C[1],lw=1.2,label=r'$2M/\sqrt{(l-1)(l+2)}$')
b.plot(ls,qc,'o',color='#222',ms=4,mfc='white',mew=1,label='exact (ECS)')
b.set_xlabel('$l$'); b.set_ylabel('$Q_c/M$'); b.set_title('(b) EM-led / pseudoscalar crossing',fontsize=9,loc='left'); b.legend(fontsize=7,frameon=False); b.grid(alpha=0.25,lw=0.4)
c=ax[2]
for i,l in enumerate(['2','3','4']):
    d=np.array(R['em'][l]); y=d[:,3].copy(); k=np.where(np.sign(y[:-1])!=np.sign(y[1:]))[0]; x=np.insert(d[:,0],k+1,np.nan); y=np.insert(y,k+1,np.nan); c.plot(x,y,'-',color=C[i],marker=MK[i],ms=3,mfc='white',mew=0.9,lw=1.4,label=f'$l={l}$')
    c.axvline(R['Qc'][l],color=C[i],lw=0.7,ls=':')
c.set_yscale('symlog',linthresh=1); c.axhline(0,color='#999',lw=0.5)
c.set_xlabel('$Q/M$'); c.set_ylabel(r'$\mathrm{Re}\,\beta\,\delta\omega_2$'); c.set_title('(c) EM-led branch',fontsize=9,loc='left'); c.legend(fontsize=7,frameon=False,loc='center right'); c.grid(alpha=0.25,lw=0.4)
fig.savefig(ROOT/'figures'/'fig_branches.pdf')
