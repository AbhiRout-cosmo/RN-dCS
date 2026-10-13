"""Fig. 3 (avoided crossing at Q_c): data from scripts 07 and 15."""
import json, pathlib, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
ROOT = pathlib.Path(__file__).resolve().parents[1]; DATA = ROOT / "data"; FIG = ROOT / "figures"
plt.rcParams.update({'font.size':9,'font.family':'serif','mathtext.fontset':'cm'})
B=json.load(open(DATA/'avoided_branches.json')); E=json.load(open(DATA/'avoided.json'))
Q=np.array(B['Q']); dec=np.array([[complex(*a),complex(*b)] for a,b in B['dec']])
Qc=2*np.sqrt(7)/8
fig,ax=plt.subplots(1,3,figsize=(7.0,2.4))
ax[0].plot(Q,dec[:,0].real,'k--',lw=1,label='EM-led ($\\beta\\to\\infty$)'); ax[0].plot(Q,dec[:,1].real,'k:',lw=1,label='scalar ($\\beta\\to\\infty$)')
cols={'10000.0':'C0','100000.0':'C1','1000000.0':'C2'}
for k,rows in B['cpl'].items():
    w=np.array([[complex(*r[0][0]),complex(*r[1][0])] for r in rows]); fr=np.array([[r[0][1],r[1][1]] for r in rows])
    lo=np.minimum(w[:,0].real,w[:,1].real); hi=np.maximum(w[:,0].real,w[:,1].real)
    lab='$\\beta=10^{%d}$'%round(np.log10(float(k)))
    ax[0].plot(Q,lo,color=cols[k],lw=1.2,label=lab); ax[0].plot(Q,hi,color=cols[k],lw=1.2)
    # Theta fraction of the upper (higher Re) mode
    up=np.where(w[:,0].real>w[:,1].real,fr[:,0],fr[:,1])
    ax[1].plot(Q,up,color=cols[k],lw=1.2,label=lab)
for a in ax[:2]: a.axvline(Qc,color='gray',lw=0.6,ls='-.')
ax[0].set_xlabel('$Q/M$'); ax[0].set_ylabel('$M\\,\\mathrm{Re}\\,\\omega$'); ax[0].legend(fontsize=6,frameon=False); ax[0].set_title('(a) $l=2$ branches',fontsize=9)
ax[1].set_xlabel('$Q/M$'); ax[1].set_ylabel('$\\Theta$ weight of upper mode'); ax[1].set_ylim(-0.02,1.02); ax[1].set_title('(b) character exchange',fontsize=9)
for l,mk in (('2','o'),('3','s')):
    rows=E[l]['rows']; b=np.array([r[0] for r in rows]); sp=np.array([r[3] for r in rows])/np.sqrt(b)
    ax[2].loglog(b,sp,mk,ms=3.5,label=f'$l={l}$, $Q=Q_c$')
    ax[2].loglog(b,sp[2]*np.sqrt(b[2]/b),'k-',lw=0.6)
ax[2].set_xlabel('$\\beta$ ($M=\\alpha=1$)'); ax[2].set_ylabel('$M|\\omega_+-\\omega_-|$'); ax[2].legend(fontsize=7,frameon=False); ax[2].set_title('(c) splitting $\\propto\\beta^{-1/2}$',fontsize=9)
fig.tight_layout(); fig.savefig(FIG/'fig_R1_avoided.pdf')
