"""Fig. 3: gravitational-led shift vs charge, (a) l=2..6 normalized to Q=0, (b) l=2 with/without chi. Reads data/sweep.json; writes figures/fig_shift_ecs.pdf.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
import json, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.linewidth':0.6,'xtick.major.width':0.6,'ytick.major.width':0.6,
    'axes.edgecolor':'#555','axes.labelcolor':'#222','xtick.color':'#444','ytick.color':'#444','font.family':'serif'})
col={2:'#2a78d6',3:'#eb6834',4:'#1baf7a',5:'#eda100',6:'#e87ba4'}; mk={2:'o',3:'s',4:'^',5:'D',6:'v'}
S=json.load(open(DATA/'sweep.json'))
fig,ax=plt.subplots(1,2,figsize=(7.0,2.7),constrained_layout=True)
# (a) grav-led, normalized Re
a=ax[0]
for l in [2,3,4,5,6]:
    d=S[str(l)]; Q=[x['Q'] for x in d]; c=np.array([complex(*x['c']) for x in d])
    a.plot(Q,c.real/c[0].real,'-',color=col[l],lw=1.5,marker=mk[l],ms=3.5,mfc='white',mew=1)
    a.annotate(f'$l={l}$',(Q[-1],c[-1].real/c[0].real),xytext=(3,0),textcoords='offset points',fontsize=7.5,color='#222',va='center')
a.set_xlabel('$Q/M$'); a.set_ylabel(r'$\mathrm{Re}\,\delta\omega(Q)\,/\,\mathrm{Re}\,\delta\omega(0)$')
a.set_title('(a) gravitational-led branch',fontsize=9,loc='left'); a.set_xlim(0,1.08); a.grid(alpha=0.25,lw=0.4)
# (b) l=2 Re & Im with/without chi
b=ax[1]; d=S['2']; Q=[x['Q'] for x in d]
c=np.array([complex(*x['c']) for x in d]); c0=np.array([complex(*x['c_nochi']) for x in d])
b.plot(Q,c.real,'-',color=col[2],lw=1.5,marker='o',ms=3.5,mfc='white',mew=1,label=r'Re, full')
b.plot(Q,c0.real,'--',color=col[2],lw=1.2,label=r'Re, $\chi=0$')
b.plot(Q,c.imag,'-',color=col[3],lw=1.5,marker='s',ms=3.5,mfc='white',mew=1,label=r'Im, full')
b.plot(Q,c0.imag,'--',color=col[3],lw=1.2,label=r'Im, $\chi=0$')
b.axhline(0,color='#999',lw=0.5)
b.set_xlabel('$Q/M$'); b.set_ylabel(r'$\beta\,\delta\omega$   ($M=\alpha=1$)')
b.set_title('(b) $l=2$, with and without $\\chi$',fontsize=9,loc='left'); b.legend(fontsize=6.5,frameon=False,loc='lower left'); b.grid(alpha=0.25,lw=0.4)
fig.savefig(ROOT/'figures'/'fig_shift_ecs.pdf')
