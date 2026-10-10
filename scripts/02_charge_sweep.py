"""Sec. XIV C, Table III, Fig. 3: linear coefficient beta*d(omega) of the gravitational-led fundamental, l=2..6, Q/M in [0,0.98],
with and without the electromagnetic channel chi.  Writes data/sweep.json.  Runtime ~10 min.
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from coef import *
import json
Ws={2:0.3736716844-0.0889623157j,3:0.5994432884-0.0927030479j,4:0.8091783775-0.0941640312j,5:1.0123815197-0.0948788764j,6:1.2120108787-0.0952804414j}
Qs=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.85,0.9,0.95,0.98]
res={}
for l in [2,3,4,5,6]:
    wg=Ws[l]; res[l]=[]
    for Q in Qs:
        Nn=300 if Q<0.9 else 400; S=60 if Q<0.9 else 90
        wGR,c,_=coef(1.,Q,l,wg,th=0.6,S=S,N=Nn)
        _,c0,_=coef(1.,Q,l,wg,chi_on=False,th=0.6,S=S,N=Nn)
        wg=wGR
        res[l].append(dict(Q=Q,wGR=[wGR.real,wGR.imag],c=[c.real,c.imag],c_nochi=[c0.real,c0.imag]))
        print(f'l={l} Q={Q:4.2f} wGR={wGR:.7f}  beta*dw={c:.5f}  (no chi: {c0:.5f})',flush=True)
json.dump(res,open(DATA/'sweep.json','w'))
