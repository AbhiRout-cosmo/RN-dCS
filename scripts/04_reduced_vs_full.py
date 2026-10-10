"""Sec. VIII A / App. E: the branch-reduced system {Z_g, Theta}, Eq. (59), reproduces the full three-field O(zeta) coefficient
(1e-7 at l=2, 1e-4 relative at l=6).  Needs data/sweep.json (script 02).
"""
import sys, pathlib
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DATA = ROOT / "data"
from branch import *
import json
S=json.load(open(DATA/'sweep.json'))
print('reduced {Z1,Theta} vs full 3x3 (grav-led linear coefficient):')
for l in ['2','6']:
    for x in S[l]:
        if x['Q'] in (0.3,0.7,0.95):
            N=300 if x['Q']<0.9 else 400; Sb=60 if x['Q']<0.9 else 90
            w0,c=coef_b(1.,x['Q'],int(l),1,complex(*x['wGR']),N=N,S=Sb)
            print(f"  l={l} Q={x['Q']}: reduced {c:.5f}  full {complex(*x['c']):.5f}  |diff|={abs(c-complex(*x['c'])):.1e}",flush=True)
