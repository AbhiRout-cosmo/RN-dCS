"""Sec. IX C, Eq. (Qc). At Q_c^2 = 4(L+1)M^2/(L+2)^2 the pseudoscalar potential equals the Chandrasekhar polar partner of the EM-led
branch, V_Theta == V_e^(+), for every l  =>  omega_Theta = omega_e^(+) = omega_e^(-) exactly (GR isospectrality)."""
import sympy as sp
r,M,L=sp.symbols('r M L',positive=True)
Q2=4*(L+1)*M**2/(L+2)**2; mu2=L-2
s=sp.sqrt(sp.factor(9*M**2+4*Q2*mu2)); s=sp.simplify(s); print('s =',s)
Del=r**2-2*M*r+Q2; f=Del/r**2; q2=3*M-s
F=Del/(r**3*(mu2*r+q2)); Vp=q2*f*sp.diff(F,r)+q2**2*F**2+mu2*(mu2+2)*F; Vm=-q2*f*sp.diff(F,r)+q2**2*F**2+mu2*(mu2+2)*F
VT=Del/r**5*(L*r+2*M-2*Q2/r); Ve=Del/r**5*(L*r+4*Q2/r-q2)
print('V_e^(-) (Chandrasekhar) == axial EM-led potential:',sp.simplify(Vm-Ve)==0)
print('V_Theta == V_e^(+) for generic L:',sp.simplify(VT-Vp)==0)
