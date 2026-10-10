"""Linear-in-1/beta coefficient beta*d(omega) by Richardson extrapolation (beta = B, 2B) on the ECS operator."""
from ecs import *
import numpy as np
def coef(M,Q,l,wguess,chi_on=True,**kw):
    """linear coefficient beta*dw as beta->inf: Richardson on beta=B,2B"""
    wGR=qnm(M,Q,l,1e300,wguess,coupled=False,**kw)
    B=2e4
    c1=B*(qnm(M,Q,l,B,wGR,chi_on=chi_on,**kw)-wGR)
    c2=2*B*(qnm(M,Q,l,2*B,wGR,chi_on=chi_on,**kw)-wGR)
    return wGR, 2*c2-c1, abs(c2-c1)
