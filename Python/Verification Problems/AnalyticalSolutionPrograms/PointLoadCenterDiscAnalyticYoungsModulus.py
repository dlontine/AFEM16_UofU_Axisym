# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 23:12:34 2016
@author: Rando
"""
#Formula 7.6, Verification Problems Fellipa, pg. 7-18 (electronic page, 18)
# Displacement in the Z direction. 
# Variation in Youngs Modulus ratio.

import numpy as np
import math

def PointLoadCenterDiscAnalyticYoungsModulus(P,E1,E2,elements,R,r,h,v):
    E = np.linspace(E1,E2,elements)
    D = E*h**3/(12*(1-v**2))
    u_z = -P/(16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_z