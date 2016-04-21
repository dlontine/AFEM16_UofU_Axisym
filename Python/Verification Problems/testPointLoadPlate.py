# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 23:12:34 2016
@author: Rando
"""

import numpy as np
import math

def testPointLoadPlate(P,v1,v2,elements,R,r,h,E):
    v = np.linspace(v1,v2,elements)
    D = E*h**3/(12*(1-v**2))
    u_z = -P/(16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_z