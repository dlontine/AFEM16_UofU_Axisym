# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:14:51 2016

@author: Alex
"""
#Formula 7.6, Verification Problems Fellipa, pg. 7-18 (electronic page, 18)
# Displacement in the Z direction. 

def PointLoadCenterDiscAnalytic(r,z,E,v,P,R,h):
    D = E*h**3/(12*(1-v**2))
    u_z = -P/(16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_z