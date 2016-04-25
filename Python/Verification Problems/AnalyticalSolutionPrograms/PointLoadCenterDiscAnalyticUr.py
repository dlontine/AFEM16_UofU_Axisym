# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
# Formula 7.7, Verification Problems Fellipa, pg. 7-18 (electronic page, 18)
# Displacement in the Z direction.
# Breaks down as r approaches zero. 

def PointLoadCenterDiscAnalyticUr(r,z,E,v,P,R,h):
    D = E*h**3/(12*(1-v**2))
    u_r = P/(8*math.pi*D) * ((3+v)/(1+v) - 1 - 2*math.log(r/R))*r*z
    return u_r