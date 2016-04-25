# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
# Formula 7.2, Verification Problems Fellipa, pg. 7-4 (electronic page, 4)
# Displacement in the Z direction.
# Breaks down as r approaches zero. 

def PointLoadCenterDiscAnalyticUr(r,z,E,v,P,R,h,Ri):
    D = E*h**3/(12*(1-v**2))
    a = Ri
    b = R
    u_r = P * a**2 * (1+v) * (b**2 + r**2 * (1 - 2*v)) / (E * (b**2 - a**2) * r)
    return u_r