# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
#Roy Mech Circular Plate, Uniform Load, Edges simply supported
#Vertical displacement at the center
import math

def CenterLoadRingEdgesSimplySupportedCenterAnalyticUz(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.0111*c**6 - 0.1724*c**5 + 1.0195*c**4 - 2.7879*c**3 + 3.1547*c**2 -1.1484
    u_z = k * P * a**2 / E * t**3
    return u_z