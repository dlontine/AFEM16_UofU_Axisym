# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
#Roy Mech Circular Plate, Uniform Load, Edges simply supported
#Vertical displacement at the center
import math

def UniformLoadRingEdgesClampedCenterAnalyticUz(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -0.0015*c**6 + 0.0230*c**5 + -0.1289*c**4 + .3166*c**3 + -0.2812*c**2 + 0.0733
    u_z = k * P * a**4 / (E * t**3)
    return u_z