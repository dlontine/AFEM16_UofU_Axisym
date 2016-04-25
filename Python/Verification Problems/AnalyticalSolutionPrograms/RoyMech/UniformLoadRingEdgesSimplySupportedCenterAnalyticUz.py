# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
#Roy Mech Circular Plate, Uniform Load, Edges simply supported
#Vertical displacement at the center
import math

def UniformLoadRingEdgesSimplySupportedCenterAnalyticUz(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.01*c**6 + -.1585*c**5 + .9563*c**4 + -2.6988*c**3 + 3.2063*c**2 + -1.4443
    u_z = k * P * a**4 / (E * t**3)
    return u_z