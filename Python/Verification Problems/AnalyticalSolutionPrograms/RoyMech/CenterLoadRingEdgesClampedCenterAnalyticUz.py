# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
#Roy Mech Circular Plate, Uniform Load, Edges simply supported
#Vertical displacement at the center
import math

def CenterLoadRingEdgesClampedCenterAnalyticUz(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -.0016*c**6 + .0233*c**5 + -.1285*c**4 + .3072*c**3 - .2544*c**2 + .051
    u_z = k * P * a**2 / E * t**3
    return u_z