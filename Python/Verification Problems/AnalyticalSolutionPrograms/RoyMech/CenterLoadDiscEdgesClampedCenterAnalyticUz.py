# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 14:07:28 2016

@author: Alex
"""
#Roy Mech Circular Plate, Uniform Load, Edges simply supported
#Vertical displacement at the center
import math

def CenterLoadDiscEdgesClampedCenterAnalyticUr(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = P*r**2/(16*math.pi*D)
    return u_z