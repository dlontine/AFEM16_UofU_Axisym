# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 23:12:34 2016

@author: Rando
"""
def testPointLoadPlate(P,v1,v2,elements,R,r):
    D = 2*R
    v = np.linspace(v1,v2,elements)
    u_z = -P/(16*pi*D) * ((3+v)/(1+v)*(R^2-r^2) + 2*r^2*log(r/R))
    return u_z