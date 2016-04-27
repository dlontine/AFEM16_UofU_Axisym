import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
import math
import numpy as np

# DICTIONARY:
# V = Finite element model object
# E = Young's Modulus
# v = Poisson's ratio
# P = Load applied (pressure or force- context driven)
# OD= Outside diameter of axisymmetric model
# inD= Inside diameter of axisymmetric model
# h = Thickness of model
# NinX = Number of elements in I (for Rectilinear Mesh)
# NinY = Number of elements in J (for Rectilinear Mesh)
# eletyp = Axisymmetric element type 
# X = Radial position from line of axisymmetry (also known as r)
# Y = Axial position from zero reference (aligned with axis of symmetry) (also known as z)

# All required inputs for each function should be defined as per parameters above
# All functions should have **kwargs as inputs so that we may implement dictionaries
#----------------------------------------------------------------------------#
# ---------------------- Analytical Slolutions ------------------------------#
#----------------------------------------------------------------------------#
def A_Plate_Point_Fellipa(E,v,P,OD,h,z=None,r=None,inD=None,**kwargs):
    D = E*h**3/(12*(1-v**2))
    R = OD/2
    u_r = P / (8*math.pi*D) * ((3+v)/(1+v) - 1 - 2*math.log(r/R))*r*z
    u_z = -P / (16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_r,u_z
    
def A_Thick_Infinite_Cyl(E,v,P,OD,h,X=None,Y=None,inD=None,**kwargs):
    #Function robustness items:    
    if Xcoord is None:
        Xcoord=0.0
    if Ycoord is None:
        Ycoord=0.0 
    u_z=0.0 #All z displacement is fixed due to boundary conditions (plane strain)
    #From Fellipa "Verification Problems" eq 7.2:
    a=OD/2.0
    b=inD/2.0
    num=a**2*(1+v)*(b**2+r**2*(1-2*v))
    den=E*(b**2-a**2)*r
    u_r=P*num/den
    return u_r,u_z

#Roymech solutions:

def A_Plate_Point_Clamped(E,v,P,OD,h,z=None,r=None,inD=None,**kwargs):
    D = E*h**3/(12*(1-v**2))
    u_z = P*(OD/2)**2/(16*math.pi*D)
    return u_z

#The diagram does not show a pinned support, rather a roller. Agree?
def A_Plate_Point_Pinned(E,v,P,OD,h,z=None,r=None,inD=None,**kwargs):
    RO=OD/2
    D = E*h**3/(12*(1-v**2))
    u_z = (5+v)*P*RO**4 / (64*(1+v)*D)
    return u_z

def A_Plate_Pressure_Clamped(E,v,P,OD,h,**kwargs):
    RO=OD/2.0
    D = E*h**3/(12*(1-v**2))
    u_z = P*RO**4/(64*D)
    return u_z

def A_Plate_Pressure_Pinned(E,v,P,OD,h,**kwargs):
    r=OD/2
    D = E*h**3/(12*(1-v**2))
    u_z = -(5+v)*P*r**4/(64*(1+v)*D)
    return u_z

def A_Washer_Point_Clamped(E,v,P,OD,h,inD,**kwargs):
    a = OD/2
    b = inD/2
    c = a/b
    t = h
    k = -.0016*c**6 + .0233*c**5 + -.1285*c**4 + .3072*c**3 - .2544*c**2 + .051
    u_z = k * P * a**2 / (E * t**3)
    return u_z

def A_Washer_Point_Pinned(E,v,P,OD,h,inD,**kwargs):
    a = OD/2
    b = inD/2
    c = a/b
    t = h
    k = 0.0111*c**6 - 0.1724*c**5 + 1.0195*c**4 - 2.7879*c**3 + 3.1547*c**2 -1.1484
    u_z = k * P * a**2 / E * t**3
    return u_z

def A_Washer_Pressure_Clamped(E,v,P,OD,h,inD,**kwargs):
    a = OD/2
    b = inD/2
    c = a/b
    t = h
    k = -0.0015*c**6 + 0.0230*c**5 + -0.1289*c**4 + .3166*c**3 + -0.2812*c**2 + 0.0733
    u_z = k * P * a**4 / (E * t**3)
    return u_z

def A_Washer_Pressure_Pinned(E,v,P,OD,h,inD,**kwargs):
    a = OD/2
    b = inD/2
    c = a/b
    t = h
    k = 0.01*c**6 + -.1585*c**5 + .9563*c**4 + -2.6988*c**3 + 3.2063*c**2 + -1.4443
    u_z = k * P * a**4 / (E * t**3)
    return u_z