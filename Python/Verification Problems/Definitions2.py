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

def A_Plate_Point_Clamped(E,v,P,OD,h,**kwargs):
    #Done by Derek
    #From Roark's handbook p492
    W=P
    a=OD/2
    D = E*h**3/(12*(1-v**2))
    y_max=-W*a**2/(16*pi*D)
    return y_max
    
def A_Plate_Point_Pinned(E,v,P,OD,h,**kwargs):
    #Done by Derek    
    #This is from the Roark's handbook p491
    W=P
    a=OD/2
    D = E*h**3/(12*(1-v**2))
    y_max= -W*a**2*(3+v)/(16*pi*D*(1+v))
    return y_max

def A_Plate_Pressure_Clamped(E,v,P,OD,h,**kwargs):
    #Done by Derek
    #From Roark's p488, 10b
    ###### FEM AND ANALYTICAL MATCH #####
    D = E*h**3/(12*(1-v**2))
    a = OD/2
    y_max=-P*a**4/(64*D)
    return y_max
    
def A_Plate_Pressure_Pinned(E,v,P,OD,h,**kwargs):
    #Done by Derek
    #Roark's handbook p488, 10a
    D = E*h**3/(12*(1-v**2))
    a = OD/2
    y_max=-P*a**4*(5+v)/(64*D*(1+v))
    return y_max
    

def A_Washer_Point_Clamped(E,v,P,OD,h,inD,**kwargs):
    #Done by Derek
    #Roark's handbook 461, 1e
    a=OD/2
    b=inD/2
    ro=b
    D = E*h**3/(12*(1-v**2))
    C1=(1+v)/2*(b/a)*log(a/b)+(1-v)/4*(a/b-b/a)
    C4=1/2*((1+v)*(b/a)+(1-v)*(a/b))
    L3=(ro/(4*a))*(((ro/a)**2+1)*log(a/ro)+(ro/a)**2-1)
    L6=(ro/(4*a))*((ro/a)**2-1+2*log(a/ro))
    y_max = -P*a**3/D * (C1*L6/C4-L3)
    return y_max


def A_Washer_Point_Pinned(E,v,P,OD,h,inD,**kwargs):
    #Done by Derek
    #Roark's handbook 459, 1
    a=OD/2
    b=inD/2
    ro=b
    D = E*h**3/(12*(1-v**2))
    C4=1/2*((1+v)*(b/a)+(1-v)*(a/b))
    C7=1/2*(1-v**2)*(a/b-b/a)
    L6=(ro/(4*a))*((ro/a)**2-1+2*log(a/ro))
    L9=(ro/a)*((1+v)/2*log(a/ro)+(1+v)/4*(1-(ro/a)**2))
    y_max = (-P*a**2/D)*(C4*L9/C7-L6)
    return y_max


def A_Washer_Pressure_Clamped(E,v,P,OD,h,inD,**kwargs):
    #Done by Derek
    #Roark's handbook p465,2e
    a=OD/2
    b=inD/2
    ro=b
    D = E*h**3/(12*(1-v**2))
    C1 =(1+v)/2*(b/a)*log(a/b)+(1+v)/4*(a/b-b/a)
    C4 =1/2*((1+v)*(b/a)+(1-v)*(a/b))
    L11=1/64*(1+4*(ro/a)**2-5*(ro/a)**4-4*(ro/a)**2*(2+(ro/a)**2)*log(a/ro))
    L14=1/16*(1-(ro/a)**4-4*(ro/a)**2*log(a/ro))
    y_max = (-P*a**4)/D*(C1*L14/C4-L11)
    return y_max

def A_Washer_Pressure_Pinned(E,v,P,OD,h,inD,**kwargs):
    #Done by Derek
    #Roark's handbook p464, 2a
    a=OD/2
    b=inD/2
    ro=b
    C1 =(1+v)/2*(b/a)*log(a/b)+(1+v)/4*(a/b-b/a)
    C7 =1/2*(1-v**2)*(a/b-b/a)
    L11=1/64*(1+4*(ro/a)**2-5*(ro/a)**4-4*(ro/a)**2*(2+(ro/a)**2)*log(a/ro))
    L17=1/4*(1-(1-v)/4*(1-(ro/a)**4)-(ro/a)**2*(1+(1+v)*log(a/ro)))
    y_max=-P*a**4/D*(C1*L17/C7-L11)
    return y_max