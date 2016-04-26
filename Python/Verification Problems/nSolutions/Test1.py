# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:31:54 2016

@author: dlontine
"""
import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions import * 



E = 1e8
v = 0.3
OD= 23
r = OD/2.0 	# Radius of plate
t = h = 0.4 	# Thickness of plate
P = E/1e5 	# Pressure applied to the plate
    
####-----FEM-----####
V=Thick_Infinite_Cyl(E,v,P,OD,t,OD/2)
zFEM=get_max_disp(V)

D = E*t**3/(12*(1-v**2)) #Flexural rigidity
zmax=(5+v)*P*r**4/(64*(1+v)*D)

####-----Exact-----####
Amax=A_Plate_Pressure_Pinned(E,v,P,r,t,0,0.1)


#E = 1e8
#v = 0.3
#OD= 23
#r = OD/2.0 	# Radius of plate
#t = h = 0.4 	# Thickness of plate
#P = E/1e5 	# Pressure applied to the plate
#    
#####-----FEM-----####
#V=Plate_Pressure_Pinned(E,v,P,OD,t)
#zFEM=get_max_disp(V)
#
#D = E*t**3/(12*(1-v**2)) #Flexural rigidity
#zmax=(5+v)*P*r**4/(64*(1+v)*D)
#
#####-----Exact-----####
#Amax=A_Plate_Pressure_Pinned(E,v,P,r,t,0,0.1)

