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


p1=dict({'E':1e6,'P':100,'OD':23.0,'h':.4,'inD':5,'X':6,'Z':0,'v':.2,})
p2=dict({'E':1e6,'P':10,'OD':23.0,'h':.4,'inD':5,'X':6,'Z':0,'v':.2,})
    
####-----FEM-----####
print('CP1')
#V=Thick_Infinite_Cyl(**p1)
#ur,uz=A_Thick_Infinite_Cyl(**p1)

V=Plate_Point_Pinned(**p1)
V=Plate_Point_Clamped(**p1)
V=Plate_Pressure_Pinned(**p2)
V=Plate_Pressure_Clamped(**p2)
V=Washer_Point_Pinned(**p1)
V=Washer_Point_Clamped(**p1)
V=Washer_Pressure_Pinned(**p2)
V=Washer_Pressure_Clamped(**p2)

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

