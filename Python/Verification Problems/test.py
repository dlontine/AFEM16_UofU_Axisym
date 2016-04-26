import os
import numpy as np
import subprocess
import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions import *
from Definitions2 import *

#This code utilizes the testing framework within python to run verification tests on the FEM solver that we have developed.

def test_1():
    # TEST CASE 1
    # Geometry: Flat circular plate with no holes
    # Supports: Simply supported at radius r
    # Loads:    Uniform pressure load across entire plate
    # ELEMENT TYPE: AxiSymmetricQuad4
    #
    # Objective of test: 
    # Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    problem = dict{'E':1e6,'v':0.3,}
    E = 1e8
    v = 0.3
    OD= 23
    r = OD/2.0 	# Radius of plate
    t = 0.4 	# Thickness of plate
    P = E/1e5 	# Pressure applied to the plate
    
    ####-----FEM-----####
    V=Plate_Pressure_Pinned(E,v,P,OD,h)
    zFEM=get_max_disp(V)
    
    ####-----Exact-----####
    Amax=A_Plate_Pressure_Pinned(E,v,P,r,t,0,0)
	
	
    print(Amax)
    print(zFEM)
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    assert np.allclose(zFEM,zmax,atol=1e-5)


def test_2():
    # TEST CASE 2
    # Geometry: Flat circular plate with no holes
    # Supports: Cantilever or clamped edges at radius r
    # Loads: 	Uniform pressure load across entire plate
    #
    # Objective of test: 
    # Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    E = 1e6
    v = 0.3
    
    r = 23/2.0 	# Radius of plate
    t = 0.4 	# Thickness of plate
    p = 1000 	# Pressure applied to the plate
    ####-----FEM-----####
    #Construct the FEM solution
    #Mesh
    #Loading Conditions
    #Evaluation of key solution point
    
    zFEM = 0.0001
    
    ####-----Exact-----####
    D = E*t**3/(12*(1-v**2)) #Flexural rigidity
    zmax=p*r**4/(64*D)		 # Analytical solution for this case
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zmax,atol=1e-10)

def test_3():
    # TEST CASE 3: ASTM Standard
    # Geometry: Flat circular plate with no holes
    # Supports: Simply supported at Ds
    # Loads: 	Force per circumference applied at Dl
    #
    # Objective of test: 
# Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    E = 1e6
    v = 0.3
    
    D = 23		# Diameter of the plate
    h = 0.4 	# Thickness of plate
    F = 1000 	# Force applied to the plate
    Ds= 11
    Dl= 5
    
    ####-----FEM-----####
    #Construct the FEM solution
    #Mesh- ABAQUS GENERATED
    #Loading Conditions
    #Evaluation of key solution point
    
    zFEM = 0.0001
    
    ####-----Exact-----####
    T1=3*F*(1-v**2)*Dl**2/(8*np.pi*E*h**3)
    T2=Ds**2/Dl**2
    T3=1+((1-v)*(Ds**2-Dl**2))/(2*(1+v)*D**2)
    T4=1+np.log(Ds/Dl)
    zmax=T1*(T2*T3-T4)		 # Analytical solution for this case
    
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zmax,atol=1e-10)

def test_4():
    # TEST CASE 4 ***********
    # Geometry: Flat circular plate with no holes
    # Supports: Cantilever or clamped edges at radius r
    # Loads: 	Point load at center of the plate
    #
    # Objective of test: 
    # Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    E = 1e6
    v = 0.3
    
    R = 23/2.0 	# Radius of plate
    t = 0.4 	# Thickness of plate
    P = 1000 	# Pressure applied to the plate
    
    ####-----FEM-----####
    #Construct the FEM solution
    #Mesh
    #Loading Conditions
    #Evaluation of key solution point
    
    zFEM = 0.0001
    
    ####-----Exact-----####
    # Solution of the deflection in z as a function of radius r
    # Solution of the deflection in r as a function of radius r
    D = E*t**3/(12*(1-v**2)) #Flexural rigidity
    r=1
    z=0
    ur = P/(8*np.pi*D)*((3+v)/(1+v)-1-2*np.log(r/R))*r*z
    uz = -P/(16*np.pi*D)*((3+v)/(1+v)*(R**2-r**2)+2*r**2*np.log(r/R))
    # Solution for the maximum deflection in z ***Special Case***
    zmax = -P/(16*np.pi*D)*((3+v)/(1+v)*R**2)
    
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zmax,atol=1e-10)
	
def test_5():
    # TEST CASE 4: WITH HOLE
    # Geometry: Flat circular plate with no holes
    # Supports: Cantilever or clamped edges at radius r
    # Loads: 	Point load at center of the plate
    #
    # Objective of test: 
    # Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    E = 1e6
    v = 0.3
    
    R = 23/2.0 	# Radius of plate
    Ri= 3		# Inside radius (HOLE)
    t = 0.4 	# Thickness of plate
    p = 1000 	# Pressure applied to the plate
    
    ####-----FEM-----####
    #Construct the FEM solution
    #Mesh
    #Loading Conditions
    #Evaluation of key solution point
    
    zFEM = 0.0001
    
    ####-----Exact-----####
    # Solution of the deflection in z as a function of radius r
    # Solution of the deflection in r as a function of radius r
    #ur = P/(8*pi*D)*((3+v)/(1+v)-1-2*log(r/R))*r*z
    #uz = -P/(16*pi*D)*((3+v)/(1+v)*(R**2-r**2)+2*r**2*log(r/R))
    # Solution for the maximum deflection in z ***Special Case***
    #zmax = -P/(16*pi*D)*((3+v)/(1+v)*R**2)
    
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zmax,atol=1e-10)

test_1()
test_2()
test_3()
test_4()
test_5()