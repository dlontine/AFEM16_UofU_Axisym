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
    # Verify that element in development produces appropriate maximum 
    # deflection at the center of the plate. 
    # See schematic in documentation for more information.
    
    #### Problem Setup ####
    problem = dict({'E':1e6,
                    'v':0.3,
                    'P': 10,
                    'OD':23,
                    'h' : .4,
                    'eletyp':AxiSymmetricQuad4,
                    'formula':1})
    
    ####-----FEM-----####
    V=Plate_Pressure_Pinned(**problem)
    zFEM=get_max_disp(V)
    ####-----Exact-----####
    zANA=A_Plate_Pressure_Pinned(**problem)
    
    print(zFEM)
    print(zANA)
    err=(zFEM-zANA)/zANA*100.
    print(err)
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zANA,atol=1e-5)

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
    problem = dict({'E':1e6,
                    'v':0.3,
                    'P': 10,
                    'OD':23,
                    'h' : .4,
                    'eletyp':AxiSymmetricQuad4,
                    'formula':1})
    
    ####-----FEM-----####
    V=Plate_Pressure_Clamped(**problem)
    zFEM=get_max_disp(V)
    ####-----Exact-----####
    zANA=-A_Plate_Pressure_Clamped(**problem)
    
    print(zFEM)
    print(zANA)
    
    err=(zFEM-zANA)/zANA*100.
    print(err)
    #Compare the solutions
    #Assert a tolerable error to pass/fail test
    #assert np.allclose(zFEM,zmax,atol=1e-10)

def test_3():
    # TEST CASE 3: Plate Point Pinned
    # Geometry: Flat circular plate with no holes
    # Supports: Simply supported at Ds
    # Loads: 	Force per circumference applied at Dl
    #
    # Objective of test: 
    # Verify that element in development produces appropriate maximum deflection at the center of the plate. 
    # See schematic in documentation (XXXX) for more information.
    
    #### Problem Setup ####
    problem = dict({'E':1e6,
                    'v':0.3,
                    'P': 100,
                    'OD':23,
                    'h' : .4,
                    'eletyp':AxiSymmetricQuad4,
                    'formula':1})
    
    ####-----FEM-----####
    V = Plate_Point_Pinned(**problem)
    zFEM = get_max_disp(V)
    
    ####-----Exact-----####
    zANA = -A_Plate_Point_Pinned(**problem)
    
    print(zFEM)
    print(zANA)
    err=(zFEM-zANA)/zANA*100.
    print(err)
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