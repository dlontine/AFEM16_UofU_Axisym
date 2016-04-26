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
    problem = dict({'E':1e6,
                    'v':0.3,
                    'P': 100,
                    'OD':23,
                    'h' : .4,
                    'eletyp':AxiSymmetricQuad4,
                    'formula':1})
    
    ####-----FEM-----####
    V = Plate_Point_Clamped(**problem)
    zFEM = get_max_disp(V)
    
    ####-----Exact-----####
    zANA = A_Plate_Point_Clamped(**problem)
    
    ####-----Error-----####
    err=(zFEM-zANA)/zANA*100.
    #assert np.allclose(zFEM,zANA,atol=2)
	
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
    problem = dict({'E':1e6,
                    'v':0.3,
                    'P': 100,
                    'OD':23,
                    'h' : .4,
                    'eletyp':AxiSymmetricQuad4,
                    'formula':1})
    
    ####-----FEM-----####
    V = Washer_Point_Pinned(**problem)
    zFEM = get_max_disp(V)
    
    ####-----Exact-----####
    zANA = A_Washer_Point_Pinned(**problem)
    
    ####-----Error-----####
    err=(zFEM-zANA)/zANA*100.
    #assert np.allclose(zFEM,zANA,atol=2)

test_1()
test_2()
test_3()
test_4()
test_5()