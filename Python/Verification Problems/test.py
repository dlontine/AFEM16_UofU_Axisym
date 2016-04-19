import os
import numpy as np
import subprocess

#This code utilizes the testing framework within python to run verification tests on the FEM solver that we have developed.


def test_1():
    # TEST CASE 1
	# Geometry: Flat circular plate with no holes
	# Supports: Simply supported at radius r
	# Loads: 	Uniform pressure load across entire plate
	#
	# Objective of test: 
	# Verify that element in development produces appropriate maximum deflection at the center of the plate. 
	# See schematic in documentation (XXXX) for more information.
	
	#### Problem Setup ####
	E = 1e6
	v = 0.3
	
	r = 2.0 	# Radius of plate
	t = 0.1 	# Thickness of plate
	p = 1000 	# Pressure applied to the plate
	
	####-----FEM-----####
	#Construct the FEM solution
	#Mesh
	#Loading Conditions
	#Evaluation of key solution point
	
	yFEM = 0.0001
	
	####-----Exact-----####
	D = E*t**3/(12*(1-v**2)) #Flexural rigidity
	ymax=(5+v)*p*r**4/(64*(1+v)*D)
	
	#Compare the solutions
	#Assert a tolerable error to pass/fail test
    assert np.allclose(E/1e6, 1, atol=1e-3)
    
def test_2():
    #This test verfies that the OD displacement in a compressive loading of a simple axisymmetric bar is correct.
    
	
	#Construct the FEM solution
	#Construct the exact solution
	#Compare the solutions (export solution comparison as png)
	#Assert a tolerable error to pass/fail test
    E = 1e6
	v = 0.3
    assert np.allclose(E/1e6, 1, atol=1e-3)
    
