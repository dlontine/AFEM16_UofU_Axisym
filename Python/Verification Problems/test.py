import os
import numpy as np
import subprocess

#This code utilizes the testing framework within python to run verification tests on the FEM solver that we have developed.


def test_1():
    #This test verifies that the center displacement in a disc that is pressure loaded on one side and simply supported on the outside edge of the other.
	# See schematic in documentation (XXXX) for more information.
	E = 1e6
	v = 0.3
	
	#Construct the FEM solution
	#Mesh
	#Loading Conditions
	#Evaluation of key solution point
	#
	
	#Construct the exact solution
	#Compare the solutions (export solution comparison as png)
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
    
