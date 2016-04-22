# -*- coding: utf-8 -*-
"""
CPoint_Plate.py
Created on Fri Apr 22 16:05:20 2016

@author: dlontine

This code makes a rectilinear mesh and applies an axisymmetric quad element.
It bends the plate according to example problem ____ where there is a point
load in the center of the plate and the plate is simply supported at the end.
"""

import sys
sys.path.insert(0, '../')
from pyfem2 import *

mesh = RectilinearMesh2D(nx=20, ny=4, lx=11.5, ly=0.4)
mat = Material('Material-1', elastic={'E':1e6, 'Nu':0.2})

V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
V.ElementBlock('ElementBlock1', ALL)
V.AssignProperties('ElementBlock1', AxiSymmetricQuad4Reduced, mat)

step = V.StaticStep()
step.FixDOF(IHI)
step = V.StaticStep()
step.ConcentratedLoad(ILO, Y, -100)
step.run()
V.WriteResults()
if not os.environ.get('NOGRAPHICS'):
    V.Plot2D(show=1, deformed=1)