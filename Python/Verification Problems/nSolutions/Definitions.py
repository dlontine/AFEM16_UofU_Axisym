# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:22:45 2016

@author: dlontine
"""
import sys
sys.path.insert(0, '../')
from pyfem2 import *

def Plate_Point_Clamped(E,v,P,eletyp=None):
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4)
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':0.2})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

Plate_Point_Clamped()

def Plate_Point_Pinned():
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4)
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':0.2})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

Plate_Point_Pinned()

def Plate_Pressure_Pinned():
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    
    kp1=nei+1         #Bottom outside edge
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4)
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':0.2})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp1)
    step = V.StaticStep()
    step.Pressure(JHI, 100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Plate_Pressure_Pinned()

def Plate_Pressure_Clamped():
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4)
    mat = Material('Material-1', elastic={'E':1e6, 'Nu':0.2})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.Pressure(JHI, 100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Plate_Pressure_Clamped()

def Washer_Point_Pinned(E,v):
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    #Define key points for use in applying loads
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    insided=5    
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4,shift=[insided,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -1000)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Washer_Point_Pinned(1e6,.1)

def Washer_Point_Clamped(E,v):
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    #Define key points for use in applying loads
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    insided=5    
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4,shift=[insided,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -1000)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Washer_Point_Clamped(1e6,.1)

def Washer_Pressure_Pinned(E,v):
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    #Define key points for use in applying loads
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    insided=5    
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4,shift=[insided,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', AxiSymmetricQuad4, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.Pressure(JHI, 100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Washer_Pressure_Pinned(1e6,.1)

def Washer_Pressure_Clamped(E,v,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    #Define key points for use in applying loads
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    insided=5    
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4,shift=[insided,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.Pressure(JHI, 100)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Washer_Pressure_Clamped(1e6,.1,AxiSymmetricQuad4Reduced)
Washer_Pressure_Clamped(1e6,.1,AxiSymmetricQuad4)
Washer_Pressure_Clamped(1e6,.1,AxiSymmetricQuad4SelectiveReduced)


def Thick_Infinite_Cyl(E,v,P,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    nej=5   #Number elements in J (Thickness)
    nei=21  #Number of elements in I (Diameter)
    #Define key points for use in applying loads
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    insided=5    
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=11.5, ly=0.4,shift=[insided,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='QuadPlate')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    V.PrescribedBC(JLO,Y)
    V.PrescribedBC(JHI,Y)
    
    #step = V.StaticStep()
    #step.FixNodes(IHI)
    step = V.StaticStep()
    step.Pressure(ILO, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
Thick_Infinite_Cyl(1e4,.3,100)