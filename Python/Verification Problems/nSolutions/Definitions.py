# -*- coding: utf-8 -*-
"""
Created on Sat Apr 23 13:22:45 2016

@author: dlontine
"""
import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *

def get_max_disp(V):
    uy_max=0
    uy_max2=0
    Xi=array(V.mesh.coord)
    ui=array(V.steps.last.dofs.reshape(Xi.shape))
    for ym in ui[:,1]:
        aym=abs(ym)
        if uy_max < aym:
            uy_max=aym
            uy_max2=ym
    return uy_max2

def get_disp_pos(GP,V):
    Xi=array(V.mesh.coord)
    ui=array(V.steps.last.dofs.reshape(Xi.shape))
    ux = ui[GP,0]
    uy = ui[GP,1]
    X  = Xi[GP,0]
    Y  = Xi[GP,1]
    return ux,uy,X,Y

def Plate_Point_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 41 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    R=OD/2.0
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePointPinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    #F=File('PlatePointPinned.exo')
    ym=get_max_disp(V)
    return V,ym

def Plate_Point_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei=41 #Number of elements in I (Diameter)
    if NinY is None:
        nej=4 #Number elements in J (Thickness)
    R=OD/2.0
    #Find Keypoints of interest
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePointClamped')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

def Plate_Pressure_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 41 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    R=OD/2.0
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePressurePinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)


def Plate_Pressure_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 41 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    R=OD/2.0
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePointPinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

def Washer_Point_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,InsideD=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    if InsideD is None:
        InsideD = OD/2.0 
    R  = OD/2.0
    Ri = InsideD/2.0
    w  = R-Ri
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPointPinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

def Washer_Point_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,InsideD=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    if InsideD is None:
        InsideD = OD/2.0 
    R  = OD/2.0
    Ri = InsideD/2.0
    w  = R-Ri
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPointClamped')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

def Washer_Pressure_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,InsideD=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    if InsideD is None:
        InsideD = OD/2.0 
    R  = OD/2.0
    Ri = InsideD/2.0
    w  = R-Ri
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPressurePinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.PinNodes(kp2)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)

def Washer_Pressure_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,InsideD=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 4 #Number elements in J (Thickness)
    if InsideD is None:
        InsideD = OD/2.0 
    R  = OD/2.0
    Ri = InsideD/2.0
    w  = R-Ri
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPressurePinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)


def Thick_Infinite_Cyl(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,InsideD=None):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 10 #Number elements in J (Thickness)
    if InsideD is None:
        InsideD = OD/2.0 
    R  = OD/2.0
    Ri = InsideD/2.0
    w  = R-Ri
    kp1=nej*(nei+1)+1 #Central point
    kp2=nei+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=nei, ny=nej, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='InfiniteCylinder')
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


#E,v,P,OD,h,NinX=None,NinY=None,eletyp=None
V,ym=Plate_Point_Pinned(1e6,.2,100,23,.4)
print(ym)
a=array(V.mesh.coord)
b=array(V.steps.last.dofs.reshape(a.shape))
Plate_Point_Clamped(1e6,.2,100,23,.4)
Plate_Pressure_Pinned(1e6,.2,100,23,.4)
Plate_Pressure_Clamped(1e6,.2,100,23,.4)
Washer_Point_Pinned(1e6,.2,100,23,.4)
Washer_Point_Clamped(1e6,.2,100,23,.4)
Washer_Pressure_Pinned(1e6,.2,100,23,.4)
Washer_Pressure_Clamped(1e6,.2,100,23,.4)
Thick_Infinite_Cyl(1e6,.2,100,23,10)