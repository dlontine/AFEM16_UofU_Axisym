import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions2 import *


# DICTIONARY:
# V = Finite element model object
# E = Young's Modulus
# v = Poisson's ratio
# P = Load applied (pressure or force- context driven)
# OD= Outside diameter of axisymmetric model
# inD= Inside diameter of axisymmetric model
# h = Thickness of model
# NinX = Number of elements in I (for Rectilinear Mesh)
# NinY = Number of elements in J (for Rectilinear Mesh)
# eletyp = Axisymmetric element type 
# X = Radial position from line of axisymmetry (also known as r)
# Y = Axial position from zero reference (aligned with axis of symmetry) (also known as z)

# All required inputs for each function should be defined as per parameters above
# All functions should have **kwargs as inputs so that we may implement dictionaries

#----------------------------------------------------------------------------#
# ---------------------- Post-Processing Fxns  ------------------------------#
#----------------------------------------------------------------------------#
def get_max_disp(V,**kwargs):
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

def get_disp_pos(GP,V,**kwargs):
    Xi=array(V.mesh.coord)
    ui=array(V.steps.last.dofs.reshape(Xi.shape))
    ux = ui[GP,0]
    uy = ui[GP,1]
    X  = Xi[GP,0]
    Y  = Xi[GP,1]
    return ux,uy,X,Y

def get_all_disp_pos(V,**kwargs):
    nnodes=V.numnod
    data = array(zeros((nnodes,4)))
    for ii in range(0,nnodes):
        data[ii,:]=get_disp_pos(ii,V)
    return data

#----------------------------------------------------------------------------#
# ------------------------- FEM Problems Setup ------------------------------#
#----------------------------------------------------------------------------#
def Plate_Point_Pinned(E,v,P,OD,h,
                       NinX=None,NinY=None,eletyp=None,formula=None,
                       **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 100 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if formula is None:
        formula=1
    R=OD/2.0
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePointPinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat, formulation=formula)
    
    #V.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    #step.PinNodes(kp2)
    #step.FixNodes(kp2)
    #step.PrescribedBC(kp2,Y)
    step.PrescribedBC(kp2,Zr)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V

def Plate_Point_Clamped(E,v,P,OD,h,
                        NinX=None,NinY=None,eletyp=None,formula=None,
                        **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX=100 #Number of elements in I (Diameter)
    if NinY is None:
        NinY=4 #Number elements in J (Thickness)
    if formula is None:
        formula=1
    R=OD/2.0
    #Find Keypoints of interest
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePointClamped')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat, formulation=formula)
    
    step = V.StaticStep()
    step.FixNodes(IHI)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V

def Plate_Pressure_Pinned(E,v,P,OD,h,
                          NinX=None,NinY=None,eletyp=None,formula=None,
                          **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 100 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if formula is None:
        formula=1
    R=OD/2.0
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=R, ly=h)
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='PlatePressurePinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat, formulation=formula)
    
    #V.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.PinNodes(kp2)
    #step.FixNodes(kp2)
    #step.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V


def Plate_Pressure_Clamped(E,v,P,OD,h,
                           NinX=None,NinY=None,eletyp=None,formula=None,
                           **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 100 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if formula is None:
        formula=1    
    R=OD/2.0
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=R, ly=h)
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
    return V

def Washer_Point_Pinned(E,v,P,OD,h,
                        NinX=None,NinY=None,eletyp=None,inD=None,formula=None,
                        **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 50 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPointPinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    #V.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.PinNodes(kp2)
    #step.FixNodes(kp2)
    #step.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.ConcentratedLoad(kp1, Y, -P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V

def Washer_Point_Clamped(E,v,P,OD,h,
                         NinX=None,NinY=None,eletyp=None,inD=None,formula=None,
                         **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 50 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
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
    return V

def Washer_Pressure_Pinned(E,v,P,OD,h,
                           NinX=None,NinY=None,eletyp=None,inD=None,formula=None,
                           **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 50 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    V = FiniteElementModel(mesh=mesh, jobid='WasherPressurePinned')
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    
    #V.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.PinNodes(kp2)
    #step.FixNodes(kp2)
    #step.PrescribedBC(kp2,Y)
    step = V.StaticStep()
    step.Pressure(JHI, P)
    step.run()
    V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V

def Washer_Pressure_Clamped(E,v,P,OD,h,
                            NinX=None,NinY=None,eletyp=None,inD=None,formula=None,
                            **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 50 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
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
    return V
    
def Universal_Plate(E,v,P,OD,h,Plate=None,Pin=None,Point=None,
                            NinX=None,NinY=None,eletyp=None,inD=None,
                            formula=None,job=None,**kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 50 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 4 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    if Plate is None:
        Plate = True
    if Pin is None:
        Pin = True
    if Point is None:
        Point = True
    if job is None:
        job = 'AxiSymPlate'
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mat = Material('Material-1', elastic={'E':E, 'Nu':v})
    
    #Case Structure for Plate vs Washer
    if Plate:
        #If Plate        
        mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=R, ly=h)
    else:
        #If Washer
        mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
    
    V = FiniteElementModel(mesh=mesh, jobid=job)
    V.ElementBlock('ElementBlock1', ALL)
    V.AssignProperties('ElementBlock1', eletyp, mat)
    step = V.StaticStep()
    # Case Structure for Pined vs Fixed Outside Diameter
    if Pin:
        step.PrescribedBC(kp2,Y)
    else:
        step.FixNodes(IHI)
    
    step = V.StaticStep()
    # Case Structure for Pressure vs Point loading
    if Point:
        step.ConcentratedLoad(kp1, Y, -P)
    else:
        step.Pressure(JHI, P)
    
    step.run()
    #V.WriteResults()
    if not os.environ.get('NOGRAPHICS'):
        V.Plot2D(show=1, deformed=1)
    return V


def Thick_Infinite_Cyl(E,v,P,OD,h,
                       NinX=None,NinY=None,eletyp=None,inD=None,formula=None,
                       **kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        NinX = 60 #Number of elements in I (Diameter)
    if NinY is None:
        NinY = 10 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    if formula is None:
        formula=1
    R  = OD/2.0
    Ri = inD/2.0
    w  = R-Ri
    kp1=NinY*(NinX+1)+1 #Central point
    kp2=NinX+1         #Bottom outside edge
    mesh = RectilinearMesh2D(nx=NinX, ny=NinY, lx=w, ly=h,shift=[Ri,0])
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
    return V

#----------------------------------------------------------------------------#
# ------------------------- Element Evaluation ------------------------------#
#----------------------------------------------------------------------------#

def Comp_Analysis(E,v,P,OD,h,inD=None,
                  NinX=None,NinY=None,eletyp=None,formula=None,
                  Model_Comparison_Function=None,A_Mod_Comp_Fun=None,**kwargs):
    d=dict({'E':E,'v':v,'P':P,'OD':OD,'h':h,'inD':inD,'NinX':NinX,'NinY':NinY,
            'eletyp':eletyp,'formula':formula})
    V=Model_Comparison_Function(**d)
    zFEM = get_max_disp(V)
    zANA = A_Mod_Comp_Fun(**d)
    err=(zFEM-zANA)#/zANA*100.
    err2=(zFEM-zANA)/zANA*100.
    print(Model_Comparison_Function)
    print(eletyp)
    print('Numerical Sol is: ',zFEM)
    print('Analytical Sol is:',zANA)
    print('Absolute Error is:',err)
    print('Percent Error is: ',err2)
    nele=V.numele
    return err,err2,nele


#def C_Plate_Point_Pinned(E,v,P,OD,h,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Plate_Point_Pinned(E,v,P,OD,h,NinX,NinY,eletyp,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Plate_Point_Pinned(E,v,P,OD,h)
#    
#    err=(zFEM-zANA)#/zANA*100.
#    err2=(zFEM-zANA)/zANA*100.
#    print('Numerical Sol is: ',zFEM)
#    print('Analytical Sol is:',zANA)
#    print('Absolute Error is:',err)
#    print('Percent Error is: ',err2)
#    nele=V.numele
#    return err,err2,nele
#    
#def C_Plate_Point_Clamped(E,v,P,OD,h,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Plate_Point_Clamped(E,v,P,OD,h,NinX,NinY,eletyp,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Plate_Point_Clamped(E,v,P,OD,h)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#
#def C_Plate_Pressure_Pinned(E,v,P,OD,h,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Plate_Pressure_Pinned(E,v,P,OD,h,NinX,NinY,eletyp,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Plate_Pressure_Pinned(E,v,P,OD,h)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#    
#def C_Plate_Pressure_Clamped(E,v,P,OD,h,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Plate_Pressure_Clamped(E,v,P,OD,h,NinX,NinY,eletyp,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Plate_Pressure_Clamped(E,v,P,OD,h)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#
#def C_Washer_Point_Pinned(E,v,P,OD,h,inD,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Washer_Point_Pinned(E,v,P,OD,h,NinX,NinY,eletyp,inD,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Washer_Point_Pinned(E,v,P,OD,h,inD)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#    
#def C_Washer_Point_Clamped(E,v,P,OD,h,inD,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Washer_Point_Clamped(E,v,P,OD,h,NinX,NinY,eletyp,inD,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Washer_Point_Clamped(E,v,P,OD,h,inD)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#
#def C_Washer_Pressure_Pinned(E,v,P,OD,h,inD,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Washer_Pressure_Pinned(E,v,P,OD,h,NinX,NinY,eletyp,inD,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Washer_Pressure_Pinned(E,v,P,OD,h,inD)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele
#    
#def C_Washer_Pressure_Clamped(E,v,P,OD,h,inD,
#                         NinX=None,NinY=None,eletyp=None,
#                         formula=None,**kwargs):
#    V = Washer_Pressure_Clamped(E,v,P,OD,h,NinX,NinY,eletyp,inD,formula)
#    zFEM = get_max_disp(V)
#    zANA = A_Washer_Pressure_Clamped(E,v,P,OD,h,inD)
#    print(zFEM)
#    print(zANA)
#    err=(zFEM-zANA)#/zANA*100.
#    print(err)
#    nele=V.numele
#    return err,nele