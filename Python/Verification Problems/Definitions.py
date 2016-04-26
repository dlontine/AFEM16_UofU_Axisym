import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *


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
def Plate_Point_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,**kwargs):
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
    return V

def Plate_Point_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,**kwargs):
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
    return V

def Plate_Pressure_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,**kwargs):
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
    return V


def Plate_Pressure_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,**kwargs):
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
    return V

def Washer_Point_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,inD=None,**kwargs):
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
    return V

def Washer_Point_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,inD=None,**kwargs):
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
    return V

def Washer_Pressure_Pinned(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,inD=None,**kwargs):
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
    return V

def Washer_Pressure_Clamped(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,inD=None,**kwargs):
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
    return V


def Thick_Infinite_Cyl(E,v,P,OD,h,NinX=None,NinY=None,eletyp=None,inD=None,**kwargs):
    if eletyp is None:
        eletyp = AxiSymmetricQuad4
    if NinX is None:
        nei = 20 #Number of elements in I (Diameter)
    if NinY is None:
        nej = 10 #Number elements in J (Thickness)
    if inD is None:
        inD = OD/2.0 
    R  = OD/2.0
    Ri = inD/2.0
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
    return V


#----------------------------------------------------------------------------#
# ---------------------- Analytical Slolutions ------------------------------#
#----------------------------------------------------------------------------#
def A_Plate_Point_Fellipa_R(r,z,E,v,P,R,h,Ri):
    D = E*h**3/(12*(1-v**2))
    a = Ri
    b = R
    u_r = P * a**2 * (1+v) * (b**2 + r**2 * (1 - 2*v)) / (E * (b**2 - a**2) * r)
    return u_r
    
def A_Plate_Point_Fellipa_Z(r,z,E,v,P,R,h):
    D = E*h**3/(12*(1-v**2))
    u_z = -P/(16*math.pi*D) * ((3+v)/(1+v)*(R**2-r**2) + 2*r**2*math.log(r/R))
    return u_z
    

def PointLoadCenterDiscAnalyticUr(r,z,E,v,P,R,h,Ri):
    D = E*h**3/(12*(1-v**2))
    a = Ri
    b = R
    u_r = P * a**2 * (1+v) * (b**2 + r**2 * (1 - 2*v)) / (E * (b**2 - a**2) * r)
    return u_r
    
def A_Thick_Infinite_Cyl(E,v,P,OD,h,X=None,Y=None,inD=None,**kwargs):
    #Function robustness items:    
    if Xcoord is None:
        Xcoord=0.0
    if Ycoord is None:
        Ycoord=0.0 
    u_z=0.0 #All z displacement is fixed due to boundary conditions (plane strain)
    #From Fellipa "Verification Problems" eq 7.2:
    a=OD/2.0
    b=inD/2.0
    num=a**2*(1+v)*(b**2+r**2*(1-2*v))
    den=E*(b**2-a**2)*r
    u_r=P*num/den
    return u_r,u_z

#Roymech solutions:

def A_Plate_Point_Clamped(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = P*r**2/(16*math.pi*D)
    return u_z

def A_Plate_Point_Pinned(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = (5+v)*P*RO**4 / (64*(1+v)*D)
    return u_z

def A_Plate_Pressure_Clamped(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = P*RO**4/(64*D)
    return u_z

def A_Plate_Pressure_Pinned(E,v,P,RO,h,z,r,RI):
    D = E*h**3/(12*(1-v**2))
    u_z = (5+v)*P*r**4/(64*(1+v)*D)
    return u_z

def A_Washer_Point_Clamped(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -.0016*c**6 + .0233*c**5 + -.1285*c**4 + .3072*c**3 - .2544*c**2 + .051
    u_z = k * P * a**2 / E * t**3
    return u_z

def A_Washer_Point_Pinned(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.0111*c**6 - 0.1724*c**5 + 1.0195*c**4 - 2.7879*c**3 + 3.1547*c**2 -1.1484
    u_z = k * P * a**2 / E * t**3
    return u_z

def A_Washer_Pressure_Clamped(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = -0.0015*c**6 + 0.0230*c**5 + -0.1289*c**4 + .3166*c**3 + -0.2812*c**2 + 0.0733
    u_z = k * P * a**4 / (E * t**3)
    return u_z

def A_Washer_Pressure_Pinned(E,v,P,RO,h,z,r,RI):
    a = RO
    b = RI
    c = a/b
    t = h
    k = 0.01*c**6 + -.1585*c**5 + .9563*c**4 + -2.6988*c**3 + 3.2063*c**2 + -1.4443
    u_z = k * P * a**4 / (E * t**3)
    return u_z





































#EOF