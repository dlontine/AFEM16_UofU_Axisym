import os
import numpy as np
import subprocess
import bisect #needed for test7
from matmodlab import *



def test_1():
    #This test checks that the modulus of elasticity is correct
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(0.01, 0, 0), descriptors='ESS', frames=500)

    sxx,exx=mps.get('S.XX','E.XX')

    E = sxx[1:] / exx[1:]

    E = sxx[-1]/exx[-1]
    assert np.allclose(E/1e6, 1, atol=1e-3)
    
def test_2():
    #This test checks that the modulus of elasticity is equal in compression and tension
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+5
    Y1 = 0.00
    m = 1.00
    #Run for tension (ensuring not to reach yield strength)
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(0.01, 0, 0), descriptors='ESS', frames=500)
    
    sxx = mps.get('S.XX')
    exx = mps.get('E.XX')
    
    Et=sxx[1:]/exx[1:]
    Et=sxx[-1]/exx[-1]
    
    #Run for compression (ensuring not to reach yield strength)
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(-0.01, 0, 0), descriptors='ESS', frames=500)
    
    sxx = mps.get('S.XX')
    exx = mps.get('E.XX')
    
    Ec=sxx[1:]/exx[1:]
    Ec=sxx[-1]/exx[-1]
    
    assert np.allclose(Et/1e6,Ec/1e6,atol=1e-4)
    
def test_3():
    #This test ensures that the von-mises model (perfectly plastic) yields at the appropriate stress in uniaxial compression
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(-0.10, 0.0, 0), descriptors='ESS', frames=500)
    
    sxx=mps.get('S.XX')
    exx=mps.get('E.XX')
    sxxmin=np.amin(sxx)
    assert np.allclose(sxxmin/1e4,-Y0/1e4,atol=1e-3)
    
def test_4():
    #This test ensures that the von-mises model (perfectly plastic) yields at the appropriate stress in uniaxial tension
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(0.10, 0, 0), descriptors='ESS', frames=500)
    
    sxx=mps.get('S.XX')
    exx=mps.get('E.XX')
    sxxmax=np.amax(sxx)
    sfail=Y0
    
    assert np.allclose(sxxmax/1e4,Y0/1e4,atol=1e-3)
    
def test_5():
    #This test ensures that the von-mises model (perfectly plastic) yields at the appropriate stress in biaxial compression
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(-0.10, -0.10, 0), descriptors='EES', frames=500)
        
    sxx=mps.get('S.XX')
    exx=mps.get('E.XX')
    sxxmin=np.amin(sxx)
    
    assert np.allclose(sxxmin/1e4,-Y0/1e4,atol=1e-3)
    
def test_6():
     #This test ensures that the von-mises model (perfectly plastic) yields at the appropriate stress in biaxial tension
    E = 1e+6
    nu = 0.2
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(0.10, 0.10, 0), descriptors='EES', frames=500)

    sxx=mps.get('S.XX')
    syy=mps.get('S.YY')
    szz=mps.get('S.ZZ')
    sxxm=np.amax(sxx)
    assert np.allclose(szz/1e4,0,atol=1e-3)
    assert np.allclose(sxx/1e4,syy/1e4,atol=1e-3)
    assert np.allclose(sxxm/1e4,Y0/1e4,atol=1e-3)
    
def test_7():
    #This test ensures that the von-mises model (perfectly plastic) yields at the appropriate stress in uniaxial strain (compression)
    E = 1e+6
    nu = 0.3
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(-0.10, 0, 0), descriptors='EEE', frames=500)
    
    sxx=mps.get('S.XX')
    szz=mps.get('S.ZZ')
    
    K=E/(3*(1-2*nu))	#Bulk modulus
    G=E/(2*(1+nu))	#Shear modulus
    alpha=(3*K-2*G)/(4*G+3*K)	#Parameter to help determine yield stress
    beta=nu/(1-nu)	#Linear elastic ratio between lateral stress and axial stress
    betap=beta+1e-6	#Tolerance above the ratio
    
    beta2=szz/sxx #Ratio between the lateral stress and axial stress (starts off as beta)
    index=bisect.bisect(beta2,betap)-1 #grab the first index that is above beta_p
    k=(1-alpha)*sxx[index] #Yield stress for SXX
    
    assert np.allclose(k/1e4,-Y0/1e4,atol=1e-3)
    

def test_8():
    #Does not yield in hydrostatic compression at very high pressure (on the order of 30% strain)
    E = 1e+6
    nu = 0.33
    Y0 = 1.0e+4
    Y1 = 0.00
    m = 1.00
    mps = MaterialPointSimulator('job')
    mps.Material('uplastic', parameters={'E':E, 'nu':nu,'Y0':Y0, 'Y1':Y1, 'm':m})
    mps.MixedStep(components=(-0.30, -0.30, -0.30), descriptors='EEE', frames=500)
    
    sxx=mps.get('S.XX')
    syy=mps.get('S.YY')
    szz=mps.get('S.ZZ')
    sxxmax=np.amin(sxx)
    assert np.allclose(sxx,syy,atol=10)#is sxx same as syy
    assert np.allclose(sxx,szz,atol=10)#is sxx same as szz
    assert np.allclose(szz,syy,atol=10)#is szz same as syy
    assert np.greater(-sxxmax,Y0) #assert that sxx is greater that Y0