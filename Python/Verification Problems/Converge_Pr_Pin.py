import sys
import pdb
import matplotlib
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions import *
from Definitions2 import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt

pconv=dict({'P':10,
            'OD':23,
            'h':.4,
            'formula':1,
            'E':5e7,
            'v':0.4})
nstep=7
xmax=200
ymax=int(pconv['h']*xmax/pconv['OD'])*2
ninx=linspace(50,xmax,nstep)
niny=linspace(2,ymax,nstep)

er_f=zeros(nstep)
ne_f=zeros(nstep)
for ii,nx in enumerate(ninx):
    pconv['NinX']=int(ninx[ii])
    pconv['NinY']=int(niny[ii])
    pconv['eletyp']=AxiSymmetricQuad4
    er_f[ii],ne_f[ii]=C_Plate_Pressure_Pinned(**pconv)
    
er_r=zeros(nstep)
ne_r=zeros(nstep)
for ii,nx in enumerate(ninx):
    pconv['NinX']=int(ninx[ii])
    pconv['NinY']=int(niny[ii])
    pconv['eletyp']=AxiSymmetricQuad4Reduced
    er_r[ii],ne_r[ii]=C_Plate_Pressure_Pinned(**pconv)
    
er_s=zeros(nstep)
ne_s=zeros(nstep)
for ii,nx in enumerate(ninx):
    pconv['NinX']=int(ninx[ii])
    pconv['NinY']=int(niny[ii])
    pconv['eletyp']=AxiSymmetricQuad4SelectiveReduced
    er_s[ii],ne_s[ii]=C_Plate_Pressure_Pinned(**pconv)

plt.plot(ne_f,er_f,'r',label='Full Int')
plt.plot(ne_r,er_r,'b',label='Red Int')
plt.plot(ne_s,er_s,'g',label='Sel Red Int')
plt.xlabel('Number of Elements')
plt.ylabel('% Error From Analytical')
plt.title('Convergence on Clamped Pressure Plate')
plt.legend()
plt.grid(True)
plt.show()

plt.plot(ne_f,er_f,'r',label='Full Int')
plt.plot(ne_r,er_r,'b',label='Red Int')
#plt.plot(ne_s,er_s,'g',label='Sel Red Int')
plt.xlabel('Number of Elements')
plt.ylabel('% Error From Analytical')
plt.title('Convergence on Clamped Pressure Plate')
plt.legend()
plt.grid(True)
plt.show()