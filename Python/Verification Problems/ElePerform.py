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

Espan=linspace(1e7,5e7,5)
vspan=linspace(0.01,0.49,5)
nE=len(Espan)
nv=len(vspan)
er=zeros([nE,nv])
problem=dict({'P':10,
              'OD':23,
              'h':.5,
              'eletyp':AxiSymmetricQuad4Reduced,
              'formula':1})
problem['E']=5
print(problem)
for i,E in enumerate(Espan):
    for j,v in enumerate(vspan):
        
        problem['E']=E
        problem['v']=v
        er[i,j]=C_Plate_Pressure_Clamped(**problem)
        

Espan2=zeros([nE,nv])
vspan2=zeros([nE,nv])
Espan2[:,0::1]=Espan
vspan2[:,0::1]=vspan
vspan2=vspan2.transpose()
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(Espan2, vspan2, er)
#ax.set_zlim3d(0, 1)
ax.set_xlabel('younsg')
ax.set_ylabel('poisson')
ax.set_zlabel('err')
plt.show()