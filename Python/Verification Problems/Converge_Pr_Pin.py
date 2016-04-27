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

def find_convergence(Model_Comparison_Function,
                     nstep=None,xmax=None,title=None,
                     saveas1=None,saveas2=None,
                     saveas3=None,saveas4=None,**kwargs):
    if nstep is None:
        nstep=7
    if xmax is None:
        xmax=200
    if title is None:
        title='Convergence on Model'
    if saveas1 is None:
        saveas1='Convergence_A.png'
    if saveas2 is None:
        saveas2='Convergence_F.png'
    if saveas2 is None:
        saveas2='Convergence_R.png'
    if saveas2 is None:
        saveas2='Convergence_S.png'
    pconv=dict({'P':40,
                'OD':23,
                'h':.4,
                'formula':1,
                'inD':23/2,
                'E':5e7,
                'v':0.4})
    ymax=int(pconv['h']*xmax/pconv['OD'])*2
    ninx=linspace(50,xmax,nstep)
    niny=linspace(2,ymax,nstep)
    
    er_f=zeros(nstep)
    ne_f=zeros(nstep)
    for ii,nx in enumerate(ninx):
        pconv['NinX']=int(ninx[ii])
        pconv['NinY']=int(niny[ii])
        pconv['eletyp']=AxiSymmetricQuad4
        er_f[ii],ne_f[ii]=Model_Comparison_Function(**pconv)
        
    er_r=zeros(nstep)
    ne_r=zeros(nstep)
    for ii,nx in enumerate(ninx):
        pconv['NinX']=int(ninx[ii])
        pconv['NinY']=int(niny[ii])
        pconv['eletyp']=AxiSymmetricQuad4Reduced
        er_r[ii],ne_r[ii]=Model_Comparison_Function(**pconv)
        
    er_s=zeros(nstep)
    ne_s=zeros(nstep)
    for ii,nx in enumerate(ninx):
        pconv['NinX']=int(ninx[ii])
        pconv['NinY']=int(niny[ii])
        pconv['eletyp']=AxiSymmetricQuad4SelectiveReduced
        er_s[ii],ne_s[ii]=Model_Comparison_Function(**pconv)
    
    plt.plot(ne_f,er_f,marker='o', linestyle='-', color='r',label='Full Int')
    plt.plot(ne_r,er_r,marker='o', linestyle='-', color='b',label='Red Int')
    #plt.plot(ne_s,er_s,marker='o', linestyle='-', color='g',label='Sel Red Int')
    plt.xlabel('Number of Elements')
    plt.ylabel('Error From Analytical')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(saveas1)
    plt.show()
    
    
    plt.plot(ne_f,er_f,marker='o', linestyle='-', color='r',label='Full Int')
    #plt.plot(ne_r,er_r,marker='o', linestyle='-', color='b',label='Red Int')
    #plt.plot(ne_s,er_s,marker='o', linestyle='-', color='g',label='Sel Red Int')
    plt.xlabel('Number of Elements')
    plt.ylabel('Error From Analytical')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(saveas2)    
    plt.show()
    
    #plt.plot(ne_f,er_f,marker='o', linestyle='-', color='r',label='Full Int')
    plt.plot(ne_r,er_r,marker='o', linestyle='-', color='b',label='Red Int')
    #plt.plot(ne_s,er_s,marker='o', linestyle='-', color='g',label='Sel Red Int')
    plt.xlabel('Number of Elements')
    plt.ylabel('Error From Analytical')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(saveas3)    
    plt.show()
    
    #plt.plot(ne_f,er_f,marker='o', linestyle='-', color='r',label='Full Int')
    #plt.plot(ne_r,er_r,marker='o', linestyle='-', color='b',label='Red Int')
    plt.plot(ne_s,er_s,marker='o', linestyle='-', color='g',label='Sel Red Int')
    plt.xlabel('Number of Elements')
    plt.ylabel('Error From Analytical')
    plt.title(title)
    plt.legend()
    plt.grid(True)
    plt.savefig(saveas4)    
    plt.show()

nsteps=10
pdict=dict({'Model_Comparison_Function':C_Washer_Point_Clamped,
            'nstep':nsteps,
            'xmax': 350/2,
            'title':'Convergence of Clamped Point Washer',
            'saveas1':'Conv_WaPoCl_1.png',
            'saveas2':'Conv_WaPoCl_2.png',
            'saveas3':'Conv_WaPoCl_3.png',
            'saveas4':'Conv_WaPoCl_4.png'})
find_convergence(**pdict)

pdict=dict({'Model_Comparison_Function':C_Washer_Point_Pinned,
            'nstep':nsteps,
            'xmax': 350/2,
            'title':'Convergence of Pinned Point Washer',
            'saveas1':'Conv_WaPoPi_1.png',
            'saveas2':'Conv_WaPoPi_2.png',
            'saveas3':'Conv_WaPoPi_3.png',
            'saveas4':'Conv_WaPoPi_4.png'})
find_convergence(**pdict)

pdict=dict({'Model_Comparison_Function':C_Washer_Pressure_Pinned,
            'nstep':nsteps,
            'xmax': 350/2,
            'title':'Convergence of Pinned Pressure Washer',
            'saveas1':'Conv_WaPrPi_1.png',
            'saveas2':'Conv_WaPrPi_2.png',
            'saveas3':'Conv_WaPoPi_3.png',
            'saveas4':'Conv_WaPoPi_4.png'})
find_convergence(**pdict)

pdict=dict({'Model_Comparison_Function':C_Washer_Pressure_Clamped,
            'nstep':nsteps,
            'xmax': 350/2,
            'title':'Convergence of Clamped Pressure Washer',
            'saveas1':'Conv_WaPrCl_1.png',
            'saveas2':'Conv_WaPrCl_2.png',
            'saveas3':'Conv_WaPrCl_3.png',
            'saveas4':'Conv_WaPrCl_4.png'})
find_convergence(**pdict)


#####
pdict=dict({'Model_Comparison_Function':C_Plate_Point_Clamped,
            'nstep':nsteps,
            'xmax': 350,
            'title':'Convergence of Clamped Point Plate',
            'saveas1':'Conv_PlPoCl_1.png',
            'saveas2':'Conv_PlPoCl_2.png',
            'saveas3':'Conv_PlPoCl_3.png',
            'saveas4':'Conv_PlPoCl_4.png',})
find_convergence(**pdict)
#
pdict=dict({'Model_Comparison_Function':C_Plate_Point_Pinned,
            'nstep':nsteps,
            'xmax': 350,
            'title':'Convergence of Pinned Point Plate',
            'saveas1':'Conv_PlPoPi_1.png',
            'saveas2':'Conv_PlPoPi_2.png',
            'saveas3':'Conv_PlPoPi_3.png',
            'saveas4':'Conv_PlPoPi_4.png'})
find_convergence(**pdict)
#
pdict=dict({'Model_Comparison_Function':C_Plate_Pressure_Pinned,
            'nstep':nsteps,
            'xmax': 350,
            'title':'Convergence of Pinned Pressure Plate',
            'saveas1':'Conv_PlPrPi_1.png',
            'saveas2':'Conv_PlPrPi_2.png',
            'saveas3':'Conv_PlPrPi_3.png',
            'saveas4':'Conv_PlPrPi_4.png'})
find_convergence(**pdict)

pdict=dict({'Model_Comparison_Function':C_Plate_Pressure_Clamped,
            'nstep':nsteps,
            'xmax': 350,
            'title':'Convergence of Clamped Pressure Plate',
            'saveas1':'Conv_PlPrCl_1.png',
            'saveas2':'Conv_PlPrCl_2.png',
            'saveas3':'Conv_PlPrCl_3.png',
            'saveas4':'Conv_PlPrCl_4.png'})
find_convergence(**pdict)


#####

