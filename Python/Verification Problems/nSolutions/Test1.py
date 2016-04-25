# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 16:31:54 2016

@author: dlontine
"""
import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions import * 

datF={'E':1e6,'v':.3,'P':100,'OD':23,'h':.4,'eletyp':AxiSymmetricQuad4}
datR={'E':1e6,'v':.3,'P':100,'OD':23,'h':.4,'eletyp':AxiSymmetricQuad4Reduced}
datS={'E':1e6,'v':.3,'P':100,'OD':23,'h':.4,'eletyp':AxiSymmetricQuad4SelectiveReduced}
V=Plate_Point_Clamped(**datF)
h00=get_max_disp(V)
V=Plate_Point_Clamped(**datR)
h01=get_max_disp(V)
V=Plate_Point_Clamped(**datS)
h02=get_max_disp(V)

print('hi')

V0=Plate_Point_Clamped(1e6,.3,100,23,.4)
h0=get_max_disp(V0)
V1=Plate_Point_Clamped(1e6,.3,100,23,.4,eletyp=AxiSymmetricQuad4Reduced)
h1=get_max_disp(V1)
V2=Plate_Point_Clamped(1e6,.3,100,23,.4,eletyp=AxiSymmetricQuad4SelectiveReduced)
h2=get_max_disp(V2)