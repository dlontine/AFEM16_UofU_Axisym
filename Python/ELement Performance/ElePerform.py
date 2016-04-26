import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *
from Definitions import *
from Definitions2 import *

Espan=linspace(1e3,1e8,1e3)
nE=length(Espan)
nv=length(vspan)
er=zeros([nE,nv])