# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 20:53:35 2016

@author: dlontine
"""
import sys
import pdb
sys.path.insert(0, '../')
from pyfem2 import *


def myfun(a,g,**kwargs):
    d=a*g
    return d

dd=dict({'a':5,'b':3,'c':2,'g':1,'t':4})
d=myfun(**dd)
print(d)