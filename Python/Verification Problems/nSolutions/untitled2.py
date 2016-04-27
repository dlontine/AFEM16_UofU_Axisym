# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 21:21:04 2016

@author: dlontine
"""

def adder(a,b,**kwargs):
    return a+b
def subtr(a,b,**kwargs):
    return a-b
def multr(a,c,**kwargs):
    return a*c

def funtunr(functional,a,b):
    d=functional(a,b)
    return d

print(funtunr(multr,2,3))
print(funtunr(adder,2,3))
print(funtunr(subtr,2,3))