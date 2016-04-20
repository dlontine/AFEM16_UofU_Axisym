from numpy import *
from .isop2_4 import CSDIsoParametricQuad4 as BaseElement
# --------------------------------------------------------------------------- #
# --------------------- BILINEAR PLANE STRAIN ELEMENT ----------------------- #
# --------------------------------------------------------------------------- #
class AxiSymmetricQuad4Reduced(BaseElement):
    ndir = 3
    nshr = 1
    integration = 4
    elefab = {'formulation': 1}
    gaussw = array([4.])
    gaussp = array([[0.,0.]])
    hourglass_control = True
    #HOURGLASS CONTROL PARAMETERS
    hglassp = array([[0.,0.,]])
    hglassv = array([[1.,-1.,1.,-1.]])
    #REST
    @property
    def formulation(self):
        return self.axisymmetric
    @formulation.setter
    def formulation(self, arg):
        assert arg in (0, 1, 2)
        self.axisymmetric = arg
    def bmatrix(self, dN, N, xi, *args):
        rp = dot(N, self.xc[:,0])
        B = zeros((4, 8))
        B[0, 0::2] = B[3, 1::2] = dN[0, :]
        B[1, 1::2] = B[3, 0::2] = dN[1, :]
        B[2, 0::2] = N / rp
        return B
