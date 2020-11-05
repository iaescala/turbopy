from __future__ import absolute_import, division, print_function

import os
import numpy as np

# This probably doesn't work if you install the package
data_path = os.path.join(__file__, 'data')

def load_atmosphere(Teff, logg, MH, vt,
                    aFe=None, CFe=None, NFe=None,
                    rFe=None, sFe=None):
    """
    Directly loads a MARCS model within the grid with no interpolation
    """
    raise NotImplementedError

def interp_atmosphere(Teff, logg, MH, vt,
                      aFe=None, CFe=None, NFe=None,
                      rFe=None, sFe=None):
    """
    Uses the Masseron interpolator to get a new MARCS model
    """
    raise NotImplementedError

class MARCSModel(object):
    def __init__(self, fname=None):
        super(MARCSModel, self).__init__()
        assert fname is None or os.path.exists(fname)
        self.fname = fname

    @staticmethod
    def load(fname, validate=False):
        assert os.path.exists(fname), fname
        atmo = MARCSModel(fname)
        if validate: raise NotImplementedError
        return atmo
    
    def get_fname(self):
        return self.fname
    
    @property
    def Teff(self):
        return self._Teff
    @Teff.setter
    def Teff(self, x):
        self._Teff = x
    
    @property
    def logg(self):
        return self._logg
    @logg.setter
    def logg(self, x):
        self._logg = x
    
    @property
    def MH(self):
        return self._MH
    @MH.setter
    def MH(self, x):
        self._MH = x
    
    @property
    def aFe(self):
        return self._AM
    @aFe.setter
    def aFe(self, x):
        self._AM = x    
    @property
    def AM(self):
        return self._AM
    @aFe.setter
    def AM(self, x):
        self._AM = x
    
