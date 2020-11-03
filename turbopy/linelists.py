from __future__ import absolute_import, division, print_function

import os
import numpy as np
from astropy.table import Table

# This probably doesn't work if you install the package
data_path = os.path.join(os.path.basename(__file__), 'data')

def get_default_linelist(wmin, wmax, species_to_skip=[]):
    """
    Returns default linelists between wmin and wmax (air wavelengths, angstroms)
    """
    return TSLineList()

class TSLineList(object):
    def __init__(self, fname=None):
        super(TSLineList, self).__init__()
        assert fname is None or os.path.exists(fname)
        self.fname = fname
    
    @staticmethod
    def load(fname, validate=False):
        assert os.path.exists(fname), fname
        ll = TSLineList(fname)
        if validate: raise NotImplementedError
        return ll
    
    def get_fname(self):
        return self.fname
    
