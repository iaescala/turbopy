from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np
from astropy.table import Table

from . import utils

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
    
def _get_levels(l2, l3):
    """ Solve for orbit levels """
    s2 = l2.split()
    s3 = l3.split()
    all_levels = ["s","p","d","f","g","h","i","k"]
    levels = []
    for split in [s2, s3]:
        if len(split)==0 or split[0] not in ["LS","JJ","JK","LK"]: return "X", "X"
        configuration = split[1]
        this_levels = re.findall(r"[a-z]", configuration)
        # Going to try last 2 levels of configuration
        try:
            levels.append(all_levels.index(this_levels[-1]))
        except IndexError:
            levels.append(-1)
        try:
            levels.append(all_levels.index(this_levels[-2]))
        except IndexError:
            levels.append(-1)
    llo, llo2, lhi, lhi2 = levels
    ## Some updates from vald3line-BPz-freeformat to make the transition match best
    if llo >= 3 or lhi >= 3:
        llo, lhi = 2, 3
    elif (llo < 0) & (lhi < 0):
        return "X","X"
    elif np.abs(lhi - llo) == 1:
        pass # keep these
    else:
        ixrow, ixcol = 4*llo + llo2+1, 4*lhi + lhi2+1
        llo, lhi = _all_level_map[ixrow, ixcol]
    all_levels = ["s","p","d","f"]
    return all_levels[llo], all_levels[lhi]

def read_vald_long(fname):
    """
    This code is meant to mimic vald3line-BPz-freeformat.f
    but fix issues with molecules etc
    """
    fdampdict1 = {11: 2.0, 14: 1.3, 20: 1.8, 26: 1.4} # neutral damping, the rest are 2.5
    fdampdict2 = {20: 1.4, 38: 1.8, 56: 3.0} # ionized damping, the rest are 2.5
    def _should_skip(data):
        ## Note: this is going to need to operate on the full final array, just writing this down now
        elems, ion, tspecies, wave, loggf, expot, fdamp, gu, raddmp, levlo, levup, ehi = data[0:11]
        Zs = [utils.elem_to_Z(el) for el in elems]
        for Z in Zs:
            if Z < 3: return True
        ion = int(ion)
        if ion < 1 or ion > 2: return True
        if expot > 15: return True
        if len(Zs) == 1:
            # This removes lines w/ upper level in continuum including autoionization lines
            if ion == 1 and ehi > utils.get_ionp1(Zs[0]): return True
            if ion == 2 and ehi > utils.get_ionp2(Zs[0]): return True
        return False
    class FinishedReading(Exception):
        pass
    def _parse_chunk(fp):
        """ Parse a 4-line chunk containing info for one line """
        ## Read in four lines
        l1 = fp.readline().strip() # most of the data
        l2 = fp.readline().strip()[1:-1] # for lower level angmom
        l3 = fp.readline().strip()[1:-1] # for upper level angmom
        l4 = fp.readline().strip() # for isotopes
        s1 = [x.strip() for x in l1.split(",")]
        s4 = l4[1:-1].strip().split()
        
        try:
            specion, wave, loggf, expot, jlo, ehi, jhi, lower, upper, mean, Rad, Stark, Waals = s1[:13]
        except ValueError as e:
            raise FinishedReading
        jhi = float(jhi)
        Rad = float(Rad)
        Waals = float(Waals)
        
        elems, ion, isos = utils.identify_fullspecstr(specion, s4[-1])
        Zs = [utils.elem_to_Z(el) for el in elems]
        tspecies = utils.make_tspecies(Zs, isos)
        inttspecies = int(float(tspecies))
        
        ## Damping atomic data
        if Waals != 0:
            fdamp = Waals
        elif ion == 1 and inttspecies in fdampdict1:
            fdamp = str(fdampdict1[inttspecies])
        elif ion == 1 and inttspecies in fdampdict2:
            fdamp = str(fdampdict2[inttspecies])
        else:
            fdamp = 2.5
        
        gu = 2*jhi + 1
        if Rad > 3.0: raddmp = 10**Rad
        else: raddmp = 1.e5
        
        if len(elems) == 1: # an atom
            levlo, levup = _get_levels(l2, l3)
        else: # a molecule
            fdamp = 2.500 
            levlo, levup = "X", "X"
        
        return tspecies, ion, wave, expot, loggf, fdamp, gu, raddmp, \
            levlo, levup, ehi
    
    alldata = []
    with open(fname) as fp:
        for i in range(3): line = fp.readline()
        while True:
            try:
                out = _parse_chunk(fp)
            except FinishedReading:
                break
            alldata.append(out)
            #elems, ion, specstr, wave, loggf, dum, dum, dum, dum, levlo, levhi = out[0:11]
            #print(specstr, levlo, levhi)
        #raise
    
    return out

# I just manually did the selection rule matrix. 12x12
# The rows are lower, going from 0-1 00 01 02 through 22 for llo1,llo2
# The cols are upper, going from 0-1 00 01 02 through 22 for lhi1,lhi2
# Index is 4*lo + hi+1
_all_level_map = np.array([
    [(0,0),(0,0),(0,1),(0,2), (0,1),(0,1),(0,1),(0,1), (0,2),(0,2),(0,1),(0,2)],
    [(0,0),(0,0),(0,1),(0,0), (0,1),(0,1),(0,1),(0,1), (0,2),(0,2),(0,1),(0,2)],
    [(1,0),(1,0),(0,1),(1,2), (0,1),(0,1),(0,1),(0,1), (0,2),(0,2),(0,1),(0,2)],
    [(2,0),(0,0),(0,1),(0,0), (0,1),(0,1),(0,1),(0,1), (0,2),(0,2),(0,1),(0,2)],
    
    [(1,0),(1,0),(1,0),(1,0), (1,1),(1,0),(1,1),(1,2), (1,2),(1,2),(1,2),(1,2)],
    [(1,0),(1,0),(1,0),(1,0), (0,1),(1,0),(0,1),(1,2), (1,2),(1,2),(1,2),(1,2)],
    [(1,0),(1,0),(1,0),(1,0), (1,1),(1,0),(1,1),(1,2), (1,2),(1,2),(1,2),(1,2)],
    [(1,0),(1,0),(1,0),(1,0), (2,1),(1,0),(2,1),(1,2), (1,2),(1,2),(1,2),(1,2)],

    [(2,0),(1,0),(2,1),(2,0), (2,1),(2,1),(2,1),(2,1), (2,2),(2,0),(2,1),(2,2)],
    [(2,0),(2,0),(2,1),(2,0), (2,1),(2,1),(2,1),(2,1), (0,2),(2,2),(2,1),(2,2)],
    [(1,0),(1,0),(2,1),(1,2), (2,1),(2,1),(2,1),(2,1), (1,2),(1,2),(2,1),(1,2)],
    [(2,0),(2,0),(2,1),(2,0), (2,1),(2,1),(2,1),(2,1), (2,2),(2,2),(2,1),(2,2)]
])
