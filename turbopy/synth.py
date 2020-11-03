from __future__ import absolute_import, division, print_function

import os, sys
import subprocess
import numpy as np
from .linelists import TSLineList, get_default_linelist
from .marcs import MARCSModel, interp_atmosphere

_lpoint_max = 100000 # hardcoded into turbospectrum, we might change this

def run_synth(wmin, wmax, dwl,
              abundances=[], atmosphere=None, linelist=None,
              Teff=None, logg=None, MH=None, vt=None,
              aFe=None, CFe=None, NFe=None, rFe=None, sFe=None,
              babsma_outfname=None, bsyn_outfname=None,
              babsma_infname=None,
):
    """
    Run a turbospectrum synthesis.
    
    wmin, wmax, dwl: which wavelength range to run
    abundances: a list of Z, XFe pairs, e.g. [(6, 0.0), (8, 1.0), (20, 0.4)]
    atmosphere: either a string pointing to the MARCS model, or a MARCS model object
    linelist: 
    
    Teff, logg, MH, vt: can specify all of these 4 parameters to call interp_atmosphere
    aFe, CFe, NFe, rFe, sFe: can add these to interp_atmosphere if desired
    
    babsma_outfname: where to save the output of babsma_lu (continous opacity)
    bsyn_outfname: where to save the output of bsyn_lu (spectrum)
    
    babsma_infname: read in a precomputed continuous opacity
    """
    
    Nwl = np.ceil((wmax-wmin)/dwl)
    if Nwl > _lpoint_max:
        raise ValueError(f"Trying to synthesize {Nwl} > {_lpoint_max} wavelength points")
    
    if linelist is None:
        linelist = get_default_linelist(wmin, wmax)
    else:
        assert isinstance(linelist, TSLineList)
    
    abundances = validate_abundances(abundances)

    if atmosphere is not None:
        if isinstance(atmosphere, str):
            atmosphere = MARCSModel.load(atmosphere)
        assert isinstance(atmosphere, MARCSModel)
    else:
        assert Teff is not None, Teff
        assert logg is not None, logg
        assert MH is not None, MH
        assert vt is not None, vt
        atmosphere = interp_atmosphere(Teff, logg, MH, vt,
                                       aFe, CFe, NFe, rFe, sFe)
    
    param_file = make_parameter_file()
    
    ## Expecting to adapt from jobovy/apogee/modelspec/turbospec.py
    print("Running Turbospectrum babsma_lu")
    try:
        p = subprocess.Popen(["babsma_lu"],
                             cwd=twd,
                             stdin=subprocess.PIPE,
                             stdout=None, stderr=None)
        stdout, stderr = p.communicate()
    except subprocess.CalledProcessError:
        # Could delete the relevant files
        raise RuntimeError("Running babsma_lu failed...")
    finally:
        raise NotImplementedError
    raise NotImplementedError
        
def make_parameter_file():
    raise NotImplementedError

def validate_abundances(abundances):
    """ Input is format [(Z1, XFe1), (Z2, XFe2), ...] """
    assert isinstance(abundances, list), abundances
    Zs = []
    XFes = []
    for Z,XFe in abundances:
        Z = int(Z)
        assert (Z >= 3) & (Z <= 92), Z
        Zs.append(Z)
        
        XFe = np.round(float(XFe),3)
        assert (XFe >= -9) & (XFe <= 9), XFe
        XFes.append(XFe)
    new_abundances = []
    for Z, XFe in zip(Zs, XFes):
        new_abundances.append([Z, XFe])
    
