from __future__ import absolute_import, division, print_function

import os, sys, shutil
import tempfile
import subprocess

import numpy as np
from .linelists import TSLineList, get_default_linelist
from .marcs import MARCSModel, interp_atmosphere

from . import utils

_lpoint_max = 100000 # hardcoded into turbospectrum, we might change this
_ERASESTR= "                                                                                "
def run_synth(wmin, wmax, dwl, *args,
              linelist=None,
              atmosphere=None,
              Teff=None, logg=None, MH=None, vt=None,
              aFe=None, CFe=None, NFe=None, rFe=None, sFe=None,
              modelopac=None,
              outfname=None, twd=None, verbose=False,
              costheta=1.0,isotopes={}
):
    """
    Run a turbospectrum synthesis.
    Based heavily on https://github.com/jobovy/apogee/blob/master/apogee/modelspec/turbospec.py
    From Nov 2020
    
    INPUT ARGUMENTS:
       wmin, wmax, dwl: which wavelength range and step to run
       lists with abundances: e.g. (6, 0.0), (8, 1.0), (20, 0.4)
          [Atomic number1,diff1]
          [Atomic number2,diff2]
          ...
          [Atomic numberM,diffM]
    
    SYNTHEIS KEYWORDS:
       isotopes= ('solar') use 'solar' or 'arcturus' isotope ratios; 
          can also be a dictionary with isotope ratios (e.g., isotopes= {'6.012':'0.9375','6.013':'0.0625'})
       costheta= (1.) cosine of the viewing angle
    
    ATMOSPHERE KEYWORDS:
       atmosphere: either a string pointing to the MARCS model file, or a MARCS model object
       
       Teff, logg, MH, vt: can specify all of these 4 parameters to call interp_atmosphere
       aFe, CFe, NFe, rFe, sFe: can add these to interp_atmosphere if desired
       
       modelopac= (None) 
                  (a) if set to an existing filename: assume babsma_lu has already been run and use this continuous opacity in bsyn_lu
                  (b) if set to a non-existing filename: store the continuous opacity in this file

    LINELIST KEYWORDS:
          air= (True) if True, perform the synthesis in air wavelengths (affects the default Hlinelist, nothing else; output is in air if air, vacuum otherwise); set to False at your own risk, as Turbospectrum expects the linelist in air wavelengths!)
          Hlinelist= (None) Hydrogen linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist; if None, then we first search for the Hlinedata.vac in the APOGEE linelist directory (if air=False) or we use the internal Turbospectrum Hlinelist (if air=True)
       linelist= (None) molecular and atomic linelists to use; can be set to the path of a linelist file or to the name of an APOGEE linelist, or lists of such files; if a single filename is given, the code will first search for files with extensions '.atoms', '.molec' or that start with 'turboatoms.' and 'turbomolec.'
    
    
    OUTPUT:
       (wavelengths,cont-norm. spectrum, spectrum (nwave))
       if keyword outfname is set to a path:
          save the output of bsyn_lu (spectrum) to outfname
    
    """
    
    Nwl = np.ceil((wmax-wmin)/dwl)
    if Nwl > _lpoint_max:
        raise ValueError(f"Trying to synthesize {Nwl} > {_lpoint_max} wavelength points")
    
    ## working directory
    if twd is not None:
        twd = tempfile.mkdtemp(dir=os.getcwd()+"/tmp")
    
    ## Linelist
    if linelist is None:
        linelist = get_default_linelist(wmin, wmax)
    else:
        assert isinstance(linelist, TSLineList)
    linelistfilenames = [linelist.get_fname()]
    rmLinelists = False
    # Link the Turbospectrum DATA directory
    os.symlink(os.getenv('TURBODATA'),os.path.join(twd,'DATA'))
    
    ## TODO: HLinelist
    
    ## Isotopes: TODO
    #if isinstance(isotopes,str) and isotopes.lower() == 'solar':
    #    isotopes= {}
    #elif isinstance(isotopes,str) and isotopes.lower() == 'arcturus':
    #    isotopes= {'6.012':'0.9375',
    #               '6.013':'0.0625'}
    #elif not isinstance(isotopes,dict):
    #    raise ValueError("'isotopes=' input not understood, should be 'solar', 'arcturus', or a dictionary")
    
    ## Stellar atmosphere
    if atmosphere is not None:
        # The MARCS models need you to set vt separately
        assert vt is not None, vt
        if isinstance(atmosphere, str):
            atmosphere = MARCSModel.load(atmosphere)
        assert isinstance(atmosphere, MARCSModel)
        atmosphere.vt = vt
        Teff, logg, MH, aFe = atmosphere.Teff, atmosphere.logg, \
            atmosphere.MH, atmosphere.AM
    else:
        assert Teff is not None, Teff
        assert logg is not None, logg
        assert MH is not None, MH
        assert vt is not None, vt
        atmosphere = interp_atmosphere(Teff, logg, MH, vt,
                                       aFe, CFe, NFe, rFe, sFe)
        atmosphere.writeto(os.path.join(twd, 'atm.mod'))
    modelfilename = atmosphere.get_fname()
    
    ## Abundances
    abundances = validate_abundances(list(args), atmosphere.MH)
    
    ## TODO
    ## THIS PART IS DIRECTLY COPIED from Jo's code
    ## Needs to be updated and tested etc
    
    if modelopac is None or \
            (isinstance(modelopac,str) and not os.path.exists(modelopac)):
        # Now write the script file for babsma_lu
        scriptfilename= os.path.join(twd,'babsma.par')
        modelopacname= os.path.join(twd,'mopac')
        _write_script(scriptfilename,
                      wmin,wmax,dwl,
                      None,
                      modelfilename,
                      1,
                      modelopacname,
                      atmosphere.MH,
                      atmosphere.AM,
                      abundances,
                      atmosphere.vt,
                      None,None,None,bsyn=False)
        # Run babsma
        sys.stdout.write('\r'+"Running Turbospectrum babsma_lu ...\r")
        sys.stdout.flush()
        if verbose:
            stdout= None
            stderr= None
        else:
            stdout= open('/dev/null', 'w')
            stderr= subprocess.STDOUT
        try:
            p= subprocess.Popen(['babsma_lu'],
                                cwd=twd,
                                stdin=subprocess.PIPE,
                                stdout=stdout,
                                stderr=stderr)
            with open(os.path.join(twd,'babsma.par'),'r') as parfile:
                for line in parfile:
                    p.stdin.write(line.encode('utf-8'))
            stdout, stderr= p.communicate()
        except subprocess.CalledProcessError:
            for linelistfilename in linelistfilenames:
                os.remove(linelistfilename,twd)
            if os.path.exists(os.path.join(twd,'DATA')):
                os.remove(os.path.join(twd,'DATA'))
            raise RuntimeError("Running babsma_lu failed ...")
        finally:
            #if os.path.exists(os.path.join(twd,'babsma.par')) \
            #   and outfname is None: #not 'saveTurboInput' in kwargs:
            #    os.remove(os.path.join(twd,'babsma.par'))
            sys.stdout.write('\r'+_ERASESTR+'\r')
            sys.stdout.flush()
        if isinstance(modelopac,str):
            shutil.copy(modelopacname,modelopac)
    else:
        shutil.copy(modelopac,twd)
        modelopacname= os.path.join(twd,os.path.basename(modelopac))

    # Now write the script file for bsyn_lu
    scriptfilename= os.path.join(twd,'bsyn.par')
    outfilename= os.path.join(twd,'bsyn.out')
    _write_script(scriptfilename,
                  wmin,wmax,dwl,
                  costheta,
                  modelfilename,
                  1,
                  modelopacname,
                  atmosphere.MH,
                  atmosphere.AM,
                  abundances, #indiv_abu,
                  None,
                  outfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=True)
    # Run bsyn
    sys.stdout.write('\r'+"Running Turbospectrum bsyn_lu ...\r")
    sys.stdout.flush()
    if verbose:
        stdout= None
        stderr= None
    else:
        stdout= open('/dev/null', 'w')
        stderr= subprocess.STDOUT
    try:
        p= subprocess.Popen(['bsyn_lu'],
                            cwd=twd,
                            stdin=subprocess.PIPE,
                            stdout=stdout,
                            stderr=stderr)
        with open(os.path.join(twd,'bsyn.par'),'r') as parfile:
            for line in parfile:
                p.stdin.write(line.encode('utf-8'))
        stdout, stderr= p.communicate()
    except subprocess.CalledProcessError:
        raise RuntimeError("Running bsyn_lu failed ...")
    finally:
        if outfname is not None:
            turbosavefilename= outfname
            if os.path.dirname(turbosavefilename) == '':
                turbosavefilename= os.path.join(os.getcwd(),turbosavefilename)
            try:
                subprocess.check_call(['tar','cvzf',turbosavefilename,
                                       os.path.basename(os.path.normpath(twd))])
            except subprocess.CalledProcessError:
                raise RuntimeError("Tar-zipping the Turbospectrum input and output failed; you will have to manually delete the temporary directory ...")
        #    # Need to remove babsma.par, bc not removed above
        #    if os.path.exists(os.path.join(twd,'babsma.par')):
        #        os.remove(os.path.join(twd,'babsma.par'))
        #if os.path.exists(os.path.join(twd,'bsyn.par')):
        #    os.remove(os.path.join(twd,'bsyn.par'))
        #if os.path.exists(modelopacname):
        #    os.remove(modelopacname)
        #if os.path.exists(modelopacname+'.mod'):
        #    os.remove(modelopacname+'.mod')
        #if os.path.exists(os.path.join(twd,'DATA')):
        #    os.remove(os.path.join(twd,'DATA'))
        #if os.path.exists(os.path.join(twd,'dummy-output.dat')):
        #    os.remove(os.path.join(twd,'dummy-output.dat'))
        #if os.path.exists(modelfilename):
        #    os.remove(modelfilename)
        #if rmLinelists:
        #    for linelistfilename in linelistfilenames[1:]:
        #        os.remove(linelistfilename)
        sys.stdout.write('\r'+_ERASESTR+'\r')
        sys.stdout.flush()
    
    # Now read the output
    turboOut= np.loadtxt(outfilename)
    # Clean up
    #os.remove(outfilename)
    #os.rmdir(twd)
    # Return wav, cont-norm, full spectrum
    return (turboOut[:,0],turboOut[:,1],turboOut[:,2])

def validate_abundances(abundances, MH):
    """ Input is format [(Z1, XFe1), (Z2, XFe2), ...] """
    assert isinstance(abundances, list), abundances
    if len(abundances) == 0: return {26: 0.0}
    Zs = []
    XFes = []
    for Z,XFe in abundances:
        Z = int(Z)
        assert (Z >= 3) & (Z <= 92), Z
        Zs.append(Z)
        
        XFe = np.round(float(XFe),3)
        assert (XFe >= -9) & (XFe <= 9), XFe
        XFes.append(XFe)
    # From Jo: Make sure to adjust C for the atmosphere's C value, by definitely including it (#61)
    if 6 not in Zs: abundances.append([6,0.])
    # Add N as well I guess
    if 7 not in Zs: abundances.append([7,0.])
    
    new_abundances = {}
    for Z, XFe in zip(Zs, XFes):
        new_abundances[Z] = XFe + atmosphere.MH + utils.get_solar(Z)
    
    return new_abundances

def _write_script(scriptfilename,
                  wmin,wmax,dw,
                  costheta,
                  modelfilename,
                  marcsfile,
                  modelopacname,
                  metals,
                  alphafe,
                  indiv_abu, # dictionary with atomic number, abundance
                  vmicro,
                  resultfilename,
                  isotopes,
                  linelistfilenames,
                  bsyn=False):
    """Write the script file for babsma and bsyn"""
    with open(scriptfilename,'w') as scriptfile:
        scriptfile.write("'LAMBDA_MIN:'  '%.3f'\n" % wmin)
        scriptfile.write("'LAMBDA_MAX:'  '%.3f'\n" % wmax)
        scriptfile.write("'LAMBDA_STEP:' '%.3f'\n" % dw)
        if bsyn:
            scriptfile.write("'INTENSITY/FLUX:' 'Flux'\n")
            scriptfile.write("'COS(THETA)    :' '%.3f'\n" % costheta)
            scriptfile.write("'ABFIND        :' '.false.'\n")
        scriptfile.write("'MODELINPUT:' '%s'\n" % modelfilename)
        if marcsfile is None:
            scriptfile.write("'MARCS-FILE:' '.false.'\n")
        scriptfile.write("'MODELOPAC:' '%s'\n" % modelopacname)
        if bsyn:
            scriptfile.write("'RESULTFILE :' '%s'\n" 
                             % resultfilename)
        scriptfile.write("'METALLICITY:'    '%.3f'\n" % metals)
        scriptfile.write("'ALPHA/Fe   :'    '%.3f'\n" % alphafe)
        scriptfile.write("'HELIUM     :'    '0.00'\n")
        scriptfile.write("'R-PROCESS  :'    '0.00'\n")
        scriptfile.write("'S-PROCESS  :'    '0.00'\n")
        # Individual abundances
        nabu= len(indiv_abu)
        if nabu > 0:
            scriptfile.write("'INDIVIDUAL ABUNDANCES:'   '%i'\n" % nabu)
            for abu in indiv_abu:
                scriptfile.write("%i %.3f\n" % (abu,indiv_abu[abu]))
        if bsyn:
            niso= len(isotopes)
            if niso > 0:
                scriptfile.write("'ISOTOPES : ' '%i'\n" % niso)
                for iso in isotopes:
                    scriptfile.write('%s %s\n' % (iso,isotopes[iso]))
            # Linelists
            nlines= len(linelistfilenames)
            scriptfile.write("'NFILES   :' '%i'\n" % nlines)
            for linelistfilename in linelistfilenames:
                scriptfile.write("%s\n" % linelistfilename)
            scriptfile.write("'SPHERICAL:'  'F'\n")
            scriptfile.write("30\n")
            scriptfile.write("300.00\n")
            scriptfile.write("15\n")
            scriptfile.write("1.30\n")
        else:
            scriptfile.write("'XIFIX:' 'T'\n")
            scriptfile.write("%.3f\n" % vmicro)
    return None
