from __future__ import absolute_import, division, print_function
import os
import numpy as np
import numpy.testing as npt
import turbopy
import tempfile

data_path = os.path.join(turbopy.__path__[0], 'data')

def test_synth_simple():
    """
    Run a default synthesis
    """
    
    twd = tempfile.mkdtemp(dir=os.getcwd()+"/tmp")
    print(twd)
    
    wmin, wmax, dwl = 6700, 6720, 0.1
    ll = turbopy.TSLineList(os.path.join(data_path, "vald-6700-6720.list"))
    atmo = turbopy.MARCSModel.load(os.path.join(data_path, "sun.mod"))
    atmo.Teff = 5777
    atmo.logg = 4.44
    atmo.MH = 0.0
    atmo.AM = 0.0
    wave, norm, flux = turbopy.run_synth(wmin, wmax, dwl,
                                         atmosphere=atmo, vt=1.0,
                                         linelist=ll, twd=twd)
    
def test_synth_linelist():
    """
    Run a default synthesis
    """
    
    twd = tempfile.mkdtemp(dir=os.getcwd()+"/tmp")
    print(twd)
    
    wmin, wmax, dwl = 6700, 6720, 0.1
    ll1 = turbopy.TSLineList(os.path.join(data_path, "vald-6700-6720.list"))
    ll2 = turbopy.TSLineList(os.path.join(data_path, "converted_BertrandPlez.002060"))
    atmo = turbopy.MARCSModel.load(os.path.join(data_path, "sun.mod"))
    atmo.Teff = 5777
    atmo.logg = 4.44
    atmo.MH = 0.0
    atmo.AM = 0.0
    wave1, norm1, flux1 = turbopy.run_synth(wmin, wmax, dwl,
                                            atmosphere=atmo, vt=1.0,
                                            linelist=ll1, twd=twd)
    wave2, norm2, flux2 = turbopy.run_synth(wmin, wmax, dwl,
                                            atmosphere=atmo, vt=1.0,
                                            linelist=ll2, twd=twd)
    
    #import matplotlib.pyplot as plt
    #fig = plt.figure()
    #plt.plot(wave1, norm1, 'k-', lw=3)
    #plt.plot(wave2, norm2, 'r-', lw=1)
    #fig.savefig("test.pdf")
    
    #np.save("test.npy", np.vstack([norm1,norm2]))

    npt.assert_almost_equal(wave1, wave2)
    npt.assert_almost_equal(norm1, norm2)
    npt.assert_almost_equal(flux1, flux2)
