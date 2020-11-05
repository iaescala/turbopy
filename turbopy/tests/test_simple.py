from __future__ import absolute_import, division, print_function
import os
import numpy as np
import numpy.testing as npt
import turbopy

data_path = os.path.join(turbopy.__path__[0], 'data')

def test_nothing_useful():
    """
    Test testing
    """
    npt.assert_equal(0.0, 0.0)
    npt.assert_almost_equal(np.pi, 4*np.arctan(1.0))
    
def test_simple():
    """
    Run a default synthesis
    """
    
    wmin, wmax, dwl = 6700, 6720, 0.1
    ll = turbopy.TSLineList(os.path.join(data_path, "vald-6700-6720.list"))
    atmo = turbopy.MARCSModel.load(os.path.join(data_path, "sun.mod"))
    wave, norm, flux = turbopy.run_synth(wmin, wmax, dwl,
                                         atmosphere=atmo,
                                         linelist=ll)
    
