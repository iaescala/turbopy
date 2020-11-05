from __future__ import absolute_import, division, print_function
import os
import numpy as np
import numpy.testing as npt
import turbopy
from turbopy import utils

things_to_test = [
    ("Ba 2", "(137)Ba+", ["Ba"], 2, [137], '56.137'),
    ("TiO 1", "(48)TiO", ["Ti","O"], 1, [48, 0], '822.000048'),
    ("TiO 1", "(49)TiO", ["Ti","O"], 1, [49, 0], '822.000049'),
    ("CH 1", "(12)CH", ["C","H"], 1, [12, 0], '106.000012'),
    ("CH 1", "CH", ["C","H"], 1, [0, 0], '106.000000'),
    ("OH 1", "(16)OH", ["O","H"], 1, [16, 0], '108.000016'),
    ("CN 1", "(13)C(14)N", ["C","N"], 1, [13, 14], '607.013014'),
    ("CN 1", "(12)C(14)N", ["C","N"], 1, [12, 14], '607.012014'),
    ("C2 1", "(12)C(12)C", ["C","C"], 1, [12, 12], '606.012012'),
    ("MgH 1", "(24)MgH", ["Mg","H"], 1, [24, 0], '112.000024'),
]

def test_identify_fullspecstr():
    """
    Tests some common cases for identifying from full VALD linelist
    """
    for specstr, fullspecstr, elems, ion, isos, tspecies in things_to_test:
        out = utils.identify_fullspecstr(specstr, fullspecstr)
        npt.assert_equal(elems, out[0])
        npt.assert_equal(ion, out[1])
        npt.assert_equal(isos, out[2])

def test_make_tspecies():
    """
    Tests some common cases for identifying from full VALD linelist
    """
    for specstr, fullspecstr, elems, ion, isos, tspecies in things_to_test:
        out = utils.identify_fullspecstr(specstr, fullspecstr)
        Zs = [utils.elem_to_Z(el) for el in out[0]]
        new_tspecies = utils.make_tspecies(Zs, out[2]).strip()
        npt.assert_equal(tspecies, new_tspecies)
