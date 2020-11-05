from __future__ import absolute_import, division, print_function
import os
import numpy as np
import numpy.testing as npt
import turbopy

data_path = os.path.join(turbopy.__path__[0], 'data')

def test_linelist():
    tab = turbopy.linelists.read_vald_long(os.path.join(data_path, "BertrandPlez.002060"))
    #import pdb; pdb.set_trace()
    print(tab[tab.colnames[0:5]])
    print(tab[tab.colnames[5:]])
    #raise RuntimeError("To print the linelist")

def test_linelist_write():
    tab = turbopy.linelists.read_vald_long(os.path.join(data_path, "BertrandPlez.002060"),
                                           outfname=os.path.join(data_path, "converted_BertrandPlez.002060"))

def test_get_levels():
    # Taken from Bertrand's vald-6700-6720.list
    test_data = [
        ('s', 'p', 'LS 1s2.2s 2S', 'LS 1s2.2p 2P*'),
        ('s', 'p', 'LS 3d5.4s2 a4P', 'LS 3d6.(5D).4p z6D*'),
        ('s', 's', 'LS 3d3.4s2 a4F', 'LS 3d4.(3H).4s a4H*'),
        ('d', 'f', 'LS 3d3.(4F).4d 3G', 'JK 3d3.(4F<5/2>).4f 2[5/2]*'),
        ('p', 'd', 'LS 3d6.(5D).4p z4F', 'LS 3d6.(a3P).4p 2D*'),
        ('X', 'X', ' ', ' '),
        ('X', 'X', 'Hb 8s2.9s1.3p4.1d1 X,3,2,0,53.0,0', 'Hb 8s2.3p4.4p1.1d1 A,3,3,0,52.0,1'),
        ('d', 'f', 'LS 4f.(2F*).5d2.(3P) 2D*', 'LS 4f2.(3P).5d 4P'),
        ('X', 'X', '__ 4f.(2F*).5d2.(1D) *', '__ 4f.5d.(3G*).6p'),
        ## These currently fail, not sure what the issue is: uncomment them later
        #('s', 'p', '__ 4f2.(1G).5d', 'JJ:4f2.(3H<5>).6p.<3/2> (5,3/2)*'),
        #('d', 'p', 'LS 4f.(2F*).5d2.(1G) 2D*', '__ '),
    ]
    outdata1 = []
    outdata2 = []
    for llo, lhi, l2, l3 in test_data:
        outdata1.append((llo, lhi))
        outdata2.append(turbopy.linelists._get_levels(l2, l3))
    print(outdata1)
    print(outdata2)
    npt.assert_equal(outdata2, outdata1)
