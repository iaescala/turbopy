from __future__ import absolute_import, division, print_function

import os
import re
import numpy as np
from astropy.table import Table

_all_elems = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
              "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
              "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
              "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
              "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
              "Fr","Ra","Ac","Th","Pa","U"]
## These are from Turbospectrum vald3line-BPz-freeformat.f
#* Solar abundances ref: Grevesse N., Asplund A., Sauval A.J. 2007,
#* Space Science Review 130, 105, DOI 10.1007/s11214-007-9173-7
_sunabund = [
     12.00,10.9275, 1.05,  1.38,  2.70,  8.39,  7.78,  8.66,  4.56,   #  1 -  9
      7.84,  6.17,  7.53,  6.37,  7.51,  5.36,  7.14,  5.50,  6.18,   # 10 - 18
      5.08,  6.31,  3.17,  4.90,  4.00,  5.64,  5.39,  7.45,  4.92,   # 19 - 27
      6.23,  4.21,  4.60,  2.88,  3.58,  2.29,  3.33,  2.56,  3.25,   # 28 - 36
      2.60,  2.92,  2.21,  2.58,  1.42,  1.92, -99.0,  1.84,  1.12,   # 37 - 45
      1.66,  0.94,  1.77,  1.60,  2.00,  1.00,  2.19,  1.51,  2.24,   # 46 - 54
      1.07,  2.17,  1.13,  1.70,  0.58,  1.45, -99.0,  1.00,  0.52,   # 55 - 63
      1.11,  0.28,  1.14,  0.51,  0.93,  0.00,  1.08,  0.06,  0.88,   # 64 - 72
     -0.17,  1.11,  0.23,  1.25,  1.38,  1.64,  1.01,  1.13,  0.90,   # 73 - 81
      2.00,  0.65, -99.0, -99.0, -99.0, -99.0, -99.0, -99.0,  0.06,   # 82 - 90
     -99.0, -0.52]                                                    # 91 - 92
#* Ionization potentials (eV), from file "atomdata"
_ionp1 = [
     13.60, 24.59,  5.39,  9.32,  8.30, 11.26, 14.53, 13.62, 17.40,
     21.56,  5.14,  7.64,  5.99,  8.15, 10.48, 10.36, 12.97, 15.76,
      4.34,  6.11,  6.54,  6.82,  6.74,  6.77,  7.44,  7.87,  7.86,
      7.64,  7.73,  9.39,  6.00,  7.88,  9.81,  9.75, 11.81, 14.00,
      4.18,  5.69,  6.38,  6.84,  6.88,  7.10,  7.28,  7.36,  7.46,
      8.33,  7.57,  8.99,  5.79,  7.34,  8.64,  9.01, 10.45, 12.13,
      3.89,  5.21,  5.58,  5.65,  5.42,  5.49,  5.55,  5.63,  5.68,
      6.16,  5.85,  5.93,  6.02,  6.10,  6.18,  6.25,  5.43,  7.00,
      7.88,  7.98,  7.87,  8.50,  9.10,  9.00,  9.22, 10.44,  6.11,
      7.42,  7.29,  8.42,  9.30, 10.75,  4.00,  5.28,  6.90,  6.08,
      9.99,  6.00
]
_ionp2 = [
       .00, 54.42, 75.64, 18.21, 25.15, 24.38, 29.60, 35.12, 35.00,
     40.96, 47.29, 15.03, 18.83, 16.35, 19.72, 23.33, 23.81, 27.62,
     31.63, 11.87, 12.80, 13.58, 14.65, 16.50, 15.64, 16.18, 17.06,
     18.17, 20.29, 17.96, 20.51, 15.93, 18.63, 21.19, 21.60, 24.36,
     27.50, 11.03, 12.24, 13.13, 14.32, 16.15, 15.26, 16.76, 18.07,
     19.42, 21.49, 16.90, 18.87, 14.63, 16.53, 18.60, 19.13, 21.21,
     25.10, 10.00, 11.06, 10.85, 10.55, 10.73, 10.90, 11.07, 11.25,
     12.10, 11.52, 11.67, 11.80, 11.93, 12.05, 12.17, 13.90, 14.90,
     16.20, 17.70, 16.60, 17.50, 20.00, 18.56, 20.50, 18.76, 20.43,
     15.03, 16.68, 19.00, 20.00, 21.00, 22.00, 10.14, 12.10, 11.50,
     99.99, 12.00
]

def Z_to_elem(Z):
    return _all_elems[Z-1]
def elem_to_Z(elem):
    try:
        ix = _all_elems.index(elem)
    except ValueError as e:
        raise ValueError(f"'{elem}' is not a valid element")
    return ix + 1
def get_solar(Z_or_elem):
    if isinstance(Z_or_elem, str): Z_or_elem = elem_to_Z(Z_or_elem)
    return _sunabund[Z_or_elem-1]
def get_ionp1(Z_or_elem):
    if isinstance(Z_or_elem, str): Z_or_elem = elem_to_Z(Z_or_elem)
    return _ionp1[Z_or_elem-1]
def get_ionp2(Z_or_elem):
    if isinstance(Z_or_elem, str): Z_or_elem = elem_to_Z(Z_or_elem)
    return _ionp2[Z_or_elem-1]

def identify_specstr(specstr):
    """ Figure out what is the species number, isotope, etc for """
    if specstr.startswith("'"): specstr = specstr[1:-1]
    species, ion = specstr.split()
    upperspecies = species.upper()
    ix_capital, = np.where(np.array([x==y for x,y in zip(species, upperspecies)]))
    Nelem = len(ix_capital)
    ix = list(ix_capital) + [len(species)]
    elems = [species[ix[i]:ix[i+1]] for i in range(Nelem)]
    if species == "C2": elems = ["C", "C"]
    Zs = [elem_to_Z(el) for el in elems]
    return elems, int(ion)
def identify_fullspecstr(specstr, fullspecstr):
    """
    Uses specstr to find the info, then fullspecstr to get isotopes
    fullspecstr examples: (48)TiO, Fe+, Mn, (12)C(14)N, (137)Ba+, ...
    """
    elems, ion = identify_specstr(specstr)
    assert len(elems) <= 2, "Doesn't work yet for 3-atom molecules"
    isos = [0 for el in elems]
    # Get any isotopes in the string
    #all_isostarts, = np.where(np.array(list(fullspecstr)) == "(")
    isotopes = re.findall(r"\((\d+)\)",fullspecstr)
    if len(isotopes) == len(elems):
        isos = [int(iso) for iso in isotopes]
    elif len(isotopes) == 1: # this works only for 2 elements
        if fullspecstr[0]=="(": isos[0] = int(isotopes[0])
        else: isos[1] = int(isotopes[0])
    return elems, ion, isos

def make_tspecies(Zs, isos=None):
    if len(Zs) == 1:
        return f"{Zs[0]:>2}.{isos[0]:03}"
    if len(Zs) == 2:
        if Zs[0] > Zs[1]:
            # swap in place, yay python
            Zs[0], Zs[1] = Zs[1], Zs[0]
            isos[0], isos[1] = isos[1], isos[0]
        if Zs[0] == Zs[1] and isos[0] > isos[1]:
            isos[0], isos[1] = isos[1], isos[0]
        return f"{Zs[0]:>2}{Zs[1]:02}.{isos[0]:03}{isos[1]:03}"
    raise ValueError("make_tspecies only works up to to-atom molecules")
