# turbopy

Python interface to Turbospectrum.



Authors
-------
 - Alex Ji (University of Chicago)

With major contributions from Jean Somalwar and Ivanna Escala

Installation
------------
This only works with python 3 (needs format strings)

### Python Dependencies
* numpy
* astropy

### Turbospectrum
* Install Turbospectrum: https://github.com/bertrandplez/Turbospectrum2019
* You will need to have `babsma_lu` and `bsyn_lu` in your `$PATH`
* Define the environment variable `$TURBODATA=/path/to/Turbospectrum2019/DATA`

Usage
-----

This is an example of how this should work once it's all going:
```
import turbopy

wmin, wmax, dwl = 6700, 6720, 0.1
ll = turbopy.TSLineList("vald-6700-6720.list")
atmo = turbopy.MARCSModel.load("sun.mod")

wave, norm, flux = turbopy.run_synth(wmin, wmax, dwl,
                                     [12.0, 0.4], [6.0, 0.5], [8.0, 0.5],
                                     atmosphere=atmo, vt=1.0,
                                     linelist=ll, outfname="sun-6700-6720.tar.gz")
```

Right now if you have a linelist and model atmosphere that you like, this will work
(based on Jo Bovy's APOGEE code).

I am now making utilities to deal with linelists, model atmospheres, and such.

Citation and Acknowledgments
----------------------------
* Please cite Plez, B., 2012, Astrophysics Source Code Library, record ascl:1205.004
  see: http://adsabs.harvard.edu/abs/2012ascl.soft05004P
* A software citation for this package will come once it is stable and releasable.
* This code has benefited greatly from Jo Bovy's `apogee` package.
  You should follow that citation for now.
  https://github.com/jobovy/apogee#citing-this-code
* The structure for this package was based on https://github.com/uwescience/shablona
