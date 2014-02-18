from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
from Cython.Distutils import build_ext

import numpy

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [
        Extension("podi_cython", 
                  sources=['cython_src/podi_cython.pyx', 
                           "cython_src/sigma_clip_mean.c",
                           "cython_src/sigma_clip_median.c",
                           "cython_src/lacosmics.c",
                       ],
                  include_dirs=["cython_src", numpy.get_include()],
                  libraries=['gslcblas', "gsl", "m"]
                  )
    ]
)
