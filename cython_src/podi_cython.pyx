
cimport cython

import numpy
cimport numpy

cdef extern void sigma_clip_mean__cy (double* pixels, int n_pixels, int n_images, double* output)


@cython.boundscheck(False)
@cython.wraparound(False)        
def sigma_clip_mean(
    numpy.ndarray[double, ndim=2, mode="c"] pixels not None,
    numpy.ndarray[double, ndim=1, mode="c"] returned not None,
):

    cdef int x, n
    # print pixels.shape

    x, n = pixels.shape[0], pixels.shape[1]
    # print "n_pixels=",x,"    n_images=",n

    sigma_clip_mean__cy(&pixels[0,0], x, n, &returned[0])
                                  
    return
