
cimport cython

import numpy
cimport numpy

cdef extern void sigma_clip_mean__cy (double* pixels, int n_pixels, int n_images, double* output,
                                      double nsigma, int max_repeat)
cdef extern void sigma_clip_median__cy (double* pixels, int n_pixels, int n_images, double* output,
                                      double nsigma, int max_repeat)
cdef extern void lacosmics__cy(double* pixels, double* output, double* output2,
                   int sx, int sy,
                   double gain, double readnoise,
                   double sigclip, double sigthres, double objlim,
                   int niter)


@cython.boundscheck(False)
@cython.wraparound(False)        
def sigma_clip_mean(
        numpy.ndarray[double, ndim=2, mode="c"] pixels not None,
        numpy.ndarray[double, ndim=1, mode="c"] returned not None,
        double nsigma = 3,
        int max_repeat = 3,
):

    cdef int x, n
    # print pixels.shape

    x, n = pixels.shape[0], pixels.shape[1]
    # print "n_pixels=",x,"    n_images=",n

    sigma_clip_mean__cy(&pixels[0,0], x, n, &returned[0], nsigma, max_repeat)
                                  
    return




@cython.boundscheck(False)
@cython.wraparound(False)        
def sigma_clip_median(
        numpy.ndarray[double, ndim=2, mode="c"] pixels not None,
        numpy.ndarray[double, ndim=1, mode="c"] returned not None,
        double nsigma = 3,
        int max_repeat = 3,
):

    cdef int x, n
    # print pixels.shape

    x, n = pixels.shape[0], pixels.shape[1]
    # print "n_pixels=",x,"    n_images=",n

    sigma_clip_median__cy(&pixels[0,0], x, n, &returned[0], nsigma, max_repeat)
                                  
    return



@cython.boundscheck(False)
@cython.wraparound(False)        
def lacosmics(
        numpy.ndarray[double, ndim=2, mode="c"] pixels not None,
        numpy.ndarray[double, ndim=2, mode="c"] returned not None,
        numpy.ndarray[double, ndim=2, mode="c"] returned2 not None,
        double gain = 0,
        double readnoise = 0,
        double sigclip = 4.5,
        double sigthres = 0.5,
        double objlim = 1.0,
        int niter = 4,
):

    cdef int x, n
    # print pixels.shape

    x, n = pixels.shape[0], pixels.shape[1]
    # print "n_pixels=",x,"    n_images=",n

    lacosmics__cy(&pixels[0,0], &returned[0,0], &returned2[0,0],
                  pixels.shape[0], pixels.shape[1],
                  gain, readnoise,
                  sigclip, sigthres, objlim,
                  niter)
                                  
    return
