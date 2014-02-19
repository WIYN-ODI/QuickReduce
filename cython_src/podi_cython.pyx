
cimport cython

import numpy
cimport numpy

cdef extern void sigma_clip_mean__cy (double* pixels, int n_pixels, int n_images, double* output,
                                      double nsigma, int max_repeat)
cdef extern void sigma_clip_median__cy (double* pixels, int n_pixels, int n_images, double* output,
                                      double nsigma, int max_repeat)
cdef extern void lacosmics__cy(double* data,
                               double* out_cleaned, int* out_mask, int* out_saturated,
                               int sx, int sy,
                               double gain, double readnoise,
                               double sigclip, double sigfrac, double objlim,
                               double saturation_limit, int verbose,
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
        numpy.ndarray[double, ndim=2, mode="c"] data_in not None,
        numpy.ndarray[double, ndim=2, mode="c"] cleaned = None,
        numpy.ndarray[int, ndim=2, mode="c"] mask = None,
        numpy.ndarray[int, ndim=2, mode="c"] saturated = None,
        double gain = 1.,
        double readnoise = 0.,
        double sigclip = 4.5,
        double sigfrac = 0.5,
        double objlim = 1.0,
        int niter = 4,
        double saturation_limit = 60000,
        int verbose = True,
):

    cdef int x, n
    # print pixels.shape

    x, n = data_in.shape[0], data_in.shape[1]
    # print "n_pixels=",x,"    n_images=",n

    if (cleaned == None):
        cleaned = numpy.ndarray(shape=(data_in.shape[0], data_in.shape[1]), dtype=numpy.float64)
    if (mask == None):
        mask = numpy.ndarray(shape=(data_in.shape[0], data_in.shape[1]), dtype=numpy.int32)
    if (saturated == None):
        saturated = numpy.ndarray(shape=(data_in.shape[0], data_in.shape[1]), dtype=numpy.int32)

    lacosmics__cy(&data_in[0,0], 
                  &cleaned[0,0], &mask[0,0], &saturated[0,0],
                  data_in.shape[0], data_in.shape[1],
                  gain, readnoise,
                  sigclip, sigfrac, objlim,
                  saturation_limit,
                  verbose,
                  niter)
                                  
    return cleaned, mask, saturated

