#!/usr/local/bin/python

import pyfits
import math
import scipy

import scipy.ndimage
import scipy.stats
import numpy
import bottleneck

def rebin_image(data, binfac):
    
    if (binfac < 1):
        stdout_write("Rebinning at the moment only supports binning to larger pixels with binfac>1\n")
        return None
    elif (binfac == 1):
        return data
    
    out_size_x, out_size_y = int(math.ceil(data.shape[0]*1.0/binfac)), int(math.ceil(data.shape[1]*1.0/binfac))
    
    if (out_size_x*binfac != data.shape[0] or out_size_y*binfac != data.shape[1]):
        # The input array size is not a multiple of the new binning
        # Create a slightly larger array to hold the data to be rebinned
        container = numpy.zeros(shape=(out_size_x*binfac, out_size_y*binfac))
        
        # And insert the original data
        container[0:data.shape[0], 0:data.shape[1]] = data[:,:]
    else:
        container = data 
    
    rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).mean(axis=-1).mean(axis=1)
    #rebinned = numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)).nm(axis=-1).nm(axis=1)
    
    #rs = numpy.array(numpy.reshape(container, (out_size_x, binfac, out_size_y, binfac)), dtype=numpy.float32)
    #rb1 = bottleneck.nanmean(rs, axis=-1)
    #rb2 = bottleneck.nanmean(rb1, axis=1)
    #rebinned = rb2
    
    return rebinned





def find_center(hdu, lx, ly, prebin=8, r_minmax=[100,200], x_minmax=[400,600], y_minmax=[400,600], dx=2, dy=2, dr=2):


    ot33b = rebin_image(hdu.data, prebin)

    ot33_orig = ot33b.copy()

    ot33b[numpy.isnan(ot33b)] = 0

    x33 = scipy.ndimage.sobel(ot33b, axis=0, mode='constant')
    y33 = scipy.ndimage.sobel(ot33b, axis=1, mode='constant')
    abs33 = numpy.hypot(x33, y33)

    mask = numpy.array(numpy.isnan(ot33_orig), dtype=numpy.float32)
    mask[3,3] = 1.
    #print mask[:8,:8]
    numpy.savetxt("mask", mask)
    kernel = numpy.ones(shape=(5,5))
    kernel_norm = kernel

    print "number valid pixels before growing",numpy.sum((mask == False))

    mask_grown_float = scipy.ndimage.filters.convolve(mask, kernel)
    mask_grown = mask_grown_float > 0
    numpy.savetxt("mask.g", mask_grown)
    #print mask_grown_float[:8,:8]

    print "number valid pixels after growing",numpy.sum((mask_grown == False))

    abs33[numpy.isnan(ot33_orig)] = numpy.NaN

    #pyfits.HDUList([pyfits.PrimaryHDU(data=ot33b)]).writeto("/scratch/ot33b.fits", clobber=True)
    #pyfits.HDUList([pyfits.PrimaryHDU(data=abs33)]).writeto("/scratch/ot33_edge.fits", clobber=True)

    abs33_binary = abs33.copy()
    abs33_binary[abs33 < 0.1] = 0
    abs33_binary[abs33 >= 0.1] = 1

    abs33[mask_grown] = numpy.NaN

    #pyfits.HDUList([pyfits.PrimaryHDU(data=abs33_binary)]).writeto("/scratch/ot33_binary.fits", clobber=True)

    # Figure out what contrast we need
    valid_pixels = numpy.isfinite(abs33)
    #print "number valid pixels",numpy.sum(valid_pixels)
    top10percent = scipy.stats.scoreatpercentile(abs33[valid_pixels].ravel(), 90)
    strong_values = abs33 > top10percent

    print "Only using pixels >",top10percent

    all_y, all_x = numpy.indices(abs33.shape)

    pixel_x = all_x[strong_values]
    pixel_y = all_y[strong_values]

    pixel_value = abs33[strong_values]

    #print pixel_x

    print numpy.sum(strong_values),"pixels with enough signal left"

    dx = 2
    dy = 2
    dr = 2
    #center_x = numpy.linspace(400,600,101)
    #center_y = numpy.linspace(400,600,101)

    lr = 100
    radius = numpy.linspace(100,200,51)

    #print center_x


    bincount = numpy.zeros(shape=(100,100,50))

    #
    # Now do the hough transformation
    #
    for i_cx in range(100):
        for i_cy in range(100):

            cx = lx + i_cx * dx
            cy = ly + i_cy * dy

            pixel_radius = numpy.sqrt( (cx-pixel_x)**2 + (cy-pixel_y)**2 )
            #print pixel_radius

            #count,edges = numpy.histogram(pixel_radius, bins=radius)
            count,edges = numpy.histogram(pixel_radius, bins=radius, weights=pixel_value)

            #print count.shape, edges.shape

            bincount[i_cx, i_cy, :] = count[:]


    fixed_radius = True

    if (fixed_radius):
        ir = 31
        center_only = bincount[:,:,ir]
        
        index = numpy.argmax(center_only)
        ix, iy = numpy.unravel_index(index, center_only.shape)

    else:
        print bincount.shape
        #numpy.savetxt("edges.txt", bincount[:,50,:])

        index = numpy.argmax(bincount)

        ix, iy, ir = numpy.unravel_index(index, bincount.shape)

        print index, numpy.unravel_index(index, bincount.shape)

    return ix*dx+lx, iy*dy+ly, ir*dr+lr, bincount, abs33



if __name__ == "__main__":

    #hdu = pyfits.open("/nas/wiyn/pupilghost_template/diff_l2-l1/diff___odi_g____-045.fits")
    hdu = pyfits.open("/nas/wiyn/pupilghost_template/frames/20130320T144234.0___dflat_ghost_tempalte_g_filter_layer_2.fits")

    prebin=8
#    lx = [0,400,400,-100,-100]
#    ly = [0, 400,-100,-100,400]

    lx = numpy.array([0, 3200, 3200, -800, -800]) / prebin
    ly = numpy.array([0, 3200, -800, -800, 3200]) / prebin

    dx, dy, dr = 16,16,16

    for i in range(1,5):
        print hdu[i].header["EXTNAME"]

        r_minmax=[800,1600] 
        x_minmax=[lx[i], lx[i]+1600]
        y_minmax=[ly[i], ly[i]+1600]

        x, y, r, bincount, edge_frame = find_center(hdu[i], lx[i], ly[i], 
                                                    x_minmax=x_minmax, y_minmax=y_minmax, r_minmax=r_minmax,
                                                    dx=dx, dy=dy, dr=dr)
        numpy.savetxt("bincount"+hdu[i].header["EXTNAME"], numpy.sum(numpy.sum(bincount, axis=0), axis=0))
        print x,y,r, " ---> ", x*prebin+1, y*prebin+1, r*8

        hdu[i].data = edge_frame
        print

    hdu.writeto("/scratch/edges.fits", clobber=True)

