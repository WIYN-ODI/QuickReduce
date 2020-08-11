#!/usr/bin/env python3

from __future__ import print_function

import os
import sys
import astropy.io.fits as pyfits
import numpy
import scipy
import scipy.optimize

def model(p, times):
    return p[0] * numpy.exp(-times/p[1]) + p[2] + p[3]*numpy.exp(-times/p[4])

def error(p, times, data, scatter):
    _model = model(p, times)
    return (data - _model)/scatter


def measure_frame(fn, boxes):

    hdulist = pyfits.open(fn)
    amplitude = numpy.zeros((boxes.shape[0],2))

    timestamp = hdulist[0].header['MJD-OBS'] if ('MJD-OBS' in hdulist[0].header) else 0.

    for i,box in enumerate(boxes):
        cutout = hdulist[1].data[box[0]:box[1], box[2]:box[3]]
        #print(cutout.shape)

        amplitude[i,:] = cutout.size, numpy.nanmedian(cutout)

    return amplitude, timestamp

if __name__ == "__main__":

    boxes = []

    ds9_fn = sys.argv[1]
    with open(ds9_fn, "r") as f:
        lines = f.readlines()[3:]
        for line in lines:
            c = [int(numpy.round(float(c),0)) for c in line.split("box(")[1].split(")")[0].split(",")[:4]]
            # print(line.split("box(")[1].split(")")[0].split(","))

            boxes.append([c[0]-c[2]//2, c[0]+c[2]//2, c[1]-c[3]//2, c[1]+c[3]//2])
    boxes = numpy.array(boxes)

    print(boxes.shape)
    # numpy.savetxt(sys.stdout, boxes.astype(numpy.int))

    # get model amplitudes
    absolute, _ = measure_frame(sys.argv[2], boxes)
    print(absolute)

    filelist = sys.argv[3:]
    timestamps = numpy.zeros((len(filelist)))
    amplitudes = numpy.zeros((len(filelist), boxes.shape[0],2))
    for i_fn, fn in enumerate(filelist):
        print(fn)
        # hdulist = pyfits.open(fn)
        # timestamps[i_fn] = hdulist[0].header['MJD-OBS']
        # for i,box in enumerate(boxes):
        #     cutout = hdulist[1].data[box[0]:box[1], box[2]:box[3]]
        #     print(cutout.shape)

        #     amplitudes[i_fn, i,:] = cutout.size, numpy.nanmedian(cutout)
        amp, ts = measure_frame(fn, boxes)
        amplitudes[i_fn, :, :] = amp
        timestamps[i_fn] = ts

    # Sort all data by timestamps
    time_sort = numpy.argsort(timestamps)
    
    start_time = numpy.min(timestamps)
    timestamps = (timestamps - start_time) * 86400.

    timestamps = timestamps[time_sort]
    amplitudes = amplitudes[time_sort, :, :]
    print(amplitudes.shape)

    with open("testdump", "wb") as dmp:
        # 
        for box in range(amplitudes.shape[1]):
            print("BOX", box)
            x = numpy.array([timestamps[:], amplitudes[:,box,1]]).T
            print(x.shape)
            numpy.savetxt(dmp, x)
            dmp.write(b"\n"*5)
            print(x)

    # Now normalize all measurements with the amplitude of the model frame
    normalized = amplitudes / absolute
    with open("testdump.norm", "wb") as dmp:
        # 
        for box in range(amplitudes.shape[1]):
            print("BOX", box)
            x = numpy.array([timestamps[:], normalized[:,box,1]]).T
            numpy.savetxt(dmp, x, fmt='%f')
            dmp.write(b"\n"*5) #, file=dmp)
            print(x)

    # Now average over all different boxes
    averages = numpy.mean(normalized[:,:,1], axis=1)
    scatter = numpy.std(normalized[:,:,1], axis=1)
    scatter[scatter<=0] = numpy.mean(scatter[scatter>0])
    print(averages)
    print(scatter)
    numpy.savetxt("testdump.avg", numpy.array([
        timestamps[:], averages[:]]).T)
    
    bestfits = []
    for _try in range(10):
        # Now we have something we can fit to get the time-dependent amplitude
        # parameters are amplitude, tau, intensity_0
        p_init = [1.0, 980., 0, numpy.random.random(), 210] #numpy.random.random()*0.5, numpy.random.random()*0.5]

        # results = numpy.optimize.
        fit = scipy.optimize.leastsq(error, p_init, (timestamps, averages, scatter), full_output=1)
        print(fit)

        best_fit = fit[0]
        bestfits.append(best_fit)

    bestfits = numpy.array(bestfits)
    print("BEST FIT:\n",bestfits)

    # uncert = numpy.sqrt(numpy.diag(fit[1]))
    uncert = fit[1]
    print("UNCERTAINTIES:\n", uncert)

    _model = numpy.array([timestamps, model(best_fit, timestamps)]).T
    numpy.savetxt("testdump.bestfit", _model)

        #amplitudes[i_fn] = data
        #numpy.savetxt(fn+".data", data)
