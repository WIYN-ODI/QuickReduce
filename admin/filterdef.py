#!/usr/bin/env python

import numpy
import os
import sys
dir, _ = os.path.split(os.path.abspath(sys.argv[0]))
sys.path.append("%s/../" % (dir))
from podi_commandline import *

from scipy.interpolate import interp1d
import matplotlib.pyplot

import podi_logging
import logging


if __name__ == "__main__":

    options = read_options_from_commandline(None)
    podi_logging.setup_logging(options)
    logger = logging.getLogger("FilterDefTool")

    # Loop over all commandline parameters and load entries starting with +
    # these are additional files that need to be multiplied with the filter-
    # curve to compute the total throughput

    convolve_file = []
    convolve_data = []
    input_filelist = []
    for x in get_clean_cmdline():
        if (x.startswith("+")):
            filename = x[1:]
            if (os.path.isfile(filename)):
                convolve_file.append(filename)
                convolve_data.append(numpy.loadtxt(filename))
        elif (os.path.isfile(x)):
            input_filelist.append(x)

    cutoff = float(cmdline_arg_set_or_default("-cutoff", 0))
    wlmin = float(cmdline_arg_set_or_default("-min", -1e9))
    wlmax = float(cmdline_arg_set_or_default("-max", 1e9))

    fig = matplotlib.pyplot.figure()
    ax = fig.add_subplot(111)

    fig2 = matplotlib.pyplot.figure()
    ax2 = fig2.add_subplot(111)
    
    pycode = []
    pycode.append("# filter-name: mean_pos, centerpos, fwhm, left_fwhm, right_fwhm, max_amplitude, mean_throughput, total_area, left-5%, right-5%")
    for filterfile in input_filelist:

        logger.info("Working on file %s" % (filterfile))

        try:
            data = numpy.loadtxt(filterfile)
        except:
            logger.error("Problem while loading %s" % (filterfile))
            continue

        #print data

        # Apply response cutoff to keep noise in the transmission curve from affecting the filter definition
        response = data[:,1]
        data[:,1][response < cutoff] = 0.
        numpy.savetxt(filterfile+".trunc1", data)

        # select only the range with response > 0 and within valid wavelength range
        good_data = (response > 0) & (data[:,0]>wlmin) & (data[:,0] < wlmax)
        index = numpy.arange(data.shape[0])
        left_edge = numpy.max([0, numpy.min(index[good_data])-1])
        right_edge = numpy.min([data.shape[0], numpy.max(index[good_data])+1])
        print left_edge, right_edge
        data = data[left_edge:right_edge+1]
        numpy.savetxt(filterfile+".trunc", data)


        # Now find, for each data point, the width of this datapoint
        width = numpy.zeros(shape=(data.shape[0]))

        #print data.shape[0]
        #print data[1:].shape[0]
        #print data[1:,0].shape[0]
        #print data[:-1,0].shape[0]

        #width = numpy.array([(data[i+1,0] - data[i-1])/2 if (i > 0 and i < data.shape[0]-1) else 0) for i in range(data.shape[0])])
        #width[1:] = (data[1:,0] - data[:-1,0]) / 2.
        for i in range(1, data.shape[0]-1):
            width[i] = (data[i+1,0] - data[i-1,0])/2

        numpy.savetxt('width', width)
        #print width

        # Check if max > 1 - if so, divide by 100
        if (numpy.max(data[:,1]) > 1):
            data[:,1] /= 100.
        data[:,1][data[:,1] < 0] = 0

        data_raw = data[:,1].copy()
        data = numpy.append(data, data_raw.reshape(data.shape[0],1), axis=1)

        _, fn = os.path.split(filterfile)
        ax.plot(data[:,0], data[:,1], label=fn)

        # Now put all additional files onto the same wavelength grid
        for i_component in range(len(convolve_file)):
            component = convolve_data[i_component]
            f = interp1d(component[:,0], component[:,1], bounds_error=False, fill_value=0)

            comp_matched = f(data[:,0])

            data[:,1] *= comp_matched
            _, bn = os.path.split(convolve_file[i_component])
            # ax.plot(data[:,0], data[:,1], label="after %s" % bn)


        max_amplitude = numpy.max(data[:,1])
        # print max_amplitude

        # now restrict the filter to only points with >0.1% throughput
        data_limited = data[data[:,1] > 0.001]#*max_amplitude]
        step_limited = width[data[:,1] > 0.001]#*max_amplitude]
        mean_throughput = numpy.sum(data_limited[:,1]*step_limited) / numpy.sum(step_limited)

        centerpos = numpy.sum(data_limited[:,0]*step_limited) / numpy.sum(step_limited)

        # Now multiply each wavelength point with its width and sum the data
        mean_pos = numpy.sum(data[:,0]*data[:,1]*width) / (numpy.sum(width*data[:,1]))
        # print mean_pos

        total_area = numpy.sum(data[:,1]*width)

        # interpolate the final filter curve to get a more accurate fwhm measurement
        fwhm_interp =  interp1d(data[:,0], data[:,1], bounds_error=False, fill_value=0)
        fine_grid = numpy.linspace(data[0,0], data[-1,0], 10*(data[-1,0]-data[0,0]))
        fwhm_fine =  fwhm_interp(fine_grid)

        left_fwhm = numpy.min(fine_grid[fwhm_fine >= 0.5*max_amplitude])
        right_fwhm = numpy.max(fine_grid[fwhm_fine >= 0.5*max_amplitude])

        left_5 = numpy.min(fine_grid[fwhm_fine >= 0.05*max_amplitude])
        right_5 = numpy.max(fine_grid[fwhm_fine >= 0.05*max_amplitude])

        fwhm = right_fwhm - left_fwhm
        ax2.plot(fine_grid, fwhm_fine, label=fn)

        # print left_fwhm, right_fwhm

        print "%30s: %7.2f %7.1f (%7.2f -- %7.2f) max=%0.4f mean=%.3f center=%.1f area=%7.2f   5%%=(%6.1f %6.1f)" % (
            fn, mean_pos, fwhm, left_fwhm, right_fwhm, max_amplitude, mean_throughput, centerpos, total_area, left_5, right_5)
        pycode.append(
            ' %30s: (%30s, %7.2f, %7.2f, %7.1f, %7.1f, %7.1f, %0.4f, %0.4f, %7.2f, %7.2f, %7.2f),' % (
                ('"%s"' % fn), ('"%s"' % fn), mean_pos, centerpos, fwhm, left_fwhm, right_fwhm, max_amplitude, mean_throughput, total_area, left_5, right_5)
        )
        
        # Save the corrected filter curve
        numpy.savetxt("%s.final" % (fn), data)

    print "\n\n\n"
    print "\n".join(pycode)+"\n"

    ax.legend(fontsize=8, loc='best')
    ax2.legend(fontsize=8, loc='best')
    #fig.show()


#    ax.set_xlim((mean_pos - fwhm, mean_pos+fwhm))
    fig.savefig("filter.png")
    fig2.savefig("filter2.png")


    podi_logging.shutdown_logging(options)
