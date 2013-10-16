#!/usr/local/bin/python
#
# Copyright 2012-2013 Ralf Kotulla
#                     kotulla@uwm.edu
#
# This file is part of the ODI QuickReduce pipeline package.
#
# If you find this program or parts thereof please make sure to
# cite it appropriately (please contact the author for the most
# up-to-date reference to use). Also if you find any problems 
# or have suggestiosn on how to improve the code or its 
# functionality please let me know. Comments and questions are 
# always welcome. 
#
# The code is made publicly available. Feel free to share the link
# with whoever might be interested. However, I do ask you to not 
# publish additional copies on your own website or other sources. 
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
#

"""
How-to use:

./podi_collectcells.py 


"""

import sys
import os
import pyfits
import numpy
import scipy
import scipy.optimize

from podi_plotting import *

gain_correct_frames = False
from podi_definitions import *



import podi_plotting




colors = ['', '#900000', '#00a000', '#0000a0']


bordersize = 75
max_polyorder=0


def create_nonlinearity_data(inputfiles):

    # Open all files, and loop over all cells in all extensions:

    all_data = []

    for filename in inputfiles:

        try:
            hdulist = pyfits.open(filename)
            exptime = hdulist[0].header['EXPTIME']
            expmeas = hdulist[0].header['EXPMEAS']
            print "Working on",filename, "exptime=",exptime
        except:
            print "#####"
            print "#####"
            print "Error opening file",filename
            print "#####"
            print "#####"
            continue

        for ext in range(1, len(hdulist)):
            if (not type(hdulist[ext]) == pyfits.hdu.image.ImageHDU):
                continue

            extname = hdulist[ext].header['EXTNAME']
            data = hdulist[ext].data
            ota = int(extname[3:5])
            otax = int(extname[3])
            otay = int(extname[4])

            for cellx in range(8):
                for celly in range(8):
                
                    cell_area = cell2ota__get_target_region(cellx, celly)
                    x1, x2, y1, y2 = cell_area

                    cell_data = data[y1:y2, x1:x2]

                    cell_center = cell_data[bordersize:-bordersize,bordersize:-bordersize]

                    median_int = numpy.median(cell_center)
                    mean_int = numpy.mean(cell_center)
                    std_int = numpy.std(cell_center)
                    #print exptime, extname, ota, cellx, celly, median_int

                    thiscell = [ota, otax, otay, cellx, celly, exptime, expmeas, median_int, mean_int, std_int]
                    # if (ext == 1 and cellx == 0 and celly == 0):
                    #     print thiscell,
                    # else:
                    #     stdout_write(".")

                    if (not numpy.isnan(median_int)):
                        stdout_write("\rOTA %02d, cell %d,%d: median=%6.0f, avg=%6.0f, std=%6.0f" % (
                                ota, cellx, celly, median_int, mean_int, std_int) )
                        
                    all_data.append(thiscell)

        hdulist.close()

        numpy.savetxt("alldata.tmp", all_data)
        stdout_write("\n\n")

    return all_data




def fit_nonlinearity_sequence(pinit, args):

    def fit_fct(p, x):
        y = numpy.zeros(x.shape)
        for i in range(p.shape[0]):
            y += p[i] * x**(i+1)
        return y
    def err_fct(p,x,y,err, fitrange_x, fitrange_y):
        yfit = fit_fct(p,x)
        in_fit_range = numpy.isfinite(x) & numpy.isfinite(y)
        if (not fitrange_x == None):
            in_fit_range = in_fit_range & (x >= fitrange_x[0]) & (x <= fitrange_x[1])
        if (not fitrange_y == None):
            in_fit_range = in_fit_range & (y >= fitrange_y[0]) & (y <= fitrange_y[1])
        if (err == None):
            return ((y-yfit))[in_fit_range]
        return ((y-yfit)/err)[in_fit_range]

    # (medlevel, exptime, None, intensity_range, exptime_range) = args

    fit = scipy.optimize.leastsq(err_fct, pinit, args=args, full_output=1)

    pfit = fit[0]
    uncert = numpy.sqrt(numpy.diag(fit[1]))

    return pfit, uncert

                          
def create_nonlinearity_fits(data, outputfits, polyorder=3, 
                             exptime_range=[0.1,2.5], intensity_range=[100,59000],
                             verbose=False):

    otas = set(data[:,0])
    #print otas

    result_count = 0 
    result_ota = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_cellx = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_celly = numpy.zeros(shape=(data.shape[0]), dtype=numpy.int)
    result_coeffs = numpy.zeros(shape=(data.shape[0],polyorder-1))
    result_coeffuncert = numpy.zeros(shape=(data.shape[0],polyorder-1))
    
    for ota in otas:
        for cellx in range(8):
            for celly in range(8):
                stdout_write("\rFitting OTA %02d, cell %1d,%1d ..." % (ota, cellx, celly))

                this_cell = (data[:,0] == ota) & (data[:,3] == cellx) & (data[:,4] == celly)

                subset = data[this_cell]
                #print ota, cellx, celly,":\n",subset

                not_nans = numpy.isfinite(subset[:,7]) & numpy.isfinite(subset[:,9])

                if (numpy.sum(not_nans) <= 0):
                    continue

                exptime = subset[:,5][not_nans]
                medlevel = subset[:,7][not_nans]
                stdlevel = subset[:,9][not_nans]

                pinit = numpy.zeros(polyorder)

                # fit = scipy.optimize.leastsq(err_fct, pinit,
                #                              args=(exptime, medlevel, stdlevel, exptime_ranges), 
                #                              full_output=1)

                # fit = scipy.optimize.leastsq(err_fct, pinit,
                #                              args=(medlevel, exptime, None, intensity_range, exptime_range), 
                #                              full_output=1)

                # pfit = fit[0]
                # uncert = numpy.sqrt(numpy.diag(fit[1]))
                args = (medlevel, exptime, None, intensity_range, exptime_range)
                pfit, uncert = fit_nonlinearity_sequence(pinit, args)



                if (verbose):
                    print ota, cellx, celly, pfit, uncert

                linear_factor = pfit[0]

                coefficients_normalized = (pfit / linear_factor)[1:]
                coefficient_errors_normalized = (uncert / linear_factor)[1:]
                
                result_ota[result_count] = ota
                result_cellx[result_count] = cellx
                result_celly[result_count] = celly
                result_coeffs[result_count] = coefficients_normalized
                result_coeffuncert[result_count] = coefficient_errors_normalized

                
                result_count += 1

    print "completed",result_count,"fits"
    stdout_write(" done!\n")

    # Prepare all data to be written to a fits file
    #print result_coeffs[:10,:]
    #print result_coeffuncert[:10,:]

    if (outputfits == None):
        return

    columns = [\
        pyfits.Column(name='OTA',    format='B', array=result_ota[:result_count], disp='ota'),
        pyfits.Column(name='CELLX',  format='B', array=result_cellx[:result_count], disp='cell-x'),
        pyfits.Column(name='CELLY',  format='B', array=result_celly[:result_count], disp='cell-y'),
        ]

    for i in range(polyorder-1):
        col = pyfits.Column(name="COEFF_X**%d" % (i+2), 
                            format='E',
                            array=result_coeffs[:result_count,i],
                            disp="polynomial coeff x^%d" % (i+2)
                            )
        columns.append(col)
    for i in range(polyorder-1):
        col = pyfits.Column(name="UNCERT_COEFF_X**%d" % (i+2), 
                            format='E',
                            array=result_coeffuncert[:result_count,i],
                            disp="uncertainty of polynomial coeff x^%d" % (i+2)
                            )
        columns.append(col)

    coldefs = pyfits.ColDefs(columns)
    tbhdu = pyfits.new_table(coldefs, tbtype='BinTableHDU')

    primhdu = pyfits.PrimaryHDU()
    primhdu.header.update("POLYORDR", polyorder)

    hdulist = pyfits.HDUList([primhdu, tbhdu])
    clobberfile(outputfits)
    hdulist.writeto(outputfits, clobber=True)
    
    return 



def load_nonlinearity_correction_table(filename, search_ota):
    
    # Load the catalog file
    hdulist = pyfits.open(filename)

    # Determine what the fittting order was
    polyorder = hdulist[0].header['POLYORDR']

    # Create an array holding all coefficients
    nonlinearity_coeffs = numpy.zeros(shape=(8,8,polyorder-1))

    # Now load the full catalog and sort the coefficients 
    # into the coefficient matrix
    ota = hdulist[1].data.field('OTA')
    cellx = hdulist[1].data.field('CELLX')
    celly = hdulist[1].data.field('CELLY')

    all_coeffs = numpy.zeros(shape=(ota.shape[0],polyorder-1))
    for order in range(polyorder-1):
        columnname = "COEFF_X**%d" % (order+2)
        all_coeffs[:,order] = hdulist[1].data.field(columnname)[:]

    in_this_ota = ota == search_ota

    cellx = cellx[in_this_ota]
    celly = celly[in_this_ota]
    all_coeffs = all_coeffs[in_this_ota]

    for i in range(cellx.shape[0]):
        nonlinearity_coeffs[cellx[i], celly[i], :] = all_coeffs[i]

    return nonlinearity_coeffs


def compute_nonlinearity_correction(data, coeffs):
    correction = numpy.zeros(data.shape)
    for i in range(coeffs.shape[0]):
        correction += coeffs[i] * data**(i+2)
    return correction
   
def compute_cell_nonlinearity_correction(data, cellx, celly, all_coeffs):
    coeffs = all_coeffs[cellx, celly, :]
    return compute_nonlinearity_correction(data, coeffs)


def create_data_fit_plot(data, fitfile, ota, cellx, celly, outputfile):

    import podi_plotting

    this_cell = (data[:,0] == ota) & (data[:,3] == cellx) & (data[:,4] == celly)
    subset = data[this_cell]

    # Determine the min and max exposure times
    exptime_min = 0
    exptime_max = 1.05 * numpy.max(subset[:,6])
    flux_min = 0
    flux_max = 70000

    medlevel = subset[:,8]
    exptimes = subset[:,6]
    #print medlevel

    fluxscaling = 1000

    # Load the fit parameters
    fittable = load_nonlinearity_correction_table(fitfile, ota)

    fig = matplotlib.pyplot.figure()
    left = 0.1
    width = 0.87

    ax1 = fig.add_axes([left, 0.3,width,0.64])
    ax2 = fig.add_axes([left, 0.1,width,0.2])

    ax1.set_xlim([flux_min, flux_max/fluxscaling])
    ax2.set_xlim([flux_min, flux_max/fluxscaling])
    ax1.set_ylim([exptime_min, exptime_max])

    ax1.scatter(medlevel/fluxscaling, exptimes, label="data")

    fit_x = numpy.linspace(flux_min, flux_max, 1000)
    delta_fit_y = compute_cell_nonlinearity_correction(fit_x, cellx, celly, fittable)
    fit_y = fit_x + delta_fit_y

    intensity_range = [0, 63000]
    exptime_range = [0, 1e9]
    args = (medlevel, exptimes, None, intensity_range, exptime_range)
    #print args


    def evaluate_poly(x, pol):
        y = numpy.zeros(x.shape)
        for i in range(pol.shape[0]):
            y += x**(i+1) * pol[i]
        return y

    colors = ('red', 'green', 'blue', 'grey')
    poly_fits = [None] * 3 #len(colors)

    error_range = [-0.6, 0.6]

    ax1.set_title("OTA %02d, cell %1d,%1d" % (ota, cellx, celly))
    for i in range(len(poly_fits)):
    
        fit, error = fit_nonlinearity_sequence(numpy.zeros(shape=(i+1)), args)
        poly_fits[i] = fit
                          
        label = "fit-order: %d" % (i+1)
        ax1.plot(fit_x/fluxscaling, evaluate_poly(fit_x, fit), label=label, c=colors[i])
        
        timediff = exptimes - evaluate_poly(medlevel, fit)
        within_errors = (timediff < error_range[1]) & (timediff > error_range[0])

        ax2.scatter(medlevel[within_errors]/fluxscaling, timediff[within_errors], c=colors[i])
        above = timediff > error_range[1]
        if (numpy.sum(above) > 0):
            med_above = medlevel[above]
            y_values = numpy.ones(shape=med_above.shape) * error_range[1]
            ax2.scatter(med_above/fluxscaling, y_values, c=colors[i], marker="^")
        below = timediff < error_range[0]
        if (numpy.sum(below) > 0):
            med_below = medlevel[below]
            y_values = numpy.ones(shape=med_below.shape) * error_range[0]
            ax2.scatter(med_below/fluxscaling, y_values, c=colors[i], marker="v")

        print "Order",i+1,":", fit

        if (i==2):
#            ax1.set_label()
            fig.text(0.93, 0.36, '3rd-order polynomial fit:\ny = x + %.4e*x^2 + %.4e*x^3' % (fit[1]/fit[0], fit[2]/fit[0]), 
                     horizontalalignment='right', verticalalignment='bottom')

    # Set maximum exposure time
    max_exptime_fit = 70000 * poly_fits[0][0]
    max_exptime_plot = numpy.max([max_exptime_fit, numpy.max(exptimes)])
    ax1.set_ylim([exptime_min, max_exptime_plot])


    ax1.legend(loc='upper left', borderaxespad=1)
    ax1.get_xaxis().set_ticklabels([]) #set_visible(False)
    ax1.set_ylabel("exposure time t_exp (~ true flux)")
    
    ax2.set_ylabel("delta t_exp")
    ax2.set_xlabel("observed flux level (x1000 cts)")
    ax2.set_ylim([-0.77,0.77])
    fig.savefig(outputfile)

    return




def create_nonlinearity_map(fitfile, outputfile, fluxlevel, minmax, labels=True):

    import podi_plotting

    # Load the catalog file
    hdulist = pyfits.open(fitfile)

    # Determine what the fittting order was
    polyorder = hdulist[0].header['POLYORDR']

    # Create an array holding all coefficients
    nonlinearity_coeffs = numpy.zeros(shape=(polyorder-1))

    # Now load the full catalog and sort the coefficients 
    # into the coefficient matrix
    ota = hdulist[1].data.field('OTA')
    cellx = hdulist[1].data.field('CELLX')
    celly = hdulist[1].data.field('CELLY')

    #print ota

    all_coeffs = numpy.zeros(shape=(ota.shape[0],polyorder-1))
    for order in range(polyorder-1):
        columnname = "COEFF_X**%d" % (order+2)
        all_coeffs[:,order] = hdulist[1].data.field(columnname)[:]


    cellsize = 0.12
    # Now loop over all OTAs and all cells and compute the corners of the cells

    all_corners = []
    all_intensity = []

    fig, ax = matplotlib.pyplot.subplots()

    data = numpy.array([fluxlevel])
    for cell in (range(ota.shape[0])):
        _ota_x = int(math.floor(ota[cell] / 10))
        _ota_y = int(math.fmod(ota[cell], 10))

        x1 = _ota_x + cellx[cell] * cellsize
        y1 = _ota_y + (7-celly[cell]) * cellsize

        corners = [[x1,y1], [x1+cellsize,y1], [x1+cellsize,y1+cellsize], [x1, y1+cellsize]]
        all_corners.append(corners)

        intensity = compute_nonlinearity_correction(data, all_coeffs[cell])[0] / fluxlevel
        #print ota[cell], _ota_x, _ota_y, cellx[cell], celly[cell], x1, y1, intensity
        all_intensity.append(intensity)

        #poly = Polygon(corners,facecolor='blue',edgecolor='none')
        #plt.gca().add_patch(poly)

        if (labels):
            intensity_text = "%.1f" % (math.fabs(intensity)*100)
            label_x = x1 + 0.5 * cellsize
            label_y = y1 + 0.5 * cellsize
            ax.text(label_x, label_y, intensity_text, fontsize=2,
                    horizontalalignment='center',
                    verticalalignment='center')

    #print all_corners
    #print all_intensity

    #cbar = matplotlib.pyplot.colorbar()
    #cbar.solids.set_edgecolor("face")
    #cbar.draw()

    cmap = matplotlib.pyplot.cm.get_cmap('spectral')

    #ax = fig.add_axes([0, 0, 1., 1.])
    corners = numpy.array(all_corners)
    
    ax.set_xlim([0,8])
    ax.set_ylim([0,8])
    
    #converter = matplotlib.colors.ColorConverter
    #colorvalues = cmap.to_rgb(all_intensity)
    #colorvalues = cmap(0.1)
    #colorvalues = [cm.jet(x) for x in np.random.rand(20)]
    #colorvalues = [matplotlib.pyplot.cm.jet(x) for x in all_intensity]
    
    nl_min = numpy.min(all_intensity[numpy.isfinite(all_intensity)]) if minmax[0] == None else float(minmax[0])
    nl_max = numpy.max(all_intensity[numpy.isfinite(all_intensity)]) if minmax[1] == None else float(minmax[1])
    
        
    colorvalues = cmap((numpy.array(all_intensity)-nl_min)/(nl_max-nl_min)) #[matplotlib.pyplot.cm.jet(x) for x in all_intensity]

    
    coll = matplotlib.collections.PolyCollection(corners, #facecolor='#505050', #
                                                 facecolor=colorvalues,
                                                 edgecolor='black', linestyle='-', linewidth=0.2,
                                                 cmap=matplotlib.pyplot.cm.get_cmap('spectral'),
                                                 )
    
    img = matplotlib.pyplot.imshow([[1e9],[1e9]], vmin=nl_min, vmax=nl_max, cmap=cmap, extent=(0,0,0,0), origin='lower')
    #colorbar = matplotlib.pyplot.colorbar(cmap=cmap)
    fig.colorbar(img) #, text="non-linearity")
    ax.set_title("Non-linearity @ %d counts" % (int(fluxlevel)))

    ax.set_xlim(-0.1,8.1)
    ax.set_ylim(-0.1,8.1)
    ax.add_collection(coll)
    #colorbar = matplotlib.pyplot.colorbar(cmap=cmap)
    fig.savefig(outputfile)

    return




def plot_cellbycell_map(fitfile, outputfile, minmax, labels=True, fontsize=2):

    import podi_plotting

    # Load the catalog file
    hdulist = pyfits.open(fitfile)

    cellsize = 0.12
    all_corners = []
    all_intensity = []
    # Now loop over all OTAs and all cells and compute the corners of the cells

    fig, ax = matplotlib.pyplot.subplots()

    for ext in range(1, len(hdulist)):
        for cellx in range(8):
            for celly in range(8):

                _ota_x = int(hdulist[ext].header['EXTNAME'][3])
                _ota_y = int(hdulist[ext].header['EXTNAME'][4])
                stdout_write("\rMeasuring OTA %d%d, cell %d,%d..." % (_ota_x, _ota_y, cellx, celly))

                x1 = _ota_x + cellx * cellsize
                y1 = _ota_y + (7-celly) * cellsize


                corners = [[x1,y1], [x1+cellsize,y1], [x1+cellsize,y1+cellsize], [x1, y1+cellsize]]
                all_corners.append(corners)

                cell_area = cell2ota__get_target_region(cellx, celly)
                cx1, cx2, cy1, cy2 = cell_area
                cell_data = hdulist[ext].data[cy1:cy2, cx1:cx2]

                cell_center = cell_data[bordersize:-bordersize,bordersize:-bordersize]

                intensity = numpy.median(cell_center)
                all_intensity.append(intensity)

                if (labels):
                    stdout_write(" %7.3f" % (intensity))
                    intensity_text = "%.2f" % (intensity)
                    label_x = x1 + 0.5 * cellsize
                    label_y = y1 + 0.5 * cellsize
                    ax.text(label_x, label_y, intensity_text, fontsize=1,
                            horizontalalignment='center',
                            color='white',
                            verticalalignment='center', zorder=99)

    cmap = matplotlib.pyplot.cm.get_cmap('spectral')

    corners = numpy.array(all_corners)
    
    ax.set_xlim([0,8])
    ax.set_ylim([0,8])
    
    nl_min = numpy.min(all_intensity[numpy.isfinite(all_intensity)]) if minmax[0] == None else float(minmax[0])
    nl_max = numpy.max(all_intensity[numpy.isfinite(all_intensity)]) if minmax[1] == None else float(minmax[1])
    
    colorvalues = cmap((numpy.array(all_intensity)-nl_min)/(nl_max-nl_min))

    coll = matplotlib.collections.PolyCollection(corners, #facecolor='#505050', #
                                                 facecolor=colorvalues,
                                                 edgecolor='black', linestyle='-', linewidth=0.2,
                                                 cmap=matplotlib.pyplot.cm.get_cmap('spectral'),
                                                 )
    
    img = matplotlib.pyplot.imshow([[1e9],[1e9]], vmin=nl_min, vmax=nl_max, cmap=cmap, extent=(0,0,0,0), origin='lower')
    fig.colorbar(img) 
    ax.set_title("Cell median intensity level")

    ax.set_xlim(-0.1,8.1)
    ax.set_ylim(-0.1,8.1)
    ax.add_collection(coll)
    fig.savefig(outputfile)

    return


if __name__ == "__main__":


    if (cmdline_arg_isset("-fit")):
        datafile = get_clean_cmdline()[1]
        
        outputfile = get_clean_cmdline()[2]

        data = numpy.loadtxt(datafile)

        min_exptime = float(cmdline_arg_set_or_default("-mint",   0.0))
        max_exptime = float(cmdline_arg_set_or_default("-maxt", 100.0))
        
        min_level = float(cmdline_arg_set_or_default("-minf",    10.0))
        max_level = float(cmdline_arg_set_or_default("-maxf", 59000.0))

        polyorder = int(cmdline_arg_set_or_default("-order", 3))

        verbose = cmdline_arg_isset("-verbose")

        stdout_write("""
Creating all fits
-----------------------------------------------------
  Polynomial order: %d
  Fit range exposure time [sec] = %.3f - %.3f
  Fit range intensity [counts]  = %d - %d
-----------------------------------------------------
""" % (polyorder, min_exptime, max_exptime, min_level, max_level))

        create_nonlinearity_fits(data, outputfile, polyorder=polyorder,
                                 intensity_range=[min_level,max_level],
                                 exptime_range=[min_exptime, max_exptime],
                                 verbose=verbose
                                 )
        stdout_write("\n")
        sys.exit(0)

    if (cmdline_arg_isset("-load")):
        fitsfile = get_clean_cmdline()[1]
        ota = int(get_clean_cmdline()[2])
        coeffs = load_nonlinearity_correction_table(fitsfile, ota)
        print coeffs
        sys.exit(0)

    if (cmdline_arg_isset("-correctraw")):
        infile = get_clean_cmdline()[1]
        hdulist = pyfits.open(infile)
        ota = int(hdulist[0].header['FPPOS'][2:4])
        catfile = get_clean_cmdline()[2]
        ota_coeffs = load_nonlinearity_correction_table(catfile, ota)
        
        for i in range(1, 65):
            stdout_write("\rcorrecting cell %d ..." % (i))
            data = hdulist[i].data
            overscan_level = numpy.median(data[:,494:])
            data -= overscan_level
            cellx = hdulist[i].header['WN_CELLX']
            celly = hdulist[i].header['WN_CELLY']
            correction = compute_cell_nonlinearity_correction(data, cellx, celly, ota_coeffs)

            #data += correction
            hdulist[i].data = data

        stdout_write(" writing ...")
        outfile = get_clean_cmdline()[3]
        hdulist.writeto(outfile, clobber=True)

        stdout_write(" done!\n\n")
        sys.exit(0)

    if (cmdline_arg_isset("-correct")):
        infile = get_clean_cmdline()[1]
        hdulist = pyfits.open(infile)

        catfile = get_clean_cmdline()[2]
        
        for ext in range(1, len(hdulist)):
            ota = int(hdulist[ext].header['EXTNAME'][3:5])
            ota_coeffs = load_nonlinearity_correction_table(catfile, ota)

            data = hdulist[ext].data
            for cx in range(0,8):
                for cy in range(0,8):

                    stdout_write("\rcorrecting ota %02d, cell %d,%d ..." % (ota, cx, cy))
                    x1, x2, y1, y2 = cell2ota__get_target_region(cx, cy)

                    data[y1:y2, x1:x2] += compute_cell_nonlinearity_correction(data[y1:y2, x1:x2], cx, cy, ota_coeffs)
                    # data -= overscan_level
                    # cellx = hdulist[i].header['WN_CELLX']
                    # celly = hdulist[i].header['WN_CELLY']
                    # correction = compute_cell_nonlinearity_correction(data, cellx, celly, ota_coeffs)

                    #data += correction
                    # hdulist[i].data = data
                    
            hdulist[ext].data = data
        stdout_write(" writing ...")
        outfile = get_clean_cmdline()[3]
        hdulist.writeto(outfile, clobber=True)

        stdout_write(" done!\n\n")
        sys.exit(0)

    if (cmdline_arg_isset("-plotdatafit")):
        datafile = get_clean_cmdline()[1]
        ota = int(get_clean_cmdline()[2])
        cellx = int(get_clean_cmdline()[3])
        celly = int(get_clean_cmdline()[4])
        fitfile = get_clean_cmdline()[5]
        outputfile = get_clean_cmdline()[6]

        data = numpy.loadtxt(datafile)
        create_data_fit_plot(data, fitfile, ota, cellx, celly, outputfile)
        sys.exit(0)

    if (cmdline_arg_isset("-nonlinmap")):
        fitfile = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]
        fluxlevel = float(get_clean_cmdline()[3])
        min_value = cmdline_arg_set_or_default("-min", None)
        max_value = cmdline_arg_set_or_default("-max", None)
        minmax = [min_value, max_value]
        labels = cmdline_arg_isset("-labels")
        create_nonlinearity_map(fitfile, outputfile, fluxlevel, minmax, labels=labels)
        sys.exit(0)

    if (cmdline_arg_isset("-cellbycellmap")):
        fitfile = get_clean_cmdline()[1]
        outputfile = get_clean_cmdline()[2]
        min_value = cmdline_arg_set_or_default("-min", None)
        max_value = cmdline_arg_set_or_default("-max", None)
        fontsize = float(cmdline_arg_set_or_default("-fontsize", 2))
        minmax = [min_value, max_value]
        labels = True #cmdline_arg_isset("-labels")
        plot_cellbycell_map(fitfile, outputfile, minmax, labels=labels, fontsize=fontsize)
        sys.exit(0)

 
    # Read the input directory that contains the individual OTA files
    inputfiles = get_clean_cmdline()[1:]
    all_data = create_nonlinearity_data(inputfiles)
    numpy.savetxt("nonlin.data", all_data)
