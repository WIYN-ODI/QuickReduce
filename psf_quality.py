#!/usr/bin/env python


import os
import sys
import pyfits
import scipy
import scipy.optimize
import numpy
from matplotlib import pyplot
import logging

import podi_logging
from podi_readfitscat import read_fits_catalog
from podi_definitions import three_sigma_clip

from podi_definitions import SXcolumn


def make_psf_plot(ota_listing, title=None):

    fig = pyplot.figure()
    pixelscale = 0.11

    nx = 5
    ny = 6
    axes = fig.subplots(6,5,sharex=True,sharey=True)

    _y,_x = numpy.indices(ota_listing[33].data.shape, dtype=numpy.float)
    print _y.shape, _x.shape
    r = numpy.hypot((_x-32.), (_y-32.)) * pixelscale

    xmin, xmax, ymin, ymax = 0, 3.5, 5e-4, 2 #-6, -1.5
    axes[0,0].set_xlim((0, xmax))
    axes[0,0].set_ylim((ymin, ymax))

    plot_x = numpy.linspace(xmin, xmax, 100)

    for ota in ota_listing:

        otax = int(numpy.floor(ota/10.))
        otay = int(numpy.fmod(ota,10))
        print ota, otax, otay

        ax = axes[ny-otay,otax-1]
        # ax.imshow(ota_listing[ota],
        #           vmin=-0.001, vmax=0.01,
        #           extent=(xmin,xmax,ymin,ymax),
        #           alpha=0.2,
        #           )

        psf = ota_listing[ota]
        data = ota_listing[ota].data
        fwhm = ota_listing[ota].fwhm

        corrected_data = (data - psf.moffat_background) / psf.moffat_intensity

        #ax.scatter(r, numpy.log10(ota_listing[ota]), s=1)
        ax.semilogy(r, corrected_data, c='grey', marker=".",
                    markersize=3, alpha=0.3,
                    markeredgecolor='none', markerfacecolor='grey',
                    linestyle='None') #, size=1)

        ax.grid(color = '#f4f4f4', linestyle = 'solid', linewidth = 1)

        ax.semilogy(plot_x, psf.gaussprofile(plot_x)/psf.intensity)

        ax.text(0.9*xmax, 1., "%.2f" % (psf.fwhm),
                horizontalalignment='right',
                verticalalignment='top',
                multialignment='center'
                )

        print psf.moffat_fit
        ax.semilogy(plot_x, (psf.moffat(plot_x)-psf.moffat_background)/psf.moffat_intensity)

        #
        # matplotlib.pyplot.imshow(X, cmap=None, norm=None, aspect=None,
        #                          interpolation=None, alpha=None, vmin=None, vmax=None,
        #                          origin=None, extent=None, shape=None, filternorm=1, filterrad=4.0, imlim=None, resample=None, url=None, hold=None, data=None, **kwargs)

    # ax.set_xlim((0, 7.))
    # ax.set_ylim((0., 1.))

    ax.set_yticks([1e-3,1e-2,1e-1,1])
    ax.set_xticks(numpy.arange(0,3.5,0.5))

    if (title is not None):
        axes[0,2].set_title(title)

    fig.subplots_adjust(wspace=0, hspace=0)
    # fig.set_tight_layout(True)
    fig.set_size_inches(6,7)
    fig.show()
    fig.savefig("test.png", dpi=200, bbox_inches='tight')

    return


def moffat_model(p, r):
    model = p[0] + p[1] * (p[3] - 1) / (numpy.pi * p[2] ** 2) \
            * numpy.power(1. + r ** 2 / p[2] ** 2, -p[3])
    return model

def moffat_error(p, data, r):
    model = moffat_model(p, r)
    diff = data - model
    return diff.ravel()


def create_safe_cutout(
        image_data, x, y, wx, wy
):
    x = int(numpy.round(x, 0)) - 1
    y = int(numpy.round(y, 0)) - 1

    corner_min = numpy.array([(y-wy), (x-wx)])
    corner_max = numpy.array([(y+wy), (x+wx)])
    dimension = corner_max - corner_min
    cutout = numpy.zeros((dimension[0], dimension[1]))

    # cutout[:,:] = numpy.NaN
    # print "image dimension:", dimension

    trunc_min = numpy.max([corner_min, [0, 0]], axis=0)
    trunc_max = numpy.min([corner_max, image_data.shape], axis=0)
    # print "limited to valid area:", trunc_min, trunc_max

    # Now extract the data and insert it into the cutout
    insert_min = trunc_min - corner_min
    insert_max = insert_min + (trunc_max - trunc_min)

    # print "truncated area:", trunc_max - trunc_min

    # print "insert region:", insert_min, insert_max
    cutout[insert_min[0]:insert_max[0], insert_min[1]:insert_max[1]] = \
        image_data[trunc_min[0]:trunc_max[0], trunc_min[1]:trunc_max[1]]

    return cutout


class PSFquality (object):

    def __init__(self, catalog_filename, pixelscale=None,
                 catalog=None,
                 image_data = None, image_extension=None,
                 use_vignets=True):

        self.logger = logging.getLogger("compPSFmodel")

        self.catalog_filename = catalog_filename
        self.data = None
        self.fwhm = 0

        self.window_x = 32
        self.window_y = 32

        self.pixelscale = 0.11 if pixelscale is None else pixelscale

        self.use_vignets_from_catalog = use_vignets

        self.catalog_in = catalog
        self.cat = None
        self.cat_mag = None
        self.cat_magerr = None
        self.cat_elongation = None
        self.cat_fwhm = None
        self.cat_x = None
        self.cat_y = None
        self.cat_flags = None
        self.cat_background = None
        self.read_catalog()


        if (image_data is not None and
                type(image_data) is not numpy.ndarray):
            hdulist = pyfits.open(image_data)
            if (image_extension is None):
                image_data = hdulist[0].data
            else:
                image_data = hdulist[image_extension].data
        self.image_data = image_data

        self.compute()

    def read_catalog(self):

        if (self.catalog_in is None):
            self.logger.info("Reading catalog from %s" % (self.catalog_filename))
            if (self.use_vignets_from_catalog):
                self.cat = read_fits_catalog(self.catalog_filename, 2, flatten=False)
            else:
                self.cat = read_fits_catalog(self.catalog_filename, 2, flatten=True)
        else:
            self.logger.info("Using user-supplied catalog from memory")
            self.cat = self.catalog_in

        print self.cat

        if (self.cat_x is None):
            self.cat_x = self.cat[:, SXcolumn['x']]
        if (self.cat_y is None):
            self.cat_y = self.cat[:, SXcolumn['y']]
        if (self.cat_mag is None):
            self.cat_mag = self.cat[:, SXcolumn['mag_auto']]
        if (self.cat_magerr is None):
            self.cat_magerr = self.cat[:, SXcolumn['mag_err_auto']]
        if (self.cat_elongation is None):
            self.cat_elongation = self.cat[:, SXcolumn['elongation']]
        if (self.cat_fwhm is None):
            self.cat_fwhm = self.cat[:, SXcolumn['fwhm_image']]
        if (self.cat_flags is None):
            self.cat_flags = self.cat[:, SXcolumn['flags']]
        if (self.cat_background is None):
            self.cat_background = self.cat[:, SXcolumn['background']]

    def compute(self):
        self.calculate_composite_PSF()
        self.fit_gauss()
        self.fit_moffat()

    def calculate_composite_PSF(self):


        cat = self.cat
        print len(cat)

        #
        # Read the relevant columns from the catalog
        #
        # TODO: Change column numbers
        mag = self.cat_mag #$cat[4]
        mag_err = self.cat_magerr #cat[5]
        flags = self.cat_flags #cat[14]
        background = self.cat_background #cat[15]
        fwhm = self.cat_fwhm #cat[18]
        elongation = self.cat_elongation #cat[13]

        #
        # Select suitable stars that are not blended, not too compact, and
        # that have good photometry
        #
        flux = numpy.power(10., -0.4 * mag)
        print mag.shape, mag_err.shape

        good = (mag < 0) & (mag_err < 0.1) & (flags == 0) & (fwhm >= 3)

        median_elongation = numpy.median(elongation)
        numpy.savetxt("elongation", elongation[good])
        _, good_elongation = three_sigma_clip(elongation[good], return_mask=True)

        # numpy.savetxt("ellipticity", cat[20])
        print "median elongation", median_elongation

        valid_fwhm = fwhm[good]
        # print valid_fwhm
        numpy.savetxt("fwhms", valid_fwhm)

        clipped, starlike = three_sigma_clip(input=valid_fwhm, return_mask=True)
        valid_fwhm[~starlike] = 0
        numpy.savetxt("fwhms2", valid_fwhm)

        all_good = good.copy()
        all_good[all_good] &= (starlike & good_elongation)

        merged = numpy.array([
            mag, mag_err, fwhm, elongation,
        ]).T
        numpy.savetxt("data_all", merged)
        merged[~all_good, :] = numpy.NaN
        numpy.savetxt("data_good", merged)

        print numpy.sum(good)

        #
        # Now we know what stars to include in the composite PSF
        # Prepare all vignet cutouts
        #

        if (self.use_vignets_from_catalog):
            vignets = cat[6]
            vignets[vignets < -1e29] = numpy.NaN

            psfs = (vignets / flux.reshape((-1, 1, 1)))[all_good]

        else:
            # self.logger.critical("The mode extracting cutouts from image is not implemented yet")
            pos_x = self.cat_x[all_good]  #cat[7][all_good]
            pos_y = self.cat_y[all_good]  #cat[8][all_good]
            bg = background[all_good]

            print pos_x.shape, pos_y.shape, flux.shape
            n_vignets = pos_x.shape[0]

            vignets = numpy.empty((n_vignets, 2*wy, 2*wx))
            for i in range(n_vignets):
                print pos_x, pos_y
                x,y = pos_x[i], pos_y[i]
                cutout = create_safe_cutout(
                    image_data=self.image_data,
                    x=x, y=y,
                    wx=self.window_x, wy=self.window_y,
                )
                print vignets.shape, cutout.shape
                vignets[i, :, :] = cutout[:,:] - bg[i]

            psfs = vignets / flux[all_good].reshape((-1,1,1))

        out_hdu = [pyfits.PrimaryHDU()]
        for _i, i in enumerate(psfs):
            img = pyfits.ImageHDU(data=i)
            img.header['OBJECT'] = "M=%.3f +/- %.3f / FWHM=%.1f / Elong=%.2f" % (
                mag[all_good][_i], mag_err[all_good][_i],
                fwhm[all_good][_i], elongation[all_good][_i]
            )
            out_hdu.append(img)

        pyfits.HDUList(out_hdu).writeto("psfs.fits", clobber=True)

        #
        # combine the remaining good PSFs
        #
        good_psf = numpy.ones((psfs.shape[0]), dtype=numpy.bool)
        psf_flux = flux[all_good]
        psf_weights = numpy.ones_like(psfs) * psf_flux.reshape((-1, 1, 1))
        psf_weights[~numpy.isfinite(psfs)] = 0

        # print good_psf
        for iter in range(1):
            combined_psf = numpy.nanmedian(psfs[good_psf], axis=0)

            # weighted = numpy.sum(psfs*psf_flux.reshape((-1,1,1)), axis=0) / \
            #         numpy.sum(psf_flux)
            weighted = numpy.nansum(psfs * psf_weights, axis=0) / \
                       numpy.sum(psf_weights, axis=0)

        print combined_psf.shape
        pyfits.PrimaryHDU(data=combined_psf).writeto("median_psf.fits", clobber=True)
        pyfits.PrimaryHDU(data=weighted).writeto("weighted_psf.fits", clobber=True)

        self.data = combined_psf

        self.y, self.x = numpy.indices(combined_psf.shape, dtype=numpy.float)
        self.x -= 0.5 * combined_psf.shape[1]
        self.y -= 0.5 * combined_psf.shape[0]
        self.r = numpy.hypot(self.x, self.y) * self.pixelscale


    def fit_gauss(self):

        def gauss_model(p, r):
            model = p[0] + p[1]*numpy.exp(-r**2/(2*p[2]**2))
            return model

        def gauss_error(p, data, r):
            model = gauss_model(p, r)
            diff = data - model
            return diff.ravel()

        p_init = [0., 0.02, 0.3]
        fit = scipy.optimize.leastsq(gauss_error, p_init,
                                     args=(self.data, self.r),
                                     full_output=True)
        #print fit

        best_fit = fit[0]
        self.background = best_fit[0]
        self.intensity = best_fit[1]
        self.gauss_sigma = best_fit[2]

        self.fwhm = self.gauss_sigma * 2.35482
        print "GAUSS:", self.fwhm, self.gauss_sigma


        return

    def fit_moffat(self):

        p_init = [0., 0.02, 1., 1.]
        fit = scipy.optimize.leastsq(moffat_error, p_init,
                                     args=(self.data, self.r),
                                     full_output=True)
        #print fit

        best_fit = fit[0]
        self.moffat_background = best_fit[0]
        self.moffat_intensity = best_fit[1]
        self.moffat_alpha = best_fit[2]
        self.moffat_beta = best_fit[3]
        self.moffat_fit = best_fit

        return


    def gaussprofile(self, r):
        return self.intensity * numpy.exp(-r**2/(2*self.gauss_sigma**2))

    def moffat(self, r):
        return moffat_model(self.moffat_fit, r)

    def info(self, logger=None):
        if (logger is None):
            logger = self.logger
        logger.info("PSF-quality: size: %dx%d, #frames=%d, FWHM=%.2f" % (
            self.window_x, self.window_y, -1, self.fwhm,
        ))

    def save2fits(self, fn):
        gauss = self.gaussprofile(self.r)
        moffat = self.moffat(self.r)
        out_list = [
            pyfits.PrimaryHDU(),
            pyfits.ImageHDU(data=self.data, name="DATA"),
            pyfits.ImageHDU(data=gauss, name="GAUSS"),
            pyfits.ImageHDU(data=moffat, name="MOFFAT"),
            pyfits.ImageHDU(data=(self.data-gauss), name="GAUSS_RESIDUALS"),
            pyfits.ImageHDU(data=(self.data-moffat), name="MOFFAT_RESIDUALS"),
        ]
        hdulist = pyfits.HDUList(out_list)
        hdulist.writeto(fn, clobber=True)


if __name__ == "__main__":

    options = {}
    podi_logging.setup_logging(options)

    cat_fn = sys.argv[1]

    try:
        image_fn = sys.argv[2]
        image_hdu = pyfits.open(image_fn)
        image_data = image_hdu[1].data
        use_vignets = False
    except:
        use_vignets = True
        image_data = None

    psf = PSFquality(cat_fn, image_data=image_data,
                     use_vignets=use_vignets)


    psf.save2fits("psfmodels.fits")

    #
    # Now we have a proper PSF, let's make some plots
    #
    ota_data = {
        33: psf,
        11: psf,
        22: psf,
        16: psf,
    }
    make_psf_plot(ota_data, title="demo plot")

    podi_logging.shutdown_logging(options)



