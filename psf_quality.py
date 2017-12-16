#!/usr/bin/env python


import os
import sys
import pyfits
import scipy
import scipy.optimize

import numpy
from matplotlib import pyplot

import podi_logging
from podi_collectcells import read_fits_catalog
from podi_definitions import three_sigma_clip



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

class PSFdata ( object ):

    def __init__(self, filename, pixelscale=None):
        self.filename = filename
        self.data = None
        self.fwhm = 0
        self.pixelscale = 0.11 if pixelscale is None else pixelscale


        self.calculate()
        self.fit_gauss()
        self.fit_moffat()


    def calculate(self):

        cat = read_fits_catalog(self.filename, 2, flatten=False)
        print len(cat)

        # TODO: Change column numbers
        mag = cat[4]
        mag_err = cat[5]
        vignets = cat[6]
        flags = cat[14]
        background = cat[15]
        fwhm = cat[18]
        elongation = cat[13]

        vignets[vignets < -1e29] = numpy.NaN

        flux = numpy.power(10., -0.4 * mag)
        print mag.shape, mag_err.shape, vignets.shape

        good = (mag < 0) & (mag_err < 0.1) & (flags == 0) & (fwhm >= 3)

        median_elongation = numpy.median(elongation)
        numpy.savetxt("elongation", elongation[good])
        _, good_elongation = three_sigma_clip(elongation[good], return_mask=True)

        numpy.savetxt("ellipticity", cat[20])
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
        # Prepare all vignet cutouts
        #
        # psfs = ((vignets - background.reshape((-1,1,1))) / flux.reshape((-1,1,1)))[good]
        psfs = (vignets / flux.reshape((-1, 1, 1)))[all_good]

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

    psf = PSFdata(cat_fn)
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



