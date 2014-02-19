/**
 *
 * (c) Ralf Kotulla, kotulla@uwm.edu
 *
 * This module implements a iterative sigma-clipping routine to speed
 * up the corresponding functionality in podi_imcombine.
 *
 */

#define true 0
#define false 1
typedef int bool;

#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define nan 0./0.

// #define printf //


double* blkrep(double* data, int sx, int sy, int facx, int facy)
{
    // allocate memory for the larger frame
    double* larger = (double*) malloc(sx*sy*facx*facy*sizeof(double));

    // Now copy the data, duplicating each pixel
    int newx = sx * facx;
    int newy = sy * facy;
    int lx,ly, _x, _y;
    
    for (lx=0; lx<newx; lx++) {
        for (ly=0; ly<newy; ly++) {
            _x = (int)( (double)lx / (double)facx );
            _y = (int)( (double)ly / (double)facy );
            
            larger[ly * newx + lx] = data[_y * sx + _x];
        }
    }

    return larger;
}



void convolve(double* input, int sx, int sy,
              double* output,
              double* kernel, int ksize)
{
    int kernel_center = (ksize-1)/2;

    int dx, dy, kx, ky, i;

    /* printf("kernel center: %d\n",kernel_center); */
    /* printf("Kernel\n%d %d %d\n%d %d %d\n%d %d %d\n", */
    /*        kernel[ 1+kernel_center + (-1+kernel_center)*ksize], kernel[ 1+kernel_center + ( 0+kernel_center)*ksize], kernel[ 1+kernel_center + (+1+kernel_center)*ksize], */
    /*        kernel[ 0+kernel_center + (-1+kernel_center)*ksize], kernel[ 0+kernel_center + ( 0+kernel_center)*ksize], kernel[ 0+kernel_center + (+1+kernel_center)*ksize], */
    /*        kernel[-1+kernel_center + (-1+kernel_center)*ksize], kernel[-1+kernel_center + ( 0+kernel_center)*ksize], kernel[-1+kernel_center + (+1+kernel_center)*ksize]); */
    /* printf("kernel_indices=\n %d %d %d\n %d %d %d\n %d %d %d\n", */
    /*        ( 1+kernel_center + (-1+kernel_center)*ksize), ( 1+kernel_center + ( 0+kernel_center)*ksize), ( 1+kernel_center + (+1+kernel_center)*ksize), */
    /*        ( 0+kernel_center + (-1+kernel_center)*ksize), ( 0+kernel_center + ( 0+kernel_center)*ksize), ( 0+kernel_center + (+1+kernel_center)*ksize), */
    /*        (-1+kernel_center + (-1+kernel_center)*ksize), (-1+kernel_center + ( 0+kernel_center)*ksize), (-1+kernel_center + (+1+kernel_center)*ksize)); */

    /* int idx; */
    /* for (kx = -1*kernel_center; kx<= kernel_center; kx++) { */
    /*     for (ky = -1*kernel_center; ky<= kernel_center; ky++) { */
    /*         idx = ky+kernel_center + (kx+kernel_center)*ksize; */
    /*         printf("kernel-idx: %d, %d --> %d = %f\n", kx, ky, idx, kernel[idx]); */
    /*     } */
    /* } */
    
    // Execute the convolution
    for (dx=0; dx<sx; dx++) {
        for (dy=0; dy<sy; dy++) {

            // Set the output pixel to 0 to start with
            output[dy + dx*sy] = 0.0;

            //
            // Compute the result for this pixel convolved with the kernel
            //
            for (kx = -1*kernel_center; kx<= kernel_center; kx++) {
                for (ky = -1*kernel_center; ky<= kernel_center; ky++) {
            /* for (kx = 0; kx<1; kx++) { */
            /*     for (ky = 0; ky<1; ky++) { */

                    if (dx+kx < 0 || dx+kx >= sx || dy+ky < 0 || dy+ky > sy) {
                        continue;
                    }

                    output[dy + dx*sy] += input[dy+ky + (dx+kx)*sy]
                        * kernel[ky+kernel_center + (kx+kernel_center)*ksize];
                    
                }
            }
            
        }
    }

    return;
}


void imreplace(double* data, int sx, int sy, double min, double max, double value)
{
    //printf("imreplace\n")

    int i;
    for (i=0; i<sx*sy; i++) {
        if (data[i] >= min && data[i] <= max) {
            data[i] = value;
        }
    }

    return;
}

    
double* blkavg(double* data, int sx, int sy, double* smaller, int facx, int facy)
{

    // Allocate memory for the smaller/binned output frame
    int n_pixels_new = sx*sy/(facx*facy), i;
    int* pixelcount = (int*)malloc(n_pixels_new*sizeof(int));
    
    // Set all values to 0
    for (i=0; i<n_pixels_new; i++) {
        smaller[i] = 0;
        pixelcount[i] = 0;
    }
    
    // Now copy the data
    int newx = sx / facx;
    int newy = sy / facy;
    int lx,ly, _x, _y;
    
    for (lx=0; lx<sx; lx++) {
        for (ly=0; ly<sy; ly++) {
            _x = (int)( (double)lx / (double)facx );
            _y = (int)( (double)ly / (double)facy );
            
            smaller[_x * newy + _y] += data[lx * sy + ly];
            pixelcount[_x * newy + _y]++;
        }
    }

    // Compute the average value by diving the sum by the
    // number of contribing pixels
    for (i=0; i<n_pixels_new; i++) {
        if (pixelcount[i] > 0) {
            smaller[i] /= pixelcount[i];
        }
    }

    // free the memory of the pixelcount array
    free(pixelcount);
    
    return smaller;
}

double* median(double* data, int sx, int sy, int wx, int wy)
{
    double* output = (double*)malloc(sx*sy*sizeof(double));

    double* neighbors = (double*)malloc(wx*wy*sizeof(double));
    
    int x, y, i, n, _x, _y, dx, dy, valid_pixel;

    dx = (wx-1)/2;
    dy = (wy-1)/2;
    
    for (x=0; x<sx; x++) {
        for (y=0; y<sy; y++) {

            // Collect all pixels in the neighborhood of this pixel
            n = 0;
            for (_x=x-dx; _x<=x+dx; _x++) {
                if (_x < 0 || _x >= sx) {
                    continue;
                }
                for (_y=y-dy; _y <= y+dy; _y++) {
                    if (_y < 0 || _y >= sy) {
                        continue;
                    }
                    neighbors[n++] = data[_y + _x*sy];
                }
            }

            // Now compute the median
            gsl_sort(neighbors, 1, n);
            output[y + x*sy] = gsl_stats_median_from_sorted_data(neighbors, 1, n);
        }
    }
    
    return output;
}

double* median_filtered(double* data, int sx, int sy, int wx, int wy, double* mask)
{
    double* output = (double*)malloc(sx*sy*sizeof(double));

    double* neighbors = (double*)malloc(wx*wy*sizeof(double));
    
    int x, y, i, n, _x, _y, dx, dy, valid_pixel;

    dx = (wx-1)/2;
    dy = (wy-1)/2;
    
    for (x=0; x<sx; x++) {
        for (y=0; y<sy; y++) {

            // Collect all pixels in the neighborhood of this pixel
            n = 0;
            for (_x=x-dx; _x<=x+dx; _x++) {
                if (_x < 0 || _x >= sx) {
                    continue;
                }
                for (_y=y-dy; _y <= y+dy; _y++) {
                    if (_y < 0 || _y >= sy) {
                        continue;
                    }
                    if (mask[_y + _x*sy] == 0) {
                        neighbors[n++] = data[_y + _x*sy];
                    }
                }
            }

            // Now compute the median
            gsl_sort(neighbors, 1, n);
            output[y + x*sy] = gsl_stats_median_from_sorted_data(neighbors, 1, n);
        }
    }
    
    return output;
}

double* get_noise_model(double* data, int sx, int sy, double gain, double readnoise)
{
    int i;
    int n_pixels = sx*sy;
    double *noise = (double*)malloc(n_pixels*sizeof(double));
    
    for (i=0; i<n_pixels; i++) {
        noise[i] = sqrt( data[i]*gain + readnoise*readnoise);
    }

    return noise;
}




double* imarith(double* a, char op, double* b, int sx, int sy)
{
    int n_pixels = sx * sy, i;
    double* output = (double*)malloc(n_pixels*sizeof(double));

    if (op == '*') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] * b[i];
        }
    }

    else if (op == '/') {
        for (i=0;i<n_pixels; i++) {
            output[i] = b[i] != 0 ? a[i] / b[i] : 0;
        }
    }
    
    else if (op == '+') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] + b[i];
        }
    }
       
    else if (op == '-') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] - b[i];
        }
    }
       
    else {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i];
        }

    }
    
    return output;
}

double* imarith_const(double* a, char op, double b, int sx, int sy)
{
    int n_pixels = sx * sy, i;
    double* output = (double*)malloc(n_pixels*sizeof(double));

    if (op == '*') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] * b;
        }
    }

    else if (op == '/') {
        for (i=0;i<n_pixels; i++) {
            output[i] = b != 0 ? a[i] / b : 0;
        }
    }
    
    else if (op == '+') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] + b;
        }
    }
       
    else if (op == '-') {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i] - b;
        }
    }
       
    else {
        for (i=0;i<n_pixels; i++) {
            output[i] = a[i];
        }

    }
    
    return output;
}

void imdel(double* data) 
{
    free(data);
    return;
}

#define MAXMEDIAN 7

void lacosmics__cy(double* data, double* output, double* output2,
                   int sx, int sy,
                   double gain, double readnoise,
                   double sigclip, double sigthres, double objlim,
                   int niter)
{

    double sigfrac = 0.25;
    
    printf("test\n");
    printf("Gain=%f\n",gain);
    printf("readnoise=%f\n",readnoise);

    int lx, ly, _x, _y, i, wx, wy, n, x, y;
    int kernel_center = 1, ksize=3;
    int dx, dy, kx, ky;
    double tmpd;
    int sx2 = sx*2, sy2 = sy*2;

    
    double* larger_2x2 = (double*)malloc(sx*2*sy*2*sizeof(double));
    double* lapla_convolved = (double*)malloc(sx*2*sy*2*sizeof(double));
    double* deriv2 = (double*)malloc(sx*sy*sizeof(double));
    double* med5 = (double*)malloc(sx*sy*sizeof(double));
    double* noise = (double*)malloc(sx*sy*sizeof(double));
    double* sigmap = (double*)malloc(sx*sy*sizeof(double));
    double* sigmap_med5 = (double*)malloc(sx*sy*sizeof(double));
    double* firstsel = malloc(sx*sy*sizeof(double));
    double* data_med3 = (double*)malloc(sx*sy*sizeof(double));
    double* data_med7 = (double*)malloc(sx*sy*sizeof(double));
    
    int* blkavg_pixelcount = (int*)malloc(sx*sy*sizeof(int));
    int* pixel_changed = (int*)malloc(sx*sy*sizeof(int));
    
    double* neighbors = (double*)malloc(MAXMEDIAN*MAXMEDIAN*sizeof(double));
    
 /*     # create Laplacian kernel */
 /* print("0 -1 0;",>> kernel) */
 /* print("-1 4 -1;",>> kernel) */
 /* print("0 -1 0;",>> kernel) */

    double laplace_kernel[3*3] = {  
          0., -1.,  0. ,
         -1. , 4., -1. ,
          0., -1.,  0. 
    };
    
 /*  # create growth kernel */
 /* print("1 1 1;",>> gkernel) */
 /* print("1 1 1;",>> gkernel) */
 /* print("1 1 1;",>> gkernel) */

    double growth_kernel[3*3] = {
         1., 1., 1. ,
         1., 1., 1. ,
         1., 1., 1. ,
    };
    
    printf("laplace-kernel test%f\n", laplace_kernel[1+1*3]);
 

    //
    // Add here: determine gain
    //


  /*   # take second-order derivative (Laplacian) of input image */
  /* # kernel is convolved with subsampled image, in order to remove negative */
  /* # pattern around high pixels */

  /* if (verbose) { */
  /*  print("Convolving image with Laplacian kernel") */
  /*  print("") */
  /*  } */
  /* blkrep(oldoutput,blk,2,2) */
  /* convolve(blk,lapla,kernel) */
  /* imreplace(lapla,0,upper=0,lower=INDEF,radius=0) */
  /* blkavg(lapla,deriv2,2,2,option="average") */


    //
    // Perform the laplace filtering
    //
    printf("Working on frame with dimensions %d x %d\n", sx, sy);

    int iteration = 0;
    niter = 1;

    
    for (iteration = 0; iteration < niter; iteration++) {
        
        // duplicate all pixels 2x2
        printf("Computing larger 2x2 blkrep array\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i] || iteration == 0) {
                    larger_2x2[ _y*2   + (_x*2  )*sy*2 ] = data[i];
                    larger_2x2[ _y*2+1 + (_x*2  )*sy*2 ] = data[i];
                    larger_2x2[ _y*2   + (_x*2+1)*sy*2 ] = data[i];
                    larger_2x2[ _y*2+1 + (_x*2+1)*sy*2 ] = data[i];
                }
                
            }
        }

        //
        //
        //
        printf("Convolving with 3x3 laplace kernel!\n");
        /* convolve(larger_2x2, 2*sx, 2*sy, lapla_convolved,  laplace_kernel, 3); */
        for (dx=0; dx<sx2; dx++) {
            for (dy=0; dy<sy2; dy++) {
                i = dy + dx*sy2;
                if (pixel_changed[dy/2 + dx/2*sy] || iteration == 0) {
                    
                    lapla_convolved[i] = 0.0;

                    //
                    // Compute the result for this pixel convolved with the kernel
                    //
                    for (kx = -1*kernel_center; kx<= kernel_center; kx++) {
                        for (ky = -1*kernel_center; ky<= kernel_center; ky++) {
                            if (!(dx+kx < 0 || dx+kx >= sx2 || dy+ky < 0 || dy+ky > sy2)) {

                                lapla_convolved[i] += larger_2x2[dy+ky + (dx+kx)*sy2]
                                    * laplace_kernel[ky+kernel_center + (kx+kernel_center)*ksize];
                            }
                        }
                    }
                    /* printf("replacing negative pixel values\n"); */
                    /* imreplace(lapla_convolved, sx*2, sy*2, -1e99, 0., 0.0); */
                    lapla_convolved[i] = lapla_convolved[i] < 0 ? 0. : lapla_convolved[i];
                }
                
            }
        }
        for (i=0; i<sx*sy*4; i++) {
            output2[i] = lapla_convolved[i];
        }
   

        printf("Running blkavg\n");
        /* blkavg(lapla_convolved, sx*2, sy*2, deriv2, 2, 2); */
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i] || iteration == 0) {
                
                    deriv2[i] = 0.25 * (
                        lapla_convolved[ _y*2   + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2   + (_x*2+1)*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2+1)*sy*2 ]
                                        );
                }
            }
        }
    
        for (i=0; i<sx*sy; i++) {
            output[i] = deriv2[i];
        }

        

        /* # create model of background flux - 5x5 box should exclude all CRs */
        /* median(oldoutput,med5,5,5,zlo=INDEF,zhi=INDEF,verb-) */
        /* imreplace(med5,0.0001,upper=0,lower=INDEF,radius=0) */
        printf("Median-filtering the data\n");
        /* double* med5 = median(data, sx, sy, 5, 5); */
        /* printf("Replacing all negative values with a small number\n"); */
        /* imreplace(med5, sx, sy, -1e99, 0.0, 0.0001); */

        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (x=0; x<sx; x++) {
            for (y=0; y<sy; y++) {
                i = y + x*sy;

                // only do work if necessary
                if (pixel_changed[i] || iteration == 0) {
                
                    // Collect all pixels in the neighborhood of this pixel
                    n = 0;
                    for (_x=x-dx; _x<=x+dx; _x++) {
                        if (_x < 0 || _x >= sx) continue;
                        for (_y=y-dy; _y <= y+dy; _y++) {
                            if (_y < 0 || _y >= sy) continue;
                            neighbors[n++] = data[_y + _x*sy];
                        }
                    }
                    // Now compute the median
                    gsl_sort(neighbors, 1, n);
                    tmpd = gsl_stats_median_from_sorted_data(neighbors, 1, n);
                    med5[i] = tmpd < 0 ? 0.0 : tmpd;

                    // Compute noise estimate
                    noise[i] = sqrt(med5[i]*gain + readnoise*readnoise);

                    // Compute significance of pixel
                    sigmap[i] = (deriv2[i] / noise[i]) / 2.;
                }
            }
        }

        /* /\* # create noise model *\/ */
        /* /\* imcalc(med5,noise,"sqrt(im1*"//usegain//" + "//readn//"**2)/"//usegain,verb-) *\/ */
        /* printf("computing noise model\n"); */
        /* double* noise = get_noise_model(med5, sx, sy, gain, readnoise); */

        /* # divide Laplacian by noise model */
        /* imarith(deriv2,"/",noise,sigmap) */
        /* printf("divide Laplacian by noise model\n"); */
        /* double* sigmap_x2 = imarith(deriv2, '/', noise, sx, sy); */

        /* # Laplacian of blkreplicated image counts edges twice: */
        /* imarith(sigmap,"/",2.,sigmap) */
        /* printf("Laplacian of blkreplicated image counts edges twice\n"); */
        /* double* sigmap = imarith_const(sigmap_x2, '/', 2., sx, sy); */
        /* free(sigmap_x2); */



    
    
        /* # removal of large structure (bright, extended objects) */
        /* imdel(med5) */
        /* median(sigmap,med5,5,5,zlo=INDEF,zhi=INDEF,verb-) */
        /* imarith(sigmap,"-",med5,sigmap) */
        printf("removing large structure\n");
    
        /* imdel(med5); */
        /* double * sigmap_med5 = median(sigmap, sx, sy, 5, 5); */
        /* double * tmp = imarith(sigmap, '-', sigmap_med5, sx, sy); */
        /* free(sigmap); */
        /* sigmap = tmp; */
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (x=0; x<sx; x++) {
            for (y=0; y<sy; y++) {
                i = y + x*sy;
                
                // only do work if necessary
                if (pixel_changed[i] || iteration == 0) {
                    n = 0;
                    for (_x=x-dx; _x<=x+dx; _x++) {
                        if (_x < 0 || _x >= sx) continue;
                        for (_y=y-dy; _y <= y+dy; _y++) {
                            if (_y < 0 || _y >= sy) continue;
                            neighbors[n++] = sigmap[_y + _x*sy];
                        }
                    }
                    // Now compute the median
                    gsl_sort(neighbors, 1, n);
                    sigmap_med5[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n);
                    // Subtract the smoothed significance map from the pixel significance map
                    sigmap[i] -= sigmap_med5[i];
                }
            }
        }
        

                
        /* # find all candidate cosmic rays */
        /* # this selection includes sharp features such as stars and HII regions */

        /* if (verbose) { */
        /*  print("Selecting candidate cosmic rays") */
        /*  print("  sigma limit = "//sigclip) */
        /*  print("") */
        /*  } */
        /* imcopy(sigmap,firstsel,verb-) */
        /* imreplace(firstsel,0,upper=sigclip,lower=INDEF,radius=0) */
        /* imreplace(firstsel,1,lower=0.1,upper=INDEF,radius=0) */
        printf("Selecting candidate CRs\n");
        for (i=0; i<sx*sy; i++) {
            firstsel[i] = sigmap[i] > sigclip ? 1 : 0;

        }
    
        /* for (i=0; i<sx*sy; i++) { */
        /*     output[i] = firstsel[i]; */
        /* } */
        /* return; */

        /* # compare candidate CRs to median filtered image */
        /* # this step rejects bright, compact sources from the initial CR list */

        /* if (verbose) { */
        /*  print("Removing suspected compact bright objects (e.g. stars)") */
        /*  print("  selecting cosmic rays > "//objlim//" times object flux") */
        /*  print("") */
        /*  } */

        /* # subtract background and smooth component of objects */
        /* median(oldoutput,med3,3,3,zlo=INDEF,zhi=INDEF,verb-) */
        /* median(med3,med7,7,7,zlo=INDEF,zhi=INDEF,verb-) */
        /* imarith(med3,"-",med7,med3) */
        /* imarith(med3,"/",noise,med3) */
        /* imreplace(med3,0.01,upper=0.01,lower=INDEF,radius=0) */
        printf("subtract background and smooth component of objects\n");

        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;

                if (pixel_changed[i] || iteration == 0) {
                    //
                    // do 7x7 median filtering
                    //
                    n=0;
                    for ( x = ( (_x-3) <  0 ?  0 : (_x-3) );
                          x < ( (_x+4) > sx ? sx : (_x+4) );
                          x++) {
                        for ( y = ( (_y-3) <  0 ?  0 : (_y-3) );
                              y < ( (_y+4) > sy ? sy : (_y+4) );
                              y++) {
                            neighbors[n++] = data[y + x*sy];
                        }
                    }
                    gsl_sort(neighbors, 1, n);
                    data_med7[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n);


                    //
                    // do 3x3 median filtering
                    //
                    n=0;
                    for ( x = ( (_x-1) <  0 ?  0 : (_x-1) );
                          x < ( (_x+2) > sx ? sx : (_x+2) );
                          x++) {
                        for ( y = ( (_y-1) <  0 ?  0 : (_y-2) );
                              y < ( (_y+2) > sy ? sy : (_y+3) );
                              y++) {
                            neighbors[n++] = data[y + x*sy];
                        }
                    }
                    gsl_sort(neighbors, 1, n);
                    data_med3[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n);

                } // end if pixel_changed

            }
        }


        /* double* data_med3 = median(data, sx, sy, 3, 3); */
        /* double* data_med7 = median(data_med3, sx, sy, 7, 7); */
        for (i=0; i<sx*sy; i++){
            tmpd = (data_med3[i] - data_med7[i]) / noise[i];
            data_med3[i] = tmpd < 0.01 ? 0.01 : tmpd;

            /* # compare CR flux to object flux */
            /* imcalc(firstsel//","//sigmap//","//med3,starreject,"(im1*im2)/im3",verb-) */

            /* # discard if CR flux <= objlim * object flux */
            /* imreplace(starreject,0,upper=objlim,lower=INDEF,radius=0) */
            /* imreplace(starreject,1,lower=0.5,upper=INDEF,radius=0) */
            /* imarith(firstsel,"*",starreject,firstsel) */
 
            if (firstsel[i] > 0) {
                firstsel[i] = sigmap[i] > data_med3[i]  * objlim ? 1 : 0;
            }

            // Also reset the mask of CR pixels to 0
            pixel_changed[i] = 0;
        }

        /* # grow CRs by one pixel and check in original sigma map */
        /* convolve(firstsel,gfirstsel,gkernel) */
        /* imreplace(gfirstsel,1,lower=0.5,upper=INDEF,radius=0) */
        /* imarith(gfirstsel,"*",sigmap,gfirstsel) */
        /* imreplace(gfirstsel,0,upper=sigclip,lower=INDEF,radius=0) */
        /* imreplace(gfirstsel,1,lower=0.1,upper=INDEF,radius=0) */
        printf("Growing mask and checking neighboring pixels\n");
        double* gfirstsel = (double*)malloc(sx*sy*sizeof(double));
        convolve(firstsel, sx, sy, gfirstsel, growth_kernel, 3);
        for (i=0; i<sx*sy; i++) {
            gfirstsel[i] = sigmap[i] > sigclip && gfirstsel[i] > 0.5 ? 1. : 0.;
        }
    
        double sigcliplow = sigfrac * sigclip;
    
        /* # grow CRs by one pixel and lower detection limit */

        /* sigcliplow = sigfrac * sigclip */

        /* if (verbose) { */
        /*  print("Finding neighbouring pixels affected by cosmic rays") */
        /*  print("  sigma limit = "//sigcliplow) */
        /*  print("") */
        /*  } */

        /* convolve(gfirstsel,finalsel,gkernel) */
        /* imreplace(finalsel,1,lower=0.5,upper=INDEF,radius=0) */
        /* imarith(finalsel,"*",sigmap,finalsel) */
        /* imreplace(finalsel,0,upper=sigcliplow,lower=INDEF,radius=0) */
        /* imreplace(finalsel,1,lower=0.1,upper=INDEF,radius=0) */

        printf("Growing mask again and checking for weaker neighboring pixels\n");
        double* finalsel = (double*)malloc(sx*sy*sizeof(double));
        convolve(gfirstsel, sx, sy, finalsel, growth_kernel, 3);
        for (i=0; i<sx*sy; i++) {
            finalsel[i] = sigmap[i] > sigcliplow && finalsel[i] > 0.5 ? 1. : 0.;
        }

        /* # determine number of CRs found in this iteration */
        /* imdel(gfirstsel) */
        /* imcalc(finalsel//","//outmask,gfirstsel,"(1-im2)*im1",verb-) */
        /* imstat(gfirstsel,fields="npix",lower=0.5,upper=INDEF,for-) | scan(npix) */

        imdel(gfirstsel);
        int crpix_found = 0;
        for (i=0; i<sx*sy; i++) {
            crpix_found += finalsel[i];
        }
        printf("Found a total of %d cosmic-ray affected pixels\n", crpix_found);
    
        /* # create cleaned output image; use 3x3 median with CRs excluded */
        /* if (verbose) { */
        /*  print("Creating output:") */
        /*  print("  bad pixel mask: "//outmask) */
        /*  print("  cleaned image: "//output) */
        /*  print("") */
        /*  } */
        /* imdel(med5) */
        /* imarith(outmask,"+",finalsel,outmask) */
        /* imreplace(outmask,1,lower=1,upper=INDEF,radius=0) */
        /* imcalc(outmask,inputmask,"(1.-10000.*im1)",verb-) */
        /* imarith(oldoutput,"*",inputmask,inputmask) */
        /* median(inputmask,med5,5,5,zloreject=-9999,zhi=INDEF,verb-) */
        /* imarith(outmask,"*",med5,med5) */
        /* if (i>1) imdel(output) */
        /* imcalc(oldoutput//","//outmask//","//med5,output,"(1.-im2)*im1+im3",verb-) */

        printf("create cleaned output image\n");
        
        double* data_filtered = median_filtered(data, sx, sy, 5, 5, finalsel);

        //
        // Replace all CR-affected pixels with the median of the surrounding
        //

        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                
                // If this is a cosmic pixel
                if (finalsel[i] > 0.5) {
                    data[i] = data_filtered[i];
                    
                    // Now mark all pixels in a 7 pixel box to be affected by the CR
                    for ( x = ( (_x-3) <  0 ?  0 : (_x-3) );
                          x < ( (_x+4) > sx ? sx : (_x+4) );
                          x++) {
                        for ( y = ( (_y-3) <  0 ?  0 : (_y-3) );
                              y < ( (_y+4) > sy ? sy : (_y+4) );
                              y++) {

                            pixel_changed[y + x*sy] = 1;
                        }
                    }
                    
                } // end if finalsel > 0.5
                
            }
        }
        
                    
        /* for (i=0; i<sx*sy; i++) { */
        /*     data[i] = finalsel[i] > 0.5 ? data_filtered[i] : data[i]; */
        /* } */
   
        for(_x=0; _x<sx; _x++){
            for(_y=0; _y<sy; _y++) {
                //output[_y + _x*sy] = pixel_changed[_y + _x*sy];
                output[_y + _x*sy] = data[_y + _x*sy];
            }
        }
        
    }
    
    
    printf("done!\n");
    return;
}
    

#ifdef __STANDALONE__

void main()
{

    printf("test\n");

    int sx=512, sy=512;

    double* data = (double*)malloc(sx*sy*sizeof(double));
    double* retval = (double*)malloc(sx*sy*sizeof(double));

    
    return;
    
}

#endif
