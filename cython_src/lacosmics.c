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

void lacosmics__cy(double* data,
                   double* out_cleaned, int* out_mask, int* out_saturated,
                   int sx, int sy,
                   double gain, double readnoise,
                   double sigclip, double sigfrac, double objlim,
                   double saturation_limit, int verbose,
                   int niter)
{
    
    printf("test\n");
    printf("Gain=%f\n",gain);
    printf("readnoise=%f\n",readnoise);

    int lx, ly, _x, _y, i, wx, wy, n, x, y;
    int kernel_center = 1, ksize=3;
    int dx, dy, kx, ky;
    double tmpd;
    int sx2 = sx*2, sy2 = sy*2, ssm, ssp;

    // memory demand calculation
    // assume sx*sy = 16M
    double* larger_2x2 = (double*)malloc(sx*2*sy*2*sizeof(double));      // 64
    double* lapla_convolved = (double*)malloc(sx*2*sy*2*sizeof(double)); // 64
    double* deriv2 = (double*)malloc(sx*sy*sizeof(double));              // 16
    double* data_med5 = (double*)malloc(sx*sy*sizeof(double));                // 16
    double* noise = (double*)malloc(sx*sy*sizeof(double));               // 16
    double* sigmap = (double*)malloc(sx*sy*sizeof(double));              // 16
    double* sigmap_med5 = (double*)malloc(sx*sy*sizeof(double));         // 16
    double* firstsel = malloc(sx*sy*sizeof(double));                     // 16
    double* data_med3 = malloc(sx*sy*sizeof(double));                     // 16
    double data_med7;
    double* gfirstsel = (double*)malloc(sx*sy*sizeof(double));           // 16
    double* finalsel = (double*)malloc(sx*sy*sizeof(double));            // 16
    double* data_filtered = (double*)malloc(sx*sy*sizeof(double));       // 16
    //                                                               total: 272 Mpixel * 8 bytes ~ 2.2 GB
    
    int* blkavg_pixelcount = (int*)malloc(sx*sy*sizeof(int));            // 16
    int* pixel_changed = (int*)malloc(sx*sy*sizeof(int));                // 16
    int* crj_iteration = (int*)malloc(sx*sy*sizeof(int));                // 16
    int* saturated = (int*)malloc(sx*sy*sizeof(int));                    // 16
    //                                                               total: 32 Mpixel * 4 bytes ~ 128 MB
    
    double* neighbors = (double*)malloc(MAXMEDIAN*MAXMEDIAN*sizeof(double));
    
    double laplace_kernel[3*3] = {  
          0., -1.,  0. ,
         -1. , 4., -1. ,
          0., -1.,  0. 
    };
    
    double growth_kernel[3*3] = {
         1., 1., 1. ,
         1., 1., 1. ,
         1., 1., 1. ,
    };
    

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
    if (verbose) printf("Working on frame with dimensions %d x %d\n", sx, sy);

    int iteration = 0;

    for (i=0; i<sx*sy; i++) {
        crj_iteration[i] = 0;
        saturated[i] = 0;
    }
    
    for (iteration = 0; iteration < niter; iteration++) {
        if (verbose) printf("\nStarting iteration %d (of %d)...\n\n", iteration+1, niter);
        
        // duplicate all pixels 2x2
        if (verbose) printf("Computing larger 2x2 blkrep array\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i]==1 || iteration == 0) {
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
        if (verbose) printf("Convolving with 3x3 laplace kernel!\n");
        /* convolve(larger_2x2, 2*sx, 2*sy, lapla_convolved,  laplace_kernel, 3); */
        for (dx=0; dx<sx2; dx++) {
            for (dy=0; dy<sy2; dy++) {
                i = dy + dx*sy2;
                if (pixel_changed[dy/2 + dx/2*sy]==1 || iteration == 0) {
                    
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

        if (verbose) printf("Running blkavg\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i]==1 || iteration == 0) {
                
                    deriv2[i] = 0.25 * (
                        lapla_convolved[ _y*2   + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2   + (_x*2+1)*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2+1)*sy*2 ]
                    );
                }
            }
        }
    
        /* for (i=0; i<sx*sy; i++) { */
        /*     output[i] = deriv2[i]; */
        /* } */

        
        if (verbose) printf("Median-filtering the data, computing noise and significance\n");
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (_x=0; _x<sx; _x++) {
            for (_y=0; _y<sy; _y++) {
                i = _y + _x*sy;
                
                // only do work if necessary
                if (pixel_changed[i]==1 || iteration == 0) {
                    n = 0;
                    for ( x = ( (_x-2) <  0 ?  0 : (_x-2) );
                          x < ( (_x+3) > sx ? sx : (_x+3) );
                          x++) {
                        for ( y = ( (_y-2) <  0 ?  0 : (_y-2) );
                              y < ( (_y+3) > sy ? sy : (_y+3) );
                              y++) {
                            neighbors[n++] = data[y + x*sy];
                        }
                    }
                    // Now compute the median
                    gsl_sort(neighbors, 1, n);
                    tmpd = gsl_stats_median_from_sorted_data(neighbors, 1, n);
                    data_med5[i] = tmpd < 0 ? 0.0 : tmpd;

                    /* // During the first iteration, create a mask for saturated stars */
                    /* if (iteration == 0 && saturation_limit > 0) { */
                    /*     if (med5[i] > saturation_limit) { */
                    /*         saturated[i] = 1; */
                    /*     } */
                    /* } */
                            
                    // Compute noise estimate
                    noise[i] = sqrt(data_med5[i]*gain + readnoise*readnoise);

                    // Compute significance of pixel
                    sigmap[i] = (deriv2[i] / noise[i]) / 2.;
                }
            }
        }

        //
        // If masking saturated pixels was requested, grow the masked area by +/- 2 pixels
        //
        if (iteration == 0 && saturation_limit > 0) {
            if (verbose) printf("Creating saturated pixel mask\n");
            ssm = -3; ssp = 4;
            for (_x=0; _x<sx; _x++) {
                for (_y=0; _y<sy; _y++) {
                    if (data_med5[_y + _x*sy] > saturation_limit) {
                        
                        for ( x = ( (_x+ssm) <  0 ?  0 : (_x+ssm) );
                              x < ( (_x+ssp) > sx ? sx : (_x+ssp) );
                              x++) {
                            for ( y = ( (_y+ssm) <  0 ?  0 : (_y+ssm) );
                                  y < ( (_y+ssp) > sy ? sy : (_y+ssp) );
                                  y++) {
                                saturated[y + x*sy] = 1;
                            }
                        }
                    }
                }
            }
        }
        
        if (verbose) printf("removing large structure\n");
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (_x=0; _x<sx; _x++) {
            for (_y=0; _y<sy; _y++) {
                i = _y + _x*sy;
                
                // only do work if necessary
                if (pixel_changed[i]==1 || iteration == 0) {
                    n = 0;
                    for ( x = ( (_x-2) <  0 ?  0 : (_x-2) );
                          x < ( (_x+3) > sx ? sx : (_x+3) );
                          x++) {
                        for ( y = ( (_y-2) <  0 ?  0 : (_y-2) );
                              y < ( (_y+3) > sy ? sy : (_y+3) );
                              y++) {

                            neighbors[n++] = sigmap[y + x*sy];
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
        

                
        if (verbose) printf("Selecting candidate CRs\n");
        for (i=0; i<sx*sy; i++) {
            firstsel[i] = ((sigmap[i] > sigclip) && (saturated[i] == 0)) ? 1 : 0;
        }
    
        if (verbose) printf("subtract background and smooth component of objects\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;

                if ((pixel_changed[i] || iteration == 0)) {
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
        
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;

                if ((pixel_changed[i] || iteration == 0) && firstsel[i] > 0) {
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
                            neighbors[n++] = data_med3[y + x*sy];
                        }
                    }
                    gsl_sort(neighbors, 1, n);
                    data_med7 = gsl_stats_median_from_sorted_data(neighbors, 1, n);
                } // end if pixel_changed


                tmpd = (data_med3[i] - data_med7) / noise[i];
                tmpd = tmpd < 0.01 ? 0.01 : tmpd;

                if (firstsel[i] > 0) {
                    firstsel[i] = sigmap[i] > (tmpd * objlim) ? 1 : 0;
                }

                // Also reset the mask of CR pixels to 0
                pixel_changed[i] = 0;
               
                
            }
        }


        if (verbose) printf("Growing mask and checking neighboring pixels\n");
        convolve(firstsel, sx, sy, gfirstsel, growth_kernel, 3);
        for (i=0; i<sx*sy; i++) {
            gfirstsel[i] = sigmap[i] > sigclip && gfirstsel[i] > 0.5 ? 1. : 0.;
        }
    
        double sigcliplow = sigfrac * sigclip;
    
        if (verbose) printf("Growing mask again and checking for weaker neighboring pixels\n");
        convolve(gfirstsel, sx, sy, finalsel, growth_kernel, 3);
        for (i=0; i<sx*sy; i++) {
            finalsel[i] = sigmap[i] > sigcliplow && finalsel[i] > 0.5 ? 1. : 0.;
        }

        int crpix_found = 0;
        for (i=0; i<sx*sy; i++) {
            crpix_found += finalsel[i];
        }
        if (verbose) printf("Found a total of %d cosmic-ray affected pixels\n", crpix_found);
    
        if (verbose) printf("create cleaned output image\n");
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (_x=0; _x<sx; _x++) {
            for (_y=0; _y<sy; _y++) {
                i = _y + _x*sy;

                // only compute the median of neighbors if we need to replace this pixel
                if (finalsel[i] > 0) {
                
                    // Collect all pixels in the neighborhood of this pixel
                    n = 0;
                    for ( x = ( (_x-1) <  0 ?  0 : (_x-1) );
                          x < ( (_x+2) > sx ? sx : (_x+2) );
                          x++) {
                        for ( y = ( (_y-1) <  0 ?  0 : (_y-2) );
                              y < ( (_y+2) > sy ? sy : (_y+3) );
                              y++) {
                            // Filter out pixels labeled as cosmics
                            if (finalsel[y + x*sy] == 0) {
                                neighbors[n++] = data[y + x*sy];
                            }
                        }
                    }
                    // Now compute the median
                    gsl_sort(neighbors, 1, n);
                    tmpd = gsl_stats_median_from_sorted_data(neighbors, 1, n);

                    // Replace this cosmic affected pixel with the median of its neighbors
                    data_filtered[i] = tmpd;

                    crj_iteration[i] = iteration + 1;

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
                    
                } else {
                    data_filtered[i] = data[i];
                }
                
            }
        }


        // If necessary, prepare for the next iteration
        if (iteration < niter-1) {
            for (i=0; i<sx*sy; i++) {
                data[i] = data_filtered[i];
            }
        }

        if (verbose) printf("Done with iteration %d (of %d)...\n\n", iteration+1, niter);

    }

    // Copy cleaned image to output buffer
    if (verbose) printf("Preparing output...\n");
    for(_x=0; _x<sx; _x++){
        for(_y=0; _y<sy; _y++) {
            out_cleaned[_y + _x*sy] = data_filtered[_y + _x*sy];
            out_mask[_y + _x*sy] = crj_iteration[_y + _x*sy];
            out_saturated[_y + _x*sy] = saturated[_y + _x*sy];
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
