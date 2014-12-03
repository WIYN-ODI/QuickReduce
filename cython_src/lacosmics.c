/**
 *
 * (c) Ralf Kotulla, kotulla@uwm.edu
 *
 * This module implements a iterative sigma-clipping routine to speed
 * up the corresponding functionality in podi_imcombine.
 *
 */

//#define _HEAPSORT_

#define true 0
#define false 1
typedef int bool;

#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <stdlib.h>
#include <time.h>

#define nan 0./0.
#define MAXMEDIAN 7
#define HEAPSIZE 2500

#define True 1
#define False 0

// #define printf //

#define tracepx //printf


void convolve(double* input, int sx, int sy,
              double* output,
              double* kernel, int ksize)
{
    int kernel_center = (ksize-1)/2;

    int dx, dy, kx, ky, i;

    // set all pixels to 0
    for (i=0; i<sx*sy; i++) {
        output[i] = 0;
    }
    
    // Execute the convolution
    for (dx=0; dx<sx; dx++) {
        for (dy=0; dy<sy; dy++) {


            // 
            // result for this pixel convolved with the kernel
            //
            /* for (kx = -1*kernel_center; kx<= kernel_center; kx++) { */
            /*     for (ky = -1*kernel_center; ky<= kernel_center; ky++) { */
            /*         if (dx+kx < 0 || dx+kx >= sx || dy+ky < 0 || dy+ky > sy) { */
            /*             continue; */
            /*         } */

            /*         output[dy + dx*sy] += input[dy+ky + (dx+kx)*sy] */
            /*             * kernel[ky+kernel_center + (kx+kernel_center)*ksize]; */
                    
            /*     } */
            /* } */

            //
            // Compute the contribtion to all neighboring pixels from this pixel
            //
            if (input[dy + dx*sy] != 0) {
                for (kx = -1*kernel_center; kx<= kernel_center; kx++) {
                    for (ky = -1*kernel_center; ky<= kernel_center; ky++) {

                        if (dx+kx < 0 || dx+kx >= sx || dy+ky < 0 || dy+ky > sy) {
                            continue;
                        }
                        output[dy+ky + (dx+kx)*sy] += input[dy + dx*sy]
                            * kernel[kernel_center-ky + (kernel_center-kx)*ksize];
                    }
                }
            }

        }
    }

    return;
}


/*
 * This routine was found at
 * http://stackoverflow.com/questions/810657/fastest-code-c-c-to-select-the-median-in-a-set-of-27-floating-point-values
 * written by user chmike
 * 
 * Note that the routine only gives correct answers if the number of input pixels is odd
 */ 
double heapMedian3( double *a, int n_full )
{
   double left[HEAPSIZE], right[HEAPSIZE], median, *p;
   int nLeft, nRight;
   int n_half = (int) ((float)n_full / 2.0 + 0.5);
   int nVal;
   //printf("nhalf=%d, n_full=%d\n", n_half, n_full);
   
   // pick first value as median candidate
   p = a;
   median = *p++;
   nLeft = nRight = 1;

   for(;;) {
       // get next value
       double val = *p++;

       // if value is smaller than median, append to left heap
       if( val < median ) {
           // move biggest value to the heap top
           int child = nLeft++, parent = (child - 1) / 2;
           while( parent && val > left[parent] ) {
               left[child] = left[parent];
               child = parent;
               parent = (parent - 1) / 2;
           }
           left[child] = val;

           // if left heap is full
           if( nLeft == n_half ) {
               // for each remaining value
               for( nVal = n_full - (p - a); nVal; --nVal ) {
                   // get next value
                   val = *p++;

                   // if value is to be inserted in the left heap
                   if( val < median ) {
                       child = left[2] > left[1] ? 2 : 1;
                       if( val >= left[child] )
                           median = val;
                       else {
                           median = left[child];
                           parent = child;
                           child = parent*2 + 1;
                           while( child < n_half ) {
                               if( child < n_half-1 && left[child+1] > left[child] )
                                   ++child;
                               if( val >= left[child] )
                                   break;
                               left[parent] = left[child];
                               parent = child;
                               child = parent*2 + 1;
                           }
                           left[parent] = val;
                       }
                   }
               }
               return median;
           }
       }

       // else append to right heap
       else  {
           // move smallest value to the heap top
           int child = nRight++, parent = (child - 1) / 2;
           while( parent && val < right[parent] ) {
               right[child] = right[parent];
               child = parent;
               parent = (parent - 1) / 2;
           }
           right[child] = val;

           // if right heap is full
           if( nRight == n_half ) {
               // for each remaining value
               for( nVal = n_full - (p - a); nVal; --nVal ) {
                   // get next value
                   val = *p++;

                   // if value is to be inserted in the right heap
                   if( val > median ) {
                       child = right[2] < right[1] ? 2 : 1;
                       if( val <= right[child] )
                           median = val;
                       else {
                           median = right[child];
                           parent = child;
                           child = parent*2 + 1;
                           while( child < n_half ) {
                               if( child < n_half-1 && right[child+1] < right[child] )
                                   ++child;
                               if( val <= right[child] )
                                   break;
                               right[parent] = right[child];
                               parent = child;
                               child = parent*2 + 1;
                           }
                           right[parent] = val;
                       }
                   }
               }
               return median;
           }
       }
   }
}


double gsl_find_median( double *neighbors, int n )
{
    gsl_sort(neighbors, 1, n);
    return gsl_stats_median_from_sorted_data(neighbors, 1, n);
}

double find_median(double *neighbors, int n) 
{
    if (n%2 == 1) {
        // Even number, use the faster heapMedian
        return heapMedian3(neighbors, n);
    }
    return gsl_find_median(neighbors, n);
}


int dumpbuffertofile(double* buf, int sx, int sy, char* filename)
{
    FILE* fp = fopen(filename, "w");
    int x,y;
    
    for (y=0; y<sy; y++) {
        for (x=0; x<sx; x++) {
            fprintf(fp, "% 4d %4d %10.4f\n", x, y, buf[y + x*sy]);
        }
    }
    fclose(fp);
    return sx*sy;
}
int dumpbuffertofile_int(int* buf, int sx, int sy, char* filename)
{
    FILE* fp = fopen(filename, "w");
    int x,y;
    
    for (y=0; y<sy; y++) {
        for (x=0; x<sx; x++) {
            fprintf(fp, "% 4d %4d % 8d\n", x, y, buf[y + x*sy]);
        }
    }
    fclose(fp);
    return sx*sy;
}
    
    
void lacosmics__cy(double* data,
                   double* out_cleaned, int* out_mask, int* out_saturated,
                   int sx, int sy,
                   double gain, double readnoise,
                   double sigclip, double sigfrac, double objlim,
                   double saturation_limit, int verbose,
                   int niter)
{
    
    if (verbose) {
        printf("Gain=%f\n",gain);
        printf("readnoise=%f\n",readnoise);
    }
    if (gain <= 0) gain = 1.;
    
    int tracepixel = 233+484*sy; //484 + 233*sy;

    
    int lx, ly, _x, _y, i, j, wx, wy, n, x, y, ix, iy;
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
    double* sigmap_prime = (double*)malloc(sx*sy*sizeof(double));        // 16
    double* firstsel = malloc(sx*sy*sizeof(double));                     // 16
    double* data_med3 = malloc(sx*sy*sizeof(double));                     // 16
//    double* data_med7 = malloc(sx*sy*sizeof(double));                     // 16
    double* data_med7 /*RK*/= malloc(sx*sy*sizeof(double));
    double* data_med7_tmpd /*RK*/= malloc(sx*sy*sizeof(double));
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
    char filename[100];
    
    int rerun_entirely = False;
    
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
    int crpix_found = 1;

    for (i=0; i<sx*sy; i++) {
        crj_iteration[i] = 0;
        saturated[i] = 0;
    }
    
    for (iteration = 0; iteration < niter && crpix_found > 0; iteration++) {
        if (verbose) printf("\nStarting iteration %d (of %d)...\n\n", iteration+1, niter);
        
        if (verbose) {
            sprintf(filename, "data_input_%d.cat", iteration);
            dumpbuffertofile(data, sx, sy, filename);
            sprintf(filename, "pixelchanged_input_%d.cat", iteration);
            dumpbuffertofile_int(pixel_changed, sx, sy, filename);
        }

        // duplicate all pixels 2x2
        tracepx("### Trace pixel: input data = %f\n", data[tracepixel]);
        if (verbose) printf("Computing larger 2x2 blkrep array\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i]==1 || iteration == 0 || rerun_entirely) {
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
                if (pixel_changed[dy/2 + dx/2*sy]==1 || iteration == 0 || rerun_entirely) {
                    
                    lapla_convolved[i] = 0.0;

                    //
                    // Compute the result for this pixel convolved with the kernel
                    //
                    for (kx = -1*kernel_center; kx<= kernel_center; kx++) {
                        for (ky = -1*kernel_center; ky<= kernel_center; ky++) {
                            /* if (!(dx+kx < 0 || dx+kx >= sx2 || dy+ky < 0 || dy+ky > sy2)) { */

                            /*     lapla_convolved[i] += larger_2x2[dy+ky + (dx+kx)*sy2] */
                            /*         * laplace_kernel[ky+kernel_center + (kx+kernel_center)*ksize]; */
                            /* } */
                            ix = ( dx+kx + sx2 ) % sx2;
                            iy = ( dy+ky + sy2 ) % sy2;
                            
                            lapla_convolved[i] += larger_2x2[iy + ix*sy2]
                                * laplace_kernel[ky+kernel_center + (kx+kernel_center)*ksize];
                        }
                    }
                    /* printf("replacing negative pixel values\n"); */
                    /* imreplace(lapla_convolved, sx*2, sy*2, -1e99, 0., 0.0); */
                    lapla_convolved[i] = lapla_convolved[i] < 0 ? 0. : lapla_convolved[i];
                }
                
            }
        }
        if (verbose) {
            sprintf(filename, "laplace_%d.cat", iteration);
            dumpbuffertofile(lapla_convolved, sx2, sy2, filename);
        }
        
            
        if (verbose) printf("Running blkavg\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;
                if (pixel_changed[i]==1 || iteration == 0 || rerun_entirely) {
                
                    deriv2[i] = 0.25 * (
                        lapla_convolved[ _y*2   + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2  )*sy*2 ] +
                        lapla_convolved[ _y*2   + (_x*2+1)*sy*2 ] +
                        lapla_convolved[ _y*2+1 + (_x*2+1)*sy*2 ]
                    );
                }
            }
        }
        if (verbose) {
            sprintf(filename, "deriv2_%d.cat", iteration);
            dumpbuffertofile(deriv2, sx, sy, filename);
        }
        tracepx("### Trace pixel: deriv2 = %f\n", deriv2[tracepixel]);
        /* for (i=0; i<sx*sy; i++) { */
        /*     output[i] = deriv2[i]; */
        /* } */

        
        if (verbose) printf("Median-filtering the data, computing noise and significance\n");
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (_x=0; _x<sx; _x++) {
            for (_y=0; _y<sy; _y++) {
                i = _y + _x*sy;
                
                // only do work if necessary
                if (pixel_changed[i]==1 || iteration == 0 || rerun_entirely) {
                    n = 0;
                    
                    for ( x = (_x-2); x < (_x+3); x++) {
                        for ( y = (_y-2); y < (_y+3); y++) {
                            ix = ( x+sx ) % sx;
                            iy = ( y+sy ) % sy;

                            neighbors[n++] = data[iy + ix*sy];
                        }
                    }

                    tmpd = find_median(neighbors, n);
                    data_med5[i] = tmpd < 0.00001 ? 0.00001 : tmpd;

                    /* // During the first iteration, create a mask for saturated stars */
                    /* if (iteration == 0 && saturation_limit > 0) { */
                    /*     if (med5[i] > saturation_limit) { */
                    /*         saturated[i] = 1; */
                    /*     } */
                    /* } */
                            
                    // Compute noise estimate
                    noise[i] = sqrt(data_med5[i]*gain + readnoise*readnoise)/gain;
                    // Compute significance of pixel
                    sigmap[i] = (deriv2[i] / noise[i]) / 2.;
                    
                }
            }
        }
        tracepx("### Trace pixel: data_med5 = %f\n", data_med5[tracepixel]);
        tracepx("### Trace pixel: sigmap = %f\n", sigmap[tracepixel]);
        tracepx("### Trace pixel: noise = %f\n", noise[tracepixel]);
        if (verbose) {
            sprintf(filename, "data_med5_%d.cat", iteration);
            dumpbuffertofile(data_med5, sx, sy, filename);
            sprintf(filename, "noise_%d.cat", iteration);
            dumpbuffertofile(noise, sx, sy, filename);
            sprintf(filename, "sigmap_%d.cat", iteration);
            dumpbuffertofile(sigmap, sx, sy, filename);
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
        if (verbose) {
            sprintf(filename, "saturated_%d.cat", iteration);
            dumpbuffertofile(saturated, sx, sy, filename);
        }
        
        if (verbose) printf("removing large structure\n");
        wx = wy = 5; dx = (wx-1)/2; dy = (wy-1)/2;
        for (_x=0; _x<sx; _x++) {
            for (_y=0; _y<sy; _y++) {
                i = _y + _x*sy;
                
                // only do work if necessary
                if (pixel_changed[i]==1 || iteration == 0 || rerun_entirely) {
                    n = 0;
                    for ( x = (_x-2); x < (_x+3); x++) {
                        for ( y = (_y-2); y < (_y+3); y++) {
                            ix = ( x+sx ) % sx;
                            iy = ( y+sy ) % sy;

                            neighbors[n++] = sigmap[iy + ix*sy];
                        }
                    }
                    // Now compute the median
                    /* gsl_sort(neighbors, 1, n); */
                    /* sigmap_med5[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n); */
                    sigmap_med5[i] = find_median(neighbors, n);
                    // Subtract the smoothed significance map from the pixel significance map
                    sigmap_prime[i] = sigmap[i] - sigmap_med5[i];
                }
            }
        }
        tracepx("### Trace pixel: sigmap_med5 = %f\n", sigmap_med5[tracepixel]);
        tracepx("### Trace pixel: sigmap_prime = %f\n", sigmap_prime[tracepixel]);
        if (verbose) {
            sprintf(filename, "sigmap_med5_%d.cat", iteration);
            dumpbuffertofile(sigmap_med5, sx, sy, filename);
            sprintf(filename, "sigmap_prime_%d.cat", iteration);
            dumpbuffertofile(sigmap_prime, sx, sy, filename);
        }
        
                
        if (verbose) printf("Selecting candidate CRs\n");
        for (i=0; i<sx*sy; i++) {
            firstsel[i] = ((sigmap_prime[i] > sigclip) && (saturated[i] == 0)) ? 1 : 0;
        }
        if (verbose) {
            sprintf(filename, "firstsel_x1_%d.cat", iteration);
            dumpbuffertofile(firstsel, sx, sy, filename);
        }
   
        if (verbose) printf("subtract background and smooth component of objects\n");
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;

                if ((pixel_changed[i]==1 || iteration == 0 || rerun_entirely)) { // RK
                    //
                    // do 3x3 median filtering
                    //
                    n=0;
                    for ( x = (_x-1); x < (_x+2); x++) {
                        for ( y = (_y-1); y < (_y+2); y++) {
                            ix = ( x+sx ) % sx;
                            iy = ( y+sy ) % sy;
                            neighbors[n++] = data[iy + ix*sy];
                        }
                    }
                    /* gsl_sort(neighbors, 1, n); */
                    /* data_med3[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n); */
                    data_med3[i] = find_median(neighbors, n);
                    
                } // end if pixel_changed
            }
        }
        tracepx("### Trace pixel: data_med3 = %f\n", data_med3[tracepixel]);
        if (verbose) {
            sprintf(filename, "data_med3_%d.cat", iteration);
            dumpbuffertofile(data_med3, sx, sy, filename);
        }
        
        /* for (i=0; i<sx*sy; i++) { */
        /*     out_cleaned[i] = data_med3[i]; */
        /* } */
        
        
        if (verbose) {
            sprintf(filename, "firstsel_xxxx_%d.cat", iteration);
            dumpbuffertofile(firstsel, sx, sy, filename);
        }
        for (_x = 0; _x<sx; _x++) {
            for (_y = 0; _y < sy; _y++) {
                i = _y + _x*sy;

                if ((pixel_changed[i]==1 || iteration == 0 || rerun_entirely) && firstsel[i] > 0) { // RK
                    //
                    // do 7x7 median filtering
                    //
                    n=0;
                    for ( x = (_x-3); x < (_x+4); x++) {
                        for ( y = (_y-3); y < (_y+4); y++) {
                            ix = ( x+sx ) % sx;
                            iy = ( y+sy ) % sy;
                            if (x<0 || x>=sx || y<0 || y >= sy) continue;
                            // without this we have periodic boundary conditions
                            
                            neighbors[n++] = data_med3[y + x*sy];
                            //neighbors[n++] = data_med3[iy + ix*sy];
                        }
                    }
                    /* gsl_sort(neighbors, 1, n); */
                    /* data_med7[i] = gsl_stats_median_from_sorted_data(neighbors, 1, n); */
                    data_med7[i] = find_median(neighbors, n);

                    data_med7_tmpd[i] = (data_med3[i] - data_med7[i]) / noise[i];
                    data_med7_tmpd[i] = data_med7_tmpd[i] < 0.01 ? 0.01 : data_med7_tmpd[i];
                    
                    // out_cleaned[i] = tmpd; // this is f in the python version
                    
                    // if (firstsel[i] > 0) {
                    if (verbose) {
                        //printf("%4d/%4d --> firstsel=%.0f &&  sigmap_prime=%10.4f > data_med7_tmpd=%10.4f * objlim=%4.1f",
                        printf("%4d/%4d --> %.0f &&  %10.4f > %10.4f * %4.1f",
                               _x, _y, firstsel[i], sigmap_prime[i], data_med7_tmpd[i], objlim);
                    }
                    
                    firstsel[i] = firstsel[i] > 0 && sigmap_prime[i] > (data_med7_tmpd[i] * objlim) ? 1 : 0;

                    if (verbose) {
                        printf("  ===> %.0f\n", firstsel[i]);
                    }
                    
                } else {
                    // Mask as not a pixel
                    firstsel[i] = 0;
                }
                
                   

                // Also reset the mask of CR pixels to 0
                pixel_changed[i] = 0;
               
                
            }
        }
        if (verbose) {
            sprintf(filename, "firstsel_xxxy_%d.cat", iteration);
            dumpbuffertofile(firstsel, sx, sy, filename);
        }
        tracepx("### Trace pixel: firstsel = %f\n", firstsel[tracepixel]);
        if (verbose) {
            sprintf(filename, "firstsel_x2_%d.cat", iteration);
            dumpbuffertofile(firstsel, sx, sy, filename);
            sprintf(filename, "data_med7_%d.cat", iteration);
            dumpbuffertofile(data_med7, sx, sy, filename);
            sprintf(filename, "data_med7_tmpd_%d.cat", iteration);
            dumpbuffertofile(data_med7_tmpd, sx, sy, filename);
        }

        if (verbose) printf("Growing mask and checking neighboring pixels\n");
        /* n=0; for(i=0; i<sx*sy; i++) if (firstsel[i] > 0.5) n++; printf("Initial #CRs: %d (>%f sigma)\n", n, sigclip); */
        convolve(firstsel, sx, sy, gfirstsel, growth_kernel, 3);
        /* n=0; for(i=0; i<sx*sy; i++) if (gfirstsel[i] > 0.5) n++; printf("Initial grown #CRs: %d\n", n); */
        for (i=0; i<sx*sy; i++) {
            gfirstsel[i] = sigmap_prime[i] > sigclip && gfirstsel[i] > 0.5 && saturated[i] == 0 ? 1. : 0.;
        }
        /* n=0; for(i=0; i<sx*sy; i++) if (gfirstsel[i] > 0.5) n++; printf("remaining grown #CRs: %d\n", n); */
        tracepx("### Trace pixel: gfirstsel = %f\n", gfirstsel[tracepixel]);
        if (verbose) {
            sprintf(filename, "gfirstsel_%d.cat", iteration);
            dumpbuffertofile(gfirstsel, sx, sy, filename);
        }
       
    
        double sigcliplow = sigfrac * sigclip;
    
        if (verbose) printf("Growing mask again and checking for weaker neighboring pixels\n");
        convolve(gfirstsel, sx, sy, finalsel, growth_kernel, 3);
        /* n=0; for(i=0; i<sx*sy; i++) if (finalsel[i] > 0.5) n++; printf("2nd grown #CRs: %d\n", n); */
        for (i=0; i<sx*sy; i++) {
            finalsel[i] = sigmap_prime[i] > sigcliplow && finalsel[i] > 0.5 && saturated[i] == 0 ? 1. : 0.;
        }
        /* n=0; for(i=0; i<sx*sy; i++) if (finalsel[i] > 0.5) n++; printf("final #CRs: %d (>%f sigma)\n", n, sigcliplow); */
        tracepx("### Trace pixel: finalsel = %f\n", finalsel[tracepixel]);
        if (verbose) {
            sprintf(filename, "finalsel_%d.cat", iteration);
            dumpbuffertofile(finalsel, sx, sy, filename);
        }

        /* for (i=0; i<sx*sy; i++) { */
        /*     out_cleaned[i] = finalsel[i]; */
        /* } */

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
                
                    crj_iteration[i] = iteration + 1;

                    // Collect all pixels in the neighborhood of this pixel
                    n = 0;
                    for ( x = (_x-2); x < (_x+3); x++) {
                        for ( y = (_y-2); y < (_y+3); y++) {
                            ix = ( x+sx ) % sx;
                            iy = ( y+sy ) % sy;
                            j = iy + ix*sy;
                            
                            // Filter out pixels labeled as cosmics
                            // Ignore all pixels masked as cosmic rays in this
                            // or any of the past iterations
                            if (crj_iteration[j] == 0 && finalsel[j] == 0) {
                                neighbors[n++] = data[j];
                            }
                        }
                    }
                    // Now compute the median
                    /* gsl_sort(neighbors, 1, n); */
                    /* tmpd = gsl_stats_median_from_sorted_data(neighbors, 1, n); */
                    if (n>1) {
                        tmpd = find_median(neighbors, n);
                    } else if (n==1) {
                        if (verbose) printf("found only a single pixel (x/y=%d,%d)!\n", _x, _y);
                        tmpd = neighbors[0];
                    } else {
                        if (verbose) printf("No valid pixels found nearby (x/y=%d,%d)!\n", _x, _y);
                        tmpd = nan;
                    }
                    
                    
                    // Replace this cosmic affected pixel with the median of its neighbors
                    data_filtered[i] = tmpd;

                    // Now mark all pixels in a 7 pixel box to be affected by the CR
                    for ( x = ( (_x-3) <  0 ?  0 : (_x-3) );
                          x < ( (_x+4) > sx ? sx : (_x+4) );
                          x++) {
                        for ( y = ( (_y-3) <  0 ?  0 : (_y-3) );
                              y < ( (_y+4) > sy ? sy : (_y+4) );
                              y++) {

                            pixel_changed[y + x*sy] = 1;
                            //printf("Marking pixel as changed (%d): %d %d\n", iteration, x, y);
                            
                        }
                    }
                    
                } else {
                    data_filtered[i] = data[i];
                }
                
            }
        }
        if (verbose) printf("Done cleaning!\n");
        if (verbose) {
            sprintf(filename, "data_filtered_%d.cat", iteration);
            dumpbuffertofile(data_filtered, sx, sy, filename);
            sprintf(filename, "pixel_changed_%d.cat", iteration);
            dumpbuffertofile_int(pixel_changed, sx, sy, filename);
        }
        
        tracepx("### Trace pixel: data_filtered = %f\n", data_filtered[tracepixel]);


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


    // Once we are done, free all memory allocated
    free(larger_2x2);
    free(lapla_convolved);
    free(deriv2);
    free(data_med5);
    free(noise);
    free(sigmap);
    free(sigmap_med5);
    free(sigmap_prime);
    free(firstsel);
    free(data_med3);
    free(gfirstsel);
    free(finalsel);
    free(data_filtered);
    free(blkavg_pixelcount);
    free(pixel_changed);
    free(crj_iteration);
    free(saturated);

    if (verbose) printf("done!\n");
    return;
}
    

#ifdef __STANDALONE__
#define NMAX 2000000

void main()
{

    printf("test\n");

    int sx=512, sy=512, i, j;

    //double* data = (double*)malloc(sx*sy*sizeof(double));
    //double* retval = (double*)malloc(sx*sy*sizeof(double));

    double neighbors_a[250], neighbors_b[250], neighbors_c[250];
    double median1, median2, median3;
    clock_t c1, c2, c3, c4;
    double time_a = 0, time_b = 0, time_c=0;

    
    printf("N=%d\n", NMAX);
    printf(" -2 % 6 = %d\n", -2%6);
    printf("  8 % 6 = %d\n",  8%6);
    printf("nan=%d\n", 0./0.);
    printf("nan==nan = %d\n", nan==nan);
    printf("isnan(nan)=%d\n", isnan(nan));
    printf("5==5=%d\n", 5==5);
    
    printf("nan!=nan = %d\n", nan!=nan);
    
        
    int random_numbers[NMAX];
    for (i=0; i<NMAX; i++) {
        // Create random numbers in the range -200 ... +800
        random_numbers[i] = rand()%1000 - 200;
    }
    c1 = clock();
    for (i=0; i<NMAX; i++) {
        j = random_numbers[i] % 600;
    }
    c2 = clock();
    for (i=0; i<NMAX; i++) {
        j = random_numbers[i] >= 600 ? random_numbers[i] - 600 : (random_numbers[i] < 0 ? random_numbers[i] + 600 : random_numbers[i]);
    }
    c3 = clock();
    time_a = (double)(c2 - c1)/CLOCKS_PER_SEC;
    time_b = (double)(c3 - c2)/CLOCKS_PER_SEC;
    printf("Timing a=%f, b=%f, a/b=%f\n", time_a, time_b, time_a/time_b);
        
    return;
    
    int n = 25;
    for (i=0; i<5000000; i++) {
        // Create some random numbers

        for (j=0; j<n; j++) {
            neighbors_a[j] = (double)(rand()%1000);
            neighbors_b[j] = neighbors_a[j];
            neighbors_c[j] = neighbors_a[j];
            // printf("% 4d%s", (int)neighbors[j], (j<n-1 ? ", " : " --> "));
        }
        c1 = clock();
        median1 = heapMedian3(neighbors_c, n);
        c2 = clock();
        gsl_sort(neighbors_b, 1, n);
        median2 = gsl_stats_median_from_sorted_data(neighbors_b, 1, n);
        c3 = clock();
        median3 = find_median(neighbors_a, n);
        c4 = clock();
        
        time_a += (double)(c2 - c1)/CLOCKS_PER_SEC;
        time_b += (double)(c3 - c2)/CLOCKS_PER_SEC;
        time_c += (double)(c4 - c3)/CLOCKS_PER_SEC;
        
        /* for (j=0; j<n; j++) { */
        /*     printf("% 4d%s", (int)neighbors_b[j], (j<n-1 ? ", " : " --> ")); */
        /* } */
        /* printf("%.1f / %.1f / %.1f\n", median1, median2, median3); */
    }
    printf("Timing a=%f, b=%f, c=%f, a/b=%f\n", time_a, time_b, time_c, time_a/time_b);
    

    c1 = clock();
    for (i=0; i<5000000; i++) {
        for (j=0; j<n; j++) {
            neighbors_a[j] = ((double)rand()/1000.);
        }
        median1 = heapMedian3(neighbors_c, n);
    }
    c2 = clock();
    time_a = (double)(c2 - c1)/CLOCKS_PER_SEC;
    

    c1 = clock();
    for (i=0; i<5000000; i++) {
        for (j=0; j<n; j++) {
            neighbors_b[j] = ((double)rand()/1000.);
        }
        gsl_sort(neighbors_b, 1, n);
        median2 = gsl_stats_median_from_sorted_data(neighbors_b, 1, n);
    }
    c2 = clock();
    time_b = (double)(c2 - c1)/CLOCKS_PER_SEC;
    

    c1 = clock();
    for (i=0; i<5000000; i++) {
        for (j=0; j<n; j++) {
            neighbors_c[j] = ((double)rand()/1000.);
        }
        median3 = find_median(neighbors_c, n);
    }
    c2 = clock();
    time_c = (double)(c2 - c1)/CLOCKS_PER_SEC;
    
    printf("Timing a=%f, b=%f, c=%f, a/b=%f\n", time_a, time_b, time_c, time_a/time_b);




    
    return;
    
}

#endif
