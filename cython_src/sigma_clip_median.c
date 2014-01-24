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

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

#define nan 0./0.

// #define printf //


void sigma_clip_median__cy(double* pixels, int n_pixels, int n_images, double* output,
                           double nsigma, int max_repeat)
{
    int x,y,l, repeat;

    int n_bytes = n_images*sizeof(bool);
    // printf("n bytes to hold mask: %d\n", n_bytes);
        
    bool *good_value = (bool*)malloc(n_images*sizeof(bool));
    // printf("\n\n\ngood-value=%x\n\n\n",good_value);

    double *pixelvalue = (double*)malloc(n_images*sizeof(double));

    
    for (x=0; x<n_pixels ; x++) {
        // set output to NaN, this is the default
        output[x] = nan; 
        
        /* printf("\n\n\n\n\n\n\rColumn %d", x+1); */
        
        //
        // Figure out how many pixels are not NaN
        //
        
        int n_good_pixels = 0;
        for (y=0; y<n_images; y++) {
            if (!isnan(pixels[x*n_images+y])) {
                good_value[y] = true;
                pixelvalue[n_good_pixels] = pixels[x*n_images+y];
                n_good_pixels++;
            }
        }
        if (n_good_pixels < 1) {
            // all pixels appear to be NaNs, nothing left to do
            continue;
        }
        
            
        //
        // Now sort all values
        //
        gsl_sort(pixelvalue, 1, n_good_pixels);
        /* for (l=0; l<n_good_pixels; l++) { */
        /*     printf("%2d -> %.0f\n", l, pixelvalue[l]); */
        /* } */
        /* printf("\n"); */
        
        
        double median, sigma_plus, sigma_minus, upper, lower, sigma;

        int start=0;
        int end=n_good_pixels;
        int i, new_start, new_end, n_valid_values;

        new_start = start;
        new_end = end;
        
        /* printf("Starting iterating: start/end=%d,%d, new-start/end=%d,%d\n", start, end, new_start, new_end); */
        /* for (i=0; i<end; i++) { */
        /*     printf("%2d: %f\n", i, pixelvalue[i]); */
        /*     } */
 
        for (repeat=0; repeat<max_repeat && (end-start)>=3; repeat++) {
            end = new_end;
            start = new_start;
            
            n_valid_values = end - start;

            /* printf("Iteration %d: start/end=%d,%d #=%d\n", repeat, start, end, n_valid_values); */

            // Compute median and the sigma-widths
            median = gsl_stats_median_from_sorted_data(&pixelvalue[start], 1, n_valid_values);
            sigma_minus = median - gsl_stats_quantile_from_sorted_data(&pixelvalue[start], 1, n_valid_values, 0.16);
            sigma_plus  = gsl_stats_quantile_from_sorted_data(&pixelvalue[start], 1, n_valid_values, 0.84) - median;

            sigma = 0.5 * (gsl_stats_quantile_from_sorted_data(&pixelvalue[start], 1, n_valid_values, 0.84) -
                           gsl_stats_quantile_from_sorted_data(&pixelvalue[start], 1, n_valid_values, 0.16));
            
            // Compute the valid range of pixels
            lower = median - nsigma * sigma;
            upper = median + nsigma * sigma;
            /* lower = median - nsigma * sigma_minus; */
            /* upper = median + nsigma * sigma_plus; */
            
            // Now figure out the start and end range of pixels within the range
            /* printf("Lower limit:\n"); */
            for (new_start=start; new_start<end; new_start++) {
                if (pixelvalue[new_start] >= lower) {
                    /* printf("Value %d is %f >= %f, taking as new lower limit\n", */
                    /*        new_start, pixelvalue[new_start], lower); */
                    break;
                }
                /* printf("Value %d is %f < %f, skipping\n", */
                /*        new_start, pixelvalue[new_start], lower); */
            }
            for (new_end = end; new_end > new_start;) {
                new_end--;
                if (pixelvalue[new_end] <= upper) {
                    /* printf("Value %d is %f <= %f, taking as new upper limit\n", */
                    /*        new_end, pixelvalue[new_end], upper); */
                    // Make stick to nomenclature that end means the index of the first
                    // entry that is no longer contained in the list
                    new_end++;
                    break;
                }
                /* printf("Value %d is %f > %f, skipping\n", */
                /*        new_end, pixelvalue[new_end], upper); */
            }
            
            
            /* printf("Was: %d/%d now is %d,%d\n", */
            /*        start, end, new_start, new_end); */
            /* printf("Iteration %d: median=%f (sigma= %f / %f) lower/upper: %f %f #=%d\n", */
            /*        repeat, median, sigma_minus, sigma_plus, lower, upper, n_valid_values); */
            /* printf("%d --> %d D=%d (%d)\n", */
            /*        pixelvalue, &pixelvalue[start], &pixelvalue[start]-pixelvalue, start); */
            /* printf("Remaining valid values:\n"); */
            /* for (i=new_start; i<new_end; i++) { */
            /*     printf("%f\n", pixelvalue[i]); */
            /* } */

            if (new_start == start && new_end == end) {
                // This iteration hasn't changed anything, so future iteration won't change
                // anything either, so we can stop right here
                /* printf("No change in this iteration skipping rest of work.\n"); */
                break;
            }
            
        }
        if ((new_end - new_start) >= 3) {
            end = new_end;
            start = new_start;
        }
        
        //
        // Now we have all pixels limited down to iteratively select only the
        // 3-sigma range for each pixel;  Now compute the mean value.
        //
        
        output[x] = gsl_stats_median_from_sorted_data(&pixelvalue[start], 1, (end-start));
        // output[x] = (double)(end-start); //gsl_stats_mean(&pixelvalue[start], 1, (end-start));
        // printf("filtered output: %f (%d, %d)\n", output[x], start, end);
        
    }

    // printf("\n\ngood-value=%x\n\n",good_value);
    /* free(good_value); */
        
    // printf("Freeing memory\n");
    free((void*)good_value);
    // printf("done freeing\n");
    
    // printf("Work done in c_sigclip.c\n");
    
    return;
}
    

#ifdef __STANDALONE__

void main()
{

    printf("test\n");

    int n_pix=4096, n_images=9;

    double* data = (double*)malloc(n_pix*n_images*sizeof(double));
    double* retval = (double*)malloc(n_pix*sizeof(double));

    int l;
    for (l=0; l<n_pix*n_images; l++) {
        data[l] = rand()/1e7;
        
    }
    
   
    sigma_clip_mean__cy(data, n_pix, n_images, retval);
    
    
    return;
    
}

#endif
