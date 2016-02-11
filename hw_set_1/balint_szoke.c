/*
* Homework set 1
*
* Balint Szoke, UniID bs2574
*
* Shows downward bias in the OLS estimate of the
* correlation coefficient in an AR(1) regression.
*
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>



double *ar1_sample_path (double *params, int n, unsigned long seed){
    /* ar1_sample_path returns an array x so we define the function as pointer 
    * and the array as static variable. It is tricky, though, because we want its length 
    * to be given by n, so it is a VLA (variable length) and cannot be defined with      
    * static. To determine the size in run-time, allocate it on the heap by malloc */
    
    /* Pull out the three parameters from params */
    double beta, alpha, sigma;
    beta = params[0], alpha = params[1], sigma = params[2];
    
    /* Seed the generator from GSL_RNG_TYPE 
        * gsl_rng_env_setup() reads GSL_RNG_TYPE and sets gsl_rng_default 
        * gsl_rng_alloc() returns an instance of a rng
        * gsl_rng_set() seeds the rng */
    gsl_rng *r;             // define a generator
    const gsl_rng_type *T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    /* Simulate an ar(1) sample path */
    static double *x;
    x = (double*)malloc(n * sizeof(double));
    double y = beta / (1-alpha); // always start it at the mean of the stat dist

    int i;
    for (i=0; i<n; i++) {
        x[i] = y;
        y = beta + alpha*y + gsl_ran_gaussian(r, sigma);
    }

    gsl_rng_free(r); // frees all the memory associated with generator r 

    return x;
}



double ols_alpha(double *data, int n ){
/* Calculates the OLS estimators for alpha
*     data  - array containing a sample path for x
*     n     - size of the sample, i.e. if sample size has length T, then n=T-1 
*/


    /* Need: sum(xy), sum(xx), sum(x) and sum(y), where y and x are the LHS and RHS vars */
    double sum_xy, sum_xx, sum_x, sum_y;
    sum_xy = 0, sum_xx = 0, sum_x = 0, sum_y = 0;

    int i;    
    for (i=0; i<n; i++ ){
        /* x[i] = data[i] and y[i]=data[i+1] */
        sum_xy += data[i]*data[i+1];
        sum_xx += data[i]*data[i];
        sum_x += data[i];
        sum_y += data[i+1];
    }
    
    double alpha = ( sum_xy - (sum_x*sum_y)/n )/( sum_xx - (sum_x*sum_x)/n );
    
    return alpha;
}


int main(void){
    
    double beta, alpha, sigma;
    beta = 1.0, alpha = 0.99, sigma=1;
    double params[3] = {beta, alpha, sigma};

    int N1, N2, N3, N4, N5;
    N1 = 11, N2 = 51, N3 = 101, N4 = 501, N5 = 1001;

    double *sample_path = ar1_sample_path(params, N5, 10021987);
    
    double alpha_hat_N1 = ols_alpha(sample_path, N1-1);
    double alpha_hat_N2 = ols_alpha(sample_path, N2-1);
    double alpha_hat_N3 = ols_alpha(sample_path, N3-1);
    double alpha_hat_N4 = ols_alpha(sample_path, N4-1);
    double alpha_hat_N5 = ols_alpha(sample_path, N5-1);
         
    
    printf("==============================================================\n");
    printf(" DOWNWARD BIAS IN THE OLS ESTIMATE OF THE AR(1) COEFFICIENT:  \n");
    printf("==============================================================\n");
    printf("True value of alpha = %g\n", alpha);
    printf("--------------------------------------------------------------\n");
    printf("OLS estimates for alpha with \n");
    printf("sample size = %d   |  %g\n", N1-1, alpha_hat_N1);
    printf("sample size = %d   |  %g\n", N2-1, alpha_hat_N2);
    printf("sample size = %d  |  %g\n", N3-1, alpha_hat_N3);
    printf("sample size = %d  |  %g\n", N4-1, alpha_hat_N4);
    printf("sample size = %d |  %g\n", N5-1, alpha_hat_N5);
    printf("==============================================================\n");

    return 0;
}



