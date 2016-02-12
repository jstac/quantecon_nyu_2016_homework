/************************************************
* Homework set 1
* 
* Daniel Csaba, UniID dc2730
*
* Shows downward bias in the OLS estimate of the
* correlation coefficient in an AR(1) regression.
*************************************************/


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

/* Estimate the coefficient on the independent variable in an AR(1) model */
double ar1_ols_est (double * params, int n, unsigned long int seed){
    double beta = params[0];
    double alpha = params[1];
    double s = params[2];
    /* Create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    
    // Generate the sample path
    int i;
    double x[n+1]; // Create an array of size n+1 for the data
    x[0] = beta / (1 - alpha); // Start at mean of stationary dist
    for (i = 1; i <= n; i++) {
        x[i] = beta + alpha * x[i-1] + gsl_ran_gaussian(r,s);
    }
    gsl_rng_free (r); // Setting the random number generator free

    // Compute auxiliary variables for estimators
    double sum_x = 0;
    double sum_y = 0;
    double sum_xy = 0;
    double sum_xx = 0;
    for (i = 1; i <= n; i++) {
        sum_x += x[i-1];
        sum_y += x[i];
        sum_xy += x[i-1] * x[i];
        sum_xx += x[i-1] * x[i-1];
    }

    // Compute the OLS estimate for alpha
    double alpha_ols =  (sum_xy - ( sum_x * sum_y)/n) / (sum_xx - (sum_x *
    sum_x)/n);
    return alpha_ols;
}

int main(void)
{   
    clock_t start, end;
    double cpu_time_used;

    int N = 100;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double alpha_hat = ar1_ols_est(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
   

    printf("=====================================================\n");
    printf("Downward bias of OLS estimate in an AR(1) regression.\n");
    printf("=====================================================\n");
    printf("Sample size          = %i\n", N);
    printf("Original parameter   = %g\n", alpha);
    printf("Estimated parameter  = %g\n", alpha_hat);
    printf("Realized bias        = %g\n", alpha_hat - alpha);
    printf("Time elapsed         = %g\n", cpu_time_used);
    printf("=====================================================\n");
    return 0;
}
