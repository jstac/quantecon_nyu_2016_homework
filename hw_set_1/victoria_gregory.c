/*
 * Homework Set 1
 * 
 * Victoria Gregory, N14207660
 *
 * Shows downward bias in the OLS estimate of the
 * correlation coefficient in an AR(1) regression.
 *
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_fit.h>
     
/* This function returns the OLS estimate of alpha */
double ar1_ts (double * params, int n,  unsigned long int seed)
{
    double beta = params[0];
    double alpha = params[1];
    double s = params[2];

    /* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    int i;
    double x = beta / (1 - alpha);  // Start at mean of stationary dist

    /* create arrays to fill with x values */
    double *xveclag = (double*)malloc( n * sizeof(double) );
    double *xvec = (double*)malloc( n * sizeof(double) );
    for (i = 1; i <= n; i++) {
        xveclag[i-1] = x;
        x = beta + alpha * x + gsl_ran_gaussian(r, s);
        xvec[i-1] = x;
     }
    gsl_rng_free (r);

    /* Run OLS on X and its lag */
    double alphahat, betahat, cov00, cov01, cov11, sumsq;
    gsl_fit_linear (xveclag, 1, xvec, 1, n, &betahat, &alphahat, &cov00, &cov01, &cov11, &sumsq);
    return alphahat;
}

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N = 100;   // use a smaller sample size
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double alpha_est = ar1_ts(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("estimated alpha = %g\n", alpha_est);
    printf("true alpha = %g\n", alpha);
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
}
