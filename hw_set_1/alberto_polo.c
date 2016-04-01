/*
* Homework set 1
*
* Alberto Polo, NYUID ap3562
*
* Shows downward bias in the OLS estimate of the correlation coefficient in an AR(1) regression
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

double regr_on_ar1 (double * params, int n, unsigned long int seed)
{
    double beta = params[0]; // location parameter
    double alpha = params[1]; // persistence parameter
    double s = params[2]; // standard deviation of Gaussian innovations
    /* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    int i;
    double x = beta / (1 - alpha);  // Start at mean of stationary dist
    double xnext = x; // Same as above
    double sum = 0; // Keep track of sum of realizations to compute the mean
    double sum2 = 0;  // Keep track of sum of squared realizations to compute the 2nd moment
    double prod = 0;  // Keep track of cross-products to compute covariance
    for (i = 1; i <= n; i++) {
        sum += xnext;
        sum2 += xnext * xnext;
        xnext = beta + alpha * x + gsl_ran_gaussian(r, s);
        prod += xnext * x;
        x = xnext;
     }

    gsl_rng_free (r);
    double mean = sum / n;
    double var = sum2 / n - mean * mean;
    double cov = prod / n - mean * mean;
    return cov / var; // Return the OLS estimate along the simulated path of the process
}

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N = 10000;  // Set length of each simulated path
    double beta = 1.0;  // Set location parameter
    double alpha = 0.99;  // Set persistence parameter
    double s = 1; // Set standard deviation of innovations
    double params[3] = {beta, alpha, s};

    start = clock();
    int j;
    double sum_regr = 0; // Keep track of sum of OLS estimates to compute their mean afterwards
    int simul = 10000; // Set number of simulated OLS coefficients
    for (j = 1; j <= simul; j++) {
        sum_regr += regr_on_ar1(params, N, j);
     }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    double mean_regr = sum_regr / simul;

    printf("mean of OLS estimator = %f\n", mean_regr);
    printf("true value of autocorrelation = %f\n", alpha);
    if(alpha - mean_regr > 0.00001) {
      printf("downward biased!\n");
    }
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
}
