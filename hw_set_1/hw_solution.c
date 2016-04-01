/* 
 * Solution to the hw exercise on downward bias of AR(1) coefficient.
 *
 * John Stachurski
 *
 * The plan is to set a fixed value of alpha and then simulate m time
 * series using that value.  For each such time series we estimate alpha
 * and take the mean of these estimates.  Since the individual time 
 * series are independent, this mean gives an estimate of the expectation
 * of the sampling distribution of the OLS estimator of alpha.
 *
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_fit.h>

int ar1_ts (double * x, double * params, int n, unsigned long int seed)
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
    x[0] = beta / (1 - alpha);  // Start at mean of stationary dist
    for (i = 1; i < n; i++) 
     {
       x[i] = beta + alpha * x[i-1] + gsl_ran_gaussian(r, s);
     }

    gsl_rng_free (r);
    return 0;
}


int fit_params(double *x, int n)
{
    double alpha_hat, beta_hat, cov00, cov01, cov11, chisq;
    double * y = x + 1;
    gsl_fit_linear (x, 1, y, 1, n-1, 
                &beta_hat, &alpha_hat, &cov00, &cov01, &cov11, &chisq);
    return alpha_hat;
}

int main(void) 
{
    /* length of time series and parameter values */
    int n = 100;  
    double alpha = 0.9;
    double beta = 1.0;
    double s = 0.1;
    double params[] = {beta, alpha, s};

    int num_reps = 1000; // number of replications
    double x[n];  // allocate memory for time series
    int i;

    double alpha_hat;
    double mean = 0.0;

    for (i=0; i < num_reps; i++) 
    {
        ar1_ts(x, params, n, i);
        alpha_hat = fit_params(x, n);
        mean += alpha_hat;
    }
    mean = mean / (double) num_reps;

    printf("True value of alpha = %g\n", alpha);
    printf("Mean alpha_hat over %d simulations = %g\n", num_reps, mean);

    return 0;
}



