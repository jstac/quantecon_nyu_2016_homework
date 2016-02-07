/*
 * Homework Assignment 1
 *
 * Chase Coleman, N10827183
 *
 * Shows downward bias in the OLS estimate of the
 * correlation coefficient in an AR(1) regression
 *
 * Relies on GNU-GSL library: http://www.gnu.org/software/gsl/
 *
 *
 */
#include <stdio.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>

void simulatex(double alpha, double beta, double sig, int nsim, double *x)
{
    /* Prepare random number generator */
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng * r;
    r = gsl_rng_alloc(T);

    /* Declare types */
    double w;
    int t = 1;

    /* Begin at stationary mean */
    x[0] = beta / (1-alpha);
    for(int t=1; t<nsim; t++){
        w = gsl_ran_gaussian(r, 1.0);
        x[t] = alpha*x[t-1] + beta + sig*w;
    }

}

int main(void)
{
    /* Initialize clock parameters */
    clock_t start, end;
    double time_used;

    /* Initialize process parameters */
    double alpha = 0.9;
    double beta = 0.0;
    double sig = 1.0;
    int nsim = 10;

    /* Allocate space for arrays */
    double sim[nsim];
    double x[nsim-1];
    double y[nsim-1];

    /* Allocate for linear regression output */
    double c0, c1, cov00, cov01, cov11, sumsq;

    /* Print problem description */
    printf("Consider process described by (alpha, beta, sigma) = (%0.2f, %0.2f, %0.2f)\n", alpha, beta, sig);
    printf("We will show that the OLS estimate of alpha is biased downwards using %i simulations.\n", nsim);

    /* Start clock and run simulation */
    start = clock();
    simulatex(alpha, beta, sig, nsim, sim);

    /* Fill vectors to pass to linear regression */
    int i = 0;
    for (i=0; i<nsim-1; i++){
        /* Set x and y values for linear regression */
        x[i] = sim[i];
        y[i] = sim[i+1];
    }

    /* Linear regression */
    gsl_fit_linear(x, 1, y, 1, nsim-1, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    /* End timer */
    end = clock();
    time_used = ((double) end - start) / CLOCKS_PER_SEC;

    /* Print finishing messages */
    printf("\n\n\n\n\n");
    printf("The true coefficients are [%0.4f, %0.4f]\n", beta, alpha);
    printf("The estimated coefficients are [%0.4f, %0.4f]\n", c0, c1);
    printf("Computation took %0.4f", time_used);

    return 0;
}

