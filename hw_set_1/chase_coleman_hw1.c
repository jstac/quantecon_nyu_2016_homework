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
 * Example usage:
 * > make all
 * For default parameters
 * > ./chase_coleman_hw1
 * For chosen parameters
 * > ./chase_coleman_hw1 alpha beta sig
 */
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <math.h>
#include <time.h>

void simulatex(double alpha, double beta, double sig, int simlength, double *x)
{
    /* Prepare random number generator */
    unsigned long seed = 61089;
    const gsl_rng_type * T;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    gsl_rng * r;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    /* Declare types */
    double w;
    int t = 1;

    /* Begin at stationary mean */
    x[0] = beta / (1-alpha);
    for(int t=1; t<simlength; t++){
        w = gsl_ran_gaussian(r, 1.0);
        x[t] = alpha*x[t-1] + beta + sig*w;
    }

}

int main(int argc, char *argv[])
{
    /* Initialize clock parameters */
    clock_t start, end;
    double time_used;

    /* Initialize process parameters */
    double alpha, beta, sig;
    if (argc == 1){
        alpha = 0.7;
        beta = 0.0;
        sig = 0.5;
    }
    else if(argc == 4){
        alpha = atof(argv[1]);
        beta = atof(argv[2]);
        sig = atof(argv[3]);
    }
    else {
        printf("Need either 0 or 3 arguments from command line!!!");
        exit(EXIT_FAILURE);
    }

    int maxsim = 1505;
    int nsims[5] = {10, 25, 250, 500, 1500};

    /* Allocate space for arrays */
    double sim[maxsim];
    double x[maxsim-1];
    double y[maxsim-1];

    /* Allocate for linear regression output */
    double c0, c1, cov00, cov01, cov11, sumsq;

    /* Print problem description */
    printf("Consider process described by (alpha, beta, sigma) = (%0.2f, %0.2f, %0.2f)\n", alpha, beta, sig);
    printf("We will show that the OLS estimate of alpha is biased downward\n");
    printf("by generating a table of regressions for various numbers of simulations.\n");

    /* Start clock and run simulation */
    start = clock();
    simulatex(alpha, beta, sig, maxsim, sim);

    /* Fill vectors to pass to linear regression */
    int i = 0;
    for (i=0; i<maxsim-1; i++){
        /* Set x and y values for linear regression */
        x[i] = sim[i];
        y[i] = sim[i+1];
    }

    /* Linear regression for each n*/
    int n;
    printf("  N   |    a   |  ahat  |   diff |\n");
    printf("------|--------|--------|--------|\n");
    for (i=0; i<5; i++){
        n = nsims[i];
        gsl_fit_linear(x, 1, y, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

        printf("%*d|%8.4f|%8.4f|%8.4f|\n", 6, n, alpha, c1, alpha-c1);
    }

    /* End timer */
    end = clock();
    time_used = ((double) end - start) / CLOCKS_PER_SEC;

    /* Print finishing messages */
    printf("\n\n\n");
    printf("Computation took %0.4f", time_used);

    return 0;
}

