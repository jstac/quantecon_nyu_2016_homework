/*
    Homework 1 for ECON-GA 3002, John Stachurski.
    Arnav Sood, asood@nyu.edu, NYUID as8554.
 
    Setup: AR(1) process, standard normal innovation. This program shows by simulation a downward bias in the OLS estimate of the correlation parameter.
 
    Disclaimer: This code is very environmentally friendly; i.e., it is liberally recycled from the examples John posted to the course GitHub.
 
    Credit goes to Timothy Petliar for independently suggesting that I use a scalar formulation instead of a matrix one, finding the right estimator and how to use it (i.e. separability), and then pointing out the requisite change in RNG.
 
    Credit goes to StackExchange for things like how to use rand().
 */


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>

// A method to provide the standard OLS estimate of alpha, given parameters. We run this a few times on different parameters from range to confirm the fact.

double ar1_estimate (double * params, int n, unsigned long int seed) {
    
    // Bind the parameter values to local variables in our stack.
    
    double beta = params[0];
    double alpha = params[1];
    double s = params[2];
    
    // Create a random number generator, with our seed, and the Gaussian spec.
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    
    // Create another random number generator, based on system time. We'll use this when we don't care about the distribution. Thanks StackExchange.
    
    srand(time(NULL));

    // Run the random number generator. The solution we use as our estimator is given: https://are.berkeley.edu/courses/EEP118/current/derive_ols.pdf, equation 8. Specifically, we compute the product terms first, but keep track of sums so we can compute means.
    
    int i;
    double x = rand() % 100; // Let's see what happens when we start at a random integer between 0 and 100.
    double sumx = 0;
    double sumy = 0;
    double sumxy = 0;
    double sumxsq = 0;
    
    for (i = 1; i < n; i++) { // We use a strict inequality, because n-1 terms are predicting n-1 terms.
        sumx += x;
        sumxsq += pow(x,2);
        // Update
        double y = beta + alpha * x + gsl_ran_gaussian(r, s);
        // Calculate dependent sums.
        sumxy += x*y;
        sumy += y;
        x = y;
    }
    
    // Now, calculate the estimate.
    double estimate = (sumxy - (n-1)*(sumx/(n-1))*(sumy/(n-1)))/(sumxsq - (n-1)*pow((sumx/(n-1)),2));

    // Reset the Gaussian estimator.
    gsl_rng_free (r);
    
    // Return the estimate.
    return estimate;
    
}

int main(void) {
    
    // Set up clock.
    
    clock_t start, end;
    double cpu_time_used;
    
    // Set up params.
    
    int N = 10000000;
    double beta = 1.0;
    double alpha1 = 0.9;
    double alpha2 = 0.5;
    double alpha3 = 0.27182818284;
    double s = 1;
    
    // We keep beta the same because it doesn't matter.
    
    double params1[3] = {beta, alpha1, s};
    double params2[3] = {beta, alpha2, s};
    double params3[3] = {beta, alpha3, s};

    // Run and time code.
    
    start = clock();
    double estimate1 = ar1_estimate(params1,N,1);
    double estimate2 = ar1_estimate(params2,N,2);
    double estimate3 = ar1_estimate(params3,N,3);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    // Print results.
    
    printf("estimate1 = %g\n", estimate1);
    printf("parameter1 = %g\n", alpha1);
    
    printf("estimate2 = %g\n", estimate2);
    printf("parameter2 = %g\n", alpha2);
    
    printf("estimate3 = %g\n", estimate3);
    printf("parameter3 = %g\n", alpha3);
    
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
    
}
