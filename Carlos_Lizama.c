/* Homework Set 1
 * Carlos Lizama, N17866309
 * Shows downward bias in the OLS estimate of
 * the correlation coeffiecient in an AR(1) regression
 */

// OLS estimator
// simulate ar(1), pretty much the same as in the example in class

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


double ols_estimator (double * params, int n, unsigned long int seed){
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
    double X[n];        // define array to save data
    X[0] = beta/(1- alpha); // start at mean of stationary dist
    double sum = X[0];     // sum n terms of the time series
    
    for (i = 1; i <= n-1; i++){
        X[i] = beta + alpha * X[i-1] + gsl_ran_gaussian(r,s);
        sum += X[i];
    }
    gsl_rng_free(r);
    double mean1 = (sum-X[n-1])/(n-1);   // mean of terms x_t
    double mean2 = (sum-X[0])/(n-1); // mean of terms x{t+1}

    // compute alpha by OLS using formula for simple linear reg.
    double num = 0;     // numerator
    double den = 0;     // denominator
    
    for (i=1; i<= n-1; i++){
        num += (X[i-1]-mean1)*(X[i]-mean2);
        den += (X[i-1]-mean1)*(X[i-1]-mean1);
    }
    return num/den;
}

// main function
int main(void){
    clock_t start, end;
    double cpu_time_used;

    int N = 100;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double alphahat = ols_estimator(params, N, 1);
    end = clock();
    cpu_time_used = ((double)(end - start))/ CLOCKS_PER_SEC;

    printf("real alpha = %g\n", alpha);
    printf("estimated alpha = %g\n", alphahat);
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
}
