/*
* Homework Set 1
*
* James Graham, jag912, N11252710
*
* This file simulates an AR(1) process and then runs an OLS regression to 
* calculate the AR(1) coefficients in a small sample. It shows that there
* is a downward bias in the OLS computed AR(1) correlation coefficient.
*
* Note the need to use 'malloc' to pre-allocate memory to the variable 
* length arrays/vectors.
*/



#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_fit.h>



double ar1_ts (double * params, int n, unsigned long int seed)
{
    double beta = params[0];
    double alpha = params[1];
    double s = params[2];

    /* Note: need to allocate size for the vectors using 'malloc' */
    double* x;
    x = malloc(n * sizeof(double));

    /* Note: create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
   
    // Simulate the AR(1) process 
    int i;
    x[0] = beta / (1 - alpha);  /* Start at mean of stationary dist */
    for (i = 1; i <= n; i++) {
        x[i] = beta + alpha * x[i-1] + gsl_ran_gaussian(r, s);
     }

    // Create the data vectors for OLS
    double* X;
    X = malloc(n * sizeof(double));
    double* Xlag;
    Xlag = malloc(n * sizeof(double));

    int j;
    for (j = 0; j<=n-1; j++)
    {
    X[j] = x[j+1];
    Xlag[j] = x[j];
    }

    // Run OLS using GSL linear projection toolbox 
   double c0;
   double c1; 
   double cov00;
   double cov01;
   double cov11;
   double sumsq;
   gsl_fit_linear (Xlag, 1, X, 1, n, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

   printf("RESULTS OF OLS ON THE AR(1) SIMULATION\n");
   printf("---------------------------------------\n");
   printf("\t \t \t OLS \t | Actual\n");
   printf("Constant \t |  %g \t | %g\n", c0, beta);
   printf("AR(1) coeff \t |  %g \t | %g\n", c1, alpha);
   printf("---------------------------------------\n");
   printf("Sample Size = %d\n", n);
   printf("---------------------------------------\n");

    gsl_rng_free (r);
    return 0;
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
    ar1_ts(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
   
    printf("Time elapsed = %g\n", cpu_time_used);
    return 0;
}




