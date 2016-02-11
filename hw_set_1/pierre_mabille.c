/*
* Homework set 1
* Pierre Mabille, N10017621
*
* Shows downward bias in the OLS estimate of the 
* correlation coefficient in an AR(1) regression.

*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
 
// function returning 2 pointers to arrays of simulated Y=X_{t+1} and X=X_t
// takes as arguments: pointer to array of parameters, integer sample size,
// and seed for random number generator
double *ar1_ts (double * params, int n,  unsigned long int seed, double
X[],double Y[])
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
    for (i = 1; i <= n; i++) {
    X[i-1] = x;
    x = beta + alpha * x + gsl_ran_gaussian(r, s);
    Y[i-1] = x;
     }
    return 0;
    gsl_rng_free (r);
}

int gsl_fit_linear(); // declare OLS function used in main()

int main() // main() calls ar1_ts() to simulate the Y and X series and
// gsl_fit_linear() to compute OLS estimator
{
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    // loop through N_array to simulate AR(1) for various sample sizes N
    static int N_array[6] = {10, 100, 1000, 10000, 100000, 500000};
    int j; 
    for (j = 0; j <= 5; j++)
    {
        int N = N_array[j] ;
        double beta = 1.0;
        double alpha = 0.9;
        double s = 1;
        double params[3] = {beta, alpha, s};
    
        int k;
        int Nsim = 100;
        double sum_beta_hat = 0;
        // for given N, simulate for Nsim different seeds
        for (k=1; k<= Nsim; k++)
        {        
            double X[N], Y[N]; // declare series simulated
            ar1_ts(params, N,k,X, Y); // call function simulating AR(1)
            // declare outcomes of OLS estimation and apply function for OLS
            double alpha_hat, beta_hat, cov00, cov01, cov11, sumsq;
            gsl_fit_linear(X,1,Y,1,N, &alpha_hat, &beta_hat, &cov00, &cov01, &cov11, &sumsq);
            sum_beta_hat += beta_hat; // sum the estimated beta coeffs
        }
    // compute the average of estimated betas for given N and print it
    double mean_beta_hat = sum_beta_hat/Nsim; 
    printf("avg. beta_hat for N=%i: %g\n",N,mean_beta_hat);
    }
    // print total computing time
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;   
    printf("time elapsed = %g\n",cpu_time_used);
    return 0;
}

