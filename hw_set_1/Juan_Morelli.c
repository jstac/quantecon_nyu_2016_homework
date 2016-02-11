#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
     

double  ar1_ts (double * params, int n, unsigned long int seed)
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

    /* create vector AR(1) process */
    int i,j;
    double  x = beta / (1 - alpha);  // Start at mean of stationary dist
    double  A[n];
    double  Apl1[n-1], Ami1[n-1];
    double  c0, c1, cov00, cov01, cov11, sumsq;
    A[0] = x;
    for (i = 1; i < n; i++) {
        x = beta + alpha * x + gsl_ran_gaussian(r, s);
        A[i]=x;
     }

    /* create vectors with time t and t-1 values for input */
     for (j=0;j<n-1; j++) {
     Apl1[j] = A[j+1];
     Ami1[j] = A[j];
     }
 
    /* use function to estimate OLS coefficients */
    gsl_fit_linear(Ami1, 1, Apl1, 1, n-1, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

    return c1;
}
 

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N[3] = {50, 100, 1000};
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();

    /* Will generate 1000 simulations for 3 different sample sizes */
    /* Then I compute the mean for each sample size and print result */
    /* I call the function defined above in each loop */ 
     int j,k;
    double C[3];
    for (j=0;j<3;j++){
        double xx=0;
        double sum = 0;
        for (k=0;k<1001;k++){
        
        sum += xx;
        xx = ar1_ts(params, N[j], 1);
        }
    
        C[j] = sum/1000;
        
        printf("alpha = %g\n", C[j]);
        }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    printf("time elapsed = %g\n", cpu_time_used);

    return 0;
}

