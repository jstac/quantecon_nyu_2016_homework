#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
     

double ar1_ts (double * params, int n, unsigned long int seed)
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
    double x = beta + (alpha * (beta / (1 - alpha))) + gsl_ran_gaussian(r, s);
    double x_minus = beta / (1 - alpha);  // Start at mean of stationary dist
    double sum_x = 0;  
    double sum_x_minus = 0;
    double x_minus_sqr = x_minus * x_minus;
    double xtx_minus = x * x_minus;
    double sum_x_minus_sqr = 0;
    double sum_xtx_minus = 0;
    double mean_x;
    double mean_x_minus;
    double mean_x_minus_sqr;
    for (i = 1; i <= n; i++) {
        
         x = beta + alpha * x_minus + gsl_ran_gaussian(r, s);
         x = x_minus;        
         sum_x += x;
         sum_x_minus += x_minus;
         sum_x_minus_sqr = x_minus_sqr;
         sum_xtx_minus += xtx_minus;
         sum_x_minus_sqr += x_minus_sqr;
         mean_x = sum_x/n;
         mean_x_minus = sum_x_minus/n;
         mean_x_minus_sqr = mean_x_minus * mean_x_minus;
    
         
    }

    gsl_rng_free (r);
    return (sum_xtx_minus - (n * mean_x * mean_x_minus)) / (sum_x_minus_sqr - (n * mean_x_minus_sqr));
}

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N = 10;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double sample_mean = ar1_ts(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("mean = %g\n", sample_mean);
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
}
