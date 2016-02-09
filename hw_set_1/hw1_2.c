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
    double y_i = beta / (1 - alpha);  // Start at mean of stationary dist
    double sum_x_i = 0; //We are going to set y_i=x_t+1 and x_i=x_t for each t
    double sum_y_i=0; // sum_x_i+1
    double sum_xi_yi=0; //sum of x_i*x_i+1
    double sum_xi_sqrd=0; // sum of x_i^2
    double x_i=0; // x_i
    double mean_x_i; //mean x_i
    double mean_y_i; //mean x_i+1
    double alpha_est; //ols estimate for alpha
    for (i = 1; i <= n; i++) {
        x_i=y_i; //Lag x_i+1 by 1 step
              
        y_i = beta + alpha * y_i + gsl_ran_gaussian(r, s);

         sum_x_i +=x_i;
         sum_xi_sqrd += x_i * x_i;
         sum_xi_yi += x_i * y_i;
         sum_y_i += y_i;
     }
    mean_x_i=sum_x_i/n; 
    mean_y_i=sum_y_i /n;
    alpha_est=(sum_xi_yi - n * mean_x_i *mean_y_i) / (sum_xi_sqrd - n * mean_x_i * mean_x_i);
    gsl_rng_free (r);
    return alpha_est;
}

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N = 10000000;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double alpha_ols = ar1_ts(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    printf("alpha = %g\n", alpha_ols);
    printf("time elapsed = %g\n", cpu_time_used);
    return 0;
}


