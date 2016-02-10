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
    double xl = 0;
    double x = beta / (1 - alpha);
    double sum_x = 0;  
    double sum_xl = 0;
    double sum_xl_sqr = 0;
    double sum_xtxl = 0;
    double mean_x;
    double mean_xl;
    double mean_xl_sqr;
    
    for (i = 1; i <= n; i++) {
        
         xl = x;
         x = beta + alpha * x + gsl_ran_gaussian(r, s); 
         sum_x += x;
         sum_xl += xl;
         sum_xl_sqr += xl * xl;
         sum_xtxl += x * xl;      
        
         mean_x = sum_x/n;
         mean_xl = sum_xl/n;
         mean_xl_sqr = mean_xl * mean_xl;
    }            
 
    gsl_rng_free (r);
    return (sum_xtxl - (n * mean_x * mean_xl)) / (sum_xl_sqr - (n * mean_xl_sqr));
}

int main()
 
 {
  double sum_alpha_hat = 0;
  double mean_alpha_hat;
  int i;
  int n = 20;
  for (i = 1; i <= n; i++) {
    clock_t start, end;
    double cpu_time_used;

    int N = 1000;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};
    start = clock();
    
    double alpha_hat = ar1_ts(params, N, i);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    sum_alpha_hat += alpha_hat; 
    
    printf("estimated alpha = %g\n", alpha_hat);
    printf("time elapsed = %g\n", cpu_time_used);
}
    mean_alpha_hat = sum_alpha_hat / (double) n; 
    printf("average alpha estimate = %g\n", mean_alpha_hat);

    return mean_alpha_hat;
}
