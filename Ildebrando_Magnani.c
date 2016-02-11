/* 
*
* Homework set 1
*
* Ildebrando Magnani, im975
*
* Show downward bias in OLS estimate of the correlation
* coefficient in a AR(1) regression.
*
*
* The program runs 100 different AR1 simulations and estimates the alpha
* coefficient (called alpha_hat) 100 different times. Hence, it takes the average
* of these 100 estimates of alpha and displays it as "evidence" of the
* downward bias (which shows up on average).
*/


#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
     

double ar1_ts (double * params, int n, unsigned long int seed)   /* function that simulates AR1 and estimates the correlation coefficient */
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
    double xl = 0;                        /* declares var lagged x */
    double x = beta / (1 - alpha);        /* declares var x and initializes it, sets initial value = mean of stationary distribution */
    double sum_x = 0;                     /* declares var sum of x variable observations */  
    double sum_xl = 0;                    /* declares sum of x lagged variable observations */
    double sum_xl_sqr = 0;           /* sum of (x lagged squared) */
    double sum_xtxl = 0;             /* sum of (x * x lagged) */
    double mean_x;                   /* sample mean of x observations */
    double mean_xl;                  /* sample mean of x lagged observations */
    double mean_xl_sqr;              /* (sample mean of x lagged)^2 */ 
    
    for (i = 1; i <= n; i++) {              /* loop that generates time series */
        
         xl = x;                      /* lags the x var */
         x = beta + alpha * x + gsl_ran_gaussian(r, s); /* AR1 */ 
         sum_x += x;                                    /* sum_x = new observation of x + previous sum */
         sum_xl += xl;                                  /* same as line above but with x lagged */
         sum_xl_sqr += xl * xl;                         /* same with each x lagged squared */                                                 
         sum_xtxl += x * xl;                            /* sums each x * x lagged */
         mean_x = sum_x/n;                 /* computes sample mean of x */
         mean_xl = sum_xl/n;               /* computes sample mean of x lagged */
         mean_xl_sqr = mean_xl * mean_xl;  /* squares the sample mean of x lagged */
    }            
 
    gsl_rng_free (r);
    return (sum_xtxl - (n * mean_x * mean_xl)) / (sum_xl_sqr - (n * mean_xl_sqr)); /* returns the estimate for alpha (correlation coeff) 
                                                                                      the standard OLS slope estimator Cov(x,xl)/Var(xl) */
}

int main(void)          /* main function that runs the AR1 simulation and correlation coefficient                          
                           estimation 100 different times. Then takes the
                           average of the 100 alpha_hat and displays it */
 {
  double sum_alpha_hat = 0;      /* declares and initializes var for the sum of each estimates for alpha */
  double mean_alpha_hat;       /* declares var for the mean of the 100 alpha estimates */
  double total_cpu_time;        
  double beta = 1.0;
  double alpha = 0.9;
  double s = 1;
  double params[3] = {beta, alpha, s};
  int i;
  int n = 100;                  /* number of AR1 simulations and estimations for correlation coeff alpha_hat */
  
  for (i = 1; i <= n; i++) {   /* loop that repeats the estimation process for alpha until i <= n */
    clock_t start, end;       
    double cpu_time_used;
    
    int N = 200;             /* number of AR1 observations */ 
    start = clock();
    
    double alpha_hat = ar1_ts(params, N, i);                  
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    total_cpu_time += cpu_time_used;         /* sums the cpu time used for each simulation/estimation */
    sum_alpha_hat += alpha_hat;              /* sums all the 100 alpha hats estimated */
    
}
    mean_alpha_hat = sum_alpha_hat / (double) n; /* computes the average estimated alpha. Will be displayed as AVERAGE ALPHA ESTIMATE below */ 
    printf("true alpha = %g\n", alpha);        /* displays true alpha */
    puts("AVERAGE ALPHA ESTIMATE is the mean of 100 alpha_hats estimated from 100 different AR(1) simulations.");
    puts("In this way, we show how the downward bias shows up, on average.");
    printf("AVERAGE ALPHA ESTIMATE = %g\n", mean_alpha_hat); /* displays AVERAGE ALPHA ESTIMATE, which is the mean of 100 alpha_hats
                                                                estimated from 100 different AR1 simulations */
    printf("time elapsed = %g\n", total_cpu_time);
   
   return 0;
}
