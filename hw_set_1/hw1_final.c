/* Homework set 1
*
*
* Timothy Petliar, N12592901 , NetID: tdp259
*
*Shows downward bias in the OLS estimate of the
*correlation coefficient in the AR(1) regression
*
* This program uses different seeds in the random generator. This
*modification is critical because the bias may be upward for certain specific
*paths, which is an issue if one seed is used, generating only one path. */


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
    int k;
    double alpha_est_sum=0; //sum of larger sample estimate over all seeds
    double l_alpha_est_sum=0;//sum of smaller sample estimate over all seeds
    double l_alpha_est;
    double alpha_est;
    for (k=1; k<=20; k++) {
    /* create a generator chosen by the environment variable GSL_RNG_TYPE with 20 different seeds */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r,k*10 ); 

    int i;
    
    double y_i = beta / (1 - alpha);  // Start at mean of stationary dist
    double sum_x_i = 0; //We are going to set y_i=x_t+1 and x_i=x_t for each t
    double sum_y_i=0; // sum_x_i+1
    double sum_xi_yi=0; //sum of x_i*x_i+1
    double sum_xi_sqrd=0; // sum of x_i^2
    double x_i=0; // x_i
    double mean_x_i; //mean x_i
    double mean_y_i; //mean x_i+1

    for (i = 1; i <= n; i++) {
        x_i=y_i; //Lag x_i+1 by 1 step
              
        y_i = beta + alpha * y_i + gsl_ran_gaussian(r, s);

         sum_x_i +=x_i;
         sum_xi_sqrd += x_i * x_i;
         sum_xi_yi += x_i * y_i;
         sum_y_i += y_i;
        
        if (i==100) { //Compute ols for small sample
            mean_x_i=sum_x_i/i;
            mean_y_i=sum_y_i /i;
            l_alpha_est=(sum_xi_yi - i* mean_x_i *mean_y_i) / (sum_xi_sqrd - i * mean_x_i * mean_x_i); 
         
} 
         
           
}

         mean_x_i=sum_x_i/n; //Computed after end of i loop for each k
         mean_y_i=sum_y_i /n;
         alpha_est=(sum_xi_yi -n * mean_x_i *mean_y_i) / (sum_xi_sqrd -n * mean_x_i * mean_x_i); //compute ols for total sample
      if (k==2)//print example of upward biased estimate
        { printf("For some paths the estimate is biased upward: %g\n",alpha_est);} 
        alpha_est_sum += alpha_est; //sum total sample estimate over seeds
        l_alpha_est_sum += l_alpha_est; //sum small sample estimate over seeds
        gsl_rng_free (r);       
    

 }
    
    printf("Average ols estimate of alpha for N=100 is %g\n", l_alpha_est_sum/20); //dislplay avg over all seeds
   
   printf("For N=10,000 estimate error is very small but still downward biased");
   
   printf(". Average ols estimate for N=10,000 is %g\n", alpha_est_sum/20);   


    
    return 0;
}

int main(void)
{
    
 
    clock_t start, end;
    double cpu_time_used;
    int N =10000;
   
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
   printf("alpha is actually %g\n", alpha);
   ar1_ts(params, N, 1);
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
   
    printf("time elapsed = %g\n", cpu_time_used);
    
    
    return 0;
}


