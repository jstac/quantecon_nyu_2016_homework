/*Homework 1 
* Laszlo Tetenyi n14167755
* Bias in the OLS estimate
*/ 
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
double ar1_ts (double * params, int n, unsigned long int seed) {
    double beta = params[0];
    double alpha = params[1];
    double s = params[2];/* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);
    int i;
    double x = beta / (1-alpha);  // Start at mean of stationary dist
    double xprev = beta / (1-alpha); // Previous observation - starts at x
    double xmean = beta / (1-alpha); // mean of x
    double xprevmean =  beta / (1-alpha); // mean of xprev
    double sumcov = 0; // sum of covariance term
    double sumvar = 0; //sum of variance term
    for (i = 1; i <= n; i++) {
        sumcov += (x-xmean)*(xprev-xprevmean);
        sumvar += (xprev-xprevmean)* (xprev-xprevmean);
        xprev = x;
        x=beta+alpha*x+gsl_ran_gaussian(r,s);
        xmean =( xmean*i + x) /(i+1);
        xprevmean=xmean;
    }
    gsl_rng_free (r);
    return sumcov/sumvar;
}


int main(void) {
    clock_t start, end;
    double cpu_time_used;
    int N1=10000000;
    int N2=100;
    double beta=1.0;
    double alpha=0.9;
    double s =1;
    double params[3]={beta,alpha,s}; 
    double OLS_consi=ar1_ts(params,N1,1);
    start=clock();
    double OLS_bias=ar1_ts(params,N2,1);
    end=clock();
    cpu_time_used=((double)(end-start))/CLOCKS_PER_SEC;
    printf("Ols is consistent as for large samplesize the OLS estimate is %g\n",OLS_consi);
    printf("OLS is biased: for small samplesize the estimate is %g\n",OLS_bias);
    printf("the OLS bias is  %g\n",OLS_bias - alpha);
    printf("time elapsed for small sample %g\n",cpu_time_used);
    return 0;
}
