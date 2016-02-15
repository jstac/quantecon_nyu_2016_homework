/*
 * Homework set 1
 *
 * Ryuichiro Izumi, UniID ri505
 * Shows downward bias in the OLD estimate of the
 * correlation coefficient in an AR(1) regression.
 *
 */

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#include <gsl/gsl_randist.h>
#include <time.h>

int main(void)
{
    clock_t start, end;
    double cpu_time_used;

    int N=10000;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double seed=1;
    /* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    int i;
    double x[N];
    double a = beta / (1 - alpha);  // Start at mean of stationary dist
    x[1]=a;
    for (i = 1; i <= N; i++) {
        x[i+1] = beta + alpha * x[i] + gsl_ran_gaussian(r, s); // Constructing AR(1) process
     }

    gsl_rng_free (r);
/*
    double params[3] = {beta, alpha, s};
*/
    start = clock();

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    double y[N-1];
    for (i=0; i <N; i++){
       y[i+1]=x[i+2];		// Explained variable X_t expressed as y here
    }

    double x2[N-1];
    double xy[N-1];
    double sum1=0;
    double sum2=0;
    int w=N-1;
    for (i=0; i<w; i++){
       x2[i+1]=x[i+1]*x[i+1];	// X'X part of OLS estimator
       xy[i+1]=x[i+1]*y[i+1];	// X'Y part
       sum1=sum1+x2[i+1];
       sum2=sum2+xy[i+1];
    }
    double betahat;
    betahat=sum1/sum2;		// (X'X)^(-1) (X'Y)
    /*printf("mean = %g\n", sample_mean);
*/
    printf("time elapsed = %g\n", cpu_time_used);
    printf("betahat = %g\n", betahat);
    return 0;

}

