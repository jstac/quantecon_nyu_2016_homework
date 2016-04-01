vi/*
 * Homework set 1
 *
 * Felipe Alves, N14713445
 *
 * Shows downward bias in the OLS estimate of the correlation coefficient
 * in AR(1) regression
 *
*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include <gsl/gsl_fit.h>

void sim_ar1_ts (double *params, double *x, double *x1, double *eps, int n, unsigned long int seed)
{
    double beta  = params[0];
    double alpha = params[1];
    double s     = params[2];

    /* create a generator chosen by the environment variable GSL_RNG_TYPE */
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(T);
    gsl_rng_set(r, seed);

    int i;
    double ee;

    *(x1+0)      = beta / (1 - alpha);  // Start at mean of stationary dist
    for (i = 0; i <= n; i++) {
        ee = gsl_ran_gaussian(r, s);
        *(x+i)   = beta + alpha * x1[i] + ee;
        *(eps+i) = ee;
        i<n ? *(x1+(i+1)) = x[i] : 0;
     }

}

void read_array (double *x)
{
    printf("X[0] = %g\n", x[0]);
}

int main(void)
{
    /* Set sample size and repetitions */
    int Nsizes[4] = {40, 100, 200, 500}, rep = 10000; 

    double Ealpha[4];
    
    /* Parameters */
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s}; //declaring a matrix I guess
    for (int i = 0; i < 4; ++i){
    
    int N = Nsizes[i];
    double X[N], X1[N], eps[N], sum_alpha=0;
    double c0, c1, cov00, cov01, cov11, chisq;

        for (int j = 0; j < rep; ++j){    
        /* Create my AR series */
        sim_ar1_ts (params, X, X1, eps, N, j);

        /* OLS estimator */
        gsl_fit_linear (X1, 1, X, 1, N, 
                        &c0, &c1, &cov00, &cov01, &cov11, &chisq);

        sum_alpha += c1;
         }
    
    Ealpha[i] = sum_alpha/rep;
     }
    printf("     \n");
    printf("  N    E[α]\n");
    printf("---- | ----\n");
    for (int i = 0; i < 4; ++i) {
    printf("%3.d    %.3g\n", Nsizes[i], Ealpha[i]);
    }
    
    return 0;
}

