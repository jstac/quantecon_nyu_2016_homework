/* Homework Set 1

* Ruixue Gong, N17593858

* Shows downward bias in the OLS estimate of the
* correlation coefficient in an AR(1) regression

* Relies on GNU-GSL library: http://www.gnu.org/software/gsl/

*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>


double OLS_estimate (double * params, int n, unsigned long int seed){

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
    
    /*Define a variable x to store a sequence of observations,which is also
    time series*/
    double x[n];
    x[0] = beta/(1- alpha); // start at mean of stationary dist


    /* Calculate the first n-term sum of time series x_t*/
    double sum = x[0];   //Define a new variable sum with initial value=x_0      
    for (i = 1; i <= n; i++){
        x[i] = beta + alpha * x[i-1] + gsl_ran_gaussian(r,s);
        sum += x[i];
    }
    gsl_rng_free(r);
   

    /* Calculate sample means of time series x_t and x_{t+1}*/
    double mean1 = (sum-x[n])/n;   // mean of terms x_t
    double mean2 = (sum-x[0])/n; // mean of terms x{t+1}

 
    /* Calculate alpha by OLS regression.*/
    double numerator = 0;     //(X'X)^{-1}
    double denominator = 0;     // XY
    
    for (i=1; i<= n; i++){
        numerator += (x[i-1]-mean1)*(x[i]-mean2);
        denominator += (x[i-1]-mean1)*(x[i-1]-mean1);
    }
   
    return numerator / denominator;  //This is the estimate of alpha
}

/*Result*/

int main(void){
    clock_t start, end;
    double cpu_time_used;
    double bias; 

    int N = 100000;
    double beta = 1.0;
    double alpha = 0.9;
    double s = 1;
    double params[3] = {beta, alpha, s};

    start = clock();
    double alpha_hat = OLS_estimate(params, N, 1); //The estimate of alpha
    end = clock();
    cpu_time_used = ((double)(end - start))/ CLOCKS_PER_SEC;
    bias=alpha_hat-alpha;   //Calculate the bias of alpha estimate


    printf("Real alpha = %g\n", alpha);
    printf("Estimated alpha = %g\n", alpha_hat);
    printf("time elapsed = %g\n", cpu_time_used);
    printf("The bias of alpha estimate =%g\n", bias);
    
    return 0;

}
