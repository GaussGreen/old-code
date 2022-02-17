#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "random_utils.h"
#include "proba_utils.h"
#include "gaussian.h"
#include "error2.h"

/* ---------------------------------------------------------------------------
// GaussianMultivariateDeviates
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Gaussian beta-correlated deviates
*/
int GaussianMultivariateDeviates(double *gaussianSequence,
                                  long nbNames,
                                  long nbPaths,
                                  const double *beta)
{
    static char routine[] = "GaussianMultivariateDeviates";
    int status = FAILURE;
    long seed = -7;
    long i = 0;
    long j = 0;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *sqrt_beta = NULL;

    Z = malloc(nbNames * nbPaths * sizeof(double));
    if(Z==NULL) goto RETURN;
    M = malloc(nbPaths * sizeof(double));
    if(M==NULL) goto RETURN;
    sqrt_beta = malloc(nbNames*sizeof(double));
    if(sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] = sqrt(1. - beta[i]*beta[i]);
    }

    /* generate the variables Z_i with gaussian random sequence */
    seed = -7;
    status = CreateGaussianRandomSequence(Z,seed,nbNames,nbPaths);
    if(status == FAILURE) goto RETURN;   
    
    /* generate the market variable with gaussian random variable */
    seed = -3;
    status = CreateGaussianRandomSequence(M,seed,1,nbPaths);
    if(status == FAILURE) goto RETURN;   

    /* get the Students variable */
    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            gaussianSequence[i+j*nbNames] =
                (beta[i]*M[j] + sqrt_beta[i]*Z[i+j*nbNames]);
        }
    }
    status = SUCCESS;
RETURN:
        if(Z) free(Z);
        if(M) free(M);
        if(sqrt_beta) free(sqrt_beta);
        if(status == FAILURE)
        {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        }
        return status;
}

/* ---------------------------------------------------------------------------
// GaussianCopulatedUniformDeviates
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Gaussian copulatd uniform deviates
*/
int GaussianCopulatedUniformDeviates(double *copulatedUniformDeviates,
                                      long nbNames,
                                      long nbPaths,
                                      const double *beta)
{
    static char routine[] = "GaussianCopulatedUniformDeviates";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    status = GaussianMultivariateDeviates(   copulatedUniformDeviates,
                                    nbNames,
                                    nbPaths,
                                    beta);
    if(status == FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedUniformDeviates[i+j*nbNames] =
                NormalCum(copulatedUniformDeviates[i+j*nbNames]);
        }
    }
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}

/* -------------------------------------------------------------------------------
// GaussianCopulaLinearCorrelation
// This function returns the actual linear correlation given
// the correlation parameter of the copula
*/
double GaussianCopulaLinearCorrelation(double beta)
{
    long i,j;
    double x_1 = -3.0;
    double x_2 = -3.0;
    double stepSize = 0.05;
    long nbSteps = 120;
    double integral = 0;
    double sqrt_beta = sqrt(1-beta*beta);
    for(i=0;i<nbSteps;i++)
    {
        x_2 = -3.0;
        for(j=0;j<nbSteps;j++)
        {
            integral +=
                NormalCum(x_1)*NormalCum(beta*x_1+sqrt_beta*x_2)*
                NormalDensity(x_1)*NormalDensity(x_2)*stepSize*stepSize;
            x_2 += stepSize;
        }
        x_1 += stepSize;
    }
    return 12*integral-3;
}


/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/

/* ---------------------------------------------------------------------------
// GaussianMultivariateDeviates_tp
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Gaussian beta-correlated deviates
// and the weights associatd to each path
*/
int GaussianMultivariateDeviates_tp(double *gaussianSequence,
                                  double *weight,
                                  long nbNames,
                                  long nbPaths,
                                  const double *beta,
                                  long seed)
{
    static char routine[] = "GaussianMultivariateDeviates_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *sqrt_beta = NULL;

    Z = malloc(nbNames * nbPaths * sizeof(double));
    if(Z==NULL) goto RETURN;
    M = malloc(nbPaths * sizeof(double));
    if(M==NULL) goto RETURN;
    sqrt_beta = malloc(nbNames*sizeof(double));
    if(sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] = sqrt(1-beta[i]*beta[i]);
    }

    /* generate the variables Z_i with gaussian random sequence */
    status = CreateGaussianRandomSequence(Z,seed,nbNames,nbPaths);
    if(status == FAILURE) goto RETURN;   
    
    /* generate the market variable with gaussian random variable */
    seed -= 3;
    status = CreateGaussianRandomSequence(M,seed,1,nbPaths);
    if(status == FAILURE) goto RETURN;   

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            gaussianSequence[i+j*nbNames] =
                (beta[i]*M[j] + sqrt_beta[i]*Z[i+j*nbNames]);
        }

        // WEIGHTS
        weight[j] = 1./nbPaths;
    }
    status = SUCCESS;

RETURN:
    if(Z) free(Z);
    if(M) free(M);
    if(sqrt_beta) free(sqrt_beta);

        if(status == FAILURE)
        {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        }
     return status;
}


/* ---------------------------------------------------------------------------
// GaussianMultivariateDeviatesSampling_tp
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Gaussian beta-correlated deviates
// and the weights associatd to each path
*/
int GaussianMultivariateDeviatesSampling_tp(double *gaussianSequence,
                                  double *weight,
                                  long nbNames,
                                  long nbPaths,
                                  const double *beta,
                                  long seed,
                                  const double *M,
                                  long nbSample)
{
    static char routine[] = "GaussianMultivariateDeviatesSampling_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;
    long k = 0;

    /* nb of path for each sample of M */
    long nbPaths_by_sample = nbPaths / nbSample;
    long nbPaths_remaining = nbPaths % nbSample;
    long nbPaths_sample;

    double totalWeight = 0.0;
    double totalWeightInverse = 0.0;
    double w = 0.0;
    long sample = 0;

    /* allocation of the Sequence */
    double *Z = NULL;
    double *sqrt_beta = NULL;
    

    Z = malloc(nbNames * (nbPaths_by_sample) * sizeof(double));
    if(Z==NULL) goto RETURN;

    sqrt_beta = malloc(nbNames*sizeof(double));
    if(sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] = sqrt(1-beta[i]*beta[i]);
    }

    /* generate the variables Z_i with gaussian random sequence */
    seed = -7;
    status = CreateGaussianRandomSequence(Z,seed,nbNames,nbPaths_by_sample);
    if(status == FAILURE) goto RETURN;   

    for(k=0;k<nbSample;k++)
    {

        nbPaths_sample = nbPaths_by_sample;
        if(k<nbPaths_remaining)
        {
            nbPaths_sample += 1;
        }

        if(k==0)
        {
            w = NormalCum((M[k+1] + M[k])*0.5);
        }
        else if(k==nbSample-1)
        {
            w = (1.0 - NormalCum((M[k-1] + M[k])*0.5));
        }
        else
        {
            w = (NormalCum((M[k+1] + M[k])*0.5) - NormalCum((M[k-1] + M[k])*0.5));
        }

        w /= nbPaths_sample;

        for(j=0;j<nbPaths_sample;j++)
        {
            for(i=0;i<nbNames;i++)
            {
                gaussianSequence[i+sample*nbNames] =
                    (beta[i] * M[k] + sqrt_beta[i] * Z[i+j*nbNames]);
            }

            // WEIGHTS
            weight[sample] = w;
            totalWeight += weight[sample];
            sample +=1;
        }
    }

    totalWeightInverse = 1./totalWeight;
    for(j=0;j<nbPaths;j++)
    {
        weight[j] *= totalWeightInverse;
    }
    status = SUCCESS;

RETURN:
    if(Z) free(Z);
    if(sqrt_beta) free(sqrt_beta);

    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}

/* ---------------------------------------------------------------------------
// GaussianCopulatedIndicator
// this function returns the survival indicator and the weight of each path
//
*/
int GaussianCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    long seed,
    const double *M,
    long nbSample)
{
    static char routine[] = "GaussianCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *gaussianSequence = NULL;
    double *T = NULL;

    gaussianSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(gaussianSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    if(M==NULL || nbSample == 0)
    {
        status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, nbNames, nbPaths, beta,seed);
        if(status ==FAILURE) goto RETURN;
    }
    else
    {
        status = GaussianMultivariateDeviatesSampling_tp(gaussianSequence, weight, nbNames, nbPaths, beta, seed, M, nbSample);
        if(status ==FAILURE) goto RETURN;
    }

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(gaussianSequence[i+j*nbNames] > T[i])
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 0;
            }
            else
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 1;
            }
        }
    }
    status = SUCCESS;
RETURN:
    if(gaussianSequence) free(gaussianSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

/* ---------------------------------------------------------------------------
// GaussianCopulatedIndicator_mtp
// this function returns the survival indicator and the weight of each path
//
*/
int GaussianCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    long seed,
    const double *M,
    long nbSample)
{
    static char routine[] = "GaussianCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double *gaussianSequence = NULL;
    double *T = NULL;

    gaussianSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(gaussianSequence == NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames*nbTimes;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    if(M==NULL || nbSample == 0)
    {
        status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, nbNames, nbPaths, beta,seed);
        if(status ==FAILURE) goto RETURN;
    }
    else
    {
        status = GaussianMultivariateDeviatesSampling_tp(gaussianSequence, weight, nbNames, nbPaths, beta, seed, M, nbSample);
        if(status ==FAILURE) goto RETURN;
    }

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
            for(k=0;k<nbTimes;k++)
            {
                if(gaussianSequence[i+j*nbNames] > T[k+i*nbTimes])
                {
                    copulatedSurvivalIndicator[i+j*nbNames] = k;
                    break;
                }
            }
        }
    }
    status = SUCCESS;
RETURN:
    if(gaussianSequence) free(gaussianSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}



/***********************************************************************/
/***********************************************************************/
/*              STOCHASTIC RECOVERY IMPLEMENTATION                     */
/***********************************************************************/

/* ---------------------------------------------------------------------------
// GaussianCopulatedIndicatorRecovery
// this function returns the survival indicator and the weight of each path
//
*/
int GaussianCopulatedIndicatorRecovery(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    double *X,                      /* (O) recovery stoch term */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    const double *alpha,
    long seed,
    const double *M,
    long nbSample)
{
    static char routine[] = "GaussianCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *beta_temp = NULL;
    double *gaussianSequence = NULL;
    double *T = NULL;

    beta_temp = malloc(2*nbNames*sizeof(double));
    if(beta_temp==NULL) goto RETURN;
    gaussianSequence = malloc(2*nbNames*nbPaths*sizeof(double));
    if(gaussianSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T==NULL) goto RETURN;
    
    memcpy(beta_temp,beta,nbNames*sizeof(double));
    memcpy(&(beta_temp[nbNames]),alpha,nbNames*sizeof(double));

    for(i=0;i<nbNames;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    if(M==NULL || nbSample == 0)
    {
        status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, 2*nbNames, nbPaths, beta_temp, seed);
        if(status ==FAILURE) goto RETURN;
    }
    else
    {
        status = GaussianMultivariateDeviatesSampling_tp(gaussianSequence, weight, 2*nbNames, nbPaths, beta_temp, seed, M, nbSample);
        if(status ==FAILURE) goto RETURN;
    }

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(gaussianSequence[i+2*j*nbNames] > T[i])
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 0;
            }
            else
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 1;
            }

            X[i+j*nbNames] = gaussianSequence[i+(2*j+1)*nbNames];
        }
    }
    status = SUCCESS;
RETURN:
    if(gaussianSequence) free(gaussianSequence);
    if(T) free(T);
    if(beta_temp) free(beta_temp);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}
