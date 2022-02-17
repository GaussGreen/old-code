#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random_utils.h"
#include "proba_utils.h"
#include "gamma_random.h"
#include "gumbel.h"
#include "error2.h"

int ClaytonCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha,
    long seed)             /* (I) */
{
    static char routine[] = "ClaytonCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double log_prob;
    double w;

    double *X = NULL;
    double *gammaSequence = NULL;
    void *randomGamma = NULL;
    

    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    gammaSequence = malloc(nbPaths*sizeof(double));
    if(gammaSequence==NULL) goto RETURN;


    status = CreateUniformSequence( X,
                                    seed,
                                    nbNames,
                                    nbPaths);
    if(status == FAILURE) goto RETURN;

    randomGamma = CreateRandomGenerator(seed-3);
    status = RandomGamma(       randomGamma,
                                gammaSequence,
                                nbPaths,
                                1. / alpha,
                                1.);
    RandomGeneratorFree(randomGamma);

    if(status == FAILURE) goto RETURN;


    for(i=0;i<nbNames;i++)
    {
        log_prob = pow(survivalProba[i], - alpha);
        for(j=0;j<nbPaths;j++)
        {

            if( log_prob >= (1.0 - log(X[i+j*nbNames]) / gammaSequence[j] ))
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 0;    
            }
            else
            {
                copulatedSurvivalIndicator[i+j*nbNames] = 1;
            }
        }
    }

    w = 1./nbPaths;
    for(j=0;j<nbPaths;j++)
    {
        weight[j] = w; 
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(gammaSequence) free(gammaSequence);
    if(X) free(X);
    return status;
}


int ClaytonCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) 0 if name defaulted before T */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha)             /* (I) */
{
    static char routine[] = "ClaytonCopulatedUniformDeviates";
    int status = FAILURE;
    int i,j;

    double *X = NULL;
    double *gammaSequence = NULL;
    void *randomGamma = NULL;
    

    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    gammaSequence = malloc(nbPaths*sizeof(double));
    if(gammaSequence==NULL) goto RETURN;


    status = CreateUniformSequence( X,
                                    -9,
                                    nbNames,
                                    nbPaths);
    if(status == FAILURE) goto RETURN;

    randomGamma = CreateRandomGenerator(-3);
    status = RandomGamma(       randomGamma,
                                gammaSequence,
                                nbPaths,
                                1. / alpha,
                                1.);
    RandomGeneratorFree(randomGamma);

    if(status == FAILURE) goto RETURN;


    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            copulatedUniformDeviates[i+j*nbNames] =
                pow ((1.0 - log(X[i+j*nbNames]) / gammaSequence[j] ) , -1./alpha);    
        }
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(gammaSequence) free(gammaSequence);
    if(X) free(X);
    return status;
}



int ClaytonCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha,
    long seed)             /* (I) */
{
    static char routine[] = "GumbelCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double log_prob;
    double w;

    double *X = NULL;
    double *gammaSequence = NULL;
    void *randomGamma = NULL;
    

    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    gammaSequence = malloc(nbPaths*sizeof(double));
    if(gammaSequence==NULL) goto RETURN;


    status = CreateUniformSequence( X,
                                    seed,
                                    nbNames,
                                    nbPaths);
    if(status == FAILURE) goto RETURN;

    randomGamma = CreateRandomGenerator(seed-3);
    status = RandomGamma(       randomGamma,
                                gammaSequence,
                                nbPaths,
                                1. / alpha,
                                1.);
    RandomGeneratorFree(randomGamma);

    if(status == FAILURE) goto RETURN;


    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
            for(k=0;k<nbTimes;k++)
            {
                log_prob = pow(survivalProba[k+i*nbTimes], - alpha);
                if( log_prob >= (1.0 - log(X[i+j*nbNames]) / gammaSequence[j] ))
                {
                    copulatedSurvivalIndicator[i+j*nbNames] = k;
                    break;
                }
            }
        }
    }

    w = 1./nbPaths;
    for(j=0;j<nbPaths;j++)
    {
        weight[j] = w; 
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(gammaSequence) free(gammaSequence);
    if(X) free(X);
    return status;
}

