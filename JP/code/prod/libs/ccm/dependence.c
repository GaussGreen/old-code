#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random_utils.h"
#include "proba_utils.h"
#include "gaussian.h"
#include "error2.h"

int DependenceCopulatedUniformDeviates(
    double *copulatedUniformDeviates, 
    long nbNames, 
    long nbPaths)
{
    static char *routine = "DependenceCopulatedUniformDeviates";
    long i = 0;
    long j = 0;
    int status = FAILURE;
    double *uniformDeviates = malloc(nbPaths*sizeof(double));
    if(uniformDeviates==NULL) goto RETURN;

    status = CreateSobolSequence(   uniformDeviates,
                                    -7,
                                    1,
                                    nbPaths);
    if(status == FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedUniformDeviates[i+j*nbNames] = uniformDeviates[j];
        }
    }
    status = SUCCESS;

RETURN:
    if(uniformDeviates) free(uniformDeviates);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}


/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/

/**
 * Calculate survival indicator function for this copula.
 */
int DependenceCopulatedIndicator(
    int* copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) */
    const double *survivalProba,          /* (O) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,
    long seed)                   /* (I) */
{
    static char routine[] = "DependenceCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double beta = 0.0;
    double *gaussianSequence = NULL;
    double *T = NULL;

    gaussianSequence = malloc(nbPaths*sizeof(double));
    if(gaussianSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T == NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, 1, nbPaths, &beta, seed);
    if(status==FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(gaussianSequence[j] > T[i])
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
    if(status == FAILURE)
    {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

/**
 * Calculate survival indicator function for this copula.
 */
int DependenceCopulatedIndicator_mtp(
    int* copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) */
    const double *survivalProba,          /* (O) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,
    long seed)                   /* (I) */
{
    static char routine[] = "DependenceCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double beta = 0.0;
    double *gaussianSequence = NULL;
    double *T = NULL;

    gaussianSequence = malloc(nbPaths*sizeof(double));
    if(gaussianSequence == NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T == NULL) goto RETURN;

    for(i=0;i<nbNames*nbTimes;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, 1, nbPaths, &beta, seed);
    if(status==FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            for(k=0;k<nbTimes;k++)
            {
                copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
                if(gaussianSequence[j] > T[k+i*nbTimes])
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
    if(status == FAILURE)
    {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

#if 0
void DependenceCopulatedTime(
    int* copulatedDefaultTimeIndex, /* (O) p if name defaulted between Tp-1 and Tp */
    double *survivalProba,          /* (O) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths)                   /* (I) */
#endif
