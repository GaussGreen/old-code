#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "gaussian.h"
#include "random_utils.h"
#include "proba_utils.h"
#include "error2.h"

int IndependenceCopulatedUniformDeviates(double *copulatedUniformDeviates,
                                          long nbNames,
                                          long nbPaths)
{
    long i = 0;
    long j = 0;
    int status = FAILURE;
    static char routine[] = "IndependenceCopulatedUniformDeviates";
    double *uniformDeviates = malloc(nbPaths*nbNames*sizeof(double));
    if(uniformDeviates==NULL) goto RETURN;

    status = CreateSobolSequence(   uniformDeviates,
                                    -7,
                                    nbNames,
                                    nbPaths);
    if(status == FAILURE)
    {
        goto RETURN;
    }
    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedUniformDeviates[i+j*nbNames] = uniformDeviates[i+j*nbNames];
        }
    }
    status = SUCCESS;
RETURN:
    if(copulatedUniformDeviates) free(copulatedUniformDeviates);
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
int IndependenceCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,
    long nbPaths,
    long seed)
{
    static char routine[] = "IndependenceCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double *beta = NULL;
    double *gaussianSequence = NULL;
    double *T = NULL;
    
    beta = malloc(nbNames*sizeof(double));
    if(beta==NULL) goto RETURN;
    gaussianSequence = malloc(nbPaths*nbNames*sizeof(double));
    if(gaussianSequence==NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        beta[i] = 0.0;
    }

    for(i=0;i<nbNames*nbTimes;i++)
    {
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, nbNames, nbPaths, beta, seed);
    if(status == FAILURE) goto RETURN;

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
    if(beta) free(beta);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}


int IndependenceCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbNames,
    long nbPaths,
    long seed)
{
    static char routine[] = "IndependenceCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *beta = NULL;
    double *gaussianSequence = NULL;
    double *T = NULL;
    
    beta = malloc(nbNames*sizeof(double));
    if(beta==NULL) goto RETURN;
    gaussianSequence = malloc(nbPaths*nbNames*sizeof(double));
    if(gaussianSequence==NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        beta[i] = 0.0;
        T[i] = NormalCumInverse(survivalProba[i]);
    }

    /* allocate the gaussian sequence */
    status = GaussianMultivariateDeviates_tp(gaussianSequence, weight, nbNames, nbPaths, beta, seed);
    if(status == FAILURE) goto RETURN;

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
    if(beta) free(beta);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}
