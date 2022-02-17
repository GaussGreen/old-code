#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random_utils.h"
#include "proba_utils.h"
#include "alpha_stable_random.h"
#include "gumbel.h"
#include "error2.h"

int GumbelCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta,
    long seed)                   /* (I) */
{
    static char routine[] = "GumbelCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double log_prob;
    double w;

    double *X = NULL;
    double *alphaSequence = NULL;
    void *randomAlpha = NULL;


    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    alphaSequence = malloc(nbPaths*sizeof(double));
    if(alphaSequence==NULL) goto RETURN;

    status = CreateUniformSequence( X,
                                    seed,
                                    nbNames,
	                                nbPaths);
    if(status == FAILURE) goto RETURN;

    randomAlpha = CreateRandomGenerator(seed-3);   
    status = RandomAlphaStableSkew( randomAlpha,
                                alphaSequence,
                                nbPaths,
                                1.0,
                                1./theta,
                                1.0);
    RandomGeneratorFree(randomAlpha);

    if(status == FAILURE) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        log_prob = pow(-log(survivalProba[i])*cos(M_PI_2/theta), theta);
        for(j=0;j<nbPaths;j++)
        {

            if( log_prob >  - log(X[i+j*nbNames]) / alphaSequence[j])
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
    if(alphaSequence) free(alphaSequence);
    if(X) free(X);
    return status;
}


int GumbelCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) 0 if name defaulted before T */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta)                   /* (I) */
{
    static char routine[] = "GumbelCopulatedUniformDeviates";
    int status = FAILURE;
    int i,j;

    double *X = NULL;
    double *alphaSequence = NULL;
    void *randomAlpha = NULL;


    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    alphaSequence = malloc(nbPaths*sizeof(double));
    if(alphaSequence==NULL) goto RETURN;

    status = CreateUniformSequence( X,
                                    -9,
                                    nbNames,
	                                nbPaths);
    if(status == FAILURE) goto RETURN;

    randomAlpha = CreateRandomGenerator(-3);   
    status = RandomAlphaStableSkew( randomAlpha,
                                alphaSequence,
                                nbPaths,
                                1.0,
                                1./theta,
                                1.0);
    RandomGeneratorFree(randomAlpha);

    if(status == FAILURE) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            copulatedUniformDeviates[i+j*nbNames] = 
                  exp(- 1./ cos(M_PI_2/theta) * pow(-log(X[i+j*nbNames]) / alphaSequence[j] , 1./theta ));    
        }
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(alphaSequence) free(alphaSequence);
    if(X) free(X);
    return status;
}



int GumbelCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta,
    long seed)                   /* (I) */
{
    static char routine[] = "GumbelCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double log_prob;
    double w;

    double *X = NULL;
    double *alphaSequence = NULL;
    void *randomAlpha = NULL;


    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    alphaSequence = malloc(nbPaths*sizeof(double));
    if(alphaSequence==NULL) goto RETURN;

    status = CreateUniformSequence( X,
                                    seed,
                                    nbNames,
	                                nbPaths);
    if(status == FAILURE) goto RETURN;

    randomAlpha = CreateRandomGenerator(seed-3);   
    status = RandomAlphaStableSkew( randomAlpha,
                                alphaSequence,
                                nbPaths,
                                1.0,
                                1./theta,
                                1.0);
    RandomGeneratorFree(randomAlpha);

    if(status == FAILURE) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;  
            for(k=0;k<nbTimes;k++)
            {
                log_prob = pow(-log(survivalProba[k+i*nbTimes])*cos(M_PI_2/theta), theta);
                if( log_prob >  - log(X[i+j*nbNames]) / alphaSequence[j])
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
    if(alphaSequence) free(alphaSequence);
    if(X) free(X);
    return status;
}


int GumbelTest( long nbNames,                   /* (I) */
                long nbPaths,                   /* (I) */
                double theta,
                double *prob)                   /* (I) */
{
    static char routine[] = "GumbelTest";
    int status = FAILURE;
    int i,j;
    double w;

    double *X = NULL;
    double *alphaSequence = NULL;
    void *randomAlpha = NULL;
    int in = 1;
    double jointProba = 0.0;
    double copulaProba = 0.0;

    w = 1./nbPaths;
    X = malloc(nbNames*nbPaths*sizeof(double));
    if(X==NULL) goto RETURN;
    alphaSequence = malloc(nbPaths*sizeof(double));
    if(alphaSequence==NULL) goto RETURN;

    status = CreateUniformSequence( X,
                                    -7,
                                    nbNames,
	                                nbPaths);
    if(status == FAILURE) goto RETURN;

    randomAlpha = CreateRandomGenerator(-12);
    status = RandomAlphaStableSkew( randomAlpha,
                                alphaSequence,
                                nbPaths,
                                1.0,
                                1./theta,
                                1.0);
    RandomGeneratorFree(randomAlpha);

    if(status == FAILURE) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            X[i+j*nbNames] = exp(-pow(-log(X[i+j*nbNames])/alphaSequence[j],1./theta)/cos(M_PI_2/theta));    
        }
    }

    for(j=0;j<nbPaths;j++)
    {
        in = 1;
        for(i=0;i<nbNames;i++)
        {
            if(X[i+j*nbNames]>prob[i])
            {
                in = 0;
                break;
            }
        }
        if(in == 1)
        {
            jointProba += w;
        }
    }

    for(i=0;i<nbNames;i++)
    {
        copulaProba += pow(-log(prob[i]),theta);
    }
    copulaProba = exp(-pow(copulaProba,1.0/theta));

    printf("Th joint proba = %lf , actual = %lf \n",copulaProba,jointProba);
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    if(alphaSequence) free(alphaSequence);
    if(X) free(X);
    return status;
}
