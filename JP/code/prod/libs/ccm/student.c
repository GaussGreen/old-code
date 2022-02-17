#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "random_utils.h"
#include "proba_utils.h"
#include "student.h"
#include "error2.h"

/*--------------------------------------------------------------
// Chi2Deviates
// this function returns an array[nbPaths]
// of chi2 deviates with FREEDOMDEGREE 
// using the Sobol sequence
*/
int Chi2DeviatesSobol(double *chi2,
                  long nbPaths,
                  long freedomDegree)
{
    static char routine[] = "Chi2Deviates"; 
    int status = FAILURE;
    long seed = 0; // seed has no effect for Sobol
    long i,j;
    double x;
    double *sobolSequence = malloc(freedomDegree * nbPaths * sizeof(double));
    if (sobolSequence==NULL) goto RETURN;

    status = CreateSobolSequence(   sobolSequence,
                                    seed,
                                    freedomDegree,
                                    nbPaths);
    if(status == FAILURE) goto RETURN;

    /* transform the uniform sobol sequence to normal */
    for(i=0;i<freedomDegree;i++)
    {
        for(j=0;j<nbPaths;j++)
        {
            sobolSequence[i+j*freedomDegree] =
                NormalCumInverse(sobolSequence[i+j*freedomDegree]);
        }
    }
    /* compute the chi2 variables */
    for(j=0;j<nbPaths;j ++)
    {
        chi2[j]=0;
        for(i=0;i<freedomDegree;i++)
        {
            x = sobolSequence[i+j*freedomDegree];
            chi2[j] += x*x;
        }
    }
    status = SUCCESS;
RETURN:
        if(sobolSequence) free(sobolSequence);
        if(status == FAILURE)
        {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        }
        return status;
}

/*--------------------------------------------------------------
// Chi2Deviates
// this function returns an array[nbPaths]
// of chi2 deviates with FREEDOMDEGREE 
// using the Sobol sequence
*/
int Chi2Deviates(double *chi2,
                  long nbPaths,
                  long freedomDegree,
                  long seed)
{
    static char routine[] = "Chi2Deviates"; 
    int status = FAILURE;
    long i,j;
    double x;
    double *sequence = malloc(freedomDegree * nbPaths * sizeof(double));
    if (sequence==NULL) goto RETURN;

    status = CreateGaussianRandomSequence(  sequence,
                                            seed,
                                            freedomDegree,
                                            nbPaths);
    if(status == FAILURE) goto RETURN;

    /* compute the chi2 variables */
    for(j=0;j<nbPaths;j ++)
    {
        chi2[j]=0;
        for(i=0;i<freedomDegree;i++)
        {
            x = sequence[i+j*freedomDegree];
            chi2[j] += x*x;
        }
    }
    status = SUCCESS;
RETURN:
        if(sequence) free(sequence);
        if(status == FAILURE)
        {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        }
        return status;
}

/*------------------------------------------------------
// StudentMultivariateDeviates
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Student deviates of the form
// (beta_i*M + sqrt(1-beta_i^2) Z_i) / sqrt(chi2/freedomDegree) 
*/
int StudentMultivariateDeviates(double *studentSequence,
                                 long nbNames,
                                 long nbPaths,
                                 const double *beta,
                                 long freedomDegree)
{
    static char routine[] = "StudentMultivariateDeviates";
    int status = FAILURE;
    long seed = -7;
    long i = 0;
    long j = 0;
    double sqrt_chi2;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *chi2 = NULL;
    double *sqrt_beta = NULL;

    Z = malloc(nbNames * nbPaths * sizeof(double));
    if(Z==NULL) goto RETURN;
    M = malloc(nbPaths * sizeof(double));
    if(M==NULL) goto RETURN;
    chi2 = malloc(nbPaths * sizeof(double));
    if(chi2==NULL) goto RETURN;
    sqrt_beta = malloc(nbNames*sizeof(double));
    if(sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] = sqrt(1. - beta[i]*beta[i]);
    }

    /* generate a chi^2 variable with Sobol sequence */
    status = Chi2DeviatesSobol(chi2, nbPaths, freedomDegree);
    if (status==FAILURE) goto RETURN;

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
        sqrt_chi2 = sqrt(chi2[j]/freedomDegree);
        for(i=0;i<nbNames;i++)
        {
            studentSequence[i+j*nbNames] =
                (beta[i]*M[j] + sqrt_beta[i]*Z[i+j*nbNames])/
                                                            sqrt_chi2;
        }
    }
    status = SUCCESS;

RETURN:
    if(M) free(M);
    if(Z) free(Z);
    if(chi2) free(chi2);
    if(sqrt_beta) free(sqrt_beta);
    if(status == FAILURE)
    {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return status;
}

/* -----------------------------------------------------------------------------
// StudentCopulatedUniformDeviates
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// of Student copulatd uniform deviates
*/
int StudentCopulatedUniformDeviates(double *copulatedUniformDeviates,
                                     long nbNames,
                                     long nbPaths,
                                     const double *beta,
                                     long freedomDegree)
{
    static char routine[] = "StudentCopulatedUniformDeviates";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    status = StudentMultivariateDeviates(copulatedUniformDeviates,
                                nbNames,
                                nbPaths,
                                beta,
                                freedomDegree);
    if(status==FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedUniformDeviates[i+j*nbNames] =
                StudentCum(copulatedUniformDeviates[i+j*nbNames],
                            freedomDegree);
        }
    }
    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}



/* -------------------------------------------------------------------------------
// StudentCopulaLinearCorrelation
// This function returns the actual linear correlation given
// the correlation parameter of the copula
// DOES NOT WORK
*/
double StudentCopulaLinearCorrelation(double beta,
                                      long freedomDegree)
{
    long i,j,k;
    double x_1 = -3.0;
    double x_2 = -3.0;
    double stepSize = 0.05;
    double stepSize_x3 = 0.1;
    double x_3 = stepSize;
    long nbSteps = 120;
    long nbSteps_x3 = 30*freedomDegree;
    double integral = 0;
    double sqrt_beta = sqrt(1-beta*beta);
    for(i=0;i<nbSteps;i++)
    {
        x_2 = -3.0;
        for(j=0;j<nbSteps;j++)
        {
            x_3 = stepSize_x3;
            for(k=0;k<nbSteps_x3;k++)
            {
                integral += StudentCum(x_1/sqrt(x_3/freedomDegree),
                    freedomDegree)*StudentCum((beta*x_1+sqrt_beta*x_2)/
                    sqrt(x_3/freedomDegree), 
                    freedomDegree)*
                    NormalDensity(x_1)*NormalDensity(x_2)*
                    Chi2Density(x_3,freedomDegree)*
                    stepSize*stepSize*stepSize_x3;
                x_3 += stepSize_x3;
            }
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

/*------------------------------------------------------
// StudentMultivariateDeviates_tp
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// and the weight[nbPaths]
// of Student deviates of the form
// (beta_i*M + sqrt(1-beta_i^2) Z_i) / sqrt(chi2/freedomDegree) 
*/
int StudentMultivariateDeviates_tp(double *studentSequence,
                                 double *weight,
                                 long nbNames,
                                 long nbPaths,
                                 const double *beta,
                                 long freedomDegree,
                                 long seed)
{
    static char routine[] = "StudentMultivariateDeviates_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;
    double sqrt_chi2;
    double inverse_freedomDegree = 1./freedomDegree;


    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *chi2 = NULL;
    double *sqrt_beta = NULL;

    Z = malloc(nbNames * nbPaths * sizeof(double));
    if (Z==NULL) goto RETURN;
    M = malloc(nbPaths * sizeof(double));
    if (M==NULL) goto RETURN;
    chi2 = malloc(nbPaths * sizeof(double));
    if (chi2==NULL) goto RETURN;
    sqrt_beta = malloc(nbNames*sizeof(double));
    if (sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] =sqrt(1. - beta[i]*beta[i]); 
    }

    /* generate a chi^2 variable with Sobol sequence */
    status = Chi2Deviates(chi2, nbPaths, freedomDegree,seed);
    if(status == FAILURE) goto RETURN;

    /* generate the variables Z_i with gaussian random sequence */
    seed -= 12;
    status = CreateGaussianRandomSequence(Z,seed,nbNames,nbPaths);
    if(status == FAILURE) goto RETURN;   
    
    /* generate the market variable with gaussian random variable */
    seed -= 3;
    status = CreateGaussianRandomSequence(M,seed,1,nbPaths);
    if(status == FAILURE) goto RETURN;   

    /* get the Students variable */
    for(j=0;j<nbPaths;j++)
    {
        sqrt_chi2 = sqrt(chi2[j]*inverse_freedomDegree);
        for(i=0;i<nbNames;i++)
        {
            studentSequence[i+j*nbNames] =
                (beta[i]*M[j] + sqrt_beta[i]*Z[i+j*nbNames])/
                                                            sqrt_chi2;
        }

        // WEIGHTS
        weight[j] = 1./ nbPaths;
    }
    status = SUCCESS;

    RETURN:
        if(Z) free(Z);
        if(M) free(M);
        if(chi2) free(chi2);
        if(sqrt_beta) free(sqrt_beta);
        if(status == FAILURE)
        {
            DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        }
        return status;
}

/*------------------------------------------------------
// StudentMultivariateDeviatesSampling_tp
// This function returns a matrix[nbNames][nbPaths]
// (stored as an array of the colums)
// and the weight[nbPaths]
// of Student deviates of the form
// (beta_i*M + sqrt(1-beta_i^2) Z_i) / sqrt(chi2/freedomDegree) 
*/
int StudentMultivariateDeviatesSampling_tp(double *studentSequence,
                                 double *weight,
                                 long nbNames,
                                 long nbPaths,
                                 const double *beta,
                                 long freedomDegree,
                                 long seed,
                                 const double *M,
                                 long M_nbSample,
                                 const double *chi2,
                                 long chi2_nbSample)
{
    static char routine[] = "StudentMultivariateDeviatesSampling_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;
    long k = 0;
    long l = 0;
    double inverse_freedomDegree = 1./freedomDegree;
    double sqrt_chi2;
    double M_k, chi2_l;

    long nbSample = chi2_nbSample * M_nbSample;

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
    if (Z==NULL) goto RETURN;
    sqrt_beta = malloc(nbNames*sizeof(double));
    if(sqrt_beta==NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        sqrt_beta[i] = sqrt(1. - beta[i]*beta[i]);   
    }

    /* generate the variables Z_i with gaussian random sequence */
    status = CreateGaussianRandomSequence(Z,seed,nbNames,nbPaths_by_sample);
    if(status == FAILURE) goto RETURN;   
    
    /* get the Students variable */
    for(k=0;k<M_nbSample;k++)
    {
        M_k = M[k];
        for(l=0;l<chi2_nbSample;l++)
        {
            chi2_l = chi2[l];
            sqrt_chi2 = sqrt(chi2_l * inverse_freedomDegree);
            
            nbPaths_sample = nbPaths_by_sample;
            if(nbPaths_remaining>0)
            {
                nbPaths_sample += 1;
                nbPaths_remaining -= 1;
            }


            if(k==0)
            {
                w = NormalCum((M[k+1] + M[k])*0.5);
            }
            else if(k==M_nbSample-1)
            {
                w = (1.0 - NormalCum((M[k-1] + M[k])*0.5));
            }
            else
            {
                w = (NormalCum((M[k+1] + M[k])*0.5) - NormalCum((M[k-1] + M[k])*0.5));
            }

            if(l==0)
            {
                w *= Chi2Cum((chi2[l+1] + chi2[l])*0.5, freedomDegree)
                        /nbPaths_sample;
            }
            else if(l==chi2_nbSample-1)
            {
                w *= (1.0 - Chi2Cum((chi2[l-1] + chi2[l])*0.5, freedomDegree))
                        /nbPaths_sample;
            }
            else
            {
                w *= (  Chi2Cum((chi2[l+1] + chi2[l])*0.5, freedomDegree) -
                        Chi2Cum((chi2[l-1] + chi2[l])*0.5, freedomDegree) )
                        /nbPaths_sample;
            }



            
            for(j=0;j<nbPaths_sample;j++)
            {
                for(i=0;i<nbNames;i++)
                {
                    studentSequence[i+sample*nbNames] =
                        (beta[i]*M_k + sqrt_beta[i]*Z[i+j*nbNames])/
                                                            sqrt_chi2;
                }

                // WEIGHTS
                weight[sample] = w;
                totalWeight += weight[sample];
                sample +=1;            
            }
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
// StudentCopulatedIndicator
// this function returns the survival indicator and the weight of each path
//
*/
int StudentCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) [nbPaths*nbNames] 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    long freedomDegree,             /* (I) */
    long seed,
    const double *M,
    long M_nbSample,
    const double *chi2,
    long chi2_nbSample)
{
    static char routine[] = "StudentCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *studentSequence = NULL;
    double *T = NULL;

    studentSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(studentSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T == NULL) goto RETURN;

    for(i=0;i<nbNames;i++)
    {
        T[i] = StudentCumInverse(survivalProba[i], freedomDegree);
    }

    /* allocate the gaussian sequence */
    if(M==NULL||chi2==NULL||M_nbSample==0||chi2_nbSample==0)
    {
        status = StudentMultivariateDeviates_tp(studentSequence, weight, nbNames, nbPaths, beta, freedomDegree, seed);
        if (status==FAILURE) goto RETURN;
    }
    else
    {
        status = StudentMultivariateDeviatesSampling_tp(studentSequence, weight, nbNames, nbPaths, beta, freedomDegree, seed, M, M_nbSample, chi2, chi2_nbSample);
        if (status==FAILURE) goto RETURN;
    }

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(studentSequence[i+j*nbNames] > T[i])
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
    if(T) free(T);
    if(studentSequence) free(studentSequence);
    if(status == FAILURE)
    {
         DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);   
    }
    return status;
}


/* ---------------------------------------------------------------------------
// StudentCopulatedIndicator_mtp
// this function returns the survival indicator and the weight of each path
//
*/
int StudentCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) [nbPaths*nbNames] 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    long freedomDegree,             /* (I) */
    long seed,
    const double *M,
    long M_nbSample,
    const double *chi2,
    long chi2_nbSample)
{
    static char routine[] = "StudentCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double *studentSequence = NULL;
    double *T = NULL;

    studentSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(studentSequence == NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T == NULL) goto RETURN;

    for(i=0;i<nbNames*nbTimes;i++)
    {
        T[i] = StudentCumInverse(survivalProba[i], freedomDegree);
    }

    /* allocate the gaussian sequence */
    if(M==NULL||chi2==NULL||M_nbSample==0||chi2_nbSample==0)
    {
        status = StudentMultivariateDeviates_tp(studentSequence, weight, nbNames, nbPaths, beta, freedomDegree, seed);
        if (status==FAILURE) goto RETURN;
    }
    else
    {
        status = StudentMultivariateDeviatesSampling_tp(studentSequence, weight, nbNames, nbPaths, beta, freedomDegree, seed, M, M_nbSample, chi2, chi2_nbSample);
        if (status==FAILURE) goto RETURN;
    }

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
            for(k=0;k<nbTimes;k++)
            {
                if(studentSequence[i+j*nbNames] > T[k+i*nbTimes])
                {
                    copulatedSurvivalIndicator[i+j*nbNames] = k;
                    break;
                }
            }
        }
    }

    status = SUCCESS;
RETURN:
    if(T) free(T);
    if(studentSequence) free(studentSequence);
    if(status == FAILURE)
    {
         DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);   
    }
    return status;
}

