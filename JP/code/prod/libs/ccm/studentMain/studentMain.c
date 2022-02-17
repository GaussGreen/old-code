#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include "../payoff_tp.h"
#include "../student.h"
#include "../gaussian.h"
#include "../gumbel.h"
#include "../dependence.h"
#include "../independence.h"
#include "../proba_utils.h"
#include "../convolution.h"
#define MALLOC(TYPE,N) (TYPE*)malloc(sizeof(TYPE)*N)
#define FREE(P) {if (P) free(P);}

/*****************************************************************/
/* Test Convolution */
static void testConvolution(int copulaType, int hasDep)
{
    long nbNames = 100;
    long freedomDegree = 10;
    COPULA copula, dep, composite;
    double *beta =  NULL;
    double *survivalProba = NULL;
    int    *lossAmt       = NULL;
    double *density       = NULL;
    double *alpha         = NULL;
    int maxLoss;
    double qM,qZ;
    int i;
    int status;
    const double *dummy = 0;
    qM=  0.5;
    qZ= -0.1;

    alpha         = malloc(nbNames*sizeof(double));
    beta          = malloc(nbNames*sizeof(double));
    survivalProba = malloc(nbNames*sizeof(double));
    lossAmt       = malloc(nbNames*sizeof(int));

    maxLoss = 0;
    for(i=0;i<nbNames;i++)
    {
        lossAmt[i]       = (i * 97) % 5;
        survivalProba[i] = 0.8 + 0.1 * cos(i * 89325);
        beta[i]          = 0.6 + 0.1 * cos(i * 23457);
        alpha[i]         = 1. - log(0.9980 - 0.001 * (i%2)) 
                              / log(survivalProba[i]);
        maxLoss         += lossAmt[i];
    }
    density = malloc((maxLoss+1)*sizeof(double));
    copula = CopulaCreate(
        copulaType, nbNames, beta, freedomDegree, 
        0., 0, 0., 0., 0., 0.,
        qM,qZ,
        NULL,0L, NULL,0L, &status);
    dep = CopulaCreate(
        2, nbNames, beta, freedomDegree, /* dependence type */
        0., 0, 0., 0., 0., 0.,
        qM,qZ,
        NULL,0L, NULL,0L, &status);
    composite = CopulaProductCreate(&copula,&dep,alpha,nbNames);
    
    if (!calcLossDensity(
        lossAmt,nbNames,maxLoss,
        !hasDep ? &copula : &composite,
        survivalProba,density))
        printf("density calc with type=%d, hasDep=%d succeeded\n",copulaType, hasDep);

    FREE(alpha);
    FREE(beta);
    FREE(survivalProba);
    FREE(density);
    FREE(lossAmt);
    CopulaFree(&copula);
    CopulaFree(&dep);
    CopulaFree(&composite);
}

/*****************************************************************/
/* main */
void main()
{
    int copulaType = 4;
    long nbPaths = 10000;
    long nbNames = 2;
    long freedomDegree = 10;

    int i = 0;
    int j = 0;
    int status = FAILURE;

    CREDIT_PORTFOLIO creditPortfolio;
    INDICATOR_SIM indicatorSim;
    COPULA copula;
    double strike1 = 0.;
    double strike2 = 100.;
    /*allocation of the Sequences */
    double *beta =  NULL;
    double *survivalProba = NULL;
    double *recovery = NULL;
    double *notional = NULL;
    double *nameMaturity = NULL;
    double *prob = NULL;
    double expectedPayoff;
    double theta = 2;
    double qM,qZ;
    double *M_sample = NULL;
    double *chi2_sample = NULL;
    long M_nbSample = 0;
    long chi2_nbSample = 0;
    long seed = -7;

    qM= 0.0;
    qZ= 0.0;

    testConvolution(0,0);
    testConvolution(7,0);
    testConvolution(0,1);
    testConvolution(7,1);
    return;
    memset(&indicatorSim,0,sizeof(INDICATOR_SIM));
    memset(&creditPortfolio,0,sizeof(CREDIT_PORTFOLIO));
    beta =  malloc(nbNames * sizeof(double));
    survivalProba = malloc(nbNames*sizeof(double));
    recovery = malloc(nbNames*sizeof(double));
    notional = malloc(nbNames*sizeof(double));
    nameMaturity = malloc(nbNames*sizeof(double));

    prob = malloc(nbNames*sizeof(double));

    /* fills beta */
    for(i=0;i<nbNames;i++)
    {
        recovery[i] = 0.5;
        notional[i] = 10000.;
        nameMaturity[i] = 100.;
        beta[i] = 0.6;
        survivalProba[i] = 0.9;
    }

    creditPortfolio = CreditPortfolioCreate(notional,
                                            recovery,
                                            nameMaturity,
                                            nbNames,
                                            strike1,
                                            strike2);

    //copula = CopulaCreate(copulaType,nbNames,beta,freedomDegree, theta, seed,qM,qZ,M_sample,M_nbSample, chi2_sample, chi2_nbSample, &status);

    indicatorSim = IndicatorSimCreate(nbPaths,nbNames,survivalProba,&copula);
    expectedPayoff = ExpectedPayoff_tp(&creditPortfolio, &indicatorSim);

    printf("tranche loss: %lf\n", expectedPayoff);

    /*for(i=0;i<300;i++)
    {
        printf("x=%lf \t %lf \t",i/10.0, Chi2Cum(i/10.,10));
    }
    */
    for(i=1;i<11;i++)
    {
        for(j=1;j<11;j++)
        {
            prob[0] = i/10.;
            prob[1] = j/10.;
            printf("%lf ",prob[0]);
            printf("%lf ",prob[1]);
            status =  GumbelTest(nbNames,                   /* (I) */
                    nbPaths,                   /* (I) */
                    theta,
                    prob);
        }
    }

    if(beta) free(beta);
    if(survivalProba) free(survivalProba);
    CreditPortfolioFree(&creditPortfolio);
    IndicatorSimFree(&indicatorSim);
    CopulaFree(&copula);
}