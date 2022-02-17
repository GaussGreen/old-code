#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "proba_utils.h"
#include "error2.h"
#include "payoff_r.h"
#include "gaussian.h"
#include "student.h"
#include "dependence.h"
#include "independence.h"
#include "gumbel.h"
#include "clayton.h"

typedef double(*PAYOFF_FUNCTION)(CREDIT_PORTFOLIO_R *, int *, double *);


/* ---------------------------------------------------------------------------
// CreditPortfolioRCreate
//
*/
CREDIT_PORTFOLIO_R CreditPortfolioRCreate(  double *notional,
                                            double *recovery,
                                            double *sigma_recovery,
                                            double *nameMaturity,
                                            long nbNames,
                                            double strike1,
                                            double strike2)
{
    int status = FAILURE;
    static char routine[] = "CreditPortfolioRCreate";
    CREDIT_PORTFOLIO_R cp;
    memset(&cp,0,sizeof(CREDIT_PORTFOLIO_R));

    // malloc
    cp.notional = malloc(nbNames*sizeof(double));
    if(cp.notional == NULL) goto RETURN;
    cp.recovery = malloc(nbNames*sizeof(double));
    if(cp.recovery == NULL) goto RETURN;
    cp.sigma_recovery = malloc(nbNames*sizeof(double));
    if(cp.sigma_recovery == NULL) goto RETURN;
    cp.nameMaturity = malloc(nbNames*sizeof(double));
    if(cp.nameMaturity == NULL) goto RETURN;

    // memcpy
    memcpy(cp.notional, notional, nbNames*sizeof(double));
    memcpy(cp.nameMaturity, nameMaturity, nbNames*sizeof(double));
    memcpy(cp.recovery, recovery, nbNames*sizeof(double));
    memcpy(cp.sigma_recovery, sigma_recovery, nbNames*sizeof(double));
    cp.strike1 = strike1;
    cp.strike2 = strike2;
    cp.nbNames = nbNames;
    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CreditPortfolioRFree(&cp);
    }
    return cp;
}

/* ---------------------------------------------------------------------------
// CreditPortfolioRFree
//
*/
void CreditPortfolioRFree(CREDIT_PORTFOLIO_R *cp)
{
    if(cp->nameMaturity) free(cp->nameMaturity);
    if(cp->notional) free (cp->notional);
    if(cp->recovery) free(cp->recovery);
    if(cp->sigma_recovery) free(cp->sigma_recovery);
}


/* ---------------------------------------------------------------------------
// IndicatorSimCreate
//
*/
INDICATOR_RECOVERY_SIM IndicatorRecoverySimCreate(   long nbPaths,
                                    long nbNames,
                                    const double *survivalProba,
                                    const COPULA_R *copula,
                                    const double *M,
                                    long M_nbSample)
{
    static char routine[] = "IndicatorSimCreate";
    int status = FAILURE;
    GAUSSIAN_PARAM *gaussian_param;
    double *beta = NULL;
    double *alpha;
    long seed;
    INDICATOR_RECOVERY_SIM is;
    memset(&is,0,sizeof(INDICATOR_RECOVERY_SIM));
    is.nbNames = nbNames;
    is.nbPaths = nbPaths;

    is.weight = malloc(nbPaths*sizeof(double));
    if(is.weight==NULL) goto RETURN;
    is.indicator = malloc(nbPaths*nbNames*sizeof(int));
    if(is.indicator==NULL) goto RETURN;
    is.X = malloc(nbPaths*nbNames*sizeof(double));
    if(is.X==NULL) goto RETURN;

    switch(copula->type)
    {
    case CR_GAUSS: /* Gaussian */
        gaussian_param  = (GAUSSIAN_PARAM *) copula->param;
        beta = gaussian_param->beta;
        alpha = gaussian_param->alpha;
        seed = gaussian_param->seed;

        status = GaussianCopulatedIndicatorRecovery( is.indicator,
                                    is.weight,
                                    is.X,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    alpha,
                                    seed,
                                    M,
                                    M_nbSample);
        if(status == FAILURE) goto RETURN;
        break;

    default:
        goto RETURN;
    }
    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        IndicatorRecoverySimFree(&is);
    }
    return is;
}

/* ---------------------------------------------------------------------------
// IndicatorRecoverySimFree
//
*/
void IndicatorRecoverySimFree(INDICATOR_RECOVERY_SIM *is)
{
    if (!is) return;
    if(is->indicator) free(is->indicator);
    if(is->weight) free(is->weight);
    if(is->X) free(is->X);
    is->weight = NULL;
    is->indicator = NULL;
    is->X = NULL;
}


/* ---------------------------------------------------------------------------
// CopulaRCreate
//
*/
COPULA_R CopulaRCreate( long type,
                        long nbNames,
                        const double *beta,
                        const double *alpha,
                        int *status)
{
    static char routine[] = "CopulaRCreate";
    GAUSSIAN_PARAM *gaussian_param = NULL;
    double *beta_param = NULL;
    double *alpha_param = NULL;
    COPULA_R c;
    memset(&c,0,sizeof(COPULA_R));
    *status = FAILURE;

    c.type = type;
    switch(type){
    case CR_GAUSS: /* Gaussian */
        gaussian_param = malloc(sizeof(GAUSSIAN_PARAM));
        memset(gaussian_param,0,sizeof(GAUSSIAN_PARAM));
        c.param = gaussian_param;
        if(gaussian_param==NULL) goto RETURN;

        beta_param = (double *)malloc(nbNames*sizeof(double));
        gaussian_param->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, beta, nbNames*sizeof(double));

        alpha_param = (double *)malloc(nbNames*sizeof(double));
        gaussian_param->alpha = alpha_param;
        if(alpha_param == NULL) goto RETURN;
        memcpy(alpha_param, alpha, nbNames*sizeof(double));
        break;

    default:
        goto RETURN;
    }
    *status = SUCCESS;
RETURN:
    if(*status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CopulaRFree(&c);                
    }
    return c;
}


/* ---------------------------------------------------------------------------
// CopulaRFree
//
*/
void CopulaRFree(COPULA_R *copula)
{
    if (!copula) return;
    switch(copula->type) {
    case CR_GAUSS: {
        GAUSSIAN_PARAM *gaussian_param = (GAUSSIAN_PARAM*) copula->param;
        if(gaussian_param->beta) free(gaussian_param->beta);
        gaussian_param->beta = NULL;
        if(gaussian_param->alpha) free(gaussian_param->alpha);
        gaussian_param->alpha = NULL;
        if(gaussian_param) free(gaussian_param);
        copula->param = NULL;
        break; }

    default:
        //assert(copula->param == NULL)
        ;
    }
}


/* --------------------------------------------------------------------------
// TranchePayoff_R
// returns the payoff of a tranche at a single time point
// given indicators
*/
double TranchePayoff_R(
    CREDIT_PORTFOLIO_R *port,/* (I) portfolio spec */
    int *indicator,        /* (I) payoff at steps */
    double *X)
{
    double totalLoss = 0.0;
    double trancheLoss = 0.0;
    long i = 0;
    double totalNotional = 0.0;
    double recovery = 0.0;
    for(i=0;i<port->nbNames;i++)
    {
        totalNotional += port->notional[i];
    }

    for(i=0;i<port->nbNames;i++)
    {
        if(indicator[i] == 0)
        {
            recovery = max(min(
                          port->recovery[i]+port->sigma_recovery[i]*X[i],
                              1.0),0.0);
            totalLoss += port->notional[i]*(1-recovery);
        }
    }

    trancheLoss =   min(
                    max(    totalLoss - port->strike1*totalNotional, 0),
                        (port->strike2 - port->strike1)*totalNotional);
    return trancheLoss;
}

/* ---------------------------------------------------------------------------
// ExpectedPayoff_R
// this function returns the expected trancheLoss
// at a single time point given indicators of default
// and weights for each path
*/
double ExpectedPayoff_R(   CREDIT_PORTFOLIO_R *port,
                           INDICATOR_RECOVERY_SIM *is)
{
    static PAYOFF_FUNCTION payoffFunction = &TranchePayoff_R;
    long j;
    long nbNames = port->nbNames;
    double expectedLoss = 0.0;

    for(j=0;j<is->nbPaths;j++)
    {
        expectedLoss += payoffFunction(port,&(is->indicator[j*nbNames]),&is->X[j*nbNames])
                        * is->weight[j];    
    }

    return expectedLoss;
}

/* ---------------------------------------------------------------------------
// PayoffDistribution_R
// this function returns the trancheLoss distribution
// at a single time point given indicators of default
// and weights for each path
*/
int PayoffDistribution_R(  CREDIT_PORTFOLIO_R *port,
                            INDICATOR_RECOVERY_SIM *is,
                            const double *sampleLoss,
                            long nbSampleLoss,
                            double *loss)
{
    static PAYOFF_FUNCTION payoffFunction = &TranchePayoff_R;
    int status = FAILURE;
    long i,j;
    long nbNames = port->nbNames;
    double l = 0.0;
    double totalNotional = 0.0;

    for(i=0;i<nbNames;i++)
    {
        totalNotional += port->notional[i];
    }

    for(i=0;i<nbSampleLoss;i++)
    {
        loss[i] = 0.0;
    }

    for(j=0;j<is->nbPaths;j++)
    {
        l = payoffFunction(port,&(is->indicator[j*nbNames]),&is->X[j*nbNames]);
        l /= totalNotional;

        if(l<=sampleLoss[0])
        {
            loss[0] += is->weight[j];
        }
 
        for(i=1;i<nbSampleLoss;i++)
        {
            if(l>sampleLoss[i-1] && l<=sampleLoss[i])
            {
                loss[i] += is->weight[j];
            }
        }
    }

    status = SUCCESS;
    return status;
}


/* ---------------------------------------------------------------------------
// NbDefaultNameDistribution_R
// this function returns the discrete probability
// of k name defaulted
// distribution[k] = proba(k name defaulted)
// k = 0 .. nbNames
*/
int NbDefaultNameDistribution_R(INDICATOR_RECOVERY_SIM i_s, double *distribution)
{
    static char routine[] = "NbDefaultNameDistribution_R";
    int status = FAILURE;
    long nbNames = i_s.nbNames;
    long nbPaths = i_s.nbPaths;
    long i = 0;
    long j = 0;
    long nbDefault = 0;
    double totalWeight = 0.0;
    for(j=0;j<nbPaths;j++)
    {
        nbDefault = 0;
        for(i=0;i<nbNames;i++)
        {
            if(i_s.indicator[i+j*nbNames] == 0)
            {
                nbDefault += 1;
            }
        }
        distribution[nbDefault] += i_s.weight[j];
        totalWeight += i_s.weight[j];
    }
    
    for(i=0;i<=nbNames;i++)
    {
        distribution[i] /= totalWeight;
    }
    status = SUCCESS;
//RETURN:
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}


