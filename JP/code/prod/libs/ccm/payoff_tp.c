#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "error2.h"
#include "payoff_tp.h"
#include "binorm.c"
#include "gaussian.h"
#include "student.h"
#include "dependence.h"
#include "independence.h"
#include "gumbel.h"
#include "clayton.h"
#include "skew.h"
#include "delta.h"
#include "rootbrent.h"

typedef double(*PAYOFF_FUNCTION)(CREDIT_PORTFOLIO *, int *);

/* ---------------------------------------------------------------------------
// IndicatorSimCreate
//
*/
INDICATOR_SIM IndicatorSimCreate(   long nbPaths,
                                    long nbNames,
                                    const double *survivalProba,
                                    const COPULA *copula)
{
    static char routine[] = "IndicatorSimCreate";
    int status = FAILURE;
    long seed;
    long *ad_seed;
    INDICATOR_SIM is;
    memset(&is,0,sizeof(INDICATOR_SIM));
    is.nbNames = nbNames;
    is.nbPaths = nbPaths;

    is.weight = malloc(nbPaths*sizeof(double));
    if(is.weight==NULL) goto RETURN;
    is.indicator = malloc(nbPaths*nbNames*sizeof(int));
    if(is.indicator==NULL) goto RETURN;

    switch(copula->type)
    {
    case C_GAUSS: {/* Gaussian */ 
        GAUSSIAN_PARAM *param = (GAUSSIAN_PARAM *) copula->param;
        double *beta = NULL;
        long M_nbSample = param->M_nbSample;
        double *M_sample = param->M_sample;
        beta = param->beta;
        seed = param->seed;
        status = GaussianCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    seed,
                                    M_sample,
                                    M_nbSample);
        if(status == FAILURE) goto RETURN;
        break;}
    case C_STUDENT: {    /* Student */
        long freedomDegree;
        STUDENT_PARAM *param = (STUDENT_PARAM *) copula->param;
        double *beta = param->beta;
        long M_nbSample = param->M_nbSample;
        double *M_sample = param->M_sample;
        long chi2_nbSample = param->chi2_nbSample;
        double *chi2_sample = param->chi2_sample;
        seed = param->seed;
        freedomDegree = (long)param->freedomDegree;
        status = StudentCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    freedomDegree,
                                    seed,
                                    M_sample,
                                    M_nbSample,
                                    chi2_sample,
                                    chi2_nbSample);
        if(status == FAILURE) goto RETURN;

        break; }
    case C_DEPENDENCE:    /* Dependence */
        ad_seed = (long *)copula->param;
        status = DependenceCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,
                                    *ad_seed); 
        if(status == FAILURE) goto RETURN;        
        break;
    case C_INDEPENDENCE:                     /* Independence */
        ad_seed = (long *)copula->param;
        status = IndependenceCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,
                                    *ad_seed);
        if(status == FAILURE) goto RETURN;        
        break;
    case C_GUMBEL: {
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)copula->param;  
        double theta = param->theta;
        seed = param->seed;
        status = GumbelCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,
                                    theta,
                                    seed
                                    );
        if(status == FAILURE) goto RETURN;        
        break; }
    case C_CLAYTON: {
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)copula->param;  
        double theta = param->theta;
        seed = param->seed;
        status = ClaytonCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,
                                    theta,
                                    seed
                                    );
        if(status == FAILURE) goto RETURN;        
        break; }
    case C_COMPOSITE : {
        INDICATOR_SIM is2;
        PRODUCT_PARAM *pp = (PRODUCT_PARAM *) copula->param;
        is2 = IndicatorProductSimCreate(nbPaths,nbNames,survivalProba,pp->c_1,pp->c_2,pp->alpha);
        memcpy(is.indicator,is2.indicator,nbNames*nbPaths*sizeof(int));
        memcpy(is.weight,is2.weight,nbPaths*sizeof(double));
        IndicatorSimFree(&is2);
        break; }
    case C_SKEW:{
        double qM,qZ;
        Q_PARAM *param = (Q_PARAM *) copula->param;
        double *beta = param->beta;
        qM = param->qM;
        qZ = param->qZ;
        seed = param->seed;
        status = QCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    qM,
                                    qZ,
                                    10000,
                                    1e-7,
                                    seed);
        if(status == FAILURE) goto RETURN;

        break; }
    case C_DELTA:{
        double xa,xb,xc,delta;
        DELTA_PARAM *param = (DELTA_PARAM *) copula->param;
        double *beta = param->beta;
        xa = param->xa;
        xb = param->xb;
        xc = param->xc;
        delta = param->delta;
        seed = param->seed;
        status = DeltaCopulatedIndicator( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    xa,
                                    xb,
                                    xc,
                                    delta,
                                    seed);
        if(status == FAILURE) goto RETURN;

        break; }
    default:
        goto RETURN;
    }
    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        IndicatorSimFree(&is);
    }
    return is;
}


/* ---------------------------------------------------------------------------
// IndicatorSimCreate_mtp
//
*/
INDICATOR_SIM IndicatorSimCreate_mtp(   long nbPaths,
                                        long nbNames,
                                        long nbTimes,
                                        const double *survivalProba,
                                        const COPULA *copula)
{
    static char routine[] = "IndicatorSimCreate_mtp";
    int status = FAILURE;
    long seed;
    long *ad_seed;
    INDICATOR_SIM is;
    memset(&is,0,sizeof(INDICATOR_SIM));
    is.nbNames = nbNames;
    is.nbPaths = nbPaths;

    is.weight = malloc(nbPaths*sizeof(double));
    if(is.weight==NULL) goto RETURN;
    is.indicator = malloc(nbPaths*nbNames*sizeof(int));
    if(is.indicator==NULL) goto RETURN;

    switch(copula->type)
    {
    case C_GAUSS: {/* Gaussian */ 
        GAUSSIAN_PARAM *param = (GAUSSIAN_PARAM *) copula->param;
        double *beta = NULL;
        long M_nbSample = param->M_nbSample;
        double *M_sample = param->M_sample;
        beta = param->beta;
        seed = param->seed;
        status = GaussianCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    seed,
                                    M_sample,
                                    M_nbSample);
        if(status == FAILURE) goto RETURN;
        break;}
    case C_STUDENT: {    /* Student */
        long freedomDegree;
        STUDENT_PARAM *param = (STUDENT_PARAM *) copula->param;
        double *beta = param->beta;
        long M_nbSample = param->M_nbSample;
        double *M_sample = param->M_sample;
        long chi2_nbSample = param->chi2_nbSample;
        double *chi2_sample = param->chi2_sample;
        seed = param->seed;
        freedomDegree = (long)param->freedomDegree;
        status = StudentCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    freedomDegree,
                                    seed,
                                    M_sample,
                                    M_nbSample,
                                    chi2_sample,
                                    chi2_nbSample);
        if(status == FAILURE) goto RETURN;

        break; }
    case C_DEPENDENCE:    /* Dependence */
        ad_seed = (long *)copula->param;
        status = DependenceCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,
                                    *ad_seed); 
        if(status == FAILURE) goto RETURN;        
        break;
    case C_INDEPENDENCE:                     /* Independence */
        ad_seed = (long *)copula->param;
        status = IndependenceCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,
                                    *ad_seed);
        if(status == FAILURE) goto RETURN;        
        break;
    case C_GUMBEL: {
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)copula->param;  
        double theta = param->theta;
        seed = param->seed;
        status = GumbelCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,
                                    theta,
                                    seed
                                    );
        if(status == FAILURE) goto RETURN;        
        break; }
    case C_CLAYTON: {
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)copula->param;  
        double theta = param->theta;
        seed = param->seed;
        status = ClaytonCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,
                                    theta,
                                    seed
                                    );
        if(status == FAILURE) goto RETURN;        
        break; }
    case C_COMPOSITE : {
        INDICATOR_SIM is2;
        PRODUCT_PARAM *pp = (PRODUCT_PARAM *) copula->param;
        is2 = IndicatorProductSimCreate_mtp(nbPaths,nbNames,nbTimes,survivalProba,pp->c_1,pp->c_2,pp->alpha);
        memcpy(is.indicator,is2.indicator,nbNames*nbPaths*sizeof(int));
        memcpy(is.weight,is2.weight,nbPaths*sizeof(double));
        IndicatorSimFree(&is2);
        break; }
    case C_SKEW:{
        double qM,qZ;
        Q_PARAM *param = (Q_PARAM *) copula->param;
        double *beta = param->beta;
        qM = param->qM;
        qZ = param->qZ;
        seed = param->seed;
        status = QCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    qM,
                                    qZ,
                                    100,
                                    10e-5,
                                    seed);
        if(status == FAILURE) goto RETURN;
        break; }
     case C_DELTA:{
        double xa,xb,xc,delta;
        DELTA_PARAM *param = (DELTA_PARAM *) copula->param;
        double *beta = param->beta;
        xa = param->xa;
        xb = param->xb;
        xc = param->xc;
        delta = param->delta;
        seed = param->seed;
        status = DeltaCopulatedIndicator_mtp( is.indicator,
                                    is.weight,
                                    survivalProba,
                                    nbTimes,
                                    nbNames,                  
                                    nbPaths,                   
                                    beta,
                                    xa,
                                    xb,
                                    xc,
                                    delta,
                                    seed);
        if(status == FAILURE) goto RETURN;
        break; }
    default:
        goto RETURN;
    }
    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        IndicatorSimFree(&is);
    }
    return is;    
}


/* ---------------------------------------------------------------------------
// IndicatorSimFree
//
*/
void IndicatorSimFree(INDICATOR_SIM *is)
{
    if (!is) return;
    if(is->indicator) free(is->indicator);
    if(is->weight) free(is->weight);
    is->weight = NULL;
    is->indicator = NULL;
}


/* ---------------------------------------------------------------------------
// CopulaCreate
//
*/
COPULA CopulaCreate(    CopulaType type,
                        long nbNames,
                        const double *beta,
                        long freedomDegree,
                        double theta,
                        long seed,
                        double qM,
                        double qZ,
                        double xa,
                        double xb,
                        double xc,
                        double delta,
                        const double *M_sample,
                        long M_nbSample,
                        const double  *chi2_sample,
                        long chi2_nbSample,
                        int *status)
{
    static char routine[] = "CopulaCreate";
    double *beta_param = NULL;
    double *theta_param = NULL;
    double *M = NULL;
    double *chi2 = NULL;
    COPULA c;
    memset(&c,0,sizeof(COPULA));
    *status = FAILURE;

    c.type = type;
    switch(type){
    case C_GAUSS: {/* Gaussian */
        GAUSSIAN_PARAM *gaussian_param = (GAUSSIAN_PARAM *)malloc(sizeof(GAUSSIAN_PARAM));
        if(gaussian_param==NULL) goto RETURN;
        memset(gaussian_param,0,sizeof(GAUSSIAN_PARAM));
        c.param = gaussian_param;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        gaussian_param->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, beta, nbNames*sizeof(double));
        M = (double *)malloc(M_nbSample*sizeof(double));
        gaussian_param->M_sample = M;
        memcpy(M,M_sample,M_nbSample*sizeof(double));
        gaussian_param->M_nbSample = M_nbSample;
        gaussian_param->seed = seed;
        break; }
    case C_STUDENT: {/* Student */
        STUDENT_PARAM *student_param = (STUDENT_PARAM *)malloc(sizeof(STUDENT_PARAM));
        if(student_param==NULL) goto RETURN;
        memset(student_param,0,sizeof(STUDENT_PARAM));
        c.param = student_param;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        student_param->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, beta, nbNames*sizeof(double));
        M = (double *)malloc(M_nbSample*sizeof(double));
        student_param->M_sample = M;
        memcpy(M,M_sample,M_nbSample*sizeof(double));
        chi2 = (double *)malloc(chi2_nbSample*sizeof(double));
        student_param->chi2_sample = chi2;
        memcpy(chi2,chi2_sample,chi2_nbSample*sizeof(double));
        student_param->M_nbSample = M_nbSample;
        student_param->chi2_nbSample = chi2_nbSample;
        student_param->freedomDegree = freedomDegree;
        student_param->seed = seed;
        break; }
    case C_DEPENDENCE: {/* dependence */
        long *param = malloc(sizeof(long));
        if(param==NULL) goto RETURN;
        *param = seed;
        c.param = param;
        break; }
    case C_INDEPENDENCE: {/* independence */
        long *param = malloc(sizeof(long));
        if(param==NULL) goto RETURN;
        *param = seed;
        c.param = param;
        break; }
    case C_GUMBEL: {/* Gumbel */
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)malloc(sizeof(ARCHIMEDEAN_PARAM));
        if(param==NULL) goto RETURN;
        memset(param,0,sizeof(ARCHIMEDEAN_PARAM));
        c.param = param;
        if(theta<=1) DR_Error("theta must be > 1 for Gumbel");
        theta_param = malloc(sizeof(double));
        if(theta_param == NULL) goto RETURN;
        c.param = theta_param;
        *theta_param = theta;
        break; }
    case C_CLAYTON: {/* Clayton */
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *)malloc(sizeof(ARCHIMEDEAN_PARAM));
        if(param==NULL) goto RETURN;
        memset(param,0,sizeof(ARCHIMEDEAN_PARAM));
        c.param = param;
        if(theta<0) DR_Error("theta must be >= 0 for Clayton");
        theta_param = malloc(sizeof(double));
        c.param = theta_param;
        if(theta_param == NULL) goto RETURN;
        *theta_param = theta;
        break; }
    case C_SKEW: {/* Q */      
        Q_PARAM *q_param = (Q_PARAM *)malloc(sizeof(Q_PARAM));
        if(q_param==NULL) goto RETURN;
        memset(q_param,0,sizeof(Q_PARAM));
        c.param = q_param;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        q_param->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, beta, nbNames*sizeof(double));
        q_param->qM = qM;
        q_param->qZ = qZ;
        q_param->seed = seed;
        break; }
    case C_DELTA: {/* delta */
        DELTA_PARAM *delta_param = (DELTA_PARAM *)malloc(sizeof(DELTA_PARAM));
        if(delta_param==NULL) goto RETURN;
        memset(delta_param,0,sizeof(DELTA_PARAM));
        c.param = delta_param;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        delta_param->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, beta, nbNames*sizeof(double));
        delta_param->xa = xa;
        delta_param->xb = xb;
        delta_param->xc = xc;
        delta_param->delta = delta;
        delta_param->seed = seed;
        break; }
    default:
        goto RETURN;
    }
    *status = SUCCESS;
RETURN:
    if(*status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CopulaFree(&c);                
    }
    return c;
}

/* --------------------------------------------------------------------------
   CopulaProductCreate
   To create a product of copulas, the two copula parameters of this function need
   to be independent. Hence these copula need to be generated with different independent
   random sequences, i.e. different seeds. Using the same random sequences for the two copula
   will make the expected tranche loss lower than its theoretical value by introducing dependence
   between the copulas.
*/
COPULA CopulaProductCreate(const COPULA *c_1, const COPULA *c_2, const double *alpha, long nbNames)
{
    static char routine[] = "CopulaProductCreate";
    int status = FAILURE;
    COPULA c;
    PRODUCT_PARAM *param = NULL;
    memset(&c,0,sizeof(COPULA));
    c.type = C_COMPOSITE;

    param = malloc(sizeof(PRODUCT_PARAM));
    if(param == NULL) goto RETURN;
    memset(param,0,sizeof(PRODUCT_PARAM));
    c.param = param;

    param->c_1 = malloc(sizeof(COPULA));
    if(param->c_1==NULL) goto RETURN;
    memset(param->c_1,0,sizeof(COPULA));
    status = CopulaCopy(c_1,param->c_1,nbNames);
    if(status == FAILURE) goto RETURN;
    
    param->c_2 = malloc(sizeof(COPULA));
    if(param->c_2==NULL) goto RETURN;
    memset(param->c_2,0,sizeof(COPULA));
    status = CopulaCopy(c_2,param->c_2,nbNames);
    if(status == FAILURE) goto RETURN;

    param->alpha = malloc(nbNames*sizeof(double));
    if(param->alpha==NULL) goto RETURN;
    memcpy(param->alpha,alpha,nbNames*sizeof(double));


    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CopulaFree(&c);
    }
    return c;
}

COPULA Copula3ProductCreate(const COPULA *c_1,
                            const COPULA *c_2,
                            const COPULA *c_3,
                            const double *alpha1,
                            const double *alpha2,
                            long nbNames)
{

    static char routine[] = "Copula3ProductCreate";
    int status = FAILURE;
    int i = 0;
    COPULA c1;
    COPULA c;
    PRODUCT_PARAM *param1 = NULL;
    PRODUCT_PARAM *param = NULL;

    memset(&c1,0,sizeof(COPULA));
    memset(&c,0,sizeof(COPULA));
    c1.type = C_COMPOSITE;

    /* first product between c_2 and c_3 */
    param1 = malloc(sizeof(PRODUCT_PARAM));
    if(param1 == NULL) goto RETURN;
    memset(param1,0,sizeof(PRODUCT_PARAM));
    c1.param = param1;

    param1->c_1 = malloc(sizeof(COPULA));
    if(param1->c_1==NULL) goto RETURN;
    memset(param1->c_1,0,sizeof(COPULA));
    status = CopulaCopy(c_2,param1->c_1,nbNames);
    if(status == FAILURE) goto RETURN;
    
    param1->c_2 = malloc(sizeof(COPULA));
    if(param1->c_2==NULL) goto RETURN;
    memset(param1->c_2,0,sizeof(COPULA));
    status = CopulaCopy(c_3,param1->c_2,nbNames);
    if(status == FAILURE) goto RETURN;

    param1->alpha = malloc(nbNames*sizeof(double));
    if(param1->alpha==NULL) goto RETURN;
    for(i=0;i<nbNames;i++)
    {
        param1->alpha[i] = 1 - alpha2[i]/(1-alpha1[i]);
    }

    /* product with c_1*/
    c.type = C_COMPOSITE;
    param = malloc(sizeof(PRODUCT_PARAM));
    if(param == NULL) goto RETURN;
    memset(param,0,sizeof(PRODUCT_PARAM));
    c.param = param;

    param->c_2 = malloc(sizeof(COPULA));
    if(param->c_2==NULL) goto RETURN;
    memset(param->c_2,0,sizeof(COPULA));
    status = CopulaCopy(&c1,param->c_2,nbNames);
    if(status == FAILURE) goto RETURN;
    
    param->c_1 = malloc(sizeof(COPULA));
    if(param->c_1==NULL) goto RETURN;
    memset(param->c_1,0,sizeof(COPULA));
    status = CopulaCopy(c_1,param->c_1,nbNames);
    if(status == FAILURE) goto RETURN;

    param->alpha = malloc(nbNames*sizeof(double));
    if(param->alpha==NULL) goto RETURN;
    for(i=0;i<nbNames;i++)
    {
        param->alpha[i] = alpha1[i];
    }

    status = SUCCESS;
RETURN:
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CopulaFree(&c);
    }
    return c;
}

/* ---------------------------------------------------------------------------
// CopulaCopy
// deep copy c_1 in c_2
*/
int CopulaCopy(const COPULA *c_1, COPULA *c_2, long nbNames)
{
    static char routine[] = "CopulaCopy";
    int status = FAILURE;
    memset(c_2,0,sizeof(COPULA));

    c_2->type = c_1->type;
    switch(c_2->type){
    case C_GAUSS: /* Gaussian */ {
        double *beta_param2 = NULL;
        double *beta_param1 = NULL;
        double *M_sample2 = NULL;
        GAUSSIAN_PARAM *param1 = c_1->param;
        GAUSSIAN_PARAM *param2 = (GAUSSIAN_PARAM *)malloc(sizeof(GAUSSIAN_PARAM));
        if(param2==NULL) goto RETURN;
        memset(param2,0,sizeof(GAUSSIAN_PARAM));
        c_2->param = param2;
        beta_param2 = (double *)malloc(nbNames*sizeof(double));
        param2->beta = beta_param2;
        if(beta_param2 == NULL) goto RETURN;
        beta_param1 = param1->beta;
        memcpy(beta_param2, beta_param1, nbNames*sizeof(double));
        M_sample2 = (double *)malloc((param1->M_nbSample)*sizeof(double));
        param2->M_sample = M_sample2;
        if(M_sample2 == NULL) goto RETURN;
        memcpy(param2->M_sample,param1->M_sample,(param1->M_nbSample)*sizeof(double));
        param2->M_nbSample = param1->M_nbSample;
        param2->seed = param1->seed;
        break; }
    case C_STUDENT: /* Student */ {
        double *beta_param = NULL;
        double *M_sample2 = NULL;
        double *chi2_sample2 = NULL;
        STUDENT_PARAM *student_param1 = c_1->param;
        STUDENT_PARAM *student_param2 = (STUDENT_PARAM *)malloc(sizeof(STUDENT_PARAM));
        if(student_param2 == NULL) goto RETURN;
        memset(student_param2,0,sizeof(STUDENT_PARAM));
        c_2->param = student_param2;
        student_param2->freedomDegree = student_param1->freedomDegree;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        student_param2->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, student_param1->beta, nbNames*sizeof(double));
        M_sample2 = (double *)malloc((student_param1->M_nbSample)*sizeof(double));
        student_param2->M_sample = M_sample2;
        if(M_sample2 == NULL) goto RETURN;
        memcpy(student_param2->M_sample,student_param1->M_sample,(student_param1->M_nbSample)*sizeof(double));
        chi2_sample2 = (double *)malloc((student_param1->chi2_nbSample)*sizeof(double));
        student_param2->chi2_sample = chi2_sample2;
        if(chi2_sample2 == NULL) goto RETURN;
        memcpy(student_param2->chi2_sample,student_param1->chi2_sample,(student_param1->chi2_nbSample)*sizeof(double));
        student_param2->M_nbSample = student_param1->M_nbSample;
        student_param2->chi2_nbSample = student_param1->chi2_nbSample;
        student_param2->seed = student_param1->seed;
        break; }
    case C_DEPENDENCE: {/* dependence */
        long *seed1 = c_1->param;
        long *seed2 = malloc(sizeof(long));
        c_2->param = seed2;
        if(seed2==NULL) goto RETURN;
        *seed2 = *seed1;
        break; }
    case C_INDEPENDENCE: {/* independence */
        long *seed1 = c_1->param;
        long *seed2 = malloc(sizeof(long));
        c_2->param = seed2;
        if(seed2==NULL) goto RETURN;
        *seed2 = *seed1;
        break; }
    case C_GUMBEL: /* Gumbel */ {
        ARCHIMEDEAN_PARAM *param1 = c_1->param;
        ARCHIMEDEAN_PARAM *param2 = (ARCHIMEDEAN_PARAM *)malloc(sizeof(ARCHIMEDEAN_PARAM));
        if(param2 == NULL) goto RETURN;
        memset(param2,0,sizeof(ARCHIMEDEAN_PARAM));
        c_2->param = param2;
        param2->theta = param1->theta;
        param2->seed = param1->seed;
        break; }
    case C_CLAYTON: /* Clayton */ {
        ARCHIMEDEAN_PARAM *param1 = c_1->param;
        ARCHIMEDEAN_PARAM *param2 = (ARCHIMEDEAN_PARAM *)malloc(sizeof(ARCHIMEDEAN_PARAM));
        if(param2 == NULL) goto RETURN;
        memset(param2,0,sizeof(ARCHIMEDEAN_PARAM));
        c_2->param = param2;
        param2->theta = param1->theta;
        param2->seed = param1->seed;
        break; }
    case C_COMPOSITE: /* Product */ {
        double *alpha = NULL;
        PRODUCT_PARAM *product_param1 = (PRODUCT_PARAM *)c_1->param;
        PRODUCT_PARAM *product_param2 = malloc(sizeof(PRODUCT_PARAM));
        if(product_param2==NULL) goto RETURN;
        memset(product_param2,0,sizeof(PRODUCT_PARAM));
        c_2->param = product_param2;
        alpha = malloc(nbNames*sizeof(double));
        product_param2->alpha = alpha;
        if(alpha==NULL) goto RETURN;
        product_param2->c_1 = malloc(sizeof(COPULA));
        if(product_param2->c_1==NULL) goto RETURN;
        memset(product_param2->c_1,0,sizeof(COPULA));
        product_param2->c_2 = malloc(sizeof(COPULA));
        if(product_param2->c_2==NULL) goto RETURN;
        memset(product_param2->c_2,0,sizeof(COPULA));

        memcpy(alpha,product_param1->alpha,nbNames*sizeof(double));
        status = CopulaCopy(product_param1->c_1,product_param2->c_1,nbNames);
        if(status == FAILURE) goto RETURN;
        status = CopulaCopy(product_param1->c_2,product_param2->c_2,nbNames);
        if(status == FAILURE) goto RETURN;
        break; }
    case C_SKEW: /* q */ {
        double *beta_param = NULL;
        Q_PARAM *q_param1 = c_1->param;
        Q_PARAM *q_param2 = (Q_PARAM *)malloc(sizeof(Q_PARAM));
        if(q_param2 == NULL) goto RETURN;
        memset(q_param2,0,sizeof(Q_PARAM));
        c_2->param = q_param2;
        q_param2->qM = q_param1->qM;
        q_param2->qZ = q_param1->qZ;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        q_param2->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, q_param1->beta, nbNames*sizeof(double));
        q_param2->seed = q_param1->seed;
        break; }
    case C_DELTA: /* delta */ {
        double *beta_param = NULL;
        DELTA_PARAM *delta_param1 = c_1->param;
        DELTA_PARAM *delta_param2 = (DELTA_PARAM *)malloc(sizeof(DELTA_PARAM));
        if(delta_param2 == NULL) goto RETURN;
        memset(delta_param2,0,sizeof(DELTA_PARAM));
        c_2->param = delta_param2;
        delta_param2->xa = delta_param1->xa;
        delta_param2->xb = delta_param1->xb;
        delta_param2->xc = delta_param1->xc;
        delta_param2->delta = delta_param1->delta;
        beta_param = (double *)malloc(nbNames*sizeof(double));
        delta_param2->beta = beta_param;
        if(beta_param == NULL) goto RETURN;
        memcpy(beta_param, delta_param1->beta, nbNames*sizeof(double));
        delta_param2->seed = delta_param1->seed;
        break; }
    default:
        goto RETURN;
    }
    status = SUCCESS;
RETURN:
    if(status ==  FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        CopulaFree(c_2);
    }
    return status;
}

/* ---------------------------------------------------------------------------
// CopulaFree
//
*/
void CopulaFree(COPULA *copula)
{
    if (!copula) return;
    switch(copula->type) {
    case C_GAUSS: {
        GAUSSIAN_PARAM *gp = (GAUSSIAN_PARAM *)copula->param;
        if(gp && gp->beta) free(gp->beta);
        if(gp && gp->M_sample) free(gp->M_sample);
        if(gp) free(gp);
        break; }
    case C_STUDENT: {
        STUDENT_PARAM *sp = (STUDENT_PARAM*)copula->param;
        if(sp && sp->beta) free(sp->beta);
        if(sp && sp->M_sample) free(sp->M_sample);
        if(sp && sp->chi2_sample) free(sp->chi2_sample);
        if (sp) free(sp);
        break; }
    case C_DEPENDENCE:
    case C_INDEPENDENCE: {
        long *seed = (long *) copula->param;
        if(seed) free(seed);
        copula->param = NULL;
        break; }
    case C_GUMBEL:
    case C_CLAYTON:{
        ARCHIMEDEAN_PARAM *param = (ARCHIMEDEAN_PARAM *) copula->param;
        if(param) free(param);
        break; }
    case C_COMPOSITE: {
        PRODUCT_PARAM *pp = (PRODUCT_PARAM*)copula->param;
        if(pp && pp->alpha) {free(pp->alpha); pp->alpha = NULL;}
        if(pp && pp->c_1) {CopulaFree(pp->c_1); free(pp->c_1); pp->c_1 = NULL;}
        if(pp && pp->c_2) {CopulaFree(pp->c_2); free(pp->c_2); pp->c_2 = NULL;}
        if(pp) free(pp);
        break; }
    case C_SKEW: {
        Q_PARAM *sp = (Q_PARAM*)copula->param;
        if(sp && sp->beta) free(sp->beta);
        if (sp) free(sp);
        break; }
    case C_DELTA: {
        DELTA_PARAM *sp = (DELTA_PARAM*)copula->param;
        if(sp && sp->beta) free(sp->beta);
        if (sp) free(sp);
        break; }
    case C_COPULA_FREE: break;
    default:
        assert(copula->param == NULL);
    }
    copula->param = NULL;
    copula->type  = C_COPULA_FREE;
}

/* --------------------------------------------------------------------------
// IndicatorProductSimCreate
//
*/
INDICATOR_SIM IndicatorProductSimCreate(long nbPaths,
                                        long nbNames,
                                        const double *survivalProba,
                                        const COPULA *c1,
                                        const COPULA *c2,
                                        const double *alpha)
{
    static char routine[] = "IndicatorProductSimCreate";
    int status = FAILURE;
    long i = 0;
    long j = 0;
    double totalWeight = 0.0;
    INDICATOR_SIM is1;
    INDICATOR_SIM is2;
    INDICATOR_SIM is;
    double *survivalProba1 = NULL;
    double *survivalProba2 = NULL;

    memset(&is,0,sizeof(INDICATOR_SIM));
    memset(&is1,0,sizeof(INDICATOR_SIM));
    memset(&is2,0,sizeof(INDICATOR_SIM));

    survivalProba1 = malloc(nbNames*sizeof(double));
    if(survivalProba1==NULL) goto RETURN;
    survivalProba2 = malloc(nbNames*sizeof(double));
    if(survivalProba2==NULL) goto RETURN;
    is.nbNames = nbNames;
    is.nbPaths = nbPaths;

    is.weight = malloc(nbPaths*sizeof(double));
    if(is.weight==NULL) goto RETURN;
    is.indicator = malloc(nbPaths*nbNames*sizeof(int));
    if(is.indicator==NULL) goto RETURN;


    for(i=0;i<nbNames;i++)
    {
        survivalProba1[i] = pow(survivalProba[i],alpha[i]);
        survivalProba2[i] = survivalProba[i] / survivalProba1[i];
    }

    if(c1->type==C_COMPOSITE)
    {
        PRODUCT_PARAM *pp1 = (PRODUCT_PARAM *) c1->param;
        is1 = IndicatorProductSimCreate(nbPaths,nbNames,survivalProba1,pp1->c_1,pp1->c_2,pp1->alpha);
    }
    else
    {
        is1 = IndicatorSimCreate(nbPaths,nbNames,survivalProba1,c1);
    }

    if(c2->type==C_COMPOSITE)
    {
        PRODUCT_PARAM *pp2 = (PRODUCT_PARAM *) c2->param;
        is2 = IndicatorProductSimCreate(nbPaths,nbNames,survivalProba2,pp2->c_1,pp2->c_2,pp2->alpha);
    }
    else
    {
        is2 = IndicatorSimCreate(nbPaths,nbNames,survivalProba2,c2);    
    }

    for(j=0;j<nbPaths;j++)
    {
        totalWeight += is1.weight[j]*is2.weight[j];
        is.weight[j] = is1.weight[j]*is2.weight[j];
        for(i=0;i<nbNames;i++)
        {
            is.indicator[i+j*nbNames] =   min( is1.indicator[i+j*nbNames],
                                               is2.indicator[i+j*nbNames]);
        }
    }

    for(j=0;j<nbPaths;j++)
    {
        is.weight[j] /= totalWeight;
    }

    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        IndicatorSimFree(&is);
    }
    IndicatorSimFree(&is1);
    IndicatorSimFree(&is2);

    if(survivalProba1) free(survivalProba1);
    if(survivalProba2) free(survivalProba2);
    return is;
}

/* --------------------------------------------------------------------------
// IndicatorProductSimCreate_mtp
//
*/
INDICATOR_SIM IndicatorProductSimCreate_mtp(long nbPaths,
                                        long nbNames,
                                        long nbTimes,
                                        const double *survivalProba,
                                        const COPULA *c1,
                                        const COPULA *c2,
                                        const double *alpha)
{
    static char routine[] = "IndicatorProductSimCreate_mtp";
    int status = FAILURE;
    long i = 0;
    long j = 0;
    double totalWeight = 0.0;
    INDICATOR_SIM is1;
    INDICATOR_SIM is2;
    INDICATOR_SIM is;
    double *survivalProba1 = NULL;
    double *survivalProba2 = NULL;

    memset(&is,0,sizeof(INDICATOR_SIM));
    memset(&is1,0,sizeof(INDICATOR_SIM));
    memset(&is2,0,sizeof(INDICATOR_SIM));

    survivalProba1 = malloc(nbNames*nbTimes*sizeof(double));
    if(survivalProba1==NULL) goto RETURN;
    survivalProba2 = malloc(nbNames*nbTimes*sizeof(double));
    if(survivalProba2==NULL) goto RETURN;
    is.nbNames = nbNames;
    is.nbPaths = nbPaths;

    is.weight = malloc(nbPaths*sizeof(double));
    if(is.weight==NULL) goto RETURN;
    is.indicator = malloc(nbPaths*nbNames*sizeof(int));
    if(is.indicator==NULL) goto RETURN;


    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbTimes;j++)
        {
            survivalProba1[i+j*nbNames] = pow(survivalProba[i+j*nbNames],alpha[i]);
            survivalProba2[i+j*nbNames] = survivalProba[i+j*nbNames] / survivalProba1[i+j*nbNames];
        }
    }

    if(c1->type==C_COMPOSITE)
    {
        PRODUCT_PARAM *pp1 = (PRODUCT_PARAM *) c1->param;
        is1 = IndicatorProductSimCreate_mtp(nbPaths,nbNames,nbTimes,survivalProba1,pp1->c_1,pp1->c_2,pp1->alpha);
    }
    else
    {
        is1 = IndicatorSimCreate_mtp(nbPaths,nbNames,nbTimes,survivalProba1,c1);
    }

    if(c2->type==C_COMPOSITE)
    {
        PRODUCT_PARAM *pp2 = (PRODUCT_PARAM *) c2->param;
        is2 = IndicatorProductSimCreate_mtp(nbPaths,nbNames,nbTimes,survivalProba2,pp2->c_1,pp2->c_2,pp2->alpha);
    }
    else
    {
        is2 = IndicatorSimCreate_mtp(nbPaths,nbNames,nbTimes,survivalProba2,c2);    
    }

    for(j=0;j<nbPaths;j++)
    {
        totalWeight += is1.weight[j]*is2.weight[j];
        is.weight[j] = is1.weight[j]*is2.weight[j];
        for(i=0;i<nbNames;i++)
        {
            is.indicator[i+j*nbNames] =   min( is1.indicator[i+j*nbNames],
                                               is2.indicator[i+j*nbNames]);
        }
    }

    for(j=0;j<nbPaths;j++)
    {
        is.weight[j] /= totalWeight;
    }

    status = SUCCESS;
RETURN:
    if(status==FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
        IndicatorSimFree(&is);
    }
    IndicatorSimFree(&is1);
    IndicatorSimFree(&is2);

    if(survivalProba1) free(survivalProba1);
    if(survivalProba2) free(survivalProba2);
    return is;
}


/* --------------------------------------------------------------------------
// TranchePayoff_tp
// returns the payoff of a tranche at a single time point
// given indicators
*/
double TranchePayoff_tp(
    CREDIT_PORTFOLIO *port,/* (I) portfolio spec */
    int *indicator)     /* (I) payoff at steps */
{
    double totalLoss = 0.0;
    double trancheLoss = 0.0;
    long i = 0;
    double totalNotional = 0.0;
    for(i=0;i<port->nbNames;i++)
    {
        totalNotional += port->notional[i];
    }

    for(i=0;i<port->nbNames;i++)
    {
        if(indicator[i] == 0)
        {
            totalLoss += port->notional[i]*(1-port->recovery[i]);
        }
    }

    trancheLoss =   min(
                    max(    totalLoss - port->strike1*totalNotional, 0),
                        (port->strike2 - port->strike1)*totalNotional);
    return trancheLoss;
}

/* ---------------------------------------------------------------------------
// ExpectedPayoff_tp
// this function returns the expected trancheLoss
// at a single time point given indicators of default
// and weights for each path
*/
double ExpectedPayoff_tp(   CREDIT_PORTFOLIO *port,
                            INDICATOR_SIM *is)
{
    static PAYOFF_FUNCTION payoffFunction = &TranchePayoff_tp;
    long j;
    long nbNames = port->nbNames;
    double expectedLoss = 0.0;

    for(j=0;j<is->nbPaths;j++)
    {
        expectedLoss += payoffFunction(port,&(is->indicator[j*nbNames]))
                        * is->weight[j];    
    }

    return expectedLoss;
}

/* ---------------------------------------------------------------------------
// PayoffDistribution_tp
// this function returns the trancheLoss distribution
// at a single time point given indicators of default
// and weights for each path
*/
int PayoffDistribution_tp(  CREDIT_PORTFOLIO *port,
                            INDICATOR_SIM *is,
                            const double *sampleLoss,
                            long nbSampleLoss,
                            double *loss)
{
    static PAYOFF_FUNCTION payoffFunction = &TranchePayoff_tp;
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
        l = payoffFunction(port,&(is->indicator[j*nbNames]));
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
// NbDefaultNameDistribution_tp
// this function returns the discrete probability
// of k name defaulted
// distribution[k] = proba(k name defaulted)
// k = 0 .. nbNames
*/
int NbDefaultNameDistribution(INDICATOR_SIM i_s, double *distribution)
{
    static char routine[] = "NbDefaultNameDistribution";
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


/* ---------------------------------------------------------------------------
// NbConditionalDefaultNameDistribution_mtp
// this function returns the discrete probability
// of k name defaulted
// distribution[k] = proba(k name defaulted)
// k = 0 .. nbNames
// [t0,t1] is the range where we observe the first default
// [t0,T] is the range where we observe the defaults conditional a default
// between t0 and t1
// idx_to, idx_t1, idx_T are all <= nbTimes
*/
int NbConditionalDefaultNameDistribution_mtp(INDICATOR_SIM i_s,
                                  double *distribution,
                                  long idx_t0,
                                  long idx_t1,
                                  long idx_T0,
                                  long idx_T1)
{
    static char routine[] = "NbConditionalDefaultNameDistribution_mtp";
    int status = FAILURE;
    long nbNames = i_s.nbNames;
    long nbPaths = i_s.nbPaths;
    long i = 0;
    long j = 0;
    long nbDefault = 0;
    double totalWeight = 0.0;
    long nbUsedPaths = 0;
    long noDefaultBefore_t0;
    long default_t0_t1;
    long nameDefaulted_t0_t1;
    for(j=0;j<nbPaths;j++)
    {
        noDefaultBefore_t0 = TRUE;
        default_t0_t1 = FALSE;

        for(i=0;i<nbNames;i++)
        {
            /* check if there si no default before t0 */
            if(i_s.indicator[i+j*nbNames] <= idx_t0)
            {
                noDefaultBefore_t0 = FALSE;
            }

            /* check if there is at least 1 default between t0 and t1 */
            if(i_s.indicator[i+j*nbNames]>idx_t0 &&
               i_s.indicator[i+j*nbNames]<=idx_t1 &&
               default_t0_t1 == FALSE)
            {
                default_t0_t1 = TRUE;
                nameDefaulted_t0_t1 = i;

            }
        }
        
        /* if no default before t0 and at least 1 between t0 and t1 */
        if(noDefaultBefore_t0 && default_t0_t1)
        {
            nbUsedPaths++;
            nbDefault = 0;
            for(i=0;i<nbNames;i++)
            {
                if(i_s.indicator[i+j*nbNames]>idx_T0 &&
                   i_s.indicator[i+j*nbNames]<=idx_T1)
                {
                    nbDefault += 1;
                }
            }
            distribution[nbDefault] += i_s.weight[j];
            totalWeight += i_s.weight[j];
        }
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

/* Convolution
*/
int Convolution(const double *v1, const double *v2, double *v, long n)
{

    long i,j;
    for(i=0;i<n;i++)
    {
        v[i] = 0.0;
        for(j=0;j<=i;j++)
        {
            v[i] += v1[j]*v2[i-j];
        }
    }
    return SUCCESS;
}

/* ------------------------------------------------------------------------
// ComplexMatrix
// returns a Matrix representation of a list of complex numbers
*/
int ComplexMatrix(double *realPart,
                  double *imagPart,
                  long nbValues,
                  double *complexMatrix)
{
    static char routine[] = "ComplexMatrix";
    int status = FAILURE;
    int i;
    for(i=0;i<nbValues;i++)
    {
        if(fabs(imagPart[i])<3e-100)
        {
        //http://www.netlib.org/lapack/lug/node50.html
        //A second basic task is to compute the Schur factorization of a matrix A. If A is complex, then its Schur factorization is A=ZTZH, where Z is unitary and T is upper triangular. If A is real, its Schur factorization is A=ZTZT, where Z is orthogonal. and T is upper quasi-triangular (1-by-1 and 2-by-2 blocks on its diagonal). The columns of Z are called the Schur vectors of A. The eigenvalues of A appear on the diagonal of T; complex conjugate eigenvalues of a real A correspond to 2-by-2 blocks on the diagonal of T.     
            complexMatrix[nbValues*i+i] = realPart[i];
        }
        else
        {
            //if(imagPart[i]!=-imagPart[i+1] || realPart[i]!=realPart[i+1])
             //   goto RETURN;
            
            complexMatrix[nbValues*i+i] = realPart[i];
            complexMatrix[nbValues*(i+1)+(i+1)] = realPart[i+1];
            complexMatrix[nbValues*i+(i+1)] = imagPart[i];
            complexMatrix[nbValues*(i+1)+i] = imagPart[i+1];
            i++;
        }
    }

    status = SUCCESS;
//RETURN:
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}


/* ------------------------------------------------------------------------
// MatrixProduct
// returns a product of 2 Matrices
*/
int MatrixProduct(double *M1,
                  double *M2,
                  long n,
                  double *M)
{
    static char routine[] = "MatrixProduct";
    int status = FAILURE;
    int i,j,k;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<n;k++)
            {
                M[j+n*i] += M1[k+i*n]*M2[j+k*n];
            }
        }
    }

    status = SUCCESS;
//RETURN:
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

void Integer(long n, double *out)
{
    long k;
    for(k=0;k<n;k++)
    {
        out[k] = k;
    }
}

void MatInteger(long n, double *out)
{
    long i,j;
    for(i=0;i<n;i++)
    {
        for(j=0;j<n;j++)
        {
            out[i+n*j] = i-j;
        }
    }
}

/***************************************************************************/
/*  LARGE POOL INDEP-GAUSS-DEP copula
/***************************************************************************/

double TrancheletLoss3C(double u,
                        double beta,
                        double a,
                        double c,
                        double K
                       )
{
    double ua = pow(u,a);
    double ub = pow(u,1.-a-c);
    double uc = pow(u,c);
    double sqr = sqrt(1. - beta*beta);
    double Ninvub = NormalCumInverse(ub);
    double NinvK  = NormalCumInverse( (1.-K) / ua);
    return 1. - uc * NormalCum( (Ninvub - sqr * NinvK)/beta );
}

double TrancheLoss3C(   double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1)
{
    double ua = pow(u,a);
    double uc = pow(u,c);
    double ub = u / (ua*uc);
    double beta2 = beta*beta;
    double sqr = sqrt(1.-beta2);
    double Ninvub = NormalCumInverse(ub);
    double Ninvk1 = NormalCumInverse((1.-K1)/ua);
    double Ninvk0 = NormalCumInverse((1.-K0)/ua);
    double N20, N21;
    BiNormalCum( -Ninvk1, Ninvub, -sqr, &N21);
    BiNormalCum( -Ninvk0, Ninvub, -sqr, &N20);
    return 1. - ua*uc / (K1-K0) * (N21 - N20);
}

double DTrancheletLoss3CDbeta(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       )
{
    double tl1 = TrancheletLoss3C(u,beta+h,a,c,K);
    double tl0 = TrancheletLoss3C(u,beta,a,c,K);
    return (tl1-tl0)/h;
}

double DTrancheletLoss3CDa(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       )
{
    double tl1 = TrancheletLoss3C(u,beta,a+h,c,K);
    double tl0 = TrancheletLoss3C(u,beta,a,c,K);
    return (tl1-tl0)/h;
}

double DTrancheletLoss3CDc(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       )
{
    double tl1 = TrancheletLoss3C(u,beta,a,c+h,K);
    double tl0 = TrancheletLoss3C(u,beta,a,c,K);
    return (tl1-tl0)/h;
}

double DTrancheletLoss3CDcAAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double Kmatch,
                        double h
                       )
{
    double tl1 = TrancheletLoss3C(u,beta,Aeq(u,beta,beta,c+h,Kmatch),c+h,K);
    double tl0 = TrancheletLoss3C(u,beta,a,c,K);
    return (tl1-tl0)/h;
}

double DTrancheletLoss3CDcATAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double Kmatch0,
                        double Kmatch1,
                        double h
                       )
{
    double tl1 = TrancheletLoss3C(u,beta,ATeq(u,beta,beta,c+h,Kmatch0,Kmatch1),c+h,K);
    double tl0 = TrancheletLoss3C(u,beta,a,c,K);
    return (tl1-tl0)/h;
}

double DTrancheLoss3CDbeta(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       )
{
    double tl1 = TrancheLoss3C(u,beta+h,a,c,K0,K1);
    double tl0 = TrancheLoss3C(u,beta,a,c,K0,K1);
    return (tl1-tl0)/h;
}

double DTrancheLoss3CDa(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       )
{
    double tl1 = TrancheLoss3C(u,beta,a+h,c,K0,K1);
    double tl0 = TrancheLoss3C(u,beta,a,c,K0,K1);
    return (tl1-tl0)/h;
}

double DTrancheLoss3CDc(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       )
{
    double tl1 = TrancheLoss3C(u,beta,a,c+h,K0,K1);
    double tl0 = TrancheLoss3C(u,beta,a,c,K0,K1);
    return (tl1-tl0)/h;
}

double DTrancheLoss3CDcATAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double Kmatch0,
                        double Kmatch1,
                        double h
                       )
{
    double tl1 = TrancheLoss3C(u,beta,ATeq(u,beta,beta,c+h,Kmatch0,Kmatch1),c+h,K0,K1);
    double tl0 = TrancheLoss3C(u,beta,a,c,K0,K1);
    return (tl1-tl0)/h;
}

double DTrancheLoss3CDcAAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double Kmatch,
                        double h
                       )
{
    double tl1 = TrancheLoss3C(u,beta,Aeq(u,beta,beta,c+h,Kmatch),c+h,K0,K1);
    double tl0 = TrancheLoss3C(u,beta,a,c,K0,K1);
    return (tl1-tl0)/h;
}

/************************************************************************/
/*  this function returns the a to put in the 3C model to
/*  offset the effect of c and get the same
/*  tranchelet EL% for strike K as CreditMetrics
/************************************************************************/
typedef struct {double u; double beta_star; double beta; double c; double K;} AParam;
static int ABrent(double x, void *data, double *out)
{   
    AParam *p           = (AParam*) data;
    double u            = p->u;
    double beta_star    = p->beta_star;
    double beta         = p->beta;
    double c            = p->c;
    double K            = p->K;
    double sqr_star     = sqrt(1. - beta_star * beta_star);
    double sqr          = sqrt(1. - beta*beta);
    double d_cm = (NormalCumInverse(u) - sqr_star* NormalCumInverse(1.-K)) / beta_star;
    double d_3c = (NormalCumInverse(pow(u,1.-x-c)) - sqr * NormalCumInverse((1.-K)/pow(u,x))) / beta;
    d_cm        = NormalCum(d_cm);
    d_3c        = pow(u,c)*NormalCum(d_3c);  
    *out = d_3c - d_cm;
    return 0;
}

double Aeq( double u,
            double beta_star,
            double beta,
            double c,
            double K)
{
    AParam p;
    double x;
    double guess = c;
    p.beta         = beta;
    p.beta_star    = beta_star;
    p.c            = c;
    p.u            = u;
    p.K            = K;
    
    RootFindBrent(&ABrent, &p, 0., 1.-c, 
            100, guess, 1e-2, 0., 1e-6, 1e-6, &x);       
    return x;
}

/************************************************************************/
/*  this function returns the a to put in the 3C model to
/*  offset the effect of c and get the same
/*  tranche EL% for strike K as CreditMetrics
/************************************************************************/
typedef struct {double u; double beta_star; double beta; double c; double K0; double K1;} ATParam;
static int ATBrent(double x, void *data, double *out)
{   
    ATParam *p           = (ATParam*) data;
    double u            = p->u;
    double beta_star    = p->beta_star;
    double beta         = p->beta;
    double c            = p->c;
    double K0            = p->K0;
    double K1           = p->K1;
    double d_cm = TrancheLoss3C(u,beta_star,0.,0.,K0,K1);
    double d_3c = TrancheLoss3C(u,beta,x,c,K0,K1);
    *out = d_3c - d_cm;
    return 0;
}

double ATeq( double u,
            double beta_star,
            double beta,
            double c,
            double K0,
            double K1)
{
    ATParam p;
    double x;
    double guess = c;
    p.beta         = beta;
    p.beta_star    = beta_star;
    p.c            = c;
    p.u            = u;
    p.K0           = K0;
    p.K1           = K1;

    RootFindBrent(&ATBrent, &p, 0., 1.-c, 
            100, guess, 1e-2, 0., 1e-6, 1e-6, &x);       
    return x;
}


double LossDensity3C(double u, double beta, double a, double c, double K)
{
    double sqr = sqrt(1.-beta*beta);
    double ua = pow(u,a);
    double ub = pow(u,1.-a-c);
    double uc = pow(u,c);
    double Ninvub = NormalCumInverse(ub);
    double NinvK  = NormalCumInverse((1.-K)/ua);
    double x = NormalDensity(1./beta * (Ninvub - sqr * NinvK));
    x *= uc * sqr / (ua * beta *NormalDensity(NinvK));
    if(fabs(K)<3e-15) return 0.0;
    if(fabs(K-1.)<3e-15) return 0.0;
    if( (1.-K)/ua >= 1.) return 0.0;
    return x;
}


/* function to be fed to the GSL adaptive integration method */
typedef struct {double u; double beta; double a; double c; double K; long n;} PARAM3C;
int derivs_gsl_moment3C(double x, double *y, void *param)
{
    PARAM3C *p = (PARAM3C*) param;
    long n = p->n;
    double r;
    switch(n)
    {
    case 0:
        r = 1;
        break;
    case 1:
        r = x;
        break;
    case 2:
        r = x*x;
        break;
    default:
        r = pow(x,n);
        break;
    }
    *y = r*LossDensity3C(p->u,p->beta,p->a, p->c,x);
    return 0;
}

double Moment3C( double u, double beta, double a, double c, long n, long nbPoint, double eps)
{
    static char routine[] = "Moment3C";
    int status      = FAILURE;
    double result   = 0.0;
    PARAM3C *p = NULL;
    p = malloc(sizeof(PARAM3C));
    if(p==NULL) goto RETURN;

    p->u     = u;
    p->a     = a;
    p->c     = c;
    p->beta  = beta;
    p->n     = n;

    /* general case, beta != 0 and beta != 1 */
    /********************************************************************/
    /* adaptive GSL method
    */
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nbPoint);
        gsl_integration_qag (&derivs_gsl_moment3C,0.,1., eps, 0, nbPoint, 6, workspace, &result, &error, p);
        gsl_integration_workspace_free(workspace);
    }
    
    result += (1-pow(u,c));
    status = SUCCESS;

RETURN:
    if(p) free(p);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return result;    
}



typedef struct {double u; double beta_star; double a; double c; long nbPoint; double eps;} BParam;
static int BBrent(double x, void *data, double *out)
{   
    BParam *p           = (BParam*) data;
    double u            = p->u;
    double beta_star    = p->beta_star;
    double c            = p->c;
    double a            = p->a;
    double eps          = p->eps;
    long   nbPoint      = p->nbPoint;
    double EL2  =   Moment3C(u,x,a/(1.-c),0,2,nbPoint,eps);
    double EL20 =   Moment3C(u,beta_star,0.,0.,2,nbPoint,eps);

    *out = EL2 - EL20;
    return 0;
}

double BetaEq(double u, double beta_star, double a, double c, long nbPoint, double eps)
{
    BParam p;
    double x;
    double guess   = beta_star;
    p.beta_star    = beta_star;
    p.a            = a;
    p.c            = c;
    p.u            = u;
    p.eps          = eps;
    p.nbPoint      = nbPoint;
    
    RootFindBrent(&BBrent, &p, 0., 1., 
            100, guess, 1e-2, 0., 1e-6, 1e-6, &x);
    return x;    
}

