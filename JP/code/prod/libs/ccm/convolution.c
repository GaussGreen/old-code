#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include "error2.h"
#include "proba_utils.h"
#include "skew.h"
#include "convolution.h"

#define MALLOC(TYPE,N) (TYPE*)malloc(sizeof(TYPE)*(N))
#define FREE(P) {if (P) free(P);}
#define SKEW_FCUMUL_NB_POINT 200L
#define SKEW_INTEG_LOW       -7.
#define SKEW_INTEG_UP         7.
#define SKEW_INTEG_NB       101
#define SKEW_DENSITY_EPSI     1e-4
#define TOWER_PROBA_EPSI      1e-12

/* ---------------------------------------------------------------------------
 * Independent Variable Loss Distribution Calculation
 */

/**
 * convolution
 * [0,1]^n       --> [0,1]^{n+1}
 * Pr(tau_i > t) |-> Pr(L = k), i \in [0,n] k \in [0, n+1]
 * calculate loss distribution density at time horizon t for independent
 * default events indicator 1_{tau_i > t}.
 */
static void convolution(
    const int    *lossAmt,  /* (I) array[n] of loss amount in units (must >= 0) */
    const double *survival, /* (I) array[n] name survival probabilities */
    int     n,        /* (I) number of names */
    int     maxLoss,  /* (I) sum of lossAmt */
    double *work,     /* (X) array[maxLoss+1] work area */
    double *out)      /* (O) array[maxLoss+1] of densities (must be allocated) */
{
    int i,j;
    double p,q;
    int umax, l;
    assert(out && lossAmt && survival && work);

    /* init density */
    memset(out, 0, sizeof(double)*(maxLoss+1));
    out[0] = 1.;

    /* loop over name */
    for (i=0, umax=1; i<n; ++i)
    {
        p = survival[i];
        q = 1.-p;
        l = lossAmt[i];
        memcpy(work,out,sizeof(double)*umax);
        for (j=0; j<umax; ++j) /* p-scale */
            out[j] *= p;

        assert(umax+l <= maxLoss+1);
        for (j=0; j<umax; ++j) /* add q-scale shift */
            out[j+l] += q * work[j];
        umax += l;
    }
    assert(umax == maxLoss+1);
}
   
/* ---------------------------------------------------------------------------
 * Frailty Copula conditional Param
 */
typedef struct {
    double *T;        /* array[n] of unconditional survival threshold */
    double *b;        /* array[n] of beta */
    double *a;        /* array[n] of (1-beta^2)^{-1/2} */
    double  qM,qZ;    /* skew of M and Z */
    int     n;        /* number of name n */
} GaussCondParam;             

/** free cond proba data struct */
static void freeCondProba(void *param)
{
    GaussCondParam *out = (GaussCondParam *)param;
    if (!out) return;
    FREE(out->T);
    FREE(out->a);
    FREE(out->b);
    memset(out,0,sizeof(GaussCondParam));
    free(out);
}

/** init param for conditional calculation (copula dependent preprocessing) */
void *initCondProba(
    const double *p, /* (I) array[n] of unconditional survival proba */
    double *b, /* (I) array[n] of beta */
    double  qM,/* (I) skew of M */
    double  qZ,/* (I) skew of Z */
    int     n) /* (I) number of names */
{
    static char routine[] = "initCondProba";
    int status = -1;
    int i;
    GaussCondParam *out = MALLOC(GaussCondParam,1);
    if (!out) goto fail;
    memset(out,0,sizeof(GaussCondParam));
    out->T = MALLOC(double,n);
    out->b = MALLOC(double,n);
    out->a = MALLOC(double,n);
    if (!out->T || !out->b || !out->a) goto fail;
    out->qM = qM;
    out->qZ = qZ;
    out->n = n;
    for (i=0; i<n; ++i) 
    {
        if (fabs(p[i]-.5) > .5) {status=-2; goto fail;}
        if (fabs(b[i]) >= 1.)    {status=-3; goto fail;}
        out->T[i] = Finv(p[i],b[i],qM,qZ,SKEW_FCUMUL_NB_POINT, 10e-8);
        out->b[i] = b[i];
        out->a[i] = pow(1.-b[i]*b[i], -.5);
    }   

    return out;
fail:
    switch (status){
    case -2: DR_Error("%s: all survival probabilities have to be in ]0,1[", routine); break;
    case -3: DR_Error("%s: all beta have to be in ]-1,1[", routine); break;
    default: DR_Error("%s: failed", routine);
    }
    freeCondProba(out);
    return NULL;
}

/** calculate survival proba conditional on M */
static int calcCondProba(
    void     *p, /* (I) ouptut of initCondProba */
    double    M, /* (I) M */
    double   *pm)/* (O) survival proba conditional on M */
{
    GaussCondParam *out = p;
    int i;
    assert(out != NULL && pm != NULL);
    for (i=0; i<out->n; ++i)
        pm[i] = NormalCum(NormMapinv((out->T[i] - out->b[i]*M) * out->a[i],
                        out->qZ));
    return 0;
}

/* ---------------------------------------------------------------------------
 * Integ Method
 */
typedef struct {
    double l; /* lower bound */
    double u; /* upper bound */
    int    n; /* number of steps */
    double step; /* step size */
} GaussIntegMethod;

static void *initIntegMethod()
{
    static char routine[] = "initIntegMethod";
    GaussIntegMethod *g = MALLOC(GaussIntegMethod,1);
    if (!g) goto fail;
    g->l = SKEW_INTEG_LOW;
    g->u = SKEW_INTEG_UP;
    g->n = SKEW_INTEG_NB;
    g->step = (g->u - g->l)/(g->n - 1);
    return g;
fail:
    return NULL;
}

static void freeIntegMethod(void *p)
{
    FREE(p);
}

static double calcIntegWeight(void *condParam, void *integMethod, double M)
{
    GaussIntegMethod *i = (GaussIntegMethod *)integMethod;
    GaussCondParam   *c = (GaussCondParam *)condParam;
    double den = i->step 
        * NormalDensity(NormMapinv(M,c->qM)) 
        * DNormMapinv(M,c->qM);
    return den;
}

/* ---------------------------------------------------------------------------
 * Loss Density
 */

/* calculate loss density for skew and gauss copula */
static int calcSkewLossDensity(
    const int    *lossAmt,   /* (I) array[n] of loss amt              */
    int           n,         /* (I) number of names                   */
    int           maxLoss,   /* (I) sum(loss amt)                     */
    const COPULA *c,         /* (I) copula 0=gauss, 1=skew            */
    const double *p,         /* (I) array[n] of uncond survival proba */
    double       *density    /* (O) array[maxLoss+1] of loss density  */
    )
{
    GaussCondParam   *condParam = NULL;
    GaussIntegMethod *integMethod = NULL;
    double *out = NULL, *work = NULL;
    double *pm  = NULL;
    int     i,j;
    int     status = -10;
    /* determine beta and skew param */
    double M,w,sumW;
    double qM,qZ;
    double *b;
    if (c->type == C_SKEW)
    {
        Q_PARAM *param = (Q_PARAM *) c->param;
        b  = param->beta;
        qM = param->qM;
        qZ = param->qZ;
    }
    else
    {
        GAUSSIAN_PARAM *param = (GAUSSIAN_PARAM *) c->param;
        if (c->type != C_GAUSS) {status = -1; goto fail;}
        b  = param->beta;
        qM = 0.;
        qZ = 0.;
    }

    /* allocate */
    out  = MALLOC(double,maxLoss+1); /* cond loss density */
    work = MALLOC(double,maxLoss+1); /* work array */
    pm   = MALLOC(double,n);         /* cond survival proba */
    if (!out || !work || !pm) goto fail;

    /* init data */
    condParam = initCondProba(p,b,qM,qZ,n);
    if (!condParam) goto fail;
    integMethod = initIntegMethod();
    if (!integMethod) goto fail;
    memset(density,0,sizeof(double)*(maxLoss+1));

    /* integrate */
    for (i=0, M=integMethod->l, sumW=0.; i<integMethod->n; ++i)
    {
        calcCondProba(condParam, M, pm);
        convolution(lossAmt,pm,n,maxLoss,work,out);

        w     = calcIntegWeight(condParam,integMethod,M);
        sumW += w;
        for (j=0; j<maxLoss; ++j)
            density[j] += w * out[j];

        M += integMethod->step;
    }
    assert(fabs(sumW-1.) < SKEW_DENSITY_EPSI);

    status = 0;        
fail:
    FREE(out);
    FREE(work);
    FREE(pm);
    freeCondProba(condParam);
    freeIntegMethod(integMethod);
    return status;
}

/* qsort utils for ordering pdep in ascending order */
#ifndef MAX
#define MAX(a,b) ((a)>(b) ? (a) : (b))
#endif
#define SIGN(a) ((a)>0. ? 1 : ((a)<0. ? -1 : 0))

typedef struct {
    int    index;
    double pdep;
} IntDoublePair;

static int IntDoublePairCompare(const void *e1, const void *e2) {
    return SIGN(((IntDoublePair *)e2)->pdep - ((IntDoublePair *)e1)->pdep);
}

/*
 * Look into the composite copula, find which one is the dependence
 * and which one is the skew or gaussian, and get probabilities.
 * note: error status gets passed to the caller for nice error msg
 */
static int calcCompositeCopulaLoss(
    const int    *lossAmt,   /* (I) array[n] of loss amt              */
    int           n,         /* (I) number of names                   */
    int           maxLoss,   /* (I) sum(loss amt)                     */
    const COPULA *c,         /* (I) copula 0=gauss, 1=skew            */
    const double *p,         /* (I) array[n] of uncond survival proba */
    double       *density    /* (O) array[maxLoss+1] of loss density  */
    )
{
    static char routine[] = "calcCompositeCopulaLoss";
    PRODUCT_PARAM *pp = (PRODUCT_PARAM *) c->param;
    int status = -10;
    double *b;
    IntDoublePair *idp = NULL;
    double *p1 = NULL, *p2 = NULL;
    double *work = NULL;
    COPULA dep, skew;
    double *pdep, *pskew;
    int    *lossAmt2 = NULL;       /* reordered loss amt */
    int     i,j;
    double sumW;
    int    shift;
    int    isSkewFirst = (pp->c_2->type == C_DEPENDENCE) ? 1 : 0;

    /* needed to ensure CopulaFree does not try to free unallocated stuff */
    dep.type  = C_COPULA_FREE;
    skew.type = C_COPULA_FREE;

    if (isSkewFirst)
    {
        if (CopulaCopy(pp->c_2,&dep,n))  goto fail;
        if (CopulaCopy(pp->c_1,&skew,n)) goto fail;
    }
    else
    {
        if (CopulaCopy(pp->c_1,&dep,n))  goto fail;
        if (CopulaCopy(pp->c_2,&skew,n)) goto fail;
    }
    if (dep.type  != C_DEPENDENCE)  {status =-3; goto fail;}
    if (skew.type == C_DEPENDENCE)  {status =-3; goto fail;}

    p1 = MALLOC(double,n);        
    p2 = MALLOC(double,n); 
    for (i=0; i<n; ++i)
    {
        if (fabs(p[i]-.5)>.5)         {status=-5; goto fail;}
        if (fabs(pp->alpha[i]-.5)>.5) {status=-6; goto fail;}
        p1[i] = pow(p[i],pp->alpha[i]);
        p2[i] = p[i] / p1[i];
        assert(fabs(p1[i]-.5) <= .5 && fabs(p2[i]-.5) <= .5);
    }
    
    if (isSkewFirst)
    {
        pdep  = p2;
        pskew = p1;
    }
    else
    {
        pdep  = p1;
        pskew = p2;
    }

    /* need to sort pdep, pskew, beta vector by ascending pdep */
    idp = MALLOC(IntDoublePair,n);
    if (!idp) goto fail;
    lossAmt2 = MALLOC(int,n);
    if (!lossAmt2) goto fail;

    for (i=0; i<n; ++i)
    {
        idp[i].index = i;
        idp[i].pdep  = pdep[i];
    }
    qsort(idp,n,sizeof(IntDoublePair),IntDoublePairCompare);

    switch (skew.type) {
    case C_SKEW:  b  = ((Q_PARAM *)skew.param)->beta; break;
    case C_GAUSS: b  = ((GAUSSIAN_PARAM *)skew.param)->beta; break;
    default: status = -1; goto fail;
    }

    work = MALLOC(double,MAX(maxLoss+1,n));
    if (!work) goto fail;   
    /* reorder pdep, pskew, beta, lossAmt */
    for (i=0; i<n; ++i) work[i] = pdep[idp[i].index];
    memcpy(pdep,work,sizeof(double)*n);
    for (i=0; i<n; ++i) work[i] = pskew[idp[i].index];
    memcpy(pskew,work,sizeof(double)*n);
    for (i=0; i<n; ++i) work[i] = b[idp[i].index];
    memcpy(b,work,sizeof(double)*n);
    for (i=0; i<n; ++i) lossAmt2[i] = lossAmt[idp[i].index];
    /* end of sorting */


    /* 
     * do a convolution algorithm using pskew
     * conditional on a dependent default scenario.
     * The loss is then shifted by the sure loss amount and multiplied by 
     * pdep. We look at n+1 states, whose proba is 
     * pdep[j-1] - pdep[j] where j in [0..n]
     *      pdep[j] < pdep[j-1] for all j
     *      pdep[-1] = 1, pdep[n] = 0 as extreme conditions
     */
    /* take care of case j=0 (all defaults) separately */
    shift = maxLoss;
    memset(density,0,sizeof(double)*maxLoss);
    density[shift] = sumW = 1.0 - pdep[0];

    for (j=1; j<=n; ++j)
    {
        double weight = j==n ? pdep[n-1] : pdep[j-1] - pdep[j];
        sumW  += weight;
        shift -= lossAmt2[j-1];
        if (fabs(weight) < TOWER_PROBA_EPSI) 
            continue;

        assert(weight>=0.);

        if (calcLossDensity(lossAmt2,j, maxLoss-shift, &skew, pskew, work+shift))
            goto fail;

        for (i=shift; i<maxLoss+1; ++i)
            density[i] += weight * work[i];
    }
    assert(shift == 0);
    assert(fabs(sumW-1.) < TOWER_PROBA_EPSI * n);
    status = 0;
fail:
    FREE(p1);
    FREE(p2);
    FREE(lossAmt2);
    FREE(idp);    
    FREE(work);
    CopulaFree(&dep);
    CopulaFree(&skew);
    return status;
}

/** calculate loss density */
int calcLossDensity(
    const int    *lossAmt,   /* (I) array[n] of loss amt              */
    int           n,         /* (I) number of names                   */
    int           maxLoss,   /* (I) sum(loss amt)                     */
    const COPULA *c,         /* (I) copula 0=gauss, 1=skew            */
    const double *p,         /* (I) array[n] of uncond survival proba */
    double       *density    /* (O) array[maxLoss+1] of loss density  */
    )
{
    static char routine[] = "calcLossDensity";
    int     status = -10;

    switch (c->type) {
    case C_SKEW: 
    case C_GAUSS:
        status = calcSkewLossDensity(lossAmt,n,maxLoss,c,p,density);
        if (status) goto fail;
        break;
    case C_COMPOSITE : 
        status = calcCompositeCopulaLoss(lossAmt,n,maxLoss,c,p,density);
        if (status) goto fail;
        break;
    default:
        status = -1;
        goto fail;
    }

    status = 0;        
fail:
    switch (status) {
    case 0: break;
    case -1: DR_Error("%s: only skew and gauss copula and their composition with dependence are implemented", routine); break;
    case -2: DR_Error("%s: loss unit must be greater than 1e-4", routine); break;
    case -3: DR_Error("%s: composite copula must have a dependence copula", routine); break;
    case -4: DR_Error("%s: names must be ordered with ascending dep survival proba", routine); break;
    case -5: DR_Error("%s: names must have probabilities in [0,1]", routine); break;
    case -6: DR_Error("%s: names must have alpha in [0,1]", routine); break;
    default: DR_Error("%s: failed", routine);
    }
    return status;
}
