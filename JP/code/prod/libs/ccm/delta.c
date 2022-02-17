#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "rootbrent.h"
#include "delta.h"
#include "random_utils.h"
#include "proba_utils.h"
#include "gaussian.h"
#include "error2.h"

#ifndef M_PI
#define M_PI	   3.14159265358979323846264338328      /* pi */
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif


/* -------------------------------------------------------------------------
** DeltaNormalisation
**
*/
double DeltaNormalisation(double a, double b, double c, double delta)
{
    
    double result = 0.0;
    result = 1. + delta * (c*c - 1. +  1./sqrt(2.*M_PI) * (b*exp(-b*b*0.5)-a*exp(-a*a*0.5)) / ( NormalCum(b)-NormalCum(a) )  );  
    result -= delta*delta *(c - 1./sqrt(2.*M_PI) * fabs(exp(-b*b*0.5)-exp(-a*a*0.5))/(NormalCum(b)-NormalCum(a)) ) * (c - 1./sqrt(2.*M_PI) *fabs(exp(-b*b*0.5) - exp(-a*a*0.5))/(NormalCum(b) - NormalCum(a))) ;
    result = 1./sqrt(result);
    return result;
} 

/* -------------------------------------------------------------------------
** FdeltaC
** Fdelta without the Dirac
*/
double FdeltaC(double x, double a, double b, double delta)
{
    static char routine[] = "FdeltaC";
    int status = FAILURE;
    double result = 0.0;
    double alpha = delta / (NormalCum(b)-NormalCum(a));

    if(x<=a)
    {
        result = NormalCum(x);
    }
    else if(x>b)
    {
        result = NormalCum(a)
                    +(1.-alpha)*(NormalCum(b)-NormalCum(a))
                    +(NormalCum(x)-NormalCum(b));
    }
    else
    {
        result = NormalCum(a) + (1.-alpha)*(NormalCum(x)-NormalCum(a));
    }
    status = SUCCESS;
    return result;
}

/* -------------------------------------------------------------------------
** FdeltaCinv
** inverse of the previous function
*/
double FdeltaCinv(double y, double a, double b, double delta)
{
    double x = 0.0;
    double alpha = delta / (NormalCum(b)-NormalCum(a));

    if(y<FdeltaC(a,a,b,delta))
    {
        x = NormalCumInverse(y);
    }
    else if(y>FdeltaC(b,a,b,delta))
    {
        x = NormalCumInverse(y + alpha*(NormalCum(b)-NormalCum(a)));
    }
    else
    {
        x = NormalCumInverse(NormalCum(a)+(y-NormalCum(a))/(1.-alpha));
    }
    return x;
}


/* -------------------------------------------------------------------------
** Fdelta
**
*/
double Fdelta(double x, double a, double b, double c, double delta)
{
    double result = 0.0;
    result = FdeltaC(x,a,b,delta);

    if(x>=c)
    {
        result += delta;
    }

    return result;
}

/* -------------------------------------------------------------------------
** Fdeltainv
** this function solves the equation Fdelta(u) = p
*/
double Fdeltainv(   double y,
                    double a,
                    double b,
                    double c,
                    double delta)
{
    double x = 0.0;

    if( y<FdeltaC(c,a,b,delta) )
    {
        x = FdeltaCinv(y,a,b,delta);
    }
    else if( y>Fdelta(c,a,b,c,delta) )
    {
        x = FdeltaCinv(y-delta,a,b,delta);
    }
    else
    {
        x = c;
    }

    return x;
}

/* -------------------------------------------------------------------------
** DeltaMap
**
*/
double DeltaMap(double x, double a, double b, double c, double delta)
{
    if(x < a && x < c)
    {
        return x;
    }
    else
    {
        return Fdeltainv(NormalCum(x),a,b,c,delta);
    }
}


/* -------------------------------------------------------------------------
** Integrandd
** the integral of this function for Zi=-infinity to Zi=+infinity
** is equal to Fd(u) = P(Xi<u)
*/
double Integrandd(double u, double Zi, double beta, double a, double b, double c, double delta)
{
    return  Fdelta((u-sqrt(1.-beta*beta)*Zi)/ beta ,a,b,c,delta) * 1./sqrt(2*M_PI) * exp(-Zi*Zi*0.5);
}

/* -------------------------------------------------------------------------
** Fd
** this function computes the integral of the previous function 
**
*/
double Fd(double u, double beta, double a, double b, double c, double delta, long nbPoints)
{
    static char routine[] = "Fd";
    int status = FAILURE;
    double result = 0.0;
    double *x = NULL;
    long k;
    double x_inf, x_sup;

    x = malloc(nbPoints*sizeof(double)); //these are the abscissas 
    if(x==NULL) goto RETURN;

    x_inf = -7.;
    x_sup = 7.;
    
    /* if beta = 1, Xi = M */
    if(fabs(beta-1.)<3e-16)
    {
        return NormalCum(u);
    }

    /* if beta = 0, Xi = Zi */
    if(fabs(beta)<3e-16)
    {
        return Fdelta(u,a,b,c,delta);
    }

    /* general case, beta != 0 and beta != 1 */
    for(k=0;k<nbPoints;k++)
    {
        x[k] = x_inf + k * (x_sup - x_inf) /nbPoints;
    }

    result = Integrandd(u,x[0],beta,a,b,c,delta)*(x[1]-x[0]);
    for (k=1;k<nbPoints-1;k++)
    {
        result += Integrandd(u,x[k],beta,a,b,c,delta)*(x[k+1] - x[k-1]);
    }
    result += Integrandd(u,x[nbPoints-1],beta,a,b,c,delta)*(x[nbPoints-1]-x[nbPoints-2]);
    result *= 0.5;

    status = SUCCESS;

RETURN:
    if(x) free(x);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return MAX(MIN(result,1.),0.);
}


/** root solving implementation */
typedef struct {double y; double beta; double a; double b; double c; double delta; long nbPoints;} ParamStruct;
static int FBrent(double x, void *data, double *out)
{
    ParamStruct *p = (ParamStruct*)data;
    *out = Fd(x, p->beta,p->a,p->b,p->c,p->delta,p->nbPoints) - p->y;
    return 0;
}

/* -------------------------------------------------------------------------
** Fdinv
** this function solves the equation F(u) = p with input p
*/
double Fdinv(double y,
            double beta,
            double a,
            double b,
            double c,
            double delta,
            long nbPoints)
{
    double x;
    ParamStruct p;
    if (y < 3e-16) 
    {
            return -3e100;
    }
    else
    {
        if (1.-y < 3e-16) return 3e100;
        p.y = y;
        p.beta = beta;
        p.a = a;
        p.b = b;
        p.c = c;
        p.delta = delta;
        p.nbPoints = nbPoints;

        RootFindBrent(&FBrent, &p, -1e10, 1e10, 
            50, 0., 1., 0., 1e-8, 1e-8, &x);
    }
    return x;
}


/* -------------------------------------------------------------------------
** DeltaMultivariateDeviates_tp
** returns beta correlated q random variables
*/
int DeltaMultivariateDeviates_tp(
    double *deltaSequence,              /* (O) [nbPaths*nbNames] */
    double *weight,                 /* (O) weight[nbPaths] */ 
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double a,                       /* (I) */
    double b,                       /* (I) */
    double c,                       /* (I) */
    double delta,                   /* (I) */
    long seed)                      /* (I) */
{
    static char routine[] = "DeltaMultivariateDeviates_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *sqrt_beta = NULL;
    double e,e2;
    e = 0.0;
    e2  =0.0;

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
    status = CreateUniformSequence(M,seed,1,nbPaths);
    if(status == FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            deltaSequence[i+j*nbNames] =
                beta[i]*Fdeltainv(M[j],a,b,c,delta) + sqrt_beta[i]*Z[i+j*nbNames];
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


/* -------------------------------------------------------------------------
** DeltaCopulatedIndicator
** this function returns the indicators of default using a q-copula
*/
int DeltaCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double a,                       /* (I) */
    double b,                       /* (I) */
    double c,
    double delta,
    long seed)
{
    static char routine[] = "DeltaCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *deltaSequence = NULL;
    double *T = NULL;

    deltaSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(deltaSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames;i++)
    {
        T[i] = Fdinv(survivalProba[i],beta[i],a,b,c,delta,200);
    }

    /* allocate the q sequence */
    status = DeltaMultivariateDeviates_tp(deltaSequence, weight, nbNames, nbPaths, beta, a,b,c,delta, seed);
    if(status ==FAILURE) goto RETURN;
    

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(deltaSequence[i+j*nbNames] > T[i])
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
    if(deltaSequence) free(deltaSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

/* ---------------------------------------------------------------------------
// DeltaCopulatedIndicator_mtp
// this function returns the survival indicator and the weight of each path
//
*/
int DeltaCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    double a,
    double b,
    double c,
    double delta,
    long seed)
{
    static char routine[] = "DeltaCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double *deltaSequence = NULL;
    double *T = NULL;

    deltaSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(deltaSequence == NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbTimes;j++)
        {
            T[i+j*nbNames] = Fdeltainv(survivalProba[i+j*nbNames],a,b,c,delta);
        }
    }

    /* allocate the gaussian sequence */
    status = DeltaMultivariateDeviates_tp(deltaSequence, weight, nbNames, nbPaths, beta, a, b, c, delta, seed);
    if(status ==FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
            for(k=0;k<nbTimes;k++)
            {
                if(deltaSequence[i+j*nbNames] > T[k+i*nbTimes])
                {
                    copulatedSurvivalIndicator[i+j*nbNames] = k;
                    break;
                }
            }
        }
    }
    status = SUCCESS;
RETURN:
    if(deltaSequence) free(deltaSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

