#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "rootbrent.h"
#include "skew.h"
#include "random_utils.h"
#include "proba_utils.h"
#include "gaussian.h"
#include "error2.h"


#ifndef M_PI
#define M_PI	   3.14159265358979323846264338328      /* pi */
#endif

#ifndef EPS
#define EPS 0.00000000000003
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif

#ifndef SIGN
#define SIGN(a) (((a) < 0) ? (-1) : (1))
#endif

typedef struct
{
    double u;
    double qM;
    double qZ;
    double beta;
} FPARAM;


/* -------------------------------------------------------------------------
** Fq
** returns the cumulative proba of a skewed var
*/
double Fq(double x, double q)
{
    return NormalCum(NormMapinv(x,q));
}

/* -------------------------------------------------------------------------
** Fqinv
** returns the cumulative proba inverse of a skewed var
*/
double Fqinv(double y, double q)
{
    return NormMap(NormalCumInverse(y),q);
}

/* --------------------------------------------------------------------------
** fq
** return the density of a skewed variable 
*/
double fq(double x, double q)
{
    return NormalDensity(NormMapinv(x,q))*DNormMapinv(x,q);
}

/* --------------------------------------------------------------------------
** f1q
** return the derivative of the density of a skewed variable 
*/
double f1q(double x, double q)
{
    double nmi  = NormMapinv(x,q);
    double Dnmi = DNormMapinv(x,q);
    return -NormalDensity(nmi)*Dnmi*Dnmi * (q + nmi);
}


/* -------------------------------------------------------------------------
** B
** this function is needed in the bounds of M (see below)
*/
double B(double q)
{
    if(fabs(q)<3e-16)
    {
        return -1e300;
    }
    else
    {
        return (q - 1./q);
    }
}


/* -------------------------------------------------------------------------
** Normalisation
** normalisation of the q-mapping to keep the variance of the variables
** equal to 1
*/
double Normalisation(double q)
{
    if(fabs(q)<3e-16)
    {
        return 1.;
    }
    else
    {
        return fabs(q) / sqrt(exp(q*q)*(exp(q*q)-1.));
    }
}

/* -------------------------------------------------------------------------
** Map
**  q mapping from a normal variable
** we normalise the transformed variable to keep the variance equal to 1
*/
double Map(double x, double q)
{
    if(fabs(q)<3e-16)
    {
        return x;
    }
    else if(fabs(q-1.)<3e-16)
    {
        return exp(x);
    }
    else
    {
        return (q + (exp(q*x) - 1.) / q);
    }
}

/* -------------------------------------------------------------------------
** Mapinv
** inverse q mapping 
*
double Mapinv(double y, double q)
{
    if(fabs(q)<3e-16)
    {
        return y;
    }
    else if(fabs(q-1.)<3e-16)
    {
        if(y>0)
        {
            return log(y);
        }
        else
        {
            /* we should never go there *
            return -1e300;
        }
    }
    else
    {
        if(y> B(q))
        {
            return log(q*(y - q) + 1.)/q;
        }
        else
        {
            /* we should never go there *
            return -1e300;
        }
    }
}
*/


/* -------------------------------------------------------------------------
** DMapinv
** derivative of inverse q mapping 
** called by Integrand
*
double DMapinv(double y, double q)
{
    if(fabs(q)<3e-16)
    {
        return 1.;
    }
    else if(fabs(q-1.)<3e-16)
    {
        return 1./y;
    }
    else
    {
        return 1 / (q*(y-q) + 1);
    }
}
*/

/* -------------------------------------------------------------------------
** NormB
** this function is needed in the bounds of M (see below)
** determines a bound of the integral where Integrand is not = 0
*/
double NormB(double q)
{
    if(fabs(q)<3e-16)
    {
        return -1e300;
    }
    else
    {
        return (q*q-1.)/sqrt(exp(q*q)*(exp(q*q)-1.));
    }
}

/* -------------------------------------------------------------------------
** NormMap
**  q mapping from a normal variable
** we normalise the transformed variable to keep the variance equal to 1
*/
double NormMap(double x, double q)
{
    if(fabs(q)<3e-16)
    {
        return x;
    }
    else if(fabs(q-1.)<3e-16)
    {
        return exp(x)/sqrt(exp(1.)*(exp(1.)-1.));
    }
    else
    {
        return SIGN(q) * (q*q+ (exp(q*x) - 1.)) /sqrt( exp(q*q) * (exp(q*q) - 1.));
    }
}

/* -------------------------------------------------------------------------
** NormMapinv
** inverse q mapping 
** normalised to keep the variance equal to 1
*/
double NormMapinv(double y, double q)
{
    if(fabs(q)<3e-16)
    {
        return y;
    }
    else if(fabs(q-1.)<3e-16)
    {
        if(y*sqrt(exp(1.)*(exp(1.)-1.)) >0)
        {
            return log(y*sqrt(exp(1.)*(exp(1.)-1.)));
        }
        else
        {
            /* we should never go there */
            return -1e300;
        }
    }
    else
    {
        if(SIGN(q)*y > NormB(q))
        {
            return log(y*SIGN(q)*sqrt(exp(q*q)*(exp(q*q)-1.)) + 1. - q*q)/q;
        }
        else
        {
            /* we should never go there */
            return -1e300;
        }
    }
}

/* -------------------------------------------------------------------------
** DNormMapinv
** derivative of inverse normalised q mapping 
** called by Integrand
*/
double DNormMapinv(double y, double q)
{
    if(fabs(q)<3e-16)
    {
        return 1.;
    }
    else if(fabs(q-1.)<3e-16)
    {
        return 1./y;
    }
    else
    {
        return 1. / (q * (y + (1.-q*q) * SIGN(q) / sqrt(exp(q*q)*(exp(q*q)-1.))) );
    }
}

/* -------------------------------------------------------------------------
** Integrand
** the integral of this function for M=-infinity to M=+infinity
** is equal to Fi(u) = P(Xi<u)
** Fi(u) = int_{-\infty}^{+\infty} F_Z ( (u-beta_i*M)/sqrt(1-beta_i^2) ) dF_M (M) 
**       = int_{-\infty}^{+\infty} N(NormMapInv( (u-beta_i*M)/sqrt(1-beta_i^2) , qZ))
**                                  * n(NormMapInv(M,qM) * d(NormMapInv(M,qM))/dM
*/
double Integrand(double u, double M, double beta, double qM, double qZ)
{
    double res;
    double z;


    /*  the mapping chosen (q mapping) implies that
        the inverse mapping is not defined everywhere because
        of the log function it contains.
        We need that mapinv(M) and mapinv(z) (see below) are defined.
        The definition interval for M is then
        [B(qM)..1/beta*(u-B(qZ)*sqrt(1-beta^2))] (see B above)
        if beta > 0,
        [B(qM).. +infinity]
        otherwise (in the case qM>=0,qZ>=0
    */
        
    if(M*SIGN(qM) <= NormB(qM))
    {
        return 0.0;
    }

    if((u-beta*M)*SIGN(qZ) <= NormB(qZ)*sqrt(1.-beta*beta))
    {
        if(qZ>=0)
        {
            return 0.0;
        }
        else
        {
                return NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
        }

    }

    if(fabs(beta-1.)<3e-16)
    {
        /* we should never go there */
        z = (u - beta*M)*1e300;
    }
    else
    {
        z = (u - beta*M)/sqrt(1-beta*beta);
    }

    res = NormalCum(NormMapinv(z,qZ))
           *NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
    return res;
}

/* function to be fed to the NR adaptative RK integration method */
int derivs(double x, double y[], double dxdy[], void *param)
{
    FPARAM *Fparam = (FPARAM*) param;
    dxdy[0] = Integrand(Fparam->u,x,Fparam->beta,Fparam->qM,Fparam->qZ);
    return 0;
}

/* function to be fed to the GSL adaptive integration method */
int derivs_gsl(double x, double *y, void *param)
{
    FPARAM *Fparam = (FPARAM*) param;
    *y = Integrand(Fparam->u,x,Fparam->beta,Fparam->qM,Fparam->qZ);
    return 0;
}


/* -------------------------------------------------------------------------
** F
** this function computes the integral of the Integrand function 
**
*/
double F(   double u,
            double beta,
            double qM,
            double qZ,
            long NbPoints,
            double eps,
            DEBUGINFO *debugInfo
            )
{
    static char routine[] = "F";
    int status = FAILURE;
    double result = 0.0;
    //double *x = NULL;
    //long k;

#ifdef NRINT
    int nok,nbad;
    double b1,b2;
    double  x_inf, x_sup;
    double y = 0;
#endif

    FPARAM *Fparam = NULL;
    Fparam = malloc(sizeof(FPARAM));
    if(Fparam==NULL) goto RETURN;


    Fparam->beta = beta;
    Fparam->qM = qM;
    Fparam->qZ = qZ;
    Fparam->u = u;

    //x = malloc(nbPoints*sizeof(double)); //these are the abscissas 
    //if(x==NULL) goto RETURN;
    
    /* if beta = 1, Xi = M */
    if(fabs(beta-1.)<3e-16)
    {
        if(u*SIGN(qM)>NormB(qM))
        {
            return NormalCum(NormMapinv(u,qM));
        }
        else
        {
            return 0.0;
        }
    }

    /* if beta = 0, Xi = Zi */
    if(fabs(beta)<3e-16)
    {
        if(u*SIGN(qZ)>NormB(qZ))
        {
            return NormalCum(NormMapinv(u,qZ));
        }
        else
        {
            return 0.0;    
        }
    }

    /* general case, beta != 0 and beta != 1 */
    /* determination of bounds of integration */
#ifdef NRINT    
    //x_inf = NormMap(-7.,qM);
    //x_sup = NormMap(7.,qM);
    
    /*
    if(fabs(qM)<3e-16)
    {
        if(qZ>3e-16 && beta<-3e-16)
        {
            x_inf = 1./beta*(u-NormB(qZ)*sqrt(1.-beta*beta));
        }
        if(qZ>3e-16 &&beta>3e-16)
        {
            x_sup = 1./beta*(u-NormB(qZ)*sqrt(1.-beta*beta));
        }
    }

    if(fabs(qZ)<3e-16)
    {
        if(qM>3e-16)
        {
            x_inf = NormB(qM);
        }
        if(qM<-3e-16)
        {
            x_sup = -NormB(qM);
        }
    }

    if(qM>3e-16 && qZ>3e-16)
    {
        b1 = NormB(qM);
        b2 = 1./beta*(u - NormB(qZ)*sqrt(1.-beta*beta));
        if(beta>3e-16)
        {
            x_inf = b1;
            x_sup = b2;
        }
        if(beta<-3e-16)
        {
            x_inf = MAX(b1,b2);
        }

    }

    if(qM<-3e-16 && qZ<-3e-16)
    {
            x_sup = - NormB(qM);
    }

    if(qM>3e-16 && qZ<-3e-16)
    {
        x_inf = NormB(qM);
    }

    if(qM<-3e-16 && qZ>3e-16)
    {
        b1 = - NormB(qM);
        b2 = 1./beta*(u - NormB(qZ)*sqrt(1.-beta*beta));
        if(beta<-3e-16)
        {
            x_inf = b2;
            x_sup = b1;
        }
        if(beta>3e-16)
        {
            x_sup = MIN(b1,b2);
        }
    }
    */
    x_inf = -7.0;
    x_sup = 7.0;

    if(x_sup<=x_inf)
    {
        result = 0.0;
        status = SUCCESS;
        goto RETURN;
    }
    
    /********************************************************************/
    /* adaptive Runge-Kutta method (see NR) */
    
    status = odeint(    &y,
                        1 ,
                        x_inf,
                        x_sup,
                        eps,
                        0.1,
                        1e-12,
                        NbPoints,
                        &nok,
                        &nbad,
                        &derivs,
                        &rkqs,
                        Fparam,
                        debugInfo);
    if(status==FAILURE) goto RETURN;
    
#else

    /********************************************************************/
    /* adaptive GSL method
    */
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);
        //if(x_sup>=10 && x_inf<= -10)
        //{
            gsl_integration_qagi (&derivs_gsl, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        /*}
        else if(x_sup < 10 && x_inf <= -10 )
        {
            gsl_integration_qagiu (&derivs_gsl,x_inf, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        }
        else if(x_sup >= 10 && x_inf > -10)
        {
            gsl_integration_qagil (&derivs_gsl,x_sup, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        }
        else 
        {
            gsl_integration_qag (&derivs_gsl, x_inf, x_sup, eps, 0, NbPoints, 6, workspace, &result, &error, Fparam);        
        }*/
        gsl_integration_workspace_free(workspace);
    }
    
#endif


    /********************************************************************/
    /* GOOD OLD INTEGRATION */
    /*
    if(qM>3e-16 && qZ>3e-16 && 1./beta*(u-NormB(qZ)*sqrt(1.-beta*beta))>NormB(qM))
    {
        x_inf = -7.;
        x_sup = MIN(NormMapinv(1./beta*(u-NormB(qZ)*sqrt(1.-beta*beta)) , qM), 7.);
    }
    else
    {
        x_inf = -7;
        x_sup = 7;
    }

    for(k=0;k<nbPoints;k++)
    {
        x[k] = x_inf + k * (x_sup - x_inf) /nbPoints;
        x[k] = NormMap(x[k],qM);
    }


    // START INTEGRATION
    result = Integrand(u,x[0],beta,qM,qZ)*(x[1]-x[0]);
    for (k=1;k<nbPoints-1;k++)
    {
        result += Integrand(u,x[k],beta,qM,qZ)*(x[k+1] - x[k-1]);
    }
    result += Integrand(u,x[nbPoints-1],beta,qM,qZ)*(x[nbPoints-1]-x[nbPoints-2]);
    result *= 0.5;
    // END INTEGRATION
    */
    //result = y ;

    status = SUCCESS;

RETURN:
    //if(x) free(x);
    if(Fparam) free(Fparam);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return MAX(MIN(result,1.),0.);
}


/** root solving implementation */
typedef struct {double y; double beta; double qM; double qZ; long nbPoints; double eps;} ParamStruct;
static int FBrent(double x, void *data, double *out)
{
    ParamStruct *p = (ParamStruct*)data;
    *out = F(x, p->beta,p->qM,p->qZ,p->nbPoints,p->eps,NULL) - p->y;
    return 0;
}



/* -------------------------------------------------------------------------
** Finv
** this function solves the equation F(u) = p with input p
*/
double Finv(double y,
            double beta,
            double qM,
            double qZ,
            long nbPoints,
            double eps)
{
    double x;
    ParamStruct p;
    if (y < 3e-16) 
    {
        if(fabs(qM)>3e-16)
        {
            return NormB(qM);
        }
        else
        {
            return -3e100;
        }
    }
    else
    {
        if (1.-y < 3e-16) return 3e100;
        p.y = y;
        p.beta = beta;
        p.qM = qM;
        p.qZ = qZ;
        p.nbPoints = nbPoints;
        p.eps = eps;

        RootFindBrent(&FBrent, &p, -1e6, 1e6, 
            1000, 0., 1., 0., 1e-8, 1e-8, &x);
    }
    return x;
}


/* -------------------------------------------------------------------------
** QMultivariateDeviates_tp
** returns beta correlated q random variables
*/
int QMultivariateDeviates_tp(
    double *qSequence,              /* (O) [nbPaths*nbNames] */
    double *weight,                 /* (O) weight[nbPaths] */ 
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double qM,                      /* (I) */
    double qZ,                      /* (I) */
    long seed)                      /* (I) */
{
    static char routine[] = "QMultivariateDeviates_tp";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *sqrt_beta = NULL;
    double Mnormalisation = Normalisation(qM);
    double Znormalisation = Normalisation(qZ);

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
    status = CreateGaussianRandomSequence(M,seed,1,nbPaths);
    if(status == FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            qSequence[i+j*nbNames] =
                (beta[i]*Map(M[j],qM)*Mnormalisation + sqrt_beta[i]*Map(Z[i+j*nbNames],qZ)*Znormalisation);
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
** QCopulatedUniformDeviates
** returns beta correlated q unif random variables
*/
int QCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) [nbPaths*nbNames] */
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double qM,                      /* (I) */
    double qZ)                      /* (I) */
{
    static char routine[] = "QCopulatedUniformDeviates";
    int status = FAILURE;
    long i = 0;
    long j = 0;

    /* allocation of the Sequences */
    double *Z = NULL;
    double *M = NULL;
    double *sqrt_beta = NULL;
    double Mnormalisation = Normalisation(qM);
    double Znormalisation = Normalisation(qZ);

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
    status = CreateGaussianRandomSequence(Z,-7,nbNames,nbPaths);
    if(status == FAILURE) goto RETURN;   
    
    /* generate the market variable with gaussian random variable */
    status = CreateGaussianRandomSequence(M,-15,1,nbPaths);
    if(status == FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedUniformDeviates[i+j*nbNames] =
                F((beta[i]*Map(M[j],qM)*Mnormalisation + sqrt_beta[i]*Map(Z[i+j*nbNames],qZ)*Znormalisation), beta[i], qM, qZ, 100, 1e-10, NULL);
        }
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
** QCopulatedIndicator
** this function returns the indicators of default using a q-copula
*/
int QCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double qM,                      /* (I) */
    double qZ,                      /* (I) */
    long nbPoints,
    double eps,
    long seed)
{
    static char routine[] = "QCopulatedIndicator";
    int status = FAILURE;
    int i,j;
    double *qSequence = NULL;
    double *T = NULL;

    qSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(qSequence == NULL) goto RETURN;
    T = malloc(nbNames*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames;i++)
    {
        T[i] = Finv(survivalProba[i],beta[i],qM,qZ,nbPoints, eps);
    }

    /* allocate the q sequence */
    status = QMultivariateDeviates_tp(qSequence, weight, nbNames, nbPaths, beta, qM,qZ, seed);
    if(status ==FAILURE) goto RETURN;
    

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            if(qSequence[i+j*nbNames] > T[i])
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
    if(qSequence) free(qSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

/* ---------------------------------------------------------------------------
// QCopulatedIndicator_mtp
// this function returns the survival indicator and the weight of each path
//
*/
int QCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    double qM,
    double qZ,
    long NbPoints,
    double eps,
    long seed)
{
    static char routine[] = "QCopulatedIndicator_mtp";
    int status = FAILURE;
    int i,j,k;
    double *qSequence = NULL;
    double *T = NULL;

    qSequence = malloc(nbNames*nbPaths*sizeof(double));
    if(qSequence == NULL) goto RETURN;
    T = malloc(nbNames*nbTimes*sizeof(double));
    if(T==NULL) goto RETURN;
    
    for(i=0;i<nbNames;i++)
    {
        for(j=0;j<nbTimes;j++)
        {
            T[i+j*nbNames] = Finv(survivalProba[i+j*nbNames],beta[i],qM,qZ,NbPoints, eps);
        }
    }

    /* allocate the gaussian sequence */
    status = QMultivariateDeviates_tp(qSequence, weight, nbNames, nbPaths, beta, qM,qZ, seed);
    if(status ==FAILURE) goto RETURN;

    for(j=0;j<nbPaths;j++)
    {
        for(i=0;i<nbNames;i++)
        {
            copulatedSurvivalIndicator[i+j*nbNames] = nbTimes;
            for(k=0;k<nbTimes;k++)
            {
                if(qSequence[i+j*nbNames] > T[k+i*nbTimes])
                {
                    copulatedSurvivalIndicator[i+j*nbNames] = k;
                    break;
                }
            }
        }
    }
    status = SUCCESS;
RETURN:
    if(qSequence) free(qSequence);
    if(T) free(T);
    if (status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);    
    }
    return status;
}

