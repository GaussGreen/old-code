#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include "rootbrent.h"
#include "skew_util.h"
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

typedef struct
{
    double u1;
    double u2;
    double qM;
    double qZ1;
    double qZ2;
    double beta1;
    double beta2;
} FPARAM2;

typedef struct
{
    double u1;
    double u2;
    double u3;
    double qM;
    double qZ1;
    double qZ2;
    double qZ3;
    double beta1;
    double beta2;
    double beta3;
} FPARAM3;

typedef struct
{
    double u;
    double qM;
    double qZ;
    double beta;
    long   n;
} MOMENT_PARAM;

/*****************************************************************************/
/*****************************************************************************/
/*  MODE  IN LARGE POOL       ************************************************/
/*****************************************************************************/

typedef struct {double y; double beta; double qM; double qZ;} ParamStruct;

double Loss( double K,
             double u,
             double beta,
             double qM,
             double qZ,
             double alpha,
             long nbPoints,
             double eps)
{
    double Finvu = Finv(pow(u,1.-alpha),beta,qM,qZ,nbPoints,eps);
    return (1-alpha)*IntegrandLn(Finvu,beta,qM,qZ,K,0) + ((K==1)?alpha:0);
}

double DLoss(double K,
             double u,
             double beta,
             double qM,
             double qZ,
             double alpha,
             long nbPoints,
             double eps)
{
    double y = Finv(pow(u,1.-alpha),beta,qM,qZ,nbPoints,eps);
    double sqr = sqrt(1.-beta*beta);
    double x   = Fqinv(1.-K,qZ);
    double ck = (y - sqr * x ) / beta;
    double temp1 = beta * f1q(x,qZ) *fq(ck,qM);
    double temp2 = sqr * f1q(ck,qM) *fq(x,qZ);
    return temp1 + temp2;
}
                

static int GBrent(double x, void *data, double *out)
{
    ParamStruct *p = (ParamStruct*)data;
    double qM,qZ, beta, y;
    double temp1, temp2;
    double sqr;
    double ck;
    qM = p->qM;
    qZ = p->qZ;
    beta = p->beta;
    y = p->y; /* = Finv(u) */
    sqr = sqrt(1.-beta*beta);
    ck = (y - sqr * x ) / beta;
    temp1 = beta * f1q(x,qZ) *fq(ck,qM);
    temp2 = sqr * f1q(ck,qM) *fq(x,qZ);
    *out = temp1 + temp2;
    return 0;
}
/* OLD CODE (a little) less clear
{
    ParamStruct *p = (ParamStruct*)data;
    double AqM, BqM, CqM, AqZ, BqZ, CqZ;
    double qM,qZ, beta, eps, y;
    long nbPoints;
    double temp1, temp2;
    qM = p->qM;
    qZ = p->qZ;
    beta = p->beta;
    nbPoints = p->nbPoints;
    eps = p->eps;
    y = p->y;
    
    BqM = SIGN(qM)*sqrt( exp(qM*qM)*(exp(qM*qM) - 1.));
    BqZ = SIGN(qZ)*sqrt( exp(qZ*qZ)*(exp(qZ*qZ) - 1.));
    CqM = 1.- qM*qM;
    CqZ = 1.- qZ*qZ;
    AqM = CqM / BqM;
    AqZ = CqZ / BqZ;

    temp1 = (Finv(y,beta,qM,qZ,nbPoints,eps) - sqrt(1.-beta*beta)*x)/beta;
    temp1 = 1./(temp1 + AqM) * (1 + 1./(qM*qM) * log( temp1*BqM + CqM ));
    temp2 = beta / sqrt(1.-beta*beta) * 1. / (x + AqZ) 
                * (1. + 1./(qZ*qZ) *log( x * BqZ + CqZ ));
    *out = temp1 + temp2;
    return 0;
}
*/

/* -------------------------------------------------------------------------
** this function returns the mode of the distribution
**
*/
double Mode(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps)
{
    double x;
    double res;
    ParamStruct p;

    //double infB = NormalCum(NormMapinv(0.5,qZ));
    //double supB =  1e300;
    double infB = Fqinv(0.5,qZ);
    //double infB = -1e300;
    double supB = Fqinv(0.9999,qZ);
    //double BqM, BqZ, CqM, CqZ;
    //double sqrt_beta2 = sqrt(1.-beta*beta);
    double ua = pow(u,1.-alpha);
    double Finvu = Finv(ua,beta,qM,qZ,nbPoints,eps);
    double guess = Fqinv(ua,qZ);
    p.y = Finvu;
    p.beta = beta;
    p.qM = qM;
    p.qZ = qZ;
/*
    BqM = SIGN(qM)*sqrt( exp(qM*qM)*(exp(qM*qM) - 1.));
    BqZ = SIGN(qZ)*sqrt( exp(qZ*qZ)*(exp(qZ*qZ) - 1.));
    CqM = 1.- qM*qM;
    CqZ = 1.- qZ*qZ;
    
    if(BqZ>0)
    {
        infB = MAX(infB, - CqZ/BqZ);
    }
    else if(BqZ<0)
    {
        supB = MIN(supB, -CqZ/BqZ);
    }
    
    if(BqM>0)
    {
        if(beta>0)
        {
            supB = MIN(supB, (Finvu +beta*CqM/BqM)/sqrt_beta2);
        }
        else if(beta<0)
        {
            infB = MAX(infB, (Finvu +beta*CqM/BqM)/sqrt_beta2);        
        }
    }
    else if(BqM<0)
    {
        if(beta>0)
        {
            infB = MAX(infB, (Finvu +beta*CqM/BqM)/sqrt_beta2);
        }
        else if(beta <0)
        {
            supB = MIN(supB, (Finvu +beta*CqM/BqM)/sqrt_beta2);        
        }
    }
    if(infB<0 && supB >0)
    {
        guess = 0.0;
    }
    else
    {
        if(supB<=0)
        {
            guess = supB-0.01;
        }
        else if(infB>=0)
        {
            guess = infB+0.01;
        }
    }
*/
    RootFindBrent(&GBrent, &p, /*infB+1e-6*/infB, /*supB-1e-6*/supB, 
            100, guess, 1e-1, 0., 1e-4, 1e-4, &x);
    
    res = 1. - Fq(x,qZ); /* actually equals K0/(1-R) */

    if(fabs(res-1.)<3e-3)
    {
        return 0.0;
    }
    else
    {
        return res;
    }
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of the mode
** with respect to qM, qZ, beta
*/
double DModeDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double mode0 = Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double mode1 = Mode(u,beta,qM+h,qZ,alpha,nbPoints,eps);
    return (mode1-mode0) / h;
}

double DModeDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double mode0 = Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double mode1 = Mode(u,beta,qM,qZ+h,alpha,nbPoints,eps);
    return (mode1-mode0) / h;
}

double DModeDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double mode0 = Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double mode1 = Mode(u,beta+h,qM,qZ,alpha,nbPoints,eps);
    return (mode1-mode0) / h;
}

double DModeDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double mode0 = Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double mode1 = Mode(u,beta,qM,qZ,alpha+h,nbPoints,eps);
    return (mode1-mode0) / h;
}

/* -------------------------------------------------------------------------
** this function returns the height of the mode of the distribution
**
*/
double FMode(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps)
{
    double mode = Mode(u,beta,qM,qZ,alpha,nbPoints,eps);
    if(fabs(mode)<3e-10) return 0.0;
    return Loss(mode,u,beta,qM,qZ,alpha,nbPoints,eps);
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of height of the the moment
** with respect to qM, qZ, beta
*/
double DFModeDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double fmode0 = FMode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double fmode1 = FMode(u,beta,qM+h,qZ,alpha,nbPoints,eps);
    return (fmode1-fmode0) / h;
}

double DFModeDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double fmode0 = FMode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double fmode1 = FMode(u,beta,qM,qZ+h,alpha,nbPoints,eps);
    return (fmode1-fmode0) / h;
}

double DFModeDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double fmode0 = FMode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double fmode1 = FMode(u,beta+h,qM,qZ,alpha,nbPoints,eps);
    return (fmode1-fmode0) / h;
}

double DFModeDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double fmode0 = FMode(u,beta,qM,qZ,alpha,nbPoints,eps);
    double fmode1 = FMode(u,beta,qM,qZ,alpha+h,nbPoints,eps);
    return (fmode1-fmode0) / h;
}

/***************************************************************************/
/***************************************************************************/
/* JOINT CUMULATIVE of 2 skewed variables                                  */
/***************************************************************************/
/* -------------------------------------------------------------------------
** Integrand2
** the integral of this function for M=-infinity to M=+infinity
** is equal to Fi(u,v) = P(Xi<u,Xj<v)
*/
double Integrand2(double u1,double u2, double M, double beta1, double beta2, double qM, double qZ1, double qZ2)
{
    double result;
    /*double res;*/
    double z1,z2;


    /*  the mapping chosen (q mapping) implies that
        the inverse mapping is not defined everywhere because
        of the log function it contains.
        We need that mapinv(M) and mapinv(z) (see below) are defined.
        The definition interval for M is then
        [B(qM)..1/beta*(u-B(qZ)*sqrt(1-beta^2))] (see B above)
        if beta > 0,
        [B(qM).. +infinity]
        otherwise
    */
        
    if(M*SIGN(qM) <= NormB(qM))
    {
        result = 0.0;
    }
    else
    {
        result = NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
    }

    if((u1-beta1*M)*SIGN(qZ1) <= NormB(qZ1) * sqrt(1.-beta1*beta1))
    {
        if(qZ1>=0)
        {
            result *= 0.0;
        }
        else
        {
            result *=1.0;
        }
    }
    else
    {
        z1 = (u1-beta1*M)/sqrt(1.-beta1*beta1);
        result *= NormalCum(NormMapinv(z1,qZ1));
    }

    if((u2-beta2*M)*SIGN(qZ2) <= NormB(qZ2) * sqrt(1.-beta2*beta2))
    {
        if(qZ2>=0)
        {
            result *= 0.0;
        }
        else
        {
            result *=1.0;
        }
    }
    else
    {
        z2 = (u2-beta2*M)/sqrt(1.-beta2*beta2);
        result *= NormalCum(NormMapinv(z2,qZ2));
    }

    /*
    if((u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai) ||
       (v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj) )
    {
        if(qZ1>=0 || qZ2>=0)
        {
            return 0.0;
        }
        else
        {
            if( (u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai) &&
                (v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj) )
            {
                return NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else if((u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai))
            {
                    if(fabs(betaj-1.)<3e-16)
                    {
                        /* we should never go there *
                        zj = (v - betaj*M)*1e300;
                    }
                    else
                    {
                        zj = (v - betaj*M)/sqrt(1-betaj*betaj);
                    }
                return  NormalCum(NormMapinv(zj,qZ2))*
                        NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else if((v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj))
            {
                    if(fabs(betai-1.)<3e-16)
                    {
                        /* we should never go there *
                        zi = (u - betai*M)*1e300;
                    }
                    else
                    {
                        zi = (u - betai*M)/sqrt(1-betai*betai);
                    }


                return  NormalCum(NormMapinv(zi,qZ1))*
                        NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else
            {
                return 0.0;
                /* never ever *
            }
        }

    }

    if(fabs(betai-1.)<3e-16)
    {
        /* we should never go there *
        zi = (u - betai*M)*1e300;
    }
    else
    {
        zi = (u - betai*M)/sqrt(1-betai*betai);
    }
    if(fabs(betaj-1.)<3e-16)
    {
        /* we should never go there *
        zj = (v - betaj*M)*1e300;
    }
    else
    {
        zj = (v - betaj*M)/sqrt(1-betaj*betaj);
    }

    res = NormalCum(NormMapinv(zi,qZ1))*NormalCum(NormMapinv(zj,qZ2))
           *NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
           */
    return result;
}

/* function to be fed to the GSL adaptive integration method */
int derivs_gsl2(double x, double *y, void *param)
{
    FPARAM2 *Fparam2 = (FPARAM2*) param;
    *y = Integrand2(Fparam2->u1,Fparam2->u2,x,Fparam2->beta1,Fparam2->beta2,Fparam2->qM,Fparam2->qZ1,Fparam2->qZ2);
    return 0;
}


/* -------------------------------------------------------------------------
** Fij
** this function computes the integral of the previous function 
** naive trapezoidal integration
**
*
double Fij( double u,
            double v,
            double betai,
            double betaj,
            double qM,
            double qZ,
            long nbPoints)
{
    static char routine[] = "Fij";
    int status = FAILURE;
    double result = 0.0;
    double *x = NULL;
    long k;
    double x_inf, x_sup;

    x = malloc(nbPoints*sizeof(double)); //these are the abscissas 
    if(x==NULL) goto RETURN;
    
    /* if beta = 1, Xi = M */
    /*
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
    */
    /* if beta = 0, Xi = Zi */
    /*
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
    */


    /* general case, beta != 0 and beta != 1 *
    if(qM>3e-16 && qZ>3e-16 && (1./betai*(u-NormB(qZ)*sqrt(1.-betai*betai))>NormB(qM)) && (1./betaj*(u-NormB(qZ)*sqrt(1.-betaj*betaj))>NormB(qM)))
    {
        x_inf = -7.;
        x_sup = MIN(MIN(NormMapinv(1./betai*(u-NormB(qZ)*sqrt(1.-betai*betai)) , qM),NormMapinv(1./betaj*(u-NormB(qZ)*sqrt(1.-betaj*betaj)) , qM)), 7.);
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

    result = IntegrandDouble(u,v,x[0],betai,betaj,qM,qZ)*(x[1]-x[0]);
    for (k=1;k<nbPoints-1;k++)
    {
        result += IntegrandDouble(u,v,x[k],betai,betaj,qM,qZ)*(x[k+1] - x[k-1]);
    }
    result += IntegrandDouble(u,v,x[nbPoints-1],betai,betaj,qM,qZ)*(x[nbPoints-1]-x[nbPoints-2]);
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

/* -------------------------------------------------------------------------
** F2
** this function computes the integral over M of the Integrand2 function 
**
*/
double F2(  double u1,
            double u2,
            double beta1,
            double beta2,
            double qM,
            double qZ1,
            double qZ2,
            long NbPoints,
            double eps
            )
{
    static char routine[] = "F2";
    int status = FAILURE;
    double result = 0.0;

    FPARAM2 *Fparam2 = NULL;
    Fparam2 = malloc(sizeof(FPARAM2));
    if(Fparam2==NULL) goto RETURN;


    Fparam2->beta1 = beta1;
    Fparam2->beta2 = beta2;
    Fparam2->qM = qM;
    Fparam2->qZ1 = qZ1;
    Fparam2->qZ2 = qZ2;
    Fparam2->u1 = u1;
    Fparam2->u2 = u2;

    /* if beta1 = 1 and beta2 = 1 */
    if(fabs(beta1-1.)<3e-16 && fabs(beta2-1.)<3e-16 )
    {
        result = NormalCum(NormMapinv(MIN(u1,u2),qM));    
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 1 and beta2 != 1 */
    if(fabs(beta1-1.)<3e-16 && fabs(beta2-1.)>3e-16 )
    {    
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);

        FPARAM *Fparam = NULL;
        Fparam = malloc(sizeof(FPARAM));
        Fparam->beta = beta2;
        Fparam->qM = qM;
        Fparam->qZ = qZ2;
        Fparam->u = u2;

        gsl_integration_qagil (&derivs_gsl,u1, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        gsl_integration_workspace_free(workspace);        
        if(Fparam) free(Fparam);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 != 1 and beta2 = 1 */
    if(fabs(beta1-1.)>3e-16 && fabs(beta2-1.)<3e-16 )
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);

        FPARAM *Fparam = NULL;
        Fparam = malloc(sizeof(Fparam));
        Fparam->beta = beta1;
        Fparam->qM = qM;
        Fparam->qZ = qZ1;
        Fparam->u = u1;

        gsl_integration_qagil (&derivs_gsl,u2, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        gsl_integration_workspace_free(workspace);        
        if(Fparam) free(Fparam);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 0 and beta2 != 0 */
    if(fabs(beta1)<3e-16 && fabs(beta2)>3e-16)
    {
        result = NormalCum(NormMapinv(u1,qZ1))*F(u2,beta2,qM,qZ2,NbPoints, eps, NULL);    
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 != 0 and beta2 = 0 */
    if(fabs(beta1)>3e-16 && fabs(beta2)<3e-16)
    {
        result = NormalCum(NormMapinv(u2,qZ2))*F(u1,beta1,qM,qZ1,NbPoints, eps, NULL);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 0 and beta2 = 0 */
    if(fabs(beta1)<3e-16 && fabs(beta2)<3e-16)
    {
        return NormalCum(NormMapinv(u1,qZ1))*NormalCum(NormMapinv(u2,qZ2));
        status = SUCCESS;
        goto RETURN;
    }



    /* general case, beta != 0 and beta != 1 */
    /********************************************************************/
    /* adaptive GSL method
    */
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);
        //if(x_sup>=10 && x_inf<= -10)
        //{
            gsl_integration_qagi (&derivs_gsl2, eps, 0, NbPoints, workspace, &result, &error, Fparam2);
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
    
    status = SUCCESS;

RETURN:
    if(Fparam2) free(Fparam2);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return MAX(MIN(result,1.),0.);
}


/***************************************************************************/
/***************************************************************************/
/* JOINT CUMULATIVE of 3 skewed variables                                    */
/***************************************************************************/
/* -------------------------------------------------------------------------
** Integrand3
** the integral of this function for M=-infinity to M=+infinity
** is equal to Fi(u,v) = P(Xi<u,Xj<v)
*/
double Integrand3(double u1,double u2,double u3, double M, double beta1, double beta2, double beta3, double qM, double qZ1, double qZ2, double qZ3)
{
    double result;
    /*double res;*/
    double z1,z2,z3;


    /*  the mapping chosen (q mapping) implies that
        the inverse mapping is not defined everywhere because
        of the log function it contains.
        We need that mapinv(M) and mapinv(z) (see below) are defined.
        The definition interval for M is then
        [B(qM)..1/beta*(u-B(qZ)*sqrt(1-beta^2))] (see B above)
        if beta > 0,
        [B(qM).. +infinity]
        otherwise
    */
        
    if(M*SIGN(qM) <= NormB(qM))
    {
        result = 0.0;
    }
    else
    {
        result = NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
    }

    if((u1-beta1*M)*SIGN(qZ1) <= NormB(qZ1) * sqrt(1.-beta1*beta1))
    {
        if(qZ1>=0)
        {
            result *= 0.0;
        }
        else
        {
            result *=1.0;
        }
    }
    else
    {
        z1 = (u1-beta1*M)/sqrt(1.-beta1*beta1);
        result *= NormalCum(NormMapinv(z1,qZ1));
    }

    if((u2-beta2*M)*SIGN(qZ2) <= NormB(qZ2) * sqrt(1.-beta2*beta2))
    {
        if(qZ2>=0)
        {
            result *= 0.0;
        }
        else
        {
            result *=1.0;
        }
    }
    else
    {
        z2 = (u2-beta2*M)/sqrt(1.-beta2*beta2);
        result *= NormalCum(NormMapinv(z2,qZ2));
    }


    if((u3-beta3*M)*SIGN(qZ3) <= NormB(qZ3) * sqrt(1.-beta3*beta3))
    {
        if(qZ3>=0)
        {
            result *= 0.0;
        }
        else
        {
            result *=1.0;
        }
    }
    else
    {
        z3 = (u3-beta3*M)/sqrt(1.-beta3*beta3);
        result *= NormalCum(NormMapinv(z3,qZ3));
    }

    /*
    if((u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai) ||
       (v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj) )
    {
        if(qZ1>=0 || qZ2>=0)
        {
            return 0.0;
        }
        else
        {
            if( (u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai) &&
                (v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj) )
            {
                return NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else if((u-betai*M)*SIGN(qZ1) <= NormB(qZ1)*sqrt(1.-betai*betai))
            {
                    if(fabs(betaj-1.)<3e-16)
                    {
                        /* we should never go there *
                        zj = (v - betaj*M)*1e300;
                    }
                    else
                    {
                        zj = (v - betaj*M)/sqrt(1-betaj*betaj);
                    }
                return  NormalCum(NormMapinv(zj,qZ2))*
                        NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else if((v-betaj*M)*SIGN(qZ2) <= NormB(qZ2)*sqrt(1.-betaj*betaj))
            {
                    if(fabs(betai-1.)<3e-16)
                    {
                        /* we should never go there *
                        zi = (u - betai*M)*1e300;
                    }
                    else
                    {
                        zi = (u - betai*M)/sqrt(1-betai*betai);
                    }


                return  NormalCum(NormMapinv(zi,qZ1))*
                        NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
            }
            else
            {
                return 0.0;
                /* never ever *
            }
        }

    }

    if(fabs(betai-1.)<3e-16)
    {
        /* we should never go there *
        zi = (u - betai*M)*1e300;
    }
    else
    {
        zi = (u - betai*M)/sqrt(1-betai*betai);
    }
    if(fabs(betaj-1.)<3e-16)
    {
        /* we should never go there *
        zj = (v - betaj*M)*1e300;
    }
    else
    {
        zj = (v - betaj*M)/sqrt(1-betaj*betaj);
    }

    res = NormalCum(NormMapinv(zi,qZ1))*NormalCum(NormMapinv(zj,qZ2))
           *NormalDensity(NormMapinv(M,qM)) * DNormMapinv(M,qM);
           */
    return result;
}

/* function to be fed to the GSL adaptive integration method */
int derivs_gsl3(double x, double *y, void *param)
{
    FPARAM3 *Fparam3 = (FPARAM3*) param;
    *y = Integrand3(Fparam3->u1,Fparam3->u2,Fparam3->u3,x,Fparam3->beta1,Fparam3->beta2,Fparam3->beta3,Fparam3->qM,Fparam3->qZ1,Fparam3->qZ2, Fparam3->qZ3);
    return 0;
}

/* -------------------------------------------------------------------------
** F3
** this function computes the integral of the Integrand3 function 
**
*/
double F3(  double u1,
            double u2,
            double u3,
            double beta1,
            double beta2,
            double beta3,
            double qM,
            double qZ1,
            double qZ2,
            double qZ3,
            long NbPoints,
            double eps
            )
{
    static char routine[] = "F3";
    int status = FAILURE;
    double result = 0.0;

    FPARAM3 *Fparam3 = NULL;
    Fparam3 = malloc(sizeof(FPARAM3));
    if(Fparam3==NULL) goto RETURN;


    Fparam3->beta1 = beta1;
    Fparam3->beta2 = beta2;
    Fparam3->beta3 = beta3;
    Fparam3->qM = qM;
    Fparam3->qZ1 = qZ1;
    Fparam3->qZ2 = qZ2;
    Fparam3->qZ3 = qZ3;
    Fparam3->u1 = u1;
    Fparam3->u2 = u2;
    Fparam3->u3 = u3;

    /* if beta1 = 1 and beta2 = 1 */
    if(fabs(beta1-1.)<3e-16 && fabs(beta2-1.)<3e-16 )
    {
        result = NormalCum(NormMapinv(MIN(u1,u2),qM));    
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 1 and beta2 != 1 */
    if(fabs(beta1-1.)<3e-16 && fabs(beta2-1.)>3e-16 )
    {    
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);

        FPARAM *Fparam = NULL;
        Fparam = malloc(sizeof(FPARAM));
        Fparam->beta = beta2;
        Fparam->qM = qM;
        Fparam->qZ = qZ2;
        Fparam->u = u2;

        gsl_integration_qagil (&derivs_gsl,u1, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        gsl_integration_workspace_free(workspace);        
        if(Fparam) free(Fparam);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 != 1 and beta2 = 1 */
    if(fabs(beta1-1.)>3e-16 && fabs(beta2-1.)<3e-16 )
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);

        FPARAM *Fparam = NULL;
        Fparam = malloc(sizeof(Fparam));
        Fparam->beta = beta1;
        Fparam->qM = qM;
        Fparam->qZ = qZ1;
        Fparam->u = u1;

        gsl_integration_qagil (&derivs_gsl,u2, eps, 0, NbPoints, workspace, &result, &error, Fparam);
        gsl_integration_workspace_free(workspace);        
        if(Fparam) free(Fparam);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 0 and beta2 != 0 */
    if(fabs(beta1)<3e-16 && fabs(beta2)>3e-16)
    {
        result = NormalCum(NormMapinv(u1,qZ1))*F(u2,beta2,qM,qZ2,NbPoints, eps, NULL);    
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 != 0 and beta2 = 0 */
    if(fabs(beta1)>3e-16 && fabs(beta2)<3e-16)
    {
        result = NormalCum(NormMapinv(u2,qZ2))*F(u1,beta1,qM,qZ1,NbPoints, eps, NULL);
        status = SUCCESS;
        goto RETURN;
    }

    /* if beta1 = 0 and beta2 = 0 */
    if(fabs(beta1)<3e-16 && fabs(beta2)<3e-16)
    {
        return NormalCum(NormMapinv(u1,qZ1))*NormalCum(NormMapinv(u2,qZ2));
        status = SUCCESS;
        goto RETURN;
    }



    /* general case, beta != 0 and beta != 1 */
    /********************************************************************/
    /* adaptive GSL method
    */
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(NbPoints);
        //if(x_sup>=10 && x_inf<= -10)
        //{
            gsl_integration_qagi (&derivs_gsl3, eps, 0, NbPoints, workspace, &result, &error, Fparam3);
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
    
    status = SUCCESS;

RETURN:
    if(Fparam3) free(Fparam3);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return MAX(MIN(result,1.),0.);
}

/*****************************************************************************/
/*****************************************************************************/
/*  n th moment LARGE POOL                                            */
/*****************************************************************************/
/*******************************************/
/* IntegrandLn                             */
/* this functions returns K^n * Pr(L in dK)*/
/*******************************************/
double IntegrandLn( double u, /* = Finv(u) */
                    double beta,
                    double qM,
                    double qZ,
                    double K,
                    long n
                    )
{
    double sqr      = sqrt(1.-beta*beta);
    double Fzinv    = Fqinv(1.-K,qZ);
    double z        = (u - sqr * Fzinv) / beta;
    double fm       = fq( z, qM);
    double res      = sqr*fm / (beta*fq(Fzinv,qZ));

    if(K==0) return 0.0;
    if(K==1) return 0.0;

    switch(n)
    {
    case 0:
        return res;
        break;
    case 1:
        return K*res;
        break;
    case 2:
        return K*K*res;
        break;
    case 3:
        return K*K*K*res;
        break;
    default:
        return pow(K,n)*res;
        break;
    }
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** Skew
** with respect to qM, qZ, beta
*/
double DIntegrandLnDqM(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h)
{
    double Finvu = Finv(u,beta,qM,qZ,nbPoints,eps);
    double int0 = IntegrandLn(Finvu,beta,qM,qZ,K,n);
    double int1 = IntegrandLn(Finvu,beta,qM+h,qZ,K,n);
    return (int1-int0) / h;
}

double DIntegrandLnDqZ(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h)
{
    double Finvu = Finv(u,beta,qM,qZ,nbPoints,eps);
    double int0 = IntegrandLn(Finvu,beta,qM,qZ,K,n);
    double int1 = IntegrandLn(Finvu,beta,qM,qZ+h,K,n);
    return (int1-int0) / h;
}

double DIntegrandLnDbeta(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h)
{
    double Finvu = Finv(u,beta,qM,qZ,nbPoints,eps);
    double int0 = IntegrandLn(Finvu,beta,qM,qZ,K,n);
    double int1 = IntegrandLn(Finvu,beta+h,qM,qZ,K,n);
    return (int1-int0) / h;
}

/* function to be fed to the GSL adaptive integration method */
int derivs_gsl_moment(double x, double *y, void *param)
{
    MOMENT_PARAM *moment_param = (MOMENT_PARAM*) param;
    *y = IntegrandLn(moment_param->u,moment_param->beta,moment_param->qM,moment_param->qZ,x, moment_param->n);
    return 0;
}

/********************************************/
/* n th moment (integral of the previous fct*/
/* over [0,1] )                             */ 
/********************************************/

double Moment(      double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long n,
                    long nbPoint,
                    double eps
                    )
{
    static char routine[] = "Moment";
    int status      = FAILURE;
    double result   = 0.0;
    MOMENT_PARAM *moment_param = NULL;
    moment_param = malloc(sizeof(MOMENT_PARAM));
    if(moment_param==NULL) goto RETURN;

    moment_param->u     = Finv(pow(u,1.-alpha),beta,qM,qZ,nbPoint,eps);
    moment_param->qM    = qM;
    moment_param->qZ    = qZ;
    moment_param->beta  = beta;
    moment_param->n     = n;

    /* general case, beta != 0 and beta != 1 */
    /********************************************************************/
    /* adaptive GSL method
    */
    {
        double error;
        gsl_integration_workspace *workspace = gsl_integration_workspace_alloc(nbPoint);
        gsl_integration_qag (&derivs_gsl_moment,0.,1., eps, 0, nbPoint, 6, workspace, &result, &error, moment_param);
        gsl_integration_workspace_free(workspace);
    }
    
    status = SUCCESS;

RETURN:
    if(moment_param) free(moment_param);
    if(status == FAILURE)
    {
        DR_Error("%s:%d %s failed", __FILE__, __LINE__, routine);
    }
    return result*(1.-alpha) + alpha;
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** the moments
** with respect to qM, qZ, beta
*/
double DMomentDqM(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double moment0 = Moment(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double moment1 = Moment(u,beta,qM+h,qZ,alpha,n,nbPoints,eps);
    return (moment1-moment0) / h;
}

double DMomentDqZ(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double moment0 = Moment(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double moment1 = Moment(u,beta,qM,qZ+h,alpha,n,nbPoints,eps);
    return (moment1-moment0) / h;
}

double DMomentDbeta(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double moment0 = Moment(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double moment1 = Moment(u,beta+h,qM,qZ,alpha,n,nbPoints,eps);
    return (moment1-moment0) / h;
}

double DMomentDalpha(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double moment0 = Moment(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double moment1 = Moment(u,beta,qM,qZ,alpha+h,n,nbPoints,eps);
    return (moment1-moment0) / h;
}

/****************************************************************************/
/*  Kn      
*/
double Kn(          double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long n,
                    long nbPoint,
                    double eps
                    )
{
    double ELn = Moment(u,beta,qM,qZ,alpha,n,nbPoint,eps);
    double EL  = 1.-u;
    double EL2,EL3;
    switch(n)
    {
    case 0:
        return ELn - 1.;
        break;
    case 1:
        return ELn - EL;
        break;
    case 2:
        return ELn - EL*EL;
        break;
    case 3:
        EL2 = Moment(u,beta,qM,qZ,alpha,2,nbPoint,eps);
        return ELn - 3*EL*EL2 + 2*EL*EL*EL;
        break;
    case 4:
        EL3 = Moment(u,beta,qM,qZ,alpha,3,nbPoint,eps);
        EL2 = Moment(u,beta,qM,qZ,alpha,2,nbPoint,eps);
        return ELn - 4*EL3*EL + 6*EL2*EL*EL - 3*EL*EL*EL*EL;
        break;
    default:
        /*  not the real Kn, a dodgy normalisation it is
            to substract the power of the expected value.
            The formula is actually 
            \sum_{k=0}^n C_n^k E[X^k] E[X]^{n-k} */
        return ELn - pow(EL,n);
        break;
    }
}


/* -------------------------------------------------------------------------
** Skew
**
**
*/
double Skew(        double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    )
{
    double EL3 = Moment(u,beta,qM,qZ,alpha,3,nbPoint,eps);
    double EL2 = Moment(u,beta,qM,qZ,alpha,2,nbPoint,eps);
    double EL  = 1.-u;
    double k3  = EL3 - 3*EL*EL2 + 2*EL*EL*EL;
    double k2  = EL2 - EL*EL;
    return k3 / (k2*sqrt(k2));
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** Skew
** with respect to qM, qZ, beta
*/
double DSkewDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double skew0 = Skew(u,beta,qM,qZ,alpha,nbPoints,eps);
    double skew1 = Skew(u,beta,qM+h,qZ,alpha,nbPoints,eps);
    return (skew1-skew0) / h;
}

double DSkewDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double skew0 = Skew(u,beta,qM,qZ,alpha,nbPoints,eps);
    double skew1 = Skew(u,beta,qM,qZ+h,alpha,nbPoints,eps);
    return (skew1-skew0) / h;
}

double DSkewDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double skew0 = Skew(u,beta,qM,qZ,alpha,nbPoints,eps);
    double skew1 = Skew(u,beta+h,qM,qZ,alpha,nbPoints,eps);
    return (skew1-skew0) / h;
}

double DSkewDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double skew0 = Skew(u,beta,qM,qZ,alpha,nbPoints,eps);
    double skew1 = Skew(u,beta,qM,qZ,alpha+h,nbPoints,eps);
    return (skew1-skew0) / h;
}

/* -------------------------------------------------------------------------
** Kurtosis
**
**
*/
double Kurtosis(    double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    )
{
    double EL4 = Moment(u,beta,qM,qZ,alpha,4,nbPoint,eps);
    double EL3 = Moment(u,beta,qM,qZ,alpha,3,nbPoint,eps);
    double EL2 = Moment(u,beta,qM,qZ,alpha,2,nbPoint,eps);
    double EL  = 1.-u;
    double k4  = EL4 - 4*EL3*EL + 6*EL2*EL*EL - 3*EL*EL*EL*EL;
    double k2  = EL2 - EL*EL;
    return k4 / (k2*k2);
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** Sigma
** with respect to qM, qZ, beta
*/
double DKurtosisDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double k0 = Kurtosis(u,beta,qM,qZ,alpha,nbPoints,eps);
    double k1 = Kurtosis(u,beta,qM+h,qZ,alpha,nbPoints,eps);
    return (k1-k0) / h;
}

double DKurtosisDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double k0 = Kurtosis(u,beta,qM,qZ,alpha,nbPoints,eps);
    double k1 = Kurtosis(u,beta,qM,qZ+h,alpha,nbPoints,eps);
    return (k1-k0) / h;
}

double DKurtosisDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double k0 = Kurtosis(u,beta,qM,qZ,alpha,nbPoints,eps);
    double k1 = Kurtosis(u,beta+h,qM,qZ,alpha,nbPoints,eps);
    return (k1-k0) / h;
}

double DKurtosisDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double k0 = Kurtosis(u,beta,qM,qZ,alpha,nbPoints,eps);
    double k1 = Kurtosis(u,beta,qM,qZ,alpha+h,nbPoints,eps);
    return (k1-k0) / h;
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** Kn
** with respect to qM, qZ, beta
*/
double DKnDqM(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double kn0 = Kn(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double kn1 = Kn(u,beta,qM+h,qZ,alpha,n,nbPoints,eps);
    return (kn1-kn0) / h;
}

double DKnDqZ(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double kn0 = Kn(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double kn1 = Kn(u,beta,qM,qZ+h,alpha,n,nbPoints,eps);
    return (kn1-kn0) / h;
}

double DKnDbeta(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double kn0 = Kn(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double kn1 = Kn(u,beta+h,qM,qZ,alpha,n,nbPoints,eps);
    return (kn1-kn0) / h;
}

double DKnDalpha(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h)
{
    double kn0 = Kn(u,beta,qM,qZ,alpha,n,nbPoints,eps);
    double kn1 = Kn(u,beta,qM,qZ,alpha+h,n,nbPoints,eps);
    return (kn1-kn0) / h;
}


/****************************************************************************/
/*  Sigma
*/
double Sigma(       double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    )
{
    double EL2 = Moment(u,beta,qM,qZ,alpha,2,nbPoint,eps);
    double EL  = 1. - u;
    return sqrt(EL2- EL*EL);
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** Sigma
** with respect to qM, qZ, beta
*/
double DSigmaDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double sigma0 = Sigma(u,beta,qM,qZ,alpha,nbPoints,eps);
    double sigma1 = Sigma(u,beta,qM+h,qZ,alpha,nbPoints,eps);
    return (sigma1-sigma0) / h;
}

double DSigmaDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double sigma0 = Sigma(u,beta,qM,qZ,alpha,nbPoints,eps);
    double sigma1 = Sigma(u,beta,qM,qZ+h,alpha,nbPoints,eps);
    return (sigma1-sigma0) / h;
}

double DSigmaDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double sigma0 = Sigma(u,beta,qM,qZ,alpha,nbPoints,eps);
    double sigma1 = Sigma(u,beta+h,qM,qZ,alpha,nbPoints,eps);
    return (sigma1-sigma0) / h;
}

double DSigmaDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h)
{
    double sigma0 = Sigma(u,beta,qM,qZ,alpha,nbPoints,eps);
    double sigma1 = Sigma(u,beta,qM,qZ,alpha+h,nbPoints,eps);
    return (sigma1-sigma0) / h;
}

/* -------------------------------------------------------------------------
** quantile functions
** 
*/
double SeniorProba(double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    double K,
                    long nbPoint,
                    double eps
                    )
{
    double ua = pow(u,1.-alpha);
    double Finvu = Finv(ua,beta,qM,qZ,nbPoint,eps);
    return (1.-alpha)*(1.- Fq((Finvu - sqrt(1.-beta*beta)*Fqinv(1-K,qZ))/beta,qM)) + alpha;
}

double Quantile(    double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    double proba,
                    long nbPoint,
                    double eps
                    )
{
    double ua = pow(u,1.-alpha);
    double Finvu = Finv(ua,beta,qM,qZ,nbPoint,eps);
    return 1. - Fq((Finvu-beta*Fqinv(1.+(alpha-proba)/(1.-alpha),qM))/sqrt(1.-beta*beta),qZ);
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** SeniorProba and Quantile
** with respect to qM, qZ, beta
*/
double DSeniorProbaDbeta(   double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double sp0 = SeniorProba(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double sp1 = SeniorProba(u,beta+h,qM,qZ,alpha,K,nbPoint,eps);
    return (sp1 - sp0) / h;
}

double DSeniorProbaDqM(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double sp0 = SeniorProba(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double sp1 = SeniorProba(u,beta,qM+h,qZ,alpha,K,nbPoint,eps);
    return (sp1 - sp0) / h;
}

double DSeniorProbaDqZ(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double sp0 = SeniorProba(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double sp1 = SeniorProba(u,beta,qM,qZ+h,alpha,K,nbPoint,eps);
    return (sp1 - sp0) / h;
}

double DSeniorProbaDalpha(  double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double sp0 = SeniorProba(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double sp1 = SeniorProba(u,beta,qM,qZ,alpha+h,K,nbPoint,eps);
    return (sp1 - sp0) / h;
}

double DQuantileDbeta(      double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double q0 = Quantile(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double q1 = Quantile(u,beta+h,qM,qZ,alpha,K,nbPoint,eps);
    return (q1 - q0) / h;
}

double DQuantileDalpha(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double q0 = Quantile(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double q1 = Quantile(u,beta,qM,qZ,alpha+h,K,nbPoint,eps);
    return (q1 - q0) / h;
}

double DQuantileDqM(        double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double q0 = Quantile(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double q1 = Quantile(u,beta,qM+h,qZ,alpha,K,nbPoint,eps);
    return (q1 - q0) / h;
}

double DQuantileDqZ(        double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        )
{
    double q0 = Quantile(u,beta,qM,qZ,alpha,K,nbPoint,eps);
    double q1 = Quantile(u,beta,qM,qZ+h,alpha,K,nbPoint,eps);
    return (q1 - q0) / h;
}
/* -------------------------------------------------------------------------
** Tranchelet price
**
*
double TrancheletPrice(double u,
                        double beta,
                        double qM,
                        double qZ,
                        double alpha,
                        double K,
                        long nbPoint,
                        double eps
                        )
{
    return SeniorProba(u,beta,qM,qZ,alpha,K,nbPoint,eps);
}

/* -------------------------------------------------------------------------
** these functions return an approximation of the derivative of
** the tranchelet price
** with respect to qM, qZ, beta
*
double DTrancheletPriceDbeta(   double u,
                                double beta,
                                double qM,
                                double qZ,
                                double K,
                                long nbPoint,
                                double eps,
                                double h
                            )
{
    double tp0 = TrancheletPrice(u,beta,qM,qZ,K,nbPoint,eps);
    double tp1 = TrancheletPrice(u,beta+h,qM,qZ,K,nbPoint,eps);
    return (tp1-tp0) / h;
}

double DTrancheletPriceDqM(   double u,
                                double beta,
                                double qM,
                                double qZ,
                                double K,
                                long nbPoint,
                                double eps,
                                double h
                            )
{
    double tp0 = TrancheletPrice(u,beta,qM,qZ,K,nbPoint,eps);
    double tp1 = TrancheletPrice(u,beta,qM+h,qZ,K,nbPoint,eps);
    return (tp1-tp0) / h;
}

double DTrancheletPriceDqZ(   double u,
                                double beta,
                                double qM,
                                double qZ,
                                double K,
                                long nbPoint,
                                double eps,
                                double h
                            )
{
    double tp0 = TrancheletPrice(u,beta,qM,qZ,K,nbPoint,eps);
    double tp1 = TrancheletPrice(u,beta,qM,qZ+h,K,nbPoint,eps);
    return (tp1-tp0) / h;
}
*/