#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "proba_utils.h"
#include "rootbrent.h"
#include "morphing.h"

#ifndef TINY 
#define TINY 3e-15
#endif

/****************************************************************************/
/* Dirac utils : DiracDensity                                               */
/****************************************************************************/
double DiracDensity(double x, double mu)
{
    return (x==mu)?1e99:0;
}

/****************************************************************************/
/* Dirac utils : DiracCum                                                   */
/****************************************************************************/
double DiracCum(double x, double mu)
{
    return ((x>=mu)?1:0);
}

/****************************************************************************/
/* variance of factor variable                                              */
/****************************************************************************/
double FactorSigma( FACTOR_DENSITY *f)
{
    long n = f->n;
    long i;
    double ex   = 0.0;
    double ex2  = 0.0;
    for(i=0;i<n;i++)
    {
        double w    = f->t[i].weight;
        double mu   = f->t[i].mu;
        double s    = f->t[i].sigma;
        ex += w*mu;
        ex2 += w*mu*mu;
        if(s>TINY)
        {
            ex2 += w*s*s;
        }        
    }
    return sqrt(ex2 - ex*ex);
}

/****************************************************************************/
/* F_X from the product of a single component M and a single componenet Z   */
/****************************************************************************/
double F_DiracM_x_DiracZ(   double x,
                            double beta,
                            double mu_fm,
                            double mu_fz)
{
    return ((x-beta*mu_fm)/sqrt(1.-beta*beta) > mu_fz)?1:0;
}

double F_DiracM_x_NormalZ(  double x,
                            double beta,
                            double mu_fm,
                            double mu_fz,
                            double sigma_fz)
{
    return NormalCum( ((x - beta*mu_fm)/sqrt(1.-beta*beta) - mu_fz )/sigma_fz);
}

double F_NormalM_x_DiracZ(  double x,
                            double beta,
                            double mu_fm,
                            double sigma_fm,
                            double mu_fz)
{
    return NormalCum( ((x - sqrt(1.-beta*beta)*mu_fz)/beta - mu_fm )/sigma_fm);
}

double F_NormalM_x_NormalZ( double x,
                            double beta,
                            double mu_fm,
                            double sigma_fm,
                            double mu_fz,
                            double sigma_fz)
{
    double beta2 = beta*beta;
    double sigma_fm2 = sigma_fm*sigma_fm;
    double sigma_fz2 = sigma_fz*sigma_fz;
    return NormalCum(   (x - beta*mu_fm - sqrt(1.-beta2)*mu_fz)
                        / sqrt((1.-beta2)*sigma_fz2 + beta2*sigma_fm2) );
}

/****************************************************************************/
/* F_X : Cumulative of X = beta * M + sqrt(1-beta2)^(1/2) * Z               */
/* multi_componenent M and Z                                                */
/****************************************************************************/
double F_X(double x, X_DENSITY *fx)
{
    long i,j;
    double beta = fx->beta;
    double result = 0.0;
    long n_fm = fx->fm.n;
    long n_fz = fx->fz.n;
    double sigmaM;
    double sigmaZ;
    if(fx->var_adjust==1)
    {
        sigmaM = FactorSigma(&(fx->fm));
        sigmaZ = FactorSigma(&(fx->fz));
    }
    else
    {
        sigmaM = 1.;
        sigmaZ = 1.;    
    }

    for(i=0;i<n_fm;i++)
    {
        for(j=0;j<n_fz;j++)
        {
            double mu_fm = fx->fm.t[i].mu;
            double mu_fz = fx->fz.t[j].mu;
            double sigma_fm = fx->fm.t[i].sigma;
            double sigma_fz = fx->fz.t[j].sigma;
            double w_fm = fx->fm.t[i].weight;
            double w_fz = fx->fz.t[j].weight;

            if(sigma_fm<TINY && sigma_fz<TINY) /* diracM*diracZ */
            {
                result += w_fm*w_fz * F_DiracM_x_DiracZ(x,beta,mu_fm/sigmaM,mu_fz/sigmaZ);        
            }
            else if(sigma_fm<TINY && sigma_fz>TINY) /* diracM*normalZ*/
            {
                result += w_fm*w_fz * F_DiracM_x_NormalZ(   x,beta,
                                                            mu_fm/sigmaM,
                                                            mu_fz/sigmaZ,
                                                            sigma_fz/sigmaZ);        
            }
            else if(sigma_fm>TINY && sigma_fz<TINY) /* normalM*diracZ */
            {
                result += w_fm*w_fz * F_NormalM_x_DiracZ(   x,beta,
                                                            mu_fm/sigmaM,
                                                            sigma_fm/sigmaM,
                                                            mu_fz)/sigmaZ;        
            }
            else if(sigma_fm>TINY && sigma_fz>TINY) /*normalM*normalZ */
            {
                result += w_fm*w_fz * F_NormalM_x_NormalZ(  x,beta,
                                                            mu_fm/sigmaM,
                                                            sigma_fm/sigmaM,
                                                            mu_fz/sigmaZ,
                                                            sigma_fz/sigmaZ);        
            }
            else
            {
                goto RETURN;
            }
        }
    }

RETURN:
return result;
}

/***************************************************************************/
/* inverse cumulative of X                                                 */
/***************************************************************************/
/** root solving implementation */
typedef struct{X_DENSITY *fx; double p;} XCUMINVPARAM;
static int XCumBrent(double x, void *data, double *out)
{
    XCUMINVPARAM *p = (XCUMINVPARAM*) data;
    *out = F_X(x, p->fx) - p->p;
    return 0;
}

double F_Xinv(double p, X_DENSITY *fx)
{
    double x;
    XCUMINVPARAM param;
    param.p = p;
    param.fx = fx;
    RootFindBrent(&XCumBrent, &param, -1e6, 1e6, 
            1000, 0., 1., 0., 1e-8, 1e-8, &x);
    return x;
}

/****************************************************************************/
/* density of a multi component factor M or Z                               */
/****************************************************************************/
double FactorDensity(double x, FACTOR_DENSITY *f)
{
    long i;
    double result = 0.0;
    double w,mu, sigma;
    long n = f->n;
    for (i=0;i<n;i++)
    {
        w = f->t[i].weight;
        mu = f->t[i].mu;
        sigma = f->t[i].sigma;
        if(sigma<TINY)
        {
            result += w*DiracDensity(x,mu);
        }
        else
        {
            double sigma_inv = 1./sigma;
            result += w*sigma_inv*NormalDensity((x-mu)*sigma_inv);
        }
    }
    return result;
}

/****************************************************************************/
/* cumulative of a multi component factor M or Z                            */
/****************************************************************************/
double FactorCum(double x, FACTOR_DENSITY *f)
{
    long i;
    double result = 0.0;
    double w,mu, sigma;
    long n = f->n;
    for (i=0;i<n;i++)
    {
        w = f->t[i].weight;
        mu = f->t[i].mu;
        sigma = f->t[i].sigma;
        if(sigma<TINY)
        {
            result += w*DiracCum(x,mu);
        }
        else
        {
            double sigma_inv = 1./sigma;
            result += w*NormalCum((x-mu)*sigma_inv);
        }
    }
    return result;
}

/****************************************************************************/
/* inverse cumulative of a multi component factor M or Z                    */
/****************************************************************************/
/** root solving implementation */
typedef struct{FACTOR_DENSITY *f; double p;} FACTCUMINVPARAM;
static int FactCumBrent(double x, void *data, double *out)
{
    FACTCUMINVPARAM *p = (FACTCUMINVPARAM*) data;
    *out = FactorCum(x, p->f) - p->p;
    return 0;
}

double FactorCuminv(double p, FACTOR_DENSITY *f)
{
    double x;
    double y = 0.0;
    FACTCUMINVPARAM param;
    param.p = p;
    param.f = f;
    if(RootFindBrent(&FactCumBrent, &param, -1e300, 1e300, 
            1000, 0., 1., 0., 1e-12, 1e-12, &x) == FAILURE)
    {
        RootFindBrent(&FactCumBrent, &param, -1e300, 1e300, 
            1000, 0., 1., 0., 1e-8, 1e-1, &x);

        y = FactorCum(x,f);
        y = FactorCum(x+1e-7,f);
        y = FactorCum(x-1e-7,f);
    }

    return x;

}

/****************************************************************************/
/* large pool loss density                                                  */
/****************************************************************************/
double LossDensity(double K, double u, X_DENSITY *fx)
{
    long var_adjust = fx->var_adjust;
    double sigmaM   = (var_adjust==0)?1.:FactorSigma(&(fx->fm));
    double sigmaZ   = (var_adjust==0)?1.:FactorSigma(&(fx->fz));
    double beta     = fx->beta;
    double Fzinv    = FactorCuminv(1.-K, &(fx->fz)) / sigmaZ;
    double Finvu    = F_Xinv(u,fx);
    double beta2    = beta*beta;
    double sqr      = sqrt(1.-beta2);
    double result   = sqr / (beta * sigmaZ * FactorDensity(Fzinv * sigmaZ ,&(fx->fz)));
    result         *= sigmaM * FactorDensity((Finvu - sqr*Fzinv)*sigmaM/beta,&(fx->fm));
    return result;
}

/****************************************************************************/
/* large pool tranchelet loss                                               */
/****************************************************************************/
double TrancheletLoss( double K, double u, X_DENSITY *fx)
{
    long var_adjust = fx->var_adjust;
    double sigmaM   = (var_adjust==0)?1.:FactorSigma(&(fx->fm));
    double sigmaZ   = (var_adjust==0)?1.:FactorSigma(&(fx->fz));
    double beta     = fx->beta;
    double Fzinv    = FactorCuminv(1.-K, &(fx->fz)) / sigmaZ;
    double Finvu    = F_Xinv(u,fx);
    double sqr      = sqrt(1.-beta*beta);
    return 1. - FactorCum( (Finvu - sqr*Fzinv)*sigmaM / beta,&(fx->fm));
}

double DpiDwM( double K0, double K1, double u, double beta, double muM, double sigmaM)
{
    double Ninvu = NormalCumInverse(u);
    double beta2 = beta*beta;
    double sqr = sqrt(1.-beta2);
    double a =  NormalCum( (Ninvu - beta*muM)/sqrt((1.-beta2)+beta2*sigmaM*sigmaM));
    double b0 = NormalDensity( (Ninvu - sqr*NormalCumInverse(1-K0))/beta);
    double b1 = NormalDensity( (Ninvu - sqr*NormalCumInverse(1-K1))/beta);
    return 1./(beta*a)*(b1-b0);
}

double DpiDwZ( double K0, double K1, double u, double beta, double muZ, double sigmaZ)
{
    double Ninvu = NormalCumInverse(u);
    double beta2 = beta*beta;
    double sqr = sqrt(1.-beta2);
    double a =  NormalCum( (Ninvu - sqr*muZ)/sqrt((1.-beta2)*sigmaZ*sigmaZ+beta2));
    double b0;
    double b1;
    double c0 = NormalDensity( (Ninvu - sqr*NormalCumInverse(1.-K1))/beta);
    double c1 = NormalDensity( (Ninvu - sqr*NormalCumInverse(1.-K0))/beta);
    if(sigmaZ>TINY)
    {
        b0 = NormalDensity( (NormalCumInverse(1-K0) - muZ)/sigmaZ);
        b1 = NormalDensity( (NormalCumInverse(1-K1) - muZ)/sigmaZ);
    }
    else
    {
        b0 = (NormalCumInverse(1-K0)>muZ)?1:0;
        b1 = (NormalCumInverse(1-K1)>muZ)?1:0;
    }
    return      (1/(beta*a) - sqr*sigmaZ / (b1*beta)) * c1
            -   (1/(beta*a) - sqr*sigmaZ / (b0*beta)) * c0;
}

double DpiDbeta(double K0, double K1, double u, double beta)
{
    double Ninvu = NormalCumInverse(u);
    double Ninv0 = NormalCumInverse(1.-K0);
    double Ninv1 = NormalCumInverse(1.-K1);
    double beta2 = beta*beta;
    double beta3 = beta2*beta;
    double sqr = sqrt(1.-beta2);
    double a0 = -Ninvu/beta2 + 1./ (beta3*sqrt(1./beta2-1.)) * Ninv0;
    double a1 = -Ninvu/beta2 + 1./ (beta3*sqrt(1./beta2-1.)) * Ninv1;
    double b0 = NormalDensity((Ninvu - sqr*Ninv0)/beta);
    double b1 = NormalDensity((Ninvu - sqr*Ninv1)/beta);
    return a1*b1-a0*b0;
}

/*
double F_MDiracAdjust( double T, double beta, double d)
{
    return NormalCum((T-beta*d)/sqrt(1.-beta*beta));
}

double F_MNormalAdjust( double T, double beta, double mu, double sigma)
{
    return NormalCum( (T-mu*beta) / sqrt(1. + beta*beta*(sigma*sigma-1.)) );
}

double F_MNormalCutAdjust( double T, double beta, double mu, double sigma, double M1, double M2)
{
    double *res1 = NULL;
    double *res2 = NULL;
    double d = sqrt( 1. + beta*beta*(sigma*sigma-1.) );
    double x1 = (mu - M1) / sigma;
    double x2 = (mu - M2) / sigma;
    double y = (T - mu*beta) / d;
    double r = -sigma*beta / d;
    BiNormalCum(x1,y,r,res1);
    BiNormalCum(x2,y,r,res2);
    return res1 - res2;
}

double F_ZDiracAdjust( double T, double beta, double d)
{
    double sqr_beta = sqrt(1.-beta*beta);
    return F_MDiracAdjust( T, sqr_beta, d);
}

double F_ZNormalAdjust( double T, double beta, double mu, double sigma)
{
    double sqr_beta = sqrt(1.-beta*beta);
    return F_MNormalAdjust( T, sqr_beta, mu, sigma);
}

double F_ZNormalCutAdjust( double T, double beta, double mu, double sigma, double Z1, double Z2)
{
    double sqr_beta = sqrt(1.-beta*beta);
    return F_MNormalCutAdjust( T, sqr_beta, mu, sigma, Z1, Z2);
}

double Fx_DiracAdjust(double m, double d)
{
    if(m>=d) return 1.0;
    return 0.0;
}

double Fx_NormalAdjust(double m, double mu, double sigma)
{
    return NormalCum((m-mu)/sigma);
}

double Fx_NormalCutAdjust(double m, double mu, double sigma, double M1, double M2)
{
    double *d = NULL;
    if(m<M1) return 0;
    if(m>=M2) return NormalCum((M2-mu)/sigma) - NormalCum((M1-mu)/sigma);
    
    //return NormalCum((m-mu)/sigma) - NormalCum((M1-mu)/sigma);
    BiNormalCum(m,m,mu,d);
    return *d;
}
*/