/* $Header$ */
#include "integration_nr.h"
#include "integration_gsl.h"

/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/


int QMultivariateDeviates_tp(
    double *qSequence,        /* (O) [nbPaths*nbNames] */
    double *weight,                 /* (O) weight[nbPaths] */ 
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double qM,                      /* (I) */
    double qZ,                      /* (I) */
    long seed);                     /* (I) */


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
    long seed);

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
    long nbPoints,
    double eps,
    long seed);


double Finv(double y,
            double beta,
            double qM,
            double qZ,
            long nbPoints,
            double eps);

double F(double u,
         double beta,
         double qM,
         double qZ,
         long NbPoints,
         double eps,
         DEBUGINFO *debugInfo
         );



double Integrand(double u,
                 double M,
                 double beta,
                 double qM,
                 double qZ);

double Mapinv(double y, double q);

double Map(double x, double q);

double DMapinv(double y, double q);

double NormMapinv(double y, double q);

double NormMap(double x, double q);

double DNormMapinv(double y, double q);

double NormB(double q);

int derivs_gsl(double x, double *y, void *param);

double Fq(double x, double q);

double Fqinv(double y, double q);

double fq(double x, double q);

double f1q(double x, double q);


int QCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) [nbPaths*nbNames] */
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double qM,                      /* (I) */
    double qZ);                      /* (I) */
