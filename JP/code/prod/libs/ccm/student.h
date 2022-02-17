/* $Header$ */
int Chi2Deviates(
                  double *chi2,             /* (0) [nbPaths] */
                  long nbPaths,             /* (I) */
                  long freedomDegree,
                  long seed);      /* (I) */

int Chi2DeviatesSobol(
                  double *chi2,             /* (0) [nbPaths] */
                  long nbPaths,             /* (I) */
                  long freedomDegree);      /* (I) */


int StudentMultivariateDeviates(
           double *studentSequence,         /* (O) [nbPaths*nbNames] */
           long nbNames,                    /* (I) */
           long nbPaths,                    /* (I) */
           const double *beta,              /* (I) */
           long freedomDegree);             /* (I) */

int StudentCopulatedUniformDeviates(
        double *copulatedUniformDeviates,   /* (O) [nbPaths*nbNames] */
        long nbNames,                       /* (I) */
        long nbPaths,                       /* (I) */
        const double *beta,                 /* (I) */
        long freedomDegree);                /* (I) */

double StudentCopulaLinearCorrelation(
        double beta,                        /* (I) */
        long freedomDegree);                /* (I) */

/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/
int StudentCopulatedIndicator(
        int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
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
        long chi2_nbSample);            

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
    long chi2_nbSample);

int StudentMultivariateDeviates_tp(       
        double *studentSequence,        /* (O) [nbPaths*nbNames] */ 
        double *weight,                 /* (I) [nbPaths] */
        long nbNames,                   /* (I) */
        long nbPaths,                   /* (I) */
        const double *beta,             /* (I) */
        long freedomDegree,
        long seed);                     /* (I) */


