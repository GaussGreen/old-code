/* $Header$ */
int GaussianMultivariateDeviates(
    double *studentSequence,        /* (O) [nbPaths*nbNames] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta);            /* (I) [nbNames] */

int GaussianCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) [nbPaths*nbNames] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta);            /* (I) [nbNames] */

double GaussianCopulaLinearCorrelation(double beta); /* (I) */

/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/
int GaussianMultivariateDeviates_tp(
    double *studentSequence,        /* (O) [nbPaths*nbNames] */
    double *weight,                 /* (O) weight[nbPaths] */ 
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,
    long seed);            /* (I) [nbNames] */

int GaussianMultivariateDeviatesSampling_tp(double *gaussianSequence,
                                  double *weight,
                                  long nbNames,
                                  long nbPaths,
                                  const double *beta,
                                  long seed,
                                  const double *M,
                                  long nbSample);

int GaussianCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,              /* (I) [nbNames] */
    long seed,
    const double *M,
    long nbSample);

int GaussianCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    long seed,
    const double *M,
    long nbSample);

/***********************************************************************/
/***********************************************************************/
/*              STOCHASTIC RECOVERY IMPLEMENTATION                     */
/***********************************************************************/


int GaussianCopulatedIndicatorRecovery(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    double *X,                      /* (O) recovery stoch term */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    const double *alpha,
    long seed,
    const double *M,
    long nbSample);
