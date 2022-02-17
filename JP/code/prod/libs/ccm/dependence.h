/* $Header$ */
int DependenceCopulatedUniformDeviates(    
          double *copulatedUniformDeviates, /* (O) [nbPaths*nbName] */
          long nbNames,                     /* (I) */
          long nbPaths);                    /* (I) */

/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/
int DependenceCopulatedIndicator(
    int* copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) [nbPaths] */
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,
    long seed);                  /* (I) */


/**
 * Calculate survival indicator function for this copula.
 */
int DependenceCopulatedIndicator_mtp(
    int* copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) */
    const double *survivalProba,          /* (O) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,
    long seed);           
