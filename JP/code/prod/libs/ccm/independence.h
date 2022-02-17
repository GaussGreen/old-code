/* $Header$ */
int IndependenceCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) [nbPaths*nbNames] */ 
    long nbNames,                   /* (I) */
    long nbPaths);                  /* (I) */

/***********************************************************************/
/***********************************************************************/
/*              SINGLE TIME POINT IMPLEMENTATION                       */
/***********************************************************************/
int IndependenceCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
                                    /*     [nbPaths*nbNames]            */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,
    long seed);                  /* (I) */


int IndependenceCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,
    long nbPaths,
    long seed);
