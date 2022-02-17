int ClaytonCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha,
    long seed);                  /* (I) */

int ClaytonCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,          /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha,
    long seed);

int ClaytonCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O)  */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double alpha);             /* (I) */
