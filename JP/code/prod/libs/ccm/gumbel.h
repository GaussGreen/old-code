int GumbelCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta,
    long seed);                  /* (I) */

int GumbelCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta,
    long seed);       

int GumbelTest( long nbNames,                   /* (I) */
                long nbPaths,                   /* (I) */
                double theta,
                double *prob);


int GumbelCopulatedUniformDeviates(
    double *copulatedUniformDeviates,/* (O) 0 if name defaulted before T */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    double theta);                   /* (I) */
