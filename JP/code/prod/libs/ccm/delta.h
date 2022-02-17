double FdeltaC(double x, double a, double b, double delta);

double FdeltaCinv(double y, double a, double b, double delta);

double Fdelta(double x, double a, double b, double c, double delta);

double Fdeltainv(   double y,
                    double a,
                    double b,
                    double c,
                    double delta);

double DeltaMap(double x, double a, double b, double c, double delta);

int DeltaMultivariateDeviates_tp(
    double *deltaSequence,              /* (O) [nbPaths*nbNames] */
    double *weight,                 /* (O) weight[nbPaths] */ 
    long nbNames,                   /* (I) proba[nbName] */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double a,                       /* (I) */
    double b,                       /* (I) */
    double c,                       /* (I) */
    double delta,                   /* (I) */
    long seed);


int DeltaCopulatedIndicator(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) [nbNames] */
    double a,                      /* (I) */
    double b,                      /* (I) */
    double c,
    double delta,
    long seed);

int DeltaCopulatedIndicator_mtp(
    int *copulatedSurvivalIndicator,/* (O) 0 if name defaulted before T */
    double *weight,                 /* (O) weight[nbPaths] */ 
    const double *survivalProba,    /* (I) proba[nbName] */
    long nbTimes,                   /* (I) */
    long nbNames,                   /* (I) */
    long nbPaths,                   /* (I) */
    const double *beta,             /* (I) */
    double a,
    double b,
    double c,
    double delta,
    long seed);
