static double gamma_large (void *random,
                           const double a);
static double gamma_frac (void *random,
                          const double a);
double gsl_ran_gamma_int (void *random,
                          const unsigned int a);
int RandomGamma (void *random,
                 double *gammaSequence,
                 long nbPaths,
                 const double a,
                 const double b);
double gsl_ran_gamma (void *random,
                      const double a,
                      const double b);
int GammaDeviates (   double *alphaSequence,
                      long nbPaths,
                      const double a, 
                      const double b);
