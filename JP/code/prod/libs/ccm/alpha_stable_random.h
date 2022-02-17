#ifndef M_PI
#define M_PI	   3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2	   1.57079632679489661923132169164      /* pi/2 */
#endif

int RandomAlphaStable (     void *randomSequence,
                            double *levySequence,
                            long nbPaths,
                            const double c, 
                            const double alpha);

int AlphaStableDeviates (   double *alphaSequence,
                            long nbPaths,
                            const double c, 
                            const double alpha);

int RandomAlphaStableSkew ( void *randomSequence,
                            double *levySequence,
                            long nbPaths,
                            const double c, 
                            const double alpha,
                            const double beta);

int AlphaStableSkewDeviates(double *alphaSequence,
                            long nbPaths,
                            const double c, 
                            const double alpha,
                            const double beta);
