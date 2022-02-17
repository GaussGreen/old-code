/* $Header$ */
#define  MAXBUFF       250
#define  SUCCESS       0
#define  FAILURE       -1
#define  TRUE          1
#define  FALSE         0

double StudentCum(
                  double x,
                  long freedomDegree);

double StudentCumInverse(
                  double x,
                  long freedomDegree);

double StableDensity(double x,
                     double alpha, double beta);

double StableCum(
                  double x,
                  double alpha, double beta);

double StableCumInverse(
                  double x,
                  double alpha, double beta);

double NormalDensity(
                  double x);

double Chi2Density(
                  double x,
                  long freedomDegree);

double Chi2Cum(double x,
               long freedomDegree);


double NormalCum(
                  double x);

double NormalCumInverse(
                  double x);

void DensityBarChart(
                  double *bar,
                  const double *sample,
                  long sampleSize,
                  double start,
                  double stepSize,
                  long nbSteps);

int Moments(const double *probas,
            const double *values,
            long nbValues,
            double *moments,
            long nbMoments);

