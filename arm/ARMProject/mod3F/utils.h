#include "DKMaille.h"
#include "DKMaille2D.h"

#define Pi 3.1415926535897932385
#define MAX(x,y)    (x >= y ? x : y)


static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1)>(maxarg2) ?\
(maxarg1) : (maxarg2)) 
#define FMIN(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1)<(maxarg2) ?\
(maxarg1) : (maxarg2)) 

double CumNormDist(double Z);

double NormDist(double x, 
				double x0, 
				double Sigma);

double Integration(double (*func)(double x, double *Y), 
				   double *Y, 
				   double LB, 
				   double HB, 
				   int NbSteps);

void Create_Strip_Analytics(DKMaille<double> &dStrip,
                            DKMaille<double> &dStripIntervalles,
                            DKMaille<double> NoticeDates,
                            int N1,
                            int N2,
                            int N3,
                            double dOptimalDate,
                            double LastDate,
                            int* NbPasTotal,
                            int* NbPasBeforeLastNotice);

double ZC_interpole(double T1, 
					double T2, 
					double ZC1, 
					double ZC2, 
					double T);


double FromRateToX(double r, 
                   double q, 
                   double K);

double FromXToRate(double X, 
                   double q,
                   double K);

double FromXToRate_deriv(   double X, 
                            double q, 
                            double K);

double BondPrice(   double A, 
                    double MeanReversion, 
                    double t, 
                    double T, 
                    double r);
