#include "skew.h"

double DLoss(double K,
             double u,
             double beta,
             double qM,
             double qZ,
             double alpha,
             long nbPoints,
             double eps);

double  Loss(double K,
             double u,
             double beta,
             double qM,
             double qZ,
             double alpha,
             long nbPoints,
             double eps);

double Mode(double u,
            double beta,
            double qM,
            double qZ,
            double alpha,
            long nbPoints,
            double eps);

double FMode(double u,
             double beta,
             double qM,
             double qZ,
             double alpha,
             long nbPoints,
             double eps);

double Fij(double u,
           double v,
           double betai,
           double betaj,
           double qM,
           double qZ,
           long nbPoints);

double F2(  double u1,
            double u2,
            double beta1,
            double beta2,
            double qM,
            double qZ1,
            double qZ2,
            long NbPoints,
            double eps);

double Moment(      double u,
                    double qM,
                    double qZ,
                    double alpha,
                    double beta,
                    long n,
                    long nbPoint,
                    double eps);

double IntegrandLn( double u,
                    double beta,
                    double qM,
                    double qZ,
                    double K,
                    long n
                    );

double Kn(          double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long n,
                    long nbPoint,
                    double eps
                    );

double Sigma(       double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    );

double Skew(        double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    );

double Kurtosis(    double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    long nbPoint,
                    double eps
                    );

double SeniorProba(double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    double K,
                    long nbPoint,
                    double eps
                    );

double Quantile(    double u,
                    double beta,
                    double qM,
                    double qZ,
                    double alpha,
                    double proba,
                    long nbPoint,
                    double eps
                    );

double DModeDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DModeDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DModeDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DModeDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);

double DFModeDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DFModeDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DFModeDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DFModeDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);

double DMomentDqM(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DMomentDqZ(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DMomentDbeta(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DMomentDalpha(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);

double DKnDqM(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DKnDqZ(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DKnDbeta(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);
double DKnDalpha(double u, double beta, double qM, double qZ, double alpha, long n, long nbPoints, double eps, double h);

double DSigmaDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSigmaDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSigmaDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSigmaDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);

double DSkewDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSkewDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSkewDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DSkewDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);

double DKurtosisDqM(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DKurtosisDqZ(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DKurtosisDbeta(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);
double DKurtosisDalpha(double u, double beta, double qM, double qZ, double alpha, long nbPoints, double eps, double h);

double DSeniorProbaDbeta(   double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DSeniorProbaDqM(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DSeniorProbaDqZ(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DSeniorProbaDalpha(  double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DQuantileDbeta(      double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DQuantileDqM(        double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DQuantileDqZ(        double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DQuantileDalpha(     double u,
                            double beta,
                            double qM,
                            double qZ,
                            double alpha,
                            double K,
                            long nbPoint,
                            double eps,
                            double h
                        );

double DIntegrandLnDqM(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h);
double DIntegrandLnDqZ(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h);
double DIntegrandLnDbeta(double u, double beta, double qM, double qZ, double K, long n, long nbPoints, double eps, double h);
