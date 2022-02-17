#ifndef _PAYOFF_TP_
#define _PAYOFF_TP_
/* $Header$ */
#include "creditPortfolio.h"

typedef enum {
    C_GAUSS        = 0L,
    C_STUDENT      = 1L,
    C_DEPENDENCE   = 2L,
    C_INDEPENDENCE = 3L,
    C_GUMBEL       = 4L,
    C_CLAYTON      = 5L,
    C_COMPOSITE    = 6L,
    C_SKEW         = 7L,
    C_DELTA        = 8L,
    C_COPULA_FREE  = 9L
} CopulaType;

typedef struct
{
    CopulaType type;
    void *param;
} COPULA;

typedef struct
{
    double *beta;
    long seed;
    double qM;
    double qZ;
} Q_PARAM;

typedef struct
{
    double *beta;
    long seed;
    double xa;
    double xb;
    double xc;
    double delta;
} DELTA_PARAM;

typedef struct
{
    double *beta;
    long seed;
    long freedomDegree;
    double *M_sample;
    double *chi2_sample;
    long M_nbSample;
    long chi2_nbSample;
} STUDENT_PARAM;

typedef struct
{
    double *beta;
    long seed;
    double *M_sample;
    long M_nbSample;
} GAUSSIAN_PARAM;

typedef struct
{
    double theta;
    long seed;
} ARCHIMEDEAN_PARAM;

typedef struct
{
  COPULA *c_1;
  COPULA *c_2;
  double *alpha;
} PRODUCT_PARAM;

typedef struct
{
    long nbNames;
    long nbPaths;
    int *indicator;
    double *weight;
} INDICATOR_SIM;

INDICATOR_SIM IndicatorSimCreate(   long nbPaths,
                                    long nbNames,
                                    const double *survivalProba,
                                    const COPULA *copula);

INDICATOR_SIM IndicatorSimCreate_mtp(   long nbPaths,
                                        long nbNames,
                                        long nbTimes,
                                        const double *survivalProba,
                                        const COPULA *copula);


void IndicatorSimFree(INDICATOR_SIM *is);

COPULA CopulaCreate(    CopulaType type,
                        long nbNames,
                        const double *beta,
                        long freedomDegree,
                        double theta,
                        long seed,
                        double qM,
                        double qZ,
                        double a,
                        double b,
                        double c,
                        double delta,
                        const double *M,
                        long M_nbSample,
                        const double  *chi2,
                        long chi2_nbSample,
                        int *status);

COPULA CopulaProductCreate(const COPULA *c_1,
                           const COPULA *c_2,
                           const double *alpha,
                           long nbNames);

COPULA Copula3ProductCreate(const COPULA *c_1,
                            const COPULA *c_2,
                            const COPULA *c_3,
                            const double *alpha1,
                            const double *alpha2,
                            long nbNames);

int CopulaCopy(const COPULA *c_1,
               COPULA *c_2,
               long nbNames);

void CopulaFree(COPULA *copula);

INDICATOR_SIM IndicatorProductSimCreate(long nbPaths,
                                        long nbNames,
                                        const double *survivalProba,
                                        const COPULA *c1,
                                        const COPULA *c2,
                                        const double *alpha);

INDICATOR_SIM IndicatorProductSimCreate_mtp(long nbPaths,
                                            long nbNames,
                                            long nbTimes,
                                            const double *survivalProba,
                                            const COPULA *c1,
                                            const COPULA *c2,
                                            const double *alpha);

double ExpectedPayoff_tp(       
            CREDIT_PORTFOLIO *port,     /* (I) */
            INDICATOR_SIM    *is);      /* (I) */

int NbDefaultNameDistribution(INDICATOR_SIM i_s, double *distribution);

int NbConditionalDefaultNameDistribution_mtp(INDICATOR_SIM i_s,
                                  double *distribution,
                                  long idx_t0,
                                  long idx_t1,
                                  long idx_T0,
                                  long idx_T1);

int PayoffDistribution_tp(  CREDIT_PORTFOLIO *port,
                            INDICATOR_SIM *is,
                            const double *sampleLoss,
                            long nbSampleLoss,
                            double *loss);

int Convolution(const double *v1, const double *v2, double *v, long n);

int ComplexMatrix(double *realPart,
                  double *imagPart,
                  long nbValues,
                  double *complexMatrix);

int MatrixProduct(double *M1,
                  double *M2,
                  long n,
                  double *M);

void Integer(long n, double *out);
void MatInteger(long n, double *out);

double TrancheletLoss3C(double u,
                        double beta,
                        double a,
                        double c,
                        double K
                       );

double TrancheLoss3C(   double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1);

double DTrancheletLoss3CDbeta(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       );

double DTrancheletLoss3CDa(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       );

double DTrancheletLoss3CDc(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double h
                       );

double DTrancheLoss3CDbeta(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       );

double DTrancheLoss3CDa(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       );

double DTrancheLoss3CDc(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double h
                       );

double Aeq( double u,
            double beta_star,
            double beta,
            double c,
            double K);

double ATeq( double u,
            double beta_star,
            double beta,
            double c,
            double K0,
            double K1);

double LossDensity3C(double u, double beta, double a, double c, double K);

double Moment3C( double u, double beta, double a, double c, long n, long nbPoint, double eps);

double BetaEq(double u, double beta_star, double a, double c, long nbPoint, double eps);

double DTrancheletLoss3CDcAAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double Kmatch,
                        double h
                       );

double DTrancheLoss3CDcAAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K0,
                        double K1,
                        double Kmatch,
                        double h
                       );

double DTrancheletLoss3CDcATAdjust(double u,
                        double beta,
                        double a,
                        double c,
                        double K,
                        double Kmatch0,
                        double Kmatch1,
                        double h
                       );

#endif

