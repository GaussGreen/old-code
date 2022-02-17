
typedef struct
{
    long    nbNames;
    double *notional;               /* [nbNames] */
    double *recovery;               /* [nbNames] */
    double *sigma_recovery;
    double *nameMaturity;           /* [nbNames] */
    double  strike1;
    double  strike2;
} CREDIT_PORTFOLIO_R;

typedef enum {
    CR_GAUSS        = 0L
} CopulaRType;

typedef struct
{
    CopulaRType type;
    void *param;
} COPULA_R;

typedef struct
{
    long nbNames;
    long nbPaths;
    int *indicator;
    double *X;
    double *weight;
} INDICATOR_RECOVERY_SIM;

typedef struct
{
    double *beta;
    double *alpha;
    long seed;
} GAUSSIAN_PARAM;

CREDIT_PORTFOLIO_R CreditPortfolioRCreate(
        double *notional,           /* (I) [nbNames] */
        double *recovery,           /* (I) [nbNames] */
        double *sigma_recovery,
        double *nameMaturity,       /* (I) [nbNames] */
        long nbNames,               
        double strike1,             
        double strike2);

void IndicatorRecoverySimFree(INDICATOR_RECOVERY_SIM *is);

INDICATOR_RECOVERY_SIM IndicatorRecoverySimCreate(   long nbPaths,
                                    long nbNames,
                                    const double *survivalProba,
                                    const COPULA_R *copula,
                                    const double *M,
                                    long M_nbSample);

void CreditPortfolioRFree(CREDIT_PORTFOLIO_R *cp);

COPULA_R CopulaRCreate( long type,
                        long nbNames,
                        const double *beta,
                        const double *alpha,
                        int *status);

void CopulaRFree(COPULA_R *copula);

double ExpectedPayoff_R(   CREDIT_PORTFOLIO_R *port,
                           INDICATOR_RECOVERY_SIM *is);

int NbDefaultNameDistribution_R(INDICATOR_RECOVERY_SIM i_s, double *distribution);

int PayoffDistribution_R(   CREDIT_PORTFOLIO_R *port,
                            INDICATOR_RECOVERY_SIM *is,
                            const double *sampleLoss,
                            long nbSampleLoss,
                            double *loss);

