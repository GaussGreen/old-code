
typedef struct {double mu, sigma, weight;} DENSITY_COMPONENT;

typedef struct
{
    long n; /* number of componenets in the density (>=1) */
    DENSITY_COMPONENT *t;
} FACTOR_DENSITY;

typedef struct
{
    double beta;
    long   var_adjust;
    FACTOR_DENSITY fm;  /* density of factor M */
    FACTOR_DENSITY fz;  /* density of factor Z */
} X_DENSITY; /* density of X=bM + (1-b2)^1/2 Z */

double DiracDensity(double x, double mu);
double DiracCum(double x, double mu);
double F_DiracM_x_DiracZ(   double x,
                            double beta,
                            double mu_fm,
                            double mu_fz);

double F_DiracM_x_NormalZ(  double x,
                            double beta,
                            double mu_fm,
                            double mu_fz,
                            double sigma_fz);

double F_NormalM_x_DiracZ(  double x,
                            double beta,
                            double mu_fm,
                            double sigma_fm,
                            double mu_fz);

double F_NormalM_x_NormalZ( double x,
                            double beta,
                            double mu_fm,
                            double sigma_fm,
                            double mu_fz,
                            double sigma_fz);

double LossDensity(         double K,
                            double u,
                            X_DENSITY *fx);

double TrancheletLoss(      double K,
                            double u,
                            X_DENSITY *fx);

double FactorCum(double x, FACTOR_DENSITY *f);

double FactorDensity(double x, FACTOR_DENSITY *f);

double F_X(double x, X_DENSITY *fx);


double DpiDwM( double K0, double K1, double u, double beta, double muM, double sigmaM);
double DpiDwZ( double K0, double K1, double u, double beta, double muZ, double sigmaZ);
double DpiDbeta( double K0, double K1, double u, double beta);
