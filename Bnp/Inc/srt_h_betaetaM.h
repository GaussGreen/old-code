#define SMAX 100
#define THETAMAX 100
#define LAMBDAMAX 5
#define XMAX 7
#define TMAX 10

static double integrand_1(double y_bar);
static double integrand_2(double y_bar);
static double integrand_3(double y_bar);
static double integrand_4(double y_bar);

Err make_M_matrix(TermStruct *l);

/*Err calc_bigM(double **M  , double beta2  , double eta2  , double lambda2  ,
   double deltau2  , double x_k2  , int i  , int j);
*/
double print_M(double beta2, double eta2, double lambda2, double deltau2,
               double x_k2);

Err srt_f_betaeta(SrtUndPtr und, SrtMdlDim mdl_dim, SwapDP *sdp, double strike,
                  double bond_strike, SrtReceiverType rec_pay, StructType type,
                  String ref_rate_code, double *answer);