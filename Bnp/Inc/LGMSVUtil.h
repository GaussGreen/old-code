#ifndef LGMSVUTILH
#define LGMSVUTILH

#include "math.h"
#include "srt_h_all.h"

/* Number max of function in multiple Integral */
#define LGMSV_NbMaxFunction 10
#define LGMSV_MINALPHA 0.005

/* -----------------------------------------------------------------------------------------------
        Definition of the structure LGMSVSolFunc
         = the general solution of a differential equation of type :
                dXt = (a*ft+b*gt+c*ht - dLambda*Xt)*dt

        ft and gt are either solution of the same type of diff eq
        or equal to the function 1
        a, b and c are three constants
   ------------------------------------------------------------------------------------------------
 */
typedef struct
{
    double dXt1;    /* X(t1) */
    double dLambda; /* Lambda */

    int    bIsft1; /* =1 if ft is the function 1, = 0 either */
    double a;      /* constant */
    void*  pft;    /* pointeur on function ft */

    int    bIsgt1; /* =1 if gt is the function 1, = 0 either */
    double b;      /* constant */
    void*  pgt;    /* pointeur on function gt */

    int    bIsht1; /* =1 if gt is the function 1, = 0 either */
    double c;      /* constant */
    void*  pht;    /* pointeur on function gt */

} LGMSVSolFunc;

/* -----------------------------------------------------------------------------------------------
        Definition of the structure LGMSVSolFunc
         = the general solution of a differential equation of type :
                dXt = (Sum(a[i]*ft[i]) - dLambda*Xt)*dt

        ft[i] are solution of the same type of diff eq
        or equal to the function 1
        a[i] are constants
   ------------------------------------------------------------------------------------------------
 */
typedef struct
{
    double dXt1;    /* X(t1) */
    double dLambda; /* Lambda */

    int    iNbFunction;                 /* Number of function ft */
    int    bIsft1[LGMSV_NbMaxFunction]; /* =1 if ft is the function 1, = 0 either */
    double a[LGMSV_NbMaxFunction];      /* constant */
    void*  pft[LGMSV_NbMaxFunction];    /* pointeur on function ft */
} LGMSVSolFunc1;

/* -----------------------------------------------------------------------------------------------
        LGMSVFuncValue
                Evaluation of
                        Xt  = Xt1*exp(-lambda*(t-t1))+integral between t1 and t of
   fs*exp(-lambda*(s-t1))

                where fs is of the same type as Xt
   ------------------------------------------------------------------------------------------------
 */
double LGMSVFuncValue1(
    /* Inputs */
    LGMSVSolFunc1* Xt, /*  Function to be valued */

    double  dDt,    /* = t-t1 */
    double* Lambda, /* Array of lambda */
    int     iLevelIntegration /* 0 for valuation */);

/* -----------------------------------------------------------------------------------------------
        LGMSVFuncValue
                Evaluation of
                        Xt  = Xt1*exp(-lambda*(t-t1))+integral between t1 and t of
   fs*exp(-lambda*(s-t1))

                where fs is of the same type as Xt
   ------------------------------------------------------------------------------------------------
 */
double LGMSVFuncValue(
    /* Inputs */
    LGMSVSolFunc Xt, /*  Function to be valued */

    double  dDt,    /* = t-t1 */
    double* Lambda, /* Array of lambda */
    int     iLevelIntegration /* 0 for valuation */);

void LGMSVFuncValue2(
    /* Inputs */
    LGMSVSolFunc* Xt, /*  Function to be valued */

    double  dDt,    /* = t-t1 */
    double* Lambda, /* Array of lambda */
    int     iLevelIntegration /* 0 for valuation */,
    double* res);

Err Get_LGMSV_TermStructure(
    char*    underlying,
    double** sigma_time,
    double** sigma,
    double** alpha,
    double** rho,
    double** lambdaeps,
    double*  Tstar,
    long*    sigma_n,
    double*  fixed_tau,
    /* extra 2F */
    int*     one2F,
    double** alpha2F_sv,
    double*  gamma2F_sv,
    double** rho2F_sv);

void ConvertAlphaRho_LGM_to_LGMSV(
    long    sigma_n,
    double* sigma_time,
    double  tstar,
    double  lambda,
    double  alpha,
    double  gamma,
    double  rho,
    double* new_alpha,
    double* new_rho);

void ConvertTS_LGMSV_to_LGM(
    long sigma_n, double* sigma_time, double* sigma, double lambda, double Tstar);

void ConvertTS_LGM_to_LGMSV(
    long    sigma_n,
    double* sigma_time,
    double* sigma,
    double  lambda,
    double  Tstar,
    /* extra for 2 Factor */
    int     one2F,
    double  alpha,
    double  gamma,
    double  rho,
    double* new_alpha,
    double* new_rho);

void Convert_Tstar_LGMSV(
    double  lambda,
    long    sigma_n,
    double* sigma_time,
    double* sigma,
    double* alphaeps,
    double* lameps,
    double* rhoeps,
    double  init_tstar,
    double  new_tstar);

typedef struct
{
    /* Market ID */
    long today;

    int one2F;

    int     num_sigma;
    long*   sigma_dates;
    double* sigma_times;
    double* sigma;

    int     num_tau;
    long*   tau_dates;
    double* tau_times;
    double* tau;

    int     num_smile;
    long*   smile_dates;
    double* smile_times;

    double alpha;
    double gamma;
    double rho;

    double* alphaSV;
    double* lamSV;
    double* rhoSV;
    double* rho2SV;

    double tstar;

} irm_sv, *IRM_SV;

void irm_sv_free_struct(IRM_SV irm_sv);

Err irm_sv_free_und_struct(SrtUndPtr pUndDesc);

Err irm_sv_get_struct_from_und(char* und, irm_sv** irmsv);

Err SrtInitIRMSVUnd(
    char* undName, /* und name */
    char* ycname,  /* mkt name */

    /* dimension */
    int one2F,

    /* volatility and tau */
    double** sigma_datas,
    int      num_sigma,
    int      sigma_col,
    double** tau_datas,
    int      num_tau,
    int      tau_col,

    /* LGM 2F parameters */
    double alpha,
    double gamma,
    double rho,

    /* Stoch vol parameters */
    double** smile_datas,
    int      num_smile,
    int      smile_col,

    /* T* of the model */
    double tstar);

Err irm_sv_get_term_struct(
    char*    und,
    int*     one2F,
    int*     num_sigma,
    double** sigma_times,
    double** sigma,
    int*     num_tau,
    double** tau_times,
    double** tau,
    double*  alpha,
    double*  gamma,
    double*  rho,
    int*     num_smile,
    double** smile_times,
    double** alphaSV,
    double** lamSV,
    double** rhoSV,
    double** rho2SV,
    double*  tstar);

Err irm_sv_get_term_struct_date(
    char*    und,
    int*     one2F,
    int*     num_sigma,
    long**   sigma_dates,
    double** sigma,
    int*     num_tau,
    long**   tau_dates,
    double** tau,
    double*  alpha,
    double*  gamma,
    double*  rho,
    int*     num_smile,
    long**   smile_dates,
    double** alphaSV,
    double** lamSV,
    double** rhoSV,
    double** rho2SV,
    double*  tstar);

/*	LGMSV model */
typedef struct
{
    double dLambdaX;
    double dTau;

    double dLambdaX2;
    double dTau2;

    int     iNbPWTime;
    double* dPWTime;
    double* dSigma;
    double* dAlpha;
    double* dLambdaEps;
    double* dLvlEps;
    double* dRho;
    double* dRho2;
    double  dTStar;
    double  dInitTStar;

    double dInitLambdaEps;

    int     iOne2F;
    double  dInitLGMAlpha;
    double  dInitLGMRho;
    double* dLGMAlpha;
    double  dLGMGamma;
    double* dLGMRho;

    long lToday;
    long lTStarDate;

    double dFlatOneFactorRho;

} LGMSV_model, *LGMSV_MODEL;

void init_NULL_LGMSV_model(LGMSV_MODEL model);

Err init_LGMSV_model(
    LGMSV_MODEL model,
    long        lToday,
    int         iOne2F,
    int         iNbPWTime,
    double      dLambdaX,
    double*     dPWTime,
    double*     dSigma,
    double*     dAlpha,
    double*     dLambdaEps,
    double*     dLvlEps,
    double*     dRho,
    double      dTStar,
    double      dLGMAlpha,
    double      dLGMGamma,
    double      dLGMRho,
    double*     dRho2);

void free_LGMSV_model(LGMSV_MODEL model);

Err fill_LGMSV_model_from_irm_sv(irm_sv* irmsv, LGMSV_MODEL model);

Err Get_LGMSV_model(char* undname, LGMSV_MODEL model);

Err Check_Model_Input_Corr_Matrix(
    int iOne2F, int iNbPWTime, double* dLGMRho, double* dRho, double* dRho2);

Err Check_Model_Corr_Matrix(LGMSV_MODEL model);

void Convert_Tstar_model(LGMSV_MODEL model, double new_tstar);

#endif