
#ifndef LGMSVCLOSEDFORM_H
#define LGMSVCLOSEDFORM_H

#include "LGMSVUtil.h"
#include "cpdcalib.h"

LGMSVSolFunc1 ***LGMSVSolFunc3tensor(long nrl, long nrh, long ncl, long nch,
                                     long ndl, long ndh);
void free_LGMSVSolFunc3tensor(LGMSVSolFunc1 ***t, long nrl, long nrh, long ncl,
                              long nch, long ndl, long ndh);

void LGMSVFillUsefullMoment(/* Input */
                            int i, int j, int k,

                            /* Output */
                            long ***IsUsefullMoment);

void LGMSVMomentCalculation(
    /* Inputs */
    int iNbCumulant, /* Number of cumulant for the approximation of the (f
                        ,Phi) density */

    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX,

    double dAlpha, double dLambdaEps,

    double dRho,

    double dTExercice, /* Exercice of the swaption in years from today */
    double dTStar,     /* Tstar */

    int iNbSigTime, /* Term Structure of g(t) */
    double *SigTime, double *Sig,

    double *lambdaArray,

    /* Outputs */
    long ***IsUsefullMoment, LGMSVSolFunc1 ***MomentFunc, double ***MomentValue,
    double *pPhitMean);

void LGMSVMomentInit2(
    /* Inputs */
    double dLambdaX, double dAlpha, double dLambdaEps, double dRho,

    LGMSVSolFunc *FuncPhi, LGMSVSolFunc *FuncPhi2, LGMSVSolFunc *FuncV2,
    LGMSVSolFunc *FuncPhiV);

void LGMSVMomentCalculation2(

    double Sig2, double dt,

    double InitPhi, double InitPhi2, double InitV2, double InitPhiV,

    LGMSVSolFunc *FuncPhi, LGMSVSolFunc *FuncPhi2, LGMSVSolFunc *FuncV2,
    LGMSVSolFunc *FuncPhiV, double *lambda,

    double *ResPhi, double *ResPhi2, double *ResV2, double *ResPhiV);

void rlft3(double ***data, double **speq, unsigned long nn1, unsigned long nn2,
           unsigned long nn3, int isign);

void LGMSVUpdateADClosedForm(double dt, double coefInt, double coefARe,
                             double coefBRe, double coefBIm, double coefCRe,
                             double coefCIm,

                             double *DRe, double *DIm, double *ARe,
                             double *AIm);

void LGMSVFillPayoff(/* Input */
                     int iNbPhi, int iNbft, double dLambdaX, double dPhitMean,
                     int iIndexPhiMean, double dPhiStep, int iIndexft0,
                     double dftStep, double dExTime, /* Exercice of the swaption
                                                        in years from today */
                     double dTStar,                  /* Tstar */
                     int iNbCoupon, /* Description of the cashflows */
                     double *CouponTime, double *Coupon,

                     /* Outputs */
                     double **Payoff);

void LGMSVSaveDensity(/* Inputs */
                      int iNbPhi, int iNbft,

                      int iIndexPhiMean, int iIndexft0,

                      double dPhiStep, double dftStep, double dPhitMean,
                      double dPhitStd,

                      double ***Density);

void LGMSVSaveMoment(/* Inputs */
                     int iNbCumulant, double ***MomentValue);

void LGMSVSavePayOff(/* Inputs */
                     int iNbPhi, int iNbft,

                     int iIndexPhiMean, int iIndexft0,

                     double dPhiStep, double dftStep, double dPhitMean,

                     double **PayOff);

void LGMSVSaveTFDensity(/* Inputs */
                        int iNbPhi, int iNbft,

                        double ***Density);

void LGMSVCalculateDensityClosedForm(
    /* Input */
    int iNbCumulant, int iNbPhi, int iNbft,

    double dPhiFreqStep, double dftFreqStep,

    /* Parameters */
    double sig0, double alpha, double rho, double T,

    /* Outputs */
    double ***Density, double **DensitySpeq);

/* -------------------------------------------------------------------------------------------------------------
        LGMSVCalib

  --------------------------------------------------------------------------------------------------------------
*/
Err LGMSVCalib(
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long ExDate, long end_date, double dStrike,

    double dTau, double dg, double dAlpha, double dRho, double dLambdaEps,
    long iNbCumulant, long iNbPhi, long iNbft,

    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaXGridLeft, double iNbSigmaXGridRight,
    /* Output */
    double *pSwaptionPrice);

Err LGMSVOption(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, double dStrike,

    double iMaxTime, double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, int iSaveFile, long iNbPhi, long iNbft,

    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaXGridLeft, double iNbSigmaXGridRight,

    /* Output */
    double *pSwaptionPrice);

Err LGMSVOptionNew(
    char *und_name,      /*	Name of the underlying */
    char *yc_name,       /*	Name of the yield curve */
    char *ref_rate_name, /*	Name of the reference rate */
    char *swaption_freq, /*	Frequency and basis of underlying swaptions */
    char *swaption_basis,

    long lExDate, long lStartDate, long lEndDate, double dStrike,
    int pay_rec, /*	pay:1 rec:-1 */

    double iMaxTime, double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, int iSaveFile, long iNbPhi, long iNbft,

    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaXGridLeft, double iNbSigmaXGridRight,

    /* Output */
    double *pSwaptionPrice);

void LGMSVFindFreq(int iNbPhi, int iNbft,

                   double AlphaEq, double Rho,

                   int endi, double expt2init, double *Coef1ReT,
                   double *Coef2ReT, double *Coef2ImT, double *Coef3ReT,
                   double *Coef3ImT, double *CoefAT, double *CoefexpT,
                   int *nbpT,

                   double dPhitMean, double dPhitStd, double dftStd,

                   double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
                   double iNbSigmaftLeft, double iNbSigmaftRight,
                   double iLimitPhi, double iLimitft, int iPriorityFreqPhi,
                   int iPriorityFreqFt,

                   /* Outputs */
                   double *dPhiFreqStep, double *dPhiStep, int *iIndexPhiMean,
                   double *dftFreqStep, double *dftStep, int *iIndexft0);

void LGMSVFindFreqNew(int iNbPhi, int iNbft,

                      double AlphaEq, double Rho,

                      int endi, double *dt, double Coef1Re, double Coef2Re,
                      double *Coef2ImT, double *Coef3ReT, double *Coef3ImT,

                      double dPhitMean, double dPhitStd, double dftStd,

                      double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
                      double iNbSigmaftLeft, double iNbSigmaftRight,
                      double iLimitPhi, double iLimitft, int iPriorityFreqPhi,
                      int iPriorityFreqFt,

                      /* Outputs */
                      double *dPhiFreqStep, double *dPhiStep,
                      int *iIndexPhiMean, double *dftFreqStep, double *dftStep,
                      int *iIndexft0);

void LGMSVCalculateDensityTPPoint(double dPhifreq, double dftFreq,

                                  int endi, double expt2init, double *Coef1ReT,
                                  double *Coef2ReT, double *Coef2ImT,
                                  double *Coef3ReT, double *Coef3ImT,
                                  double *CoefAT, double *CoefexpT, int *nbpT,

                                  double PhiMean,

                                  /* Outputs */
                                  double *LogTFRe, double *LogTFIm);

void LGMSVCalculateDensityPrecalc(int iNbPhi, int iNbft, double dPhiFreqStep,
                                  double dftFreqStep,

                                  int endi, double expt2init, double *Coef1ReT,
                                  double *Coef2ReT, double *Coef2ImT,
                                  double *Coef3ReT, double *Coef3ImT,
                                  double *CoefAT, double *CoefexpT, int *nbpT,

                                  double PhiMean, int SaveFile,

                                  /* Outputs */
                                  double ***Density, double **DensitySpeq);

void LGMSVCalculateDensityPrecalcNew(int iNbPhi, int iNbft, double dPhiFreqStep,
                                     double dftFreqStep,

                                     int endi, double *dt, double Coef1Re,
                                     double Coef2Re, double *Coef2ImT,
                                     double *Coef3ReT, double *Coef3ImT,

                                     double PhiMean, int SaveFile,

                                     /* Outputs */
                                     double ***Density, double **DensitySpeq);

void LGMSVOptionPricePrecalc(
    int iNbPhi, int iNbft,

    double dLambdaX, double dTStar, int endi, double expt2init,
    double *Coef1ReT, double *Coef2ReT, double *Coef2ImT, double *Coef3ReT,
    double *Coef3ImT, double *CoefAT, double *CoefexpT, int *nbpT,

    double dPhitMean, int SaveFile, double ***Density, double **DensitySpeq,
    double **PayOff,

    double dExTime, int iNbCoupon, double *CouponTime, double *Coupon,

    double dPhiFreqStep, double dPhiStep, int iIndexPhiMean, double dftFreqStep,
    double dftStep, int iIndexft0,

    double *Price);

void LGMSVClosedFormNew(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX, int iNbSigTime,             /* Term Structure of g(t) */
    double *SigTime, double *Sig, double dTStar, /* Tstar in years from today */
    double dAlpha, double dLambdaEps, double dRho,

    /* Product description */
    long lExDate,   /* Exercice date of the swaption  */
    double dExTime, /* Exercice of the swaption in years from today */
    int iNbCoupon,  /* Description of the cashflows */
    double *CouponTime, long *CouponDate, double *Coupon,
    char *cYieldCurve, /* Yield Curve */

    /* Parameter of grids */
    int iNbPhi, /* Number of Phi : Should be a power of two */
    int iNbft,  /* Number of ft : Should be a power of two */
    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, double iMaxTime, int SaveFile,

    /* Outputs */
    double *Price);

void LGMSVClosedFormOnVol(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX, int iNbSigTime,             /* Term Structure of g(t) */
    double *SigTime, double *Sig, double dTStar, /* Tstar in years from today */
    double dAlpha, double dLambdaEps, double dRho,

    /* Product description */
    long lExDate,   /* Exercice date of the swaption  */
    double dExTime, /* Exercice of the swaption in years from today */
    int iNbCoupon,  /* Description of the cashflows */
    double *CouponTime, long *CouponDate, double *Coupon,
    char *cYieldCurve, /* Yield Curve */

    /* Parameter of grids */
    int iNbPsi,   /* Number of Psi : Should be a power of two */
    int iNbDrift, /* Number of ft : Should be a power of two */
    double iNbSigmaPsiGridLeft, double iNbSigmaPsiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPsi, double iRatioFt, int iPriorityFreqPsi,
    int iPriorityFreqFt, double iMaxTime, int SaveFile,

    /* Outputs */
    double *Price);

Err LGMSVOptionGRFNPrice(
    int iNbPhi, int iNbft,

    double dLambdaX, double dTStar, int endi,

    double *dt, double Coef1Re, double Coef2Re, double *Coef2ImT,
    double *Coef3ReT, double *Coef3ImT,

    double dPhitMean, int SaveFile,

    double dExTime, long dExDate,

    double dPhiFreqStep, double dPhiStep, int iIndexPhiMean, double dftFreqStep,
    double dftStep, int iIndexft0,

    /*	Product data */
    char *yc, void **func_parm_tab,

    /* GRFN function */
    Err (*payoff_func)(/* Time Information */
                       double dDate, double dTime, void *func_parm,

                       /* Market data	*/
                       void *cYieldCurve,

                       /*	Model data Information	*/
                       double dLambdaX,

                       /* Martingale Mesure QTstar information */
                       double dTstar, /* In time */

                       /* Grid data Information	*/
                       int iNbPhi, int iNbft,

                       int iIndexPhiMean, int iIndexft0, double dPhiStep,
                       double dftStep, double dPhitMean,

                       /* Tensor of results to be updated		*/
                       /* 4 dimensions : Phit  ,Xt  ,Epst  ,Product	*/
                       int iNbProduct, double ***PayoffTensor),

    int iNbProduct, double *Price);

void LGMSVCalib2(
    /* Parameter of diffusion */
    /* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
    double dLambdaX, double dTStar, double dAlpha, double dLambdaEps,
    double dRho,

    /* Product description */
    char *cYieldCurve, double *dExTime, int nbEx, int *iNbCoupon,
    double **CouponTime, double **Coupon, double *TargetPrice,
    double *TargetVol,

    /*	Additional informations to estimate Vega */
    double *Fwd, double *Strike, double *Level,

    /* Parameter of grids */
    int iNbPhi, /* Number of Phi : Should be a power of two */
    int iNbft,  /* Number of ft : Should be a power of two */
    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, double iMaxTime,

    /* Newton parameters */
    double Precision, int NbIterMax,

    /* Outputs */
    double *sigma);

Err cpd_calibSV_diagonal_2(
    /*	Market */
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *instr_freq, /*	Frequency and basis of instruments */
    char *instr_basis,
    /*	If ex_date is NULL  ,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates  ,
                                                      all supposed to be on or
                         after today */
    long *ex_date,    /*	Supposed to be sorted */
    char **end_tenor, /*	Tenors of the underlying instruments
                                                              or "DIAG" */
    long end_date,    /*	End date for diagonal */
    double *strike,   /*	Strikes
                                              0: ATM */
    /*	Model */
    double dLambdaX, double dAlpha, double dLambdaEps, double dRho,
    double dTStar,

    /* Parameter of grids */
    int iNbPhi, /* Number of Phi : Should be a power of two */
    int iNbft,  /* Number of ft : Should be a power of two */
    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, double iMaxTime,

    /* Newton parameters */
    double Precision, int NbIterMax,
    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data);

Err cpd_calibSV_diagonal_new(
    /*	Market */
    char *yc_name,        /*	Name of the yield curve */
    char *vol_curve_name, /*	Name of the market vol curve */
    char *ref_rate_name,  /*	Name of the reference rate */
    Err (*get_cash_vol)(  /*	Function to get cash vol from the market */
                        char *vol_curve_name, double start_date,
                        double end_date, double cash_strike, int zero,
                        char *ref_rate_name, double *vol, double *power),
    char *instr_freq, /*	Frequency and basis of instruments */
    char *instr_basis,
    /*	If ex_date is NULL  ,
    exercise dates will be generated 2bd before start */
    /*	Structure */
    int num_ex_dates, /*	Exercise dates  ,
                                                      all supposed to be on or
                         after today */
    long *ex_date,    /*	Supposed to be sorted */
    char **end_tenor, /*	Tenors of the underlying instruments
                                                              or "DIAG" */
    long end_date,    /*	End date for diagonal */
    double *strike,   /*	Strikes
                                              0: ATM */
    /*	Model */
    double dLambdaX, double dAlpha, double dLambdaEps, double dRho,
    double dTStar,

    /* Parameter of grids */
    int iNbPhi, /* Number of Phi : Should be a power of two */
    int iNbft,  /* Number of ft : Should be a power of two */
    double iNbSigmaPhiGridLeft, double iNbSigmaPhiGridRight,
    double iNbSigmaftLeft, double iNbSigmaftRight,

    double iRatioPhi, double iRatioFt, int iPriorityFreqPhi,
    int iPriorityFreqFt, double iMaxTime,

    /* Newton parameters */
    double Precision, int NbIterMax,
    /*	Output */
    int *num_sig, /*	Answer */
    double **sig_time, double **sig,
    /*	Parameters */
    CPD_DIAG_CALIB_PARAM param,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data);

void LGMSV_prod_comp(double re1, double im1, double re2, double im2,
                     double *re3, double *im3);

void LGMSV_div_comp(double re1, double im1, double re2, double im2, double *re3,
                    double *im3);

void LGMSV_sqr_comp(double re1, double im1, double *re3, double *im3);

static void LGMSV_solve_poly2_comp(double Are, double Aim, double Bre,
                                   double Bim, double Cre, double Cim,
                                   double *x1re, double *x1im, double *x2re,
                                   double *x2im, double *DeltaRe,
                                   double *DeltaIm);

#endif