#ifndef __OTCUTILS_H
#define __OTCUTILS_H

#define MAXCPN 512
#include "CPDProdstruct.h"
#include "CPDVol.h"

/* ==================================
Structure for OTC
=================================== */
typedef struct {
  int smileModel; // 0 = SSL; 1 = Log Mix; 2 = Beta Mix
  int Options_number;

  /* Fwd Smile method */
  int FwdVolMethod;
  char *FSMsigma;
  char *FSMalpha;
  char *FSMbeta;
  char *FSMrho;

  /* Smile Calibration */
  double smile_sigma;
  double smile_sigmaBeta;
  double smile_alpha;
  double smile_beta;
  double smile_rho;

  /* SSL */
  double SSLstd;
  int SSLnIter;
  double SSLshift;
  double SSLvolUp;
  double SSLvolDown;

  /* Beta Mix */
  double BMnStd;
  double BMsigma;
  double BMalpha;
  double BMbeta;
  double BMlambda;
  double BMpi;
  double BMfwd1;
  double BMfwd2;
  double BMsig1;
  double BMsig2;

  /* CUMULATIVE */
  long CUMULnPoints;
  int CUMULnStd;
  double CUMULPrecision;
  int CUMULlinear;
  double **CUMUL;

  /* COPULA */
  int COPULAdo_pecs;
  int COPULAnSimul;
  double **COPULAgauss;
  double **COPULAmatrix;

  /* OTC */
  int OTCsmileAndFlat;
  int OTCpayoff;
  double OTCfwdSmileVisuStd; /* Gets the used smile at the call dates for three
                                strikes calculated with OTCfwdSmileStd */
  int OTCsmile;
  int OTCsmileFee;
  int OTC;

  /* KO */
  int KO_do_optim;
  double **KO_FxGauss;
  double **KO_FxMatrix;
  double **KO_RatesGauss;
  double **KO_RatesMatrix;
  double **KO_Barrier;     // First row is the barrier and the second row is the
                           // adjustment
  int KO_FMC;              // Fast Monte Carlo
  double KO_FMC_precision; // Fast Monte Carlo Precision
  int KO_FMC_min_paths;    // Fast Monte Carlo Min Paths

  /* Random Numbers */
  long Seed;

  /* CORREL */
  int CORRELfirstIndex;
  int CORRELsecondIndex;
  int CORRELfirstLong;
  int CORRELsecondLong;
  long CORRELtStar;

  /* Funding */
  int SpeedFunding;

  /* Which otc are calculated */
  int *OTCcalculate; // 1: calculate ; 2: let GMA interpolate

  /* Change the calibration strategy */
  int use_GMA;
  int LongShort;
} otcpd_params, *OTCPDPARAMS;

void cpd_init_otc_params(
    otcpd_params *otc_params, CPD_STR cpd, int ncall, int nSimul, int Do_pecs,
    int nPoints, int nStd, double CummulPrecision, int CummulLinear, double std,
    int nIter, int PayoffFunction, double fwdSmileVisuNstd,
    /* otc */
    int fwdVolMethod,        /*	0=Sliding 3F vol; 1=Converging 3F vol; 2=Sliding
                            Cvg Sbeta; 3=Cvg Cvg Sbeta */
    int smileOtc,            /*	use smile parameters for the OTC */
    int smileFee,            /*	use smile parameters for the OTC Fees */
    int smileModel, int otc, /*	which call to keep */
    double BMpi,
    /* Correl */
    int firstIndex, int secondIndex, int firstLong, int secondLong,
    long correlTstar, char *CPDsigma, char *CPDalpha, char *CPDbeta,
    char *CPDrho,
    /* Change the Funding for speed */
    int FundingSpeedUp,
    /* Do not calculate all OTC*/
    int nStart, int oneOutOfN,
    /*Fast MC*/
    int FMC_do, double FMC_precision, int FMC_min_paths,
    /*Which GMA*/
    int use_GMA,
    /* Change the calibration strategy */
    int LongShort);

Err cpd_free_otc_params(otcpd_params *otc_params);

typedef struct _otcpd_precalc {
  double **df_const_dom_fx_val;
  double **df_lin_dom_fx_val;
  double **df_const_for_fx_val;
  double **df_lin_for_fx_val;
  double **df_const_pd_pay;
  double **df_lin_pd_pay;
  double **df_const_fund_pay;
  double **df_lin_fund_pay;
  double *df_const_fund_start;
  double *df_lin_fund_start;

  long *fx_mkt_vol_date;
  double *fx_mkt_vol_time;
  double *fx_mkt_smile_alpha;
  double *fx_mkt_smile_beta;
  double *fx_mkt_smile_rho;
  int num_fx_mkt_vol;

  double **pd_fwd3F_vol;
  double **pd_fwdsmile_vol;
  double **pd_fwdsmile_alpha;
  double **pd_fwdsmile_beta;
  double **pd_fwdsmile_rho;
  double **pd_fwdsmile_pi;

  double **FxAdj;
  double *FeeFxAdj;

  double **NotFxAdj;
  double *NotFeeFxAdj;

  double **KO_barrier;
} otcpd_precalc, *OTCPDPRECALC;

void cpd_init_otc_precalc(otcpd_precalc *otc_precalc);

void cpd_free_otc_precalc(otcpd_precalc *otc_precalc, int nCall,
                          int nCoupons_pd, int nCoupons_fund);

Err cpd_alloc_otc_precalc(otcpd_precalc *otc_precalc, int nCall, int nCoupon_pd,
                          int nCoupon_fund, int nMktDates);

Err cpd_otc_precalc(otcpd_precalc *otc_precalc, otcpd_params *OTC_params,
                    CPD_STR cpd, CPD_UND und, SMILE_VOL_MARKET smile_mkt,
                    double *mergeTimes, int mergeNtimes, double *sigDom,
                    double *sigFor, double *sigFx, double *corDF,
                    double *corDFx, double *corFFx);

/***********************************************************************************************************************************
 Calibration of the sum of shifted lognormal parameters to SABR to obtain a
copula based estimation of OTC
***********************************************************************************************************************************/

Err CalibSSLtoSABR(double t,   /* Maturity in years */
                   double fwd, /* Forward */
                   double std, /* # std to calibrate skew */
                   /* SABR Parameters */
                   double sigma, double alpha, double beta, double rho,
                   /*	Results */
                   double *Volup, double *Voldown, double *Fwd_Shift,
                   double *error,
                   /*	Parameters */
                   int nIter);

Err SSLGradient(double strike, double *ParamSSL, double *price,
                double *gradient, int nbParamsSSL);

Err SSLCumulative(double Forward, double Maturity, double VolUp, double VolDown,
                  double FwdShift, int Npoints, int Nstd, double **Cumulative);

Err BMCumulative(double Forward, double Maturity, double BMbeta, double BMfwd1,
                 double BMfwd2, double BMsig1, double BMsig2, double BMpi,
                 int nPoints, int nStd, double **Cumulative);

Err SSLcumul(double Forward, double Maturity, double Strike, double VolUp,
             double VolDown, double Shift, double *Cumul);

Err BMcumul(double Maturity, double BMfwd1, double BMfwd2, double BMsig1,
            double BMsig2, double BMpi, double BMbeta, double strike,
            double strikeBump, double *cumul);

Err SSLprice(double Forward, double Maturity, double Strike, double VolUp,
             double VolDown, double Shift, int IsCall, double *Price);

/* ==================================
Functions for OTC
=================================== */
Err GenerateRealisations(double Mean1, double Mean2, double **cumulative,
                         int Nsimul, int Npoints,
                         int GetGaussOrMarginal, // 0:Gauss / 1:Marginal
                         int FxIndex, int isLinear, double **matrix);

Err OTC_balsam_generation(long nbPaths, /* Must be odd */
                          long nbSteps, double **matrix, long *seed);

Err OTCdfBis(
    double realisation, double lambda, double Tstart, double Tend,
    double Tstar, /* measure Date always equals the last date of the deal */
    double *volTimes, double *vol, int nVolTimes, double *df);

Err OTCdfPhi(double lambda, double Tstart, double *volTimes, double *vol,
             int nVolTimes, double *Phi);

Err OTCdf(double realisation, double lambda, double Tstart, double Tend,
          double Tstar, /* measure Date always equals the OTC Date */
          double Phi, double *df);

Err OTCgetMoments(double endTime, double ldaDom, double ldaFor, double TStar,
                  double *mergeTimes, int mergeNtimes, double *SigmaDom,
                  double *SigmaFor, double *SigmaFx, double *CorrDomFor,
                  double *CorrDomFx, double *CorrForfx,

                  double *dom_std, double *for_fwd, double *for_std,

                  double *dom_for_cov, double *dom_fx_cov, double *for_fx_cov);

Err OTCcvxtAdj(double t, double Tfix, double Tpay, double ldaDom, double ldaFor,
               double *mergeTimes, int mergeNtimes, double *SigmaDom,
               double *SigmaFor, double *SigmaFx, double *CorrDomFor,
               double *CorrDomFx, double *CorrForFx,
               /* Results */
               double *adjustment);

Err OTCcorrelateVariables(double domStd, double forStd, double corDF,
                          double corDFx, double corFFx, long nbPaths,
                          double **inputMatrix, double **outputMatrix);

Err OTCgetFwdSABR(double *forward, CPD_STR cpd, CPD_UND und,
                  SMILE_VOL_MARKET smile_mkt, int otc, double *optionsSigma,
                  double *optionsAlpha, double *optionsBeta, double *optionsRho,
                  double *optionsPi, int noptions, double *mergeTimes,
                  int mergeNtimes, double *sigDom, double ldaDom,
                  double *sigFor, double ldaFor, double *sigFx, double *corDF,
                  double *corDFx, double *corFFx, char *CPDsigma,
                  char *CPDalpha, char *CPDbeta, char *CPDrho, int method);

Err OTCfwdSABR(double forward, double fwdStartTime, double optionsTimeFix,
               double optionsTimeSettle, SMILE_VOL_MARKET smile_mkt,
               SMILE_PARAMETERS smile_params, double *mergeTimes,
               int mergeNtimes, double *sigDom, double ldaDom, double *sigFor,
               double ldaFor, double *sigFx, double *corDF, double *corDFx,
               double *corFFx, int method);

Err OTCgetSMILEparams(CPD_STR cpd, CPD_UND und, double Tfix, double Tval,
                      long DateVal, SMILE_VOL_MARKET smile_mkt,
                      SMILE_PARAMETERS smile_params, int get_full_smile);

Err OTCgetSABRparams(CPD_STR cpd, CPD_UND und, double Tfix, double Tval,
                     long DateVal, long *fx_mkt_vol_date,
                     double *fx_mkt_smile_alpha, double *fx_mkt_smile_beta,
                     double *fx_mkt_smile_rho, int num_fx_mkt_vol,
                     double *smileVol, double *smileAlpha, double *smileBeta,
                     double *smileRho, int smileFee);

Err OTCcalibSmileModel(double Tex,     // Maturity in years
                       double forward, // Forward
                       otcpd_params *OTC_params, SMILE_PARAMETERS smile_params,
                       double *error);

Err OTCgetCumulative(double forward, double Tex, otcpd_params *OTC_params);

Err OTCmatrixRealloc(double ***mat, long nrl, long nhl, long ncl, long nch,
                     long new_nrl, long new_nhl, long new_ncl, long new_nch);

Err OTCgetChangedPayoff(CPD_STR cpd, CPD_UND und, double *NewCoupons);

/* =================================================================================================

                PAYOFF FUNCTIONS

=================================================================================================
*/
Err OTCfwdSABRCalib(CPD_STR cpd, CPD_UND und, double *pd_fwdsmile_vol,
                    double *pd_fwdsmile_alpha, double *pd_fwdsmile_beta,
                    double *pd_fwdsmile_rho, int otc, int nSimul,
                    double **matrix, double tstarTime, int tstarDate,
                    double fwdSmileVisuNstd, double *FxAdj, double **Results);

Err OTCPDpayoffFeeCombinedSmileAndFlat(
    CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
    otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
    SMILE_PARAMETERS smile_params, double pd_not, double **matrix,
    double tstarTime, int tstarDate, double **Results);

Err OTCfundingValuation(double domReal, double forReal, double fxReal,
                        CPD_STR cpd, CPD_UND und, otcpd_precalc *OTC_precalc,
                        double tstarTime, double phiDom, double phiFor, int otc,
                        double dfNotExchge, double dfNotExchgeShort,
                        double dfFirstFunding, double *fundDf, double *fund_leg,
                        double *fund_legShort);

Err OTCNewfundingValuation(double domReal, double forReal, double fxReal,
                           CPD_STR cpd, CPD_UND und, otcpd_precalc *OTC_precalc,
                           double tstarTime, double phiDom, double phiFor,
                           int otc, double dfNotExchge, double dfNotExchgeShort,
                           double dfFirstFunding, double *NewCoupons,
                           double *fundDfDom, double *fundDfFor,
                           double *fund_leg, double *fund_legShort);

Err OTCPDcorrel(CPD_STR cpd, CPD_UND und, double pd_not,
                double *pd_fwdsmile_vol, double *pd_fwdsmile_alpha,
                double *pd_fwdsmile_beta, double *pd_fwdsmile_rho, int smileFee,
                int otc, int nSimul, double **matrix, double tstarTime,
                int tstarDate, double *FxAdj, double *FeeFxAdj, int firstIndex,
                int secondIndex, int firstLong, int secondLong,
                long correlTstar, double **Results);

Err OTCtestPayoff(CPD_UND und, double **matrix, int nSimul, double tstarTime,
                  int tstarDate, int nStrikes, double **Results);

Err fwdSABRpayoff(long today, long forwardFixDate, double forwardFixTime,
                  long valDate, double fixTime, double valTime, double forward,
                  double fwdsmile_vol, double fwdsmile_alpha,
                  double fwdsmile_beta, double fwdsmile_rho, double **matrix,
                  int nSimul, double tstarTime, int tstarDate, char *dom_yc,
                  char *for_yc, double *sigma_time_rates, int sigma_n_rates,
                  double *sigma_dom, double lda_dom, double *sigma_for,
                  double lda_for, double *sigma_time_fx, double *sigma_fx,
                  int sigma_n_fx, double *corr_times, double *correl_dom_for,
                  double *correl_dom_fx, double *correl_for_fx,
                  int corr_n_times, int nStrikes, double *Strikes, int isCall,
                  int OutputVol, double *Results);

Err fwdSABRpayoffAlphaFudge(
    long today, long forwardFixDate, double forwardFixTime, long valDate,
    double fixTime, double valTime, double forward, double fwdsmile_vol,
    double fwdsmile_alpha, double fwdsmile_beta, double fwdsmile_rho,
    double **matrix, int nSimul, double tstarTime, int tstarDate, char *dom_yc,
    char *for_yc, double *sigma_time_rates, int sigma_n_rates,
    double *sigma_dom, double lda_dom, double *sigma_for, double lda_for,
    double *sigma_time_fx, double *sigma_fx, int sigma_n_fx, double *corr_times,
    double *correl_dom_for, double *correl_dom_fx, double *correl_for_fx,
    int corr_n_times, int nStrikes, double *Strikes, int isCall, int OutputVol,
    double SwitchLevel, double Proba, double AlphaFudge, double *Results);

Err KOPDpayoffFee(CPD_STR cpd, CPD_UND und, otcpd_params *OTC_params,
                  otcpd_precalc *OTC_precalc, SMILE_VOL_MARKET smile_mkt,
                  SMILE_PARAMETERS smile_params, long start_date, double pd_not,
                  double tstarTime, int tstarDate, double **Results,
                  double ***save_values);

/* =================================================================================================

                UTILS FUNCTIONS

=================================================================================================
*/
Err OTCgetOptValue(int IsOption, int IsCap, double cap, double Strike,
                   double Notional, double Forward, double Maturity,
                   SMILE_PARAMETERS smile_params, int IsBetaVol, double *Value);

Err OTCgetCopula(double forward, double forwardFixTime, double tstarTime,
                 double *mergeTimes, int mergeNtimes, double *sigDom,
                 double lda_dom, double *sigFor, double lda_for, double *sigFx,
                 double *corDF, double *corDFx, double *corFFx, double smilevol,
                 double smilealpha, double smilebeta, double smilerho,
                 int nSimul, int nPoints, int nIter, int nStd, double std,
                 int smileModel, double **matrix);

Err OTCutilsInterpolate(double X, double XX, double XXX, double F, double FF,
                        double FFF, double Fseek, int method, double *XX1,
                        double *XX2);

Err OTCutilsGetIndex(double *probas, int Nprobas, double prob, int isLinear,
                     int *index1, int *index2, int *index3, int *method);

Err OTCutilsGetFxRealisations(otcpd_params *otc_params);

Err OTCutilsCorrelateFxAndRates(otcpd_params *OTC_params, double meanDom,
                                double stdDom, double meanFor, double stdFor,
                                double corDF, double corDFx, double corFFx);

typedef struct {

  double dForward;
  double dMaturity;
  double dSqMaturity;
  double dProba1;
  double dProba2;

  int iNbStrikes;
  double *dStrikes;
  double *dTargetVols;
  double *dTargetPrices;

  double dTolerance;

  SrtDiffusionType eVolType;
  int iOptimiseOnVol;

  long iMaxIter;
  int iCalibVersion;

} CalibSL_Params, *CALIBSL_PARAMS;

Err alloc_calib_sl_params(int iNbStrikes, CALIBSL_PARAMS sParams);

void free_calib_sl_params(CALIBSL_PARAMS sParams);

Err Calib_MixedSL_ToSmile(CALIBSL_PARAMS sParams, double *dShift,
                          double *dVolDown, double *dVolUp,
                          double *dFittingError);

void OTCutilsSaveMatrix(int iNRows, int iNCols, double **Matrix);

void OTCutilsSaveDouble(char *fileName, double toSave);

Err Srt_smooth_matrix(int method, int num_adj, int NRow, int NCol,
                      double **in_matrix, double **out_matrix);
#endif
