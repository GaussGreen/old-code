
#include "utallhdr.h"
#include <NUM_H_ALLHDR.H>

/* ========================================================================== */

#ifndef SRT_H_RESETABLE_H
#define SRT_H_RESETABLE_H

/* structure made to pass elements from the wrapper to the c funtion
for the RangeSwapMinCoupon*/
typedef struct {
  double *Strikes;
  int NumStrikesInVol;
  SRT_Boolean AdjForSpread;
} elemforCms;

Err srt_f_mertonfwdvols(int seedo, double smiledate, int n_resetdates,
                        double *resetdates, int n_strikes, double *strikes,
                        double *initparam, char *UseDRS, char *cRefRateCode,
                        char *cMarketId, char *szYieldCurveName,
                        char *szVolCurveName,
                        Err (*GetVol)(double dStart, double dEnd,
                                      double dStrike, SRT_Boolean bAdjForSpread,
                                      double dForward, double dSpread,
                                      double *pdBsVol),
                        elemforCms passCms, double *impvols);
/* =============================================================================
  FUNCTION     : srt_f_resetablejumps(...)
 ============================================================================ */

Err srt_f_resetable(int n_resetdates, double *resetdates, double *initparam,
                    char *FloatCoupon, double margin, double width,
                    char *cRefRateCode, char *cMarketId, char *szYieldCurveName,
                    char *szVolCurveName, double *pv);

/* =============================================================================
  FUNCTION     : NewtonJumps
 ============================================================================ */

Err NewtonJumps(
    int nmat,    /* = 1+nbre de digitales dans la periode de reset */
    double *DRS, /* term structure des DRS a la date de reset: [1..nmat] */
    double *maturites, /* [1..nmat] */
    double *isFriday,  /* [1..nmat] flag qui indique si chaque	maturite  est
                          un Friday */
    double *
        *param,   /*parametres de la dynamique de chaque DRS [1..nmat][1..5] */
    double width, /* Largeur du corridor */
    double callspread, double *eps, double *prix);

/* =============================================================================
  FUNCTION     : srt_f_timeswapjumps
 ============================================================================ */

Err srt_f_timeswapjumps(int n_resetdates, double *resetdates, double *initparam,
                        char *FloatCoupon, double margin, double upperbarrier,
                        double lowerbarrier, char *cRefRateCode,
                        char *cMarketId, char *szYieldCurveName,
                        char *szVolCurveName, double *pv);
/* =============================================================================
  FUNCTION     : PrixTimeSwap
 ============================================================================ */

/* Computes the price of a Time Swap with a fixed range*/
double PrixTimeSwap(
    double eps, double *DRS, /* Term structure of the DRS. tab[1..nmat] */
    int nmat,         /* Number of dates +1 (the first date of the period)*/
    double *mat,      /* tab[1..nmat] */
    double *IsFriday, /*tab[1..nmat] */
    double *
        *param, /*tab[1..nmat][1..5] *. term structure of the jump parameters */
    double width, /* width of the range */
    double callspread);

/* =============================================================================
  FUNCTION     : srt_f_FlooredTimeSwap(...)
 ============================================================================ */

Err srt_f_FlooredTimeSwap(
    long Startdate, long Enddate, double coupon, double minimum, double UpBa,
    double LoBa, double epsilon, int lag, int endlag, double lowvolshift,
    double upvolshift, char *FloatingLeg, char *FloatingMinimum,
    char *AccRefRateCode, double corrup, double corrdown, char *Approx,
    char *TimeSwap, char *dailyorperiod, char *WeekendRule, char *RecPay,
    String sYieldCurveName, String sVolCurveName, SrtDiffusionType CMSlognorm,
    String cMarketId, String cRefRateCode,
    char *(*getLogVol)(char *szVolCurve, double dStartDate, double dEndDate,
                       double dStrike, double *dVol, double *dPower),
    elemforCms passCms, double *pv);

Err srt_f_sabrmontecarlo(double forw, double vovol, double beta, double rho,
                         double num_paths, double num_steps, double sigma,
                         double maturity, double *strike, int num_strikes,
                         int modeltype, int sampletype, double *impvol);
/* =============================================================================
  FUNCTION     : optmertonpremiumtimedependent(...)
 ============================================================================ */
Err optmertonpremiumtimedependent(double dFwd, double Strike, int n_periods,
                                  double *sigma, double *U1, double *lambda1,
                                  double *U2, double *lambda2, double *T,
                                  char *Logornorm, double *answer);

Err optmertonsmiletimedependent(double dFwd, int n_strikes, double *strikes,
                                int n_periods, /* The 6 following variables are
                                                  vectors of length nperiods */
                                double *param, double *mat, char *Logornorm,
                                double *impvols);

Err optcalibmertontimedep(double dFwd_loc, int n_strikes_loc, int n_periods_loc,
                          double *strikes_loc, double *market_vols,
                          double *param, double *maturity, char *Logornorm,
                          double *chisq, long *ia, double *calibvols);

Err srt_f_mertonfwdvols_td(int seedo, double smiledate, int n_resetdates,
                           double *resetdates, int n_strikes, double *strikes,
                           double *initparam, char *cRefRateCode,
                           char *cMarketId, char *szYieldCurveName,
                           char *szVolCurveName, double *res);
/**************************************************************************************
FUNCTION:	srt_f_resetable_td(...)
***************************************************************************************/
Err srt_f_resetable_td(int n_resetdates, double *resetdates, double *initparam,
                       char *FloatCoupon, double margin, double width,
                       double callspread, char *ResetType, char *cRefRateCode,
                       char *cMarketId, char *szVolCurveName,
                       char *szYieldCurveName, double **pv);

Err srt_f_simpleresetable(
    double StartDate, double EndDate,
    char *VolType, /*Flat or Sliding or Input*/
    double *
        VolsInput, /*If above is "Input"  , this contains 4 vols:lowercallspread
                      , lower  , upper  , uppercallspread*/
    char *FloatCoupon, SrtCompounding ResetFrequency, double margin,
    double width, double callspread, double VolShift, char *ResetType,
    double lowerbarrier, char *cRefRateCode, char *cMarketId,
    char *szVolCurveName, char *szYieldCurveName, double *pv);

Err srt_f_simpleresetcap(
    long StrikeDate, long StartDate, long EndDate,
    char *VolType,   /*Flat or NormSliding or LogSliding or Input*/
    double VolInput, /*Contains a vol input buy the user used only if above is
                        input*/
    char *cRefRateCode, char *cMarketId, char *szVolCurveName,
    char *szYieldCurveName, SrtCallPutType CallPut, double *pv);

#endif