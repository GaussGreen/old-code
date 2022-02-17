/* ===================================================================================
   FILENAME:      swp_h_swaption.h

   PURPOSE:       Compute swaption prices and greeks  , and implied volatilities
   ===================================================================================
 */

#ifndef SWP_H_SWAPTION_H
#define SWP_H_SWAPTION_H

/* -------------------  SWAPTION PRICE AND DERIVATIVES   -----------------------
 */

/*  ------------------------------------------------------------------------------
 */
/* Compute the price or greeks of a swaption (for use outside SORT) */
Err swp_f_Swaption(long start, long end_nfp, String compStr, String basisStr,
                   double vol, double strike, String recPayStr,
                   String refRateCode, String ycName, String greekStr,
                   String normLogStr, double *result);

/* -----------------------------------------------------------------------------------
 */
/* Compute the price or greeks of a swaption (for use inside SORT)
   Caution: swapdp->spot_lag should be set */

Err swp_f_Swaption_SwapDP(SwapDP *swapdp, double vol, double strike,
                          SrtReceiverType recPay, String refRateCode,
                          String ycName, SrtGreekType greek,
                          SrtDiffusionType logNorm, double *result);

/* ----------------------  IMPLIED VOLATILITY FOR SWAPTIONS
 * ------------------------ */

/* Compute the implied volatility given the price of a swaption (for use outside
 * SROT) */

Err swp_f_SwaptionImpliedVol(double premium, long start, long end_nfp,
                             String compStr, String basisStr, double strike,
                             String recPayStr, String refRateCode,
                             String ycName, String normLogStr, double *impvol);

/* Compute the implied volatility given the price of a swaption (for use inside
 * SROT) */

Err swp_f_SwaptionImpliedVol_SwapDP(double premium, SwapDP *swapdp,
                                    double strike, SrtReceiverType recPay,
                                    String refRateCode, String ycName,
                                    SrtDiffusionType logNorm, double *impvol);

/* -------------------------  CASH SETTLED SWAPTIONS
 * -------------------------- */

/* Compute the premium of a cash settled swaption (for use outside SORT) */
Err swp_f_CashSettledSwaption(long start, long nfp, String strComp,
                              String strBasis, double strike, double vol,
                              String strRecPay, String refRateCode,
                              String ycName, String strGreek, String strNormLog,
                              double precision, double *answer);

/* Compute the premium of a cash settled swaption (for use inside SORT) */
Err swp_f_CashSettledSwaption_SwapDP(SwapDP *swapdp, double vol, double strike,
                                     SrtReceiverType rec_pay,
                                     String refRateCode, String ycName,
                                     SrtGreekType greek,
                                     SrtDiffusionType norm_or_log,
                                     double precision, double *answer);

/* -------------------------  CASH SETTLED SWAPTIONS
 * -------------------------- */

/* Compute the premium of a cash settled swaption (for use outside SORT) */
Err swp_f_QuantoCashSettledSwaption(long start, long nfp, String strComp,
                                    String strBasis, double strike, double vol,
                                    double correlationfxcms, double volfx,
                                    String strRecPay, String refRateCode,
                                    String ycName, String strGreek,
                                    String strNormLog, double precision,
                                    double *answer);

/* Compute the premium of a cash settled swaption (for use inside SORT) */
Err swp_f_QuantoCashSettledSwaption_SwapDP(
    SwapDP *swapdp, double vol, double strike, double correlationfxcms,
    double volfx, SrtReceiverType rec_pay, String refRateCode, String ycName,
    SrtGreekType greek, SrtDiffusionType norm_or_log, double precision,
    double *answer);

/*	Alan -- Compute the price of a quanto swaption  ,
        Mad spliting -- to include vector of end dates */

Err swp_f_QuantoSwaptionPrice(
    double dFwdSwapRate, double dStrike, SrtCallPutType SrtPayRec, Date dToday,
    Date dDomExerciceDate, Date dForFixingDate, long dStartDate,
    Date *pdEndDates, double *pdCvg, int iNumDates, double dCmsNumOfPeriods,
    SrtCompounding SrtForFrequency, SrtDiffusionType SrtDomVolType,
    SrtDiffusionType SrtForVolType, double *pdForCmsVol, double *pdVolForwardFx,
    double *pdCorrelationCmsForwardFx, long lNumOfVolsAndCorrels,
    char *sDomYCName, SRT_Boolean bAdjustCmsVol, double *dQuantoSwaption);

/*	Alan -- Compute the price of a quanto swaption  ,
        West version -- Compute all details */

Err swp_f_QuantoSwaption(
    double dFwdSwapRate, double dStrike, char *sPayRec, long dStartDate,
    long dEndDate, double dCmsNumOfPeriods, char *sDomCompounding,
    char *sDomBasis, char *sForCompounding, SrtDiffusionType SrtDomVolType,
    SrtDiffusionType SrtForVolType, double *pdForCmsVol, double *pdVolForwardFx,
    double *pdCorrelationCmsForwardFx, long lNumOfVolsAndCorrels,
    char *sDomYCName, char *sForYCName, SRT_Boolean bAdjustCmsVol,
    double *dQuantoSwaption);

#endif