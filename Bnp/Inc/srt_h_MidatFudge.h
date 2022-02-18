#ifndef SRT_H_MIDATFUDGE_H
#define SRT_H_MIDATFUDGE_H

#include "uterror.h"
#include "utstring.h"
#include "uttypes.h"

/*   Main Function of Pricing */
Err srt_f_MidatFudge(
    long    nNumExercise,
    long*   plExercise,
    long*   plExercisePayDates,
    double* pdExercisePremium,
    long    lSwapStart,
    long    lSwapTheoEnd,
    char*   szSwapFreq,
    char*   szSwapBasis,
    double  dCouponRate,
    double  dFundingMargin,
    String  szRefRate,
    String  szPayRec,
    String  szMarket,
    String  szYieldcurve,
    Err (*pfGetVol)(
        double  dStart,
        double  dEnd,
        double  dStrike,
        double  dForward,
        double  dSpread,
        double* pdBsVol),
    String szVolType,
    String szModel,
    /* Outputs from the addin */
    double*** pdAnswer);

/* Type of Vol */
typedef enum
{
    VOLLOGNORMAL,
    VOLNORMAL,
    VOLBETA,
    VOLSABR
} VOL_TYPE;

/* Type of Model */
typedef enum
{
    BASIC,     /* CONVERGING LOGNORMAL without optimisation of the Ex Boundary */
    SLIDDING,  /* For a given dynamic of the smile forward we gonna */
    CONVERGING /* optimise the exercise boundary */
} MIDATFUDGE_MODEL;

/* Swaption in Midat Structure */
typedef struct
{
    /* time from today / exercise date */
    double exp_time;
    /* Start Date */
    long start;
    /* TheoEnd */
    long theoend;
    /* Strike of the Swaption */
    double coupon_rate;
    /* Volatility which would preaval at today / exercise date */
    double vol;
    double beta;
    double alpha;
    double rho;
    /* Swap Rate */
    double fwdswap;
    /* vol type */
    VOL_TYPE voltype;
} SWAPTION_MIDAT;

/* Midat Structure */
typedef struct
{
    /* Number of Exercise Dates */
    long ntEx;
    /* Exercise Dates */
    long* tEx;
    /* Exercise Pay Dates */
    long* tExPay;
    /* Exercise Premiums */
    double* tExPremium;
    /* Number of swaption per Exercise */
    long* nSwaption;
    /* Swaption at exercise Date */
    SWAPTION_MIDAT** Swaption;
    /* Most expensive Swaption (on average) at exercise date */
    long* MostExp;
    /* Exercise boundary (in term of swaprate) */
    double* ExBound;
    /* RefRate */
    char* refrate;
    /* Swap Freq */
    SrtCompounding SwapFreq;
    /* Swap Basis */
    SrtBasisCode SwapBasis;
    /* Model */
    MIDATFUDGE_MODEL model;
} MIDAT_FUDGE;

Err srt_f_initMIDAT_FUDGE(
    long    nNumEx,
    long*   Ex,
    long*   ExPay,
    double* ExPremium,
    long    lStart,
    long    lTheoEnd,
    double  strike,
    String  szSwapBasis,
    String  szSwapFreq,
    String  szRefRate,
    String  szYieldcurve,
    Err (*pfGetVol)(
        double  dStart,
        double  dEnd,
        double  dStrike,
        double  dForward,
        double  dSpread,
        double* pdBsVol),
    String       szVolType,
    String       szModel,
    MIDAT_FUDGE* midat);

Err srt_f_freeMIDAT_FUDGE(MIDAT_FUDGE* midat);

Err interp_MIDAFUDGE_model(char* str, MIDATFUDGE_MODEL* val);

#endif