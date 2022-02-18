#ifndef SWP_H_CMS_H
#define SWP_H_CMS_H

/* These are the original version of the CMS rate using Flat vol and swaption replication
        profile */

Err swp_f_cmswp(
    SwapDP*          swapdp,
    String           refRatecode,
    double           strike,
    double           volatility,
    SrtReceiverType  pay_or_receive,
    Date             pay_date,
    double           delta,
    int              nswaps,
    String           ycName,
    SrtDiffusionType lognrm_nrm,
    double*          ans);

Err swp_f_cmsrate(
    double           forward,
    double           num_periods,
    SrtCompounding   frequency,
    double           volatility,
    double           maturity,
    double           delay,
    double           delta,
    int              num_swaps,
    SrtDiffusionType lognrm_nrm,
    double*          ans);

Err swp_f_cmsratefwd(
    SwapDP*          swapdp,
    String           refRatecode,
    double           volatility,
    Date             pay_date,
    double           delta,
    int              nswaps,
    String           ycName,
    SrtDiffusionType lognrm_nrm,
    double*          ans);

Err swp_f_cmsratedp(
    SwapDP*          swapdp,
    Date             today,
    double           forward,
    double           volatility,
    Date             pay_date,
    double           delta,
    int              nswaps,
    SrtDiffusionType lognrm_nrm,
    double*          ans);

Err swp_f_cmsratedp_clsdfrm(
    SwapDP*          swapdp,
    Date             today,
    double           forward,
    double           volatility,
    Date             pay_date,
    SrtDiffusionType lognrm_nrm,
    double*          ans);

Err swp_f_Cms_Rate(
    double           dFwdSwapRate, /* Forward Swap Rate */
    double           dMaturity,    /* Cms Maturity as double */
    double           dNumPeriods,
    double           dFrequency, /* Cms Frequency */
    double           dDelay,
    double           dRateConv, /* date adjustment */
    SrtDiffusionType VolType,   /* Vol type Lognormal or Normal */
    double           dFlatVol,  /* Flat Vol if used */
    int              iMethod,   /* 0: Use Flat Vol, 1: linear interpolation, 2: FullSmile*/
    Date             dStart,    /* The following parameters are used in the GetVol function */
    Date             dEnd,      /* and are useless in the rest of the code */
    SRT_Boolean      bAdjForSpread,
    double           dSpread,
    char*            szVolCurveName,
    long             lNumStrikesInVol,
    double*          pdStrikesVol,
    double*          dCmsRateValue);

Err swp_f_Tec_Rate(
    double           dFwdSwapRate, /*Forward Swap Rate */
    double           dMaturity,
    double           dNumPeriods,
    double           dAtmStrike,
    double           dMargin,
    double           dFrequency,    /* Tec Frequency */
    double           dPaymentPower, /* Tec Power */
    double           dDelay,
    double           dRateConv, /* date adjustment */
    SrtDiffusionType VolType,
    double           dFlatVol,
    int              iMethod, /* 0: Use Flat Vol, 1: linear interpolation, 2: FullSmile*/
    Date             dStart,  /* The following parameters are used in the GetVol function */
    Date             dEnd,    /* and are useless in the rest of the code */
    SRT_Boolean      bAdjForSpread,
    double           dSpread,
    char*            szVolCurveName,
    long             lNumStrikesInVol,
    double*          pdStrikesVol,
    double*          dTecRateValue);

double swp_f_cmsrateNewYork(
    double maturity,
    double underlying,
    double forward,
    double ATMvol,
    double alpha,
    double beta,
    double rho,
    double numperiod,
    double precision);
#endif
