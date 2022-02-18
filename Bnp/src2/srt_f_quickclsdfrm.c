/* ===================================================================================
   FILENAME:      srt_f_quickclsdfrm.c

   PURPOSE:       Compute quick prices for saptions, cap floors for a generic Cheyette
                  Beta model
                                  the formula is a complete approximation and offers absolutely
                                  no accuracy.
                                  It is just used in the calibration with FIXED POINT as a way of
                                  getting a Beta sensitivity...
   =================================================================================== */
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

/* --------------------------------------------------------------------------------- */

double srt_f_quick_beta_bond_option(
    int             n,
    double          bond_strike,
    TermStruct*     ts,
    double          fixing_time,
    double*         coupon,
    double*         pay_time,
    double*         df,
    SrtReceiverType rec_pay,
    SrtMdlType      mdl_type,
    SrtMdlDim       eModelDim,
    String          szYieldCurveName)
{
    Err         err = NULL;
    double**    ppdSigmaCurve;
    long        lNumSigmas;
    long        lNumSigmaCols;
    double**    ppdTauCurve;
    long        lNumTaus;
    long        lNumTauCols;
    double      dBeta;
    double      dBeta1;
    double      dBeta2;
    double      dFra;
    TermStruct* psLgmTermStruct;
    SrtCurvePtr psYieldCurve;
    double      dLgmPrice;
    double      dLevel;
    double      dForward;
    double      dStrike;
    double      dImpliedBsnVol;
    double      dEqBetaVol;
    double      dEqBsVol;
    double      dIntrinsic;
    long        lToday;
    double      dPrice;
    double      dMaturity;
    int         i;

    /* Get the original Term Structure */
    if (eModelDim == TWO_FAC)
    {
        lNumSigmaCols = 6;
        ppdSigmaCurve = (double**)malloc(lNumSigmaCols * sizeof(double*));
        lNumTauCols   = 3;
        ppdTauCurve   = (double**)malloc(lNumTauCols * sizeof(double*));

        err = srt_f_display_IRM_TwoFac_TermStruct(
            ts,
            &(ppdSigmaCurve[0]),
            &(ppdSigmaCurve[1]),
            &(ppdSigmaCurve[2]),
            &(ppdSigmaCurve[3]),
            &(ppdSigmaCurve[4]),
            &(ppdSigmaCurve[5]),
            &lNumSigmas,
            &(ppdTauCurve[0]),
            &(ppdTauCurve[1]),
            &(ppdTauCurve[3]),
            &lNumTaus);
    }
    else
    {
        lNumSigmaCols = 6;
        ppdSigmaCurve = (double**)malloc(lNumSigmaCols * sizeof(double*));
        lNumTauCols   = 2;
        ppdTauCurve   = (double**)malloc(lNumTauCols * sizeof(double*));

        err = srt_f_display_IRM_OneFac_TermStruct(
            ts,
            &(ppdSigmaCurve[0]),
            &(ppdSigmaCurve[1]),
            &(ppdSigmaCurve[2]),
            &(ppdSigmaCurve[3]),
            &(ppdSigmaCurve[4]),
            &(ppdSigmaCurve[5]),
            &lNumSigmas,
            &(ppdTauCurve[0]),
            &(ppdTauCurve[1]),
            &lNumTaus);
    }

    /* Rebuilds a fake LGM TS by a quick Normal approximation of the Volatility */
    for (i = 0; i < lNumSigmas; i++)
    {
        /* Computes the FRA at the volatility date */
        dFra = swp_f_zr(
            (Ddate)(ppdSigmaCurve[0][i]), (Ddate)(ppdSigmaCurve[0][i]) + 14.0, szYieldCurveName);

        /* Gets the Beta */
        dBeta = (Ddate)(ppdSigmaCurve[2][i]);

        /* The rough normal vol is the Beta vol times the FRA power Beta */
        ppdSigmaCurve[1][i] *= pow(dFra, dBeta);

        if (eModelDim == TWO_FAC)
        {
            dBeta = (Ddate)(ppdSigmaCurve[4][i]);
            ppdSigmaCurve[3][i] *= pow(dFra, dBeta);
        }
    }

    /* Initialises the corresponding LGM Term Struct */
    psYieldCurve = lookup_curve(szYieldCurveName);
    lToday       = get_clcndate_from_yldcrv(psYieldCurve);
    err          = srt_f_init_IRM_TermStruct(
        lToday,
        ppdSigmaCurve,
        lNumSigmaCols,
        lNumSigmas,
        ppdTauCurve,
        lNumTauCols,
        lNumTaus,
        mdl_type,
        eModelDim,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0, /* vasicek parms */
        0,
        0,
        NULL,
        &psLgmTermStruct);
    if (err)
    {
        free(ppdSigmaCurve);
        free(ppdTauCurve);
    }

    /* Price the deal with this LGM Term Structure */
    dLgmPrice = srt_f_lgm_coupon_bond_option(
        n, bond_strike, psLgmTermStruct, fixing_time, coupon, pay_time, df, rec_pay, eModelDim);

    /* Compute a rough estimate of the FOrward and the Strike */
    dLevel  = 0.0;
    dStrike = 0.0;
    for (i = 1; i <= n; i++)
    {
        dLevel += df[i];
        dStrike += coupon[i];
    }
    dForward = (df[0] - df[n]) / dLevel;
    dStrike /= n;

    /* Gets a rough estimate of the equivalent BSN vol */
    dMaturity = (fixing_time - lToday) * YEARS_IN_DAY;
    dLgmPrice /= dLevel;
    dIntrinsic = dForward - dStrike;
    if (rec_pay == SRT_RECEIVER)
        dIntrinsic *= -1.0;
    if (dLgmPrice < dIntrinsic)
    {
        dImpliedBsnVol = 0.0;
    }
    else
    {
        err = srt_f_optimpvol(
            dLgmPrice,
            dForward,
            dStrike,
            dMaturity,
            1.0,
            (rec_pay == SRT_RECEIVER ? SRT_PUT : SRT_CALL),
            SRT_NORMAL,
            &dImpliedBsnVol);
        if (err)
            return 0.0;
    }

    /* Transforms this vol into an equivalent Beta Vol (through BS) for ATM option */
    if (eModelDim == ONE_FAC)
        dBeta = find_beta(fixing_time, ts);
    else
    {
        find_tf_beta(fixing_time, ts, &dBeta1, &dBeta2);
        dBeta = 0.5 * (dBeta1 + dBeta2);
    }

    err = srt_f_optbetavoltoblkvol(dForward, dForward, dImpliedBsnVol, dMaturity, 0.0, &dEqBsVol);

    err = srt_f_optblkvoltobetavol(dForward, dForward, dEqBsVol, dMaturity, dBeta, &dEqBetaVol);

    /* Reprices the Deal with this BEta Vol (with the right strike) */
    dPrice = srt_f_optblkschbetaquick(
        dForward,
        dStrike,
        dEqBetaVol,
        dMaturity,
        dBeta,
        1.0,
        (rec_pay == SRT_RECEIVER ? SRT_PUT : SRT_CALL),
        PREMIUM);

    /* That's it */

    return dPrice;
}

/* -------------------------------------------------------------------------------- */

double srt_f_quick_beta_capfloor(
    TermStruct*     ts,
    double*         fixing_time,
    double*         period_time, /* period_time[0] is the strike payment time */
    double*         df,
    double*         payment, /* This is cvg * ( cash_fwd + spread ) */
    int             num_caplets,
    SrtReceiverType rec_pay,
    SrtMdlType      mdl_type,
    SrtMdlDim       eModelDim,
    String          yc_name)
{
    double  price;
    double  caplet;
    int     i;
    double* coupon;

    /*	We think of a cap as a sum of zero coupon bond options */
    price     = 0.0;
    coupon    = (double*)srt_malloc(2 * sizeof(double));
    coupon[0] = coupon[1] = 0.0;
    for (i = 0; i < num_caplets; i++)
    {
        /* Skip caplets that fixed today or before */
        if (fixing_time[i] > 0)
        {
            caplet = srt_f_quick_beta_bond_option(
                1,
                1.0 / (1.0 + payment[i + 1]),
                ts,
                fixing_time[i],
                coupon,
                &(period_time[i]),
                &(df[i]),
                rec_pay,
                mdl_type,
                eModelDim,
                yc_name);
            caplet *= (1 + payment[i + 1]);
            price += caplet;
        }
    }

    return price;
}
