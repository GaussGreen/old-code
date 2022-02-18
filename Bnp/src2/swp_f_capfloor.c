/* =============================================================================

   FILENAME:    swp_f_capfloor.c

  PURPOSE:     Calculation for caps/floors: premium, implied volatilities...

   ============================================================================= */

#include "math.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_capfloor.h"
#include "swp_h_swap_generic.h"

#define MAXITERATION 100
#define PREM_TOL_PCT2 0.001
#define PREM_TOL2 (double)1e-20

/* -------------------  CAP/FLOOR PRICE AND DERIVATIVES   ---------------------- */

/* ----------------------------------------------------------------------------- */
/* Compute the price or Greeks of a cap/floor on a Reference Rate (use outside SORT)*/
Err swp_f_CapFloor(
    long   start,
    long   end_nfp,
    double strike,
    Err (*GetBSVol)(
        double start, double end, double strike, double dForward, double dSpread, double* vol),
    String  capFloorStr,
    String  szRefRateCode,
    String  ycName,
    String  greekStr,
    String  logNormStr,
    double* result)
{
    Err              err = NULL;
    SwapDP           swapdp;
    SrtReceiverType  capFloor;
    SrtGreekType     greek;
    SrtDiffusionType logNorm;

    /* Get the compounding and the basis of the cap from the reference rate */
    err = swp_f_get_ref_rate_details(szRefRateCode, &swapdp.basis_code, &swapdp.compd);
    if (err)
        return err;

    /* Set the SwapDP for the cap/floor */
    err = swp_f_setSwapDP(start, end_nfp, swapdp.compd, swapdp.basis_code, &swapdp);
    if (err)
        return err;

    /* Modify the capFloorStr into a type */
    err = interp_rec_pay(capFloorStr, &capFloor);
    if (err)
        return err;

    /* Modify the GreekString into a type */
    greek = PREMIUM;
    if (greekStr)
    {
        err = interp_greeks(greekStr, &greek);
        if (err)
            return err;
    }

    /* Get the Black Scholes type: normal or lognormal */
    logNorm = SRT_LOGNORMAL;
    if (logNormStr)
    {
        err = interp_diffusion_type(logNormStr, &logNorm);
        if (err)
            return err;
    }

    /* Main call to the pricing function for caps/floors */
    err = swp_f_CapFloor_SwapDP(
        &swapdp, strike, GetBSVol, capFloor, szRefRateCode, ycName, greek, logNorm, result);
    if (err)
        return err;

    /* Return a success message */
    return NULL;

} /* END swp_f_CapFloor(...) */

/* ----------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------------- */
/* Compute the price or Greeks of a cap/floor on a Reference Rate (use inside SORT)*/
Err swp_f_CapFloor_SwapDP(
    SwapDP* swapdp,
    double  strike,
    Err (*GetBSVol)(
        double start, double end, double strike, double dForward, double dSpread, double* vol),
    SrtReceiverType  capFloor,
    String           szRefRateCode,
    String           ycName,
    SrtGreekType     greek,
    SrtDiffusionType logNorm,
    double*          result)
{
    Err         err = NULL;
    Date        today;
    GenSwapLeg  floatleg;
    double      forward;
    double      vol;
    double      mat;
    double      df;
    double      cvg;
    double      caplet;
    int         i;
    SrtCurvePtr yccrv;

    /* Get the curve */
    yccrv = lookup_curve(ycName);
    if (!yccrv)
        return serror("Could not find yc %s in swp_f_Swaption", ycName);

    /* Get today from Yield Curve */
    today = get_today_from_curve(yccrv);

    /* Make a floating leg corresponding to the cap */
    err =
        swp_f_make_FloatingLeg(swapdp->start, swapdp->end, today, szRefRateCode, yccrv, &floatleg);
    if (err)
        return err;

    /* Computes the discount factors for each payment date */
    err = df_list(floatleg.pay_date, ycName, &(floatleg.df));

    /* Sets the maturities from today to fixing date (spot lag before spot) */
    err = maturity_list(floatleg.fixing_date, today, &(floatleg.mat));
    if (err)
        return err;

    /* Loop on all the caplets */
    *result = 0.0;
    for (i = 0; i < floatleg.leg_length - 1; i++)
    {
        /* If fixing is before today, discard the caplet (no fixings yet) */
        mat = floatleg.mat.d[i];
        if (mat <= 0.0)
            continue;

        /* The forward rate incorporates the fwd_cash and the spread */
        forward = floatleg.fwd.d[i] + floatleg.spread.d[i];

        /* Get the volatility for the period (start, end) and the strike */
        err = GetBSVol(
            floatleg.start_date.date[i],
            floatleg.end_date.date[i],
            strike,
            forward,
            floatleg.spread.d[i],
            &vol);

        if (err)
            return err;

        /* Get coverage and discount factor */
        cvg = floatleg.cvg.d[i];
        df  = floatleg.df.d[i + 1];

        /* Call to the BlackScholes type formulas */
        err = swp_f_Caplet_Struct(forward, strike, vol, mat, capFloor, greek, logNorm, &caplet);
        if (err)
            return err;

        /* Increase the result by the result for currant caplet */
        *result += caplet * cvg * df;

    } /* END for (i=0; i < floatleg.leg_length) loop on caplets */

    /* Frees the allocated memory */
    err = swp_f_freein_GenSwapLeg(&floatleg);
    if (err)
        return err;

    /* Return a success message */
    return NULL;

} /* END Err swp_f_CapFloor_SwapDP(...) */

/* ----------------------------------------------------------------------------- */
/* Call to the relevant BlackScholes formula for the caplet pricing */
Err swp_f_Caplet_Struct(
    double           forward,
    double           strike,
    double           vol,
    double           mat,
    SrtReceiverType  capFloor,
    SrtGreekType     greek,
    SrtDiffusionType logNorm,
    double*          premium)
{
    Err            err = NULL;
    SrtCallPutType callPut;

    /* Set call_put according to cap or floor: cap is a call on rate */
    callPut = (capFloor == SRT_PAYER ? SRT_CALL : SRT_PUT);

    /* Call the relevant BS pricing function */
    if (logNorm == SRT_LOGNORMAL)
    {
        *premium = srt_f_optblksch(forward, strike, vol, mat, 1.0, callPut, greek);
    }
    else
    {
        *premium = srt_f_optblknrm(forward, strike, vol, mat, 1.0, callPut, greek);
    }

    /* Return a success message */
    return NULL;

} /* END Err swp_f_Caplet_Struct (...) */

/* ----------------------------------------------------------------------------- */

/* -----------------  IMPLIED VOLATILITY FOR CAP/FLOOR   ----------------------- */

/* ----------------------------------------------------------------------------- */
/* Computes the implied vol of a cap/floor using Newton (for use outside SORT) */
Err swp_f_CapFloorImpliedVol(
    double  premium,
    long    start,
    long    end_nfp,
    double  strike,
    String  capFloorStr,
    String  refRateCode,
    String  ycName,
    String  logNormStr,
    double* impvol)
{
    Err              err = NULL;
    SwapDP           swapdp;
    SrtReceiverType  capFloor;
    SrtDiffusionType logNorm;

    /* Get the compounding and the basis of the cap from the reference rate */
    err = swp_f_get_ref_rate_details(refRateCode, &swapdp.basis_code, &swapdp.compd);
    if (err)
        return err;

    /* Set the SwapDP for the cap/floor */
    err = swp_f_setSwapDP(start, end_nfp, swapdp.compd, swapdp.basis_code, &swapdp);
    if (err)
        return err;

    /* Modify the capFloorStr into a type */
    err = interp_rec_pay(capFloorStr, &capFloor);
    if (err)
        return err;

    /* Get the Black Scholes type: normal or lognormal */
    logNorm = SRT_LOGNORMAL;
    if (logNormStr)
    {
        err = interp_diffusion_type(logNormStr, &logNorm);
        if (err)
            return err;
    }

    /* Main call to the implied volatility function for caps/floors */
    err = swp_f_CapFloorImpVol_SwapDP(
        premium, &swapdp, strike, capFloor, refRateCode, ycName, logNorm, impvol);
    if (err)
        return err;

    /* Return a success message */
    return NULL;

} /* END swp_f_CapFloorImpliedVol(...) */

/* ----------------------------------------------------------------------------- */

/* This function cheets, and returns a null vol in *vol , no matter
   the start date, the end date and the strike */
static Err NullVol(
    double start, double end, double strike, double dForward, double dSpread, double* vol)
{
    *vol = 0.0;
    return NULL;
}

static double _capfloorconstantvol;

/* This function cheets, and returns the same vol in *vol , no matter
   the start date, the end date and the strike ( the vol stored in _cafloorconstantvol )*/
static Err ConstantVol(
    double start, double end, double strike, double dForward, double dSpread, double* vol)
{
    *vol = _capfloorconstantvol;
    return NULL;
}
/* ----------------------------------------------------------------------------- */

/* Computes the implied vol of a cap/floor using Newton (for use inside SORT) */
Err swp_f_CapFloorImpVol_SwapDP(
    double           premium,
    SwapDP*          swapdp,
    double           strike,
    SrtReceiverType  capFloor,
    String           refRateCode,
    String           ycName,
    SrtDiffusionType logNorm,
    double*          impvol)
{
    Err         err = NULL;
    Date        today;
    double      intrinsic;
    double      vol_guess;
    double      vol_shift;
    double      tmp_vol;
    double      prem_new;
    double      deriv;
    double      prem_shift;
    int         i;
    SrtCurvePtr yccrv;

    /* Get the yield curve */
    yccrv = lookup_curve(ycName);
    if (!yccrv)
        return serror("Could not find yc %s in swp_f_Swaption", ycName);

    /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
    swapdp->spot_lag = get_spotlag_from_curve(yccrv);

    /* Get today from Yield Curve */
    today = get_today_from_curve(yccrv);

    /* Compute the intrinsic value of the option, with a null volatility */
    err = swp_f_CapFloor_SwapDP(
        swapdp, strike, NullVol, capFloor, refRateCode, ycName, PREMIUM, logNorm, &intrinsic);
    if (err)
        return err;

    /* Make sure the premium is above the intrinsic */
    if (premium < intrinsic)
        return serror("Premium below intrinsic in swp_f_CapFloorImpVol_SwapDP");

    /* Sets initial guess for vol */
    if (logNorm == SRT_NORMAL)
    {
        vol_guess = .01;
        vol_shift = 0.00000001;
    }
    else if (logNorm == SRT_LOGNORMAL)
    {
        vol_guess = .10;
        vol_shift = 0.0000001;
    }

    /* Start Newton like iterations (maximum of 20 iterations) */
    for (i = 0; i <= 20; i++)
    {
        /* Set the static volatility to the volatility guess */
        _capfloorconstantvol = vol_guess;

        /* Compute the premium with this new volatility */
        err = swp_f_CapFloor_SwapDP(
            swapdp,
            strike,
            ConstantVol,
            capFloor,
            refRateCode,
            ycName,
            PREMIUM,
            logNorm,
            &prem_new);

        /* Checks if relative difference in premium is lower than 1.0e-06: if it is, stop */
        if (fabs(prem_new / premium - 1.0) < 1.0e-06)
        {
            *impvol = vol_guess;
            return NULL;
        }

        /* Set the static volatility to a shifted volatility guess */
        _capfloorconstantvol = vol_guess + vol_shift;

        /* Gets a new point with a slighly shifted vol */
        err = swp_f_CapFloor_SwapDP(
            swapdp,
            strike,
            ConstantVol,
            capFloor,
            refRateCode,
            ycName,
            PREMIUM,
            logNorm,
            &prem_shift);

        /* Computes the derivative of price wrt volatility */
        deriv = (prem_shift - prem_new) / vol_shift;

        if (deriv != 0)
        {
            /* Gets the new point via Newton method */
            tmp_vol = vol_guess - ((prem_new - premium) / deriv);

            /* Make sure we do not go too far */
            if (tmp_vol > (vol_guess * 10))
                vol_guess = vol_guess * 2.0;
            else if (tmp_vol < (vol_guess / 10))
                vol_guess = vol_guess / 2.0;
            else
                vol_guess = tmp_vol;
        }
        else
        {
            /* We are on a point where the deriv is flat: we try to get out of it */
            if (prem_new < premium)
                vol_guess = vol_guess * 2.0;
            else
                vol_guess = vol_guess / 1.5;
        }

    } /* END for(i=0;i<20;i++) loop on */

    /* Return an error message: no convergence so far */
    return serror("Newton iterations failed  in swp_f_CapFloorImpVol_SwapDP");

} /* END swp_f_CapFloorImpVol_SwapDP(...) */

/* ============================================================================== */
