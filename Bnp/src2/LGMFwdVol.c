#include "math.h"
#include "srt_h_all.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"

/*	Compute prices of options on forward swaps in LGM and AUTOCAL	*/

/*	1.-	LGM	*/
Err LGMFwdSwaptionLGM(
    long    exp,
    long    start,
    long    nfp_or_end,
    char*   cpdStr,
    char*   basisStr,
    char*   recPayStr,
    double  strike,
    char*   refRateCodeStr,
    char*   sUndPtrName,
    double* price)
{
    Err             err = NULL;
    SrtUndPtr       sUndPtr;
    SwapDP          sdp;
    SrtReceiverType rec_pay;
    SrtMdlDim       mdl_dim;
    Date            today;
    GenSwapLeg      fixedleg;
    GenSwapLeg      spreadleg;
    GenSwapLeg      bigleg;
    TermStruct*     ts;
    String          yc_name;

    /*	Gets the sUndPtrerlying through its name */
    if (sUndPtrName)
    {
        sUndPtr = lookup_und(sUndPtrName);
    }
    else
    {
        return serror("Heinous error: No sUndPtrerlying passed to LGMFwdSwaptionLGM");
    }

    /*	Checks the sUndPtrerlying type */
    if (sUndPtr)
    {
        if (!ISUNDTYPE(sUndPtr, INTEREST_RATE_UND) ||
            ((SrtIrDesc*)(sUndPtr->spec_desc))->mdl_type != LGM)
        {
            return serror("sUndPtrerlying must be of type Interest Rate, Model must be LGM");
        }
    }
    else
    {
        return serror("Could not find sUndPtrerlying in market list");
    }

    /*	Gets the Number of Factors in the Model */
    if (err = get_underlying_mdldim(sUndPtr, &mdl_dim))
    {
        return err;
    }

    /*	Gets the curve	*/
    yc_name = get_ycname_from_irund(sUndPtr);

    /*	Gets today from the underlying */
    today = (Date)get_clcndate_from_curve(lookup_curve(yc_name));

    /*	Populates the SwapDp structire with Spot Lag, Start Date,...*/
    sdp.spot_lag = get_spotlag_from_underlying(sUndPtr);

    if (err = swp_f_initSwapDP(start, nfp_or_end, cpdStr, basisStr, &sdp))
    {
        return err;
    }

    /*	Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr, &rec_pay))
    {
        return err;
    }

    /*	Make the fixed leg (coupon of strike) with initial and final */
    if (err = swp_f_make_FixedAndNotionalsLeg(&sdp, strike, 1.0, 1.0, today, &fixedleg))
    {
        return err;
    }

    /*	For the moment: assume receive fixed leg */
    fixedleg.rec_pay = SRT_RECEIVER;

    /*	Make the spread leg: dates, times, cvg, spreads cash vs Livor...*/
    if (err = swp_f_make_SpreadLeg(sdp.start, sdp.end, today, refRateCodeStr, &spreadleg))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        return err;
    }

    /*	For the moment: assume pay floating */
    spreadleg.rec_pay = SRT_PAYER;

    /*	Merge fixed and floating in one single big leg */
    if (err = swp_f_merge_SwapLegs(&fixedleg, &spreadleg, &bigleg))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        return err;
    }
    bigleg.sdp.spot_lag = sdp.spot_lag;

    /*	Make the times from today to the swap dates */
    if (err = time_list(bigleg.pay_date, bigleg.today, &(bigleg.time)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    /*	Make the maturities from today to the fixing dates */
    if (err = maturity_list(bigleg.fixing_date, bigleg.today, &(bigleg.mat)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    /*	Compute the discount factors for each cash flow date in the big list */
    if (err = df_list(bigleg.pay_date, yc_name, &(bigleg.df)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    /*	Gets the Term Structure from the underlying */
    if (err = get_underlying_ts(sUndPtr, &ts))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    *price = srt_f_lgm_coupon_bond_option(
        bigleg.leg_length - 1,
        -bigleg.payment.d[0],
        ts,
        (exp - today) * YEARS_IN_DAY,
        bigleg.payment.d,
        bigleg.time.d,
        bigleg.df.d,
        rec_pay,
        mdl_dim);

    swp_f_freein_GenSwapLeg(&fixedleg);
    swp_f_freein_GenSwapLeg(&spreadleg);
    swp_f_freein_GenSwapLeg(&bigleg);
    return NULL;
}

/*	2.-	Autocal */
Err LGMFwdSwaptionATC(
    long    exp,
    long    start,
    long    nfp_or_end,
    char*   cpdStr,
    char*   basisStr,
    char*   recPayStr,
    double  strike,
    char*   refRateCodeStr,
    char*   yc_name,
    LGM_TS* tsPtr,
    double* price)
{
    Err             err = NULL;
    SwapDP          sdp;
    SrtReceiverType rec_pay;
    Date            today;
    GenSwapLeg      fixedleg;
    GenSwapLeg      spreadleg;
    GenSwapLeg      bigleg;
    double          ystar;
    int             i;
    int             n;
    double*         a;
    double          dst;
    double          sqz;
    double*         g;
    double          gst;

    /*	Gets today from the underlying */
    if (!lookup_curve(yc_name))
    {
        return serror("Cannot find yeild curve %s", yc_name);
    }
    today = (Date)get_clcndate_from_curve(lookup_curve(yc_name));

    /*	Populates the SwapDp structire with Spot Lag, Start Date,...*/
    sdp.spot_lag = get_ccyparam_from_yldcrv(lookup_curve(yc_name))->spot_lag;

    if (err = swp_f_initSwapDP(start, nfp_or_end, cpdStr, basisStr, &sdp))
    {
        return err;
    }

    /*	Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr, &rec_pay))
    {
        return err;
    }

    /*	Make the fixed leg (coupon of strike) with initial and final */
    if (err = swp_f_make_FixedAndNotionalsLeg(&sdp, strike, 0.0, 0.0, today, &fixedleg))
    {
        return err;
    }

    /*	For the moment: assume receive fixed leg */
    fixedleg.rec_pay = SRT_RECEIVER;

    /*	Make the spread leg: dates, times, cvg, spreads cash vs Livor...*/
    if (err = swp_f_make_SpreadLeg(sdp.start, sdp.end, today, refRateCodeStr, &spreadleg))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        return err;
    }

    /*	For the moment: assume pay floating */
    spreadleg.rec_pay = SRT_PAYER;

    /*	Merge fixed and floating in one single big leg */
    if (err = swp_f_merge_SwapLegs(&fixedleg, &spreadleg, &bigleg))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        return err;
    }
    bigleg.sdp.spot_lag = sdp.spot_lag;

    /*	Make the times from today to the swap dates */
    if (err = time_list(bigleg.pay_date, bigleg.today, &(bigleg.time)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    /*	Make the maturities from today to the fixing dates */
    if (err = maturity_list(bigleg.fixing_date, bigleg.today, &(bigleg.mat)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    /*	Compute the discount factors for each cash flow date in the big list */
    if (err = df_list(bigleg.pay_date, yc_name, &(bigleg.df)))
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }

    n = bigleg.leg_length - 1;
    a = bigleg.payment.d + 1;
    for (i = 0; i < n; i++)
    {
        a[i] *= bigleg.df.d[i + 1];
    }
    a[n - 1] += bigleg.df.d[n];
    dst = bigleg.df.d[0];

    sqz = sqrt(LGMZetaFromTS(exp, today, tsPtr));

    g = (double*)calloc(n, sizeof(double));
    if (!g)
    {
        swp_f_freein_GenSwapLeg(&fixedleg);
        swp_f_freein_GenSwapLeg(&spreadleg);
        swp_f_freein_GenSwapLeg(&bigleg);
        return err;
    }
    for (i = 0; i < n; i++)
    {
        g[i] = LGMGFromTS(bigleg.pay_date.date[i + 1], tsPtr);
    }
    gst = LGMGFromTS(start, tsPtr);

    *price = LGMRecVal(&ystar, n, a - 1, dst, sqz, g - 1, gst);

    if (rec_pay == SRT_PAYER)
    {
        for (i = 0; i < n; i++)
        {
            *price -= a[i];
        }
        *price += dst;
    }

    swp_f_freein_GenSwapLeg(&fixedleg);
    swp_f_freein_GenSwapLeg(&spreadleg);
    swp_f_freein_GenSwapLeg(&bigleg);
    free(g);
    return NULL;
}