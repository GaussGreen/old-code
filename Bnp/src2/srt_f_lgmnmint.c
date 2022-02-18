/*------------------------------------------------------------------------*
 *
 * Module: srt_f_lgmnmint.c
 * Author: E.Auld
 * Jan 96
 * This module prices an all or nothing style rangefloater in LGM
 * by performing a double numerical integral.
 *
 * band is set at t1 around rate at t1;
 * libor of t2 plus a margin is payed at t3 if within the band.
 *
 * How this works: in order to use the numerical integration function
 * sm_qsimp, you need lots of global variables.  This method isn't actually
 * used (it is too slow); but that's why the code is structured that way
 *
 *
 *------------------------------------------------------------------------*/

#include "math.h"
#include "srt_h_all.h"
#include "utallhdr.h"

/* this was 32 points but ten have weight < 0.0000000001, these have been
   removed */
static double gauher_out[] = {
    6.0889643091,  0.0000000023, 5.4500332736,  0.0000000888, 4.8326046132,  0.0000020596,
    4.2320211100,  0.0000305598, 3.6447812499,  0.0003025570, 3.0681351690,  0.0020620511,
    2.4998404152,  0.0099034612, 1.9380049059,  0.0341098475, 1.3809801993,  0.0853448082,
    0.8272849038,  0.1565389938, 0.2755464192,  0.2117055699, -0.2755464192, 0.2117055699,
    -0.8272849038, 0.1565389938, -1.3809801993, 0.0853448082, -1.9380049059, 0.0341098475,
    -2.4998404152, 0.0099034612, -3.0681351690, 0.0020620511, -3.6447812499, 0.0003025570,
    -4.2320211100, 0.0000305598, -4.8326046132, 0.0000020596, -5.4500332736, 0.0000000888,
    -6.0889643091, 0.0000000023};

#define GAUHER_MAX 22

#define LGM_SIMP_TOL (1e-8)

#define LGMCMS_TODO 200 /* num of evaluations of inner integral */

/**************************************************************
  global variables for numerical integration
***************************************************************/

static SrtSample   NI_sam[3];
static SrtIRMTmInf NI_tminf[3];
static double      NI_strike;
static SrtRtFnc    NI_fwd_rtf, NI_lvl_rtf;
static int         NI_index = 0;
static SrtErr      global_err;

/*
 * more global variables for inner integral
 */
static SrtRtFnc CMS_RTF;
static SrtRtFnc CMSDF_RTF;
static double   CMS_STRIKE, CMS_STRIKE2, CMS_MARGIN;

/*
 * more global variables for outer integral
 */
static SrtRtFnc RFLTR_RT1_RTF;
static SrtRtFnc RFLTR_DF1_RTF;
static double   RFLTR_STRIKE;
static double   RFLTR_STRIKE2;

/*------------------------------------------------------------------------*
 *
 * stuff for inverting the relationship between r and a given rate
 * i.e. given a certain swap rate, what value of r (in LGM) corresponds
 * to that rate.  Of course there is not always a solution...
 * We need this (for example) to find limits of integration.
 *------------------------------------------------------------------------*/

/*
 * find the value of r such that rt evaluates to val
 * and put that r in sam. Assumes that other values of sam ---
 * fwd rate and phi --- have already been set
 */
static SrtRtFnc   RTSAFE_RT;
static SrtSample* RTSAFE_SAM;
static double     RTSAFE_VAL;

static Err rtsafe_func(double x, double* fx, double* dfx)
{
    double short_rate, tmp, tmp2, dx = 1.0e-6;

    sam_get(RTSAFE_SAM[0], 0, STATEVAR) = x;
    /* sad fix */
    sam_get(RTSAFE_SAM[0], 0, SHORT_RATE) = x;

    short_rate = x + sam_get(NI_tminf[NI_index + 1].fwd_sam, 0, F_0_t);
    /*  short_rate = x + sam_get(RTSAFE_SAM[0],0,F_0_t); */
    srt_f_rtevalyp(RTSAFE_RT, RTSAFE_SAM, 0, ONE_FAC, LGM);
    srt_f_rtevl(RTSAFE_RT, short_rate, &tmp);
    *fx = tmp - RTSAFE_VAL;

    x += dx;
    sam_get(RTSAFE_SAM[0], 0, STATEVAR) = x;
    /* sad fix */
    sam_get(RTSAFE_SAM[0], 0, SHORT_RATE) = x;

    short_rate = x + sam_get(NI_tminf[NI_index + 1].fwd_sam, 0, F_0_t);
    /*  short_rate = x + sam_get(RTSAFE_SAM[0],0,F_0_t); */
    srt_f_rtevalyp(RTSAFE_RT, RTSAFE_SAM, 0, ONE_FAC, LGM);
    srt_f_rtevl(RTSAFE_RT, short_rate, &tmp2);

    *dfx = (tmp2 - tmp) / dx;

    return NULL;
}

static SrtErr srt_f_lgmfindr(SrtRtFnc rt, SrtSample* sam, double val, double* answer)
{
    double tmp = 0.0;
    SrtErr err;

    RTSAFE_RT                    = rt;
    RTSAFE_SAM                   = sam;
    RTSAFE_VAL                   = val;
    err                          = rtsafe(rtsafe_func, -40.0, 40.0, 1.0e-7, 10, &tmp);
    sam_get(sam[0], 0, STATEVAR) = *answer = tmp;
    return err;
}

/*------------------------------------------------------------------------*
 * get stuff off the TermStructure
 *------------------------------------------------------------------------*/

static void lgm_populate_tminf(
    TermStruct* ts, double prev_time, double cur_time, SrtIRMTmInf* tminf)
{
    double prev_G, prev_H;
    tminf->ev.onef.F   = F_func(cur_time, ts);
    tminf->ev.onef.Psi = Psi_func(cur_time, ts);
    G_H_func(cur_time, ts, &tminf->rf.onef.G, &tminf->rf.onef.H);
    sam_get(tminf->fwd_sam, 0, PHI) = tminf->rf.onef.G * tminf->ev.onef.F * tminf->ev.onef.F;

    G_H_func(prev_time, ts, &prev_G, &prev_H);
    tminf->rf.onef.stdev_x = tminf->ev.onef.F * sqrt(tminf->rf.onef.G - prev_G);
}

/*--------------------------------------------------------------
 * price cms type double digital option in LGM;
 * (range floater once strike has been set);
 * payoff is (cms(t) + margin)*DF(paydt)*if(cms_in_band,1,0)
 --------------------------------------------------------------*/
/*
 * dr is supposed to be Normal (0,1)
 */

static double lgmcms_integrand(double dr)
{
    double cms, df, short_rate, payoff;

    sam_get(NI_sam[NI_index + 1], 0, STATEVAR) =
        sam_get(NI_sam[NI_index], 0, STATEVAR) * NI_tminf[NI_index + 1].ev.onef.F /
            NI_tminf[NI_index].ev.onef.F +
        NI_tminf[NI_index + 1].ev.onef.F * NI_tminf[NI_index].rf.onef.G *
            (NI_tminf[NI_index + 1].ev.onef.Psi - NI_tminf[NI_index].ev.onef.Psi) +
        dr * NI_tminf[NI_index + 1].rf.onef.stdev_x;

    sam_get(NI_sam[NI_index + 1], 0, SHORT_RATE) = short_rate =
        sam_get(NI_sam[NI_index + 1], 0, STATEVAR) +
        sam_get(NI_tminf[NI_index + 1].fwd_sam, 0, F_0_t);

    /* sad fix */
    sam_get(NI_sam[NI_index + 1], 0, SHORT_RATE) = sam_get(NI_sam[NI_index + 1], 0, STATEVAR);

    srt_f_rtevalyp(CMS_RTF, &NI_sam[NI_index + 1], 0, ONE_FAC, LGM);
    srt_f_rtevalyp(CMSDF_RTF, &NI_sam[NI_index + 1], 0, ONE_FAC, LGM);
    srt_f_rtevl(CMS_RTF, short_rate, &cms);
    srt_f_rtevl(CMSDF_RTF, short_rate, &df);

    payoff =
        df * gauss(dr) * (cms + CMS_MARGIN) * ((cms < CMS_STRIKE && cms > CMS_STRIKE2) ? 1.0 : 0.0);

    return payoff;
}

/*--------------------------------------------------------------
 * price cms type double digital option in LGM;
 * (range floater once strike has been set);
 * payoff is (cms(t2) + margin)*DF(paydt)*if(cms_in_band,1,0)
 * where band = rate(t1)+- strike
 --------------------------------------------------------------*/

static double lgmrfltr_integrand(double dr)
{
    int    i;
    double df, rate;
    double lower = -10.0, upper = 10.0;
    double (*f)(double);
    double short_rate, payoff;
    double mu;
    double x, dx;
    SrtErr err;

    sam_get(NI_sam[NI_index + 1], 0, STATEVAR) =
        sam_get(NI_sam[NI_index], 0, STATEVAR) * NI_tminf[NI_index + 1].ev.onef.F /
            NI_tminf[NI_index].ev.onef.F +
        NI_tminf[NI_index + 1].ev.onef.F * NI_tminf[NI_index].rf.onef.G *
            (NI_tminf[NI_index + 1].ev.onef.Psi - NI_tminf[NI_index].ev.onef.Psi) +
        dr * NI_tminf[NI_index + 1].rf.onef.stdev_x;

    sam_get(NI_sam[NI_index + 1], 0, SHORT_RATE) = short_rate =
        sam_get(NI_sam[NI_index + 1], 0, STATEVAR) +
        sam_get(NI_tminf[NI_index + 1].fwd_sam, 0, F_0_t);

    /* sad fix */
    sam_get(NI_sam[NI_index + 1], 0, SHORT_RATE) = sam_get(NI_sam[NI_index + 1], 0, STATEVAR);

    srt_f_rtevalyp(RFLTR_RT1_RTF, &NI_sam[NI_index + 1], 0, ONE_FAC, LGM);
    srt_f_rtevalyp(RFLTR_DF1_RTF, &NI_sam[NI_index + 1], 0, ONE_FAC, LGM);
    srt_f_rtevl(RFLTR_RT1_RTF, short_rate, &rate);
    srt_f_rtevl(RFLTR_DF1_RTF, short_rate, &df);

    CMS_STRIKE  = RFLTR_STRIKE + rate;
    CMS_STRIKE2 = -RFLTR_STRIKE2 + rate;

    NI_index++;
    /* find limits of integration; note that answers here are returned
     in terms of STATEVAR, so we have to transform to get variance 1 */

    err = srt_f_lgmfindr(CMS_RTF, &NI_sam[NI_index + 1], CMS_STRIKE, &upper);
    if (err)
        global_err = err;
    err = srt_f_lgmfindr(CMS_RTF, &NI_sam[NI_index + 1], CMS_STRIKE2, &lower);
    if (err)
        global_err = err;

    mu = sam_get(NI_sam[NI_index], 0, STATEVAR) * NI_tminf[NI_index + 1].ev.onef.F /
             NI_tminf[NI_index].ev.onef.F +
         NI_tminf[NI_index + 1].ev.onef.F * NI_tminf[NI_index].rf.onef.G *
             (NI_tminf[NI_index + 1].ev.onef.Psi - NI_tminf[NI_index].ev.onef.Psi);

    upper = (upper - mu) / NI_tminf[NI_index + 1].rf.onef.stdev_x;
    lower = (lower - mu) / NI_tminf[NI_index + 1].rf.onef.stdev_x;

    if (fabs(upper - lower) < 1.0e-6)
        payoff = 0.0;
    else
    {
        dx     = 1.0;
        f      = lgmcms_integrand;
        payoff = 0.0;
        dx     = (upper - lower) / ((double)LGMCMS_TODO);
        for (i = 0, x = lower + dx / 2; i < LGMCMS_TODO; i++, x += dx)
        {
            payoff += lgmcms_integrand(x);
        }

        payoff *= dx * df;
    }
    NI_index--;
    return payoff;
}

/*
 *  payoff for a call is:
 *  MAX(rate2 - (rate1+strike),0.0) +
 *  margin * (rate2 - (rate1+strike) > 0.0)
 *  payed at paydt.  Rate1 and rate2 are fixed at fix1 and fix2.
 */

/*
         err = srt_f_lgmrfltr(start,end,freq,bas,start2,end2,fix1,fix2,
          paydt,strike,strike2,margin,und,&answer);
*/

SrtErr srt_f_lgmrfltr(
    Date      start1,
    Date      enfp1,
    char*     freq,
    char*     bas,
    Date      start2,
    Date      enfp2,
    Date      fix1,
    Date      fix2,
    Date      paydt,
    double    strike,
    double    strike2,
    double    margin,
    SrtUndPtr und,
    double    rate1,
    double    rate2,
    double*   answer)
{
    SrtCurvePtr crv;
    SrtRtFnc    r1_rtf, r2_rtf, df_f1_f2_rtf, df_f2_pdt_rtf;
    SrtErr      err;
    char*       ycname;
    TermStruct* ts;
    Date        today;
    double (*f)(double), t[3], mu, lower, upper, price;
    int i;

    today  = get_today_from_underlying(und);
    ycname = get_ycname_from_irund(und);

    /*
     * deal with stupid cases
     */
    if (fix2 <= fix1 || start1 < fix1 || start2 < fix2)
        return serror("LGM rfltr bad fixing dates");
    if (paydt < fix2 || paydt < today)
        return serror("LGM rfltr bad pay date");
    if (fix2 <= today)
    {
        /*    und = lookup_und(ycname);
            price = ((rate2 > rate1 - strike2) && (rate2 < rate1 + strike) ?
                margin + rate2 : 0);
            price *= df_mkt(today,paydt,und);
            *answer = price;
            return NULL;
        */
        /* shouldn't be using this program any more */
        return serror("Range Floater both fixings already happened.");
    }

    /*
     * rates we will need
     */
    err = srt_f_rtmk("SWAP", start1, NULL, enfp1, NULL, freq, bas, &r1_rtf);
    err = srt_f_rtmk("SWAP", start2, NULL, enfp2, NULL, freq, bas, &r2_rtf);
    err = srt_f_rtmk("DF", fix1, NULL, fix2, NULL, freq, bas, &df_f1_f2_rtf);
    err = srt_f_rtmk("DF", fix2, NULL, paydt, NULL, freq, bas, &df_f2_pdt_rtf);
    err = srt_f_rtinityp(r1_rtf, und, (Ddate)fix1);
    err = srt_f_rtinityp(r2_rtf, und, (Ddate)fix2);
    err = srt_f_rtinityp(df_f1_f2_rtf, und, (Ddate)fix1);
    err = srt_f_rtinityp(df_f2_pdt_rtf, und, (Ddate)fix2);

    /*
     * get ts
     */

    err = get_underlying_ts(und, &ts);
    if (err)
        return err;
    crv = lookup_curve(ycname);

    /*
     * set global variables
     */
    RFLTR_RT1_RTF = r1_rtf;
    RFLTR_DF1_RTF = df_f1_f2_rtf;
    RFLTR_STRIKE  = strike;
    RFLTR_STRIKE2 = strike2;
    CMS_RTF       = r2_rtf;
    CMSDF_RTF     = df_f2_pdt_rtf;
    CMS_MARGIN    = margin;
    NI_index      = 0;

    /*
     * deal with case where strike has already been set
     */
    if (fix1 <= today)
    {
        t[0] = 0.0;
        t[1] = YEARS_IN_DAY * (fix2 - today);
        lgm_populate_tminf(ts, t[0], t[0], &NI_tminf[0]);
        for (i = 1; i < 2; i++)
        {
            lgm_populate_tminf(ts, t[i - 1], t[i], &NI_tminf[i]);
        }
        sam_get(NI_tminf[0].fwd_sam, 0, F_0_t) = swp_f_zr(today, fix2, ycname);
        sam_get(NI_tminf[1].fwd_sam, 0, F_0_t) = swp_f_zr(fix2, fix2 + 30.0, ycname);
        sam_get(NI_sam[0], 0, PHI)             = sam_get(NI_tminf[0].fwd_sam, 0, PHI);
        sam_get(NI_sam[1], 0, PHI)             = sam_get(NI_tminf[1].fwd_sam, 0, PHI);
        sam_get(NI_sam[0], 0, STATEVAR)        = 0.0;

        CMS_STRIKE  = RFLTR_STRIKE + rate1;
        CMS_STRIKE2 = -RFLTR_STRIKE2 + rate1;

        /* find limits of integration; note that answers here are returned
         in terms of STATEVAR, so we have to transform to get variance 1 */

        srt_f_lgmfindr(CMS_RTF, &NI_sam[NI_index + 1], CMS_STRIKE, &upper);
        srt_f_lgmfindr(CMS_RTF, &NI_sam[NI_index + 1], CMS_STRIKE2, &lower);

        mu = sam_get(NI_sam[NI_index], 0, STATEVAR) * NI_tminf[NI_index + 1].ev.onef.F /
                 NI_tminf[NI_index].ev.onef.F +
             NI_tminf[NI_index + 1].ev.onef.F * NI_tminf[NI_index].rf.onef.G *
                 (NI_tminf[NI_index + 1].ev.onef.Psi - NI_tminf[NI_index].ev.onef.Psi);

        upper = (upper - mu) / NI_tminf[NI_index + 1].rf.onef.stdev_x;
        lower = (lower - mu) / NI_tminf[NI_index + 1].rf.onef.stdev_x;

        if (fabs(upper - lower) < 1.0e-6)
            price = 0.0;
        else
        {
            f     = lgmcms_integrand;
            price = sm_qsimp(f, lower, upper, LGM_SIMP_TOL);
        }
        price *= swp_f_df(today, fix2, ycname);
    }
    else
    {
        /*
         * get drift, variance of r at interesting times;
         * and get other important bits of information off the term structure
         */
        t[0] = 0.0;
        t[1] = YEARS_IN_DAY * (fix1 - today);
        t[2] = YEARS_IN_DAY * (fix2 - today);

        lgm_populate_tminf(ts, t[0], t[0], &NI_tminf[0]);
        for (i = 1; i < 3; i++)
        {
            lgm_populate_tminf(ts, t[i - 1], t[i], &NI_tminf[i]);
        }
        sam_get(NI_tminf[0].fwd_sam, 0, F_0_t) = swp_f_zr(today, fix1, ycname);
        sam_get(NI_tminf[1].fwd_sam, 0, F_0_t) = swp_f_zr(fix1, fix2, ycname);
        sam_get(NI_tminf[2].fwd_sam, 0, F_0_t) = swp_f_zr(fix2, 2 * fix2 - fix1, ycname);
        sam_get(NI_sam[0], 0, PHI)             = sam_get(NI_tminf[0].fwd_sam, 0, PHI);
        sam_get(NI_sam[1], 0, PHI)             = sam_get(NI_tminf[1].fwd_sam, 0, PHI);
        sam_get(NI_sam[2], 0, PHI)             = sam_get(NI_tminf[2].fwd_sam, 0, PHI);
        sam_get(NI_sam[0], 0, STATEVAR)        = 0.0;

        /*
         integrate by Gauss Hermite quadrature
         we can't do this for the other integral because the limits of integration
         aren't +- infinity (although we could find some other fancy thing)
         */
        price      = 0;
        global_err = NULL;
        for (i = 0; i < GAUHER_MAX; i++)
        {
            price += lgmrfltr_integrand(gauher_out[2 * i]) * gauher_out[2 * i + 1];
            if (global_err)
                return global_err;
        }
        price *= swp_f_df(today, fix1, ycname);
    }
    *answer = price;

    srt_f_rtfre(r1_rtf);
    srt_f_rtfre(r2_rtf);
    srt_f_rtfre(df_f1_f2_rtf);
    srt_f_rtfre(df_f2_pdt_rtf);
    return err;
}
