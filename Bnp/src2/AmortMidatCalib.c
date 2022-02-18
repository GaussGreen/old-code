
#include "AmortMidatCalib.h"

#include "AmortSwapUtils.h"
#include "CPDCalib.h"
#include "DiagCalibDLM.h"
#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "srt_h_lgmprotos.h"
#include "srt_h_lgmtypes.h"
#include "swp_h_swap_pricing.h"

#define MAX_CPN 600
#define MAX_INST 600
#define NUM_HERMITE 6
#define CALPRESAMORTMIDAT 1.0e-06
#define NITER 10
#define BUMP_LAM_LM 5.0E-04
#define LAMBDA_REF 0.05

extern double x[NUM_HERMITE + 1], w[NUM_HERMITE + 1];

static double MAX_FACT = 0.001;

static int USE_JUMPS = 0;

static double CALPRESNOT = 1.0e-04;

/*Price of the sum of diagonal swaptions within LGM1F*/
double amortMidat_lgmsumsopval1F(
    int ncpn,         /*	Number of cash-flow dates, including
                                              start and end date */
    int    nex,       /*	Number of exercise dates */
    int    ex_cpn[],  /*	Index of the first cash-flow to be exercised */
    double df[],      /*	Df to cash flow dates */
    double cvg[],     /*	cvg from i-1 to i */
    double cpn_G[],   /*	G at cash-flow dates */
    double ex_zeta[], /*	Z at exercise dates */
    double ex_G[],    /*	G at exercise dates */
    double strike[])  /*	Strikes */
{
    int    i;
    double price;

    price = 0;

    for (i = 0; i < nex; ++i)
    {
        price += lgmsopval1F(
            ncpn - ex_cpn[i],
            df + ex_cpn[i],
            cvg + ex_cpn[i],
            cpn_G + ex_cpn[i],
            ex_zeta[i],
            ex_G[i],
            strike[i]);
    }

    return price;
}

/*	Precalculate exponential factors in the IV, and decide what variable to use
        in integration, and which one to use in y* solving  */
static void amort_static_lgm2Fcalcexpfact(
    int    ncpn,
    double cpn[],
    double cpn_G1[],
    double cpn_G2[],
    double ex_zeta1,
    double sqz1,
    double ex_zeta2,
    double sqz2,
    double ex_zeta12,
    double c1,
    double c2,
    double ex_G1,
    double ex_G2,
    double cstfact[],
    double n1fact[],
    double n2fact[],
    int*   highest) /*	Index of the variable to be used for solving
                                          - the one with the highest sensitivity */
{
    int    i;
    double s1 = 0.0, s2 = 0.0;
    double half_z1 = 0.5 * ex_zeta1, half_z2 = 0.5 * ex_zeta2;
    double c1sqz2 = c1 * sqz2, c2sqz2 = c2 * sqz2;

    for (i = 0; i < ncpn; i++)
    {
        cstfact[i] = half_z1 * (cpn_G1[i] - ex_G1) * (cpn_G1[i] - ex_G1);
        cstfact[i] += half_z2 * (cpn_G2[i] - ex_G2) * (cpn_G2[i] - ex_G2);
        cstfact[i] += ex_zeta12 * (cpn_G1[i] - ex_G1) * (cpn_G2[i] - ex_G2);

        n1fact[i] = (cpn_G1[i] - ex_G1) * sqz1 + (cpn_G2[i] - ex_G2) * c1sqz2;

        n2fact[i] = (cpn_G2[i] - ex_G2) * c2sqz2;

        s1 -= cpn[i] * n1fact[i] * exp(-cstfact[i]);
        s2 -= cpn[i] * n2fact[i] * exp(-cstfact[i]);
    }

    if (fabs(s1) > fabs(s2))
    {
        *highest = 1;
    }
    else
    {
        *highest = 2;
    }
}

/*	In order to help static_lgmystar: PV of coupons when reconstruction is
                df (t, Ti) = df (0, t, Ti) * exp (a[i] + b[i] * Gaussian)
                and Gaussian is y */
static double amort_static_lgmiv(
    double y,     /*	The y */
    int    ncpn,  /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    double a[],   /*	The ai */
    double b[])   /*	The bi */
{
    int    i;
    double val;

    val = 0.0;
    for (i = 0; i < ncpn; i++)
    {
        val += cpn[i] * exp(a[i] + b[i] * y);
    }

    return val;
}

#define YPREC 1.0e-08
/*	Solve y / static_lgmiv (y) = 0 with Newton */
static double amort_static_solvelgmy(
    int    ncpn,  /*	Num coupons */
    double cpn[], /*	Discounted coupons */
    double a[],   /*	The ai */
    double b[],   /*	The bi */
    int*   dir)     /*	Output, whether iv is increasing in y
                                          1: decreasing, -1: increasing */
{
    double y1, y2, iv1, iv2, temp;
    int    it;

    y1 = 0.0;

    iv1 = amort_static_lgmiv(y1, ncpn, cpn, a, b);

    y2 = 1.0e-02;

    iv2 = amort_static_lgmiv(y2, ncpn, cpn, a, b);

    if (iv2 < iv1)
    {
        *dir = 1;
    }
    else
    {
        *dir = -1;
    }

    if (fabs(iv1) < YPREC || fabs(iv2) < YPREC)
    {
        if (fabs(iv1) < fabs(iv2))
        {
            return y1;
        }
        else
        {
            return y2;
        }
    }

    for (it = 1; it < 25; it++)
    {
        if (fabs(y2 - y1) < YPREC || fabs(iv2 - iv1) < YPREC)
            break;

        temp = y2;
        y2 -= iv2 * (y2 - y1) / (iv2 - iv1);
        y2  = DMIN(20.0, DMAX(-20.0, y2));
        y1  = temp;
        iv1 = iv2;

        iv2 = amort_static_lgmiv(y2, ncpn, cpn, a, b);

        if (fabs(iv2) < YPREC)
            break;
    }

    if (fabs(iv2) > YPREC)
    {
        /*	We should return an error here */
    }

    return y2;
}

/*	Solve y* in the model: 2F case */
static double amort_static_lgmystar2F(
    int    ncpn,  /*	Num coupons */
    double cpn[], /*	Discounted coupons */
                  /*	As output from static_lgm2Fcalcexpfact */
    double cstfact[],
    double n1fact[],
    double n2fact[],
    double highest,
    /*	Condition on Intergral variable = n */
    double n,
    int*   dir) /*	Output, whether iv is increasing in y
                                      1: decreasing, -1: increasing */
{
    int i;

    double a[MAX_CPN], b[MAX_CPN];

    for (i = 0; i < ncpn; i++)
    {
        a[i] = -cstfact[i];

        if (highest == 1)
        {
            a[i] -= n2fact[i] * n;
            b[i] = -n1fact[i];
        }
        else
        {
            a[i] -= n1fact[i] * n;
            b[i] = -n2fact[i];
        }
    }

    /*	Solver */
    return amort_static_solvelgmy(ncpn, cpn, a, b, dir);
}

double amort_lgmopval2F(
    int    ncpn,      /*	Number of cash-flows */
    double cpn[],     /*	Discounted Cash-Flows */
    double cpn_G1[],  /*	G1 at cash-flow dates */
    double cpn_G2[],  /*	G2 at cash-flow dates */
    double ex_zeta1,  /*	Z1 at exercise date */
    double ex_zeta2,  /*	Z2 at exercise date */
    double ex_zeta12, /*	Z12 at exercise date */
    double ex_G1,     /*	G1 at exercise date */
    double ex_G2)     /*	G2 at exercise date */
{
    double cstfact[MAX_CPN], n1fact[MAX_CPN], n2fact[MAX_CPN];
    double sqz1 = sqrt(ex_zeta1), sqz2 = sqrt(ex_zeta2);
    double c1 = ex_zeta12 / sqz1 / sqz2, c1sqz2 = c1 * sqz2, c2 = sqrt(1.0 - c1 * c1),
           c2sqz2 = c2 * sqz2;
    double *X = &(x[1]), *W = &(w[1]);
    double  ystar;
    int     i, j;
    int     highest, dir;
    double  temp, val = 0.0;

    if ((W[0] == 0) && (X[0] == 0))
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    /*	Check for intrinsic */
    if (ex_zeta1 < 1.0e-10 && ex_zeta2 < 1.0e-10)
    {
        for (i = 0; i < ncpn; i++)
        {
            val += cpn[i];
        }

        return max(0.0, val);
    }

    /*	Precalculate exponential factors in the IV */
    amort_static_lgm2Fcalcexpfact(
        ncpn,
        cpn,
        cpn_G1,
        cpn_G2,
        ex_zeta1,
        sqz1,
        ex_zeta2,
        sqz2,
        ex_zeta12,
        c1,
        c2,
        ex_G1,
        ex_G2,
        cstfact,
        n1fact,
        n2fact,
        &highest);

    /*	Do Hermite integration */
    for (i = 0; i < NUM_HERMITE; i++)
    {
        /*	Disregard points of little weight */
        //		if (fabs (X[i]) > 5.0) continue;

        temp = 0.0;
        for (j = 0; j < ncpn; j++)
        {
            ystar = amort_static_lgmystar2F(
                ncpn,
                cpn,
                cstfact,
                n1fact,
                n2fact,
                highest,
                X[i] - (highest == 1 ? n2fact[j] : n1fact[j]),
                &dir);

            temp +=
                cpn[j] * static_lgmsafenorm(dir * (ystar + (highest == 1 ? n1fact[j] : n2fact[j])));
        }

        val += W[i] * temp;
    }

    return val;
}

/*	Useful functions designed to get term structures */

Err Get_LGM_TermStructureOneOrTwoFact(
    char*    underlying,
    double** sigma_date,
    double** sigma,
    long*    sigma_n,
    double** tau_date,
    double** tau,
    long*    tau_n,
    double*  alpha,
    double*  gamma,
    double*  rho)

{
    SrtUndPtr   und;
    TermStruct* ts = NULL;
    SrtMdlDim   mdl_dim;
    double *    beta = NULL, *beta2 = NULL, *vovol = NULL, *rhots = NULL, *meanvol = NULL;
    double*     sigma2  = NULL;
    double*     rhovect = NULL;
    double*     tau2    = NULL;
    long        today;
    int         i;
    Err         err;

    und = lookup_und(underlying);
    if (!und)
    {
        return serror("Couldn't find underlying named %s", underlying);
    }

    if (get_underlying_type(und) != INTEREST_RATE_UND)
    {
        return serror("Underlying %s is not of type IR", underlying);
    }

    if (get_mdltype_from_irund(und) != LGM)
    {
        return serror("Underlying %s is not of type LGM", underlying);
    }

    ts    = get_ts_from_irund(und);
    today = get_today_from_underlying(und);

    err = get_underlying_mdldim(und, &mdl_dim);

    if (mdl_dim == ONE_FAC)
    {
        err = srt_f_display_IRM_OneFac_TermStruct(
            ts, sigma_date, sigma, &beta, &vovol, &rhots, &meanvol, sigma_n, tau_date, tau, tau_n);
        *alpha = 1.0;
        *gamma = 0.0;
        *rho   = 1.0;
    }
    else
    {
        err = srt_f_display_IRM_TwoFac_TermStruct(
            ts,
            sigma_date,
            sigma,
            &beta,
            &sigma2,
            &beta2,
            &rhovect,
            sigma_n,
            tau_date,
            tau,
            &tau2,
            tau_n);
        *alpha = sigma2[0] / (*sigma)[0];
        *gamma = 1.0 / tau2[0] - 1.0 / (*tau)[0];
        *rho   = rhovect[0];
    }

    for (i = 0; i < *sigma_n; i++)
    {
        (*sigma_date)[i] = ((*sigma_date)[i] - today) / 365.0;
    }
    for (i = 0; i < *tau_n; i++)
    {
        (*tau_date)[i] = ((*tau_date)[i] - today) / 365.0;
    }

    if (beta)
        free(beta);
    if (vovol)
        free(vovol);
    if (rhots)
        free(rhots);
    if (meanvol)
        free(meanvol);

    if (mdl_dim == TWO_FAC)
    {
        free(beta2);
        free(sigma2);
        free(tau2);
        free(rhovect);
    }

    return err;
}

double amortMidat_lgmamortsopval1F(
    int ncpn,        /*	Number of cash-flow dates, including
                                             start and end date */
    double cpn[],    /*	Notional */
    double df[],     /*	Df to cash flow dates */
    double cvg[],    /*	cvg from i-1 to i */
    double cpn_G1[], /*	G1 at cash-flow dates */
    double ex_zeta1, /*	Z1 at exercise date */
    double ex_G1)    /*	G1 at exercise date */
{
    int           i;
    static double coupon[MAX_CPN];

    for (i = 0; i < ncpn; i++)
    {
        coupon[i] = cpn[i];
    }

    return lgmopval1F(ncpn, coupon, cpn_G1, ex_zeta1, ex_G1);
}

double amortMidat_lgmamortsopval2F(
    int ncpn,         /*	Number of cash-flow dates, including
                                              start and end date */
    double cpn[],     /*	Notional */
    double df[],      /*	Df to cash flow dates */
    double cvg[],     /*	cvg from i-1 to i */
    double cpn_G1[],  /*	G1 at cash-flow dates */
    double cpn_G2[],  /*	G2 at cash-flow dates */
    double ex_zeta1,  /*	Z1 at exercise date */
    double ex_zeta2,  /*	Z2 at exercise date */
    double ex_zeta12, /*	Z12 at exercise date */
    double ex_G1,     /*	G1 at exercise date */
    double ex_G2)     /*	G2 at exercise date */
//	double		strike)								/*	Strike
//*/
{
    int           i;
    static double coupon[MAX_CPN];

    for (i = 0; i < ncpn; i++)
    {
        coupon[i] = cpn[i];
    }

    return amort_lgmopval2F(
        ncpn, coupon, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12, ex_G1, ex_G2);
}

/*	Value of European Swap option within LGM 2F */
double amortMidat_lgmsopval2F(
    int ncpn,         /*	Number of cash-flow dates, including
                                              start and end date */
    double df[],      /*	Df to cash flow dates */
    double cvg[],     /*	cvg from i-1 to i */
    double cpn_G1[],  /*	G1 at cash-flow dates */
    double cpn_G2[],  /*	G2 at cash-flow dates */
    double ex_zeta1,  /*	Z1 at exercise date */
    double ex_zeta2,  /*	Z2 at exercise date */
    double ex_zeta12, /*	Z12 at exercise date */
    double ex_G1,     /*	G1 at exercise date */
    double ex_G2,     /*	G2 at exercise date */
    double strike)    /*	Strike */
{
    int           i;
    static double cpn[MAX_CPN];

    cpn[0] = -df[0];

    for (i = 1; i < ncpn; i++)
    {
        cpn[i] = df[i] * cvg[i] * strike;
    }

    cpn[ncpn - 1] += df[ncpn - 1];

    return amort_lgmopval2F(ncpn, cpn, cpn_G1, cpn_G2, ex_zeta1, ex_zeta2, ex_zeta12, ex_G1, ex_G2);
}

/*	Value of Sum of diagonal swaptions within LGM 2F */
double amortMidat_lgmsumsopval2F(
    int ncpn,                 /*	Number of cash-flow dates, including
                                                      start and end date */
    int    nex,               /*	Number of exercise dates */
    int    ex_cpn[],          /*	Index of the first cash-flow to be exercised */
    double df[],              /*	Df to cash flow dates */
    double cvg[],             /*	cvg from i-1 to i */
    double cpn_G1[],          /*	G1 at cash-flow dates */
    double cpn_G2[],          /*	G2 at cash-flow dates */
    double ex_zeta1[],        /*	Z1 at exercise dates */
    double ex_zeta2[],        /*	Z2 at exercise dates */
    double ex_zeta12[],       /*	Z12 at exercise dates */
    double ex_G1[],           /*	G1 at exercise dates */
    double ex_G2[],           /*	G2 at exercise dates */
    double swaption_strike[]) /*	Standard Swaption Coupons */
{
    int    i;
    double price;

    price = 0;

    for (i = 0; i < nex; ++i)
    {
        price += amortMidat_lgmsopval2F(
            ncpn - ex_cpn[i],
            df + ex_cpn[i],
            cvg + ex_cpn[i],
            cpn_G1 + ex_cpn[i],
            cpn_G2 + ex_cpn[i],
            ex_zeta1[i],
            ex_zeta2[i],
            ex_zeta12[i],
            ex_G1[i],
            ex_G2[i],
            swaption_strike[i]);
    }

    return price;
}

double amortMidat_lgmsumsopval2F_b(
    int ncpn,              /*	Number of cash-flow dates, including
                                                   start and end date */
    int    nex,            /*	Number of exercise dates */
    int    ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double df[],           /*	Df to cash flow dates */
    double cpn_G1[],       /*	G1 at cash-flow dates */
    double cpn_G2[],       /*	G2 at cash-flow dates */
    double ex_zeta1[],     /*	Z1 at exercise dates */
    double ex_zeta2[],     /*	Z2 at exercise dates */
    double ex_zeta12[],    /*	Z12 at exercise dates */
    double ex_G1[],        /*	G1 at exercise dates */
    double ex_G2[],        /*	G2 at exercise dates */
    double swaption_cpn[]) /*	Standard Swaption Coupons */
{
    int    i, j;
    double price;
    double tempo;

    price = 0;

    for (i = 0; i < nex; ++i)
    {
        j = ex_cpn[i];

        tempo           = swaption_cpn[j];
        swaption_cpn[j] = -df[j];

        price += amort_lgmopval2F(
            ncpn - j,
            swaption_cpn + j,
            cpn_G1 + j,
            cpn_G2 + j,
            ex_zeta1[i],
            ex_zeta2[i],
            ex_zeta12[i],
            ex_G1[i],
            ex_G2[i]);

        swaption_cpn[j] = tempo;
    }

    return price;
}

double amortMidat_lgmsumsopval1F_b(
    int ncpn,              /*	Number of cash-flow dates, including
                                                   start and end date */
    int    nex,            /*	Number of exercise dates */
    int    ex_cpn[],       /*	Index of the first cash-flow to be exercised */
    double df[],           /*	Df to cash flow dates */
    double cpn_G[],        /*	G at cash-flow dates */
    double ex_zeta[],      /*	Z at exercise dates */
    double ex_G[],         /*	G at exercise dates */
    double swaption_cpn[]) /*	Standard Swaption Coupons */
{
    int    i, j;
    double price;
    double tempo;

    price = 0;

    for (i = 0; i < nex; ++i)
    {
        j = ex_cpn[i];

        tempo           = swaption_cpn[j];
        swaption_cpn[j] = -df[j];

        price += lgmopval1F(ncpn - j, swaption_cpn + j, cpn_G, ex_zeta[i], ex_G[i]);

        swaption_cpn[j] = tempo;
    }

    return price;
}

/*	Value of European Cap within LGM 1F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one exercise date to the
   following */
double amortMidat_lgmcapval1F(
    int ncpn,          /*	Number of cash-flow dates, including
                                               start and end date */
    double cpn[],      /*	Notional */
    double df[],       /*	Df to cash flow dates */
    double cvg[],      /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    int    nex,        /*	Number of exercise dates */
    int    ex_cpn[],   /*	For each exercise date, first coupon
                                               to be exercised */
    int ex_ncpn[],     /*	For each exercise date, number of coupons
                                               to be exercised */
    double ex_zeta1[], /*	Z1 at exercise date */
    double ex_G1[])    /*	G1 at exercise date */

{
    double val;
    double temp;
    int    i;

    val = 0.0;

    for (i = 0; i < nex; i++)
    {
        temp           = cpn[ex_cpn[i]];
        cpn[ex_cpn[i]] = -df[ex_cpn[i]];
        if (i < nex - 1)
        {
            cpn[ex_cpn[i] + ex_ncpn[i] - 1] =
                cpn[ex_cpn[i] + ex_ncpn[i] - 1] + df[ex_cpn[i] + ex_ncpn[i] - 1];
        }

        val += amortMidat_lgmamortsopval1F(
            ex_ncpn[i],
            &(cpn[ex_cpn[i]]),
            &(df[ex_cpn[i]]),
            &(cvg[ex_cpn[i]]),
            &(cpn_G1[ex_cpn[i]]),
            ex_zeta1[i],
            ex_G1[i]);

        cpn[ex_cpn[i]] = temp;
        if (i < nex - 1)
        {
            cpn[ex_cpn[i] + ex_ncpn[i] - 1] =
                cpn[ex_cpn[i] + ex_ncpn[i] - 1] - df[ex_cpn[i] + ex_ncpn[i] - 1];
        }
    }

    return val;
}

/*	Value of European Cap within LGM 1F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one exercise date to the
   following */
double amortMidat_lgmcapval1F_b(
    int ncpn,             /*	Number of cash-flow dates, including
                                                  start and end date */
    double cpn[],         /*	Notional */
    double df[],          /*	Df to cash flow dates */
    double cvg[],         /*	cvg from i-1 to i */
    double cpn_G1[],      /*	G1 at cash-flow dates */
    int    nex,           /*	Number of exercise dates */
    double ex_sstrike[],  /*	Strikes for cap */
    double floatcoupon[], /*	Float Coupon*/
    int    ex_cpn[],      /*	For each exercise date, first coupon
                                                  to be exercised */
    int ex_ncpn[],        /*	For each exercise date, number of coupons
                                                  to be exercised */
    double ex_zeta1[],    /*	Z1 at exercise date */
    double ex_G1[])       /*	G1 at exercise date */

{
    double val;
    double coupon[MAX_CPN];
    int    i, j;

    val = 0.0;

    for (i = 0; i < nex; i++)
    {
        coupon[ex_cpn[i]] = -df[ex_cpn[i]];
        for (j = ex_cpn[i] + 1; j < ex_cpn[i] + ex_ncpn[i] - 1; ++j)
        {
            coupon[j] = floatcoupon[j];
        }
        coupon[ex_cpn[i] + ex_ncpn[i] - 1] = floatcoupon[j - 1] + df[ex_cpn[i] + ex_ncpn[i] - 1] +
                                             ex_sstrike[i] * df[ex_cpn[i] + ex_ncpn[i] - 1];

        val += amortMidat_lgmamortsopval1F(
            ex_ncpn[i],
            &(coupon[ex_cpn[i]]),
            &(df[ex_cpn[i]]),
            &(cvg[ex_cpn[i]]),
            &(cpn_G1[ex_cpn[i]]),
            ex_zeta1[i],
            ex_G1[i]);
    }

    return val;
}

/*	Value of European Cap within LGM 1F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one exercise date to the
   following */
double amortMidat_lgmcapval2F_b(
    int ncpn,             /*	Number of cash-flow dates, including
                                                  start and end date */
    double cpn[],         /*	Notional */
    double df[],          /*	Df to cash flow dates */
    double cvg[],         /*	cvg from i-1 to i */
    double cpn_G1[],      /*	G1 at cash-flow dates */
    double cpn_G2[],      /*	G2 at cash-flow dates */
    int    nex,           /*	Number of exercise dates */
    double ex_sstrike[],  /*	Strikes for cap */
    double floatcoupon[], /*	Float Coupon*/
    int    ex_cpn[],      /*	For each exercise date, first coupon
                                                  to be exercised */
    int ex_ncpn[],        /*	For each exercise date, number of coupons
                                                  to be exercised */
    double ex_zeta1[],    /*	Z1 at exercise date */
    double ex_zeta2[],    /*	Z2 at exercise date */
    double ex_zeta12[],   /*	Z12 at exercise date */
    double ex_G1[],       /*	G1 at exercise date */
    double ex_G2[])       /*	G2 at exercise date */

{
    double val;
    double coupon[MAX_CPN];
    int    i, j;

    val = 0.0;

    for (i = 0; i < nex; i++)
    {
        coupon[ex_cpn[i]] = -df[ex_cpn[i]];
        for (j = ex_cpn[i] + 1; j < ex_cpn[i] + ex_ncpn[i] - 1; ++j)
        {
            coupon[j] = floatcoupon[j];
        }
        coupon[ex_cpn[i] + ex_ncpn[i] - 1] = floatcoupon[j - 1] + df[ex_cpn[i] + ex_ncpn[i] - 1] +
                                             ex_sstrike[i] * df[ex_cpn[i] + ex_ncpn[i] - 1];

        val += amortMidat_lgmamortsopval2F(
            ex_ncpn[i],
            &(coupon[ex_cpn[i]]),
            &(df[ex_cpn[i]]),
            &(cvg[ex_cpn[i]]),
            &(cpn_G1[ex_cpn[i]]),
            &(cpn_G2[ex_cpn[i]]),
            ex_zeta1[i],
            ex_zeta2[i],
            ex_zeta12[i],
            ex_G1[i],
            ex_G2[i]);
    }

    return val;
}

/*	Value of European Cap within LGM 2F */
/*	NOTE: cap = sum of swaptions
                each swaption has for underlyings all the coupons from one exercise date to the
   following */
double amortMidat_lgmcapval2F(
    int ncpn,           /*	Number of cash-flow dates, including
                                                start and end date */
    double cpn[],       /*	Notional */
    double df[],        /*	Df to cash flow dates */
    double cvg[],       /*	cvg from i-1 to i */
    double cpn_G1[],    /*	G1 at cash-flow dates */
    double cpn_G2[],    /*	G2 at cash-flow dates */
    int    nex,         /*	Number of exercise dates */
    int    ex_cpn[],    /*	For each exercise date, first coupon
                                                to be exercised */
    int ex_ncpn[],      /*	For each exercise date, number of coupons
                                                to be exercised */
    double ex_zeta1[],  /*	Z1 at exercise date */
    double ex_zeta2[],  /*	Z2 at exercise date */
    double ex_zeta12[], /*	Z12 at exercise date */
    double ex_G1[],     /*	G1 at exercise date */
    double ex_G2[])     /*	G2 at exercise date */
{
    double val;
    double temp;
    int    i;

    val = 0.0;

    for (i = 0; i < nex; i++)
    {
        temp           = cpn[ex_cpn[i]];
        cpn[ex_cpn[i]] = -df[ex_cpn[i]];
        if (i < nex - 1)
        {
            cpn[ex_cpn[i] + ex_ncpn[i] - 1] =
                cpn[ex_cpn[i] + ex_ncpn[i] - 1] + df[ex_cpn[i] + ex_ncpn[i] - 1];
        }

        val += amortMidat_lgmamortsopval2F(
            ex_ncpn[i],
            &(cpn[ex_cpn[i]]),
            &(df[ex_cpn[i]]),
            &(cvg[ex_cpn[i]]),
            &(cpn_G1[ex_cpn[i]]),
            &(cpn_G1[ex_cpn[i]]),
            ex_zeta1[i],
            ex_zeta2[i],
            ex_zeta12[i],
            ex_G1[i],
            ex_G2[i]);

        cpn[ex_cpn[i]] = temp;
        if (i < nex - 1)
        {
            cpn[ex_cpn[i] + ex_ncpn[i] - 1] =
                cpn[ex_cpn[i] + ex_ncpn[i] - 1] - df[ex_cpn[i] + ex_ncpn[i] - 1];
        }
    }

    return val;
}

Err computeZeta(
    double  exer_time,
    double  lambda,
    double* sigma_date,
    double* sigma,
    int     sigma_n,
    double* exer_zeta)
{
    Err    err = NULL;
    int    i;
    double zeta;

    if (exer_time <= sigma_date[0])
    {
        zeta = sigma[0] * sigma[0] * (exp(2 * lambda * exer_time) - 1.0) / (2 * lambda);
    }
    else
    {
        zeta = sigma[0] * sigma[0] * (exp(2 * lambda * sigma_date[0]) - 1.0) / (2 * lambda);

        i = 1;
        while ((exer_time > sigma_date[i]) && (i < sigma_n))
        {
            zeta += sigma[i] * sigma[i] *
                    (exp(2 * lambda * sigma_date[i]) - exp(2 * lambda * sigma_date[i - 1])) /
                    (2 * lambda);
            ++i;
        }

        if (i < sigma_n)
        {
            zeta += sigma[i] * sigma[i] *
                    (exp(2 * lambda * exer_time) - exp(2 * lambda * sigma_date[i - 1])) /
                    (2 * lambda);
        }
        else
        {
            zeta += sigma[sigma_n - 1] * sigma[sigma_n - 1] *
                    (exp(2 * lambda * exer_time) - exp(2 * lambda * sigma_date[sigma_n - 1])) /
                    (2 * lambda);
        }
    }
    *exer_zeta = zeta;

    return err;
}

Err computeZeta_ts(
    double exer_time,

    double* lambda_time,
    double* lambda,
    int     lambda_n,
    double  gamma,

    double* sigma_date,
    double* sigma,
    int     sigma_n,

    double* exer_zeta)
{
    Err err = NULL;
    int i;

    for (i = 0; i < lambda_n; i++)
    {
        lambda[i] += gamma;
    }

    export_lgmcalczeta1_tauts(
        sigma_n,
        sigma_date,
        sigma,
        lambda_n,
        lambda_time,
        lambda,
        0.0,
        gamma,
        0.0,
        1,
        &exer_time,
        exer_zeta);

    for (i = 0; i < lambda_n; i++)
    {
        lambda[i] -= gamma;
    }

    return err;
}

Err europSwaption_clsdfrm(
    char*     yc_name,
    SrtUndPtr und,
    int       notperiod,
    long      lStartDate,
    long      lEndDate,
    char*     fixFreq,
    char*     fixBasis,
    char*     refRateName,
    double    exer_fee,
    double    strike,
    double    margin,
    char*     recPayStr,
    double*   price)
{
    Err         err = NULL;
    Date        lToday;
    String      und_name;
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts = NULL;

    SwapDP swapdp;
    long   float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    long    fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    int k, n;

    int float_index, fix_index;

    int     nbcoupon;
    double* coupon_dates;
    double  coupon[MAX_CPN];
    double  coupon_time[MAX_CPN];
    double  coupon_df[MAX_CPN];
    double  coupon_cvg[MAX_CPN];

    double cpn_G[MAX_CPN];
    double cpn_G2[MAX_CPN];

    double ex_G[MAX_CPN];
    double ex_G2[MAX_CPN];
    double exer_time[MAX_CPN];

    double lambda;
    double alpha;
    double gamma;
    double rho;

    //	double sig1;
    //	double sig2;
    //	double tau1;
    //	double tau2;

    double ex_zeta;
    double ex_zeta2;
    double ex_zeta12;

    double *sigma_date = NULL, *sigma = NULL;
    int     sigma_n;
    double *tau_date = NULL, *tau = NULL;
    int     tau_n;

    double ivalue;

    SrtReceiverType rec_pay;
    double          rec;

    /* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr, &rec_pay))
    {
        return err;
    }

    if (rec_pay == SRT_RECEIVER)
    {
        rec = 1;
    }
    else
    {
        rec = -1;
    }

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* Gets Today from underlying */
    und_name = get_underlying_name(und);
    lToday   = get_today_from_underlying(und);

    /* IF LGM */
    if (mdl_type == LGM || mdl_type == NEWLGM)
    {
        //----------------------------------------------------------------------------
        //---------SWAPDP Build Swap Schedule-----------------------------------------
        //----------------------------------------------------------------------------

        err = swp_f_initSwapDP(lStartDate, lEndDate, fixFreq, fixBasis, &swapdp);
        if (err)
        {
            goto FREE_RETURN;
        }

        err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
            &swapdp,
            lToday,
            refRateName,
            &float_pay_dates,
            &float_nb_pay_dates,
            &float_fixing_dates,
            &float_start_dates,
            &float_end_dates,
            &float_cvgs,
            &float_spreads,
            &float_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        float_pay_times = dvector(0, float_nb_pay_dates - 1);
        for (n = 0; n < float_nb_pay_dates; ++n)
        {
            float_pay_times[n] = (float_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        err = swp_f_make_FixedLegDatesAndCoverages(
            &swapdp,
            lToday,
            &fix_pay_dates,
            &fix_nb_pay_dates,
            &fix_start_dates,
            &fix_end_dates,
            &fix_cvgs,
            &fix_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
        for (n = 0; n < fix_nb_pay_dates; ++n)
        {
            fix_pay_times[n] = (fix_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        nbcoupon = float_nb_pay_dates;

        coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
        memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
        num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates, fix_pay_times);
        num_f_sort_vector(nbcoupon, coupon_dates);
        num_f_unique_vector(&nbcoupon, coupon_dates);

        for (k = 0; k < nbcoupon; ++k)
        {
            coupon[k] = 0;
        }

        ivalue         = 0;
        float_index    = 1;
        fix_index      = 1;
        ivalue         = -rec * (1 + exer_fee) * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon[0]      = -rec * (1 + exer_fee) * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_time[0] = coupon_dates[0];
        coupon_df[0]   = swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_cvg[0]  = 0;
        for (k = 1; k < nbcoupon; ++k)
        {
            coupon_time[k] = coupon_dates[k];
            coupon_df[k]   = swp_f_df(lToday, float_pay_dates[k], yc_name);
            coupon_cvg[k]  = float_cvgs[k - 1];
            if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
                (coupon_dates[k] != float_pay_times[float_index]))
            {
                ivalue = ivalue + rec * strike * fix_cvgs[fix_index - 1] *
                                      swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                coupon[k] += rec * strike * fix_cvgs[fix_index - 1] *
                             swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_dates[k] != fix_pay_times[fix_index]) &&
                (coupon_dates[k] == float_pay_times[float_index]))
            {
                ivalue = ivalue + -rec *
                                      ((float_spreads[float_index - 1] + margin) *
                                       float_cvgs[float_index - 1]) *
                                      swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                coupon[k] +=
                    -rec *
                    ((float_spreads[float_index - 1] + margin) * float_cvgs[float_index - 1]) *
                    swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < float_nb_pay_dates - 1)
                {
                    ivalue = ivalue + rec * strike * fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    coupon[k] += rec * strike * fix_cvgs[fix_index - 1] *
                                 swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue = ivalue - rec *
                                          ((float_spreads[float_index - 1] + margin) *
                                           float_cvgs[float_index - 1]) *
                                          swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] +=
                        -rec *
                        ((float_spreads[float_index - 1] + margin) * float_cvgs[float_index - 1]) *
                        swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
                else
                {
                    ivalue = ivalue + rec * strike * fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    coupon[k] += rec * strike * fix_cvgs[fix_index - 1] *
                                 swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue = ivalue - rec *
                                          (-1 + (float_spreads[float_index - 1] + margin) *
                                                    float_cvgs[float_index - 1]) *
                                          swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] += -rec *
                                 (-1 + (float_spreads[float_index - 1] + margin) *
                                           float_cvgs[float_index - 1]) *
                                 swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
            }
        }

        //-----------------------------------------------------------------
        //----------------End Of Build Swap Schedule-----------------------
        //-----------------------------------------------------------------
        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            return err;
        }

        if (tau_n > 1)
        {
            err = "Constant Tau required";
            return err;
        }

        lambda = 1.0 / tau[0];

        exer_time[0] =
            (add_unit(lStartDate, -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) / 365.0;

        if (mdl_dim == ONE_FAC)
        {
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            *price = lgmopval1F(nbcoupon, coupon, cpn_G, ex_zeta, ex_G[0]);

            if (*price < DMAX(0.0, ivalue))
            {
                *price = DMAX(0.0, ivalue);
            }
        }
        else
        {
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            export_lgmsetupG2(lambda, gamma, nbcoupon, coupon_time, cpn_G2, 1, exer_time, ex_G2);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda + gamma, sigma_date, sigma, sigma_n, &ex_zeta2);
            if (err)
            {
                return err;
            }
            ex_zeta2 = ex_zeta2 * alpha * alpha;

            err = computeZeta(
                exer_time[0], lambda + 0.5 * gamma, sigma_date, sigma, sigma_n, &ex_zeta12);
            if (err)
            {
                return err;
            }
            ex_zeta12 = ex_zeta12 * alpha * rho;

            *price = amort_lgmopval2F(
                nbcoupon, coupon, cpn_G, cpn_G2, ex_zeta, ex_zeta2, ex_zeta12, ex_G[0], ex_G2[0]);

            if (*price < DMAX(0.0, ivalue))
            {
                *price = DMAX(0.0, ivalue);
            }
        }
    }
    else
    {
        err = "Not a LGM underlying ";
    }

FREE_RETURN:

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);

    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);
    if (fix_cvgs)
        free(fix_cvgs);

    if (sigma_date)
        free(sigma_date);
    if (sigma)
        free(sigma);

    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);

    if (coupon_dates)
        free(coupon_dates);

    return err;
}

Err europAmortSwaption_clsdfrm_new(
    char*     yc_name,
    SrtUndPtr und,
    int       notperiod,

    char*   fixbasis,
    int     NFix,
    long*   FixStartDates,
    long*   FixEndDates,
    double* FixRates,
    double* FixNotionals,

    char*   refrate,
    int     NFloat,
    long*   FloatStartDates,
    long*   FloatEndDates,
    double* FloatMargins,
    double* FloatSpreads,
    double* FloatNotionals,

    double  exer_fee,
    char*   recPayStr,
    double* price)
{
    Err         err = NULL;
    Date        lToday;
    String      und_name;
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts = NULL;

    double* float_pay_times = NULL;
    double* fix_pay_times   = NULL;

    int k, n;

    int float_index, fix_index;

    int     nbcoupon;
    double* coupon_dates = NULL;
    double  coupon[MAX_CPN];
    double  coupon_time[MAX_CPN];
    double  coupon_df[MAX_CPN];
    double  coupon_cvg[MAX_CPN];

    double cpn_G[MAX_CPN];
    double cpn_G2[MAX_CPN];

    double ex_G[MAX_CPN];
    double ex_G2[MAX_CPN];
    double exer_time[MAX_CPN];

    double  lambda;
    double* pdlambda;
    double* pdlambda_time;
    int     lambda_n;

    double alpha;
    double gamma;
    double rho;

    double ex_zeta;
    double ex_zeta2;
    double ex_zeta12;

    double *sigma_date = NULL, *sigma = NULL;
    int     sigma_n;
    double *tau_date = NULL, *tau = NULL;
    int     tau_n;

    double ivalue;

    SrtReceiverType rec_pay;
    double          rec;

    SrtBasisCode fixbas, floatbas;
    int          floatfreq;

    double* FixCoverages   = NULL;
    double* FloatCoverages = NULL;

    err = swp_f_get_ref_rate_details(refrate, &floatbas, &floatfreq);

    /* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr, &rec_pay))
    {
        goto FREE_RETURN;
        return err;
    }

    if (rec_pay == SRT_RECEIVER)
    {
        rec = 1;
    }
    else
    {
        for (k = 0; k < NFix; ++k)
        {
            FixNotionals[k] = -FixNotionals[k];
        }

        for (k = 0; k < NFloat; ++k)
        {
            FloatNotionals[k] = -FloatNotionals[k];
        }
    }

    FixCoverages   = (double*)calloc(NFix, sizeof(double));
    FloatCoverages = (double*)calloc(NFloat, sizeof(double));

    err = interp_basis(fixbasis, &fixbas);
    if (err)
    {
        goto FREE_RETURN;
        return err;
    }
    for (k = 0; k < NFix; ++k)
    {
        FixCoverages[k] = coverage(FixStartDates[k], FixEndDates[k], fixbas);
    }

    if (err)
    {
        goto FREE_RETURN;
        return err;
    }
    for (k = 0; k < NFloat; ++k)
    {
        FloatCoverages[k] = coverage(FloatStartDates[k], FloatEndDates[k], floatbas);
    }

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* Gets Today from underlying */
    und_name = get_underlying_name(und);
    lToday   = get_today_from_underlying(und);

    /* IF LGM */
    if (mdl_type == LGM || mdl_type == NEWLGM)
    {
        //----------------------------------------------------------------------------
        //---------SWAPDP Build Swap Schedule-----------------------------------------
        //----------------------------------------------------------------------------

        float_pay_times    = dvector(0, NFloat);
        float_pay_times[0] = (FloatStartDates[0] - lToday) * YEARS_IN_DAY;
        for (n = 1; n <= NFloat; ++n)
        {
            float_pay_times[n] = (FloatEndDates[n - 1] - lToday) * YEARS_IN_DAY;
        }

        fix_pay_times    = dvector(0, NFix);
        fix_pay_times[0] = (FixStartDates[0] - lToday) * YEARS_IN_DAY;
        for (n = 1; n <= NFix; ++n)
        {
            fix_pay_times[n] = (FixEndDates[n - 1] - lToday) * YEARS_IN_DAY;
        }

        nbcoupon = NFloat + 1;

        coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
        memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
        num_f_concat_vector(&nbcoupon, &coupon_dates, NFix + 1, fix_pay_times);
        num_f_sort_vector(nbcoupon, coupon_dates);
        num_f_unique_vector(&nbcoupon, coupon_dates);

        for (k = 0; k < nbcoupon; ++k)
        {
            coupon[k] = 0;
        }

        ivalue      = 0;
        float_index = 1;
        fix_index   = 1;
        ivalue    = -(FloatNotionals[0] + exer_fee) * swp_f_df(lToday, FloatStartDates[0], yc_name);
        coupon[0] = -(FloatNotionals[0] + exer_fee) * swp_f_df(lToday, FloatStartDates[0], yc_name);
        coupon_time[0] = coupon_dates[0];
        coupon_df[0]   = swp_f_df(lToday, FloatStartDates[0], yc_name);
        coupon_cvg[0]  = 0;
        for (k = 1; k < nbcoupon; ++k)
        {
            coupon_time[k] = coupon_dates[k];
            coupon_df[k]   = swp_f_df(lToday, FloatEndDates[k - 1], yc_name);
            coupon_cvg[k]  = FloatCoverages[k - 1];
            if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
                (coupon_dates[k] != float_pay_times[float_index]))
            {
                ivalue = ivalue + FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                                      FixCoverages[fix_index - 1] *
                                      swp_f_df(lToday, FixEndDates[fix_index - 1], yc_name);
                coupon[k] += FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                             FixCoverages[fix_index - 1] *
                             swp_f_df(lToday, FixEndDates[fix_index], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_dates[k] != fix_pay_times[fix_index]) &&
                (coupon_dates[k] == float_pay_times[float_index]))
            {
                ivalue =
                    ivalue + -(FloatNotionals[float_index] - FloatNotionals[float_index - 1] +
                               FloatNotionals[float_index - 1] *
                                   (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                                   FloatCoverages[float_index - 1]) *
                                 swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                coupon[k] += -(FloatNotionals[float_index] - FloatNotionals[float_index - 1] +
                               FloatNotionals[float_index - 1] *
                                   (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                                   FloatCoverages[float_index - 1]) *
                             swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < NFloat)
                {
                    ivalue = ivalue + FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                                          FixCoverages[fix_index - 1] *
                                          swp_f_df(lToday, FixEndDates[fix_index - 1], yc_name);
                    coupon[k] += FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                                 FixCoverages[fix_index - 1] *
                                 swp_f_df(lToday, FixEndDates[fix_index - 1], yc_name);
                    fix_index = fix_index + 1;

                    ivalue = ivalue -
                             (FloatNotionals[float_index] - FloatNotionals[float_index - 1] +
                              FloatNotionals[float_index - 1] *
                                  (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                                  FloatCoverages[float_index - 1]) *
                                 swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                    coupon[k] +=
                        -(FloatNotionals[float_index] - FloatNotionals[float_index - 1] +
                          FloatNotionals[float_index - 1] *
                              (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                              FloatCoverages[float_index - 1]) *
                        swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                    float_index = float_index + 1;
                }
                else
                {
                    ivalue = ivalue + FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                                          FixCoverages[fix_index - 1] *
                                          swp_f_df(lToday, FixEndDates[fix_index - 1], yc_name);
                    coupon[k] += FixRates[fix_index - 1] * FixNotionals[fix_index - 1] *
                                 FixCoverages[fix_index - 1] *
                                 swp_f_df(lToday, FixEndDates[fix_index - 1], yc_name);
                    fix_index = fix_index + 1;

                    ivalue = ivalue +
                             -(-FloatNotionals[float_index - 1] +
                               FloatNotionals[float_index - 1] *
                                   (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                                   FloatCoverages[float_index - 1]) *
                                 swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                    coupon[k] +=
                        -(-FloatNotionals[float_index - 1] +
                          FloatNotionals[float_index - 1] *
                              (FloatSpreads[float_index - 1] + FloatMargins[float_index - 1]) *
                              FloatCoverages[float_index - 1]) *
                        swp_f_df(lToday, FloatEndDates[float_index - 1], yc_name);
                    float_index = float_index + 1;
                }
            }
        }

        //-----------------------------------------------------------------
        //----------------End Of Build Swap Schedule-----------------------
        //-----------------------------------------------------------------
        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            goto FREE_RETURN;
            return err;
        }

        //// previously ...
#if 0
		if(tau_n>1)
		{
			err = "Constant Tau required";
			goto FREE_RETURN;
			return err;
		}
#endif

        if (tau_n == 1)
        {
            lambda = 1.0 / tau[0];

            exer_time[0] =
                (add_unit(FixStartDates[0], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) /
                365.0;

            if (mdl_dim == ONE_FAC)
            {
                export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                *price = lgmopval1F(nbcoupon, coupon, cpn_G, ex_zeta, ex_G[0]);

                if (*price < DMAX(0.0, ivalue))
                {
                    *price = DMAX(0.0, ivalue);
                }
            }
            else
            {
                export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                export_lgmsetupG2(
                    lambda, gamma, nbcoupon, coupon_time, cpn_G2, 1, exer_time, ex_G2);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta(
                    exer_time[0], lambda + gamma, sigma_date, sigma, sigma_n, &ex_zeta2);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }
                ex_zeta2 = ex_zeta2 * alpha * alpha;

                err = computeZeta(
                    exer_time[0], lambda + 0.5 * gamma, sigma_date, sigma, sigma_n, &ex_zeta12);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }
                ex_zeta12 = ex_zeta12 * alpha * rho;

                *price = amort_lgmopval2F(
                    nbcoupon,
                    coupon,
                    cpn_G,
                    cpn_G2,
                    ex_zeta,
                    ex_zeta2,
                    ex_zeta12,
                    ex_G[0],
                    ex_G2[0]);

                if (*price < DMAX(0.0, ivalue))
                {
                    *price = DMAX(0.0, ivalue);
                }
            }

        }  ///

        else  // nTau > 1
        {
            lambda_n      = tau_n;
            pdlambda      = dvector(0, lambda_n);
            pdlambda_time = dvector(0, lambda_n);

            for (k = 0; k < lambda_n; ++k)
            {
                pdlambda_time[k] = tau_date[k];
                pdlambda[k]      = 1. / tau[k];
            }

            exer_time[0] =
                (add_unit(FixStartDates[0], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) /
                365.0;

            if (mdl_dim == ONE_FAC)
            {
                export_lgmsetupG_ts(
                    lambda_n,
                    pdlambda_time,
                    pdlambda,
                    nbcoupon,
                    coupon_time,
                    cpn_G,
                    1,
                    exer_time,
                    ex_G);

                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta_ts(
                    exer_time[0],
                    pdlambda_time,
                    pdlambda,
                    lambda_n,
                    0.,
                    sigma_date,
                    sigma,
                    sigma_n,
                    &ex_zeta);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                *price = lgmopval1F(nbcoupon, coupon, cpn_G, ex_zeta, ex_G[0]);

                if (*price < DMAX(0.0, ivalue))
                {
                    *price = DMAX(0.0, ivalue);
                }
            }
            else
            {
                export_lgmsetupG_ts(
                    lambda_n,
                    pdlambda_time,
                    pdlambda,
                    nbcoupon,
                    coupon_time,
                    cpn_G,
                    1,
                    exer_time,
                    ex_G);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                export_lgmsetupG2_ts(
                    lambda_n,
                    pdlambda_time,
                    pdlambda,
                    gamma,
                    nbcoupon,
                    coupon_time,
                    cpn_G2,
                    1,
                    exer_time,
                    ex_G2);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta_ts(
                    exer_time[0],
                    pdlambda_time,
                    pdlambda,
                    lambda_n,
                    0.,
                    sigma_date,
                    sigma,
                    sigma_n,
                    &ex_zeta);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }

                err = computeZeta_ts(
                    exer_time[0],
                    pdlambda_time,
                    pdlambda,
                    lambda_n,
                    gamma,
                    sigma_date,
                    sigma,
                    sigma_n,
                    &ex_zeta2);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }
                ex_zeta2 = ex_zeta2 * alpha * alpha;

                err = computeZeta_ts(
                    exer_time[0],
                    pdlambda_time,
                    pdlambda,
                    lambda_n,
                    0.5 * gamma,
                    sigma_date,
                    sigma,
                    sigma_n,
                    &ex_zeta12);
                if (err)
                {
                    goto FREE_RETURN;
                    return err;
                }
                ex_zeta12 = ex_zeta12 * alpha * rho;

                *price = amort_lgmopval2F(
                    nbcoupon,
                    coupon,
                    cpn_G,
                    cpn_G2,
                    ex_zeta,
                    ex_zeta2,
                    ex_zeta12,
                    ex_G[0],
                    ex_G2[0]);

                if (*price < DMAX(0.0, ivalue))
                {
                    *price = DMAX(0.0, ivalue);
                }
            }
        }
    }
    else
    {
        err = "Not a LGM underlying ";
    }

FREE_RETURN:

    if (FixCoverages)
        free(FixCoverages);
    if (FloatCoverages)
        free(FloatCoverages);

    if (float_pay_times)
        free_dvector(float_pay_times, 0, NFloat);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, NFix);

    if (sigma_date)
        free(sigma_date);
    if (sigma)
        free(sigma);

    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);

    if (coupon_dates)
        free(coupon_dates);

    return err;
}

Err europAmortSwaption_clsdfrm(
    char*     yc_name,
    SrtUndPtr und,
    int       notperiod,
    long      lStartDate,
    long      lEndDate,
    char*     fixFreq,
    char*     fixBasis,
    char*     refRateName,
    double    exer_fee,
    int       Nfix,
    double*   fixNotional,
    double*   fixRates,
    int       Nfloat,
    double*   floatNotional,
    double*   margin,
    char*     recPayStr,
    double*   price)
{
    Err         err = NULL;
    Date        lToday;
    String      und_name;
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts = NULL;

    SwapDP swapdp;
    long   float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    long    fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    int k, n;

    int float_index, fix_index;

    int     nbcoupon;
    double* coupon_dates;
    double  coupon[MAX_CPN];
    double  coupon_time[MAX_CPN];
    double  coupon_df[MAX_CPN];
    double  coupon_cvg[MAX_CPN];

    double cpn_G[MAX_CPN];
    double cpn_G2[MAX_CPN];

    double ex_G[MAX_CPN];
    double ex_G2[MAX_CPN];
    double exer_time[MAX_CPN];

    double lambda;
    double alpha;
    double gamma;
    double rho;

    //	double sig1;
    //	double sig2;
    //	double tau1;
    //	double tau2;

    double ex_zeta;
    double ex_zeta2;
    double ex_zeta12;

    double *sigma_date, *sigma;
    int     sigma_n;
    double *tau_date, *tau;
    int     tau_n;

    double ivalue;

    SrtReceiverType rec_pay;
    double          rec;

    /* Interprets the Rec Pay string into a type  */
    if (err = interp_rec_pay(recPayStr, &rec_pay))
    {
        return err;
    }

    if (rec_pay == SRT_RECEIVER)
    {
        rec = 1;
    }
    else
    {
        for (k = 0; k < Nfix; ++k)
        {
            fixNotional[k] = -fixNotional[k];
        }

        for (k = 0; k < Nfloat; ++k)
        {
            floatNotional[k] = -floatNotional[k];
        }
    }

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* Gets Today from underlying */
    und_name = get_underlying_name(und);
    lToday   = get_today_from_underlying(und);

    /* IF LGM */
    if (mdl_type == LGM || mdl_type == NEWLGM)
    {
        //----------------------------------------------------------------------------
        //---------SWAPDP Build Swap Schedule-----------------------------------------
        //----------------------------------------------------------------------------

        err = swp_f_initSwapDP(lStartDate, lEndDate, fixFreq, fixBasis, &swapdp);
        if (err)
        {
            goto FREE_RETURN;
        }

        err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
            &swapdp,
            lToday,
            refRateName,
            &float_pay_dates,
            &float_nb_pay_dates,
            &float_fixing_dates,
            &float_start_dates,
            &float_end_dates,
            &float_cvgs,
            &float_spreads,
            &float_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        float_pay_times = dvector(0, float_nb_pay_dates - 1);
        for (n = 0; n < float_nb_pay_dates; ++n)
        {
            float_pay_times[n] = (float_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        err = swp_f_make_FixedLegDatesAndCoverages(
            &swapdp,
            lToday,
            &fix_pay_dates,
            &fix_nb_pay_dates,
            &fix_start_dates,
            &fix_end_dates,
            &fix_cvgs,
            &fix_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
        for (n = 0; n < fix_nb_pay_dates; ++n)
        {
            fix_pay_times[n] = (fix_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        nbcoupon = float_nb_pay_dates;

        coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
        memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
        num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates, fix_pay_times);
        num_f_sort_vector(nbcoupon, coupon_dates);
        num_f_unique_vector(&nbcoupon, coupon_dates);

        for (k = 0; k < nbcoupon; ++k)
        {
            coupon[k] = 0;
        }

        ivalue      = 0;
        float_index = 1;
        fix_index   = 1;
        ivalue    = -(floatNotional[0] + exer_fee) * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon[0] = -(floatNotional[0] + exer_fee) * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_time[0] = coupon_dates[0];
        coupon_df[0]   = swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_cvg[0]  = 0;
        for (k = 1; k < nbcoupon; ++k)
        {
            coupon_time[k] = coupon_dates[k];
            coupon_df[k]   = swp_f_df(lToday, float_pay_dates[k], yc_name);
            coupon_cvg[k]  = float_cvgs[k - 1];
            if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
                (coupon_dates[k] != float_pay_times[float_index]))
            {
                ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                      fix_cvgs[fix_index - 1] *
                                      swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                coupon[k] += fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_dates[k] != fix_pay_times[fix_index]) &&
                (coupon_dates[k] == float_pay_times[float_index]))
            {
                ivalue = ivalue + -(floatNotional[float_index] - floatNotional[float_index - 1] +
                                    floatNotional[float_index - 1] *
                                        (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                        float_cvgs[float_index - 1]) *
                                      swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < float_nb_pay_dates - 1)
                {
                    ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                          fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    coupon[k] += fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                 fix_cvgs[fix_index - 1] *
                                 swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue =
                        ivalue - (floatNotional[float_index] - floatNotional[float_index - 1] +
                                  floatNotional[float_index - 1] *
                                      (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                      float_cvgs[float_index - 1]) *
                                     swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                 swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
                else
                {
                    ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                          fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    coupon[k] += fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                 fix_cvgs[fix_index - 1] *
                                 swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue =
                        ivalue + -(-floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                     swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] += -(-floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                 swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
            }
        }

        //-----------------------------------------------------------------
        //----------------End Of Build Swap Schedule-----------------------
        //-----------------------------------------------------------------
        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            return err;
        }

        if (tau_n > 1)
        {
            err = "Constant Tau required";
            return err;
        }

        lambda = 1.0 / tau[0];

        exer_time[0] =
            (add_unit(lStartDate, -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) / 365.0;

        if (mdl_dim == ONE_FAC)
        {
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            *price = lgmopval1F(nbcoupon, coupon, cpn_G, ex_zeta, ex_G[0]);

            if (*price < DMAX(0.0, ivalue))
            {
                *price = DMAX(0.0, ivalue);
            }
        }
        else
        {
            /*			err = get_underlying_ts(und, &ts);
                                    if(err)
                                    {
                                            return err;
                                    }

                                    err = find_2f_rho (exer_time[0], ts, &rho);
                                    if(err)
                                    {
                                            return err;
                                    }

                                    err = find_2f_sig (exer_time[0], ts, &sig1, &sig2);
                                    if(err)
                                    {
                                            return err;
                                    }
                                    alpha = sig2/sig1;

                                    err = find_2f_tau (exer_time[0], ts, &tau1, &tau2);
                                    if(err)
                                    {
                                            return err;
                                    }
                                    gamma = 1/tau2 - 1/tau1;
            */
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            export_lgmsetupG2(lambda, gamma, nbcoupon, coupon_time, cpn_G2, 1, exer_time, ex_G2);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda + gamma, sigma_date, sigma, sigma_n, &ex_zeta2);
            if (err)
            {
                return err;
            }
            ex_zeta2 = ex_zeta2 * alpha * alpha;

            err = computeZeta(
                exer_time[0], lambda + 0.5 * gamma, sigma_date, sigma, sigma_n, &ex_zeta12);
            if (err)
            {
                return err;
            }
            ex_zeta12 = ex_zeta12 * alpha * rho;

            *price = amort_lgmopval2F(
                nbcoupon, coupon, cpn_G, cpn_G2, ex_zeta, ex_zeta2, ex_zeta12, ex_G[0], ex_G2[0]);

            if (*price < DMAX(0.0, ivalue))
            {
                *price = DMAX(0.0, ivalue);
            }
        }
    }
    else
    {
        err = "Not a LGM underlying ";
    }

FREE_RETURN:

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);

    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);
    if (fix_cvgs)
        free(fix_cvgs);

    if (sigma_date)
        free(sigma_date);
    if (sigma)
        free(sigma);

    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);

    if (coupon_dates)
        free(coupon_dates);

    return err;
}

Err europZCSwaption_clsdfrm(
    char*     yc_name,
    SrtUndPtr und,
    int       notperiod,
    long      lStartDate,
    long      lEndDate,
    char*     fixFreq,
    char*     fixBasis,
    char*     refRateName,
    int       Nfix,
    double*   fixNotional,
    double*   fixRates,
    int       Nfloat,
    double*   floatNotional,
    double*   margin,
    double*   price)
{
    Err         err = NULL;
    Date        lToday;
    String      und_name;
    SrtMdlType  mdl_type;
    SrtMdlDim   mdl_dim;
    TermStruct* ts = NULL;

    SwapDP swapdp;
    long   float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    long    fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    int k, n;

    int float_index, fix_index;

    int     nbcoupon;
    double  finalFixCoupon;
    double  accruedTerm;
    double* coupon_dates;
    double  coupon[MAX_CPN];
    double  coupon_time[MAX_CPN];
    double  coupon_df[MAX_CPN];
    double  coupon_cvg[MAX_CPN];

    double cpn_G[MAX_CPN];
    double cpn_G2[MAX_CPN];

    double ex_G[MAX_CPN];
    double ex_G2[MAX_CPN];
    double exer_time[MAX_CPN];

    double lambda;
    double alpha;
    double gamma;
    double rho;

    double ex_zeta;
    double ex_zeta2;
    double ex_zeta12;

    double *sigma_date, *sigma;
    int     sigma_n;
    double *tau_date, *tau;
    int     tau_n;

    double ivalue;

    /* Gets the Model type (LGM, Cheyette,...) */
    err = get_underlying_mdltype(und, &mdl_type);

    /* Gets the Number of Factors in the Model */
    err = get_underlying_mdldim(und, &mdl_dim);

    /* Gets Today from underlying */
    und_name = get_underlying_name(und);
    lToday   = get_today_from_underlying(und);

    /* IF LGM */
    if (mdl_type == LGM || mdl_type == NEWLGM)
    {
        //----------------------------------------------------------------------------
        //---------SWAPDP Build Swap Schedule-----------------------------------------
        //----------------------------------------------------------------------------

        err = swp_f_initSwapDP(lStartDate, lEndDate, fixFreq, fixBasis, &swapdp);
        if (err)
        {
            goto FREE_RETURN;
        }

        err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
            &swapdp,
            lToday,
            refRateName,
            &float_pay_dates,
            &float_nb_pay_dates,
            &float_fixing_dates,
            &float_start_dates,
            &float_end_dates,
            &float_cvgs,
            &float_spreads,
            &float_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        float_pay_times = dvector(0, float_nb_pay_dates - 1);
        for (n = 0; n < float_nb_pay_dates; ++n)
        {
            float_pay_times[n] = (float_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        err = swp_f_make_FixedLegDatesAndCoverages(
            &swapdp,
            lToday,
            &fix_pay_dates,
            &fix_nb_pay_dates,
            &fix_start_dates,
            &fix_end_dates,
            &fix_cvgs,
            &fix_nb_dates);
        if (err)
        {
            goto FREE_RETURN;
        }
        fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
        for (n = 0; n < fix_nb_pay_dates; ++n)
        {
            fix_pay_times[n] = (fix_pay_dates[n] - lToday) * YEARS_IN_DAY;
        }

        nbcoupon = float_nb_pay_dates;

        coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
        memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
        num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates, fix_pay_times);
        num_f_sort_vector(nbcoupon, coupon_dates);
        num_f_unique_vector(&nbcoupon, coupon_dates);

        for (k = 0; k < nbcoupon; ++k)
        {
            coupon[k] = 0;
        }

        ivalue         = 0;
        float_index    = 1;
        fix_index      = 1;
        ivalue         = -floatNotional[0] * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon[0]      = -floatNotional[0] * swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_time[0] = coupon_dates[0];
        coupon_df[0]   = swp_f_df(lToday, float_pay_dates[0], yc_name);
        coupon_cvg[0]  = 0;
        for (k = 1; k < nbcoupon; ++k)
        {
            coupon_time[k] = coupon_dates[k];
            coupon_df[k]   = swp_f_df(lToday, float_pay_dates[k], yc_name);
            coupon_cvg[k]  = float_cvgs[k - 1];
            if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
                (coupon_dates[k] != float_pay_times[float_index]))
            {
                ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                      fix_cvgs[fix_index - 1] *
                                      swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                //				coupon[k] += fixRates[fix_index-1] *
                //fixNotional[fix_index-1] * fix_cvgs[fix_index-1] * swp_f_df (lToday,
                //fix_pay_dates[fix_index], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_dates[k] != fix_pay_times[fix_index]) &&
                (coupon_dates[k] == float_pay_times[float_index]))
            {
                ivalue = ivalue + -(floatNotional[float_index] - floatNotional[float_index - 1] +
                                    floatNotional[float_index - 1] *
                                        (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                        float_cvgs[float_index - 1]) *
                                      swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < float_nb_pay_dates - 1)
                {
                    ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                          fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    //					coupon[k] += fixRates[fix_index-1] * fixNotional[fix_index-1] *
                    //fix_cvgs[fix_index-1] * swp_f_df (lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue =
                        ivalue - (floatNotional[float_index] - floatNotional[float_index - 1] +
                                  floatNotional[float_index - 1] *
                                      (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                      float_cvgs[float_index - 1]) *
                                     swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                 swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
                else
                {
                    ivalue = ivalue + fixRates[fix_index - 1] * fixNotional[fix_index - 1] *
                                          fix_cvgs[fix_index - 1] *
                                          swp_f_df(lToday, fix_pay_dates[fix_index], yc_name);
                    //					coupon[k] += fixRates[fix_index-1] * fixNotional[fix_index-1] *
                    //fix_cvgs[fix_index-1] * swp_f_df (lToday, fix_pay_dates[fix_index], yc_name);
                    fix_index = fix_index + 1;

                    ivalue =
                        ivalue + -(-floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                     swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    coupon[k] += -(-floatNotional[float_index - 1] +
                                   floatNotional[float_index - 1] *
                                       (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                       float_cvgs[float_index - 1]) *
                                 swp_f_df(lToday, float_pay_dates[float_index], yc_name);
                    float_index = float_index + 1;
                }
            }
        }

        finalFixCoupon = 0.0;
        accruedTerm    = 1.0;
        for (k = fix_nb_dates - 1; k >= 0; --k)
        {
            finalFixCoupon += accruedTerm * fixRates[k] * fixNotional[k] * fix_cvgs[k];
            accruedTerm = accruedTerm * (1 + fixRates[k] * fix_cvgs[k]);
        }

        coupon[nbcoupon - 1] +=
            finalFixCoupon * swp_f_df(lToday, fix_end_dates[fix_nb_dates - 1], yc_name);

        //-----------------------------------------------------------------
        //----------------End Of Build Swap Schedule-----------------------
        //-----------------------------------------------------------------
        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            return err;
        }

        if (tau_n > 1)
        {
            err = "Constant Tau required";
            return err;
        }

        lambda = 1.0 / tau[0];

        exer_time[0] =
            (add_unit(lStartDate, -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) - lToday) / 365.0;

        if (mdl_dim == ONE_FAC)
        {
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            *price = lgmopval1F(nbcoupon, coupon, cpn_G, ex_zeta, ex_G[0]);

            if (*price < ivalue)
            {
                *price = ivalue;
            }
        }
        else
        {
            export_lgmsetupG(lambda, nbcoupon, coupon_time, cpn_G, 1, exer_time, ex_G);
            if (err)
            {
                return err;
            }

            export_lgmsetupG2(lambda, gamma, nbcoupon, coupon_time, cpn_G2, 1, exer_time, ex_G2);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda, sigma_date, sigma, sigma_n, &ex_zeta);
            if (err)
            {
                return err;
            }

            err = computeZeta(exer_time[0], lambda + gamma, sigma_date, sigma, sigma_n, &ex_zeta2);
            if (err)
            {
                return err;
            }
            ex_zeta2 = ex_zeta2 * alpha * alpha;

            err = computeZeta(
                exer_time[0], lambda + 0.5 * gamma, sigma_date, sigma, sigma_n, &ex_zeta12);
            if (err)
            {
                return err;
            }
            ex_zeta12 = ex_zeta12 * alpha * rho;

            *price = amort_lgmopval2F(
                nbcoupon, coupon, cpn_G, cpn_G2, ex_zeta, ex_zeta2, ex_zeta12, ex_G[0], ex_G2[0]);

            if (*price < ivalue)
            {
                *price = ivalue;
            }
        }
    }
    else
    {
        err = "Not a LGM underlying ";
    }

FREE_RETURN:

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);

    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);

    if (sigma_date)
        free(sigma_date);
    if (sigma)
        free(sigma);

    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);

    if (coupon_dates)
        free(coupon_dates);

    return err;
}

/*	Calibrate zeta to amortized diagonal given G: 2F case */
Err AmortMidat_lgmcalibzeta2F(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn[],      /*	Discounted Cash-Flow */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    int    ex_bool[],  /*	Exercise or not */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    double ex_G1[],    /*	G1 at exercise date */
    double ex_G2[],    /*	G2 at exercise date */
    // double			strike[],							/*	Strikes
    // */
    double mkt_price[],     /*	Market prices */
    double ex_fee[],        /*	Exercise Fees for Diagonal */
    double ex_zeta[],       /*	Output: zetas (1) */
    double floatnotional[], /*	Float Notional	*/
    /*	Lambda, Alpha, gamma, rho */
    double lambda,
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    double s1, s2, s_last, ds;
    int    i, j, it;
    double l1 = lambda, l2 = lambda + gamma;
    double z1, z2, z12;
    double q1, q2, q12;
    double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
    double pr1, pr2;
    double exp_fact1, exp_fact2, exp_fact12;
    double temp;
    double minvar, maxvar, d;
    double quad_var;
    int    niter;
    int    firstcalibdone;
    double convert_ratio;
    int    k, l, m;

    static double coupon[MAX_CPN];

    double** res_iter = NULL;

    for (i = 0; i < ncpn; ++i)
    {
        coupon[i] = cpn[i];
    }

    i = 0;
    while (i < nex && !ex_bool[i])
    {
        i++;
    }

    if (i == nex)
    {
        return "No option to calibrate";
    }

    res_iter = dmatrix(0, 5 * NITER, 0, 1);

    if (!res_iter)
    {
        return "Memory allocation faillure in lgmprcapgiventauts_momentum_dlm";
    }

    convert_ratio = sqrt(
        lambda / LAMBDA_REF * (1.0 - exp(-2.0 * LAMBDA_REF * cpn_time[i])) /
        (1.0 - exp(-2.0 * lambda * cpn_time[i])));

    /*	Initial guess: total vol of 0.005 - 0.07 bps */
    /*	s = local variance */
    s1 = 0.02;
    s1 *= convert_ratio / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s1 *= s1;

    s2 = 0.04;
    s2 *= convert_ratio / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s2 *= s2;

    /*	Initialisation */
    quad_var = 0.0;
    s_last   = s1;
    ds       = s2 / s1;
    zeta1 = zeta2 = zeta12 = 0.0;
    t                      = 0.0;

    /*	If only one option, no skipping last */
    if (nex == 1)
    {
        skip_last = 0;
    }

    firstcalibdone = 0;
    for (i = 0; i < nex - (skip_last > 0); i++)
    {
        j         = ex_cpn[i];
        exp_fact1 = ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
        exp_fact2 = ((exp(2 * l2 * ex_time[i]) - exp(2 * l2 * t)) / 2 / l2) /
                    ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
        exp_fact12 = ((exp((l1 + l2) * ex_time[i]) - exp((l1 + l2) * t)) / (l1 + l2)) /
                     ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);

        if ((ex_bool[i] == 1) || (i == 0))
        {
            /*	Initial guess: last vol */
            s1 = s_last;
            /*	Initial guess 2: last vol * ds */
            s2 = ds * s1;

            /*	First price */
            dz1 = s1 * exp_fact1;
            z1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            z2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            z12  = zeta12 + dz12;

            coupon[j] = -floatnotional[j] * cpn_df[j] - ex_fee[i] * cpn_df[j];

            pr1 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                z1,
                z2,
                z12,
                ex_G1[i],
                ex_G2[i]);

            if (fabs(mkt_price[i] - pr1) < CALPRESNOT)
            {
                s_last = s1;
                quad_var += s1 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = z1;
                zeta2              = z2;
                zeta12             = z12;
                t                  = ex_time[i];
                continue;
            }

            /*	Second price */
            dz1 = s2 * exp_fact1;
            q1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            q2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            q12  = zeta12 + dz12;

            pr2 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                q1,
                q2,
                q12,
                ex_G1[i],
                ex_G2[i]);

            if (fabs(mkt_price[i] - pr2) < CALPRESNOT)
            {
                s_last = s2;
                quad_var += s2 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = q1;
                zeta2              = q2;
                zeta12             = q12;
                t                  = ex_time[i];
                continue;
            }

            /*	First option: no limits, 15 iterations */
            if (firstcalibdone == 0)
            {
                niter  = 3 * NITER;
                minvar = 1.0e-16;
                maxvar = 100.0;

                /* make sure we are between the bounds */
                it = 1;

                while (pr1 > mkt_price[i] && it < 6)
                {
                    s1 *= 0.25;
                    dz1  = s1 * exp_fact1;
                    z1   = zeta1 + dz1;
                    dz2  = dz1 * alpha * alpha * exp_fact2;
                    z2   = zeta2 + dz2;
                    dz12 = dz1 * alpha * rho * exp_fact12;
                    z12  = zeta12 + dz12;

                    pr1 = amortMidat_lgmamortsopval2F(
                        ncpn - j,
                        coupon + j,
                        cpn_df + j,
                        cpn_cvg + j,
                        cpn_G1 + j,
                        cpn_G2 + j,
                        z1,
                        z2,
                        z12,
                        ex_G1[i],
                        ex_G2[i]);

                    it++;
                }

                it = 1;
                while (pr2 < mkt_price[i] && it < 6)
                {
                    s2 *= 4.0;
                    dz1  = s2 * exp_fact1;
                    z1   = zeta1 + dz1;
                    dz2  = dz1 * alpha * alpha * exp_fact2;
                    z2   = zeta2 + dz2;
                    dz12 = dz1 * alpha * rho * exp_fact12;
                    z12  = zeta12 + dz12;

                    pr2 = amortMidat_lgmamortsopval2F(
                        ncpn - j,
                        coupon + j,
                        cpn_df + j,
                        cpn_cvg + j,
                        cpn_G1 + j,
                        cpn_G2 + j,
                        z1,
                        z2,
                        z12,
                        ex_G1[i],
                        ex_G2[i]);

                    it++;
                }
            }
            else
            /*	Next ones: limited vol variation and only 5 iterations */
            {
                niter  = NITER;
                minvar = quad_var / t * MAX_FACT;
                maxvar = 4.0 * quad_var / t / MAX_FACT;
            }

            d = 0.0;

            if (s1 < s2)
            {
                res_iter[0][0] = s1;
                res_iter[0][1] = pr1;

                res_iter[1][0] = s2;
                res_iter[1][1] = pr2;
            }
            else
            {
                res_iter[0][0] = s2;
                res_iter[0][1] = pr2;

                res_iter[1][0] = s1;
                res_iter[1][1] = pr1;
            }

            k = 2;

            for (it = 1; it < niter; it++)
            {
                temp = s2;

                /*	Newton iteration */

                //				s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 -
                //pr1);

                s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                /*	Out of lower bound */
                if (s2 < minvar)
                {
                    /*	Calibrate to market price + 1bp */
                    d  = CALPRESNOT;
                    s2 = temp;

                    /*
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    */

                    s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                    /*	Still out of bounds */
                    if (s2 < minvar)
                    {
                        s2 = minvar;
                    }
                }

                /*	Out of upper bound */
                if (s2 > maxvar)
                {
                    /*	Calibrate to market price - 1bp */
                    d  = -CALPRESNOT;
                    s2 = temp;

                    /*
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    */

                    s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);
                    /*	Out of lower bound */
                    if (s2 < minvar)
                    {
                        /*	Calibrate to market price + 1bp */
                        d  = CALPRESNOT;
                        s2 = temp;

                        /*
                        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                        */

                        s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                        /*	Still out of bounds */
                        if (s2 < minvar)
                        {
                            s2 = minvar;
                        }
                    }
                }

                if (fabs(s2 - temp) < 1.0e-16)
                    break;

                s1  = temp;
                pr1 = pr2;

                /*	Reprice with new vol */
                dz1 = s2 * exp_fact1;
                q1  = zeta1 + dz1;

                dz2 = dz1 * alpha * alpha * exp_fact2;
                q2  = zeta2 + dz2;

                dz12 = dz1 * alpha * rho * exp_fact12;
                q12  = zeta12 + dz12;

                pr2 = amortMidat_lgmamortsopval2F(
                    ncpn - j,
                    coupon + j,
                    cpn_df + j,
                    cpn_cvg + j,
                    cpn_G1 + j,
                    cpn_G2 + j,
                    q1,
                    q2,
                    q12,
                    ex_G1[i],
                    ex_G2[i]);

                k++;

                /* Save Res */
                l = 0;
                while (l < k - 1 && res_iter[l][1] < pr2)
                {
                    l++;
                }

                if (l < k - 1)
                {
                    for (m = k - 2; m >= l; m--)
                    {
                        res_iter[m + 1][0] = res_iter[m][0];
                        res_iter[m + 1][1] = res_iter[m][1];
                    }
                }

                res_iter[l][0] = s2;
                res_iter[l][1] = pr2;

                if (fabs((mkt_price[i] + d) - pr2) < CALPRESNOT)
                    break;
            }

            /*	If failed, keep last vol */
            if (fabs((mkt_price[i] + d) - pr2) > CALPRESNOT)
            {
                if (!USE_JUMPS)
                {
                    smessage("Failed to calibrate for exercise date %d", i + 1);
                    s2 = s_last;
                }
                else
                {
                    smessage("Jump too large for exercise date %d", i + 1);
                }

                dz1  = s2 * exp_fact1;
                dz2  = dz1 * alpha * alpha * exp_fact2;
                dz12 = dz1 * alpha * rho * exp_fact12;

                q1  = zeta1 + dz1;
                q2  = zeta2 + dz2;
                q12 = zeta12 + dz12;
            }

            /* Back populate the previous Zeta */
            l = i - 1;
            while (l >= 0 && !ex_bool[l])
            {
                exp_fact1  = ((exp(2 * l1 * ex_time[l]) - exp(2 * l1 * t)) / 2 / l1);
                ex_zeta[l] = zeta1 + s2 * exp_fact1;
                l--;
            }

            s_last = s2;
            quad_var += s2 * (ex_time[i] - t);
            ex_zeta[i] = zeta1 = q1;
            zeta2              = q2;
            zeta12             = q12;
            t                  = ex_time[i];
            firstcalibdone     = 1;
        }
        //		else
        //		{
        //			s_last = s2;
        //			quad_var += s2 * (ex_time[i] - t);
        //			ex_zeta[i] = zeta1 + s2 * exp_fact1;
        //			zeta2 = zeta2 + s2 * exp_fact1 * alpha * alpha * exp_fact2;
        //			zeta12 = zeta12 + s2 * exp_fact1 * alpha * rho * exp_fact12;
        //			t = ex_time[i];
        //		}
    }

    /* Back populate the previous Zeta */
    l = nex - 1;
    while (l >= 0 && !ex_bool[l])
    {
        exp_fact1  = ((exp(2 * l1 * ex_time[l]) - exp(2 * l1 * t)) / 2 / l1);
        ex_zeta[l] = zeta1 + s2 * exp_fact1;
        l--;
    }

    if (res_iter)
        free_dmatrix(res_iter, 0, 5 * NITER, 0, 1);

    return NULL;
}

/*	Calibrate zeta to amortized diagonal given G: 2F case */
Err AmortMidat_lgmcalibzeta2F_ts(
    int    ncpn,       /*	Total number of cash-flow dates */
    double cpn[],      /*	Discounted Cash-Flow */
    double cpn_time[], /*	Cash-Flow times */
    double cpn_df[],   /*	Df to cash-flow dates */
    double cpn_cvg[],  /*	cvg from i-1 to i */
    double cpn_G1[],   /*	G1 at cash-flow dates */
    double cpn_G2[],   /*	G2 at cash-flow dates */
    int    nex,        /*	Total number of exercise dates */
    int    ex_bool[],  /*	Exercise or not */
    double ex_time[],  /*	Exercise times */
    int    ex_cpn[],   /*	Index of the first cash-flow to be exercised */
    double ex_G1[],    /*	G1 at exercise date */
    double ex_G2[],    /*	G2 at exercise date */
    // double			strike[],							/*	Strikes
    // */
    double mkt_price[],     /*	Market prices */
    double ex_fee[],        /*	Exercise Fees for Diagonal */
    double ex_zeta[],       /*	Output: zetas (1) */
    double floatnotional[], /*	Float Notional	*/
    /*	Lambda, Alpha, gamma, rho */
    int    nlambda,
    double lambda_time[],
    double lambda[],
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    double s1, s2, s_last, ds;
    int    i, j, it;
    // double		l1 = lambda, l2 = lambda + gamma;
    // double		l1 , l2;
    double z1, z2, z12;
    double q1, q2, q12;
    double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
    double pr1, pr2;
    double exp_fact1, exp_fact2, exp_fact12;
    double temp;
    double minvar, maxvar, d;
    double quad_var;
    int    niter;
    int    firstcalibdone;
    double convert_ratio;
    int    k, l, m;

    static double coupon[MAX_CPN];

    double** res_iter = NULL;

    if (!lambda_time || !lambda)
    {
        return "Invalid lambda times/values ";
    }

    for (i = 0; i < ncpn; ++i)
    {
        coupon[i] = cpn[i];
    }

    i = 0;
    while (i < nex && !ex_bool[i])
    {
        i++;
    }

    if (i == nex)
    {
        return "No option to calibrate";
    }

    res_iter = dmatrix(0, 5 * NITER, 0, 1);

    if (!res_iter)
    {
        return "Memory allocation faillure in lgmprcapgiventauts_momentum_dlm";
    }

    convert_ratio = sqrt(
        lambda[0] / LAMBDA_REF * (1.0 - exp(-2.0 * LAMBDA_REF * cpn_time[i])) /
        (1.0 - exp(-2.0 * lambda[0] * cpn_time[i])));

    /*	Initial guess: total vol of 0.005 - 0.07 bps */
    /*	s = local variance */
    s1 = 0.02;
    s1 *= convert_ratio / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s1 *= s1;

    s2 = 0.04;
    s2 *= convert_ratio / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s2 *= s2;

    /*	Initialisation */
    quad_var = 0.0;
    s_last   = s1;
    ds       = s2 / s1;
    zeta1 = zeta2 = zeta12 = 0.0;
    t                      = 0.0;

    /*	If only one option, no skipping last */
    if (nex == 1)
    {
        skip_last = 0;
    }

    firstcalibdone = 0;
    for (i = 0; i < nex - (skip_last > 0); i++)
    {
        j = ex_cpn[i];

        exp_fact1 = export_lgmcalcexpfact_tauts(t, ex_time[i], nlambda, lambda_time, lambda, 0.0);
        exp_fact2 =
            export_lgmcalcexpfact_tauts(t, ex_time[i], nlambda, lambda_time, lambda, gamma) /
            exp_fact1;
        exp_fact12 =
            export_lgmcalcexpfact_tauts(t, ex_time[i], nlambda, lambda_time, lambda, 0.5 * gamma) /
            exp_fact1;

        if ((ex_bool[i] == 1) || (i == 0))
        {
            /*	Initial guess: last vol */
            s1 = s_last;
            /*	Initial guess 2: last vol * ds */
            s2 = ds * s1;

            /*	First price */
            dz1 = s1 * exp_fact1;
            z1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            z2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            z12  = zeta12 + dz12;

            coupon[j] = -floatnotional[j] * cpn_df[j] - ex_fee[i] * cpn_df[j];

            pr1 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                z1,
                z2,
                z12,
                ex_G1[i],
                ex_G2[i]);

            if (fabs(mkt_price[i] - pr1) < CALPRESNOT)
            {
                s_last = s1;
                quad_var += s1 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = z1;
                zeta2              = z2;
                zeta12             = z12;
                t                  = ex_time[i];
                continue;
            }

            /*	Second price */
            dz1 = s2 * exp_fact1;
            q1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            q2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            q12  = zeta12 + dz12;

            pr2 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                q1,
                q2,
                q12,
                ex_G1[i],
                ex_G2[i]);

            if (fabs(mkt_price[i] - pr2) < CALPRESNOT)
            {
                s_last = s2;
                quad_var += s2 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = q1;
                zeta2              = q2;
                zeta12             = q12;
                t                  = ex_time[i];
                continue;
            }

            /*	First option: no limits, 15 iterations */
            if (firstcalibdone == 0)
            {
                niter  = 3 * NITER;
                minvar = 1.0e-16;
                maxvar = 100.0;

                /* make sure we are between the bounds */
                it = 1;

                while (pr1 > mkt_price[i] && it < 6)
                {
                    s1 *= 0.25;
                    dz1  = s1 * exp_fact1;
                    z1   = zeta1 + dz1;
                    dz2  = dz1 * alpha * alpha * exp_fact2;
                    z2   = zeta2 + dz2;
                    dz12 = dz1 * alpha * rho * exp_fact12;
                    z12  = zeta12 + dz12;

                    pr1 = amortMidat_lgmamortsopval2F(
                        ncpn - j,
                        coupon + j,
                        cpn_df + j,
                        cpn_cvg + j,
                        cpn_G1 + j,
                        cpn_G2 + j,
                        z1,
                        z2,
                        z12,
                        ex_G1[i],
                        ex_G2[i]);

                    it++;
                }

                it = 1;
                while (pr2 < mkt_price[i] && it < 6)
                {
                    s2 *= 4.0;
                    dz1  = s2 * exp_fact1;
                    z1   = zeta1 + dz1;
                    dz2  = dz1 * alpha * alpha * exp_fact2;
                    z2   = zeta2 + dz2;
                    dz12 = dz1 * alpha * rho * exp_fact12;
                    z12  = zeta12 + dz12;

                    pr2 = amortMidat_lgmamortsopval2F(
                        ncpn - j,
                        coupon + j,
                        cpn_df + j,
                        cpn_cvg + j,
                        cpn_G1 + j,
                        cpn_G2 + j,
                        z1,
                        z2,
                        z12,
                        ex_G1[i],
                        ex_G2[i]);

                    it++;
                }
            }
            else
            /*	Next ones: limited vol variation and only 5 iterations */
            {
                niter  = NITER;
                minvar = quad_var / t * MAX_FACT;
                maxvar = 4.0 * quad_var / t / MAX_FACT;
            }

            d = 0.0;

            if (s1 < s2)
            {
                res_iter[0][0] = s1;
                res_iter[0][1] = pr1;

                res_iter[1][0] = s2;
                res_iter[1][1] = pr2;
            }
            else
            {
                res_iter[0][0] = s2;
                res_iter[0][1] = pr2;

                res_iter[1][0] = s1;
                res_iter[1][1] = pr1;
            }

            k = 2;

            for (it = 1; it < niter; it++)
            {
                temp = s2;

                /*	Newton iteration */

                //				s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 -
                //pr1);

                s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                /*	Out of lower bound */
                if (s2 < minvar)
                {
                    /*	Calibrate to market price + 1bp */
                    d  = CALPRESNOT;
                    s2 = temp;

                    /*
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    */

                    s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                    /*	Still out of bounds */
                    if (s2 < minvar)
                    {
                        s2 = minvar;
                    }
                }

                /*	Out of upper bound */
                if (s2 > maxvar)
                {
                    /*	Calibrate to market price - 1bp */
                    d  = -CALPRESNOT;
                    s2 = temp;

                    /*
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    */

                    s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);
                    /*	Out of lower bound */
                    if (s2 < minvar)
                    {
                        /*	Calibrate to market price + 1bp */
                        d  = CALPRESNOT;
                        s2 = temp;

                        /*
                        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                        */

                        s2 = solve_for_next_coef_dlm(res_iter, k, mkt_price[i] + d, 2);

                        /*	Still out of bounds */
                        if (s2 < minvar)
                        {
                            s2 = minvar;
                        }
                    }
                }

                if (fabs(s2 - temp) < 1.0e-16)
                    break;

                s1  = temp;
                pr1 = pr2;

                /*	Reprice with new vol */
                dz1 = s2 * exp_fact1;
                q1  = zeta1 + dz1;

                dz2 = dz1 * alpha * alpha * exp_fact2;
                q2  = zeta2 + dz2;

                dz12 = dz1 * alpha * rho * exp_fact12;
                q12  = zeta12 + dz12;

                pr2 = amortMidat_lgmamortsopval2F(
                    ncpn - j,
                    coupon + j,
                    cpn_df + j,
                    cpn_cvg + j,
                    cpn_G1 + j,
                    cpn_G2 + j,
                    q1,
                    q2,
                    q12,
                    ex_G1[i],
                    ex_G2[i]);

                k++;

                /* Save Res */
                l = 0;
                while (l < k - 1 && res_iter[l][1] < pr2)
                {
                    l++;
                }

                if (l < k - 1)
                {
                    for (m = k - 2; m >= l; m--)
                    {
                        res_iter[m + 1][0] = res_iter[m][0];
                        res_iter[m + 1][1] = res_iter[m][1];
                    }
                }

                res_iter[l][0] = s2;
                res_iter[l][1] = pr2;

                if (fabs((mkt_price[i] + d) - pr2) < CALPRESNOT)
                    break;
            }

            /*	If failed, keep last vol */
            if (fabs((mkt_price[i] + d) - pr2) > CALPRESNOT)
            {
                if (!USE_JUMPS)
                {
                    smessage("Failed to calibrate for exercise date %d", i + 1);
                    s2 = s_last;
                }
                else
                {
                    smessage("Jump too large for exercise date %d", i + 1);
                }

                dz1  = s2 * exp_fact1;
                dz2  = dz1 * alpha * alpha * exp_fact2;
                dz12 = dz1 * alpha * rho * exp_fact12;

                q1  = zeta1 + dz1;
                q2  = zeta2 + dz2;
                q12 = zeta12 + dz12;
            }

            s_last = s2;
            quad_var += s2 * (ex_time[i] - t);
            ex_zeta[i] = zeta1 = q1;
            zeta2              = q2;
            zeta12             = q12;
            t                  = ex_time[i];
            firstcalibdone     = 1;
        }
        //		else
        //		{
        //			s_last = s2;
        //			quad_var += s2 * (ex_time[i] - t);
        //			ex_zeta[i] = zeta1 + s2 * exp_fact1;
        //			zeta2 = zeta2 + s2 * exp_fact1 * alpha * alpha * exp_fact2;
        //			zeta12 = zeta12 + s2 * exp_fact1 * alpha * rho * exp_fact12;
        //			t = ex_time[i];
        //		}
    }

    /*	Replace last vol by last-1 if relevant */
    if (skip_last)
    {
        exp_fact1 =
            export_lgmcalcexpfact_tauts(t, ex_time[nex - 1], nlambda, lambda_time, lambda, 0.0);

        dz1 = s_last * exp_fact1;

        // dz1 = s_last * (( exp (2 * l1 * ex_time[nex-1]) - exp (2 * l1 * t)) / 2 / l1);
        ex_zeta[nex - 1] = zeta1 + dz1;
    }

    if (res_iter)
        free_dmatrix(res_iter, 0, 5 * NITER, 0, 1);

    return NULL;
}

/*	Calibrate zeta to amortized diagonal given G: 2F case */
Err zcAmortMidat_lgmcalibzeta2F(
    int    ncpn,          /*	Total number of cash-flow dates */
    double cpn[],         /*	Discounted Cash-Flow */
    double zc_last_cpn[], /*	Discounted Last Fix CF */
    double cpn_time[],    /*	Cash-Flow times */
    double cpn_df[],      /*	Df to cash-flow dates */
    double cpn_cvg[],     /*	cvg from i-1 to i */
    double cpn_G1[],      /*	G1 at cash-flow dates */
    double cpn_G2[],      /*	G2 at cash-flow dates */
    int    nex,           /*	Total number of exercise dates */
    int    ex_bool[],     /*	Exercise or not */
    double ex_time[],     /*	Exercise times */
    int    ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    double ex_G1[],       /*	G1 at exercise date */
    double ex_G2[],       /*	G2 at exercise date */
    // double			strike[],							/*	Strikes
    // */
    double mkt_price[],     /*	Market prices */
    double ex_zeta[],       /*	Output: zetas (1) */
    double floatnotional[], /*	Float Notional	*/
    /*	Lambda, Alpha, gamma, rho */
    double lambda,
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    double s1, s2, s_last, ds;
    int    i, j, it;
    double l1 = lambda, l2 = lambda + gamma;
    double z1, z2, z12;
    double q1, q2, q12;
    double t, zeta1, zeta2, zeta12, dz1, dz2, dz12;
    double pr1, pr2;
    double exp_fact1, exp_fact2, exp_fact12;
    double temp;
    double minvar, maxvar, d;
    double quad_var;
    int    niter;
    int    firstcalibdone;
    double minvol, maxvol, ds_long;

    static double coupon[MAX_CPN];

    for (i = 0; i < ncpn; ++i)
    {
        coupon[i] = cpn[i];
    }

    /*	Initial guess: total vol of 85 - 115 bps */
    /*	s = local variance of first factor */
    minvol = 0.0085;
    maxvol = 0.03;

    s1 = minvol / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s1 *= s1;
    s2 = maxvol / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s2 *= s2;
    ds_long = s2 / s1;

    if (ex_time[0] < 1.0 / 12.0)
    {
        minvol = 0.0085;
        maxvol = 0.05;
    }

    s1 = minvol / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s1 *= s1;
    s2 = maxvol / sqrt(1.0 + alpha * alpha + 2.0 * rho * alpha);
    s2 *= s2;

    /*	Initialisation */
    quad_var = 0.0;
    s_last   = s1;
    ds       = s2 / s1;
    zeta1 = zeta2 = zeta12 = 0.0;
    t                      = 0.0;

    /*	If only one option, no skipping last */
    if (nex == 1)
    {
        skip_last = 0;
    }

    firstcalibdone = 0;
    for (i = 0; i < nex - (skip_last > 0); i++)
    {
        if (ex_time[i] > 0.25)
        {
            ds = ds_long;
        }

        j         = ex_cpn[i];
        exp_fact1 = ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
        exp_fact2 = ((exp(2 * l2 * ex_time[i]) - exp(2 * l2 * t)) / 2 / l2) /
                    ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);
        exp_fact12 = ((exp((l1 + l2) * ex_time[i]) - exp((l1 + l2) * t)) / (l1 + l2)) /
                     ((exp(2 * l1 * ex_time[i]) - exp(2 * l1 * t)) / 2 / l1);

        if ((ex_bool[i] == 1) || (i == 0))
        {
            /*	Initial guess: last vol */
            s1 = s_last;
            /*	Initial guess 2: last vol * ds */
            s2 = ds * s1;

            /*	First price */
            dz1 = s1 * exp_fact1;
            z1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            z2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            z12  = zeta12 + dz12;

            coupon[j] = -floatnotional[j] * cpn_df[j];
            coupon[ncpn - 1] += zc_last_cpn[i];

            pr1 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                z1,
                z2,
                z12,
                ex_G1[i],
                ex_G2[i]);
            //				strike[i]);

            if (fabs(mkt_price[i] - pr1) < CALPRESNOT)
            {
                s_last = s1;
                quad_var += s1 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = z1;
                zeta2              = z2;
                zeta12             = z12;
                t                  = ex_time[i];
                continue;
            }

            /*	Second price */
            dz1 = s2 * exp_fact1;
            q1  = zeta1 + dz1;

            dz2 = dz1 * alpha * alpha * exp_fact2;
            q2  = zeta2 + dz2;

            dz12 = dz1 * alpha * rho * exp_fact12;
            q12  = zeta12 + dz12;

            pr2 = amortMidat_lgmamortsopval2F(
                ncpn - j,
                coupon + j,
                cpn_df + j,
                cpn_cvg + j,
                cpn_G1 + j,
                cpn_G2 + j,
                q1,
                q2,
                q12,
                ex_G1[i],
                ex_G2[i]);

            if (fabs(mkt_price[i] - pr2) < CALPRESNOT)
            {
                s_last = s2;
                quad_var += s2 * (ex_time[i] - t);
                ex_zeta[i] = zeta1 = q1;
                zeta2              = q2;
                zeta12             = q12;
                t                  = ex_time[i];
                continue;
            }

            /*	First option: no limits, 15 iterations */
            if (firstcalibdone == 0)
            {
                niter  = 3 * NITER;
                minvar = 1.0e-16;
                maxvar = 100.0;
            }
            else
            /*	Next ones: limited vol variation and only 5 iterations */
            {
                niter  = NITER;
                minvar = quad_var / t * MAX_FACT;
                maxvar = 4.0 * quad_var / t / MAX_FACT;
            }

            d = 0.0;

            for (it = 1; it < niter; it++)
            {
                temp = s2;

                /*	Newton iteration */
                s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

                /*	Out of lower bound */
                if (s2 < minvar)
                {
                    /*	Calibrate to market price + 1bp */
                    d  = CALPRESNOT;
                    s2 = temp;
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    /*	Still out of bounds */
                    if (s2 < minvar)
                    {
                        s2 = minvar;
                    }
                }

                /*	Out of upper bound */
                if (s2 > maxvar)
                {
                    /*	Calibrate to market price - 1bp */
                    d  = -CALPRESNOT;
                    s2 = temp;
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    /*	Still out of bounds */
                    if (s2 > maxvar)
                    {
                        s2 = maxvar;
                    }
                }

                if (fabs(s2 - temp) < 1.0e-16)
                    break;

                s1  = temp;
                pr1 = pr2;

                /*	Reprice with new vol */
                dz1 = s2 * exp_fact1;
                q1  = zeta1 + dz1;

                dz2 = dz1 * alpha * alpha * exp_fact2;
                q2  = zeta2 + dz2;

                dz12 = dz1 * alpha * rho * exp_fact12;
                q12  = zeta12 + dz12;

                pr2 = amortMidat_lgmamortsopval2F(
                    ncpn - j,
                    coupon + j,
                    cpn_df + j,
                    cpn_cvg + j,
                    cpn_G1 + j,
                    cpn_G2 + j,
                    q1,
                    q2,
                    q12,
                    ex_G1[i],
                    ex_G2[i]);

                if (fabs((mkt_price[i] + d) - pr2) < CALPRESNOT)
                    break;
            }

            /*	If failed, keep last vol */
            if (fabs((mkt_price[i] + d) - pr2) > CALPRESNOT)
            {
                if (s2 < 0)
                {
                    smessage("Negative Variance for exercise date %d", i + 1);
                }
                else
                {
                    smessage("Jump too large for exercise date %d", i + 1);
                }

                s2 = s_last;

                dz1  = s2 * exp_fact1;
                dz2  = dz1 * alpha * alpha * exp_fact2;
                dz12 = dz1 * alpha * rho * exp_fact12;

                q1  = zeta1 + dz1;
                q2  = zeta2 + dz2;
                q12 = zeta12 + dz12;
            }

            s_last = s2;
            quad_var += s2 * (ex_time[i] - t);
            ex_zeta[i] = zeta1 = q1;
            zeta2              = q2;
            zeta12             = q12;
            t                  = ex_time[i];
            firstcalibdone     = 1;

            coupon[ncpn - 1] -= zc_last_cpn[i];
        }
        else
        {
            //			s_last = s2;
            //			quad_var += s2 * (ex_time[i] - t);
            //			ex_zeta[i] = zeta1 + s2 * exp_fact1;
            //			zeta2 = zeta2 + s2 * exp_fact1 * alpha * alpha * exp_fact2;
            //			zeta12 = zeta12 + s2 * exp_fact1 * alpha * rho * exp_fact12;
            //			t = ex_time[i];
        }
    }

    /*	Replace last vol by last-1 if relevant */
    if (skip_last)
    {
        dz1              = s_last * ((exp(2 * l1 * ex_time[nex - 1]) - exp(2 * l1 * t)) / 2 / l1);
        ex_zeta[nex - 1] = zeta1 + dz1;
    }

    return NULL;
}

/*	Calibrate zeta to diagonal given G: 1F case */
Err AmortMidat_lgmcalibzeta1F(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flow */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    double cpn_G[],         /*	G at cash-flow dates */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    double ex_G[],          /*	G at exercise date */
    double mkt_price[],     /*	Market prices */
    double ex_fee[],        /*	Exercise Fees for Diagonal */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double lambda,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    double        s1, s2, s_last, ds;
    double        t, zeta;
    int           i, j, nb_iter, pos, m, it;
    double        dz1, dz2, z1, z2;
    double        pr1, pr2;
    double        exp_fact;
    double        temp;
    double        minvar, maxvar, d;
    double        quad_var;
    double        convert_ratio;
    int           niter;
    int           skip_calib;
    int           firstcalibdone;
    static double coupon[MAX_CPN];

    double** res_iter = NULL;

    for (i = 0; i < ncpn; ++i)
    {
        coupon[i] = cpn[i];
    }

    i = 0;
    while (i < nex && !ex_bool[i])
    {
        i++;
    }

    if (i == nex)
    {
        return "No option to calibrate";
    }

    res_iter = dmatrix(0, 5 * NITER, 0, 1);

    if (!res_iter)
    {
        return "Memory allocation faillure in lgmprcapgiventauts_momentum_dlm";
    }

    convert_ratio = sqrt(
        lambda / LAMBDA_REF * (1.0 - exp(-2.0 * LAMBDA_REF * cpn_time[i])) /
        (1.0 - exp(-2.0 * lambda * cpn_time[i])));

    /*	Initial guess: total vol of 0.005 - 0.07 bps */
    /*	s = local variance */
    s1 = 0.02;
    s1 *= convert_ratio;
    s1 *= s1;

    s2 = 0.03;
    s2 *= convert_ratio;
    s2 *= s2;

    /*	Initialisation */
    quad_var = 0.0;
    s_last   = s1;
    //	ds = s2 / s1;

    ds   = 0.75;
    zeta = 0.0;
    t    = 0.0;

    /*	If only one option, no skipping last */
    if (nex == 1)
    {
        skip_last = 0;
    }
    firstcalibdone = 0;
    for (i = 0; i < nex - (skip_last > 0); i++)
    {
        j        = ex_cpn[i];
        exp_fact = ((exp(2 * lambda * ex_time[i]) - exp(2 * lambda * t)) / 2 / lambda);

        if ((ex_bool[i] == 1) || (i == 0))
        {
            /*	Initial guess: last vol */
            s1 = s_last;
            /*	Initial guess 2: last vol * ds */
            if (i > 0)
            {
                s2 = ds * s1;
            }

            /*	First price */
            dz1 = s1 * exp_fact;
            z1  = zeta + dz1;

            coupon[j] = -floatnotional[j] * cpn_df[j] - ex_fee[i] * cpn_df[j];

            pr1 = amortMidat_lgmamortsopval1F(
                ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z1, ex_G[i]);

            if (fabs(mkt_price[i] - pr1) < CALPRESNOT)
            {
                s_last = s1;
                quad_var += s1 * (ex_time[i] - t);
                ex_zeta[i] = zeta = z1;
                t                 = ex_time[i];
                continue;
            }

            /*	Second price */
            dz2 = s2 * exp_fact;
            z2  = zeta + dz2;

            pr2 = amortMidat_lgmamortsopval1F(
                ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

            if (fabs(mkt_price[i] - pr2) < CALPRESNOT)
            {
                s_last = s2;
                quad_var += s2 * (ex_time[i] - t);
                ex_zeta[i] = zeta = z2;
                t                 = ex_time[i];
                continue;
            }

            /*	First option: no limits, 15 iterations */
            if (firstcalibdone == 0)
            {
                skip_calib = 0;
                niter      = 3 * NITER;
                minvar     = 1.0e-16;
                maxvar     = 100.0;

                /* make sure we are between the bounds */
                it = 1;

                while (pr1 > mkt_price[i] && it < 6)
                {
                    s2  = s1;
                    pr2 = pr1;

                    s1 *= 0.25;
                    dz1 = s1 * exp_fact;
                    z1  = zeta + dz1;

                    pr1 = amortMidat_lgmamortsopval1F(
                        ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z1, ex_G[i]);

                    it++;
                }

                it = 1;
                while (pr2 < mkt_price[i] && it < 10 && skip_calib == 0)
                {
                    s1  = s2;
                    pr1 = pr2;

                    s2 *= 1.5;
                    dz2 = s2 * exp_fact;
                    z2  = zeta + dz2;

                    pr2 = amortMidat_lgmamortsopval1F(
                        ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

                    if (pr2 < pr1)
                    {
                        /* we reached the limit */
                        skip_calib = 1;
                        s2         = s1;
                        pr2        = pr1;
                    }

                    it++;
                }
            }
            else
            /*	Next ones: limited vol variation and only 5 iterations */
            {
                niter  = NITER;
                minvar = quad_var / t * MAX_FACT;
                maxvar = 4.0 * quad_var / t / MAX_FACT;
            }

            d = 0.0;

            if (s1 < s2)
            {
                res_iter[0][0] = s1;
                res_iter[0][1] = pr1;

                res_iter[1][0] = s2;
                res_iter[1][1] = pr2;
            }
            else
            {
                res_iter[0][0] = s2;
                res_iter[0][1] = pr2;

                res_iter[1][0] = s1;
                res_iter[1][1] = pr1;
            }

            nb_iter = 2;

            if (skip_calib == 0)
            {
                for (it = 1; it < niter; it++)
                {
                    temp = s2;

                    /*	Newton iteration */

                    s2 = solve_for_next_coef_dlm(res_iter, nb_iter, mkt_price[i] + d, 2);

                    /*	Out of lower bound */
                    if (s2 < minvar)
                    {
                        /*	Calibrate to market price + 1bp */
                        d  = CALPRESNOT;
                        s2 = temp;

                        /*
                        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                        */

                        s2 = solve_for_next_coef_dlm(res_iter, nb_iter, mkt_price[i] + d, 2);

                        /*	Still out of bounds */
                        if (s2 < minvar)
                        {
                            s2 = minvar;
                        }
                    }

                    /*	Out of upper bound */
                    if (s2 > maxvar)
                    {
                        /*	Calibrate to market price - 1bp */
                        d  = -CALPRESNOT;
                        s2 = temp;

                        /*
                        s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                        */

                        s2 = solve_for_next_coef_dlm(res_iter, nb_iter, mkt_price[i] + d, 2);
                        /*	Out of lower bound */
                        if (s2 < minvar)
                        {
                            /*	Calibrate to market price + 1bp */
                            d  = CALPRESNOT;
                            s2 = temp;

                            /*
                            s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                            */

                            s2 = solve_for_next_coef_dlm(res_iter, nb_iter, mkt_price[i] + d, 2);

                            /*	Still out of bounds */
                            if (s2 < minvar)
                            {
                                s2 = minvar;
                            }
                        }
                    }

                    if (fabs(s2 - temp) < 1.0e-16)
                        break;

                    s1  = temp;
                    pr1 = pr2;

                    /*	Reprice with new vol */
                    dz2 = s2 * exp_fact;
                    z2  = zeta + dz2;

                    pr2 = amortMidat_lgmamortsopval1F(
                        ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

                    nb_iter++;

                    /* Save Res */
                    pos = 0;
                    while (pos < nb_iter - 1 && res_iter[pos][1] < pr2)
                    {
                        pos++;
                    }

                    if (pos < nb_iter - 1)
                    {
                        for (m = nb_iter - 2; m >= pos; m--)
                        {
                            res_iter[m + 1][0] = res_iter[m][0];
                            res_iter[m + 1][1] = res_iter[m][1];
                        }
                    }

                    res_iter[pos][0] = s2;
                    res_iter[pos][1] = pr2;

                    if (fabs((mkt_price[i] + d) - pr2) < CALPRESNOT)
                        break;
                }
            }

            /*	If failed, keep last vol */
            if (fabs((mkt_price[i] + d) - pr2) > CALPRESNOT)
            {
                if (!USE_JUMPS)
                {
                    smessage("Failed to calibrate for exercise date %d", i + 1);
                    s2 = s_last;
                }
                else
                {
                    smessage("Jump too large for exercise date %d", i + 1);
                }

                dz2 = s2 * exp_fact;
                z2  = zeta + dz2;
            }

            s_last = s2;
            quad_var += s2 * (ex_time[i] - t);
            ex_zeta[i] = zeta = z2;
            t                 = ex_time[i];
            firstcalibdone    = 1;
            skip_calib        = 0;
        }
        //		else
        //		{
        //			s_last = s2;
        //			quad_var += s2 * (ex_time[i] - t);
        //			temp_zeta = zeta;
        //			zeta = temp_zeta + s2 * exp_fact;
        //			ex_zeta[i] = temp_zeta + s2 * exp_fact;
        //			t = ex_time[i];
        //		}
    }

    /*	Replace last vol by last-1 if relevant */
    if (skip_last)
    {
        dz2 = s_last * ((exp(2 * lambda * ex_time[nex - 1]) - exp(2 * lambda * t)) / 2 / lambda);
        ex_zeta[nex - 1] = zeta + dz2;
    }

    if (res_iter)
        free_dmatrix(res_iter, 0, 5 * NITER, 0, 1);

    return NULL;
}

/*	Calibrate zeta to diagonal given G: 1F case */
Err zcAmortMidat_lgmcalibzeta1F(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flow */
    double zc_last_cpn[],   /*	Discounted Last Fix CF */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    double cpn_G[],         /*	G at cash-flow dates */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    double ex_G[],          /*	G at exercise date */
    double mkt_price[],     /*	Market prices */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double lambda,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    double        s1, s2, s_last, ds;
    double        t, zeta;
    int           i, j, it;
    double        dz1, dz2, z1, z2;
    double        pr1, pr2;
    double        exp_fact;
    double        temp;
    double        minvar, maxvar, d;
    double        quad_var;
    int           niter;
    int           firstcalibdone;
    static double coupon[MAX_CPN];

    for (i = 0; i < ncpn; ++i)
    {
        coupon[i] = cpn[i];
    }

    /*	Initial guess: total vol of 85 - 115 bps */
    /*	s = local variance */
    s1 = 0.0085;
    s1 *= s1;
    s2 = 0.0115;
    s2 *= s2;

    /*	Initialisation */
    quad_var = 0.0;
    s_last   = s1;
    ds       = s2 / s1;
    zeta     = 0.0;
    t        = 0.0;

    /*	If only one option, no skipping last */
    if (nex == 1)
    {
        skip_last = 0;
    }
    firstcalibdone = 0;
    for (i = 0; i < nex - (skip_last > 0); i++)
    {
        j        = ex_cpn[i];
        exp_fact = ((exp(2 * lambda * ex_time[i]) - exp(2 * lambda * t)) / 2 / lambda);

        if ((ex_bool[i] == 1) || (i == 0))
        {
            /*	Initial guess: last vol */
            s1 = s_last;
            /*	Initial guess 2: last vol * ds */
            s2 = ds * s1;

            /*	First price */
            dz1 = s1 * exp_fact;
            z1  = zeta + dz1;

            coupon[j] = -floatnotional[j] * cpn_df[j];
            coupon[ncpn - 1] += zc_last_cpn[i];

            pr1 = amortMidat_lgmamortsopval1F(
                ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z1, ex_G[i]);

            if (fabs(mkt_price[i] - pr1) < CALPRESNOT)
            {
                s_last = s1;
                quad_var += s1 * (ex_time[i] - t);
                ex_zeta[i] = zeta = z1;
                t                 = ex_time[i];
                continue;
            }

            /*	Second price */
            dz2 = s2 * exp_fact;
            z2  = zeta + dz2;

            pr2 = amortMidat_lgmamortsopval1F(
                ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

            if (fabs(mkt_price[i] - pr2) < CALPRESNOT)
            {
                s_last = s2;
                quad_var += s2 * (ex_time[i] - t);
                ex_zeta[i] = zeta = z2;
                t                 = ex_time[i];
                continue;
            }

            /*	First option: no limits, 15 iterations */
            if (firstcalibdone == 0)
            {
                niter  = 3 * NITER;
                minvar = 1.0e-16;
                maxvar = 100.0;
            }
            else
            /*	Next ones: limited vol variation and only 5 iterations */
            {
                niter  = NITER;
                minvar = quad_var / t * MAX_FACT;
                maxvar = 4.0 * quad_var / t / MAX_FACT;
            }

            d = 0.0;

            for (it = 1; it < niter; it++)
            {
                temp = s2;

                /*	Newton iteration */
                s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);

                /*	Out of lower bound */
                if (s2 < minvar)
                {
                    /*	Calibrate to market price + 1bp */
                    d  = CALPRESNOT;
                    s2 = temp;
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    /*	Still out of bounds */
                    if (s2 < minvar)
                    {
                        s2 = minvar;
                    }
                }

                /*	Out of upper bound */
                if (s2 > maxvar)
                {
                    /*	Calibrate to market price - 1bp */
                    d  = -CALPRESNOT;
                    s2 = temp;
                    s2 -= (pr2 - (mkt_price[i] + d)) * (s2 - s1) / (pr2 - pr1);
                    /*	Still out of bounds */
                    if (s2 > maxvar)
                    {
                        s2 = maxvar;
                    }
                }

                if (fabs(s2 - temp) < 1.0e-16)
                    break;

                s1  = temp;
                pr1 = pr2;

                /*	Reprice with new vol */
                dz2 = s2 * exp_fact;
                z2  = zeta + dz2;

                pr2 = amortMidat_lgmamortsopval1F(
                    ncpn - j, coupon + j, cpn_df + j, cpn_cvg + j, cpn_G + j, z2, ex_G[i]);

                if (fabs((mkt_price[i] + d) - pr2) < CALPRESNOT)
                    break;
            }

            /*	If failed, keep last vol */
            if (fabs((mkt_price[i] + d) - pr2) > CALPRESNOT)
            {
                if (s2 < 0)
                {
                    smessage("Negative Variance for exercise date %d", i + 1);
                }
                else
                {
                    smessage("Jump too large for exercise date %d", i + 1);
                }
                //				smessage ("Failed to calibrate for exercise date %d",
                //i+1);
                s2  = s_last;
                dz2 = s2 * exp_fact;
                z2  = zeta + dz2;
            }

            s_last = s2;
            quad_var += s2 * (ex_time[i] - t);
            ex_zeta[i] = zeta = z2;
            t                 = ex_time[i];
            firstcalibdone    = 1;

            coupon[ncpn - 1] -= zc_last_cpn[i];
        }
        else
        {
            //			s_last = s2;
            //			quad_var += s2 * (ex_time[i] - t);
            //			temp_zeta = zeta;
            //			zeta = temp_zeta + s2 * exp_fact;
            //			ex_zeta[i] = temp_zeta + s2 * exp_fact;
            //			t = ex_time[i];
        }
    }

    /*	Replace last vol by last-1 if relevant */
    if (skip_last)
    {
        dz2 = s_last * ((exp(2 * lambda * ex_time[nex - 1]) - exp(2 * lambda * t)) / 2 / lambda);
        ex_zeta[nex - 1] = zeta + dz2;
    }

    return NULL;
}

/*	Calibrate and price cap given lambda */
Err AmortMidat_lgmprcapgivenlambda(
    int    ncpn,          /*	Total number of cash-flow dates */
    double cpn[],         /*	Discounted Cash-Flows */
    double cpn_time[],    /*	Cash-Flow times */
    double cpn_df[],      /*	Df to cash-flow dates */
    double cpn_cvg[],     /*	cvg from i-1 to i */
    int    nex,           /*	Total number of exercise dates */
    int    ex_bool[],     /*	Exercise or not */
    double ex_time[],     /*	Exercise times */
    int    ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],    /*	Number of coupons in each caplet */
    double ex_lprice[],   /*	Market prices for diagonal */
    double ex_fee[],      /*	Exercise Fee for diagonal */
    double ex_sstrike[],  /*	Strikes for cap */
    double floatcoupon[], /*	Float Coupon */
    // double			swaption_strikes[],					/*	Strikes for standard swaptions
    // */
    double swaption_cpn[], /*	Discounted Cash-Flows */
    // double			cpn_cvg_standard[],					/*	standard cvg from i-1 to i
    // */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double lambda,          /*	Lambda */
    int    one2F,           /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last,  /*	If 1, the last option is disregarded
                                               and the forward volatility is flat from option
                                               n-1 */
    int     price_cap, /*	0: just calibrate */
    double* ex_sprice) /*	Cap price as output */
{
    static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN], zeta2[MAX_CPN],
        zeta12[MAX_CPN];
    Err err;

    /*	Setup G function */
    export_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

    if (one2F == 2)
    {
        export_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time, ex_G2);
    }

    if (one2F == 1)
    {
        err = AmortMidat_lgmcalibzeta1F(
            ncpn,
            cpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            cpn_G,
            nex,
            ex_bool,
            ex_time,
            ex_cpn,
            ex_G,
            ex_lprice,
            ex_fee,
            ex_zeta,
            floatnotional,
            lambda,
            skip_last);

        if (err)
        {
            return err;
        }

        if (price_cap == 1)
        {
            *ex_sprice = amortMidat_lgmcapval1F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                ex_G);
        }
        else if (price_cap == 2)
        {
            *ex_sprice = amortMidat_lgmsumsopval1F_b(
                ncpn, nex, ex_cpn, cpn_df, cpn_G, ex_zeta, ex_G, swaption_cpn);
        }
    }
    else
    {
        err = AmortMidat_lgmcalibzeta2F(
            ncpn,
            cpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            cpn_G,
            cpn_G2,
            nex,
            ex_bool,
            ex_time,
            ex_cpn,
            ex_G,
            ex_G2,
            ex_lprice,
            ex_fee,
            ex_zeta,
            floatnotional,
            lambda,
            alpha,
            gamma,
            rho,
            skip_last);

        if (err)
        {
            return err;
        }

        export_lgmcalczeta2zeta12(nex, ex_time, ex_zeta, lambda, alpha, gamma, rho, zeta2, zeta12);

        if (price_cap == 1)
        {
            *ex_sprice = amortMidat_lgmcapval2F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                cpn_G2,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                zeta2,
                zeta12,
                ex_G,
                ex_G2);
        }
        else if (price_cap == 2)
        {
            *ex_sprice = amortMidat_lgmsumsopval2F_b(
                ncpn,
                nex,
                ex_cpn,
                cpn_df,
                cpn_G,
                cpn_G2,
                ex_zeta,
                zeta2,
                zeta12,
                ex_G,
                ex_G2,
                swaption_cpn);
        }
    }

    return NULL;
}

Err AmortMidat_lgmprcapgivenlambda_ts(
    int    ncpn,          /*	Total number of cash-flow dates */
    double cpn[],         /*	Discounted Cash-Flows */
    double cpn_time[],    /*	Cash-Flow times */
    double cpn_df[],      /*	Df to cash-flow dates */
    double cpn_cvg[],     /*	cvg from i-1 to i */
    int    nex,           /*	Total number of exercise dates */
    int    ex_bool[],     /*	Exercise or not */
    double ex_time[],     /*	Exercise times */
    int    ex_cpn[],      /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],    /*	Number of coupons in each caplet */
    double ex_lprice[],   /*	Market prices for diagonal */
    double ex_fee[],      /*	Exercise Fee for diagonal */
    double ex_sstrike[],  /*	Strikes for cap */
    double floatcoupon[], /*	Float Coupon */
    // double			swaption_strikes[],					/*	Strikes for standard swaptions
    // */
    double swaption_cpn[], /*	Discounted Cash-Flows */
    // double			cpn_cvg_standard[],					/*	standard cvg from i-1 to i
    // */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/

    int     nNumLambdas,
    double* pdLambdaValue,
    double* pdLambdaTime,

    int one2F, /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last,  /*	If 1, the last option is disregarded
                                               and the forward volatility is flat from option
                                               n-1 */
    int     price_cap, /*	0: just calibrate */
    double* ex_sprice) /*	Cap price as output */
{
    static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN], zeta2[MAX_CPN],
        zeta12[MAX_CPN];
    Err err;

    //// does not allow 1 factor at this point
    if (one2F == 1)
    {
        return serror(
            "AmortMidat_lgmprcapgivenlambda_ts(..) :: only 2 factors are allowed at the moment !");
    }

    /*	Setup G function */
    export_lgmsetupG_ts(
        nNumLambdas, pdLambdaTime, pdLambdaValue, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

    if (one2F == 2)
    {
        export_lgmsetupG2_ts(
            nNumLambdas,
            pdLambdaTime,
            pdLambdaValue,
            gamma,
            ncpn,
            cpn_time,
            cpn_G2,
            nex,
            ex_time,
            ex_G2);
    }

    if (one2F == 1)
    {
#if 0
		err = AmortMidat_lgmcalibzeta1F(
			ncpn,
			cpn,
			cpn_time,
			cpn_df,
			cpn_cvg,
			cpn_G,
			nex,
			ex_bool,
			ex_time,
			ex_cpn,
			ex_G,
			ex_lprice,
			ex_fee,
			ex_zeta,
			floatnotional,
			lambda,
			skip_last);


		if (err)
		{
			return err;
		}
#endif

        if (price_cap == 1)
        {
            *ex_sprice = amortMidat_lgmcapval1F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                ex_G);
        }
        else if (price_cap == 2)
        {
            *ex_sprice = amortMidat_lgmsumsopval1F_b(
                ncpn, nex, ex_cpn, cpn_df, cpn_G, ex_zeta, ex_G, swaption_cpn);
        }
    }
    else
    {
        err = AmortMidat_lgmcalibzeta2F_ts(
            ncpn,
            cpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            cpn_G,
            cpn_G2,
            nex,
            ex_bool,
            ex_time,
            ex_cpn,
            ex_G,
            ex_G2,
            ex_lprice,
            ex_fee,
            ex_zeta,
            floatnotional,
            nNumLambdas,
            pdLambdaTime,
            pdLambdaValue,
            alpha,
            gamma,
            rho,
            skip_last);

        if (err)
        {
            return err;
        }

        export_lgmcalczeta2zeta12_ts(
            nex,
            ex_time,
            ex_zeta,
            nNumLambdas,
            pdLambdaTime,
            pdLambdaValue,
            alpha,
            gamma,
            rho,
            zeta2,
            zeta12);

        if (price_cap == 1)
        {
            *ex_sprice = amortMidat_lgmcapval2F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                cpn_G2,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                zeta2,
                zeta12,
                ex_G,
                ex_G2);
        }
        else if (price_cap == 2)
        {
            *ex_sprice = amortMidat_lgmsumsopval2F_b(
                ncpn,
                nex,
                ex_cpn,
                cpn_df,
                cpn_G,
                cpn_G2,
                ex_zeta,
                zeta2,
                zeta12,
                ex_G,
                ex_G2,
                swaption_cpn);
        }
    }

    return NULL;
}

/*	Calibrate and price cap given lambda */
Err zcAmortMidat_lgmprcapgivenlambda(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flows */
    double zc_last_cpn[],   /*	Discounted Last Fix CF */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lprice[],     /*	Market prices for diagonal */
    double ex_sstrike[],    /*	Strikes for cap */
    double floatcoupon[],   /*	Float Coupon */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double lambda,          /*	Lambda */
    int    one2F,           /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last,  /*	If 1, the last option is disregarded
                                               and the forward volatility is flat from option
                                               n-1 */
    int     price_cap, /*	0: just calibrate */
    double* ex_sprice) /*	Cap price as output */
{
    static double cpn_G[MAX_CPN], cpn_G2[MAX_CPN], ex_G[MAX_CPN], ex_G2[MAX_CPN], zeta2[MAX_CPN],
        zeta12[MAX_CPN];
    Err err;

    /*	Setup G function */
    export_lgmsetupG(lambda, ncpn, cpn_time, cpn_G, nex, ex_time, ex_G);

    if (one2F == 2)
    {
        export_lgmsetupG2(lambda, gamma, ncpn, cpn_time, cpn_G2, nex, ex_time, ex_G2);
    }

    if (one2F == 1)
    {
        err = zcAmortMidat_lgmcalibzeta1F(
            ncpn,
            cpn,
            zc_last_cpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            cpn_G,
            nex,
            ex_bool,
            ex_time,
            ex_cpn,
            ex_G,
            ex_lprice,
            ex_zeta,
            floatnotional,
            lambda,
            skip_last);

        if (err)
        {
            return err;
        }

        if (price_cap)
        {
            *ex_sprice = amortMidat_lgmcapval1F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                ex_G);
        }
    }
    else
    {
        err = zcAmortMidat_lgmcalibzeta2F(
            ncpn,
            cpn,
            zc_last_cpn,
            cpn_time,
            cpn_df,
            cpn_cvg,
            cpn_G,
            cpn_G2,
            nex,
            ex_bool,
            ex_time,
            ex_cpn,
            ex_G,
            ex_G2,
            ex_lprice,
            ex_zeta,
            floatnotional,
            lambda,
            alpha,
            gamma,
            rho,
            skip_last);

        if (err)
        {
            return err;
        }

        export_lgmcalczeta2zeta12(nex, ex_time, ex_zeta, lambda, alpha, gamma, rho, zeta2, zeta12);

        if (price_cap)
        {
            *ex_sprice = amortMidat_lgmcapval2F_b(
                ncpn,
                cpn,
                cpn_df,
                cpn_cvg,
                cpn_G,
                cpn_G2,
                nex,
                ex_sstrike,
                floatcoupon,
                ex_cpn,
                ex_sncpn,
                ex_zeta,
                zeta2,
                zeta12,
                ex_G,
                ex_G2);
        }
    }

    return NULL;
}

static int NCPN_LM, NEX_LM, *EX_CPN_LM, *EX_LENDCPN_LM, *EX_SENDCPN_LM, ONE2F_LM, SKIP_LAST_LM,
    NB_MOM_LM, FREQ_SHORT_LM, SHIFT_FREQ_LM;

static double *CPN_TIME_LM, *CPN_DF_LM, *CPN_CVG_LM, *EX_TIME_LM, *EX_LSTRIKE_LM, *EX_LPRICE_LM,
    *EX_SSTRIKE_LM, *EX_SPRICE_LM, *EX_ZETA_LM, *LAM_TIME_LM, LAM_LM[MAX_CPN], ALPHA_LM, GAMMA_LM,
    RHO_LM, CPN_G_LM[MAX_CPN], CPN_G2_LM[MAX_CPN], EX_G_LM[MAX_CPN], EX_G2_LM[MAX_CPN],
    ZETA2_LM[MAX_CPN], ZETA12_LM[MAX_CPN], PRICE_LM[MAX_CPN], SENSI_LM[MAX_CPN][MAX_CPN];

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err amortMidat_lgmcalibzetalambda(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flows */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lprice[],     /*	Market prices for diagonal */
    double ex_fee[],        /*	Exercise Fee for diagonal */
    double ex_sstrike[],    /*	Strikes for cap */
    double floatcoupon[],   /*	Float Coupon */
    double ex_sprice,       /*	Market price for cap */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double swaption_cpn[],  /*	Standard Swaption Coupons	*/
    // double			standard_cvg[],						/*	Standard Coverages
    // */
    int fix_lambda,          /*	0: calib lambda to cap, 1: fix lambda calib
                                                     to diagonal */
    int     cap_or_swaption, /*	0: calib lambda to cap, 1: calib lambda to sum of diag swaptions */
    double* lambda,          /*	Lambda: may be changed in the process */
    int     one2F,           /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    int    it, compt;
    double lam1, lam2;
    double pr1, pr2;
    double fact, temp, temppr;
    int    ifact;
    Err    err;

    if (!fix_lambda && nex < 2)
    {
        return serror(
            "Cannot calibrate lambda and zeta with less than 2 exercises - choose fix lambda");
    }

    err = AmortMidat_lgmprcapgivenlambda(
        ncpn,
        cpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        nex,
        ex_bool,
        ex_time,
        ex_cpn,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        ex_sstrike,
        floatcoupon,
        swaption_cpn,
        ex_zeta,
        floatnotional,
        *lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last,
        (fix_lambda ? 0 : 1) + cap_or_swaption,
        &pr1);

    if (err || fix_lambda)
    {
        return err;
    }

    if (fabs(ex_sprice - pr1) < CALPRESAMORTMIDAT)
    {
        return NULL;
    }

    if (pr1 > ex_sprice)
    {
        ifact = -1;
    }
    else
    {
        ifact = 1;
    }

    if ((floatnotional[0] > floatnotional[ncpn - 2]) && (cap_or_swaption == 1))
    {
        ifact = -ifact;
    }
    //	else
    //	{
    //		ifact = - ifact;
    //	}

    lam1 = *lambda;
    fact = 0.25 * fabs(lam1);
    if (fact < 0.01)
    {
        fact = 0.01;
    }

    compt = 0;
    do
    {
        pr2 = pr1;

        fact *= 1.5;
        lam1 += ifact * fact;

        if (err = AmortMidat_lgmprcapgivenlambda(
                ncpn,
                cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_fee,
                ex_sstrike,
                floatcoupon,
                swaption_cpn,
                ex_zeta,
                floatnotional,
                lam1,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                1 + cap_or_swaption,
                &pr1))
        {
            return err;
        }

        if (fabs(ex_sprice - pr1) < CALPRESAMORTMIDAT)
        {
            *lambda = lam1;
            return NULL;
        }
        compt += 1;
    } while ((ifact * pr1 < ifact * ex_sprice) && (compt < 10));

    lam2 = lam1 - ifact * fact;

    if (lam1 > lam2)
    {
        temp   = lam1;
        temppr = pr1;
        lam1   = lam2;
        pr1    = pr2;
        lam2   = temp;
        pr2    = temppr;
    }

    for (it = 1; it < 20; it++)
    {
        temp   = lam2;
        temppr = pr2;
        lam2   = 0.5 * (lam1 + lam2);

        if (err = AmortMidat_lgmprcapgivenlambda(
                ncpn,
                cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_fee,
                ex_sstrike,
                floatcoupon,
                swaption_cpn,
                ex_zeta,
                floatnotional,
                lam2,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                1 + cap_or_swaption,
                &pr2))
        {
            return err;
        }

        if (fabs(ex_sprice - pr2) < CALPRESAMORTMIDAT)
        {
            *lambda = lam2;
            return NULL;
        }

        if (pr2 < ex_sprice)
        {
            pr1  = pr2;
            lam1 = lam2;
            lam2 = temp;
            pr2  = temppr;
        }

        if (fabs(lam2 - lam1) < CALPRESAMORTMIDAT || fabs(pr2 - pr1) < CALPRESAMORTMIDAT)
            break;

        while (pr1 > pr2 && it < 10)
        {
            it++;

            temp = lam1;
            lam1 = lam2;
            pr1  = pr2;
            lam2 += lam2 - temp;

            if (err = AmortMidat_lgmprcapgivenlambda(
                    ncpn,
                    cpn,
                    cpn_time,
                    cpn_df,
                    cpn_cvg,
                    nex,
                    ex_bool,
                    ex_time,
                    ex_cpn,
                    ex_sncpn,
                    ex_lprice,
                    ex_fee,
                    ex_sstrike,
                    floatcoupon,
                    swaption_cpn,
                    ex_zeta,
                    floatnotional,
                    lam2,
                    one2F,
                    alpha,
                    gamma,
                    rho,
                    skip_last,
                    1 + cap_or_swaption,
                    &pr2))
            {
                return err;
            }

            if (fabs(ex_sprice - pr2) < CALPRESAMORTMIDAT)
            {
                *lambda = lam2;
                return NULL;
            }
        }
    }

    smessage(
        "Calibration of lambda, error in bp: %.00f",
        10000 * min(fabs(ex_sprice - pr1), fabs(ex_sprice - pr2)));

    *lambda = lam2;

    return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err amortMidat_lgmcalibzetalambda_ts(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flows */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lprice[],     /*	Market prices for diagonal */
    double ex_fee[],        /*	Exercise Fee for diagonal */
    double ex_sstrike[],    /*	Strikes for cap */
    double floatcoupon[],   /*	Float Coupon */
    double ex_sprice,       /*	Market price for cap */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    double swaption_cpn[],  /*	Standard Swaption Coupons	*/
    // double			standard_cvg[],						/*	Standard Coverages
    // */
    int fix_lambda,      /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    int cap_or_swaption, /*	0: calib lambda to cap, 1: calib lambda to sum of diag swaptions */

    int     nNumLambda,
    double* pdLambdaValue,
    double* pdLambdaTime,

    int one2F, /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    int    it;
    double lam1, lam2;
    double pr1, pr2;
    double fact, temp, temppr;
    int    ifact;
    Err    err;

    //// only works for fixed lambda
    if (fix_lambda != 1)
    {
        return serror(
            "amortMidat_lgmcalibzetalambda_ts(...) :: this version only works for fixed lambda !");
    }

    if (!fix_lambda && nex < 2)
    {
        return serror(
            "Cannot calibrate lambda and zeta with less than 2 exercises - choose fix lambda");
    }

    err = AmortMidat_lgmprcapgivenlambda_ts(
        ncpn,
        cpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        nex,
        ex_bool,
        ex_time,
        ex_cpn,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        ex_sstrike,
        floatcoupon,
        swaption_cpn,
        ex_zeta,
        floatnotional,
        nNumLambda,
        pdLambdaValue,
        pdLambdaTime,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last,
        (fix_lambda ? 0 : 1),
        &pr1);

    if (err || fix_lambda)
    {
        return err;
    }

    if (fabs(ex_sprice - pr1) < CALPRESNOT)
    {
        return NULL;
    }

    if (pr1 > ex_sprice)
    {
        ifact = -1;
    }
    else
    {
        ifact = 1;
    }

    // lam1 = *lambda;
    lam1 = pdLambdaValue[0];
    fact = 0.25 * fabs(lam1);
    if (fact < 0.01)
    {
        fact = 0.01;
    }

    do
    {
        pr2 = pr1;

        fact *= 1.5;
        lam1 += ifact * fact;

        if (err = AmortMidat_lgmprcapgivenlambda_ts(
                ncpn,
                cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_fee,
                ex_sstrike,
                floatcoupon,
                swaption_cpn,
                ex_zeta,
                floatnotional,
                nNumLambda,
                pdLambdaValue,
                pdLambdaTime,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                (fix_lambda ? 0 : 1),
                &pr1))
        {
            return err;
        }

        if (fabs(ex_sprice - pr1) < CALPRESNOT)
        {
            //*lambda = lam1;
            return NULL;
        }
    } while (ifact * pr1 < ifact * ex_sprice);

    lam2 = lam1 - ifact * fact;

    if (lam1 > lam2)
    {
        temp   = lam1;
        temppr = pr1;
        lam1   = lam2;
        pr1    = pr2;
        lam2   = temp;
        pr2    = temppr;
    }

    for (it = 1; it < 20; it++)
    {
        temp   = lam2;
        temppr = pr2;
        lam2   = 0.5 * (lam1 + lam2);

        if (err = AmortMidat_lgmprcapgivenlambda_ts(
                ncpn,
                cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_fee,
                ex_sstrike,
                floatcoupon,
                swaption_cpn,
                ex_zeta,
                floatnotional,
                nNumLambda,
                pdLambdaValue,
                pdLambdaTime,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                (fix_lambda ? 0 : 1),
                &pr1))
        {
            return err;
        }

        if (fabs(ex_sprice - pr2) < CALPRESNOT)
        {
            //*lambda = lam2;
            return NULL;
        }

        if (pr2 < ex_sprice)
        {
            pr1  = pr2;
            lam1 = lam2;
            lam2 = temp;
            pr2  = temppr;
        }

        if (fabs(lam2 - lam1) < CALPRESNOT || fabs(pr2 - pr1) < CALPRESNOT)
            break;

        while (pr1 > pr2 && it < 10)
        {
            it++;

            temp = lam1;
            lam1 = lam2;
            pr1  = pr2;
            lam2 += lam2 - temp;

            if (err = AmortMidat_lgmprcapgivenlambda_ts(
                    ncpn,
                    cpn,
                    cpn_time,
                    cpn_df,
                    cpn_cvg,
                    nex,
                    ex_bool,
                    ex_time,
                    ex_cpn,
                    ex_sncpn,
                    ex_lprice,
                    ex_fee,
                    ex_sstrike,
                    floatcoupon,
                    swaption_cpn,
                    ex_zeta,
                    floatnotional,
                    nNumLambda,
                    pdLambdaValue,
                    pdLambdaTime,
                    one2F,
                    alpha,
                    gamma,
                    rho,
                    skip_last,
                    (fix_lambda ? 0 : 1),
                    &pr1))
            {
                return err;
            }

            if (fabs(ex_sprice - pr2) < CALPRESNOT)
            {
                //	*lambda = lam2;
                return NULL;
            }
        }
    }

    smessage(
        "Calibration of lambda, error in bp: %.00f",
        10000 * min(fabs(ex_sprice - pr1), fabs(ex_sprice - pr2)));

    //*lambda = lam2;
    return NULL;
}

/*	Calibrate zeta to diagonal and lambda to cap: both 1F and 2F */
Err zcamortMidat_lgmcalibzetalambda(
    int    ncpn,            /*	Total number of cash-flow dates */
    double cpn[],           /*	Discounted Cash-Flows */
    double zc_last_cpn[],   /*	Discounted Last Fix CF */
    double cpn_time[],      /*	Cash-Flow times */
    double cpn_df[],        /*	Df to cash-flow dates */
    double cpn_cvg[],       /*	cvg from i-1 to i */
    int    nex,             /*	Total number of exercise dates */
    int    ex_bool[],       /*	Exercise or not */
    double ex_time[],       /*	Exercise times */
    int    ex_cpn[],        /*	Index of the first cash-flow to be exercised */
    int    ex_sncpn[],      /*	Number of coupons in each caplet */
    double ex_lprice[],     /*	Market prices for diagonal */
    double ex_sstrike[],    /*	Strikes for cap */
    double floatcoupon[],   /*	Float Coupon */
    double ex_sprice,       /*	Market price for cap */
    double ex_zeta[],       /*	Output: zetas */
    double floatnotional[], /*	Float Notional	*/
    int    fix_lambda,      /*	0: calib lambda to cap, 1: fix lambda calib
                                                    to diagonal */
    double* lambda,         /*	Lambda: may be changed in the process */
    int     one2F,          /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,
    int    skip_last) /*	If 1, the last option is disregarded
                                           and the forward volatility is flat from option
                                           n-1 */
{
    int    it;
    double lam1, lam2;
    double pr1, pr2;
    double fact, temp, temppr;
    int    ifact;
    Err    err;

    if (!fix_lambda && nex < 2)
    {
        return serror(
            "Cannot calibrate lambda and zeta with less than 2 exercises - choose fix lambda");
    }

    err = zcAmortMidat_lgmprcapgivenlambda(
        ncpn,
        cpn,
        zc_last_cpn,
        cpn_time,
        cpn_df,
        cpn_cvg,
        nex,
        ex_bool,
        ex_time,
        ex_cpn,
        ex_sncpn,
        ex_lprice,
        ex_sstrike,
        floatcoupon,
        ex_zeta,
        floatnotional,
        *lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last,
        (fix_lambda ? 0 : 1),
        &pr1);

    if (err || fix_lambda)
    {
        return err;
    }

    if (fabs(ex_sprice - pr1) < CALPRESAMORTMIDAT)
    {
        return NULL;
    }

    if (pr1 > ex_sprice)
    {
        ifact = -1;
    }
    else
    {
        ifact = 1;
    }

    lam1 = *lambda;
    fact = 0.25 * fabs(lam1);
    if (fact < 0.01)
    {
        fact = 0.01;
    }

    do
    {
        pr2 = pr1;

        fact *= 1.5;
        lam1 += ifact * fact;

        if (err = zcAmortMidat_lgmprcapgivenlambda(
                ncpn,
                cpn,
                zc_last_cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_sstrike,
                floatcoupon,
                ex_zeta,
                floatnotional,
                lam1,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                1,
                &pr1))
        {
            return err;
        }

        if (fabs(ex_sprice - pr1) < CALPRESAMORTMIDAT)
        {
            *lambda = lam1;
            return NULL;
        }
    } while (ifact * pr1 < ifact * ex_sprice);

    lam2 = lam1 - ifact * fact;

    if (lam1 > lam2)
    {
        temp   = lam1;
        temppr = pr1;
        lam1   = lam2;
        pr1    = pr2;
        lam2   = temp;
        pr2    = temppr;
    }

    for (it = 1; it < 20; it++)
    {
        temp   = lam2;
        temppr = pr2;
        lam2   = 0.5 * (lam1 + lam2);

        if (err = zcAmortMidat_lgmprcapgivenlambda(
                ncpn,
                cpn,
                zc_last_cpn,
                cpn_time,
                cpn_df,
                cpn_cvg,
                nex,
                ex_bool,
                ex_time,
                ex_cpn,
                ex_sncpn,
                ex_lprice,
                ex_sstrike,
                floatcoupon,
                ex_zeta,
                floatnotional,
                lam2,
                one2F,
                alpha,
                gamma,
                rho,
                skip_last,
                1,
                &pr2))
        {
            return err;
        }

        if (fabs(ex_sprice - pr2) < CALPRESAMORTMIDAT)
        {
            *lambda = lam2;
            return NULL;
        }

        if (pr2 < ex_sprice)
        {
            pr1  = pr2;
            lam1 = lam2;
            lam2 = temp;
            pr2  = temppr;
        }

        if (fabs(lam2 - lam1) < CALPRESAMORTMIDAT || fabs(pr2 - pr1) < CALPRESAMORTMIDAT)
            break;

        while (pr1 > pr2 && it < 10)
        {
            it++;

            temp = lam1;
            lam1 = lam2;
            pr1  = pr2;
            lam2 += lam2 - temp;

            if (err = zcAmortMidat_lgmprcapgivenlambda(
                    ncpn,
                    cpn,
                    zc_last_cpn,
                    cpn_time,
                    cpn_df,
                    cpn_cvg,
                    nex,
                    ex_bool,
                    ex_time,
                    ex_cpn,
                    ex_sncpn,
                    ex_lprice,
                    ex_sstrike,
                    floatcoupon,
                    ex_zeta,
                    floatnotional,
                    lam2,
                    one2F,
                    alpha,
                    gamma,
                    rho,
                    skip_last,
                    1,
                    &pr2))
            {
                return err;
            }

            if (fabs(ex_sprice - pr2) < CALPRESAMORTMIDAT)
            {
                *lambda = lam2;
                return NULL;
            }
        }
    }

    smessage(
        "Calibration of lambda, error in bp: %.00f",
        10000 * min(fabs(ex_sprice - pr1), fabs(ex_sprice - pr2)));

    *lambda = lam2;

    return NULL;
}

Err amortMidat_modify_dates(
    int   Nfix,
    long* fix_start_datesIn, /*	Start date of the amortised swap */
    long* fix_end_datesIn,   /*	End date for the amortised swap */

    int   Nfloat,
    long* float_start_datesIn, /*	Start date of the amortised swap */
    long* float_end_datesIn,

    long* fix_start_datesOut, /*	Start date of the amortised swap */
    long* fix_end_datesOut,   /*	End date for the amortised swap */

    long* float_start_datesOut, /*	Start date of the amortised swap */
    long* float_end_datesOut)   /*	End date for the amortised swap */
{
    Err    err = NULL;
    int    i;
    int    floatindex;
    double lag1, lag2;

    if (Nfix <= Nfloat)
    {
        floatindex = 0;
        for (i = 0; i < Nfix; ++i)
        {
            while ((float_end_datesIn[floatindex] < fix_end_datesIn[i]) && floatindex < Nfloat)
            {
                float_start_datesOut[floatindex] = float_start_datesIn[floatindex];
                float_end_datesOut[floatindex]   = float_end_datesIn[floatindex];
                floatindex                       = floatindex + 1;
            }

            if ((floatindex > 0) && (floatindex < Nfloat))
            {
                lag1 = fabs(float_end_datesIn[floatindex] - fix_end_datesIn[i]);
                lag2 = fabs(float_end_datesIn[floatindex - 1] - fix_end_datesIn[i]);

                if (lag1 < lag2)
                {
                    if (floatindex < Nfloat - 1)
                    {
                        float_start_datesOut[floatindex + 1] = fix_end_datesIn[i];
                    }
                    float_end_datesOut[floatindex] = fix_end_datesIn[i];
                }
                else
                {
                    float_start_datesOut[floatindex]   = fix_end_datesIn[i];
                    float_end_datesOut[floatindex - 1] = fix_end_datesIn[i];
                }
            }
            else if (floatindex == 0)
            {
                if (Nfix == Nfloat)
                {
                    float_start_datesOut[floatindex] = fix_start_datesIn[i];
                    float_end_datesOut[floatindex]   = fix_end_datesIn[i];
                }
                else
                {
                    float_start_datesOut[floatindex] = float_start_datesIn[floatindex];
                    float_end_datesOut[floatindex]   = float_end_datesIn[floatindex];
                }
                floatindex = floatindex + 1;
            }
            else if (floatindex == Nfloat)
            {
                float_end_datesOut[floatindex - 1] = fix_end_datesIn[i];
            }

            fix_start_datesOut[i] = fix_start_datesIn[i];
            fix_end_datesOut[i]   = fix_end_datesIn[i];
        }
    }
    else
    {
        err = "Fix freq should be larger than Float Freq";
    }

    return err;
}

Err amortMidat_modify_dates2(
    int   Nfix,
    long* fix_pay_datesIn, /*	Pay date for the amortised swap */

    int   Nfloat,
    long* float_pay_datesIn,

    long* fix_pay_datesOut, /*	Pay date for the amortised swap */

    long* float_pay_datesOut) /*	Pay date of the amortised swap */
{
    Err    err = NULL;
    int    i;
    int    floatindex;
    double lag1, lag2;

    if (Nfix <= Nfloat)
    {
        floatindex = 0;
        for (i = 0; i < Nfix; ++i)
        {
            while ((float_pay_datesIn[floatindex] < fix_pay_datesIn[i]) && floatindex < Nfloat)
            {
                float_pay_datesOut[floatindex] = float_pay_datesIn[floatindex];
                floatindex                     = floatindex + 1;
            }

            if ((floatindex > 0) && (floatindex < Nfloat))
            {
                lag1 = fabs(float_pay_datesIn[floatindex] - fix_pay_datesIn[i]);
                lag2 = fabs(float_pay_datesIn[floatindex - 1] - fix_pay_datesIn[i]);

                if (lag1 < lag2)
                {
                    //					if(floatindex<Nfloat-1)
                    //					{
                    //						float_pay_datesOut[floatindex+1] =
                    //fix_pay_datesIn[i];
                    //					}
                    float_pay_datesOut[floatindex] = fix_pay_datesIn[i];
                }
                else
                {
                    float_pay_datesOut[floatindex - 1] = fix_pay_datesIn[i];
                }
            }
            else if (floatindex == 0)
            {
                if (Nfix == Nfloat)
                {
                    float_pay_datesOut[floatindex] = fix_pay_datesIn[i];
                }
                else
                {
                    float_pay_datesOut[floatindex] = float_pay_datesIn[floatindex];
                }
                floatindex = floatindex + 1;
            }
            else if (floatindex == Nfloat)
            {
                float_pay_datesOut[floatindex - 1] = fix_pay_datesIn[i];
            }

            fix_pay_datesOut[i] = fix_pay_datesIn[i];
            fix_pay_datesOut[i] = fix_pay_datesIn[i];
        }
    }
    else
    {
        err = "Fix freq should be larger than Float Freq";
    }

    return err;
}

//	Calibrate lgm: main function
Err amortMidat_cpd_calib_diagonal_new(
    int   notperiod,
    char* yc_name,
    char* vol_curve_name,

    char* default_ref_rate_name,
    char* swaption_basis,
    char* swaption_freq,
    Err (*get_cash_vol)(
        char*   vol_curve_name,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),
    double vol_shift,
    int    shift_type,

    int*    ex_bool,
    double* diag_prices,
    double* ex_fee,

    char* fix_freq,
    char* fix_basis,
    int   Nfix,
    long* fix_start_dates,
    long* fix_end_dates,
    long* fixpaydates,

    double* fix_rates,
    double* fix_notionals,

    char* float_freq,
    char* float_basis,
    int   Nfloat,
    long* float_start_dates,
    long* float_end_dates,
    long* floatpaydates,

    double* float_margins,
    double* float_spreads,
    double* float_notionals,

    // Calibration Parameters
    double* short_strike,
    int     strike_type,
    double  max_std_short,

    int fix_lambda,
    int one_f_equi,
    int skip_last,

    int    use_jump,
    double max_var_jump,

    double* lambda,
    int     one2F,
    //	Alpha, Gamma, Rho (2F only)
    double alpha,
    double gamma,
    double rho,

    int*     num_sig,
    double** sig_time,
    double** sig)
{  /// start of amortMidat_cpd_calib_diagonal_new
    int            i, j, k, n, float_index, fix_index;
    int            nbcoupon;
    SrtCompounding fixfreq, floatfreq;
    SrtBasisCode   fixbasis, floatbasis, swapbasis;
    double         ex_lstrike[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_sncpn[MAX_CPN];
    double      cpn_G[MAX_CPN];
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;
    double      coupon[MAX_CPN];
    double      swaption_cpn[MAX_CPN];
    double      floatcoupon[MAX_CPN];
    double      shortcoupon[MAX_CPN];
    double      coupon_time[MAX_CPN];
    double      coupon_df[MAX_CPN];
    double      coupon_cvg[MAX_CPN];
    double      coupon_cvg_standard[MAX_CPN];
    int         ex_coupon[MAX_CPN];

    double temp;

    int fix_float_mult;

    int last_i;

    double* dIRRs       = NULL;
    double* dIRRstrikes = NULL;

    double *float_cvgs = NULL, *float_pay_times = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;
    long *  float_pay_dates = NULL, *fix_pay_dates = NULL;

    double* coupon_times = NULL;

    double maxNotional;

    int    nbexercise;
    double exer_time[MAX_CPN];
    long   exer_date[MAX_CPN];
    int    exer_bool[MAX_CPN];

    int firstDoCalib, firstDoCalibFloat;

    int cap_or_swaption;

    int compt;

    double calibstrike;

    int UseVol;

    //-----If no lambda calibration : change strike_type---------------
    if (fix_lambda)
    {
        strike_type = 0;
    }

    //-----For IRR equivalent strike : Use vol or not ?---------
    UseVol = 0;
    if ((strike_type == 5) || (strike_type == 7))
    {
        UseVol = 1;
    }

    //-----Secondary instrument : Short or Diag Swaption ?---------
    cap_or_swaption = 0;
    if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
    {
        cap_or_swaption = 1;
    }

    //-----Parameter controlling vol jump size---------
    MAX_FACT  = max_var_jump;
    USE_JUMPS = use_jump;

    //-----Initialize outpu pointers-----------------
    *sig_time = NULL;
    *sig      = NULL;

    //-----Retreive yield curve from its name---------
    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    //-----If LGM2F compute Hermite quadrature weights and ponits---------
    if (one2F == 2)
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    //-----convert char in Sort structure---------
    err = interp_basis(fix_basis, &fixbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(fix_freq, &fixfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(float_basis, &floatbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(float_freq, &floatfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(swaption_basis, &swapbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----get freq and basis corresponding to ref_rate---------
    //	err = swp_f_get_ref_rate_details(ref_rate, &floatbasis, &floatfreq);

    //-------------Compute Coverages--------------
    fix_cvgs = dvector(0, Nfix - 1);
    for (i = 0; i < Nfix; ++i)
    {
        fix_cvgs[i] = coverage(fix_start_dates[i], fix_end_dates[i], fixbasis);
    }
    float_cvgs = dvector(0, Nfloat - 1);
    for (i = 0; i < Nfloat; ++i)
    {
        float_cvgs[i] = coverage(float_start_dates[i], float_end_dates[i], floatbasis);
    }
    fix_float_mult = floatfreq / fixfreq;

    //-------------Find First Exercise Date----------------
    if (add_unit(fix_start_dates[0], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) <= today)
    {
        ex_bool[0] = 0;
    }

    i = 0;
    while ((ex_bool[i] == 0) && (i < Nfix))
    {
        ++i;
    }

    if (i == Nfix)
    {
        err = "At least one maturity must be calibrated : Check DoCalib Value";
        goto FREE_RETURN;
    }

    firstDoCalib = i;

    //--------------Shift all data---------------------------
    for (i = firstDoCalib; i < Nfix; ++i)
    {
        ex_bool[i - firstDoCalib]     = ex_bool[i];
        exer_bool[i - firstDoCalib]   = ex_bool[i];
        diag_prices[i - firstDoCalib] = diag_prices[i];
        ex_fee[i - firstDoCalib]      = ex_fee[i];

        fix_start_dates[i - firstDoCalib] = fix_start_dates[i];
        fix_end_dates[i - firstDoCalib]   = fix_end_dates[i];
        fixpaydates[i - firstDoCalib]     = fixpaydates[i];

        fix_rates[i - firstDoCalib]     = fix_rates[i];
        fix_notionals[i - firstDoCalib] = fix_notionals[i];
        fix_cvgs[i - firstDoCalib]      = fix_cvgs[i];

        if (short_strike)
        {
            short_strike[i - firstDoCalib] = short_strike[i];
        }

        coupon_cvg_standard[i - firstDoCalib] =
            coverage(fix_start_dates[i - firstDoCalib], fix_end_dates[i - firstDoCalib], swapbasis);
    }

    i = 0;
    while (float_start_dates[i] < fix_start_dates[0] - 10)
    {
        ++i;
    }

    firstDoCalibFloat = i;

    for (i = firstDoCalibFloat; i < Nfloat; ++i)
    {
        float_start_dates[i - firstDoCalibFloat] = float_start_dates[i];
        float_end_dates[i - firstDoCalibFloat]   = float_end_dates[i];
        floatpaydates[i - firstDoCalibFloat]     = floatpaydates[i];

        float_margins[i - firstDoCalibFloat]   = float_margins[i];
        float_spreads[i - firstDoCalibFloat]   = float_spreads[i];
        float_notionals[i - firstDoCalibFloat] = float_notionals[i];
        float_cvgs[i - firstDoCalibFloat]      = float_cvgs[i];
    }

    float_pay_times    = dvector(0, Nfloat);
    float_pay_dates    = lngvector(0, Nfloat);
    float_pay_times[0] = (float_start_dates[0] - today) * YEARS_IN_DAY;
    float_pay_dates[0] = float_start_dates[0];
    for (n = 1; n <= Nfloat; ++n)
    {
        float_pay_times[n] = (floatpaydates[n - 1] - today) * YEARS_IN_DAY;
        float_pay_dates[n] = floatpaydates[n - 1];
    }

    fix_pay_times    = dvector(0, Nfix);
    fix_pay_dates    = lngvector(0, Nfix);
    fix_pay_times[0] = (fix_start_dates[0] - today) * YEARS_IN_DAY;
    fix_pay_dates[0] = fix_start_dates[0];
    for (n = 1; n <= Nfix; ++n)
    {
        fix_pay_times[n] = (fixpaydates[n - 1] - today) * YEARS_IN_DAY;
        fix_pay_dates[n] = fixpaydates[n - 1];
    }

    //----------Rescale Notionals---------------------
    maxNotional = fabs(fix_notionals[0]);
    for (i = 1; i < Nfix - firstDoCalib; ++i)
    {
        if (maxNotional < fabs(fix_notionals[i]))
        {
            maxNotional = fabs(fix_notionals[i]);
        }
    }
    for (i = 0; i < Nfloat - firstDoCalibFloat; ++i)
    {
        if (maxNotional < fabs(float_notionals[i]))
        {
            maxNotional = fabs(float_notionals[i]);
        }
    }

    CALPRESNOT = maxNotional * CALPRESAMORTMIDAT;

    //-----------------------------------------------------
    //-------------Build Amort Swap Coupons----------------
    //-----------------------------------------------------

    nbcoupon = Nfloat + 1 - firstDoCalibFloat;

    coupon_times = (double*)calloc(nbcoupon, sizeof(double));
    memcpy(coupon_times, float_pay_times, nbcoupon * sizeof(double));

    for (k = 0; k < nbcoupon; ++k)
    {
        coupon[k]       = 0;
        floatcoupon[k]  = 0;
        swaption_cpn[k] = 0;
    }

    nbexercise  = Nfix - firstDoCalib;
    float_index = 1;
    fix_index   = 1;
    coupon[0]   = -float_notionals[0] * swp_f_df(today, float_start_dates[0], yc_name);
    //	swaption_cpn[0] = - swp_f_df (today, float_start_dates[0], yc_name);
    floatcoupon[0] = -swp_f_df(today, float_start_dates[0], yc_name);
    coupon_time[0] = coupon_times[0];
    coupon_df[0]   = swp_f_df(today, float_start_dates[0], yc_name);
    coupon_cvg[0]  = 0;
    for (k = 1; k < nbcoupon; ++k)
    {
        coupon_time[k] = coupon_times[k];
        coupon_df[k]   = swp_f_df(today, floatpaydates[k - 1], yc_name);
        coupon_cvg[k]  = float_cvgs[k - 1];
        if ((fabs(coupon_times[k] - fix_pay_times[fix_index]) < 10.0 / 365.0) &&
            (fabs(coupon_times[k] - float_pay_times[float_index]) > 10.0 / 365.0))
        {
            coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                         fix_cvgs[fix_index - 1] *
                         swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
            ex_coupon[fix_index - 1] = k - fix_float_mult;
            exer_date[fix_index - 1] =
                add_unit(fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
            exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
            fix_index                = fix_index + 1;
        }
        else if (
            (fabs(coupon_times[k] - fix_pay_times[fix_index]) > 10.0 / 365.0) &&
            (fabs(coupon_times[k] - float_pay_times[float_index]) < 10.0 / 365.0))
        {
            temp = swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            coupon[k] += -(float_notionals[float_index] - float_notionals[float_index - 1] +
                           float_notionals[float_index - 1] *
                               (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                         swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            floatcoupon[k] += -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                              float_cvgs[float_index - 1] *
                              swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            float_index = float_index + 1;
        }
        else
        {
            if (float_index < nbcoupon - 1)
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(float_notionals[float_index] - float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                    float_cvgs[float_index - 1] *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(-float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(-1.0 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
        }
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Amort Swap Coupon-------------------
    //-----------------------------------------------------------------

    //	Underlyings

    //	Long

    dff = swp_f_df(today, floatpaydates[Nfloat - 1 - firstDoCalibFloat], yc_name);
    for (i = 0; i < nbexercise; i++)
    {
        ex_lstrike[i] = fix_rates[i];
        ex_lprice[i]  = diag_prices[i];
    }

    dIRRs = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRs)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //-----------------------------------------------------------------
    //----------------Compute Diagonal IRRs of the amort Swap----------
    //-----------------------------------------------------------------
    if ((strike_type == 4) || (strike_type == 5) || (strike_type == 6) || (strike_type == 7))
    {
        err = ComputeAmortSwapDiagonalIRRsNew2(
            yc_name,
            vol_curve_name,
            default_ref_rate_name,
            fixfreq,
            fixbasis,

            Nfix - firstDoCalib,
            fix_start_dates,
            fix_end_dates,
            fix_cvgs,
            fix_rates,
            fix_notionals,
            ex_fee,

            Nfloat - firstDoCalibFloat,
            float_start_dates,
            float_end_dates,
            float_cvgs,
            float_margins,
            float_spreads,
            float_notionals,

            dIRRs,
            UseVol);
        if (err)
        {
            goto FREE_RETURN;
        }
    }

    dIRRstrikes = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRstrikes)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //--------------------------------------------------------------
    //----------Secondary instruments for Lambda calibration --------
    //--------------------------------------------------------------
    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute equivalent swaption strike using diagonal irrs-----
            else if ((strike_type == 6) || (strike_type == 7))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    fix_end_dates[Nfix - firstDoCalib - 1],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute diagonal swap, level and vol---------------------
            if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
            {
                lvl = 0.0;
                for (k = i; k < Nfix - firstDoCalib; k++)
                {
                    lvl += coupon_cvg_standard[k] * swp_f_df(today, fix_end_dates[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_end_dates[k - 1], yc_name);

                ex_slvl[i] = lvl;

                ex_sfwd[i] = (dfi - dff) / lvl;

                //	ATM std
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_end_dates[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute short swap, level and vol---------------------
            else
            {
                lvl = 0.0;
                for (k = i; k < i + 1; k++)
                {
                    lvl += fix_cvgs[k] * swp_f_df(today, fix_end_dates[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_end_dates[k - 1], yc_name);

                ex_slvl[i] = lvl;
                ex_sfwd[i] = (dfi - dff) / lvl;

                //	ATM std
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_end_dates[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            //---------------------------------------------------------------
            //---------Compute effective strike of secondary instrument -----
            //---------------------------------------------------------------

            //-----------------------Compute std-----------------------------
            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(exer_time[i]);

            //-----------ATM-----------------------------
            if (strike_type == 0)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            //-----------Custom Cash Strike--------------
            else if ((strike_type == 1) || (strike_type == 8))
            {
                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i];
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom Swap Strike -------------
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        fix_start_dates[i],
                        fix_end_dates[k - 1],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        default_ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i] - spr;
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom STD Strike --------------
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }
            //-----------IRR Equivalent Strike ----------
            else if (
                (strike_type == 4) || (strike_type == 5) || (strike_type == 6) ||
                (strike_type == 7))
            {
                ex_sstrike[i] = dIRRstrikes[i];
            }

            //----------Apply max std--------------------
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            //-----	Make sure strikes are positive (actually more than 1bp)--
            //-----	otherwise use ATM	-------------------------------------
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            calibstrike = ex_sstrike[i];

            //--------Price secondary instrument---------------
            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_end_dates[k - 1],
                calibstrike,
                0,
                default_ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            cap_price += ex_sprice[i];
            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }
    }

    //------------------------------------------------------------
    //--------- Fill Standard Swaption Coupons--------------------
    //------------------------------------------------------------
    if ((strike_type >= 6) && (!fix_lambda))
    {
        swaption_cpn[0] = -swp_f_df(today, float_start_dates[0], yc_name);
        float_index     = 1;
        fix_index       = 1;
        for (k = 1; k < nbcoupon; ++k)
        {
            if ((coupon_times[k] == fix_pay_times[fix_index]) &&
                (coupon_times[k] != float_pay_times[float_index]))
            {
                swaption_cpn[k] += ex_sstrike[fix_index - 1] * coupon_cvg_standard[fix_index - 1] *
                                   swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_times[k] != fix_pay_times[fix_index]) &&
                (coupon_times[k] == float_pay_times[float_index]))
            {
                swaption_cpn[k] +=
                    -((float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                      float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < nbcoupon - 1)
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] *
                                       coupon_cvg_standard[fix_index - 1] *
                                       swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                    fix_index   = fix_index + 1;
                    float_index = float_index + 1;
                }
                else
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                                       swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                    fix_index = fix_index + 1;
                    swaption_cpn[k] +=
                        -(-1 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                        swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                    float_index = float_index + 1;
                }
            }
        }
    }
    //------------------------------------------------------------
    //--------- End Of Fill Standard Swaption Coupons-------------
    //------------------------------------------------------------

    //	The 1F equivalent case
    if (one2F == 2 && fix_lambda && one_f_equi && (strike_type < 6))
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute short swap rate, level and vol----------
            lvl = 0.0;
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_end_dates[k + 1], yc_name);
            }
            dfi = swp_f_df(today, fix_end_dates[i], yc_name);
            dff = swp_f_df(today, fix_end_dates[k], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            //--------Compute effective short strike--------------

            //--------Equivalent IRR strike-----------------------
            if ((strike_type == 4) || (strike_type == 5))
            {
                ex_sstrike[i] = dIRRstrikes[i];
                if (dIRRstrikes[i] < 0.0001)
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //--------ATM-----------------------------------------
            else
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }

        //--------Calibrate LGM1F vol TS with given lambda ---------
        err = amortMidat_lgmcalibzetalambda(
            nbcoupon,
            coupon,
            coupon_time,
            coupon_df,
            coupon_cvg,
            nbexercise,
            exer_bool,
            exer_time,
            ex_coupon,
            ex_sncpn,
            ex_lprice,
            ex_fee,
            shortcoupon,
            floatcoupon,
            0.0,
            ex_zeta,
            float_notionals,
            swaption_cpn,
            1,
            0,
            lambda,
            1,
            0.0,
            0.0,
            0.0,
            skip_last);

        if (err)
        {
            goto FREE_RETURN;
        }

        //--------Price the cap within the calibrated LGM1F model----------------
        export_lgmsetupG(*lambda, nbcoupon, coupon_time, cpn_G, nbexercise, exer_time, ex_G);

        cap_price = amortMidat_lgmcapval1F_b(
            nbcoupon,
            coupon,
            coupon_df,
            coupon_cvg,
            cpn_G,
            nbexercise,
            shortcoupon,
            floatcoupon,
            ex_coupon,
            ex_sncpn,
            ex_zeta,
            ex_G);

        fix_lambda = 0;
    }

    //--------Calibrate lambda and zeta----------------
    err = amortMidat_lgmcalibzetalambda(
        nbcoupon,
        coupon,
        coupon_time,
        coupon_df,
        coupon_cvg,
        nbexercise,
        exer_bool,
        exer_time,
        ex_coupon,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        shortcoupon,
        floatcoupon,
        cap_price,
        ex_zeta,
        float_notionals,
        swaption_cpn,
        fix_lambda,
        cap_or_swaption,
        lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    //	3.)	Transform into sigma

    compt = 0;
    for (i = 0; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            compt += 1;
        }
    }

    *num_sig  = compt;
    *sig_time = (double*)calloc(*num_sig, sizeof(double));
    *sig      = (double*)calloc(*num_sig, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    (*sig_time)[0] = exer_time[0];
    (*sig)[0]      = sqrt(ex_zeta[0] * 2 * (*lambda) / (exp(2 * (*lambda) * exer_time[0]) - 1.0));

    last_i = 0;
    k      = 1;
    for (i = 1; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            (*sig_time)[k] = exer_time[i];
            if (ex_zeta[i] > ex_zeta[last_i])
            {
                (*sig)[k] = sqrt(
                    (ex_zeta[i] - ex_zeta[last_i]) * 2 * (*lambda) /
                    (exp(2 * (*lambda) * exer_time[i]) - exp(2 * (*lambda) * exer_time[last_i])));
            }
            else
            {
                smessage(
                    "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                    exer_time[i]);
                for (j = k; j < *num_sig; j++)
                {
                    (*sig_time)[j] = exer_time[k];
                    (*sig)[j]      = (*sig)[k - 1];
                }
                i = nbexercise;
            }

            last_i = i;
            k      = k + 1;
        }
    }

FREE_RETURN:

    if (fix_cvgs)
        free_dvector(fix_cvgs, 0, Nfix - 1);
    if (float_cvgs)
        free_dvector(float_cvgs, 0, Nfloat - 1);

    if (float_pay_dates)
        free_lngvector(float_pay_dates, 0, Nfloat);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, Nfloat);

    if (fix_pay_dates)
        free_lngvector(fix_pay_dates, 0, Nfix);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, Nfix);

    if (coupon_times)
        free(coupon_times);

    if (dIRRs)
        free(dIRRs);
    if (dIRRstrikes)
        free(dIRRstrikes);

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;
    }

    return err;
}  //// end of amortMidat_cpd_calib_diagonal_new

//	Calibrate lgm: main function with lambda ts enabled
Err amortMidat_cpd_calib_diagonal_new_ts(
    int   notperiod,
    char* yc_name,
    char* vol_curve_name,

    char* default_ref_rate_name,
    char* swaption_basis,
    char* swaption_freq,
    Err (*get_cash_vol)(
        char*   vol_curve_name,
        double  start_date,
        double  end_date,
        double  cash_strike,
        int     zero,
        char*   ref_rate_name,
        double* vol,
        double* power),
    double vol_shift,
    int    shift_type,

    int*    ex_bool,
    double* diag_prices,
    double* ex_fee,

    char* fix_freq,
    char* fix_basis,
    int   Nfix,
    long* fix_start_dates,
    long* fix_end_dates,
    long* fixpaydates,

    double* fix_rates,
    double* fix_notionals,

    char* float_freq,
    char* float_basis,
    int   Nfloat,
    long* float_start_dates,
    long* float_end_dates,
    long* floatpaydates,

    double* float_margins,
    double* float_spreads,
    double* float_notionals,

    // Calibration Parameters
    double* short_strike,
    int     strike_type,
    double  max_std_short,

    int fix_lambda,
    int one_f_equi,
    int skip_last,

    int    use_jump,
    double max_var_jump,

    int     nNum_lambda,
    double* pdLambdaValue,
    double* pdLambdaDate,

    int one2F,
    //	Alpha, Gamma, Rho (2F only)
    double alpha,
    double gamma,
    double rho,

    int*     num_sig,
    double** sig_time,
    double** sig)
{  /// start of amortMidat_cpd_calib_diagonal_new_ts
    int            i, j, k, n, float_index, fix_index;
    int            nbcoupon;
    SrtCompounding fixfreq, floatfreq;
    SrtBasisCode   fixbasis, floatbasis, swapbasis;
    double         ex_lstrike[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_sncpn[MAX_CPN];
    double      cpn_G[MAX_CPN];
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;
    double      coupon[MAX_CPN];
    double      swaption_cpn[MAX_CPN];
    double      floatcoupon[MAX_CPN];
    double      shortcoupon[MAX_CPN];
    double      coupon_time[MAX_CPN];
    double      coupon_df[MAX_CPN];
    double      coupon_cvg[MAX_CPN];
    double      coupon_cvg_standard[MAX_CPN];
    int         ex_coupon[MAX_CPN];

    double temp;

    int fix_float_mult;

    int last_i;

    double* dIRRs       = NULL;
    double* dIRRstrikes = NULL;

    double *float_cvgs = NULL, *float_pay_times = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;
    long *  float_pay_dates = NULL, *fix_pay_dates = NULL;

    double* coupon_times = NULL;

    double maxNotional;

    int    nbexercise;
    double exer_time[MAX_CPN];
    long   exer_date[MAX_CPN];
    int    exer_bool[MAX_CPN];

    int firstDoCalib, firstDoCalibFloat;

    int cap_or_swaption;

    int compt;

    double calibstrike;

    int UseVol;

    double exp_fact;

    //-----If no lambda calibration : change strike_type---------------
    if (fix_lambda)
    {
        strike_type = 0;
    }

    //-----For IRR equivalent strike : Use vol or not ?---------
    UseVol = 0;
    if ((strike_type == 5) || (strike_type == 7))
    {
        UseVol = 1;
    }

    //-----Secondary instrument : Short or Diag Swaption ?---------
    cap_or_swaption = 0;
    if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
    {
        cap_or_swaption = 1;
    }

    //-----Parameter controlling vol jump size---------
    MAX_FACT  = max_var_jump;
    USE_JUMPS = use_jump;

    //-----Initialize outpu pointers-----------------
    *sig_time = NULL;
    *sig      = NULL;

    //-----Retreive yield curve from its name---------
    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    //-----If LGM2F compute Hermite quadrature weights and ponits---------
    if (one2F == 2)
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    //-----convert char in Sort structure---------
    err = interp_basis(fix_basis, &fixbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(fix_freq, &fixfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(float_basis, &floatbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(float_freq, &floatfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(swaption_basis, &swapbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----get freq and basis corresponding to ref_rate---------
    //	err = swp_f_get_ref_rate_details(ref_rate, &floatbasis, &floatfreq);

    //-------------Compute Coverages--------------
    fix_cvgs = dvector(0, Nfix - 1);
    for (i = 0; i < Nfix; ++i)
    {
        fix_cvgs[i] = coverage(fix_start_dates[i], fix_end_dates[i], fixbasis);
    }
    float_cvgs = dvector(0, Nfloat - 1);
    for (i = 0; i < Nfloat; ++i)
    {
        float_cvgs[i] = coverage(float_start_dates[i], float_end_dates[i], floatbasis);
    }
    fix_float_mult = floatfreq / fixfreq;

    //-------------Find First Exercise Date----------------
    if (add_unit(fix_start_dates[0], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) <= today)
    {
        ex_bool[0] = 0;
    }

    i = 0;
    while ((ex_bool[i] == 0) && (i < Nfix))
    {
        ++i;
    }

    if (i == Nfix)
    {
        err = "At least one maturity must be calibrated : Check DoCalib Value";
        goto FREE_RETURN;
    }

    firstDoCalib = i;

    //--------------Shift all data---------------------------
    for (i = firstDoCalib; i < Nfix; ++i)
    {
        ex_bool[i - firstDoCalib]     = ex_bool[i];
        exer_bool[i - firstDoCalib]   = ex_bool[i];
        diag_prices[i - firstDoCalib] = diag_prices[i];
        ex_fee[i - firstDoCalib]      = ex_fee[i];

        fix_start_dates[i - firstDoCalib] = fix_start_dates[i];
        fix_end_dates[i - firstDoCalib]   = fix_end_dates[i];
        fixpaydates[i - firstDoCalib]     = fixpaydates[i];

        fix_rates[i - firstDoCalib]     = fix_rates[i];
        fix_notionals[i - firstDoCalib] = fix_notionals[i];
        fix_cvgs[i - firstDoCalib]      = fix_cvgs[i];

        if (short_strike)
        {
            short_strike[i - firstDoCalib] = short_strike[i];
        }

        coupon_cvg_standard[i - firstDoCalib] =
            coverage(fix_start_dates[i - firstDoCalib], fix_end_dates[i - firstDoCalib], swapbasis);
    }

    i = 0;
    while (float_start_dates[i] < fix_start_dates[0] - 10)
    {
        ++i;
    }

    firstDoCalibFloat = i;

    for (i = firstDoCalibFloat; i < Nfloat; ++i)
    {
        float_start_dates[i - firstDoCalibFloat] = float_start_dates[i];
        float_end_dates[i - firstDoCalibFloat]   = float_end_dates[i];
        floatpaydates[i - firstDoCalibFloat]     = floatpaydates[i];

        float_margins[i - firstDoCalibFloat]   = float_margins[i];
        float_spreads[i - firstDoCalibFloat]   = float_spreads[i];
        float_notionals[i - firstDoCalibFloat] = float_notionals[i];
        float_cvgs[i - firstDoCalibFloat]      = float_cvgs[i];
    }

    float_pay_times    = dvector(0, Nfloat);
    float_pay_dates    = lngvector(0, Nfloat);
    float_pay_times[0] = (float_start_dates[0] - today) * YEARS_IN_DAY;
    float_pay_dates[0] = float_start_dates[0];
    for (n = 1; n <= Nfloat; ++n)
    {
        float_pay_times[n] = (floatpaydates[n - 1] - today) * YEARS_IN_DAY;
        float_pay_dates[n] = floatpaydates[n - 1];
    }

    fix_pay_times    = dvector(0, Nfix);
    fix_pay_dates    = lngvector(0, Nfix);
    fix_pay_times[0] = (fix_start_dates[0] - today) * YEARS_IN_DAY;
    fix_pay_dates[0] = fix_start_dates[0];
    for (n = 1; n <= Nfix; ++n)
    {
        fix_pay_times[n] = (fixpaydates[n - 1] - today) * YEARS_IN_DAY;
        fix_pay_dates[n] = fixpaydates[n - 1];
    }

    //----------Rescale Notionals---------------------
    maxNotional = fabs(fix_notionals[0]);
    for (i = 1; i < Nfix - firstDoCalib; ++i)
    {
        if (maxNotional < fabs(fix_notionals[i]))
        {
            maxNotional = fabs(fix_notionals[i]);
        }
    }
    for (i = 0; i < Nfloat - firstDoCalibFloat; ++i)
    {
        if (maxNotional < fabs(float_notionals[i]))
        {
            maxNotional = fabs(float_notionals[i]);
        }
    }

    CALPRESNOT = maxNotional * CALPRESAMORTMIDAT;

    //-----------------------------------------------------
    //-------------Build Amort Swap Coupons----------------
    //-----------------------------------------------------

    nbcoupon = Nfloat + 1 - firstDoCalibFloat;

    coupon_times = (double*)calloc(nbcoupon, sizeof(double));
    memcpy(coupon_times, float_pay_times, nbcoupon * sizeof(double));

    for (k = 0; k < nbcoupon; ++k)
    {
        coupon[k]       = 0;
        floatcoupon[k]  = 0;
        swaption_cpn[k] = 0;
    }

    nbexercise  = Nfix - firstDoCalib;
    float_index = 1;
    fix_index   = 1;
    coupon[0]   = -float_notionals[0] * swp_f_df(today, float_start_dates[0], yc_name);
    //	swaption_cpn[0] = - swp_f_df (today, float_start_dates[0], yc_name);
    floatcoupon[0] = -swp_f_df(today, float_start_dates[0], yc_name);
    coupon_time[0] = coupon_times[0];
    coupon_df[0]   = swp_f_df(today, float_start_dates[0], yc_name);
    coupon_cvg[0]  = 0;
    for (k = 1; k < nbcoupon; ++k)
    {
        coupon_time[k] = coupon_times[k];
        coupon_df[k]   = swp_f_df(today, floatpaydates[k - 1], yc_name);
        coupon_cvg[k]  = float_cvgs[k - 1];
        if ((fabs(coupon_times[k] - fix_pay_times[fix_index]) < 10.0 / 365.0) &&
            (fabs(coupon_times[k] - float_pay_times[float_index]) > 10.0 / 365.0))
        {
            coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                         fix_cvgs[fix_index - 1] *
                         swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
            ex_coupon[fix_index - 1] = k - fix_float_mult;
            exer_date[fix_index - 1] =
                add_unit(fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
            exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
            fix_index                = fix_index + 1;
        }
        else if (
            (fabs(coupon_times[k] - fix_pay_times[fix_index]) > 10.0 / 365.0) &&
            (fabs(coupon_times[k] - float_pay_times[float_index]) < 10.0 / 365.0))
        {
            temp = swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            coupon[k] += -(float_notionals[float_index] - float_notionals[float_index - 1] +
                           float_notionals[float_index - 1] *
                               (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                         swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            floatcoupon[k] += -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                              float_cvgs[float_index - 1] *
                              swp_f_df(today, floatpaydates[float_index - 1], yc_name);
            float_index = float_index + 1;
        }
        else
        {
            if (float_index < nbcoupon - 1)
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(float_notionals[float_index] - float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                    float_cvgs[float_index - 1] *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(-float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(-1.0 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
        }
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Amort Swap Coupon-------------------
    //-----------------------------------------------------------------

    //	Underlyings

    //	Long

    dff = swp_f_df(today, floatpaydates[Nfloat - 1 - firstDoCalibFloat], yc_name);
    for (i = 0; i < nbexercise; i++)
    {
        ex_lstrike[i] = fix_rates[i];
        ex_lprice[i]  = diag_prices[i];
    }

    dIRRs = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRs)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //-----------------------------------------------------------------
    //----------------Compute Diagonal IRRs of the amort Swap----------
    //-----------------------------------------------------------------
    if ((strike_type == 4) || (strike_type == 5) || (strike_type == 6) || (strike_type == 7))
    {
        err = ComputeAmortSwapDiagonalIRRsNew2(
            yc_name,
            vol_curve_name,
            default_ref_rate_name,
            fixfreq,
            fixbasis,

            Nfix - firstDoCalib,
            fix_start_dates,
            fix_end_dates,
            fix_cvgs,
            fix_rates,
            fix_notionals,
            ex_fee,

            Nfloat - firstDoCalibFloat,
            float_start_dates,
            float_end_dates,
            float_cvgs,
            float_margins,
            float_spreads,
            float_notionals,

            dIRRs,
            UseVol);
        if (err)
        {
            goto FREE_RETURN;
        }
    }

    dIRRstrikes = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRstrikes)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //--------------------------------------------------------------
    //----------Secondary instruments for Lambda calibration --------
    //--------------------------------------------------------------
    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute equivalent swaption strike using diagonal irrs-----
            else if ((strike_type == 6) || (strike_type == 7))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    fix_end_dates[Nfix - firstDoCalib - 1],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute diagonal swap, level and vol---------------------
            if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
            {
                lvl = 0.0;
                for (k = i; k < Nfix - firstDoCalib; k++)
                {
                    lvl += coupon_cvg_standard[k] * swp_f_df(today, fix_end_dates[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_end_dates[k - 1], yc_name);

                ex_slvl[i] = lvl;

                ex_sfwd[i] = (dfi - dff) / lvl;

                //	ATM std
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_end_dates[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute short swap, level and vol---------------------
            else
            {
                lvl = 0.0;
                for (k = i; k < i + 1; k++)
                {
                    lvl += fix_cvgs[k] * swp_f_df(today, fix_end_dates[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_end_dates[k - 1], yc_name);

                ex_slvl[i] = lvl;
                ex_sfwd[i] = (dfi - dff) / lvl;

                //	ATM std
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_end_dates[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            //---------------------------------------------------------------
            //---------Compute effective strike of secondary instrument -----
            //---------------------------------------------------------------

            //-----------------------Compute std-----------------------------
            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(exer_time[i]);

            //-----------ATM-----------------------------
            if (strike_type == 0)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            //-----------Custom Cash Strike--------------
            else if ((strike_type == 1) || (strike_type == 8))
            {
                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i];
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom Swap Strike -------------
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        fix_start_dates[i],
                        fix_end_dates[k - 1],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        default_ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i] - spr;
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom STD Strike --------------
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }
            //-----------IRR Equivalent Strike ----------
            else if (
                (strike_type == 4) || (strike_type == 5) || (strike_type == 6) ||
                (strike_type == 7))
            {
                ex_sstrike[i] = dIRRstrikes[i];
            }

            //----------Apply max std--------------------
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            //-----	Make sure strikes are positive (actually more than 1bp)--
            //-----	otherwise use ATM	-------------------------------------
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            calibstrike = ex_sstrike[i];

            //--------Price secondary instrument---------------
            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_end_dates[k - 1],
                calibstrike,
                0,
                default_ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            cap_price += ex_sprice[i];
            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }
    }

    //------------------------------------------------------------
    //--------- Fill Standard Swaption Coupons--------------------
    //------------------------------------------------------------
    if ((strike_type >= 6) && (!fix_lambda))
    {
        swaption_cpn[0] = -swp_f_df(today, float_start_dates[0], yc_name);
        float_index     = 1;
        fix_index       = 1;
        for (k = 1; k < nbcoupon; ++k)
        {
            if ((coupon_times[k] == fix_pay_times[fix_index]) &&
                (coupon_times[k] != float_pay_times[float_index]))
            {
                swaption_cpn[k] += ex_sstrike[fix_index - 1] * coupon_cvg_standard[fix_index - 1] *
                                   swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_times[k] != fix_pay_times[fix_index]) &&
                (coupon_times[k] == float_pay_times[float_index]))
            {
                swaption_cpn[k] +=
                    -((float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                      float_cvgs[float_index - 1]) *
                    swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < nbcoupon - 1)
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] *
                                       coupon_cvg_standard[fix_index - 1] *
                                       swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                    fix_index   = fix_index + 1;
                    float_index = float_index + 1;
                }
                else
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                                       swp_f_df(today, fixpaydates[fix_index - 1], yc_name);
                    fix_index = fix_index + 1;
                    swaption_cpn[k] +=
                        -(-1 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                        swp_f_df(today, floatpaydates[float_index - 1], yc_name);
                    float_index = float_index + 1;
                }
            }
        }
    }
    //------------------------------------------------------------
    //--------- End Of Fill Standard Swaption Coupons-------------
    //------------------------------------------------------------

    //	The 1F equivalent case
    if (one2F == 2 && fix_lambda && one_f_equi && (strike_type < 6))
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute short swap rate, level and vol----------
            lvl = 0.0;
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_end_dates[k + 1], yc_name);
            }
            dfi = swp_f_df(today, fix_end_dates[i], yc_name);
            dff = swp_f_df(today, fix_end_dates[k], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            //--------Compute effective short strike--------------

            //--------Equivalent IRR strike-----------------------
            if ((strike_type == 4) || (strike_type == 5))
            {
                ex_sstrike[i] = dIRRstrikes[i];
                if (dIRRstrikes[i] < 0.0001)
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //--------ATM-----------------------------------------
            else
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }

        //--------Calibrate LGM1F vol TS with given lambda ---------
        err = amortMidat_lgmcalibzetalambda_ts(
            nbcoupon,
            coupon,
            coupon_time,
            coupon_df,
            coupon_cvg,
            nbexercise,
            exer_bool,
            exer_time,
            ex_coupon,
            ex_sncpn,
            ex_lprice,
            ex_fee,
            shortcoupon,
            floatcoupon,
            0.0,
            ex_zeta,
            float_notionals,
            swaption_cpn,
            1,
            0,
            nNum_lambda,
            pdLambdaValue,
            pdLambdaDate,
            1,
            0.0,
            0.0,
            0.0,
            skip_last);

        if (err)
        {
            goto FREE_RETURN;
        }

        //--------Price the cap within the calibrated LGM1F model----------------
        export_lgmsetupG_ts(
            nNum_lambda,
            pdLambdaValue,
            pdLambdaDate,
            nbcoupon,
            coupon_time,
            cpn_G,
            nbexercise,
            exer_time,
            ex_G);

        cap_price = amortMidat_lgmcapval1F_b(
            nbcoupon,
            coupon,
            coupon_df,
            coupon_cvg,
            cpn_G,
            nbexercise,
            shortcoupon,
            floatcoupon,
            ex_coupon,
            ex_sncpn,
            ex_zeta,
            ex_G);

        fix_lambda = 0;
    }

    //--------Calibrate lambda and zeta----------------
    err = amortMidat_lgmcalibzetalambda_ts(
        nbcoupon,
        coupon,
        coupon_time,
        coupon_df,
        coupon_cvg,
        nbexercise,
        exer_bool,
        exer_time,
        ex_coupon,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        shortcoupon,
        floatcoupon,
        cap_price,
        ex_zeta,
        float_notionals,
        swaption_cpn,
        fix_lambda,
        cap_or_swaption,
        nNum_lambda,
        pdLambdaValue,
        pdLambdaDate,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    //	3.)	Transform into sigma

    compt = 0;
    for (i = 0; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            compt += 1;
        }
    }

    *num_sig  = compt;
    *sig_time = (double*)calloc(*num_sig, sizeof(double));
    *sig      = (double*)calloc(*num_sig, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    (*sig_time)[0] = exer_time[0];

    //// previously
    //// (*sig)[0] = sqrt ( ex_zeta[0] * 2 * (*lambda) / (exp (2 * (*lambda) * exer_time[0]) - 1.0)
    ///);

    exp_fact =
        static_lgmcalcexpfact_tauts(0.0, exer_time[0], nNum_lambda, pdLambdaDate, pdLambdaValue);

    (*sig)[0] = sqrt(ex_zeta[0] / exp_fact);

    last_i = 0;
    k      = 1;
    for (i = 1; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            (*sig_time)[k] = exer_time[i];
            if (ex_zeta[i] > ex_zeta[last_i])
            {
                /// previously
                ///(*sig)[k] = sqrt ((ex_zeta[i] - ex_zeta[last_i]) * 2 * (*lambda)
                ///				/ (exp (2 * (*lambda) * exer_time[i]) - exp ( 2 * (*lambda) *
                ///exer_time[last_i])));

                exp_fact = static_lgmcalcexpfact_tauts(
                    exer_time[last_i], exer_time[i], nNum_lambda, pdLambdaDate, pdLambdaValue);

                (*sig)[k] = sqrt((ex_zeta[i] - ex_zeta[last_i]) / exp_fact);
            }
            else
            {
                smessage(
                    "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                    exer_time[i]);
                for (j = k; j < *num_sig; j++)
                {
                    (*sig_time)[j] = exer_time[k];
                    (*sig)[j]      = (*sig)[k - 1];
                }
                i = nbexercise;
            }

            last_i = i;
            k      = k + 1;
        }
    }

FREE_RETURN:

    if (fix_cvgs)
        free_dvector(fix_cvgs, 0, Nfix - 1);
    if (float_cvgs)
        free_dvector(float_cvgs, 0, Nfloat - 1);

    if (float_pay_dates)
        free_lngvector(float_pay_dates, 0, Nfloat);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, Nfloat);

    if (fix_pay_dates)
        free_lngvector(fix_pay_dates, 0, Nfix);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, Nfix);

    if (coupon_times)
        free(coupon_times);

    if (dIRRs)
        free(dIRRs);
    if (dIRRstrikes)
        free(dIRRstrikes);

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;
    }

    return err;
}  //// end of amortMidat_cpd_calib_diagonal_new_ts

/*	Calibrate lgm: main function */
Err amortMidat_cpd_calib_diagonal_new2(
    int   notperiod,
    char* yc_name,        /*	Name of the yield curve */
    char* vol_curve_name, /*	Name of the market vol curve */

    char* default_ref_rate_name, /*	Name of the default reference rate */
    char* swaption_basis,        /*	 */
    char* swaption_freq,         /*	 */
    Err (*get_cash_vol)(         /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    double vol_shift,
    int    shift_type, /*	0:	Additive
                                       1:	Multiplicative */

    int*    ex_bool,     /*	Exercise or not : NULL = True */
    double* diag_prices, /* Diagonal Prices */
    double* ex_fee,      /* Exercise Fees of Diagonal Swaptions */

    char*   fix_freq,
    char*   fix_basis,
    int     Nfix,
    long*   fix_start_dates, /*	Start date of the amortised swap */
    long*   fix_end_dates,   /*	End date for the amortised swap */
    long*   fix_pay_dates,   /*	Pay date for the amortised swap */
    double* fix_rates,       /*	Diagonal swaption strikes */
    double* fix_notionals,   /* Amortised Notional */

    char*   float_freq,
    char*   float_basis,
    int     Nfloat,
    long*   float_start_dates, /*	Start date of the amortised swap */
    long*   float_end_dates,   /*	End date for the amortised swap */
    long*   float_pay_dates,   /*	Pay date for the amortised swap */
    double* float_margins,     /*	Margins */
    double* float_spreads,     /*	Spreads */
    double* float_notionals,   /*	Amortised Notional */

    /* Calibration Parameters*/
    double* short_strike, /*	For Calibration : Short swaption strikes NULL = ATM */
    int     strike_type,  /*	0: ATM
                                                  1: CASH
                                                  2: SWAP
                                                  3: STD
                                                  4: IRR Flat
                                                  5: IRR with vol
                                                  6: SWAPTION IRR strike
                                                  7: SWAPTION IRR vol strike
                                                  8: SWAPTION short_strike*/

    double max_std_short,

    int fix_lambda,      /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    int one_f_equi,      /*	1F equivalent flag:
                                                 if set to 1, then 2F lambda will calibrate
                                                 to the cap priced within calibrated 1F
                                                 with the given lambda */
    int skip_last,       /*	If 1, the last option is disregarded
                                                 and the forward volatility is flat from option
                                                 n-1 */
    int    use_jump,     /*	Allow Jump (and vol = 0) in vol TS */
    double max_var_jump, /*	Maximum multiplicative variation in variance */

    double* lambda, /*	Lambda: may be changed in the process */
    int     one2F,  /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double alpha,
    double gamma,
    double rho,

    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig)
{
    int            i, j, k, n, float_index, fix_index;
    int            nbcoupon;
    SrtCompounding fixfreq, floatfreq;
    SrtBasisCode   fixbasis, floatbasis, swapbasis;
    double         ex_lstrike[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_sncpn[MAX_CPN];
    double      cpn_G[MAX_CPN];
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;
    double      coupon[MAX_CPN];
    double      swaption_cpn[MAX_CPN];
    double      floatcoupon[MAX_CPN];
    double      shortcoupon[MAX_CPN];
    double      coupon_time[MAX_CPN];
    double      coupon_df[MAX_CPN];
    double      coupon_cvg[MAX_CPN];
    double      coupon_cvg_standard[MAX_CPN];
    int         ex_coupon[MAX_CPN];

    double temp;

    int fix_float_mult;

    int last_i;

    double* dIRRs       = NULL;
    double* dIRRstrikes = NULL;

    double *float_cvgs = NULL, *float_pay_times = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;
    long *  float_pay_datesOut = NULL, *fix_pay_datesOut = NULL;

    double* coupon_times = NULL;

    double maxNotional;

    int    nbexercise;
    double exer_time[MAX_CPN];
    long   exer_date[MAX_CPN];
    int    exer_bool[MAX_CPN];

    int firstDoCalib, firstDoCalibFloat;

    int cap_or_swaption;

    int compt;

    double calibstrike;

    int UseVol;

    //---------------------------------------
    //---------------------------------------

    fix_pay_datesOut = (long*)calloc(Nfix, sizeof(long));
    if (!fix_pay_datesOut)
    {
        err = "memory allocation failed in AmortMidatCalib";
        smessage("memory allocation failed in AmortMidatCalib");
        goto FREE_RETURN;
    }

    float_pay_datesOut = (long*)calloc(Nfloat, sizeof(long));
    if (!float_pay_datesOut)
    {
        err = "memory allocation failed in AmortMidatCalib";
        smessage("memory allocation failed in AmortMidatCalib");
        goto FREE_RETURN;
    }

    err = amortMidat_modify_dates2(
        Nfix, fix_pay_dates, Nfloat, float_pay_dates, fix_pay_datesOut, float_pay_datesOut);

    //-----If no lambda calibration : change strike_type---------------
    if (fix_lambda)
    {
        strike_type = 0;
    }

    //-----For IRR equivalent strike : Use vol or not ?---------
    UseVol = 0;
    if ((strike_type == 5) || (strike_type == 7))
    {
        UseVol = 1;
    }

    //-----Secondary instrument : Short or Diag Swaption ?---------
    cap_or_swaption = 0;
    if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
    {
        cap_or_swaption = 1;
    }

    //-----Parameter controlling vol jump size---------
    MAX_FACT  = max_var_jump;
    USE_JUMPS = use_jump;

    //-----Initialize outpu pointers-----------------
    *sig_time = NULL;
    *sig      = NULL;

    //-----Retreive yield curve from its name---------
    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    //-----If LGM2F compute Hermite quadrature weights and ponits---------
    if (one2F == 2)
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    //-----convert char in Sort structure---------
    err = interp_basis(fix_basis, &fixbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(fix_freq, &fixfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(float_basis, &floatbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_compounding(float_freq, &floatfreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----convert char in Sort structure---------
    err = interp_basis(swaption_basis, &swapbasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    //-----get freq and basis corresponding to ref_rate---------
    //	err = swp_f_get_ref_rate_details(ref_rate, &floatbasis, &floatfreq);

    //-------------Compute Coverages--------------
    fix_cvgs = dvector(0, Nfix - 1);
    for (i = 0; i < Nfix; ++i)
    {
        fix_cvgs[i] = coverage(fix_start_dates[i], fix_end_dates[i], fixbasis);
    }
    float_cvgs = dvector(0, Nfloat - 1);
    for (i = 0; i < Nfloat; ++i)
    {
        float_cvgs[i] = coverage(float_start_dates[i], float_end_dates[i], floatbasis);
    }
    fix_float_mult = floatfreq / fixfreq;

    //-------------Find First Exercise Date----------------
    if (add_unit(fix_start_dates[0], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) <= today)
    {
        ex_bool[0] = 0;
    }

    i = 0;
    while ((ex_bool[i] == 0) && (i < Nfix))
    {
        ++i;
    }

    if (i == Nfix)
    {
        err = "At least one maturity must be calibrated : Check DoCalib Value";
        goto FREE_RETURN;
    }

    firstDoCalib = i;

    //--------------Shift all data---------------------------
    for (i = firstDoCalib; i < Nfix; ++i)
    {
        ex_bool[i - firstDoCalib]     = ex_bool[i];
        exer_bool[i - firstDoCalib]   = ex_bool[i];
        diag_prices[i - firstDoCalib] = diag_prices[i];
        ex_fee[i - firstDoCalib]      = ex_fee[i];

        fix_start_dates[i - firstDoCalib]  = fix_start_dates[i];
        fix_end_dates[i - firstDoCalib]    = fix_end_dates[i];
        fix_pay_datesOut[i - firstDoCalib] = fix_pay_datesOut[i];
        fix_rates[i - firstDoCalib]        = fix_rates[i];
        fix_notionals[i - firstDoCalib]    = fix_notionals[i];
        fix_cvgs[i - firstDoCalib]         = fix_cvgs[i];

        if (short_strike)
        {
            short_strike[i - firstDoCalib] = short_strike[i];
        }

        coupon_cvg_standard[i - firstDoCalib] =
            coverage(fix_start_dates[i - firstDoCalib], fix_end_dates[i - firstDoCalib], swapbasis);
    }

    i = 0;
    while (float_start_dates[i] < fix_start_dates[0] - 10)
    {
        ++i;
    }

    firstDoCalibFloat = i;

    for (i = firstDoCalibFloat; i < Nfloat; ++i)
    {
        float_start_dates[i - firstDoCalibFloat] = float_start_dates[i];
        float_end_dates[i - firstDoCalibFloat]   = float_end_dates[i];

        float_pay_datesOut[i - firstDoCalibFloat] = float_pay_datesOut[i];
        float_margins[i - firstDoCalibFloat]      = float_margins[i];
        float_spreads[i - firstDoCalibFloat]      = float_spreads[i];
        float_notionals[i - firstDoCalibFloat]    = float_notionals[i];
        float_cvgs[i - firstDoCalibFloat]         = float_cvgs[i];
    }

    float_pay_times    = dvector(0, Nfloat);
    float_pay_times[0] = (float_start_dates[0] - today) * YEARS_IN_DAY;
    for (n = 1; n <= Nfloat; ++n)
    {
        float_pay_times[n] = (float_pay_datesOut[n - 1] - today) * YEARS_IN_DAY;
    }

    fix_pay_times    = dvector(0, Nfix);
    fix_pay_times[0] = (fix_start_dates[0] - today) * YEARS_IN_DAY;
    for (n = 1; n <= Nfix; ++n)
    {
        fix_pay_times[n] = (fix_pay_datesOut[n - 1] - today) * YEARS_IN_DAY;
    }

    //----------Rescale Notionals---------------------
    maxNotional = fabs(fix_notionals[0]);
    for (i = 1; i < Nfix - firstDoCalib; ++i)
    {
        if (maxNotional < fabs(fix_notionals[i]))
        {
            maxNotional = fabs(fix_notionals[i]);
        }
    }
    for (i = 0; i < Nfloat - firstDoCalibFloat; ++i)
    {
        if (maxNotional < fabs(float_notionals[i]))
        {
            maxNotional = fabs(float_notionals[i]);
        }
    }

    CALPRESNOT = maxNotional * CALPRESAMORTMIDAT;

    //-----------------------------------------------------
    //-------------Build Amort Swap Coupons----------------
    //-----------------------------------------------------

    nbcoupon = Nfloat + 1 - firstDoCalibFloat;

    coupon_times = (double*)calloc(nbcoupon, sizeof(double));
    memcpy(coupon_times, float_pay_times, nbcoupon * sizeof(double));

    for (k = 0; k < nbcoupon; ++k)
    {
        coupon[k]       = 0;
        floatcoupon[k]  = 0;
        swaption_cpn[k] = 0;
    }

    nbexercise     = Nfix - firstDoCalib;
    float_index    = 1;
    fix_index      = 1;
    coupon[0]      = -float_notionals[0] * swp_f_df(today, float_start_dates[0], yc_name);
    floatcoupon[0] = -swp_f_df(today, float_start_dates[0], yc_name);
    coupon_time[0] = coupon_times[0];
    coupon_df[0]   = swp_f_df(today, float_start_dates[0], yc_name);
    coupon_cvg[0]  = 0;
    for (k = 1; k < nbcoupon; ++k)
    {
        coupon_time[k] = coupon_times[k];
        coupon_df[k]   = swp_f_df(today, float_pay_datesOut[k - 1], yc_name);
        coupon_cvg[k]  = float_cvgs[k - 1];
        if ((coupon_times[k] == fix_pay_times[fix_index]) &&
            (coupon_times[k] != float_pay_times[float_index]))
        {
            coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                         fix_cvgs[fix_index - 1] *
                         swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
            ex_coupon[fix_index - 1] = k - fix_float_mult;
            exer_date[fix_index - 1] =
                add_unit(fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
            exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
            fix_index                = fix_index + 1;
        }
        else if (
            (coupon_times[k] != fix_pay_times[fix_index]) &&
            (coupon_times[k] == float_pay_times[float_index]))
        {
            temp = swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
            coupon[k] += -(float_notionals[float_index] - float_notionals[float_index - 1] +
                           float_notionals[float_index - 1] *
                               (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                         swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
            floatcoupon[k] += -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                              float_cvgs[float_index - 1] *
                              swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
            float_index = float_index + 1;
        }
        else
        {
            if (float_index < nbcoupon - 1)
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(float_notionals[float_index] - float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                    float_cvgs[float_index - 1] *
                    swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                coupon[k] += fix_rates[fix_index - 1] * fix_notionals[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_start_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] +=
                    -(-float_notionals[float_index - 1] +
                      float_notionals[float_index - 1] *
                          (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                          float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                floatcoupon[k] +=
                    -(-1.0 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
        }
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Amort Swap Coupon-------------------
    //-----------------------------------------------------------------

    /*	Underlyings */

    /*	Long */

    dff = swp_f_df(today, float_pay_datesOut[Nfloat - 1 - firstDoCalibFloat], yc_name);
    for (i = 0; i < nbexercise; i++)
    {
        ex_lstrike[i] = fix_rates[i];
        ex_lprice[i]  = diag_prices[i];
    }

    dIRRs = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRs)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //-----------------------------------------------------------------
    //----------------Compute Diagonal IRRs of the amort Swap----------
    //-----------------------------------------------------------------
    if ((strike_type == 4) || (strike_type == 5) || (strike_type == 6) || (strike_type == 7))
    {
        err = ComputeAmortSwapDiagonalIRRsNew2(
            yc_name,
            vol_curve_name,
            default_ref_rate_name,
            fixfreq,
            fixbasis,

            Nfix - firstDoCalib,
            fix_start_dates,
            fix_end_dates,
            fix_cvgs,
            fix_rates,
            fix_notionals,
            ex_fee,

            Nfloat - firstDoCalibFloat,
            float_start_dates,
            float_end_dates,
            float_cvgs,
            float_margins,
            float_spreads,
            float_notionals,

            dIRRs,
            UseVol);
        if (err)
        {
            goto FREE_RETURN;
        }
    }

    dIRRstrikes = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRstrikes)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    //--------------------------------------------------------------
    //----------Secondary instruments for Lambda calibration --------
    //--------------------------------------------------------------
    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute equivalent swaption strike using diagonal irrs-----
            else if ((strike_type == 6) || (strike_type == 7))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    fix_end_dates[Nfix - firstDoCalib - 1],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute diagonal swap, level and vol---------------------
            if ((strike_type == 6) || (strike_type == 7) || (strike_type == 8))
            {
                lvl = 0.0;
                for (k = i; k < Nfix - firstDoCalib; k++)
                {
                    lvl += coupon_cvg_standard[k] * swp_f_df(today, fix_pay_datesOut[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_pay_datesOut[k - 1], yc_name);

                ex_slvl[i] = lvl;

                ex_sfwd[i] = (dfi - dff) / lvl;

                /*	ATM std */
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_pay_datesOut[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            //---------Compute short swap, level and vol---------------------
            else
            {
                lvl = 0.0;
                for (k = i; k < i + 1; k++)
                {
                    lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_datesOut[k], yc_name);
                }
                dfi = swp_f_df(today, fix_start_dates[i], yc_name);
                dff = swp_f_df(today, fix_pay_datesOut[k - 1], yc_name);

                ex_slvl[i] = lvl;
                ex_sfwd[i] = (dfi - dff) / lvl;

                /*	ATM std */
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_pay_datesOut[k - 1],
                    ex_sfwd[i],
                    0,
                    default_ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            //---------------------------------------------------------------
            //---------Compute effective strike of secondary instrument -----
            //---------------------------------------------------------------

            //-----------------------Compute std-----------------------------
            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(exer_time[i]);

            //-----------ATM-----------------------------
            if (strike_type == 0)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            //-----------Custom Cash Strike--------------
            else if ((strike_type == 1) || (strike_type == 8))
            {
                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i];
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom Swap Strike -------------
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        fix_start_dates[i],
                        fix_pay_datesOut[k - 1],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        default_ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                if (short_strike)
                {
                    ex_sstrike[i] = short_strike[i] - spr;
                }
                else
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //-----------Custom STD Strike --------------
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }
            //-----------IRR Equivalent Strike ----------
            else if (
                (strike_type == 4) || (strike_type == 5) || (strike_type == 6) ||
                (strike_type == 7))
            {
                ex_sstrike[i] = dIRRstrikes[i];
            }

            //----------Apply max std--------------------
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            //-----	Make sure strikes are positive (actually more than 1bp)--
            //-----	otherwise use ATM	-------------------------------------
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            calibstrike = ex_sstrike[i];

            //--------Price secondary instrument---------------
            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_pay_datesOut[k - 1],
                calibstrike,
                0,
                default_ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            cap_price += ex_sprice[i];
            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }
    }

    //------------------------------------------------------------
    //--------- Fill Standard Swaption Coupons--------------------
    //------------------------------------------------------------
    if ((strike_type >= 6) && (!fix_lambda))
    {
        swaption_cpn[0] = -swp_f_df(today, float_start_dates[0], yc_name);
        float_index     = 1;
        fix_index       = 1;
        for (k = 1; k < nbcoupon; ++k)
        {
            if ((coupon_times[k] == fix_pay_times[fix_index]) &&
                (coupon_times[k] != float_pay_times[float_index]))
            {
                swaption_cpn[k] += ex_sstrike[fix_index - 1] * coupon_cvg_standard[fix_index - 1] *
                                   swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
                fix_index = fix_index + 1;
            }
            else if (
                (coupon_times[k] != fix_pay_times[fix_index]) &&
                (coupon_times[k] == float_pay_times[float_index]))
            {
                swaption_cpn[k] +=
                    -((float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                      float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                if (float_index < nbcoupon - 1)
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] *
                                       coupon_cvg_standard[fix_index - 1] *
                                       swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
                    fix_index   = fix_index + 1;
                    float_index = float_index + 1;
                }
                else
                {
                    swaption_cpn[k] += ex_sstrike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                                       swp_f_df(today, fix_pay_datesOut[fix_index - 1], yc_name);
                    fix_index = fix_index + 1;
                    swaption_cpn[k] +=
                        -(-1 + (float_spreads[float_index - 1] + float_margins[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                        swp_f_df(today, float_pay_datesOut[float_index - 1], yc_name);
                    float_index = float_index + 1;
                }
            }
        }
    }
    //------------------------------------------------------------
    //--------- End Of Fill Standard Swaption Coupons-------------
    //------------------------------------------------------------

    /*	The 1F equivalent case */
    if (one2F == 2 && fix_lambda && one_f_equi && (strike_type < 6))
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            //---------Compute equivalent caplet strike using diagonal irrs-----
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    fixfreq,
                    fixbasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = Nfloat - firstDoCalibFloat - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            //---------Compute short swap rate, level and vol----------
            lvl = 0.0;
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_datesOut[k], yc_name);
            }
            dfi = swp_f_df(today, fix_pay_dates[i - 1], yc_name);
            dff = swp_f_df(today, fix_pay_dates[k - 1], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            //--------Compute effective short strike--------------

            //--------Equivalent IRR strike-----------------------
            if ((strike_type == 4) || (strike_type == 5))
            {
                ex_sstrike[i] = dIRRstrikes[i];
                if (dIRRstrikes[i] < 0.0001)
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            //--------ATM-----------------------------------------
            else
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }

        //--------Calibrate LGM1F vol TS with given lambda ---------
        err = amortMidat_lgmcalibzetalambda(
            nbcoupon,
            coupon,
            coupon_time,
            coupon_df,
            coupon_cvg,
            nbexercise,
            exer_bool,
            exer_time,
            ex_coupon,
            ex_sncpn,
            ex_lprice,
            ex_fee,
            shortcoupon,
            floatcoupon,
            0.0,
            ex_zeta,
            float_notionals,
            swaption_cpn,
            1,
            0,
            lambda,
            1,
            0.0,
            0.0,
            0.0,
            skip_last);

        if (err)
        {
            goto FREE_RETURN;
        }

        //--------Price the cap within the calibrated LGM1F model----------------
        export_lgmsetupG(*lambda, nbcoupon, coupon_time, cpn_G, nbexercise, exer_time, ex_G);

        cap_price = amortMidat_lgmcapval1F_b(
            nbcoupon,
            coupon,
            coupon_df,
            coupon_cvg,
            cpn_G,
            nbexercise,
            shortcoupon,
            floatcoupon,
            ex_coupon,
            ex_sncpn,
            ex_zeta,
            ex_G);

        fix_lambda = 0;
    }

    //--------Calibrate lambda and zeta----------------
    err = amortMidat_lgmcalibzetalambda(
        nbcoupon,
        coupon,
        coupon_time,
        coupon_df,
        coupon_cvg,
        nbexercise,
        exer_bool,
        exer_time,
        ex_coupon,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        shortcoupon,
        floatcoupon,
        cap_price,
        ex_zeta,
        float_notionals,
        swaption_cpn,
        fix_lambda,
        cap_or_swaption,
        lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	3.)	Transform into sigma */

    compt = 0;
    for (i = 0; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            compt += 1;
        }
    }

    *num_sig  = compt;
    *sig_time = (double*)calloc(*num_sig, sizeof(double));
    *sig      = (double*)calloc(*num_sig, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    (*sig_time)[0] = exer_time[0];
    (*sig)[0]      = sqrt(ex_zeta[0] * 2 * (*lambda) / (exp(2 * (*lambda) * exer_time[0]) - 1.0));

    last_i = 0;
    k      = 1;
    for (i = 1; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            (*sig_time)[k] = exer_time[i];
            if (ex_zeta[i] > ex_zeta[last_i])
            {
                (*sig)[k] = sqrt(
                    (ex_zeta[i] - ex_zeta[last_i]) * 2 * (*lambda) /
                    (exp(2 * (*lambda) * exer_time[i]) - exp(2 * (*lambda) * exer_time[last_i])));
            }
            else
            {
                smessage(
                    "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                    exer_time[i]);
                for (j = k; j < *num_sig; j++)
                {
                    (*sig_time)[j] = exer_time[k];
                    (*sig)[j]      = (*sig)[k - 1];
                }
                i = nbexercise;
            }

            last_i = i;
            k      = k + 1;
        }
    }

FREE_RETURN:

    if (fix_cvgs)
        free_dvector(fix_cvgs, 0, Nfix - 1);
    if (float_cvgs)
        free_dvector(float_cvgs, 0, Nfloat - 1);

    if (float_pay_datesOut)
        free(float_pay_datesOut);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, Nfloat);

    if (fix_pay_datesOut)
        free(fix_pay_datesOut);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, Nfix);

    if (coupon_times)
        free(coupon_times);

    if (dIRRs)
        free(dIRRs);
    if (dIRRstrikes)
        free(dIRRstrikes);

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;
    }

    return err;
}

/*	Calibrate lgm: main function */
Err amortMidat_cpd_calib_diagonal(
    int   notperiod,
    char* yc_name,               /*	Name of the yield curve */
    char* vol_curve_name,        /*	Name of the market vol curve */
    char* default_ref_rate_name, /*	Name of the reference rate */
    Err (*get_cash_vol)(         /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    double vol_shift,
    int    shift_type,     /*	0:	Additive
                                           1:	Multiplicative */
    int* ex_bool,          /*	Exercise or not
                                                   NULL = True */
    long    start_date,    /*	Start date of the amortised swap */
    long    end_date,      /*	End date for the amortised swap */
    double* long_strike,   /*	Diagonal swaption strikes
                                                                   NULL = ATM */
    double* short_strike,  /*	Short swaption strikes
                                                                   NULL = ATM */
    int strike_type,       /*	0: ATM
                                                   1: CASH
                                                   2: SWAP
                                                   3: STD
                                                   4: IRR Flat
                                                   5: IRR with vol*/
    double* diag_prices,   /* Diagonal Prices */
    double* ex_fee,        /* Exercise Fees of Diagonal Swaptions */
    double* fixNotional,   /* Amortised Notional
                                           NULL = No amortisation */
    double* floatNotional, /* Amortised Notional
                                           NULL = No amortisation */
    double* margin,

    double max_std_short,
    char*  ref_rate_name,
    char*  swaption_freq, /*	Frequency, basis and ref. rate of calibrated swaptions */
    char*  swaption_basis,
    int    fix_lambda,   /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    int one_f_equi,      /*	1F equivalent flag:
                                                 if set to 1, then 2F lambda will calibrate
                                                 to the cap priced within calibrated 1F
                                                 with the given lambda */
    int skip_last,       /*	If 1, the last option is disregarded
                                                 and the forward volatility is flat from option
                                                 n-1 */
    double max_var_jump, /*	Maximum multiplicative variation in variance */

    double* lambda, /*	Lambda: may be changed in the process */
    int     one2F,  /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   gamma,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data) /*	NULL = don't save calibration instrument data */
{
    int            i, j, k, ncpn, n, float_index, fix_index;
    int            nbcoupon;
    SrtCompounding ifreq;
    SrtBasisCode   ibasis;
    double         swaption_strikes[MAX_CPN];
    double         ex_lstrike[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_sncpn[MAX_CPN];
    double      cpn_G[MAX_CPN];
    long        theo_date, act_date;
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;
    double      coupon[MAX_CPN];
    double      swaption_cpn[MAX_CPN];
    double      floatcoupon[MAX_CPN];
    double      shortcoupon[MAX_CPN];
    double      coupon_time[MAX_CPN];
    double      coupon_df[MAX_CPN];
    double      coupon_cvg[MAX_CPN];
    int         ex_coupon[MAX_CPN];

    double temp;

    int tem;

    int fix_float_mult;

    int last_i;

    double* dIRRs       = NULL;
    double* dIRRstrikes = NULL;

    //---------------------------------------
    //------------SWAPDP---------------------
    //---------------------------------------
    SwapDP swapdp;
    long   float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    long    fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    double* coupon_dates = NULL;

    double maxNotional;

    int    nbexercise;
    double exer_time[MAX_CPN];
    long   exer_date[MAX_CPN];
    int    exer_bool[MAX_CPN];

    int firstDoCalib;

    int cap_or_swaption;

    int compt;

    double dtest;

    double calibstrike;

    int UseVol;
    //---------------------------------------
    //---------------------------------------

    if (fix_lambda)
    {
        strike_type = 0;
    }

    UseVol = 0;
    if ((strike_type == 5) || (strike_type == 6))
    {
        UseVol = 1;
    }
    MAX_FACT = max_var_jump;

    *sig_time = NULL;
    *sig      = NULL;

    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    if (one2F == 2)
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    /*	1.)	Setup the bond schedule and its coupons */

    /*	Coupons */

    err = interp_compounding(swaption_freq, &ifreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = interp_basis(swaption_basis, &ibasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    theo_date = end_date;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn      = 1;

    while (act_date > today)
    {
        theo_date = add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
        act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        ncpn++;
    }
    ncpn--;

    if (ncpn < 2)
    {
        err = "Not enough coupons in amortMidat_cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    theo_date = end_date;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    i         = ncpn - 1;

    //--------------------------------------------
    //-------------Find Start Date----------------
    //--------------------------------------------
    if (add_unit(start_date, -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) <= today)
    {
        ex_bool[0] = 0;
    }

    i = 0;
    while ((ex_bool[i] == 0) && (i < ncpn))
    {
        ++i;
    }

    if (i == ncpn)
    {
        err = "At least one maturity must be calibrated : Check DoCalib Value";
        goto FREE_RETURN;
    }

    firstDoCalib = i;

    //--------------------------------------------
    //--------- If Start = today -----------------
    //--------------------------------------------

    err = swp_f_initSwapDP(start_date, theo_date, swaption_freq, swaption_basis, &swapdp);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp,
        today,
        ref_rate_name,
        &float_pay_dates,
        &float_nb_pay_dates,
        &float_fixing_dates,
        &float_start_dates,
        &float_end_dates,
        &float_cvgs,
        &float_spreads,
        &float_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp,
        today,
        &fix_pay_dates,
        &fix_nb_pay_dates,
        &fix_start_dates,
        &fix_end_dates,
        &fix_cvgs,
        &fix_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }

    start_date = fix_start_dates[firstDoCalib];

    tem = firstDoCalib;

    for (i = firstDoCalib; i < fix_nb_dates; ++i)
    {
        ex_bool[i - firstDoCalib]     = ex_bool[i];
        exer_bool[i - firstDoCalib]   = ex_bool[i];
        fixNotional[i - firstDoCalib] = fixNotional[i];
        diag_prices[i - firstDoCalib] = diag_prices[i];
        long_strike[i - firstDoCalib] = long_strike[i];
        ex_fee[i - firstDoCalib]      = ex_fee[i];
        if (short_strike)
        {
            short_strike[i - firstDoCalib] = short_strike[i];
        }
    }

    for (i = (int)(0.5 + fix_cvgs[1] / float_cvgs[1]) * firstDoCalib; i < float_nb_dates; ++i)
    {
        floatNotional[i - (int)(0.5 + fix_cvgs[1] / float_cvgs[1]) * firstDoCalib] =
            floatNotional[i];
    }

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_cvgs)
        free(fix_cvgs);

    float_fixing_dates = NULL;
    float_start_dates  = NULL;
    float_end_dates    = NULL;
    float_pay_dates    = NULL;
    float_cvgs         = NULL;
    float_spreads      = NULL;
    float_pay_times    = NULL;

    fix_start_dates = NULL;
    fix_end_dates   = NULL;
    fix_pay_dates   = NULL;
    fix_cvgs        = NULL;
    fix_pay_times   = NULL;

    //--------------------------------------------
    //---------SWAPDP-----------------------------
    //---------Build Swap Schedule----------------
    //--------------------------------------------

    err = swp_f_initSwapDP(start_date, theo_date, swaption_freq, swaption_basis, &swapdp);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp,
        today,
        ref_rate_name,
        &float_pay_dates,
        &float_nb_pay_dates,
        &float_fixing_dates,
        &float_start_dates,
        &float_end_dates,
        &float_cvgs,
        &float_spreads,
        &float_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }
    float_pay_times = dvector(0, float_nb_pay_dates - 1);
    for (n = 0; n < float_nb_pay_dates; ++n)
    {
        float_pay_times[n] = (float_pay_dates[n] - today) * YEARS_IN_DAY;
    }

    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp,
        today,
        &fix_pay_dates,
        &fix_nb_pay_dates,
        &fix_start_dates,
        &fix_end_dates,
        &fix_cvgs,
        &fix_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }
    fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
    for (n = 0; n < fix_nb_pay_dates; ++n)
    {
        fix_pay_times[n] = (fix_pay_dates[n] - today) * YEARS_IN_DAY;
    }

    fix_float_mult = (int)(fix_cvgs[fix_nb_dates - 1] / float_cvgs[float_nb_dates - 1] + 0.5);

    //------------------------------------------------
    //----------Rescale Notionals---------------------
    //------------------------------------------------
    maxNotional = fabs(fixNotional[0]);
    for (i = 1; i < fix_nb_dates; ++i)
    {
        if (maxNotional < fabs(fixNotional[i]))
        {
            maxNotional = fabs(fixNotional[i]);
        }
    }
    for (i = 0; i < float_nb_dates; ++i)
    {
        if (maxNotional < fabs(floatNotional[i]))
        {
            maxNotional = fabs(floatNotional[i]);
        }
    }

    CALPRESNOT = maxNotional * CALPRESAMORTMIDAT;

    //-------------------------------------------
    //-------------------------------------------

    nbcoupon = float_nb_pay_dates;

    coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
    memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
    num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates, fix_pay_times);
    num_f_sort_vector(nbcoupon, coupon_dates);
    num_f_unique_vector(&nbcoupon, coupon_dates);

    for (k = 0; k < nbcoupon; ++k)
    {
        coupon[k]       = 0;
        floatcoupon[k]  = 0;
        swaption_cpn[k] = 0;
    }

    nbexercise      = fix_nb_pay_dates - 1;
    float_index     = 1;
    fix_index       = 1;
    coupon[0]       = -floatNotional[0] * swp_f_df(today, float_pay_dates[0], yc_name);
    swaption_cpn[0] = -swp_f_df(today, float_pay_dates[0], yc_name);
    floatcoupon[0]  = -swp_f_df(today, float_pay_dates[0], yc_name);
    coupon_time[0]  = coupon_dates[0];
    coupon_df[0]    = swp_f_df(today, float_pay_dates[0], yc_name);
    coupon_cvg[0]   = 0;
    for (k = 1; k < nbcoupon; ++k)
    {
        coupon_time[k] = coupon_dates[k];
        coupon_df[k]   = swp_f_df(today, float_pay_dates[k], yc_name);
        coupon_cvg[k]  = float_cvgs[k - 1];
        if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
            (coupon_dates[k] != float_pay_times[float_index]))
        {
            coupon[k] += long_strike[fix_index - 1] * fixNotional[fix_index - 1] *
                         fix_cvgs[fix_index - 1] *
                         swp_f_df(today, fix_pay_dates[fix_index], yc_name);
            swaption_cpn[k] += long_strike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                               swp_f_df(today, fix_pay_dates[fix_index], yc_name);
            ex_coupon[fix_index - 1] = k - fix_float_mult;
            exer_date[fix_index - 1] =
                add_unit(fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
            exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
            fix_index                = fix_index + 1;
        }
        else if (
            (coupon_dates[k] != fix_pay_times[fix_index]) &&
            (coupon_dates[k] == float_pay_times[float_index]))
        {
            temp = swp_f_df(today, float_pay_dates[float_index], yc_name);
            coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                           floatNotional[float_index - 1] *
                               (float_spreads[float_index - 1] + margin[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                         swp_f_df(today, float_pay_dates[float_index], yc_name);
            swaption_cpn[k] += -((float_spreads[float_index - 1] + margin[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                               swp_f_df(today, float_pay_dates[float_index], yc_name);
            floatcoupon[k] += -(float_spreads[float_index - 1] + margin[float_index - 1]) *
                              float_cvgs[float_index - 1] *
                              swp_f_df(today, float_pay_dates[float_index], yc_name);
            float_index = float_index + 1;
        }
        else
        {
            if (float_index < float_nb_pay_dates - 1)
            {
                coupon[k] += long_strike[fix_index - 1] * fixNotional[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fix_pay_dates[fix_index], yc_name);
                swaption_cpn[k] += long_strike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                                   swp_f_df(today, fix_pay_dates[fix_index], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(today, float_pay_dates[float_index], yc_name);
                floatcoupon[k] += -(float_spreads[float_index - 1] + margin[float_index - 1]) *
                                  float_cvgs[float_index - 1] *
                                  swp_f_df(today, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                coupon[k] += long_strike[fix_index - 1] * fixNotional[fix_index - 1] *
                             fix_cvgs[fix_index - 1] *
                             swp_f_df(today, fix_pay_dates[fix_index], yc_name);
                swaption_cpn[k] += long_strike[fix_index - 1] * fix_cvgs[fix_index - 1] *
                                   swp_f_df(today, fix_pay_dates[fix_index], yc_name);
                ex_coupon[fix_index - 1] = k - fix_float_mult;
                exer_date[fix_index - 1] = add_unit(
                    fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] += -(-floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(today, float_pay_dates[float_index], yc_name);
                swaption_cpn[k] +=
                    -(-1 + (float_spreads[float_index - 1] + margin[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_dates[float_index], yc_name);
                floatcoupon[k] +=
                    -(-1.0 + (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
        }
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Swap Schedule-----------------------
    //-----------------------------------------------------------------

    /*	Underlyings */

    /*	Long */

    dff = swp_f_df(today, float_pay_dates[nbcoupon - 1], yc_name);
    for (i = 0; i < nbexercise; i++)
    {
        ex_lstrike[i] = long_strike[i];
        ex_lprice[i]  = diag_prices[i];
    }

    dIRRs = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRs)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    if ((strike_type == 4) || (strike_type == 5) || (strike_type == 6))
    {
        err = ComputeAmortSwapDiagonalIRRs(
            yc_name,
            vol_curve_name,
            default_ref_rate_name,
            start_date,
            end_date,
            ifreq,
            ibasis,
            fix_nb_dates,
            fixNotional,
            long_strike,
            ex_fee,
            float_nb_dates,
            floatNotional,
            margin,
            dIRRs,
            UseVol);
        if (err)
        {
            goto FREE_RETURN;
        }
    }

    dIRRstrikes = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRstrikes)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    if ((strike_type == 4) || (strike_type == 5) || (strike_type == 6))
    {
        err = ComputeSwapRateFromIRR(
            yc_name,
            vol_curve_name,
            today,
            default_ref_rate_name,
            start_date,
            end_date,
            ifreq,
            ibasis,
            dIRRs[0],
            UseVol,
            &dtest);
        if (err)
        {
            goto FREE_RETURN;
        }
    }

    /*	Short */

    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    ifreq,
                    ibasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            else if (strike_type == 6)
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    end_date,
                    ifreq,
                    ibasis,
                    dIRRs[i],
                    UseVol,
                    &swaption_strikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = fix_float_mult * nbexercise - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            if (strike_type == 6)
            {
                lvl = 0.0;
                for (k = i; k < fix_nb_pay_dates - 1; k++)
                {
                    lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k + 1], yc_name);
                }
                dfi = swp_f_df(today, fix_pay_dates[i], yc_name);
                dff = swp_f_df(today, fix_pay_dates[k], yc_name);

                ex_slvl[i] = lvl;

                ex_sfwd[i] = (dfi - dff) / lvl;

                /*	ATM std */
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_pay_dates[k],
                    ex_sfwd[i],
                    0,
                    ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }
            else
            {
                lvl = 0.0;
                for (k = i; k < i + 1; k++)
                {
                    lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k + 1], yc_name);
                }
                dfi = swp_f_df(today, fix_pay_dates[i], yc_name);
                dff = swp_f_df(today, fix_pay_dates[k], yc_name);

                ex_slvl[i] = lvl;
                ex_sfwd[i] = (dfi - dff) / lvl;

                /*	ATM std */
                err = get_cash_vol(
                    vol_curve_name,
                    add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                    fix_pay_dates[k],
                    ex_sfwd[i],
                    0,
                    ref_rate_name,
                    &std,
                    &power);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(exer_time[i]);

            /*	Strike */
            if ((!short_strike) || (!strike_type))
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            else if (strike_type == 1)
            {
                ex_sstrike[i] = short_strike[i];
            }
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        float_pay_dates[ex_coupon[i]],
                        fix_pay_dates[k],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                ex_sstrike[i] = short_strike[i] - spr;
            }
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }
            else if ((strike_type == 4) || (strike_type == 5))
            {
                //				ex_sstrike[i] = ex_sfwd[i];
                ex_sstrike[i] = dIRRstrikes[i];
            }

            /*	Apply max std */
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            /*	Make sure strikes are positive (actually more than 1bp)
                            otherwise use ATM	*/
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            calibstrike = ex_sstrike[i];
            if (strike_type == 6)
            {
                calibstrike = long_strike[i];  // swaption_strikes[i];
            }

            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_pay_dates[k],
                //				ex_sstrike[i],
                calibstrike,
                0,
                ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    //					ex_sstrike[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    //					ex_sstrike[i],
                    calibstrike,
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            if (strike_type == 6)
            {
                cap_price += ex_sprice[i];
            }
            else
            {
                cap_price += ex_sprice[i];
            }
            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }
    }

    /*	The 1F equivalent case */
    if (one2F == 2 && fix_lambda && one_f_equi)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            if ((strike_type == 4) || (strike_type == 5))
            {
                err = ComputeSwapRateFromIRR(
                    yc_name,
                    vol_curve_name,
                    today,
                    default_ref_rate_name,
                    float_start_dates[ex_coupon[i]],
                    float_end_dates[ex_coupon[i]],
                    ifreq,
                    ibasis,
                    dIRRs[i],
                    UseVol,
                    &dIRRstrikes[i]);
                if (err)
                {
                    goto FREE_RETURN;
                }
            }

            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = fix_float_mult * nbexercise - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            lvl = 0.0;
            //			for (k=(int)(ex_coupon[i]*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1));
            //k<(ex_coupon[i]+ex_sncpn[i])*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1); k++)
            //			{
            //				lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k+1],
            //yc_name);
            //			}
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k + 1], yc_name);
            }
            //			dfi = swp_f_df (today, float_pay_dates[ex_coupon[i]], yc_name);
            dfi = swp_f_df(today, fix_pay_dates[i], yc_name);
            dff = swp_f_df(today, fix_pay_dates[k], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            if ((strike_type == 4) || (strike_type == 5))
            {
                //				ex_sstrike[i] = ex_sfwd[i];
                ex_sstrike[i] = dIRRstrikes[i];
                if (dIRRstrikes[i] < 0.0001)
                {
                    ex_sstrike[i] = ex_sfwd[i];
                }
            }
            else
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }

        err = amortMidat_lgmcalibzetalambda(
            nbcoupon,
            coupon,
            coupon_time,
            coupon_df,
            coupon_cvg,
            nbexercise,
            exer_bool,
            exer_time,
            ex_coupon,
            ex_sncpn,
            ex_lprice,
            ex_fee,
            shortcoupon,
            floatcoupon,
            0.0,
            ex_zeta,
            floatNotional,
            swaption_cpn,
            1,
            0,
            lambda,
            1,
            0.0,
            0.0,
            0.0,
            skip_last);

        if (err)
        {
            goto FREE_RETURN;
        }

        export_lgmsetupG(*lambda, nbcoupon, coupon_time, cpn_G, nbexercise, exer_time, ex_G);

        cap_price = amortMidat_lgmcapval1F_b(
            nbcoupon,
            coupon,
            coupon_df,
            coupon_cvg,
            cpn_G,
            nbexercise,
            shortcoupon,
            floatcoupon,
            ex_coupon,
            ex_sncpn,
            ex_zeta,
            ex_G);

        fix_lambda = 0;
    }

    /*	2.)	Calibrate lambda and zeta */

    cap_or_swaption = 0;
    if (strike_type == 6)
    {
        cap_or_swaption = 1;
    }

    err = amortMidat_lgmcalibzetalambda(
        nbcoupon,
        coupon,
        coupon_time,
        coupon_df,
        coupon_cvg,
        nbexercise,
        exer_bool,
        exer_time,
        ex_coupon,
        ex_sncpn,
        ex_lprice,
        ex_fee,
        shortcoupon,
        floatcoupon,
        cap_price,
        ex_zeta,
        floatNotional,
        swaption_cpn,
        fix_lambda,
        cap_or_swaption,
        lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	3.)	Transform into sigma */

    compt = 0;
    for (i = 0; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            compt += 1;
        }
    }

    *num_sig  = compt;
    *sig_time = (double*)calloc(*num_sig, sizeof(double));
    *sig      = (double*)calloc(*num_sig, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    (*sig_time)[0] = exer_time[0];
    (*sig)[0]      = sqrt(ex_zeta[0] * 2 * (*lambda) / (exp(2 * (*lambda) * exer_time[0]) - 1.0));

    last_i = 0;
    k      = 1;
    for (i = 1; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            (*sig_time)[k] = exer_time[i];
            if (ex_zeta[i] > ex_zeta[last_i])
            {
                (*sig)[k] = sqrt(
                    (ex_zeta[i] - ex_zeta[last_i]) * 2 * (*lambda) /
                    (exp(2 * (*lambda) * exer_time[i]) - exp(2 * (*lambda) * exer_time[last_i])));
            }
            else
            {
                smessage(
                    "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                    exer_time[i]);
                for (j = k; j < *num_sig; j++)
                {
                    (*sig)[j] = (*sig)[k - 1];
                }
                i = nbexercise;
            }

            last_i = i;
            k      = k + 1;
        }
    }

    /*	4.)	Save instrument data if required */
    if (inst_data)
    {
        inst_data->num_inst      = nbexercise;
        inst_data->num_insts     = nbexercise;
        inst_data->start_dates   = (long*)calloc(nbexercise, sizeof(long));
        inst_data->end_dates     = (long*)calloc(nbexercise, sizeof(long));
        inst_data->short_strikes = (double*)calloc(nbexercise, sizeof(double));
        inst_data->long_strikes  = (double*)calloc(nbexercise, sizeof(double));

        if (!inst_data->start_dates || !inst_data->end_dates || !inst_data->short_strikes ||
            !inst_data->long_strikes)
        {
            err = "Allocation error (4) in cpd_calib_diagonal";
            goto FREE_RETURN;
        }

        for (i = 0; i < nbexercise; i++)
        {
            inst_data->start_dates[i] = float_pay_dates[ex_coupon[i]];
            inst_data->end_dates[i]   = float_pay_dates[nbcoupon - 1];
            if (!fix_lambda)
            {
                inst_data->short_strikes[i] = ex_sstrike[i];
            }
            else
            {
                inst_data->short_strikes[i] = 0.0;
            }
            inst_data->long_strikes[i] = ex_lstrike[i];
        }
    }

FREE_RETURN:

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_cvgs)
        free(fix_cvgs);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);

    if (coupon_dates)
        free(coupon_dates);

    if (dIRRs)
        free(dIRRs);
    if (dIRRstrikes)
        free(dIRRstrikes);

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;

        /*
        if (inst_data)
        {
                cpd_free_calib_inst_data (inst_data);
        }
        */
    }

    return err;
}

/*	Calibrate lgm: main function */
Err ZCamortMidat_cpd_calib_diagonal(
    int   notperiod,
    char* yc_name,               /*	Name of the yield curve */
    char* vol_curve_name,        /*	Name of the market vol curve */
    char* default_ref_rate_name, /*	Name of the reference rate */
    Err (*get_cash_vol)(         /*	Function to get cash vol from the market */
                        char*   vol_curve_name,
                        double  start_date,
                        double  end_date,
                        double  cash_strike,
                        int     zero,
                        char*   ref_rate_name,
                        double* vol,
                        double* power),
    double vol_shift,
    int    shift_type,     /*	0:	Additive
                                           1:	Multiplicative */
    int* ex_bool,          /*	Exercise or not
                                                   NULL = True */
    long    start_date,    /*	Start date of the amortised swap */
    long    end_date,      /*	End date for the amortised swap */
    double* long_strike,   /*	Diagonal swaption strikes
                                                                   NULL = ATM */
    double* short_strike,  /*	Short swaption strikes
                                                                   NULL = ATM */
    int strike_type,       /*	0: ATM
                                                   1: CASH
                                                   2: SWAP
                                                   3: STD
                                                   4: IRR */
    double* diag_prices,   /* Diagonal Prices */
    double* fixNotional,   /* Amortised Notional
                                           NULL = No amortisation */
    double* floatNotional, /* Amortised Notional
                                           NULL = No amortisation */
    double* margin,

    double max_std_short,
    char*  ref_rate_name,
    char*  swaption_freq, /*	Frequency, basis and ref. rate of calibrated swaptions */
    char*  swaption_basis,
    int    fix_lambda,   /*	0: calib lambda to cap, 1: fix lambda calib
                                                 to diagonal */
    int one_f_equi,      /*	1F equivalent flag:
                                                 if set to 1, then 2F lambda will calibrate
                                                 to the cap priced within calibrated 1F
                                                 with the given lambda */
    int skip_last,       /*	If 1, the last option is disregarded
                                                 and the forward volatility is flat from option
                                                 n-1 */
    double max_var_jump, /*	Maximum multiplicative variation in variance */

    double* lambda, /*	Lambda: may be changed in the process */
    int     one2F,  /*	Number of factors */
    /*	Alpha, Gamma, Rho (2F only) */
    double   alpha,
    double   gamma,
    double   rho,
    int*     num_sig, /*	Answer */
    double** sig_time,
    double** sig,
    /*	Calibration instrument data */
    CPD_CALIB_INST_DATA inst_data) /*	NULL = don't save calibration instrument data */
{
    int            i, j, k, ncpn, n, float_index, fix_index;
    int            nbcoupon;
    SrtCompounding ifreq;
    SrtBasisCode   ibasis;
    double         ex_lstrike[MAX_CPN], ex_lprice[MAX_CPN], ex_sfwd[MAX_CPN], ex_slvl[MAX_CPN],
        ex_sstrike[MAX_CPN], ex_svol[MAX_CPN], ex_sprice[MAX_CPN], ex_zeta[MAX_CPN], ex_G[MAX_CPN];
    int         ex_sncpn[MAX_CPN];
    double      cpn_G[MAX_CPN];
    long        theo_date, act_date;
    long        today;
    double      lvl, dfi, dff;
    double      power;
    double      swp_rte, spr;
    double      cap_price;
    double      std;
    SrtCurvePtr yc_ptr;
    Err         err = NULL;
    double      coupon[MAX_CPN];
    double      floatcoupon[MAX_CPN];
    double      shortcoupon[MAX_CPN];
    double      coupon_time[MAX_CPN];
    double      coupon_df[MAX_CPN];
    double      coupon_cvg[MAX_CPN];
    int         ex_coupon[MAX_CPN];

    double temp;

    int tem;

    int fix_float_mult;

    int last_i;

    double* dIRRs       = NULL;
    double* dIRRstrikes = NULL;

    double  CumCoupon;
    double* lastFixCoupon = NULL;

    //---------------------------------------
    //------------SWAPDP---------------------
    //---------------------------------------
    SwapDP swapdp;
    long   float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    long    fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    double* coupon_dates = NULL;

    double maxNotional;

    int    nbexercise;
    double exer_time[MAX_CPN];
    long   exer_date[MAX_CPN];
    int    exer_bool[MAX_CPN];

    int firstDoCalib;

    int compt;

    //	double dtest;
    //---------------------------------------
    //---------------------------------------

    MAX_FACT = max_var_jump;

    *sig_time = NULL;
    *sig      = NULL;

    yc_ptr = lookup_curve(yc_name);
    if (!yc_ptr)
    {
        err = "Yield Curve not found";
        goto FREE_RETURN;
    }
    today = get_today_from_curve(yc_ptr);

    if (one2F == 2)
    {
        HermiteStandard(x, w, NUM_HERMITE);
    }

    /*	1.)	Setup the bond schedule and its coupons */

    /*	Coupons */

    err = interp_compounding(swaption_freq, &ifreq);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = interp_basis(swaption_basis, &ibasis);
    if (err)
    {
        goto FREE_RETURN;
    }

    theo_date = end_date;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    ncpn      = 1;

    while (act_date > today)
    {
        theo_date = add_unit(theo_date, -12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
        act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
        ncpn++;
    }
    ncpn--;

    if (ncpn < 2)
    {
        err = "Not enough coupons in amortMidat_cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    theo_date = end_date;
    act_date  = bus_date_method(theo_date, MODIFIED_SUCCEEDING);
    i         = ncpn - 1;

    //--------------------------------------------
    //-------------Find Start Date----------------
    //--------------------------------------------
    if (add_unit(start_date, -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING) <= today)
    {
        ex_bool[0] = 0;
    }

    i = 0;
    while ((ex_bool[i] == 0) && (i < ncpn))
    {
        ++i;
    }

    if (i == ncpn)
    {
        err = "At least one maturity must be calibrated : Check DoCalib Value";
        goto FREE_RETURN;
    }

    firstDoCalib = i;

    //--------------------------------------------
    //--------- If Start = today -----------------
    //--------------------------------------------

    err = swp_f_initSwapDP(start_date, theo_date, swaption_freq, swaption_basis, &swapdp);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp,
        today,
        ref_rate_name,
        &float_pay_dates,
        &float_nb_pay_dates,
        &float_fixing_dates,
        &float_start_dates,
        &float_end_dates,
        &float_cvgs,
        &float_spreads,
        &float_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp,
        today,
        &fix_pay_dates,
        &fix_nb_pay_dates,
        &fix_start_dates,
        &fix_end_dates,
        &fix_cvgs,
        &fix_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }

    start_date = fix_start_dates[firstDoCalib];

    tem = firstDoCalib;

    for (i = firstDoCalib; i < fix_nb_dates; ++i)
    {
        ex_bool[i - firstDoCalib]     = ex_bool[i];
        exer_bool[i - firstDoCalib]   = ex_bool[i];
        fixNotional[i - firstDoCalib] = fixNotional[i];
        diag_prices[i - firstDoCalib] = diag_prices[i];
        long_strike[i - firstDoCalib] = long_strike[i];
        if (short_strike)
        {
            short_strike[i - firstDoCalib] = short_strike[i];
        }
    }

    for (i = (int)(0.5 + fix_cvgs[1] / float_cvgs[1]) * firstDoCalib; i < float_nb_dates; ++i)
    {
        floatNotional[i - (int)(0.5 + fix_cvgs[1] / float_cvgs[1]) * firstDoCalib] =
            floatNotional[i];
    }

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_cvgs)
        free(fix_cvgs);

    float_fixing_dates = NULL;
    float_start_dates  = NULL;
    float_end_dates    = NULL;
    float_pay_dates    = NULL;
    float_cvgs         = NULL;
    float_spreads      = NULL;
    float_pay_times    = NULL;

    fix_start_dates = NULL;
    fix_end_dates   = NULL;
    fix_pay_dates   = NULL;
    fix_cvgs        = NULL;
    fix_pay_times   = NULL;

    //--------------------------------------------
    //---------SWAPDP-----------------------------
    //---------Build Swap Schedule----------------
    //--------------------------------------------

    err = swp_f_initSwapDP(start_date, theo_date, swaption_freq, swaption_basis, &swapdp);
    if (err)
    {
        goto FREE_RETURN;
    }

    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp,
        today,
        ref_rate_name,
        &float_pay_dates,
        &float_nb_pay_dates,
        &float_fixing_dates,
        &float_start_dates,
        &float_end_dates,
        &float_cvgs,
        &float_spreads,
        &float_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }
    float_pay_times = dvector(0, float_nb_pay_dates - 1);
    for (n = 0; n < float_nb_pay_dates; ++n)
    {
        float_pay_times[n] = (float_pay_dates[n] - today) * YEARS_IN_DAY;
    }

    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp,
        today,
        &fix_pay_dates,
        &fix_nb_pay_dates,
        &fix_start_dates,
        &fix_end_dates,
        &fix_cvgs,
        &fix_nb_dates);
    if (err)
    {
        goto FREE_RETURN;
    }
    fix_pay_times = dvector(0, fix_nb_pay_dates - 1);
    for (n = 0; n < fix_nb_pay_dates; ++n)
    {
        fix_pay_times[n] = (fix_pay_dates[n] - today) * YEARS_IN_DAY;
    }

    fix_float_mult = (int)(fix_cvgs[fix_nb_dates - 1] / float_cvgs[float_nb_dates - 1] + 0.5);

    //------------------------------------------------
    //----------Rescale Notionals---------------------
    //------------------------------------------------
    maxNotional = fabs(fixNotional[0]);
    for (i = 1; i < fix_nb_dates; ++i)
    {
        if (maxNotional < fabs(fixNotional[i]))
        {
            maxNotional = fabs(fixNotional[i]);
        }
    }
    for (i = 0; i < float_nb_dates; ++i)
    {
        if (maxNotional < fabs(floatNotional[i]))
        {
            maxNotional = fabs(floatNotional[i]);
        }
    }

    CALPRESNOT = maxNotional * CALPRESAMORTMIDAT;

    //-------------------------------------------
    //-------------------------------------------

    nbcoupon = float_nb_pay_dates;

    coupon_dates = (double*)calloc(nbcoupon, sizeof(double));
    memcpy(coupon_dates, float_pay_times, nbcoupon * sizeof(double));
    num_f_concat_vector(&nbcoupon, &coupon_dates, fix_nb_pay_dates, fix_pay_times);
    num_f_sort_vector(nbcoupon, coupon_dates);
    num_f_unique_vector(&nbcoupon, coupon_dates);

    for (k = 0; k < nbcoupon; ++k)
    {
        coupon[k]      = 0;
        floatcoupon[k] = 0;
    }

    nbexercise     = fix_nb_pay_dates - 1;
    float_index    = 1;
    fix_index      = 1;
    coupon[0]      = -floatNotional[0] * swp_f_df(today, float_pay_dates[0], yc_name);
    floatcoupon[0] = -swp_f_df(today, float_pay_dates[0], yc_name);
    coupon_time[0] = coupon_dates[0];
    coupon_df[0]   = swp_f_df(today, float_pay_dates[0], yc_name);
    coupon_cvg[0]  = 0;
    for (k = 1; k < nbcoupon; ++k)
    {
        coupon_time[k] = coupon_dates[k];
        coupon_df[k]   = swp_f_df(today, float_pay_dates[k], yc_name);
        coupon_cvg[k]  = float_cvgs[k - 1];
        if ((coupon_dates[k] == fix_pay_times[fix_index]) &&
            (coupon_dates[k] != float_pay_times[float_index]))
        {
            //			coupon[k] += long_strike[fix_index-1] * fixNotional[fix_index-1] *
            //fix_cvgs[fix_index-1] * swp_f_df (today, fix_pay_dates[fix_index], yc_name);
            ex_coupon[fix_index - 1] = k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
            exer_date[fix_index - 1] =
                add_unit(fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
            exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
            fix_index                = fix_index + 1;
        }
        else if (
            (coupon_dates[k] != fix_pay_times[fix_index]) &&
            (coupon_dates[k] == float_pay_times[float_index]))
        {
            temp = swp_f_df(today, float_pay_dates[float_index], yc_name);
            coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                           floatNotional[float_index - 1] *
                               (float_spreads[float_index - 1] + margin[float_index - 1]) *
                               float_cvgs[float_index - 1]) *
                         swp_f_df(today, float_pay_dates[float_index], yc_name);
            floatcoupon[k] += -(float_spreads[float_index - 1] + margin[float_index - 1]) *
                              float_cvgs[float_index - 1] *
                              swp_f_df(today, float_pay_dates[float_index], yc_name);
            float_index = float_index + 1;
        }
        else
        {
            if (float_index < float_nb_pay_dates - 1)
            {
                //				coupon[k] += long_strike[fix_index-1] *
                //fixNotional[fix_index-1] * fix_cvgs[fix_index-1] * swp_f_df (today,
                //fix_pay_dates[fix_index], yc_name);
                ex_coupon[fix_index - 1] = k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
                exer_date[fix_index - 1] = add_unit(
                    fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] += -(floatNotional[float_index] - floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(today, float_pay_dates[float_index], yc_name);
                floatcoupon[k] += -(float_spreads[float_index - 1] + margin[float_index - 1]) *
                                  float_cvgs[float_index - 1] *
                                  swp_f_df(today, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
            else
            {
                //				coupon[k] += long_strike[fix_index-1] *
                //fixNotional[fix_index-1] * fix_cvgs[fix_index-1] * swp_f_df (today,
                //fix_pay_dates[fix_index], yc_name);
                ex_coupon[fix_index - 1] = k - (float_nb_pay_dates - 1) / (fix_nb_pay_dates - 1);
                exer_date[fix_index - 1] = add_unit(
                    fix_pay_dates[fix_index - 1], -notperiod, SRT_BDAY, MODIFIED_SUCCEEDING);
                exer_time[fix_index - 1] = (exer_date[fix_index - 1] - today) * YEARS_IN_DAY;
                fix_index                = fix_index + 1;
                coupon[k] += -(-floatNotional[float_index - 1] +
                               floatNotional[float_index - 1] *
                                   (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                   float_cvgs[float_index - 1]) *
                             swp_f_df(today, float_pay_dates[float_index], yc_name);
                floatcoupon[k] +=
                    -(-1.0 + (float_spreads[float_index - 1] + margin[float_index - 1]) *
                                 float_cvgs[float_index - 1]) *
                    swp_f_df(today, float_pay_dates[float_index], yc_name);
                float_index = float_index + 1;
            }
        }
    }

    lastFixCoupon = (double*)calloc(fix_nb_pay_dates - 1, sizeof(double));
    if (!lastFixCoupon)
    {
        err = "Allocation failed in zcMidat Calib";
        goto FREE_RETURN;
    }

    temp      = 1;
    CumCoupon = 0.0;
    dff       = swp_f_df(today, fix_end_dates[fix_nb_dates - 1], yc_name);
    for (k = fix_nb_dates - 1; k >= 0; --k)
    {
        CumCoupon += long_strike[k] * fix_cvgs[k] * fixNotional[k] * temp;
        temp             = temp * (1 + long_strike[k] * fix_cvgs[k]);
        lastFixCoupon[k] = CumCoupon * dff;
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Swap Schedule-----------------------
    //-----------------------------------------------------------------

    /*	Underlyings */

    /*	Long */

    dff = swp_f_df(today, float_pay_dates[nbcoupon - 1], yc_name);
    for (i = 0; i < nbexercise; i++)
    {
        ex_lstrike[i] = long_strike[i];
        ex_lprice[i]  = diag_prices[i];
    }

    dIRRs = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRs)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    /*	if ((strike_type == 4)||(strike_type == 5))
            {
                    err = ComputeAmortSwapDiagonalIRRs(yc_name,
                                                    default_ref_rate_name,
                                                    start_date,
                                                    end_date,
                                                    ifreq,
                                                    ibasis,
                                                    fix_nb_dates,
                                                    fixNotional,
                                                    long_strike,
                                                    float_nb_dates,
                                                    floatNotional,
                                                    margin,
                                                    dIRRs);
                    if(err)
                    {
                            goto FREE_RETURN;
                    }
            }
    */
    dIRRstrikes = (double*)calloc(nbexercise, sizeof(double));
    if (!dIRRstrikes)
    {
        err = "memory allocation failed in amort_midat_calib";
        goto FREE_RETURN;
    }

    /*	if ((strike_type == 4)||(strike_type == 5))
            {
                    err = ComputeSwapRateFromIRR(yc_name,
                                                                                    today,
                                                                                    default_ref_rate_name,
                                                                                    start_date,
                                                                                    end_date,
                                                                                    ifreq,
                                                                                    ibasis,
                                                                                    dIRRs[0],
                                                                                    &dtest);
                    if(err)
                    {
                            goto FREE_RETURN;
                    }
            }
    */
    /*	Short */

    if (!fix_lambda)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            /*			if ((strike_type == 4)||(strike_type == 5))
                                    {
                                            err = ComputeSwapRateFromIRR(yc_name,
                                                                                            today,
                                                                                            default_ref_rate_name,
                                                                                            float_start_dates[ex_coupon[i]],
                                                                                            float_end_dates[ex_coupon[i]],
                                                                                            ifreq,
                                                                                            ibasis,
                                                                                            dIRRs[i],
                                                                                            &dIRRstrikes[i]);
                                            if(err)
                                            {
                                                    goto FREE_RETURN;
                                            }
                                    }
            */
            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = fix_float_mult * nbexercise - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            lvl = 0.0;
            //			for (k=(int)(ex_coupon[i]*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1));
            //k<(ex_coupon[i]+ex_sncpn[i])*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1); k++)
            //			{
            //				lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k+1],
            //yc_name);
            //			}
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k + 1], yc_name);
            }
            //			dfi = swp_f_df (today, float_pay_dates[ex_coupon[i]], yc_name);
            dfi = swp_f_df(today, fix_pay_dates[i], yc_name);
            dff = swp_f_df(today, fix_pay_dates[k], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            /*	ATM std */
            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_pay_dates[k],
                ex_sfwd[i],
                0,
                ref_rate_name,
                &std,
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            std += (shift_type == 1 ? std * vol_shift : vol_shift);
            if (power > 0.5)
            {
                power = srt_f_optblksch(
                    ex_sfwd[i], ex_sfwd[i], std, exer_time[i], 1.0, SRT_CALL, PREMIUM);
                err = srt_f_optimpvol(
                    power, ex_sfwd[i], ex_sfwd[i], exer_time[i], 1.0, SRT_CALL, SRT_NORMAL, &std);
            }
            std *= sqrt(exer_time[i]);

            /*	Strike */
            if ((!short_strike) || (!strike_type))
            {
                ex_sstrike[i] = ex_sfwd[i];
            }
            else if (strike_type == 1)
            {
                ex_sstrike[i] = short_strike[i];
            }
            else if (strike_type == 2)
            {
                if (err = swp_f_ForwardRate(
                        float_pay_dates[ex_coupon[i]],
                        fix_pay_dates[k],
                        swaption_freq,
                        swaption_basis,
                        yc_name,
                        ref_rate_name,
                        &swp_rte))
                {
                    goto FREE_RETURN;
                }

                spr = swp_rte - ex_sfwd[i];

                ex_sstrike[i] = short_strike[i] - spr;
            }
            else if (strike_type == 3)
            {
                ex_sstrike[i] = ex_sfwd[i] + short_strike[i] * std;
            }
            else if ((strike_type == 4) || (strike_type == 5))
            {
                ex_sstrike[i] = ex_sfwd[i];
                //				ex_sstrike[i] = dIRRstrikes[i];
            }

            /*	Apply max std */
            if (ex_sstrike[i] > ex_sfwd[i] + max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] + max_std_short * std;
            }
            else if (ex_sstrike[i] < ex_sfwd[i] - max_std_short * std)
            {
                ex_sstrike[i] = ex_sfwd[i] - max_std_short * std;
            }

            /*	Make sure strikes are positive (actually more than 1bp)
                            otherwise use ATM	*/
            if (ex_sstrike[i] < 1.0e-04)
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            err = get_cash_vol(
                vol_curve_name,
                add_unit(exer_date[i], notperiod, SRT_BDAY, MODIFIED_SUCCEEDING),
                fix_pay_dates[k],
                ex_sstrike[i],
                0,
                ref_rate_name,
                &(ex_svol[i]),
                &power);
            if (err)
            {
                goto FREE_RETURN;
            }
            ex_svol[i] += (shift_type == 1 ? ex_svol[i] * vol_shift : vol_shift);

            if (power > 0.5)
            {
                ex_sprice[i] = srt_f_optblksch(
                    ex_sfwd[i],
                    ex_sstrike[i],
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }
            else
            {
                ex_sprice[i] = srt_f_optblknrm(
                    ex_sfwd[i],
                    ex_sstrike[i],
                    ex_svol[i],
                    exer_time[i],
                    ex_slvl[i],
                    SRT_PUT,
                    PREMIUM);
            }

            cap_price += ex_sprice[i];
            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }
    }

    /*	The 1F equivalent case */
    if (one2F == 2 && fix_lambda && one_f_equi)
    {
        cap_price = 0.0;
        for (i = 0; i < nbexercise; i++)
        {
            /*			if ((strike_type == 4)||(strike_type == 5))
                                    {
                                            err = ComputeSwapRateFromIRR(yc_name,
                                                                                            today,
                                                                                            default_ref_rate_name,
                                                                                            float_start_dates[ex_coupon[i]],
                                                                                            float_end_dates[ex_coupon[i]],
                                                                                            ifreq,
                                                                                            ibasis,
                                                                                            dIRRs[i],
                                                                                            &dIRRstrikes[i]);
                                            if(err)
                                            {
                                                    goto FREE_RETURN;
                                            }
                                    }
            */
            if (i < nbexercise - 1)
            {
                ex_sncpn[i] = ex_coupon[i + 1] - ex_coupon[i] + 1;
            }
            else
            {
                ex_sncpn[i] = fix_float_mult * nbexercise - ex_coupon[i] + 1;
            }

            if (ex_sncpn[i] < 2)
            {
                err = "One exercise date controls less than 2 coupons in cpd_calib_diagonal";
                goto FREE_RETURN;
            }

            lvl = 0.0;
            //			for (k=(int)(ex_coupon[i]*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1));
            //k<(ex_coupon[i]+ex_sncpn[i])*(fix_nb_pay_dates-1)/(float_nb_pay_dates-1); k++)
            //			{
            //				lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k+1],
            //yc_name);
            //			}
            for (k = i; k < i + 1; k++)
            {
                lvl += fix_cvgs[k] * swp_f_df(today, fix_pay_dates[k + 1], yc_name);
            }
            //			dfi = swp_f_df (today, float_pay_dates[ex_coupon[i]], yc_name);
            dfi = swp_f_df(today, fix_pay_dates[i], yc_name);
            dff = swp_f_df(today, fix_pay_dates[k], yc_name);

            ex_slvl[i] = lvl;
            ex_sfwd[i] = (dfi - dff) / lvl;

            if ((strike_type == 4) || (strike_type == 5))
            {
                //				ex_sstrike[i] = ex_sfwd[i];
                ex_sstrike[i] = dIRRstrikes[i];
            }
            else
            {
                ex_sstrike[i] = ex_sfwd[i];
            }

            shortcoupon[i] = ex_sstrike[i] * fix_cvgs[i];
        }

        err = zcamortMidat_lgmcalibzetalambda(
            nbcoupon,
            coupon,
            lastFixCoupon,
            coupon_time,
            coupon_df,
            coupon_cvg,
            nbexercise,
            exer_bool,
            exer_time,
            ex_coupon,
            ex_sncpn,
            ex_lprice,
            shortcoupon,
            floatcoupon,
            0.0,
            ex_zeta,
            floatNotional,
            1,
            lambda,
            1,
            0.0,
            0.0,
            0.0,
            skip_last);

        if (err)
        {
            goto FREE_RETURN;
        }

        export_lgmsetupG(*lambda, nbcoupon, coupon_time, cpn_G, nbexercise, exer_time, ex_G);

        cap_price = amortMidat_lgmcapval1F_b(
            nbcoupon,
            coupon,
            coupon_df,
            coupon_cvg,
            cpn_G,
            nbexercise,
            shortcoupon,
            floatcoupon,
            ex_coupon,
            ex_sncpn,
            ex_zeta,
            ex_G);

        fix_lambda = 0;
    }

    /*	2.)	Calibrate lambda and zeta */

    err = zcamortMidat_lgmcalibzetalambda(
        nbcoupon,
        coupon,
        lastFixCoupon,
        coupon_time,
        coupon_df,
        coupon_cvg,
        nbexercise,
        exer_bool,
        exer_time,
        ex_coupon,
        ex_sncpn,
        ex_lprice,
        shortcoupon,
        floatcoupon,
        cap_price,
        ex_zeta,
        floatNotional,
        fix_lambda,
        lambda,
        one2F,
        alpha,
        gamma,
        rho,
        skip_last);

    if (err)
    {
        goto FREE_RETURN;
    }

    /*	3.)	Transform into sigma */

    compt = 0;
    for (i = 0; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            compt += 1;
        }
    }

    *num_sig  = compt;
    *sig_time = (double*)calloc(*num_sig, sizeof(double));
    *sig      = (double*)calloc(*num_sig, sizeof(double));

    if (!sig_time || !sig)
    {
        err = "Allocation error (3) in cpd_calib_diagonal";
        goto FREE_RETURN;
    }

    (*sig_time)[0] = exer_time[0];
    (*sig)[0]      = sqrt(ex_zeta[0] * 2 * (*lambda) / (exp(2 * (*lambda) * exer_time[0]) - 1.0));

    last_i = 0;
    k      = 1;
    for (i = 1; i < nbexercise; i++)
    {
        if (exer_bool[i] > 0)
        {
            (*sig_time)[k] = exer_time[i];
            if (ex_zeta[i] > ex_zeta[last_i])
            {
                (*sig)[k] = sqrt(
                    (ex_zeta[i] - ex_zeta[last_i]) * 2 * (*lambda) /
                    (exp(2 * (*lambda) * exer_time[i]) - exp(2 * (*lambda) * exer_time[last_i])));
            }
            else
            {
                smessage(
                    "Diagonal calibration failed at exercise year %.2f - Calibration stopped",
                    exer_time[i]);
                for (j = k; j < *num_sig; j++)
                {
                    (*sig)[j] = (*sig)[k - 1];
                }
                i = nbexercise;
            }

            last_i = i;
            k      = k + 1;
        }
    }

    /*	4.)	Save instrument data if required */
    if (inst_data)
    {
        inst_data->num_inst      = nbexercise;
        inst_data->num_insts     = nbexercise;
        inst_data->start_dates   = (long*)calloc(nbexercise, sizeof(long));
        inst_data->end_dates     = (long*)calloc(nbexercise, sizeof(long));
        inst_data->start_datess  = (long*)calloc(nbexercise, sizeof(long));
        inst_data->end_datess    = (long*)calloc(nbexercise, sizeof(long));
        inst_data->short_strikes = (double*)calloc(nbexercise, sizeof(double));
        inst_data->long_strikes  = (double*)calloc(nbexercise, sizeof(double));

        if (!inst_data->start_dates || !inst_data->end_dates || !inst_data->short_strikes ||
            !inst_data->long_strikes)
        {
            err = "Allocation error (4) in cpd_calib_diagonal";
            goto FREE_RETURN;
        }

        for (i = 0; i < nbexercise; i++)
        {
            inst_data->start_dates[i]  = float_pay_dates[ex_coupon[i]];
            inst_data->end_dates[i]    = float_pay_dates[nbcoupon - 1];
            inst_data->start_datess[i] = float_pay_dates[ex_coupon[i]];
            inst_data->end_datess[i]   = float_pay_dates[nbcoupon - 1];

            if (!fix_lambda)
            {
                inst_data->short_strikes[i] = ex_sstrike[i];
            }
            else
            {
                inst_data->short_strikes[i] = 0.0;
            }
            inst_data->long_strikes[i] = ex_lstrike[i];
        }
    }

FREE_RETURN:

    if (lastFixCoupon)
        free(lastFixCoupon);

    if (float_fixing_dates)
        free(float_fixing_dates);
    if (float_start_dates)
        free(float_start_dates);
    if (float_end_dates)
        free(float_end_dates);
    if (float_pay_dates)
        free(float_pay_dates);
    if (float_cvgs)
        free(float_cvgs);
    if (float_spreads)
        free(float_spreads);
    if (float_pay_times)
        free_dvector(float_pay_times, 0, float_nb_pay_dates - 1);

    if (fix_start_dates)
        free(fix_start_dates);
    if (fix_end_dates)
        free(fix_end_dates);
    if (fix_pay_dates)
        free(fix_pay_dates);
    if (fix_cvgs)
        free(fix_cvgs);
    if (fix_pay_times)
        free_dvector(fix_pay_times, 0, fix_nb_pay_dates - 1);

    if (coupon_dates)
        free(coupon_dates);

    if (dIRRs)
        free(dIRRs);
    if (dIRRstrikes)
        free(dIRRstrikes);

    if (err)
    {
        if (*sig_time)
            free(*sig_time);
        *sig_time = NULL;

        if (*sig)
            free(*sig);
        *sig = NULL;

        /*
        if (inst_data)
        {
                cpd_free_calib_inst_data (inst_data);
        }
        */
    }

    return err;
}
