/******************************************************************
  srt_f_bndcrv.c
  Author: E.Auld
  Math: R.Roy

  This module comprises the SrtBndCrv object, which is an interest rate
  curve estimated from a set of forward bond prices.  Once initialized,
  a SrtBndCrv can give discount factors and zero rates; forward yields and
  forward bond prices

  The model:
  f(0,t) = b_0 + b_1 * t + b_2 * t^2
    sum{(i=1 to n) b_i+2 * (t-k_i)^2 * Indicator(t > k_i)}
  Here the k_i are knot points. b_i's are chosen to minimize
  weighted sum of squared errors using Levenberg Marquart.

  b_i+2 will be meaningless if k_i >= longest maturity bond;
  in this case we only optimize over b_0..b_i+1.
  We force the following condition: f'(0,k_i_last_optimised+1) = 0;
  this is solved by forcing the value of b_i_last_optimised+2; hence we do not
  optimize over b_i_last_optimised+1.

  A note on dates/coverages:
  This code prices all bonds as if the coupon to be payed is accrued for the
  real amount of time between the coupon dates.

*******************************************************************/
#include "math.h"
#include "num_h_allhdr.h"
#include "stdarg.h"
#include "swp_h_all.h"
#include "utallhdr.h"

/******************************************************************
  internal definitions
*******************************************************************/

#define DERIV_INTERNAL (1.0e-8)

/******************************************************************
  static global variables
*******************************************************************/

static SrtBndCrv*      global_crv;
static int             global_num_bond;
static SrtBasisCode*   global_basis;
static SrtCompounding* global_compounding;
static Ddate*          global_prev_cpn_date;
static int*            global_num_coupons;
static Ddate**         global_coupon_dates;
static double**        global_coupons;
static double*         global_weight;
static double*         global_fwd_price;
static double*         global_redemption;
static int             global_algo_flag;

/*******************************************************************
 Initialisation of knot points
*******************************************************************/

/* compare to double */
static int compare_double(const void* elem1, const void* elem2)
{
    double* delem1 = (double*)elem1;
    double* delem2 = (double*)elem2;

    if (fabs(*delem1 - *delem2) < EPS)
        return 0;
    else if (*delem1 > *delem2)
        return 1;
    else
        return -1;
}

/* The knot point array will contain all date from the bond array */
static Err populate_knot_point(
    SrtBndCrv* crv, int* numknot, int* num_coupons, Ddate** coupon_dates, int numbond, Ddate today)
{
    long   modified_end;
    double times;
    int    i, j;
    int    found; /* boolean */

    /* First fill the knot point array */
    for (i = 0; i < numbond; i++)
    {
        found = 0;

        modified_end =
            bus_date_method((Date)coupon_dates[i][num_coupons[i] - 1], MODIFIED_SUCCEEDING);
        /* times in year of the bond */
        /* It needs two be at least 2 days */
        times = (modified_end - today) * YEARS_IN_DAY;
        if (times < 0)
            continue;

        /* look if the times is already in the knot point array */
        for (j = 1; j <= *numknot; j++)
        {
            if (fabs((crv->knot_point)[j] - times) < YEARS_IN_DAY / 10)
            {
                found = 1;
                break;
            }
        }

        /* the times is not in the knot array */
        /* it is added to the knot point array and the knot date array */
        if (!found)
        {
            (*numknot)++;

            (crv->knot_point) =
                (double*)realloc((crv->knot_point), (*numknot + 1) * sizeof(double));
            (crv->knot_date) = (double*)realloc((crv->knot_date), (*numknot + 1) * sizeof(double));
            if ((!crv->knot_point) || (!crv->knot_date))
            {
                srt_free(crv->knot_point);
                srt_free(crv->knot_date);
                return serror("Memory allocation error");
            }

            (crv->knot_point)[*numknot] = times;
            (crv->knot_date)[*numknot]  = modified_end;
        }
    }

    /* Order it by date */
    qsort((void*)(&crv->knot_point[1]), *numknot, sizeof(double), compare_double);
    qsort((void*)(&crv->knot_date[1]), *numknot, sizeof(double), compare_double);

    return NULL;
}

/*******************************************************************
 Constraints on length
*******************************************************************/

/* Smoothness of the curve */
static double srt_f_bndcrvlen(SrtBndCrv* crv)
{
    int    i;
    double tmp;
    double crv_len = 0;
    double lastint = (crv->lastdate - crv->today) * YEARS_IN_DAY;

    /* we integrate a specific function f''(0,t)^2 between the nodes */
    tmp = 6 * crv->zcpn_yld_coef[0][3];
    tmp *= tmp;
    tmp *= (crv->knot_point[1] - 0);

    crv_len += tmp;

    /* loop */
    for (i = 1; i < crv->lastknot - 1; i++)
    {
        tmp = 6 * crv->zcpn_yld_coef[i][3];
        tmp *= tmp;
        tmp *= (crv->knot_point[i + 1] - crv->knot_point[i]);

        crv_len += tmp;
    }

    /* last one */
    if (crv->knot_point[crv->lastknot] < lastint)
    {
        tmp = 6 * crv->zcpn_yld_coef[crv->lastknot][3];
        tmp *= tmp;
        tmp *= (lastint - crv->knot_point[crv->lastknot]);

        crv_len += tmp;
    }

    return crv_len;
}

/******************************************************************
  Calculation functions
*******************************************************************/

/*
 * calculate forward value of a bond using a SrtBndCrv
 * forward value is to cpn_date[0]
 * accrued interest between prev_cpn_date and cpn_date[0] is
 * subtracted
 */
static double fwdclnprc(
    SrtBndCrv*     crv,
    Ddate          prev_cpn_date,
    Ddate*         cpn_date,
    int            num_cpn_date,
    double*        cpn,
    double         redemption,
    SrtBasisCode   bc,
    SrtCompounding compd)
{
    double dr, df, cov, fwd_price = 0.0, acc_int;
    int    i;

    /*
     * first coupon
     */
    srt_f_bndcrv_df((void*)crv, (Ddate)cpn_date[0], (Ddate)cpn_date[1], &df);
    if (bc == BASIS_ACT_USD)
        cov = 1.0 / compd;
    else
        cov = coverage((Date)prev_cpn_date, (Date)cpn_date[1], bc);
    fwd_price = cpn[1] * cov * df;

    /*
     * other coupons
     */
    for (i = 2; i < num_cpn_date; i++)
    {
        srt_f_bndcrv_df((void*)crv, (Ddate)cpn_date[0], (Ddate)cpn_date[i], &df);
        if (bc == BASIS_ACT_USD)
            cov = 1.0 / compd;
        else
            cov = coverage((Date)cpn_date[i - 1], (Date)cpn_date[i], bc);
        fwd_price += cpn[i] * cov * df;
    }

    /*
     * notional
     */
    fwd_price += redemption * df;

    /*
     * subtract accrued interest
     */

    if (bc == BASIS_ACT_USD)
    {
        dr = day_count_date((Date)prev_cpn_date, (Date)cpn_date[0], bc);
        dr = dr / day_count_date((Date)prev_cpn_date, (Date)cpn_date[1], bc);
        dr = dr / compd;
    }
    else
    {
        dr = coverage((Date)prev_cpn_date, (Date)cpn_date[0], bc);
    }

    acc_int = cpn[1] * dr;
    fwd_price -= acc_int;

    return (fwd_price);
}

/******************************************************************
  Levenberg-Marquart functions
*******************************************************************/

/*
 * reinitialize the crv object after b0..b_n have been changed.
 */
static void reset_zcpn_yld_coefs(SrtBndCrv* crv)
{
    double tmp = 0.0;
    int    i, numpts = crv->num_knot_point;

    /**********/
    /* SPLINE */
    /**********/
    if (global_algo_flag == 2)
    {
        /*
         * choose the next to last coefficient (b_lastknot+2) to force f'(0,k_lastknot+1) = 0
         */

        for (i = 1; i <= crv->lastknot - 1; i++)
        {
            tmp += crv->fwd_rate_coef[i + 2] *
                   (crv->knot_point[crv->lastknot + 1] - crv->knot_point[i]);
        }

        tmp = (-crv->fwd_rate_coef[1] -
               2.0 * crv->fwd_rate_coef[2] * crv->knot_point[crv->lastknot + 1] - 2.0 * tmp) /
              (2.0 * (crv->knot_point[crv->lastknot + 1] - crv->knot_point[crv->lastknot]));

        crv->fwd_rate_coef[crv->lastknot + 2] = tmp;

        /*
         * between knot points, -log(df(0,t)) is a cubic; coefs are here
         * we say crv->zcpn_yld_coef[i][0..3] is the third order polynomial
         * (index 3 is the highest order coefficient) between knot point i
         * and knot point i+1.  Here we consider the zeroth knot point to always
         * have value zero.
         */

        crv->zcpn_yld_coef[0][3] = crv->fwd_rate_coef[2] / 3.0;
        crv->zcpn_yld_coef[0][2] = 0.5 * crv->fwd_rate_coef[1];
        crv->zcpn_yld_coef[0][1] = crv->fwd_rate_coef[0];
        crv->zcpn_yld_coef[0][0] = 0.0;

        for (i = 1; i <= numpts; i++)
        {
            tmp                      = crv->fwd_rate_coef[i + 2];
            crv->zcpn_yld_coef[i][3] = crv->zcpn_yld_coef[i - 1][3] + tmp / 3.0;

            tmp *= crv->knot_point[i];
            crv->zcpn_yld_coef[i][2] = crv->zcpn_yld_coef[i - 1][2] - tmp;

            tmp *= crv->knot_point[i];
            crv->zcpn_yld_coef[i][1] = crv->zcpn_yld_coef[i - 1][1] + tmp;

            tmp *= crv->knot_point[i];
            crv->zcpn_yld_coef[i][0] = crv->zcpn_yld_coef[i - 1][0] - tmp / 3.0;
        }
    }
    /**********/
    /* LINEAR */
    /**********/
    else if (global_algo_flag == 1)
    {
        /* boundary conditions */
        /*		for(i=-1 ; i <= crv->lastknot-1 ; i++)
                        {
                                tmp -= crv->fwd_rate_coef[i+2];
                        }

                        crv->fwd_rate_coef[crv->lastknot+2] = tmp;
        */
        /*
         * between knot points, -log(df(0,t)) is a second order polynomial; coefs are here
         * we say crv->zcpn_yld_coef[i][0..3] is the third order polynomial
         * (index 3 is the highest order coefficient) between knot point i
         * and knot point i+1.  Here we consider the zeroth knot point to always
         * have value zero.
         */

        crv->zcpn_yld_coef[0][3] = 0.0;
        /*		crv->zcpn_yld_coef[0][2] = 0.5*crv->fwd_rate_coef[1];
         */
        crv->zcpn_yld_coef[0][2] = 0.0;
        crv->zcpn_yld_coef[0][1] = crv->fwd_rate_coef[0];
        crv->zcpn_yld_coef[0][0] = 0.0;

        for (i = 1; i <= numpts; i++)
        {
            crv->zcpn_yld_coef[i][3] = 0.0;

            tmp                      = crv->fwd_rate_coef[i + 2];
            crv->zcpn_yld_coef[i][2] = crv->zcpn_yld_coef[i - 1][2] + tmp / 2.0;

            tmp *= crv->knot_point[i];
            crv->zcpn_yld_coef[i][1] = crv->zcpn_yld_coef[i - 1][1] - tmp;

            tmp *= crv->knot_point[i];
            crv->zcpn_yld_coef[i][0] = crv->zcpn_yld_coef[i - 1][0] + tmp / 2.0;
        }
    }
    else
    {
        /*
         * between knot points, -log(df(0,t)) is a second order polynomial; coefs are here
         * we say crv->zcpn_yld_coef[i][0..3] is the third order polynomial
         * (index 3 is the highest order coefficient) between knot point i
         * and knot point i+1.  Here we consider the zeroth knot point to always
         * have value zero.
         */

        crv->zcpn_yld_coef[0][3] = 0.0;
        /*		crv->zcpn_yld_coef[0][2] = 0.5*crv->fwd_rate_coef[1];
         */
        crv->zcpn_yld_coef[0][2] = 0.0;
        crv->zcpn_yld_coef[0][1] = crv->fwd_rate_coef[0];
        crv->zcpn_yld_coef[0][0] = 0.0;

        for (i = 1; i <= numpts; i++)
        {
            crv->zcpn_yld_coef[i][3] = 0.0;
            crv->zcpn_yld_coef[i][2] = 0.0;

            tmp                      = crv->fwd_rate_coef[i + 2];
            crv->zcpn_yld_coef[i][1] = crv->zcpn_yld_coef[i - 1][1] + tmp;

            crv->zcpn_yld_coef[i][0] = 0.0;
        }
    }
}

/*
 * calculate integral(0,d){f(0,s)ds}
 */
static double calc_cum_f_to_dt(SrtBndCrv* crv, Ddate d)
{
    double t = (d - crv->today) * YEARS_IN_DAY;
    double cum_f, *coef;
    int    i = 0;

    if (t < EPS)
        return 0.0;

    /*
     * find index of last knot point <= t
     */
    while (i < crv->num_knot_point && crv->knot_point[i + 1] < t)
        i++;

    coef  = crv->zcpn_yld_coef[i];
    cum_f = coef[0] + t * (coef[1] + t * (coef[2] + t * coef[3]));

    return cum_f;
}

/* Calculate the forward price for the jth bond */
/* or the length */
static void value_jth_bond(int j, double* answer)
{
    *answer = fwdclnprc(
        global_crv,
        global_prev_cpn_date[j],
        global_coupon_dates[j],
        global_num_coupons[j],
        global_coupons[j],
        global_redemption[j],
        global_basis[j],
        global_compounding[j]);
}

/*
 * this is the function that gets passed to the optimization function
 * It calculates f(x) and the derivatives of f at x.  In this particular
 * case x is the index of the jth bond.  a[] is the things we are optimizing
 * over; dyda is the derivative of f(x) at a w.r.t. all those things.
 */
static Err funcs(double x, double a[], double* yfit, double dyda[], int ma)
{
    double up, mid;
    int    i, j;

    j = DTOL(x);

    /* we have to be careful for the last one */
    /* The last coefficient is for the length function */
    memcpy(global_crv->fwd_rate_coef, &a[1], (ma - 1) * sizeof(double));
    reset_zcpn_yld_coefs(global_crv);

    /* For a bond */
    if (j < global_num_bond)
    {
        value_jth_bond(j, &mid);

        /* this will be a wasted step when it reaches b_n-1 */
        /* there is two in the difference of the indices of a[] and fwd_rate_coef[] */
        for (i = 0; i < ma - 1; i++)
        {
            global_crv->fwd_rate_coef[i] += DERIV_INTERNAL;
            reset_zcpn_yld_coefs(global_crv);
            value_jth_bond(j, &up);
            global_crv->fwd_rate_coef[i] -= DERIV_INTERNAL;
            dyda[i + 1] = (up - mid) / DERIV_INTERNAL;
        }
    }
    /* For the length */
    else
    {
        mid = srt_f_bndcrvlen(global_crv);

        /* this will be a wasted step when it reaches b_n-1 */
        /* there is two in the difference of the indices of a[] and fwd_rate_coef[] */
        for (i = 0; i < ma - 1; i++)
        {
            global_crv->fwd_rate_coef[i] += DERIV_INTERNAL;
            reset_zcpn_yld_coefs(global_crv);
            up = srt_f_bndcrvlen(global_crv);
            global_crv->fwd_rate_coef[i] -= DERIV_INTERNAL;
            dyda[i + 1] = (up - mid) / DERIV_INTERNAL;
        }
    }

    /* It is for the smoothness so the derivation is null */
    dyda[ma] = 0.0;

    *yfit = mid;

    return NULL;
}

/******************************************************************
  External Functions
*******************************************************************/

/*
 * function to initialize a bond curve by choosing a quadratic spline
 * on forward rates that minimizes a weighted sum of squared errors.
 * uses a nonlinear optimization routine (Levenberg Marquardt) from
 * Numerical Recipes.
 */

SrtErr srt_f_bndcrv_init(
    SrtBndCrv**     crvptr,        /* Return of the bond curve */
    SrtBasisCode*   basis,         /* array containing basis of bonds */
    SrtCompounding* compounding,   /* array containing compounding of bonds */
    Ddate*          prev_cpn_date, /* array containing previous coupon date of bonds */
    int*            num_coupons,   /* array containing number of coupon in bonds */
    Ddate**         coupon_dates,  /* array containing coupon dates of bonds*/
    double**        coupons,       /* array containing coupons of bonds*/
    double*         weight,        /* array containing weight (standart deviation) of bonds */
    double*         fwd_price,     /* array containing the forward clean price */
    double*         redemption,    /* array containing the redemption */
    int             numbond,       /* number of bonds */
    double*         knot_point,    /* array containing the knot points (in times) */
    int             numknot,       /* number of knot points */
    Ddate           today,         /* today's date */
    double          coeff_line,    /* coefficient for smoothness */
    int             num_iter,      /* Number of iterations for levenberg */
    int             algo_flag)                 /* 0 = CONSTANT, 1 = LINEAR, 2 = SPLINE */
{
    SrtBndCrv* crv;
    int        i;
    long *     ia, PARA_NUMB, PRICE_NUMB;
    double     chisq;
    double *   x, *y, *sig, *p;
    Date       lastdate, knot_date_last;
    double     knot_point_last;
    Err        err;

    /*
     * allocate memory for the curve and initialise the different arrays in the curve
     */
    crv = (SrtBndCrv*)srt_calloc(1, sizeof(SrtBndCrv));
    if (!crv)
    {
        return serror("alloc failure in srt_f_bndcrv_init");
    }

    global_crv           = crv;
    global_num_bond      = numbond;
    global_basis         = basis;
    global_compounding   = compounding;
    global_prev_cpn_date = prev_cpn_date;
    global_num_coupons   = num_coupons;
    global_coupon_dates  = coupon_dates;
    global_coupons       = coupons;
    global_fwd_price     = fwd_price;
    global_redemption    = redemption;
    global_weight        = weight;
    global_algo_flag     = algo_flag;

    lastdate = (Date)0;
    for (i = 0; i < numbond; i++)
    {
        lastdate = IMAX((Date)coupon_dates[i][num_coupons[i] - 1], lastdate);
    }

    /* Check the knots */
    if (knot_point == NULL)
    {
        /* if there is no knot points, all bond dates become knots points */
        if (populate_knot_point(crv, &numknot, num_coupons, coupon_dates, numbond, today) != NULL)
        {
            srt_free(crv);
            return serror("alloc failure in srt_f_bndcrv_init");
        }
    }
    else
    {
        /* we have to add one point if the lastdate is superior to the last knot date */
        knot_date_last = (Date)add_unit(
            DTOL(crv->spot), DTOL(12 * knot_point[numknot - 1]), SRT_MONTH, MODIFIED_SUCCEEDING);
        knot_point_last = (double)(knot_date_last - today) * YEARS_IN_DAY;
        if (lastdate > knot_point_last)
        {
            numknot++;
        }

        crv->knot_date  = dvector(1, numknot);
        crv->knot_point = dvector(1, numknot);
        if ((!crv->knot_date) || (!crv->knot_point))
        {
            free_dvector(crv->knot_date, 1, numknot);
            free_dvector(crv->knot_point, 1, numknot);
            srt_free(crv);
            return serror("alloc failure in srt_f_bndcrv_init");
        }

        for (i = 1; i <= numknot - 1; i++)
        {
            crv->knot_date[i] = (Ddate)add_unit(
                DTOL(crv->spot), DTOL(12 * knot_point[i - 1]), SRT_MONTH, MODIFIED_SUCCEEDING);
            crv->knot_point[i] = (double)(crv->knot_date[i] - today) * YEARS_IN_DAY;
        }

        /* A point should be add */
        if (lastdate > knot_point_last)
        {
            knot_date_last          = lastdate + (lastdate - knot_date_last);
            knot_point_last         = (double)(knot_date_last - today) * YEARS_IN_DAY;
            crv->knot_date[numknot] = (Ddate)add_unit(
                DTOL(crv->spot), DTOL(12 * knot_point_last), SRT_MONTH, MODIFIED_SUCCEEDING);
        }
        else
        {
            crv->knot_date[numknot] = (Ddate)add_unit(
                DTOL(crv->spot), DTOL(12 * knot_point[numknot]), SRT_MONTH, MODIFIED_SUCCEEDING);
        }

        crv->knot_point[numknot] = (double)(crv->knot_date[numknot] - today) * YEARS_IN_DAY;

        /* Order it by date */
        qsort((void*)(&crv->knot_point[1]), numknot, sizeof(double), compare_double);
        qsort((void*)(&crv->knot_date[1]), numknot, sizeof(double), compare_double);
    }

    /* Allocate memory */
    crv->fwd_rate_coef = dvector(0, numknot + 2);
    crv->zcpn_yld_coef = dmatrix(0, numknot, 0, 3);
    if ((!crv->fwd_rate_coef) || (!crv->zcpn_yld_coef))
    {
        free_dmatrix(crv->zcpn_yld_coef, 0, numknot, 0, 3);
        free_dvector(crv->fwd_rate_coef, 0, numknot + 2);
        free_dvector(crv->knot_date, 1, numknot);
        free_dvector(crv->knot_point, 1, numknot);
        srt_free(crv);
        return serror("alloc failure in srt_f_bndcrv_init");
    }

    /* initialise the different fields in the curve structure */
    crv->num_knot_point = numknot;
    crv->lastdate       = lastdate;
    crv->today          = today;
    crv->algo_flag      = algo_flag;

    PARA_NUMB  = numknot + 3 + 1; /* The last one for the length */
    PRICE_NUMB = numbond + 1;     /* The last one for the length */

    p   = dvector(1, PARA_NUMB);   /* things to optimize over */
    ia  = lngvector(1, PARA_NUMB); /* flags to optimize or not */
    x   = dvector(1, PRICE_NUMB);  /* indices of bonds */
    y   = dvector(1, PRICE_NUMB);  /* prices of bonds */
    sig = dvector(1, PRICE_NUMB);  /* weights of bonds */

    if ((!p) || (!ia) || (!x) || (!y) || (!sig))
    {
        free_dvector(p, 1, PARA_NUMB);
        free_lngvector(ia, 1, PARA_NUMB);
        free_dvector(sig, 1, PRICE_NUMB);
        free_dvector(x, 1, PRICE_NUMB);
        free_dvector(y, 1, PRICE_NUMB);

        free_dmatrix(crv->zcpn_yld_coef, 0, numknot, 0, 3);
        free_dvector(crv->fwd_rate_coef, 0, numknot + 2);
        free_dvector(crv->knot_date, 1, numknot);
        free_dvector(crv->knot_point, 1, numknot);
        srt_free(crv);
        return serror("alloc failure in srt_f_bndcrv_init");
    }

    /*
     * put knot points in bond curve, initialize coefficients we are optimizing
     * over, (these are the coefficients of the quadratic forward rate spline)
     * and calculate the corresponding coefficients of the cubic equation this
     * implies for the zero coupon yields.
     */

    /* modify the coeff */
    /* we don't want to have 0 or 1 */
    if (coeff_line <= EPS)
        coeff_line = EPS;
    else if (coeff_line >= 1.0 - EPS)
        coeff_line = 1.0 - EPS;
    else
    {
        coeff_line = log(pow(coeff_line, 1.0 / 10.0) * (exp(1) - 1) + 1);

        if (coeff_line <= EPS)
            coeff_line = EPS;
        else if (coeff_line >= 1.0 - EPS)
            coeff_line = 1.0 - EPS;
    }

    /*
     * optimize by calling mrqmin (Numerical Recipes levenberg marquart routine);
     * we want to minimize the (weighted) sum of squared errors for our bonds.
     */

    for (i = 1; i <= PARA_NUMB; i++)
        ia[i] = 1; /* SHOW THAT WE FIT ALL PARAMS */

    /* we optimise the coefficient in p */
    memset(&p[1], 0, PARA_NUMB * sizeof(double));

    /*
     * the last knot points are only meaningful if the last bond is after
     * those knot points;
     */

    i = numknot;
    while ((i > 0) && (lastdate <= crv->knot_date[i]))
    {
        ia[i + 3] = 0;
        i--;
    }

    /* initialise the last knot optimised */
    crv->lastknot = i;

    /************************************/
    /* CONSTRAINTS DUE TO THE ALGORITHM */
    /************************************/

    /* The length doesn't need to be optimised if the algo is linear */
    /* and the convexity of the first coeff is zero */
    if ((algo_flag == 1) || (algo_flag == 0)) /* LINEAR or CONSTANT */
    {
        /* coeff of smoothness */
        coeff_line = EPS;

        /* Length */
        p[PARA_NUMB] = 0.0;

        /* Convexity */
        ia[3] = 0;
        p[3]  = 0;

        /* linear  change */
        ia[2] = 0;
        p[2]  = 0;
    }
    else if (algo_flag == 2) /* SPLINE */
    {
        /* Length */
        p[PARA_NUMB] = 1.0; /* see below */
    }

    /* INITIALISE THE DIFFERENT COEFFICIENT FOR THE LEVENBERG MARQUARDT ALGORITHM */
    /* For the length */
    /* To be the smallest as possible */
    ia[PARA_NUMB] = 0;
    x[PRICE_NUMB] = PRICE_NUMB - 1;
    y[PRICE_NUMB] = 0;
    /* renormalise */
    sig[PRICE_NUMB] = ((double)(lastdate - today) * YEARS_IN_DAY) / coeff_line;

    /* For the bond prices */
    for (i = 0; i < PRICE_NUMB - 1; i++)
    {
        sig[i + 1] = weight[i] / 10.0 / (1.0 - coeff_line);
        x[i + 1]   = (double)i;
        y[i + 1]   = fwd_price[i];
    }

    /*******************************************/
    /* OK NOW WE USE THE DIFFERENT CONSTRAINTS */
    /*******************************************/

    /* LAST CONSTRAINT */
    /* f'(k_lastknot+1) = 0, fix ia[lastknot+3] */
    if (algo_flag == 2) /* SPLINE */
    {
        ia[crv->lastknot + 3] = 0;
    }

    /************************/
    /* AT LAST LE ALGORITHM */
    /************************/
    if (err = levenberg_marquardt_select(
            x, y, sig, PRICE_NUMB, p, ia, PARA_NUMB, num_iter, funcs, &chisq))
    {
        free_dmatrix(crv->zcpn_yld_coef, 0, numknot, 0, 3);
        free_dvector(crv->fwd_rate_coef, 0, numknot + 2);
        free_dvector(crv->knot_date, 1, numknot);
        free_dvector(crv->knot_point, 1, numknot);
        srt_free(crv);
        return err;
    }

    /*
    free date schedules,
    other stuff
    */
    free_dvector(p, 1, PARA_NUMB);
    free_lngvector(ia, 1, PARA_NUMB);
    free_dvector(sig, 1, PRICE_NUMB);
    free_dvector(x, 1, PRICE_NUMB);
    free_dvector(y, 1, PRICE_NUMB);

    /* Return */
    if (err)
        return err;

    *crvptr = crv;

    return NULL;
}

/*
 * instantaneous forward rate at date start
 */
SrtErr srt_f_bndcrv_fwdr(SrtBndCrv* crv, Ddate d, double* answer)
{
    double t = (d - crv->today) * YEARS_IN_DAY;
    double fr;
    int    i;

    if (d < crv->today)
        return serror("date must be >= today");

    fr = crv->fwd_rate_coef[0] + t * (crv->fwd_rate_coef[1] + t * crv->fwd_rate_coef[2]);

    /*
     * find index of last knot point >= t
     */
    for (i = 1; i <= crv->num_knot_point && crv->knot_point[i] < t; i++)
    {
        fr += pow((t - crv->knot_point[i]), crv->algo_flag) * crv->fwd_rate_coef[i + 2];
    }

    *answer = fr;

    return NULL;
}

/*
  level from SrtBndCrv
*/
SrtErr srt_f_bndcrv_level(SrtBndCrv* crv, SwapDP* sdp, double* answer)
{
    SrtErr   err = NULL;
    int      i;
    DateList dl;
    double   level = 0;
    double   cov, df;

    dl = SwapDP_to_DateList(sdp, MODIFIED_SUCCEEDING);

    for (i = 0; i < dl.len - 1; i++)
    {
        if ((err = srt_f_bndcrv_df((void*)crv, dl.date[i], dl.date[i + 1], &df)) != NULL)
            return err;

        if (sdp->basis_code == BASIS_ACT_USD)
            cov = 1.0 / sdp->compd;
        else
            cov = coverage(dl.date[i], dl.date[i + 1], sdp->basis_code);

        level += df * cov;
    }

    srt_free(dl.date);

    *answer = level;

    return NULL;
}

/*
  discount factor from SrtBndCrv
*/
SrtErr srt_f_bndcrv_df(void* crv, Ddate start, Ddate end, double* answer)
{
    double     cum_f_to_s, cum_f_to_e, df;
    SrtBndCrv* bond_crv = (SrtBndCrv*)crv;

    if (fabs(start - end) < EPS)
    {
        *answer = 1.0;
        return NULL;
    }

    cum_f_to_s = calc_cum_f_to_dt(bond_crv, start);
    cum_f_to_e = calc_cum_f_to_dt(bond_crv, end);

    df      = exp(-(cum_f_to_e - cum_f_to_s));
    *answer = df;

    return NULL;
}

/*
  zero rate from SrtBndCrv
  Impossible to do it for US (US = ACT/365)
*/
SrtErr srt_f_bndcrv_zr(SrtBndCrv* crv, Ddate start, Ddate end, SrtBasisCode basis, double* answer)
{
    double cum_f_to_s, cum_f_to_e, zr;

    if (fabs(start - end) < EPS)
    {
        *answer = 0.0;
        return NULL;
    }

    cum_f_to_s = calc_cum_f_to_dt(crv, start);
    cum_f_to_e = calc_cum_f_to_dt(crv, end);

    zr      = (cum_f_to_e - cum_f_to_s) / coverage((Date)start, (Date)end, basis);
    *answer = zr;

    return NULL;
}

/*
  fwd par bond yield
*/
SrtErr srt_f_bndcrv_fwdtr(SrtBndCrv* crv, Date s, Date e_nfp, char* f, char* b, double* answer)
{
    SrtRtFnc rt;
    SrtErr   err;
    int      i, len;
    Ddate*   d;
    double*  df;

    err = srt_f_rtmk("SWAP", s, NULL, e_nfp, NULL, f, b, &rt);
    if (err)
        return err;

    df  = srt_f_rtgetdf(rt);
    d   = srt_f_rtgetdt(rt);
    len = srt_f_rtgetlen(rt);

    for (i = 0; i < len; i++)
    {
        srt_f_bndcrv_df((void*)crv, crv->today, d[i], &df[i]);
    }

    srt_f_rtevl(rt, crv->fwd_rate_coef[0], answer);
    srt_f_rtfre(rt);

    return NULL;
}

/*
   fwd bond price
*/
SrtErr srt_f_bndcrv_fwdpr(
    SrtBndCrv* crv, SwapDP* sdp, double cpn, double redemption, double* answer)
{
    int      i;
    DateList dl;
    double*  cpn_array;
    Ddate*   date_array;

    dl = SwapDP_to_DateList(sdp, MODIFIED_SUCCEEDING);

    cpn_array  = (double*)calloc(dl.len, sizeof(double));
    date_array = (Ddate*)calloc(dl.len, sizeof(Ddate));
    if ((!cpn_array) || (!date_array))
        return serror("Memory Allocation Failure");

    for (i = 0; i < dl.len; i++)
    {
        date_array[i] = (Ddate)dl.date[i];
        cpn_array[i]  = cpn;
    }

    if (dl.type == NOTBROKEN)
        dl.prev = dl.date[0];

    *answer = fwdclnprc(
        crv, dl.prev, date_array, dl.len, cpn_array, redemption, sdp->basis_code, sdp->compd);

    srt_free(dl.date);
    srt_free(cpn_array);
    srt_free(date_array);

    return NULL;
}

/*
   fwd bond yield
*/
SrtErr srt_f_bndcrv_fwdyld(
    SrtBndCrv* crv, SwapDP* sdp, double cpn, double redemption, double* answer)
{
    SrtErr err;

    err = srt_f_bndcrv_fwdpr(crv, sdp, cpn, redemption, answer);
    if (err)
        return err;

    (*answer) = yield_fct(*sdp, cpn, (*answer), 1.0, cpn);

    return NULL;
}

/*******************************************************************
 Utility Functions
*******************************************************************/

/*
  find out lastdate
*/
Date srt_f_bndcrv_lastdate(SrtBndCrv* crv)
{
    return (Date)DTOL(crv->lastdate);
}

/*
  find out today
*/
Date srt_f_bndcrv_today(SrtBndCrv* crv)
{
    return (Date)DTOL(crv->today);
}

/*
  free a SrtBndCrv
*/
void srt_f_bndcrv_free(SrtBndCrv** crv)
{
    if ((*crv)->knot_date)
    {
        free_dvector((*crv)->knot_date, 1, (*crv)->num_knot_point);
        (*crv)->knot_date = NULL;
    }

    if ((*crv)->knot_point)
    {
        free_dvector((*crv)->knot_point, 1, (*crv)->num_knot_point);
        (*crv)->knot_point = NULL;
    }

    if ((*crv)->fwd_rate_coef)
    {
        free_dvector((*crv)->fwd_rate_coef, 0, (*crv)->num_knot_point + 2);
        (*crv)->fwd_rate_coef = NULL;
    }

    if ((*crv)->zcpn_yld_coef)
    {
        free_dmatrix((*crv)->zcpn_yld_coef, 0, (*crv)->num_knot_point, 0, 5);
        (*crv)->zcpn_yld_coef = NULL;
    }

    if (*crv)
        srt_free((*crv));
}
