/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:
        Author		:

        Created		:

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "LGM_CMS.h"

#include "AmortMidatCalib.h"
#include "math.h"
#include "num_h_gausslegendre.h"
#include "num_h_proba.h"
#include "srt_h_all.h"
#include "srt_h_allFx3F.h"
#include "swp_h_cms.h"

const double Pi = 3.14159265358979323846;

#define XL_MAX_ID_LEN 40

// A modifier
//	const int XL_MAX_ID_LEN=1234567890;

// compute l'int in sigma(i) from 0 to t
double compute_sum_sigma(double t, double* sigma, double* sigma_date, int nb_sigma, double lambda)
{
    double t1;
    double t2;
    double sigma_temp;
    double res = 0;
    int    i   = 0;
    t2         = 0;
    while ((sigma_date[i] < t) & (i < nb_sigma))
    {
        if (i > 0)
        {
            t2 = sigma_date[i];
            t1 = sigma_date[i - 1];
        }
        else  // i=0
        {
            t2 = sigma_date[i];
            t1 = 0;
        }
        res += sigma[i] * sigma[i] * (exp(2 * lambda * t2) - exp(2 * lambda * t1)) / (2 * lambda);
        i++;
    }
    t1 = t2;
    t2 = t;
    if (i == nb_sigma)
    {
        sigma_temp = sigma[i - 1];
    }
    else
    {
        sigma_temp = sigma[i];
    }
    res += sigma_temp * sigma_temp * (exp(2 * lambda * t2) - exp(2 * lambda * t1)) / (2 * lambda);
    return res;
}

double compute_sum_diagterm(
    double  t,
    double* sigma_1,
    double* sigma_2,
    double* sigma_date,
    int     nb_sigma,
    double  lambda_1,
    double  lambda_2)
{
    double t1;
    double t2;
    double sigma1_temp;
    double sigma2_temp;
    double res = 0;
    int    i   = 0;
    t2         = 0;
    while ((sigma_date[i] < t) & (i < nb_sigma))
    {
        if (i > 0)
        {
            t2 = sigma_date[i];
            t1 = sigma_date[i - 1];
        }
        else  // i=0
        {
            t2 = sigma_date[i];
            t1 = 0;
        }
        res += sigma_1[i] * sigma_2[i] *
               (exp((lambda_1 + lambda_2) * t2) - exp((lambda_1 + lambda_2) * t1)) /
               (lambda_1 + lambda_2);
        i++;
    }
    t1 = t2;
    t2 = t;
    if (i == nb_sigma)
    {
        sigma1_temp = sigma_1[i - 1];
        sigma2_temp = sigma_2[i - 1];
    }
    else
    {
        sigma1_temp = sigma_1[i];
        sigma2_temp = sigma_2[i];
    }
    res += sigma1_temp * sigma2_temp *
           (exp((lambda_1 + lambda_2) * t2) - exp((lambda_1 + lambda_2) * t1)) /
           (lambda_1 + lambda_2);
    return res;
}

double compute_A(
    double t,
    double _T,
    double sum_sigma_1,
    double sum_sigma_2,
    double sum_diagterm,
    double rho,
    double lambda_1,
    double lambda_2,
    double Tpay)
{
    double Ares     = 0;
    double int_low1 = (1 - exp(lambda_1 * (Tpay - t))) / lambda_1;
    double int_low2 = (1 - exp(lambda_2 * (Tpay - t))) / lambda_2;
    double int_up1  = (1 - exp(lambda_1 * (Tpay - _T))) / lambda_1;
    double int_up2  = (1 - exp(lambda_2 * (Tpay - _T))) / lambda_2;
    Ares += -(double)(1) / 2 * sum_sigma_1 * exp(-2 * lambda_1 * Tpay) *
            (int_up1 * int_up1 - int_low1 * int_low1);
    Ares += -(double)(1) / 2 * sum_sigma_2 * exp(-2 * lambda_2 * Tpay) *
            (int_up2 * int_up2 - int_low2 * int_low2);
    Ares += -sum_diagterm * rho * exp(-(lambda_1 + lambda_2) * Tpay) *
            (int_up1 * int_up2 - int_low1 * int_low2);
    return Ares;
}

double compute_c1(
    double  t,
    double  _T,
    double  Tpay,
    double* sigma_1,
    double* sigma_2,
    double* sigma_date,
    int     nb_sigma,
    double  rho,
    double  lambda_1,
    double  lambda_2)
{
    double t1;
    double t2;
    int    i     = 0;
    double temp1 = 0;
    double temp2 = 0;
    double temp3 = 0;
    double sigma1_temp;
    double sigma2_temp;
    double res  = 0;
    double coef = 0;
    if ((_T == t) | (t == 0))
    {
        return 0;
    }
    else
    {
        t2 = 0;
        while ((sigma_date[i] < t) & (i < nb_sigma))
        {
            if (i > 0)
            {
                t2 = sigma_date[i];
                t1 = sigma_date[i - 1];
            }
            else  // i=0
            {
                t2 = sigma_date[i];
                t1 = 0;
            }
            temp1 += sigma_1[i] * sigma_1[i] *
                     (exp(-2 * lambda_1 * (Tpay - t2)) - exp(-2 * lambda_1 * (Tpay - t1))) /
                     (2 * lambda_1);
            temp2 += sigma_2[i] * sigma_2[i] *
                     (exp(-2 * lambda_2 * (Tpay - t2)) - exp(-2 * lambda_2 * (Tpay - t1))) /
                     (2 * lambda_2);
            temp3 += sigma_1[i] * sigma_2[i] *
                     (exp(-(lambda_1 + lambda_2) * (Tpay - t2)) -
                      exp(-(lambda_1 + lambda_2) * (Tpay - t1))) /
                     (lambda_1 + lambda_2);
            i++;
        }
        t1 = t2;
        t2 = t;
        if (i == nb_sigma)
        {
            sigma1_temp = sigma_1[i - 1];
            sigma2_temp = sigma_2[i - 1];
        }
        else
        {
            sigma1_temp = sigma_1[i];
            sigma2_temp = sigma_2[i];
        }
        temp1 += sigma1_temp * sigma1_temp *
                 (exp(-2 * lambda_1 * (Tpay - t2)) - exp(-2 * lambda_1 * (Tpay - t1))) /
                 (2 * lambda_1);
        temp2 += sigma2_temp * sigma2_temp *
                 (exp(-2 * lambda_2 * (Tpay - t2)) - exp(-2 * lambda_2 * (Tpay - t1))) /
                 (2 * lambda_2);
        temp3 += sigma1_temp * sigma2_temp *
                 (exp(-(lambda_1 + lambda_2) * (Tpay - t2)) -
                  exp(-(lambda_1 + lambda_2) * (Tpay - t1))) /
                 (lambda_1 + lambda_2);

        coef = (exp(lambda_2 * (Tpay - _T)) - exp(lambda_2 * (Tpay - t))) / lambda_2;
        coef *= lambda_1 / (exp(lambda_1 * (Tpay - _T)) - exp(lambda_1 * (Tpay - t)));
        temp2 *= coef * coef * rho * rho;
        temp3 *= 2 * coef * rho;
        res = temp1 + temp2 + temp3;
        res = sqrt(res);
        res *= (exp(lambda_1 * (Tpay - _T)) - exp(lambda_1 * (Tpay - t))) / lambda_1;
        return res;
    }
}

double compute_c2(
    double  t,
    double  _T,
    double  Tpay,
    double* sigma_2,
    double* sigma_date,
    int     nb_sigma,
    double  rho,
    double  lambda_2)
{
    double t1;
    double t2;
    int    i    = 0;
    double temp = 0;
    double sigma2_temp;
    double res = 0;
    if ((_T == t) | (t == 0))
    {
        return 0;
    }
    else
    {
        t2 = 0;
        while ((sigma_date[i] < t) & (i < nb_sigma))
        {
            if (i > 0)
            {
                t2 = sigma_date[i];
                t1 = sigma_date[i - 1];
            }
            else  // i=0
            {
                t2 = sigma_date[i];
                t1 = 0;
            }
            temp += sigma_2[i] * sigma_2[i] *
                    (exp(-2 * lambda_2 * (Tpay - t2)) - exp(-2 * lambda_2 * (Tpay - t1))) /
                    (2 * lambda_2);
            i++;
        }
        t1 = t2;
        t2 = t;
        if (i == nb_sigma)
        {
            sigma2_temp = sigma_2[i - 1];
        }
        else
        {
            sigma2_temp = sigma_2[i];
        }
        temp = temp + sigma2_temp * sigma2_temp *
                          (exp(-2 * lambda_2 * (Tpay - t2)) - exp(-2 * lambda_2 * (Tpay - t1))) /
                          (2 * lambda_2);
        temp = (1 - rho * rho) * temp;
        res  = sqrt(temp);
        res *= (exp(lambda_2 * (Tpay - _T)) - exp(lambda_2 * (Tpay - t))) / lambda_2;
        return res;
    }
}

double compute_CMS_GaussLeg(
    double          Tf,
    double          Tpay,
    double*         DF_fix,
    double*         fix_coupon_date,
    double*         fix_cvg,  // cvg[0]=0, the same size as fix_coupon_date
    int             nb_fix_coupon,
    double*         DF_float,
    double          DF_TF,
    double*         float_coupon_date,
    double*         float_cvg,  // the same size as float_coupon_date
    int             nb_float_coupon,
    double*         spread,
    double*         sigma_1,
    double*         sigma_2,
    double*         sigma_date,
    int             nb_sigma,
    double          lambda_1,
    double          lambda_2,
    double          rho,
    int             nb_integration,
    double          bound,  // must be positive
    SrtReceiverType ProductType,
    double          strike)
{
    double  y1, y1r, y1m, dy1, y2, y2r, y2m, dy2, s1, s2;
    double* x            = NULL;
    double* w            = NULL;
    double* Afloat       = NULL;
    double* Afix         = NULL;
    double* c1float      = NULL;
    double* c2float      = NULL;
    double* c1fix        = NULL;
    double* c2fix        = NULL;
    double* ZC_float     = NULL;
    double* ZC_fix       = NULL;
    double* coupon_float = NULL;
    double  res          = 0;
    int     i            = 0;
    int     j            = 0;
    int     k            = 0;
    int     l            = 0;
    double  sum_sigma_1  = 0;
    double  sum_sigma_2  = 0;
    double  sum_diagterm = 0;
    double  swap;
    double  swap_num;
    double  swap_den;

    // compute the different A(Tf,coupon_date) and Ci
    Afloat       = (double*)calloc(nb_float_coupon, sizeof(double));
    Afix         = (double*)calloc(nb_fix_coupon, sizeof(double));
    c1float      = (double*)calloc(nb_float_coupon, sizeof(double));
    c2float      = (double*)calloc(nb_float_coupon, sizeof(double));
    c1fix        = (double*)calloc(nb_fix_coupon, sizeof(double));
    c2fix        = (double*)calloc(nb_fix_coupon, sizeof(double));
    ZC_float     = (double*)calloc(nb_float_coupon, sizeof(double));
    ZC_fix       = (double*)calloc(nb_fix_coupon, sizeof(double));
    coupon_float = (double*)calloc(nb_float_coupon, sizeof(double));

    // compute abscissas and weights for Gauss Legendre integration
    x = (double*)calloc(nb_integration, sizeof(double));
    w = (double*)calloc(nb_integration, sizeof(double));
    GaussLeg(-1, 1, x - 1, w - 1, nb_integration);

    // compute B(t,T)=B(0,T)/B(0,t)*exp(A(t,T)+c1(t,T)*z1+c2(t,T)*z2) under QTp
    sum_sigma_1 = compute_sum_sigma(Tf, sigma_1, sigma_date, nb_sigma, lambda_1);
    sum_sigma_2 = compute_sum_sigma(Tf, sigma_2, sigma_date, nb_sigma, lambda_2);
    sum_diagterm =
        compute_sum_diagterm(Tf, sigma_1, sigma_2, sigma_date, nb_sigma, lambda_1, lambda_2);

    for (i = 0; i < nb_float_coupon; i++)
    {
        Afloat[i] = compute_A(
            Tf,
            float_coupon_date[i],
            sum_sigma_1,
            sum_sigma_2,
            sum_diagterm,
            rho,
            lambda_1,
            lambda_2,
            Tpay);
        c1float[i] = compute_c1(
            Tf,
            float_coupon_date[i],
            Tpay,
            sigma_1,
            sigma_2,
            sigma_date,
            nb_sigma,
            rho,
            lambda_1,
            lambda_2);
        c2float[i] = compute_c2(
            Tf, float_coupon_date[i], Tpay, sigma_2, sigma_date, nb_sigma, rho, lambda_2);
    }
    for (i = 0; i < nb_fix_coupon; i++)
    {
        Afix[i] = compute_A(
            Tf,
            fix_coupon_date[i],
            sum_sigma_1,
            sum_sigma_2,
            sum_diagterm,
            rho,
            lambda_1,
            lambda_2,
            Tpay);
        c1fix[i] = compute_c1(
            Tf,
            fix_coupon_date[i],
            Tpay,
            sigma_1,
            sigma_2,
            sigma_date,
            nb_sigma,
            rho,
            lambda_1,
            lambda_2);
        c2fix[i] =
            compute_c2(Tf, fix_coupon_date[i], Tpay, sigma_2, sigma_date, nb_sigma, rho, lambda_2);
    }

    // the same bound in y1 and y2
    y1m = 0.5 * (-bound + bound);
    y1r = 0.5 * (bound + bound);
    y2m = y1m;
    y2r = y1r;
    y1  = 0;
    y2  = 0;
    s1  = 0;
    for (l = 0; l < nb_integration; l++)
    {
        s2  = 0;
        dy1 = y1r * x[l];
        y1  = y1m + dy1;
        for (k = 0; k < nb_integration; k++)
        {
            dy2 = y2r * x[k];
            y2  = y2m + dy2;

            // compute the swap for z1=y1 and z2=y2
            swap        = 0;
            swap_num    = 0;
            swap_den    = 0;
            ZC_float[0] = DF_float[0] / DF_TF * exp(Afloat[0] + c1float[0] * y1 + c2float[0] * y2);
            for (i = 0; i < nb_float_coupon - 1; i++)
            {
                ZC_float[i + 1] = DF_float[i + 1] / DF_TF *
                                  exp(Afloat[i + 1] + c1float[i + 1] * y1 + c2float[i + 1] * y2);
                coupon_float[i] =
                    ((double)(1) / float_cvg[i + 1] * (ZC_float[i] / ZC_float[i + 1] - 1) +
                     spread[i]) *
                    float_cvg[i + 1];
                swap_num += coupon_float[i] * ZC_float[i + 1];
            }
            j = 0;
            i = 0;
            // to avoid computing twice the zero coupons
            while ((i < nb_float_coupon - 1) & (j < nb_fix_coupon))
            {
                if (float_coupon_date[i] > fix_coupon_date[j])
                {
                    ZC_fix[j] = DF_fix[j] / DF_TF * exp(Afix[j] + c1fix[j] * y1 + c2fix[j] * y2);
                    j++;
                }
                else if (
                    (float_coupon_date[i] < fix_coupon_date[j]) &
                    (fix_coupon_date[j] < float_coupon_date[i + 1]))
                {
                    ZC_fix[j] = DF_fix[j] / DF_TF * exp(Afix[j] + c1fix[j] * y1 + c2fix[j] * y2);
                    j++;
                }
                else if (float_coupon_date[i] == fix_coupon_date[j])
                {
                    ZC_fix[j] = ZC_float[i];
                    j++;
                }
                else if (fix_coupon_date[j] == float_coupon_date[i + 1])
                {
                    ZC_fix[j] = ZC_float[i + 1];
                    j++;
                }
                else
                {
                    i++;
                }
            }
            while (j < nb_fix_coupon)
            {
                ZC_fix[j] = DF_fix[j] / DF_TF * exp(Afix[j] + c1fix[j] * y1 + c2fix[j] * y2);
                j++;
            }
            for (j = 0; j < nb_fix_coupon; j++)
            {
                swap_den += fix_cvg[j] * ZC_fix[j];
            }
            // select the good payoff: CMS or option on CMS
            if (ProductType == SRT_RECEIVER)
            {
                swap = strike - swap_num / swap_den;
                if (swap < 0)
                {
                    swap = 0;
                }
            }
            else if (ProductType == SRT_PAYER)
            {
                swap = swap_num / swap_den - strike;
                if (swap < 0)
                {
                    swap = 0;
                }
            }
            else
            {
                swap = swap_num / swap_den;
            }
            s2 += w[k] * swap * exp(-y1 * y1 / 2 - y2 * y2 / 2) / (2 * Pi);
        }
        s1 += w[l] * s2;
    }
    return s1 * y1r * y2r;

    if (Afloat)
        free(Afloat);
    if (Afix)
        free(Afix);
    if (c1float)
        free(c1float);
    if (c2float)
        free(c2float);
    if (c1fix)
        free(c1fix);
    if (c2fix)
        free(c2fix);
    if (ZC_float)
        free(ZC_float);
    if (ZC_fix)
        free(ZC_fix);
    if (coupon_float)
        free(coupon_float);
    if (x)
        free(x);
    if (w)
        free(w);
}

double compute_CMScov_GL(
    double  Tf1,
    double  Tpay1,
    double* DF_fix1,
    double* fix_coupon_date1,
    double* fix_cvg1,  // cvg[0]=0, the same size as fix_coupon_date
    int     nb_fix_coupon1,
    double* DF_float1,
    double  DF_TF1,
    double* float_coupon_date1,
    double* float_cvg1,  // the same size as float_coupon_date
    int     nb_float_coupon1,
    double* spread1,
    double  Tf2,
    double  Tpay2,
    double* DF_fix2,
    double* fix_coupon_date2,
    double* fix_cvg2,  // cvg[0]=0, the same size as fix_coupon_date
    int     nb_fix_coupon2,
    double* DF_float2,
    double  DF_TF2,
    double* float_coupon_date2,
    double* float_cvg2,  // the same size as float_coupon_date
    int     nb_float_coupon2,
    double* spread2,
    double* sigma_11,
    double* sigma_21,
    double* sigma_date1,
    int     nb_sigma1,
    double  lambda_11,
    double  lambda_21,
    double  rho1,
    double* sigma_12,
    double* sigma_22,
    double* sigma_date2,
    int     nb_sigma2,
    double  lambda_12,
    double  lambda_22,
    double  rho2,
    int     nb_integration,
    double  bound  // must be positive
)
{
    double  y1, y1r, y1m, dy1, y2, y2r, y2m, dy2, s1, s2;
    double* x             = NULL;
    double* w             = NULL;
    double* Afloat1       = NULL;
    double* Afix1         = NULL;
    double* c1float1      = NULL;
    double* c2float1      = NULL;
    double* c1fix1        = NULL;
    double* c2fix1        = NULL;
    double* ZC_float1     = NULL;
    double* ZC_fix1       = NULL;
    double* coupon_float1 = NULL;
    double* Afloat2       = NULL;
    double* Afix2         = NULL;
    double* c1float2      = NULL;
    double* c2float2      = NULL;
    double* c1fix2        = NULL;
    double* c2fix2        = NULL;
    double* ZC_float2     = NULL;
    double* ZC_fix2       = NULL;
    double* coupon_float2 = NULL;
    double  res           = 0;
    int     i             = 0;
    int     j             = 0;
    int     k             = 0;
    int     l             = 0;
    double  sum_sigma_11  = 0;
    double  sum_sigma_21  = 0;
    double  sum_diagterm1 = 0;
    double  swap1;
    double  swap_num1;
    double  swap_den1;
    double  sum_sigma_12  = 0;
    double  sum_sigma_22  = 0;
    double  sum_diagterm2 = 0;
    double  swap2;
    double  swap_num2;
    double  swap_den2;
    double  cms1, cms1temp;
    double  cms2, cms2temp;
    double  cov;

    // compute the different A(Tf,coupon_date) and Ci
    Afloat1       = (double*)calloc(nb_float_coupon1, sizeof(double));
    Afix1         = (double*)calloc(nb_fix_coupon1, sizeof(double));
    c1float1      = (double*)calloc(nb_float_coupon1, sizeof(double));
    c2float1      = (double*)calloc(nb_float_coupon1, sizeof(double));
    c1fix1        = (double*)calloc(nb_fix_coupon1, sizeof(double));
    c2fix1        = (double*)calloc(nb_fix_coupon1, sizeof(double));
    ZC_float1     = (double*)calloc(nb_float_coupon1, sizeof(double));
    ZC_fix1       = (double*)calloc(nb_fix_coupon1, sizeof(double));
    coupon_float1 = (double*)calloc(nb_float_coupon1, sizeof(double));

    Afloat2       = (double*)calloc(nb_float_coupon2, sizeof(double));
    Afix2         = (double*)calloc(nb_fix_coupon2, sizeof(double));
    c1float2      = (double*)calloc(nb_float_coupon2, sizeof(double));
    c2float2      = (double*)calloc(nb_float_coupon2, sizeof(double));
    c1fix2        = (double*)calloc(nb_fix_coupon2, sizeof(double));
    c2fix2        = (double*)calloc(nb_fix_coupon2, sizeof(double));
    ZC_float2     = (double*)calloc(nb_float_coupon2, sizeof(double));
    ZC_fix2       = (double*)calloc(nb_fix_coupon2, sizeof(double));
    coupon_float2 = (double*)calloc(nb_float_coupon2, sizeof(double));

    // compute abscissas and weights for Gauss Legendre integration
    x = (double*)calloc(nb_integration, sizeof(double));
    w = (double*)calloc(nb_integration, sizeof(double));
    GaussLeg(-1, 1, x - 1, w - 1, nb_integration);

    // compute B(t,T)=B(0,T)/B(0,t)*exp(A(t,T)+c1(t,T)*z1+c2(t,T)*z2) under QTp
    sum_sigma_11 = compute_sum_sigma(Tf1, sigma_11, sigma_date1, nb_sigma1, lambda_11);
    sum_sigma_21 = compute_sum_sigma(Tf1, sigma_21, sigma_date1, nb_sigma1, lambda_21);
    sum_diagterm1 =
        compute_sum_diagterm(Tf1, sigma_11, sigma_21, sigma_date1, nb_sigma1, lambda_11, lambda_21);

    sum_sigma_12 = compute_sum_sigma(Tf2, sigma_12, sigma_date2, nb_sigma2, lambda_12);
    sum_sigma_22 = compute_sum_sigma(Tf2, sigma_22, sigma_date2, nb_sigma2, lambda_22);
    sum_diagterm2 =
        compute_sum_diagterm(Tf2, sigma_12, sigma_22, sigma_date2, nb_sigma2, lambda_12, lambda_22);

    for (i = 0; i < nb_float_coupon1; i++)
    {
        Afloat1[i] = compute_A(
            Tf1,
            float_coupon_date1[i],
            sum_sigma_11,
            sum_sigma_21,
            sum_diagterm1,
            rho1,
            lambda_11,
            lambda_21,
            Tpay1);
        c1float1[i] = compute_c1(
            Tf1,
            float_coupon_date1[i],
            Tpay1,
            sigma_11,
            sigma_21,
            sigma_date1,
            nb_sigma1,
            rho1,
            lambda_11,
            lambda_21);
        c2float1[i] = compute_c2(
            Tf1, float_coupon_date1[i], Tpay1, sigma_21, sigma_date1, nb_sigma1, rho1, lambda_21);
    }
    for (i = 0; i < nb_fix_coupon1; i++)
    {
        Afix1[i] = compute_A(
            Tf1,
            fix_coupon_date1[i],
            sum_sigma_11,
            sum_sigma_21,
            sum_diagterm1,
            rho1,
            lambda_11,
            lambda_21,
            Tpay1);
        c1fix1[i] = compute_c1(
            Tf1,
            fix_coupon_date1[i],
            Tpay1,
            sigma_11,
            sigma_21,
            sigma_date1,
            nb_sigma1,
            rho1,
            lambda_11,
            lambda_21);
        c2fix1[i] = compute_c2(
            Tf1, fix_coupon_date1[i], Tpay1, sigma_21, sigma_date1, nb_sigma1, rho1, lambda_21);
    }

    for (i = 0; i < nb_float_coupon2; i++)
    {
        Afloat2[i] = compute_A(
            Tf2,
            float_coupon_date2[i],
            sum_sigma_12,
            sum_sigma_22,
            sum_diagterm2,
            rho2,
            lambda_12,
            lambda_22,
            Tpay2);
        c1float2[i] = compute_c1(
            Tf2,
            float_coupon_date2[i],
            Tpay2,
            sigma_12,
            sigma_22,
            sigma_date2,
            nb_sigma2,
            rho2,
            lambda_12,
            lambda_22);
        c2float2[i] = compute_c2(
            Tf2, float_coupon_date2[i], Tpay2, sigma_22, sigma_date2, nb_sigma2, rho2, lambda_22);
    }
    for (i = 0; i < nb_fix_coupon2; i++)
    {
        Afix2[i] = compute_A(
            Tf2,
            fix_coupon_date2[i],
            sum_sigma_12,
            sum_sigma_22,
            sum_diagterm2,
            rho2,
            lambda_12,
            lambda_22,
            Tpay2);
        c1fix2[i] = compute_c1(
            Tf2,
            fix_coupon_date2[i],
            Tpay2,
            sigma_12,
            sigma_22,
            sigma_date2,
            nb_sigma2,
            rho2,
            lambda_12,
            lambda_22);
        c2fix2[i] = compute_c2(
            Tf2, fix_coupon_date2[i], Tpay2, sigma_22, sigma_date2, nb_sigma2, rho2, lambda_22);
    }

    // the same bound in y1 and y2
    y1m  = 0.5 * (-bound + bound);
    y1r  = 0.5 * (bound + bound);
    y2m  = y1m;
    y2r  = y1r;
    y1   = 0;
    y2   = 0;
    s1   = 0;
    cms1 = 0;
    cms2 = 0;
    for (l = 0; l < nb_integration; l++)
    {
        s2       = 0;
        cms1temp = 0;
        cms2temp = 0;
        dy1      = y1r * x[l];
        y1       = y1m + dy1;
        for (k = 0; k < nb_integration; k++)
        {
            dy2 = y2r * x[k];
            y2  = y2m + dy2;

            // compute the swap for z1=y1 and z2=y2
            swap1     = 0;
            swap_num1 = 0;
            swap_den1 = 0;
            swap2     = 0;
            swap_num2 = 0;
            swap_den2 = 0;
            ZC_float1[0] =
                DF_float1[0] / DF_TF1 * exp(Afloat1[0] + c1float1[0] * y1 + c2float1[0] * y2);
            ZC_float2[0] =
                DF_float2[0] / DF_TF2 * exp(Afloat2[0] + c1float2[0] * y1 + c2float2[0] * y2);

            for (i = 0; i < nb_float_coupon1 - 1; i++)
            {
                ZC_float1[i + 1] =
                    DF_float1[i + 1] / DF_TF1 *
                    exp(Afloat1[i + 1] + c1float1[i + 1] * y1 + c2float1[i + 1] * y2);
                coupon_float1[i] =
                    ((double)(1) / float_cvg1[i + 1] * (ZC_float1[i] / ZC_float1[i + 1] - 1) +
                     spread1[i]) *
                    float_cvg1[i + 1];
                swap_num1 += coupon_float1[i] * ZC_float1[i + 1];
            }
            for (i = 0; i < nb_float_coupon2 - 1; i++)
            {
                ZC_float2[i + 1] =
                    DF_float2[i + 1] / DF_TF2 *
                    exp(Afloat2[i + 1] + c1float2[i + 1] * y1 + c2float2[i + 1] * y2);
                coupon_float2[i] =
                    ((double)(1) / float_cvg2[i + 1] * (ZC_float2[i] / ZC_float2[i + 1] - 1) +
                     spread2[i]) *
                    float_cvg2[i + 1];
                swap_num2 += coupon_float2[i] * ZC_float2[i + 1];
            }
            for (j = 0; j < nb_fix_coupon1; j++)
            {
                ZC_fix1[j] = DF_fix1[j] / DF_TF1 * exp(Afix1[j] + c1fix1[j] * y1 + c2fix1[j] * y2);
                swap_den1 += fix_cvg1[j] * ZC_fix1[j];
            }
            for (j = 0; j < nb_fix_coupon2; j++)
            {
                ZC_fix2[j] = DF_fix2[j] / DF_TF2 * exp(Afix2[j] + c1fix2[j] * y1 + c2fix2[j] * y2);
                swap_den2 += fix_cvg2[j] * ZC_fix2[j];
            }
            swap1 = swap_num1 / swap_den1;
            swap2 = swap_num2 / swap_den2;
            cov   = swap1 * swap2;
            s2 += w[k] * cov * exp(-y1 * y1 / 2 - y2 * y2 / 2) / (2 * Pi);
            cms1temp += w[k] * swap1 * exp(-y1 * y1 / 2 - y2 * y2 / 2) / (2 * Pi);
            cms2temp += w[k] * swap2 * exp(-y1 * y1 / 2 - y2 * y2 / 2) / (2 * Pi);
        }
        s1 += w[l] * s2;
        cms1 += w[l] * cms1temp;
        cms2 += w[l] * cms2temp;
    }
    s1   = s1 * y1r * y2r;
    cms1 = cms1 * y1r * y2r;
    cms2 = cms2 * y1r * y2r;
    return s1 - cms1 * cms2;

    if (Afloat1)
        free(Afloat1);
    if (Afix1)
        free(Afix1);
    if (c1float1)
        free(c1float1);
    if (c2float1)
        free(c2float1);
    if (c1fix1)
        free(c1fix1);
    if (c2fix1)
        free(c2fix1);
    if (ZC_float1)
        free(ZC_float1);
    if (ZC_fix1)
        free(ZC_fix1);
    if (coupon_float1)
        free(coupon_float1);
    if (Afloat2)
        free(Afloat2);
    if (Afix2)
        free(Afix2);
    if (c1float2)
        free(c1float2);
    if (c2float2)
        free(c2float2);
    if (c1fix2)
        free(c1fix2);
    if (c2fix2)
        free(c2fix2);
    if (ZC_float2)
        free(ZC_float2);
    if (ZC_fix2)
        free(ZC_fix2);
    if (coupon_float2)
        free(coupon_float2);

    if (x)
        free(x);
    if (w)
        free(w);
}

Err compute_CMScorrel_GL(
    double   Tf1,
    double   Tpay1,
    double*  DF_fix1,
    double*  fix_coupon_date1,
    double*  fix_cvg1,  // cvg[0]=0, the same size as fix_coupon_date
    int      nb_fix_coupon1,
    double*  DF_float1,
    double   DF_TF1,
    double*  float_coupon_date1,
    double*  float_cvg1,  // the same size as float_coupon_date
    int      nb_float_coupon1,
    double*  spread1,
    double   Tf2,
    double   Tpay2,
    double*  DF_fix2,
    double*  fix_coupon_date2,
    double*  fix_cvg2,  // cvg[0]=0, the same size as fix_coupon_date
    int      nb_fix_coupon2,
    double*  DF_float2,
    double   DF_TF2,
    double*  float_coupon_date2,
    double*  float_cvg2,  // the same size as float_coupon_date
    int      nb_float_coupon2,
    double*  spread2,
    double*  sigma_11,
    double*  sigma_21,
    double*  sigma_date1,
    int      nb_sigma1,
    double   lambda_11,
    double   lambda_21,
    double   rho1,
    double*  sigma_12,
    double*  sigma_22,
    double*  sigma_date2,
    int      nb_sigma2,
    double   lambda_12,
    double   lambda_22,
    double   rho2,
    int      nb_integration,
    double   bound,  // must be positive
    double** output
    //		  5 outputs: 1-correl
    //					 2-CMS1
    //					 3-CMS2
    //					 4-Vol(CMS1)
    //					 5-Vol(CMS2)
)
{
    double  y1, y1r, y1m, dy1, y2, y2r, y2m, dy2;
    double* x             = NULL;
    double* w             = NULL;
    double* Afloat1       = NULL;
    double* Afix1         = NULL;
    double* c1float1      = NULL;
    double* c2float1      = NULL;
    double* c1fix1        = NULL;
    double* c2fix1        = NULL;
    double* ZC_float1     = NULL;
    double* ZC_fix1       = NULL;
    double* coupon_float1 = NULL;
    double* Afloat2       = NULL;
    double* Afix2         = NULL;
    double* c1float2      = NULL;
    double* c2float2      = NULL;
    double* c1fix2        = NULL;
    double* c2fix2        = NULL;
    double* ZC_float2     = NULL;
    double* ZC_fix2       = NULL;
    double* coupon_float2 = NULL;
    double  res           = 0;
    int     i             = 0;
    int     j             = 0;
    int     k             = 0;
    int     l             = 0;
    double  sum_sigma_11  = 0;
    double  sum_sigma_21  = 0;
    double  sum_diagterm1 = 0;
    double  swap1;
    double  swap_num1;
    double  swap_den1;
    double  sum_sigma_12  = 0;
    double  sum_sigma_22  = 0;
    double  sum_diagterm2 = 0;
    double  swap2;
    double  swap_num2;
    double  swap_den2;
    double  cms1, cms1temp, varcms1, varcms1tmp;
    double  cms2, cms2temp, varcms2, varcms2tmp;
    double  cov, cov_temp;
    double  integ;

    Err err = NULL;

    // compute the different A(Tf,coupon_date) and Ci
    Afloat1       = (double*)calloc(nb_float_coupon1, sizeof(double));
    Afix1         = (double*)calloc(nb_fix_coupon1, sizeof(double));
    c1float1      = (double*)calloc(nb_float_coupon1, sizeof(double));
    c2float1      = (double*)calloc(nb_float_coupon1, sizeof(double));
    c1fix1        = (double*)calloc(nb_fix_coupon1, sizeof(double));
    c2fix1        = (double*)calloc(nb_fix_coupon1, sizeof(double));
    ZC_float1     = (double*)calloc(nb_float_coupon1, sizeof(double));
    ZC_fix1       = (double*)calloc(nb_fix_coupon1, sizeof(double));
    coupon_float1 = (double*)calloc(nb_float_coupon1, sizeof(double));

    Afloat2       = (double*)calloc(nb_float_coupon2, sizeof(double));
    Afix2         = (double*)calloc(nb_fix_coupon2, sizeof(double));
    c1float2      = (double*)calloc(nb_float_coupon2, sizeof(double));
    c2float2      = (double*)calloc(nb_float_coupon2, sizeof(double));
    c1fix2        = (double*)calloc(nb_fix_coupon2, sizeof(double));
    c2fix2        = (double*)calloc(nb_fix_coupon2, sizeof(double));
    ZC_float2     = (double*)calloc(nb_float_coupon2, sizeof(double));
    ZC_fix2       = (double*)calloc(nb_fix_coupon2, sizeof(double));
    coupon_float2 = (double*)calloc(nb_float_coupon2, sizeof(double));

    // compute abscissas and weights for Gauss Legendre integration
    x = (double*)calloc(nb_integration, sizeof(double));
    w = (double*)calloc(nb_integration, sizeof(double));
    GaussLeg(-1, 1, x - 1, w - 1, nb_integration);

    // compute B(t,T)=B(0,T)/B(0,t)*exp(A(t,T)+c1(t,T)*z1+c2(t,T)*z2) under QTp
    sum_sigma_11 = compute_sum_sigma(Tf1, sigma_11, sigma_date1, nb_sigma1, lambda_11);
    sum_sigma_21 = compute_sum_sigma(Tf1, sigma_21, sigma_date1, nb_sigma1, lambda_21);
    sum_diagterm1 =
        compute_sum_diagterm(Tf1, sigma_11, sigma_21, sigma_date1, nb_sigma1, lambda_11, lambda_21);

    sum_sigma_12 = compute_sum_sigma(Tf2, sigma_12, sigma_date2, nb_sigma2, lambda_12);
    sum_sigma_22 = compute_sum_sigma(Tf2, sigma_22, sigma_date2, nb_sigma2, lambda_22);
    sum_diagterm2 =
        compute_sum_diagterm(Tf2, sigma_12, sigma_22, sigma_date2, nb_sigma2, lambda_12, lambda_22);

    for (i = 0; i < nb_float_coupon1; i++)
    {
        Afloat1[i] = compute_A(
            Tf1,
            float_coupon_date1[i],
            sum_sigma_11,
            sum_sigma_21,
            sum_diagterm1,
            rho1,
            lambda_11,
            lambda_21,
            Tpay1);
        Afloat1[i] += log(DF_float1[i]);
        c1float1[i] = compute_c1(
            Tf1,
            float_coupon_date1[i],
            Tpay1,
            sigma_11,
            sigma_21,
            sigma_date1,
            nb_sigma1,
            rho1,
            lambda_11,
            lambda_21);
        c2float1[i] = compute_c2(
            Tf1, float_coupon_date1[i], Tpay1, sigma_21, sigma_date1, nb_sigma1, rho1, lambda_21);
    }
    for (i = 1; i < nb_fix_coupon1; i++)
    {
        Afix1[i] = compute_A(
            Tf1,
            fix_coupon_date1[i],
            sum_sigma_11,
            sum_sigma_21,
            sum_diagterm1,
            rho1,
            lambda_11,
            lambda_21,
            Tpay1);
        Afix1[i] += log(fix_cvg1[i] * DF_fix1[i]);
        c1fix1[i] = compute_c1(
            Tf1,
            fix_coupon_date1[i],
            Tpay1,
            sigma_11,
            sigma_21,
            sigma_date1,
            nb_sigma1,
            rho1,
            lambda_11,
            lambda_21);
        c2fix1[i] = compute_c2(
            Tf1, fix_coupon_date1[i], Tpay1, sigma_21, sigma_date1, nb_sigma1, rho1, lambda_21);
    }

    for (i = 0; i < nb_float_coupon2; i++)
    {
        Afloat2[i] = compute_A(
            Tf2,
            float_coupon_date2[i],
            sum_sigma_12,
            sum_sigma_22,
            sum_diagterm2,
            rho2,
            lambda_12,
            lambda_22,
            Tpay2);
        Afloat2[i] += log(DF_float2[i]);
        c1float2[i] = compute_c1(
            Tf2,
            float_coupon_date2[i],
            Tpay2,
            sigma_12,
            sigma_22,
            sigma_date2,
            nb_sigma2,
            rho2,
            lambda_12,
            lambda_22);
        c2float2[i] = compute_c2(
            Tf2, float_coupon_date2[i], Tpay2, sigma_22, sigma_date2, nb_sigma2, rho2, lambda_22);
    }
    for (i = 1; i < nb_fix_coupon2; i++)
    {
        Afix2[i] = compute_A(
            Tf2,
            fix_coupon_date2[i],
            sum_sigma_12,
            sum_sigma_22,
            sum_diagterm2,
            rho2,
            lambda_12,
            lambda_22,
            Tpay2);
        Afix2[i] += log(fix_cvg2[i] * DF_fix2[i]);
        c1fix2[i] = compute_c1(
            Tf2,
            fix_coupon_date2[i],
            Tpay2,
            sigma_12,
            sigma_22,
            sigma_date2,
            nb_sigma2,
            rho2,
            lambda_12,
            lambda_22);
        c2fix2[i] = compute_c2(
            Tf2, fix_coupon_date2[i], Tpay2, sigma_22, sigma_date2, nb_sigma2, rho2, lambda_22);
    }

    // the same bound in y1 and y2
    y1m = 0.5 * (-bound + bound);
    y1r = 0.5 * (bound + bound);
    y2m = y1m;
    y2r = y1r;
    y1  = 0;
    y2  = 0;

    cms1    = 0;
    cms2    = 0;
    varcms1 = 0;
    varcms2 = 0;
    cov     = 0;

    for (l = 0; l < nb_integration; l++)
    {
        cms1temp   = 0.0;
        cms2temp   = 0.0;
        varcms1tmp = 0.0;
        varcms2tmp = 0.0;
        cov_temp   = 0.0;

        dy1 = y1r * x[l];
        y1  = y1m + dy1;

        for (k = 0; k < nb_integration; k++)
        {
            dy2 = y2r * x[k];
            y2  = y2m + dy2;

            // compute the swap for z1=y1 and z2=y2
            swap1     = 0.0;
            swap_num1 = 0.0;
            swap_den1 = 0.0;
            swap2     = 0.0;
            swap_num2 = 0.0;
            swap_den2 = 0.0;

            ZC_float1[0] = exp(Afloat1[0] + c1float1[0] * y1 + c2float1[0] * y2);
            ZC_float2[0] = exp(Afloat2[0] + c1float2[0] * y1 + c2float2[0] * y2);

            for (i = 0; i < nb_float_coupon1 - 1; i++)
            {
                ZC_float1[i + 1] =
                    exp(Afloat1[i + 1] + c1float1[i + 1] * y1 + c2float1[i + 1] * y2);
                coupon_float1[i] = ZC_float1[i] - ZC_float1[i + 1] +
                                   spread1[i] * float_cvg1[i + 1] * ZC_float1[i + 1];
                swap_num1 += coupon_float1[i];
            }

            for (i = 0; i < nb_float_coupon2 - 1; i++)
            {
                ZC_float2[i + 1] =
                    exp(Afloat2[i + 1] + c1float2[i + 1] * y1 + c2float2[i + 1] * y2);
                coupon_float2[i] = ZC_float2[i] - ZC_float2[i + 1] +
                                   spread2[i] * float_cvg2[i + 1] * ZC_float2[i + 1];
                swap_num2 += coupon_float2[i];
            }

            for (j = 1; j < nb_fix_coupon1; j++)
            {
                swap_den1 += exp(Afix1[j] + c1fix1[j] * y1 + c2fix1[j] * y2);
            }

            for (j = 1; j < nb_fix_coupon2; j++)
            {
                swap_den2 += exp(Afix2[j] + c1fix2[j] * y1 + c2fix2[j] * y2);
            }

            swap1 = swap_num1 / swap_den1;
            swap2 = swap_num2 / swap_den2;

            integ = w[k] * exp(-y1 * y1 / 2.0 - y2 * y2 / 2.0);

            cms1temp += swap1 * integ;
            cms2temp += swap2 * integ;
            varcms1tmp += swap1 * swap1 * integ;
            varcms2tmp += swap2 * swap2 * integ;
            cov_temp += swap1 * swap2 * integ;
        }

        cms1 += w[l] * cms1temp;
        cms2 += w[l] * cms2temp;
        varcms1 += w[l] * varcms1tmp;
        varcms2 += w[l] * varcms2tmp;
        cov += w[l] * cov_temp;
    }

    cms1    = cms1 * y1r * y2r / (2.0 * Pi);
    cms2    = cms2 * y1r * y2r / (2.0 * Pi);
    varcms1 = varcms1 * y1r * y2r / (2.0 * Pi);
    varcms2 = varcms2 * y1r * y2r / (2.0 * Pi);
    cov     = cov * y1r * y2r / (2.0 * Pi);

    (*output)[1] = cms1;
    (*output)[2] = cms2;
    (*output)[3] = sqrt((varcms1 - cms1 * cms1) / Tf1);
    (*output)[4] = sqrt((varcms2 - cms2 * cms2) / Tf2);
    (*output)[0] =
        (cov - cms1 * cms2) / (sqrt(varcms1 - cms1 * cms1) * sqrt(varcms2 - cms2 * cms2));

    if (Afloat1)
        free(Afloat1);
    if (Afix1)
        free(Afix1);
    if (c1float1)
        free(c1float1);
    if (c2float1)
        free(c2float1);
    if (c1fix1)
        free(c1fix1);
    if (c2fix1)
        free(c2fix1);
    if (ZC_float1)
        free(ZC_float1);
    if (ZC_fix1)
        free(ZC_fix1);
    if (coupon_float1)
        free(coupon_float1);
    if (Afloat2)
        free(Afloat2);
    if (Afix2)
        free(Afix2);
    if (c1float2)
        free(c1float2);
    if (c2float2)
        free(c2float2);
    if (c1fix2)
        free(c1fix2);
    if (c2fix2)
        free(c2fix2);
    if (ZC_float2)
        free(ZC_float2);
    if (ZC_fix2)
        free(ZC_fix2);
    if (coupon_float2)
        free(coupon_float2);

    if (x)
        free(x);
    if (w)
        free(w);

    return err;
}

Err LGM2FCMS(
    char*     yc_name,
    char*     vol_curve_name,
    SrtUndPtr und,
    double    alpha,
    double    gamma,
    double    rho,
    int       do_calib,
    long      lTfixing,
    long      lStartDate,
    long      lEndDate,
    char*     fixFreq,
    char*     fixBasis,
    char*     refRateName,
    long      lTpay,
    long      lToday,
    int       nb_integration,
    double    bound,
    char*     ProductTypeStr,
    double    strike,
    char*(
        GetCpdAutocalCashVol)(char* lpszVolCurveName, double dStart, double dEnd, double dStrike, int bIsACap, char* lpszRefRateCode, double* pdVol, double* pdPower),
    double* price)
{
    Err err = NULL;

    int     sigma_n;
    double* sigma_date = NULL;
    double* sigma      = NULL;
    double* sigma_1    = NULL;
    double* sigma_2    = NULL;

    int*    nb_UNO             = NULL;
    long*   primary_ex_date1   = NULL;
    double* primary_strike1    = NULL;
    char**  end_tenor_primary1 = NULL;

    double* secundary_strike    = NULL;
    char*   secundary_basis     = NULL;
    char*   secundary_freq      = NULL;
    long*   secundary_ex_date   = NULL;
    char**  end_tenor_secundary = NULL;

    int     nlam;
    double* lam_time1  = NULL;
    double* lam1       = NULL;
    double* lam_shift1 = NULL;

    double*    lam2     = NULL;
    double*    tau_date = NULL;
    double*    tau      = NULL;
    int        tau_n;
    String     und_name;
    SrtMdlType mdl_type;
    SrtMdlDim  mdl_dim;

    double fixed_tau = 15;
    double Tf;
    double Tp;
    double Tstart;
    double Tend;

    int                  fix_lambda = 0;
    int                  one2F      = 2;
    int                  one_f_equi = 0;
    cpd_diag_calib_param paraml1;
    cpd_diag_calib_param params1;
    diag_calib_lm_params lm_params1;

    int  tmp;
    char char_tmp[XL_MAX_ID_LEN];
    int  k, n;

    // TO BUILD
    SwapDP swapdp;
    int    float_nb_dates, float_nb_pay_dates;
    long * float_fixing_dates = NULL, *float_start_dates = NULL, *float_end_dates = NULL,
         *float_pay_dates = NULL;
    double *float_cvgs = NULL, *float_spreads = NULL, *float_pay_times = NULL;
    int     fix_nb_dates, fix_nb_pay_dates;
    long *  fix_start_dates = NULL, *fix_end_dates = NULL, *fix_pay_dates = NULL;
    double *fix_cvgs = NULL, *fix_pay_times = NULL;

    double* float_coupon_date = NULL;
    double* float_cvg         = NULL;
    double* fix_coupon_date   = NULL;
    double* fix_cvg           = NULL;
    double* DF_float          = NULL;
    double* DF_fix            = NULL;
    double* spread            = NULL;
    double  DF_Tf;
    double  lambda_1;
    double  lambda_2;
    int     nb_fix_coupon;
    int     nb_float_coupon;

    SrtReceiverType ProductType;

    // END "TO BUILD"
    if (lTfixing == 0)
    {
        lTfixing = add_unit(lStartDate, -2, 0, 0);
    }
    if (lTpay == 0)
    {
        lTpay = lStartDate;
    }

    Tf     = (lTfixing - lToday) * YEARS_IN_DAY;
    Tp     = (lTpay - lToday) * YEARS_IN_DAY;
    Tstart = (lStartDate - lToday) * YEARS_IN_DAY;
    Tend   = (lEndDate - lToday) * YEARS_IN_DAY;

    // by default bound is 10 and nb_integration is 200
    if (bound == 0)
    {
        bound = 10;
    }
    if (nb_integration == 0)
    {
        nb_integration = 200;
    }

    if (err = interp_rec_pay(ProductTypeStr, &ProductType))
    {
        return err;
    }

    nlam = 1;
    lam1 = (double*)calloc(nlam, sizeof(double));
    lam2 = (double*)calloc(nlam, sizeof(double));

    if (do_calib == 1)
    {
        if (alpha == 0)
        {
            alpha = 1.39;
        }
        if (gamma == 0)
        {
            gamma = 0.21;
        }
        if (rho == 0)
        {
            rho = -0.85;
        }

        //----------------Choice of the Instruments------------------------

        primary_ex_date1  = (long*)calloc(1, sizeof(long));
        secundary_ex_date = (long*)calloc(1, sizeof(long));
        secundary_freq    = (char*)calloc(256, sizeof(char));
        primary_strike1   = (double*)calloc(1, sizeof(double));
        secundary_strike  = (double*)calloc(1, sizeof(double));
        secundary_basis   = (char*)calloc(256, sizeof(char));

        nb_UNO                = (int*)calloc(1, sizeof(int));
        end_tenor_primary1    = (char**)calloc(1, sizeof(char*));
        end_tenor_primary1[0] = (char*)calloc(256, sizeof(char));

        end_tenor_secundary    = (char**)calloc(1, sizeof(char*));
        end_tenor_secundary[0] = (char*)calloc(256, sizeof(char));

        // Choice of the strikes
        if (ProductType == SRT_FORWARD)
        {
            primary_strike1[0]  = 0;  // ATM strike
            secundary_strike[0] = 0;
        }
        else
        {
            primary_strike1[0]  = strike;
            secundary_strike[0] = strike;
        }

        nb_UNO    = calloc(1, sizeof(int));
        nb_UNO[0] = 1;

        primary_ex_date1[0] = (long)(lTfixing);

        cpd_calib_set_default_param(&paraml1);
        cpd_calib_set_default_param(&params1);
        diag_calib_lm_params_set_default_param(&lm_params1);
        paraml1.strike_type = 2;  //"SWAP"=2
        params1.strike_type = 2;

        lam_time1  = (double*)calloc(nlam, sizeof(double));
        lam_shift1 = (double*)calloc(nlam, sizeof(double));

        lam_time1[0]  = 1.0;
        lam_shift1[0] = 0;
        lam1[0]       = 0.05;

        // Tenor of the first instruments
        if (Tend - Tstart < 1)
        {
            tmp = (int)((Tend - Tstart + 1.0e-08) * 12);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "M");
            strcpy(end_tenor_primary1[0], char_tmp);
        }
        else
        {
            tmp = (int)(Tend - Tstart + 1.0e-08);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "Y");
            strcpy(end_tenor_primary1[0], char_tmp);
        }

        // CALIBRATION OF THE SECOND INSTRUMENTS

        strcpy(secundary_basis, "MM");
        if (Tp - Tstart < 0.49)
        {
            strcpy(secundary_freq, "Q");
        }
        else if (Tp - Tstart < 0.98)
        {
            strcpy(secundary_freq, "S");
        }
        else
        {
            strcpy(secundary_freq, "A");
        }

        // Tenor of the second instruments
        if (Tp - Tstart < 1)
        {
            tmp = (int)((Tp - Tstart + 1.0e-08) * 12);
            if (tmp != 0)
            {
                sprintf(char_tmp, "%d", tmp);
                strcat(char_tmp, "M");
                strcpy(end_tenor_secundary[0], char_tmp);
            }
            else
            {
                fix_lambda = 1;
                lam1[0]    = 1 / fixed_tau;
            }
        }
        else
        {
            tmp = (int)(Tp - Tstart + 1.0e-08);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "Y");
            strcpy(end_tenor_secundary[0], char_tmp);
        }

        secundary_ex_date[0] = (long)(lTfixing);

        // Calibration
        err = cpd_calib_diagonal_dlm(
            //	Market
            yc_name,
            vol_curve_name,
            GetCpdAutocalCashVol,
            // Get vol ref
            refRateName,
            // Long Instruments //
            fixFreq,
            fixBasis,
            refRateName,
            1,
            primary_ex_date1,
            nb_UNO,
            end_tenor_primary1,
            (long)(lEndDate),
            primary_strike1,
            &paraml1,
            // Short Instruments //
            secundary_freq,
            secundary_basis,
            refRateName,
            1,
            secundary_ex_date,
            nb_UNO,
            end_tenor_secundary,
            (long)(lTpay),
            secundary_strike,
            NULL,
            &params1,
            //	Model
            fix_lambda,
            one_f_equi,
            nlam,
            lam_time1,
            lam1,
            lam_shift1,
            2,
            alpha,
            gamma,
            rho,
            //	Shift Parameters
            0,
            0,
            0,
            0,
            0,
            0,
            //	Output //
            &sigma_n,
            &sigma_date,
            &sigma,
            //	Parameters //
            &lm_params1,
            //	Calibration instrument data //
            NULL);
    }
    // END OF CALIBRATION
    else
    {
        // Gets the Model type (LGM, Cheyette,...)
        err = get_underlying_mdltype(und, &mdl_type);

        // Gets the Number of Factors in the Model
        err = get_underlying_mdldim(und, &mdl_dim);

        // Gets Today from underlying
        und_name = get_underlying_name(und);

        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            return err;
        }

        if (tau_n > 1)
        {
            err = "Constant Tau required";
            goto FREE_RETURN;
        }
        lam1[0] = 1.0 / tau[0];
    }

    lam2[0] = lam1[0] + gamma;
    sigma_1 = (double*)calloc(sigma_n, sizeof(double));
    sigma_2 = (double*)calloc(sigma_n, sizeof(double));
    for (k = 0; k < sigma_n; k++)
    {
        sigma_1[k] = sigma[k];
        sigma_2[k] = sigma_1[k] * alpha;
    }

    //------------------- - BUILD - ----------------------
    DF_Tf = swp_f_df(lToday, lTfixing, yc_name);

    lambda_1 = lam1[0];
    lambda_2 = lam2[0];

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
    nb_fix_coupon   = fix_nb_pay_dates;
    nb_float_coupon = float_nb_pay_dates;

    float_coupon_date = (double*)calloc(nb_float_coupon, sizeof(double));
    DF_float          = (double*)calloc(nb_float_coupon, sizeof(double));
    float_cvg         = (double*)calloc(nb_float_coupon, sizeof(double));
    spread            = (double*)calloc(nb_float_coupon - 1, sizeof(double));
    fix_coupon_date   = (double*)calloc(nb_fix_coupon, sizeof(double));
    DF_fix            = (double*)calloc(nb_fix_coupon, sizeof(double));
    fix_cvg           = (double*)calloc(nb_fix_coupon, sizeof(double));

    float_coupon_date[0] = (float_pay_dates[0] - lToday) * YEARS_IN_DAY;
    fix_coupon_date[0]   = (fix_pay_dates[0] - lToday) * YEARS_IN_DAY;
    DF_float[0]          = swp_f_df(lToday, float_pay_dates[0], yc_name);
    DF_fix[0]            = swp_f_df(lToday, fix_pay_dates[0], yc_name);
    float_cvg[0]         = 0;
    fix_cvg[0]           = 0;
    for (k = 1; k < nb_float_coupon; k++)
    {
        float_coupon_date[k] = (float_pay_dates[k] - lToday) * YEARS_IN_DAY;
        DF_float[k]          = swp_f_df(lToday, float_pay_dates[k], yc_name);
        float_cvg[k]         = float_cvgs[k - 1];
        spread[k - 1]        = float_spreads[k - 1];
    }
    for (k = 1; k < nb_fix_coupon; k++)
    {
        fix_coupon_date[k] = (fix_pay_dates[k] - lToday) * YEARS_IN_DAY;
        DF_fix[k]          = swp_f_df(lToday, fix_pay_dates[k], yc_name);
        fix_cvg[k]         = fix_cvgs[k - 1];
    }

    //-----------------------------------------------------------------
    //----------------End Of Build Swap Schedule-----------------------
    //-----------------------------------------------------------------

    *price = compute_CMS_GaussLeg(
        Tf,
        Tp,
        DF_fix,
        fix_coupon_date,
        fix_cvg,  // cvg[0]=0, the same size as fix_coupon_date
        nb_fix_coupon,
        DF_float,
        DF_Tf,
        float_coupon_date,
        float_cvg,  // the same size as float_coupon_date
        nb_float_coupon,
        float_spreads,
        sigma_1,
        sigma_2,
        sigma_date,
        sigma_n,
        lambda_1,
        lambda_2,
        rho,
        nb_integration,
        bound,
        ProductType,
        strike);

FREE_RETURN:

    if (float_coupon_date)
        free(float_coupon_date);
    if (fix_coupon_date)
        free(fix_coupon_date);
    if (DF_float)
        free(DF_float);
    if (DF_fix)
        free(DF_fix);
    if (float_cvg)
        free(float_cvg);
    if (spread)
        free(spread);
    if (fix_cvg)
        free(fix_cvg);
    if (sigma_1)
        free(sigma_1);
    if (sigma_2)
        free(sigma_2);

    if (sigma_date)
        free(sigma_date);
    if (sigma)
        free(sigma);
    if (nb_UNO)
        free(nb_UNO);
    if (primary_ex_date1)
        free(primary_ex_date1);
    if (primary_strike1)
        free(primary_strike1);
    if (end_tenor_primary1)
        free(end_tenor_primary1);
    if (secundary_strike)
        free(secundary_strike);
    if (secundary_ex_date)
        free(secundary_ex_date);
    if (end_tenor_secundary)
        free(end_tenor_secundary);
    if (lam_time1)
        free(lam_time1);
    if (lam1)
        free(lam1);
    if (lam_shift1)
        free(lam_shift1);
    if (lam2)
        free(lam2);
    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);

    if (secundary_freq)
        free(secundary_freq);
    if (secundary_basis)
        free(secundary_basis);

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

    return err;
}

Err LGM2FDRS(
    char*     yc_name,
    char*     vol_curve_name,
    SrtUndPtr und,
    double    alpha,
    double    gamma,
    double    rho,
    int       do_calib,
    long      lTfixing,
    long      lStartDate,
    long      lEndDate,
    char*     fixFreq,
    char*     fixBasis,
    char*     refRateName,
    long      lTpay,
    long      lToday,
    int       nb_integration,
    double    bound,
    char*(
        GetCpdAutocalCashVol)(char* lpszVolCurveName, double dStart, double dEnd, double dStrike, int bIsACap, char* lpszRefRateCode, double* pdVol, double* pdPower),
    double** output)
// 3 outputs:
// 1-CMS(Tp)-CMS(Ts) / 2-CMS(Tp) / 3-CMS(Ts)
{
    char*  ProductType_loc = NULL;
    double strike          = 0;
    double price;
    Err    err = NULL;

    ProductType_loc = (char*)calloc(256, sizeof(char));
    strcpy(ProductType_loc, "RATE");

    err = LGM2FCMS(
        yc_name,
        vol_curve_name,
        und,
        alpha,
        gamma,
        rho,
        do_calib,
        lTfixing,
        lStartDate,
        lEndDate,
        fixFreq,
        fixBasis,
        refRateName,
        lTpay,
        lToday,
        nb_integration,
        bound,
        ProductType_loc,
        strike,
        GetCpdAutocalCashVol,
        &price);
    (*output)[1] = price;
    err          = LGM2FCMS(
        yc_name,
        vol_curve_name,
        und,
        alpha,
        gamma,
        rho,
        do_calib,
        lTfixing,
        lStartDate,
        lEndDate,
        fixFreq,
        fixBasis,
        refRateName,
        lStartDate,  // instead of lTpay
        lToday,
        nb_integration,
        bound,
        ProductType_loc,
        strike,
        GetCpdAutocalCashVol,
        &price);
    (*output)[2] = price;
    (*output)[0] = (*output)[1] - (*output)[2];
    return err;
    if (ProductType_loc)
        free(ProductType_loc);
}

Err LGM2FCMSCorrel(
    char*     yc_name,
    char*     vol_curve_name,
    SrtUndPtr und,
    double    alpha,
    double    gamma,
    double    rho,
    int       do_calib,
    long      lTfixing,
    long      lStartDate,
    long      lEndDate1,
    char*     fixFreq1,
    char*     fixBasis1,
    char*     refRateName,
    long      lTpay,
    long      lToday,
    long      lEndDate2,
    char*     fixFreq2,
    char*     fixBasis2,
    int       nb_integration,
    double    bound,
    char*(
        GetCpdAutocalCashVol)(char* lpszVolCurveName, double dStart, double dEnd, double dStrike, int bIsACap, char* lpszRefRateCode, double* pdVol, double* pdPower),
    double** output)  // dimension=7: OUTPUT:1-Correl / 2-CMS1 / 3-CMS2 / 4-VAR(CMS1) / 5-VAR(CMS2)
                      // / 6-sigma1 / 7-lambda1
{
    Err     err = NULL;
    int     sigma_n;
    double* sigma_date = NULL;
    double* sigma      = NULL;
    double* sigma_1    = NULL;
    double* sigma_2    = NULL;

    int*    nb_UNO             = NULL;
    long*   primary_ex_date1   = NULL;
    double* primary_strike1    = NULL;
    char**  end_tenor_primary1 = NULL;

    double* secundary_strike    = NULL;
    long*   secundary_ex_date   = NULL;
    char**  end_tenor_secundary = NULL;

    int     nlam;
    double* lam_time1  = NULL;
    double* lam1       = NULL;
    double* lam_shift1 = NULL;

    double*    lam2     = NULL;
    double*    tau_date = NULL;
    double*    tau      = NULL;
    int        tau_n;
    String     und_name;
    SrtMdlType mdl_type;
    SrtMdlDim  mdl_dim;

    double fixed_tau = 15;
    double Tf;
    double Tp;
    double Tstart;
    double Tend1;
    double Tend2;

    int                  fix_lambda = 0;
    int                  one2F      = 2;
    int                  one_f_equi = 0;
    cpd_diag_calib_param paraml1;
    cpd_diag_calib_param params1;
    diag_calib_lm_params lm_params1;

    int  tmp;
    char char_tmp[XL_MAX_ID_LEN];
    int  k;

    double* output_tmp = NULL;
    if (lTfixing == 0)
    {
        lTfixing = add_unit(lStartDate, -2, 0, 0);
    }
    if (lTpay == 0)
    {
        lTpay = lStartDate;
    }

    // By default, the second CMS is the same as the first one
    if (strlen(fixFreq2) == 0)
    {
        fixFreq2 = fixFreq1;
    }
    if (strlen(fixBasis2) == 0)
    {
        fixBasis2 = fixBasis1;
    }
    if (lEndDate2 == 0)
    {
        lEndDate2 = lEndDate1;
    }

    Tf     = (lTfixing - lToday) * YEARS_IN_DAY;
    Tp     = (lTpay - lToday) * YEARS_IN_DAY;
    Tstart = (lStartDate - lToday) * YEARS_IN_DAY;

    // by default bound is 10 and nb_integration is 200
    if (bound == 0)
    {
        bound = 10;
    }
    if (nb_integration == 0)
    {
        nb_integration = 200;
    }

    output_tmp = (double*)calloc(5, sizeof(double));
    // OUTPUT:1-Correl / 2-CMS1 / 3-CMS2 / 4-VAR(CMS1) / 5-VAR(CMS2)

    nlam = 1;
    lam1 = (double*)calloc(nlam, sizeof(double));
    lam2 = (double*)calloc(nlam, sizeof(double));

    if (do_calib == 1)
    {
        if (alpha == 0)
        {
            alpha = 1.39;
        }
        if (gamma == 0)
        {
            gamma = 0.21;
        }
        if (rho == 0)
        {
            rho = -0.85;
        }

        //----------------Choice of the Instruments------------------------

        primary_ex_date1  = (long*)calloc(1, sizeof(long));
        secundary_ex_date = (long*)calloc(1, sizeof(long));
        primary_strike1   = (double*)calloc(1, sizeof(double));
        secundary_strike  = (double*)calloc(1, sizeof(double));

        nb_UNO                = (int*)calloc(1, sizeof(int));
        end_tenor_primary1    = (char**)calloc(1, sizeof(char*));
        end_tenor_primary1[0] = (char*)calloc(256, sizeof(char));

        end_tenor_secundary    = (char**)calloc(1, sizeof(char*));
        end_tenor_secundary[0] = (char*)calloc(256, sizeof(char));

        // Choice of the strikes
        primary_strike1[0]  = 0;
        secundary_strike[0] = 0;

        nb_UNO    = calloc(1, sizeof(int));
        nb_UNO[0] = 1;

        // Tenor of the first instruments
        Tend1 = (lEndDate1 - lToday) * YEARS_IN_DAY;
        if (Tend1 - Tstart < 1)
        {
            tmp = (int)((Tend1 - Tstart + 1.0e-08) * 12);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "M");
            strcpy(end_tenor_primary1[0], char_tmp);
        }
        else
        {
            tmp = (int)(Tend1 - Tstart + 1.0e-08);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "Y");
            strcpy(end_tenor_primary1[0], char_tmp);
        }

        primary_ex_date1[0] = (long)(lTfixing);

        cpd_calib_set_default_param(&paraml1);
        cpd_calib_set_default_param(&params1);
        diag_calib_lm_params_set_default_param(&lm_params1);
        paraml1.strike_type = 2;  //"SWAP"=2
        params1.strike_type = 2;

        lam_time1  = (double*)calloc(nlam, sizeof(double));
        lam_shift1 = (double*)calloc(nlam, sizeof(double));

        lam_time1[0]  = 1.0;
        lam_shift1[0] = 0;
        lam1[0]       = 0.05;

        // CALIBRATION OF THE SECOND CMS

        // Tenor of the second instruments
        Tend2 = (lEndDate2 - lToday) * YEARS_IN_DAY;
        if (Tend2 - Tstart < 1)
        {
            tmp = (int)((Tend2 - Tstart + 1.0e-08) * 12);
            if (tmp != 0)
            {
                sprintf(char_tmp, "%d", tmp);
                strcat(char_tmp, "M");
                strcpy(end_tenor_secundary[0], char_tmp);
            }
            else
            {
                fix_lambda = 1;
                lam1[0]    = 1 / fixed_tau;
            }
        }
        else
        {
            tmp = (int)(Tend2 - Tstart + 1.0e-08);
            sprintf(char_tmp, "%d", tmp);
            strcat(char_tmp, "Y");
            strcpy(end_tenor_secundary[0], char_tmp);
        }

        secundary_ex_date[0] = (long)(lTfixing);

        // Calibration
        err = cpd_calib_diagonal_dlm(
            //	Market
            yc_name,
            vol_curve_name,
            GetCpdAutocalCashVol,
            // Get vol ref
            refRateName,
            // Long Instruments //
            fixFreq1,
            fixBasis1,
            refRateName,
            1,
            primary_ex_date1,
            nb_UNO,
            end_tenor_primary1,
            (long)(lEndDate1),
            primary_strike1,
            &paraml1,
            // Short Instruments //
            fixFreq2,
            fixBasis2,
            refRateName,
            1,
            secundary_ex_date,
            nb_UNO,
            end_tenor_secundary,
            (long)(lEndDate2),
            secundary_strike,
            NULL,
            &params1,
            //	Model
            fix_lambda,
            one_f_equi,
            nlam,
            lam_time1,
            lam1,
            lam_shift1,
            2,
            alpha,
            gamma,
            rho,
            //	Shift Parameters
            0,
            0,
            0,
            0,
            0,
            0,
            //	Output //
            &sigma_n,
            &sigma_date,
            &sigma,
            //	Parameters //
            &lm_params1,
            //	Calibration instrument data //
            NULL);
    }
    // END OF CALIBRATION
    else
    {
        /* Gets the Model type (LGM, Cheyette,...) */
        err = get_underlying_mdltype(und, &mdl_type);

        /* Gets the Number of Factors in the Model */
        err = get_underlying_mdldim(und, &mdl_dim);

        /* Gets Today from underlying */
        und_name = get_underlying_name(und);

        err = Get_LGM_TermStructureOneOrTwoFact(
            und_name, &sigma_date, &sigma, &sigma_n, &tau_date, &tau, &tau_n, &alpha, &gamma, &rho);
        if (err)
        {
            return err;
        }

        if (tau_n > 1)
        {
            err = "Constant Tau required";
            goto FREE_RETURN;
        }
        lam1[0] = 1.0 / tau[0];
    }

    lam2[0] = lam1[0] + gamma;
    sigma_1 = (double*)calloc(sigma_n, sizeof(double));
    sigma_2 = (double*)calloc(sigma_n, sizeof(double));
    for (k = 0; k < sigma_n; k++)
    {
        sigma_1[k] = sigma[k];
        sigma_2[k] = sigma_1[k] * alpha;
    }

    err = LGM2FCMSCorrelBuild(
        yc_name,
        und_name,
        lTfixing,
        lStartDate,
        lEndDate1,
        lEndDate2,
        lTpay,
        lToday,
        fixFreq1,
        fixBasis1,
        refRateName,
        fixFreq2,
        fixBasis2,
        sigma_1,
        sigma_2,
        sigma_date,
        sigma_n,
        lam1[0],
        lam2[0],
        alpha,
        gamma,
        rho,
        nb_integration,
        bound,
        &output_tmp);

    if (err)
        goto FREE_RETURN;

    (*output)[0] = output_tmp[0];
    (*output)[1] = output_tmp[1];
    (*output)[2] = output_tmp[2];
    (*output)[3] = output_tmp[3];
    (*output)[4] = output_tmp[4];
    (*output)[5] = sigma_1[0];
    (*output)[6] = lam1[0];

FREE_RETURN:

    if (primary_ex_date1)
        free(primary_ex_date1);
    if (secundary_ex_date)
        free(secundary_ex_date);
    if (primary_strike1)
        free(primary_strike1);
    if (secundary_strike)
        free(secundary_strike);
    if (nb_UNO)
        free(nb_UNO);
    if (sigma)
        free(sigma);
    if (sigma_1)
        free(sigma_1);
    if (sigma_2)
        free(sigma_2);
    if (sigma_date)
        free(sigma_date);
    if (lam1)
        free(lam1);
    if (lam_time1)
        free(lam_time1);
    if (lam_shift1)
        free(lam_shift1);
    if (end_tenor_primary1)
        free(end_tenor_primary1);
    if (lam2)
        free(lam2);
    if (end_tenor_secundary)
        free(end_tenor_secundary);
    if (tau_date)
        free(tau_date);
    if (tau)
        free(tau);
    if (output_tmp)
        free(output_tmp);

    return err;
}

Err LGM2FCMSCorrelBuild(
    char*    yc_name,
    char*    und,
    long     lTfixing,
    long     lStartDate,
    long     lEndDate1,
    long     lEndDate2,
    long     lTpay,
    long     lToday,
    char*    fixFreq1,
    char*    fixBasis1,
    char*    refRateName,
    char*    fixFreq2,
    char*    fixBasis2,
    double*  sigma_1,
    double*  sigma_2,
    double*  sigma_date,
    int      sigma_n,
    double   lambda1,
    double   lambda2,
    double   alpha,
    double   gamma,
    double   rho,
    int      nb_integration,
    double   bound,
    double** output)
{
    Err err = NULL;

    SwapDP swapdp1;
    int    float_nb_dates1, float_nb_pay_dates1;
    long * float_fixing_dates1 = NULL, *float_start_dates1 = NULL, *float_end_dates1 = NULL,
         *float_pay_dates1 = NULL;
    double *float_cvgs1 = NULL, *float_spreads1 = NULL, *float_pay_times1 = NULL;
    int     fix_nb_dates1, fix_nb_pay_dates1;
    long *  fix_start_dates1 = NULL, *fix_end_dates1 = NULL, *fix_pay_dates1 = NULL;
    double *fix_cvgs1 = NULL, *fix_pay_times1 = NULL;

    SwapDP swapdp2;
    int    float_nb_dates2, float_nb_pay_dates2;
    long * float_fixing_dates2 = NULL, *float_start_dates2 = NULL, *float_end_dates2 = NULL,
         *float_pay_dates2 = NULL;
    double *float_cvgs2 = NULL, *float_spreads2 = NULL, *float_pay_times2 = NULL;
    int     fix_nb_dates2, fix_nb_pay_dates2;
    long *  fix_start_dates2 = NULL, *fix_end_dates2 = NULL, *fix_pay_dates2 = NULL;
    double *fix_cvgs2 = NULL, *fix_pay_times2 = NULL;

    int     k, n;
    double* float_coupon_date1 = NULL;
    double* float_cvg1         = NULL;
    double* fix_coupon_date1   = NULL;
    double* fix_cvg1           = NULL;
    double* DF_float1          = NULL;
    double* DF_fix1            = NULL;
    double* spread1            = NULL;
    double  Tf;
    double  DF_Tf1;
    double  Tp;
    int     nb_fix_coupon1;
    int     nb_float_coupon1;

    double* float_coupon_date2 = NULL;
    double* float_cvg2         = NULL;
    double* fix_coupon_date2   = NULL;
    double* fix_cvg2           = NULL;
    double* DF_float2          = NULL;
    double* DF_fix2            = NULL;
    double* spread2            = NULL;
    double  DF_Tf2;
    int     nb_fix_coupon2;
    int     nb_float_coupon2;

    double* output_tmp = NULL;

    Tf     = (lTfixing - lToday) * YEARS_IN_DAY;
    DF_Tf1 = swp_f_df(lToday, lTfixing, yc_name);
    Tp     = (lTpay - lToday) * YEARS_IN_DAY;

    DF_Tf2 = swp_f_df(lToday, lTfixing, yc_name);

    //----------------------------------------------------------------------------
    //---------SWAPDP Build Swap Schedule-----------------------------------------
    //----------------------------------------------------------------------------

    err = swp_f_initSwapDP(lStartDate, lEndDate1, fixFreq1, fixBasis1, &swapdp1);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_initSwapDP(lStartDate, lEndDate2, fixFreq2, fixBasis2, &swapdp2);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp1,
        lToday,
        refRateName,
        &float_pay_dates1,
        &float_nb_pay_dates1,
        &float_fixing_dates1,
        &float_start_dates1,
        &float_end_dates1,
        &float_cvgs1,
        &float_spreads1,
        &float_nb_dates1);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FloatLegDatesCoveragesAndSpreads(
        &swapdp2,
        lToday,
        refRateName,
        &float_pay_dates2,
        &float_nb_pay_dates2,
        &float_fixing_dates2,
        &float_start_dates2,
        &float_end_dates2,
        &float_cvgs2,
        &float_spreads2,
        &float_nb_dates2);
    if (err)
    {
        goto FREE_RETURN;
    }
    float_pay_times1 = dvector(0, float_nb_pay_dates1 - 1);
    float_pay_times2 = dvector(0, float_nb_pay_dates2 - 1);
    for (n = 0; n < float_nb_pay_dates1; ++n)
    {
        float_pay_times1[n] = (float_pay_dates1[n] - lToday) * YEARS_IN_DAY;
    }
    for (n = 0; n < float_nb_pay_dates2; ++n)
    {
        float_pay_times2[n] = (float_pay_dates2[n] - lToday) * YEARS_IN_DAY;
    }
    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp1,
        lToday,
        &fix_pay_dates1,
        &fix_nb_pay_dates1,
        &fix_start_dates1,
        &fix_end_dates1,
        &fix_cvgs1,
        &fix_nb_dates1);
    if (err)
    {
        goto FREE_RETURN;
    }
    err = swp_f_make_FixedLegDatesAndCoverages(
        &swapdp2,
        lToday,
        &fix_pay_dates2,
        &fix_nb_pay_dates2,
        &fix_start_dates2,
        &fix_end_dates2,
        &fix_cvgs2,
        &fix_nb_dates2);
    if (err)
    {
        goto FREE_RETURN;
    }
    fix_pay_times1 = dvector(0, fix_nb_pay_dates1 - 1);
    fix_pay_times2 = dvector(0, fix_nb_pay_dates2 - 1);
    for (n = 0; n < fix_nb_pay_dates1; ++n)
    {
        fix_pay_times1[n] = (fix_pay_dates1[n] - lToday) * YEARS_IN_DAY;
    }
    for (n = 0; n < fix_nb_pay_dates2; ++n)
    {
        fix_pay_times2[n] = (fix_pay_dates2[n] - lToday) * YEARS_IN_DAY;
    }
    nb_fix_coupon1   = fix_nb_pay_dates1;
    nb_float_coupon1 = float_nb_pay_dates1;
    nb_fix_coupon2   = fix_nb_pay_dates2;
    nb_float_coupon2 = float_nb_pay_dates2;

    float_coupon_date1 = (double*)calloc(nb_float_coupon1, sizeof(double));
    DF_float1          = (double*)calloc(nb_float_coupon1, sizeof(double));
    float_cvg1         = (double*)calloc(nb_float_coupon1, sizeof(double));
    spread1            = (double*)calloc(nb_float_coupon1 - 1, sizeof(double));
    fix_coupon_date1   = (double*)calloc(nb_fix_coupon1, sizeof(double));
    DF_fix1            = (double*)calloc(nb_fix_coupon1, sizeof(double));
    fix_cvg1           = (double*)calloc(nb_fix_coupon1, sizeof(double));

    float_coupon_date2 = (double*)calloc(nb_float_coupon2, sizeof(double));
    DF_float2          = (double*)calloc(nb_float_coupon2, sizeof(double));
    float_cvg2         = (double*)calloc(nb_float_coupon2, sizeof(double));
    spread2            = (double*)calloc(nb_float_coupon2 - 1, sizeof(double));
    fix_coupon_date2   = (double*)calloc(nb_fix_coupon2, sizeof(double));
    DF_fix2            = (double*)calloc(nb_fix_coupon2, sizeof(double));
    fix_cvg2           = (double*)calloc(nb_fix_coupon2, sizeof(double));

    float_coupon_date1[0] = (float_pay_dates1[0] - lToday) * YEARS_IN_DAY;
    fix_coupon_date1[0]   = (fix_pay_dates1[0] - lToday) * YEARS_IN_DAY;
    DF_float1[0]          = swp_f_df(lToday, float_pay_dates1[0], yc_name);
    DF_fix1[0]            = swp_f_df(lToday, fix_pay_dates1[0], yc_name);
    float_cvg1[0]         = 0;
    fix_cvg1[0]           = 0;
    for (k = 1; k < nb_float_coupon1; k++)
    {
        float_coupon_date1[k] = (float_pay_dates1[k] - lToday) * YEARS_IN_DAY;
        DF_float1[k]          = swp_f_df(lToday, float_pay_dates1[k], yc_name);
        float_cvg1[k]         = float_cvgs1[k - 1];
        spread1[k - 1]        = float_spreads1[k - 1];
    }
    for (k = 1; k < nb_fix_coupon1; k++)
    {
        fix_coupon_date1[k] = (fix_pay_dates1[k] - lToday) * YEARS_IN_DAY;
        DF_fix1[k]          = swp_f_df(lToday, fix_pay_dates1[k], yc_name);
        fix_cvg1[k]         = fix_cvgs1[k - 1];
    }

    float_coupon_date2[0] = (float_pay_dates2[0] - lToday) * YEARS_IN_DAY;
    fix_coupon_date2[0]   = (fix_pay_dates2[0] - lToday) * YEARS_IN_DAY;
    DF_float2[0]          = swp_f_df(lToday, float_pay_dates2[0], yc_name);
    DF_fix2[0]            = swp_f_df(lToday, fix_pay_dates2[0], yc_name);
    float_cvg2[0]         = 0;
    fix_cvg2[0]           = 0;
    for (k = 1; k < nb_float_coupon2; k++)
    {
        float_coupon_date2[k] = (float_pay_dates2[k] - lToday) * YEARS_IN_DAY;
        DF_float2[k]          = swp_f_df(lToday, float_pay_dates2[k], yc_name);
        float_cvg2[k]         = float_cvgs2[k - 1];
        spread2[k - 1]        = float_spreads2[k - 1];
    }
    for (k = 1; k < nb_fix_coupon2; k++)
    {
        fix_coupon_date2[k] = (fix_pay_dates2[k] - lToday) * YEARS_IN_DAY;
        DF_fix2[k]          = swp_f_df(lToday, fix_pay_dates2[k], yc_name);
        fix_cvg2[k]         = fix_cvgs2[k - 1];
    }
    //-----------------------------------------------------------------
    //----------------End Of Build Swap Schedule-----------------------
    //-----------------------------------------------------------------

    output_tmp = (double*)calloc(5, sizeof(double));
    err        = compute_CMScorrel_GL(
        Tf,
        Tp,
        DF_fix1,
        fix_coupon_date1,
        fix_cvg1,
        nb_fix_coupon1,
        DF_float1,
        DF_Tf1,
        float_coupon_date1,
        float_cvg1,
        nb_float_coupon1,
        spread1,
        Tf,
        Tp,
        DF_fix2,
        fix_coupon_date2,
        fix_cvg2,
        nb_fix_coupon2,
        DF_float2,
        DF_Tf2,
        float_coupon_date2,
        float_cvg2,
        nb_float_coupon2,
        spread2,
        sigma_1,
        sigma_2,
        sigma_date,
        sigma_n,
        lambda1,
        lambda2,
        rho,
        sigma_1,
        sigma_2,
        sigma_date,
        sigma_n,
        lambda1,
        lambda2,
        rho,
        nb_integration,
        bound,
        &output_tmp);
    (*output)[0] = output_tmp[0];
    (*output)[1] = output_tmp[1];
    (*output)[2] = output_tmp[2];
    (*output)[3] = output_tmp[3];
    (*output)[4] = output_tmp[4];

FREE_RETURN:

    if (float_coupon_date1)
        free(float_coupon_date1);
    if (fix_coupon_date1)
        free(fix_coupon_date1);
    if (DF_float1)
        free(DF_float1);
    if (DF_fix1)
        free(DF_fix1);
    if (float_cvg1)
        free(float_cvg1);
    if (spread1)
        free(spread1);
    if (fix_cvg1)
        free(fix_cvg1);

    if (float_fixing_dates1)
        free(float_fixing_dates1);
    if (float_start_dates1)
        free(float_start_dates1);
    if (float_end_dates1)
        free(float_end_dates1);
    if (float_pay_dates1)
        free(float_pay_dates1);

    if (float_cvgs1)
        free(float_cvgs1);
    if (float_spreads1)
        free(float_spreads1);
    if (float_pay_times1)
        free_dvector(float_pay_times1, 0, float_nb_pay_dates1 - 1);

    if (fix_start_dates1)
        free(fix_start_dates1);
    if (fix_end_dates1)
        free(fix_end_dates1);
    if (fix_pay_dates1)
        free(fix_pay_dates1);
    if (fix_pay_times1)
        free_dvector(fix_pay_times1, 0, fix_nb_pay_dates1 - 1);
    if (fix_cvgs1)
        free(fix_cvgs1);

    if (float_coupon_date2)
        free(float_coupon_date2);
    if (fix_coupon_date2)
        free(fix_coupon_date2);
    if (DF_float2)
        free(DF_float2);
    if (DF_fix2)
        free(DF_fix2);
    if (float_cvg2)
        free(float_cvg2);
    if (spread2)
        free(spread2);
    if (fix_cvg2)
        free(fix_cvg2);

    if (float_fixing_dates2)
        free(float_fixing_dates2);
    if (float_start_dates2)
        free(float_start_dates2);
    if (float_end_dates2)
        free(float_end_dates2);
    if (float_pay_dates2)
        free(float_pay_dates2);

    if (float_cvgs2)
        free(float_cvgs2);
    if (float_spreads2)
        free(float_spreads2);
    if (float_pay_times2)
        free_dvector(float_pay_times2, 0, float_nb_pay_dates2 - 1);

    if (fix_start_dates2)
        free(fix_start_dates2);
    if (fix_end_dates2)
        free(fix_end_dates2);
    if (fix_pay_dates2)
        free(fix_pay_dates2);
    if (fix_pay_times2)
        free_dvector(fix_pay_times2, 0, fix_nb_pay_dates2 - 1);
    if (fix_cvgs2)
        free(fix_cvgs2);

    if (output_tmp)
        free(output_tmp);

    return err;
}
