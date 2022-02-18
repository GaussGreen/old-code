/*****************************************************************************************/
/*                                                                                       */
/* FUNCTION     	: */
/*														                                 */
/*                                                                                       */
/* PURPOSE      	:
 */
/*                                                                                       */
/* DESCRIPTION  	:
 */
/*																						 */
/*															                             */
/*                                                                                       */
/*		                                                                                 */
/* PARAMETERS                                                                            */
/*	INPUT	    : */
/*              : */
/*              : */
/*              : */
/*				: */
/*												                                         */
/*												                                         */
/* RETURNS      :                               	      	                             */
/*                                                                                       */
/*****************************************************************************************/

/* -------------------------------------------------------------------------- */

#include "math.h"
#include "opfnctns.h"
#include "srt_h_all.h"

/* -------------------------------------------------------------------------- */

void MomentsComputing(
    int      iDimension,
    double   dMaturity,
    double*  dForwards,
    double*  dVols,
    double** dCorrels,
    double*  dMoments)
{
    double dMoment1, dMoment2, dMoment3;
    double dMoment2Container;
    double dMoment3Container1, dMoment3Container2;
    long   i, j, k;

    dMoment1 = 0.;
    dMoment2 = 0.;
    dMoment3 = 0.;

    for (i = 0; i < iDimension; i++)
    {
        dMoment1 = dMoment1 + dForwards[i];

        dMoment2Container = 0.;
        for (j = 0; j < iDimension; j++)
        {
            dMoment2Container = dMoment2Container + dForwards[j] * exp(dVols[i] * dVols[j] *
                                                                       dCorrels[i][j] * dMaturity);
        }
        dMoment2 = dMoment2 + dForwards[i] * dMoment2Container;

        dMoment3Container2 = 0.;
        for (j = 0; j < iDimension; j++)
        {
            dMoment3Container1 = 0.;
            for (k = 0; k < iDimension; k++)
            {
                dMoment3Container1 =
                    dMoment3Container1 + dForwards[k] * exp((dVols[i] * dVols[j] * dCorrels[i][j] +
                                                             dVols[i] * dVols[k] * dCorrels[i][k] +
                                                             dVols[k] * dVols[j] * dCorrels[k][j]) *
                                                            dMaturity);
            }
            dMoment3Container2 = dMoment3Container2 + dForwards[j] * dMoment3Container1;
        }
        dMoment3 = dMoment3 + dForwards[i] * dMoment3Container2;
    }

    dMoments[0] = dMoment1;
    dMoments[1] = dMoment2 - dMoment1 * dMoment1;
    dMoments[2] = dMoment3 - 3. * dMoment1 * dMoment2 + 2. * dMoment1 * dMoment1 * dMoment1;
}

void LogNormalParametersComputing(double* dMoments, double* dParameters)
{
    double dBigX, dU, dZ, dX, dY;
    double dLambda1, dLambda2;

    dBigX = (dMoments[2] + sqrt(pow(dMoments[2], 2.) + 4. * pow(dMoments[1], 3.))) / 2.;
    dU    = pow(dBigX, 1. / 3.);
    dZ    = dU - dMoments[1] / dU;
    dX    = dMoments[1] / dZ;
    dY    = dZ + dX;

    dLambda1 = dX;
    dLambda2 = dX * dY;

    dParameters[0] = 2. * log(dLambda1) - 0.5 * log(dLambda2);
    dParameters[1] = log(dLambda2) - 2. * log(dLambda1);
    dParameters[2] = dMoments[0] - dLambda1;
}

Err srt_f_BasketCapPrice(
    int              iNbRates,
    double           dStartDate,
    double           dStrike,
    double*          dForwards,
    double*          dVols,
    double**         dCorrels,
    SrtGreekType     SrtGreek,
    SrtCallPutType   SrtCallPut,
    SrtDiffusionType SrtVolType,
    double*          dAnswer)
{
    Err     err = NULL;
    int     i   = 0;
    double  dForward;
    double  dSigma;
    double* dMoments;
    double* dParameters;

    dMoments    = (double*)calloc(3, sizeof(double));
    dParameters = (double*)calloc(3, sizeof(double));

    /* Greek Switch */
    switch (SrtGreek)
    {
    case PREMIUM:

        MomentsComputing(iNbRates, dStartDate, dForwards, dVols, dCorrels, dMoments);

        LogNormalParametersComputing(dMoments, dParameters);

        dStrike -= dParameters[2];
        dForward = exp(dParameters[0] + 0.5 * dParameters[1]);
        dSigma   = sqrt(dParameters[1]);

        if (SrtVolType == SRT_NORMAL)
            *dAnswer = 0.;
        else
            *dAnswer = srt_f_optblksch(dForward, dStrike, dSigma, 1.0, 1.0, SrtCallPut, SrtGreek);
        break;

    default:

        *dAnswer = 0.;
        break;
    }
    return err;
}

/*	Low level lognormal */
static double inflation_multiput_ll(
    double    T,
    int       nidx,
    double*   ai,
    double*   fyi,
    double*   ki,
    double*   si,
    double    fx,
    double    kx,
    double    sx,
    double**  rho, /*	index 0..nidx-1: xi, index nidx: x */
    int       npth,
    double*** g) /*	Sobol cube */
{
    int     i, j, k;
    double *ex = NULL, **cov = NULL, **h = NULL, **chol = NULL, *yti = NULL;
    double  sum;
    double  payoff, prem = -99999.99;

    /*	Calculate cov */
    cov = dmatrix(0, nidx - 1, 0, nidx - 1);
    if (!cov)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i < nidx; i++)
    {
        cov[i][i] = si[i] * si[i] * T;
        for (j = 0; j < i; j++)
        {
            cov[i][j] = cov[j][i] = si[i] * si[j] * rho[i][j] * T;
        }
    }

    /*	Calculate ex under Qx */
    ex = dvector(0, nidx - 1);
    if (!exp)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i < nidx; i++)
    {
        ex[i] = (-0.50 * si[i] * si[i] + rho[i][nidx] * sx * si[i]) * T;
    }

    /*	Allocate memory for Sobol points */
    h    = dmatrix(0, npth - 1, 0, nidx - 1);
    chol = dmatrix(0, nidx - 1, 0, nidx - 1);

    if (!h || !chol)
    {
        goto FREE_RETURN;
    }

    /*	Calculate cholesky */
    nr_choldc(nidx, cov, chol);

    /*	Incorporate covariance */
    for (i = 0; i < npth; i++)
    {
        for (j = 0; j < nidx; j++)
        {
            sum = 0.0;
            for (k = 0; k <= j; k++)
            {
                sum += chol[j][k] * g[i][0][k];
            }
            h[i][j] = sum;
        }
    }

    yti = dvector(0, nidx - 1);
    if (!yti)
    {
        goto FREE_RETURN;
    }

    prem = 0.0;
    /*	Calculate price */
    for (i = 0; i < npth; i++)
    {
        /*	Reconstruct indices */
        for (j = 0; j < nidx; j++)
        {
            yti[j] = fyi[j] * exp(ex[j] + h[i][j]);
        }

        /*	Compute payoff */
        payoff = 0.0;

        for (j = 0; j < nidx; j++)
        {
            payoff += ai[j] * min(yti[j], ki[j]);
        }

        payoff = max(0, payoff - kx);

        /*	Update price */
        prem += payoff / npth;
    }
    prem *= fx;

FREE_RETURN:

    if (ex)
        free_dvector(ex, 0, nidx - 1);
    if (cov)
        free_dmatrix(cov, 0, nidx - 1, 0, nidx - 1);
    if (h)
        free_dmatrix(h, 0, npth - 1, 0, nidx - 1);
    if (chol)
        free_dmatrix(chol, 0, nidx - 1, 0, nidx - 1);
    if (yti)
        free_dvector(yti, 0, nidx - 1);

    return prem;
}

/*	Low level shifted */
static double inflation_multiput_sl_ll(
    double    T,
    int       nidx,
    double*   ai,
    double*   fyi,
    double*   ki,
    double*   shifti,
    double*   voli,
    double    fx,
    double    kx,
    double    sx,
    double**  rho, /*	index 0..nidx-1: xi, index nidx: x */
    int       npth,
    double*** g) /*	Sobol cube */
{
    int     i;
    double *fyi_ = NULL, *ki_ = NULL;
    double  result;

    fyi_ = dvector(0, nidx - 1);
    ki_  = dvector(0, nidx - 1);
    if (!fyi_ || !ki_)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i < nidx; i++)
    {
        fyi_[i] = fyi[i] + shifti[i];
        ki_[i]  = ki[i] + shifti[i];
        kx += ai[i] * shifti[i];
    }

    result = inflation_multiput_ll(T, nidx, ai, fyi_, ki_, voli, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (fyi_)
        free_dvector(fyi_, 0, nidx - 1);
    if (ki_)
        free_dvector(ki_, 0, nidx - 1);

    return result;
}

static void atm_beta_2_sl(double fwd, double atm, double beta, double* shift, double* vol)
{
    *shift = (1.0 - beta) / beta * fwd;
    *vol   = beta * atm;
}

static void sl_2_atm_beta(double fwd, double shift, double vol, double* atm, double* beta)
{
    *beta = fwd / (fwd + shift);
    *atm  = vol / (*beta);
}

/*	Price X * max ( 0 , sum [ai * min (Yi , Ki)] - Kx ) */
double inflation_multiput(
    double   T,
    int      nidx,
    double*  ai,
    double*  fyi,
    double*  ki,
    double*  si,
    double   fx,
    double   kx,
    double   sx,
    double** rho, /* index 0..nidx-1: xi, index nidx: x */
    int      npth)
{
    double*** g = NULL;
    double    result;

    g = f3tensor(0, npth - 1, 0, 0, 0, nidx - 1);
    if (!g)
    {
        goto FREE_RETURN;
    }

    /*	Initialise Sobol sequences */
    sobol_init(0, npth - 1, 0, 0, 0, nidx - 1);
    sobol_cube(g, 0, npth - 1, 0, 0, 0, nidx - 1);

    result = inflation_multiput_ll(T, nidx, ai, fyi, ki, si, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (g)
        free_f3tensor(g, 0, npth - 1, 0, 0, 0, nidx);
    sobol_free();
    return result;
}

/*	Price X * max ( 0 , sum [ai * min (Yi , Ki)] - Kx ) in a shifted-log model */
double inflation_multiput_sl(
    double   T,
    int      nidx,
    double*  ai,
    double*  fyi,
    double*  ki,
    double*  shifti,
    double*  voli,
    double   fx,
    double   kx,
    double   sx,
    double** rho, /*	index 0..nidx-1: xi, index nidx: x */
    int      npth,
    int      atmbumptype,  /*	0: no bump, 1: add, 2: mult */
    int      betabumptype, /*	0: no bump, 1: add, 2: mult */
    double*  atmbump,
    double*  betabump)
{
    double*** g = NULL;
    double    atm, beta, result;
    int       i;

    g = f3tensor(0, npth - 1, 0, 0, 0, nidx - 1);
    if (!g)
    {
        goto FREE_RETURN;
    }

    /*	Initialise Sobol sequences */
    sobol_init(0, npth - 1, 0, 0, 0, nidx - 1);
    sobol_cube(g, 0, npth - 1, 0, 0, 0, nidx - 1);

    /*	Bump */
    if (atmbumptype || betabumptype)
    {
        for (i = 0; i < nidx; i++)
        {
            sl_2_atm_beta(fyi[i], shifti[i], voli[i], &atm, &beta);

            switch (atmbumptype)
            {
            case 0:
            default:
                break;

            case 1:
                atm += atmbump[i];
                break;

            case 2:
                atm *= 1.0 + atmbump[i];
                break;
            }

            switch (betabumptype)
            {
            case 0:
            default:
                break;

            case 1:
                beta += betabump[i];
                break;

            case 2:
                beta *= 1.0 + betabump[i];
                break;
            }

            atm_beta_2_sl(fyi[i], atm, beta, shifti + i, voli + i);
        }
    }

    result = inflation_multiput_sl_ll(T, nidx, ai, fyi, ki, shifti, voli, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (g)
        free_f3tensor(g, 0, npth - 1, 0, 0, 0, nidx - 1);
    sobol_free();
    return result;
}

/*	Find shifted-log parameters by calibration to 2 strikes */
/*	The model is:
. F = L - shift or L = F + shift where L is a lognormal martingale with a certain vol
. If the vol of L is negative in the case of a sub-normal skew, then
        F = -L - shift or L = -F - shift where L is a lognormal martingale with vol abs (vol)
        and led by a reversed Brownian Motion from F
. The implied beta is also calculated */

static double price_given_vol_and_shift(
    double t,     /*	Maturity in years */
    double fwd,   /*	Fwd */
    double k,     /*	Strike */
    double vol,   /*	Vol of L */
    double shift) /*	Shift */
{
    if (vol > 0.0)
    {
        /*	Supernormal case */
        return srt_f_optblksch(fwd + shift, k + shift, vol, t, 1.0, SRT_CALL, PREMIUM);
    }
    else
    {
        /*	Subnormal case */
        return srt_f_optblksch(-fwd - shift, -k - shift, -vol, t, 1.0, SRT_PUT, PREMIUM);
    }
}

static double price_given_shift(
    double  t,    /*	Maturity in years */
    double  fwd,  /*	Fwd */
    double  k1,   /*	Strike 1 */
    double  p1,   /*	Price 1: to deduce vol */
    double  k2,   /*	Strike 2: price 2 will be returned */
    double* vol,  /*	Vol: will be changed on exit */
    double  shift) /*	Shift */
{
    if (*vol > 0.0)
    {
        /*	Supernormal case */
        srt_f_optimpvol(p1, fwd + shift, k1 + shift, t, 1.0, SRT_CALL, SRT_LOGNORMAL, vol);
    }
    else
    {
        /*	Subnormal case */
        srt_f_optimpvol(p1, -fwd - shift, -k1 - shift, t, 1.0, SRT_PUT, SRT_LOGNORMAL, vol);
        *vol *= -1;
    }

    return price_given_vol_and_shift(t, fwd, k2, *vol, shift);
}

static void implement_shift_limits(
    double  fwd,
    double  k1,
    double  k2,
    double  p1,
    double  p2,
    int     shift_sign,
    int     vol_sign,
    double  s1,
    double  s2,
    double* s)
{
    double slim, lim;

    if (*s * shift_sign < 0)
    {
        if (shift_sign > 0)
        {
            slim = min(s1, s2);
        }
        else
        {
            slim = max(s1, s2);
        }
        *s = 0.5 * slim;
    }

    if (vol_sign > 0)
    {
        lim = -min(fwd, min(k1, k2));
        if (*s < lim)
        {
            if (fabs(s1) > 1.0e-08)
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = min(s1, s2);
                }
                else
                {
                    slim = s1;
                }
            }
            else
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = s2;
                }
                else
                {
                    slim = lim + 0.50 * fabs(lim);
                }
            }
            *s = 0.5 * (lim + slim);
        }
    }
    else
    {
        lim = -max(fwd, max(k1, k2));
        if (*s > lim)
        {
            if (fabs(s1) > 1.0e-08)
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = max(s1, s2);
                }
                else
                {
                    slim = s1;
                }
            }
            else
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = s2;
                }
                else
                {
                    slim = lim - 0.50 * fabs(lim);
                }
            }
            *s = 0.5 * (lim + slim);
        }
    }

    if (vol_sign > 0)
    {
        lim = max(p1, p2) - fwd;
        if (*s < lim)
        {
            if (fabs(s1) > 1.0e-08)
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = min(s1, s2);
                }
                else
                {
                    slim = s1;
                }
            }
            else
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = s2;
                }
                else
                {
                    slim = lim + 0.50 * fabs(lim);
                }
            }
            *s = 0.5 * (lim + slim);
        }
    }
    else
    {
        lim = -max(p1 + k1, p2 + k2);
        if (*s > lim)
        {
            if (fabs(s1) > 1.0e-08)
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = max(s1, s2);
                }
                else
                {
                    slim = s1;
                }
            }
            else
            {
                if (fabs(s2) > 1.0e-08)
                {
                    slim = s2;
                }
                else
                {
                    slim = lim - 0.50 * fabs(lim);
                }
            }
            *s = 0.5 * (lim + slim);
        }
    }
}

Err find_shifted_log_params(
    double t,   /*	Maturity in years */
    double fwd, /*	Forward */
    double k1,  /*	1st strike */
    double v1,  /*	BS implied vol for 1st strike */
    double k2,  /*	2nd strike */
    double v2,  /*	BS implied vol for 2nd strike */
    /*	Results */
    double* shift,
    double* vol,
    double* beta)
{
    double p1, p2, nv1, nv2, l1, l2, l, s1, s2, s, x1, x2, x, temp1, temp2;
    double vega1, vega2, tempk, tempv, tempp;
    int    shift_sign, vol_sign;
    int    it;

    p1 = srt_f_optblksch(fwd, k1, v1, t, 1.0, SRT_CALL, PREMIUM);
    p2 = srt_f_optblksch(fwd, k2, v2, t, 1.0, SRT_CALL, PREMIUM);

    vega1 = srt_f_optblksch(fwd, k1, v1, t, 1.0, SRT_CALL, VEGA);
    vega2 = srt_f_optblksch(fwd, k2, v2, t, 1.0, SRT_CALL, VEGA);

    if (vega2 > vega1)
    {
        tempk = k1;
        tempv = v1;
        tempp = p1;

        k1 = k2;
        v1 = v2;
        p1 = p2;

        k2 = tempk;
        v2 = tempv;
        p2 = tempp;
    }

    srt_f_optimpvol(p1, fwd, k1, t, 1.0, SRT_CALL, SRT_NORMAL, &nv1);
    srt_f_optimpvol(p2, fwd, k2, t, 1.0, SRT_CALL, SRT_NORMAL, &nv2);

    if (fabs(v2 - v1) < 1.0e-05)
    {
        *shift = 0;
        *vol   = v1;
        *beta  = 1;
        return NULL;
    }

    if (fabs(nv2 - nv1) < 1.0e-06)
    {
        *shift = 1.e06;
        *vol   = nv1 * 1.e-06;
        *beta  = 0;
        return NULL;
    }

    if (fabs(k2 - k1) < 1.0e-06)
    {
        return "Error: different vols for same strike";
    }

    l1 = l2 = l = 1.0;
    shift_sign = vol_sign = 1;
    if ((v2 - v1) / (k2 - k1) > 0)
    {
        shift_sign = -1;
    }
    if ((nv2 - nv1) / (k2 - k1) < 0)
    {
        shift_sign = -1;
        vol_sign   = -1;
    }
    l1 *= vol_sign;
    l2 *= vol_sign;
    l *= vol_sign;

    temp1 = 0.5 * (fwd + k1);
    temp2 = 0.5 * (fwd + k2);
    s     = temp1 * temp2 * (v1 - v2) / (temp2 * v2 - temp1 * v1);
    implement_shift_limits(fwd, k1, k2, p1, p2, shift_sign, vol_sign, 0.0, 0.0, &s);

    s1 = 0.8 * s;
    s2 = 2 * s - s1;
    implement_shift_limits(fwd, k1, k2, p1, p2, shift_sign, vol_sign, 0.0, s, &s1);
    implement_shift_limits(fwd, k1, k2, p1, p2, shift_sign, vol_sign, s, s1, &s2);
    if (fabs(s2 - s1) < 1.0e-04)
        s2 += 1.0e-04;

    x1 = price_given_shift(t, fwd, k1, p1, k2, &l1, s1);
    x2 = price_given_shift(t, fwd, k1, p1, k2, &l2, s2);

    it = 1;
    s  = s1 + (p2 - x1) / (x2 - x1) * (s2 - s1);
    implement_shift_limits(fwd, k1, k2, p1, p2, shift_sign, vol_sign, s1, s2, &s);

    x = price_given_shift(t, fwd, k1, p1, k2, &l, s);

    while (fabs(x - p2) > 1.0e-08 && it < 15)
    {
        if (fabs(x1 - p2) < fabs(x2 - p2))
        {
            s2 = s;
            x2 = x;
        }
        else
        {
            s1 = s;
            x1 = x;
        }

        s = s1 + (p2 - x1) / (x2 - x1) * (s2 - s1);
        implement_shift_limits(fwd, k1, k2, p1, p2, shift_sign, vol_sign, s1, s2, &s);

        x = price_given_shift(t, fwd, k1, p1, k2, &l, s);

        it++;
    }

    *shift = s;
    *vol   = l;
    *beta  = fwd / (fwd + s);

    return NULL;
}

char* i_stellar_cpn(
    /*	The product */
    int     prod_type,     /*	0: Stellar, 1: I-Stellar */
    double  fix_mat_years, /*	(Fixing date - max (Ref date, Today)) / 365 */
    double  vol_mat_years, /*	(Fixing date - Today) / 365 */
    int     nidx,          /*	Number of eqd indices in play */
    double* idx_weights,   /*	Weights of the indices in the payoff e.g. 0.45 */
    double* idx_ref,       /*	The base values for indices
                                                           historical fixing if Ref Date <= Today
                                                           or interpolated from forward curve if Ref date
                              > Today */
    double* idx_str,       /*	Stellar strikes in % of ref values, e.g. 1.065 */
    int*    idx_qto,       /*	Wether eqd index is quantoed */
    double  cpi_ref,       /*	The base value for CPI
                                                           historical fixing if Ref Date <= Today
                                                           or interpolated from forward curve if Ref date >
                              Today */
    double global_str,     /*	Global strike for the payoff */
    /*	The IR market */
    double pay_df, /*	DF to payment date times coverage */
    /*	The EQD market */
    double* idx_spots,                      /*	Spots of indices */
    double* idx_fwd,                        /*	Forwards of indices for fixing date interpolated from
                                                                            forward curves */
    char** idx_vol_name,                    /*	Names of the vol curve of the indices */
    double (*get_eqd_vol)(                  /*	GetVol function for eqd indices */
                          double mat_years, /*	maturity in years */
                          double strike_spot, /*	strike in % of spot, e.g. 1.065 */
                          char*  name),        /*	name of the vol curve */
    /*	Smile parameters
            0: ATMS
            1: ATMF
            2: ISTR
            3: Custom */
    int     strike_1_type,
    int     strike_2_type,
    double* strike_1_custom,
    double* strike_2_custom,
    /*	The CPI market */
    double cpi_fwd, /*	Forward CPI for fixing date interpolated from forward curve */
    double cpi_vol, /*	CPI vol for fixing date interpolated from vol curve
                                                    WE NEED CUMULATIVE CPI VOL */
    /*	The Fx market */
    double* fx_vol,    /*	For each eqd index, the volatility of the corresponding Fx
                                                       interpolated from the implied BS volatility
                          curve of Fx    for the fixing date */
    double* fx_correl, /*	Correl between the eqd index and the fx
                                                       EXPRESSED IN EQD CCY / TRADE(cpi) CCY
                                                       interpolated from correlation curve by tenor
                                                       (not by fixed maturity) */
    /*	The basket correl as interpolated for the fixing date from the relevant
            term structures BY TENOR */
    double** rho, /*	0..nidx-1: eqd indices, nidx: CPI */
    /*	Parameters */
    int integ_points, /*	10 to 15, default 12 */
    /*	MAD: 0 everywhere */
    int     atmf_bump_type, /*	0: none, 1: add, 2: mult */
    int     beta_bump_type,
    double* atmf_bump,
    double* beta_bump,
    /*	Result */
    double*         pv,
    int             store_info,
    i_stellar_info* info)
{
    int     i, j;
    double  str1, str2, vol1, vol2;
    double *shift = NULL, *vol = NULL, *beta = NULL;
    Err     err = NULL;

    if (!prod_type)
    {
        cpi_ref = cpi_fwd = 1.0;
        cpi_vol           = 0.0;
    }

    if (store_info)
    {
        info->mat     = fix_mat_years;
        info->df      = pay_df;
        info->num_eqd = nidx;
        for (i = 0; i < nidx; i++)
        {
            info->fwd[i] = idx_fwd[i];
            for (j = 0; j < nidx; j++)
            {
                info->rho[i][j] = rho[i][j];
            }
        }
    }

    shift = dvector(0, nidx - 1);
    vol   = dvector(0, nidx - 1);
    beta  = dvector(0, nidx - 1);

    if (!shift || !vol || !beta)
    {
        err = "Allocation error in i_stellar_cpn";
        goto FREE_RETURN;
    }

    /*	For each eqd index, calibrate shifted-log dynamics */
    for (i = 0; i < nidx; i++)
    {
        switch (strike_1_type)
        {
        case 0:
        default:
            str1 = idx_spots[i];
            break;
        case 1:
            str1 = idx_fwd[i];
            break;
        case 2:
            str1 = idx_str[i] * idx_ref[i];
            break;
        case 3:
            str1 = strike_1_custom[i] * idx_spots[i];
            break;
        }

        switch (strike_2_type)
        {
        case 0:
        default:
            str2 = idx_spots[i];
            break;
        case 1:
            str2 = idx_fwd[i];
            break;
        case 2:
            str2 = idx_str[i] * idx_ref[i];
            break;
        case 3:
            str2 = strike_2_custom[i] * idx_spots[i];
            break;
        }

        vol1 = get_eqd_vol(vol_mat_years, str1 / idx_spots[i], idx_vol_name[i]);
        vol2 = get_eqd_vol(vol_mat_years, str2 / idx_spots[i], idx_vol_name[i]);

        if (store_info)
        {
            info->strike1[i] = str1;
            info->strike2[i] = str2;
            info->vol1[i]    = vol1;
            info->vol2[i]    = vol2;
        }

        err = find_shifted_log_params(
            fix_mat_years, idx_fwd[i], str1, vol1, str2, vol2, shift + i, vol + i, beta + i);

        if (err)
            goto FREE_RETURN;

        /*	For each eqd index, adjust for quanto */
        if (idx_qto[i])
        {
            idx_fwd[i] =
                (idx_fwd[i] + shift[i]) * exp(fx_correl[i] * fx_vol[i] * vol[i] * fix_mat_years) -
                shift[i];
        }

        if (store_info)
        {
            info->qto[i]   = idx_fwd[i];
            info->shift[i] = shift[i];
        }

        /*	For each eqd index, normalise by base values */
        idx_fwd[i] /= idx_ref[i];
        shift[i] /= idx_ref[i];

        /*	Store info */
        if (store_info)
        {
            info->ref[i]        = idx_ref[i];
            info->spot[i]       = idx_spots[i];
            info->vol[i]        = vol[i];
            info->beta[i]       = beta[i];
            info->used_fwd[i]   = idx_fwd[i];
            info->used_shift[i] = shift[i];
        }

    } /*	End of loop on eqd indices */

    /*	Normalise CPI by base values */
    if (store_info)
    {
        info->cpi_fwd = cpi_fwd;
        info->cpi_vol = cpi_vol;
    }

    cpi_fwd /= cpi_ref;

    if (store_info)
    {
        info->used_cpi_fwd = cpi_fwd;
    }

    /*	Value coupon */
    *pv = inflation_multiput_sl(
        fix_mat_years,
        nidx,
        idx_weights,
        idx_fwd,
        idx_str,
        shift,
        vol,
        cpi_fwd,
        global_str,
        cpi_vol,
        rho,
        (int)(pow(2.0, (double)integ_points) + 1.0e-08),
        atmf_bump_type,
        beta_bump_type,
        atmf_bump,
        beta_bump);

    /*	Multiply by DF */
    *pv *= pay_df;

    /*	Free memory */
FREE_RETURN:
    if (shift)
        free_dvector(shift, 0, nidx - 1);
    if (vol)
        free_dvector(vol, 0, nidx - 1);
    if (beta)
        free_dvector(beta, 0, nidx - 1);

    /*	Return */
    return err;
}

/*	CMT Stellar */
/*	XXXXXXXXXXX	*/

/*	Low level lognormal */
static double cmt_stellar_ll(
    double    T,
    int       nidx,
    double*   ai,
    double*   fyi,
    double*   ki,
    double*   si,
    double    ax,
    double    fx,
    double    kx,
    double    sx,
    double**  rho, /*	index 0..nidx-1: xi, index nidx: x */
    int       npth,
    double*** g) /*	Sobol cube */
{
    int     i, j, k;
    double *ex = NULL, **cov = NULL, **h = NULL, **chol = NULL, *yti = NULL;
    double  xt;
    double  sum;
    double  payoff, prem = -99999.99;

    /*	Calculate cov */
    cov = dmatrix(0, nidx, 0, nidx);
    if (!cov)
    {
        goto FREE_RETURN;
    }

    cov[nidx][nidx] = sx * sx * T;
    for (i = 0; i < nidx; i++)
    {
        cov[i][i]    = si[i] * si[i] * T;
        cov[i][nidx] = cov[nidx][i] = si[i] * sx * rho[i][nidx] * T;
        for (j = 0; j < i; j++)
        {
            cov[i][j] = cov[j][i] = si[i] * si[j] * rho[i][j] * T;
        }
    }

    /*	Calculate ex */
    ex = dvector(0, nidx);
    if (!exp)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i <= nidx; i++)
    {
        ex[i] = -0.50 * cov[i][i];
    }

    /*	Allocate memory for Sobol points */
    h    = dmatrix(0, npth - 1, 0, nidx);
    chol = dmatrix(0, nidx, 0, nidx);

    if (!h || !chol)
    {
        goto FREE_RETURN;
    }

    /*	Calculate cholesky */
    nr_choldc(nidx + 1, cov, chol);

    /*	Incorporate covariance */
    for (i = 0; i < npth; i++)
    {
        for (j = 0; j <= nidx; j++)
        {
            sum = 0.0;
            for (k = 0; k <= j; k++)
            {
                sum += chol[j][k] * g[i][0][k];
            }
            h[i][j] = sum;
        }
    }

    yti = dvector(0, nidx - 1);
    if (!yti)
    {
        goto FREE_RETURN;
    }

    prem = 0.0;
    /*	Calculate price */
    for (i = 0; i < npth; i++)
    {
        /*	Reconstruct indices */
        for (j = 0; j < nidx; j++)
        {
            yti[j] = fyi[j] * exp(ex[j] + h[i][j]);
        }
        xt = fx * exp(ex[nidx] + h[i][nidx]);

        /*	Compute payoff */
        payoff = 0.0;

        for (j = 0; j < nidx; j++)
        {
            payoff += ai[j] * min(yti[j], ki[j]);
        }

        payoff = max(ax * xt, payoff - kx);

        /*	Update price */
        prem += payoff / npth;
    }

FREE_RETURN:

    if (ex)
        free_dvector(ex, 0, nidx);
    if (cov)
        free_dmatrix(cov, 0, nidx, 0, nidx);
    if (h)
        free_dmatrix(h, 0, npth - 1, 0, nidx);
    if (chol)
        free_dmatrix(chol, 0, nidx, 0, nidx);
    if (yti)
        free_dvector(yti, 0, nidx - 1);

    return prem;
}

/*	Low level shifted */
static double cmt_stellar_sl_ll(
    double    T,
    int       nidx,
    double*   ai,
    double*   fyi,
    double*   ki,
    double*   shifti,
    double*   voli,
    double    ax,
    double    fx,
    double    kx,
    double    sx,
    double**  rho, /*	index 0..nidx-1: xi, index nidx: x */
    int       npth,
    double*** g) /*	Sobol cube */
{
    int     i;
    double *fyi_ = NULL, *ki_ = NULL;
    double  result;

    fyi_ = dvector(0, nidx - 1);
    ki_  = dvector(0, nidx - 1);
    if (!fyi_ || !ki_)
    {
        goto FREE_RETURN;
    }

    for (i = 0; i < nidx; i++)
    {
        fyi_[i] = fyi[i] + shifti[i];
        ki_[i]  = ki[i] + shifti[i];
        kx += ai[i] * shifti[i];
    }

    result = cmt_stellar_ll(T, nidx, ai, fyi_, ki_, voli, ax, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (fyi_)
        free_dvector(fyi_, 0, nidx - 1);
    if (ki_)
        free_dvector(ki_, 0, nidx - 1);

    return result;
}

/*	Price max ( ax * X , sum [ai * min (Yi , Ki)] - Kx ) */
double cmt_stellar(
    double   T,
    int      nidx,
    double*  ai,
    double*  fyi,
    double*  ki,
    double*  si,
    double   ax,
    double   fx,
    double   kx,
    double   sx,
    double** rho, /* index 0..nidx-1: xi, index nidx: x */
    int      npth)
{
    double*** g = NULL;
    double    result;

    g = f3tensor(0, npth - 1, 0, 0, 0, nidx);
    if (!g)
    {
        goto FREE_RETURN;
    }

    /*	Initialise Sobol sequences */
    sobol_init(0, npth - 1, 0, 0, 0, nidx);
    sobol_cube(g, 0, npth - 1, 0, 0, 0, nidx);

    result = cmt_stellar_ll(T, nidx, ai, fyi, ki, si, ax, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (g)
        free_f3tensor(g, 0, npth - 1, 0, 0, 0, nidx);
    sobol_free();
    return result;
}

/*	Price max ( ax * X , sum [ai * min (Yi , Ki)] - Kx ) in a shifted-log model */
double cmt_stellar_sl(
    double   T,
    int      nidx,
    double*  ai,
    double*  fyi,
    double*  ki,
    double*  shifti,
    double*  voli,
    double   ax,
    double   fx,
    double   kx,
    double   sx,
    double** rho, /*	index 0..nidx-1: xi, index nidx: x */
    int      npth,
    int      atmbumptype,  /*	0: no bump, 1: add, 2: mult */
    int      betabumptype, /*	0: no bump, 1: add, 2: mult */
    double*  atmbump,
    double*  betabump)
{
    double*** g = NULL;
    double    atm, beta, result;
    int       i;

    g = f3tensor(0, npth - 1, 0, 0, 0, nidx);
    if (!g)
    {
        goto FREE_RETURN;
    }

    /*	Initialise Sobol sequences */
    sobol_init(0, npth - 1, 0, 0, 0, nidx);
    sobol_cube(g, 0, npth - 1, 0, 0, 0, nidx);

    /*	Bump */
    if (atmbumptype || betabumptype)
    {
        for (i = 0; i < nidx; i++)
        {
            sl_2_atm_beta(fyi[i], shifti[i], voli[i], &atm, &beta);

            switch (atmbumptype)
            {
            case 0:
            default:
                break;

            case 1:
                atm += atmbump[i];
                break;

            case 2:
                atm *= 1.0 + atmbump[i];
                break;
            }

            switch (betabumptype)
            {
            case 0:
            default:
                break;

            case 1:
                beta += betabump[i];
                break;

            case 2:
                beta *= 1.0 + betabump[i];
                break;
            }

            atm_beta_2_sl(fyi[i], atm, beta, shifti + i, voli + i);
        }
    }

    result = cmt_stellar_sl_ll(T, nidx, ai, fyi, ki, shifti, voli, ax, fx, kx, sx, rho, npth, g);

FREE_RETURN:

    if (g)
        free_f3tensor(g, 0, npth - 1, 0, 0, 0, nidx);
    sobol_free();
    return result;
}

char* cmt_stellar_cpn(
    /*	The product */
    int     prod_type,     /*	0: Stellar, 1: CMT-Stellar */
    double  fix_mat_years, /*	(Fixing date - max (Ref date, Today)) / 365 */
    double  vol_mat_years, /*	(Fixing date - Today) / 365 */
    int     nidx,          /*	Number of eqd indices in play */
    double* idx_weights,   /*	Weights of the indices in the payoff e.g. 0.45 */
    double* idx_ref,       /*	The base values for indices
                                                           historical fixing if Ref Date <= Today
                                                           or interpolated from forward curve if Ref date
                              > Today */
    double* idx_str,       /*	Stellar strikes in % of ref values, e.g. 1.065 */
    int*    idx_qto,       /*	Wether eqd index is quantoed */
    double  global_str,    /*	Global strike for the payoff */
    /*	The IR market */
    double pay_df, /*	DF to payment date times coverage */
    /*	The EQD market */
    double* idx_spots,                      /*	Spots of indices */
    double* idx_fwd,                        /*	Forwards of indices for fixing date interpolated from
                                                                            forward curves */
    char** idx_vol_name,                    /*	Names of the vol curve of the indices */
    double (*get_eqd_vol)(                  /*	GetVol function for eqd indices */
                          double mat_years, /*	maturity in years */
                          double strike_spot, /*	strike in % of spot, e.g. 1.065 */
                          char*  name),        /*	name of the vol curve */
    /*	Smile parameters
            0: ATMS
            1: ATMF
            2: ISTR
            3: Custom */
    int     strike_1_type,
    int     strike_2_type,
    double* strike_1_custom,
    double* strike_2_custom,
    /*	The CMT market */
    double cmt_weight,
    double cmt_fwd, /*	Forward CMT for fixing date interpolated from forward curve */
    double cmt_vol, /*	CMT vol for fixing date interpolated from vol curve
                                                    WE NEED CUMULATIVE CMT VOL */
    /*	The Fx market */
    double* fx_vol,    /*	For each eqd index, the volatility of the corresponding Fx
                                                       interpolated from the implied BS volatility
                          curve of Fx    for the fixing date */
    double* fx_correl, /*	Correl between the eqd index and the fx
                                                       EXPRESSED IN EQD CCY / TRADE(cpi) CCY
                                                       interpolated from correlation curve by tenor
                                                       (not by fixed maturity) */
    /*	The basket correl as interpolated for the fixing date from the relevant
            term structures BY TENOR */
    double** rho, /*	0..nidx-1: eqd indices, nidx: CMT */
    /*	Parameters */
    int integ_points, /*	10 to 15, default 12 */
    /*	MAD: 0 everywhere */
    int     atmf_bump_type, /*	0: none, 1: add, 2: mult */
    int     beta_bump_type,
    double* atmf_bump,
    double* beta_bump,
    /*	Result */
    double*           pv,
    int               store_info,
    cmt_stellar_info* info)
{
    int     i, j;
    double  str1, str2, vol1, vol2;
    double *shift = NULL, *vol = NULL, *beta = NULL;
    Err     err = NULL;

    if (!prod_type)
    {
        cmt_weight = 0.0;
        cmt_fwd    = 0.0;
        cmt_vol    = 0.0;
    }

    if (store_info)
    {
        info->mat     = fix_mat_years;
        info->df      = pay_df;
        info->num_eqd = nidx;
        for (i = 0; i < nidx; i++)
        {
            info->fwd[i] = idx_fwd[i];
            for (j = 0; j < nidx; j++)
            {
                info->rho[i][j] = rho[i][j];
            }
        }
    }

    shift = dvector(0, nidx - 1);
    vol   = dvector(0, nidx - 1);
    beta  = dvector(0, nidx - 1);

    if (!shift || !vol || !beta)
    {
        err = "Allocation error in i_stellar_cpn";
        goto FREE_RETURN;
    }

    /*	For each eqd index, calibrate shifted-log dynamics */
    for (i = 0; i < nidx; i++)
    {
        switch (strike_1_type)
        {
        case 0:
        default:
            str1 = idx_spots[i];
            break;
        case 1:
            str1 = idx_fwd[i];
            break;
        case 2:
            str1 = idx_str[i] * idx_ref[i];
            break;
        case 3:
            str1 = strike_1_custom[i] * idx_spots[i];
            break;
        }

        switch (strike_2_type)
        {
        case 0:
        default:
            str2 = idx_spots[i];
            break;
        case 1:
            str2 = idx_fwd[i];
            break;
        case 2:
            str2 = idx_str[i] * idx_ref[i];
            break;
        case 3:
            str2 = strike_2_custom[i] * idx_spots[i];
            break;
        }

        vol1 = get_eqd_vol(vol_mat_years, str1 / idx_spots[i], idx_vol_name[i]);
        vol2 = get_eqd_vol(vol_mat_years, str2 / idx_spots[i], idx_vol_name[i]);

        if (store_info)
        {
            info->strike1[i] = str1;
            info->strike2[i] = str2;
            info->vol1[i]    = vol1;
            info->vol2[i]    = vol2;
        }

        err = find_shifted_log_params(
            fix_mat_years, idx_fwd[i], str1, vol1, str2, vol2, shift + i, vol + i, beta + i);

        if (err)
            goto FREE_RETURN;

        /*	For each eqd index, adjust for quanto */
        if (idx_qto[i])
        {
            idx_fwd[i] =
                (idx_fwd[i] + shift[i]) * exp(fx_correl[i] * fx_vol[i] * vol[i] * fix_mat_years) -
                shift[i];
        }

        if (store_info)
        {
            info->qto[i]   = idx_fwd[i];
            info->shift[i] = shift[i];
        }

        /*	For each eqd index, normalise by base values */
        idx_fwd[i] /= idx_ref[i];
        shift[i] /= idx_ref[i];

        /*	Store info */
        if (store_info)
        {
            info->ref[i]        = idx_ref[i];
            info->spot[i]       = idx_spots[i];
            info->vol[i]        = vol[i];
            info->beta[i]       = beta[i];
            info->used_fwd[i]   = idx_fwd[i];
            info->used_shift[i] = shift[i];
        }

    } /*	End of loop on eqd indices */

    if (store_info)
    {
        info->cmt_fwd = cmt_fwd;
        info->cmt_vol = cmt_vol;
    }

    /*	Value coupon */
    *pv = cmt_stellar_sl(
        fix_mat_years,
        nidx,
        idx_weights,
        idx_fwd,
        idx_str,
        shift,
        vol,
        cmt_weight,
        cmt_fwd,
        global_str,
        cmt_vol,
        rho,
        (int)(pow(2.0, (double)integ_points) + 1.0e-08),
        atmf_bump_type,
        beta_bump_type,
        atmf_bump,
        beta_bump);

    /*	Multiply by DF */
    *pv *= pay_df;

    /*	Free memory */
FREE_RETURN:
    if (shift)
        free_dvector(shift, 0, nidx - 1);
    if (vol)
        free_dvector(vol, 0, nidx - 1);
    if (beta)
        free_dvector(beta, 0, nidx - 1);

    /*	Return */
    return err;
}