/* ===================================================================================
   FILENAME:      BVM_Perturbation.c

   PURPOSE:       Computes Option Prices via Perturbation Method in BVM
   =================================================================================== */

#pragma warning(disable : 4786)  // Disable long name warnings

#include "Fx3FUtils.h"
#include "math.h"
#include "num_h_GaussianIntegral.h"
#include "num_h_allhdr.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_vol.h"
#include "utconst.h"

#define PI 3.141592653589793238462643383279502884197169399375105820974944592308

static double BlackDerivative(double F, double sigma, double T, double K, int iF, int iSigma)
{
    double o1, o2, o3, o4, o5, o6, o7, o8, o9;
    double res;

    switch (iF)
    {
    case 4:
        switch (iSigma)
        {
        case 0:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 6.);
            o4  = pow(sigma, 2.);
            o5  = log(o2);
            o6  = pow(o5, 2.);
            o7  = o4 * T;
            res = (0.0249338925250895 * K * exp((-0.5 * o6) / (o4 * T) - 0.125 * o4 * T) *
                   (-8. * o4 * o6 * (6. + o7) * T + 16. * pow(o5, 4.) +
                    o3 * (-4. + o7) * pow(T, 3.)) *
                   pow(T, -3. / 2.) * sqrt(o2)) /
                  o3;
            break;
        }
        break;

    case 3:
        switch (iSigma)
        {
        case 0:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            res = -0.199471140200716 * K * (2. * o4 + 3. * o3 * T) *
                  exp(-0.125 * o3 * T - (0.5 * pow(o4, 2.)) / (o3 * T)) * pow(F, -3.) *
                  pow(sigma, -3.) * pow(T, -3. / 2.) * sqrt(o2);
            break;
        case 1:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            o5  = o3 * T;
            res = 0.0498677850501791 * K * (2. * o4 + o5) *
                  (-4. * o4 * (o4 + o5) + 3. * o3 * (4. + o5) * T) *
                  exp(-0.125 * o3 * T - (0.5 * pow(o4, 2.)) / (o3 * T)) * pow(F, -3.) *
                  pow(sigma, -6.) * pow(T, -5. / 2.) * sqrt(o2);
            break;
        }
        break;

    case 2:
        switch (iSigma)
        {
        case 0:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, GAMMA);
            break;
        case 1:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            o5  = pow(o4, 2.);
            res = -0.0997355701003582 * K * (-4. * o5 + o3 * T * (4. + o3 * T)) *
                  exp((-0.5 * o5) / (o3 * T) - 0.125 * o3 * T) * pow(F, -2.) * pow(sigma, -4.) *
                  pow(T, -3. / 2.) * sqrt(o2);
            break;
        case 2:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            o5  = pow(o4, 2.);
            o6  = o3 * T;
            res = 0.0249338925250895 * K * exp((-0.5 * o5) / (o3 * T) - 0.125 * o3 * T) *
                  pow(F, -2.) * pow(sigma, -7.) *
                  (-8. * o3 * o5 * (10. + o6) * T + 16. * pow(o4, 4.) +
                   (32. + o3 * (4. + o6) * T) * pow(sigma, 4.) * pow(T, 2.)) *
                  pow(T, -5. / 2.) * sqrt(o2);
            break;
        }
        break;

    case 1:
        switch (iSigma)
        {
        case 0:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, DELTA);
            break;
        case 1:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, VANNA);
            break;
        case 2:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            o5  = pow(o4, 2.);
            res = (-0.04986778505017909 * exp((-0.5 * o5) / (o3 * T) - 0.125 * o3 * T) *
                   pow(sigma, -5.) *
                   (-4. * o3 * o5 * T - 2. * o3 * o4 * T * (8. + o3 * T) + 8. * pow(o4, 3.) +
                    pow(sigma, 6.) * pow(T, 3.)) *
                   pow(T, -3. / 2.)) /
                  sqrt(o2);
            break;
        case 3:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 8.);
            o4  = pow(sigma, 2.);
            o5  = log(o2);
            o6  = pow(o5, 2.);
            o7  = o4 * T;
            o8  = pow(sigma, 4.);
            o9  = pow(T, 2.);
            res = (0.012466946262544772 * exp((-0.5 * o6) / (o4 * T) - 0.125 * o4 * T) *
                   (2. * o5 *
                        (-(o8 * o9 * (96. + o4 * (12. + o7) * T)) +
                         4. * o5 *
                             (-((6. + o7) * o8 * o9) +
                              2. * o5 * (-2. * o6 + o4 * o5 * T + o4 * (14. + o7) * T))) +
                    o3 * (-4. + o7) * pow(T, 4.)) *
                   pow(T, -5. / 2.)) /
                  (o3 * sqrt(o2));
            break;
        }
        break;

    case 0:
        switch (iSigma)
        {
        case 0:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, PREMIUM);
            break;
        case 1:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, VEGA);
            break;
        case 2:
            res = srt_f_optblksch(F, K, sigma, T, 1.0, SRT_CALL, VOLGA);
            break;
        case 3:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 6.);
            o4  = pow(sigma, 2.);
            o5  = log(o2);
            o6  = pow(o5, 2.);
            o7  = o4 * T;
            res = (0.0249338925250895 * K * exp((-0.5 * o6) / (o4 * T) - 0.125 * o4 * T) *
                   (-8. * o4 * o6 * (6. + o7) * T + 16. * pow(o5, 4.) +
                    o3 * (-4. + o7) * pow(T, 3.)) *
                   pow(T, -3. / 2.) * sqrt(o2)) /
                  o3;
            break;
        case 4:
            o1  = 1 / K;
            o2  = F * o1;
            o3  = pow(sigma, 2.);
            o4  = log(o2);
            o5  = pow(o4, 2.);
            o6  = o3 * T;
            res = -0.00623347313127239 * K * exp((-0.5 * o5) / (o3 * T) - 0.125 * o3 * T) *
                  pow(sigma, -9.) *
                  (48. * o3 * (12. + o6) * T * pow(o4, 4.) - 64. * pow(o4, 6.) -
                   12. * o5 * (64. + o3 * (8. + o6) * T) * pow(sigma, 4.) * pow(T, 2.) +
                   (-12. + o6) * pow(sigma, 10.) * pow(T, 5.)) *
                  pow(T, -5. / 2.) * sqrt(o2);
            break;
        }
        break;
    }

    return res;
}

double bvm_local_vol(double x, double a, double b, double c, int type)
{
    if (type == 0)
    {
        return (1 - exp(-x / a));
    }
    else
    {
        return (1 - exp(-x / a)) / x;
    }
}

static double First_Order_Integrand_function(
    double date,
    double FwdInit,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    L_or_N,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    double* f1,
    double* f2)
{
    double res;
    double d1, d2;
    double Gamma, Vega, Volga, Vanna;
    double tau, sqrt_tau;

    tau      = Maturity - date;
    sqrt_tau = sqrt(tau);

    d1 = (log(Fwd / Strike) + 0.5 * Vol * Vol * tau) / (Vol * sqrt_tau);
    d2 = d1 - Vol * sqrt_tau;

    if (L_or_N == 1)
    {
        Gamma = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, GAMMA);

        Vega = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VEGA);

        Volga = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VOLGA);

        Vanna = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VANNA);

        *f1 =
            Vol * Vol * Fwd * Fwd *
            (local_vol(Fwd, gamma, 0, 0, L_or_N) / local_vol(FwdInit, gamma, 0, 0, L_or_N) - 1.0) *
            Gamma;

        *f2 = Lambda * (VolInfinity - Vol) * Vega;

        *f2 += 0.5 * alpha * alpha * Vol * Vol * Volga;

        *f2 += rho * alpha * Vol * Fwd * Vol * Vanna;

        res = *f1 + *f2;
    }
    else
    {
        Gamma = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, GAMMA);

        Vega = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VEGA);

        Volga = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VOLGA);

        Vanna = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VANNA);

        *f1 =
            Vol * Vol *
            (local_vol(Fwd, gamma, 0, 0, L_or_N) / local_vol(FwdInit, gamma, 0, 0, L_or_N) - 1.0) *
            Gamma;

        *f2 = Lambda * (VolInfinity - Vol) * Vega;

        *f2 += 0.5 * alpha * alpha * Vol * Vol * Volga;

        *f2 += rho * alpha * Vol * Vol * Vanna;

        res = *f1 + *f2;
    }

    //	Gamma = exp(-0.5 * d1 * d1) / sqrt(2 * PI) / (Fwd * Vol * sqrt_tau);
    //	Vega = Fwd * sqrt_tau * exp(-0.5 * d1 * d1) / sqrt(2*PI);
    //	Volga = Fwd * sqrt_tau * exp(-0.5 * d1 * d1) / sqrt(2*PI)
    //			* d1 * d2 / Vol;
    //	Vanna = - d2 * exp(-0.5 * d1 * d1) / sqrt(2*PI) / Vol;

    return res;
}

static double Expectation_In_Fwd(
    double today,
    double date,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    log_or_norm,
    double (*local_vol)(double x, double a, double b, double c, int lgo_or_norm),
    int     NQuadrature,
    double* w,
    double* x,
    double* f1,
    double* f2)
{
    int    i;
    double res;
    double Forward;
    double tempf1, tempf2;

    res = 0.0;
    *f1 = 0.0;
    *f2 = 0.0;
    for (i = 1; i <= NQuadrature; ++i)
    {
        if (log_or_norm == 1)
        {
            Forward =
                Fwd * exp(-0.5 * (date - today) * Vol * Vol) * exp(Vol * sqrt(date - today) * x[i]);
        }
        else
        {
            Forward = Fwd + Vol * sqrt(date - today) * x[i];
        }

        res += w[i] * First_Order_Integrand_function(
                          date,
                          Fwd,
                          Forward,
                          Vol,
                          Strike,
                          Maturity,
                          gamma,
                          alpha,
                          rho,
                          VolInfinity,
                          Lambda,
                          log_or_norm,
                          local_vol,
                          &tempf1,
                          &tempf2);

        *f1 += w[i] * tempf1;
        *f2 += w[i] * tempf2;
    }

    return res;
}

static double First_Order_Integral_In_Time(
    double today,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    log_or_norm,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg,
    double* f1,
    double* f2)
{
    int    i;
    double res;
    double tempf1, tempf2;

    res = 0.0;
    *f1 = 0.0;
    *f2 = 0.0;
    for (i = 1; i <= NQuad_Leg; ++i)
    {
        if (today < x_leg[i])
        {
            res += w_leg[i] * Expectation_In_Fwd(
                                  today,
                                  x_leg[i],
                                  Fwd,
                                  Vol,
                                  Strike,
                                  Maturity,
                                  gamma,
                                  alpha,
                                  rho,
                                  VolInfinity,
                                  Lambda,
                                  log_or_norm,
                                  local_vol,
                                  NQuad_Her,
                                  w_her,
                                  x_her,
                                  &tempf1,
                                  &tempf2);

            *f1 += w_leg[i] * tempf1;
            *f2 += w_leg[i] * tempf2;
        }
    }

    return res;
}

static double Second_Order_Integrand_function(
    double date,
    double FwdInit,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    L_or_N,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg)
{
    double res;
    double Gamma, Vanna, Vega, Volga, Gamma2, Vanna12;
    double LocVolFwd, LocVolFwdInit;

    double dFwd, dVol;
    double f, f1, f2;
    double f_f_up, f1_f_up, f2_f_up;
    double f_f_down, f1_f_down, f2_f_down;
    double f_vol_up, f1_vol_up, f2_vol_up;
    double f_vol_down, f1_vol_down, f2_vol_down;
    double f_f_up_vol_up, f1_f_up_vol_up, f2_f_up_vol_up;

    LocVolFwd     = local_vol(Fwd, gamma, 0, 0, L_or_N);
    LocVolFwdInit = local_vol(FwdInit, gamma, 0, 0, L_or_N);

    dFwd = Fwd / 100.0;
    dVol = Vol / 100.0;

    f = First_Order_Integral_In_Time(
        date,
        Fwd,
        Vol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1,
        &f2);

    f_f_up = First_Order_Integral_In_Time(
        date,
        Fwd + dFwd,
        Vol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1_f_up,
        &f2_f_up);

    f_f_down = First_Order_Integral_In_Time(
        date,
        Fwd - dFwd,
        Vol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1_f_down,
        &f2_f_down);

    f_vol_up = First_Order_Integral_In_Time(
        date,
        Fwd,
        Vol + dVol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1_vol_up,
        &f2_vol_up);

    f_vol_down = First_Order_Integral_In_Time(
        date,
        Fwd,
        Vol - dVol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1_vol_down,
        &f2_vol_down);

    f_f_up_vol_up = First_Order_Integral_In_Time(
        date,
        Fwd + dFwd,
        Vol + dVol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        L_or_N,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1_f_up_vol_up,
        &f2_f_up_vol_up);

    if (L_or_N == 1)
    {
        Gamma = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, GAMMA);

        Vanna = srt_f_optblksch(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VANNA);

        Vega    = (f_vol_up - f) / dVol;
        Volga   = (f_vol_up - 2 * f + f_vol_down) / (dVol * dVol);
        Gamma2  = (f2_f_up - 2 * f2 + f2_f_down) / (dFwd * dFwd);
        Vanna12 = (f_f_up_vol_up - f_vol_up - f_f_up + f) / (dFwd * dVol);

        res = 0.5 * Vol * Vol * Fwd * Fwd * (LocVolFwd / LocVolFwdInit - 1.0) *
              (LocVolFwd / LocVolFwdInit - 1.0) * Gamma;

        res += rho * Fwd * (LocVolFwd / LocVolFwdInit - 1.0) * alpha * Vol * Vol * Vanna;

        res += Vol * Vol * Fwd * Fwd * (LocVolFwd / LocVolFwdInit - 1.0) * Gamma2;

        res += Lambda * (VolInfinity - Vol) * Vega;

        res += 0.5 * alpha * alpha * Vol * Vol * Volga;

        res += rho * alpha * Vol * Fwd * Vol * Vanna12;
    }
    else
    {
        Gamma = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, GAMMA);

        Vanna = srt_f_optblknrm(Fwd, Strike, Vol, Maturity - date, 1.0, SRT_CALL, VANNA);

        Vega    = (f_vol_up - f) / dVol;
        Volga   = (f_vol_up - 2 * f + f_vol_down) / (dVol * dVol);
        Gamma2  = (f2_f_up - 2 * f2 + f2_f_down) / (dFwd * dFwd);
        Vanna12 = (f_f_up_vol_up - f_vol_up - f_f_up + f) / (dFwd * dVol);

        res = 0.5 * Vol * Vol * (LocVolFwd / LocVolFwdInit - 1.0) *
              (LocVolFwd / LocVolFwdInit - 1.0) * Gamma;

        res += rho * (LocVolFwd / LocVolFwdInit - 1.0) * alpha * Vol * Vol * Vanna;

        res += Vol * Vol * (LocVolFwd / LocVolFwdInit - 1.0) * Gamma2;

        res += Lambda * (VolInfinity - Vol) * Vega;

        res += 0.5 * alpha * alpha * Vol * Vol * Volga;

        res += rho * alpha * Vol * Vol * Vanna12;
    }

    return res;
}

static double Second_Order_Expectation_In_Fwd(
    double today,
    double date,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    log_or_norm,
    double (*local_vol)(double x, double a, double b, double c, int lgo_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg)
{
    int    i;
    double res;
    double Forward;

    res = 0.0;
    for (i = 1; i <= NQuad_Her; ++i)
    {
        if (log_or_norm == 1)
        {
            Forward = Fwd * exp(-0.5 * (date - today) * Vol * Vol) *
                      exp(Vol * sqrt(date - today) * x_her[i]);
        }
        else
        {
            Forward = Fwd + Vol * sqrt(date - today) * x_her[i];
        }

        res += w_her[i] * Second_Order_Integrand_function(
                              date,
                              Fwd,
                              Forward,
                              Vol,
                              Strike,
                              Maturity,
                              gamma,
                              alpha,
                              rho,
                              VolInfinity,
                              Lambda,
                              log_or_norm,
                              local_vol,
                              NQuad_Her,
                              w_her,
                              x_her,
                              NQuad_Leg,
                              w_leg,
                              x_leg);
    }

    return res;
}

static double Second_Order_Integral_In_Time(
    double today,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    int    log_or_norm,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg)
{
    int    i;
    double res;

    res = 0.0;
    for (i = 1; i <= NQuad_Leg; ++i)
    {
        res += w_leg[i] * Second_Order_Expectation_In_Fwd(
                              today,
                              x_leg[i],
                              Fwd,
                              Vol,
                              Strike,
                              Maturity,
                              gamma,
                              alpha,
                              rho,
                              VolInfinity,
                              Lambda,
                              log_or_norm,
                              local_vol,
                              NQuad_Her,
                              w_her,
                              x_her,
                              NQuad_Leg,
                              w_leg,
                              x_leg);
    }

    return res;
}

double BVM_First_Order_Perturbation_Price(
    double today,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg,
    int     LOG_NORM)
{
    double res;
    double black, pertvalue, iv;
    double locvol;
    int    log_or_norm;
    double f1, f2;

    if (Strike >= Fwd)
    {
        log_or_norm = 1;
    }
    else
    {
        log_or_norm = 0;
    }

    if ((LOG_NORM == 0) || (LOG_NORM == 1))
    {
        log_or_norm = LOG_NORM;
    }
    locvol = local_vol(Fwd, gamma, 0, 0, log_or_norm);

    if (log_or_norm == 1)
    {
        black =
            srt_f_optblksch(Fwd, Strike, Vol * locvol, Maturity - today, 1.0, SRT_CALL, PREMIUM);
    }
    else
    {
        black =
            srt_f_optblknrm(Fwd, Strike, Vol * locvol, Maturity - today, 1.0, SRT_CALL, PREMIUM);
    }

    pertvalue = First_Order_Integral_In_Time(
        today,
        Fwd,
        Vol * locvol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        log_or_norm,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1,
        &f2);

    res = black + pertvalue;
    iv  = max(0, Fwd - Strike);

    if (res > iv)
    {
        return res;
    }
    else
    {
        return black;
    }
}

double BVM_Second_Order_Perturbation_Price(
    double today,
    double Fwd,
    double Vol,
    double Strike,
    double Maturity,
    double gamma,
    double alpha,
    double rho,
    double VolInfinity,
    double Lambda,
    double (*local_vol)(double x, double a, double b, double c, int log_or_norm),
    int     NQuad_Her,
    double* w_her,
    double* x_her,
    int     NQuad_Leg,
    double* w_leg,
    double* x_leg)
{
    double res;
    double black, firstpertvalue, secondpertvalue, iv;
    double locvol;
    int    log_or_norm;
    double f1, f2;

    if (Strike >= Fwd)
    {
        log_or_norm = 1;
    }
    else
    {
        log_or_norm = 0;
    }

    //	log_or_norm = 1;
    locvol = local_vol(Fwd, gamma, 0, 0, log_or_norm);

    if (log_or_norm == 1)
    {
        black =
            srt_f_optblksch(Fwd, Strike, Vol * locvol, Maturity - today, 1.0, SRT_CALL, PREMIUM);
    }
    else
    {
        black =
            srt_f_optblknrm(Fwd, Strike, Vol * locvol, Maturity - today, 1.0, SRT_CALL, PREMIUM);
    }

    firstpertvalue = First_Order_Integral_In_Time(
        today,
        Fwd,
        Vol * locvol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        log_or_norm,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg,
        &f1,
        &f2);

    secondpertvalue = Second_Order_Integral_In_Time(
        today,
        Fwd,
        Vol * locvol,
        Strike,
        Maturity,
        gamma,
        alpha,
        rho,
        VolInfinity,
        Lambda,
        log_or_norm,
        local_vol,
        NQuad_Her,
        w_her,
        x_her,
        NQuad_Leg,
        w_leg,
        x_leg);

    res = black + firstpertvalue + secondpertvalue;
    iv  = max(0, Fwd - Strike);

    if (res > iv)
    {
        return res;
    }
    else
    {
        return black;
    }
}

Err monte_carlo_sabrbvm(
    double** rndm_mat,
    double   forw,
    double   vovol,
    double   beta,
    double   rho,
    double   sigma,
    double   num_paths,
    double   maturity,
    double   num_steps,
    int      path_num,
    double*  Forward,
    int      modeltype,
    int      sampletype,
    int*     pathnogood)
{
    double g[2];
    double a1;
    double a2;
    double Eta;
    double delta_t;
    double sqrt_delta_t;
    double delta;
    double logEta, logEta0;
    int    step;
    Err    err = NULL;
    long   idum;
    double rand1, rand2;
    double statevar;  //,statevar0,statevar1;
    idum = -path_num;
    uniform(&idum);

    *pathnogood = 0;

    /* use a1 and a2 as symmetrical correlated values for the Brownians */

    a1 = sqrt((1 + rho) / 2);
    a2 = sqrt((1 - rho) / 2);

    step         = 1;
    delta_t      = maturity / num_steps;
    sqrt_delta_t = sqrt(delta_t);
    Eta          = sigma;
    logEta       = log(sigma);
    logEta0      = log(sigma);
    Forward[0]   = forw;
    Forward[1]   = forw;
    if (beta != 1)
    {
        statevar = pow(forw, 2 - 2 * beta) / pow(sigma * (1 - beta), 2);
        delta    = (1 - 2 * beta) / (1 - beta);
    }

    while (step < num_steps)

    {
        if (sampletype == 0) /*Randsam*/
        {
            rand1 = inv_cumnorm_fast(uniform(&idum));
            rand2 = inv_cumnorm_fast(uniform(&idum));
            g[0]  = a1 * rand1 + a2 * rand2;
            g[1]  = a1 * rand1 - a2 * rand2;
        }

        else
        {
            g[0] = a1 * rndm_mat[1][step] + a2 * rndm_mat[2][step];
            g[1] = a1 * rndm_mat[1][step] - a2 * rndm_mat[2][step];
        }

        switch (modeltype)

        {
        case 0:
            /* Euler Discretization*/

            Forward[0] += Eta * (1 - exp(-beta * Forward[0])) * g[0] * sqrt_delta_t;
            Forward[1] += Eta * (1 - exp(-beta * Forward[1])) * g[0] * sqrt_delta_t;

            //			Forward[0]+=Eta*pow(Forward[0],beta)*g[0]*sqrt_delta_t;
            //			Forward[1]+=sigma*pow(Forward[1],beta)*g[0]*sqrt_delta_t;
            break;

        case 1:

            /* Mihlstein discretization */

            Forward[0] +=
                (1 - exp(-beta * Forward[0])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;
            Forward[1] +=
                (1 - exp(-beta * Forward[1])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;
            //+jumps[step]/sigma+(U2*lambda2-U1*lambda1)*delta_t);
            break;
        }

        //			logEta = log(Eta);
        //			logEta = logEta + (-0.5 * vovol*vovol + lambda * ( logEta0 - logEta ) )*
        //delta_t
        //							+ vovol*g[1]*sqrt_delta_t;
        //			Eta = exp(logEta);

        Eta = Eta * exp(vovol * g[1] * sqrt_delta_t - 0.5 * vovol * vovol * delta_t);

        if (Forward[0] < 0)
        {
            /*				*pathnogood = 1;*/
            Forward[0] = 0;
        }

        step++;
    }

    /* free memory*/

    return NULL;
}

Err monte_carlo_sabrbvm3(
    double  forw,
    double  vovol,
    double  beta,
    double  rho,
    double  sigma,
    double  lambda,
    double  num_paths,
    double  maturity,
    double  num_steps,
    int     path_num,
    double* Forward,
    int     modeltype,
    int     sampletype,
    int*    pathnogood)
{
    double g[2];
    double a1;
    double a2;
    double Eta;
    double delta_t;
    double sqrt_delta_t;
    double delta;
    double logEta, logEta0;
    int    step;
    Err    err = NULL;
    long   idum;
    double rand1, rand2;
    double statevar;  //,statevar0,statevar1;
    idum = -path_num;
    uniform(&idum);

    *pathnogood = 0;

    /* use a1 and a2 as symmetrical correlated values for the Brownians */

    a1 = sqrt((1 + rho) / 2);
    a2 = sqrt((1 - rho) / 2);

    step         = 1;
    delta_t      = maturity / num_steps;
    sqrt_delta_t = sqrt(delta_t);
    Eta          = sigma;
    logEta       = log(sigma);
    logEta0      = log(sigma);
    Forward[0]   = forw;
    Forward[1]   = forw;
    if (beta != 1)
    {
        statevar = pow(forw, 2 - 2 * beta) / pow(sigma * (1 - beta), 2);
        delta    = (1 - 2 * beta) / (1 - beta);
    }

    while (step < num_steps)

    {
        rand1 = inv_cumnorm_fast(uniform(&idum));
        rand2 = inv_cumnorm_fast(uniform(&idum));
        g[0]  = a1 * rand1 + a2 * rand2;
        g[1]  = a1 * rand1 - a2 * rand2;

        switch (modeltype)

        {
        case 0:
            /* Euler Discretization*/

            Forward[0] += Eta * (1 - exp(-beta * Forward[0])) * g[0] * sqrt_delta_t;
            Forward[1] += Eta * (1 - exp(-beta * Forward[1])) * g[0] * sqrt_delta_t;

            break;

        case 1:

            /* Mihlstein discretization */

            Forward[0] +=
                (1 - exp(-beta * Forward[0])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;
            Forward[1] +=
                (1 - exp(-beta * Forward[1])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;

            break;
        }

        logEta = log(Eta);
        logEta = logEta + (-0.5 * vovol * vovol + lambda * (logEta0 - logEta)) * delta_t +
                 vovol * g[1] * sqrt_delta_t;
        Eta = exp(logEta);

        //			Eta = Eta*exp(vovol*g[1]*sqrt_delta_t-0.5*vovol*vovol*delta_t);

        if (Forward[0] < 0)
        {
            Forward[0] = 0;
        }

        step++;
    }

    /* free memory*/

    return NULL;
}

Err monte_carlo_sabrbvm2(
    double  forw,
    double  vovol,
    double  beta,
    double  rho,
    double  sigma,
    double  num_paths,
    double  maturity,
    double  num_steps,
    int     path_num,
    double* Forward,
    int     modeltype,
    int     sampletype,
    int*    pathnogood)
{
    double g[2];
    double a1;
    double a2;
    double Eta;
    double delta_t;
    double sqrt_delta_t;
    double delta;
    int    step;
    Err    err = NULL;
    long   idum;
    double rand1, rand2;
    double statevar;  //,statevar0,statevar1;
    idum = -path_num;
    uniform(&idum);

    *pathnogood = 0;

    /* use a1 and a2 as symmetrical correlated values for the Brownians */

    a1 = sqrt((1 + rho) / 2);
    a2 = sqrt((1 - rho) / 2);

    step         = 1;
    delta_t      = maturity / num_steps;
    sqrt_delta_t = sqrt(delta_t);
    Eta          = sigma;
    Forward[0]   = forw;
    Forward[1]   = forw;
    if (beta != 1)
    {
        statevar = pow(forw, 2 - 2 * beta) / pow(sigma * (1 - beta), 2);
        delta    = (1 - 2 * beta) / (1 - beta);
    }

    while (step < num_steps)

    {
        rand1 = inv_cumnorm_fast(uniform(&idum));
        rand2 = inv_cumnorm_fast(uniform(&idum));
        g[0]  = a1 * rand1 + a2 * rand2;
        g[1]  = a1 * rand1 - a2 * rand2;
        switch (modeltype)

        {
        case 0:
            /* Euler Discretization*/

            Forward[0] += Eta * (1 - exp(-beta * Forward[0])) * g[0] * sqrt_delta_t;
            Forward[1] += Eta * (1 - exp(-beta * Forward[1])) * g[0] * sqrt_delta_t;

            //			Forward[0]+=Eta*pow(Forward[0],beta)*g[0]*sqrt_delta_t;
            //			Forward[1]+=sigma*pow(Forward[1],beta)*g[0]*sqrt_delta_t;
            break;

        case 1:

            /* Mihlstein discretization */

            Forward[0] +=
                (1 - exp(-beta * Forward[0])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;
            Forward[1] +=
                (1 - exp(-beta * Forward[1])) * (Eta * g[0] * sqrt_delta_t) +
                0.5 * Eta * Eta * beta * exp(-beta * Forward[0]) * (g[0] * g[0] - 1) * delta_t;
            //+jumps[step]/sigma+(U2*lambda2-U1*lambda1)*delta_t);
            break;
        }

        Eta = Eta * exp(vovol * g[1] * sqrt_delta_t - 0.5 * vovol * vovol * delta_t);

        if (Forward[0] < 0)
        {
            /*				*pathnogood = 1;*/
            Forward[0] = 0;
        }

        step++;
    }

    /* free memory*/

    return NULL;
}

Err sabrBVMmontecarlo2(
    double   forward,
    double*  strike,
    int      nb_strike,
    double   maturity,
    double   sigma_beta,
    double   alpha,
    double   beta,
    double   rho,
    double   lambda,
    int      npaths,
    int      nsteps,
    int      do_balsam,
    double** res)
{
    double   dt, sqdt;
    long     i, j;
    double   Z0, fwdZ, Z, sig;
    double   rho1, rho2, alpha2;
    double   w1, w2;
    double   lnsig, lnsig0;
    double   var;
    double   intsig;
    double** matrix = NULL;
    double*  matrixi;

    double impliedvol;

    long seed = -123456789;
    Err  err  = NULL;

    dt   = maturity / nsteps;
    sqdt = sqrt(dt);

    Z0     = forward;
    lnsig0 = log(sigma_beta);

    rho1 = alpha * rho * sqdt;
    rho2 = alpha * (1.0 - rho * rho) * sqdt;

    alpha2 = 0.5 * alpha * alpha * dt;

    for (j = 0; j < nb_strike; j++)
    {
        res[j][0] = 0.0;
        res[j][1] = 0.0;
    }

    /* fill the Brownian matrix */
    /* npaths has to be odd */
    npaths = 2 * ((long)(npaths / 2)) + 1;

    if (do_balsam)
    {
        matrix = dmatrix(0, npaths - 1, 0, nsteps - 1);
        if (!matrix)
        {
            return NULL;
        }

        err = balsam_generation(npaths, nsteps, matrix);
    }

    for (i = 0; i < npaths; i++)
    {
        Z     = Z0;
        lnsig = lnsig0;

        var    = 0.0;
        intsig = 0.0;

        if (do_balsam)
        {
            matrixi = matrix[i];
        }

        for (j = 0; j < nsteps; j++)
        {
            if (do_balsam)
            {
                w1 = matrixi[j];
            }
            else
            {
                w1 = gauss_sample(&seed);
            }
            w2  = gauss_sample(&seed);
            sig = exp(lnsig);
            Z += sig * (1 - exp(-beta * Z)) * sqdt * w1;
            Z = DMAX(0.0001, Z);
            lnsig += lambda * (lnsig0 - lnsig) - alpha2 * dt + rho1 * w1 + rho2 * w2;
        }

        fwdZ = Z;

        for (j = 0; j < nb_strike; j++)
        {
            res[j][0] += DMAX(0.0, fwdZ - strike[j]);
            res[j][1] += DMAX(0.0, fwdZ - strike[j]) * DMAX(0.0, fwdZ - strike[j]);
        }
    }

    for (j = 0; j < nb_strike; j++)
    {
        res[j][0] = res[j][0] / npaths;

        err = srt_f_optimpvol(
            res[j][0], forward, strike[j], maturity, 1, SRT_CALL, SRT_LOGNORMAL, &(impliedvol));
        if (err)
        {
            res[j][0] = 0;
            smessage("error in implied vol");
            err = NULL;
        }
        else
        {
            res[j][0] = impliedvol;
        }
        res[j][1] = sqrt((res[j][1] - res[j][0] * res[j][0]) / npaths);
    }

    if (matrix)
        free_dmatrix(matrix, 0, npaths - 1, 0, nsteps - 1);

    return NULL;
}

Err sabrBVMmontecarlo(
    double  forw,
    double  vovol,
    double  beta,
    double  rho,
    double  num_paths,
    double  num_steps,
    double  sigma,
    double  maturity,
    double* strike,
    int     num_strikes,
    int     modeltype,
    int     sampletype, /*0: randsam, 1: abs, 2:sobol, 3: SPECTRUNC*/
    double* impvol)
{
    double*** cube = NULL;
    int       cur_path;
    double*   Forward;
    double**  price = NULL;
    double*   times_at_steps;
    int       i;
    long      seed;
    double*   blkbeta;
    int       pathnogood;
    Err       err         = NULL;
    double    actnumpaths = num_paths;
    double*   correl      = NULL;
    double**  stdev       = NULL;
    double    correlation = 0;

    /*Memory allocation*/
    price  = dmatrix(1, 2, 1, num_strikes);
    stdev  = dmatrix(1, 2, 1, num_strikes);
    correl = dvector(1, num_strikes);
    //	cube = f3tensor(1,(long) num_paths,1,2,1,(long) num_steps);
    Forward        = dvector(1, 2);
    blkbeta        = dvector(1, num_strikes);
    times_at_steps = dvector(1, (long)num_steps);
    /* Choice of a method */

    for (i = 1; i < num_steps + 1; i++)
    {
        times_at_steps[i] = i * maturity / num_steps;
    }

    seed = -123456789;

    /* Spectrunc*/
    /*	if (sampletype == 3)
            {
                    if (err = SpecTruncCube(cube,times_at_steps,0.0001,1,(long)
       num_paths,1,2,1,(long) num_steps,&seed))
                    {
                            free_dmatrix(price,1,2,1,num_strikes);
                            free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
                            free_dvector(Forward,1,2);
                            free_dvector(blkbeta,1,num_strikes);
                            smessage("error in sobol cube");
                            return err;
                    }
            }
    */
    /*Sobol*/
    /*	if (sampletype==2)
            {
                    if (err = sobol_init(1,(long) num_paths,1,2,1,(long) num_steps))
                    {
                            free_dmatrix(price,1,2,1,num_strikes);
                            free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
                            free_dvector(Forward,1,2);
                            free_dvector(blkbeta,1,num_strikes);
                            smessage("error in sobol cube");
                            return err;
                    }



                    if (err = sobol_cube(cube,1,(long) num_paths,1,2,1,(long) num_steps))
                    {
                            free_dmatrix(price,1,2,1,num_strikes);
                            free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
                            free_dvector(Forward,1,2);
                            free_dvector(blkbeta,1,num_strikes);
                            smessage("error in sobol cube");
                            return err;
                    }

            }
    */
    /*Abs*/
    /*	else if (sampletype == 1)
            {
                    if (err = ABSCube(cube,1,(long) num_paths,1,2,1,(long) num_steps, &seed))
                    {
                            free_dmatrix(price,1,2,1,num_strikes);
                            free_dvector(Forward,1,2);
                            free_f3tensor(cube,1,(long) num_paths,1,2,1,(long) num_steps);
                            free_dvector(blkbeta,1,num_strikes);
                            smessage("error in abs cube");
                            return err;
                    }

            }
    */

    for (cur_path = 1; cur_path <= num_paths; cur_path++)
    {
        if (err = monte_carlo_sabrbvm(
                cube[cur_path],
                forw,
                vovol,
                beta,
                rho,
                sigma,
                num_paths,
                maturity,
                num_steps,
                cur_path + 1,
                Forward,
                modeltype,
                sampletype,
                &pathnogood))
        {
            free_dmatrix(price, 1, 2, 1, num_strikes);
            free_f3tensor(cube, 1, (long)num_paths, 1, 2, 1, (long)num_steps);
            free_dvector(Forward, 1, 2);
            free_dvector(blkbeta, 1, num_strikes);
            smessage("error in monte carlo");
            return err;
        }

        if (pathnogood == 0)
        {
            for (i = 1; i <= num_strikes; i++)
            {
                price[1][i] += DMAX(Forward[0] - strike[i], 0);
                stdev[1][i] += DMAX(Forward[0] - strike[i], 0) * DMAX(Forward[0] - strike[i], 0);
                correl[i] += DMAX(Forward[0] - strike[i], 0) * DMAX(Forward[1] - strike[i], 0);
                price[2][i] += DMAX(Forward[1] - strike[i], 0);
                stdev[2][i] += DMAX(Forward[1] - strike[i], 0) * DMAX(Forward[1] - strike[i], 0);
            }
        }
        else
        {
            actnumpaths--;
        }
    }

    for (i = 1; i <= num_strikes; i++)
    {
        correlation = (correl[i] - price[1][i] * price[2][i] / actnumpaths) /
                      sqrt(
                          (stdev[1][i] - price[1][i] * price[1][i] / actnumpaths) *
                          (stdev[2][i] - price[2][i] * price[2][i] / actnumpaths));

        price[1][i] = price[1][i] / actnumpaths;

        if (err = srt_f_optimpvol(
                price[1][i], forw, strike[i], maturity, 1, SRT_CALL, SRT_LOGNORMAL, &(impvol[i])))
        {
            impvol[i] = 0;
            smessage("error in implied vol");
            err = NULL;
        }
    }

    /*free memory*/
    free_dmatrix(price, 1, 2, 1, num_strikes);
    free_dmatrix(stdev, 1, 2, 1, num_strikes);
    free_dvector(correl, 1, num_strikes);
    free_f3tensor(cube, 1, (long)num_paths, 1, 2, 1, (long)num_steps);
    free_dvector(Forward, 1, 2);
    free_dvector(blkbeta, 1, num_strikes);
    free_dvector(times_at_steps, 1, (long)num_steps);
    return NULL;
}

Err sabrBVMmontecarlo3(
    double  forw,
    double  vovol,
    double  beta,
    double  rho,
    double  lambda,
    double  num_paths,
    double  num_steps,
    double  sigma,
    double  maturity,
    double* strike,
    int     num_strikes,
    int     modeltype,
    int     sampletype, /*0: randsam, 1: abs, 2:sobol, 3: SPECTRUNC*/
    double* impvol)
{
    double*** cube = NULL;
    int       cur_path;
    double*   Forward;
    double**  price = NULL;
    double*   times_at_steps;
    int       i;
    long      seed;
    double*   blkbeta;
    int       pathnogood;
    Err       err         = NULL;
    double    actnumpaths = num_paths;
    double*   correl      = NULL;
    double**  stdev       = NULL;
    double    correlation = 0;

    /*Memory allocation*/
    price          = dmatrix(1, 2, 1, num_strikes);
    stdev          = dmatrix(1, 2, 1, num_strikes);
    correl         = dvector(1, num_strikes);
    Forward        = dvector(1, 2);
    blkbeta        = dvector(1, num_strikes);
    times_at_steps = dvector(1, (long)num_steps);
    /* Choice of a method */

    for (i = 1; i < num_steps + 1; i++)
    {
        times_at_steps[i] = i * maturity / num_steps;
    }

    seed = -123456789;

    for (cur_path = 1; cur_path <= num_paths; cur_path++)
    {
        if (err = monte_carlo_sabrbvm3(
                forw,
                vovol,
                beta,
                rho,
                sigma,
                lambda,
                num_paths,
                maturity,
                num_steps,
                cur_path + 1,
                Forward,
                modeltype,
                sampletype,
                &pathnogood))
        {
            free_dmatrix(price, 1, 2, 1, num_strikes);
            free_dvector(Forward, 1, 2);
            free_dvector(blkbeta, 1, num_strikes);
            smessage("error in monte carlo");
            return err;
        }

        if (pathnogood == 0)
        {
            for (i = 1; i <= num_strikes; i++)
            {
                price[1][i] += DMAX(Forward[0] - strike[i], 0);
                stdev[1][i] += DMAX(Forward[0] - strike[i], 0) * DMAX(Forward[0] - strike[i], 0);
                correl[i] += DMAX(Forward[0] - strike[i], 0) * DMAX(Forward[1] - strike[i], 0);
                price[2][i] += DMAX(Forward[1] - strike[i], 0);
                stdev[2][i] += DMAX(Forward[1] - strike[i], 0) * DMAX(Forward[1] - strike[i], 0);
            }
        }
        else
        {
            actnumpaths--;
        }
    }

    for (i = 1; i <= num_strikes; i++)
    {
        correlation = (correl[i] - price[1][i] * price[2][i] / actnumpaths) /
                      sqrt(
                          (stdev[1][i] - price[1][i] * price[1][i] / actnumpaths) *
                          (stdev[2][i] - price[2][i] * price[2][i] / actnumpaths));

        price[1][i] = price[1][i] / actnumpaths;

        if (err = srt_f_optimpvol(
                price[1][i], forw, strike[i], maturity, 1, SRT_CALL, SRT_LOGNORMAL, &(impvol[i])))
        {
            impvol[i] = 0;
            smessage("error in implied vol");
            err = NULL;
        }
    }

    /*free memory*/
    free_dmatrix(price, 1, 2, 1, num_strikes);
    free_dmatrix(stdev, 1, 2, 1, num_strikes);
    free_dvector(correl, 1, num_strikes);
    free_dvector(Forward, 1, 2);
    free_dvector(blkbeta, 1, num_strikes);
    free_dvector(times_at_steps, 1, (long)num_steps);
    return NULL;
}

//--------------------------------------------------------------------------------
//-----------------------Pat's formula applied to BVM-----------------------------
//--------------------------------------------------------------------------------
double BVMLOGVol(
    double forward,
    double strike,
    double maturity,
    double betavol,
    double vovol,
    double rho,
    double eta)
{
    double favg;
    double zeta, xhat, zeta_xhat;
    double gamma1, gamma2;
    double volavg;
    double logvol;

    favg   = 0.5 * (forward + strike);
    volavg = (1 - exp(-eta * favg)) / eta;
    zeta   = (vovol / betavol) *
           (eta * (forward - strike) + log(1 - exp(-eta * forward)) - log(1 - exp(-eta * strike)));
    xhat = log((sqrt(1 - 2 * rho * zeta + zeta * zeta) - rho + zeta) / (1 - rho));
    if (fabs(zeta) < 1.0e-10)
    {
        zeta_xhat = 1.0;
    }
    else
    {
        zeta_xhat = zeta / xhat;
    }
    gamma1 = exp(-favg * eta) / (1 - exp(-favg * eta)) * eta;
    gamma2 = -gamma1 * eta;

    if (forward == strike)
    {
        logvol = betavol * zeta_xhat * volavg / forward;
    }
    else
    {
        logvol = betavol * zeta_xhat * log(forward / strike) /
                 (eta * (forward - strike) +
                  (log(1 - exp(-forward * eta)) - log(1 - exp(-strike * eta))));
    }

    logvol = logvol * (1 + ((2 * gamma2 - gamma1 * gamma1 + 1 / (favg * favg)) / 24.0 * betavol *
                                betavol * volavg * volavg +
                            0.25 * rho * vovol * betavol * gamma1 * volavg +
                            (2 - 3 * rho * rho) * vovol * vovol / 24.0) *
                               maturity);

    return logvol;
}

//--------------------------------------------------------------------------------
//-----------------------Pat's equivalent Normal Vol------------------------------
//--------------------------------------------------------------------------------
double SABRNormalVol(
    double forward,
    double strike,
    double maturity,
    double betavol,
    double vovol,
    double rho,
    double beta)
{
    double favg;
    double zeta, xhat, zeta_xhat;
    double gamma1, gamma2;
    double volavg;
    double normvol;

    favg   = 0.5 * (forward + strike);
    volavg = pow(favg, beta);
    zeta   = (vovol / betavol) * (forward - strike) / volavg;
    xhat   = log((sqrt(1 - 2 * rho * zeta + zeta * zeta) - rho + zeta) / (1 - rho));
    if (fabs(zeta) < 1.0e-10)
    {
        zeta_xhat = 1.0;
    }
    else
    {
        zeta_xhat = zeta / xhat;
    }
    gamma1 = beta / favg;
    gamma2 = beta * (beta - 1) / (favg * favg);

    if (forward == strike)
    {
        normvol = betavol * zeta_xhat * pow(forward, beta);
    }
    else
    {
        normvol = betavol * zeta_xhat * (1 - beta) * (forward - strike) /
                  (pow(forward, 1 - beta) - pow(strike, 1 - beta));
    }

    normvol = normvol *
              (1 + ((2 * gamma2 - gamma1 * gamma1) / 24.0 * betavol * betavol * volavg * volavg +
                    0.25 * rho * vovol * betavol * gamma1 * volavg +
                    (2 - 3 * rho * rho) * vovol * vovol / 24.0) *
                       maturity);

    return normvol;
}

double BVMApproxSimulUnderlying(double X, double XStar, double GStar, double Sigma, double Lambda)
{
    double res;

    //	res = exp(X - XStar) * GStar;
    res = log(1 + exp(X - XStar) * (exp(Lambda * GStar) - 1)) / Lambda;

    return res;
}

Err BVMApproxPrice(
    double  Xinit,
    double  XStar,
    double  GStar,
    double  Sigma,
    double  Lambda,
    double  Maturity,
    long    NHermite,
    long    NStrikes,
    double* Strikes,
    double* ImpliedVol,
    double* Fwd)
{
    Err     err = NULL;
    int     i, j;
    double* x_her = NULL;
    double* w_her = NULL;
    double  X, Underlying;

    x_her = (double*)calloc(NHermite + 1, sizeof(double));
    w_her = (double*)calloc(NHermite + 1, sizeof(double));
    if ((!x_her) || (!w_her))
    {
        err = "Allocation failed in BVMApproxPrice";
        goto FREE_RETURN;
    }

    err = HermiteStandard(x_her, w_her, NHermite);
    if (err)
    {
        goto FREE_RETURN;
    }

    *Fwd = 0;
    for (j = 0; j < NStrikes; ++j)
    {
        ImpliedVol[j] = 0;
    }

    for (i = 1; i <= NHermite; ++i)
    {
        X          = Xinit + Sigma * sqrt(Maturity) * x_her[i];
        Underlying = BVMApproxSimulUnderlying(X, XStar, GStar, Sigma, Lambda);
        *Fwd += w_her[i] * Underlying;
        for (j = 0; j < NStrikes; ++j)
        {
            ImpliedVol[j] += w_her[i] * DMAX(0, Underlying - Strikes[j]);
        }
    }

    //	for(j=0;j<NStrikes;++j)
    //	{
    //		err = srt_f_optimpvol(ImpliedVol[j], *Fwd, Strikes[j], Maturity, 1, SRT_CALL,
    //SRT_LOGNORMAL, &(impliedvol)); 		ImpliedVol[j] = impliedvol;
    //	}

FREE_RETURN:

    if (x_her)
    {
        free(x_her);
        x_her = NULL;
    }

    if (w_her)
    {
        free(w_her);
        w_her = NULL;
    }

    return err;
}

Err BVMApproxPriceFromFwd(
    double  Fwd,
    double  XStar,
    double  GStar,
    double  Sigma,
    double  Lambda,
    double  Maturity,
    long    NHermite,
    long    NStrikes,
    double* Strikes,
    double* ImpliedVol,
    double* FwdOutPut)
{
    Err    err = NULL;
    int    NIter;
    double Xinit, Xinit1, Xinit2;
    double Fwd1, Fwd2;

    Xinit  = XStar;
    Xinit1 = XStar;
    Xinit2 = XStar * 1.01;

    err = BVMApproxPrice(
        Xinit1, XStar, GStar, Sigma, Lambda, Maturity, NHermite, 0, NULL, NULL, &Fwd1);

    if (fabs(Fwd - Fwd1) < 1.0e-14)
    {
        Xinit = Xinit1;
        NIter = 20;
    }

    err = BVMApproxPrice(
        Xinit2, XStar, GStar, Sigma, Lambda, Maturity, NHermite, 0, NULL, NULL, &Fwd2);

    if (fabs(Fwd - Fwd2) < 1.0e-14)
    {
        Xinit  = Xinit2;
        Xinit1 = Xinit2;
        NIter  = 20;
    }

    if (fabs(Xinit1 - Xinit2) > 1.0e-20)
    {
        Xinit1 = Xinit1 + (Fwd - Fwd1) * (Xinit2 - Xinit1) / (Fwd2 - Fwd1);
        err    = BVMApproxPrice(
            Xinit1, XStar, GStar, Sigma, Lambda, Maturity, NHermite, 0, NULL, NULL, &Fwd1);
    }
    else
    {
        NIter = 20;
    }

    NIter = 0;
    while ((NIter < 20) && (fabs(Fwd - Fwd1) > 1.0e-7))
    {
        Xinit2 = Xinit1 * 1.01;

        err = BVMApproxPrice(
            Xinit2, XStar, GStar, Sigma, Lambda, Maturity, NHermite, 0, NULL, NULL, &Fwd2);

        if (fabs(Fwd - Fwd2) < 1.0e-14)
        {
            Xinit = Xinit2;
            NIter = 20;
        }

        if (fabs(Xinit1 - Xinit2) > 1.0e-20)
        {
            Xinit1 = Xinit1 + (Fwd - Fwd1) * (Xinit2 - Xinit1) / (Fwd2 - Fwd1);
            err    = BVMApproxPrice(
                Xinit1, XStar, GStar, Sigma, Lambda, Maturity, NHermite, 0, NULL, NULL, &Fwd1);
        }
        else
        {
            NIter = 20;
        }

        NIter++;
    }

    Xinit = Xinit1;

    err = BVMApproxPrice(
        Xinit,
        XStar,
        GStar,
        Sigma,
        Lambda,
        Maturity,
        NHermite,
        NStrikes,
        Strikes,
        ImpliedVol,
        FwdOutPut);

    return err;
}

static void compute_I(double d, long n, double* I)
{
    int i;

    I[0] = 1.0 - norm_accurate(d);
    I[1] = d * exp(-0.5 * d * d) * INV_SQRT_TWO_PI - d * I[0];

    for (i = 2; i < n; ++i)
    {
        I[i] = (i - 1) * I[i - 2] - d * I[i - 1];
    }
}

double g_function(double x, double xStar, double gStar, double lambda, int nderiv)
{
    double res;
    double C;

    switch (nderiv)
    {
    case 0:
        if (x > 50.0)
        {
            C   = exp(-xStar) * (exp(lambda * gStar) - 1);
            res = log(C) / lambda + exp(-x) / (C * lambda) + 0.5 * exp(-2 * x) / (C * C * lambda);
        }
        else
        {
            res = log(1 + exp(x - xStar) * (exp(lambda * gStar) - 1)) / lambda;
        }
        break;

    case 1:
        res = -exp(x - xStar) * (exp(lambda * gStar) - 1) /
              (lambda * (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)));
        break;

    case 2:
        res = exp(x - xStar) * (exp(lambda * gStar) - 1) / lambda /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar));
        break;

    case 3:
        res = -exp(x - xStar) * (exp(lambda * gStar) - 1) *
              (1 + exp(x - xStar) - exp(x - xStar + lambda * gStar)) / lambda /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar));
        break;

    case 4:
        res = exp(x - xStar) * (exp(lambda * gStar) - 1) *
              (1 + exp(2 * (x - xStar)) + exp(2 * (x - xStar + lambda * gStar)) -
               2 * exp(2 * (x - xStar) + lambda * gStar) - 4 * exp(x - xStar + 2 * lambda * gStar) +
               4 * exp(x - xStar)) /
              lambda / (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
              (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar));
        break;

    case 5:
        res =
            -exp(x - xStar) * (exp(lambda * gStar) - 1) *
            (1 + 11 * exp(x - xStar) - 22 * exp(2 * (x - xStar) + lambda * gStar) -
             exp(3 * (x - xStar + lambda * gStar)) + 11 * exp(2 * (x - xStar)) -
             3 * exp(3 * (x - xStar) + lambda * gStar) +
             3 * exp(3 * (x - xStar) + 2 * lambda * gStar) - 11 * exp(x - xStar + lambda * gStar) +
             exp(3 * (x - xStar)) + 11 * exp(2 * (x - xStar + lambda * gStar))) /
            lambda / (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
            (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
            (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
            (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar)) /
            (exp(x - xStar) - 1 - exp(x - xStar + lambda * gStar));
        break;

    case 6:
        res = log(1 + exp(x - xStar) * (exp(lambda * gStar) - 1)) / lambda;
        break;

    case 7:
        res = log(1 + exp(x - xStar) * (exp(lambda * gStar) - 1)) / lambda;
        break;

    case 8:
        res = log(1 + exp(x - xStar) * (exp(lambda * gStar) - 1)) / lambda;
        break;
    }

    return res;
}

Err BVMApproxAccurate(
    double  Xinit,
    double  Strike,
    double  Maturity,
    double  XStar,
    double  GStar,
    double  sigma,
    double  lambda,
    long    Order,
    double* Price)
{
    Err     err = NULL;
    int     i;
    double  price;
    double  vol, powvol, facti;
    double  H;
    double* I_func = NULL;

    if (Order > 5)
    {
        smessage("Order should be lower than 5 in BVMApproxAccurate");
        Order = 5;
        // err = "Order should be lower than 5 in BVMApproxAccurate";
    }

    I_func = (double*)calloc(Order, sizeof(double));
    if (!I_func)
    {
        smessage("Memory Allocation failed in BVMApproxAccurate");
        err = "Memory Allocation failed in BVMApproxAccurate";
        goto FREE_RETURN;
    }

    H   = XStar + log((exp(lambda * Strike) - 1) / (exp(lambda * GStar) - 1));
    vol = sigma * sqrt(Maturity);
    compute_I((H - Xinit) / vol, Order, I_func);
    powvol = 1.0;
    facti  = 1.0;
    price  = 0;
    for (i = 1; i < Order; ++i)
    {
        powvol = powvol * vol;
        facti  = i * facti;
        price += g_function(H, XStar, GStar, lambda, i) * powvol / facti * I_func[i];
    }

    *Price = price;

FREE_RETURN:

    if (I_func)
        free(I_func);

    return err;
}

Err BVMApproxAccurate2(
    double  Xinit,
    double  Strike,
    double  Maturity,
    double  XStar,
    double  GStar,
    double  sigma,
    double  lambda,
    long    NQuadrature,
    double* Price)
{
    Err     err = NULL;
    int     i;
    double  price;
    double  vol;
    double  H;
    double* x = NULL;
    double* w = NULL;
    double  test;

    x = (double*)calloc(NQuadrature + 1, sizeof(double));
    w = (double*)calloc(NQuadrature + 1, sizeof(double));
    if ((!x) || (!w))
    {
        smessage("Memory Allocation failed in BVMApproxAccurate2");
        err = "Memory Allocation failed in BVMApproxAccurate2";
        goto FREE_RETURN;
    }

    H    = XStar + log((exp(lambda * Strike) - 1) / (exp(lambda * GStar) - 1));
    vol  = sigma * sqrt(Maturity);
    test = g_function(H, XStar, GStar, lambda, 0);

    err = GaussianIntegral((H - Xinit) / vol, 1000, 0, 1, NQuadrature, x, w);

    price = 0;
    for (i = 1; i <= NQuadrature; ++i)
    {
        price += w[i] * DMAX(0, g_function(Xinit + vol * x[i], XStar, GStar, lambda, 0) - Strike);
    }

    *Price = price;

FREE_RETURN:

    if (x)
        free(x);
    if (w)
        free(w);

    return err;
}
