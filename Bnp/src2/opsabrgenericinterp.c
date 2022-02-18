/* ==========================================================================
   FILE_NAME:	opsabrgenericcalibinterp.c

   PURPOSE:		MODIFIED IMPLEMENTATION OF SABR Generic with smile extrapolation using a
   beta model

   DATE:		24/03/04

   AUTHOR:		J.B. - P.T.
   ========================================================================== */

#include "opsabrgenericinterp.h"

#include "math.h"
#include "opfnctns.h"
#include "opsabrgenericcalib.h"

// Nag routines for the Bessel function
#include <nag.h>
#include <nagg01.h>
#include <nags.h>

// NAG for the optimizer
#include <nage04.h>
#include <nagx02.h>

static void NAG_CALL
BMMlsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);
static void NAG_CALL
BMM2lsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);
static void NAG_CALL
BMM3lsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);
static void NAG_CALL
BMMGenlsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);
static void NAG_CALL
BMMlsqfunqwithpi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);
static void NAG_CALL
BMM2lsqfunqwithzeta(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm);

#define UNSurSQRT2Pi 0.39894228040143267793994605993438

//****************************************************************************************************************
//************************* Aux Functions for moments of Bernouilli and Log **************************************
//****************************************************************************************************************

double FindPiFromSkew(double s)
{
    double temp1, temp2;

    temp1 = 0.5 * (1 - s / sqrt(4 + s * s));
    temp2 = (s * s - sqrt(s * s * (4 + s * s)) + 4) / (2 * s * s + 8);
    return 0.5 * (1 - s / sqrt(4 + s * s));
}

double SkewTri(double Pi0, double Pi)
{
    return (1 - 2 * Pi) / sqrt((1 - Pi) * (1 - Pi0) * Pi);
}

double KurtTri(double Pi0, double Pi)
{
    return 1 / (1 - Pi0) * (Pi * Pi / (1 - Pi) + (1 - Pi) * (1 - Pi) / Pi);
}

double SkewLog(double alpha)
{
    double alpsq;

    alpsq = alpha * alpha;
    return (exp(3 * alpsq) - 3 * exp(alpsq) + 2) / pow(exp(alpsq) - 1, 1.5);
}

double KurtLog(double alpha)
{
    double alpsq = alpha * alpha;
    return (exp(6 * alpsq) - 4 * exp(3 * alpsq) + 6 * exp(alpsq) - 3) / pow(exp(alpsq) - 1, 2);
}

double MomentLog(double alpha, double N)
{
    return exp(0.5 * N * alpha * alpha * (N - 1));
}

double D(int j, double beta, int N)
{
    double res = 1;
    int    i;

    for (i = 0; i <= N; i++)
    {
        if (i != j)
            res *= 1 / (pow(beta + j, 2) - pow(beta + i, 2));
    }
    res *= pow(2, N);
    return res;
}

double Mtmt(int N, double sigma, double alpha, double T)
{
    int    j;
    double lambda, eta, res = 0;

    lambda = 2 * alpha;
    eta    = -0.5 * alpha;

    for (j = 0; j <= N; j++)
    {
        res += D(j, eta / lambda, N) * exp((lambda * lambda * j * j * 0.5 + lambda * eta * j) * T);
    }
    res *= fact(N) / pow(lambda, 2 * N);
    return res * pow(sigma, 2 * N);
}

// moment function for integral of Sigma(t)^2
double VarianceMoment(double var, int N)
{
    double res, temp;
    double factor = fact(N) / pow(2 * var, N);
    int    i, j;
    res = 0;
    for (j = 0; j < N + 1; j++)
    {
        temp = 1;
        for (i = 0; i < N + 1; i++)
        {
            if (i != j)
            {
                temp /= (j + i - 0.5) * (j - i);
            }
        }
        res += temp * exp((2 * j - 1) * j * var);
    }
    res *= factor;
    return res;
}

double CalibPiPi0VetaSub(double Pi, double Pi0)
{
    return pow(1 - 2 * Pi, 2) / ((1 - Pi0) * Pi * (1 - Pi));
}

// Match 3 first moments with a given veta (Pi0)
// modified to match vol and not cumul variance
Err CalibPiPi0Veta(
    double  SigmaBeta,
    double  alpha,
    double  veta,
    double  Pi0,
    double* ResPi0,
    double* ResPi,
    double* ResSig0,
    double* ResSig1,
    double* ResSig2)
{
    // assumes that Pi0 and veta ranges have been checked
    double m1, m2, m3, m4;
    double Mt1, Mt2, Mt3, Mt4;
    double Pi;
    double temp, target, val, lambda;
    double magicnumber = 0.5;
    double threshold   = 0.01;  // mimimum value of Sig1 and Sig2 in percent of Sig
    int    i;
    Err    err = NULL;
    m1         = MomentLog(alpha, 1);
    m2         = MomentLog(alpha, 2);
    m3         = MomentLog(alpha, 3);
    m4         = MomentLog(alpha, 4);
    Mt1        = m1;
    Mt2        = m2 - m1 * m1;
    Mt3        = m3 - 3 * m2 * m1 + 2 * pow(m1, 3);
    Mt4        = m4 - 4 * m3 * m1 + 6 * m2 * m1 * m1 - 3 * pow(m1, 4);

    // Pi calibrate
    temp   = magicnumber / 2;
    val    = 0;
    target = veta * veta * Mt3 * Mt3 / Mt2 / Mt2 / Mt2;
    for (i = 0; i < 15; i++)
    {
        val = CalibPiPi0VetaSub(temp, Pi0);
        if (val < target)
            temp -= magicnumber / pow(2, i + 2);
        else
            temp += magicnumber / pow(2, i + 2);
    }
    Pi = temp;

    // lambda is abs(X2-X1) lambda is negative because skew always positive
    lambda = -sqrt(Mt2 / (1 - Pi0) / (1 - Pi) / Pi);

    // Test if outputs > 0
    if (Mt1 + Pi * lambda < threshold)
    {
        err = "Alpha too high in CalibPiPi0Veta";
    }
    // skew is always positive for log normal and Pi <0.5 so we know where to put which one is X1
    // and X2
    *ResSig0 = SigmaBeta * Mt1;
    *ResSig1 = SigmaBeta * (Mt1 - (1 - Pi) * lambda);
    *ResSig2 = SigmaBeta * (Mt1 + Pi * lambda);

    *ResPi  = Pi;
    *ResPi0 = Pi0;

    return err;
}

//****************************************************************************************************************
//******************************** SABR-like Mapping Functions ***************************************************
//****************************************************************************************************************

Err BMMGetStates(
    double  forward,
    double  maturity,
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  pi,
    double* pSig1,
    double* pSig2,
    double* pFwd1,
    double* pFwd2)
{
    Err    err = NULL;
    double lambda;
    double b, alphasqrmat, unsurpiunmoinspi;
    double threshold = 0.01;  // mimimum value of Sig1 and Sig2 in percent of Sig

    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((pi < 0.0) || (pi > 1.0))
    {
        err = "Error: Pi should be between 0 and 1";
        goto FREE_RETURN;
    }

    //	pi = FindPiFromSkew(pi * SkewLog(alpha));

    if (pi * (1 - pi) != 0)
    {
        unsurpiunmoinspi = 1 / (pi * (1 - pi));

        sigma *= sqrt(1 - rho * rho);  // because of pb on the FX for high rho
        lambda = rho * sqrt(unsurpiunmoinspi) * sigma * pow(forward, beta - 1) * sqrt(maturity);

        *pFwd1 = forward * (1 + (1 - pi) * lambda);
        *pFwd2 = forward * (1 - (pi)*lambda);

        if (fabs(alpha) > 1e-5)
        {
            alphasqrmat = alpha * alpha * maturity;
            b           = sigma * sqrt(((exp(alphasqrmat) - 1) / alphasqrmat - 1));
        }
        else
            b = 0;

        *pSig1 = sigma + sqrt((1 - pi) / pi) * b;
        *pSig2 = sigma - sqrt(pi / (1 - pi)) * b;
        if (*pSig2 <= threshold * sigma)
        {
            err = "Alpha too high";
            goto FREE_RETURN;
        }
        if ((*pFwd1) * (*pFwd2) <= 0)
        {
            err = "Forwards too far (Correlation to High or to Low)";
            goto FREE_RETURN;
        }
    }
    else
    {
        *pFwd1 = forward;
        *pFwd2 = forward;

        *pSig1 = sigma;
        *pSig2 = sigma;
    }

FREE_RETURN:
    return err;
}

Err BMM2GetStates(
    double  forward,
    double  maturity,
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  zeta,
    double* pSig1,
    double* pSig2,
    double* pFwd1,
    double* pFwd2,
    double* pPi)
{
    Err    err = NULL;
    double alphasqrmat, b, newpi, lambda;
    double threshold = 0.001;  // mimimum value of Sig1 and Sig2 in percent of Sig
    double temprho   = 0;

    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if (sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((zeta < 0.0))
    {
        err = "Error: Zeta should be positive";
        goto FREE_RETURN;
    }

    newpi = FindPiFromSkew(zeta * SkewLog(0.5 * alpha * sqrt(maturity)));

    if (newpi * (1 - newpi) * alpha != 0)
    {
        sigma *= sqrt(1 - rho * rho);
        lambda = rho / sqrt(newpi * (1 - newpi)) * sigma * pow(forward, beta - 1) * sqrt(maturity);

        *pFwd1 =
            max(forward * (1 + (1 - newpi) * lambda), 1e-4);  // new to handle Westminster newton
        *pFwd2 = max(forward * (1 - newpi * lambda), 1e-4);   // new to handle Westminster newton
        if ((*pFwd1) * (*pFwd2) <= 0)
        {
            err = "Forwards too far (Correlation to High or to Low)";
            goto FREE_RETURN;
        }

        alphasqrmat = alpha * alpha * maturity;
        b           = sigma * sqrt(((exp(alphasqrmat) - 1) / alphasqrmat - 1));

        *pSig1 = sigma + sqrt((1 - newpi) / newpi) * b;
        *pSig2 = sigma - sqrt(newpi / (1 - newpi)) * b;
        if (*pSig2 <= threshold * sigma)
        {
            *pSig2 = threshold * sigma;
            //	err = "Alpha too high";	goto FREE_RETURN;
        }

        *pPi = newpi;
    }
    else
    {
        *pFwd1 = forward;
        *pFwd2 = forward;
        *pSig1 = sigma;
        *pSig2 = sigma;
        *pPi   = 0.5;
    }

FREE_RETURN:
    return err;
}

// 3 states mapping
Err BMM3GetStates(
    double  forward,
    double  maturity,
    double  sigma,
    double  alpha,
    double  beta,
    double  rho,
    double  veta,
    double  pi0,
    double* pSig,
    double* pSig1,
    double* pSig2,
    double* pFwd,
    double* pFwd1,
    double* pFwd2,
    double* pPi0,
    double* pPi)
{
    Err    err = NULL;
    double pi;
    double lambda;

    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if (sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((veta < 0.0) || (veta > 1.0))
    {
        err = "Error: Veta should be between 0 and 1";
        goto FREE_RETURN;
    }
    if ((pi0 < 0.0) || (pi0 > 1.0))
    {
        err = "Error: Pi0 should be between 0 and 1";
        goto FREE_RETURN;
    }

    if (pi0 == 1)
    {
        *pFwd  = forward;
        *pFwd1 = forward;
        *pFwd2 = forward;
        *pSig  = sigma;
        *pSig1 = sigma;
        *pSig2 = sigma;
        *pPi0  = pi0;
        *pPi   = 0.5;
    }
    else
    {
        // Version with given Pi0, kurtosis of vol is not matched
        err =
            CalibPiPi0Veta(sigma, alpha * sqrt(maturity), veta, pi0, pPi0, pPi, pSig, pSig1, pSig2);
        if (err)
            goto FREE_RETURN;
        pi = *pPi;

        lambda = (rho * UNSurSQRT2Pi * (*pSig) * pow(forward, beta - 1) * sqrt(maturity)) /
                 ((1 - pi0) * pi * (1 - pi));

        *pFwd  = forward;
        *pFwd1 = forward * (1 + (1 - pi) * lambda);
        *pFwd2 = forward * (1 - (pi)*lambda);
        if ((*pFwd1) * (*pFwd2) <= 0)
        {
            err = "Forwards too far (Correlation to High or to Low)";
            goto FREE_RETURN;
        }
    }

FREE_RETURN:
    return err;
}

Err BMMGenGetStates(
    double           forward,
    double           maturity,
    double           sigma,
    double           alpha,
    double           a,
    double           b,
    double           c,
    double           rho,
    double           pi,
    SrtDiffusionType TypeVolLoc,
    double*          pSig1,
    double*          pSig2,
    double*          pFwd1,
    double*          pFwd2)
{
    double         temp, alphasqrmat, unsurpiunmoinspi, lambda;
    FuncVolLocType vol_loc;
    Err            err       = NULL;
    double         threshold = 0.01;  // mimimum value of Sig1 and Sig2 in percent of Sig

    if (vol_parameters_check_from_SrtDiffusionType(a, b, c, 30, TypeVolLoc) == 0)
    {
        err = "Generic parameters of vol loc out of range";
        goto FREE_RETURN;
    }
    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if (sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((pi < 0.0) || (pi > 1.0))
    {
        err = "Error: pi should be between 0 and 1";
        goto FREE_RETURN;
    }

    GetLocVolFromDiffusionType(TypeVolLoc, &vol_loc);
    temp = vol_loc(forward, a, b, c, 0);

    pi = FindPiFromSkew(pi * SkewLog(alpha));
    if (pi * (1 - pi) != 0)
    {
        unsurpiunmoinspi = 1 / (pi * (1 - pi));

        lambda = rho * UNSurSQRT2Pi * unsurpiunmoinspi * sigma * temp / forward * sqrt(maturity) /
                 sqrt(1 - rho * rho);

        *pFwd1 = forward * (1 + (1 - pi) * lambda);
        *pFwd2 = forward * (1 - (pi)*lambda);

        if (alpha != 0)
        {
            alphasqrmat = alpha * alpha * maturity;
            b           = sigma * sqrt(((exp(alphasqrmat) - 1) / alphasqrmat - 1));
        }
        else
            b = 0;

        *pSig1 = sigma + sqrt((1 - pi) / pi) * b;
        *pSig2 = sigma - sqrt(pi / (1 - pi)) * b;
        if (*pSig2 <= threshold * sigma)
        {
            err = "Alpha too high";
            goto FREE_RETURN;
        }
        if ((*pFwd1) * (*pFwd2) <= 0)
        {
            err = "Forwards too far (Correlation to High or to Low)";
            goto FREE_RETURN;
        }
    }
    else
    {
        *pFwd1 = forward;
        *pFwd2 = forward;

        *pSig1 = sigma;
        *pSig2 = sigma;
    }
FREE_RETURN:
    return err;
}

//****************************************************************************************************************
//************************************* BMM Pricing Functions ****************************************************
//****************************************************************************************************************

Err BMM_Option_Price(
    double         forward,
    double         strike,
    double         maturity,
    double         sigma,
    double         alpha,
    double         beta,
    double         rho,
    double         pi,
    SrtCallPutType optiontype,
    double*        price)
{
    double res = 0, Fwd1, Fwd2, Sig1, Sig2;
    Err    err = NULL;

    err = BMMGetStates(forward, maturity, sigma, alpha, beta, rho, pi, &Sig1, &Sig2, &Fwd1, &Fwd2);
    if (err)
        goto FREE_RETURN;
    if (beta >= 0)
    {
        // Fwd up   Vol Up
        res = pi * srt_f_optblkschbeta(Fwd1, strike, Sig1, maturity, beta, 1, optiontype, PREMIUM);
        // Fwd up   Vol Dn
        res += (1 - pi) *
               srt_f_optblkschbeta(Fwd2, strike, Sig2, maturity, beta, 1, optiontype, PREMIUM);
    }
    else
    {
        // Fwd up   Vol Up
        res = pi * srt_f_optblknrm(Fwd1, strike, Sig1, maturity, 1.0, optiontype, PREMIUM);
        // Fwd up   Vol Dn
        res += (1 - pi) * srt_f_optblknrm(Fwd2, strike, Sig2, maturity, 1.0, optiontype, PREMIUM);
    }
    *price = res;

FREE_RETURN:
    return err;
}

Err BMM_Option_Price_quick(
    double         forward,
    double         strike,
    double         maturity,
    double         sigma,
    double         alpha,
    double         beta,
    double         rho,
    double         pi,
    SrtCallPutType optiontype,
    double*        price)
{
    double res = 0, Fwd1, Fwd2, Sig1, Sig2;
    Err    err = NULL;

    err = BMMGetStates(forward, maturity, sigma, alpha, beta, rho, pi, &Sig1, &Sig2, &Fwd1, &Fwd2);
    if (err)
        goto FREE_RETURN;
    if (beta >= 0)
    {
        // Fwd up   Vol Up
        res = pi *
              srt_f_optblkschbetaquick(Fwd1, strike, Sig1, maturity, beta, 1, optiontype, PREMIUM);
        // Fwd up   Vol Dn
        res += (1 - pi) *
               srt_f_optblkschbetaquick(Fwd2, strike, Sig2, maturity, beta, 1, optiontype, PREMIUM);
    }
    else
    {
        // Fwd up   Vol Up
        res = pi * srt_f_optblknrm(Fwd1, strike, Sig1, maturity, 1.0, optiontype, PREMIUM);
        // Fwd up   Vol Dn
        res += (1 - pi) * srt_f_optblknrm(Fwd2, strike, Sig2, maturity, 1.0, optiontype, PREMIUM);
    }
    *price = res;

FREE_RETURN:
    return err;
}

Err BMM_Option_Price_From_States(
    double         Strike,
    double         Maturity,
    double         Beta,
    double         Fwd1,
    double         Fwd2,
    double         Sig1,
    double         Sig2,
    double         Pi,
    SrtCallPutType optiontype,
    double*        Price)
{
    // Fwd up   Vol Up
    *Price = Pi * srt_f_optblkschbeta(Fwd1, Strike, Sig1, Maturity, Beta, 1, optiontype, PREMIUM);
    // Fwd up   Vol Dn
    *Price +=
        (1 - Pi) * srt_f_optblkschbeta(Fwd2, Strike, Sig2, Maturity, Beta, 1, optiontype, PREMIUM);
    return NULL;
}

Err BMM_Option_Price2(
    double         forward,
    double         strike,
    double         maturity,
    double         sigma,
    double         alpha,
    double         beta,
    double         rho,
    double         zeta,
    SrtCallPutType optiontype,
    double*        price)
{
    double res = 0;
    double Sig1, Sig2, Fwd1, Fwd2, Pi;
    Err    err = NULL;
    err        = BMM2GetStates(
        forward, maturity, sigma, alpha, beta, rho, zeta, &Sig1, &Sig2, &Fwd1, &Fwd2, &Pi);
    if (err)
        goto FREE_RETURN;
    // Fwd up   Vol Up
    res = Pi * srt_f_optblkschbeta(Fwd1, strike, Sig1, maturity, beta, 1, optiontype, PREMIUM);
    // Fwd up   Vol Dn
    res +=
        (1 - Pi) * srt_f_optblkschbeta(Fwd2, strike, Sig2, maturity, beta, 1, optiontype, PREMIUM);
    *price = res;

FREE_RETURN:
    return err;
}

Err BMM_Option_Price3(
    double         forward,
    double         strike,
    double         maturity,
    double         sigma,
    double         alpha,
    double         beta,
    double         rho,
    double         veta,
    double         pi0,
    SrtCallPutType optiontype,
    double*        price)
{
    // assume range check have been performed
    Err    err = NULL;
    double res = 0;
    double Sig, Sig1, Sig2, Fwd, Fwd1, Fwd2, Pi0, Pi;
    err = BMM3GetStates(
        forward,
        maturity,
        sigma,
        alpha,
        beta,
        rho,
        veta,
        pi0,
        &Sig,
        &Sig1,
        &Sig2,
        &Fwd,
        &Fwd1,
        &Fwd2,
        &Pi0,
        &Pi);
    if (err)
        goto FREE_RETURN;
    // Fwd up   Vol Up
    res = (1 - Pi0) * Pi *
          srt_f_optblkschbeta(Fwd1, strike, Sig1, maturity, beta, 1, optiontype, PREMIUM);
    // Fwd up   Vol Dn
    res += (1 - Pi0) * (1 - Pi) *
           srt_f_optblkschbeta(Fwd2, strike, Sig2, maturity, beta, 1, optiontype, PREMIUM);
    // Fwd Vol
    if (Pi0 != 0)
        res += Pi0 * srt_f_optblkschbeta(Fwd, strike, Sig, maturity, beta, 1, optiontype, PREMIUM);
    *price = res;

FREE_RETURN:
    return err;
}

Err BMM_Option_Price3_From_States(
    double         forward,
    double         strike,
    double         maturity,
    double         beta,
    double         fwd,
    double         fwd1,
    double         fwd2,
    double         sig,
    double         sig1,
    double         sig2,
    double         pi0,
    double         pi,
    SrtCallPutType optiontype,
    double*        price)
{
    double res = 0;
    // Fwd up   Vol Up
    res = (1 - pi0) * pi *
          srt_f_optblkschbeta(fwd1, strike, sig1, maturity, beta, 1, optiontype, PREMIUM);
    // Fwd up   Vol Dn
    res += (1 - pi0) * (1 - pi) *
           srt_f_optblkschbeta(fwd2, strike, sig2, maturity, beta, 1, optiontype, PREMIUM);
    // Fwd Vol
    if (pi0 != 0)
        res += pi0 * srt_f_optblkschbeta(fwd, strike, sig, maturity, beta, 1, optiontype, PREMIUM);
    *price = res;
    return NULL;
}

Err BMMGen_Option_Price(
    double           forward,
    double           strike,
    double           maturity,
    double           sigma,
    double           alpha,
    double           a,
    double           b,
    double           c,
    double           rho,
    double           pi,
    SrtCallPutType   optiontype,
    double*          price,
    SrtDiffusionType TypeVolLoc)
{
    double res = 0, tempprice;
    double Fwd1, Fwd2, Sig1, Sig2;
    Err    err = NULL;

    err = BMMGenGetStates(
        forward, maturity, sigma, alpha, a, b, c, rho, pi, TypeVolLoc, &Sig1, &Sig2, &Fwd1, &Fwd2);
    if (err)
        goto FREE_RETURN;

    // Fwd 1   Vol 1
    err = BMMGenPrice(Fwd1, strike, Sig1, maturity, a, b, c, &tempprice, TypeVolLoc);
    if (err)
        goto FREE_RETURN;
    res = pi * tempprice;
    // Fwd 2   Vol 2
    err = BMMGenPrice(Fwd2, strike, Sig2, maturity, a, b, c, &tempprice, TypeVolLoc);
    if (err)
        goto FREE_RETURN;
    res += (1 - pi) * tempprice;

    *price = res;

FREE_RETURN:
    return NULL;
}

//****************************************************************************************************************
//************************************* BMM Calibration Functions ************************************************
//****************************************************************************************************************

Err BMMCalibOnPrices(
    double  forward,
    double  maturity,
    double  atmvol,
    double  alpha,
    double  beta,
    double  rho,
    double  pi,
    int     nCalibStrikes,
    double* CalibStrikes,
    double* CalibVols,
    double* CalibWeights,
    double* CalibSigmabeta,
    double* CalibBeta,
    double* CalibAlpha,
    double* CalibRho,
    double* CalibPi,
    double* Calibres)
{
    double          pistep, fjac[10][10], fvec[10], x[5];
    int             exit;
    double          res = -1;
    int             i, tdj = 10;
    static NagError fail;
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer
    Opt_params      params;
    Err             err = NULL;
    // dummy parameters to check range
    double tpFwd1, tpFwd2, tpSig1, tpSig2, Sigma, alphastep;
    double ATMVega;
    e04xxc(&options);

#ifdef _DEBUG
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    fail.print          = FALSE;
/*		options.print_level	= Nag_Soln_Iter_Full;
                options.list		= TRUE;
                fail.print			= TRUE;*/
#endif _DEBUG

#ifdef NDEBUG
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    fail.print          = FALSE;
#endif NDEBUG

    // for weigthing factor
    ATMVega = srt_f_optblksch(forward, forward, atmvol, maturity, 1.0, SRT_CALL, VEGA) * 100;

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = nCalibStrikes;
    params.Vec_marketvols = CalibVols;
    params.Vec_strikes    = CalibStrikes;
    params.Vec_weights    = CalibWeights;
    params.Beta           = beta;
    params.Pi             = pi;
    params.VegaFact       = 1 / ATMVega;
    params.ExitFlag       = 0;

    nagcomm.p = (void*)&params;

    // trick to help calibrate and check ranges
    // Move alpha down if needed
    Sigma = pow(forward, 1 - beta) * atmvol;
    BMMGetStates(
        forward, maturity, Sigma, alpha, beta, rho, pi, &tpSig1, &tpSig2, &tpFwd1, &tpFwd2);

    alphastep = alpha;
    while ((tpSig2 < 0.01 * Sigma) && (alphastep > 0.01 * alpha))
    {
        alphastep = 0.9 * alphastep;
        err       = BMMGetStates(
            forward, maturity, Sigma, alphastep, beta, rho, pi, &tpSig1, &tpSig2, &tpFwd1, &tpFwd2);
    }

    if (err)
        return err;
    alpha = alphastep;
    // End of shitty trick

    x[0] = sqrt(pow(forward, 1 - beta) * atmvol);
    x[1] = sqrt(alpha);
    x[2] = acos(rho);

    exit = 0;
    i    = 0;

    if (nCalibStrikes <= 3)
    {
        nag_opt_lsq_no_deriv(
            nCalibStrikes,
            3,
            BMMlsqfunqnopi,
            x,
            &res,
            fvec,
            (double*)fjac,
            tdj,
            &options,
            &nagcomm,
            &fail);
        *CalibBeta = beta;
        *CalibPi   = pi;
        *Calibres  = sqrt(fvec[0] * fvec[0] + fvec[1] * fvec[1] + fvec[2] * fvec[2]) / 3 / ATMVega;
    }
    else
    {
        pistep = pi;
        while ((exit == 0) && (i < 10))
        {
            x[3] = sqrt(acos(pistep));
            nag_opt_lsq_no_deriv(
                nCalibStrikes,
                4,
                BMMlsqfunqwithpi,
                x,
                &res,
                fvec,
                (double*)fjac,
                tdj,
                &options,
                &nagcomm,
                &fail);
            if (fail.errnum == 0)
                exit = 1;
            else
            {
                x[0] = sqrt(pow(forward, 1 - beta) * atmvol);
                x[1] = sqrt(alpha);
                x[2] = acos(rho);
                INIT_FAIL(fail);
                pistep = pistep * 0.75;
                i++;
            }
        }
        *CalibBeta = beta;
        *CalibPi   = pow(cos(x[3]), 2);
        *Calibres =
            sqrt(fvec[0] * fvec[0] + fvec[1] * fvec[1] + fvec[2] * fvec[2] + fvec[3] * fvec[3]) /
            4 / ATMVega;
    };

    *CalibSigmabeta = x[0] * x[0];
    *CalibAlpha     = x[1] * x[1];
    *CalibRho       = cos(x[2]);
    return NULL;
};

Err BMM2CalibOnPrices(
    double  forward,
    double  maturity,
    double  atmvol,
    double  alpha,
    double  beta,
    double  rho,
    double  zeta,
    int     nCalibStrikes,
    double* CalibStrikes,
    double* CalibVols,
    double* CalibWeights,
    double* CalibSigmabeta,
    double* CalibBeta,
    double* CalibAlpha,
    double* CalibRho,
    double* CalibZeta,
    double* Calibres)
{
    double          zetastep, fjac[10][10], fvec[10], x[5];
    int             exit;
    double          res = -1;
    int             i, tdj = 10;
    static NagError fail;
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer
    Opt_params      params;
    Err             err = NULL;
    // dummy parameters to check range
    double tpFwd1, tpFwd2, tpSig1, tpSig2, Sigma, alphastep, tpPi;
    double ATMVega;
    e04xxc(&options);

#ifdef _DEBUG
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    fail.print          = FALSE;
/*		options.print_level	= Nag_Soln_Iter_Full;
                options.list		= TRUE;
                fail.print			= TRUE;*/
#endif _DEBUG

#ifdef NDEBUG
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    fail.print          = FALSE;
#endif NDEBUG

    // for weigthing factor
    ATMVega = srt_f_optblksch(forward, forward, atmvol, maturity, 1.0, SRT_CALL, VEGA) * 100;

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = nCalibStrikes;
    params.Vec_marketvols = CalibVols;
    params.Vec_strikes    = CalibStrikes;
    params.Vec_weights    = CalibWeights;
    params.Beta           = beta;
    params.Pi             = zeta;
    params.VegaFact       = 1 / ATMVega;
    params.ExitFlag       = 0;

    nagcomm.p = (void*)&params;

    // trick to help calibrate and check ranges
    // Move alpha down if needed
    Sigma = pow(forward, 1 - beta) * atmvol;
    BMM2GetStates(
        forward,
        maturity,
        Sigma,
        alpha,
        beta,
        rho,
        zeta,
        &tpSig1,
        &tpSig2,
        &tpFwd1,
        &tpFwd2,
        &tpPi);

    alphastep = alpha;
    while ((tpSig2 < 0.01 * Sigma) && (alphastep > 0.01 * alpha))
    {
        alphastep = 0.9 * alphastep;
        err       = BMM2GetStates(
            forward,
            maturity,
            Sigma,
            alphastep,
            beta,
            rho,
            zeta,
            &tpSig1,
            &tpSig2,
            &tpFwd1,
            &tpFwd2,
            &tpPi);
    }

    if (err)
        return err;
    alpha = alphastep;
    // End of shitty trick

    x[0] = sqrt(pow(forward, 1 - beta) * atmvol);
    x[1] = sqrt(alpha);
    x[2] = acos(rho);

    exit = 0;
    i    = 0;

    if (nCalibStrikes <= 3)
    {
        nag_opt_lsq_no_deriv(
            nCalibStrikes,
            3,
            BMM2lsqfunqnopi,
            x,
            &res,
            fvec,
            (double*)fjac,
            tdj,
            &options,
            &nagcomm,
            &fail);
        *CalibBeta = beta;
        *CalibZeta = zeta;
        *Calibres  = sqrt(fvec[0] * fvec[0] + fvec[1] * fvec[1] + fvec[2] * fvec[2]) / 3 / ATMVega;
    }
    else
    {
        zetastep = zeta;
        while ((exit == 0) && (i < 10))
        {
            x[3] = sqrt(acos(zetastep));
            nag_opt_lsq_no_deriv(
                nCalibStrikes,
                4,
                BMM2lsqfunqwithzeta,
                x,
                &res,
                fvec,
                (double*)fjac,
                tdj,
                &options,
                &nagcomm,
                &fail);
            if (fail.errnum == 0)
                exit = 1;
            else
            {
                x[0] = sqrt(pow(forward, 1 - beta) * atmvol);
                x[1] = sqrt(alpha);
                x[2] = acos(rho);
                INIT_FAIL(fail);
                zetastep = zetastep * 0.75;
                i++;
            }
        }
        *CalibBeta = beta;
        *CalibZeta = cos(pow(x[3], 2));
        *Calibres =
            sqrt(fvec[0] * fvec[0] + fvec[1] * fvec[1] + fvec[2] * fvec[2] + fvec[3] * fvec[3]) /
            4 / ATMVega;
    };

    *CalibSigmabeta = x[0] * x[0];
    *CalibAlpha     = x[1] * x[1];
    *CalibRho       = cos(x[2]);
    return NULL;
};

Err BMMCalibOnSabrStates(
    double           forward,
    double           maturity,
    double           sigma,
    double           beta,
    double           alpha,
    double           rho,
    double           pi,
    double           NStd,
    double*          Fwd1,
    double*          Fwd2,
    double*          Sig1,
    double*          Sig2,
    double*          Pi,
    SrtDiffusionType TypeInput,
    double*          Calibres)
{
    int    i, switchtonormal = 0;
    double CalibSigma, CalibBeta, CalibAlpha, CalibRho, CalibPi, oldbeta;
    double CalibStrikes[3];  //	= dvector(1,3);
    double CalibPrices[3];   //		= dvector(1,3);
    double CalibWeights[3];  //		= dvector(1,3);
    Err    err = NULL;

    if ((maturity == 0) || (pi == 0))
        return NULL;
    if (beta < 0)
    {
        oldbeta        = beta;
        beta           = 0.0;
        switchtonormal = 1;
    }

    srt_f_optsarbvol(
        forward, forward, maturity, sigma, alpha, beta, rho, TypeInput, SRT_LOGNORMAL, &sigma);

    // reset forward to 1 and adjust the srikes

    // Fill 3 strikes
    CalibStrikes[0] = forward * exp(-NStd * sqrt(maturity) * sigma);
    CalibStrikes[1] = forward;
    CalibStrikes[2] = forward * exp(NStd * sqrt(maturity) * sigma);

    CalibWeights[0] = 1.0;
    CalibWeights[1] = 1.0;
    CalibWeights[2] = 1.0;

    for (i = 0; i < 3; i++)
        srt_f_optsarbvol(
            forward,
            CalibStrikes[i],
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            SRT_LOGNORMAL,
            SRT_LOGNORMAL,
            &CalibPrices[i]);  // CalibPrices = LOG vols sabr

    for (i = 0; i < 3; i++)
        CalibPrices[i] = srt_f_optblksch(
            forward, CalibStrikes[i], CalibPrices[i], maturity, 1.0, SRT_CALL, PREMIUM);

    if (switchtonormal == 1)
    {
        beta = oldbeta;
    }

    err = BMMCalibOnPrices(
        forward,
        maturity,
        sigma,
        alpha,
        beta,
        rho,
        pi,
        3,
        CalibStrikes,
        CalibPrices,
        CalibWeights,
        &CalibSigma,
        &CalibBeta,
        &CalibAlpha,
        &CalibRho,
        &CalibPi,
        Calibres);
    if (err)
        return err;
    err = BMMGetStates(
        forward,
        maturity,
        CalibSigma,
        CalibAlpha,
        CalibBeta,
        CalibRho,
        CalibPi,
        Sig1,
        Sig2,
        Fwd1,
        Fwd2);
    if (err)
        return err;

    //	*Fwd1 *= forward / newforward;
    //	*Fwd2 *= forward / newforward;

    // free_dvector(CalibStrikes,1,3);
    // free_dvector(CalibPrices,1,3);

    return NULL;
};

Err BMM2CalibOnSabrStates(
    double           forward,
    double           maturity,
    double           sigma,
    double           beta,
    double           alpha,
    double           rho,
    double           zeta,
    double           NStd,
    double*          Fwd1,
    double*          Fwd2,
    double*          Sig1,
    double*          Sig2,
    double*          Pi,
    SrtDiffusionType TypeInput,
    double*          Calibres)
{
    double CalibSigma, CalibBeta, CalibAlpha, CalibRho, CalibZeta;
    Err    err = NULL;

    if (maturity <= 0)
        return NULL;

    err = BMM2CalibOnSabr(
        forward,
        maturity,
        sigma,
        alpha,
        beta,
        rho,
        zeta,
        NStd,
        TypeInput,
        &CalibSigma,
        &CalibAlpha,
        &CalibBeta,
        &CalibRho,
        &CalibZeta,
        Calibres);
    if (err)
        return err;
    err = BMM2GetStates(
        forward,
        maturity,
        CalibSigma,
        CalibAlpha,
        CalibBeta,
        CalibRho,
        CalibZeta,
        Sig1,
        Sig2,
        Fwd1,
        Fwd2,
        Pi);
    if (err)
        return err;
    return NULL;
};

Err BMMCalibOnSabr(
    double           forward,
    double           maturity,
    double           sigma,
    double           alpha,
    double           beta,
    double           rho,
    double           pi,
    double           NStd,
    SrtDiffusionType TypeInput,
    double*          CalibSigma,
    double*          CalibAlpha,
    double*          CalibBeta,
    double*          CalibRho,
    double*          CalibPi,
    double*          Calibres)
{
    int    i, switchtonormal = 0;
    Err    err = NULL;
    double CalibStrikes[3];
    double CalibPrices[3];
    double CalibWeights[3];
    double oldbeta;

    if (beta < 0)
    {
        oldbeta        = beta;
        beta           = 0.0;
        switchtonormal = 1;
    }

    srt_f_optsarbvol(
        forward, forward, maturity, sigma, alpha, beta, rho, TypeInput, SRT_LOGNORMAL, &sigma);
    // Fill the 3 strikes
    CalibStrikes[0] = forward * exp(-NStd * sqrt(maturity) * sigma);
    CalibStrikes[1] = forward;
    CalibStrikes[2] = forward * exp(NStd * sqrt(maturity) * sigma);

    CalibWeights[0] = 1.0;
    CalibWeights[1] = 1.0;
    CalibWeights[2] = 1.0;

    for (i = 0; i < 3; i++)
        srt_f_optsarbvol(
            forward,
            CalibStrikes[i],
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            SRT_LOGNORMAL,
            SRT_LOGNORMAL,
            &CalibPrices[i]);  // CalibPrices = LOG vols sabr
    for (i = 0; i < 3; i++)
        CalibPrices[i] = srt_f_optblksch(
            forward, CalibStrikes[i], CalibPrices[i], maturity, 1.0, SRT_CALL, PREMIUM);

    if (switchtonormal == 1)
    {
        beta = oldbeta;
    }

    err = BMMCalibOnPrices(
        forward,
        maturity,
        sigma,
        alpha,
        beta,
        rho,
        pi,
        3,
        CalibStrikes,
        CalibPrices,
        CalibWeights,
        CalibSigma,
        CalibBeta,
        CalibAlpha,
        CalibRho,
        CalibPi,
        Calibres);

    return err;
};

Err BMM2CalibOnSabr(
    double           forward,
    double           maturity,
    double           atmvol,
    double           alpha,
    double           beta,
    double           rho,
    double           zeta,
    double           NSTD,
    SrtDiffusionType TypeInput,
    double*          Calibsigma,
    double*          Calibalpha,
    double*          Calibbeta,
    double*          Calibrho,
    double*          Calibzeta,
    double*          Calibres)
{
    double          fjac[10][10], fvec[10], x[5];
    double          CalibPrices[5], CalibStrikes[5], CalibWeights[5];
    double          vNStd[3] = {-1, 0, 1};
    double          res      = 0;
    int             i;
    double          sigma, sigmaLog, tempvol;
    double          pFwd1, pFwd2, pSig1, pSig2, pPi;
    Err             err = NULL;
    static NagError fail;
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer
    Opt_params      params;

    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if (atmvol < 0.0)
    {
        err = "Error: Sigma should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // Compute SigmaBeta vol and atmlogvol
    srt_f_optsarbvol(
        forward, forward, maturity, atmvol, alpha, beta, rho, TypeInput, SRT_BETAVOL, &sigma);
    srt_f_optsarbvol(
        forward, forward, maturity, atmvol, alpha, beta, rho, TypeInput, SRT_LOGNORMAL, &sigmaLog);

    // compute strikes and prices:
    for (i = 0; i < 3; i++)
    {
        CalibStrikes[i] = forward * exp(vNStd[i] * NSTD * sqrt(maturity) * sigmaLog);
        srt_f_optsarbvol(
            forward,
            CalibStrikes[i],
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            SRT_BETAVOL,
            SRT_LOGNORMAL,
            &tempvol);
        CalibPrices[i] =
            srt_f_optblksch(forward, CalibStrikes[i], tempvol, maturity, 1, SRT_CALL, SRT_PREMIUM);
        CalibWeights[i] = 1;
    }
    e04xxc(&options);
    fail.print          = FALSE;
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    // options.optim_tol		= 0.0000001;

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = 3;
    params.Vec_marketvols = CalibPrices;
    params.Vec_strikes    = CalibStrikes;
    params.Vec_weights    = CalibWeights;
    params.Beta           = beta;
    params.Pi             = zeta;
    nagcomm.p             = (void*)&params;

    // trick to help calibrate
    // Move alpha down if needed
    i   = 0;
    err = BMM2GetStates(
        forward, maturity, sigma, alpha, beta, rho, zeta, &pFwd1, &pFwd2, &pSig1, &pSig2, &pPi);
    while (err && i < 15)
    {
        alpha = 0.9 * alpha;
        err   = BMM2GetStates(
            forward, maturity, sigma, alpha, beta, rho, zeta, &pFwd1, &pFwd2, &pSig1, &pSig2, &pPi);
        i++;
    }
    if (err)
        return err;
    // End of shitty trick

    x[0] = sqrt(sigma);  // set initial values for parameters
    x[1] = sqrt(alpha);
    x[2] = acos(rho);

    nag_opt_lsq_no_deriv(
        params.nb_strikes,
        3,
        BMM2lsqfunqnopi,
        x,
        &res,
        fvec,
        (double*)fjac,
        10,
        &options,
        &nagcomm,
        &fail);

    // write results
    *Calibsigma = x[0] * x[0];
    *Calibalpha = x[1] * x[1];
    *Calibbeta  = beta;
    *Calibrho   = cos(x[2]);
    *Calibzeta  = zeta;
    *Calibres   = sqrt(fvec[0] * fvec[0] + fvec[1] * fvec[1] + fvec[2] * fvec[2]) / 3;

FREE_RETURN:
    return err;
};

Err BMM3CalibOnSabr(
    double           forward,
    double           maturity,
    double           atmvol,
    double           alpha,
    double           beta,
    double           rho,
    double           veta,
    double           pi0,
    double           NSTD,
    SrtDiffusionType TypeInput,
    double*          Calibsigma,
    double*          Calibalpha,
    double*          Calibbeta,
    double*          Calibrho,
    double*          Calibveta,
    double*          Calibpi0,
    double*          Calibres)
{
    double          fjac[10][10], fvec[10], x[5];
    double          CalibPrices[5], CalibStrikes[5];
    double          vNStd[3] = {-1, 0, 1};
    double          res      = 0;
    int             i;
    double          sigma, sigmaLog, tempvol;
    double          pFwd, pFwd1, pFwd2, pSig, pSig1, pSig2, pPi, pPi0;
    Err             err = NULL;
    static NagError fail;
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer
    Opt_params      params;

    if (forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if (atmvol < 0.0)
    {
        err = "Error: Sigma should be positive";
        goto FREE_RETURN;
    }
    if (alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if ((rho < -1.0) || (rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // Compute SigmaBeta vol and atmlogvol
    srt_f_optsarbvol(
        forward, forward, maturity, atmvol, alpha, beta, rho, TypeInput, SRT_BETAVOL, &sigma);
    srt_f_optsarbvol(
        forward, forward, maturity, atmvol, alpha, beta, rho, TypeInput, SRT_LOGNORMAL, &sigmaLog);

    // compute strikes and prices:
    for (i = 0; i < 3; i++)
    {
        CalibStrikes[i] = forward * exp(vNStd[i] * NSTD * sqrt(maturity) * sigmaLog);
        srt_f_optsarbvol(
            forward,
            CalibStrikes[i],
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            SRT_BETAVOL,
            SRT_LOGNORMAL,
            &tempvol);
        CalibPrices[i] =
            srt_f_optblksch(forward, CalibStrikes[i], tempvol, maturity, 1, SRT_CALL, SRT_PREMIUM);
    }
    e04xxc(&options);
    fail.print          = FALSE;
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    options.optim_tol   = 0.0000001;

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = 3;
    params.Vec_marketvols = CalibPrices;
    params.Vec_strikes    = CalibStrikes;
    params.Beta           = beta;
    params.Veta           = veta;
    params.Pi0            = pi0;
    nagcomm.p             = (void*)&params;

    // trick to help calibrate
    // Move alpha down if needed
    i   = 0;
    err = BMM3GetStates(
        forward,
        maturity,
        sigma,
        alpha,
        beta,
        rho,
        veta,
        pi0,
        &pFwd,
        &pFwd1,
        &pFwd2,
        &pSig,
        &pSig1,
        &pSig2,
        &pPi,
        &pPi0);
    while (err && i < 15)
    {
        alpha = 0.9 * alpha;
        err   = BMM3GetStates(
            forward,
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            veta,
            pi0,
            &pFwd,
            &pFwd1,
            &pFwd2,
            &pSig,
            &pSig1,
            &pSig2,
            &pPi,
            &pPi0);
    }
    if (err)
        return err;
    // End of shitty trick

    x[0] = sqrt(sigma);  // set initial values for parameters
    x[1] = sqrt(alpha);
    x[2] = acos(rho);

    nag_opt_lsq_no_deriv(
        params.nb_strikes,
        3,
        BMM3lsqfunqnopi,
        x,
        &res,
        fvec,
        (double*)fjac,
        10,
        &options,
        &nagcomm,
        &fail);

    // write results
    *Calibsigma = x[0] * x[0];
    *Calibalpha = x[1] * x[1];
    *Calibbeta  = beta;
    *Calibrho   = cos(x[2]);
    *Calibveta  = veta;
    *Calibpi0   = pi0;
    *Calibres   = res;

FREE_RETURN:
    return err;
};

Err BMMGenCalibOnSabrWithStrikesAndFirstGuess(
    double           forward,
    double           maturity,
    double           atmvol,
    double           alpha,
    double           beta,
    double           a,
    double           b,
    double           c,
    double           rho,
    double           pi,
    SrtDiffusionType TypeVolLoc,
    int              nCalibStrikes,
    double*          CalibStrikes,
    double*          CalibVols,
    double*          CalibSigmabeta,
    double*          CalibAlpha,
    double*          CalibA,
    double*          CalibB,
    double*          CalibC,
    double*          CalibRho,
    double*          CalibPi,
    double*          Calibres)
{
    double          fjac[10][10], fvec[10], x[5];
    int             exit;
    double          res, pistep;
    int             i, tdj = 10;
    static NagError fail;
    Nag_Comm        nagcomm;
    Nag_E04_Opt     options;  // optional parameters of the NAG optimizer
    Opt_params      params;
    FuncVolLocType  vol_loc;

    if ((maturity == 0) || (pi == 0) || (nCalibStrikes == 0))
        return NULL;

    e04xxc(&options);
    fail.print          = FALSE;
    options.print_level = Nag_NoPrint;
    options.list        = FALSE;
    options.optim_tol   = 0.0000001;

    params.forward        = forward;
    params.maturity       = maturity;
    params.nb_strikes     = nCalibStrikes;
    params.Vec_marketvols = CalibVols;
    params.Vec_strikes    = CalibStrikes;
    params.a              = a;
    params.b              = b;
    params.c              = c;
    params.TypeVolLoc     = TypeVolLoc;
    params.Pi             = pi;

    nagcomm.p = (void*)&params;

    GetLocVolFromDiffusionType(TypeVolLoc, &vol_loc);
    x[0] = sqrt(forward * atmvol / vol_loc(forward, a, b, c, 0));
    x[1] = sqrt(alpha);
    x[2] = acos(rho);
    if (nCalibStrikes <= 3)
    {
        nag_opt_lsq_no_deriv(
            nCalibStrikes,
            3,
            BMMGenlsqfunqnopi,
            x,
            &res,
            fvec,
            (double*)fjac,
            tdj,
            &options,
            &nagcomm,
            &fail);
        // modified for BMMGen
        *CalibA   = a;
        *CalibB   = b;
        *CalibC   = c;
        *CalibPi  = pi;
        *Calibres = res;
    }
    else
    {
        exit   = 0;
        i      = 0;
        pistep = pi;
        while ((exit == 0) && (i < 10))
        {
            x[3] = sqrt(acos(pistep));
            nag_opt_lsq_no_deriv(
                nCalibStrikes,
                4,
                BMMlsqfunqwithpi,
                x,
                &res,
                fvec,
                (double*)fjac,
                tdj,
                &options,
                &nagcomm,
                &fail);
            if (fail.errnum == 0)
                exit = 1;
            else
            {
                x[0] = sqrt(forward * atmvol / vol_loc(forward, a, b, c, 0));
                x[1] = sqrt(alpha);
                x[2] = acos(rho);
                INIT_FAIL(fail);
                pistep = pistep * 0.75;
                i++;
            }
        }
        *CalibA   = a;
        *CalibB   = b;
        *CalibC   = c;
        *CalibPi  = cos(pow(x[3], 2));
        *Calibres = res;
    };

    *CalibSigmabeta = x[0] * x[0];
    *CalibAlpha     = x[1] * x[1];
    *CalibRho       = cos(x[2]);
    return NULL;
};

//****************************************************************************************************************
//**************************************** Ancillary Functions (NAG) *********************************************
//****************************************************************************************************************

static void NAG_CALL BMMlsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params* pparams        = (Opt_params*)(comm->p);
    double      forward        = pparams->forward;
    double      maturity       = pparams->maturity;
    double      beta           = pparams->Beta;
    double*     Vec_strikes    = pparams->Vec_strikes;
    double*     Vec_marketvols = pparams->Vec_marketvols;
    double*     Vec_weights    = pparams->Vec_weights;
    double      Pi             = pparams->Pi;
    double      Price;
    int         i;
    double      KeepOn = 0, res = 0;

    double sig = x[0] * x[0], alpha = x[1] * x[1], rho = cos(x[2]);

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5, min(rhosuralpha, 5));
    }

    if (pparams->ExitFlag)
    {
        comm->flag = -1;
    }
    else
    {
        for (i = 0; i < m; i++)
        {
            BMM_Option_Price(
                forward, Vec_strikes[i], maturity, sig, alpha, beta, rho, Pi, SRT_CALL, &Price);

            fvec[i] = (Price - Vec_marketvols[i]) * Vec_weights[i];
            KeepOn += (fabs(fvec[i]) > 0.0000001);
        }
        if (!KeepOn)
        {
            pparams->ExitFlag = 1;
        }
    }
}

static void NAG_CALL
BMM2lsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params* pparams        = (Opt_params*)(comm->p);
    double      forward        = pparams->forward;
    double      maturity       = pparams->maturity;
    double      beta           = pparams->Beta;
    double*     Vec_strikes    = pparams->Vec_strikes;
    double*     Vec_marketvols = pparams->Vec_marketvols;
    double*     Vec_weights    = pparams->Vec_weights;
    double      pi             = pparams->Pi;
    double      Price;
    int         i;

    double sigma = x[0] * x[0];
    double alpha = x[1] * x[1];
    double rho   = cos(x[2]);

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5.0, min(rhosuralpha, 5.0));
    }

    for (i = 0; i < m; i++)
    {
        BMM_Option_Price2(
            forward, Vec_strikes[i], maturity, sigma, alpha, beta, rho, pi, SRT_CALL, &Price);
        fvec[i] = (Price - Vec_marketvols[i]) * Vec_weights[i];
    }
}

static void NAG_CALL
BMM3lsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params* pparams        = (Opt_params*)(comm->p);
    double      forward        = pparams->forward;
    double      maturity       = pparams->maturity;
    double      beta           = pparams->Beta;
    double*     Vec_strikes    = pparams->Vec_strikes;
    double*     Vec_marketvols = pparams->Vec_marketvols;
    double      veta           = pparams->Veta;
    double      pi0            = pparams->Pi0;
    double      Price;
    int         i;

    double sigma = x[0] * x[0];
    double alpha = x[1] * x[1];
    double rho   = cos(x[2]);

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5, min(rhosuralpha, 5));
    }

    for (i = 0; i < m; i++)
    {
        BMM_Option_Price3(
            forward,
            Vec_strikes[i],
            maturity,
            sigma,
            alpha,
            beta,
            rho,
            veta,
            pi0,
            SRT_CALL,
            &Price);
        fvec[i] = (Price - Vec_marketvols[i]);
    }
}

static void NAG_CALL
BMMGenlsqfunqnopi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params*      pparams        = (Opt_params*)(comm->p);
    double           forward        = pparams->forward;
    double           maturity       = pparams->maturity;
    double           a              = pparams->a;
    double           b              = pparams->b;
    double           c              = pparams->c;
    double*          Vec_strikes    = pparams->Vec_strikes;
    double*          Vec_marketvols = pparams->Vec_marketvols;
    SrtDiffusionType TypeVolLoc     = pparams->TypeVolLoc;
    double           Pi             = pparams->Pi;
    double           Price;
    int              i;

    double sig = x[0] * x[0], alpha = x[1] * x[1], rho = cos(x[2]);

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5, min(rhosuralpha, 5));
    }

    for (i = 0; i < m; i++)
    {
        BMMGen_Option_Price(
            forward,
            Vec_strikes[i],
            maturity,
            sig,
            alpha,
            a,
            b,
            c,
            rho,
            Pi,
            SRT_CALL,
            &Price,
            TypeVolLoc);

        fvec[i] = (Price - Vec_marketvols[i]);
    }
}

static void NAG_CALL
BMMlsqfunqwithpi(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params* pparams        = (Opt_params*)(comm->p);
    double      forward        = pparams->forward;
    double      maturity       = pparams->maturity;
    double      beta           = pparams->Beta;
    double*     Vec_strikes    = pparams->Vec_strikes;
    double*     Vec_marketvols = pparams->Vec_marketvols;
    double*     Vec_weights    = pparams->Vec_weights;
    double      Price;
    int         i;

    double sig = x[0] * x[0], alpha = x[1] * x[1], rho = cos(x[2]), Pi = pow(cos(x[3]), 2);

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5, min(rhosuralpha, 5));
    }

    if (alpha > 10)
    {
        comm->flag = -1;
    }
    else
    {
        for (i = 0; i < m; i++)
        {
            BMM_Option_Price(
                forward, Vec_strikes[i], maturity, sig, alpha, beta, rho, Pi, SRT_CALL, &Price);

            fvec[i] = (Price - Vec_marketvols[i]) * Vec_weights[i];
        }
    }
}

static void NAG_CALL
BMM2lsqfunqwithzeta(Integer m, Integer n, double x[], double fvec[], Nag_Comm* comm)
{
    Opt_params* pparams        = (Opt_params*)(comm->p);
    double      forward        = pparams->forward;
    double      maturity       = pparams->maturity;
    double      beta           = pparams->Beta;
    double*     Vec_strikes    = pparams->Vec_strikes;
    double*     Vec_marketvols = pparams->Vec_marketvols;
    double*     Vec_weights    = pparams->Vec_weights;
    double      Price;
    int         i;

    double sig = x[0] * x[0], alpha = x[1] * x[1], rho = cos(x[2]), Zeta = cos(pow(x[3], 2));

    double rhosuralpha;
    if (alpha == 0)
        rho = 0;
    else
    {
        rhosuralpha = rho / alpha;
        rho         = alpha * max(-5, min(rhosuralpha, 5));
    }

    if (alpha > 10)
    {
        comm->flag = -1;
    }
    else
    {
        for (i = 0; i < m; i++)
        {
            BMM_Option_Price2(
                forward, Vec_strikes[i], maturity, sig, alpha, beta, rho, Zeta, SRT_CALL, &Price);

            fvec[i] = (Price - Vec_marketvols[i]) * Vec_weights[i];
        }
    }
}

double op_bmm_calib(
    double F,
    double K,
    double T,
    double sigma,  // ATMLOG vol
    double alpha,
    double beta,
    double rho,
    double pi)
{
    long   k;
    double sigma_beta1, sigma_beta2;
    double vol1, vol2, vega;
    double pricetemp;
    double error;

    sigma_beta1 = sigma * pow(F, 1.0 - beta);

    BMM_Option_Price(F, K, T, sigma_beta1, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol1);

    error = fabs(vol1 - sigma);

    /* First run */
    if (vol1 < sigma)
    {
        sigma_beta2 = sigma_beta1 * 1.005;
    }
    else
    {
        sigma_beta2 = sigma_beta1 * 0.995;
    }

    BMM_Option_Price(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    sigma_beta2 = sigma_beta1 + (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);

    BMM_Option_Price(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    error = fabs(vol2 - sigma);
    k     = 0;

    while (error > 0.000001 && k < 50)
    {
        vega        = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
        sigma_beta1 = sigma_beta2;
        vol1        = vol2;

        sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;

        BMM_Option_Price(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
        srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

        error = fabs(vol2 - sigma);

        k++;
    }
    return sigma_beta2;
}

double op_bmm2_calib(
    double F,
    double K,
    double T,
    double sigma,  // ATMLOG vol
    double alpha,
    double beta,
    double rho,
    double pi)
{
    long   k;
    double sigma_beta1, sigma_beta2;
    double vol1, vol2, vega;
    double pricetemp;
    double error, diffvol;
    Err    err = NULL;

    sigma_beta1 = sigma * pow(F, 1.0 - beta);

    BMM_Option_Price2(F, K, T, sigma_beta1, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol1);

    error = fabs(vol1 - sigma);
    if (error < 1e-6)
        return sigma_beta1;

    /* First run */
    if (vol1 < sigma)
    {
        sigma_beta2 = sigma_beta1 * 1.1;
    }
    else
    {
        sigma_beta2 = sigma_beta1 * 0.9;
    }
    BMM_Option_Price2(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    diffvol = fabs(vol1 - vol2);
    if (diffvol < 1e-6)
    {
        return sigma_beta1;
    }
    sigma_beta2 = sigma_beta1 + (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);

    BMM_Option_Price2(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    error = fabs(vol2 - sigma);
    k     = 0;

    while (error > 1e-6 && k < 50)
    {
        vega        = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
        sigma_beta1 = sigma_beta2;
        vol1        = vol2;

        sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;

        BMM_Option_Price2(F, K, T, sigma_beta2, alpha, beta, rho, pi, SRT_CALL, &pricetemp);
        srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);
        diffvol = fabs(vol1 - vol2);
        if (diffvol < 1e-6)
        {
            return sigma_beta1;
        }
        error = fabs(vol2 - sigma);

        k++;
    }
    return sigma_beta2;
}

double op_bmm3_calib(
    double F,
    double K,
    double T,
    double sigma,  // ATMLOG vol
    double alpha,
    double beta,
    double rho,
    double veta,
    double pi0)
{
    long   k;
    double sigma_beta1, sigma_beta2;
    double vol1, vol2, vega;
    double pricetemp;
    double error;

    sigma_beta1 = sigma * pow(F, 1.0 - beta);

    BMM_Option_Price3(F, K, T, sigma_beta1, alpha, beta, rho, veta, pi0, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol1);

    error = fabs(vol1 - sigma);

    /* First run */
    if (vol1 < sigma)
    {
        sigma_beta2 = sigma_beta1 * 1.005;
    }
    else
    {
        sigma_beta2 = sigma_beta1 * 0.995;
    }
    BMM_Option_Price3(F, K, T, sigma_beta2, alpha, beta, rho, veta, pi0, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    sigma_beta2 = sigma_beta1 + (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);

    BMM_Option_Price3(F, K, T, sigma_beta2, alpha, beta, rho, veta, pi0, SRT_CALL, &pricetemp);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    error = fabs(vol2 - sigma);
    k     = 0;

    while (error > 0.000001 && k < 50)
    {
        vega        = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
        sigma_beta1 = sigma_beta2;
        vol1        = vol2;

        sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;

        BMM_Option_Price3(F, K, T, sigma_beta2, alpha, beta, rho, veta, pi0, SRT_CALL, &pricetemp);
        srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

        error = fabs(vol2 - sigma);

        k++;
    }
    return sigma_beta2;
}

double op_bmmgen_calib(
    double           F,
    double           K,
    double           T,
    double           sigma,  // ATMLOG vol
    double           alpha,
    double           a,
    double           b,
    double           c,
    double           rho,
    double           pi,
    SrtDiffusionType TypeVolLoc)
{
    long           k;
    double         sigma_beta1, sigma_beta2;
    double         vol1, vol2, vega;
    double         pricetemp;
    double         error;
    FuncVolLocType vol_loc;
    GetLocVolFromDiffusionType(TypeVolLoc, &vol_loc);

    sigma_beta1 = sigma * F / vol_loc(F, a, b, c, 0);

    BMMGen_Option_Price(
        F, K, T, sigma_beta1, alpha, a, b, c, rho, pi, SRT_CALL, &pricetemp, TypeVolLoc);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol1);

    error = fabs(vol1 - sigma);

    /* First run */
    if (vol1 < sigma)
    {
        sigma_beta2 = sigma_beta1 * 1.005;
    }
    else
    {
        sigma_beta2 = sigma_beta1 * 0.995;
    }

    BMMGen_Option_Price(
        F, K, T, sigma_beta2, alpha, a, b, c, rho, pi, SRT_CALL, &pricetemp, TypeVolLoc);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    sigma_beta2 = sigma_beta1 + (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);

    BMMGen_Option_Price(
        F, K, T, sigma_beta2, alpha, a, b, c, rho, pi, SRT_CALL, &pricetemp, TypeVolLoc);
    srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

    error = fabs(vol2 - sigma);
    k     = 0;

    while (error > 0.000001 && k < 50)
    {
        vega        = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
        sigma_beta1 = sigma_beta2;
        vol1        = vol2;

        sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;

        BMMGen_Option_Price(
            F, K, T, sigma_beta2, alpha, a, b, c, rho, pi, SRT_CALL, &pricetemp, TypeVolLoc);
        srt_f_optimpvol(pricetemp, F, K, T, 1, SRT_CALL, SRT_LOGNORMAL, &vol2);

        error = fabs(vol2 - sigma);
        k++;
    }
    return sigma_beta2;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
//	srt_f_op_bmmvol(), Returns BSlog(@K) //copy from the sabr gen version..... should be
// improved
//
//	Returns vol(K) Given one ATM Vol and the smiles parameters
//	input SigmaBeta ATMLog, ATmNorm
//	output Log or Norm
//
/////////////////////////////////////////////////////////////////////////////////////////////

Err srt_f_optbmmvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           Pi,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    double res;
    double SigmaBeta;
    double pricetemp;
    Err    err = NULL;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((Pi < -1.0) || (Pi > 1.0))
    {
        err = "Error: Zeta should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // trick to handle ATS input vols.
    if ((input == SRT_LOGNORMAL_ATS) || (input == SRT_NORMAL_ATS))
    {
        if (input == SRT_LOGNORMAL_ATS)
        {
            pricetemp = srt_f_optblksch(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        if (input == SRT_NORMAL_ATS)
        {
            pricetemp = srt_f_optblknrm(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        err = srt_f_optimpvol(
            pricetemp, Forward, Strike, Maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &Sigma);
        input = SRT_LOGNORMAL;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta

        SigmaBeta = op_bmm_calib(Forward, Forward, Maturity, Sigma, Alpha, Beta, Rho, Pi);

    else  // Case input = Beta
        SigmaBeta = Sigma;

    if (output == SRT_BETAVOL)
        res = SigmaBeta;
    else
    {  // Convert to LOG
        BMM_Option_Price(
            Forward, Strike, Maturity, SigmaBeta, Alpha, Beta, Rho, Pi, SRT_CALL, &pricetemp);
        srt_f_optimpvol(pricetemp, Forward, Strike, Maturity, 1, SRT_CALL, SRT_LOGNORMAL, &res);

        if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
        {
            err = srt_f_optsarbvol(
                Forward, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
            if (err)
                goto FREE_RETURN;
        }
    }

    *vol = res;

FREE_RETURN:

    return err;
}

Err srt_f_optbmmprice(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           Pi,
    SrtDiffusionType input,
    SrtCallPutType   call_put,
    double*          price)
{
    double SigmaBeta, pricetemp;
    Err    err = NULL;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((Pi < -1.0) || (Pi > 1.0))
    {
        err = "Error: Zeta should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // trick to handle ATS input vols.
    if ((input == SRT_LOGNORMAL_ATS) || (input == SRT_NORMAL_ATS))
    {
        if (input == SRT_LOGNORMAL_ATS)
        {
            pricetemp = srt_f_optblksch(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        if (input == SRT_NORMAL_ATS)
        {
            pricetemp = srt_f_optblknrm(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        err = srt_f_optimpvol(
            pricetemp, Forward, Strike, Maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &Sigma);
        input = SRT_LOGNORMAL;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta

        SigmaBeta = op_bmm2_calib(Forward, Forward, Maturity, Sigma, Alpha, Beta, Rho, Pi);

    else  // Case input = Beta
        SigmaBeta = Sigma;

    BMM_Option_Price2(
        Forward, Strike, Maturity, SigmaBeta, Alpha, Beta, Rho, Pi, call_put, &pricetemp);
    *price = pricetemp;

FREE_RETURN:

    return err;
}

Err srt_f_optbmmcalvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           Pi,
    double           NStd,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    Err    err = NULL;
    double CalibSigma, CalibBeta, CalibAlpha, CalibRho, CalibPi, CalibRes;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((Pi <= 0) || (Pi >= 1.0))
    {
        err = "Error: Pi should be > 0 and < 1";
        goto FREE_RETURN;
    }
    if (NStd <= 0)
    {
        err = "Error: NStd should be > 0=";
        goto FREE_RETURN;
    }

    err = BMMCalibOnSabr(
        Forward,
        Maturity,
        Sigma,
        Alpha,
        Beta,
        Rho,
        Pi,
        NStd,
        input,
        &CalibSigma,
        &CalibAlpha,
        &CalibBeta,
        &CalibRho,
        &CalibPi,
        &CalibRes);

    if (!err)
        err = srt_f_optbmmvol(
            Forward,
            Strike,
            Maturity,
            CalibSigma,
            CalibAlpha,
            CalibBeta,
            CalibRho,
            CalibPi,
            SRT_BETAVOL,
            output,
            vol);

FREE_RETURN:
    return err;
}

Err srt_f_optbmmvolfromstates(
    double           Strike,
    double           Maturity,
    double           Beta,
    double           Fwd1,
    double           Fwd2,
    double           Sig1,
    double           Sig2,
    double           Pi,
    SrtDiffusionType output,
    double*          vol)
{
    double res;
    double Fwd;
    Err    err = NULL;

    // Control of passed arguments
    if ((Fwd1 < 0.0) || (Fwd2 < 0.0))
    {
        err = "Error: Forwards should be positive";
        goto FREE_RETURN;
    }
    if ((Sig1 < 0.0) || (Sig2 < 0.0))
    {
        err = "Error: Sigmas should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Pi < 0.0) || (Pi > 1.0))
    {
        err = "Error: Pi should be > 0 and < 1";
        goto FREE_RETURN;
    }
    if ((output != SRT_NORMAL) && (output != SRT_LOGNORMAL))
    {
        err = "output should be NORMAL or LOGNORMAL";
        goto FREE_RETURN;
    }

    Fwd = Pi * Fwd1 + (1 - Pi) * Fwd2;
    BMM_Option_Price_From_States(
        Strike, Maturity, Beta, Fwd1, Fwd2, Sig1, Sig2, Pi, SRT_CALL, &res);
    srt_f_optimpvol(res, Fwd, Strike, Maturity, 1, SRT_CALL, output, &res);

    if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
    {
        err = srt_f_optsarbvol(Fwd, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
        if (err)
            goto FREE_RETURN;
    }

    *vol = res;

FREE_RETURN:

    return err;
}

Err srt_f_optbmm2vol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           Zeta,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    double res;
    double SigmaBeta;
    double pricetemp;
    Err    err = NULL;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((Zeta < 0.0))
    {
        err = "Error: Zeta should be > 0";
        goto FREE_RETURN;
    }

    // trick to handle ATS input vols.
    if ((input == SRT_LOGNORMAL_ATS) || (input == SRT_NORMAL_ATS))
    {
        if (input == SRT_LOGNORMAL_ATS)
        {
            pricetemp = srt_f_optblksch(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        if (input == SRT_NORMAL_ATS)
        {
            pricetemp = srt_f_optblknrm(Forward, Strike, Sigma, Maturity, 1.0, SRT_CALL, PREMIUM);
        }

        err = srt_f_optimpvol(
            pricetemp, Forward, Strike, Maturity, 1.0, SRT_CALL, SRT_LOGNORMAL, &Sigma);
        input = SRT_LOGNORMAL;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta

        SigmaBeta = op_bmm2_calib(Forward, Forward, Maturity, Sigma, Alpha, Beta, Rho, Zeta);

    else  // Case input = Beta
        SigmaBeta = Sigma;

    if (output == SRT_BETAVOL)
        res = SigmaBeta;
    else
    {  // Convert to LOG
        BMM_Option_Price2(
            Forward, Strike, Maturity, SigmaBeta, Alpha, Beta, Rho, Zeta, SRT_CALL, &pricetemp);
        srt_f_optimpvol(pricetemp, Forward, Strike, Maturity, 1, SRT_CALL, SRT_LOGNORMAL, &res);

        if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
        {
            err = srt_f_optsarbvol(
                Forward, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
            if (err)
                goto FREE_RETURN;
        }
    }

    *vol = res;

FREE_RETURN:

    return err;
}

Err srt_f_optbmm2calvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           Zeta,
    double           NStd,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    Err    err = NULL;
    double CalibSigma, CalibBeta, CalibAlpha, CalibRho, CalibZeta, CalibRes;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((Zeta < 0))
    {
        err = "Error: Zeta should be >= 0";
        goto FREE_RETURN;
    }
    if (NStd <= 0)
    {
        err = "Error: NStd should be > 0";
        goto FREE_RETURN;
    }

    err = BMM2CalibOnSabr(
        Forward,
        Maturity,
        Sigma,
        Alpha,
        Beta,
        Rho,
        Zeta,
        NStd,
        input,
        &CalibSigma,
        &CalibAlpha,
        &CalibBeta,
        &CalibRho,
        &CalibZeta,
        &CalibRes);
    if (!err)
        err = srt_f_optbmm2vol(
            Forward,
            Strike,
            Maturity,
            CalibSigma,
            CalibAlpha,
            CalibBeta,
            CalibRho,
            CalibZeta,
            input,
            output,
            vol);

FREE_RETURN:
    return err;
}

Err srt_f_optbmm3vol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           veta,
    double           pi0,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    double res;
    double SigmaBeta;
    double pricetemp;
    Err    err = NULL;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta

        SigmaBeta = op_bmm3_calib(Forward, Forward, Maturity, Sigma, Alpha, Beta, Rho, veta, pi0);

    else  // Case input = Beta
        SigmaBeta = Sigma;

    if (output == SRT_BETAVOL)
        res = SigmaBeta;
    else
    {  // Convert to LOG
        BMM_Option_Price3(
            Forward,
            Strike,
            Maturity,
            SigmaBeta,
            Alpha,
            Beta,
            Rho,
            veta,
            pi0,
            SRT_CALL,
            &pricetemp);
        srt_f_optimpvol(pricetemp, Forward, Strike, Maturity, 1, SRT_CALL, SRT_LOGNORMAL, &res);

        if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
        {
            err = srt_f_optsarbvol(
                Forward, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
            if (err)
                goto FREE_RETURN;
        }
    }

    *vol = res;

FREE_RETURN:

    return err;
}

Err srt_f_optbmm3calvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           Beta,
    double           Rho,
    double           veta,
    double           pi0,
    double           NStd,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol)
{
    Err    err = NULL;
    double CalibSigma, CalibBeta, CalibAlpha, CalibRho, CalibVeta, CalibPi0, CalibRes;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }
    if ((pi0 <= 0) || (pi0 >= 1.0))
    {
        err = "Error: Pi0 should be > 0 and < 1";
        goto FREE_RETURN;
    }
    if ((veta <= 0) || (veta >= 1.0))
    {
        err = "Error: Veta should be > 0 and < 1";
        goto FREE_RETURN;
    }
    if (NStd <= 0)
    {
        err = "Error: NStd should be > 0=";
        goto FREE_RETURN;
    }

    err = BMM3CalibOnSabr(
        Forward,
        Maturity,
        Sigma,
        Alpha,
        Beta,
        Rho,
        veta,
        pi0,
        NStd,
        input,
        &CalibSigma,
        &CalibAlpha,
        &CalibBeta,
        &CalibRho,
        &CalibVeta,
        &CalibPi0,
        &CalibRes);
    if (!err)
        err = srt_f_optbmm3vol(
            Forward,
            Strike,
            Maturity,
            CalibSigma,
            CalibAlpha,
            CalibBeta,
            CalibRho,
            CalibVeta,
            CalibPi0,
            input,
            output,
            vol);

FREE_RETURN:
    return err;
}

Err srt_f_optbmm3volfromstates(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Beta,
    double           Fwd,
    double           Fwd1,
    double           Fwd2,
    double           Sig,
    double           Sig1,
    double           Sig2,
    double           Pi0,
    double           Pi,
    SrtDiffusionType output,
    double*          vol)
{
    double res;
    Err    err = NULL;

    BMM_Option_Price3_From_States(
        Forward, Strike, Maturity, Beta, Fwd, Fwd1, Fwd2, Sig, Sig1, Sig2, Pi0, Pi, SRT_CALL, &res);
    srt_f_optimpvol(res, Fwd, Strike, Maturity, 1, SRT_CALL, output, &res);

    if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
    {
        err = srt_f_optsarbvol(Fwd, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
        if (err)
            goto FREE_RETURN;
    }
    *vol = res;

FREE_RETURN:

    return err;
}

Err srt_f_optbmmgenvol(
    double           Forward,
    double           Strike,
    double           Maturity,
    double           Sigma,
    double           Alpha,
    double           a,
    double           b,
    double           c,
    double           Rho,
    double           Pi,
    SrtDiffusionType input,
    SrtDiffusionType output,
    double*          vol,
    SrtDiffusionType TypeVolLoc)
{
    double res;
    double SigmaBeta;
    double pricetemp;
    Err    err = NULL;

    // Control of passed arguments
    if (Forward < 0.0)
    {
        err = "Error: Forward should be positive";
        goto FREE_RETURN;
    }
    if (Sigma < 0.0)
    {
        err = "Error: SigmaBeta should be positive";
        goto FREE_RETURN;
    }
    if (Alpha < 0.0)
    {
        err = "Error: Alpha should be positive";
        goto FREE_RETURN;
    }
    if (Maturity < 0.0)
    {
        err = "Error: Maturity should better be positive";
        goto FREE_RETURN;
    }
    if ((Rho < -1.0) || (Rho > 1.0))
    {
        err = "Error: Rho should be > -1 and < 1";
        goto FREE_RETURN;
    }

    // Process According to imput
    if (input == SRT_NORMAL)  // first convert to LOG
    {
        err = srt_f_optsarbvol(
            Forward, Forward, Maturity, Sigma, 0, 0, 0, SRT_BETAVOL, SRT_LOGNORMAL, &Sigma);
        if (err)
            goto FREE_RETURN;
    }
    if (input != SRT_BETAVOL)  // From LOG to SigmaBeta

        SigmaBeta =
            op_bmmgen_calib(Forward, Forward, Maturity, Sigma, Alpha, a, b, c, Rho, Pi, TypeVolLoc);

    else  // Case input = Beta
        SigmaBeta = Sigma;

    if (output == SRT_BETAVOL)
        res = SigmaBeta;
    else
    {  // Convert to LOG
        BMMGen_Option_Price(
            Forward,
            Strike,
            Maturity,
            SigmaBeta,
            Alpha,
            a,
            b,
            c,
            Rho,
            Pi,
            SRT_CALL,
            &pricetemp,
            TypeVolLoc);
        srt_f_optimpvol(pricetemp, Forward, Strike, Maturity, 1, SRT_CALL, SRT_LOGNORMAL, &res);

        if (output == SRT_NORMAL)  // if output = NORM, Convert to LOG
        {
            err = srt_f_optsarbvol(
                Forward, Strike, Maturity, res, 0, 1, 0, SRT_BETAVOL, SRT_NORMAL, &res);
            if (err)
                goto FREE_RETURN;
        }
    }

    *vol = res;

FREE_RETURN:

    return err;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
// Pat hagan formulae for Beta and BMMGEN tools (high order formula)
// Also implementation of NAG Bessel to price optiions under a Beta model
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

// Beta implementation using NAG does not bring anything to the story
Err BetaPrice2(
    double  Forward,
    double  Strike,
    double  Maturity,
    double  SigmaBeta,
    double  Beta,
    double  Rho,
    double* Price)
{
    NagError Fail;

    double Sigma02T = pow(SigmaBeta, 2.0) * Maturity;
    double Bessel1, Bessel2;

    // Bessel = s18eec(Strike/pow(Forward,Theta)/Sigma02T,2.0/abs(Theta),&Fail);
    // res =
    // abs(Theta)/Sigma02T*sqrt(Forward)/pow(Strike,1.5+2*Theta)*exp(-(pow(Strike,-2*Theta)+pow(Forward,-2*Theta))/(2*Sigma02T))*Bessel;

    double a = pow(Strike, 2 * (1 - Beta)) / (Sigma02T * pow(1 - Beta, 2.0));
    double c = pow(Forward, 2 * (1 - Beta)) / (Sigma02T * pow(1 - Beta, 2.0));
    double b = 1 / (1 - Beta);

    INIT_FAIL(Fail);

    Bessel1 = g01gcc(a, 2 + b, c, 1e-15, 500, &Fail);
    Bessel2 = g01gcc(c, b, a, 1e-15, 500, &Fail);

    *Price = Forward * (1 - Bessel1) - Strike * Bessel2;

    return NULL;
}

// Beta implementation using high order Pat, see BMMGenPrice() for generic version
Err BetaPrice3(
    double Forward, double Strike, double SigmaBeta, double Maturity, double Beta, double* Price)
{
    FuncVolLocType vol_loc = vol_sabr;
    double         Fav     = 0.5 * (Forward + Strike);
    double         Theta   = Forward - Strike;
    double         VolFav  = vol_loc(Fav, Beta, 0, 0, 0);
    double         Delta   = SigmaBeta * SigmaBeta * VolFav * VolFav * Maturity;
    double         Gam1    = vol_loc(Fav, Beta, 0, 0, 1) / VolFav;
    double         Gam2    = vol_loc(Fav, Beta, 0, 0, 2) / VolFav;
    double         Gam3    = Gam2 / Fav * (Beta - 2);
    double         Gam4    = Gam3 / Fav * (Beta - 3);

    double vol;
    vol = SigmaBeta * VolFav / Fav *
          (1 + Delta / 24 * (2 * Gam2 - 1 * Gam1 * Gam1 + 1 / pow(Fav, 2.0)) +
           pow(Theta, 2.0) / 24 * (1 * Gam2 - 2 * Gam1 * Gam1 + 2 / pow(Fav, 2.0))

           + pow(Delta, 2.0) / 480 *
                 (2 * Gam4 + 4 * Gam1 * Gam3 + 3 * Gam2 * Gam2 - 3 * Gam1 * Gam1 * Gam2 +
                  0.75 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                  0.5 / pow(Fav, 2.0) * (10 * Gam2 - 5 * Gam1 * Gam1 + 5 / pow(Fav, 2.0))) +
           Delta * pow(Theta, 2.0) / 2880 *
               (6 * Gam4 - 18 * Gam1 * Gam3 + 14 * Gam2 * Gam2 - 29 * Gam1 * Gam1 * Gam2 +
                11 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                1 / pow(Fav, 2.0) * (35 * Gam2 - 40 * Gam1 * Gam1 + 40 / pow(Fav, 2.0))) +
           pow(Theta, 4.0) / 1440 *
               (0.75 * Gam4 - 6 * Gam1 * Gam3 - 2 * Gam2 * Gam2 + 17 * Gam1 * Gam1 * Gam2 -
                8 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                1 / pow(Fav, 2.0) * (5 * Gam2 - 10 * Gam1 * Gam1 + 10 / pow(Fav, 2.0))));

    *Price = srt_f_optblksch(Forward, Strike, vol, Maturity, 1.0, SRT_CALL, PREMIUM);

    return NULL;
}

// Pat Hagan high order formula for generic local vol -> same trick as SABRGEN)
Err BMMGenPrice(
    double           Forward,
    double           Strike,
    double           SigmaBeta,
    double           Maturity,
    double           a,
    double           b,
    double           c,
    double*          Price,
    SrtDiffusionType TypeVolLoc)
{
    double         Fav, Theta, VolFav, Delta, Gam1, Gam2, Gam3, Gam4, vol;
    FuncVolLocType vol_loc;
    GetLocVolFromDiffusionType(TypeVolLoc, &vol_loc);

    Fav    = 0.5 * (Forward + Strike);
    Theta  = Forward - Strike;
    VolFav = vol_loc(Fav, a, b, c, 0);
    Delta  = SigmaBeta * SigmaBeta * VolFav * VolFav * Maturity;
    Gam1   = vol_loc(Fav, a, b, c, 1) / VolFav;
    Gam2   = vol_loc(Fav, a, b, c, 2) / VolFav;
    Gam3   = Gam2 / Fav * (a - 2);
    Gam4   = Gam3 / Fav * (a - 3);

    vol = SigmaBeta * VolFav / Fav *
          (1 + Delta / 24 * (2 * Gam2 - 1 * Gam1 * Gam1 + 1 / pow(Fav, 2.0)) +
           pow(Theta, 2.0) / 24 * (1 * Gam2 - 2 * Gam1 * Gam1 + 2 / pow(Fav, 2.0))

           + pow(Delta, 2.0) / 480 *
                 (2 * Gam4 + 4 * Gam1 * Gam3 + 3 * Gam2 * Gam2 - 3 * Gam1 * Gam1 * Gam2 +
                  0.75 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                  0.5 / pow(Fav, 2.0) * (10 * Gam2 - 5 * Gam1 * Gam1 + 5 / pow(Fav, 2.0))) +
           Delta * pow(Theta, 2.0) / 2880 *
               (6 * Gam4 - 18 * Gam1 * Gam3 + 14 * Gam2 * Gam2 - 29 * Gam1 * Gam1 * Gam2 +
                11 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                1 / pow(Fav, 2.0) * (35 * Gam2 - 40 * Gam1 * Gam1 + 40 / pow(Fav, 2.0))) +
           pow(Theta, 4.0) / 1440 *
               (0.75 * Gam4 - 6 * Gam1 * Gam3 - 2 * Gam2 * Gam2 + 17 * Gam1 * Gam1 * Gam2 -
                8 * (pow(Gam1, 4.0) - 1 / pow(Fav, 4.0)) +
                1 / pow(Fav, 2.0) * (5 * Gam2 - 10 * Gam1 * Gam1 + 10 / pow(Fav, 2.0))));
    *Price = srt_f_optblksch(Forward, Strike, vol, Maturity, 1.0, SRT_CALL, PREMIUM);

    return NULL;
}
