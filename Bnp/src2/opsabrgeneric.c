/* ==========================================================================
   FILE_NAME:	opsabrgeneric.c

   PURPOSE:		MODIFIED IMPLEMENTATION OF SABR (generic local vol)

   DATE:		03/03/04

   AUTHOR:		P.T.
   ========================================================================== */

#include "opsabrgeneric.h"

#include "opfnctns.h"
//#include "Fx3FUtils.h"
//#include "FxSabrGrfn.h"
#include "math.h"
#include "opsabrcalib.h"

#define DEFAULT_PERC 0.005
#define LIMIT_DOWN_P 0.0005

//
// op_sabrgen()			Returns BSlog(@K)
//
// op_sabrgen_calib()	Given the ATMLOGVol, solves for the SigmaBeta vol
//
// vol_bvm()				Specify your local vol of the form vol_mylocvol(double,
// double,double,double, int)
//

/////////////////////////////////////////////////////////////////////////////////////////////
//
// op_sabrgen(), Returns BSlog(@K)
//
/////////////////////////////////////////////////////////////////////////////////////////////

double op_sabrgen(
    double F,
    double K,
    double T,
    double sigma,  // SigmaBeta Vol
    double alpha,
    double a,
    double b,
    double c,
    double rho,
    double (*vol_local)(double x, double a, double b, double c, int type))
{
    double Af, AfP1, AfP2, AfIF, AfIK;
    double Fb, gam1, gam2, Z, XZ;
    double res;

    if (T > 1.0e-14)
    {
        Fb = 0.5 * (F + K);  // 0.5 * (F + K);  sqrt(F * K);  old version

        Af   = vol_local(Fb, a, b, c, 0);
        AfP1 = vol_local(Fb, a, b, c, 1);
        AfP2 = vol_local(Fb, a, b, c, 2);
        AfIF = vol_local(F, a, b, c, 3);
        AfIK = vol_local(K, a, b, c, 3);

        gam1 = AfP1 / Af;
        gam2 = AfP2 / Af;

        //	Z = alpha / sigma * (AfIF - AfIK); // old version seems to be worse.

        Z  = alpha / sigma * (F - K) / Af;
        XZ = log((sqrt(1.0 - 2.0 * rho * Z + Z * Z) + Z - rho) / (1.0 - rho));

        if (fabs(F - K) < 1.0E-10)
        {
            res = sigma * Af / F;
        }
        else
        {
            double explosive_term = log(F / K) / (AfIF - AfIK) * Z / XZ;
            res                   = sigma * explosive_term;
        }

        res *=
            (1.0 + ((2.0 * gam2 - gam1 * gam1 + 1.0 / Fb / Fb) * Fb * Fb / 24.0 +
                    alpha * Fb / sigma / Af *
                        (0.25 * rho * gam1 * Fb +
                         (2.0 - 3.0 * rho * rho) / 24.0 * alpha * Fb / sigma / Af)) *
                       sigma * sigma * Af * Af / Fb / Fb * T);
    }
    else
    {
        res = 0.1;
    }

    return res;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
// Given the ATMLOGVol, solves for the SigmaBeta vol
//
/////////////////////////////////////////////////////////////////////////////////////////////

double op_sabrgen_calib(
    double F,
    double K,
    double T,
    double sigma,  // ATMLOG vol
    double alpha,
    double a,
    double b,
    double c,
    double rho,
    double (*vol_local)(double x, double a, double b, double c, int type))
{
    long   k;
    double sigma_beta1, sigma_beta2;
    double vol1, vol2, vega;
    double error;

    sigma_beta1 = sigma * F / vol_local(F, a, b, c, 0);
    vol1        = op_sabrgen(F, K, T, sigma_beta1, alpha, a, b, c, rho, vol_local);

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

    vol2 = op_sabrgen(F, K, T, sigma_beta2, alpha, a, b, c, rho, vol_local);

    sigma_beta2 = sigma_beta1 + (sigma - vol1) * (sigma_beta2 - sigma_beta1) / (vol2 - vol1);
    vol2        = op_sabrgen(F, K, T, sigma_beta2, alpha, a, b, c, rho, vol_local);

    error = fabs(vol2 - sigma);
    k     = 0;

    while (error > 0.000001 && k < 50)
    {
        vega        = (vol2 - vol1) / (sigma_beta2 - sigma_beta1);
        sigma_beta1 = sigma_beta2;
        vol1        = vol2;

        sigma_beta2 = sigma_beta1 + (sigma - vol1) / vega;
        vol2        = op_sabrgen(F, K, T, sigma_beta2, alpha, a, b, c, rho, vol_local);

        error = fabs(vol2 - sigma);

        k++;
    }
    //	smessage (" nb loops op_sabrgen_calib : %d", k);
    return sigma_beta2;
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
// Returns a pointer to your favorite local vol //add yours here
//
//
/////////////////////////////////////////////////////////////////////////////////////////////
void GetLocVolFromDiffusionType(SrtDiffusionType TypeVolLoc, FuncVolLocType* vol_loc)
{
    switch (TypeVolLoc)
    {
    case SRT_SABRVOL:
        *vol_loc = vol_sabr;
        break;
    case SRT_BVMVOL:
        *vol_loc = vol_bvm;
        break;
    case SRT_BVM2VOL:
        *vol_loc = vol_bvm2;
        break;
    case SRT_BVMHVOL:
        *vol_loc = vol_bvmh;
        break;
    case SRT_BVMH2VOL:
        *vol_loc = vol_bvmh2;
        break;
    case SRT_BVMCVOL:
        *vol_loc = vol_bvmc;
        break;
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////
//
// The loc vol functions
//
// x: forward
// a: param1
// b: param2
// c: param3
// type: specify function(0), derivative(1), sec derivative(2), prime of inverse(3), eq Beta (4)
//
/////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: SABR, Vol = F^beta
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_sabr(double x, double beta, double not_used1, double not_used2, int type)
{
    static double res;

    /* Vol = F^beta */

    switch (type)
    {
    case 0:
    {
        res = pow(x, beta);

        return res;
        break;
    }
    case 1:
    {
        res = beta * pow(x, beta - 1.0);

        return res;
        break;
    }
    case 2:
    {
        res = beta * (beta - 1.0) * pow(x, beta - 2.0);

        return res;
        break;
    }
    case 3:
    {
        if (beta == 1)
            res = log(x);
        else
            res = 1.0 / (1.0 - beta) * pow(x, 1.0 - beta);

        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_shifted_log, Vol = a * F + b
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_shifted_log(double x, double a, double b, double c, int type)
{
    /* Vol = a * F + b */
    switch (type)
    {
    case 0:
    {
        return (a * x + b);
        break;
    }
    case 1:
    {
        return (a);
        break;
    }
    case 2:
    {
        return 0;
        break;
    }
    case 3:
    {
        return log(a * x + b) / a;
        break;
    }
    }

    return 0.0;
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_log_quadra, Vol = F * (a * lnF ^ 2 + b * lnF + c)
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_log_quadra(double x, double a, double b, double c, int type)
{
    static double lnx, res, x1, x2, delta;

    /* Vol = F * (a * lnF ^ 2 + b * lnF + c) */

    switch (type)
    {
    case 0:
    {
        lnx = log(x);
        res = x * (lnx * (a * lnx + b) + c);

        return res;
        break;
    }
    case 1:
    {
        lnx = log(x);
        res = (lnx * (a * lnx + b) + c);
        res += 2.0 * a * lnx + b;

        return res;
        break;
    }
    case 2:
    {
        lnx = log(x);
        res = (2.0 * a * lnx + b + 2.0 * a) / x;

        return res;
        break;
    }
    case 3:
    {
        lnx = log(x);

        if (fabs(a) < 1.0E-10)
        {
            if (fabs(b) < 1.0E-10)
            {
                res = lnx / c;
            }
            else
            {
                res = log(fabs(b * lnx + c)) / b;
            }
        }
        else
        {
            delta = b * b - 4.0 * a * c;

            if (delta > 1.0E-10)
            {
                delta = sqrt(delta);

                x1 = (-b + delta) / (2.0 * a);
                x2 = (-b - delta) / (2.0 * a);

                res = 1.0 / (x1 - x2) * log(fabs((lnx - x1) / (lnx - x2))) / a;
            }
            else if (delta < -1.0E-10)
            {
                delta = sqrt(-delta);

                res = (2.0 * a * lnx + b) / delta;
                res = 2.0 / delta * atan(res);
            }
            else
            {
                res = -1.0 / (a * lnx + b);
            }
        }

        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_quadra, Vol = a * F ^ 2 + b * F + c
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_quadra(double x, double a, double b, double c, int type)
{
    static double res, x1, x2, delta;

    /* Vol = a * F ^ 2 + b * F + c */

    switch (type)
    {
    case 0:
    {
        res = x * (a * x + b) + c;

        return res;
        break;
    }
    case 1:
    {
        res = 2.0 * a * x + b;

        return res;
        break;
    }
    case 2:
    {
        res = 2.0 * a;

        return res;
        break;
    }
    case 3:
    {
        if (fabs(a) < 1.0E-10)
        {
            if (fabs(b) < 1.0E-10)
            {
                res = x / c;
            }
            else
            {
                res = log(fabs(b * x + c)) / b;
            }
        }
        else
        {
            delta = b * b - 4.0 * a * c;

            if (delta > 1.0E-10)
            {
                delta = sqrt(delta);

                x1 = (-b + delta) / (2.0 * a);
                x2 = (-b - delta) / (2.0 * a);

                res = 1.0 / (x1 - x2) * log(fabs((x - x1) / (x - x2))) / a;
            }
            else if (delta < -1.0E-10)
            {
                delta = sqrt(-delta);

                res = (2.0 * a * x + b) / delta;
                res = 2.0 / delta * atan(res);
            }
            else
            {
                res = -1.0 / (a * x + b);
            }
        }

        return res;
        break;
    }
    default:
    {
        return 0.0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_bvm, Vol = (1 - exp(-x / delay))
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_bvm(double x, double delay, double not_used1, double not_used2, int type)
{
    static double res, x1, x2, delta;

    /* Vol = (1 - exp(-x / delay)) */

    switch (type)
    {
    case 0:
    {
        res = (1 - exp(-x / delay));

        return res;
        break;
    }
    case 1:
    {
        res = exp(-x / delay) / delay;

        return res;
        break;
    }
    case 2:
    {
        res = -exp(-x / delay) / delay / delay;

        return res;
        break;
    }
    case 3:
    {
        res = x + log(1 - exp(-x / delay)) * delay;

        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_bvm2,  Vol = x^beta (1 - exp((-x / x0)^(1-beta)))
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_bvm2(double F, double F0, double Beta, double not_used2, int type)
{
    double res;

    /* Vol = x^Beta (1 - exp((-x / x0)^(1-Beta))) */

    F = F / F0;

    if (Beta == 0)
    {
        return vol_bvm(F, F0, 0, 0, type);
    }
    switch (type)
    {
    case 0:
    {
        res = pow(F, Beta) * (1 - exp(-pow(F, 1.0 - Beta)));

        return res;
        break;
    }
    case 1:
    {
        double o2 = 1. - Beta;
        double o3 = pow(F, o2);
        res       = (exp(-o3) * (F - Beta * F + Beta * (-1. + exp(o3)) * pow(F, Beta))) / F;
        res /= F0;

        return res;
        break;
    }
    case 2:
    {
        double o1 = -1. + Beta;
        double o3 = 1. - Beta;
        double o4 = pow(F, o3);
        res       = o1 * exp(-o4) *
              (-(o1 * pow(F, 2.)) + Beta * (-1. + exp(o4)) * pow(F, 2. * Beta) -
               Beta * pow(F, 1. + Beta)) *
              pow(F, -2. - Beta);
        res /= F0 * F0;
        return res;
        break;
    }
    case 3:
    {
        double o2 = 1. - Beta;
        res       = log(-1. + exp(pow(F, o2))) / o2;
        res *= F0;

        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_bvmh
//
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_bvmh(double F, double F0, double Beta, double notused2, int type)
{
    static double res;

    // Vol = F^Beta*tanh(F^(1-Beta))
    F = F / F0;

    switch (type)
    {
    case 0:
    {
        res = pow(F, Beta) * tanh(pow(F, 1.0 - Beta));

        return res;
        break;
    }
    case 1:
    {
        double o1 = -1. + Beta;
        double o2 = -Beta;
        double o3 = 1. + o2;
        double o4 = pow(F, o3);
        res       = -(o1 * pow(1 / cosh(o4), 2.)) + Beta * pow(F, o1) * tanh(o4);
        res /= F0;

        return res;
        break;
    }
    case 2:
    {
        double o1 = -Beta;
        double o2 = -1. + Beta;
        double o3 = 1. + o1;
        double o4 = pow(F, o3);
        double o5 = tanh(o4);
        res       = o2 * pow(F, -2. + o1) *
              (o5 * Beta * pow(F, 2. * Beta) -
               F * (2. * F * o2 * o5 + Beta * pow(F, Beta)) * pow(1 / cosh(o4), 2.));
        res /= F0 * F0;

        return res;
        break;
    }
    case 3:
    {
        res = log(sinh(pow(F, 1.0 - Beta))) / (1.0 - Beta);
        res *= F0;

        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_bvmh2
//
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_bvmh2(double F, double F0, double F1, double Beta, int type)
{
    static double res;
    int           i;

    // Vol = x/delay1 * tanh[x/delay1] / tanh[x/delay2]

    switch (type)
    {
    case 0:
    {
        res = pow(F / F1, Beta) * tanh(F / F0) / tanh(pow(F / F1, Beta));

        return res;
        break;
    }
    case 1:
    {
        double o4 = pow(F / F1, Beta);
        double o6 = tanh(F / F0);
        res       = (1 / F0 * o4 *
               (-(F0 * o4 * o6 * Beta * pow(1 / sinh(o4), 2.)) +
                1 / tanh(o4) * (F0 * o6 * Beta + F * pow(1 / cosh(F / F0), 2.)))) /
              F;

        return res;
        break;
    }
    case 2:
    {
        double o1  = pow(F0, 2.);
        double o2  = 1 / F1;
        double o3  = F * o2;
        double o4  = pow(o3, Beta);
        double o5  = 1 / sinh(o4);
        double o6  = pow(o5, 2.);
        double o7  = 1 / F0;
        double o8  = F * o7;
        double o9  = 1 / sinh(o8);
        double o10 = pow(o9, 2.);
        double o11 = tanh(o8);
        res        = (o4 * pow(F, -2.) *
               (-(F0 * o4 * o6 * Beta * (2. * F * o10 + F0 * o11 * (-1. + 3. * Beta))) +
                1 / tanh(o4) *
                    (2. * F * o10 * (-(F * o11) + F0 * Beta) +
                     o1 * o11 * Beta * (-1. + Beta + 2. * o6 * Beta * pow(o3, 2. * Beta))))) /
              o1;

        return res;
        break;
    }
    case 3:
    {
        int    nbstep  = 100;
        double tempvar = 0.5 * (F0 + F1);
        double mesh    = (F - tempvar) / nbstep;
        double tempres = 0;
        for (i = 0; i < nbstep; i++)
        {
            // tempres += mesh / (tempvar / F0 * tanh(tempvar / F0) / tanh(tempvar / F1));
            tempres += mesh / vol_bvmh2(tempvar, F0, F1, Beta, 0);
            tempvar += mesh;
        }
        res = tempres;
        return res;
        break;
    }
    default:
    {
        return 0;
    }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// local vol function: vol_bvmc,  Vol = (F/F0)^(betamax+betamin) / ((F/F0)^betamax + (F/F0)^betamin)
//
/////////////////////////////////////////////////////////////////////////////////////////////
double vol_bvmc(double F, double F0, double BetaMin, double BetaMax, int type)
{
    double res;

    /* Vol = (F/F0)^(betamax+betamin) / ((F/F0)^betamax + (F/F0)^betamin)   */

    F = F / F0;  // Rescale by zero and do not forget to adjust the derivatives

    // reorder BetaMin and BetaMax so that the if's  are fine

    if (BetaMin > BetaMax)
    {
        double temp = BetaMin;
        BetaMin     = BetaMax;
        BetaMax     = temp;
    }

    if (BetaMin == BetaMax && BetaMax == 0)
    {
        switch (type)
        {
        case 0:
        {
            res = 1;
            return res;
            break;
        }
        case 1:
        {
            res = 0;
            return res;
            break;
        }
        case 2:
        {
            res = 0;
            return res;
            break;
        }
        case 3:
        {
            res = F;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMin == BetaMax && BetaMax == 1)
    {
        switch (type)
        {
        case 0:
        {
            res = F;
            return res;
            break;
        }
        case 1:
        {
            res = 1;
            res = res / F0;
            return res;
            break;
        }
        case 2:
        {
            res = 0;
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            res = log(F);
            res = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMin == BetaMax)
    {
        switch (type)
        {
        case 0:
        {
            res = pow(F, BetaMin);
            return res;
            break;
        }
        case 1:
        {
            res = BetaMin * pow(F, BetaMin - 1.0);
            res = res / F0;
            return res;
            break;
        }
        case 2:
        {
            res = BetaMin * (BetaMin - 1.0) * pow(F, BetaMin - 2.0);
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            res = pow(F, 1.0 - BetaMin) / (1.0 - BetaMin);
            res = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMin == 0 && BetaMax == 1)
    {
        switch (type)
        {
        case 0:
        {
            res = F / (1. + F);
            return res;
            break;
        }
        case 1:
        {
            res = pow(1. + F, -2.);
            res = res / F0;
            return res;
            break;
        }
        case 2:
        {
            res = -2. * pow(1. + F, -3.);
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            res = F + log(F);
            res = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMin == 0)
    {
        switch (type)
        {
        case 0:
        {
            double o2 = pow(F, BetaMax);
            res       = o2 / (1. + o2);
            return res;
            break;
        }
        case 1:
        {
            res = BetaMax * pow(F, -1. + BetaMax) * pow(1. + pow(F, BetaMax), -2.);
            res = res / F0;
            return res;
            break;
        }
        case 2:
        {
            double o1 = pow(F, BetaMax);
            res =
                -(BetaMax * (1. - BetaMax + (1. + BetaMax) * o1) * pow(F, -2. + BetaMax) *
                  pow(1. + o1, -3.));
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            double o2 = 1. - BetaMax;
            res       = F + pow(F, o2) / o2;
            res       = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMin == 1)
    {
        switch (type)
        {
        case 0:
        {
            res = pow(F, 1. + BetaMax) / (F + pow(F, BetaMax));
            return res;
            break;
        }
        case 1:
        {
            double o1 = pow(F, BetaMax);
            res       = o1 * (BetaMax * F + o1) * pow(F + o1, -2.);
            res       = res / F0;
            return res;
            break;
        }
        case 2:
        {
            double o1 = pow(F, BetaMax);
            res       = -(
                (-1. + BetaMax) * o1 * (-(BetaMax * F) + (-2. + BetaMax) * o1) * pow(F + o1, -3.));
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            double o2 = 1. - BetaMax;
            res       = log(F) + pow(F, o2) / o2;
            res       = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    else if (BetaMax == 1)
    {
        switch (type)
        {
        case 0:
        {
            res = pow(F, 1. + BetaMin) / (F + pow(F, BetaMin));
            return res;
            break;
        }
        case 1:
        {
            double o1 = pow(F, BetaMin);
            res       = o1 * (BetaMin * F + o1) * pow(F + o1, -2.);
            res       = res / F0;
            return res;
            break;
        }
        case 2:
        {
            double o1 = pow(F, BetaMin);
            res       = -(
                (-1. + BetaMin) * o1 * (-(BetaMin * F) + (-2. + BetaMin) * o1) * pow(F + o1, -3.));
            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            double o2 = 1. - BetaMin;
            res       = log(F) + pow(F, o2) / o2;
            res       = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
    // Main Case
    else
    {
        switch (type)
        {
        case 0:
        {
            res = pow(F, BetaMax) / (1. + pow(F, -BetaMin + BetaMax));

            return res;
            break;
        }
        case 1:
        {
            double o1 = pow(F, BetaMax);
            double o2 = pow(F, BetaMin);
            res =
                (BetaMin * o1 + BetaMax * o2) * pow(F, -1. + BetaMax + BetaMin) * pow(o1 + o2, -2.);

            res = res / F0;
            return res;
            break;
        }
        case 2:
        {
            res = pow(F, -2. + BetaMax + BetaMin) *
                  ((-1. + BetaMin) * BetaMin * pow(F, 2. * BetaMax) +
                   (-1. + BetaMax) * BetaMax * pow(F, 2. * BetaMin) -
                   (BetaMax + BetaMin - 4. * BetaMax * BetaMin + pow(BetaMax, 2.) +
                    pow(BetaMin, 2.)) *
                       pow(F, BetaMax + BetaMin)) *
                  pow(pow(F, BetaMax) + pow(F, BetaMin), -3.);

            res = res / (F0 * F0);
            return res;
            break;
        }
        case 3:
        {
            double o1 = 1. - BetaMax;
            double o2 = 1. - BetaMin;
            res       = pow(F, o1) / o1 + pow(F, o2) / o2;

            res = F0 * res;
            return res;
            break;
        }
        default:
        {
            return 0;
        }
        }
    }
}
/////////////////////////////////////////////////////////////////////////////////////////////
//
// op_sabr_distri,  returns the density or cumulative distribution at strike K
//
/////////////////////////////////////////////////////////////////////////////////////////////

double op_sabr_distri(
    double F,
    double K,
    double T,
    double sigma,
    double alpha,
    double a,
    double b,
    double c,
    double rho,
    double perc,
    int    type, /* 0: density, 1: cumul */
    double (*vol_local)(double x, double a, double b, double c, int type))
{
    double price1, price2, price3;
    double eps, res;

    if (fabs(perc) < 1.0E-10)
    {
        perc = DEFAULT_PERC;
    }

    eps = F * perc;

    if (K - eps < 0)
    {
        eps = K / 2.0;
    }

    price1 = op_sabrgen(F, K - eps, T, sigma, alpha, a, b, c, rho, vol_local);
    price1 = srt_f_optblksch(F, K - eps, price1, T, 1.0, SRT_CALL, PREMIUM);

    price3 = op_sabrgen(F, K + eps, T, sigma, alpha, a, b, c, rho, vol_local);
    price3 = srt_f_optblksch(F, K + eps, price3, T, 1.0, SRT_CALL, PREMIUM);

    if (type)
    {
        res = 1.0 - (price1 - price3) / (2.0 * eps);
    }
    else
    {
        price2 = op_sabrgen(F, K, T, sigma, alpha, a, b, c, rho, vol_local);
        price2 = srt_f_optblksch(F, K, price2, T, 1.0, SRT_CALL, PREMIUM);

        res = (price1 + price3 - 2.0 * price2) / (eps * eps);
    }

    return res;
}

// checking of the boundaries for a specified vol_loc
// type
// lowerbounds : 10->a, 11->b, 12->c
// upperbounds : 20->a, 21->b, 22->c
// check all : 30
double vol_parameters_check_from_SrtDiffusionType(
    double a, double b, double c, int type, SrtDiffusionType VolLocType)
{
    double alb = -DBL_MAX, aub = DBL_MAX, blb = -DBL_MAX, bub = DBL_MAX, clb = -DBL_MAX,
           cub = DBL_MAX;

    switch (VolLocType)
    {
    case SRT_SABRVOL:
    {
        alb = 0.0;
        aub = 1;
        break;
    }
    case SRT_BVMVOL:
    {
        alb = DBL_EPSILON;
        break;
    }
    case SRT_BVM2VOL:
    {
        alb = DBL_EPSILON;
        blb = 0;
        bub = 1;
        break;
    }
    case SRT_BVMHVOL:
    {
        alb = DBL_EPSILON;
        blb = 0;
        bub = 1;
        break;
    }
    case SRT_BVMH2VOL:
    {
        alb = DBL_EPSILON;
        blb = DBL_EPSILON;
        clb = 0;
        cub = 1;
        break;
    }
    case SRT_BVMCVOL:
    {
        alb = DBL_EPSILON;
        blb = 0;
        bub = 1;
        clb = 0;
        cub = 10;
        break;
    }
    default:
    {
    }
    }

    switch (type)
    {
    case 10:
        return alb;
    case 11:
        return blb;
    case 12:
        return clb;
    case 20:
        return aub;
    case 21:
        return bub;
    case 22:
        return cub;
    case 30:
    {
        if ((a < aub) && (a > alb) && (b < bub) && (b > blb) && (c < cub) && (c > clb))
            return 1;  // ok
        else
            return 0;  // out of bounds
    }
    default:
    {
        return 0;
    }
    }
}
