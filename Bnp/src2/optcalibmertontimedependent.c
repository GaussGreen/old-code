/* Calibrates with Levenberg the parameters of a simple Merton model*/

#include "utallhdr.h"

#include "math.h"
#include "num_h_allhdr.h"
#include "num_h_levenberg.h"
#include "opfnctns.h"
#include "srt_h_all.h"
#include "srt_h_resetable.h"

/* Free Memory */
void free_levenberg_memory2(double** covar, double** alpha, int nparam);

#define n_opt 5     /* number of parameters in the model */
#define shift 0.001 /* shift used to compute gradients */

/* Global Variables declaration */
static double*  strikes = NULL;
static double*  dvols   = NULL;
static double   dFwd;
static int      n_strikes, nperiods;
static double*  vols     = NULL;
static double** gradient = NULL;
static double*  mat;
static char*    Lognorm;

/* Computes the gradient of the price */
Err priceforcalibtimedep(double str, double param[], double* value, double* Deriv, int n_param)
{
    int     i, j, k;
    double* value_shiftee;
    double* param_shiftes;
    int     out_range = 0;
    Err     err       = NULL;

    if (str == strikes[1])
    {
        if (strcmp(Lognorm, "Normal") == 0)
        {
            if ((param[1] < 0) || (param[1] > 0.1) || (param[2] < 0) || (param[2] > 0.02) ||
                (param[3] < 0) || (param[3] > 10) || (param[4] > 0) || (param[4] < -0.02) ||
                (param[5] < 0) || (param[5] > 10))
                out_range = 1;

            if (out_range == 1)
            {
                param[1] = 0.01;
                param[2] = 0.0025;
                param[3] = 1;
                param[4] = -0.0025;
                param[5] = 1;
            }
        }
        else
        {
            if ((param[1] < 0) || (param[1] > 0.5) || (param[2] < 0) || (param[2] > 1) ||
                (param[3] < 0) || (param[3] > 5) || (param[4] > 0) || (param[4] < -1) ||
                (param[5] < 0) || (param[5] > 5))
                out_range = 1;

            if (out_range == 1)
            {
                param[1] = 0.1;
                param[2] = 0.15;
                param[3] = 1;
                param[4] = -0.15;
                param[5] = 1;
            }
        }
        /* Computes the implied vols for all strikes */
        err = optmertonsmiletimedependent(
            dFwd, n_strikes, strikes, nperiods, param, mat, Lognorm, vols /*implied vols*/
        );
        if (err)
        {
            smessage("Error in optmertonsmiletimedependent");
            return err;
        }

        /* Calculates the derivatives: optimization parameters */
        param_shiftes = dvector(1, n_param);
        value_shiftee = dvector(1, n_strikes);
        for (j = 1; j < (n_param + 1); ++j)
        {
            for (k = 1; k < (n_param + 1); k++)
                param_shiftes[k] = param[k];
            param_shiftes[j] += param[j] * shift;

            err = optmertonsmiletimedependent(
                dFwd,
                n_strikes,
                strikes,
                nperiods,
                param_shiftes,
                mat,
                Lognorm,
                value_shiftee /*resulting implied vols*/
            );
            if (err)
            {
                smessage("Error in optmertonsmiletimedependent");
                free_dvector(param_shiftes, 1, n_param);
                free_dvector(value_shiftee, 1, n_strikes);
                return err;
            }

            for (i = 1; i < (n_strikes + 1); ++i)
                gradient[i][j] = (value_shiftee[i] - vols[i]) / (shift * param[j]);
        }

        out_range = 0;

        if (strcmp(Lognorm, "Normal") == 0)
        {
            if ((param[1] < 0) || (param[1] > 0.1) || (param[2] < 0) || (param[2] > 0.02) ||
                (param[3] < 0) || (param[3] > 10) || (param[4] > 0) || (param[4] < -0.02) ||
                (param[5] < 0) || (param[5] > 10))
                out_range = 1;
        }

        else
        {
            if ((param[1] < 0) || (param[1] > 0.5) || (param[2] < 0) || (param[2] > 1) ||
                (param[3] < 0) || (param[3] > 5) || (param[4] > 0) || (param[4] < -1) ||
                (param[5] < 0) || (param[5] > 5))
                out_range = 1;
        }
        if (out_range == 1)
            for (i = 1; i < (n_strikes + 1); ++i)
            {
                vols[i] = -10000;
                for (j = 1; j < (n_param + 1); ++j)
                    gradient[i][j] = 10000;
            }
        free_dvector(param_shiftes, 1, n_param);
        free_dvector(value_shiftee, 1, n_strikes);
    }

    /* Returns the elements of the vector after each calls of the function*/
    for (i = 1; i < (n_strikes + 1); i++)
    {
        if (str == strikes[i])
        {
            *value = vols[i];
            for (j = 1; j < (n_param + 1); j++)
                Deriv[j] = gradient[i][j];
            break;
        }
    }
    return err;
}

/* Functions which finds the parameters which fit the market vols*/
Err optcalibmertontimedep(
    double  dFwd_loc,
    int     n_strikes_loc,
    int     n_periods_loc,
    double* strikes_loc,
    double* market_vols,
    double* param,
    double* maturity,
    char*   Logornorm,
    double* chisq,
    long*   ia,
    double* calibvols)
{
    double* sig;
    double  Strike_inf, Strike_sup;
    int     i, i1, i2;
    int     nparam;
    double  totalmat;
    Err     err = NULL;

    totalmat  = 0;
    Lognorm   = strdup(Logornorm);
    dFwd      = dFwd_loc;
    n_strikes = n_strikes_loc;
    nperiods  = n_periods_loc;
    nparam    = 5 * nperiods;
    strikes   = dvector(1, n_strikes);
    dvols     = dvector(1, n_strikes);
    vols      = dvector(1, n_strikes);
    gradient  = dmatrix(1, n_strikes, 1, nparam);
    mat       = dvector(1, nperiods);
    for (i = 1; i < (n_strikes + 1); i++)
        strikes[i] = strikes_loc[i];
    for (i = 1; i < (n_strikes + 1); i++)
        dvols[i] = market_vols[i];
    for (i = 1; i <= nperiods; i++)
    {
        mat[i] = maturity[i];
        totalmat += maturity[i];
    }
    sig = dvector(1, n_strikes);

    /* Assignment of weights to the different strikes. More weight is given to strikes
    closer to the ATM. Also, basically no weight is given to strikes further away than
    2 sd's from the money.
    */
    i = 1;

    while (strikes[i] < dFwd)
        i++;

    if (strcmp(Lognorm, "Normal") == 0)
    {
        Strike_inf = dFwd - 2 * market_vols[i] * sqrt(totalmat);
        Strike_sup = dFwd + 2 * market_vols[i] * sqrt(totalmat);
    }
    else
    {
        Strike_inf = dFwd * (1 - 2 * market_vols[i] * sqrt(totalmat));
        Strike_sup = dFwd * (1 + 2 * market_vols[i] * sqrt(totalmat));
    }
    i1 = 1;
    i2 = n_strikes;

    for (i = 1; i < n_strikes - 1; i++)
    {
        if ((strikes[i] < Strike_inf) && (strikes[i + 1] > Strike_inf))
            i1 = i;
        if ((strikes[i] < Strike_sup) && (strikes[i + 1] > Strike_sup))
            i2 = i;
    }

    for (i = 1; i < i1; i++)
        sig[i] = 100;
    for (i = i1; i < i2 + 1; i++)
        sig[i] = 10 * fabs(strikes[i] - dFwd) + 0.0001;
    for (i = i2 + 1; i < n_strikes + 1; i++)
        sig[i] = 100;

    err = levenberg_marquardt_select(
        strikes, dvols, sig, n_strikes_loc, param, ia, nparam, 10, priceforcalibtimedep, chisq);
    if (err)
    {
        smessage("Error in levenberg_marquardt_select");
    }
    else
    {
        err = optmertonsmiletimedependent(
            dFwd, n_strikes, strikes, nperiods, param, mat, Logornorm, calibvols);
        if (err)
            smessage("Error in optmertonsmile");
    }

    free_dvector(strikes, 1, n_strikes);
    free_dvector(dvols, 1, n_strikes);
    free_dvector(sig, 1, n_strikes);
    free_dvector(vols, 1, n_strikes);
    free_dmatrix(gradient, 1, n_strikes, 1, nparam);
    free_dvector(mat, 1, nperiods);
    free(Lognorm);
    return err;
}
