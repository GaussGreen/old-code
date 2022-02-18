/* -------------------------------------------------------------------------------------------------------------------------------------------------
 */
/* USDStochShift.c */
/* Contains functions:

*/

/* Include files */
#include "USDStochShift.h"

#include "math.h"
#include "opFnctns.h"

/* --------------------------------------------------------------------------------------------------------------------

Price a call option using the stochastic shift model.  The process is

  d(F+L) = sigma dW
  dL = alpha dZ
  dW*dZ = rho dt

The option is a put option with strike K and expiry tau.  The models are
iModel = 0:		double integration

---------------------------------------------------------------------------------------------------------------------
*/

static double _StochShift_h(
    double x, double rho, double f, double l, double alpha, double sigma, double t, double K)
{
    return rho * x +
           ((f + l) * exp(-.5 * sigma * sigma * t + sigma * sqrt(t) * x) - l - K) / alpha / sqrt(t);
}

char* _StochShift_Exact(
    double                   in_dF,
    double                   in_dL,
    double                   in_dalpha,
    double                   in_dsigma,
    double                   in_drho,
    double                   in_dK,
    double                   in_dTau,
    struct StochShift_Param* in_Param,
    double*                  out_dPrice)
{
    /* Variable declaration */
    int     i, N;
    double  dSum = 0, h;
    double *x = 0, *w = 0;
    char*   err = 0;
    /* Variable initialization and memory allocation */
    N           = in_Param->iIntPoints;
    x           = calloc(N + 1, sizeof(double));
    w           = calloc(N + 1, sizeof(double));
    *out_dPrice = -1.0;

    /* Generate the Gaussian points */
    if (err = HermiteStandard(x, w, N))
        goto FREE_RETURN;

    /* Perform the integration */
    dSum = 0.0;
    for (i = 1; i < N; i++)
    {
        h = _StochShift_h(x[i], in_drho, in_dF, in_dL, in_dalpha, in_dsigma, in_dTau, in_dK);

        dSum += w[i] * (h * norm(h) + gauss(h));
    }
    *out_dPrice = in_dalpha * sqrt(in_dTau) * dSum;

/* Free the memory */
FREE_RETURN:
    free_and_zero(&x);
    free_and_zero(&w);

    /* Return */
    return err;
}

/* This function expands the price as a Taylor series in sigma */
char* _StochShift_Approx(
    double                   in_dF,
    double                   in_dL,
    double                   in_dalpha,
    double                   in_dsigma,
    double                   in_drho,
    double                   in_dK,
    double                   in_dTau,
    struct StochShift_Param* in_Param,
    double*                  out_dPrice)
{
    /* Variable Declaration */
    double normStd, normStd2, normStd4, h0, Nh0, nh0, fpl, fpl2, intrin, intrin2, logStd, logStd2;
    double dTerm0, dTerm2, dTerm4;

    /* Variable Initialization */
    normStd  = in_dalpha * sqrt(in_dTau);
    normStd2 = normStd * normStd;
    normStd4 = normStd2 * normStd2;
    h0       = (in_dF - in_dK) / normStd;
    Nh0      = norm(h0);
    nh0      = gauss(h0);
    fpl      = in_dF + in_dL;
    fpl2     = fpl * fpl;
    intrin   = in_dF - in_dK;
    intrin2  = intrin * intrin;
    logStd   = in_dsigma * sqrt(in_dTau);
    logStd2  = logStd * logStd;

    /* Terms in Taylor series */
    dTerm0 = normStd * (h0 * Nh0 + nh0);
    dTerm2 = 0.5 * fpl2 * logStd2 * nh0 / normStd;
    dTerm4 = 0.125 * logStd2 * logStd2 / normStd4 / normStd * fpl2 * nh0 *
             (fpl2 * intrin2 - 4.0 * fpl * intrin * normStd2 - fpl2 * normStd2 + 2 * normStd4);

    /* Value */
    *out_dPrice = dTerm0 + dTerm2 + dTerm4;

    /* return */
    return 0;
}

/* triage */
char* StochShift_PriceCall(
    double                   in_dF,
    double                   in_dL,
    double                   in_dalpha,
    double                   in_dsigma,
    double                   in_drho,
    double                   in_dK,
    double                   in_dTau,
    struct StochShift_Param* in_Param,
    double*                  out_dPrice)
{
    /* Check that the inputs are sensible */

    /* Call the correct model */
    switch (in_Param->iModel)
    {
    case SS_EXACT:
        return _StochShift_Exact(
            in_dF, in_dL, in_dalpha, in_dsigma, in_drho, in_dK, in_dTau, in_Param, out_dPrice);
    case SS_APPROX:
        return _StochShift_Approx(
            in_dF, in_dL, in_dalpha, in_dsigma, in_drho, in_dK, in_dTau, in_Param, out_dPrice);
    }
    return "Unknown model";
}
