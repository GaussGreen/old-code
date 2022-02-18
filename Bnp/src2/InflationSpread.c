#include "InflationSpread.h"

#include "math.h"
#include "nag.h"
#include "nagd01.h"
#include "srt_h_ts.h"
#include "srt_h_ts_irm.h"
#include "srt_h_types.h"
#include "srt_h_und_struct.h"
#include "srtaccess.h"
#include "utallhdr.h"

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Local Structures */
struct sInflationArg
{
    double dCorrelation;
    //	Date dtExpiry;
    //	DateList dtMaturity;
    double  dRateMean;
    double  dRateStdDev;
    int     iNumDates;
    double *dvPsiDiff, *dvPsiDiff2, dF, dG;
    double* dvForwardDF;
    double  dExpiryDF;
    //	double dMaturityDF;
    //	char* csUndName;
};

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Local Function Declarations */
Err srt_InflationSpreadOption(
    double dRate, double dInflation, struct sInflationArg* argPtr, double* pdOption);

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/*                       */
/* Local Function Bodies */
/*                       */
/* ----------------------------------------------------------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* The integrand */
double NAG_CALL inflationIntegrand(Integer iDim, double* x, Nag_User* comm)
{
    // Variable declarations
    double                dDensity, dRate, dSpread, dDF, dOption;
    struct sInflationArg* argPtr       = comm->p;
    double                dCorrelation = argPtr->dCorrelation;

    // Given the N(0,1) variates calculate the short rate and the spread
    dRate   = argPtr->dRateMean + argPtr->dRateStdDev * x[0];
    dSpread = 0.0;

    // From the realized value of the spread, calculate the forward spreads

    // From the value of the short rate, calculate the required discount bonds
    // At the same time compute the forward
    srt_InflationSpreadOption(dRate, dSpread, argPtr, &dDF);
    dOption = dDF;

    // Calculate the uncorrelated 2D Gaussian
    dDensity = 0.5 * exp(-0.5 * (x[0] * x[0] + x[1] * x[1])) / SRT_PI;

    // Return the probability weighted value of the option
    return dOption * dDensity;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Price the option */
Err srt_InflationSpreadOption(
    double dRate, double dInflation, struct sInflationArg* argPtr, double* pdOption)
{
    int i;
    /* Calculate the level */
    *pdOption = 0.0;
    for (i = 0; i < argPtr->iNumDates; i++)
        *pdOption += argPtr->dvForwardDF[i] * exp(argPtr->dvPsiDiff[i] * dRate / argPtr->dF -
                                                  0.5 * argPtr->dvPsiDiff2[i] * argPtr->dG);

    return 0;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Set up the integrand structure */
Err setInflationStructure(
    char* csUndName, Date dtExpiry, Date* dtvMaturity, int iNumDates, struct sInflationArg* sArg)
{
    /* Variable declarations */
    TermStruct* ptsShortRate;
    SrtUndPtr   srtUndPtr;
    char*       csYCName;
    double      dExpiry, dMaturity, dH, dPsiExpiry;
    Err         err;
    Date        dtToday;
    int         i;

    /* get the und from the underlying name */
    srtUndPtr = lookup_und(csUndName);

    /* get the underlying term structure */
    if (err = get_underlying_ts(srtUndPtr, &ptsShortRate))
        return err;

    /* Get the name of the discount curve associated to the underlying */
    csYCName = get_discname_from_underlying(srtUndPtr);

    /* get today from the underlying */
    dtToday = get_today_from_underlying(srtUndPtr);

    /* Calculate the Expiry dependent functions */
    dExpiry         = (dtExpiry - dtToday) * YEARS_IN_DAY;
    sArg->dExpiryDF = swp_f_df(dtToday, dtExpiry, csYCName);
    dPsiExpiry      = Psi_func(dExpiry, ptsShortRate);
    sArg->dF        = F_func(dExpiry, ptsShortRate);
    G_H_func(dExpiry, ptsShortRate, &sArg->dG, &dH);
    sArg->dRateMean   = 0.0; /* swp_f_zr( dtExpiry-7.0, dtExpiry+7.0, csYCName ); */
    sArg->dRateStdDev = sArg->dF * sqrt(sArg->dG);

    /* Allocate the memory for the Maturity dependent variables */
    sArg->iNumDates   = iNumDates;
    sArg->dvPsiDiff   = dvector(0, sArg->iNumDates - 1);
    sArg->dvPsiDiff2  = dvector(0, sArg->iNumDates - 1);
    sArg->dvForwardDF = dvector(0, sArg->iNumDates - 1);

    /* Calculate the Maturity dependent variables */
    for (i = 0; i < sArg->iNumDates; i++)
    {
        dMaturity            = (dtvMaturity[i] - dtToday) * YEARS_IN_DAY;
        sArg->dvForwardDF[i] = swp_f_df(dtExpiry, dtvMaturity[i], csYCName);
        sArg->dvPsiDiff[i]   = dPsiExpiry - Psi_func(dMaturity, ptsShortRate);
        sArg->dvPsiDiff2[i]  = sArg->dvPsiDiff[i] * sArg->dvPsiDiff[i];
    }

    /* Other arguments */
    sArg->dCorrelation = 0.0;

    /* Return */
    return 0;
}

void deleteInflationStructure(struct sInflationArg* sArg)
{
    free_dvector(sArg->dvForwardDF, 0, sArg->iNumDates - 1);
    free_dvector(sArg->dvPsiDiff, 0, sArg->iNumDates - 1);
    free_dvector(sArg->dvPsiDiff2, 0, sArg->iNumDates - 1);
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/*                          */
/* Exported Function Bodies */
/*                          */
/* ----------------------------------------------------------------------------------------------------------------------------------
 */

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Prices a trigger option on an inflation spread swap */
Err srt_InflationSpread(
    char*   csUndName,
    Date    dtExpiry,
    Date*   dtvMaturity,
    int     iNumDates,
    double* dpIntArg,
    double* dpPrice)
{
    // Variable declarations
    double               dEps, dAcc;
    double *             dvLowerLimits, *dvUpperLimits;
    Nag_User             nagArg;
    struct sInflationArg sArg;
    Integer              iMinPts, iMaxPts, iDim;
    static NagError      neFail;
    Err                  err;

    /* Calibrate the mean reversion so that the forwards are correctly reproduced */

    /* Set up the structure for the integrand */
    nagArg.p = (Pointer)&sArg;
    if (err = setInflationStructure(csUndName, dtExpiry, dtvMaturity, iNumDates, &sArg))
        return err;

    /* Set the upper and lower limits, the number of points and the desired accuracy (or use default
     * values) */
    iDim          = 2;
    dvLowerLimits = dvector(0, 1);
    dvUpperLimits = dvector(0, 1);
    if (dpIntArg)
    {
        dvLowerLimits[0] = dpIntArg[0];
        dvLowerLimits[1] = dpIntArg[1];
        dvUpperLimits[0] = dpIntArg[2];
        dvUpperLimits[1] = dpIntArg[3];
        iMinPts          = (int)dpIntArg[4];
        iMaxPts          = (int)dpIntArg[5];
        dEps             = dpIntArg[6];
    }
    else
    {
        dvLowerLimits[0] = -10.0;
        dvLowerLimits[1] = -10.0;
        dvUpperLimits[0] = 10.0;
        dvUpperLimits[1] = 10.0;
        iMinPts          = 100;
        iMaxPts          = 1000000;
        dEps             = 1.0e-06;
    }

    /* Integrate over a 2D Gaussian */
    nag_multid_quad_adapt_1(
        iDim,
        &inflationIntegrand,
        dvLowerLimits,
        dvUpperLimits,
        &iMinPts,
        iMaxPts,
        dEps,
        dpPrice,
        &dAcc,
        &nagArg,
        &neFail);
    *dpPrice *= sArg.dExpiryDF;

    if (neFail.code != NE_NOERROR)
        err = neFail.message;

    /* Free the memory */
    free_dvector(dvLowerLimits, 0, 1);
    free_dvector(dvUpperLimits, 0, 1);
    deleteInflationStructure(&sArg);

    /* Return */
    return err;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Calculates the inflation forward */
Err srt_f_get_vasicek_var_sr(double time, TermStruct* ts, double* var_sr)
{
    Err    err = NULL;
    double G, H;

    G_H_func(time, ts, &G, &H);

    (*var_sr) = F_func(time, ts) * F_func(time, ts) * G;

    return NULL;
}

Err srt_get_vasicek_risk_neutral_sr_mean(double time, TermStruct* ts, double* vasicek_mean_sr)
{
    SrtLst*           ls;
    IrmTermStructVal *tsval, *tsval_p;
    Err               err = NULL;

    ls = ts->head;

    while ((ls != NULL) && (((IrmTermStructVal*)ls->element->val.pval)->time < time))
        ls = ls->next;

    if (ls == NULL)
        ls = ts->tail;

    tsval = (IrmTermStructVal*)ls->element->val.pval;

    if (ls != ts->head)
    {
        tsval_p = (IrmTermStructVal*)ls->previous->element->val.pval;
        time -= tsval_p->time;
    }

    (*vasicek_mean_sr) =
        tsval->vasicek_mean_sr + (tsval->mean_rev_level / tsval->F) * (exp(time / tsval->tau) - 1);

    return err;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Calculates the mean of a Vasicek underlying in a forward measure defined by an LGM */
Err srt_get_vasicek_forward_measure_und_mean(
    double      dTime,
    double      dFwdTime,
    double      dCorrelation,
    TermStruct* ptsVasicek,
    TermStruct* ptsLGM,
    double*     pdMean)
{
    /* Variable declarations */
    SrtLst *          plsVasicek, *plsLGM;
    IrmTermStructVal *pvalVasicek, *pvalLGM;
    Err               err;
    double            dIntegral, dPsiFwdMeasure, dUpperTime, dLowerTime, dTimeDiff, lambda, a;

    /* The mean is equal to the mean in the risk neutral measure, plus a correction.  Begin with the
     * risk-neutral mean */
    if (err = srt_get_vasicek_risk_neutral_sr_mean(dTime, ptsVasicek, pdMean))
        return err;

    /* Point at the start of the lists */
    plsVasicek = ptsVasicek->head;
    plsLGM     = ptsLGM->head;

    /* Move the loops until they are pointing to the element beyond the time */
    while ((plsVasicek != NULL) &&
           (((IrmTermStructVal*)plsVasicek->element->val.pval)->time < dTime))
        plsVasicek = plsVasicek->next;
    if (plsVasicek == NULL)
        plsVasicek = ptsVasicek->tail;

    while ((plsLGM != NULL) && (((IrmTermStructVal*)plsLGM->element->val.pval)->time < dTime))
        plsLGM = plsLGM->next;
    if (plsLGM == NULL)
        plsLGM = ptsLGM->tail;

    /* Variable initializations */
    dPsiFwdMeasure = Psi_func(dFwdTime, ptsLGM);
    dIntegral      = 0.0;
    dUpperTime     = dTime;

    /* Loop over until we reach the beginning of both schedules */
    while (plsVasicek != ptsVasicek->head && plsLGM != ptsLGM->head)
    {
        /* Get the current term structure element */
        pvalVasicek = (IrmTermStructVal*)plsVasicek->element->val.pval;
        pvalLGM     = (IrmTermStructVal*)plsLGM->element->val.pval;

        /* Get the time */
        dLowerTime = max(pvalVasicek->time, pvalLGM->time);
        dTimeDiff  = dUpperTime - dLowerTime;

        /* Calculate the next bit of the integral */
        a      = 1.0 / pvalVasicek->tau;
        lambda = pvalLGM->Lambda;
        dIntegral += pvalLGM->sig / pvalLGM->F * pvalVasicek->sig * pvalVasicek->tau /
                     pvalVasicek->F *
                     ((dPsiFwdMeasure - pvalLGM->Psi - pvalLGM->F / lambda) / (a + lambda) *
                          (exp((a + lambda) * dTimeDiff) - 1.0) +
                      pvalLGM->F / a / lambda * (exp(a * dTimeDiff) - 1.0));

        /* Move backwards in time if we have used up the the entire interval */
        if (pvalVasicek->time == dLowerTime)
            plsVasicek = plsVasicek->previous;
        if (pvalLGM->time == dLowerTime)
            plsLGM = plsLGM->previous;

        /* Get the new times */
        dUpperTime = dLowerTime;
    }

    (*pdMean) += dCorrelation * dIntegral;

    return err;
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Calculates the inflation forward */
Err srt_InflationSpreadForward(
    char* csVasicekUnd, char* csLGMUnd, Date dtExpiry, double dCorrelation, double* pdForward)
{
    /* Variable declarations */
    TermStruct *ptsVasicek, *ptsLGM;
    SrtUndPtr   p_undVasicek, p_undLGM;
    char*       csYCName;
    double      dExpiry;
    Err         err;
    Date        dtToday;

    /* get the und from the underlying name */
    p_undVasicek = lookup_und(csVasicekUnd);
    p_undLGM     = lookup_und(csLGMUnd);

    /* get the underlying term structure */
    if (err = get_underlying_ts(p_undVasicek, &ptsVasicek))
        return err;

    if (err = get_underlying_ts(p_undLGM, &ptsLGM))
        return err;

    /* Get the name of the discount curve associated to the underlying */
    csYCName = get_discname_from_underlying(p_undLGM);

    /* get today from the underlying */
    dtToday = get_today_from_underlying(p_undLGM);

    /* Get the expiry time */
    dExpiry = (dtExpiry - dtToday) * YEARS_IN_DAY;

    /* Calculate the mean of the inflation spread */
    return srt_get_vasicek_forward_measure_und_mean(
        dExpiry, dExpiry, dCorrelation, ptsVasicek, ptsLGM, pdForward);
}

/* ----------------------------------------------------------------------------------------------------------------------------------
 */
/* Calibrates an OU process to market data */
Err srt_CalibrateInflationSpread(
    char*   csYCname,
    double  dSpotSpread,
    double  dRevSpeed,
    double  dCorrelation,
    double  dSpreadVol,
    double  dTau,
    Date    dtToday,
    int     iNumDates,
    double* dvFwdSpread,
    double* dvLGMvol,
    Date*   dtvDates,
    double* dvCalibratedLevels)
{
    /* Variable declarations */
    int     i, j;
    double  dSum1, dSum2, dTimeM1;
    double *dvTimes, dLambda;

    /* Assign memory */
    dvTimes = dvector(0, iNumDates - 1);

    /* Loop over the dates, calculating the times and the difference of the exponentials */
    /* Convert the dates to times */
    dLambda = 1.0 / dTau;
    for (i = 0; i < iNumDates; i++)
    {
        dvTimes[i] = (dtvDates[i] - dtToday) * YEARS_IN_DAY;
    }

    /* Loop over the dates, calcualating the reversion level up to that time */
    dTimeM1 = 0.0;
    for (i = 0; i < iNumDates; i++)
    {
        dSum1 = dvLGMvol[0] *
                ((exp(dRevSpeed * dvTimes[0]) - 1.0) / dRevSpeed -
                 exp(dLambda * dvTimes[i]) * (exp((dRevSpeed + dLambda) * dvTimes[0]) - 1.0) /
                     (dRevSpeed + dLambda));
        dSum2 = 0.0;
        for (j = 0; j < i; j++)
        {
            dSum1 += dvLGMvol[j + 1] *
                     ((exp(dRevSpeed * dvTimes[j + 1]) - exp(dRevSpeed * dvTimes[j])) / dRevSpeed -
                      exp(dLambda * dvTimes[i]) *
                          (exp((dRevSpeed + dLambda) * dvTimes[j + 1]) -
                           exp((dRevSpeed + dLambda) * dvTimes[j])) /
                          (dRevSpeed + dLambda));
            dSum2 += dvCalibratedLevels[j] *
                     (exp(dRevSpeed * dvTimes[j + 1]) - exp(dRevSpeed * dvTimes[j]));
        }
        dSum1 *= dCorrelation * dSpreadVol * dTau * exp(-dRevSpeed * dvTimes[i]);
        dSum2 *= exp(-dRevSpeed * dvTimes[i]);

        dvCalibratedLevels[i] = (dvFwdSpread[i] / swp_f_df(dtToday, dtvDates[i], csYCname) -
                                 dSpotSpread * exp(-dRevSpeed * dvTimes[i]) - dSum1 - dSum2) /
                                (1 - exp(-dRevSpeed * (dvTimes[i] - dTimeM1)));
        dTimeM1 = dvTimes[i];
    }

    /* Free memory */
    free_dvector(dvTimes, 0, iNumDates - 1);

    /* return */
    return 0;
}