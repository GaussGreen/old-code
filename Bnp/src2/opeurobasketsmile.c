
/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "opfnctns.h"

#include "math.h"
#include "num_h_proba.h"
#include "num_h_sobol.h"

#include "utallhdr.h"

/*******************************************************************************
 *
 * FUNCTION      : srt_f_EuroBasketSmile(...)
 *
 * PURPOSE       : price an option with payoff
 *		 	( sum_i Wi * Xi - K )+
 *
 * DESCRIPTION   : uses a MonteCarlo type numerical integration with SOBOL
 *		  technique to reduce the variance of the result
 *
 * CALLS         : Sobol_init, Sobol_Cube
 *
 * PARAMETERS    : nFwds         - number of Forwards
 *               : Fwds          - forward price of underlying
 *               : k             - strike
 *               : mat           - maturity
 *               : disc          - discount factor
 *               : call_put      - option type
 *				: paths			- number of MC paths used for integration
 *               : greek         - greek wanted (only premium)
 *				: Points		- Num of points for the approx of density of each
 *underlying
 *
 * RETURNS       : ??            - ??
 *
 *******************************************************************************/

Err Free_memory_at_end(
    double**** pppdSobolNoise,
    double**   pdFwdsAtMaturity,
    double***  ppdCorrelCoef,
    double***  ppdPoints,
    double***  ppdCumProbaValue,
    double***  ppdSpreadSquare,
    long       paths,
    long       nFwds,
    long       nNumPoints)
{
    Err err;

    if (*pppdSobolNoise)
        free_f3tensor(*pppdSobolNoise, 1, paths, 0, nFwds - 1, 1, 1);
    *pppdSobolNoise = NULL;
    if (*pdFwdsAtMaturity)
        free_dvector(*pdFwdsAtMaturity, 0, nFwds - 1);
    *pdFwdsAtMaturity = NULL;
    if (*ppdCorrelCoef)
        free_dmatrix(*ppdCorrelCoef, 0, nFwds, 0, nFwds);
    *ppdCorrelCoef = NULL;

    err = sobol_free();
    if (*ppdPoints)
        free_dmatrix(*ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);

    if (*ppdCumProbaValue)
        free_dmatrix(*ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    *ppdPoints        = NULL;
    *ppdCumProbaValue = NULL;

    if (*ppdSpreadSquare)
        free_dmatrix(*ppdSpreadSquare, 0, nFwds - 1, 0, nFwds - 1);
    *ppdSpreadSquare = NULL;

    return err;
}

Err srt_f_EuroBasketSmile(
    int            nFwds,
    double*        Fwds,
    int            nWeights,
    double*        Weights,
    double         strike,
    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           paths,
    long           points,
    SrtGreekType   greek,
    double*        premium,
    double**       ppdSpreadVol)
{
    Err       err             = NULL;
    double*** pppdSobolNoise  = NULL;
    double**  ppdSpreadSquare = NULL;
    double**  ppdCorrelCoef;
    double*   pdFwdsAtMaturity = NULL;
    long      nNumPoints       = points;
    double**  ppdPoints        = NULL;
    double**  ppdCumProbaValue = NULL;
    double    dApproxStdev;
    double    dCoeffAdjustStrike;
    double    dCumProba;
    double    dCoeff;
    double    dAdjStrike;
    double    dTemp;
    double    dTemp1;
    double    dTemp2;
    double    dSum;
    double    dBasket;
    double    dCumVariance;
    double    dCashFlow;
    double    dSqMat;
    long      i, j, k;
    double    cp;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    /* Compute the Coeff matrix from the Correlation matrix used later to correlate the noise */
    ppdCorrelCoef = dmatrix(0, nFwds, 0, nFwds);
    err           = compute_eigen_from_correl(Correlation, nFwds, ppdCorrelCoef);
    if (err)
    {
        free_dmatrix(ppdCorrelCoef, 0, nFwds, 0, nFwds);
        ppdCorrelCoef = NULL;
        return err;
    }

    /* In order to be the more accurate the gonna fill a Sobol cube for one time step */
    pppdSobolNoise = f3tensor(1, paths, 0, nFwds - 1, 1, 1);

    /* Draw the Gaussian Noise thanks to the Correlation and to the Sobol Method */
    /* Initialisation of the Sobol Noise */
    if (err = sobol_init(1, paths, 0, nFwds - 1, 1, 1))
    {
        err = Free_memory_at_end(
            &pppdSobolNoise,
            &pdFwdsAtMaturity,
            &ppdCorrelCoef,
            &ppdPoints,
            &ppdCumProbaValue,
            &ppdSpreadSquare,
            paths,
            nFwds,
            nNumPoints);

        return err;
    }
    /* Compute the Noise */
    if (err = sobol_cube(pppdSobolNoise, 1, paths, 0, nFwds - 1, 1, 1))
    {
        err = Free_memory_at_end(
            &pppdSobolNoise,
            &pdFwdsAtMaturity,
            &ppdCorrelCoef,
            &ppdPoints,
            &ppdCumProbaValue,
            &ppdSpreadSquare,
            paths,
            nFwds,
            nNumPoints);

        return err;
    }

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* If the Spread Vol matrix is allocated then return the Vols */
    if (ppdSpreadVol)
        ppdSpreadSquare = dmatrix(0, nFwds - 1, 0, nFwds - 1);

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */
    for (i = 0; i < nFwds; i++)
    {
        /* set the Approximation for the Stdev */
        dApproxStdev = Vols[i] / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        /* For the shift of the Strike */
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        /* set the value at the boundary 10 stdev */
        ppdPoints[i][0] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev - 10.0 * dApproxStdev);
        ppdCumProbaValue[i][0] = 0.0;
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + 10.0 * dApproxStdev);
        ppdCumProbaValue[i][2 * nNumPoints] = 1.0;

        for (j = 1; j < 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            /* Compute the Cum Proba associated to this points P(S < pdPoints) */

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;
            err        = op_sabr_pricing(
                Fwds[i],
                dAdjStrike,
                mat,
                1,
                SRT_CALL,
                Alpha[i],
                Beta[i],
                Rho[i],
                Vols[i],
                SABR_ATM_BETA,
                &dTemp1);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;
            err = op_sabr_pricing(
                Fwds[i],
                dAdjStrike,
                mat,
                1,
                SRT_CALL,
                Alpha[i],
                Beta[i],
                Rho[i],
                Vols[i],
                SABR_ATM_BETA,
                &dTemp2);

            if (err)
            {
                err = Free_memory_at_end(
                    &pppdSobolNoise,
                    &pdFwdsAtMaturity,
                    &ppdCorrelCoef,
                    &ppdPoints,
                    &ppdCumProbaValue,
                    &ppdSpreadSquare,
                    paths,
                    nFwds,
                    nNumPoints);

                return err;
            }

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
            /*
            dAdjStrike = (log(Fwds[i] / ppdPoints[i][j] ) - 0.5 * Vols[i] * Vols[i] * mat) /
            (Vols[i] * sqrt(mat)); ppdCumProbaValue[i][j] = 1.0 - norm(dAdjStrike);
            */
        }
    }

    /* Compute the Price and the normal vols via a Monte Carlo Methods */
    dSum   = 0.0;
    dSqMat = sqrt(mat);
    for (i = 1; i <= paths; i++)
    {
        /* Correlate the Noise according to the Correaltion Matrix */
        err = correl_random_eigen(pppdSobolNoise[i], 1, 1, nFwds, ppdCorrelCoef);
        if (err)
        {
            err = Free_memory_at_end(
                &pppdSobolNoise,
                &pdFwdsAtMaturity,
                &ppdCorrelCoef,
                &ppdPoints,
                &ppdCumProbaValue,
                &ppdSpreadSquare,
                paths,
                nFwds,
                nNumPoints);

            return err;
        }

        dBasket = 0.0;
        for (j = 0; j < nFwds; j++)
        {
            /* Set the Value of the Forwards at Maturity */
            /* In order to do so we interpolate the value of the Cum proba */

            dCumProba = norm(pppdSobolNoise[i][j][1]);

            pdFwdsAtMaturity[j] =
                interp(ppdCumProbaValue[j], ppdPoints[j], 2 * nNumPoints + 1, dCumProba, 1, &dTemp);

            /* Sum for the Basket */
            dBasket += Weights[j] * pdFwdsAtMaturity[j];

            if (ppdSpreadVol)
            {
                for (k = 0; k < j; k++)
                    ppdSpreadSquare[j][k] += (pdFwdsAtMaturity[j] - pdFwdsAtMaturity[k]) *
                                             (pdFwdsAtMaturity[j] - pdFwdsAtMaturity[k]);
            }
        }

        dCashFlow = cp * (dBasket - strike);
        if (dCashFlow > 0)
            dSum += dCashFlow;
    }

    (*premium) = dSum / paths;
    (*premium) *= disc;

    if (ppdSpreadVol)
    {
        /* Return the Normal Spread Vols Equivalent */
        for (i = 0; i < nFwds; i++)
        {
            ppdSpreadVol[i][i] = 0.0;
            for (j = 0; j < i; j++)
            {
                dCumVariance =
                    ppdSpreadSquare[i][j] / paths - (Fwds[i] - Fwds[j]) * (Fwds[i] - Fwds[j]);
                ppdSpreadVol[i][j] = sqrt(dCumVariance / mat);
                ppdSpreadVol[j][i] = ppdSpreadVol[i][j];
            }
        }
    }

    /* Free the Memory */
    err = Free_memory_at_end(
        &pppdSobolNoise,
        &pdFwdsAtMaturity,
        &ppdCorrelCoef,
        &ppdPoints,
        &ppdCumProbaValue,
        &ppdSpreadSquare,
        paths,
        nFwds,
        nNumPoints);

    return err;
}

/*--------------------------------------------------------------------------------------------------------
This is the same function for pricing a basket under the SABR model. The difference with the
preceding one is that, for each Forward (in total we have nFdws) it draws a random variable which is
compatible with the observed smile, that is a variable distributed accrding to the implied marginal
cumulative distributions. In addition, the copula function define a joint density for the n
forwards. The parameter to input is StudDegree that corresponds to the degree of the Student Copula.
In general for StudDegree->Infinity we recover the gaussian case but in this implementation, the
gaussian copula is given by assigning StudDegree=0
----------------------------------------------------------------------------------------------------------
*/

// FILE *fp2;

Err srt_f_EuroBasketSmileCpl(
    int            nFwds,
    double*        Fwds,
    int            nWeights,
    double*        Weights,
    double         strike,
    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           nPaths,
    long           points,
    int            StudDegree,  // this is the only difference: input the degree of the
                                // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    double**         ppdSpreadVol,
    int              n_conv,
    SrtDiffusionType TypeInput,
    SrtMCSamType     MCType)
{
    Err err = NULL;

    double **ppdSpreadSquare = NULL, *pdFwdsAtMaturity = NULL, **ppdPoints = NULL,
           **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *avg = NULL, *avg2 = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum, dBasket;
    double dCashFlow, dSqMat, cp, sprd = 0.0;

    double n_StdDev, PrMin, RaccRatio;
    //	double Gap,dx,PrMin;

    SrtCallPutType   CallPut  = SRT_CALL;
    SrtDiffusionType log_norm = SRT_LOGNORMAL;

    long i, j, k, kfail, n_eff, nNumPoints = points;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* If the Spread Vol matrix is allocated then return the Vols */
    if (ppdSpreadVol)
    {
        ppdSpreadSquare = dmatrix(0, nFwds - 1, 0, nFwds - 1);
        avg             = dvector(0, nFwds - 1);
        avg2            = dvector(0, nFwds - 1);
    }

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */

    for (i = 0; i < nFwds; i++)
    {
        // set the Approximation for the Stdev //
        dApproxStdev = Vols[i] / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary n_StdDev number of standard deviations
        n_StdDev = 5. * sqrt(mat);

        ppdPoints[i][0] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);
        //		ppdCumProbaValue[i][0] = 0.0;
        //		ppdCumProbaValue[i][2*nNumPoints] = 1.0;

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }

        ppdCumProbaValue[i][0] = 0.0;

        //// find the point where eventually SABR function fails

        kfail = 0;

        if (Fwds[i] >= 0.036)
        {
            n_eff = 0;
        }
        else if (Fwds[i] >= 0.02 && Fwds[i] < 0.036)
        {
            n_eff = nNumPoints;
        }
        else if (Fwds[i] < 0.02)
        {
            n_eff = 2 * nNumPoints - 5;
        }

        for (j = 0; j < n_eff; j++)
        {
            if (ppdCumProbaValue[i][j + 1] <= ppdCumProbaValue[i][j])
            {
                kfail = j + 1;
                PrMin = ppdCumProbaValue[i][j + 1];
            }
        }

        //// then use a lower vovol

        switch (kfail)
        {
        case 0:

            break;
        default:

            for (k = 1; k <= kfail + 1; k++)
            {
                dAdjStrike = ppdPoints[i][k] * dCoeffAdjustStrike;

                dTemp1 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    Vols[i],
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    TypeInput,
                    log_norm,
                    CallPut,
                    greek);

                dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

                dTemp2 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    Vols[i],
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    TypeInput,
                    log_norm,
                    CallPut,
                    greek);

                ppdCumProbaValue[i][k] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][k];
            }

            RaccRatio = PrMin / ppdCumProbaValue[i][kfail + 1];

            for (k = 1; k <= kfail + 1; k++)
            {
                ppdCumProbaValue[i][k] *= RaccRatio;
            }

            break;
        }

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdCumProbaValue[i][j] /= ppdCumProbaValue[i][2 * nNumPoints];
        }
    }

    /*
            for (i=0; i< nFwds; i++)
            {
                    // set the Approximation for the Stdev
                    dApproxStdev = Vols[i] / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

                    // For the shift of the Strike
                    dCoeffAdjustStrike = exp( dApproxStdev / nNumPoints);
                    dCoeff = dCoeffAdjustStrike / (dCoeffAdjustStrike*dCoeffAdjustStrike - 1.0);

                    // set the value at the boundary 10 stdev
                    n_StdDev=6.;

    //		ppdPoints[i][0] = Fwds[i] * exp( - 0.5 * dApproxStdev * dApproxStdev - n_StdDev *
    dApproxStdev);
    //		ppdPoints[i][2*nNumPoints] = Fwds[i] * exp( - 0.5 * dApproxStdev * dApproxStdev +
    n_StdDev * dApproxStdev);

                    ppdPoints[i][0]=DMAX(0.00001,Fwds[i]-5.*sqrt(Fwds[i]*Fwds[i]*(exp(dApproxStdev*dApproxStdev)-1.)));
                    ppdPoints[i][2*nNumPoints]=Fwds[i]+5.*sqrt(Fwds[i]*Fwds[i]*(exp(dApproxStdev*dApproxStdev)-1.));

    //		ppdCumProbaValue[i][0] = 0.0;
    //		ppdCumProbaValue[i][2*nNumPoints] = 1.0;

                    dx=(ppdPoints[i][2*nNumPoints]-ppdPoints[i][0])/(2*nNumPoints);

                    for (j=1; j<=2*nNumPoints; j++)
                    {

                    ppdPoints[i][j]=ppdPoints[i][0]+dx*j;

                            // Compute the Cum Proba associated to this points P(S < pdPoints)

                            dAdjStrike = ppdPoints[i][j-1];
                            dTemp1 = srt_f_optblkschbetastochquick( Fwds[i],
                                                                                                            dAdjStrike,
                                                                                                            mat,Vols[i],
                                                                                                            Alpha[i],
                                                                                                            Beta[i],
                                                                                                            Rho[i],
                                                                                                            1.0,
                                                                                                            TypeInput,
                                                                                                            log_norm,
                                                                                                            CallPut,
                                                                                                            greek);

                            dAdjStrike = ppdPoints[i][j];
                            dTemp2 = srt_f_optblkschbetastochquick( Fwds[i],
                                                                                                            dAdjStrike,
                                                                                                            mat,Vols[i],
                                                                                                            Alpha[i],
                                                                                                            Beta[i],
                                                                                                            Rho[i],
                                                                                                            1.0,
                                                                                                            TypeInput,
                                                                                                            log_norm,
                                                                                                            CallPut,
                                                                                                            greek);

                            ppdCumProbaValue[i][j] = 1.0 + (dTemp2 - dTemp1) / dx;


                    }

                    ppdCumProbaValue[i][0] =
    ppdCumProbaValue[i][1]-fabs(ppdCumProbaValue[i][1]-ppdCumProbaValue[i][2])

                    for (j=0; j<=2*nNumPoints; j++)
                    {
                            ppdCumProbaValue[i][j]-=ppdCumProbaValue[i][0];

                    }

                    for (j=1; j<=2*nNumPoints; j++)
                    {
                            ppdCumProbaValue[i][j]/=ppdCumProbaValue[i][2*nNumPoints];
                    }
            }

    */

    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        avg[i]    = 0.0;
        avg2[i]   = 0.0;
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
        for (j = 0; j < nFwds; j++)
        {
            ppdSpreadSquare[i][j] = 0.0;  // initialise
        }
    }

    CorrCube = dmatrix(0, nPaths, 0, nFwds - 1);
    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        nPaths,
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        0,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum   = 0.0;
    dSqMat = sqrt(mat);

    //		output=dvector(0,2*paths);

    for (i = 3; i < nPaths; i++)
    {
        dBasket = 0.0;
        for (j = 0; j < nFwds; j++)
        {
            pdFwdsAtMaturity[j] = CorrCube[i][j];

            avg[j] += (pdFwdsAtMaturity[j] - avg[j]) / (i - 2);
            dBasket += Weights[j] * pdFwdsAtMaturity[j] * Fwds[j] / avg[j];
            avg2[j] += (pdFwdsAtMaturity[j] * Fwds[j] / avg[j] - avg2[j]) / (i - 2);

            /*
                                            if (ppdSpreadVol)
                                            {
                                                    for (k=0; k<j; k++) {

                                                            sprd_ran=pdFwdsAtMaturity[j] -
               pdFwdsAtMaturity[k]; ppdSpreadSquare[j][k] +=
               ((sprd_ran)*(sprd_ran)-ppdSpreadSquare[j][k])/(i-2);

                                                    }
                                            }
            */
        }

        dCashFlow = cp * (dBasket - strike);
        dSum += (DMAX(dCashFlow, 0.0) - dSum) / (i - 2);
    }

    /*
            SortSeries(paths-2,output);

            fp2= fopen ("C:\\WORKAREA\\Devl\\Cpl_Spread\\StudentTDis_new.txt", "w+");
            for (i=1;i<=paths-2;i++) {
                    fprintf(fp2,"{ %f,  %f}, \n", output[i],1.*i/(paths-2));
            }

            for (i=0;i<2*nNumPoints+1;i++) {
                    fprintf(fp2,"{ %f,  %f}, \n",ppdPoints[1][i], ppdCumProbaValue[1][i]);
            }
                    free_dvector(output,0,2*paths);

            fclose(fp2);


            fp2= fopen ("C:\\WORKAREA\\Devl\\Cpl_Spread\\StudentTDis_new.txt", "w+");
            for (i=1;i<2*nNumPoints+1;i++) {

                    fprintf(fp2,"{ %f,  %f}, \n",(ppdPoints[0][i]+ppdPoints[0][i-1])/2.,
                                                                         (ppdCumProbaValue[0][i]-ppdCumProbaValue[0][i-1])/
                                                                                                     (ppdPoints[0][i]-ppdPoints[0][i-1])
                                                            );
            }
            fclose(fp2);
    */

    (*premium) = disc * dSum;
    /*
            if (ppdSpreadVol)
            {
                     Return the Normal Spread Vols Equivalent
                    for (i=0; i< nFwds; i++)
                    {
                            ppdSpreadVol[i][i] = 0.0;
                            for (j=0; j<i; j++) {

                                    dCumVariance =
       ppdSpreadSquare[i][j]-(Fwds[i]-Fwds[j])*(Fwds[i]-Fwds[j]); ppdSpreadVol[i][j] =
       ppdSpreadVol[j][i]= sqrt( dCumVariance / mat ) ;

                            }
                    }
            }

      */

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdSpreadSquare, 0, nFwds - 1, 0, nFwds - 1);
    free_dvector(avg, 0, nFwds - 1);
    free_dvector(avg2, 0, nFwds - 1);
    free_dmatrix(CorrCube, 0, nPaths, 0, nFwds - 1);

    return err;
}

/*--------------------------------------------------------------------------------------------------------
This is the same function for pricing a standard option on the weighted geometrical average
under the SABR model. The copula function define a joint density for the n forwards.
The parameter to input is StudDegree that corresponds to the degree of the Student Copula. In
general for StudDegree->Infinity we recover the gaussian case but in this implementation, the
gaussian copula is given by assigning StudDegree=0
----------------------------------------------------------------------------------------------------------
*/
Err srt_f_WeightedPdtSmileCpl(
    int            nFwds,
    double*        Fwds,
    int            nWeights,
    double*        Weights,
    double         strike,
    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           nPaths,
    long           points,
    int            StudDegree,  // this is the only difference: input the degree of the
                                // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    double**         ppdSpreadVol,
    double           dCoupon,
    double           dMargin,
    int              n_conv,
    SrtDiffusionType TypeInput,
    SrtMCSamType     MCType)
{
    Err err = NULL;

    double **ppdSpreadSquare = NULL, *pdFwdsAtMaturity = NULL, **ppdPoints = NULL,
           **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *avg = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum;
    double dCashFlow, dSqMat, cp, sprd = 0.0;

    double n_StdDev, PrMin, RaccRatio, dGeomAvg;
    double SigmaBeta;

    SrtCallPutType   CallPut  = SRT_CALL;
    SrtDiffusionType log_norm = SRT_LOGNORMAL;

    long i, j, k, kfail, n_eff, nNumPoints = points;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* If the Spread Vol matrix is allocated then return the Vols */
    if (ppdSpreadVol)
    {
        ppdSpreadSquare = dmatrix(0, nFwds - 1, 0, nFwds - 1);
        avg             = dvector(0, nFwds - 1);
    }

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */

    for (i = 0; i < nFwds; i++)
    {
        /* Sigma Beta */
        srt_f_optsarbvol(
            Fwds[i],
            Fwds[i],
            mat,
            Vols[i],
            Alpha[i],
            Beta[i],
            Rho[i],
            TypeInput,
            SRT_BETAVOL,
            &SigmaBeta);

        // set the Approximation for the Stdev //
        dApproxStdev = SigmaBeta / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary n_StdDev number of standard deviations
        n_StdDev = 5. * sqrt(mat);

        ppdPoints[i][0] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);

        //		ppdCumProbaValue[i][0] = 0.0;
        //		ppdCumProbaValue[i][2*nNumPoints] = 1.0;

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }

        ppdCumProbaValue[i][0] = 0.0;

        //// find the point where eventually SABR function fails

        kfail = 0;

        if (Fwds[i] >= 0.036)
        {
            n_eff = 0;
        }
        else if (Fwds[i] >= 0.02 && Fwds[i] < 0.036)
        {
            n_eff = nNumPoints;
        }
        else if (Fwds[i] < 0.02)
        {
            n_eff = 2 * nNumPoints - 5;
        }

        for (j = 0; j < n_eff; j++)
        {
            if (ppdCumProbaValue[i][j + 1] <= ppdCumProbaValue[i][j])
            {
                kfail = j + 1;
                PrMin = ppdCumProbaValue[i][j + 1];
            }
        }

        //// then use a lower vovol

        switch (kfail)
        {
        case 0:

            break;
        default:

            for (k = 1; k <= kfail + 1; k++)
            {
                dAdjStrike = ppdPoints[i][k] * dCoeffAdjustStrike;

                dTemp1 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

                dTemp2 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                ppdCumProbaValue[i][k] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][k];
            }

            RaccRatio = PrMin / ppdCumProbaValue[i][kfail + 1];

            for (k = 1; k <= kfail + 1; k++)
            {
                ppdCumProbaValue[i][k] *= RaccRatio;
            }

            break;
        }

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdCumProbaValue[i][j] /= ppdCumProbaValue[i][2 * nNumPoints];
        }
    }

    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        avg[i]    = 0.0;
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
        for (j = 0; j < nFwds; j++)
        {
            ppdSpreadSquare[i][j] = 0.0;  // initialise
        }
    }

    CorrCube = dmatrix(0, nPaths, 0, nFwds - 1);
    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        nPaths,
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        0,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum   = 0.0;
    dSqMat = sqrt(mat);

    for (i = 3; i < nPaths; i++)
    {
        dGeomAvg = 1.0;
        for (j = 0; j < nFwds; j++)
        {
            pdFwdsAtMaturity[j] = CorrCube[i][j];
            avg[j] += (pdFwdsAtMaturity[j] - avg[j]) / (i - 2);
            dGeomAvg = dGeomAvg * pow(pdFwdsAtMaturity[j] * Fwds[j] / avg[j], Weights[j]);
        }

        if (dGeomAvg > strike)
        {
            dCashFlow = dGeomAvg + dMargin;
        }
        else
        {
            dCashFlow = dCoupon;
        }
        dSum += (dCashFlow - dSum) / (i - 2);
    }

    (*premium) = disc * dSum;

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdSpreadSquare, 0, nFwds - 1, 0, nFwds - 1);
    free_dvector(avg, 0, nFwds - 1);
    free_dmatrix(CorrCube, 0, nPaths, 0, nFwds - 1);

    return err;
};

/*--------------------------------------------------------------------------------------------------------
This is the same function for pricing a digital BestOf under the SABR model. The
copula function define a joint density for the n forwards.
The parameter to input is StudDegree that corresponds to the degree of the Student Copula. In
general for StudDegree->Infinity we recover the gaussian case but in this implementation, the
gaussian copula is given by assigning StudDegree=0
----------------------------------------------------------------------------------------------------------
*/

Err srt_f_BestOfDigitalSmileCpl(
    int            nFwds,
    double*        Fwds,
    int            nWeights,
    double*        Weights,
    double         strike,
    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           nPaths,
    long           points,
    int            StudDegree,  // this is the only difference: input the degree of the
                                // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    double**         ppdSpreadVol,
    double           dCoupon,
    double           dMargin,
    int              n_conv,
    SrtDiffusionType TypeInput,
    SrtMCSamType     MCType)
{
    Err err = NULL;

    double **ppdSpreadSquare = NULL, *pdFwdsAtMaturity = NULL, **ppdPoints = NULL,
           **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *avg = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum;
    double dCashFlow, dSqMat, cp, sprd = 0.0;

    double n_StdDev, PrMin, RaccRatio, dBestOf;
    double SigmaBeta;

    SrtCallPutType   CallPut  = SRT_CALL;
    SrtDiffusionType log_norm = SRT_LOGNORMAL;

    long i, j, k, kfail, n_eff, nNumPoints = points;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* If the Spread Vol matrix is allocated then return the Vols */
    if (ppdSpreadVol)
    {
        ppdSpreadSquare = dmatrix(0, nFwds - 1, 0, nFwds - 1);
        avg             = dvector(0, nFwds - 1);
    }

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */

    for (i = 0; i < nFwds; i++)
    {
        /* Sigma Beta */
        srt_f_optsarbvol(
            Fwds[i],
            Fwds[i],
            mat,
            Vols[i],
            Alpha[i],
            Beta[i],
            Rho[i],
            TypeInput,
            SRT_BETAVOL,
            &SigmaBeta);

        // set the Approximation for the Stdev //
        dApproxStdev = SigmaBeta / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary n_StdDev number of standard deviations
        n_StdDev = 5. * sqrt(mat);

        ppdPoints[i][0] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);

        //		ppdCumProbaValue[i][0] = 0.0;
        //		ppdCumProbaValue[i][2*nNumPoints] = 1.0;

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }

        ppdCumProbaValue[i][0] = 0.0;

        //// find the point where eventually SABR function fails

        kfail = 0;

        if (Fwds[i] >= 0.036)
        {
            n_eff = 0;
        }
        else if (Fwds[i] >= 0.02 && Fwds[i] < 0.036)
        {
            n_eff = nNumPoints;
        }
        else if (Fwds[i] < 0.02)
        {
            n_eff = 2 * nNumPoints - 5;
        }

        for (j = 0; j < n_eff; j++)
        {
            if (ppdCumProbaValue[i][j + 1] <= ppdCumProbaValue[i][j])
            {
                kfail = j + 1;
                PrMin = ppdCumProbaValue[i][j + 1];
            }
        }

        //// then use a lower vovol

        switch (kfail)
        {
        case 0:

            break;
        default:

            for (k = 1; k <= kfail + 1; k++)
            {
                dAdjStrike = ppdPoints[i][k] * dCoeffAdjustStrike;

                dTemp1 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

                dTemp2 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                ppdCumProbaValue[i][k] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][k];
            }

            RaccRatio = PrMin / ppdCumProbaValue[i][kfail + 1];

            for (k = 1; k <= kfail + 1; k++)
            {
                ppdCumProbaValue[i][k] *= RaccRatio;
            }

            break;
        }

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdCumProbaValue[i][j] /= ppdCumProbaValue[i][2 * nNumPoints];
        }
    }

    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        avg[i]    = 0.0;
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
        for (j = 0; j < nFwds; j++)
        {
            ppdSpreadSquare[i][j] = 0.0;  // initialise
        }
    }

    CorrCube = dmatrix(0, nPaths, 0, nFwds - 1);
    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        nPaths,
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        0,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum   = 0.0;
    dSqMat = sqrt(mat);

    for (i = 3; i < nPaths; i++)
    {
        dBestOf = 0.0;
        for (j = 0; j < nFwds; j++)
        {
            pdFwdsAtMaturity[j] = CorrCube[i][j];

            avg[j] += (pdFwdsAtMaturity[j] - avg[j]) / (i - 2);
            dBestOf = max(dBestOf, pdFwdsAtMaturity[j] * Fwds[j] / avg[j] + Weights[j]);
        }

        if (dBestOf > strike)
        {
            dCashFlow = dBestOf + dMargin;
        }
        else
        {
            dCashFlow = dCoupon;
        }
        dSum += (dCashFlow - dSum) / (i - 2);
    }

    (*premium) = disc * dSum;

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdSpreadSquare, 0, nFwds - 1, 0, nFwds - 1);
    free_dvector(avg, 0, nFwds - 1);
    free_dmatrix(CorrCube, 0, nPaths, 0, nFwds - 1);

    return err;
}

/*--------------------------------------------------------------------------------------------------------
This is the same function for pricing a digital BestOf under the SABR model. The
copula function define a joint density for the n forwards.
The parameter to input is StudDegree that corresponds to the degree of the Student Copula. In
general for StudDegree->Infinity we recover the gaussian case but in this implementation, the
gaussian copula is given by assigning StudDegree=0
----------------------------------------------------------------------------------------------------------
*/

Err srt_f_BestOfSmileCpl(
    int            nFwds,
    double*        Fwds,
    int            nWeights,
    double*        Weights,
    double         strike,
    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           nPaths,
    long           points,
    int            StudDegree,  // this is the only difference: input the degree of the
                                // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    double**         ppdSpreadVol,
    int              n_conv,
    SrtDiffusionType TypeInput,
    SrtMCSamType     MCType)
{
    Err err = NULL;

    double **ppdSpreadSquare = NULL, *pdFwdsAtMaturity = NULL, **ppdPoints = NULL,
           **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *avg = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum;
    double dSqMat, cp, sprd = 0.0;

    double n_StdDev, PrMin, RaccRatio, dBestOf;
    double SigmaBeta;

    SrtCallPutType   CallPut  = SRT_CALL;
    SrtDiffusionType log_norm = SRT_LOGNORMAL;

    long i, j, k, kfail, n_eff, nNumPoints = points;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* If the Spread Vol matrix is allocated then return the Vols */
    if (ppdSpreadVol)
    {
        ppdSpreadSquare = dmatrix(0, nFwds - 1, 0, nFwds - 1);
        avg             = dvector(0, nFwds - 1);
    }

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */

    for (i = 0; i < nFwds; i++)
    {
        /* Sigma Beta */
        srt_f_optsarbvol(
            Fwds[i],
            Fwds[i],
            mat,
            Vols[i],
            Alpha[i],
            Beta[i],
            Rho[i],
            TypeInput,
            SRT_BETAVOL,
            &SigmaBeta);

        // set the Approximation for the Stdev //
        dApproxStdev = SigmaBeta / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary n_StdDev number of standard deviations
        n_StdDev = 5. * sqrt(mat);

        ppdPoints[i][0] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);

        //		ppdCumProbaValue[i][0] = 0.0;
        //		ppdCumProbaValue[i][2*nNumPoints] = 1.0;

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                SigmaBeta,
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                SRT_BETAVOL,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }

        ppdCumProbaValue[i][0] = 0.0;

        //// find the point where eventually SABR function fails

        kfail = 0;

        if (Fwds[i] >= 0.036)
        {
            n_eff = 0;
        }
        else if (Fwds[i] >= 0.02 && Fwds[i] < 0.036)
        {
            n_eff = nNumPoints;
        }
        else if (Fwds[i] < 0.02)
        {
            n_eff = 2 * nNumPoints - 5;
        }

        for (j = 0; j < n_eff; j++)
        {
            if (ppdCumProbaValue[i][j + 1] <= ppdCumProbaValue[i][j])
            {
                kfail = j + 1;
                PrMin = ppdCumProbaValue[i][j + 1];
            }
        }

        //// then use a lower vovol

        switch (kfail)
        {
        case 0:

            break;
        default:

            for (k = 1; k <= kfail + 1; k++)
            {
                dAdjStrike = ppdPoints[i][k] * dCoeffAdjustStrike;

                dTemp1 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

                dTemp2 = srt_f_optblkschbetastochquick(
                    Fwds[i],
                    dAdjStrike,
                    mat,
                    SigmaBeta,
                    0.0,
                    Beta[i],
                    Rho[i],
                    1.0,
                    SRT_BETAVOL,
                    log_norm,
                    CallPut,
                    greek);

                ppdCumProbaValue[i][k] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][k];
            }

            RaccRatio = PrMin / ppdCumProbaValue[i][kfail + 1];

            for (k = 1; k <= kfail + 1; k++)
            {
                ppdCumProbaValue[i][k] *= RaccRatio;
            }

            break;
        }

        for (j = 1; j <= 2 * nNumPoints; j++)
        {
            ppdCumProbaValue[i][j] /= ppdCumProbaValue[i][2 * nNumPoints];
        }
    }

    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        avg[i]    = 0.0;
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
        for (j = 0; j < nFwds; j++)
        {
            ppdSpreadSquare[i][j] = 0.0;  // initialise
        }
    }

    CorrCube = dmatrix(0, nPaths, 0, nFwds - 1);
    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        nPaths,
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        0,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum   = 0.0;
    dSqMat = sqrt(mat);

    for (i = 3; i < nPaths; i++)
    {
        dBestOf = 0.0;
        for (j = 0; j < nFwds; j++)
        {
            pdFwdsAtMaturity[j] = CorrCube[i][j];
            avg[j] += (pdFwdsAtMaturity[j] - avg[j]) / (i - 2);
            dBestOf = max(dBestOf, max(pdFwdsAtMaturity[j] * Fwds[j] / avg[j] - Weights[j], 0.0));
        }

        dSum += (DMAX(0.0, cp * (dBestOf - strike)) - dSum) / (i - 2);
    }

    (*premium) = disc * dSum;

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdSpreadSquare, 0, nFwds - 1, 0, nFwds - 1);
    free_dvector(avg, 0, nFwds - 1);
    free_dmatrix(CorrCube, 0, nPaths, 0, nFwds - 1);

    return err;
}

/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/

Err srt_f_CallableBasketCpl(
    double** InFwds,
    double*  InWeights,
    double*  InStrikes,
    double** InVols,
    double*  InBeta,
    double*  InAlpha,
    double*  InRho,
    double*  InLvl,
    double** InCorr,
    double   mat,
    double   disc,

    SrtCallPutType call_put,

    long paths,
    long points,
    int  StudDegree,  // this is the only difference: input the degree of the
                      // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    int              n_conv,
    const long       nFwdsRows,
    const long       nFwdsCols,
    SrtDiffusionType TypeInput,
    SrtDiffusionType log_norm)
{
    Err err = NULL;

    double * pdFwdsAtMaturity = NULL, **ppdPoints = NULL, **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *Fwds = NULL, *Vols = NULL, **Correlation;
    double * Alpha = NULL, *Beta = NULL, *Rho = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum, dBasket0,
        dBasket1;
    double dCashFlow, cp, sprd = 0.0;

    double n_StdDev;

    SrtCallPutType CallPut = SRT_CALL;
    SrtMCSamType   MCType  = SOBOL;
    long           i, j, l, k, m, nNumPoints = points, nFwds = 0;

    Fwds  = dvector(0, nFwdsRows * nFwdsCols - 1);
    Vols  = dvector(0, nFwdsRows * nFwdsCols - 1);
    Alpha = dvector(0, nFwdsRows * nFwdsCols - 1);
    Beta  = dvector(0, nFwdsRows * nFwdsCols - 1);
    Rho   = dvector(0, nFwdsRows * nFwdsCols - 1);

    Correlation = dmatrix(0, nFwdsRows * nFwdsCols - 1, 0, nFwdsRows * nFwdsCols - 1);

    //////////////////////////////////////////  defines the forwards and the vols

    k = 0;
    for (i = 0; i < nFwdsCols; i++)
    {
        for (j = 0; j < nFwdsRows; j++)
        {
            Fwds[k] = InFwds[j][i];
            Vols[k] = InVols[j][i];

            Alpha[k] = InAlpha[j];
            Beta[k]  = InBeta[j];
            Rho[k]   = InRho[j];

            k++;
        }
    }
    nFwds = k;

    l = 0;
    for (i = 0; i < nFwds; i++)
    {
        Correlation[i][i] = 1.0;

        if (l != 2)
        {
            k = 0;
            for (j = 0; j < i; j++)
            {
                Correlation[i][j] = InCorr[0][1] / pow(fabs(1. * i - 1. * j), 0.1);
                Correlation[j][i] = Correlation[i][j];

                if (k == 2)
                {
                    k                 = -1;
                    Correlation[i][j] = Correlation[j][i] = InCorr[0][2];
                }
                k++;
            }

            l++;
        }
        else if (l == 2)
        {
            k = 0;
            for (j = 0; j < i; j++)
            {
                Correlation[i][j] = Correlation[j][i] = InCorr[0][2];

                if (k == 2)
                {
                    k = -1;

                    Correlation[i][j] = InCorr[0][1] / pow(fabs(1. * i - 1. * j), 0.1);
                    Correlation[j][i] = Correlation[i][j];
                }
                k++;
            }

            l = 0;
        }
    }

    ////////////////////////////////////////////////////////////////////////////

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */
    for (i = 0; i < nFwds; i++)
    {
        // set the Approximation for the Stdev //
        dApproxStdev = Vols[i] / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary 10 stdev
        n_StdDev = 6.;

        ppdPoints[i][0]        = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdCumProbaValue[i][0] = 0.0;
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);
        ppdCumProbaValue[i][2 * nNumPoints] = 1.0;

        for (j = 1; j < 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////   Deconvolve the marginals
    ///////////////////////////////////////////////////////////////////////////////////////////////

    /*
                    ConvolveMarginals(
                                                            ppdPoints,
                                                            ppdCumProbaValue,
                                                            2*nNumPoints+1,
                                                            nFwds,
                                                            n_conv,
                                                            -1
                                                     );
    */
    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
    }

    CorrCube = dmatrix(0, 2 * (paths + 1) * (n_conv + 1), 0, nFwds - 1);

    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        (paths + 1) * (n_conv + 1),
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        0,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum = 0.0;

    //		output=dvector(0,2*paths);

    for (i = 3; i < 2 * paths; i++)
    {
        m         = 0;
        dCashFlow = 0.0;

        for (k = 0; k < nFwdsCols; k++)
        {
            dBasket0 = InStrikes[0];
            dBasket1 = 0.0;

            for (j = 0; j < nFwdsRows - 1; j++)
            {
                pdFwdsAtMaturity[j] = 0.0;

                for (l = 1; l <= n_conv + 1; l++)
                {
                    pdFwdsAtMaturity[j] += CorrCube[i + (l - 1) * (paths + 1)][m];
                }

                dBasket1 += InWeights[j] * pdFwdsAtMaturity[j];
                m++;
            }

            pdFwdsAtMaturity[nFwdsRows - 1] = 0.0;
            for (l = 1; l <= n_conv + 1; l++)
            {
                pdFwdsAtMaturity[nFwdsRows - 1] += CorrCube[i + (l - 1) * (paths + 1)][m];
            }

            dCashFlow += (dBasket0 + DMAX(cp * (dBasket1 - InStrikes[1]), 0.0) -
                          DMAX(cp * (dBasket1 - InStrikes[2]), 0.0) -
                          InWeights[nFwdsRows - 1] * pdFwdsAtMaturity[nFwdsRows - 1]) *
                         InLvl[k];

            m++;  // counter for the forwards
        }
        dSum += (DMAX(dCashFlow, 0.0) - dSum) / (i - 2);
    }

    (*premium) = disc * dSum;

    free_dvector(Fwds, 0, nFwdsRows * nFwdsCols - 1);
    free_dvector(Vols, 0, nFwdsRows * nFwdsCols - 1);
    free_dvector(Alpha, 0, nFwdsRows * nFwdsCols - 1);
    free_dvector(Beta, 0, nFwdsRows * nFwdsCols - 1);
    free_dvector(Rho, 0, nFwdsRows * nFwdsCols - 1);

    free_dmatrix(Correlation, 0, nFwdsRows * nFwdsCols - 1, 0, nFwdsRows * nFwdsCols - 1);

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(CorrCube, 0, 2 * (paths + 1) * (n_conv + 1), 0, nFwds - 1);

    return err;
}

/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/
/*-------------------------------------------------------------------------------------------------------*/

Err srt_f_CompoundBasketCpl(
    int     nFwds,
    double* Fwds,
    int     nWeights,
    double* Weights,

    int     nStrikes,
    double* Strikes,

    int            nVols,
    double*        Vols,
    int            nBeta,
    double*        Beta,
    int            nAlpha,
    double*        Alpha,
    int            nRho,
    double*        Rho,
    double**       Correlation,
    double         mat,
    double         disc,
    SrtCallPutType call_put,
    long           paths,
    long           points,
    int            StudDegree,  // this is the only difference: input the degree of the
                                // Student t Copula. 0 Corresponds to the gaussian copula
    SrtGreekType     greek,
    double*          premium,
    int              n_conv,
    SrtDiffusionType TypeInput,
    SrtDiffusionType log_norm)
{
    Err err = NULL;

    double * pdFwdsAtMaturity = NULL, **ppdPoints = NULL, **ppdCumProbaValue = NULL;
    double **CorrCube = NULL, *mean_v = NULL, *CorrMatr = NULL;

    double dApproxStdev, dCoeffAdjustStrike, dCoeff, dAdjStrike, dTemp1, dTemp2, dSum, dBasket0,
        dBasket1;
    double dCashFlow, dSqMat, cp, sprd = 0.0;

    double n_StdDev;

    SrtCallPutType CallPut = SRT_CALL;
    SrtMCSamType   MCType  = SOBOL;
    long           i, j, l, nNumPoints = points, idum = -6123787;

    /* Set the coefficient for the payoff */
    cp = (call_put == SRT_CALL) ? 1 : -1;

    pdFwdsAtMaturity = dvector(0, nFwds - 1);

    /* Allocate Memory for the Inverse Cum Proba computation */
    ppdPoints        = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);
    ppdCumProbaValue = dmatrix(0, nFwds - 1, 0, 2 * nNumPoints);

    /* Define the points */
    /* Compute the cumulative proba = 1 + d_K Call(K) */

    for (i = 0; i < nFwds; i++)
    {
        // set the Approximation for the Stdev //
        dApproxStdev = Vols[i] / pow(Fwds[i], 1.0 - Beta[i]) * sqrt(mat);

        // For the shift of the Strike
        dCoeffAdjustStrike = exp(dApproxStdev / nNumPoints);
        dCoeff             = dCoeffAdjustStrike / (dCoeffAdjustStrike * dCoeffAdjustStrike - 1.0);

        // set the value at the boundary 10 stdev
        n_StdDev = 6.;

        ppdPoints[i][0]        = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev -
                                        n_StdDev * dApproxStdev);  // modifier....
        ppdCumProbaValue[i][0] = 0.0;
        ppdPoints[i][2 * nNumPoints] =
            Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev + n_StdDev * dApproxStdev);
        ppdCumProbaValue[i][2 * nNumPoints] = 1.0;

        for (j = 1; j < 2 * nNumPoints; j++)
        {
            ppdPoints[i][j] = Fwds[i] * exp(-0.5 * dApproxStdev * dApproxStdev +
                                            (j - nNumPoints) * 4.0 / nNumPoints * dApproxStdev);

            // Compute the Cum Proba associated to this points P(S < pdPoints)

            dAdjStrike = ppdPoints[i][j] * dCoeffAdjustStrike;

            dTemp1 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            dAdjStrike /= dCoeffAdjustStrike * dCoeffAdjustStrike;

            dTemp2 = srt_f_optblkschbetastochquick(
                Fwds[i],
                dAdjStrike,
                mat,
                Vols[i],
                Alpha[i],
                Beta[i],
                Rho[i],
                1.0,
                TypeInput,
                log_norm,
                CallPut,
                greek);

            ppdCumProbaValue[i][j] = 1.0 + (dTemp1 - dTemp2) * dCoeff / ppdPoints[i][j];
        }
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////   Deconvolve the marginals
    ///////////////////////////////////////////////////////////////////////////////////////////////

    ConvolveMarginals(ppdPoints, ppdCumProbaValue, 2 * nNumPoints + 1, nFwds, n_conv, -1);

    // loads the cube of correlated variables compatible with the Student Copula for n
    // degrees of freedom

    mean_v = dvector(0, nFwds - 1);
    for (i = 0; i < nFwds; i++)
    {
        mean_v[i] = 0.0;  // resets the means of the gaussian variables to 0
    }

    CorrMatr = dvector(0, nFwds - 1);

    for (i = 3; i < paths; i++)
    {
        /*
                        CorrMatr = GetStudentCplDev_Rand (
                                                                        StudDegree,
                                                                        mean_v,
                                                                        Correlation,
                                                                        paths+1,
                                                                        StudDegree+nFwds,
                                                                        ppdPoints,
                                                                        ppdCumProbaValue,
                                                                        2*nNumPoints+1,
                                                                        n_conv,
                                                                        &idum
                                                          );
        */
    }

    CorrCube = dmatrix(0, 2 * (paths + 1) * (n_conv + 1), 0, nFwds - 1);
    GetStudentCplDev(
        StudDegree,
        mean_v,
        Correlation,
        (paths + 1) * (n_conv + 1),
        StudDegree + nFwds,
        ppdPoints,
        ppdCumProbaValue,
        2 * nNumPoints + 1,
        n_conv,
        MCType,
        CorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */

    dSum   = 0.0;
    dSqMat = sqrt(mat);

    //		output=dvector(0,2*paths);

    for (i = 3; i < 2 * paths; i++)
    {
        dBasket0 = Strikes[0];
        dBasket1 = 0.0;

        for (j = 0; j < nFwds - 1; j++)
        {
            pdFwdsAtMaturity[j] = 0.0;

            for (l = 1; l <= n_conv + 1; l++)
            {
                pdFwdsAtMaturity[j] += CorrCube[i + (l - 1) * (paths + 1)][j];
            }

            dBasket1 += Weights[j] * pdFwdsAtMaturity[j];
        }

        pdFwdsAtMaturity[nFwds - 1] = 0.0;
        for (l = 1; l <= n_conv + 1; l++)
        {
            pdFwdsAtMaturity[nFwds - 1] += CorrCube[i + (l - 1) * (paths + 1)][nFwds - 1];
        }

        dCashFlow = dBasket0 + DMAX(cp * (dBasket1 - Strikes[1]), 0.0) -
                    DMAX(cp * (dBasket1 - Strikes[2]), 0.0) -
                    Weights[nFwds - 1] * pdFwdsAtMaturity[nFwds - 1];

        dSum += (DMAX(dCashFlow, 0.0) - dSum) / (i - 2);
        //			dSum += (dCashFlow-dSum)/(i-2);
    }

    (*premium) = disc * dSum;

    free_dvector(pdFwdsAtMaturity, 0, nFwds - 1);
    free_dvector(mean_v, 0, nFwds - 1);

    free_dmatrix(ppdPoints, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(ppdCumProbaValue, 0, nFwds - 1, 0, 2 * nNumPoints);
    free_dmatrix(CorrCube, 0, 2 * (paths + 1) * (n_conv + 1), 0, nFwds - 1);

    free_dvector(CorrMatr, 0, nFwds - 1);

    return err;
}

Err srt_f_Basket_MixedSL_Copula(
    double dMaturity,

    int     iNbFwd,
    double* dShiftedFwdDown,
    double* dShiftedFwdUp,
    double* dProba,
    double* dShift,
    double* dVolDown,
    double* dVolUp,

    double** dCorrelation,

    double* dWeights,
    int     iNbStrike,
    double* dStrike,

    SrtCallPutType eCallPut,
    long           lNbPaths,
    long           lNbPoints,
    double         dNbStd,
    int            iStudDegree,  // this is the only difference: input the degree of the
                                 // Student t Copula. 0 Corresponds to the gaussian copula
    SrtMCSamType MCType,
    int          iNConv,
    double*      dPremium)
{
    Err err = NULL;

    double *pdFwdsAtMaturity = NULL, **dPoints = NULL, *dDriftDown = NULL, *dDriftUp = NULL,
           **dCumProbaValue = NULL;

    double **dCorrCube = NULL, *dMeanV = NULL;

    double dApproxStdev, dApproxDrift, dApproxShiftedFwd, dx;
    double dBasket, dCashFlow, cp;
    double dConv;

    SrtCallPutType   CallPut  = SRT_CALL;
    SrtDiffusionType log_norm = SRT_LOGNORMAL;

    long i, j, lNbPoints2;

    /* Set the coefficient for the payoff */
    cp         = (eCallPut == SRT_CALL) ? 1.0 : -1.0;
    lNbPoints2 = 2 * lNbPoints;

    /* Allocate Memory for the Inverse Cum Proba computation */
    pdFwdsAtMaturity = calloc(iNbFwd, sizeof(double));
    dDriftDown       = calloc(iNbFwd, sizeof(double));
    dDriftUp         = calloc(iNbFwd, sizeof(double));
    dMeanV           = calloc(iNbFwd, sizeof(double));

    dPoints        = dmatrix(0, iNbFwd - 1, 0, lNbPoints2);
    dCumProbaValue = dmatrix(0, iNbFwd - 1, 0, lNbPoints2);
    dCorrCube      = dmatrix(0, lNbPaths, 0, iNbFwd - 1);

    if (!pdFwdsAtMaturity || !dDriftDown || !dDriftUp || !dPoints || !dCumProbaValue || !dMeanV ||
        !dCorrCube)
    {
        err = "Memory allocation faillure in srt_f_Spread_MixedSL_Copula";
        goto FREE_RETURN;
    }

    /* Define the points and the cumulative */
    for (i = 0; i < iNbFwd; i++)
    {
        // set the Approximation for the Stdev //
        dApproxStdev      = dProba[i] * dVolDown[i] + (1.0 - dProba[i]) * dVolUp[i];
        dApproxDrift      = -0.5 * (dProba[i] * dVolDown[i] * dVolDown[i] +
                               (1.0 - dProba[i]) * dVolUp[i] * dVolUp[i]);
        dApproxShiftedFwd = dProba[i] * dShiftedFwdDown[i] + (1.0 - dProba[i]) * dShiftedFwdUp[i];

        if (dApproxShiftedFwd < 0)
        {
            dConv = -1.0;
        }
        else
        {
            dConv = 1.0;
        }

        dVolDown[i] *= sqrt(dMaturity);
        dDriftDown[i] = -0.5 * dVolDown[i] * dVolDown[i];
        dDriftDown[i] -= log(dApproxShiftedFwd / dShiftedFwdDown[i]);

        dVolUp[i] *= sqrt(dMaturity);
        dDriftUp[i] = -0.5 * dVolUp[i] * dVolUp[i];
        dDriftUp[i] -= log(dApproxShiftedFwd / dShiftedFwdUp[i]);

        dx = 2.0 * dNbStd * dApproxStdev / (lNbPoints2 * 1.0);

        dPoints[i][0] = dApproxDrift - dConv * dNbStd * dApproxStdev;

        dCumProbaValue[i][0] = 0.0;

        for (j = 1; j <= lNbPoints2; j++)
        {
            dPoints[i][j] = dPoints[i][j - 1] + dConv * dx;

            dCumProbaValue[i][j] = dProba[i] * norm((dPoints[i][j] - dDriftDown[i]) / dVolDown[i]);
            dCumProbaValue[i][j] +=
                (1.0 - dProba[i]) * norm((dPoints[i][j] - dDriftUp[i]) / dVolUp[i]);

            if (dConv < 0.0)
            {
                dCumProbaValue[i][j] = 1.0 - dCumProbaValue[i][j];
            }
        }

        dCumProbaValue[i][lNbPoints2] = 1.0;
    }

    /* Calculates the true coordinates */
    for (i = 0; i < iNbFwd; i++)
    {
        dApproxShiftedFwd = dProba[i] * dShiftedFwdDown[i] + (1.0 - dProba[i]) * dShiftedFwdUp[i];

        if (dApproxShiftedFwd < 0)
        {
            dConv = -1.0;
        }
        else
        {
            dConv = 1.0;
        }

        for (j = 0; j <= lNbPoints2; j++)
        {
            dPoints[i][j] = dApproxShiftedFwd * exp(dPoints[i][j]) - dShift[i];
        }
    }

    /* loads the cube of correlated variables compatible with the Student Copula for n degrees of
     * freedom */
    for (i = 0; i < iNbFwd; i++)
    {
        dMeanV[i] = 0.0;  // resets the means of the gaussian variables to 0
    }

    GetStudentCplDev(
        iStudDegree,
        dMeanV,
        dCorrelation,
        lNbPaths,
        iStudDegree + iNbFwd,
        dPoints,
        dCumProbaValue,
        lNbPoints2 + 1,
        0,
        MCType,
        dCorrCube);

    /* Compute the Price and the normal vols via a Monte Carlo Methods */
    memset(dPremium, 0, iNbStrike * sizeof(double));

    for (i = 3; i < lNbPaths; i++)
    {
        dBasket = 0.0;

        for (j = 0; j < iNbFwd; j++)
        {
            dBasket += dWeights[j] * dCorrCube[i][j];
        }

        for (j = 0; j < iNbStrike; j++)
        {
            dCashFlow = cp * (dBasket - dStrike[j]);
            dPremium[j] += (DMAX(dCashFlow, 0.0) - dPremium[j]) / (1.0 * (i - 2));
        }
    }

FREE_RETURN:

    if (pdFwdsAtMaturity)
        free(pdFwdsAtMaturity);
    if (dDriftDown)
        free(dDriftDown);
    if (dDriftUp)
        free(dDriftUp);
    if (dMeanV)
        free(dMeanV);

    if (dPoints)
        free_dmatrix(dPoints, 0, iNbFwd - 1, 0, lNbPoints2);
    if (dCumProbaValue)
        free_dmatrix(dCumProbaValue, 0, iNbFwd - 1, 0, lNbPoints2);
    if (dCorrCube)
        free_dmatrix(dCorrCube, 0, lNbPaths, 0, iNbFwd - 1);

    return err;
}