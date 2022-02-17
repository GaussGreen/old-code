/* ===================================================================================
   FILENAME:      InflationUtils.h / P.A

   PURPOSE:
   ===================================================================================
 */
#ifndef __INFLATIONUTILS_H
#define __INFLATIONUTILS_H

#include "Fx3FCalib.h"
#include "InflationPricing.h"
#include "srt_h_und_struct.h"
#include "uterror.h"

#define DAYSINYEAR 365.0

/* ----------------------------------------------------------------------------------------------------------
        interp_assetswaptype
        Gestion of type
   -----------------------------------------------------------------------------------------------------------
 */
Err interp_assetswaptype(const char *constStr, InflationAssetSwapType *val);

/* -----------------------------------------------------------------------------------------------------------------------------
        Free function for structure Model and CPIMKT


   -----------------------------------------------------------------------------------------------------------------------------
 */
void free_InflationModel(InflationModel *Model);
void free_InflationCPIMKT(Inflation_CPIMkt *CPIMKT);

/* -----------------------------------------------------------------------------------------------------------------------------
        Inflation dfReal


   -----------------------------------------------------------------------------------------------------------------------------
 */

double dfReal(long lSettle, long lDate, WrapInfo *RealCurve);

/* Try */
double InflationTry(long lToday, long lValDate, WrapInfo *RealCurve);

/* -----------------------------------------------------------------------------------------------------------------------------
        MODEL Free : Vol functions
                for all adjustments calculation


                for the moment Inflation_VolTs contains a 3F VolTs but should be
   a generic VolTsModel inside all the following functions we should cast a void
   ptr to the right


   -----------------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------
        Inflation_GetAdjFwdCPI

        Calculation of
        E(under QTp) of (FFX(Tf  ,Tf) knowing Ft) = E(under QTf) of (FFX(Tf ,Tf)
  knowing Ft) * Cvxty_adj Cvxty_adj = exp (int(t->Tf) (<dFFX(t  ,Tf)/FFX(t  ,Tf)
  ,(G(t  ,Tp)-G(t  ,Tf))dW>

  -------------------------------------------------------------------------------------------------
*/
Err Inflation_GetAdjFwdCPI(               /* Inputs */
                           long lFixDate, /*	Fixing date of the cash CPI */
                           long lPayDate, /*	Pay date  */
                           InflationModel *Model, /*	Model used */

                           /* Outputs */
                           double *ptr_dFwdCPIAdj);

/* ------------------------------------------------------------------------------------------------
        Inflation_GetCPIImpliedVol

        calculate the cumulated vol from Tstart to Tend of FFX(t  ,Tval)
  -------------------------------------------------------------------------------------------------
*/
Err Inflation_GetCPIImpliedVol(/* Inputs */
                               long lValDate, long lStartDate, long lEndDate,
                               InflationModel *Model, /*	Model used */

                               /* Outputs */
                               double *ptr_dCPIImpliedVol);

/* ------------------------------------------------------------------------------------------------
        Inflation_GetCorrCPI1CPI2

        calculate the cumulated vol and corr from Tstart to Tend of FFX(t
  ,Tval1) and FFX(t  ,Tval2)
  -------------------------------------------------------------------------------------------------
*/
Err Inflation_GetCorrCPI1CPI2(                 /* Inputs */
                              long lStartDate, /* Start date */
                              long lEndDate,   /* End date */

                              long lValDate1, /* Val date1 */
                              long lValDate2, /* Val date2 */

                              InflationModel *Model, /*	Model used */

                              /* Outputs */
                              double *ptr_dCorr12, double *ptr_dVol1,
                              double *ptr_dVol2);

/* -----------------------------------------------------------------------------------------------------------------------------
        3F MODEL

        Specific new 3F factors functions






   -----------------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------
        Fx3DtsFwdCumCovar_corr

        Calculate cumulative covariance of F ( mat1 ) and F ( mat2 ) between T0
  and EndDate
  -------------------------------------------------------------------------------------------------
*/
Err Fx3DtsFwdCumCovar_corr(double T0,    /*	Forward start date */
                           double Tval1, /*	Value date of the 1st forward */
                           double Tval2, /*	Value date of the 2nd forward */
                           double Tfix,  /*	End Date */
                           /*	Model data */
                           double *maturity_rates, long nbrMat,
                           double *sig_curve_dom, double lda_dom,
                           double *sig_curve_for, double lda_for,
                           double *maturity_fx, double *sig_curve_fx,
                           long nbrMat_fx, double *maturity_corr,
                           double *correl_dom_for, double *correl_dom_fx,
                           double *correl_for_fx, long nbrMat_corr,
                           /*	Result */
                           double *Cov12);
/* ------------------------------------------------------------------------------------------------
        SLMapping_SumLogNormal

        Map  alpha*Xt+Beta*Yt to a SL model
                where	dXt = sx * Xt dWx
                                dYt = sy * Yt dWy
                                <dWx  ,dWy> = rho dt

        Z + shiftZ = alpha*Xt + Beta*Yt

        The mapping is made by matching the 3 order moments
  -------------------------------------------------------------------------------------------------
*/
Err SLMapping_SumLogNormal(               /* Inputs */
                           double dT,     /* Time to maturity */
                           double dAlpha, /*	coef on X */
                           double dFwdX,  /*	Forward X */
                           double dVolX,  /*  LogNormal vol X */
                           double dBeta,  /*	coef on Y */
                           double dFwdY,  /*	Forward Y */
                           double dVolY,  /*  LogNormal vol Y */
                           double dRhoXY, /*  Correlation */

                           /* Outputs */
                           double *dSLShift, /* Calibrated SL Shift */
                           double *dSLVol /* Calibrated SL Vol */);

#endif