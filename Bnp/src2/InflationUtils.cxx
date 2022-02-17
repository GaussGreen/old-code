/* ===================================================================================
   FILENAME:      InflationUtils.cxx / P.A

   PURPOSE:
   ===================================================================================
 */

#include "InflationUtils.h"
#include "Fx3FUtils.h"

/* ----------------------------------------------------------------------------------------------------------
        interp_assetswaptype
        Gestion of type
   -----------------------------------------------------------------------------------------------------------
 */
Err interp_assetswaptype(const char *constStr, InflationAssetSwapType *val) {
  char str[250 + 1];
  strncpy(str, constStr, 250);
  str[250] = '\0';

  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "FLAT_NOTIONAL")) {
    *val = FlatNotional;
    return 0;
  }
  if (!strcmp(str, "FLAT_INDEXED_NOTIONAL")) {
    *val = FlatIndexedNotional;
    return 0;
  }
  if (!strcmp(str, "INDEXED_NOTIONAL")) {
    *val = IndexedNotional;
    return 0;
  }
  if (!strcmp(str, "FLAT_DIRTY_NOTIONAL")) {
    *val = FlatDirtyNotional;
    return 0;
  }
  if (!strcmp(str, "INDEXED_DIRTY_NOTIONAL")) {
    *val = IndexedDirtyNotional;
    return 0;
  }

  return "assetswaptype code not recognized.";
}

/* -----------------------------------------------------------------------------------------------------------------------------
        Free function for structure Model and CPIMKT


   -----------------------------------------------------------------------------------------------------------------------------
 */
void free_InflationModel(InflationModel *Model) {
  Inflation_3FVolTs *VolTs = NULL;

  switch (Model->iModelType) {
  case Infl3F:
    /* Casting the parameter to a 3FVolTs */
    if (Model->ModelVolTs) {
      VolTs = (Inflation_3FVolTs *)(Model->ModelVolTs);

      /* free */
      if (VolTs->dDomSig)
        free(VolTs->dDomSig);
      if (VolTs->dForSig)
        free(VolTs->dForSig);
      if (VolTs->dFxSig)
        free(VolTs->dFxSig);
      if (VolTs->dFxSigtime)
        free(VolTs->dFxSigtime);
      if (VolTs->dRateSigTime)
        free(VolTs->dRateSigTime);
      if (VolTs->dRhoDomFor)
        free(VolTs->dRhoDomFor);
      if (VolTs->dRhoDomFx)
        free(VolTs->dRhoDomFx);
      if (VolTs->dRhoForFx)
        free(VolTs->dRhoForFx);
      if (VolTs->dRhoTime)
        free(VolTs->dRhoTime);

      free(Model->ModelVolTs);
    }

    break;
  }
}

void free_InflationCPIMKT(Inflation_CPIMkt *CPIMKT) {
  if (CPIMKT->CPIFixing.dFixValue)
    free(CPIMKT->CPIFixing.dFixValue);

  if (CPIMKT->CPIFixing.lFixingDate)
    free(CPIMKT->CPIFixing.lFixingDate);
}

/* -----------------------------------------------------------------------------------------------------------------------------
        Inflation dfReal


   -----------------------------------------------------------------------------------------------------------------------------
 */
double dfReal(long lSettle, long lDate, WrapInfo *RealCurve) {
  // Different casting if we are stripping or not
  if (RealCurve->iIsStripping) {
    return (*(RealCurve->StripDf))(lSettle, lDate, RealCurve->BetaPtr);
  } else {
    return swp_f_df(lSettle, lDate, RealCurve->CurveName);
  }
}

double InflationTry(long lToday, long lValDate, WrapInfo *RealCurve) {
  return dfReal(lToday, lValDate, RealCurve);
};

/* -----------------------------------------------------------------------------------------------------------------------------
        MODEL Free : Vol functions
                for all adjustments calculation


                for the moment Inflation_VolTs contains a 3F VolTs but should be
   a generic VolTsModel inside all the following functions we should cast a void
   ptr to the right model


   -----------------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------
        Inflation_GetAdjFwdCPI

        Calculation of
        E(under QTp) of (FFX(Tf        ,Tf) knowing Ft) = E(under QTf) of
  (FFX(Tf       ,Tf) knowing Ft) * Cvxty_adj Cvxty_adj = exp (int(t->Tf)
  (<dFFX(t        ,Tf)/FFX(t        ,Tf)
        ,(G(t        ,Tp)-G(t        ,Tf))dW>

  -------------------------------------------------------------------------------------------------
*/
Err Inflation_GetAdjFwdCPI(               /* Inputs */
                           long lFixDate, /*	Fixing date of the cash CPI */
                           long lPayDate, /*	Pay date  */
                           InflationModel *Model, /*	Model used */

                           /* Outputs */
                           double *ptr_dFwdCPIAdj) {
  long lToday;
  double dAdj;
  double dTf, dTp;
  Err err = NULL;
  Inflation_3FVolTs *VolTs = NULL;

  /* Extract today */
  lToday = Model->lToday;

  /* Case Fixing date is before today */
  if (lFixDate <= lToday) {
    dAdj = 0.0;
  } else {
    /* Time in years calculation */
    dTf = ((double)(lFixDate - lToday)) / DAYSINYEAR;
    dTp = ((double)(lPayDate - lToday)) / DAYSINYEAR;

    /* Switch on model type */
    switch (Model->iModelType) {
    case Infl3F:
      /* Casting the parameter to a 3FVolTs */
      VolTs = (Inflation_3FVolTs *)(Model->ModelVolTs);

      /* Calculation of the integral */
      err = Fx3DtsFwdPayAdjustment_corr(
          0.0, /*	Forward start date */
          dTf, /*	Value date of the forward */
          dTf, /*	Original pay date */
          dTp, /*	Pay date of the forward */
          dTf, /*	Fix date of the forward */

          /*	Model data */
          VolTs->dRateSigTime, VolTs->lRateNbDate, VolTs->dDomSig,
          VolTs->dDomLambda, VolTs->dForSig, VolTs->dForLambda,
          VolTs->dFxSigtime, VolTs->dFxSig, VolTs->lFxNbDate, VolTs->dRhoTime,
          VolTs->dRhoDomFor, VolTs->dRhoDomFx, VolTs->dRhoForFx,
          VolTs->lCorrNbDate,

          /*	Result */
          &dAdj);
      if (err)
        return err;
      break;

    default:
      return "Inflation Model not available";
      break;
    }
  }

  /* return the exp of the integral */
  *ptr_dFwdCPIAdj = exp(dAdj);
  return err;
}

/* ------------------------------------------------------------------------------------------------
        Inflation_GetCPIImpliedVol

        calculate the cumulated vol from Tstart to Tend of FFX(t        ,Tval)
  -------------------------------------------------------------------------------------------------
*/
Err Inflation_GetCPIImpliedVol(/* Inputs */
                               long lValDate, long lStartDate, long lEndDate,
                               InflationModel *Model, /*	Model used */

                               /* Outputs */
                               double *ptr_dCPIImpliedVol) {
  long lToday;
  double dVol;
  double dTval, dTstart, dTend;
  Err err = NULL;
  Inflation_3FVolTs *VolTs = NULL;

  /* Extract today */
  lToday = Model->lToday;

  /* Transform date in time */
  dTval = ((double)(lValDate - lToday)) / DAYSINYEAR;
  dTstart = ((double)(lStartDate - lToday)) / DAYSINYEAR;
  dTend = ((double)(lEndDate - lToday)) / DAYSINYEAR;

  /* Switch on model type */
  switch (Model->iModelType) {
  case Infl3F:
    /* Casting the parameter to a 3FVolTs */
    VolTs = (Inflation_3FVolTs *)Model->ModelVolTs;

    /* Calculation of the implied vol */
    err = Fx3DtsImpliedVol_corr(
        dTval, dTstart, dTend, VolTs->dRateSigTime, VolTs->lRateNbDate,
        VolTs->dDomSig, VolTs->dDomLambda, VolTs->dForSig, VolTs->dForLambda,
        VolTs->dFxSigtime, VolTs->dFxSig, VolTs->lFxNbDate, VolTs->dRhoTime,
        VolTs->dRhoDomFor, VolTs->dRhoDomFx, VolTs->dRhoForFx,
        VolTs->lCorrNbDate,

        /*	Result */
        &dVol);

    if (err)
      return err;

    break;

  default:
    return "Inflation Model not available";
    break;
  }

  /* return the exp */
  *ptr_dCPIImpliedVol = dVol;
  return err;
}

/* ------------------------------------------------------------------------------------------------
        Inflation_GetCorrCPI1CPI2

        calculate the cumulated vol and corr from Tstart to Tend of FFX(t
        ,Tval1) and FFX(t        ,Tval2)
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
                              double *ptr_dVol2) {
  long lToday;
  double dCov, dVar1, dVar2;
  double dTval1, dTval2, dTstart, dTend;
  Err err = NULL;
  Inflation_3FVolTs *VolTs = NULL;

  /* Extract today */
  lToday = Model->lToday;

  /* Transform date in time */
  dTval1 = ((double)(lValDate1 - lToday)) / DAYSINYEAR;
  dTval2 = ((double)(lValDate2 - lToday)) / DAYSINYEAR;
  dTstart = ((double)(lStartDate - lToday)) / DAYSINYEAR;
  dTend = ((double)(lEndDate - lToday)) / DAYSINYEAR;

  /* Switch on model type */
  switch (Model->iModelType) {
  case Infl3F:
    /* Casting the parameter to a 3FVolTs */
    VolTs = (Inflation_3FVolTs *)Model->ModelVolTs;

    /* Calculation of the dVar1 */
    err = Fx3DtsImpliedVol_corr(
        dTval1, min(dTstart, dTval1), min(dTend, dTval1), VolTs->dRateSigTime,
        VolTs->lRateNbDate, VolTs->dDomSig, VolTs->dDomLambda, VolTs->dForSig,
        VolTs->dForLambda, VolTs->dFxSigtime, VolTs->dFxSig, VolTs->lFxNbDate,
        VolTs->dRhoTime, VolTs->dRhoDomFor, VolTs->dRhoDomFx, VolTs->dRhoForFx,
        VolTs->lCorrNbDate,

        /*	Result */
        &dVar1);
    if (err)
      return err;
    dVar1 = dVar1 * dVar1 * (min(dTend, dTval1) - min(dTstart, dTval1));

    /* Calculation of dVar2 */
    err = Fx3DtsImpliedVol_corr(
        dTval2, min(dTstart, dTval2), min(dTend, dTval2), VolTs->dRateSigTime,
        VolTs->lRateNbDate, VolTs->dDomSig, VolTs->dDomLambda, VolTs->dForSig,
        VolTs->dForLambda, VolTs->dFxSigtime, VolTs->dFxSig, VolTs->lFxNbDate,
        VolTs->dRhoTime, VolTs->dRhoDomFor, VolTs->dRhoDomFx, VolTs->dRhoForFx,
        VolTs->lCorrNbDate,

        /*	Result */
        &dVar2);
    if (err)
      return err;

    dVar2 = dVar2 * dVar2 * (min(dTend, dTval2) - min(dTstart, dTval2));

    /* Calculation of the cov12 */
    err = Fx3DtsFwdCumCovar_corr(
        min(min(dTstart, dTval1), dTval2), /*	Forward start date */
        dTval1,                            /*	Value date of the 1st forward */
        dTval2,                            /*	Value date of the 2nd forward */
        min(min(dTend, dTval1), dTval2),   /*	Fix date of both forwards */

        /*	Model data */
        VolTs->dRateSigTime, VolTs->lRateNbDate, VolTs->dDomSig,
        VolTs->dDomLambda, VolTs->dForSig, VolTs->dForLambda, VolTs->dFxSigtime,
        VolTs->dFxSig, VolTs->lFxNbDate, VolTs->dRhoTime, VolTs->dRhoDomFor,
        VolTs->dRhoDomFx, VolTs->dRhoForFx, VolTs->lCorrNbDate,

        /*	Result */
        &dCov);
    if (err)
      return err;

    break;

  default:
    return "Inflation Model not available";
    break;
  }

  /* Vols and correlation  */
  if (fabs(dTend - dTstart) > 1.0e-08) {
    *ptr_dVol1 = sqrt(dVar1 / (min(dTend, max(dTval1, dTval2)) -
                               min(min(dTstart, dTval1), dTval2)));
    *ptr_dVol2 = sqrt(dVar2 / (min(dTend, max(dTval1, dTval2)) -
                               min(min(dTstart, dTval1), dTval2)));
  } else {
    *ptr_dVol1 = 0.0;
    *ptr_dVol2 = 0.0;
  }

  if ((*ptr_dVol1 == 0.0) || (*ptr_dVol2 == 0.0)) {
    *ptr_dCorr12 = 0.0;
  } else {
    *ptr_dCorr12 = dCov / sqrt(dVar1 * dVar2);
  }

  return err;
}

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
                           double *Cov12) {
  double sig_dom, sig_for, sig_fx;
  double T1, T2, t1, t2, U1, U2;
  double covar, covar_partial;
  int i, j, k;
  long StartIndex, EndIndex, StartIndex2, EndIndex2, StartIndex3, EndIndex3;
  Err err = NULL;

  if (T0 > Tfix) {
    err = "end_date before start_date in Fx3DtsFwdCumCovar_corr";
    return err;
  }

  if (Tfix == 0) {
    (*Cov12) = 0;
    return err;
  }

  StartIndex = Get_Index(T0, maturity_fx, nbrMat_fx);
  EndIndex = Get_Index(Tfix, maturity_fx, nbrMat_fx);

  covar = 0.0;

  for (i = StartIndex; i < EndIndex + 1; i++) {
    if (i > StartIndex) {
      T1 = maturity_fx[i - 1];
    } else {
      /* First part */
      T1 = T0;
    }

    if (i == EndIndex || StartIndex == EndIndex) {
      /* Last part */
      T2 = Tfix;
    } else {
      T2 = maturity_fx[i];
    }

    StartIndex2 = Get_Index(T1, maturity_rates, nbrMat);
    EndIndex2 = Get_Index(T2, maturity_rates, nbrMat);

    sig_fx = sig_curve_fx[i];

    for (j = StartIndex2; j < EndIndex2 + 1; j++) {
      if (j > StartIndex2) {
        t1 = maturity_rates[j - 1];
      } else {
        /* First part */
        t1 = T1;
      }

      if (j == EndIndex2 || StartIndex2 == EndIndex2) {
        /* Last part */
        t2 = T2;
      } else {
        t2 = maturity_rates[j];
      }

      sig_dom = sig_curve_dom[j];
      sig_for = sig_curve_for[j];

      StartIndex3 = Get_Index(t1, maturity_corr, nbrMat_corr);
      EndIndex3 = Get_Index(t2, maturity_corr, nbrMat_corr);

      for (k = StartIndex3; k < EndIndex3 + 1; k++) {
        if (k > StartIndex3) {
          U1 = maturity_corr[k - 1];
        } else {
          /* First part */
          U1 = t1;
        }

        if (k == EndIndex3 || StartIndex3 == EndIndex3) {
          /* Last part */
          U2 = t2;
        } else {
          U2 = maturity_corr[k];
        }

        covar_partial = sig_fx * sig_fx * (U2 - U1) +
                        sig_fx * (correl_dom_fx[k] * sig_dom *
                                      (Etha_Func(lda_dom, Tval1, U1, U2) +
                                       Etha_Func(lda_dom, Tval2, U1, U2))

                                  - correl_for_fx[k] * sig_for *
                                        (Etha_Func(lda_for, Tval1, U1, U2) +
                                         Etha_Func(lda_for, Tval2, U1, U2))) +
                        sig_dom * sig_dom *
                            Psi2_Func(lda_dom, lda_dom, Tval1, Tval2, U1, U2) +
                        sig_for * sig_for *
                            Psi2_Func(lda_for, lda_for, Tval1, Tval2, U1, U2) -
                        correl_dom_for[k] * sig_dom * sig_for *
                            (Psi2_Func(lda_dom, lda_for, Tval1, Tval2, U1, U2) +
                             Psi2_Func(lda_for, lda_dom, Tval1, Tval2, U1, U2));

        covar += covar_partial;
      }
    }
  }

  *Cov12 = covar;

  return err;
}

/* -----------------------------------------------------------------------------------------------------------------------------
        3F MODEL

        Specific Lognormal functions






   -----------------------------------------------------------------------------------------------------------------------------
 */

/* ------------------------------------------------------------------------------------------------
        CovarXYZ

        Calculation of E((XT-E(XT))*(YT-E(YT))*(ZT-E(ZT)))
                where X        , Y        , Z are lognormal
  -------------------------------------------------------------------------------------------------
*/
double CovarXYZ(double dT,
                /* Forward */
                double dFwdX, double dFwdY, double dFwdZ,
                /* LogNormal Vol */
                double dVolX, double dVolY, double dVolZ,
                /* correl */
                double dRhoXY, double dRhoXZ, double dRhoYZ) {
  double dCovar;

  dCovar = 2.0 - exp(dRhoXZ * dVolX * dVolZ * dT) -
           exp(dRhoXY * dVolX * dVolY * dT) - exp(dRhoYZ * dVolY * dVolZ * dT) +
           exp((dRhoXZ * dVolX * dVolZ + dRhoXY * dVolX * dVolY +
                dRhoYZ * dVolY * dVolZ) *
               dT);

  dCovar *= dFwdX * dFwdY * dFwdZ;

  return dCovar;
}

/* ------------------------------------------------------------------------------------------------
        CovarXY

        Calculation of E((XT-E(XT))*(YT-E(YT)))
                where X        , Y are lognormal
  -------------------------------------------------------------------------------------------------
*/
double CovarXY(double dT,
               /* Forward */
               double dFwdX, double dFwdY,
               /* LogNormal Vol */
               double dVolX, double dVolY,
               /* correl */
               double dRhoXY) {
  double dCovar;

  dCovar = dFwdX * dFwdY * (exp(dRhoXY * dVolX * dVolY * dT) - 1.0);

  return dCovar;
}

/* ------------------------------------------------------------------------------------------------
        SLMapping_SumLogNormal

        Map  alpha*Xt+Beta*Yt to a SL model
                where	dXt = sx * Xt dWx
                                dYt = sy * Yt dWy
                                <dWx        ,dWy> = rho dt

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
                           double *dSLVol /* Calibrated SL Vol */) {
  double m1, m2;
  double c, d, e;
  double dtemp, dY, dFwdZ;
  Err err = NULL;

  /* Calculation of var(alpha*XT+beta*YT */
  m1 = dAlpha * dAlpha * CovarXY(dT, dFwdX, dFwdX, dVolX, dVolX, 1.0) +
       dBeta * dBeta * CovarXY(dT, dFwdY, dFwdY, dVolY, dVolY, 1.0) +
       2.0 * dAlpha * dBeta * CovarXY(dT, dFwdX, dFwdY, dVolX, dVolY, dRhoXY);

  /* Calculation of the covar of order 3 */
  m2 =
      dAlpha * dAlpha * dAlpha *
          CovarXYZ(dT, dFwdX, dFwdX, dFwdX, dVolX, dVolX, dVolX, 1.0, 1.0,
                   1.0) +
      3.0 * dAlpha * dAlpha * dBeta *
          CovarXYZ(dT, dFwdX, dFwdX, dFwdY, dVolX, dVolX, dVolY, 1.0, dRhoXY,
                   dRhoXY) +
      3.0 * dAlpha * dBeta * dBeta *
          CovarXYZ(dT, dFwdX, dFwdY, dFwdY, dVolX, dVolY, dVolY, dRhoXY, dRhoXY,
                   1.0) +
      dBeta * dBeta * dBeta *
          CovarXYZ(dT, dFwdY, dFwdY, dFwdY, dVolY, dVolY, dVolY, 1.0, 1.0, 1.0);

  /* Find shifted vol */
  d = m1;
  c = d * d / m2;
  e = d / c / c;

  dtemp = exp(1.0 / 3.0 * log(0.5 * e + sqrt(e + 0.25 * e * e) + 1));
  dY = 1.0 / dtemp + dtemp - 1.0;
  *dSLVol = sqrt(log(dY) / dT);

  /* Find FwdZ */
  dFwdZ = c * (dY + 2.0);
  *dSLShift = dFwdZ - (dAlpha * dFwdX + dBeta * dFwdY);

  return err;
}
