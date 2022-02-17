/* =======================================================================================

   FUNCNAME        :swp_f_cmsopt

   DESCRIPTION     :Contains all versions of
                                        value of cms (rate and option) through
   integration with swaptions

   =======================================================================================
 */

#include "math.h"
#include "opfnctns.h"
#include "opsabrgeneric.h"
#include "opsabrgenericcalib.h"
#include "opsabrgenericinterp.h"
#include "swp_h_all.h"
#include "swp_h_cms.h"
#include "swp_h_cmsopt.h"
#include "swp_h_vol.h"
#include "utmemory.h"

/****** TEST FUNCTIONS *******************/

/* VOL FUNCTIONS */
static double CAPSTRIKE = 0.20;
static double BOUNDARYINT = 0.50;
static double CMS_SPACE_NUMSTEP = 600;
static int CMS_SPACE_TYPE = 0; /* 1 normal  , 0 exponential */
static int CMS_SPACE_TYPE_F = 1;
static int CMS_SPACE_TYPE_C = 0;

Err swp_f_cmsoption(double forward, double num_periods,
                    SrtCompounding frequency, double strike, double volatility,
                    double maturity, SrtReceiverType pay_or_receive,
                    double delay, double delta, int num_swaps,
                    SrtDiffusionType lognrm_nrm, double *ans) {
  Err err;

  err = swp_f_Cms_Option(forward, maturity, num_periods, strike,
                         (double)frequency, pay_or_receive, delay, 1,
                         lognrm_nrm, volatility, 0, (Date)0, (Date)0, SRT_FALSE,
                         0, NULL, 0, NULL, ans);
  return err;
}

Err swp_f_Cms_Option(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* CMS Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dFrequency,                  /* CMS Frequency */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol  , 1: linear interpolation  , 2: FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dCMSOptionValue) {
  Err err = NULL;
  double *pdCoeffApprox = NULL, *pdVolVector = NULL;
  double dFactor_n = 1.0, dVol = 0.0, dVal = 0.0, dAdjForwardRate, dDeltaStrike,
         dOpt, r, dInvLvl, dAmt, dInvLvlm1, dInvLvlm2, dSens, dSens1;
  double power = 0.0; /* only used for the call of the vols */
  int i, m = 0, iLeftRightPos = 0, dSignRecPay = 1, iLower, iHigher;
  SRT_Boolean bConverged = SRT_FALSE;
  SrtCallPutType CallPut;
  int OldType;

  if (err = conv_rec_in_call(PayRec, &CallPut))
    return err;

  /* Correct the forward and strike to bring them to an annual basis */
  dAdjForwardRate = dFwdSwapRate * dRateConv;
  dStrike *= dRateConv;

  /* Convert from delay in years to delay in periods */
  dDelay *= dFrequency;

  /* check for consistency */
  if (VolType == SRT_NORMAL) {
    OldType = CMS_SPACE_TYPE;
    CMS_SPACE_TYPE = 1;
  }

  /* Special cases first */
  if (dMaturity <= 0.0) {
    /* At maturity or zero vol  , so value
       is MAX(dForward-dStrike  ,0) for EPAY  , etc */
    dVal = dAdjForwardRate - dStrike;
    if (PayRec == SRT_RECEIVER)
      dVal *= -1.0;
    /* max (...  ,0) */
    if (dVal < 0.0)
      dVal = 0.0;
  }
  /* Maturity is non zero */
  else {
    if (PayRec == SRT_RECEIVER)
      CMS_SPACE_TYPE = CMS_SPACE_TYPE_F;
    else
      CMS_SPACE_TYPE = CMS_SPACE_TYPE_C;

    /* Compute the direction of the integration PAY : up  , REC : down */
    if (PayRec == SRT_RECEIVER)
      dSignRecPay *= -1;

    /* If Error in the vol matrix then switch to constant vol*/
    if (iMethod == 1 && lNumStrikesInVol == 0)
      iMethod = 0;

    /* Set a flat vol just in case
            otherwise get the vol for the strike */
    if (!iMethod)
      dVol = dFlatVol;
    else if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, dStrike, &dVol,
                                  &power))
      return err;

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
      dDeltaStrike = exp(log((-(dSignRecPay - 1.0) * 1.0e-10 +
                              (dSignRecPay + 1.0) * BOUNDARYINT / 2.0) /
                             dStrike) /
                         CMS_SPACE_NUMSTEP);
    else
      dDeltaStrike = (-(dSignRecPay - 1.0) * 1.0e-10 +
                      (dSignRecPay + 1.0) * BOUNDARYINT / 2.0 - dStrike) /
                     CMS_SPACE_NUMSTEP;

    /* Now compute the Smile linear interpolation */
    if (iMethod == 1) {
      pdCoeffApprox = dvector(0, lNumStrikesInVol);
      if (pdCoeffApprox == NULL)
        return "Cms allocation Error";

      pdVolVector = dvector(0, lNumStrikesInVol);
      if (pdVolVector == NULL)
        return "Cms allocation Error";

      /* Compute linear interpolation coeficients between two vols
         and get the position of the strike in the grid */
      i = 1;
      if (CMS_SPACE_TYPE == SRT_NORMAL)
        r = dStrike + 2.0 * dDeltaStrike;
      else
        r = dStrike * dDeltaStrike * dDeltaStrike;

      while (i < lNumStrikesInVol && r >= pdStrikesVol[i])
        i++;
      m = i - 1;

      if (PayRec == SRT_PAYER) {
        iLeftRightPos = 1;
        iLower = max(m - iLeftRightPos, 1);
        iHigher = lNumStrikesInVol - 1;
      } else {
        iLower = 1;
        iHigher = min(m + 1, lNumStrikesInVol - 1);
      }

      if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd,
                               pdStrikesVol[iLower - 1],
                               &pdVolVector[iLower - 1], &power))
        return err;

      for (i = iLower; i <= iHigher; i++) {
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, pdStrikesVol[i],
                                 &pdVolVector[i], &power))
          return err;
        pdCoeffApprox[i - 1] = (pdVolVector[i] - pdVolVector[i - 1]) /
                               (pdStrikesVol[i] - pdStrikesVol[i - 1]);
      }
      pdCoeffApprox[lNumStrikesInVol - 1] = .0;
    }

    /* Treat the Receiver case first: put on the rate */
    if (PayRec == SRT_RECEIVER) {
      /* If the Model is Lognormal  , we have to restrict to positive strikes */
      if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
        /* MAX(dStrike-dForward  ,0) will always be zero here */
        if (dStrike < CMS_EPSILON) {
          *dCMSOptionValue = 0.0;
          return err;
        }
        /* If there is not enough space for two swaptions  , just do the one and
         * return */
        if ((dStrike * (1.0 - dDeltaStrike * dDeltaStrike)) >
            (dStrike - CMS_EPSILON)) {
          if (VolType == SRT_NORMAL)
            *dCMSOptionValue =
                srt_f_optblknrm(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                                CallPut, PREMIUM);
          else
            *dCMSOptionValue =
                srt_f_optblksch(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                                CallPut, PREMIUM);

          dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
          dSens = (1.0 - (1.0 / dSens)) / dAdjForwardRate;
          *dCMSOptionValue = *dCMSOptionValue * dSens / dRateConv;
          return err;
        }

      } /* END if (VolType == SRTLOGNORMAL) */
    }
    /* ... now this is a payer (call on the rate) */
    else if (dStrike < CMS_EPSILON)
      dStrike = CMS_EPSILON;

    /* First term: swaption  , centered on strike */
    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                             CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                             CallPut, PREMIUM);

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
      r = dDeltaStrike * dStrike / dFrequency;
    else
      r = (dDeltaStrike + dStrike) / dFrequency;

    dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
    dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
    dVal = dInvLvl * dOpt;

    /* Second term: swaption  , centered on r * frequency */
    if (iMethod)
      if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                               &dVol, &power))
        return err;

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
      r *= dDeltaStrike;
      dFactor_n = dDeltaStrike;
    } else
      r += dDeltaStrike;

    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);

    dInvLvlm1 = dInvLvl;
    dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
    dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));

    if (CMS_SPACE_TYPE == SRT_NORMAL)
      dAmt = 2.0 * (dInvLvl - dInvLvlm1);
    else
      dAmt = (dDeltaStrike + 1.0) / dDeltaStrike * (dInvLvl - dInvLvlm1);

    dVal += dAmt * dOpt;

    /* Following Terms: Swaptions centered on Strike + i * Delta_k */
    i = 0;
    while ((bConverged == SRT_FALSE) && (i <= CMS_SPACE_NUMSTEP)) {
      /* Get the Vol at strike K + (i - 1) * delta_k */
      if (iMethod == 2)
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                                 &dVol, &power))
          return err;

      if (iMethod == 1) {
        if (dSignRecPay *
                (r * dFrequency - pdStrikesVol[min(m, lNumStrikesInVol - 1)]) >=
            0.0)
          m += dSignRecPay;
        m = min(lNumStrikesInVol - 1 + iLeftRightPos, max(m, iLeftRightPos));
        dVol = pdVolVector[m - iLeftRightPos] +
               pdCoeffApprox[m - iLeftRightPos] *
                   (r * dFrequency - pdStrikesVol[m - iLeftRightPos]);
        if (r * dFrequency <= pdStrikesVol[0])
          dVol = pdVolVector[0];
      }

      /* Move to the next strike (K + i * delta_k) or K * (delta_k)^i at
         which we intersect the profile */
      if (VolType == SRT_LOGNORMAL)
        dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);
      else
        dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);

      if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
        r *= dDeltaStrike;
        dFactor_n *= dDeltaStrike;
      } else
        r += dDeltaStrike;

      /* Coefficient from the series (recursive in 2 steps)
         to compute the weights */
      dInvLvlm2 = dInvLvlm1;
      dInvLvlm1 = dInvLvl;
      if (fabs(r) > 1.0e-10) {
        dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
        dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
      } else
        dInvLvl = 1.0 / dNumPeriods;

      if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
        dAmt = 1.0 / (dFactor_n * (dDeltaStrike - 1.0)) *
               ((dDeltaStrike * dFactor_n - 1.0) * dInvLvl -
                (dDeltaStrike + 1.0) * (dFactor_n - 1.0) * dInvLvlm1 +
                (dFactor_n - dDeltaStrike) * dInvLvlm2);
      else
        dAmt = i * dInvLvl - 2.0 * (i - 1) * dInvLvlm1 + (i - 2) * dInvLvlm2;

      dVal += dAmt * dOpt;

      /* two ways of convergence :: the integral has converged or the boundary
       * is reached */
      /* ((fabs(dAmt * dOpt)) < CMS_TOLERANCE) || */
      if ((r * dFrequency > BOUNDARYINT)) {
        /* Series has converged */
        bConverged = SRT_TRUE;
        break;
      }
      i++;
    }

    dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dDelay);
    dSens1 = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
    dSens *= (1.0 - (1.0 / dSens1)) / dAdjForwardRate;

    dVal *= dSens;
    if (iMethod == 1) {
      free_dvector(pdCoeffApprox, 0, lNumStrikesInVol);
      free_dvector(pdVolVector, 0, lNumStrikesInVol);
      pdCoeffApprox = NULL;
      pdVolVector = NULL;
    }

  } /* END if(dMaturity != 0.0 ) */

  /* retrieve the Old param */
  if (VolType == SRT_NORMAL)
    CMS_SPACE_TYPE = OldType;

  /* Return the Basis adjusted Rate */
  *dCMSOptionValue = dVal / dRateConv;
  return err;
}

/* Calculates the CMS option taking into account the fact that in NY  ,
market data are physical-delivery swaptions and not cash-settle swaptions */
Err swp_f_Cms_OptionNY(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* CMS Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dFrequency,                  /* CMS Frequency */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol  , 1: linear interpolation  , 2: FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dCMSOptionValue,
    String ycName, String refRateCode) {
  Err err = NULL;
  double *pdCoeffApprox = NULL, *pdVolVector = NULL;
  long dStartSliding, dEndSliding;
  double dFactor_n = 1.0, dVol = 0.0, dVal = 0.0, dAdjForwardRate, dDeltaStrike,
         dOpt, r, dInvLvl, dAmt, dInvLvlm1, dInvLvlm2, dSens, dSens1;
  long spot;
  double power = 0.0; /* only used for the call of the vols */
  int i, m = 0, iLeftRightPos = 0, dSignRecPay = 1, iLower, iHigher;
  SRT_Boolean bConverged = SRT_FALSE;
  SrtCallPutType CallPut;
  int OldType;
  SwapDP swapdp;
  SrtCurvePtr yccrv;
  double level, levelcash;
  double forward;
  double dSens1Sliding;

  /* Gets the yield curve */
  yccrv = lookup_curve(ycName);
  if (!yccrv)
    return serror("Could not find yc in swp_f_Swaption", ycName);

  /* Get the spot lag from the Yield Curve and sets it into the SwapDP */
  (&swapdp)->spot_lag = get_spotlag_from_curve(yccrv);

  /* Get spot and clcn date */
  spot = get_spotdate_from_yldcrv(yccrv);
  dStartSliding = spot;
  dEndSliding = spot + (dEnd - dStart);
  dEndSliding = add_unit(dEndSliding, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* Creates the SwapDP from start  ,end  , comp and basis */
  err =
      swp_f_initSwapDP(dStartSliding, dEndSliding, "SEMIANNUAL", "BB", &swapdp);
  if (err)
    return err;

  /* Computes the level payment of the swap */
  err = swp_f_Level_SwapDP(&swapdp, ycName, &level);
  if (err)
    return err;

  /* Computes the swap */
  err = swp_f_ForwardRate_SwapDP(&swapdp, ycName, refRateCode, &forward);
  if (err)
    return err;

  /* Computes the cash level */
  dSens1Sliding = pow(1.0 + (forward / dFrequency), dNumPeriods);
  levelcash = (1.0 - (1.0 / dSens1Sliding)) / forward;

  if (err = conv_rec_in_call(PayRec, &CallPut))
    return err;

  /* Correct the forward and strike to bring them to an annual basis */
  dAdjForwardRate = dFwdSwapRate * dRateConv;
  dStrike *= dRateConv;

  /* Convert from delay in years to delay in periods */
  dDelay *= dFrequency;

  /* check for consistency */
  if (VolType == SRT_NORMAL) {
    OldType = CMS_SPACE_TYPE;
    CMS_SPACE_TYPE = 1;
  }

  /* Special cases first */
  if (dMaturity <= 0.0) {
    /* At maturity or zero vol  , so value
       is MAX(dForward-dStrike  ,0) for EPAY  , etc */
    dVal = dAdjForwardRate - dStrike;
    if (PayRec == SRT_RECEIVER)
      dVal *= -1.0;
    /* max (...  ,0) */
    if (dVal < 0.0)
      dVal = 0.0;
  }

  /* Maturity is non zero */
  else {
    /* Compute the direction of the integration PAY : up  , REC : down */
    if (PayRec == SRT_RECEIVER)
      dSignRecPay *= -1;

    /* If Error in the vol matrix then switch to constant vol*/
    if (iMethod == 1 && lNumStrikesInVol == 0)
      iMethod = 0;

    /* Set a flat vol just in case
            otherwise get the vol for the strike */
    if (!iMethod)
      dVol = dFlatVol;
    else if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, dStrike, &dVol,
                                  &power))
      return err;

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
      dDeltaStrike = exp(log((-(dSignRecPay - 1.0) * 1.0e-10 +
                              (dSignRecPay + 1.0) * BOUNDARYINT / 2.0) /
                             dStrike) /
                         CMS_SPACE_NUMSTEP);
    else
      dDeltaStrike = (-(dSignRecPay - 1.0) * 1.0e-10 +
                      (dSignRecPay + 1.0) * BOUNDARYINT / 2.0 - dStrike) /
                     CMS_SPACE_NUMSTEP;

    /* Now compute the Smile linear interpolation */
    if (iMethod == 1) {
      pdCoeffApprox = dvector(0, lNumStrikesInVol);
      if (pdCoeffApprox == NULL)
        return "Cms allocation Error";

      pdVolVector = dvector(0, lNumStrikesInVol);
      if (pdVolVector == NULL)
        return "Cms allocation Error";

      /* Compute linear interpolation coeficients between two vols
         and get the position of the strike in the grid */
      i = 1;
      if (CMS_SPACE_TYPE == SRT_NORMAL)
        r = dStrike + 2.0 * dDeltaStrike;
      else
        r = dStrike * dDeltaStrike * dDeltaStrike;

      while (i < lNumStrikesInVol && r >= pdStrikesVol[i])
        i++;
      m = i - 1;

      if (PayRec == SRT_PAYER) {
        iLeftRightPos = 1;
        iLower = max(m - iLeftRightPos, 1);
        iHigher = lNumStrikesInVol - 1;
      } else {
        iLower = 1;
        iHigher = min(m + 1, lNumStrikesInVol - 1);
      }

      if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd,
                               pdStrikesVol[iLower - 1],
                               &pdVolVector[iLower - 1], &power))
        return err;

      for (i = iLower; i <= iHigher; i++) {
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, pdStrikesVol[i],
                                 &pdVolVector[i], &power))
          return err;
        pdCoeffApprox[i - 1] = (pdVolVector[i] - pdVolVector[i - 1]) /
                               (pdStrikesVol[i] - pdStrikesVol[i - 1]);
      }
      pdCoeffApprox[lNumStrikesInVol - 1] = .0;
    }

    /* Treat the Receiver case first: put on the rate */
    if (PayRec == SRT_RECEIVER) {
      /* If the Model is Lognormal  , we have to restrict to positive strikes */
      if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
        /* MAX(dStrike-dForward  ,0) will always be zero here */
        if (dStrike < CMS_EPSILON) {
          *dCMSOptionValue = 0.0;
          return err;
        }
        /* If there is not enough space for two swaptions  , just do the one and
         * return */
        if ((dStrike * (1.0 - dDeltaStrike * dDeltaStrike)) >
            (dStrike - CMS_EPSILON)) {
          if (VolType == SRT_NORMAL)
            *dCMSOptionValue =
                srt_f_optblknrm(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                                CallPut, PREMIUM);
          else
            *dCMSOptionValue =
                srt_f_optblksch(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                                CallPut, PREMIUM);

          dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
          dSens = (1.0 - (1.0 / dSens)) / dAdjForwardRate;
          *dCMSOptionValue = *dCMSOptionValue * dSens / dRateConv;
          return err;
        }

      } /* END if (VolType == SRTLOGNORMAL) */
    }
    /* ... now this is a payer (call on the rate) */
    else if (dStrike < CMS_EPSILON)
      dStrike = CMS_EPSILON;

    /* First term: swaption  , centered on strike */
    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                             CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, dStrike, dVol, dMaturity, 1.0,
                             CallPut, PREMIUM);

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
      r = dDeltaStrike * dStrike / dFrequency;
    else
      r = (dDeltaStrike + dStrike) / dFrequency;

    dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
    dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
    dVal = dInvLvl * dOpt;

    /* Second term: swaption  , centered on r * frequency */
    if (iMethod)
      if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                               &dVol, &power))
        return err;

    if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
      r *= dDeltaStrike;
      dFactor_n = dDeltaStrike;
    } else
      r += dDeltaStrike;

    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);

    dInvLvlm1 = dInvLvl;
    dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
    dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));

    if (CMS_SPACE_TYPE == SRT_NORMAL)
      dAmt = 2.0 * (dInvLvl - dInvLvlm1);
    else
      dAmt = (dDeltaStrike + 1.0) / dDeltaStrike * (dInvLvl - dInvLvlm1);

    dVal += dAmt * dOpt;

    /* Following Terms: Swaptions centered on Strike + i * Delta_k */
    i = 0;
    while ((bConverged == SRT_FALSE) && (i <= CMS_SPACE_NUMSTEP)) {
      /* Get the Vol at strike K + (i - 1) * delta_k */
      if (iMethod == 2)
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                                 &dVol, &power))
          return err;

      if (iMethod == 1) {
        if (dSignRecPay *
                (r * dFrequency - pdStrikesVol[min(m, lNumStrikesInVol - 1)]) >=
            0.0)
          m += dSignRecPay;
        m = min(lNumStrikesInVol - 1 + iLeftRightPos, max(m, iLeftRightPos));
        dVol = pdVolVector[m - iLeftRightPos] +
               pdCoeffApprox[m - iLeftRightPos] *
                   (r * dFrequency - pdStrikesVol[m - iLeftRightPos]);
        if (r * dFrequency <= pdStrikesVol[0])
          dVol = pdVolVector[0];
      }

      /* Move to the next strike (K + i * delta_k) or K * (delta_k)^i at
         which we intersect the profile */
      if (VolType == SRT_LOGNORMAL)
        dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);
      else
        dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);

      if (CMS_SPACE_TYPE == SRT_LOGNORMAL) {
        r *= dDeltaStrike;
        dFactor_n *= dDeltaStrike;
      } else
        r += dDeltaStrike;

      /* Coefficient from the series (recursive in 2 steps)
         to compute the weights */
      dInvLvlm2 = dInvLvlm1;
      dInvLvlm1 = dInvLvl;
      if (fabs(r) > 1.0e-10) {
        dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
        dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
      } else
        dInvLvl = 1.0 / dNumPeriods;

      if (CMS_SPACE_TYPE == SRT_LOGNORMAL)
        dAmt = 1.0 / (dFactor_n * (dDeltaStrike - 1.0)) *
               ((dDeltaStrike * dFactor_n - 1.0) * dInvLvl -
                (dDeltaStrike + 1.0) * (dFactor_n - 1.0) * dInvLvlm1 +
                (dFactor_n - dDeltaStrike) * dInvLvlm2);
      else
        dAmt = i * dInvLvl - 2.0 * (i - 1) * dInvLvlm1 + (i - 2) * dInvLvlm2;

      dVal += dAmt * dOpt;

      /* two ways of convergence :: the integral has converged or the boundary
       * is reached */
      /* ((fabs(dAmt * dOpt)) < CMS_TOLERANCE) || */
      if ((r * dFrequency > BOUNDARYINT)) {
        /* Series has converged */
        bConverged = SRT_TRUE;
        break;
      }
      i++;
    }

    dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dDelay);
    dSens1 = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
    dSens *= (1.0 - (1.0 / dSens1)) / dAdjForwardRate;

    dVal *= dSens;
    if (iMethod == 1) {
      free_dvector(pdCoeffApprox, 0, lNumStrikesInVol);
      free_dvector(pdVolVector, 0, lNumStrikesInVol);
      pdCoeffApprox = NULL;
      pdVolVector = NULL;
    }

  } /* END if(dMaturity != 0.0 ) */

  /* retrieve the Old param */
  if (VolType == SRT_NORMAL)
    CMS_SPACE_TYPE = OldType;

  /* Return the Basis adjusted Rate */
  *dCMSOptionValue = dVal / dRateConv;

  /* Because in NY  , there are no cash settle swaptions */
  *dCMSOptionValue *= level / levelcash;
  return err;
}

/*	This calculates the price of a TEC option (taking into account the
 second order convexity of the underlying rate an of the non linear payoff) */
Err swp_f_Tec_Option(
    double dFwdSwapRate,                /*Forward Swap Rate */
    double dMaturity,                   /* TEC Maturity as double */
    double dNumPeriods, double dStrike, /* Option Strike */
    double dMargin,                     /* the Tec Margin */
    double dFrequency,                  /* TEC Frequency */
    double dPaymentPower,               /* TEC Power */
    SrtReceiverType PayRec,             /* Pay or Rec */
    double dDelay,                      /* Delay */
    double dRateConv,                   /* date adjustment */
    SrtDiffusionType VolType,           /* Vol type Lognormal or Normal */
    double dFlatVol,                    /* Flat Vol if used */
    int iMethod, /* 0: Use Flat Vol  , 1: linear interpolation  , 2: FullSmile*/
    Date dStart, /* The following parameters are used in the GetVol function */
    Date dEnd,   /* and are useless in the rest of the code */
    SRT_Boolean bAdjForSpread, double dSpread, char *szVolCurveName,
    long lNumStrikesInVol, double *pdStrikesVol, double *dTECOptionValue) {

  Err err = NULL;

  double *pdCoeffApprox = NULL, *pdVolVector = NULL;
  double dVol = .0, dVal = .0, dDeltaStrike, dAdjForwardRate, dOpt = .0, r = .0,
         dInvLvl, dAmt = .0, dLvl, dInvLvlm1, dInvLvlm2, dSens, dSens1;
  double power = 0.0;
  int i, j, m = 0, iNumSwaps, dSignRecPay = 1, iLeftRightPos = 0, iLower,
            iHigher;
  SRT_Boolean bConverged = SRT_FALSE;
  SrtCallPutType CallPut;
  double *dTECw = NULL;

  if (err = conv_rec_in_call(PayRec, &CallPut))
    return err;

  /* Correct the forward  , strike and margin
     to bring them to an annual basis */
  dAdjForwardRate = dFwdSwapRate * dRateConv;
  dStrike *= dRateConv;
  dMargin *= dRateConv;

  /* Convert from delay in years to delay in periods */
  dDelay *= dFrequency;

  /* Special cases first */
  if (dMaturity <= 0.0) {
    /* At maturity or zero vol  , so value is
                    MAX((1+CmsForward+m)^(1/k)-(1+dStrike)^(1/k)  ,0) for EPAY
                      , etc */
    dVal = pow((1 + dFwdSwapRate + dMargin), (1 / dPaymentPower)) -
           pow((1 + dStrike), (1 / dPaymentPower));

    if (PayRec == SRT_RECEIVER)
      dVal *= -1.0;

    if (dVal < 0.0)
      dVal = 0.0;
  }
  /* ...Maturity is non zero */
  else {
    /* Compute the direction of the integration PAY : up  , REC : down */
    if (PayRec == SRT_RECEIVER)
      dSignRecPay *= -1;

    /* If Error in the vol matrix then switch to constant vol*/
    if (iMethod == 1 && lNumStrikesInVol == 0)
      iMethod = 0;

    /* Set a flat vol just in case
            otherwise get the vol for the strike */
    if (!iMethod)
      dVol = dFlatVol;
    else if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd,
                                  dStrike - dMargin, &dVol, &power))
      return err;

    /* Let's be sure that we won't use strikes < 0
    dDeltaStrike = DEFAULT_CMS_DELTA_STRIKE;
    iNumSwaps = MAX_CMS_SWAPS;
    if (PayRec == SRT_RECEIVER)
            iNumSwaps = min((int)((dStrike - dMargin - CMS_EPSILON) /
    dDeltaStrike)  , iNumSwaps);

     reverse the sign if needed
    dDeltaStrike *= dSignRecPay; */

    dDeltaStrike =
        (-(dSignRecPay - 1.0) * 1.0e-10 +
         (dSignRecPay + 1.0) * BOUNDARYINT / 2.0 - dStrike + dMargin) /
        CMS_SPACE_NUMSTEP;
    iNumSwaps = (int)CMS_SPACE_NUMSTEP;

    /* Now compute the Smile linear interpolation */
    if (iMethod == 1) {
      pdCoeffApprox = dvector(0, lNumStrikesInVol);
      if (pdCoeffApprox == NULL)
        return "Tec allocation Error";

      pdVolVector = dvector(0, lNumStrikesInVol);
      if (pdVolVector == NULL)
        return "Tec allocation Error";

      /* Compute linear interpolation coeficients between two vols
         and get the position of the strike in the grid */
      i = 1;
      r = 2.0 * dDeltaStrike + dStrike - dMargin;
      while (i < lNumStrikesInVol && r >= pdStrikesVol[i])
        i++;
      m = i - 1;

      if (PayRec == SRT_PAYER) {
        iLeftRightPos = 1;
        iLower = max(m - iLeftRightPos, 1);
        iHigher = lNumStrikesInVol - 1;
      } else {
        iLower = 1;
        iHigher = min(m + 1, lNumStrikesInVol - 1);
      }

      if (err = err = swp_f_truncvol(szVolCurveName, dStart, dEnd,
                                     pdStrikesVol[iLower - 1],
                                     &pdVolVector[iLower - 1], &power))
        return err;

      for (i = iLower; i <= iHigher; i++) {
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, pdStrikesVol[i],
                                 &pdVolVector[i], &power))
          return err;
        pdCoeffApprox[i - 1] = (pdVolVector[i] - pdVolVector[i - 1]) /
                               (pdStrikesVol[i] - pdStrikesVol[i - 1]);
      }
      pdCoeffApprox[lNumStrikesInVol - 1] = .0;
    }

    /* Treat the Receiver case first: put on the rate */
    if (PayRec == SRT_RECEIVER) {
      /* If the Model is Lognormal  , we have to restrict to positive strikes */
      if (VolType == SRT_LOGNORMAL) {
        if (dStrike - dMargin < CMS_EPSILON) {
          *dTECOptionValue = 0.0;
          return err;
        }
        /* If there is not enough space for two swaptions  , just do the one and
         * return */
        if (2.0 * dDeltaStrike > dStrike - dMargin - CMS_EPSILON) {
          if (VolType == SRT_NORMAL)
            dVal = srt_f_optblknrm(dAdjForwardRate, dStrike - dMargin, dVol,
                                   dMaturity, 1.0, CallPut, PREMIUM);
          else
            dVal = srt_f_optblksch(dAdjForwardRate, dStrike - dMargin, dVol,
                                   dMaturity, 1.0, CallPut, PREMIUM);

          r = (dDeltaStrike + dStrike - dMargin) / dFrequency;
          dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
          dSens1 = pow(1.0 + r, dNumPeriods);
          dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dDelay) *
                  (1.0 - (1.0 / dSens)) / dAdjForwardRate;
          dSens1 = pow(1.0 + r, dDelay) * (1.0 - (1.0 / dSens1)) / r;
          dTECw[0] = dFrequency * dSens1;
          dVal *= (1 / dPaymentPower) *
                  (pow((1 + dStrike), (1 / dPaymentPower) - 1) +
                   (1 / dPaymentPower - 1) * dDeltaStrike *
                       pow((1 + dStrike), (1 / dPaymentPower) - 2)) *
                  dTECw[0] * dSens;

          *dTECOptionValue = dVal / dRateConv;
          return err;
        }

      } /* END if (m_eVolType == ELOGNORMAL) */
    }
    /* ... now this is a payer (call on the rate) */
    else if (dStrike < CMS_EPSILON + dMargin)
      dStrike = CMS_EPSILON + dMargin;

    r = (dDeltaStrike + dStrike - dMargin) / dFrequency;

    /* Initialise the dVector dTECw */
    dTECw = dvector(0, (long)CMS_SPACE_NUMSTEP);

    /* First term: first swaption  , centered on option strike - m */

    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, dStrike - dMargin, dVol,
                             dMaturity, 1.0, CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, dStrike - dMargin, dVol,
                             dMaturity, 1.0, CallPut, PREMIUM);

    dLvl = pow(1.0 + r, dNumPeriods);
    dInvLvl = r / (pow(1.0 + r, dDelay) * (1.0 - (1.0 / dLvl)));
    dTECw[0] = dFrequency * dInvLvl;

    if (dPaymentPower != 1)
      dAmt = (1 / dPaymentPower) *
             (pow((1 + dStrike), (1 / dPaymentPower - 1)) +
              (1 / dPaymentPower - 1) * dDeltaStrike *
                  pow((1 + dStrike), (1 / dPaymentPower - 2))) *
             dTECw[0];
    else
      dAmt = dTECw[0];

    dVal = dAmt * dOpt;

    /* Second term: Swaption centered on option strike - m + delta_k ( pos. or
     * neg.)*/
    if (iMethod)
      if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                               &dVol, &power))
        return err;

    if (VolType == SRT_NORMAL)
      dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);
    else
      dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                             1.0, CallPut, PREMIUM);

    r += dDeltaStrike / dFrequency;

    dInvLvlm1 = dTECw[0];
    dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
    dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
    dTECw[1] = 2.0 * (dInvLvl - dInvLvlm1);

    dAmt = 0;
    /* Avoid making useless loops */
    if (dPaymentPower != 1) {
      for (j = 0; j < 2; j++)
        dAmt += pow((1 + dStrike + (1 - j) * dDeltaStrike),
                    (1 / dPaymentPower - 2)) *
                dTECw[j];
      dAmt =
          (1 / dPaymentPower) * pow((1 + dStrike), (1 / dPaymentPower - 1)) *
              dTECw[1] +
          (1 / dPaymentPower) * (1 / dPaymentPower - 1) * dDeltaStrike * dAmt;
    } else
      dAmt = dTECw[1];

    dVal += dAmt * dOpt;

    /* Following Terms: Swaptions centered on Strike + i * Delta_k  (series) */
    i = 3;
    while ((bConverged == SRT_FALSE) && (i <= iNumSwaps)) {
      /* Get the Vol at strike K + (i - 1) * delta_k */
      if (iMethod == 2)
        if (err = swp_f_truncvol(szVolCurveName, dStart, dEnd, r * dFrequency,
                                 &dVol, &power))
          return err;

      if (iMethod == 1) {
        if (dSignRecPay *
                (r * dFrequency - pdStrikesVol[min(m, lNumStrikesInVol - 1)]) >=
            0.0)
          m += dSignRecPay;
        m = min(lNumStrikesInVol - 1 + iLeftRightPos, max(m, iLeftRightPos));
        dVol = pdVolVector[m - iLeftRightPos] +
               pdCoeffApprox[m - iLeftRightPos] *
                   (r * dFrequency - pdStrikesVol[m - iLeftRightPos]);
        if (r * dFrequency <= pdStrikesVol[0])
          dVol = pdVolVector[0];
      }

      if (VolType == SRT_NORMAL)
        dOpt = srt_f_optblknrm(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);
      else
        dOpt = srt_f_optblksch(dAdjForwardRate, r * dFrequency, dVol, dMaturity,
                               1.0, CallPut, PREMIUM);

      /* Move to the next strike (K + n * delta_k)
              at which we intersect the profile */
      r += dDeltaStrike / dFrequency;

      /* Coefficient from the series (recursive in 2 steps)
              to compute the weights */

      dInvLvlm2 = dInvLvlm1;
      dInvLvlm1 = dInvLvl;
      if (fabs(r) > 1.0e-10) {
        dInvLvl = r / (1.0 - (1.0 / (pow(1.0 + r, dNumPeriods))));
        dInvLvl = dFrequency * dInvLvl / (pow(1.0 + r, dDelay));
      } else {
        dInvLvl = 1.0 / dNumPeriods;
      }
      dTECw[i - 1] =
          (i * dInvLvl) - (2.0 * (i - 1) * dInvLvlm1) + ((i - 2) * dInvLvlm2);

      dAmt = 0;

      if (dPaymentPower != 1) {
        for (j = 0; j < i; j++)
          dAmt += pow((1 + dStrike + ((i - 1) - j) * dDeltaStrike),
                      (1 / dPaymentPower - 2)) *
                  dTECw[j];

        dAmt =
            (1 / dPaymentPower) * pow((1 + dStrike), (1 / dPaymentPower - 1)) *
                dTECw[i - 1] +
            (1 / dPaymentPower) * (1 / dPaymentPower - 1) * dDeltaStrike * dAmt;
      } else
        dAmt = dTECw[i - 1];

      dVal += dAmt * dOpt;

      /* fabs(dAmt * dOpt) < CMS_TOLERANCE */
      if (r * dFrequency > BOUNDARYINT) {
        /* Series has converged */
        bConverged = SRT_TRUE;
        break;
      }
      i++;
    }

    dSens = pow(1.0 + (dAdjForwardRate / dFrequency), dDelay);
    dSens1 = pow(1.0 + (dAdjForwardRate / dFrequency), dNumPeriods);
    dSens *= (1.0 - (1.0 / dSens1)) / dAdjForwardRate;

    dVal *= dSens;

    /* Free the dTECw vector */
    if (dTECw)
      free_dvector(dTECw, 0, (long)CMS_SPACE_NUMSTEP);
    dTECw = NULL;

    if (iMethod == 1) {
      free_dvector(pdCoeffApprox, 0, lNumStrikesInVol);
      free_dvector(pdVolVector, 0, lNumStrikesInVol);
      pdCoeffApprox = NULL;
      pdVolVector = NULL;
    }
  } /* END if(dMaturity != 0.0 )*/

  /* Return the Basis adjusted Rate */
  *dTECOptionValue = dVal / dRateConv;
  return err;
}

/*	--------------------------------------------------------------------
                        Alan: the following code is used in the direct
   integration
        -------------------------------------------------------------------- */

/* Level cash Numeraire used for direct integration */
double CMSNumeraire(double dForwardRate, double dDelay, double dFrequency,
                    double dNumPeriods) {
  double Y, IL, IB;

  Y = 1.0 + dForwardRate / dFrequency;
  IL = 1.0 - 1.0 / pow(Y, dNumPeriods);
  IB = pow(Y, dDelay);
  return IB * IL / dForwardRate;
}

double CMSdNumeraire(double dForwardRate, double dDelay, double dFrequency,
                     double dNumPeriods) {
  double Y, IL, IB;

  Y = 1.0 + dForwardRate / dFrequency;
  IL = 1.0 - 1.0 / pow(Y, dNumPeriods);
  IB = pow(Y, dDelay);
  return 1.0 / dForwardRate + (dNumPeriods - dDelay) / dFrequency / Y -
         dNumPeriods / dFrequency / Y / IL;
}

/*Pat version of Hermite points & Weights for Gaussian integration */
Err PatExtendedHermite(double *x, double *w, double xMax, int n) {
  double dx = .0;
  double dCoefa, dCoefb;
  double dSum = 0, dSum2 = 0, dSum4 = 0;
  int i;
  dx = xMax / (double)n;

  /* Init the Gaussian Points */
  x[0] = -n * dx;
  x[2 * n] = -x[0];
  x[1] = (1 - n) * dx;
  x[2 * n - 1] = -x[1];
  x[2] = (2 - n) * dx;
  x[2 * n - 2] = -x[2];

  /* Init the weights */
  w[0] = w[2 * n] = INV_SQRT_TWO_PI * 11 / 24 * exp(-x[0] * x[0] / 2.0) * dx;

  w[1] = w[2 * n - 1] = INV_SQRT_TWO_PI * exp(-x[1] * x[1] / 2.0) * dx;

  w[2] = w[2 * n - 2] =
      INV_SQRT_TWO_PI * 25 / 24 * exp(-x[2] * x[2] / 2.0) * dx;

  /* Init the Coefficient adjustments of the weights */
  dSum = w[0] + w[1] + w[2];
  dSum2 = x[0] * x[0] * w[0] + x[1] * x[1] * w[1] + x[2] * x[2] * w[2];

  dSum4 = x[0] * x[0] * x[0] * x[0] * w[0] + x[1] * x[1] * x[1] * x[1] * w[1] +
          x[2] * x[2] * x[2] * x[2] * w[2];

  /*Main computation of the hermite points */
  for (i = 3; i <= 2 * n - 3; i++) {
    x[i] = (i - n) * dx;
    w[i] = INV_SQRT_TWO_PI * exp(-x[i] * x[i] / 2.0) * dx;
    dSum += w[i];
    dSum2 += x[i] * x[i] * w[i];
    dSum4 += x[i] * x[i] * x[i] * x[i] * w[i];
  }

  /* Compute the Coefficient adjustments of the weights */
  dCoefa = (dSum4 - dSum2) / (dSum4 * dSum - dSum2 * dSum2);
  dCoefb = (dSum - dSum2) / (dSum4 * dSum - dSum2 * dSum2);

  for (i = 0; i <= 2 * n; i++)
    w[i] *= dCoefa + dCoefb * x[i] * x[i];

  return NULL;
}

/* A generic Version of the Tec / Cms payoff:
        if PayRec ==	SRT_PAYER is a Cms Cap  ,
                                        SRT_RECEIVER is a Cms Floor  ,
                                        (SrtReceiverType) 3 is a Cms Rate  ,
                                        (SrtReceiverType) 4 is a Tec Cap  ,
                                        (SrtReceiverType) 5 is a Tec Floor  ,
                                        (SrtReceiverType) 6 is a Tec Rate.*/

double TecCmsGenericPayoffDensity(double S, double dStrike,
                                  SrtReceiverType PayRec, double dMargin,
                                  double dTecPower) {
  switch (PayRec) {
  /* Cms Cap Case */
  case SRT_PAYER:
    if (S >= dStrike)
      return (S - dStrike);
    break;
  /* Cms Floor Case */
  case SRT_RECEIVER:
    if (S < dStrike)
      return (dStrike - S);
    break;
  /* Cms forward Case */
  case ((SrtReceiverType)3):
    return S;
    break;
  /* Tec Cap Case */
  case ((SrtReceiverType)4):
    if (S >= dStrike - dMargin)
      return pow((1.0 + S + dMargin), (1.0 / dTecPower)) -
             pow((1.0 + dStrike), (1.0 / dTecPower));
    break;
  /* Tec Floor Case */
  case ((SrtReceiverType)5):
    if (S < dStrike - dMargin)
      return pow((1.0 + dStrike), (1.0 / dTecPower)) -
             pow((1.0 + S + dMargin), (1.0 / dTecPower));
    break;
  /* Tec Forward case */
  case ((SrtReceiverType)6):
    return pow((1.0 + S + dMargin), (1.0 / dTecPower)) - 1.0;
    break;
  }
  return 0;
}

Err swp_f_FlatTecCmsOption(
    double dFwdSwapRate, Date StartDate, Date TheoEndDate, double dMaturity,
    double dNumPeriods, double dStrike, double dFrequency,
    SrtReceiverType PayRec,   /* Pay or Rec */
    double dDelay,            /* Delay */
    SrtDiffusionType VolType, /* Vol type Lognormal or Normal */
    double dFlatVol,          /* Flat Vol if used */
    double dRateConv,         /* Basis adjustment */
    double dMargin,           /* Margin for Tec */
    double dPower,            /* Power for Tec */
    double *dGenericValue) {
  static SRT_Boolean bComputeWeight = SRT_FALSE;
  static double *pdx, *pdWeight;
  Err err = NULL;
  double S, Numeraire, dxMax = 6.0, dAdjForwardRate, dAdjNumeraire = 0.0,
                       dStdDev;
  int i, NumHermitePoints = 250;

  /* compute only one time the value of the Integration grid points and weights
   */
  if (bComputeWeight == SRT_FALSE) {
    pdx = (double *)malloc((2 * NumHermitePoints + 1) * sizeof(double));
    pdWeight = (double *)malloc((2 * NumHermitePoints + 1) * sizeof(double));
    /*HermiteStandard(pdx  , pdWeight  , NumHermitePoints);*/
    PatExtendedHermite(pdx, pdWeight, dxMax, NumHermitePoints);
    bComputeWeight = SRT_TRUE;
  }

  /* Convert from delay in years to delay in periods */
  dDelay *= dFrequency;

  /* Correct the forward and strike to bring them to an annual basis */
  dAdjForwardRate = dFwdSwapRate * dRateConv;
  dStrike *= dRateConv;

  /* Compute the Std Dev of the underlying S for the Integration */
  dStdDev = dFlatVol * sqrt(dMaturity);

  *dGenericValue = 0.0;
  for (i = 0; i <= 2 * NumHermitePoints; i++) {
    /* Switch case for Normal / Lognormal case */
    if (VolType == SRT_NORMAL)
      S = dAdjForwardRate + dStdDev * pdx[i];
    else
      S = dAdjForwardRate * exp(-dStdDev * dStdDev / 2.0 + dStdDev * pdx[i]);

    /* Compute the numeraire for the current point of the grid */
    Numeraire = CMSNumeraire(S, dDelay, dFrequency, dNumPeriods);
    dAdjNumeraire += pdWeight[i] / Numeraire;

    /* Compute the payoff value for the current point of the grid */
    *dGenericValue +=
        pdWeight[i] *
        TecCmsGenericPayoffDensity(S, dStrike, PayRec, dMargin, dPower) /
        Numeraire;
  }

  /* Adjust the value of the expectation of the numeraire */
  *dGenericValue /= dAdjNumeraire * dRateConv;

  return err;
}

/* Generic (Normal or LogNormal) Caller of Black Scholes functions for the CMS
 */
double
srt_f_OptGeneric(double Fwd, double Strike, double Vol,
                 SrtDiffusionType VolType, /* Vol type Lognormal or Normal */
                 double Mat, SrtCallPutType CallPut) {
  if (VolType == SRT_NORMAL)
    return srt_f_optblknrm(Fwd, Strike, Vol, Mat, 1.0, CallPut, PREMIUM);
  return srt_f_optblksch(Fwd, Strike, Vol, Mat, 1.0, CallPut, PREMIUM);
}

Err swp_f_BnpCmsRate(
    double dFwdSwapRate, /*Forward Swap Rate */
    double dMaturity, double dNumPeriods, double dFrequency, double dDelay,
    double dVol, SrtDiffusionType VolType, int iMethod, double dlambda,
    double dCmsSpread,
    double *dBnpCmsRateValue) { /* See Cms Bnp Approximation for the poor
                                   notation of the variable */
  double dL = 0.0, dM = 0.0, dTemp = 0.0, dStdDev, dK, Y, dRate;

  Y = 1 + dFwdSwapRate / dFrequency;

  if (iMethod == 5)
    dStdDev = exp(dVol * dVol * dMaturity) - 1.0;
  else
    dStdDev = dVol * dVol * dMaturity;

  if (VolType == SRT_LOGNORMAL)
    dStdDev *= dFwdSwapRate * dFwdSwapRate;

  if (iMethod < 3)
    dK = -(dNumPeriods + 1) / 2.0 / dFrequency * Y;
  if (iMethod == 2) {
    dL = -(dNumPeriods + 2) * Y * dK / 3.0 / dFrequency;
    dM = -(dNumPeriods + 3) * Y * dL / 4.0 / dFrequency;
    dTemp = (3.0 * dK * dL - 3 * dM - dK * dK * dK) * dStdDev * dStdDev;
  }
  if ((iMethod == 3) || (iMethod == 4) || (iMethod == 5)) {
    dK = dNumPeriods / dFrequency / Y / (pow(Y, dNumPeriods) - 1) -
         1.0 / dFwdSwapRate;
  }
  if (iMethod == 4) {
    dTemp = dCmsSpread + dlambda * dNumPeriods * dStdDev * dVol * dVol *
                             dMaturity / dFrequency;
  }

  dRate = dFwdSwapRate - dK * dStdDev + dTemp;

  /* Apply delay adjustment assuming normal vol - E. Nahume  */
  if (VolType == SRT_LOGNORMAL)
    dStdDev *= dRate * dRate / dFwdSwapRate / dFwdSwapRate;
  *dBnpCmsRateValue = dRate - dDelay * dStdDev / (1 + dRate / dFrequency);

  return NULL;
}

void Implied_SpreadLambda(double dFwdSwapRate,  /*Forward Swap Rate */
                          double dCmsFullSmile, /*Forward Cms rate */
                          double dCmsVegaFS,    /* Vega full smile ???*/
                          double dMaturity, double dNumPeriods,
                          double dFrequency, double dVol,
                          SrtDiffusionType VolType, double *dlambda,
                          double *dCmsSpread) {
  double dK, Y, dStdDev;

  dStdDev = dVol * dVol * dMaturity;

  if (VolType == SRT_LOGNORMAL)
    dStdDev *= dFwdSwapRate * dFwdSwapRate;

  Y = 1 + dFwdSwapRate / dFrequency;
  dK = dNumPeriods / dFrequency / Y / (pow(Y, dNumPeriods) - 1) -
       1.0 / dFwdSwapRate;

  *dlambda = (dCmsVegaFS + 2 * dK * dStdDev / dVol) /
             (4.0 * dNumPeriods / dFrequency * dStdDev * dVol * dMaturity);

  *dCmsSpread =
      dCmsFullSmile - dFwdSwapRate + dK * dStdDev -
      *dlambda * dNumPeriods * dStdDev * dVol * dVol * dMaturity / dFrequency;
}

/* Functions to BOUND the domain of integration and Cap the vols */
/* SET FUNCTIONS */
void set_numberpoints(double numpoints, int type) {
  CMS_SPACE_NUMSTEP = numpoints;
  CMS_SPACE_TYPE = type;
}

void set_upperboundvol(double strike) { CAPSTRIKE = strike; }

void set_upperboundIntegral(double bound) { BOUNDARYINT = bound; }

/* Function to get the static variables */
double get_upperbound(void) { return BOUNDARYINT; }

double get_upperstrike(void) { return CAPSTRIKE; }

/* NEW FUNCTIONS for the CMS/TEC smile */

/* Part 2: the trunc getvol function */
Err swp_f_truncvol(char *vol_id, Ddate start, Ddate end, double strike,
                   double *volatility, double *power) {
  Err err = NULL;
  double adjstrike;

  /* cut the smile vol at CAPSTRIKE */
  if (strike > CAPSTRIKE)
    adjstrike = CAPSTRIKE;
  else
    adjstrike = strike;

  err = swp_f_vol(vol_id, start, end, adjstrike, volatility, power);
  return err;
}

void CMSTECinitParams(CMSParams *Paramstruct) {
  Paramstruct->CAP_IntUD = 0.5;
  Paramstruct->CAP_VolUD = 0.2;
  Paramstruct->CAP_Mesh = 2;

  Paramstruct->FLOOR_IntUD = 0.00005;
  Paramstruct->FLOOR_VolUD = 0.01;
  Paramstruct->FLOOR_Mesh = 2;

  Paramstruct->Nx = 41;
}

void CMSTECinitvol(CMSVol *Volstruct) {
  // initialise
  Volstruct->voltype_input = SRT_LOGNORMAL;
  Volstruct->voltype_used = SRT_LOGNORMAL;

  Volstruct->input_method = 0;
  Volstruct->comp_method = 0;
  Volstruct->szVolCurveName = NULL;
  Volstruct->input_SABR = 0;
  Volstruct->input_SABRAF = 0;
  Volstruct->alpha = 0.0;
  Volstruct->beta = 1.0;
  Volstruct->rho = 0.0;
  Volstruct->zeta = 0.0;
  Volstruct->lambda = 0.0;
  Volstruct->vol = 0.1;
  Volstruct->volinfty = 0.0;
  Volstruct->atm = 0;
  Volstruct->bump_type = 0;
  Volstruct->bump_vol = 0.0;
  Volstruct->is_vol_bumped = 0;
}

void CMSTECfreevol(CMSVol *Volstruct) {
  if (Volstruct->szVolCurveName)
    free(Volstruct->szVolCurveName);
  Volstruct->szVolCurveName = NULL;
}

Err CMSTECsetupvol(CMSVol *Volstruct, CMSGreeks *Greekstruct, double strike,
                   double swaprate, double maturity, SrtReceiverType PayRec) {
  Err err = NULL;
  double power;
  double yesno;
  int index;
  SABR_VOL_TYPE sabrtype;
  double tempvol;
  int isSABR = 0, isSABRAF = 0;

  /* Get the inputs with no changes */
  if (Volstruct->input_method == 4) {
    double strikeused;
    if (Volstruct->atm == 1)
      strikeused = swaprate;
    else
      strikeused = strike;

    err = swp_f_vol(Volstruct->szVolCurveName, Volstruct->start, Volstruct->end,
                    strikeused, &(Volstruct->vol), &power);
    if (err)
      return err;

    Volstruct->comp_method = 0;
  } else if (Volstruct->input_method == 2 || Volstruct->input_method == 3 ||
             Volstruct->input_method == 250) // 3 for BMM in test
  {
    // Detects if the market is SABR by retrieving the sigmabetavol.
    /* old method
    err = swp_f_SABRvol(Volstruct->szVolCurveName  , Volstruct->start  ,
    Volstruct->end  , swaprate  , &yesno  ,  &power  , SABR_BETAVOL);
    */
    err = swp_f_IsSmileVol(Volstruct->szVolCurveName, &yesno);

    // TODO: Joe to move this to our wrapper function
    if (err || (fabs(yesno) < 1e-8)) {
      isSABR = 0;
      isSABRAF = 0;
    } else {
      err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                          Volstruct->end, swaprate, &yesno, &power, SABR_ZETA);
      if (err) {
        isSABR = 1;
        isSABRAF = 0;
      } else {
        isSABR = 0;
        isSABRAF = 1;
      }
    }

    if (isSABR || isSABRAF) {
      if (Volstruct->voltype_input == SRT_LOGNORMAL)
        index = 0;
      else if (Volstruct->voltype_input == SRT_NORMAL)
        index = 1;
      else /* SRT_BETAVOL */
        index = 2;

      if (err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                              Volstruct->end, swaprate, &(Volstruct->vol),
                              &power, index))
        return err;

      if (err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                              Volstruct->end, swaprate, &(Volstruct->beta),
                              &power, 6))
        return err;

      if (err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                              Volstruct->end, swaprate, &(Volstruct->alpha),
                              &power, 5))
        return err;

      if (err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                              Volstruct->end, swaprate, &(Volstruct->rho),
                              &power, 7))
        return err;

      if (isSABR) {
        Volstruct->comp_method = 0;
        Volstruct->input_SABR = 1;
      } else {
        if (err = swp_f_SABRvol(Volstruct->szVolCurveName, Volstruct->start,
                                Volstruct->end, swaprate, &(Volstruct->zeta),
                                &power, 9))
          return err;
        Volstruct->comp_method = 0;
        Volstruct->input_SABR = 0;
        Volstruct->input_SABRAF = 1;
      }
    } else

      Volstruct->comp_method = 2;
  }

  /* Apply the bump's if needed */
  Volstruct->is_vol_bumped = 0;
  if (Greekstruct) {
    if (Greekstruct->greek_type) {
      if ((fabs(Greekstruct->alpha_bump) > 0.0) &&
          ((Greekstruct->greek_type == 1) || (Greekstruct->greek_type == 3)))
        Volstruct->alpha = Volstruct->alpha + Greekstruct->alpha_bump;

      if ((fabs(Greekstruct->rho_bump) > 0.0) &&
          ((Greekstruct->greek_type == 1) || (Greekstruct->greek_type == 5)))
        Volstruct->rho = Volstruct->rho + Greekstruct->rho_bump;

      if ((fabs(Greekstruct->beta_bump) > 0.0) &&
          ((Greekstruct->greek_type == 1) || (Greekstruct->greek_type == 4)))
        Volstruct->beta = Volstruct->beta + Greekstruct->beta_bump;

      if ((fabs(Greekstruct->vol_bump) > 0.0) &&
          ((Greekstruct->greek_type == 1) || (Greekstruct->greek_type == 2))) {
        if (Volstruct->comp_method == 2) {
          /* only used for the full smile case with strike's vols */
          Volstruct->bump_type = Greekstruct->vol_bump_type;
          Volstruct->bump_vol = Greekstruct->vol_bump;
          Volstruct->is_vol_bumped = 1;
        } else {
          if (Greekstruct->vol_bump_type == 0)
            Volstruct->vol = Volstruct->vol + Greekstruct->vol_bump;
          else
            Volstruct->vol = Volstruct->vol * Greekstruct->vol_bump;
        }
      }
    }
  }

  /* Convert the inputs +bumped into more practical parameters */
  if (Volstruct->input_SABR == 1) {
    /* convert the vol into a Betavol */
    if (Volstruct->voltype_input == SRT_LOGNORMAL)
      sabrtype = SABR_ATM_LOG;
    else if (Volstruct->voltype_input == SRT_NORMAL)
      sabrtype = SABR_ATM_NORM;
    else
      sabrtype = SABR_ATM_BETA;

    err = vol_conv(Volstruct->vol, sabrtype, &tempvol, SABR_ATM_BETA, swaprate,
                   swaprate, maturity, Volstruct->alpha, Volstruct->beta,
                   Volstruct->rho);
    if (err)
      return err;

    if (Volstruct->voltype_input == SRT_BMMVOL)
      Volstruct->voltype_used = SRT_BMMVOL;
    else if (Volstruct->voltype_input == SRT_BMM2VOL)
      Volstruct->voltype_used = SRT_BMM2VOL;
    else
      Volstruct->voltype_used = SRT_BETAVOL;
    Volstruct->vol = tempvol;
  } else if (Volstruct->input_SABRAF == 1) {
    /* convert the vol into a SABRAF Betavol */
    if (Volstruct->voltype_input == SRT_LOGNORMAL)
      sabrtype = SRT_LOGNORMAL;
    else if (Volstruct->voltype_input == SRT_NORMAL)
      sabrtype = SRT_NORMAL;
    else
      sabrtype = SRT_BETAVOL;

    /*	Calculate equivalent vol of the type to be frozen	only ATM types or
     * BETA*/
    err = srt_f_optbmm2vol(swaprate, swaprate, maturity, Volstruct->vol,
                           Volstruct->alpha, Volstruct->beta, Volstruct->rho,
                           Volstruct->zeta, sabrtype, SRT_BETAVOL, &tempvol);
    if (err)
      return err;
    Volstruct->vol = tempvol;
    Volstruct->voltype_used = SRT_BMM2VOL;
  }

  return err;
}

/* 3: usefull functions for the integration */
/* 3.1 numeraire associated to the default cash settle swaption of the market */
void CMS_get_numeraire(double x, /* the swap */
                       double nfp, double invf, double *num, double *d1num,
                       double *d2num) {
  double Y, y;
  double d1Y, d2Y;

  y = 1.0 + x * invf;
  Y = pow(y, -nfp);
  y = 1.0 / y;

  d1Y = -nfp * Y * y * invf;
  d2Y = -(nfp + 1.) * d1Y * y * invf;

  Y = 1.0 / (1.0 - Y);

  /* numeraire */
  (*num) = x * Y;
  /* first derivative of the numeraire */
  (*d1num) = 1.0 + x * Y * d1Y;
  (*d1num) *= Y;
  /* second derivative of the numeraire */
  (*d2num) = (1.0 + x * Y * d1Y) * 2.0 * d1Y + x * d2Y;
  (*d2num) *= Y * Y;
}

/* 3.2 discount factor (delay) and its derivatives */
void CMS_get_discount(double x, double cvg, double spd, double *ZC,
                      double *d1ZC, double *d2ZC) {
  double df;

  // The following expression gives a better fit to the previous
  // swp_CMS_Option function when the payment delay is big
  // than the original 1 / (1+cvg*S+spd)

  df = 1.0 / pow(1.0 + x + spd, cvg);
  *ZC = df;
  *d1ZC = -cvg * df / (1.0 + x + spd);
  *d2ZC = -(cvg + 1.0) * (*d1ZC) / (1.0 + x + spd);

  // Initial setup of the function
  /*
          df = 1.0 / (1.0 + cvg * x + spd);
          // Fwd ZC
          *ZC = df;
          // first derivative
          *d1ZC = - cvg * df * df;
          // second derivative
          *d2ZC = - 2.0 * cvg * df * (*d1ZC);
  */
}

/* 3.3 the payoff and its derivatives */
void CMS_get_payoff(double x, double K, /* real strike */
                    double invr, double invp, double m, double cp, double *pay,
                    double *d1pay, double *d2pay) {
  double y;
  double Y;
  double temp;
  double tempay;

  if (invp == 1.0) {
    temp = cp * (x * invr + m - K);
    if (temp > 0) {
      *pay = temp;
      *d1pay = cp;
    } else {
      *pay = 0.0;
      *d1pay = 0.0;
    }
    *d2pay = 0.0;
  } else {
    y = 1.0 + x * invr + m;
    Y = cp * pow(y, invp);
    tempay = Y - cp * pow(1.0 + K, invp);
    if (tempay > 0.0) {
      y = 1.0 / y;
      *pay = tempay;
      temp = Y * invr * invp * y;
      *d1pay = temp;
      *d2pay = temp * invr * (invp - 1.0) * y;
    } else {
      *pay = 0.0;
      *d1pay = 0.0;
      *d2pay = 0.0;
    }
  }
}

/* 3.4 : the reduced payoff second derivative */
double CMS_get_reduced_payoff2(double x, double nfp, double invf,
                               double K, /* real strike */
                               double invp, double invr, double m, double spdZC,
                               double cp, double cvgZC) {
  double num, d1num, d2num;
  double ZC, d1ZC, d2ZC;
  double pay, d1pay, d2pay;

  CMS_get_numeraire(x, nfp, invf, &num, &d1num, &d2num);
  CMS_get_payoff(x, K, invr, invp, m, cp, &pay, &d1pay, &d2pay);

  if (cvgZC > 10e-3) {
    CMS_get_discount(x, cvgZC, spdZC, &ZC, &d1ZC, &d2ZC);
    return ZC * pay * d2num + ZC * d2pay * num + d2ZC * pay * num +
           2.0 * (ZC * d1pay * d1num + d1ZC * pay * d1num + d1ZC * d1pay * num);
  } else {
    return pay * d2num + d2pay * num + 2.0 * d1pay * d1num;
  }
}

/* 3.4 : the full payoff */
double CMS_get_pay(double x, double K, /*real strike */
                   double invp, double invr, double m, double cp) {
  double Temppay;
  double yield;

  if (invp == 1.0) {
    yield = cp * (x * invr + m - K);
    Temppay = yield > 0.0 ? yield : 0.0;
  } else {
    yield = 1.0 + x * invr + m;
    yield = cp * (pow(yield, invp) - pow(1.0 + K, invp));
    Temppay = yield > 0.0 ? yield : 0.0;
  }

  return Temppay;
}

/* 3.5 get the factor of the control variate */
double CMS_get_controlvariatefactor(double strikeadj, double nfp, double invr,
                                    double m, double invp, double invf,
                                    double cvgZC, double spdZC) {
  double tempDisc = 1.0;
  double tempLvl;
  double tempDer;
  double yield;

  if (invp == 1.0) {
    tempDer = invr;
  } else {
    yield = 1.0 + strikeadj * invr + m;
    tempDer = invr * invp * pow(yield, invp - 1.0);
  }

  yield = 1.0 + strikeadj * invf;
  tempLvl = (1.0 - pow(yield, -nfp)) / strikeadj;

  if (cvgZC > 1.0e-3) {
    // The following expression gives a better fit to the previous
    // swp_CMS_Option function when the payment delay is big
    // than the original 1 / (1+cvg*S+spd)
    tempDisc = 1.0 / pow(1.0 + strikeadj + spdZC, cvgZC);

    // Original setup of the function where the delay
    // is approximated by this linear expression
    // tempDisc = 1.0 / (1.0 + cvgZC * strikeadj + spdZC);
  }

  return tempDisc * tempDer / tempLvl;
}

void CMS_get_lastpart(double x, double forward, double vol, double maturity,
                      double cp, double *m0, double *m1, double *m2) {
  double d;
  double stdev;

  stdev = vol * sqrt(maturity);
  d = log(forward / x) / stdev - 0.5 * stdev;
  *m0 = norm(cp * d);
  *m1 = forward * norm(cp * (d + stdev));
  *m2 = forward * forward * norm(cp * (d + 2.0 * stdev)) * exp(stdev * stdev);
}

/* build the points for the integration */
void CMS_get_pdX(double strikeadj, long Nx, double IntUD, double cp, int Mesh,
                 double *pX, double *pdX) {
  long i;
  double spacing;
  double prev;
  double next;

  /* mid points */
  if (Mesh == 0) {
    spacing = (IntUD - strikeadj) / Nx;
    pX[0] = strikeadj + spacing / 2.0;
    pdX[0] = spacing;
    for (i = 1; i < Nx; i++)
      pX[i] = pX[i - 1] + spacing;
  } else
      /* exponential */
      if (Mesh == 1) {
    spacing = exp(log(IntUD / strikeadj) / Nx);
    prev = strikeadj;
    for (i = 0; i < Nx; i++) {
      next = prev * spacing;
      pdX[i] = (next - prev);
      pX[i] = (next + prev) / 2.0;
      prev = next;
    }
  } else
      /* Gauss legendre distributed */
      if (Mesh == 2) {
    GaussLeg(strikeadj, IntUD, pX - 1, pdX - 1, Nx);
  }
}

/* compute the weights for the integration */
void swp_f_FinalCms_Weights(                      /* underlying */
                            double *x,            /* point of computation */
                            long nx, double dnfp, /* number of full period of
                                                     the default swap rate */
                            double dfreq, /* frequency of the default swap rate
                                             1:A  , 2:S  , 4:Q */
                            double dRateConv, /* adjustment swap default/swap
                                                 und */
                            double dPower,    /* 1: CMS  , 4: tec  , etc.... */
                            double dMargin,
                            /* option */
                            double dRealStrike, double cp,
                            /* delay parameters */
                            double dCvgZC, /* delay in nfp */
                            double dspdZC, /* spread for ZC */
                            /* Output */
                            double *dCMSWeights) {
  double invf = 1.0 / dfreq;
  double invr = 1.0 / dRateConv;
  double invp = 1.0 / dPower;
  long n;

  for (n = 0; n < nx; n++)
    dCMSWeights[n] = CMS_get_reduced_payoff2(
        x[n], dnfp, invf, dRealStrike, invp, invr, dMargin, dspdZC, cp, dCvgZC);
}

/* compute the weitghts for the integration */
Err swp_f_FinalCms_Opts(           /* underlying */
                        double *x, /* point of computation */
                        long nx,
                        /* option */
                        double dDefaultSwap, double dStrikeAdj, double cp,
                        double dMaturity,
                        /* Flat */
                        CMSVol *Volstruct,
                        /* Numerical params */
                        double VolUD,
                        /* Output */
                        double *dCMSOpts) {
  Err err = NULL;
  long n;
  double volused = Volstruct->vol;
  double power;
  double betavol = 0.0;
  double alpha = 0.0;
  double beta = 0.0;
  double rho = 0.0;
  double zeta = 0.0;
  double strikeused = dStrikeAdj;
  // double *CalibStrikes = NULL;
  // double *CalibVols	 = NULL;
  // double atmvol  ,newbetavol  , newbeta  , newalpha  , newrho  , newprob;
  // int i;
  double Fwd1, Fwd2, Sig1, Sig2, Pi;
  double res;

  SrtCallPutType call_put;

  /* set the call or the put type function */
  if (cp > 0)
    call_put = SRT_CALL;
  else
    call_put = SRT_PUT;

  /* for speed reasons and in order to be ready for the Heston vol  , we sort by
   * method */
  if (Volstruct->comp_method == 0) {
    volused = Volstruct->vol;
    if (Volstruct->voltype_used == SRT_LOGNORMAL) {
      for (n = 0; n < nx; n++)
        dCMSOpts[n] = srt_f_optblksch(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
    } else if (Volstruct->voltype_used == SRT_NORMAL) {
      for (n = 0; n < nx; n++)
        dCMSOpts[n] = srt_f_optblknrm(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
    } else if (Volstruct->voltype_used == SRT_BETAVOL) {
      betavol = Volstruct->vol;
      beta = Volstruct->beta;
      alpha = Volstruct->alpha;
      rho = Volstruct->rho;

      /* in case of the Floor */
      err = srt_f_optbetastochvoltoblkvol(dDefaultSwap, VolUD, betavol, alpha,
                                          rho, dMaturity, beta, &volused);
      if (err)
        return err;

      for (n = 0; n < nx; n++) {
        if ((cp * (VolUD - x[n]) >= -10e-17)) {
          err =
              srt_f_optbetastochvoltoblkvol(dDefaultSwap, x[n], betavol, alpha,
                                            rho, dMaturity, beta, &volused);
          if (err)
            return err;
        }
        dCMSOpts[n] = srt_f_optblksch(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
      }
    } else if (Volstruct->voltype_used == SRT_HESTONVOL) {
      /*
      err = HestonPrice(
                      dDefaultSwap  ,
                      x  ,
                      nx  ,
                      dMaturity  ,
                      Volstruct->vol  ,
                      Volstruct->alpha  ,
                      Volstruct->volinfty  ,
                      double	b  ,
                      double  Gamma  ,
                      Volstruct->rho  ,
                      1.0  ,
                      double  UpperBound  ,
                      call_put  ,
                      SRT_PREMIUM  ,
                      SRT_Boolean isVolInfFix  ,
                      int IntegrType  ,
                      int nSteps  ,
                      dCMSOpts
                      );
      */
      if (err)
        return err;
    } else if (Volstruct->voltype_used == SRT_BMM2VOL) {
      /* beta mix still under tests*/
      betavol = Volstruct->vol;
      beta = Volstruct->beta;
      alpha = Volstruct->alpha;
      rho = Volstruct->rho;

      if (Volstruct->input_SABRAF == 1) {
        zeta = Volstruct->zeta;
        BMM2GetStates(dDefaultSwap, dMaturity, betavol, alpha, beta, rho, zeta,
                      &Sig1, &Sig2, &Fwd1, &Fwd2, &Pi);
      } else {
        BMM2CalibOnSabrStates(dDefaultSwap, dMaturity, betavol, beta, alpha,
                              rho, Volstruct->Pi, Volstruct->NstD, &Fwd1, &Fwd2,
                              &Sig1, &Sig2, &Pi, SRT_BETAVOL, &res);
        if (res > 1e-3)
          return "SABRAF could'nt calibrate";
      }

      srt_f_optbmmvolfromstates(VolUD, dMaturity, beta, Fwd1, Fwd2, Sig1, Sig2,
                                Pi, SRT_LOGNORMAL, &volused);
      if (err)
        return err;

      for (n = 0; n < nx; n++) {
        if ((cp * (VolUD - x[n]) >= -10e-17))
          srt_f_optbmmvolfromstates(x[n], dMaturity, beta, Fwd1, Fwd2, Sig1,
                                    Sig2, Pi, SRT_LOGNORMAL, &volused);
        dCMSOpts[n] = srt_f_optblksch(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
      }
    } else if (Volstruct->voltype_used == SRT_BMMVOL) {
      /* beta mix still under tests*/
      betavol = Volstruct->vol;
      beta = Volstruct->beta;
      alpha = Volstruct->alpha;
      rho = Volstruct->rho;

      BMMCalibOnSabrStates(dDefaultSwap, dMaturity, betavol, beta, alpha, rho,
                           Volstruct->Pi, Volstruct->NstD, &Fwd1, &Fwd2, &Sig1,
                           &Sig2, &Pi, SRT_BETAVOL, &res);
      if (res > 1e-3)
        return "BMM could'nt calibrate";

      srt_f_optbmmvolfromstates(VolUD, dMaturity, beta, Fwd1, Fwd2, Sig1, Sig2,
                                Volstruct->Pi, SRT_LOGNORMAL, &volused);
      if (err)
        return err;

      for (n = 0; n < nx; n++) {
        if ((cp * (VolUD - x[n]) >= -10e-17))
          srt_f_optbmmvolfromstates(x[n], dMaturity, beta, Fwd1, Fwd2, Sig1,
                                    Sig2, Volstruct->Pi, SRT_LOGNORMAL,
                                    &volused);
        dCMSOpts[n] = srt_f_optblksch(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
      }
    } else
      return "CMS :: vol type unknown";

  } else if (Volstruct->comp_method == 2) {
    err = swp_f_vol(Volstruct->szVolCurveName, Volstruct->start, Volstruct->end,
                    VolUD, &volused, &power);
    if (err)
      return err;

    if (Volstruct->is_vol_bumped == 1) {
      if (Volstruct->bump_type == 0)
        volused += Volstruct->bump_vol;
      else
        volused *= Volstruct->bump_vol;
    }

    for (n = 0; n < nx; n++) {
      /* update the vol up-down to the strikevolbound */
      if ((cp * (VolUD - x[n]) >= -10e-17)) {
        err = swp_f_vol(Volstruct->szVolCurveName, Volstruct->start,
                        Volstruct->end, x[n], &volused, &power);
        if (err)
          return err;

        if (Volstruct->is_vol_bumped == 1) {
          if (Volstruct->bump_type == 0)
            volused += Volstruct->bump_vol;
          else
            volused *= Volstruct->bump_vol;
        }
      }

      if (Volstruct->voltype_used == SRT_LOGNORMAL) {
        dCMSOpts[n] = srt_f_optblksch(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
      } else if (Volstruct->voltype_used == SRT_NORMAL) {
        dCMSOpts[n] = srt_f_optblknrm(dDefaultSwap, x[n], volused, dMaturity,
                                      1.0, call_put, SRT_PREMIUM);
      } else
        return "CMS :: vol type unknown";
    }

  } else
    return "CMS :: method unknown";

  return err;
}

/* Main function for TEC-CMS pricing */
Err swp_f_CMSTEC_Option(                     /* underlying */
                        double dFwdSwapRate, /*Forward Swap Rate */
                        double dnfp,  /* number of full period of the default
                                         swap rate */
                        double dfreq, /* frequency of the default swap rate 1:A
                                         , 2:S  , 4:Q */
                        double dRateConv, /* adjustment swap default/swap und */
                        /* Option parameters */
                        double dMaturity,       /* CMS Maturity as double */
                        double dStrike,         /* Option Strike */
                        SrtReceiverType PayRec, /* Pay or Rec */
                        double dPower, /* 1: CMS  , 4: tec  , etc.... */
                        double dMargin,
                        /* delay parameters */
                        double dCvgZC, /* cvg of the delay */
                        double dFwdZC, /* B(0  ,Tpay) / B(0  ,Tstart) */
                        /* vol struct */
                        CMSVol *Volstruct, /* Flat Vol if used */
                        /* greek struct */
                        CMSGreeks *Greekstruct,
                        /* numerical param */
                        CMSParams *params,
                        /* Output */
                        double *dCMSOptionValue) {
  Err err = NULL;
  /* Integration param */
  long i;
  double dx;
  double dIntegral;
  double left;
  double right;
  /* last part */
  double factor;
  double dTV;
  /* Temp param */
  double dDefaultSwap;
  double dNumeraireFactor, num1, num2;
  /* Option param */
  double cp = 1.0;
  double invr;
  double invp;
  double invf;
  double strikeadj;
  /* ZC param */
  double cvgZC = dCvgZC;
  double spdZC = 0.0;
  /* Common params */
  double vol = 0.0;
  double power;
  /* Grid  , weights and Options */
  double *pdX = NULL;
  double *pdDX = NULL;
  double *pdW = NULL;
  double *pdOpts = NULL;

  double VolUD;
  double IntUD;
  int Mesh;
  double Fwd1, Fwd2, Sig1, Sig2, Pi, res;

  /* Initialise param */
  invp = 1.0 / dPower;
  invr = 1.0 / dRateConv;
  invf = 1.0 / dfreq;
  strikeadj = (dStrike - dMargin) * dRateConv;
  if (PayRec == SRT_RECEIVER)
    cp = -1.0;
  dDefaultSwap = dFwdSwapRate * dRateConv;

  if (PayRec == SRT_PAYER) {
    VolUD = params->CAP_VolUD;
    IntUD = params->CAP_IntUD;
    Mesh = params->CAP_Mesh;
  } else {
    VolUD = params->FLOOR_VolUD;
    IntUD = params->FLOOR_IntUD;
    Mesh = params->FLOOR_Mesh;
  }

  err = CMSTECsetupvol(Volstruct, Greekstruct, strikeadj, dDefaultSwap,
                       dMaturity, PayRec);
  if (err)
    return err;

  /* ZC params */
  if (cvgZC > 10e-4)
    // Original setup of the function where ZCFwd = 1 / (1 + cvg*S + spd)
    // spdZC =  ( 1.0 / dFwdZC  - cvgZC * dDefaultSwap - 1.0);

    // The following expression gives a better fit to the previous
    // swp_CMS_Option function when the payment delay is big
    // than the original 1 / (1+cvg*S+spd)
    spdZC = pow(1.0 / dFwdZC, 1.0 / cvgZC) - dDefaultSwap - 1.0;

  /* compute the integration space grid */
  pdDX = (double *)malloc(2 * params->Nx * sizeof(double));
  pdX = (double *)malloc(2 * params->Nx * sizeof(double));

  if (cp == 1) {
    left = strikeadj;
    right = VolUD;
  } else {
    left = IntUD;
    right = VolUD;
  }

  CMS_get_pdX(left, params->Nx, right, cp, Mesh, pdX, pdDX);
  if (cp == 1) {
    left = VolUD;
    right = IntUD;
  } else {
    left = VolUD;
    right = strikeadj;
  }
  CMS_get_pdX(left, params->Nx, right, cp, Mesh, pdX + params->Nx,
              pdDX + params->Nx);

  /* compute the weights */
  pdW = (double *)malloc(2 * params->Nx * sizeof(double));
  swp_f_FinalCms_Weights(pdX, 2 * params->Nx, dnfp, dfreq, dRateConv, dPower,
                         dMargin, dStrike, cp, cvgZC, spdZC, pdW);

  /* compute the options */
  pdOpts = (double *)malloc(2 * params->Nx * sizeof(double));
  err = swp_f_FinalCms_Opts(pdX, 2 * params->Nx, dDefaultSwap, strikeadj, cp,
                            dMaturity, Volstruct, VolUD,
                            /* Output */
                            pdOpts);
  if (err)
    goto FREE_CMS_END;

  /* do the integration */
  dIntegral = 0.0;
  if (Mesh == 0) {
    dx = pdDX[0];
    for (i = 0; i < params->Nx; i++)
      dIntegral += pdW[i] * pdOpts[i] * dx;

    dx = pdDX[params->Nx];
    for (i = params->Nx; i < 2 * params->Nx; i++)
      dIntegral += pdW[i] * pdOpts[i] * dx;
  } else {
    for (i = 0; i < 2 * params->Nx; i++)
      dIntegral += pdDX[i] * pdW[i] * pdOpts[i];
  }

  /* add the TV of the control variate factor */
  factor = CMS_get_controlvariatefactor(strikeadj, dnfp, invr, dMargin, invp,
                                        invf, cvgZC, spdZC);

  /* Get the vol */
  if (Volstruct->comp_method == 0) {
    if ((Volstruct->voltype_used == SRT_NORMAL) ||
        (Volstruct->voltype_used == SRT_LOGNORMAL))
      vol = Volstruct->vol;
    else if (Volstruct->voltype_used == SRT_BMM2VOL) {
      if (Volstruct->input_SABRAF == 1) {
        BMM2GetStates(dDefaultSwap, dMaturity, Volstruct->vol, Volstruct->alpha,
                      Volstruct->beta, Volstruct->rho, Volstruct->zeta, &Sig1,
                      &Sig2, &Fwd1, &Fwd2, &Pi);
      } else {
        BMM2CalibOnSabrStates(dDefaultSwap, dMaturity, Volstruct->vol,
                              Volstruct->beta, Volstruct->alpha, Volstruct->rho,
                              Volstruct->Pi, Volstruct->NstD, &Fwd1, &Fwd2,
                              &Sig1, &Sig2, &Pi, SRT_BETAVOL, &res);
        if (res > 1e-3)
          return "SABRAF could'nt calibrate";
      }

      srt_f_optbmmvolfromstates(strikeadj, dMaturity, Volstruct->beta, Fwd1,
                                Fwd2, Sig1, Sig2, Pi, SRT_LOGNORMAL, &vol);
      if (err)
        return err;
    } else if (Volstruct->voltype_used == SRT_BMMVOL) {
      BMMCalibOnSabrStates(dDefaultSwap, dMaturity, Volstruct->vol,
                           Volstruct->beta, Volstruct->alpha, Volstruct->rho,
                           Volstruct->Pi, Volstruct->NstD, &Fwd1, &Fwd2, &Sig1,
                           &Sig2, &Pi, SRT_BETAVOL, &res);
      if (res > 1e-3)
        return " BMM could'nt calibrate";

      srt_f_optbmmvolfromstates(strikeadj, dMaturity, Volstruct->beta, Fwd1,
                                Fwd2, Sig1, Sig2, Volstruct->Pi, SRT_LOGNORMAL,
                                &vol);
      if (err)
        return err;
    } else {
      err = srt_f_optbetastochvoltoblkvol(
          dDefaultSwap, strikeadj, Volstruct->vol, Volstruct->alpha,
          Volstruct->rho, dMaturity, Volstruct->beta, &vol);
      if (err)
        return err;
    }
  } else {
    err = swp_f_vol(Volstruct->szVolCurveName, Volstruct->start, Volstruct->end,
                    strikeadj, &vol, &power);
    if (err)
      return err;

    if (Volstruct->is_vol_bumped == 1) {
      if (Volstruct->bump_type == 0)
        vol += Volstruct->bump_vol;
      else
        vol *= Volstruct->bump_vol;
    }
  }

  /* compute the TV according to the Vol type definition */
  if (Volstruct->voltype_used == SRT_NORMAL)
    dTV = srt_f_optblknrm(dDefaultSwap, strikeadj, vol, dMaturity, 1.0,
                          ((cp > 0.0) ? SRT_CALL : SRT_PUT), SRT_PREMIUM);
  else
    dTV = srt_f_optblksch(dDefaultSwap, strikeadj, vol, dMaturity, 1.0,
                          ((cp > 0.0) ? SRT_CALL : SRT_PUT), SRT_PREMIUM);

  dIntegral += factor * dTV;

  /* convert into the Qtpay numeraire */
  CMS_get_numeraire(dDefaultSwap, dnfp, invf, &dNumeraireFactor, &num1, &num2);
  dIntegral /= dNumeraireFactor;
  if (cvgZC > 1.0e-3)
    dIntegral /= dFwdZC;

FREE_CMS_END:

  if (pdX)
    free(pdX);
  if (pdDX)
    free(pdDX);
  if (pdW)
    free(pdW);
  if (pdOpts)
    free(pdOpts);
  pdX = NULL;
  pdDX = NULL;
  pdW = NULL;
  pdOpts = NULL;

  /* Return the Basis adjusted Rate */
  *dCMSOptionValue = dIntegral;

  return err;
}

Err swp_f_CMSTEC_rate(                     /* underlying */
                      double dFwdSwapRate, /*Forward Swap Rate */
                      double dnfp,  /* number of full period of the default swap
                                       rate */
                      double dfreq, /* frequency of the default swap rate 1:A  ,
                                       2:S  , 4:Q */
                      double dRateConv, /* adjustment swap default/swap und */
                      /* Option parameters */
                      double dMaturity, /* CMS Maturity as double */
                      double dPower,    /* 1: CMS  , 4: tec  , etc.... */
                      double dMargin,
                      /* delay parameters */
                      double dCvgZC, /* cvg of the delay */
                      double dFwdZC, /* B(0  ,Tpay) / B(0  ,Tstart) */
                      /* vol structure */
                      CMSVol *Volstruct, /* Flat Vol*/
                      /* Greek struct */
                      CMSGreeks *Greekstruct,
                      /* numerical param */
                      CMSParams *params,
                      /* Output */
                      double *dCMSTECrate) {
  Err err = NULL;
  double call;
  double put;

  err = swp_f_CMSTEC_Option(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                            dFwdSwapRate, SRT_PAYER, dPower, dMargin, dCvgZC,
                            dFwdZC, Volstruct, Greekstruct, params, &call);
  if (err)
    return err;

  err = swp_f_CMSTEC_Option(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                            dFwdSwapRate, SRT_RECEIVER, dPower, dMargin, dCvgZC,
                            dFwdZC, Volstruct, Greekstruct, params, &put);
  if (err)
    return err;

  /* return the rate inside the decompouding formula */
  *dCMSTECrate = pow(call - put + pow(1.0 + dFwdSwapRate, 1 / dPower), dPower) -
                 1.0 - dMargin;

  return err;
}

Err swp_f_CMSTEC_GREEKS(                     /* underlying */
                        double dFwdSwapRate, /*Forward Swap Rate */
                        double dnfp,  /* number of full period of the default
                                         swap rate */
                        double dfreq, /* frequency of the default swap rate 1:A
                                         , 2:S  , 4:Q */
                        double dRateConv, /* adjustment swap default/swap und */
                        /* Option parameters */
                        double dMaturity,       /* CMS Maturity as double */
                        double dStrike,         /* Option Strike */
                        SrtReceiverType PayRec, /* Pay or Rec */
                        double dPower, /* 1: CMS  , 4: tec  , etc.... */
                        double dMargin,
                        /* delay parameters */
                        double dCvgZC, /* cvg of the delay */
                        double dFwdZC, /* B(0  ,Tpay) / B(0  ,Tstart) */
                        /* vol structure */
                        CMSVol *Volstruct, /* Flat Vol*/
                        /* greek structure */
                        CMSGreeks *Greekstruct,
                        /* numerical param */
                        CMSParams *params,
                        /* Output */
                        double *dCMSTEC) {
  Err err = NULL;

  if (PayRec == SRT_FORWARD) {
    err = swp_f_CMSTEC_rate(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                            dPower, dMargin, dCvgZC, dFwdZC, Volstruct,
                            Greekstruct, params, dCMSTEC);
    return err;
  } else {
    err = swp_f_CMSTEC_Option(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                              dStrike, PayRec, dPower, dMargin, dCvgZC, dFwdZC,
                              Volstruct, Greekstruct, params, dCMSTEC);
    return err;
  }
}

Err swp_f_CMSTEC(                     /* underlying */
                 double dFwdSwapRate, /*Forward Swap Rate */
                 double
                     dnfp, /* number of full period of the default swap rate */
                 double dfreq, /* frequency of the default swap rate 1:A  , 2:S
                                  , 4:Q */
                 double dRateConv, /* adjustment swap default/swap und */
                 /* Option parameters */
                 double dMaturity,       /* CMS Maturity as double */
                 double dStrike,         /* Option Strike */
                 SrtReceiverType PayRec, /* Pay or Rec */
                 double dPower,          /* 1: CMS  , 4: tec  , etc.... */
                 double dMargin,
                 /* delay parameters */
                 double dCvgZC, /* cvg of the delay */
                 double dFwdZC, /* B(0  ,Tpay) / B(0  ,Tstart) */
                 /* vol structure */
                 CMSVol *Volstruct, /* Flat Vol*/
                 /* numerical param */
                 CMSParams *params,
                 /* Output */
                 double *dCMSTEC) {
  Err err = NULL;

  if (PayRec == SRT_FORWARD) {
    err = swp_f_CMSTEC_rate(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                            dPower, dMargin, dCvgZC, dFwdZC, Volstruct, NULL,
                            params, dCMSTEC);
    return err;
  } else {
    err = swp_f_CMSTEC_Option(dFwdSwapRate, dnfp, dfreq, dRateConv, dMaturity,
                              dStrike, PayRec, dPower, dMargin, dCvgZC, dFwdZC,
                              Volstruct, NULL, params, dCMSTEC);
    return err;
  }
}

/* Get the Weights and the points for the intregration */
/* The numeraire is not taken into account */

Err swp_f_CMSTEC_Weights(              /* underlying */
                         double dnfp,  /* number of full period of the default
                                          swap rate */
                         double dfreq, /* frequency of the default swap rate 1:A
                                          , 2:S  , 4:Q */
                         double
                             dRateConv, /* adjustment swap default/swap und */
                         /* Option parameters */
                         double dStrike,         /* Option Strike */
                         SrtReceiverType PayRec, /* Pay or Rec */
                         double dPower, /* 1: CMS  , 4: tec  , etc.... */
                         double dMargin,
                         /* delay parameters */
                         double dCvgZC, /* cvg of the delay */
                         double dspdZC, /* spread to retreive B(0  ,Tpay) / B(0
                                           ,Tstart) */
                         /* numerical param */
                         CMSParams *params,
                         /* Output */
                         double *pdW, double *pdX) {
  Err err = NULL;
  /* Integration param */
  long i;
  double dx;
  double left;
  double right;
  /* Option param */
  double cp = 1.0;
  double invr;
  double invp;
  double invf;
  double strikeadj;
  /* ZC param */
  double cvgZC = dCvgZC;
  double spdZC = dspdZC;
  /* Grid  , weights and Options */
  double *pdDX = NULL;

  double VolUD;
  double IntUD;
  int Mesh;

  /* fudge for the outputs */
  int shift = 1;
  int id = 0;

  /* Initialise param */
  invp = 1.0 / dPower;
  invr = 1.0 / dRateConv;
  invf = 1.0 / dfreq;
  strikeadj = (dStrike - dMargin) * dRateConv;
  if (PayRec == SRT_RECEIVER)
    cp = -1.0;

  if (PayRec == SRT_PAYER) {
    VolUD = params->CAP_VolUD;
    IntUD = params->CAP_IntUD;
    Mesh = params->CAP_Mesh;
  } else {
    shift = 0;
    id = (int)2.0 * params->Nx;
    VolUD = params->FLOOR_VolUD;
    IntUD = params->FLOOR_IntUD;
    Mesh = params->FLOOR_Mesh;
  }

  /* compute the integration space grid */
  pdDX = (double *)malloc((2 * params->Nx + 1) * sizeof(double));
  if (cp == 1) {
    left = strikeadj;
    right = VolUD;
  } else {
    left = IntUD;
    right = VolUD;
  }

  CMS_get_pdX(left, params->Nx, right, cp, Mesh, pdX + shift, pdDX + shift);
  if (cp == 1) {
    left = VolUD;
    right = IntUD;
  } else {
    left = VolUD;
    right = strikeadj;
  }
  CMS_get_pdX(left, params->Nx, right, cp, Mesh, pdX + params->Nx + shift,
              pdDX + params->Nx + shift);

  /* compute the weights */
  swp_f_FinalCms_Weights(pdX + shift, 2 * params->Nx, dnfp, dfreq, dRateConv,
                         dPower, dMargin, dStrike, cp, cvgZC, spdZC,
                         pdW + shift);

  /* do the multiplication */
  if (Mesh == 0) {
    dx = pdDX[shift];
    for (i = 0; i < params->Nx; i++)
      pdW[i + shift] *= dx;

    dx = pdDX[params->Nx + shift];
    for (i = params->Nx; i < 2 * params->Nx; i++)
      pdW[i + shift] *= dx;
  } else {
    for (i = 0; i < 2 * params->Nx; i++)
      pdW[i + shift] *= pdDX[i + shift];
  }

  /* free the memory */
  if (pdDX)
    free(pdDX);
  pdDX = NULL;

  /* add the first point of the control variate factor */
  pdW[id] = CMS_get_controlvariatefactor(strikeadj, dnfp, invr, dMargin, invp,
                                         invf, cvgZC, spdZC);
  pdX[id] = strikeadj;

  return err;
}