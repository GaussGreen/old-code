/*-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

        Description	:	implementation of the Heston model: Calibration

                                (C) 2002 BNP Paribas.. All rights reserved.

        Author		:	 Stefano Galluccio

        Created		:	14.11.2002

        History		:

-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-*/

#include "opfnctns.h"
#include <math.h"
#include <opHeston.h"
#include <utallhdr.h"

/*-----------------------------------------------------------------------------------------------
           This is the generic add-in to calibrate the shifted-Heston model.
Depending on the type of calibration that is specified as input  , it performs
three different algorithms

                   1) If CalibrType = "MATCH_3STRIKES"  , it calibrates on the
BS volatilties associated to an increasing set of (three) strikes that are
inpout through the pointer *CalibrDblArgs. In all cases  , the middle strike HAS
TO BE at-the-money for the routine to work. First  , the local vol is fitted to
                   match the ATM volatility  , then it is freezed and the
correlation is fitted to the volatility skew associated to the two expreme
strikes. Finally  , the first two parameters are frozem and the volvol (alpha)
is fitted to the smile convexity. Then the loop is repeated until convergence.

  `		   2) If CalibrType = "MATCH_CONV"  , it minimizes the distance
between the Heston-generated smile and a set of volatility provided by the user
in input. Here  , a global measure of convexity is used to determine the best
fit. In this case none of the input volatilities is exactly matched.

                   3) If CalibrType = "CHI2_MIN"  , it performs the same
algorithm as MATCH_CONV. The only difference is that the distance between the
two sets of BS volatilities is measured in terms of a standard CHI2 distance. In
this case none of the input volatilities is exactly matched.
-----------------------------------------------------------------------------------------------*/

#define NRANSI
#define ITMAX 100
#define HESTON_EPS 3.0e-8
#define R 0.61803399
#define C (1.0 - R)
#define SHFT2(a, b, c)                                                         \
  (a) = (b);                                                                   \
  (b) = (c);
#define SHFT3(a, b, c, d)                                                      \
  (a) = (b);                                                                   \
  (b) = (c);                                                                   \
  (c) = (d);

Err CalibrateHestonModel(const double dForward, const double dBeta,
                         const double dTau, const double dExpiry,
                         const double dUpperBound, const int iNumSteps,
                         const int iIntegerType, const SRT_Boolean isVolInfFix,
                         const SrtCalibrationType CalibrType,
                         const double *CalibrDblArgs,
                         const double *CalibrInitArgs, double *CalibratedArgs)

{
  Err err = NULL;
  double /* ShiftHeston=HESTON_CONST*(1.-dBeta)/dBeta  , */
      SigmaIftyHeston = CalibrDblArgs[0],
      LambdaHeston = 1. / dTau, SigmaHeston = CalibrInitArgs[0],
      AlphaHeston = CalibrInitArgs[1], RhoHeston = CalibrInitArgs[2],
      BetaHeston, Sigma_new = SigmaHeston, Sigma_old, Alpha_new = AlphaHeston,
      Alpha_old;
  unsigned int n_iter = 0;

  switch (CalibrType) {

  case (MATCH_3STRIKES):

    do {

      Sigma_old = Sigma_new;
      Alpha_old = Alpha_new;

      // STEP 1:  calibration of HestonSigma on the ATMVol
      //          Alphaheston and RhoHeston are frozen

      err = HestonCalibrateSigmaInit(dForward, dExpiry, CalibrDblArgs[6],
                                     AlphaHeston, SigmaIftyHeston, LambdaHeston,
                                     dBeta, RhoHeston, dUpperBound, iNumSteps,
                                     iIntegerType, isVolInfFix, &SigmaHeston);

      // STEP 2:  calibration of RhoHeston on the risk reversal
      //          Sigmaheston and AlphaHeston are frozen

      err = HestonCalibrRho(dForward, dExpiry, CalibrDblArgs, AlphaHeston,
                            SigmaIftyHeston, LambdaHeston, dBeta, SigmaHeston,
                            dUpperBound, iNumSteps, iIntegerType, isVolInfFix,
                            &RhoHeston);

      // STEP 3:  calibration of AlphaHeston on the butterfly
      //          Sigmaheston and RhoHeston are frozen

      err = HestonCalibrAlpha(dForward, dExpiry, CalibrDblArgs, SigmaHeston,
                              SigmaIftyHeston, LambdaHeston, dBeta, RhoHeston,
                              dUpperBound, iNumSteps, iIntegerType, isVolInfFix,
                              &AlphaHeston);

      Sigma_new = SigmaHeston;
      Alpha_new = AlphaHeston;
      n_iter++;

    } while ((fabs(Sigma_new - Sigma_old) > 1.e-7 ||
              fabs(Alpha_new - Alpha_old) > 1.e-3) &&
             n_iter <= ITMAX);

    CalibratedArgs[0] = SigmaHeston;
    CalibratedArgs[1] = AlphaHeston;
    CalibratedArgs[2] = RhoHeston;
    CalibratedArgs[3] = dBeta;

    break;

  default:

    err = CalibrateHestonToSABR(
        dForward, dExpiry, CalibrDblArgs[2], CalibrDblArgs[3], CalibrDblArgs[4],
        CalibrDblArgs[5], &SigmaHeston, &AlphaHeston, &SigmaIftyHeston,
        &LambdaHeston, &BetaHeston, &RhoHeston, dUpperBound, iNumSteps,
        iIntegerType, CalibrDblArgs[1], isVolInfFix, CalibrType);

    CalibratedArgs[0] = SigmaHeston;
    CalibratedArgs[1] = AlphaHeston;
    CalibratedArgs[2] = RhoHeston;
    CalibratedArgs[3] = BetaHeston;

    break;
  }

  return err;
}

/*-----------------------------------------------------------------------------------------------
            Given the set of SABR parameters  , this function returns the
equivalent Heston-shift parameter such that the slope of the two smiles ATM are
matched
-----------------------------------------------------------------------------------------------*/

Err HestonGetShift(double Forward, double Maturity, double SigmaBetaSABR,
                   double AlphaSABR, double BetaSABR, double RhoSABR,
                   double ATMVol, double *Beta) {

  double slope, volshift, shift;
  Err err = NULL;

  if (SigmaBetaSABR == 0.0) {

  } else {

    /* Calculates the ATM Vol */
    err = srt_f_optsarbvol(Forward, Forward, Maturity, SigmaBetaSABR, AlphaSABR,
                           BetaSABR, RhoSABR, SRT_BETAVOL, SRT_LOGNORMAL,
                           &ATMVol);

    if (err)
      return err;
  }

  /* Calculates the shift by matching the smile slope ATM between SABR and
   * Heston*/
  slope = (BetaSABR - 1) * ATMVol / Forward *
          (1 - 0.5 / (1 + (1 - BetaSABR) * (1 - BetaSABR) / 24 * ATMVol *
                              ATMVol * Maturity));
  volshift = 2 * slope * Forward + ATMVol;
  shift = ATMVol * Forward / volshift *
              (1. - ATMVol * ATMVol * Maturity / 24.) /
              (1. - volshift * volshift * Maturity / 24.) -
          Forward;

  // Finally  , transforms the shift into equivalent Beta

  //	(*Beta)=HESTON_CONST/(shift+HESTON_CONST);
  (*Beta) = Forward / (shift + Forward);

  return err;
}

/*-----------------------------------------------------------------------------------------------

        This function assumes that a whole SABR smile has been generated and

-----------------------------------------------------------------------------------------------*/

Err CalibrateHestonToSABR(double Forward, double Maturity, double SigmaBeta,
                          double AlphaSABR, double BetaSABR, double RhoSABR,
                          double *SigmaHeston, double *AlphaHeston,
                          double *SigmaIftyHeston, double *LambdaHeston,
                          double *BetaHeston, double *RhoHeston,
                          double UpperBound, int nSteps, int iIntegerType,
                          double nStdDev, SRT_Boolean isVolInfFix,
                          SrtCalibrationType CalibrType) {
  Err err = NULL;
  double ATMVol;
  double slope, volshift;
  double KMin, KMax, volKMin, volKMax;
  double alphamin, alphamax, alphamid, precalpha;
  double sigma0;
  double *Strikes;
  double *Prices;
  double *Vols;
  double conv, convSABR, dApproxStdev;
  double x, y, z;
  int n_iter_max = 100, n_iter_max2 = 20, i, IntegrType = 0;
  int current_iter = 0, nK = 9, trial_iter = 0;
  double f1, f2, x0, x1, x2, x3, *StrikeVec, xMINold, xMINnew, MINnew,
      tol = 1.e-5, ax, bx, cx;
  double atm_BSSABR_vol, atm_BSHESTON_vol, gap, ShiftHeston;

  Strikes = dvector(1, 2);
  Prices = dvector(1, 2);
  Vols = dvector(1, 2);

  /* Calculates the ATM Vol */
  err =
      srt_f_optsarbvol(Forward, Forward, Maturity, SigmaBeta, AlphaSABR,
                       BetaSABR, RhoSABR, SRT_BETAVOL, SRT_LOGNORMAL, &ATMVol);

  /*
  KMax = Forward*(1+nStdDev*ATMVol*sqrt(Maturity));
  if (nStdDev*Forward - KMax > 0.01)
  {
          KMin = nStdDev*Forward - KMax;
  }
  else
  {
          KMin = 0.01;
  }
  */

  dApproxStdev = SigmaBeta / pow(Forward, 1.0 - BetaSABR) * sqrt(Maturity);

  KMin = Forward * exp(-0.5 * dApproxStdev * dApproxStdev -
                       nStdDev * dApproxStdev); // modifier....
  KMax = Forward *
         exp(-0.5 * dApproxStdev * dApproxStdev + nStdDev * dApproxStdev);

  Strikes[1] = KMin;
  Strikes[2] = KMax;

  /* Calculates the vols far from the money */
  err =
      srt_f_optsarbvol(Forward, KMin, Maturity, SigmaBeta, AlphaSABR, BetaSABR,
                       RhoSABR, SRT_BETAVOL, SRT_LOGNORMAL, &volKMin);
  err =
      srt_f_optsarbvol(Forward, KMax, Maturity, SigmaBeta, AlphaSABR, BetaSABR,
                       RhoSABR, SRT_BETAVOL, SRT_LOGNORMAL, &volKMax);

  x = volKMax + volKMin;
  y = x * 0.5;
  z = y - ATMVol;
  convSABR = z;

  /* Calculates the shift */
  slope = (BetaSABR - 1) * ATMVol / Forward *
          (1 - 0.5 / (1 + (1 - BetaSABR) * (1 - BetaSABR) / 24 * ATMVol *
                              ATMVol * Maturity));
  volshift = 2 * slope * Forward + ATMVol;
  ShiftHeston = ATMVol * Forward / volshift *
                    (1 - ATMVol * ATMVol * Maturity / 24) /
                    (1 - volshift * volshift * Maturity / 24) -
                Forward;

  //. Converts Heston Shift into Heston Beta

  //	(*BetaHeston)=HESTON_CONST/(ShiftHeston+HESTON_CONST);
  //	(*BetaHeston)=Forward/(ShiftHeston+Forward);
  (*BetaHeston) = BetaSABR;

  /* Calculates the rho */
  (*RhoHeston) = RhoSABR;

  /******************************************************************************************************/
  /* Calibrates the alpha */
  /******************************************************************************************************/

  switch (CalibrType) {

  case MATCH_CONV:

    alphamin = 0.001;
    alphamax = min(2., 2 * AlphaSABR);
    precalpha = 0.000001;

    /* Dichotomy on the alpha */
    while ((alphamax - alphamin > precalpha) && (current_iter < n_iter_max)) {
      alphamid = 0.5 * (alphamin + alphamax);

      err = HestonCalibrateSigmaInit(
          Forward, Maturity, ATMVol, alphamid, *SigmaIftyHeston, *LambdaHeston,
          *BetaHeston, *RhoHeston, UpperBound, nSteps, iIntegerType,
          isVolInfFix, &sigma0);

      if (err)
        return err;

      err = HestonPrice(Forward, Strikes, 2, Maturity, sigma0, alphamid,
                        *SigmaIftyHeston, *LambdaHeston, *BetaHeston,
                        *RhoHeston, 1, UpperBound, SRT_CALL, PREMIUM,
                        isVolInfFix, IntegrType, nSteps, Prices);

      if (err)
        return err;

      err = srt_f_optimpvol(Prices[1], Forward, Strikes[1], Maturity, 1,
                            SRT_CALL, SRT_LOGNORMAL, &(Vols[1]));
      err = srt_f_optimpvol(Prices[2], Forward, Strikes[2], Maturity, 1,
                            SRT_CALL, SRT_LOGNORMAL, &(Vols[2]));

      conv = ((Vols[1] + Vols[2]) * 0.5 - ATMVol);

      if (conv < convSABR)
        alphamin = alphamid;
      else
        alphamax = alphamid;
      current_iter += 1;
    }

    (*AlphaHeston) = alphamid;

    err = HestonCalibrateSigmaInit(Forward, Maturity, ATMVol, alphamid,
                                   *SigmaIftyHeston, *LambdaHeston, *BetaHeston,
                                   *RhoHeston, UpperBound, nSteps, iIntegerType,
                                   isVolInfFix, &sigma0);

    (*SigmaHeston) = sigma0;

    break;
  case CHI2_MIN:

    StrikeVec = dvector(1, nK);

    for (i = 1; i <= nK; i++) {

      StrikeVec[i] = KMin + (i - 1.) * (KMax - KMin) / (nK - 1.);
    }

    ax = 0.000001;
    cx = max(2.0, AlphaSABR);
    bx = 0.1; // (ax+cx)/16.; //(ax+cx)/8.;
    precalpha = 0.0000000001;

    err = HestonCalibrateSigmaInit(
        Forward, Maturity, ATMVol, bx, *SigmaIftyHeston, *LambdaHeston,
        *BetaHeston, *RhoHeston, UpperBound, nSteps, iIntegerType, isVolInfFix,
        &sigma0); /* 0.09 this is the first guess */
    xMINnew = bx;

    do {

      xMINold = xMINnew;
      x0 = ax;
      x3 = cx;

      if (fabs(cx - bx) > fabs(bx - ax)) {
        x1 = bx;
        x2 = bx + C * (cx - bx);
      } else {
        x2 = bx;
        x1 = bx - C * (bx - ax);
      }

      err = srt_f_optsarbvol(Forward, Forward, Maturity, SigmaBeta, AlphaSABR,
                             BetaSABR, RhoSABR, SRT_BETAVOL, SRT_LOGNORMAL,
                             &atm_BSSABR_vol);

      /**********************  computes the Chi2 at x1 ***********************/

      err = HestonATMVol(Forward, Maturity, sigma0, x1, *SigmaIftyHeston,
                         *LambdaHeston, *BetaHeston, *RhoHeston, 1.0,
                         UpperBound, nSteps, isVolInfFix, &atm_BSHESTON_vol);

      gap = atm_BSSABR_vol - atm_BSHESTON_vol;

      f1 = SmileChi2(Forward, Maturity, sigma0, x1, *SigmaIftyHeston,
                     *LambdaHeston, *BetaHeston, *RhoHeston, UpperBound, nSteps,
                     iIntegerType, nK, StrikeVec, SigmaBeta, AlphaSABR,
                     BetaSABR, RhoSABR, gap, isVolInfFix);

      /**********************  computes the Chi2 at x1 **********************/

      err = HestonATMVol(Forward, Maturity, sigma0, x2, *SigmaIftyHeston,
                         *LambdaHeston, *BetaHeston, *RhoHeston, 1.0,
                         UpperBound, nSteps, isVolInfFix, &atm_BSHESTON_vol);

      gap = atm_BSSABR_vol - atm_BSHESTON_vol;

      f2 = SmileChi2(Forward, Maturity, sigma0, x2, *SigmaIftyHeston,
                     *LambdaHeston, *BetaHeston, *RhoHeston, UpperBound, nSteps,
                     iIntegerType, nK, StrikeVec, SigmaBeta, AlphaSABR,
                     BetaSABR, RhoSABR, gap, isVolInfFix);

      while (fabs(x3 - x0) > tol * (fabs(x1) + fabs(x2))) {

        if (f2 < f1) {

          SHFT3(x0, x1, x2, R * x1 + C * x3);
          f1 = f2;

          err =
              HestonATMVol(Forward, Maturity, sigma0, x2, *SigmaIftyHeston,
                           *LambdaHeston, *BetaHeston, *RhoHeston, 1.0,
                           UpperBound, nSteps, isVolInfFix, &atm_BSHESTON_vol);

          gap = atm_BSSABR_vol - atm_BSHESTON_vol;

          f2 = SmileChi2(Forward, Maturity, sigma0, x2, *SigmaIftyHeston,
                         *LambdaHeston, *BetaHeston, *RhoHeston, UpperBound,
                         nSteps, iIntegerType, nK, StrikeVec, SigmaBeta,
                         AlphaSABR, BetaSABR, RhoSABR, gap, isVolInfFix);

        } else {

          SHFT3(x3, x2, x1, R * x2 + C * x0);
          f2 = f1;

          err =
              HestonATMVol(Forward, Maturity, sigma0, x1, *SigmaIftyHeston,
                           *LambdaHeston, *BetaHeston, *RhoHeston, 1.0,
                           UpperBound, nSteps, isVolInfFix, &atm_BSHESTON_vol);

          gap = atm_BSSABR_vol - atm_BSHESTON_vol;

          f1 = SmileChi2(Forward, Maturity, sigma0, x1, *SigmaIftyHeston,
                         *LambdaHeston, *BetaHeston, *RhoHeston, UpperBound,
                         nSteps, iIntegerType, nK, StrikeVec, SigmaBeta,
                         AlphaSABR, BetaSABR, RhoSABR, gap, isVolInfFix);
        }
      }

      if (f1 < f2) {
        xMINnew = x1;
        MINnew = f1;
      } else {
        xMINnew = x2;
        MINnew = f2;
      }

      err = HestonCalibrateSigmaInit(
          Forward, Maturity, ATMVol, xMINnew, *SigmaIftyHeston, *LambdaHeston,
          *BetaHeston, *RhoHeston, UpperBound, nSteps, iIntegerType,
          isVolInfFix, &sigma0);

      if (sigma0 > 0.15) /* > 0.2 */ {
        /* sigma0=0.09; */

        bx = bx + pow(-1.0, trial_iter) * 0.01 * (trial_iter + 1);
        err = HestonCalibrateSigmaInit(
            Forward, Maturity, ATMVol, bx, *SigmaIftyHeston, *LambdaHeston,
            *BetaHeston, *RhoHeston, UpperBound, nSteps, iIntegerType,
            isVolInfFix, &sigma0);

        current_iter = -1;
        trial_iter++;
        if (bx < 0.0)
          current_iter = n_iter_max2;
      }

      current_iter += 1;

    } while ((xMINnew - xMINold > precalpha) &&
             (current_iter < n_iter_max2 + 1));

    (*AlphaHeston) = xMINnew;
    err = HestonCalibrateSigmaInit(Forward, Maturity, ATMVol, xMINnew,
                                   *SigmaIftyHeston, *LambdaHeston, *BetaHeston,
                                   *RhoHeston, UpperBound, nSteps, iIntegerType,
                                   isVolInfFix, &sigma0);
    *SigmaHeston = sigma0;

    free_dvector(StrikeVec, 1, 11);
    break;
  }

  free_dvector(Strikes, 1, 2);
  free_dvector(Prices, 1, 2);
  free_dvector(Vols, 1, 2);
  return err;
}

/*-----------------------------------------------------------------------------------------------
  Assuming AlphaHeston (the VolVol) and RhoHeston (Correlation) constant  , this
function determines the Heston local (or initial) volatility to match the ATM BS
volatility
-----------------------------------------------------------------------------------------------*/

Err HestonCalibrateSigmaInit(double Forward, double Maturity, double ATMVol,
                             double AlphaHeston, double SigmaIftyHeston,
                             double LambdaHeston, double dBeta,
                             double RhoHeston, double UpperBound, int nSteps,
                             int iIntegerType, SRT_Boolean isVolInfFix,
                             double *result) {

  /*	double	SigmaMin  ,SigmaMax  ,SigmaMid;

          double	VolMin  ,VolMax  ,VolMid; */
  int n_iter_max = 50;
  int current_iter;
  Err err = NULL;
  double SigmaGuess, VolGuess;
  double Sigma1, Vol1, Sigma2, Vol2, Vol1Shift;
  double errorv;
  double prec = 1.e-7;
  double shift, ShiftHeston;
  double der;

  /* First Guess */

  //	ShiftHeston=HESTON_CONST*(1.-dBeta)/dBeta;
  ShiftHeston = Forward * (1. - dBeta) / dBeta;
  //	SigmaGuess = ATMVol*Forward/(Forward+ShiftHeston);
  if (dBeta < 0.9999) {
    SigmaGuess = ATMVol * pow(Forward, 1.0 - dBeta);
  } else {
    SigmaGuess = ATMVol;
  }
  if (isVolInfFix)
    SigmaIftyHeston = SigmaGuess;

  err = HestonATMVol(Forward, Maturity, SigmaGuess, AlphaHeston,
                     SigmaIftyHeston, LambdaHeston, dBeta, RhoHeston, 1,
                     UpperBound, nSteps, isVolInfFix, &VolGuess);
  if (err) {
    return err;
  }

  errorv = ATMVol - VolGuess;
  if (fabs(errorv) < prec) {
    (*result) = VolGuess;
    return err;
  }

  /* Calculates a second point with the derivative */
  Sigma1 = SigmaGuess;
  Vol1 = VolGuess;
  shift = 0.0005;

  if (errorv < 0)
    shift *= -1;

  if (isVolInfFix)
    SigmaIftyHeston = Sigma1 + shift;

  err = HestonATMVol(Forward, Maturity, Sigma1 + shift, AlphaHeston,
                     SigmaIftyHeston, LambdaHeston, dBeta, RhoHeston, 1,
                     UpperBound, nSteps, isVolInfFix, &Vol1Shift);

  der = (Vol1Shift - Vol1) / shift;

  Sigma2 = Sigma1 + errorv / der;

  current_iter = 1;

  /* Newton */
  while (current_iter <= n_iter_max) {

    if (isVolInfFix)
      SigmaIftyHeston = Sigma2;

    err = HestonATMVol(Forward, Maturity, Sigma2, AlphaHeston, SigmaIftyHeston,
                       LambdaHeston, dBeta, RhoHeston, 1, UpperBound, nSteps,
                       isVolInfFix, &Vol2);

    if (err)
      return err;

    errorv = ATMVol - Vol2;
    if (fabs(errorv) < prec) {
      (*result) = Sigma2;
      return err;
    }
    der = (Vol2 - Vol1) / (Sigma2 - Sigma1);

    Sigma1 = Sigma2;
    Vol1 = Vol2;
    Sigma2 = Sigma1 + errorv / der;
    current_iter++;
  }

  (*result) = Sigma2;
  return err;
}

/*-----------------------------------------------------------------------------------------------
  Assuming SigmaHeston (the Local Vol) and AlphaHeston (the VolVol) constant  ,
this function determines the Heston correlation (RhoHeston) to match the skew at
two different points.
-----------------------------------------------------------------------------------------------*/

Err HestonCalibrRho(const double dForward, const double dExpiry,
                    const double *CalibrDblArgs, const double AlphaHeston,
                    const double SigmaIftyHeston, const double LambdaHeston,
                    const double dBeta, const double SigmaHeston,
                    const double dUpperBound, const int iNumSteps,
                    const int iIntegerType, const SRT_Boolean isVolInfFix,
                    double *RhoHeston) {

  Err err = NULL;

  int i, iter;
  double a = -0.99, b = 0.99, c = 0.99, d, e, min1, min2;
  double fa, fb, fc, p, q, r, s, tol1, xm, tol = 1.e-9;
  double *Strikes = NULL, *NewVols = NULL;
  static int failindex = 0;

  Strikes = dvector(1, 3);
  NewVols = dvector(1, 3);

  for (i = 1; i <= 3; i++)
    Strikes[i] = CalibrDblArgs[i + 1];

  // finds the value of the objective function at the two interval extremes
  // first

  err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                  AlphaHeston, dBeta, a, LambdaHeston, HESTON_TO_LOG,
                  dUpperBound, iNumSteps, iIntegerType, isVolInfFix, NewVols);

  // evaluates the difference between the risk aversals
  fa = (CalibrDblArgs[7] - CalibrDblArgs[5]) - (NewVols[3] - NewVols[1]);

  err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                  AlphaHeston, dBeta, b, LambdaHeston, HESTON_TO_LOG,
                  dUpperBound, iNumSteps, iIntegerType, isVolInfFix, NewVols);

  // evaluates the difference between the risk aversals
  fb = (CalibrDblArgs[7] - CalibrDblArgs[5]) - (NewVols[3] - NewVols[1]);

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    ++failindex;

    if (failindex > 2) {
      smessage("Rho Calibration failed : Change input Beta");
    }

    free_dvector(Strikes, 1, 3);
    free_dvector(NewVols, 1, 3);
    return err;
  }

  fc = fb;
  for (iter = 1; iter <= ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * HESTON_EPS * fabs(b) + 0.5 * tol;
    xm = 0.5 * (c - b);
    if (fabs(xm) <= tol1 || fb == 0.0) {

      (*RhoHeston) = c; // it returned b instead !
      free_dvector(Strikes, 1, 3);
      free_dvector(NewVols, 1, 3);
      return err;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = fabs(p);
      min1 = 3.0 * xm * q - fabs(tol1 * q);
      min2 = fabs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1, xm);

    err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                    AlphaHeston, dBeta, b, LambdaHeston, HESTON_TO_LOG,
                    dUpperBound, iNumSteps, iIntegerType, isVolInfFix, NewVols);

    // evaluates the difference between the risk aversals
    fb = (CalibrDblArgs[7] - CalibrDblArgs[5]) - (NewVols[3] - NewVols[1]);
  }

  //	it reached the maximum number of iterations
  free_dvector(Strikes, 1, 3);
  free_dvector(NewVols, 1, 3);
  return err;
}

/*-----------------------------------------------------------------------------------------------
  Assuming SigmaHeston (the VolVol) and RhoHeston (Correlation) constant  , this
function determines the Heston local (or initial) volatility to match the ATM BS
volatility
-----------------------------------------------------------------------------------------------*/

Err HestonCalibrAlpha(const double dForward, const double dExpiry,
                      const double *CalibrDblArgs, const double SigmaHeston,
                      const double SigmaIftyHeston, const double LambdaHeston,
                      const double dBeta, const double RhoHeston,
                      const double dUpperBound, const int iNumSteps,
                      const int iIntegerType, const SRT_Boolean isVolInfFix,
                      double *AlphaHeston) {

  Err err = NULL;

  int i, iter;
  double a = 1.e-4, b = 2., c = 2., d, e, min1, min2;
  double fa, fb, fc, p, q, r, s, tol1, xm, tol = 1.e-8;
  double *Strikes = NULL, *NewVols = NULL;
  static int failindex = 0;

  Strikes = dvector(1, 3);
  NewVols = dvector(1, 3);

  for (i = 1; i <= 3; i++)
    Strikes[i] = CalibrDblArgs[i + 1];

  // finds the value of the objective function at the two interval extremes
  // first

  err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                  a, dBeta, RhoHeston, LambdaHeston, HESTON_TO_LOG, dUpperBound,
                  iNumSteps, iIntegerType, isVolInfFix, NewVols);

  // evaluates the difference between the butterflies
  fa = (CalibrDblArgs[7] - 2. * CalibrDblArgs[6] + CalibrDblArgs[5]) -
       (NewVols[3] - 2. * NewVols[2] + NewVols[1]);

  err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                  b, dBeta, RhoHeston, LambdaHeston, HESTON_TO_LOG, dUpperBound,
                  iNumSteps, iIntegerType, isVolInfFix, NewVols);

  // evaluates the difference between the butterflies
  fb = (CalibrDblArgs[7] - 2. * CalibrDblArgs[6] + CalibrDblArgs[5]) -
       (NewVols[3] - 2. * NewVols[2] + NewVols[1]);

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
    ++failindex;

    if (failindex > 2) {
      smessage("Alpha Calibration failed : Change input Beta");
    }

    free_dvector(Strikes, 1, 3);
    free_dvector(NewVols, 1, 3);
    return err;
  }

  fc = fb;
  for (iter = 1; iter <= ITMAX; iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c = a;
      fc = fa;
      e = d = b - a;
    }
    if (fabs(fc) < fabs(fb)) {
      a = b;
      b = c;
      c = a;
      fa = fb;
      fb = fc;
      fc = fa;
    }
    tol1 = 2.0 * HESTON_EPS * fabs(b) + 0.5 * tol;
    xm = 0.5 * (c - b);
    if (fabs(xm) <= tol1 || fb == 0.0) {

      (*AlphaHeston) = c; // it returned b instead !
      free_dvector(Strikes, 1, 3);
      free_dvector(NewVols, 1, 3);
      return err;
    }

    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s = fb / fa;
      if (a == c) {
        p = 2.0 * xm * s;
        q = 1.0 - s;
      } else {
        q = fa / fc;
        r = fb / fc;
        p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
        q = (q - 1.0) * (r - 1.0) * (s - 1.0);
      }
      if (p > 0.0)
        q = -q;
      p = fabs(p);
      min1 = 3.0 * xm * q - fabs(tol1 * q);
      min2 = fabs(e * q);
      if (2.0 * p < (min1 < min2 ? min1 : min2)) {
        e = d;
        d = p / q;
      } else {
        d = xm;
        e = d;
      }
    } else {
      d = xm;
      e = d;
    }
    a = b;
    fa = fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1, xm);

    err = HestonVol(dForward, Strikes, 3, dExpiry, SigmaHeston, SigmaIftyHeston,
                    b, dBeta, RhoHeston, LambdaHeston, HESTON_TO_LOG,
                    dUpperBound, iNumSteps, iIntegerType, isVolInfFix, NewVols);

    // evaluates the difference between the butterflies
    fb = (CalibrDblArgs[7] - 2. * CalibrDblArgs[6] + CalibrDblArgs[5]) -
         (NewVols[3] - 2. * NewVols[2] + NewVols[1]);
  }

  //	it reached the maximum number of iterations
  free_dvector(Strikes, 1, 3);
  free_dvector(NewVols, 1, 3);

  return err;
}

/*-----------------------------------------------------------------------------------------------
                At any call  , this function returns the chi2 distance between
                SABR and Heston smiles on nK points.
                This  , in turn  ,  is used in the calibration routine under the
specifics "CHI2_MIN" where the algorithm is in fact based on the minimization of
the Chi2 distance
-----------------------------------------------------------------------------------------------*/

double SmileChi2(double Forward, double Maturity, double Sigma0,
                 double AlphaHeston, double SigmaIftyHeston,
                 double LambdaHeston, double BetaHeston, double RhoHeston,
                 double UpperBound, int nSteps, int IntegerType, int nK,
                 double *StrikeVec, double SigmaBeta, double AlphaSABR,
                 double BetaSABR, double RhoSABR, double gap,
                 SRT_Boolean isVolInfFix) {
  double BS_SABRVol, *BS_HESTONVol, Chi2 = 0.0;
  Err err = NULL;
  int i;

  BS_HESTONVol = dvector(1, nK);

  err =
      HestonVol(Forward, StrikeVec, nK, Maturity, Sigma0, SigmaIftyHeston,
                AlphaHeston, BetaHeston, RhoHeston, LambdaHeston, HESTON_TO_LOG,
                UpperBound, nSteps, IntegerType, isVolInfFix, BS_HESTONVol);

  for (i = 1; i <= nK; i++) {

    err = srt_f_optsarbvol(Forward, StrikeVec[i], Maturity, SigmaBeta,
                           AlphaSABR, BetaSABR, RhoSABR, SRT_BETAVOL,
                           SRT_LOGNORMAL, &BS_SABRVol);

    Chi2 += (BS_HESTONVol[i] + gap - BS_SABRVol) *
            (BS_HESTONVol[i] + gap - BS_SABRVol);
  }

  free_dvector(BS_HESTONVol, 1, nK);
  return Chi2;
}

#undef C
#undef R
#undef SHFT2
#undef SHFT3
#undef ITMAX
#undef NRANSI
#undef HESTON_EPS
