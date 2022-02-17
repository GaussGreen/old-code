/* ===================================================================================
   FILENAME:      swp_f_fwd_swaption.c

   PURPOSE:       Computes price of a fwd swaption  , with maturity T  , to
   enter into a Fwd swap between T-1 > T and T_2 > T. To model the fwd vol  ,
   the function uses a method based on CMS replication and Copula simulation in
   a Swap market model context
   ===================================================================================
 */

#pragma warning(disable : 4786) // Disable long name warnings

#include "math.h"
#include "num_h_allhdr.h"
#include "opHeston.h"
#include "opfnctns.h"
#include "swp_h_all.h"
#include "swp_h_fwd_swaption.h"
#include "swp_h_swap_pricing.h"
#include "swp_h_swap_simple.h"
#include "swp_h_swaption.h"
#include "swp_h_vol.h"

// IRD
// #include "IRDCommonHeaders.h"
// #include "DteDate.h"

Err PriceFwdSwaption(char *cYCname, char *cRefRname, long xlMaturityDate,
                     long xlStartDate, long xlEndDate, SrtCompounding srtFreq,
                     SrtBasisCode srtBasis, int nStrikes, double *Strikes,
                     SrtCallPutType srtCallPut, int xlStudDegree, int xlNumSim,
                     char *xlTypeVol, int xlNPts, double xlUpBound,
                     int xlnClasses, double *SigmaHeston, double *AlphaHeston,
                     double *LambdaHeston, double *BetaHeston,
                     double *RhoHeston,

                     double *SigmaHestonShort, double *AlphaHestonShort,
                     double *LambdaHestonShort, double *BetaHestonShort,
                     double *RhoHestonShort, double **CorrelLong,
                     double *dPrice, double *dPriceSwnShort,
                     double *dPriceSwnLong)

{

  ///// Objects and variables declaration

  Err err = NULL;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  SwapDP SwapShort, SwapLong;
  Date lToday, lSpotLag;
  Date *lPayDatesShort = NULL, *lStartDatesShort = NULL, *lEndDatesShort = NULL,
       *lPayDatesLong = NULL, *lStartDatesLong = NULL, *lEndDatesLong = NULL;
  double *dCoveragesShort = NULL, *dCoveragesLong = NULL, **dCopulaCube = NULL,
         *dMean_v = NULL, *dAvg = NULL;
  double **dPPDx = NULL, **dPPDy = NULL, **dCPDy = NULL, **dPPD_Shortx = NULL,
         **dPPD_Shorty = NULL, **dCPD_Shorty = NULL;
  int i, j, k, iNDatesShort, iNDatesLong, iNPayDatesShort, iNPayDatesLong,
      nClass = 1 << xlnClasses, flag;
  double *dPayoff = NULL, *dPayoffSwnShort = NULL, *dPayoffSwnLong = NULL,
         dFwdLong, dFwdShort, dLvlLong, dLvlShort, *dCoveragesAux = NULL;
  Date *lPayDatesAux = NULL, *lStartDatesAux = NULL, *lEndDatesAux = NULL;
  int iNDatesAux, iNPayDatesAux;

  char *cBasis, *cFreq;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);
  /////////////////////////////////////////////////////////////////////////
  ///// defines the pointer to the curve and retrieves today and spot dates
  /////////////////////////////////////////////////////////////////////////

  pCurve = lookup_curve(cYCname);
  ///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not
  ///found");
  lToday = (Date)get_today_from_curve(pCurve);
  lSpotLag = (Date)get_spotlag_from_curve(pCurve);

  ///////////////////////////////////////////////////////////////
  ///// defines the relevant swap objects
  ///////////////////////////////////////////////////////////////

  err = swp_f_setSwapDP(xlMaturityDate, xlStartDate, srtFreq, srtBasis,
                        &SwapShort);
  err = swp_f_setSwapDP(xlStartDate, xlEndDate, srtFreq, srtBasis, &SwapLong);

  err = swp_f_make_FixedLegDatesAndCoverages(
      &SwapShort, lToday, &lPayDatesShort, &iNPayDatesShort, &lStartDatesShort,
      &lEndDatesShort, &dCoveragesShort, &iNDatesShort);

  err = swp_f_make_FixedLegDatesAndCoverages(
      &SwapLong, lToday, &lPayDatesAux, &iNPayDatesAux, &lStartDatesAux,
      &lEndDatesAux, &dCoveragesAux, &iNDatesAux);

  iNPayDatesLong = iNPayDatesShort + iNPayDatesAux - 1;
  iNDatesLong = iNDatesShort + iNDatesAux;

  dCoveragesLong = dvector(0, iNDatesLong - 1);
  lStartDatesLong = lngvector(0, iNDatesLong - 1);
  lEndDatesLong = lngvector(0, iNDatesLong - 1);
  lPayDatesLong = lngvector(0, iNPayDatesLong - 1);

  lPayDatesLong[0] = lPayDatesShort[0];
  for (i = 0; i < iNDatesShort; i++) {

    lPayDatesLong[i + 1] = lPayDatesShort[i + 1];
    lStartDatesLong[i] = lStartDatesShort[i];
    lEndDatesLong[i] = lEndDatesShort[i];
    dCoveragesLong[i] = dCoveragesShort[i];
  }

  for (i = 0; i < iNDatesAux; i++) {

    lPayDatesLong[i + iNDatesShort + 1] = lPayDatesAux[i + 1];
    lStartDatesLong[i + iNDatesShort] = lStartDatesAux[i];
    lEndDatesLong[i + iNDatesShort] = lEndDatesAux[i];
    dCoveragesLong[i + iNDatesShort] = dCoveragesAux[i];
  }

  ///////////////////////////////////////////////////////////////
  ////////  Retrieves the densities
  ///////////////////////////////////////////////////////////////

  dPPDx = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);
  dPPDy = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);
  dCPDy = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);

  dPPD_Shortx = dmatrix(0, iNPayDatesShort - 1, 0, 2 * nClass);
  dPPD_Shorty = dmatrix(0, iNPayDatesShort - 1, 0, 2 * nClass);
  dCPD_Shorty = dmatrix(0, iNPayDatesShort - 1, 0, 2 * nClass);

  err = Fwd_Swaption_Get_Cum_Density(
      lToday, cYCname, cRefRname, srtFreq, srtBasis,

      iNPayDatesLong, lPayDatesLong, dCoveragesLong,

      SigmaHeston, AlphaHeston, LambdaHeston, BetaHeston, RhoHeston,

      xlNPts, xlUpBound, nClass,

      dPPDx, dPPDy, dCPDy);

  ///////////////////////////////////////////////////////////////
  ////////  Retrieves the Matrix of Joint distributions
  ///////////////////////////////////////////////////////////////

  dMean_v = dvector(0, iNDatesLong - 1);
  dAvg = dvector(0, iNDatesLong - 1);
  for (i = 0; i < iNDatesLong; i++) {
    dAvg[i] = 0.0;
    dMean_v[i] = 0.0;
  }

  dCopulaCube = dmatrix(0, xlNumSim, 0, iNDatesLong - 1);
  GetStudentCplDev(xlStudDegree, dMean_v, CorrelLong, xlNumSim,
                   xlStudDegree + iNDatesLong, dPPDx, dCPDy, 2 * nClass + 1, 0,
                   SOBOL, dCopulaCube);

  ///////////////////////////////////////////////////////////////
  ////////  Simulation
  ///////////////////////////////////////////////////////////////

  //   swp_f_ForwardRate(lPayDatesLong[0]  ,lPayDatesLong[1]  ,cFreq  ,cBasis
  //   ,cYCname  ,cRefRname  ,&dFwd);

  dPayoff = dvector(0, nStrikes - 1);
  dPayoffSwnShort = dvector(0, nStrikes - 1);
  dPayoffSwnLong = dvector(0, nStrikes - 1);

  for (i = 0; i < nStrikes; i++) {

    dPrice[i] = 0.0;
    dPriceSwnShort[i] = 0.0;
    dPriceSwnLong[i] = 0.0;
  }

  flag = -1;
  if (srtCallPut == SRT_CALL)
    flag = 1;

  j = 3;
  for (i = 3; i < xlNumSim; i++) {

    err = GetFwdSwaptionPayoff(
        dCopulaCube, iNDatesLong, lEndDatesLong, dCoveragesLong,

        iNDatesShort, lEndDatesShort, dCoveragesShort,

        nStrikes, Strikes, i, dPayoff, dPayoffSwnShort, dPayoffSwnLong, 1);

    for (k = 0; k < nStrikes; k++) {

      dPrice[k] += (dPayoff[k] - dPrice[k]) / (j - 2);
      dPriceSwnShort[k] += (dPayoffSwnShort[k] - dPriceSwnShort[k]) / (j - 2);
      dPriceSwnLong[k] += (dPayoffSwnLong[k] - dPriceSwnLong[k]) / (j - 2);
    }

    j++;
  }

  // prices a put by call/put parity

  if (flag != 1) {

    err = swp_f_ForwardRate(lPayDatesLong[0], lPayDatesLong[iNPayDatesLong - 1],
                            cFreq, cBasis, cYCname, cRefRname, &dFwdLong);
    err =
        swp_f_ForwardRate(lPayDatesShort[0], lPayDatesLong[iNPayDatesShort - 1],
                          cFreq, cBasis, cYCname, cRefRname, &dFwdShort);

    err =
        swp_f_LevelPayment(lPayDatesLong[0], lPayDatesLong[iNPayDatesLong - 1],
                           cFreq, cBasis, cYCname, cRefRname, &dLvlLong);
    err = swp_f_LevelPayment(lPayDatesShort[0],
                             lPayDatesLong[iNPayDatesShort - 1], cFreq, cBasis,
                             cYCname, cRefRname, &dLvlShort);

    for (k = 0; k < nStrikes; k++) {
      dPrice[k] -= dLvlLong * (dFwdLong - Strikes[k]) -
                   dLvlShort * (dFwdShort - Strikes[k]);
    }
  }

  ///////////////////////////////////////////////////////////////
  ////////  Frees the memory and return
  ///////////////////////////////////////////////////////////////

  free_dmatrix(dPPDx, 0, iNPayDatesLong - 1, 0, 2 * nClass);
  free_dmatrix(dPPDy, 0, iNPayDatesLong - 1, 0, 2 * nClass);
  free_dmatrix(dCPDy, 0, iNPayDatesLong - 1, 0, 2 * nClass);

  free_dmatrix(dPPD_Shortx, 0, iNPayDatesShort - 1, 0, 2 * nClass);
  free_dmatrix(dPPD_Shorty, 0, iNPayDatesShort - 1, 0, 2 * nClass);
  free_dmatrix(dCPD_Shorty, 0, iNPayDatesShort - 1, 0, 2 * nClass);

  free_dvector(dMean_v, 0, iNDatesLong - 1);
  free_dvector(dAvg, 0, iNDatesLong - 1);
  free_dmatrix(dCopulaCube, 0, xlNumSim, 0, iNDatesLong - 1);

  free_dvector(dPayoff, 0, nStrikes - 1);
  free_dvector(dPayoffSwnShort, 0, nStrikes - 1);
  free_dvector(dPayoffSwnLong, 0, nStrikes - 1);

  free_dvector(dCoveragesLong, 0, iNDatesLong - 1);
  free_lngvector(lStartDatesLong, 0, iNDatesLong - 1);
  free_lngvector(lEndDatesLong, 0, iNDatesLong - 1);
  free_lngvector(lPayDatesLong, 0, iNPayDatesLong - 1);

  return err;
}

/*-----------------------------------------------------------------------------------*/
/*        Function that transforms smile from the swap measure into Fwd measure
 */
/*-----------------------------------------------------------------------------------*/

Err Fwd_Swaption_Get_Cum_Density(Date lToday, char *cYCname, char *cRefRname,
                                 SrtCompounding srtFreq, SrtBasisCode srtBasis,

                                 int iNPayDatesLong, Date *lPayDatesLong,
                                 double *dCoveragesLong,

                                 double *SigmaHeston, double *AlphaHeston,
                                 double *LambdaHeston, double *BetaHeston,
                                 double *RhoHeston,

                                 int iNPts, double dUpBound, int nClass,

                                 double **dPPDx, double **dPPDy, double **dCPDy)

{
  Err err = NULL;
  int i;
  double dMat, dFwd, dLevel, dCov;
  char *cBasis, *cFreq;
  //	double shiftHeston;

  dMat = (lPayDatesLong[0] - lToday) / 365.;
  //		dDisc=swp_f_df(lToday  , lPayDatesLong[0]  , cYCname);

  //////////////////////////////////////////
  /// translate SORT types into strings
  //////////////////////////////////////////

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);

  if (cFreq == "ANNUAL") {
    dCov = 1.0;
  } else if (cFreq == "SEMIANNUAL") {
    dCov = 0.5;
  } else if (cFreq == "QUARTERLY") {
    dCov = 0.25;
  } else if (cFreq == "MONTHLY") {
    dCov = 1. / 12.;
  }

  //////////////////////////////////////////////////////////
  ////////////////// generates the cumulative densities
  //////////////////////////////////////////////////////////

  for (i = 1; i < iNPayDatesLong; i++) {

    //      computes the Fwd swap rate

    //		shiftHeston = BetaHeston[i-1];

    err = swp_f_ForwardRate(lPayDatesLong[0], lPayDatesLong[i], cFreq, cBasis,
                            cYCname, cRefRname, &dFwd);
    err = swp_f_LevelPayment(lPayDatesLong[0], lPayDatesLong[i], cFreq, cBasis,
                             cYCname, cRefRname, &dLevel);

    err = GetDensity_Lvl_Measure(
        nClass, dMat, dFwd, SigmaHeston[i - 1], AlphaHeston[i - 1], 1.0,
        LambdaHeston[i - 1], BetaHeston[i - 1], RhoHeston[i - 1], dLevel,
        dUpBound, iNPts, i - 1, dCov, dPPDx, dPPDy, dCPDy);
  }

  /*
          for (i=1;i<iNPayDatesShort;i++) {

  //      computes the Fwd swap rate

                  swp_f_ForwardRate(lPayDatesShort[0]  ,lPayDatesShort[i] ,cFreq
  ,cBasis  ,cYCname  ,cRefRname  ,&dFwd); err = GetDensity_Lvl_Measure(nClass
  ,dMat  ,dFwd  ,SigmaHestonShort[i-1]  ,AlphaHestonShort[i-1]  ,1.0
  ,LambdaHestonShort[i-1]  , ShiftHestonShort[i-1]  ,RhoHestonShort[i-1]  ,dDisc
  ,dUpBound  ,iNPts  ,i-1  ,dCov  ,dPPD_Shortx  ,dPPD_Shorty  ,dCPD_Shorty);
          }

  */

  return err;
}

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns the payoff of a fwd swpation conditional to the
 */
/*          realization of the spanning set of swap rates */
/*-----------------------------------------------------------------------------------
 */

Err GetFwdSwaptionPayoff(double **dCopulaCube, int iNDatesLong,
                         Date *lEndDatesLong, double *dCoveragesLong,

                         int iNDatesShort, Date *lEndDatesShort,
                         double *dCoveragesShort,

                         int nStrikes, double *xlStrike, int iMC,
                         double *dPayoff, double *dPayoffSwnShort,
                         double *dPayoffSwnLong, int flag) {

  Err err = NULL;
  double *DfSet = NULL, dLevelShort, dLevelLong;
  int i;

  DfSet = dvector(0, iNDatesLong - 1);

  err = GetSpanningSetDf(dCopulaCube, iNDatesLong, lEndDatesLong,
                         dCoveragesLong, iMC, DfSet);

  /*	if (!((*isAF) == SRT_TRUE)) {
                  free_dvector(DfSet  ,0  ,iNDatesLong-1);
                  return err;
          }
  */

  err = GetSimulatedLevel(DfSet, iNDatesLong, lEndDatesLong, dCoveragesLong,
                          &dLevelLong);
  err = GetSimulatedLevel(DfSet, iNDatesShort, lEndDatesShort, dCoveragesShort,
                          &dLevelShort);

  for (i = 0; i < nStrikes; i++) {

    //	(*dPayoff)= DfSet[iNDatesLong-1];

    dPayoff[i] = DMAX(
        0.0,
        flag *
            (dLevelLong * (dCopulaCube[iMC][iNDatesLong - 1] - xlStrike[i]) -
             dLevelShort * (dCopulaCube[iMC][iNDatesShort - 1] - xlStrike[i])));

    //	(*dPayoff)=dLevelLong;

    //	(*dPayoff)=dCopulaCube[iMC][iNDatesLong-1];

    //	(*dPayoff)=DMAX(0.0  , dCopulaCube[iMC][iNDatesLong-1]-xlStrike);

    dPayoffSwnShort[i] =
        DMAX(0.0, flag * (dLevelShort *
                          (dCopulaCube[iMC][iNDatesShort - 1] - xlStrike[i])));

    dPayoffSwnLong[i] =
        DMAX(0.0, flag * (dLevelLong *
                          (dCopulaCube[iMC][iNDatesLong - 1] - xlStrike[i])));

    //			dPayoffSwnShort[i]=flag*(dLevelShort*(dCopulaCube[iMC][iNDatesShort-1])
    //) ;

    //			dPayoffSwnLong[i]=flag*(dLevelLong*(dCopulaCube[iMC][iNDatesLong-1])
    //) ;
  }

  free_dvector(DfSet, 0, iNDatesLong - 1);
  return err;
}

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns the level computed at each iteration */
/*                      of the spanning set of swap rates */
/*-----------------------------------------------------------------------------------
 */

Err GetSimulatedLevel(double *DfSet, int iNDates, Date *lEndDates,
                      double *dCoverages, double *dLevel) {
  Err err = NULL;
  int i;

  (*dLevel) = 0.0;
  for (i = 0; i < iNDates; i++) {
    (*dLevel) += dCoverages[i] * DfSet[i];
  }

  return err;
}

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns a set of n probability densities for the Heston
 * model */
/*-----------------------------------------------------------------------------------
 */

//-------------------------------------------------------------------
//-----From Qlevel to QT with level cash approx with regular mesh----
//-------------------------------------------------------------------
Err GetDensity_Lvl_Measure(int nClass, double dMat, double dFwd, double dSigmaH,
                           double dAlphaH, double dSigmaIftyH, double dMeanRevH,
                           double dBetaH, double dRhoH, double dDisc,
                           double dUpBound, int nSteps, int nTen, double dCov,
                           double **dPPDx, double **dPPDy, double **dCPDy) {

  double dEqBSVol, dStdDev, dK, LvlC_Fwd, LvlC_x;
  double *dVec = NULL, *dX = NULL, *dY = NULL, dThreshold = 1.e-3;
  int i;
  double dShiftH;
  Err err = NULL;

  //	dShiftH = HESTON_CONST*(1.-dBetaH)/dBetaH;
  dShiftH = dFwd * (1. - dBetaH) / dBetaH;

  dEqBSVol = (dFwd + dShiftH) / dFwd * dSigmaH;
  dStdDev = 6.;
  /////	dCoeffAdjustStrike = exp( dEqBSVol*sqrt(dMat) / nClass);

  ///// gets the two boundaries

  dPPDx[nTen][0] = -dShiftH + dFwd * exp(-0.5 * dEqBSVol * dEqBSVol * dMat -
                                         dStdDev * dEqBSVol * sqrt(dMat));
  //	dPPDx[nTen][2*nClass]=-dShiftH + dFwd * exp( -
  //0.5*dEqBSVol*dEqBSVol*dMat + dStdDev*dEqBSVol*sqrt(dMat));

  dX = dvector(1, 1);
  dY = dvector(1, 1);

  dX[1] = 0.15;
  for (i = 1; i <= 15; i++) {

    dX[1] += 0.04;
    HestonPrice(dFwd, dX, 1, dMat, dSigmaH, dAlphaH, dSigmaIftyH, dMeanRevH,
                dBetaH, dRhoH, dDisc, dUpBound, SRT_CALL, DENSITY, SRT_TRUE, 2,
                nSteps, dY);
    if (dY[1] < dThreshold) {

      dPPDx[nTen][2 * nClass] = dX[1];
      break;
    }
  }

  free_dvector(dX, 1, 1);
  free_dvector(dY, 1, 1);

  dK = (dPPDx[nTen][2 * nClass] - dPPDx[nTen][0]) / (2. * nClass);
  ///// runs through the points

  for (i = 1; i < 2 * nClass; i++) {
    dPPDx[nTen][i] = dPPDx[nTen][0] + i * dK;
    dCPDy[nTen][i] = 0.0;
    dPPDy[nTen][i] = 0.0;
  }
  /// computes the density of the heston model on the whole domain;

  HestonPrice(dFwd, dPPDx[nTen] - 1, 2 * nClass + 1, dMat, dSigmaH, dAlphaH,
              dSigmaIftyH, dMeanRevH, dBetaH, dRhoH, dDisc, dUpBound, SRT_CALL,
              DENSITY, SRT_TRUE, 2, nSteps, dPPDy[nTen] - 1);

  /// computes the Level cash asset LvlC(Fwd)

  LvlC_Fwd = (1. - 1. / pow(1 + dFwd * dCov, nTen + 1)) / dFwd;

  /// computes the new density

  LvlC_x =
      (1. - 1. / pow(1 + dPPDx[nTen][0] * dCov, nTen + 1)) / dPPDx[nTen][0];
  dPPDy[nTen][0] *= LvlC_Fwd / LvlC_x;
  dCPDy[nTen][0] = 0.0;

  for (i = 1; i <= 2 * nClass; i++) {

    LvlC_x =
        (1. - 1. / pow(1 + dPPDx[nTen][i] * dCov, nTen + 1)) / dPPDx[nTen][i];
    dPPDy[nTen][i] *= LvlC_Fwd / LvlC_x;
    dCPDy[nTen][i] =
        dCPDy[nTen][i - 1] + (dPPDy[nTen][i] + dPPDy[nTen][i - 1]) / 2. * dK;
  }

  for (i = 1; i <= 2 * nClass; i++)
    dCPDy[nTen][i] /= dCPDy[nTen][2 * nClass];

  return err;
}

/*---------------------------------------------------------------------------------------
 */
/*        Function that returns a set n Discount factors compatible with the set
 * of swaps */
/*---------------------------------------------------------------------------------------
 */

Err GetSpanningSetDf(double **dCopulaCube, int iNDates, Date *lEndDates,
                     double *dCoverages, int iMC, double *DfSet) {

  Err err = NULL;
  double dProd, dSum;
  int i, k, l;

  for (i = 0; i < iNDates; i++) {

    dSum = 0.0;
    for (k = 0; k < i; k++) {

      dProd = 1.0;
      for (l = k; l <= i; l++)
        dProd *= (1. + dCoverages[l] * dCopulaCube[iMC][l]);

      dSum += dCoverages[k] / dProd;
    }
    DfSet[i] = -dSum * dCopulaCube[iMC][i] +
               1. / (1 + dCoverages[i] * dCopulaCube[iMC][i]);
  }

  return err;
}

/*---------------------------------------------------------------------------------------
 */
/*        Function that calibrates the smile of a caplet given */
/*---------------------------------------------------------------------------------------
 */

Err CalibrateCorrelation(char *cYCname, char *cRefRname, int xlMaturityDate,
                         int nDates, long *lDates, double *dCoverages,
                         double dCoverageShort, double dFwd,

                         SrtCompounding srtFreq, SrtBasisCode srtBasis,
                         int nStrikes, double *Strikes, double *BSVols,

                         SrtCallPutType srtCallPut, int xlStudDegree,
                         int xlNumSim, char *xlTypeVol, int xlNPts,
                         double xlUpBound, int xlnClasses,

                         double *SigmaHeston, double *AlphaHeston,
                         double *LambdaHeston, double *BetaHeston,
                         double *RhoHeston,

                         double **CorrelLong, double **CorrelOutput)

{

  ///// Objects and variables declaration

  Err err = NULL;
  SrtCurvePtr pCurve = lookup_curve(cYCname);
  Date lToday, lSpotLag;

  double **dCopulaCube = NULL, *dMean_v = NULL, *dAvg = NULL, dMaturity;
  double **dPPDx = NULL, **dPPDy = NULL, **dCPDy = NULL;
  int i, j, nClass = 1 << xlnClasses;
  double *dCoveragesLong = NULL, x1 = 0.001, x2 = 0.1, tol = 1.e-6;

  Date *lPayDatesLong = NULL;
  int iNPayDatesLong = nDates + 1;

  char *cBasis, *cFreq;
  double *dPayoff;

  err = translate_basis(&cBasis, srtBasis);
  err = translate_compounding(&cFreq, srtFreq);
  /////////////////////////////////////////////////////////////////////////
  ///// defines the pointer to the curve and retrieves today and spot dates
  /////////////////////////////////////////////////////////////////////////

  pCurve = lookup_curve(cYCname);
  ///// if (!pCurve) throw XSort("Fatal: (FwdSwaption) Yield Curve Object not
  ///found");
  lToday = (Date)get_today_from_curve(pCurve);
  lSpotLag = (Date)get_spotlag_from_curve(pCurve);

  ///////////////////////////////////////////////////////////////
  ////////  Retrieves the densities
  ///////////////////////////////////////////////////////////////

  lPayDatesLong = lngvector(0, iNPayDatesLong - 1);
  dCoveragesLong = dvector(0, iNPayDatesLong - 1);

  lPayDatesLong[0] = (long)xlMaturityDate;
  for (i = 1; i < iNPayDatesLong; i++) {

    lPayDatesLong[i] = lDates[i - 1];
    dCoveragesLong[i - 1] = dCoverages[i - 1];
  }

  dPPDx = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);
  dPPDy = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);
  dCPDy = dmatrix(0, iNPayDatesLong - 1, 0, 2 * nClass);

  err = Fwd_Swaption_Get_Cum_Density(
      lToday, cYCname, cRefRname, srtFreq, srtBasis,

      iNPayDatesLong, lPayDatesLong, dCoveragesLong,

      SigmaHeston, AlphaHeston, LambdaHeston, BetaHeston, RhoHeston,

      xlNPts, xlUpBound, nClass,

      dPPDx, dPPDy, dCPDy);

  dMaturity = (xlMaturityDate - lToday) / 365.;

  ///////////////////////////////////////////////////////////////
  ////////  Retrieves the Matrix of Joint distributions
  ///////////////////////////////////////////////////////////////

  dMean_v = dvector(0, nDates - 1);
  for (i = 0; i < nDates; i++)
    dMean_v[i] = 0.0;

  dCopulaCube = dmatrix(0, xlNumSim, 0, nDates - 1);

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  ///////////////////////////  TEST FOR PRICING ////////////////////////////////

  dPayoff = dvector(0, nStrikes - 1);

  GetStudentCplDev(xlStudDegree, dMean_v, CorrelLong, xlNumSim,
                   xlStudDegree + nDates, dPPDx, dCPDy, 2 * nClass + 1, 0,
                   SOBOL, dCopulaCube);

  /// run MC simulation and evaluates the caplet

  for (i = 3; i < xlNumSim; i++) {

    err = GetCapletPayoff(dCopulaCube, nDates, lDates, dCoverages,

                          dCoverageShort,

                          nStrikes, Strikes, i, dPayoff);

    for (j = 0; j < nDates; j++) {
      CorrelOutput[j][0] += (dPayoff[j] - CorrelOutput[j][0]) / (i - 2);
    }
  }

  free_dvector(dPayoff, 0, nStrikes - 1);

  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////

  /*
          ///////////////////////////////////////////////////////////////
          ////////  Find the decay rate
          ///////////////////////////////////////////////////////////////

          err =  FindDecayRate (
                                                          x1  ,
                                                          x2  ,
                                                          tol  ,

                                                          dMaturity  ,
                                                          nDates  ,
                                                          lDates  ,

                                                          dCoverages  ,
                                                          dCoverageShort  ,

                                                          nStrikes  ,
                                                          Strikes  ,
                                                          BSVols  ,
                                                          dFwd  ,

                                                          dCopulaCube  ,

                                                          dPPDx  ,
                                                          dCPDy  ,
                                                          xlStudDegree  ,
                                                          dMean_v  ,
                                                          nClass  ,
                                                          xlNumSim  ,

                                                          SigmaHeston  ,
                                                          AlphaHeston  ,
                                                          LambdaHeston  ,
                                                          BetaHeston  ,
                                                          RhoHeston  ,

                                                          CorrelLong  ,
                                                          &dDecay
                                           );
  */

  ///////////////////////////////////////////////////////////////
  ////////  Frees the memory and return
  ///////////////////////////////////////////////////////////////

  free_dmatrix(dPPDx, 0, nDates - 1, 0, 2 * nClass);
  free_dmatrix(dPPDy, 0, nDates - 1, 0, 2 * nClass);
  free_dmatrix(dCPDy, 0, nDates - 1, 0, 2 * nClass);

  free_dvector(dMean_v, 0, nDates - 1);
  free_dmatrix(dCopulaCube, 0, xlNumSim, 0, nDates - 1);

  free_lngvector(lPayDatesLong, 0, iNPayDatesLong - 1);
  free_dvector(dCoveragesLong, 0, iNPayDatesLong - 1);

  return err;
}

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns the payoff of a fwd swpation conditional to the
 */
/*          realization of the spanning set of swap rates */
/*-----------------------------------------------------------------------------------
 */

#define NRANSI
#define ITMAX_DECAY 100
#define EPS_DECAY 3.0e-8

Err FindDecayRate(double x1, double x2, double tol,

                  double Maturity, int iNDates, Date *lEndDates,
                  double *dCoverages, double dCoverageShort,

                  int nStrikes, double *Strikes, double *BSVols, double dFwd,

                  double **dCopulaCube,

                  double **dPPDx, double **dCPDy, int xlStudDegree,
                  double *dMean_v, int nClass, int nSim,

                  double *SigmaHestonShort, double *AlphaHestonShort,
                  double *LambdaHestonShort, double *BetaHestonShort,
                  double *RhoHestonShort,

                  double **dCorrel, double *dDecay) {

  Err err = NULL;
  int iter;
  double a = x1, b = x2, c = x2, d, e, min1, min2;
  double fc, p, q, r, s, tol1, xm, fa, fb;

  err = EvalChi2Error(
      a, Maturity, iNDates, lEndDates, dCoverages, dCoverageShort, nStrikes,
      Strikes, BSVols, dFwd, dCopulaCube, dPPDx, dCPDy, xlStudDegree, dMean_v,
      nClass, nSim, SigmaHestonShort, AlphaHestonShort, LambdaHestonShort,
      BetaHestonShort, RhoHestonShort, dCorrel, &fa);

  err = EvalChi2Error(
      b, Maturity, iNDates, lEndDates, dCoverages, dCoverageShort, nStrikes,
      Strikes, BSVols, dFwd, dCopulaCube, dPPDx, dCPDy, xlStudDegree, dMean_v,
      nClass, nSim, SigmaHestonShort, AlphaHestonShort, LambdaHestonShort,
      BetaHestonShort, RhoHestonShort, dCorrel, &fb);

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
    return err;

  fc = fb;
  for (iter = 1; iter <= ITMAX_DECAY; iter++) {
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
    tol1 = 2.0 * EPS_DECAY * fabs(b) + 0.5 * tol;
    xm = 0.5 * (c - b);

    if (fabs(xm) <= tol1 || fb == 0.0) {

      (*dDecay) = c;
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

    err = EvalChi2Error(
        b, Maturity, iNDates, lEndDates, dCoverages, dCoverageShort, nStrikes,
        Strikes, BSVols, dFwd, dCopulaCube, dPPDx, dCPDy, xlStudDegree, dMean_v,
        nClass, nSim, SigmaHestonShort, AlphaHestonShort, LambdaHestonShort,
        BetaHestonShort, RhoHestonShort, dCorrel, &fb);
  }

  return err;
}

#undef ITMAX_DECAY
#undef EPS_DECAY
#undef NRANSI

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns the payoff of a fwd swpation conditional to the
 */
/*          realization of the spanning set of swap rates */
/*-----------------------------------------------------------------------------------
 */

Err EvalChi2Error(double dTryDecay, double Maturity, int iNDates,
                  Date *lEndDates, double *dCoverages,

                  double dCoverageShort, int nStrikes, double *Strikes,
                  double *BSVols, double dFwd,

                  double **dCopulaCube, double **dPPDx, double **dCPDy,
                  int xlStudDegree, double *dMean_v, int nClass, int nSim,
                  double *SigmaHestonShort, double *AlphaHestonShort,
                  double *LambdaHestonShort, double *BetaHestonShort,
                  double *RhoHestonShort, double **Correl, double *ChiError) {

  Err err = NULL;

  int i, j;
  double *dPayoff = NULL, *dPrice = NULL, ImplVol;

  dPayoff = dvector(0, nStrikes - 1);
  dPrice = dvector(0, nStrikes - 1);
  /// redefines the new correlation matrix

  for (i = 0; i < iNDates; i++) {
    for (j = i + 1; j < iNDates; j++) {
      Correl[i][j] = Correl[j][i] = exp(-dTryDecay * fabs(i - j));
    }
    Correl[i][i] = 1.;
  }

  /// retrieves the new Copula Cube

  GetStudentCplDev(xlStudDegree, dMean_v, Correl, nSim, xlStudDegree + iNDates,
                   dPPDx, dCPDy, 2 * nClass + 1, 0, SOBOL, dCopulaCube);

  /// run MC simulation and evaluates the caplet

  for (i = 3; i < nSim; i++) {

    err = GetCapletPayoff(dCopulaCube, iNDates, lEndDates, dCoverages,

                          dCoverageShort,

                          nStrikes, Strikes, i, dPayoff);

    for (j = 0; j < nStrikes; j++) {
      dPrice[j] += (dPayoff[j] - dPrice[j]) / (i - 2);
    }
  }

  /*
          for (k=0; k< nStrikes; k++) {
                          dPrice[k] *= Df;
          }
  */

  /// computes the BS smile for the caplet

  (*ChiError) = 0.0;
  for (i = 0; i < nStrikes; i++) {

    err = srt_f_optimpvol(dPrice[i], dFwd, Strikes[i], Maturity, 1.0, SRT_CALL,
                          SRT_LOGNORMAL, &ImplVol);

    (*ChiError) += (ImplVol - BSVols[i]) * (ImplVol - BSVols[i]);
  }

  free_dvector(dPayoff, 0, nStrikes - 1);
  free_dvector(dPrice, 0, nStrikes - 1);
  return err;
}

/*-----------------------------------------------------------------------------------
 */
/*        Function that returns the payoff of a fwd swpation conditional to the
 */
/*          realization of the spanning set of swap rates */
/*-----------------------------------------------------------------------------------
 */

Err GetCapletPayoff(double **dCopulaCube, int iNDates, Date *lEndDates,
                    double *dCoverages, double dCoverageShort,

                    int nStrikes, double *Strikes, int iMC, double *dPayoff) {

  Err err = NULL;
  double *DfSet = NULL, dLibor, Df = 0.0;
  int i;

  DfSet = dvector(0, iNDates - 1);

  err =
      GetSpanningSetDf(dCopulaCube, iNDates, lEndDates, dCoverages, iMC, DfSet);

  // interpolates the simulated Discount factors

  Df = log(DfSet[0]) / dCoverages[0] * dCoverageShort;
  Df = exp(Df);
  dLibor = 1. / Df - 1.;

  for (i = 0; i < nStrikes; i++) {
    dPayoff[i] = DMAX(0.0, 1. - (1. + Strikes[i] * dCoverageShort) * Df);
  }

  free_dvector(DfSet, 0, iNDates - 1);
  return err;
}
