/**************************************************************************************************
 *																								  *
 *		Function: srt_f_rangeswap
 **
 *																								  *
 *		Calculates the PV of the fixed leg of a  range swap with minimum
 *daily coupon and		  * an approximation of the PV when there is a
 *minimum PERIOD coupon.						  *
 *																								  *
 *		Note that the daily PV is always bigger than the period PV. The
 *two prices computed here  * are upper bounds to the real PV.
 **
 *																								  *
 *		If the TimeSwap flag is set to YES the PV of the fixed leg of a
 *regular TimeSwap		  * is computed.
 **
 *																								  *
 *		see the FIRST Document : Time Swap and Range Swap with a Minimum
 *Coupon					  * for details.
 **
 *																								  *
 *		Author: Ezra Nahum
 ** Last modified Oct 99 1999.
 **
 *																								  *
 ***************************************************************************************************/
#include "OPFNCTNS.H"
#include "UTALLHDR.H"
#include "math.h"
#include "num_h_gamma.h"
#include "srt_h_all.h"
#include "srt_h_resetable.h"
#include "swp_h_cms.h"
#include <swp_h_all.h"
#include <swp_h_cmsopt.h"
#include <swp_h_external_fct.h"
#include <swp_h_vol.h"

/********************************************************************************
FraDBspread computes  , for a given index  , the price of the double barrier
digital (in fact the difference of two callspreads with strikes located at the
lower and upper barriers of the Range Swap). The Cms adjustment is made.
*********************************************************************************/
Err get_fra(long Fradate, int spotlag, char *cRefRateCode,
            char *szYieldCurveName, double *dFra);

Err FraDBspread(long today, long Startdate, long paydate, int spotlag, int lag,
                double lowerstrike, double upperstrike, double epsilon,
                double lowvolshift, double upvolshift, char *szYieldCurveName,
                char *szVolCurveName, SrtDiffusionType CMSlognorm,
                char *cMarketId, char *cRefRateCode,
                char *(*getLogVol)(char *szVolCurve, double dStartDate,
                                   double dEndDate, double dStrike,
                                   double *dVol, double *dPower),
                elemforCms passCms, double *DBspread) {

  double lowercall, lowercalleps, uppercall, uppercalleps;
  long enddate, theoenddate;
  double dFra, dCms, delay;
  double lowvol, lowvoleps, upvol, upvoleps;
  double lognorm;
  double maturity;
  double num_of_months;
  double lowshift;
  double upshift;
  double dspread;
  long FraSpot;
  long Fradate;
  SwapDP sdpFra;
  SrtBasisCode float_basis;
  SrtCompounding float_compounding;

  Err err = NULL;

  /* Calculate the vol type.  This fixes a previous bug wherein the vol type was
   * assumed to be lognormal */
  (*getLogVol)(szVolCurveName, Startdate, paydate, upperstrike, &upvol,
               &lognorm);

  /*Get the fixing date*/
  Fradate = add_unit(Startdate, -lag, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* start date of the index in question*/
  FraSpot = add_unit(Fradate, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* get info from the RefRate Code*/
  if (err = swp_f_get_ref_rate_details(cRefRateCode, &float_basis,
                                       &float_compounding))
    smessage("Error if swp_f_get_ref_rate_details");

  num_of_months = 12 / float_compounding;

  /* maturity of the fra */
  maturity = coverage(today, Fradate, BASIS_ACT_365);

  /* delay  , used for computation of the CMS rate*/
  delay = coverage(FraSpot, paydate, BASIS_ACT_365);

  /* enddate of the fra and theoenddate*/

  theoenddate =
      add_unit(FraSpot, (int)num_of_months, SRT_MONTH, NO_BUSDAY_CONVENTION);

  enddate = add_unit(theoenddate, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* Compute the Fra */
  /* first initialise the sdpFra */
  err = swp_f_setSwapDP(FraSpot, theoenddate, float_compounding, float_basis,
                        &sdpFra);
  if (err)
    return err;
  /* input the spot lag  */
  sdpFra.spot_lag = spotlag;

  /* computation of the Fra */
  err =
      swp_f_ForwardRate_SwapDP(&sdpFra, szYieldCurveName, cRefRateCode, &dFra);

  /* computes the spread: for CMS rate*/
  dspread = swp_f_spread(today, Fradate, cRefRateCode);

  /*Computes the CMSRate*/
  if (err = swp_f_Cms_Rate(
          dFra, maturity, 1, (double)float_compounding, delay, 365.0 / 360.0,
          CMSlognorm, 0.0, 1 /* Full Smile Approx */, FraSpot, theoenddate,
          passCms.AdjForSpread, passCms.AdjForSpread * dspread, szVolCurveName,
          passCms.NumStrikesInVol, passCms.Strikes, &dCms))

    smessage("Error in swp_f_cmsrate");

  /* computation of the Vols for the Call Spreads*/

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         lowerstrike, &lowvol, &lognorm))
    smessage("Error in swp_f_vol");

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         lowerstrike - epsilon, &lowvoleps, &lognorm))
    smessage("Error in swp_f_vol");

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         upperstrike, &upvol, &lognorm))
    smessage("Error in swp_f_vol");

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         upperstrike + epsilon, &upvoleps, &lognorm))
    smessage("Error in swp_f_vol");

  /* The vol is shifted so as to be conservative*/

  if ((dCms > lowerstrike) && (dCms < upperstrike)) {
    lowshift = -lowvolshift;
    upshift = -upvolshift;
  }

  else {
    lowshift = lowvolshift;
    upshift = upvolshift;
  }

  if (lowerstrike <= epsilon)
    DBspread[0] = 1;
  else {

    if (lognorm == 0.0) {
      lowercall = srt_f_optblknrm(dCms, lowerstrike, lowvol + lowshift,
                                  maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps =
          srt_f_optblknrm(dCms, lowerstrike - epsilon, lowvoleps + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);
    } else {
      lowercall = srt_f_optblksch(dCms, lowerstrike, lowvol + lowshift,
                                  maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps =
          srt_f_optblksch(dCms, lowerstrike - epsilon, lowvoleps + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);
    }

    DBspread[0] = (lowercalleps - lowercall) / epsilon;
  }

  if (lognorm == 0.0) {
    uppercall = srt_f_optblknrm(dCms, upperstrike, upvol + upshift, maturity,
                                1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblknrm(dCms, upperstrike + epsilon, upvoleps + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);
  } else {
    uppercall = srt_f_optblksch(dCms, upperstrike, upvol + upshift, maturity,
                                1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblksch(dCms, upperstrike + epsilon, upvoleps + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);
  }

  DBspread[0] = DBspread[0] - (uppercall - uppercalleps) / epsilon;

  return err;
}

Err FraDBspreadFloating(
    long today, long Startdate, long paydate, long StartDateIdx, int spotlag,
    int lag, double lowerstrike, double upperstrike, double epsilon,
    double lowvolshift, double upvolshift, char *AccRefRateCode, double corrup,
    double corrdown, char *szYieldCurveName, char *szVolCurveName,
    SrtDiffusionType CMSlognorm, char *cMarketId, char *cRefRateCode,
    char *(*getLogVol)(char *szVolCurve, double dStartDate, double dEndDate,
                       double dStrike, double *dVol, double *dPower),
    elemforCms passCms, double *DBspread) {

  double lowercall, lowercalleps, uppercall, uppercalleps;
  long enddate, theoenddate, enddateacc, theoenddateacc;
  double dFra, dFraacc, dCms;
  double lowvol, lowvoleps, upvol, upvoleps;
  double atmvolindex;
  double atmvolacc;
  double correl, delay;
  double lognorm;
  double lognormacc;
  double maturity, maturityacc;
  double num_of_months, num_of_monthsacc;
  double lowshift;
  double upshift;
  double lowvolacc, upvolacc;
  double correlcorrec;
  double dspread;
  long FraSpot;
  long Fradate;
  SwapDP sdpFra;
  SwapDP sdpFraacc;
  SrtBasisCode float_basis;
  SrtCompounding float_compounding;
  SrtBasisCode float_basisacc;
  SrtCompounding float_compoundingacc;
  Err err = NULL;

  /*get the fixing date*/
  Fradate = add_unit(StartDateIdx, -lag, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* start date of the index Fra*/
  FraSpot = add_unit(Fradate, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);

  /*Get details for the index Fra and the Accrual Fra*/
  if (err = swp_f_get_ref_rate_details(cRefRateCode, &float_basis,
                                       &float_compounding)) {
    smessage("Error in swp_f_get_ref_rate_details");
  }

  if (err = swp_f_get_ref_rate_details(AccRefRateCode, &float_basisacc,
                                       &float_compoundingacc)) {
    smessage("Error in swp_f_get_ref_rate_details");
  }

  num_of_months = 12 / float_compounding;

  num_of_monthsacc = 12 / float_compoundingacc;

  maturity = coverage(today, Fradate, BASIS_ACT_365);

  theoenddate =
      add_unit(FraSpot, (int)num_of_months, SRT_MONTH, NO_BUSDAY_CONVENTION);

  enddate = add_unit(theoenddate, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

  /* computes the delay  , for use in the CMS rate computation*/
  delay = coverage(FraSpot, paydate, BASIS_ACT_365);

  /* Compute the Fwd */
  /* first initialise the sdpFra */
  err = swp_f_setSwapDP(FraSpot, theoenddate, float_compounding, float_basis,
                        &sdpFra);
  if (err)
    return err;
  /* input the spot lag  */
  sdpFra.spot_lag = spotlag;

  /* computation of the Fwd */
  err =
      swp_f_ForwardRate_SwapDP(&sdpFra, szYieldCurveName, cRefRateCode, &dFra);

  if (err)
    return err;

  /* computes the spread  , for use in the CMS rate computation*/
  dspread = swp_f_spread(today, Fradate, cRefRateCode);

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         dFra, &atmvolindex, &lognorm)) {
    smessage("Error in swp_f_vol");
  }

  /* Computes the CMS rate */
  if (err = swp_f_Cms_Rate(
          dFra, maturity, 1, (double)float_compounding, delay, 365.0 / 360.0,
          CMSlognorm, 0.0, 1 /* Full Smile Approx */, FraSpot, theoenddate,
          passCms.AdjForSpread, passCms.AdjForSpread * dspread, szVolCurveName,
          passCms.NumStrikesInVol, passCms.Strikes, &dCms))

  {
    smessage("Error in swp_f_cmsrate");
  }

  /* start date of the Accrual fra  , and elements*/

  theoenddateacc = add_unit(Startdate, (int)num_of_monthsacc, SRT_MONTH,
                            NO_BUSDAY_CONVENTION);

  enddateacc = add_unit(theoenddateacc, 0, SRT_BDAY, MODIFIED_SUCCEEDING);

  maturityacc = coverage(
      today, add_unit(Startdate, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING),
      BASIS_ACT_365);

  /* Compute the AccFwd */
  /* first initialise the sdpFra */
  err = swp_f_setSwapDP(Startdate, theoenddateacc, float_compoundingacc,
                        float_basisacc, &sdpFraacc);
  if (err)
    return err;
  /* input the spot lag  */
  sdpFraacc.spot_lag = spotlag;

  /* computation of the AccFwd */
  err = swp_f_ForwardRate_SwapDP(&sdpFraacc, szYieldCurveName, AccRefRateCode,
                                 &dFraacc);

  /* Computes the vols necessary for the calculations of the call spreads*/

  if (err =
          (*getLogVol)(szVolCurveName, (double)(Startdate), (double)enddateacc,
                       dFraacc, &atmvolacc, &lognormacc)) {
    smessage("Error in swp_f_vol");
  }

  if (err =
          (*getLogVol)(szVolCurveName, (double)(Startdate), (double)enddateacc,
                       lowerstrike, &lowvolacc, &lognormacc)) {
    smessage("Error in swp_f_vol");
  }

  if (err =
          (*getLogVol)(szVolCurveName, (double)(Startdate), (double)enddateacc,
                       upperstrike, &upvolacc, &lognormacc)) {
    smessage("Error in swp_f_vol");
  }

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         lowerstrike, &lowvol, &lognorm)) {
    smessage("Error in swp_f_vol");
  }

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         lowerstrike - epsilon, &lowvoleps, &lognorm)) {
    smessage("Error in swp_f_vol");
  }

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         upperstrike, &upvol, &lognorm)) {
    smessage("Error in swp_f_vol");
  }

  if (err = (*getLogVol)(szVolCurveName, (double)(FraSpot), (double)enddate,
                         upperstrike + epsilon, &upvoleps, &lognorm)) {
    smessage("Error in swp_f_vol");
  }

  /* correlation interpolation*/
  correl = (corrdown - corrup) * ((double)(FraSpot - Startdate)) /
               ((double)(paydate - Startdate)) +
           corrup;

  /* The vol is shifted so as to be conservative*/

  if ((dFra > lowerstrike) && (dFra < upperstrike)) {
    lowshift = -lowvolshift;
    upshift = -upvolshift;
  }

  else {
    lowshift = lowvolshift;
    upshift = upvolshift;
  }

  /* Index Lognormal  , Accrual Lognormal case*/
  if ((lognorm == 1.0) && (lognormacc == 1.0)) {

    if (lowerstrike <= epsilon) {
      DBspread[0] = 1;
    } else {
      correlcorrec =
          exp(correl * lowvol * lowvolacc * DMIN(maturityacc, maturity));

      lowercall =
          srt_f_optblksch(dCms * correlcorrec, lowerstrike, lowvol + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps = srt_f_optblksch(dCms * correlcorrec, lowerstrike - epsilon,
                                     lowvoleps + lowshift, maturity, 1.0,
                                     SRT_CALL, PREMIUM);

      DBspread[0] = (lowercalleps - lowercall) / epsilon;
    }

    correlcorrec = exp(correl * upvol * upvolacc * DMIN(maturityacc, maturity));

    uppercall =
        srt_f_optblksch(dCms * correlcorrec, upperstrike, upvol + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblksch(dCms * correlcorrec, upperstrike + epsilon,
                        upvoleps + upshift, maturity, 1.0, SRT_CALL, PREMIUM);

    DBspread[0] = DBspread[0] - (uppercall - uppercalleps) / epsilon;
  }

  /* Index Normal  , Accrual Lognormal case*/
  if ((lognorm == 0.0) && (lognormacc == 1.0)) {

    if (lowerstrike <= epsilon) {
      DBspread[0] = 1;
    } else {
      correlcorrec = correl * lowvol * lowvolacc * DMIN(maturityacc, maturity);

      lowercall =
          srt_f_optblknrm(dCms + correlcorrec, lowerstrike, lowvol + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps = srt_f_optblknrm(dCms + correlcorrec, lowerstrike - epsilon,
                                     lowvoleps + lowshift, maturity, 1.0,
                                     SRT_CALL, PREMIUM);

      DBspread[0] = (lowercalleps - lowercall) / epsilon;
    }

    correlcorrec = correl * upvol * upvolacc * DMIN(maturityacc, maturity);

    uppercall =
        srt_f_optblknrm(dCms + correlcorrec, upperstrike, upvol + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblknrm(dCms + correlcorrec, upperstrike + epsilon,
                        upvoleps + upshift, maturity, 1.0, SRT_CALL, PREMIUM);

    DBspread[0] = DBspread[0] - (uppercall - uppercalleps) / epsilon;
  }

  /* Index Normal  , Accrual Normal case*/
  if ((lognorm == 0.0) && (lognormacc == 0.0)) {

    if (lowerstrike <= epsilon) {
      DBspread[0] = 1;
    } else {
      lowercall = srt_f_optblknrm(dCms, lowerstrike, lowvol + lowshift,
                                  maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps =
          srt_f_optblknrm(dCms, lowerstrike - epsilon, lowvoleps + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);

      DBspread[0] =
          (lowercalleps - lowercall) / epsilon +
          correl * lowvolacc * DMIN(maturity, maturityacc) *
              exp(-0.5 * ((dCms - lowerstrike) / (lowvol * sqrt(maturity))) *
                  ((dCms - lowerstrike) / (lowvol * sqrt(maturity)))) /
              (sqrt(maturity) * sqrt(2 * SRT_PI) * dFraacc);
    }

    uppercall = srt_f_optblknrm(dCms, upperstrike, upvol + upshift, maturity,
                                1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblknrm(dCms, upperstrike + epsilon, upvoleps + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);

    DBspread[0] =
        DBspread[0] - (uppercall - uppercalleps) / epsilon -
        correl * upvolacc * DMIN(maturity, maturityacc) *
            exp(-0.5 * ((dCms - upperstrike) / (upvol * sqrt(maturity))) *
                ((dCms - upperstrike) / (upvol * sqrt(maturity)))) /
            (sqrt(maturity) * sqrt(2 * SRT_PI) * dFraacc);
  }

  /* Index Lognormal  , Accrual Normal case*/
  if ((lognorm == 1.0) && (lognormacc == 0.0)) {

    if (lowerstrike <= epsilon) {
      DBspread[0] = 1;
    } else {
      lowercall = srt_f_optblksch(dCms, lowerstrike, lowvol + lowshift,
                                  maturity, 1.0, SRT_CALL, PREMIUM);
      lowercalleps =
          srt_f_optblksch(dCms, lowerstrike - epsilon, lowvoleps + lowshift,
                          maturity, 1.0, SRT_CALL, PREMIUM);

      DBspread[0] =
          (lowercalleps - lowercall) / epsilon +
          correl * lowvolacc * DMIN(maturity, maturityacc) *
              exp(-0.5 *
                  (log(dCms / lowerstrike) / (lowvol * sqrt(maturity)) -
                   lowvol * sqrt(maturity) / 2) *
                  (log(dCms / lowerstrike) / (lowvol * sqrt(maturity)) -
                   lowvol * sqrt(maturity) / 2)) /
              (sqrt(maturity) * sqrt(2 * SRT_PI) * dFraacc);
    }
    correlcorrec = exp(correl * upvol * upvolacc * DMIN(maturityacc, maturity));

    uppercall =
        srt_f_optblksch(dCms * correlcorrec, upperstrike, upvol + upshift,
                        maturity, 1.0, SRT_CALL, PREMIUM);
    uppercalleps =
        srt_f_optblksch(dCms * correlcorrec, upperstrike + epsilon,
                        upvoleps + upshift, maturity, 1.0, SRT_CALL, PREMIUM);

    DBspread[0] = DBspread[0] - (uppercall - uppercalleps) / epsilon -
                  correl * upvolacc * DMIN(maturity, maturityacc) *
                      exp(-0.5 *
                          (log(dCms / upperstrike) / (upvol * sqrt(maturity)) -
                           upvol * sqrt(maturity) / 2) *
                          (log(dCms / upperstrike) / (upvol * sqrt(maturity)) -
                           upvol * sqrt(maturity) / 2)) /
                      (sqrt(maturity) * sqrt(2 * SRT_PI) * dFraacc);
  }

  return err;
}

/***********************************************
Computation of the PV of the fixed leg
providing an upper bound when Paribas pays the fixed leg
and a lower bound when Paribas receives it in order to yield
a conservative price in all cases.
When the minimum is accrued daily  , the pv is exact and constitutes
an upper bound for the period case.
************************************************/

Err PVComputeFloating(long today, long StartDate, long EndDate, int spotlag,
                      int endlag, double NumFixingDays, double *proba,
                      double *isfriday, double coupon, double minimum,
                      double corrup, double corrdown, char *FloatingCoupon,
                      char *FloatingMinimum, char *approx, char *TimeSwap,
                      char *WeekendRule, char *dailyorperiod, char *RecPay,
                      char *AccRefRateCode, char *szYieldCurveName,
                      char *szVolCurveName, char *cMarketId, double *pv)

{

  double avep; /* average of the digitals*/
  int i, j;
  double Nstar; /* minimum number of days the index has to be in the period to
                   get more than the minimum*/
  double TotalNumDays; /*Total number of days counting towards the Time Swap */
  double variance; /* variance of the number of days the index spent within the
                      barriers */
  double varp1, varp2; /*used for intermediate calculations*/
  double dFraacc;
  double m, sigma;
  double alpha, beta, lambda; /*used for Beta Approx*/
  Err err = NULL;

  /* Compute the AccFwd */
  get_fra(add_unit(StartDate, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING), spotlag,
          AccRefRateCode, szYieldCurveName, &dFraacc);

  /*initialization*/
  avep = 0;
  varp1 = 0;
  varp2 = 0;

  /***********************************************
  Test wether the Fridays count as 3 days or 1 day:
                  WeekendRule = YES: 3 days
                  WeekendRule = NO:  1 day
  ************************************************/

  if (strcmp(WeekendRule, "YES") == 0) {

    TotalNumDays = EndDate - StartDate;

  }

  else if (strcmp(WeekendRule, "NO") == 0)

  {

    TotalNumDays = NumFixingDays + endlag;

    for (i = 0; i < NumFixingDays; i++) {
      isfriday[i] = 0;
    }
  }

  /* Computation of Nstar */
  if (strcmp(FloatingCoupon, "YES") == 0) {
    if (strcmp(FloatingMinimum, "YES") == 0) {
      Nstar = (int)floor(DMAX(dFraacc + minimum, 0) * TotalNumDays /
                         (dFraacc + coupon)) +
              1;
    }

    else {
      Nstar =
          (int)floor(DMAX(minimum, 0) * TotalNumDays / (dFraacc + coupon)) + 1;
    }
  } else {
    if (strcmp(FloatingMinimum, "YES") == 0) {
      Nstar =
          (int)floor(DMAX(dFraacc + minimum, 0) * TotalNumDays / coupon) + 1;
    }

    else {
      Nstar = (int)floor(DMAX(minimum, 0) * TotalNumDays / coupon) + 1;
    }
  }

  /* computation of the weighted (or not  , depending on WeekendRule) average of
   * the digitals*/

  for (i = 0; i < NumFixingDays; i++) {

    avep += (1 + 2 * isfriday[i]) * proba[i] / TotalNumDays;
    varp1 += (1 + 2 * isfriday[i]) * (1 + 2 * isfriday[i]) * proba[i];
  }

  avep += endlag * proba[(int)NumFixingDays - 1] / TotalNumDays;

  /*************************************************************
  Computation of the PV in different cases:
  1. If it is a minimum DAILY coupon (dailyorperiod = YES)
  2. If it is a minimum PERIOD coupon (dailyorperiod = NO):
          a. an upper bound is given in the case Paribas Pays the minimum Period
  coupon (RecPay = PAY) b. a lower bound is given in the case  Paribas receives
  the minimum Period coupon (RecPay = REC)
  **************************************************************/

  if ((strcmp(TimeSwap, "YES") == 0) || (strcmp(dailyorperiod, "DAILY") == 0)) {

    pv[0] = avep;
  }

  else

  {

    if (strcmp(RecPay, "PAY") == 0) {

      for (i = 0; i < NumFixingDays; i++) {
        for (j = i + 1; j < NumFixingDays; j++) {

          /*varp2 +=
           * sqrt(proba[i]*proba[j])*(1+2*isfriday[i])*(1+2*isfriday[j]);*/
          varp2 += sqrt(proba[i] * proba[j]) *
                   (sqrt(proba[i] * proba[j]) +
                    (corrup + (corrdown - corrup) * (j - i) / TotalNumDays) *
                        sqrt((1 - proba[i]) * (1 - proba[j]))) *
                   (1 + 2 * isfriday[i]) * (1 + 2 * isfriday[j]);
          /*varp2+=
           * (-(j-i)/TotalNumDays+1)*proba[i]*(1+2*isfriday[i])*(1+2*isfriday[j]);
           */
        }
      }

      variance = varp1 + 2 * varp2 - TotalNumDays * TotalNumDays * avep * avep;

    }

    else if (strcmp(RecPay, "REC") == 0) {

      variance = TotalNumDays * avep * (1 - avep);
    }

    /* we make sure the variance  , because of numerical errors is not
     * negative*/
    variance = DMAX(variance, 0.000000001);

    if (strcmp(approx, "Gamma") == 0) {
      lambda = avep * TotalNumDays;
      alpha = lambda * avep * TotalNumDays;

      pv[0] = exp(gammln(alpha + 1) - gammln(alpha)) *
              (1 - gammp(alpha + 1, Nstar) * exp(-gammln(alpha + 1))) / lambda;
      pv[0] = pv[0] / TotalNumDays;
    }

    else if (strcmp(approx, "Normal") == 0) {

      pv[0] = sqrt(variance / (2 * SRT_PI)) *
                  exp(-(Nstar - TotalNumDays * avep) *
                      (Nstar - TotalNumDays * avep) / (2 * variance)) +
              TotalNumDays * avep *
                  norm((TotalNumDays * avep - Nstar) / sqrt(variance));
      pv[0] = pv[0] / TotalNumDays;
    }

    else if (strcmp(approx, "Beta") == 0) {
      alpha = avep *
              (avep * (1 - avep) * TotalNumDays * TotalNumDays / variance - 1);
      beta = alpha * (1 - avep) / avep;

      pv[0] = TotalNumDays * betacmp(alpha + 1, beta) *
              (1 - betaincmp(alpha + 1, beta, Nstar / TotalNumDays) /
                       betacmp(alpha + 1, beta)) /
              betacmp(alpha, beta);
      pv[0] = pv[0] / TotalNumDays;
    }

    else {
      m = 2 * log(TotalNumDays * avep) -
          0.5 * log(variance + TotalNumDays * avep * TotalNumDays * avep);
      sigma = sqrt(log(variance + TotalNumDays * avep * TotalNumDays * avep) -
                   2 * log(TotalNumDays * avep));

      pv[0] = exp(m + sigma * sigma / 2) *
              norm((m + sigma * sigma - log(Nstar)) / sigma);
      pv[0] = pv[0] / TotalNumDays;
    }
  }

  return err;
}

Err srt_f_FlooredTimeSwap(
    long Startdate, long Enddate, double coupon,
    double minimum, /* minimum coupon or minimum margin*/
    double UpBa,    /*upper barrier*/
    double LoBa,    /* lower barrier*/
    double epsilon, /* difference between the strikes used to compute the
                       callspread*/
    int lag, int endlag, double lowvolshift, double upvolshift,
    char *FloatingLeg, char *FloatingMinimum, char *AccRefRateCode,
    double corrup,   /* maximum correlation */
    double corrdown, /* minimum correlation */
    char *approx,    /*The type of approximation made on the distribution of the
                        number of days spent within the range*/
    char *TimeSwap,  /* flag: yes if it's a timeswap  , no if it's a range swap
                        with min coupon*/
    char *dailyorperiod, /* daily coupon or period coupon*/
    char *WeekendRule,   /* count fridays as 3 days or not*/
    char *RecPay, char *szYieldCurveName, char *szVolCurveName,
    SrtDiffusionType CMSlognorm, char *cMarketId, char *cRefRateCode,
    char *(*getLogVol)(char *szVolCurve, double dStartDate, double dEndDate,
                       double dStrike, double *dVol, double *dPower),
    elemforCms passCms, double *pv)

{

  double NumFixingDays; /* the number of fixing days counting towards  the Time
                           Swap*/
  double *probafix = NULL; /* will store the values of the everyday callspreads
                              when the coupon is fixed*/
  double *probafloat = NULL; /*will store the values of the everyday callspreads
                                when the coupon is floating*/
  double *isfriday = NULL; /* for each day  , isfriday will be 1 if that day is
                              a friday and 0 otherwise*/
  long firstfriday; /* Date of the first friday after the start of the Time
                       Swap*/
  int j;
  long d;
  long today;
  int spotlag;
  double dFraacc;
  double pvfix, pvfloat;
  Err err = NULL;
  SrtCrvPtr yldcrv;

  yldcrv = lookup_curve(szYieldCurveName);
  today = get_today_from_curve(yldcrv);
  spotlag = get_spotlag_from_curve(yldcrv);

  /*first test*/

  if (Startdate >= Enddate) {
    smessage("The Start Date is bigger than the End Date");
    return err;
  }

  NumFixingDays = 0;
  d = Startdate;

  /*Getting the first friday*/

  for (j = 0; j < 5; j++) {
    if (week_day(d + j) == FRIDAY) {
      firstfriday = d + j;
    }
  }

  /* computing the number of fixing days counting towards the Time Swap*/

  while (d < Enddate - endlag) {
    if (week_day(d) == FRIDAY) {
      d = d + 3;
    }

    else {
      d++;
    }

    NumFixingDays++;
  }

  /*memory allocation*/
  probafix = dvector(0, (long)NumFixingDays - 1);
  probafloat = dvector(0, (long)NumFixingDays - 1);
  isfriday = dvector(0, (long)NumFixingDays - 1);

  /*Initialization of the vectors*/

  for (j = 0; j < NumFixingDays; j++) {
    probafix[j] = 0;
    probafloat[j] = 0;
    isfriday[j] = 0;
  }

  /*fill in isfriday*/

  j = (int)(firstfriday - Startdate);

  while (j < NumFixingDays) {

    isfriday[j] = 1;

    j = j + 5;
  }

  /* Compute the AccFwd */
  /*This function takes the fixing dates as input*/
  get_fra(add_unit(Startdate, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING), spotlag,
          AccRefRateCode, szYieldCurveName, &dFraacc);

  /* fill in proba which contains  , for each fixing day of the period  ,
  the probability that the index is within the range*/

  d = Startdate;
  j = 0;
  while (d < Enddate - endlag) {

    /* callspreads under the forward neutral measure  ,  fixed coupon*/
    if (strcmp(RecPay, "PAY") == 0) {
      if (err = FraDBspread(today, d, Enddate, spotlag, lag, LoBa, UpBa,
                            epsilon / 10000, lowvolshift, upvolshift,
                            szYieldCurveName, szVolCurveName, CMSlognorm,
                            cMarketId, cRefRateCode, getLogVol, passCms,
                            &(probafix[j]))) {
        smessage("Error in FraDBspread");
        free_dvector(probafix, 0, (long)NumFixingDays - 1);
        free_dvector(isfriday, 0, (long)NumFixingDays - 1);
        return err;
      }
    } else {
      if (err = FraDBspread(today, d, Enddate, spotlag, lag, LoBa, UpBa,
                            -epsilon / 10000, lowvolshift, upvolshift,
                            szYieldCurveName, szVolCurveName, CMSlognorm,
                            cMarketId, cRefRateCode, getLogVol, passCms,
                            &(probafix[j]))) {
        smessage("Error in FraDBspread");
        free_dvector(probafix, 0, (long)NumFixingDays - 1);
        free_dvector(isfriday, 0, (long)NumFixingDays - 1);
        return err;
      }
    }

    if (week_day(d) == FRIDAY) {
      d = d + 3;
    } else {
      d++;
    }
    j++;
  }

  if (strcmp(FloatingLeg, "YES") == 0) {

    d = Startdate;
    j = 0;
    while (d < Enddate - endlag) {

      /* callspreads under the numeraire of the accrual fra*/
      if (strcmp(RecPay, "PAY") == 0) {
        if (err = FraDBspreadFloating(
                today, Startdate, Enddate, d, spotlag, lag, LoBa, UpBa,
                epsilon / 10000, lowvolshift, upvolshift, AccRefRateCode,
                corrup, corrdown, szYieldCurveName, szVolCurveName, CMSlognorm,
                cMarketId, cRefRateCode, getLogVol, passCms,
                &(probafloat[j]))) {
          smessage("Error in FraDBspreadFloating");
          free_dvector(probafloat, 0, (long)NumFixingDays - 1);
          free_dvector(isfriday, 0, (long)NumFixingDays - 1);
          return err;
        }
      } else {

        if (err = FraDBspreadFloating(
                today, Startdate, Enddate, d, spotlag, lag, LoBa, UpBa,
                -epsilon / 10000, lowvolshift, upvolshift, AccRefRateCode,
                corrup, corrdown, szYieldCurveName, szVolCurveName, CMSlognorm,
                cMarketId, cRefRateCode, getLogVol, passCms,
                &(probafloat[j]))) {
          smessage("Error in FraDBspreadFloating");
          free_dvector(probafloat, 0, (long)NumFixingDays - 1);
          free_dvector(isfriday, 0, (long)NumFixingDays - 1);
          return err;
        }
      }

      if (week_day(d) == FRIDAY) {
        d = d + 3;
      } else {
        d++;
      }
      j++;
    }
  }

  /*Computing the expectations*/
  /*first for the expectations under the forward measure*/

  if (err = PVComputeFloating(
          today, Startdate, Enddate, spotlag, endlag, NumFixingDays, probafix,
          isfriday, coupon, minimum, corrup, corrdown, FloatingLeg,
          FloatingMinimum, approx, TimeSwap, WeekendRule, dailyorperiod, RecPay,
          AccRefRateCode, szYieldCurveName, szVolCurveName, cMarketId,
          &pvfix)) {
    smessage("Error in PVComputeFloating");
    free_dvector(probafix, 0, (long)NumFixingDays - 1);
    free_dvector(isfriday, 0, (long)NumFixingDays - 1);
    return err;
  }
  /* then the expectation in the numeraire of dFraacc*/

  if (strcmp(FloatingLeg, "YES") == 0) {
    if (err = PVComputeFloating(
            today, Startdate, Enddate, spotlag, endlag, NumFixingDays,
            probafloat, isfriday, coupon, minimum, corrup, corrdown,
            FloatingLeg, FloatingMinimum, approx, TimeSwap, WeekendRule,
            dailyorperiod, RecPay, AccRefRateCode, szYieldCurveName,
            szVolCurveName, cMarketId, &pvfloat)) {
      smessage("Error in PVComputeFloating");
      free_dvector(probafloat, 0, (long)NumFixingDays - 1);
      free_dvector(isfriday, 0, (long)NumFixingDays - 1);
      return err;
    }
  }

  if (strcmp(TimeSwap, "YES") == 0) {
    if (strcmp(FloatingLeg, "YES") == 0) {
      pv[0] = dFraacc * pvfloat + coupon * pvfix;
    } else {
      pv[0] = coupon * pvfix;
    }
  } else {
    if (strcmp(FloatingLeg, "YES") == 0) {
      if (strcmp(FloatingMinimum, "YES") == 0) {
        pv[0] = DMAX(dFraacc + minimum, 0) + (coupon - minimum) * pvfix;
      } else {
        pv[0] = minimum + dFraacc * pvfloat + (coupon - minimum) * pvfix;
      }
    } else {
      if (strcmp(FloatingMinimum, "YES") == 0) {
        pv[0] = DMAX(dFraacc + minimum, 0) + (coupon - minimum) * pvfix -
                dFraacc * pvfloat;
      } else {
        pv[0] = minimum + (coupon - minimum) * pvfix;
      }
    }
  }

  /*free memory*/

  free_dvector(probafix, 0, (long)NumFixingDays - 1);
  free_dvector(probafloat, 0, (long)NumFixingDays - 1);
  free_dvector(isfriday, 0, (long)NumFixingDays - 1);

  return err;
}
