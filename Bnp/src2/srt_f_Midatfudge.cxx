/* ===============================================================================

   FILENAME:       srt_f_MidatFudge.cxx

   PURPOSE:        Pricing of a Midat with a Fudge

   ===============================================================================
 */
#include "srt_h_Midatfudge.h"
#include "swp_h_all.h"
#include "utallhdr.h"

Err interp_MIDAFUDGE_vol(char *str, VOL_TYPE *val) {
  if (!str)
    return "Empty string in interp_MIDAFUDGE_vol";
  strupper(str);
  strip_white_space(str);

  if (!(strcmp(str, "LOGNORMAL") && strcmp(str, "LOG"))) {
    *val = VOLLOGNORMAL;
    return 0;
  }
  if (!(strcmp(str, "NORMAL") && strcmp(str, "NORMAL"))) {
    *val = VOLNORMAL;
    return 0;
  }
  if (!(strcmp(str, "BETA") && strcmp(str, "BETAVOL"))) {
    *val = VOLBETA;
    return 0;
  }
  if (!(strcmp(str, "SABR") && strcmp(str, "SABRVOL"))) {
    *val = VOLSABR;
    return 0;
  }

  return serror("unknown vol in Midatfudge. %s", str);
}

Err interp_MIDAFUDGE_model(char *str, MIDATFUDGE_MODEL *val) {
  if (!str)
    return "Empty string in interp_MIDAFUDGE_model";
  strupper(str);
  strip_white_space(str);

  if (!strcmp(str, "SLIDDING")) {
    *val = SLIDDING;
    return 0;
  }

  if (!strcmp(str, "CONVERGING")) {
    *val = SLIDDING;
    return 0;
  }

  if (!strcmp(str, "BASIC")) {
    *val = SLIDDING;
    return 0;
  }

  return serror("unknown Midatfudge model. %s", str);
}

Err srt_f_MidatFudge(long nNumExercise, long *plExercise,
                     long *plExercisePayDates, double *pdExercisePremium,
                     long lSwapStart, long lSwapTheoEnd, char *szSwapFreq,
                     char *szSwapBasis, double dCouponRate,
                     double dFundingMargin, String szRefRate, String szPayRec,
                     String szMarket, String szYieldcurve,
                     Err (*pfGetVol)(double dStart, double dEnd, double dStrike,
                                     double dForward, double dSpread,
                                     double *pdBsVol),
                     String szVolType, String szModel,
                     /* Outputs from the addin */
                     double ***pdAnswer) {
  Err err = NULL;
  MIDAT_FUDGE midat;

  err = srt_f_initMIDAT_FUDGE(
      nNumExercise, plExercise, plExercisePayDates, pdExercisePremium,
      lSwapStart, lSwapTheoEnd, dCouponRate, szSwapBasis, szSwapFreq, szRefRate,
      szYieldcurve, pfGetVol, szVolType, szModel, &midat);
  if (err)
    return err;

  return err;
}

/* Initialisation of the Big structure */
Err srt_f_initMIDAT_FUDGE(long nNumEx, long *Ex, long *ExPay, double *ExPremium,
                          long lStart, long lTheoEnd, double strike,
                          String szSwapBasis, String szSwapFreq,
                          String szRefRate, String szYieldcurve,
                          Err (*pfGetVol)(double dStart, double dEnd,
                                          double dStrike, double dForward,
                                          double dSpread, double *pdBsVol),
                          String szVolType, String szModel,
                          MIDAT_FUDGE *midat) {
  Err err = NULL;
  SrtBasisCode SwapBasis;
  SrtCompounding SwapFreq;
  SwapDP Swapdp;
  int spotlag;
  SrtDateList date_list;
  long index;
  long first_period_start;
  double forward;
  long *plFixing = NULL;
  long i, j;
  long nTempNum;
  long *tEx = NULL;
  long *tExPay = NULL;
  long *NumSwaption = NULL;
  long *MostExp = NULL;
  double *ExBound = NULL;
  double *tExPr = NULL;
  VOL_TYPE voltype;
  MIDATFUDGE_MODEL model;
  SWAPTION_MIDAT **Swaption;

  /* Interpret the Basis and the Freq */
  err = interp_basis(szSwapBasis, &SwapBasis);
  if (err)
    return err;

  err = interp_compounding(szSwapFreq, &SwapFreq);
  if (err)
    return err;

  /* Allocate Memory for the Structure */
  (*midat).ntEx = nNumEx;
  (*midat).SwapFreq = SwapFreq;
  (*midat).SwapBasis = SwapBasis;

  tEx = lvector(1, nNumEx);
  tExPay = lvector(1, nNumEx);
  NumSwaption = lvector(1, nNumEx);
  MostExp = lvector(1, nNumEx);
  tExPr = dvector(1, nNumEx);
  ExBound = dvector(1, nNumEx);

  for (i = 1; i <= nNumEx; i++) {
    tEx[i] = Ex[i];
    tExPay[i] = ExPay[i];
    tExPr[i] = ExPremium[i];
  }

  (*midat).tEx = tEx;
  (*midat).tExPay = tExPay;
  (*midat).tExPremium = tExPr;
  (*midat).nSwaption = NumSwaption;
  (*midat).MostExp = MostExp;

  strcpy((*midat).refrate, szRefRate);

  err = interp_MIDAFUDGE_vol(szVolType, &voltype);
  if (err)
    return err;

  err = interp_MIDAFUDGE_model(szModel, &model);
  if (err)
    return err;

  (*midat).model = model;

  /* Initialise the Schedule for the swap */
  err = swp_f_initSwapDP(lStart, lTheoEnd, szSwapBasis, szSwapFreq, &Swapdp);
  if (err)
    return err;

  /* Get the spot lag from the RefRate and attach it to the SwapDP */
  err = srt_f_get_spot_lag_from_refrate(szRefRate, &spotlag);
  if (err)
    return err;

  Swapdp.spot_lag = spotlag;

  /* Builds the schedule of the underlying swap */
  date_list = SwapDP_to_DateList(&Swapdp, MODIFIED_SUCCEEDING);

  /* allocation of the Fixing Array */
  plFixing = lvector(0, date_list.len - 2);

  /* Fix the start and the first fixing */
  for (i = 0; i <= date_list.len - 2; i++) {
    first_period_start = date_list.date[i];
    plFixing[i] =
        add_unit(first_period_start, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);
  }

  index = 0;

  /* Allocation of the memory */
  Swaption = (SWAPTION_MIDAT **)malloc(nNumEx * sizeof(SWAPTION_MIDAT *));

  for (i = 1; i <= nNumEx; i++) {
    while ((plFixing[index] <= tEx[i]) && (index < date_list.len - 2))
      index++;

    /* From the remaining start date and fixing date we gonna
       find the relevant swaption */
    nTempNum = date_list.len - 1 - index;
    (*midat).nSwaption[i] = nTempNum;

    if (nTempNum == 0)
      Swaption[i] = NULL;
    else {
      Swaption[i] = (SWAPTION_MIDAT *)malloc(nTempNum * sizeof(SWAPTION_MIDAT));
      for (j = 1; j <= nTempNum; j++) {
        Swaption[i][j].exp_time = (plFixing[index + j - 1] - tEx[i]) / 365.0;
        Swaption[i][j].start = date_list.date[index + j - 1];
        Swaption[i][j].theoend = date_list.date[date_list.len - 1];
        Swaption[i][j].cxxoupon_rate = strike;

        /* compute the ForwardSwapRate today */
        err = swp_f_ForwardRate(date_list.date[index + j - 1],
                                date_list.date[date_list.len - 1], szSwapFreq,
                                szSwapBasis, szYieldcurve, szRefRate, &forward);
        if (err)
          return err;

        Swaption[i][j].fwdswap = forward;
        Swaption[i][j].voltype = voltype;
      }
    }
  }

  /* Attach the MIDAT_SWAPTION information */
  (*midat).Swaption = Swaption;

  /* Extracts the relevant fixing date for the deal */

  /* Allocate Memory for the Swaptions */

  return err;
}