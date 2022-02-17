/* ==========================================================================
FILENAME: SrtInitEQUnd.C

PURPOSE: Initialise an IR underlying and stores it in the underlying list
========================================================================== */

#include "SrtAccess.h"
#include "grf_h_all.h"
#include "srt_h_all.h"

char *SrtInitEQUnd(char *undName, char *model, double spot, char *discName,

                   char *dividendName, char *repoName,

                   int volCrvRows, int volCrvCols, double **volCrvVals,

                   /* SRVG parameter */
                   double omega, double beta,
                   double gamma, /* gamma-smile parameters */

                   double voldrift, double vovol,
                   double rho /* eq. vol parameters */
) {
  Err err;
  SrtCurvePtr crv, repocrv;
  TermStruct *ts;
  SrtMdlType modelType;
  SrtMdlDim modelDim;
  long today;
  char *ccy;
  SrtUndListPtr und_list;
  /* For Eq Stoch Rates */
  SrtUndPtr dom_und;

  /* remobe extra information from the name and make it uppercase */
  strupper(undName);
  strip_white_space(undName);

  /* Check the model : BLACK_SCHOLES (defaulted if no model is passed) or
   * FX_STOCH_RATES */
  if (!strcmp(model, ""))
    model = "BS";
  if (err = srt_f_interp_model(model, &modelType, &modelDim)) {
    return serror(err);
  }

  if ((modelType == EQ_STOCH_RATES) || (modelType == EQ_STOCH_RATES_SRVGS)) {
    /* The first name corresponds to an Interest Rate Underlying */
    dom_und = lookup_und(discName);
    if (!dom_und)
      return serror(
          "With Eq_Stoch_Rates  , the first curve should be an UNDERLYING");
    if (!ISUNDTYPE(dom_und, INTEREST_RATE_UND))
      return serror(
          "With Eq_Stoch_Rates  , the first curve should be an IR UNDERLYING");

    /* Get its currency */
    ccy = get_underlying_ccy(dom_und);

    /* Get today from it */
    today = get_today_from_underlying(dom_und);
  } else {
    /* Get the Discount Curve */
    if ((crv = lookup_curve(discName)) == NULL) {
      return serror("Fatal: (SrtInitEQUnd) Cannot find discount curve %s",
                    discName);
    }
    if (!ISCURVETYPE(crv, YIELD_CURVE)) {
      return serror(
          "Fatal: (SrtInitEQUnd) discount curve must be a yield curve curve");
    }

    /* Get today and the currency from the discount curve */
    today = get_today_from_curve(crv);
    ccy = get_curve_ccy(crv);
  }

  /* Get the Growth Curve */
  if ((crv = lookup_curve(dividendName)) == NULL) {
    return serror("Fatal: (SrtInitEQUnd) Cannot find dividend curve %s",
                  dividendName);
  }
  if (!ISCURVETYPE(crv, YIELD_CURVE) && !ISCURVETYPE(crv, DVD_CURVE)) {
    return serror("Fatal: (SrtInitEQUnd) growth curve must be a yield curve or "
                  "a dividend curve");
  }

  if ((repocrv = lookup_curve(repoName)) == NULL) {
    return serror("Fatal: (SrtInitEQUnd) Cannot find repo curve %s", repoName);
  }

  if (!ISCURVETYPE(repocrv, REPO_CURVE)) {
    return serror("Fatal: (SrtInitEQUnd) need a repo curve");
  }

  /* Initialise the Underlying Term Structure of Volatilities */
  if (err = srt_f_init_EQ_TermStruct(today, volCrvVals, volCrvCols, volCrvRows,
                                     undName, discName, modelType,

                                     omega, beta, gamma, voldrift, vovol, rho,
                                     &ts)) {
    return serror(err);
  }

  /* Get the underlying list and check it has been initialised  */
  und_list = get_underlying_list();
  if (und_list == NULL)
    return serror("No Underlying list defined: call SrtInit before");

  /*  Puts the Underlying in the Market List */
  err = srt_f_addundtolist(und_list, undName, "EQUITY_UND", ccy, model,
                           discName, dividendName, repoName, ts, spot);
  if (err) {
    return serror("Fatal: (SrtInitEQUnd) Failed to add EQ underlying");
  }

  /* return a success message */
  return NULL;

} /* END char* SrtInitEQUnd(...) */
