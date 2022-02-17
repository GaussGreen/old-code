/* ==========================================================================
   FILENAME:  SrtInitFXUnd.C

   PURPOSE:   Initialise an FX underlying and stores it in the underlying list
              The Underlying name is optional: if none is passed        , the
   default name will come from the currencies        , as "FOR/DOM".

   ========================================================================== */
#include "SrtAccess.h"
#include "grf_h_all.h"
#include "srt_h_all.h"
#include "swp_h_curve_struct.h"

char *SrtInitFXUnd(char *undName, /* Optional */
                   double spot, char *model, char *domDiscName,
                   char *forDiscName, int volCrvRows, int volCrvCols,
                   double **volCrvVals,

                   /* OUTPUT */
                   char *definite_name) {
  Err err;
  SrtUndListPtr und_list;
  SrtCurvePtr crv;
  SrtMdlType modelType, dom_mdl_type;
  SrtMdlDim modelDim;
  TermStruct *ts;
  long today;
  char *dom_ccy;
  char *for_ccy;

  /* For Jumping Numeraire */
  SrtUndPtr dom_und;
  SrtUndPtr for_und;

  /* Check the model : BLACK_SCHOLES (defaulted if no model is passed) or
   * FX_STOCH_RATES */
  if (!strcmp(model, ""))
    model = "BS";
  if (err = srt_f_interp_model(model, &modelType, &modelDim)) {
    return serror(err);
  }

  /* Get the underlyings/curves currencies and check for Curve/Model consistency
   */
  if (modelType == BLACK_SCHOLES || modelType == NORMAL_BS) {
    /* The first name corresponds to a yield curve */
    crv = lookup_curve(domDiscName);
    if (!crv)
      return serror(
          "With FX Black-Scholes        , the first curve should be a YC");
    if (!ISCURVETYPE(crv, YIELD_CURVE))
      return serror(
          "With FX Black-Scholes        , the first curve should be a YC");
    /* Get its currency */
    dom_ccy = get_curve_ccy(crv);
    /* Get today from it */
    today = get_clcndate_from_yldcrv(crv);

    /* The second name corresponds to a yield curve */
    crv = lookup_curve(forDiscName);
    if (!crv)
      return serror(
          "With FX Black-Scholes        , the second curve should be a YC");
    if (!ISCURVETYPE(crv, YIELD_CURVE))
      return serror(
          "With FX Black-Scholes        , the second curve should be a YC");
    /* Get its currency */
    for_ccy = get_curve_ccy(crv);
  } else if (modelType == FX_STOCH_RATES) {
    /* The first name corresponds to an Interest Rate Underlying */
    dom_und = lookup_und(domDiscName);
    if (!dom_und)
      return serror("With Fx_Stoch_Rates        , the first curve should be an "
                    "UNDERLYING");
    if (!ISUNDTYPE(dom_und, INTEREST_RATE_UND))
      return serror("With Fx_Stoch_Rates        , the first curve should be an "
                    "UNDERLYING");

    /* Get its currency */
    dom_ccy = get_underlying_ccy(dom_und);
    err = get_underlying_mdltype(dom_und, &dom_mdl_type);

    /* Get today from it */
    today = get_today_from_underlying(dom_und);

    /* The second name corresponds to an Interest Rate Underlying */
    for_und = lookup_und(forDiscName);

    if (!for_und)
      return serror("With Fx_Stoch_Rates        , the second curve should be "
                    "an UNDERLYING");

    if (!ISUNDTYPE(for_und, INTEREST_RATE_UND))
      return serror("With Fx_Stoch_Rates        , the second curve should be "
                    "an UNDERLYING");

    /* Get its currency */
    for_ccy = get_underlying_ccy(for_und);

  } /* END if (modelType == FX_STOCH_RATES) */

  /* Build the default name: "FOR/DOM" or the UndName if passed*/
  if ((!undName) || (!strcmp(undName, ""))) {
    sprintf(definite_name, "%s/%s", for_ccy, dom_ccy);
  } else {
    strcpy(definite_name, undName);
  }

  /* Initialise the Underlying Term Structure of Volatilities */
  if (err = srt_f_init_FX_TermStruct(today, volCrvVals, volCrvCols, volCrvRows,
                                     modelType, 0.0, 0.0, definite_name,
                                     domDiscName, forDiscName, &ts))

  {
    return err;
  }

  /* Get the underlying list and check it has been initialised  */
  und_list = get_underlying_list();
  if (und_list == NULL)
    return serror("No Underlying list defined: call SrtInit before");

  /* Puts the Underlying in the Market List */
  err = srt_f_addundtolist(und_list, definite_name, "FX_UND", dom_ccy, model,
                           domDiscName, forDiscName, NULL, ts, spot);
  if (err)
    return serror("Fatal: (SrtInitFXUnd) Failed to add FX underlying");

  /* Return a success message */
  return NULL;

} /* END char *SrtInitFXUnd(...) */
