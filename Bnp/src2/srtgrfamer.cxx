/**********************************************************************
 *      Name: SrtGrfAmer.cxx                                            *
 *  Function: Entry point to GRF_AMER with raw data                   *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 09/01/96                                                *
 *--------------------------------------------------------------------*
 *    Inputs: Raw data from anywhere                                  *
 *   Returns:                                                         *
 *   Globals: Expects und and request list structures to exist        *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 09/01/96 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "SrtAccess.h"
#include "srt_h_all.h"
#include "srt_h_grfclsdfrm.h"

char *SrtGrfAmer(int numParams, char **paramStrings, char **valueStrings,
                 char *undName, double start, double end, int nfp, int delay,
                 char *cpdStr, char *basisStr, double *strikes, int num_strikes,
                 char *recPayStr, double *price)

{
  Err err = NULL;
  int status = 0;
  SrtUndPtr und;
  SrtGrfnParam grfnparam;
  SrtReceiverType rec_pay;
  int compd;
  int basis;
  char *obj_name = NULL;
  char *yc_name = NULL;
  int dts;
  long num_dates;
  double *strikes2 = NULL;
  double temp;
  int uFreeIt = 0;
  int i;

  /* Overwrite parametes with user defined parameters */

  if (err = srt_f_set_GrfnParams(numParams, paramStrings, valueStrings,
                                 &grfnparam)) {
    return err;
  }

  strupper(cpdStr);
  strupper(basisStr);

  err = interp_basis(basisStr, &basis);
  if (err) {
    return err;
  }
  err = interp_compounding(cpdStr, &compd);
  if (err) {
    return err;
  }
  err = interp_rec_pay(recPayStr, &rec_pay);
  if (err) {
    return err;
  }

  if (!(und = lookup_und(undName))) {
    return serror("gr_amer: Could not find underlying in market list");
  }

  dts = get_spotlag_from_underlying(und);

  num_dates = nfp + 1;

  if (num_strikes == 1 && num_dates > 1) {
    temp = strikes[0];
    strikes2 = calloc(num_dates, sizeof(double));
    if (!strikes2) {
      return serror("gr_amer: Memory allocation error");
    }
    num_strikes = num_dates;
    for (i = 0; i < num_strikes; i++) {
      strikes2[i] = temp;
    }
    uFreeIt = 1;
  } else {
    strikes2 = strikes;
    uFreeIt = 0;
  }

  if (num_dates != num_strikes) {
    srt_free(strikes2);
    return serror("number of dates must match number of strikes");
  }

  err = grf_ameswp_clsdfrm(start, end, nfp, delay, compd, basis, strikes2,
                           rec_pay, und, &grfnparam, price);

  if (uFreeIt) {
    srt_free(strikes2);
  }
  return err;
}
