/**********************************************************************
 *      Name: SrtGrfMidat.cxx                                           *
 *  Function: Entry point to GRF_MIDAT with raw data                  *
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

char *SrtGrfMidat(int numParams, char **paramStrings, char **valueStrings,
                  char *undName, int num_exercise_dates, long *exercise_dates,
                  long *exercise_start_dates, double *exercise_premiums,
                  int num_prod_dates, long *prod_dates, double *prod_cfs,
                  char *recPayStr, double *price) {
  Err err = NULL;
  SrtUndPtr und;
  SrtGrfnParam grfnparam;
  SrtReceiverType rec_pay;

  /* Overwrite defaults with user defined parameters */

  if (err = srt_f_set_GrfnParams(numParams, paramStrings, valueStrings,
                                 &grfnparam)) {
    return err;
  }

  err = interp_rec_pay(recPayStr, &rec_pay);
  if (err) {
    return err;
  }

  if (!(und = lookup_und(undName))) {
    return serror("gr_midat: Could not find underlying in market list");
  }

  err = grfn_midat_clsdfrm(num_exercise_dates, exercise_dates,
                           exercise_start_dates, exercise_premiums,
                           num_prod_dates, prod_dates, prod_cfs, rec_pay, und,
                           &grfnparam, price);

  return err;
}
