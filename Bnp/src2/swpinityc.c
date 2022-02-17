/*********************************************************************
 *      Name: SwpInitYC.c                                             *
 *  Function: Add a yield curve object to SrtMkt structure            *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 20/10/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 20/10/95 FOS     Created for SORT5-GRFN3 port to NT for SPG        *
 **********************************************************************/
#include "SwpAccess.h"

char *SwpInitYC(char *ycName, char *ycType, char *ccyName, long today,
                long spotdate) {
  char *err = NULL;
  SrtCcyParam *ccy_param, *temp;
  SrtCurveListPtr curve_list;

  /* Get global pointer to curve list */
  curve_list = get_curve_list();
  if (!curve_list)
    return serror("No Curve List has been defined: call SrtInit first");

  /* Gets the Currency Parameters */
  ccy_param = new_CcyParam();
  err = swp_f_get_CcyParam_from_CcyStr(ccyName, &temp);
  if (err)
    return err;
  memcpy(ccy_param, temp, sizeof(SrtCcyParam));

  /* Add a YC object to the curve list */
  err = swp_f_addcurvetolist(
      curve_list, ycName, "YIELD_CURVE", ycType, ccyName, today, spotdate,
      ccy_param, /* Can be used to pass ccy_prm or ccy_str*/
      NULL);

  return err;

} /* END SrtInitYC(...) */
