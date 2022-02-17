/**********************************************************************
 *      Name: SwpFwdRate.cxx                                            *
 *  Function: Calculates a forward rate                               *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 08/01/96                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 08/01/96 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/

#include "SwpAccess.h"

char *SwpFwdRate(long start, long nfp_or_end, char *cpdStr, char *basisStr,
                 char *ycName, char *refRateCode, double *fwdRate) {
  Err err;

  /* Compute the Swap rate taking spreads into account */
  err = swp_f_ForwardRate(start, nfp_or_end, cpdStr, basisStr, ycName,
                          refRateCode, fwdRate);
  if (err)
    return err;

  /* Return a success message */
  return NULL;

} /* END char *SwpFwdRate(...) */
