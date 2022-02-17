/**********************************************************************
 *      Name: SwpMargin.c                                             *
 *  Function: Calculates the swap margin                              *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 05/01/96                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 05/01/96 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "SwpAccess.h"
#include "Swp_h_all.h"

char *SwpMargin(double pv, long start, long nfp_or_end, char *cpdStr,
                char *basisStr, char *ycName, double *margin) {
  Err err;

  /* Computes a margin on the leg of a swap to reach a given PV */
  err =
      swp_f_SwapMargin(pv, start, nfp_or_end, cpdStr, basisStr, ycName, margin);
  if (err)
    return err;

  /* Return a success message */
  return NULL;
}
