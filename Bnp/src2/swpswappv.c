/**********************************************************************
 *      Name: SwpSwapPV.c                                             * 
 *  Function: Calculates the PV of a swap                             *
 * Copyright: (C) Paribas Capital Markets Ltd.                        *
 *--------------------------------------------------------------------*
 *    Author: Finbarr O'Sullivan                                      *
 *      Date: 12/12/95                                                *
 *--------------------------------------------------------------------*
 *    Inputs:                                                         *
 *   Returns:                                                         *
 *   Globals:                                                         *
 *--------------------------------------------------------------------*
 * Modification Record                                                *
 * Date     Inits   Comments                                          *
 * dd/mm/yy                                                           *
 * 12/12/95 FOS     Created for SORT5-GRFN3 port to NT                *
 **********************************************************************/
#include "Swp_h_all.h"
#include "SwpAccess.h"


char *SwpSwapPV(long start, long  nfp_or_end,
                char *cpdStr, char *basisStr, double strike, 
                char *dcntName, char *refRateCode, double *pv)
{
Err 	err;

/* Computes the PV of a Swap, taking spreads into account */
	err = swp_f_SwapPv(
		strike,
		start,
		nfp_or_end,
		cpdStr,
		basisStr,
		dcntName,
		refRateCode,
		pv);
	if (err)
		return err;

/* Return a success message */
	return NULL;
}	
