/****************************************************************
 * Module:	DRL
 * Submodule:	TS - Curve Data Structure
 * File:	drical.h
 * Function:	Curve routines.
 * Author:	Christian Daher
 * Revision:	$Header$
 *****************************************************************/
#ifndef	_drical_H
#define	_drical_H
#include "drlstd.h"

#include <stdio.h>

#include "dritkwrp.h"


/*-------------------------------------------------------------
 *
 */


DLL_EXPORT(int)
DriCalibrateVnfmParams(
	const char *logFnam,		/* (I) log file (or NULL)*/
	const char *logLink,		/* (I) log file link (or NULL)*/
	const char *curCode,		/* (I) currency code (or NULL)*/
	const char *tagName,		/* (I) calibration name (or NULL) */
	FILE *fpIn,			/* (I) input */
	TDrWrapperData *drwData);	/* (I) mkt env */






#endif	/* _drical_H */

