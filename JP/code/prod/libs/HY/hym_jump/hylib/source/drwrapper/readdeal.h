/*------------------------------------------------------------------
HEADER FILE:    readdeal.h

CREATED BY:     Neil Yang - Feb 2000
 

PURPOSE:        Data conversion functions to be used in the 
                Kapital DR wrapper in connection with  ALIB
                based executables.
---------------------------------------------------------------------- */
#ifndef _READDEAL_H
#define _READDEAL_H


#include <stdio.h>
#include "bastypes.h"
#include "tcurve.h"
#include "kwrapgen.h"
/*#include "cfileio.h"*/

/*----------------------------------------------------------------------------
FUNCTION:    DrKapWrapGetCreditDefaultSwap

CREATED BY:  Neil Yang -- Feb 2000

DESCRIPTION: Kapital wrapper function to read the deal data for a credit default swap
			 hycds_t.dat file (in DR format) into the necessary market
             parameters.
             
             Restrictions are applied concerning the format of the data
             file.
----------------------------------------------------------------------------*/
GTO_EXPORT(int) DrKapWrapGetCreditDefaultSwap(
    char            *filename,              /* (I) Usually "resetccs_t.dat"     */
    double          *ppy,                   /* (O) period per year */
	double          *beta,                  /* (O) */
	double          *dps,                   /* (O) */
	double          *x,                     /* (O) */
	double          *lim,                   /* (O) */
	double          *lim1,                  /* (O) */
	double          *lim2,                  /* (O) */
	double          *recovery,              /* (O) */
	double          *notional,              /* (O) */
	TDate           *maturityDate,          /* (O) */
	char            *fee);           /* (O) */


                
#endif     
