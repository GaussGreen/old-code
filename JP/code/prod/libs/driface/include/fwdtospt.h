/*-------------------------------------------------------------------------

C SOURCE FILE:   fwdtospt.h

CREATED BY:      Julia Chislenko September 1998

CONTAINS:        Object wrapper DRIFwdToSpotRatesOL

$Header$
----------------------------------------------------------------------------*/
#ifndef DRI_FWD_TO_SPT_H
#define DRI_FWD_TO_SPT_H

#include <ctype.h>                 /* toupper */
#include "bastypes.h"
#include "cgeneral.h"
#include "ldate.h"                 /* GTO_ACT_365F */
#include "macros.h"
#include "tcurve.h"                /* GtoNewTCurve */
#include "zcurveo.h"               /* TZeroCurve */
#include "dcconvo.h"               /* TDayCountConv */
#include "dtivlo.h"                /* TDateInterval */

/*----------------------------------------------------------------------
  FUNCTION:       DriFwdToSpotRatesOL
  
  CREATED BY:     Julia Chislenko  September, 1998
  
  PURPOSE:        Given n fwd rates compute n spot rates that generate
                  the zero curve producing original fwd rates. The 
		  interface except for "solveOrder" array is the same
                  as in Doug's GtoZeroCurveFwdToSpotRates, but instead
                  of using one-dimensional solver n times with predetermined 
		  order of matching fwds, this routine makes one call
		  to the n-dim solver (using IMSL).
 */

GTO_EXPORT(int)  DriFwdToSpotRatesOL
  (TDate           *valueDate,       /* (I) Used only if stubZC is empty */
   TZeroCurve      *stubCurve,       /* (I) Stub zero curve */
   TDate           *fwdStartDates,   /* (I) Fwd swap start dates */
   TDate           *fwdMatDates,     /* (I) Fwd swap maturity dates 
				      *      Have to ascend */
   TDateInterval  **fwdPayIntervals, /* (I) Fwd frequencies */
   TDayCountConv  **fwdDCCs,         /* (I) Fwd day count convs */
   double          *fwdSwapRates,    /* (I) Fwd starting swap rates */
   TDate           *reqSpotMatDates, /* (I) Requested spot maturity dates
				      *      Have to ascend */
   long            *reqSwapFlags,    /* (I) 0=CashRate, 1=SwapRate
				       *      Swap rates should follow cash */
   long            *reqSwapFreq,     /* (I) Requested output swap freq */
   TDayCountConv   *reqCashDCC,      /* (I) Req. money market day count conv*/
   TDayCountConv   *reqSwapFixedDCC, /* (I) Req. swap fixed day count conv */
   long            *couponInterpType,/* (I) Zero curve coupon interp type */
   char            *zeroInterpStr,   /* (I) Zero curve zero interp type */
   char            *badDayConv,      /* (I) Bad day convention */
   char            *holidayFile,     /* (I) Holidays for bad day adjustment */
   TMatrix2D      **spotRates,       /* (O) Vector of spot cash/swap rates */
   TMatrix2D      **tweakMatrix);    /* (O) Matrix dSpot/dFwd */

#endif
