/************************************************************************
 * Module:      driface
 * File:        cpiswap.h
 * Function:    CPI linked swap pricer
 * Author:      Julia Chislenko Sep 1999

$ Header:$
 ************************************************************************/
#ifndef DRI_CPI_SWAP
#define DRI_CPI_SWAP

#include "drlstd.h"

#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"             /* MAX */

/*f---------------------------------------------------------------------
 * Pricing routine for CPI linked swaps.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCPISwap(   
TDate       today,              /* (I) spot date */
char        instrType,   	/* (I) (T)IP/(C)urrentPay/(A)ccretingCurrentPay */
double      origNotl,           /* (I) original notional */
double      coupon,             /* (I) coupon */
TDayCount   dayCountConv,       /* (I) for coupon payments */
double     *cpiFixings,         /* (I) [CPI_NUM_FIXINGS] Spot, last, first(for TIPs only) */
TCurve     *nomZC,              /* (I) nominal zero curve */
TCurve     *realZC,             /* (I) real zero curve */
TCurve     *discZC,             /* (I) discount zero curve */
long        numPayDates,        /* (I) number of payment dates */
TDate      *accStartDates,      /* (I) [numPayDates] accrue start dates */
TDate      *accEndDates,        /* (I) [numPayDates] accrue end dates */
TDate      *payDates,           /* (I) [numPayDates] payment dates */
double      nomBpVol,           /* (I) nominal rate bp volatility */
double      realBpVol,          /* (I) real rate bp volatility */
double      nomMR,              /* (I) nominal (normal) mean reversion */
double      realMR,             /* (I) real (normal) mean reversion */
double      correlation,        /* (I) between nom and real interest rates */    
double     *pv);		/* (O) present value */

/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriCPISwap}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "cpiswap_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCpiSwapW(char *dataFnam);

#endif
