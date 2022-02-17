/************************************************************************
 * Module:      driface
 * File:        cpis.h
 * Function:    CPI linked security
 * Author: vb

$ Header:$
 ************************************************************************/
#ifndef DRI_CPIS
#define DRI_CPIS

#include "drlstd.h"
#include "bastypes.h"

typedef struct{
	/* market data */
	TDate       today;              /* (I) spot date */
	TCurve     *nomZC;              /* (I) nominal zero curve */
	TCurve     *realZC;             /* (I) real zero curve */
	TCurve     *discZC;             /* (I) discount zero curve */
	double		spotCPI;			/* (I) CPI at the spot date */
	/* deal data */
	TDate		maturDate;	/* (I) date of the deal maturity */
	double      origNotl;           /* (I) original notional */
	double      firstCPI;             /* (I) theCPI at deal intiation*/
	double		leverage;		/* (I) leverage factor */

	/* option data */
	double		fee;				/* (I) redemption fee added to tips yield */
	double		finSpread;		/* Morgan credit spread for financing */
	long          numExDates;        /* (I) number of exercise dates */
	TDate      *exDates;      /* (I) [numExDates] exersice dates */
	TDate      *ntDates;        /* (I) [numExDates] notify dates */
	double		trsCoupon;		/* coupon of pegged treasury */
	long		trsNumDates; /* number of treasury dates */
	TDate		*trsDates;		/* pegged date teasury */
	double		tipsCoupon;		/* coupon of pegged tips */
	long		tipsNumDates; /* number of tips dates */
	TDate		*tipsDates;		/* pegged date tips */
	/* vol market data */
	double      cpiBpVol;           /* (I) nominal rate bp volatility */
	double      cpiMR;              /* (I) nominal (normal) mean reversion */

	double     priceSwap; 		/* (O) underlying price */
	double     priceOption; 	/* (O) Option price */
	double     price; 				/* (O) total price */
} TCPISecurity;

#ifdef __cplusplus
extern "C"
{
#endif

/*f---------------------------------------------------------------------
 * Pricing routine for CPI linked swaps.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCPISecurity(   
TCPISecurity	*cpis
);

/*f---------------------------------------------------------------------
 * Dr-Wrapper for {\tt DriCPISwap}.
 * The argument {\tt dataFnam} specifies the name of the file
 * containing the data (if it is NULL, the default name
 * "cpiswap_w.dat" is used.
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriCpiSecurityW(char *dataFnam);

#ifdef __cplusplus
}
#endif

#endif
