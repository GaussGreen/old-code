/****************************************************************************/
/*      Declaration file for structure of data.                             */
/****************************************************************************/
/*      drieq.h                                                             */
/****************************************************************************/
#ifndef DRIEQ_H
#define DRIEQ_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"
#include "tcurve.h"
#include "ldate.h"
#include "cerror.h"
#include "eqsettle.h"
#include "eqdiv.h"
#include "gtobf.h"

#include "dritkwrp.h"           /* TDrWrapperData routines */

#define         MAXNBDATE       500          /* Maximum number of elements in input date array     */



/*t-@CDOC(idxn=TEqStatData)-------------------------------------------------
 * Equaty static file data structure.
 *
 * \index{TEqStatData}
 */

typedef struct _TEqStatData{

	TDividendList     *divList;
	char              settleType; 
	TEquitySettlement *stm;

}TEqStatData;

/*e*/	


/*t-@CDOC(idxn=TEqStmPrivate)-------------------------------------------------
 * Equaty settlement data structure.
 * Copied from ALIB (eqsettle.c)
 * Structure with easier working
 * data for determination of settlement dates
 * for stock.
 *
 * \index{TEqStmPrivate}
 */

typedef struct _TEqStmPrivate {
 
    /*
     * For rolling settlement, ie T+n days
     * (US system, S&P500)
     *
     * if tradeDate > stmPeriodDates[i], use period[i]
     * if tradeDate > stmPeriodDates[last], use period[last+1]
     *
     * 'settlementHolName' is used to decide settlement dates with T+n
     */
    TDateList   *stmPeriodDates;  /* Date when stmPeriod switches */
    long        *stmPeriods;      /* settle periods for stmPeriodDates
                                   * Must be 1 more stmPeriods than
                                   * stmPeriodDates
                                   */
    long         rollingIndex;    /* Index into stmPeriodDates[] and
                                   * stmPeriods[] */
 
    char        *settlementHolName;
 
    /*
     * Allows more efficient searching
     * GTO_EQ_STM_DIRN_INC    Increasing dates
     * GTO_EQ_STM_DIRN_NONE   No pre-arranged direction
     * GTO_EQ_STM_DIRN_DEC    Decreasing dates
     *
     */
    long         direction;
 
 
    /*
     * For fixed settlement dates as published
     * by an exchange.
     * (French system, CAC40)
     * Before first ltd use rolling stm at stmPeriodBefore
     * After last ltd use rolling stm at stmPeriodAfter
     */
    TDateList   *lastTradingDates;
    TDateList   *settlementDates;
    long         settleIndex;
    long         stmPeriodBefore;
    long         stmPeriodAfter;
 
 
} TEqStmPrivate ;

/*e*/	


GTO_EXPORT(int)
DriTEqStatGet(char 	   *pathdir,      /* (I) tmp direcory name (or NULL) */
	      char 	   *eqStaFnam,    /* (I) file name (equity.sta)  */
	      long         busDayConv,    /* (I) Business Day Conv    */
	      char         *holidayFile,  /* (I) Holiday File       */
	      TEqStatData  **eqStatData); /* (O) equaty static data */


GTO_EXPORT(extern TEqStatData *)
DrNewTEqStatData(TDividendList     *divList,
		 char		   settleType,
                 TEquitySettlement *stm);


GTO_EXPORT(extern void)
DrFreeTEqStatData(TEqStatData *eqStatData);


/*---------------------------------------------------------------------
 * Given a swap zero curve, build a funding curve for a given index.
 * - read funding spreads from equity DR wrapper environment
 * - recover MM and par swap rates from swap curve
 * - add spreads to swap rates, and generate the corresponding zero curve
 *
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriTFundingCurveGet(   
char *pathdir,			/* (I) tmp direcory name (or NULL) */
TDrWrapperData *drWrapper,	/* (I) dr wrapper data */
TCurve        **indxZC);        /* (O) index zero curve with funding spread */


/*---------------------------------------------------------------------
 * Read data from equity.dyn in London format
 * Basis of vol curve is set to 4L (quarterly)
 *
 * Returns SUCCESS/FAILURE.
 */

DLL_EXPORT(int)
DriEqDynDataGet(   
char    *pathdir,	       	/* (I) tmp direcory name (or NULL) */
char    *eqDynFnam,    		/* (I) file name (equity.dyn)  */
double  *spotValue,            	/* (O)  */
double  *correlation,          	/* (O)  */
TCurve **indxVolCurve);        	/* (O) index volatility curve */


#endif  /* DRIEQ_H */
