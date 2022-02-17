/*-----------------------------------------------------------------------
C SOURCE FILE:  fwdtospttt.c
  
CREATED BY:     Julia Chislenko Sptember, 1998

PURPOSE:        Test driver for DriFwdToSpotRatesOL

$Header: /home/drfiprod/cvsroot/home/drfidev/driface/test/src/fwdtosptt.c,v 1.1 1999/05/20 18:54:22 dliu Exp $
-------------------------------------------------------------------------*/

#include "cgeneral.h"           /* General stuff */
#include "cerror.h"             /* GtoErrMsg */
#include "date_sup.h"           /* GtoFreq2TDateInterval */
#include "yearfrac.h"           /* GtoStringToDayCountConv */
#include "convert.h"            /* GtoStringToDateInterval */
#include "cerrsup.h"            /* GtoLoggingSet */
#include "tcurve.h"             /* GtoPrintTCurve */
#include "interp.h"             /* GtoStringToInterpType */
#include "zcurveo.h"            /* TZeroCurve */
#include "zr2coup.h"            /* GtoZerosToCouponsPoint */
#include "zr2simp.h"            /* GtoZerosToSimplePoint */
#include "gtomat.h"             /* GtoMatrix2D */
#include "gtozc.h"              /* GtoZCCash, GtoZCSwaps */
#include "mainlot.inc"
#include "testpub.h"

#include "fwdtospt.h"

ARG_DEF g_ArgInfoArray[] =
{
    {"ValueDate:",         AT_DATE},
    {"StubZC:",            AT_OBJECT, "ZC"},
    {"FwdStartDates:",     AT_ARRAY(AT_DATE)},
    {"FwdMatDates:",       AT_ARRAY(AT_DATE)},
    {"FwdPayIntervals:",   AT_ARRAY(AT_OBJECT),"IVL"},
    {"FwdDayCountConvs:",  AT_ARRAY(AT_OBJECT),"DCC"},
    {"FwdSwapRates:",      AT_ARRAY(AT_DOUBLE)},
    {"ReqSpotMatDates:",   AT_ARRAY(AT_DATE)},
    {"ReqSwapFlags:",      AT_ARRAY(AT_LONG)},  
    {"ReqSwapFreq:",       AT_LONG},
    {"ReqCashDCC:",        AT_OBJECT, "DCC"},
    {"ReqSwapDCC:",        AT_OBJECT, "DCC"},
    {"CouponInterp:",      AT_LONG},
    {"ZeroInterp:",        AT_CHAR_BLOCK},
    {"BadDayConv:",        AT_CHAR_BLOCK},
    {"HolidayFile:",       AT_CHAR_BLOCK}
};


long g_NumArgs = sizeof(g_ArgInfoArray) / sizeof(ARG_DEF);


/*
 * callRoutine(): Called by driver.
 */

void callRoutine(FILE *fp)

{
    static char routine[] = "FWD_TO_SPOT_RATES_O";
    int status = FAILURE;

    TCurve        *cashZC = NULL;
    TCurve        *fullZC = NULL;
    double        *prices = NULL;
    double        *spotRates = NULL;
    TMatrix2D     *spotRatesM = NULL;
    TMatrix2D     *tweaks = NULL;
    long          idx;
    long          numRates;
    long          numCashRates;
    long          zeroInterp;
    static long   count = 1;

    TDate         *valueDate = GET_ARRAY_DATA(0, TDate);
    TZeroCurve    *stubZC = GET_SCALAR_PTR(1, TZeroCurve);
    TDate         *fwdStartDates = GET_ARRAY_DATA(2, TDate);
    TDate         *fwdMatDates = GET_ARRAY_DATA(3, TDate);
    TDateInterval **fwdPayIntervals = GET_ARRAY_PTR(4, TDateInterval);
    TDayCountConv **fwdDCCs = GET_ARRAY_PTR(5, TDayCountConv);
    double        *fwdSwapRates = GET_ARRAY_DATA(6, double);
    TDate         *reqSpotMatDates = GET_ARRAY_DATA(7, TDate);
    long          *reqSwapFlags = GET_ARRAY_DATA(8, long);
    long          *reqSwapFreq = GET_ARRAY_DATA(9,long);
    TDayCountConv *reqCashDCC = GET_SCALAR_PTR(10, TDayCountConv);
    TDayCountConv *reqSwapDCC = GET_SCALAR_PTR(11, TDayCountConv);
    long          *couponInterp = GET_ARRAY_DATA(12,long);
    char          *zeroInterpStr = GET_ARRAY_DATA(13, char);
    char          *badDayConv = GET_ARRAY_DATA(14, char);
    char          *holidayFile = GET_ARRAY_DATA(15, char);


    printf ("%ld: Calling DriFwdToSpotRatesOL...\n", count++);


    
#ifdef LOGGING
    GtoLoggingSet(1);
#endif


    numRates = fwdStartDates[0];

    /* Call the routine
     */
    if (DriFwdToSpotRatesOL
        (valueDate, stubZC, 
         fwdStartDates, fwdMatDates, fwdPayIntervals, fwdDCCs,
         fwdSwapRates, 
         reqSpotMatDates, reqSwapFlags, 
         reqSwapFreq, reqCashDCC, reqSwapDCC, 
         couponInterp, zeroInterpStr,
         badDayConv, holidayFile,
         &spotRatesM,
	 &tweaks) IS FAILURE)
    {
        goto done;
    }

    /* Generate zero curve and check fwd rates
     */
    prices = NEW_ARRAY(double, numRates);
    spotRates = NEW_ARRAY(double, numRates);
    if (prices IS NULL ||
	spotRates == NULL)
        goto done;

    for (idx=0; idx < numRates; idx++)
    {
	if(GtoMatrixGetValue(spotRatesM,
			     (int)idx,
			     0,
			     spotRates+idx) == FAILURE)
	    goto done;
        prices[idx] = 1.;
	
    }


    if (GtoStringToInterpType(zeroInterpStr+1, &zeroInterp) IS FAILURE)
        goto done;              /* Failed */

     /* Count number of swaps
     */
    for (idx=1; idx<=numRates; idx++)
    {
	if(reqSwapFlags[idx] == 1)
	    break;
    }
    
    numCashRates = idx-1;
   
    if((cashZC = GtoZCCash(stubZC->curve, 
			reqSpotMatDates+1, 
			spotRates, 
			numCashRates,
			*reqCashDCC)) == NULL)
        goto done;   
   
    if((fullZC = GtoZCSwaps (cashZC, 
		     NULL,
		     reqSpotMatDates+1+numCashRates, 
		     spotRates+numCashRates, 
		     prices, 
		     numRates-numCashRates,
		     reqSwapFreq[1], 
		     reqSwapFreq[1], /* FloatSide=N/A */
		     *reqSwapDCC,      
		     *reqSwapDCC,  /* FloatSide=N/A */
		     couponInterp[1], 
		     zeroInterp,
		     0/* FloatSide=N/A */, FALSE,/* Dont value floating side */
		     badDayConv[1], 
		     holidayFile+1))
       == NULL) goto done;


    /* Print out the results 
     */
    fprintf(fp, "\nOutputs:\n");
    fprintf(fp, "Date        Spot          OrigFwd     ImpliedFwd    Bp diff\n");
    for (idx=1; idx <= numRates; idx++)
    {
        double newRate;

        if (GtoZerosToCouponsPoint/*Adj*/
            (fullZC,
	     zeroInterp, 
             fwdStartDates[idx], fwdPayIntervals[idx], 
             fwdMatDates[idx], *(fwdDCCs[idx]),  GTO_STUB_SIMPLE,
             FALSE, /*badDayConv[1], badDayConv[1], holidayFile+1,*/ 
             &newRate) IS FAILURE)
            goto done;

        fprintf(fp, "%10.10s  %3.9f%c  %3.9f %3.9f %2.9f\n",
                GtoFormatDate(reqSpotMatDates[idx]), 
                spotRates[idx-1], 
                reqSwapFlags[idx] ? 'S' : 'M',
                fwdSwapRates[idx], newRate,(newRate-fwdSwapRates[idx])*10000.);


    }
    GtoErrMsg("\n\n");
    GtoMatrixPrint(tweaks, "dSpot/dFwd");

    status = SUCCESS;

  done:
    GtoFreeTCurve(cashZC);
    GtoFreeTCurve(fullZC);
    FREE_ARRAY(prices);
    FREE_ARRAY(spotRates);
    GtoMatrixFree(spotRatesM);
    GtoMatrixFree(tweaks);

    if (status IS  FAILURE)
        GtoErrMsg("%s: Failed.\n", routine);

    return;
}



int testInit (TTestInitData *initData)
{
    return (GtoRegisterClasses(initData->om));
}
