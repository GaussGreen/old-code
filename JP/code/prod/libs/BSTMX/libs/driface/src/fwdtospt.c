/*-------------------------------------------------------------------------

C SOURCE FILE:   fwdtospt.c

CREATED BY:      Julia Chislenko September 1998

CONTAINS:        Object wrapper for DRIFwdToSpotRates

$Header$
----------------------------------------------------------------------------*/
#include "drlmatrixo.h"               /* MatrixNew */

#include "fwdtospt.h"               /* Prototype consistency */

#include <ctype.h>                 /* toupper */
#include "bastypes.h"
#include "cgeneral.h"
#include "ldate.h"                 /* GTO_ACT_365F */
#include "date_sup.h"              /* GtoCountDatesArrayWrap */
#include "cerror.h"
#include "convert.h"
#include "date_sup.h"              
#include "interp.h"                /* GTO_LINEAR_INTERP */
#include "macros.h"
#include "tcurve.h"                /* GtoNewTCurve */
#include "zcurveo.h"               /* TZeroCurve */
#include "zr2coup.h"               /* GtoZerosToCouponsPoint */
#include "zr2simp.h"               /* GtoZerosToSimplePoint */
#include "gtomat.h"                /* Matrix2D */
#include "gtozc.h"                 /* GtoZCCash, GtoZCSwaps */
#include "zcswpadj.h"              /* GtoZeroCurveParSwapsAdjusted */
#include "duration.h"              /* MMmodDuration & BondModDuration */

#include "imsl.h"
#include "alimsl.h"

#define DRI_MAX_FWD_ERR    1e-5
#define DRI_TWEAK_SIZE     1e-5


#define ZC_DEFAULT_BASIS 1.
#define ZC_DEFAULT_DAYCNT GTO_ACT_365F

#define xDRI_ZC_DEBUG
#define xDRI_MATRIX_DEBUG
#define xDRI_DEBUG_STEP

/* ZC generation parameters 
 * used by the solver
 */
typedef struct
{
    TCurve        *stubZC;          /* Stub zero curve */

    TDate         *spotDates;       /* [numRates] */
    int            numCashRates;
    long           cashDCC;

    TDate         *spotSwapDates;   /* Spot starting swap dates */
    double        *prices;          /* all 1. for par */
    long           numSwaps;        /* Num spot starting swaps */
    int            swapFreq;
    long           swapDCC;
    long           couponInterpType;
    long           zeroInterpType;  /* Zero curve zero interp type */

    long           badDayConv;      /* Bad day convention */
    char          *holidayFile;     /* Holiday file */

    TDate           *fwdStartDates;   /* (I) Fwd swap start dates */
    TDate           *fwdMatDates;     /* (I) Fwd swap maturity dates */
    TDateInterval  **fwdPayIntervals; /* (I) Fwd frequencies */
    long            *fwdDCCArray;     /* (I) Fwd day count convs */
    double          *fwdSwapRates;    /* (I) Fwd starting swap rates */

    double          *actFwdRates;     /* (O) Solver's solution */
    
    int             status;
    

} TSwapZcParams;

/* Global parameter structure - to be used by the solver
 */
TSwapZcParams p;

/* Match forward rates  - function called by the solver
 * and by tweak function
 */
void matchFwdRates (
int              n,                    /* (I) */
double           *spotRates,           /* (I) */ 
double           *fwdErrors);          /* (O) */ 


/* Allocate n x n and compute tweak matrix - 
 * first tweak spot rates to get dFwd/dSpot,
 * then invert the matrix using IMSL
 */
static double *driDSpotOverDFwdMatrix(long    numRates, 
				     double *spotRates);

static void driLogInputs
  (TCurve          *stubZC,          /* (I) Stub zero curve */

   TDate           *fwdStartDates,   /* (I) Fwd swap start dates(incr. order)*/
   TDate           *fwdMatDates,     /* (I) Fwd swap maturity dates */
   TDateInterval  **fwdPayIntervals, /* (I) Fwd frequencies */
   long            *fwdDCCs,         /* (I) Fwd day count convs */
   double          *fwdSwapRates,    /* (I) Fwd starting swap rates */

   TDate           *reqSpotMatDates, /* (I) Req spot maturity dates */
   long            *reqSwapFlags,    /* (I) 1=SwapRate, 0=CashRate */
   long             numRates,        /* (I) Len of all arrays */

   long             reqSwapFreq,     /* (I) Requested swap frequency */
   long             reqCashDCC,      /* (I) Requested cash day count conv */
   long             reqSwapDCC,      /* (I) Requested swap day count conv */
   long             couponInterpType,/* (I) Zero curve coupon interp type */
   long             zeroInterpType,  /* (I) Zero curve zero interp type */
   char             badDayConv,      /* (I) Bad day convention */
   char            *holidayFile);    /* (I) Holidays for bad day adjustment */

/**** Temporarily moved here until the libs move to
 *    ALIB 9.3
 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */ 
static char * tempGtoFormatInterpType();

/*f----------------------------------------------------------------------
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
   TDayCountConv   *reqSwapDCC,      /* (I) Req. swap fixed day count conv */
   long            *couponInterpType,/* (I) Zero curve coupon interp type */
   char            *zeroInterpStr,   /* (I) Zero curve zero interp type */
   char            *badDayConv,      /* (I) Bad day convention */
   char            *holidayFile,     /* (I) Holidays for bad day adjustment */
   TMatrix2D      **spotRates,       /* (O) Vector of spot cash/swap rates */
   TMatrix2D      **tweakMatrix)     /* (O) Matrix dSpot/dFwd */
{
    static char      routine[]="DriFwdToSpotRatesOL";
    int              status=FAILURE;  

    long      *fwdDCCArray=NULL;
    long       idx;
    long       zeroInterpType;
    long       numRates = 0;            /* Total number rates */
    long       numSwapRates;            /* Swap rates returned */
    long       numCashRates;            /* Cash rates returned */
    long       numFwdStartDates = GtoCountDatesArrayWrap(fwdStartDates);
    long       numFwdMatDates = GtoCountDatesArrayWrap(fwdMatDates);
    long       numReqSpotMatDates = GtoCountDatesArrayWrap(reqSpotMatDates);

    double *solveRates = NULL;
    double *tweaks = NULL;
    double *prices = NULL;
    long   *flags = NULL; 
    double *actFwdRates = NULL;

    TCurve *stubZC = NULL;
    TCurve *guessZC = NULL;
    TCurve *guessZCcash = NULL;
    TInterpInfo *interpInfo = NULL;
    double *guessRates= NULL;
    TDateInterval intvl;
    TDate prevDate;
    TBoolean   newStubZC = FALSE; 


#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Checking inputs\n");
#endif

    /* Initialize stubZC
     */
    if (stubCurve IS NULL) 
    {
        stubZC = GtoNewTCurve (valueDate[1],
                               0,                 /* # points */
                               ZC_DEFAULT_BASIS,
                               ZC_DEFAULT_DAYCNT);
        if (stubZC == NULL)
            goto done;        

        newStubZC = TRUE; 
    }
    else 
    {
        stubZC = stubCurve->curve; 
    }
    
   /* Check for array lengths
     */
    if (numFwdStartDates ISNT numFwdMatDates ||
        numFwdStartDates ISNT numReqSpotMatDates)
    {
        GtoErrMsg("%s: # fwdStartDates (%ld), # fwdMatDates (%ld)\n"
                  "\t and #reqSpotMatDates (%ld) not equal.\n",
                  routine, numFwdStartDates, numFwdMatDates,
                  numReqSpotMatDates);
        goto done;
    }

    numRates = numFwdStartDates;

#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Converting wrapper arguments\n");
#endif

    /* Convert day count convs from pointers to longs
     */
    fwdDCCArray = NEW_ARRAY(long, numRates);
    if (fwdDCCArray IS NULL)
        goto done;
    for (idx=0; idx < numRates; idx++)
        fwdDCCArray[idx] = *(fwdDCCs[idx+1]);
    
    /* Convert zero interp type.
     */
    if (GtoStringToInterpType (zeroInterpStr+1, 
                               &zeroInterpType) IS FAILURE)
    {
        GtoErrMsg("%s: Failed converting zero interp type.\n", routine);
        goto done;
    }

    GTO_IF_LOGGING
        (driLogInputs
         (stubZC, fwdStartDates+1, fwdMatDates+1, 
          fwdPayIntervals+1, fwdDCCArray, fwdSwapRates+1, 
          reqSpotMatDates+1, reqSwapFlags+1, numRates,
          reqSwapFreq[1], *reqCashDCC, *reqSwapDCC, 
          couponInterpType[1], zeroInterpType,
          badDayConv[1], holidayFile+1))
         
	numRates = numFwdStartDates;

    /* Count number of swaps
     */
    numCashRates = numRates;
    prevDate = stubZC->fBaseDate;
    
    for (idx=1; idx<=numRates; idx++)
    {
	long type = reqSwapFlags[idx];
	TBoolean isSwap = (type==1);
	
	if(!isSwap && type != 0)
	{
	    GtoErrMsg("%s: Unknown rate type %d", routine, type);
	    goto done;
	}
	if(reqSpotMatDates[idx]<=prevDate)
	{
	    GtoErrMsg("%s: Requested maturity[%d] (%s) <= maturity[%d] (%s)\n",
		      routine, idx, GtoFormatDate(reqSpotMatDates[idx]),
		      idx-1, GtoFormatDate(prevDate));
	    goto done;
	}
	if(numCashRates<numRates && !isSwap)
	{
	    GtoErrMsg("%s: All swap rates must follow cash rates, and "
		      "maturities must ascend.\n",
		      routine, idx, GtoFormatDate(reqSpotMatDates[idx]),
		      idx-1, GtoFormatDate(prevDate));
	    goto done;
	}
	if(numCashRates==numRates && isSwap)
	{
	    numCashRates = idx-1;
	}
	prevDate = reqSpotMatDates[idx];
    }

    numSwapRates = numRates - numCashRates;

    /* Prices - fill in with 1. since only dealing with par rates
     */
    if(numSwapRates > 0)
    {
	prices = NEW_ARRAY(double, numSwapRates);
	if (prices == NULL)
	    goto done;
    }
    for (idx=0; idx < numSwapRates; idx++)
        prices[idx] = 1.;
    
    /* Make sure that stub zero curve and requested 
     * maturities don't overlap.
     */
    if (stubZC->fNumItems > 0 &&
	reqSpotMatDates[idx] <= stubZC->fArray[stubZC->fNumItems-1].fDate)
    {
	GtoErrMsg("%s: Last stub zero curve date (%s) comes on or before\n"
		  "\trequested spot maturity date(%s).\n",
		  routine, 
		  GtoFormatDate(stubZC->fArray[stubZC->fNumItems-1].fDate),
		  GtoFormatDate(reqSpotMatDates[idx]));
	status = FAILURE;
    }
    actFwdRates = NEW_ARRAY(double, numRates);
    if(actFwdRates == NULL)
	goto done;

#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Filling in solver params\n");
#endif

    /* Fill in solver parameters
     */
    p.stubZC = stubZC; 

    p.spotDates = reqSpotMatDates+1;
    p.numCashRates = numCashRates;
    p.cashDCC = *reqCashDCC;
    p.spotSwapDates = reqSpotMatDates+1+numCashRates;
    p.prices = prices;
    p.numSwaps = numSwapRates;
    p.swapFreq = reqSwapFreq[1];
    p.swapDCC = *reqSwapDCC;
    p.couponInterpType = couponInterpType[1];
    p.zeroInterpType = zeroInterpType;
    p.badDayConv = (long)(toupper(badDayConv[1]));
    p.holidayFile = holidayFile+1;
    p.fwdStartDates = fwdStartDates+1;
    p.fwdMatDates = fwdMatDates+1;
    p.fwdPayIntervals = fwdPayIntervals+1;
    p.fwdDCCArray = fwdDCCArray;
    p.fwdSwapRates = fwdSwapRates+1;   /* target rates */
    p.status = SUCCESS;
    p.actFwdRates = actFwdRates;
    
#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Computing init guess\n");
#endif

    /* Compute initial guess by generating boot-strapped zero curve
     * and then recovering swap rates
     */
    guessRates = NEW_ARRAY(double, numRates);
    flags = NEW_ARRAY(long, numRates);

    if (guessRates == NULL ||
	flags == NULL )
	goto done;
    for (idx=0; idx < numRates; idx++)
        flags[idx] = 1;
    
    interpInfo = GtoInterpInfoNew(p.zeroInterpType,
				  NULL, NULL);
    if(interpInfo == NULL)
	goto done;

    if(p.stubZC == NULL)
    {
	prevDate = p.stubZC->fBaseDate + 1;
	guessZCcash = GtoZCCash(NULL, 
				&prevDate, 
				p.fwdSwapRates, 
				1, /* numCashRates */
				p.cashDCC);
	if(guessZCcash == NULL)
	    goto done;
    }
    if((guessZC = GtoZeroCurveParSwapsAdjusted 
	(p.stubZC == NULL? guessZCcash: p.stubZC,
	 numRates,
	 p.fwdStartDates,
	 p.fwdMatDates,
	 p.fwdSwapRates,
	 flags,
	 NULL, /* adjustments */
	 fwdDCCArray[0],
	 p.fwdPayIntervals[0],
	 GTO_STUB_SIMPLE,
	 FALSE, /* stubAtFront */
	 p.badDayConv,
	 p.badDayConv,
	 p.holidayFile,
	 interpInfo)) == NULL)
	goto done;;
					    
#ifdef DRI_ZC_DEBUG
   GtoPrintTCurve(guessZC,"Guess Zero Curve");
#endif

#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Computing init rates\n");
#endif

    for (idx=0; idx < numCashRates; idx++)
    {
	if(GtoZerosToSimplePoint(guessZC,
				 p.zeroInterpType,
				 guessZC->fBaseDate,
				 p.spotDates[idx],
				 p.cashDCC,
				 &guessRates[idx]) == FAILURE)
	    goto done;
    }    
    if(GtoFreq2TDateInterval((long)p.swapFreq, &intvl) == FAILURE)
	goto done;

    for (idx=0; idx < p.numSwaps; idx++)
    {
	if(GtoZerosToCouponsPoint/*Adj*/(guessZC,
				     p.zeroInterpType,
				     guessZC->fBaseDate,
				     &intvl,
				     p.spotSwapDates[idx],
				     p.swapDCC,
				     GTO_STUB_SIMPLE,
				     FALSE, /*StubAtFront*/
/*p.badDayConv,p.badDayConv,p.holidayFile,*/
				     &guessRates[numCashRates+idx]) == FAILURE)
	    goto done;
    }
    
#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Solving with IMSL\n");
#endif

    /* Convert fwd rates to spot rates using IMSL 
     * multidimensional solver
     */
    GtoImslCmathInit ();
    solveRates = imsl_d_zeros_sys_eqn(matchFwdRates, 
				      (int)numRates, 
				      IMSL_XGUESS, guessRates,
				      IMSL_ERR_REL, 1.0e-8, /* sqrtEps */
				      IMSL_MAX_ITN, 100,    /* maxIter */
				      0);
    
    for (idx=0; idx < numRates; idx++)
    {
	if(solveRates != NULL && solveRates[idx]<=0.)
	{
	    GtoErrMsg("Negative spot rate[%d] %f%%\n", 
		      idx, solveRates[idx]*100.);
	    p.status=FAILURE;
	}
	if(ABS(p.actFwdRates[idx]-p.fwdSwapRates[idx])>DRI_MAX_FWD_ERR)
	{
	    GtoErrMsg("Error %fbp for fwdSwapRate[%d] \n", 
		      (p.actFwdRates[idx]-p.fwdSwapRates[idx])*10000.,
		      idx);
	    p.status=FAILURE;
	}
    }
    if(solveRates == NULL  ||
       /*       GtoImslCmathErrorStatus () == FAILURE || */
       p.status == FAILURE)
    {
	GtoErrMsg("%s: Can not solve for spot rates.\n", routine);
	goto done;
    }

    
#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Computing tweak matrix\n");
#endif

    /* Compute tweak matrix 
     */
    if((tweaks = driDSpotOverDFwdMatrix(numRates, solveRates)) == NULL)
	goto done;

    /* Fill in outputs
     */
    if((*spotRates = MatrixNew(numRates, 1, solveRates)) == NULL)
    {
	GtoErrMsg ("%s: Failed to construct spotRates object.\n", routine);
	goto done;
    }
    if((*tweakMatrix = MatrixNew(numRates, numRates, tweaks)) == NULL)
    {
	GtoErrMsg ("%s: Failed to construct tweakMatrix object.\n", routine);
	goto done;
    }
    
    status = SUCCESS;

  done: 
    
#ifdef DRI_DEBUG_STEP
    GtoErrMsg("Freeing memory\n");
#endif

    if(newStubZC)
    {
	GtoFreeTCurve(stubZC);
    }	
    FREE_ARRAY(solveRates);
    FREE_ARRAY(tweaks);
    FREE_ARRAY(actFwdRates);
    FREE_ARRAY(prices);
    FREE_ARRAY(flags);
    FREE_ARRAY(guessRates);
    FREE_ARRAY(fwdDCCArray);
    GtoFreeTCurve(guessZCcash);
    GtoFreeTCurve(guessZC);
    GtoInterpInfoDelete(interpInfo);
    if (status == FAILURE) 
        GtoErrMsg ("%s: Failed.\n", routine); 
    
    return (status); 
}



/* Match forward rates  - function called by the solver
   and tweaking routine
 */
void matchFwdRates (
int              n,                    /* (I) */
double           *spotRates,           /* (I) */ 
double           *fwdErrors)           /* (O) */ 
{
    int status = FAILURE;
    TCurve *cashZC = NULL;
    TCurve *zc = NULL;
    int i;

    double *currRate;

#ifdef DRI_ZC_DEBUG
   static int solveRun = 0;
#endif

    /* Generate ZC using spot rates - have to use ZCSwaps 
     * to handle holidays
     */
    cashZC = GtoZCCash(p.stubZC, 
		       p.spotDates, 
		       spotRates, 
		       p.numCashRates,
		       p.cashDCC);
    if (cashZC IS NULL)
        goto done;              /* Failed */
   
    zc = GtoZCSwaps (cashZC, 
		     NULL,
		     p.spotSwapDates, 
		     spotRates+p.numCashRates, 
		     p.prices, 
		     p.numSwaps,
		     p.swapFreq, 
		     p.swapFreq, /* FloatSide=N/A */
		     p.swapDCC,      
		     p.swapDCC,  /* FloatSide=N/A */
		     p.couponInterpType, 
		     p.zeroInterpType,
		     0/* FloatSide=N/A */, FALSE,/* Dont value floating side */
		     p.badDayConv, 
		     p.holidayFile);
    if (zc IS NULL)
        goto done; 
					    
#ifdef DRI_ZC_DEBUG
    solveRun++;
    GtoErrMsg("\n\nSolver run %d\n\n", solveRun);
    GtoPrintTCurve(zc,"Solve Zero Curve");
#endif

    /* Compute implied fwd rates and the errors
     */
    for (i=0, currRate=p.actFwdRates ; i<n; i++, currRate++)
    {
	if(GtoZerosToCouponsPoint/*Adj*/(zc,
				     p.zeroInterpType,
				     p.fwdStartDates[i],
				     p.fwdPayIntervals[i],
				     p.fwdMatDates[i],
				     p.fwdDCCArray[i],
				     GTO_STUB_SIMPLE,
				     FALSE, /*StubAtFront*/
/*p.badDayConv,p.badDayConv,p.holidayFile,*/
				     currRate) == FAILURE)
	    goto done;

	fwdErrors[i] = *currRate-p.fwdSwapRates[i];
					    
#ifdef DRI_ZC_DEBUG
    GtoErrMsg("\nNew fwd[%d] %f error %f\n", i, temp, fwdErrors[i]);
#endif
	/* penalize negative rates */
	if(*currRate<0.)
	{
	    fwdErrors[i] = 1.; 
	}
    }

    status = SUCCESS;

  done:  

    GtoFreeTCurve(cashZC);
    GtoFreeTCurve(zc);

    if(status == FAILURE)
    {
	for (i=0; i<n; i++)
	{
	    fwdErrors[i] = 0.;
	}
    }

    p.status = status;
    return;
}


static void driLogInputs
  (TCurve          *stubZC,          /* (I) Stub zero curve */

   TDate           *fwdStartDates,   /* (I) Fwd swap start dates(incr. order)*/
   TDate           *fwdMatDates,     /* (I) Fwd swap maturity dates */
   TDateInterval  **fwdPayIntervals, /* (I) Fwd frequencies */
   long            *fwdDCCs,         /* (I) Fwd day count convs */
   double          *fwdSwapRates,    /* (I) Fwd starting swap rates */

   TDate           *reqSpotMatDates, /* (I) Req spot maturity dates */
   long            *reqSwapFlags,    /* (I) 1=SwapRate, 0=CashRate */
   long             numRates,        /* (I) Len of all arrays */

   long             reqSwapFreq,     /* (I) Requested swap frequency */
   long             reqCashDCC,      /* (I) Requested cash day count conv */
   long             reqSwapDCC,      /* (I) Requested swap day count conv */
   long             couponInterpType,/* (I) Zero curve coupon interp type */
   long             zeroInterpType,  /* (I) Zero curve zero interp type */
   char             badDayConv,      /* (I) Bad day convention */
   char            *holidayFile)     /* (I) Holidays for bad day adjustment */
{
    long  idx;

    GtoErrMsg("\nLogged Inputs for DriFwdToSpotRatesOL:\n");

    GtoErrMsg("Req output swap freq:           %ld\n", reqSwapFreq);
    GtoErrMsg("Req output cash day count conv: %s\n", 
              GtoFormatDayCountConv(reqCashDCC));
    GtoErrMsg("Req output swap day count conv: %s\n", 
              GtoFormatDayCountConv(reqSwapDCC));
    GtoErrMsg("Zero curve coupon interp type:  %ld\n", couponInterpType);
    GtoErrMsg("Zero curve zero interp type:    %s\n", 
              tempGtoFormatInterpType(zeroInterpType));
    GtoErrMsg("Bad day convention:             %c\n", (char)badDayConv);
    GtoErrMsg("Holiday file:                   %s\n", holidayFile);

    GtoErrMsg(
    "  FwdStart     FwdMat     ReqMat  FwdPrd  FwdDCC  FwdRates\n");
    for (idx=0; idx < numRates; idx++)
    {
        GtoErrMsg("%10.10s %10.10s %10.10s%c %3s  %8s  %3.6f\n",
                  GtoFormatDate(fwdStartDates[idx]),
                  GtoFormatDate(fwdMatDates[idx]),
                  GtoFormatDate(reqSpotMatDates[idx]),
                  reqSwapFlags[idx] ? 'S' : 'M',
                  GtoFormatDateInterval(fwdPayIntervals[idx]),
                  GtoFormatDayCountConv(fwdDCCs[idx]),
                  fwdSwapRates[idx]);
    }

   GtoPrintTCurve(stubZC,"Stub Zero Curve");

   return;
}
 

/* Allocate n x n and compute tweak matrix - 
 * first tweak spot rates to get dFwd/dSpot,
 * invert the matrix using IMSL,
 * then multiply by the ratio of durations to
 * get positions
 * 10/9/98 - removed multiplication by durations
 */
static double *driDSpotOverDFwdMatrix(long    numRates, 
				     double *spotRates)
{
    static char      routine[]="driDSpotOverDFwdMatrix";
    int              status=FAILURE;  

    double *d = NULL;
    double *inv = NULL;
    long size = numRates*numRates;
    double *fwdDur=NULL, *spotDur=NULL;

    /*    double *ptr;
    double freq;
    double years;
    TDate valueDate = p.stubZC->fBaseDate;
    */
    long i,j;

    spotDur = NEW_ARRAY(double, numRates);
    fwdDur = NEW_ARRAY(double, numRates);
    if(spotDur==NULL || fwdDur==NULL)
	goto done;

    d = NEW_ARRAY(double, size);
    inv = NEW_ARRAY(double, size);
    if(d==NULL || inv==NULL)
	goto done;

    /* Precompute spot durations 
     */
    /*
    for(j=0; j<p.numCashRates; j++)
    {
	if(GtoMoneyMarketModDuration(spotRates[j],
				     p.spotDates[j]-valueDate,
				     (p.cashDCC==GTO_ACT_360? 360:365),
				     &spotDur[j]) == FAILURE)
	    goto done;
    }
    for(; j<numRates; j++)  
    {
	if(GtoDayCountFraction(valueDate,
			       p.spotSwapDates[j-p.numCashRates],
			       p.swapDCC,
			       &years) ==  FAILURE ||
	   GtoBondModDuration(spotRates[j],
			      spotRates[j],
			      (long)p.swapFreq,
			      years,
			      GTO_STUB_SIMPLE,
			      &spotDur[j]) == FAILURE)
	    goto done;
    }
    */
    /* Precompute fwd durations 
     * removed by traders' request
     */
    /*
    for(j=0; j<numRates; j++)  
    {
	if(GtoDayCountFraction(p.fwdStartDates[j],
			       p.fwdMatDates[j],
			       p.fwdDCCArray[j],
			       &years) ==  FAILURE ||
	   GtoDateIntervalToFreq(p.fwdPayIntervals[j],
				 &freq) ==  FAILURE ||
	   GtoBondModDuration(p.actFwdRates[j],
			      p.actFwdRates[j],
			      (long)freq,
			      years,
			      GTO_STUB_SIMPLE,
			      &fwdDur[j]) == FAILURE)
	    goto done;
    }
    */
    /* j-th column are the tweaks of all fwds w/respect to j-th spot
     */
    for (j=0; j<numRates; j++)
    {
	/* tweak up 
	 */
	spotRates[j] += DRI_TWEAK_SIZE; 
	matchFwdRates(numRates,
		   spotRates,
		   inv);
	if(p.status == FAILURE)
	    goto done;
	for (i=0; i<numRates; i++) /* fill in j-th column */
	{
	    d[j+i*numRates] = inv[i];
	}
	/* tweak down
	 */
	spotRates[j] -= 2.*DRI_TWEAK_SIZE; 
	matchFwdRates(numRates,
		   spotRates,
		   inv);
	if(p.status == FAILURE)
	    goto done;

	for (i=0; i<numRates; i++) /* fill in j-th column */
	{
	    double *temp = d+(j+i*numRates);
	    *temp = (*temp - inv[i])/(2.*DRI_TWEAK_SIZE);

	}
	/* Recover spot rate
	 */
	spotRates[j] += DRI_TWEAK_SIZE; 
    }

    /* Invert matrix
     */
    GtoImslCmathInit ();
    imsl_d_lin_sol_gen(numRates,
		       d,
		       NULL,
		       IMSL_INVERSE_USER, inv,
		       IMSL_INVERSE_ONLY,
		       0);
    if(GtoImslCmathErrorStatus () == FAILURE)
	goto done;

 /* remove adjustment for durations */
    /* for(i=0, ptr=inv; i<numRates; i++)
       {
       for(j=0; j<numRates; j++, ptr++)
       {
	    *ptr *= spotDur[i]/fwdDur[j];
	    }
	    }*/

#ifdef DRI_MATRIX_DEBUG
    for(i=0; i<numRates; i++)
    {
	for(j=0; j<numRates; j++)
	{
	    GtoErrMsg("%f\t", d[j+i*numRates]);
	}
    }

    GtoErrMsg("\n\nInverse matrix");

    for(i=0; i<numRates; i++)
    {
	GtoErrMsg("\n");

	for(j=0; j<numRates; j++)
	{
	    GtoErrMsg("%f\t", inv[j+i*numRates]);
	}
    }
    GtoErrMsg("\n");

#endif

    status = SUCCESS;

  done: 
    FREE_ARRAY(d);
    FREE(spotDur);
    FREE(fwdDur);

    if (status == FAILURE) 
    {
        GtoErrMsg ("%s: Failed.\n", routine); 
	FREE(inv);
	return (NULL);
    }

    return (inv); 
}


/* 
 **** Temporarily moved here until the libs move to
 *    ALIB 9.3
 *!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */ 

/*----------------------------------------------------------------------------
FUNCTION: GtoFormatInterpType
 
CREATED: 3/98 by D. Gallager
 
DESCRIPTION: Formats an interp type. Note that since these are used
             to encode objects, these strings must be understandable
             by routine GtoStringToInterpType above.
---------------------------------------------------------------------------*/
static char * tempGtoFormatInterpType
    (long interpType)              /* (I) The interp. type */
{
    switch(interpType)
    {
      case GTO_LINEAR_INTERP:
        return("Linear");
 
      case GTO_SPLINE_INTERP:
        return("Spline");
 
      case GTO_FLAT_FORWARDS:
        return("Flat_Forwards");
 
      case GTO_INSTANEOUS_LINFWDS:
        return("Instantaneous_LinFwds");
 
      case GTO_PARABOLIC_FORWARDS:
        return("Parabolic_Forwards");
 
      case GTO_CONST_SPOT_VOL_INTERP:
        return("Vol_Interp");
 
      case GTO_POLY_INTERP_2ND_ORDER:
        return("P2");
 
      case GTO_POLY_INTERP_3RD_ORDER:
        return("P3");
 
      case GTO_POLY_INTERP_4TH_ORDER:
        return("P4");
 
      case GTO_POLY_INTERP_5TH_ORDER:
        return("P5");
 
      case GTO_POLY_INTERP_6TH_ORDER:
        return("P6");
 
      case GTO_POLY_INTERP_7TH_ORDER:
        return("P7");
 
      case GTO_POLY_INTERP_8TH_ORDER:
        return("P8");
 
      case GTO_POLY_INTERP_9TH_ORDER:
        return("P9");
 
      case GTO_POLY_INTERP_10TH_ORDER:
        return("P10");
      default:
        return("UNKNOWN");
 
    } /* switch */
      /* NOTREACHED */
 
}
 
