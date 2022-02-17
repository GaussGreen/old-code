/****************************************************************************/
/*      Calculation of forward price. 	                            	    */
/* 	Author:	David Liu,	Oct. 1998				    */
/*      fwdprice.c                                                          */
/****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cgeneral.h"
#include "bastypes.h"
#include "macros.h"
#include "tcurve.h"
#include "ldate.h"
#include "convert.h"
#include "datelist.h"
#include "cerror.h"
#include "eqfwdpr.h"
#include "drieq.h"
#include "eqfwd.h"

#include "drlstd.h"
#include "drlvtype.h"


#define __LONDON__
/*#undef __LONDON__ */

#define DRI_INTERP_TYPE GTO_LINEAR_INTERP

void    dlinterp (long, double *, long, long, double, double);


/*f-@CDOC(catn="DR Wrapper")------------------------------------------------
 * Compute forward price given the equity static data and funding curve.
 * Currently two underlying routines are available: 
 * 1) the ALIB GtoEqForwardPrice function
 * 2) based on London code with proper modification for settlement and dividend
 *    dates adjustments.   
 */
GTO_EXPORT(int)
DrForwardPriceGen(
	double		spotPrice,	/* (I) spot asset price   */
	TEqStatData	*eqStatData,	/* (I) Equity static data */
	TCurve		*zcCoF,		/* (I) Funding zero curve */
	long	        busDayConv,     /* (I) Business Day Conv    */
	char	        *holidayFile,   /* (I) Holiday File       */
	long		numFwd,		/* (I) Number of fwd points */
	TDate		*fwdDate,	/* (I) Forward dates to compute fwdPrice
					 *     fwdDate[0] must be the spotDate
					 *     where spotPrice is known  */
	double		*fwdPrice)	/* (O) Forward price      */
{
	static  char routine[] = "DrForwardPriceGen"; 
	int     status = FAILURE; 

	TDividendList     *divList   = eqStatData->divList; 

#ifdef __LONDON__
	TEqStmPrivate     *stm       = (TEqStmPrivate *)eqStatData->stm;  
	char		  settleType = eqStatData->settleType; 
	int	i;

	long	numSettle; 
	long	numDiv = divList->fNumItems;

	TDate 	*exDivDates  = NULL;
	TDate 	*payDivDates = NULL;
	double	*divAmounts  = NULL;
	long	*divTypes    = NULL;

	if (settleType == 'R')
		numSettle = stm->stmPeriods[0];
	else
		numSettle = stm->lastTradingDates->fNumItems;

	if ((exDivDates  = NEW_ARRAY(TDate,  numDiv)) == NULL ||
	    (payDivDates = NEW_ARRAY(TDate,  numDiv)) == NULL ||
	    (divAmounts  = NEW_ARRAY(double, numDiv)) == NULL ||
	    (divTypes    = NEW_ARRAY(long,   numDiv)) == NULL )
		goto done;

	for (i=0; i<=divList->fNumItems-1; i++){
		exDivDates[i] = divList->fArray[i].exDivDate;
		payDivDates[i] = divList->fArray[i].payDivDate;
		divAmounts[i] = divList->fArray[i].divAmount;
		divTypes[i]   = divList->fArray[i].divType;
	}

	if (DrLondonForwardPrice(spotPrice,
				 numDiv,
				 divAmounts,
				 exDivDates, 
				 payDivDates, 
				 divTypes, 
				 numSettle,
				 stm->lastTradingDates->fArray,
				 stm->settlementDates->fArray,
				 settleType,
				 busDayConv,
				 holidayFile,
				 zcCoF,			 	
				 fwdDate,
				 numFwd,
				 fwdPrice) == FAILURE)
		goto done;

#else	/* use GtoEqForwardPrice */
	TDateList *fwdDateList = NULL;
	TCurve	  *zcDummy = NULL;
	TDate	  zcDummyDates[2];
	double	  zcDummyRates[2];

	if ((fwdDateList = GtoNewDateListFromDates(fwdDate,
					           numFwd)) == NULL)
		goto done;

	/* 
         * Dummy zero curve with 0 rates for cost of funding to be 
	 * consistent with EDG convention.
         */
	
	zcDummyDates[0] = zcCoF->fArray[0].fDate;
	zcDummyDates[1] = zcCoF->fArray[zcCoF->fNumItems - 1].fDate;
	
	zcDummyRates[0] = 0e0;
	zcDummyRates[1] = 0e0;

	if ((zcDummy = GtoMakeTCurve(zcCoF->fBaseDate,
				     zcDummyDates,
				     zcDummyRates,
			    	     2,
				     zcCoF->fBasis,
				     zcCoF->fDayCountConv)) == NULL)
		goto done;

	if (GtoEqForwardPrice(fwdDate[0],
			      spotPrice,
			      divList,
			      eqStatData->stm,
			      holidayFile,
			      busDayConv,
			      holidayFile,
			      busDayConv,
			      zcDummy,
			      zcCoF,
			      fwdDateList,
			      fwdPrice) == FAILURE)
		goto done;
	
#endif

        status = SUCCESS;

done:
	if (status != SUCCESS)
	    GtoErrMsg("%s: failed.\n", routine);

#ifdef __LONDON__
	FREE(exDivDates);
	FREE(payDivDates);
	FREE(divAmounts);
	FREE(divTypes);
#else
	GtoFreeDateList(fwdDateList);
	GtoFreeTCurve(zcDummy);
#endif

        return (status);

#undef __LONDON__
}



/*****  DrForwardPrice  ******************************************************/
/*
*       Compute forward price at each forward date given the funding
*       curve and dividend curve.  The first point of forward dates
*	fwdDate[0] should always be the spot date at which spot price 
*	is given.
*       interpolation to cope with commodities and equity.
*/
GTO_EXPORT(int)
DrLondonForwardPrice(
        double   spotPrice,       /* (I) spot price                          */
        int      NbDiv,           /* (I) Number of dividend / forward prices */
        double   *dividends,      /* (I) Dividends / Forward prices          */
        TDate    *divDates,       /* (I) ex-Dividend dates                   */
        TDate    *payDivDates,    /* (I) Dividend pay dates                  */
        long     *divTypes,       /* (I) Dividend / forward type             */
        int      NbSettle,        /* (I) Nb of settlement dates(<0 if rolling*/
        TDate    *LastTrading,    /* (I) Last trading date of settl period   */
        TDate    *SettleDate,     /* (I) Corresponding settlement date       */
        char     settleType,      /* (I) Settlement type:'F'ixed or 'R'olling*/
	long	 busDayConv,	  /* (I) Business day conv  		     */
	char	 *holidayFile,    /* (I) "NONE" for weekends only
                                   *     "No_Weekends" for no adjustments    */
	TCurve	 *zcCoF,	  /* (I) Funding Zero Curve	             */ 
        TDate    *fwdDates,       /* (I) Forward dates to compute fwd price  */
        int      NbFwd,           /* (I) Total number of time points         */
        double   *fwdPrices)      /* (O) Fwd price at eachtime forward date  */
{
#define __DEBUG__
#undef __DEBUG__
	static  char routine[] = "DrLondonForwardPrice";
	int	status = FAILURE; 

        double  y, dt; 	    

        int     i1,        /* Node index of the previous volatility point */
                i2,        /* Node index of the current volatility point  */
                i3,        /* Node index of the next volatility point     */
                i, j, k;  

	int 	numMerge;  /* Total number of combined dates from divDates
			      and fwdDates  */

	TDate	*adjFwdDates = NULL;
	double 	*fullDateZero = NULL;
	double 	*fullSettleDateZero = NULL;
	double 	*fullPrices = NULL;
	TDate 	*fullDates = NULL;
	TDate 	*fullSettleDates = NULL;

	TDate	spotDate;
	double	spotZero;

	double 	divDateZero, payDateZero;

	/* 
	 * Adjust the dividend for the delay of payment.
	 * If a dividend is paid after the ex-date, its value at 
 	 * ex-date is the PV over the period pay-date to ex-date.
	 * each discrete div is adjusted by this amount here.
	 */
	for (i=0; i <= NbDiv-1; i++)
	{
	    if (payDivDates[i] > divDates[i])
	    {
	      	if ((GtoDiscountDate(divDates[i],
                                     zcCoF,
                                     DRI_INTERP_TYPE,
                                     &divDateZero) == FAILURE) ||
	    	    (GtoDiscountDate(payDivDates[i],
                                     zcCoF,
                                     DRI_INTERP_TYPE,
                                     &payDateZero) == FAILURE))
                     goto done;
	
		dividends[i] *= payDateZero/divDateZero;
	    }
 	    else
	    if (payDivDates[i] < divDates[i])
	    {
		GtoErrMsg("%s: Dividend pay date (%s) before ex-date (%s).\n",
			  routine, GtoFormatDate(payDivDates[i]),
			  GtoFormatDate(divDates[i]));
		goto done;
	    }
	}

	/* holiday adjustment for the input forward dates */

	if ((adjFwdDates = NEW_ARRAY(TDate,  NbFwd)) == NULL)
	{
	    GtoErrMsg("%s: Error allocating memory for adjFwdDates.\n",
                      routine);
            goto done;
        }

	for (i = 0; i <= NbFwd-1; i++)
	{
	    if (GtoBusinessDay( fwdDates[i],
				busDayConv,
				holidayFile,
				&adjFwdDates[i]) == FAILURE)
		goto done;
	}

	spotDate = adjFwdDates[0];

	/*
	 *  Merge the adjFwdDates and divDates into one single array
	 *  fullDates and sort out in ascending order and remove
	 *  duplicate dates if exist.
	 */

	if ((fullDates = NEW_ARRAY(TDate,  NbFwd + NbDiv)) == NULL)
	{
	    GtoErrMsg("%s: Error allocating memory for fullDates.\n",
                      routine);
            goto done;
        }
	
	j = 0;

	/* only the divDates after spotDate will be included */
	for (i = 0; i <= NbDiv-1; i++)
	{
	    if (divDates[i] >= spotDate)
	    {
	    	fullDates[j] = divDates[i];
		j++;
	    }
	}

	for (i = 0; i<=NbFwd-1; i++, j++)
	    fullDates[j] = adjFwdDates[i];
	
	numMerge = j;

#ifdef  __DEBUG__
	GtoErrMsg("%s: Merged forward dates before sorting\n", routine);
	for (i=0; i<=numMerge-1; i++)
	    GtoErrMsg("\t%d\t%s\n", i+1, GtoFormatDate(fullDates[i]));
#endif	

	if (DrlVTypeVectSort(fullDates,
			     &numMerge,
			     DRL_TDATE_T,
			     TRUE) == FAILURE)
	    goto done;

#ifdef  __DEBUG__
	GtoErrMsg("%s: Merged forward dates after sorting\n", routine);
	for (i=0; i<=numMerge-1; i++)
	    GtoErrMsg("\t%d\t%s\n", i+1, GtoFormatDate(fullDates[i]));
#endif	

        /* 
         *  Calculate settlement date for each forward date.
         */

	if ((fullSettleDates  = NEW_ARRAY(TDate,  numMerge)) == NULL )
	{
            GtoErrMsg("%s: Error allocating memory for fullSettleDates.\n",
                      routine);
            goto done;
        }

        if (settleType == 'R')     /* Rolling settlement */
	{
	    for (i = 0; i <= numMerge-1; i++)                                  
            {        
		if(GtoDateFromBusDaysOffset(fullDates[i],
                                	    (long)NbSettle,
                                	    holidayFile,
                               		    &fullSettleDates[i]) == FAILURE)
            		goto done;
            }  
	}
        else    		   /* Fixed settlement */
        {
	    i1 = i2 = 0;
	    for (k = 0; k < NbSettle; k++)
            {        
		/* 
		 * All forward dates before current LastTrading[k]
		 * will settle on corresponding SettleDate[k].
		 */
		while ((fullDates[i2] <= LastTrading[k]) && (i2 <= numMerge-1))	
                    i2++; 

                for (i = i1; i < i2; i++)
                    fullSettleDates[i] = SettleDate[k];
                        
                i1 = i2;
            }  
                
	    /*
	     * If there are some forward date left after the last 
             * settlement date, just assume it settles immediately.
	     */
            for (i = i1; i <= numMerge-1; i++)        
                fullSettleDates[i] = fullDates[i];
        }  


	/* 
	 * Compute the zero coupon price at each forward date
	 * and forward settlement date. 
	 */

	if ((fullDateZero       = NEW_ARRAY(double, numMerge)) == NULL ||
	    (fullSettleDateZero = NEW_ARRAY(double, numMerge)) == NULL )
	{
	    GtoErrMsg("%s: Error allocating memory for fullDateZero.\n",
                      routine);
            goto done;
        }

	for (i=0; i <= numMerge-1; i++)
	{
	    if ((GtoDiscountDate(fullDates[i],
                                 zcCoF,
                                 DRI_INTERP_TYPE,
                                 &fullDateZero[i]) == FAILURE) ||
	        (GtoDiscountDate(fullSettleDates[i],
                                 zcCoF,
                                 DRI_INTERP_TYPE,
                                 &fullSettleDateZero[i]) == FAILURE))
                goto done;
		
	    /* discount between fwd date and spotDate
	     * This allows to handle the case where spot price
 	     * is known at forward starting date
	     * or the settlement date of spot trade is different
 	     * from that of value date of zero curves, 
	     * e.g. bond that follows T+1 settlement convention.
	     */
	    if (i==0) spotZero = fullDateZero[0];
	    fullDateZero[i] /= spotZero;
	    fullSettleDateZero[i] /= spotZero;
	}


	/* forward prices at each date of the full list */

	if ((fullPrices  = NEW_ARRAY(double, numMerge)) == NULL) 
	{
	    GtoErrMsg("%s: Error allocating memory for fullPrices.\n",
                      routine);
            goto done;
        }

        /* 
         *  Adjust quoted spot level for settlemnt delay and pay dividend.
	 *  Get zero corresponding to spot settlement date
         */
        fullPrices[0] = spotPrice * fullSettleDateZero[0];

#ifdef SPOTDIVIDEND
        /* Pay today's dividend if discrete */ 
        for (k = 0; k < NbDiv; k++)         
        {        
            if (fullDates[0] == divDates[k])
            {
		switch (divTypes[k]){
		    case GtoDIVIDENDRATE_AMOUNT:    /* Fixed $ dividend */
		    {
		    	fullPrices[0] -= dividends[k];
                    	break;
                    }
                    case GtoDIVIDENDRATE_PERCENT:   /* Fixed dividend yield */
                    {
		    	/* 
		     	 * Dividend yield dividend[k] is applied to the
		     	 * quoted stock price. 
		     	 */
                    	fullPrices[0] -= spotPrice * dividends[k];
                    	break;                                          
                    }
                    default:
		    {
			break;
	            }
                }  /* switch */
            } 
        } 
#endif

	
        /* 
         *	Calculate forward curve.
         */
	
        i1 = i2 = 0;
                          
        for (k = 0; k <= NbDiv-1; k++)
        {        
	   /* 
	    * If current type is a continuous dividend then 
	    * the following are as well 
	    */
           if (divTypes[k] == GtoDIVIDENDRATE_CONTINUOUS)    
		break;                                                  
        
           /* 
	    *  All forward dates before the current dividend date
	    */        
	   while ((fullDates[i2] < divDates[k]) && (i2 < numMerge-1))    
		i2++;     /* It exists as dividend dates are critical dates */

           if (i2 == i1)
           	continue;

           switch (divTypes[k])
           {
           case GtoDIVIDENDRATE_AMOUNT:          /* Fixed $ dividend */
           {
		/*
		 * Forward price is flat between two dividend dates
		 */
                for (i = i1 + 1; i <= i2; i++)
                    fullPrices[i] = fullPrices[i1];      

                /* 
		 * If fullDates[i2] = divDates[k] 
		 * Subtract PV of $ dividend from previous ex-dividend price
		 */
                if (fullDates[i2] == divDates[k])               
                    fullPrices[i2] -= dividends[k] * fullDateZero[i2]; 

                break;
           }
           case GtoDIVIDENDRATE_PERCENT:         /* Fixed dividend yield */
           {
                for (i = i1 + 1; i <= i2; i++)
                    fullPrices[i] = fullPrices[i1];
                        
                if (fullDates[i2] == divDates[k])
                    fullPrices[i2] -= (fullPrices[i2]/fullSettleDateZero[i2])  
				      * dividends[k] * fullDateZero[i2];	

                break;
           }
           default:          /* Linear interpolation for commodity */
           {
                for (i = i1 + 1; i <= i2; i++)
                {        
                    dlinterp ( fullDates[i], 
                               &y, 
                               fullDates[i1], 
                               fullDates[i2], 
                               fullPrices[i1], 
                               dividends[k]);
                                                        
		    /*
		     * Forward prices have to be discounted
		     */
                    fullPrices[i] = y * fullDateZero[i];   
                }  
           }

           }  /* switch */
        
           i1 = i2;

        }  /* for k */


        if (k < NbDiv)      /* There are some continuous dividend left */
        {
           while ((fullDates[i2] < divDates[k]) && (i2 < numMerge-1))
                i2++;

           /* Flat forward from last discrete dividend date */
           for (i = i1 + 1; i <= i2; i++)        
                fullPrices[i] = fullPrices[i1];
        }  

        i3 = 0;

        for (; k < NbDiv; k++)
        {
           if (k == NbDiv - 1) /* Only point left: we fill dates until last */
                i3 = numMerge - 1;
           else                /* Some continuous dividends left afterwards */
           {
                while ((fullDates[i3] < divDates[k+1]) && (i3 < numMerge-1))
                     i3++;    
           }  
                        
	   /* Continuous dividends are not input in % */
           for (i = i2 + 1; i <= i3; i++) {
		if (GtoDayCountFraction(fullDates[i-1],
					fullDates[i],
					GTO_ACT_365F,
					&dt) == FAILURE)
			goto done;
                fullPrices[i] = fullPrices[i-1] / pow (1 + dividends[k], dt);   
	   }

           i2 = i3;
        }  

        /* 
	 * Flat forward stock after the last discrete dividend 
	 * or forward price 
	 */
        for (i = i2 + 1; i <= numMerge-1; i++)     
            fullPrices[i] = fullPrices[i1];
                                   
        /* 
         * Forward value forward curve.
         */

        for (i = 0; i <= numMerge-1; i++)
            fullPrices[i] /= fullSettleDateZero[i];     

	/*
	 * Select the forward prices at desired input forward dates
	 */
	for (k = 0; k <= NbFwd-1; k++)
	for (i = 0; i <= numMerge-1; i++)
	{ 
	    if (adjFwdDates[k] == fullDates[i])
	    {
		fwdPrices[k] = fullPrices[i];
		break;
	    }
	}	

        status = SUCCESS;

done:
	if (status != SUCCESS)
	    GtoErrMsg("%s: failed.\n", routine);

	FREE(fullDates);
	FREE(fullSettleDates);
        FREE(fullPrices);
        FREE(fullDateZero);
        FREE(fullSettleDateZero);
	FREE(adjFwdDates);

        return (status);

#undef __DEBUG__

} 



/*****  DrForwardDiscZero  **************************************************/
/*
*       Compute discount factor at each forward date given the discount
*       curve.  The first point of forward dates
*	fwdDate[0] should always be the base date. 
*/
GTO_EXPORT(int)
DriForwardDiscZero(
	TEqStatData	*eqStatData,	/* (I) Equity static data */
	TCurve		*discZcCurve,	/* (I) Discount zero curve */
	char	        *holidayFile,   /* (I) Holiday File       */
	long		numFwd,		/* (I) Number of fwd points */
	TDate		*fwdDate,	/* (I) Forward dates */
        double   	*fwdDiscZero)   /* (O) disc zero at each fwd date  */
{
	static  char routine[] = "DriForwardDiscZero";
	int	status = FAILURE; 

	TEqStmPrivate     *stm       = (TEqStmPrivate *)eqStatData->stm;  
	char		  settleType = eqStatData->settleType; 

	TDate		  *LastTrading = NULL;
        TDate             *SettleDate  = NULL;

	TDate   *fwdSettleDate = NULL;

	int	i, i1, i2, k;

	long	numSettle; 

	LastTrading = stm->lastTradingDates->fArray;
        SettleDate  = stm->settlementDates->fArray;

	if ((fwdSettleDate  = NEW_ARRAY(TDate,  numFwd)) == NULL)
		goto done;

	if (settleType == 'R')
		numSettle = stm->stmPeriods[0];
	else
		numSettle = stm->lastTradingDates->fNumItems;

        if (settleType == 'R')     /* Rolling settlement */
	{
	    for (i = 0; i <= numFwd-1; i++)                                  
            {        
		if(GtoDateFromBusDaysOffset(fwdDate[i],
                                	    numSettle,
                                	    holidayFile,
                               		    &fwdSettleDate[i]) == FAILURE)
            		goto done;
            }  
	}
        else    		   /* Fixed settlement */
        {
	    i1 = i2 = 0;
	    for (k = 0; k <= numSettle-1; k++)
            {        
		/* 
		 * All forward dates before current LastTrading[k]
		 * will settle on corresponding SettleDate[k].
		 */
		while ((fwdDate[i2] <= LastTrading[k]) && (i2 <= numFwd-1))	
                    i2++; 

                for (i = i1; i < i2; i++)
                    fwdSettleDate[i] = SettleDate[k];
                        
                i1 = i2;
            }  
                
	    /*
	     * If there are some forward date left after the last 
             * settlement date, just assume it settles immediately.
	     */
            for (i = i1; i <= numFwd-1; i++)        
                fwdSettleDate[i] = fwdDate[i];
        }  

	for (i = 0; i <= numFwd-1; i++)                                  
        {        
            if(GtoDiscountDate(fwdSettleDate[i],
                               discZcCurve,
                               DRI_INTERP_TYPE,
			       &fwdDiscZero[i]) == FAILURE)
		goto done;
	}
                           
        status = SUCCESS;

done:
	if (status != SUCCESS)
	    GtoErrMsg("%s: failed.\n", routine);
	
	FREE(fwdSettleDate);

        return (status);

} 
	


/*****  dlinterp  ***********************************************************/
/*
*       Linear interpolation between 2 Dates.
*/
void    dlinterp (  long d,
                    double *y,
                    long d0,
                    long d1,
                    double y0,
                    double y1)
{
	long	a, b;
        double  l0, l1;

	GtoDaysDiff(d, d1, GTO_ACT_365, &a);
	GtoDaysDiff(d0, d1, GTO_ACT_365, &b);
        l0 = ((double) a)/((double) b);

	GtoDaysDiff(d0, d, GTO_ACT_365, &a); 
	GtoDaysDiff(d0, d1, GTO_ACT_365, &b);
        l1 = ((double) a)/((double) b);

        *y=l0*y0+l1*y1;

        return;

}  /* dlinterp */

