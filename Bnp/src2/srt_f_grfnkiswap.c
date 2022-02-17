/* ===============================================================================
   
   FILENAME:       srt_f_grfnkiswap.c

   PURPOSE:        Pricing of a Knock In Swap with Grfn for a given model
                   Builds a tableau and prices the deal

   =============================================================================== */
#include "grf_h_all.h"
#include "srt_h_all.h"
#include "grf_h_public.h"
#include "srt_h_grfnkiswap.h"
#include "srtaccess.h"


/* ------------------------------------------------------------------------------------------- */
/* Cell in the first column (c[0,...]) of the knock in swap tableau : the KI Index FRA (a[0],a[1]) */
static void GrfnKiSwapIndexCell(
			String	    szIndexBasis,
			long         lIndexNumMonths,
			String      szUndName,
			String      szIndexRefRateCode,
			String      szIndexCell)
{
	sprintf(szIndexCell,
			"FRA(a[0,i],%d,\"%s\",\"%s\",\"%s\")",
			lIndexNumMonths,
			szIndexBasis,
			szUndName,
			szIndexRefRateCode);

}

/* ------------------------------------------------------------------------------------------ */
/* Cell in the second column (c[1,...]) of the knock in swap tableau : the Call Spread for PV > 0 */
static void GrfnKiSwapPositivePvCallSpreadCell(
			SrtBarrierType  eKnockUpOrDown,
			double          dCallSpreadWidth,
			String      szPositivePvCallSpreadCell)
{
	if (eKnockUpOrDown == SRT_UP)
	{
		if ( dCallSpreadWidth > 0.0)
		{
			sprintf(szPositivePvCallSpreadCell,
				"(MAX(c[0,i]-a[1,i],0)-MAX(c[0,i]-a[1,i]-%f,0))/%f",
				dCallSpreadWidth,
				dCallSpreadWidth);
		}
		else
		{
			sprintf(szPositivePvCallSpreadCell,
				"if ( c[0,i]> a[1,i], 1, 0)");
		}
	}
	else
	{
		if ( dCallSpreadWidth > 0.0)
		{
			sprintf(szPositivePvCallSpreadCell,
				"(MAX(a[1,i]-c[0,i],0)-MAX(a[1,i]-%f-c[0,i],0))/%f",
				dCallSpreadWidth,
				dCallSpreadWidth);
		}
		else
		{
			sprintf(szPositivePvCallSpreadCell,
				"if ( c[0,i] < a[1,i], 1, 0)");
		}
	}

}

/* ----------------------------------------------------------------------------------------- */
/* Cell in the third column (c[2,...]) of the knock in swap tableau : the Call Spread for PV < 0 */
static void GrfnKiSwapNegativePvCallSpreadCell(
			SrtBarrierType  eKnockUpOrDown,
			double          dCallSpreadWidth,
			String      szNegativePvCallSpreadCell)
{
	if (eKnockUpOrDown == SRT_UP)
	{
		if ( dCallSpreadWidth > 0.0)
		{
		
			sprintf(szNegativePvCallSpreadCell,
				"(MAX(c[0,i]-a[1,i]+%f,0)-MAX(c[0,i]-a[1,i],0))/%f",
				dCallSpreadWidth,
				dCallSpreadWidth);
		}
		else
		{
			sprintf(szNegativePvCallSpreadCell,
				"if ( c[0,i] > a[1,i], 1, 0)");
		}
	}
	else
	{
		if ( dCallSpreadWidth > 0.0)
		{
		
			sprintf(szNegativePvCallSpreadCell,
				"(MAX(a[1,i] + %f-c[0,i],0)-MAX(a[1,i]-c[0,i],0))/%f",
				dCallSpreadWidth,
				dCallSpreadWidth);
		}
		else
		{
			sprintf(szNegativePvCallSpreadCell,
				"if ( c[0,i] < a[1,i], 1, 0)");
		}
	}

}

/* ----------------------------------------------------------------------------------------- */
/* Cell in the fourth column (c[3,...]) of the knock in swap tableau : the Fixed Leg Value */
static void GrfnKiSwapFixedLegPvCell(
			SRT_Boolean      bStubPaidAtStart,
			String      szFixedLegPvCell)
{
	if (bStubPaidAtStart == SRT_NO)
		sprintf(szFixedLegPvCell,
			"PVRNG(6,7,a[2,i])-a[5,i]*df(now,a[3,i])");
	else
	if (bStubPaidAtStart == SRT_YES)
		sprintf(szFixedLegPvCell,
			"PVRNG(6,7,a[2,i])-a[5,i]");

}

/* ----------------------------------------------------------------------------------------- */
/* Cell in the fifth column (c[4,...]) of the knock in swap tableau : the Float Leg Value */
static void GrfnKiSwapFloatLegPvCell(
			SrtReceiverType     eRecPay,
			String             szFloatLegPvCell)
{
	if (eRecPay == SRT_RECEIVER)
		sprintf(szFloatLegPvCell,
			"PVRNG(8,9,a[2,i])-df(now,a[2,i])*a[4,i]");
	else
	if (eRecPay == SRT_PAYER)
		sprintf(szFloatLegPvCell,
			"PVRNG(8,9,a[2,i])+df(now,a[2,i])*a[4,i]");

}

/* ----------------------------------------------------------------------------------------- */
/* Cell in the sixth column (c[5,...]) of the knock in swap tableau : the Full Swap Value */
static void GrfnKiSwapFullSwapPvCell(
			String             szFullSwapPvCell)
{
	sprintf(szFullSwapPvCell,
		"c[3,i]+c[4,i]");
	
}

/* ----------------------------------------------------------------------------------------- */
/* Cell in the seventh column (c[6,...]) of the knock in swap tableau : have we knocked in ? */
static void GrfnKiSwapKnockedInCell(
			String             szKnockedInCell)
{
	sprintf(szKnockedInCell,
			"if(va>0,va:=2,va:= if(c[5,i]>0,c[1,i],c[2,i]))");
	
}

/* ----------------------------------------------------------------------------------------- */
/* First Cell in the seventh column (c[6,...]) of the knock in swap tableau : va initialisation */
static void GrfnKiSwapKnockedInFirstRowCell(
			String             szKnockedInFirstRowCell)
{
	sprintf(szKnockedInFirstRowCell,
			"va:= if(c[5,i]>0,c[1,i],c[2,i])");
	
}
/* ----------------------------------------------------------------------------------------- */
/* Cell in the eight column (c[7,...]) of the knock in swap tableau : the Cash Flow */
static void GrfnKiSwapCashFlowCell(
			String             szCashFlowCell)
{
	sprintf(szCashFlowCell,
			"if(va>1.5,0,va*c[5,i])");
	
}

/* ----------------------------------------------------------------------------------------- */


/* Function that builds the full grfn tableau for a Knock In Swap */
static Err GrfnKnockInSwapMakeTableau(
			long            lNumEventDates,
			String	       szIndexBasis,
			long            lIndexNumMonths,
			String         szUndName,
			String         szIndexRefRateCode,
			SrtReceiverType eRecPay,
			SrtBarrierType  eKnockUpOrDown,
			double          dCallSpreadWidth,
			SRT_Boolean         bStubPaidAtStart,
			long          *plNumRows,
			long          *plNumCols,
			GrfnCell  ***pppsTableau
		)
{
Err 	err=NULL;
int     i;

/* Sets the tableau dimensions : 8 columns */
	*plNumCols 		= 8;
	*plNumRows 		= lNumEventDates;

/* Allocate space for the Tableau (including strings) */	
	*pppsTableau 	= GrfnCellmatrix(*plNumRows,*plNumCols,GRFN_DEF_ARGBUFSZ);

/* Fills in the tableau row by row */
	for (i=0; i < *plNumRows; i++)
	{
	/* The first column : the Kknock In Index (a FRA ) */
		GrfnKiSwapIndexCell(szIndexBasis, lIndexNumMonths, szUndName, szIndexRefRateCode, (*pppsTableau)[i][0].sval);
		(*pppsTableau)[i][0].type = GRFNSCELL;

	/* The second column: the call spread value to be used for a positive PV */
		GrfnKiSwapPositivePvCallSpreadCell(eKnockUpOrDown,dCallSpreadWidth, (*pppsTableau)[i][1].sval);
		(*pppsTableau)[i][1].type = GRFNSCELL;
	
	/* The third column: the call spread value to be used for a negative PV */
		GrfnKiSwapNegativePvCallSpreadCell(eKnockUpOrDown,dCallSpreadWidth, (*pppsTableau)[i][2].sval);
		(*pppsTableau)[i][2].type = GRFNSCELL;
	
	/* The fourth column: the fixed leg pv (as a range of CF and dates) */
		GrfnKiSwapFixedLegPvCell(bStubPaidAtStart, (*pppsTableau)[i][3].sval);
		(*pppsTableau)[i][3].type = GRFNSCELL;

	/* The fifth column: the float leg pv (as exchange of notionals and spreads) */
		GrfnKiSwapFloatLegPvCell(eRecPay, (*pppsTableau)[i][4].sval);
		(*pppsTableau)[i][4].type = GRFNSCELL;
	
	/* The sixth column: the full swap PV */
		GrfnKiSwapFullSwapPvCell( (*pppsTableau)[i][5].sval);
		(*pppsTableau)[i][5].type = GRFNSCELL;
	
	/* The seventh column: a variable to know if we have knocked in  */
		if (i==0)
			GrfnKiSwapKnockedInFirstRowCell((*pppsTableau)[i][6].sval);
		else
			GrfnKiSwapKnockedInCell( (*pppsTableau)[i][6].sval);
		(*pppsTableau)[i][6].type = GRFNSCELL;

	/* The eight and last column: the cash flow */
		GrfnKiSwapCashFlowCell( (*pppsTableau)[i][7].sval);
		(*pppsTableau)[i][7].type = GRFNSCELL;
	
	} /* END loop on all tableau rows */
                                             
	return err;
}
					
/* ---------------------------------------------------------------------------- */

#define FREE_GRFN_KNOCK_IN_SWAP_AUXILIARY_MEMORY {\
		for (i =0; i< (*plNumAuxColumns);i++){\
			if ((*pppdAuxRanges)[i])\
				free_dvector( (*pppdAuxRanges)[i], 0, (*pplAuxRangesLength)[i] -1 );\
			(*pppdAuxRanges)[i] = NULL;}\
		srt_free((*pppdAuxRanges));\
		srt_free((*pplAuxRangesLength));\
	}

/* Function that builds the Auxiliary Ranges for the Knock In Swap */
static Err GrfnKnockInSwapMakeAuxiliary(
		double  ***pppdAuxRanges,
		long      **pplAuxRangesLength,
		long        *plNumAuxColumns,
		long        *plFixedStartDates,
		long        *plFixedPayDates,
		long          lNumFixedDates,
		double      *pdFixedCoupons,
		double      *pdNotionals,
		SrtBasisCode  eFixedBasis,
		SrtReceiverType eRecPay,
		long        *plFloatStartDates,
		long        *plFloatPayDates,
		long          lNumFloatDates,
		double        dFloatMargin,
		String       szFloatRefRate,
		SrtBasisCode  eFloatBasis,
		double      *pdKnockInLevels,
		long        *plKnockInIndexStartDates,
		long        *plKnockInSwapStartDates,
		long          lNumKnockInDates
		)
{
double           dThisPeriodNotional;
double           dNextPeriodNotional;
double           dSpread;
double           dCoverage;
Err              err = NULL;
int              i,j;
SRT_Boolean          bExitFlag;

/* Allocate Memory for the Auxiliary ranges */
	(*plNumAuxColumns) = 10;
	(*pplAuxRangesLength) = (long *) srt_malloc ( (*plNumAuxColumns) * sizeof( long));
	(*pppdAuxRanges) = (double **) srt_malloc( (*plNumAuxColumns) * sizeof(double**));
	for (i =0; i< (*plNumAuxColumns);i++)
		(*pppdAuxRanges)[i] = NULL;
	
/* Aux 0: the Index Start Dates */
	(*pplAuxRangesLength)[0] = lNumKnockInDates;
	(*pppdAuxRanges)[0] = dvector( 0, (*pplAuxRangesLength)[0] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[0]; i++)
	{
		(*pppdAuxRanges)[0][i] = (double) plKnockInIndexStartDates[i];
	}	
	
/* Aux 1: the Knock In Levels */
	(*pplAuxRangesLength)[1] = lNumKnockInDates;
	(*pppdAuxRanges)[1] = dvector( 0, (*pplAuxRangesLength)[1] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[1]; i++)
	{
		(*pppdAuxRanges)[1][i] = (double) pdKnockInLevels[i];
	}

/* Aux 2: the Knock In Swap Start Dates */
	(*pplAuxRangesLength)[2] = lNumKnockInDates;
	(*pppdAuxRanges)[2] = dvector( 0, (*pplAuxRangesLength)[2] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[2]; i++)
	{
		(*pppdAuxRanges)[2][i] = (double) plKnockInSwapStartDates[i];
	}

/* Aux 3: the Swap Corresponding Next Pay Dates */
	(*pplAuxRangesLength)[3] = lNumKnockInDates;
	(*pppdAuxRanges)[3] = dvector( 0, (*pplAuxRangesLength)[3] -1 );
	j=0;
	for ( i =0; i < (*pplAuxRangesLength)[3]; i++)
	{
		while (plFixedPayDates[j] <= plKnockInSwapStartDates[i] )
		{
			j++;
			if (j==lNumFixedDates)
			{
				FREE_GRFN_KNOCK_IN_SWAP_AUXILIARY_MEMORY;
				return serror("Some KI dates are after last Swap Pay Date");
			}
		}
		(*pppdAuxRanges)[3][i] = (double) plFixedPayDates[j];
	}

/* Aux 4: the Corresponding Swap Notional (based on the notional for the payment) */
	(*pplAuxRangesLength)[4] = lNumKnockInDates;
	(*pppdAuxRanges)[4] = dvector( 0, (*pplAuxRangesLength)[4] -1 );
	j=0;
	for ( i =0; i < (*pplAuxRangesLength)[4]; i++)
	{
		while (plFixedPayDates[j] <= plKnockInSwapStartDates[i] )
		{
			j++;
			if (j==lNumFixedDates)
			{	
				FREE_GRFN_KNOCK_IN_SWAP_AUXILIARY_MEMORY;
				return serror("Some KI dates are after last Swap Pay Date");
			}
		}
		(*pppdAuxRanges)[4][i] = pdNotionals[j];
	}

/* Aux 5: the Coupon Stub for the Fixed Payment (with Notional - from SwapStart to KiSwapStart) */
	(*pplAuxRangesLength)[5] = lNumKnockInDates;
	(*pppdAuxRanges)[5] = dvector( 0, (*pplAuxRangesLength)[5] -1 );
	j=0;
	for ( i =0; i < (*pplAuxRangesLength)[4]; i++)
	{
		while (plFixedPayDates[j] <= plKnockInSwapStartDates[i] )
		{
			j++;
			if (j==lNumFixedDates)
			{
				return serror("Some KI dates are after last Swap Pay Date");
				FREE_GRFN_KNOCK_IN_SWAP_AUXILIARY_MEMORY;
			}
		}
		(*pppdAuxRanges)[5][i] = coverage(plFixedStartDates[j],plKnockInSwapStartDates[i],eFixedBasis)
				* pdNotionals[j] * pdFixedCoupons[j];
		if (eRecPay == SRT_PAYER)
			(*pppdAuxRanges)[5][i] *= -1.0;
	}

/* Aux 6: the Fixed Leg Payment Dates */
	(*pplAuxRangesLength)[6] = lNumFixedDates;
	(*pppdAuxRanges)[6] = dvector( 0, (*pplAuxRangesLength)[6] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[6]; i++)
	{
		(*pppdAuxRanges)[6][i] = (double) plFixedPayDates[i];
	}

/* Aux 7: the Fixed Leg Payments (Full Coupon * Coverage * Notional) */
	(*pplAuxRangesLength)[7] = lNumFixedDates;
	(*pppdAuxRanges)[7] = dvector( 0, (*pplAuxRangesLength)[7] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[7]; i++)
	{
		(*pppdAuxRanges)[7][i] =  pdFixedCoupons[i] * pdNotionals[i] * 
			 coverage(plFixedStartDates[i],plFixedPayDates[i],eFixedBasis);
		if (eRecPay == SRT_PAYER)
			(*pppdAuxRanges)[7][i] *= -1.0;
	}
	
/* Aux 8: the Floating Leg Payment Dates */
	(*pplAuxRangesLength)[8] = lNumFloatDates;
	(*pppdAuxRanges)[8] = dvector( 0, (*pplAuxRangesLength)[8] -1 );
	for ( i =0; i < (*pplAuxRangesLength)[8]; i++)
	{
		(*pppdAuxRanges)[8][i] = (double) plFloatPayDates[i];
	}

/* Aux 9: the equivalent Float Payments (Spreads + Margin + Notional Exchange for each period) */
	(*pplAuxRangesLength)[9] = lNumFloatDates;
	(*pppdAuxRanges)[9] = dvector( 0, (*pplAuxRangesLength)[9] -1 );
	if (plFixedPayDates[lNumFixedDates-1] != plFloatPayDates[lNumFloatDates-1])
	{
		FREE_GRFN_KNOCK_IN_SWAP_AUXILIARY_MEMORY;
		return serror("Last Payment Dates do not match between Fixed and Floating");
	}
	dNextPeriodNotional = 0.0;
	j=lNumFixedDates -1;
	bExitFlag = SRT_NO;
	for ( i =(*pplAuxRangesLength)[9] - 1; i >= 0; i--)
	{
	/* Finds the Notional for this period (according to fixed pay dates) */
		while (!bExitFlag)
		{
			if (plFixedPayDates[j] == plFloatPayDates[i] )
			{
				bExitFlag = SRT_YES;
			}
			else
			if (j == 0)
			{
				bExitFlag = SRT_YES;
			}
			else
			if (plFixedPayDates[j-1] < plFloatPayDates[i] )
			{
				bExitFlag = SRT_YES;
			}
			else
			{
				j--;
			}
		}
		dThisPeriodNotional = pdNotionals[j];

	/* The Spread and Margin on the Floating Rate */
		dSpread = swp_f_spread( plFloatStartDates[i], plFloatPayDates[i], szFloatRefRate); 
		dCoverage = coverage(plFloatStartDates[i],plFloatPayDates[i],eFloatBasis);
	
	/* The Floating Leg (REC) Representation: Notional at Beginning + Margin at End - Notional at End */
		(*pppdAuxRanges)[9][i] =  - dThisPeriodNotional 
				+ dCoverage * (dSpread + dFloatMargin) * dThisPeriodNotional
				+ dNextPeriodNotional;
		if (eRecPay == SRT_RECEIVER)
			(*pppdAuxRanges)[9][i] *= -1.0;
	/* Keep track of the Notional */
		dNextPeriodNotional = dThisPeriodNotional;
	
	}

	return NULL;
}

/* ---------------------------------------------------------------------------

                                 MAIN FUNCTION
   
	Prices a KnockIn Swap in Grfn, building the Tableau internally
   --------------------------------------------------------------------------- */

#define FREE_GRFN_KNOCK_IN_SWAP_MEMORY {\
		if(*pppsGrfnTableau)\
			grfn_free_GrfnCellmatrix(*pppsGrfnTableau,*plNumRows,*plNumCols);\
		*pppsGrfnTableau = NULL;\
		for (i =0; i< *plNumAuxColumns;i++){\
			if ((*pppdAuxRanges)[i])\
				free_dvector( (*pppdAuxRanges)[i], 0, (*pplAuxRangesLength)[i] -1 );\
			(*pppdAuxRanges)[i] = NULL;}\
		srt_free((*pppdAuxRanges));\
		srt_free((*pplAuxRangesLength));\
	}


Err SrtGrfnKnockInSwapPrice(
	long            *plKnockInFixingDates,
	long            *plKnockInIndexStartDates,
	long            *plKnockInSwapStartDates,
	double          *pdKnockInLevels,
	long              lNumKnockInDates,
	
	String           szIndexRefRate,
	String           szKnockUpOrDown,
	double            dCallSpreadWidth,
	SRT_Boolean           bStubPaidAtStart,

	long            *plFixedStartDates,
	long            *plFixedPayDates,
	double          *pdFixedCoupons,
	double          *pdNotionals,
	long              lNumFixedDates,
	String           szFixedBasis,
	String           szRecPay,
		
	long            *plFloatStartDates,
	long            *plFloatPayDates,
	long              lNumFloatDates,
	double            dFloatMargin,
	String           szFloatRefRate,
	
	String           szUndName,

	String         *pszGrfnParamNames,
	String         *pszGrfnParamValues,
	int               iNumGrfnParams, 
	
/* Outputs from Grfn */	
	double 	        *pdKnockInSwapPrice,
	long            *plNumRows,
	long            *plNumCols,
	GrfnCell 	***pppsGrfnTableau,
	double        ***pppdLastPath,
	long            *plNumAuxColumns,
	long          **pplAuxRangesLength,
	double      ***pppdAuxRanges
	)			
{
Err                err		= NULL;
SrtGrfnParam     sGrfnParams;

SrtUndPtr         sUndPtr;
long             lNumEventdates	= 0;
SrtIOStruct    *psIOList;

SrtReceiverType   eRecPay;
SrtBarrierType    eKnockUpOrDown;
SrtBasisCode      eFixedBasis;
SrtBasisCode      eFloatBasis;
SrtBasisCode      eIndexBasis;
SrtCompounding    eFloatFrequency;
SrtCompounding    eIndexFrequency;
String           szIndexBasis;
int               i;
long              lToday;
long             lIndexNumMonths;

		
/* Gets the Underlying from its name */
	sUndPtr = lookup_und(szUndName);
		
/* Gets today from the Underlying */
	lToday = get_today_from_underlying(sUndPtr);

/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(
				iNumGrfnParams,
				pszGrfnParamNames,
				pszGrfnParamValues,
				&sGrfnParams))
	{
		return err;
	}

/* Sets the Receiver/Payer type */
	err = interp_rec_pay(szRecPay,&eRecPay);
	if (err)
		return err;

/* Sets the Basis type for the fixed leg*/
	err = interp_basis(szFixedBasis,&eFixedBasis);
	if (err)
		return err;

/* Sets the Up/Down type */
	err = interp_barrier_type(szKnockUpOrDown,&eKnockUpOrDown);
	if (err)
		return err;

/* Gets the Floating Reference Rate information : basis, frequency */
	err = swp_f_get_ref_rate_details(szFloatRefRate, &eFloatBasis, &eFloatFrequency);
	if(err)
		return err;

/* Gets the Index Reference Rate information : basis, frequency */
	err = swp_f_get_ref_rate_details(szIndexRefRate, &eIndexBasis, &eIndexFrequency);
	if(err)
		return err;
	lIndexNumMonths = (long) ( 12 / eIndexFrequency);

/* Transforms the Index basis back into a string */
	translate_basis( &szIndexBasis, eIndexBasis);

/* Create a I/O list for the prices */
	err = srt_f_IOstructcreate(&psIOList,"");
	if (err)
		return(err);


/* Remove all the KI dates that are before Today (including the e.o.day flag) */
	while (	plKnockInFixingDates[0] < lToday + sGrfnParams.end_of_day_flg)
	{
		lNumKnockInDates--;
		plKnockInFixingDates++;
		plKnockInIndexStartDates++;
		plKnockInSwapStartDates++;
		pdKnockInLevels++;
	}

/* Construct a Grfn Tableau that describes the Flexi Cap  */
	*pppsGrfnTableau	= NULL;
	err = GrfnKnockInSwapMakeTableau(
			lNumKnockInDates,
			szIndexBasis,
			lIndexNumMonths,
			szUndName,
			szIndexRefRate,
			eRecPay,
			eKnockUpOrDown,
			dCallSpreadWidth,
			bStubPaidAtStart,
			plNumRows,
			plNumCols,
			pppsGrfnTableau
		);
	if (err)
	{
		srt_f_IOstructfree(&psIOList);
		return err;
	}

/* Builds the Auxiliary ranges used by the Tableau */
	err = GrfnKnockInSwapMakeAuxiliary(
				pppdAuxRanges,
				pplAuxRangesLength,
				plNumAuxColumns,
				plFixedStartDates,
				plFixedPayDates,
				lNumFixedDates,
				pdFixedCoupons, 
				pdNotionals,
				eFixedBasis,
				eRecPay,
				plFloatStartDates,
				plFloatPayDates,
				lNumFloatDates,
				dFloatMargin,
				szFloatRefRate,
				eFloatBasis,
				pdKnockInLevels,
				plKnockInIndexStartDates,
				plKnockInSwapStartDates,
				lNumKnockInDates
				);
	if (err)
	{
		FREE_GRFN_KNOCK_IN_SWAP_MEMORY;	
		return err;
	}

/* Allocate space for the Grfn Cells (results of the last path) */
	(*pppdLastPath) = dmatrix(0, (*plNumRows) -1 , 0, (*plNumCols) -1 );

/* Call GRFN to value the deal */
	err = srt_f_grfn(
				sUndPtr,
				&sGrfnParams,
				lNumKnockInDates,
				&plKnockInFixingDates,
				plNumRows,
				plNumCols,
				pppsGrfnTableau, 
				0,
				NULL,
				(*plNumAuxColumns),
				(*pplAuxRangesLength),
				(*pppdAuxRanges),
				psIOList,
				(*pppdLastPath),
				0);
	if (err)
	{
		FREE_GRFN_KNOCK_IN_SWAP_MEMORY;	
		return err;
	}

/* Gets the PV from the I/O list */
	err = srt_f_IOstructgetpremiumval(*psIOList, pdKnockInSwapPrice);
	if (err)
	{
		FREE_GRFN_KNOCK_IN_SWAP_MEMORY;	
		return err;
	}


/* Returns a success message */	
	return err;

} /* END Err SrtGrfnKnockInSwapPrice(...) */


/* ------------------------------------------------------------------------------ */
		
/* ---------------------------------------------------------------------------

                                 MAIN FUNCTION
   
	Calibrates on a set of Caps/Floors for a Flexi Cap
   --------------------------------------------------------------------------- */


Err SrtGrfnKnockInSwapCalibrate(
	long            *plIndexFixingDates,
	long            *plIndexStartDates,
	double          *pdKnockInLevels,
	long              lNumDates,
	String           szIndexRefRateCode,
	String           szKnockUpOrDown,
	
	String           szModelName,
	double            dTau,
	SRT_Boolean           bUseTwoFactor,
	double            dAlpha,
	double            dGamma,
	double            dRho,
	String         *pszGrfnParamNames,
	String         *pszGrfnParamValues,
	int               iNumGrfnParams, 
	String           szYieldCurveName,
	Err              (*pfGetVol)(double dStart, double dEnd, double dStrike, 
						double dForward, double dSpread, double *pdBsVol),
	String           szVolType,
/* Outputs from this calibration */	
	String           szUndName,
	double 	      **ppdCapletPrices,
	long            *plNumCaplets
	)
{
Err err = NULL;
int         iNumSigmaRows;
int         iNumSigmaCols;
double  **ppdSigmaCurve;
int         iNumTauRows;
int         iNumTauCols;
double  **ppdTauCurve;
char        szRealModelName[64];
SrtBasisCode     eBasis;
SrtCompounding   eFrequency;
String          szBasis;
String          szFrequency;
Date             lToday;			
SrtCurvePtr      sCurve;
int              iNumInstruments;
long           *plNfp;
String        *pszFrequency;
String        *pszBasis;
double         *pdStrike;
double         *pdBondStrike;
String        *pszProductType;
String        *pszRecPay;
String        *pszRefRate;
double         *pdOptionPrice;
double		   *pdATMOptionPrice;
double			forward;
SrtGrfnParam     sGrfnParams;
int              i;
char            szCapFloor[64];
SrtBarrierType   eKnockUpOrDown;
char            szPremium[32];

/* Builds the Full Model Name */
	if (szModelName == NULL)
	{
		sprintf(szRealModelName, "LGM");
	}
	else
	{
		sprintf(szRealModelName, szModelName);
	}
	
	if (bUseTwoFactor == SRT_YES)
	{
		strcat(szRealModelName, "2F");
	}

/* Creates a fake initial Term Structure with just one value for Sigma and Tau */
	iNumSigmaRows = 1;
	iNumSigmaCols = 1;
	ppdSigmaCurve = dmatrix(0, iNumSigmaCols-1, 0, iNumSigmaRows -1);
	ppdSigmaCurve[0][0] = 0.02;
	iNumTauRows = 1;
	iNumTauCols = 1;
	ppdTauCurve = dmatrix(0, iNumTauCols-1, 0, iNumTauRows -1);
	ppdTauCurve[0][0] = dTau;
	
/* Sets a default name for the underlying */
	sprintf(szUndName,"FlexiCap%sUnd",szModelName);

/* Creates an initial underlying with the input values for Tau, Alpha, Gamma and Rho */
	err = SrtInitIRUnd(
			szUndName, 
			szYieldCurveName, 
			szRealModelName, 
			iNumSigmaRows,
			iNumSigmaCols,
			ppdSigmaCurve, 
			iNumTauRows,
			iNumTauCols,
			ppdTauCurve,
			0.0, 
			dAlpha,
			dGamma,
			dRho,
			0.0,
			0.0,
			0.0, /* vasicek parms */
			0,
			0,
			NULL);
	
/* Free the memory allocated for the Sigma and Tau Cruves */
	free_dmatrix(ppdSigmaCurve  , 0, iNumSigmaCols-1, 0, iNumSigmaRows -1);
	free_dmatrix(ppdTauCurve  , 0, iNumTauCols-1, 0, iNumTauRows -1);

	if (err)
		return err;


/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(
				iNumGrfnParams,
				pszGrfnParamNames,
				pszGrfnParamValues,
				&sGrfnParams))
	{
		return err;
	}

/* Sets the Cap/Floor type */
	err = interp_barrier_type(szKnockUpOrDown,&eKnockUpOrDown);
	if (err)
		return err;
	if (eKnockUpOrDown == SRT_UP)
		sprintf(szCapFloor,"%s","CAP");
	else
		sprintf(szCapFloor,"%s","FLOOR");

/* Gets the Reference Rate information : basis, frequency */
	err = swp_f_get_ref_rate_details(szIndexRefRateCode, &eBasis, &eFrequency);
	if(err)
		return err;

/* Transforms the basis and frequency back into a string */
	translate_basis( &szBasis, eBasis);
	translate_compounding( &szFrequency, eFrequency);

/* Gets the Yield Curve (as well as today) */
	sCurve = lookup_curve(szYieldCurveName);
	lToday = get_today_from_curve(sCurve);

/* Select all the Fixing Dates on Today or after */
	iNumInstruments = lNumDates;
	while ( plIndexFixingDates[0] < lToday +  sGrfnParams.end_of_day_flg)
	{
		if (iNumInstruments <= 1)
		{
			return serror("All caplets have expired in SrtGrfnFlexiCapCalibrate");
		}
		else
		{
			iNumInstruments--;
			plIndexFixingDates ++;
			plIndexStartDates ++;
		}	
	}

/* Allocate memory for the calibration instruments (caplets in this case) */	
	plNfp = (long *) srt_malloc(iNumInstruments * sizeof(long));
	pszFrequency = svector_size(0, iNumInstruments-1,64);
	pszBasis = svector_size(0, iNumInstruments-1,64);
	pdStrike = (double *) srt_malloc(iNumInstruments * sizeof(double));
	pdBondStrike = (double *) srt_malloc(iNumInstruments * sizeof(double));
	pszRecPay = svector_size(0, iNumInstruments-1,64);
	pszRefRate = svector_size(0, iNumInstruments-1,64);
	pszProductType = svector_size(0, iNumInstruments-1,64);
	pdOptionPrice = (double *) srt_malloc(iNumInstruments * sizeof(double));
	pdATMOptionPrice = (double *) srt_malloc(iNumInstruments * sizeof(double));
	
/* Fills in the Fields for the instruments */
	strcpy(szPremium,"PREMIUM");
	for (i=0;i <iNumInstruments;i++)
	{
		forward = 0.0;
		plNfp[i] = 1;
		sprintf(pszFrequency[i], szFrequency );
		sprintf(pszBasis[i], szBasis) ;
		pdStrike[i] = pdKnockInLevels[i];
		pdBondStrike[i] = 1.0;
		sprintf(pszRecPay[i], szCapFloor );
		sprintf(pszRefRate[i], szIndexRefRateCode) ;
		sprintf(pszProductType[i], "CAPFLOOR") ;
		err = swp_f_CapFloor(plIndexStartDates[i], plNfp[i],pdStrike[i], pfGetVol,
					pszRecPay[i],pszRefRate[i], szYieldCurveName, 
					szPremium,szVolType, &pdOptionPrice[i]);
		err =  swp_f_ForwardRate(plIndexStartDates[i], plNfp[i],pszFrequency[i],pszBasis[i],
								szYieldCurveName,pszRefRate[i],&forward);
		err = swp_f_CapFloor(plIndexStartDates[i], plNfp[i],forward, pfGetVol,
					pszRecPay[i],pszRefRate[i], szYieldCurveName, szPremium,
					szVolType, &pdATMOptionPrice[i]);
	}
	
/* Calls the Calibration Routine */
	err = SrtBootstrap(
				iNumGrfnParams, 
				pszGrfnParamNames,
				pszGrfnParamValues,
				&iNumInstruments,
				plIndexStartDates, 
				plNfp,
				pszFrequency,
				pszBasis,
				pdStrike, 
				pdBondStrike, 
				pszProductType, 
				pszRecPay,       
				pszRefRate,
				pdOptionPrice,
				pdATMOptionPrice,
				dTau, 
				dAlpha,
				dGamma,
				dRho,
				szUndName);

/* Transfers the Caplets Prices information */
	*plNumCaplets = iNumInstruments;
	(*ppdCapletPrices) = dvector(0, *plNumCaplets -1);
	for (i=0;i <*plNumCaplets;i++)
		(*ppdCapletPrices)[i] = pdOptionPrice[i];

/* Free all */
	srt_free(plNfp);
	free_svector_size(pszFrequency , 0, iNumInstruments-1,64);
	free_svector_size(pszBasis , 0, iNumInstruments-1,64);
	srt_free(pdStrike);
	srt_free(pdBondStrike);
	free_svector_size(pszRecPay, 0, iNumInstruments-1,64);
	free_svector_size(pszRefRate, 0, iNumInstruments-1,64);
	free_svector_size(pszProductType , 0, iNumInstruments-1,64);
	srt_free(pdOptionPrice);
	srt_free(pdATMOptionPrice);

/* Return the error message if any */
	return err;
	
	
} /* Err Err SrtGrfnKnockInSwapCalibrate(...) */
