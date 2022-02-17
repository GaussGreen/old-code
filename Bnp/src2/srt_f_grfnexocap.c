/* ===============================================================================
   
   FILENAME:       srt_f_grfnexocap.c

   PURPOSE:        Pricing of a Exotics Cap with Grfn for a given model
                   Builds tableaux and prices deals

   =============================================================================== */

#include "grf_h_all.h"
#include "srt_h_all.h"
#include "grf_h_public.h"
#include "srt_h_grfnflexicap.h"
#include "srt_h_grfnexocap.h"
#include "srtaccess.h"


#define FREE_GRFN_EXOTIC_CAP_MEMORY {\
		if(*pppsGrfnTableau)\
			grfn_free_GrfnCellmatrix((*pppsGrfnTableau),(*plNumRows),(*plNumCols));\
		*pppsGrfnTableau = NULL;\
		if (*pplAuxRangesLength)\
			srt_free((*pplAuxRangesLength));\
		*pplAuxRangesLength = NULL;\
		free_dmatrix((*pppdAuxRanges) ,0, (*plNumAuxColumns) - 1, 0, lNumDates -1 );\
		*pppdAuxRanges = NULL;\
		free_dmatrix((*pppdLastPath) , 0, (*plNumRows) -1 , 0, (*plNumCols) -1 );\
		*pppdLastPath = NULL;\
	}



/* Cell in the first column of the chooser cap tableau : caplet */

static void GrfnChooserCapletCell(SrtReceiverType		eCapFloor,
								  double				dStrike,
								  String				szCashBasis,
								  String				szUndName,
								  String				szRefRateCode,
								  String				szCapletCell)
{

	if (eCapFloor == SRT_RECEIVER)
	{
		sprintf(szCapletCell,
			"max(%f - FRA(a[0,i],a[1,i],\"%s\",\"%s\",\"%s\"),0)*a[2,i]*df(now,a[1,i])",
			dStrike, 
			szCashBasis, 
			szUndName,
			szRefRateCode);
	}
	else
	if (eCapFloor == SRT_PAYER)
	{
		sprintf(szCapletCell,
			"max(FRA(a[0,i],a[1,i],\"%s\",\"%s\",\"%s\") - %f,0)*a[2,i]*df(now,a[1,i])",
			szCashBasis, 
			szUndName,
			szRefRateCode, 
			dStrike);
	}
} /* End of GrfnChooserCapletCell */


/* Cell in the last row of the chooser cap tableau : c[0,i] */

static void GrfnChooserEndCell(String	szEndCell)
{
	sprintf(szEndCell, "c[0,i]");
}


/* Cell used in Midat term */

static void GrfnChooserRowCell(String	szEndCell)
{
	sprintf(szEndCell, "c[0,i] + PV[j-1]");
}


/* Midat cell : if Startcolumn the max(PV[j],c[0,i]) otherwise max(PV[j],c[0,i]+PV[j-1]) */

static void GrfnChooserMidatCapCell(SRT_Boolean StartMidatColumn,
									String	szMidatCell)
{
char szTemp[GRFN_DEF_ARGBUFSZ];

	if (StartMidatColumn)
		GrfnChooserEndCell(szTemp);
	else
		GrfnChooserRowCell(szTemp);

	sprintf(szMidatCell, "max(PV[j],%s)", szTemp);
}

/* Barrier cell final cash flow: c[0,i] * c[1,i] (condition * caplet) */

static void GrfnBarrierCapCell(String	szFinalCashflowCell)
{
	sprintf(szFinalCashflowCell, "c[1,i] * c[2,i]");
}


static void GrfnBarrierCapletCell(SrtReceiverType   eCapFloor,
								  double            dStrike,
								  String			szCapletCell)
{
	if (eCapFloor == SRT_RECEIVER)
	{
		sprintf(szCapletCell,
			"max(%f - c[0,i], 0) * a[2,i] * df(now, a[1,i])",
			dStrike);
	}
	else
	if (eCapFloor == SRT_PAYER)
	{
		sprintf(szCapletCell,
			"max(c[0,i] - %f, 0) * a[2,i] * df(now, a[1,i])",
			dStrike);
	
	}
}

static void GrfnDigitalConditionCell(SRT_Boolean	bIsFirstCell,
									 double		dBarrier,
									 SRT_Boolean	bIsInOrOut, /* Out by default SRT_NO */
									 SRT_Boolean	bIsUpOrDown, /* Down by default SRT_NO */
								     String		szFinalCashflowCell)
{
	
	if(bIsUpOrDown) /* Up case */
		sprintf(szFinalCashflowCell,
				"if(c[0,i]>%f,%d,%d)",
				dBarrier,
				(bIsInOrOut==SRT_YES?1:0), 
				(bIsInOrOut==SRT_YES?0:1));
	else
		sprintf(szFinalCashflowCell,
				"if(c[0,i]<%f,%d,%d)",
				dBarrier,
				(bIsInOrOut==SRT_YES?1:0), 
				(bIsInOrOut==SRT_YES?0:1));

	if (bIsFirstCell==SRT_NO)
		sprintf(szFinalCashflowCell,
		"if(c[1,i-1]=%d,%d,%s)",
		(bIsInOrOut==SRT_YES?1:0),
		(bIsInOrOut==SRT_YES?1:0),
		szFinalCashflowCell);
}


/* Function that build the full grfn tableau for a Chooser cap */
static Err GrfnChooserMakeTableau(long              lNumEventDates,
								  double            dStrike,
								  long              lNumExercises,
								  SrtReceiverType   eCapFloor,
								  String			szCashBasis,
								  String			szUndName,
								  String			szRefRateCode,
								  long				*plNumRows,
								  long				*plNumCols,
								  GrfnCell			***pppsTableau)
{
Err 	err=NULL;
int     i, j;

	/* Sets the tableau dimensions : lNumExercices + 1  columns (c[0] -> N[lNumExercices] */
	
	*plNumCols 		= max(lNumExercises + 1, 2);
	*plNumRows 		= lNumEventDates;

	/* Allocate space for the Tableau (including strings) */	
	*pppsTableau 	= GrfnCellmatrix(*plNumRows,
									 *plNumCols,
									 GRFN_DEF_ARGBUFSZ);

	/* Fills in the tableau row by row */
	
	for (i=0; i < *plNumRows; i++)
	{
		/* The first column : the Caplet */
		GrfnChooserCapletCell(eCapFloor, 
							  dStrike, 
							  szCashBasis, 
							  szUndName,
							  szRefRateCode,
							  (*pppsTableau)[i][0].sval);
		(*pppsTableau)[i][0].type = GRFNSCELL;
		
		/* The second column: the first midat on caplet */
		
		if(i == *plNumRows - 1)
			GrfnChooserEndCell((*pppsTableau)[i][1].sval);
		else
			GrfnChooserMidatCapCell(SRT_YES, (*pppsTableau)[i][1].sval);		
		(*pppsTableau)[i][1].type = GRFNSCELL;
	
		/* The third column: the Midat term,
		   only if the number of column > 2 thanks Colin */
		if (*plNumCols > 2)
		{
			if(i == *plNumRows - 1)
				GrfnChooserEndCell((*pppsTableau)[i][2].sval);
			else
				GrfnChooserMidatCapCell(SRT_NO, (*pppsTableau)[i][2].sval);		
			(*pppsTableau)[i][2].type = GRFNSCELL;
		}
		/* Last columns are just copies of the third column */
		for(j=3; j<*plNumCols; j++)
		{
			sprintf((*pppsTableau)[i][j].sval, (*pppsTableau)[i][2].sval);
			(*pppsTableau)[i][j].type = GRFNSCELL;
		}
	} 
                      
	return err;
} /* End of GrfnChooserMakeTableau */


/* Function that build the full grfn tableau for a Barrier cap */
static Err GrfnBarrierMakeTableau(long              lNumEventDates,
								  double            dStrike,
								  double			dBarrier,
								  SrtReceiverType   eCapFloor,
								  SRT_Boolean			bIsInOrOut,
								  SRT_Boolean			bIsUpOrDown,
								  String			szCashBasis,
								  String			szUndName,
								  String			szRefRateCode,
								  long				*plNumRows,
								  long				*plNumCols,
								  GrfnCell			***pppsTableau)
{
Err 	err=NULL;
int     i;

	/* Sets the tableau dimensions : number of columns = 4 */
	
	*plNumCols 		= 5;
	*plNumRows 		= lNumEventDates;

	/* Allocate space for the Tableau (including strings) */	
	*pppsTableau 	= GrfnCellmatrix(*plNumRows,
									 *plNumCols,
									 GRFN_DEF_ARGBUFSZ);

	/* Fills in the tableau row by row */
	
	for (i=0; i < *plNumRows; i++)
	{
		/* The first column : c[0] - FRA */
		GrfnFlexiCapFraCell(szCashBasis, 
							szUndName,
							szRefRateCode,
							(*pppsTableau)[i][0].sval);
		(*pppsTableau)[i][0].type = GRFNSCELL;

		/* The second column : c[1] - Condition */
		GrfnDigitalConditionCell((i==0?SRT_TRUE:SRT_NO),
								 dBarrier,
								 bIsInOrOut,
								 bIsUpOrDown,
								 (*pppsTableau)[i][1].sval);
		(*pppsTableau)[i][1].type = GRFNSCELL;
		
		/* The third column: c[2] - Caplet at strike = dStrike */
		GrfnBarrierCapletCell(eCapFloor, 
							  dStrike, 
							  (*pppsTableau)[i][2].sval);
		(*pppsTableau)[i][2].type = GRFNSCELL;

		/* The fourth column: c[3] - Caplet at strike = dBarrier */
		GrfnBarrierCapletCell(eCapFloor, 
							  dBarrier, 
							  (*pppsTableau)[i][3].sval);
		(*pppsTableau)[i][3].type = GRFNSCELL;

		/* The last column: c[4] - final cashflow */
		GrfnBarrierCapCell((*pppsTableau)[i][4].sval);
		(*pppsTableau)[i][4].type = GRFNSCELL;
	} 
                      
	return err;
} /* End of GrfnBarrierMakeTableau */



/* the following function is based on the FlexiCap original code and price the ChooserCap */
Err SrtGrfnChooserPrice(long		*plFixingDates,
						long		*plStartDates,
						long		*plPayDates,
						double		*pdCoverages,
						long		lNumDates,
						String		szRefRateCode,
						double		dStrike,
						long		lMaxNumExercise,
						long		lNumAlreadyExercised,
						String		szCapFloor,
						String		szUndName,
						String		*pszGrfnParamNames,
						String		*pszGrfnParamValues,
						int			iNumGrfnParams,
						double		dFullCapRealPrice,
						/* Outputs from Grfn */
						double		*pdChooserPrice,
						double		*pdGrfnFullCapFloorPrice,
						long		*plNumRows,
						long		*plNumCols,
						GrfnCell 	***pppsGrfnTableau,
						double		***pppdLastPath,
						long		*plNumAuxColumns,
						long		**pplAuxRangesLength,
						double		***pppdAuxRanges)			
{
Err					err		= NULL;
SrtReceiverType		eCapFloor;
long				lNumEventdates	= 0;
SrtGrfnParam		sGrfnParams;
SrtUndPtr			sUndPtr;
SrtIOStruct			*psIOList;
SrtBasisCode		eBasis;
SrtCompounding		eFrequency;
String				szCashBasis;
double				*pdColumnPvs;
long				lNumPvs;
int					i;
long				lToday;

	/* Gets the Underlying from its name */
	sUndPtr = lookup_und(szUndName);
		
	/* Gets today from the Underlying */
	lToday = get_today_from_underlying(sUndPtr);

	/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(iNumGrfnParams,
								   pszGrfnParamNames,
								   pszGrfnParamValues,
								   &sGrfnParams))
		return err;

	/* Sets the Receiver/Payer type */
	if (err = interp_rec_pay(szCapFloor,&eCapFloor))
		return err;

	/* Gets the Reference Rate information : basis, frequency */
	if(err = swp_f_get_ref_rate_details(szRefRateCode, &eBasis, &eFrequency))
		return err;

	/* Transforms the basis back into a string */
	translate_basis( &szCashBasis, eBasis);

	/* Create a I/O list for the prices */
	if (err = srt_f_IOstructcreate(&psIOList,""))
		return(err);

	/* Remove everything that is before Today */
	while (	plFixingDates[0] < lToday + sGrfnParams.end_of_day_flg)
	{
		lNumDates --;
		plFixingDates ++;
		plStartDates ++;
		plPayDates ++;
		pdCoverages ++;	
	}

	/* Construct a Grfn Tableau that describes the Chooser  */
	(*pppsGrfnTableau)	= NULL;

	if (err = GrfnChooserMakeTableau(lNumDates,
									 dStrike,
									 lMaxNumExercise,
									 eCapFloor,
									 szCashBasis,
									 szUndName,
									 szRefRateCode,
									 plNumRows,
									 plNumCols,
									 pppsGrfnTableau))
	{
		srt_f_IOstructfree(&psIOList);
		return err;
	}

	/* Sets the Auxiliary ranges ; A[0] = start, A[1] = end, A[2] = coverage */
	(*plNumAuxColumns) = 3;
	(*pplAuxRangesLength) = (long *) srt_malloc ( (*plNumAuxColumns) * sizeof( long));
	(*pplAuxRangesLength)[0] = lNumDates;
	(*pplAuxRangesLength)[1] = lNumDates;
	(*pplAuxRangesLength)[2] = lNumDates;
	(*pppdAuxRanges) = dmatrix(0, (*plNumAuxColumns) - 1, 0, lNumDates -1 );
	for ( i =0; i < lNumDates; i++)
	{
		(*pppdAuxRanges)[0][i] = (double) plStartDates[i];
		(*pppdAuxRanges)[1][i] = (double) plPayDates[i];
		(*pppdAuxRanges)[2][i] = (double) pdCoverages[i];
	}


	/* Allocate space for the Grfn Cells (results of the last path) */
	(*pppdLastPath) = dmatrix(0, (*plNumRows) -1 , 0, (*plNumCols) -1 );

	/* Call GRFN to value the deal */
	if(err = srt_f_grfn(sUndPtr,
						&sGrfnParams,
						lNumDates,
						&plFixingDates,
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
						0))
	{
		FREE_GRFN_EXOTIC_CAP_MEMORY;
		return err;
	}

	/* Gets all the columns PV from the I/O list */
	if (err = srt_f_IOstructgetcolpvs(*psIOList, &pdColumnPvs, &lNumPvs))
	{
		FREE_GRFN_EXOTIC_CAP_MEMORY;
		return err;
	}

	/* The Flexi Cap PV is in the last column */
	*pdChooserPrice = pdColumnPvs[(*plNumCols) -1];
	
	/* The Full Cap PV is in the first column */
	*pdGrfnFullCapFloorPrice =  pdColumnPvs[0];

	/* Returns the success message */
	/* FREE_GRFN_EXOTIC_CAP_MEMORY; */
	return err;

} /* End of chooser price */

/* the following function is based on the FlexiCap original code and price the BarrierCap */
Err SrtGrfnBarrierPrice(long		*plFixingDates,
						long		*plStartDates,
						long		*plPayDates,
						double		*pdCoverages,
						long		lNumDates,
						String		szRefRateCode,
						double		dStrike,
						double		dBarrier,
						String		szCapFloor,
						SRT_Boolean		bIsInOrOut,
						SRT_Boolean		bIsUpOrDown,
						String		szUndName,
						String		*pszGrfnParamNames,
						String		*pszGrfnParamValues,
						int			iNumGrfnParams,
						/* Outputs from Grfn */
						double		*pdBarrierPrice,
						double		*pdGrfnFullCapFloorPrice,
						double		*pdGrfnCapFloorAtBarrier,
						long		*plNumRows,
						long		*plNumCols,
						GrfnCell 	***pppsGrfnTableau,
						double		***pppdLastPath,
						long		*plNumAuxColumns,
						long		**pplAuxRangesLength,
						double		***pppdAuxRanges)			
{
Err					err		= NULL;
SrtReceiverType		eCapFloor;
long				lNumEventdates	= 0;
SrtGrfnParam		sGrfnParams;
SrtUndPtr			sUndPtr;
SrtIOStruct			*psIOList;
SrtBasisCode		eBasis;
SrtCompounding		eFrequency;
String				szCashBasis;
double				*pdColumnPvs;
long				lNumPvs;
int					i;
long				lToday;

	/* Gets the Underlying from its name */
	sUndPtr = lookup_und(szUndName);
		
	/* Gets today from the Underlying */
	lToday = get_today_from_underlying(sUndPtr);

	/* Set and Overwrite defaults with user defined parameters */
	if (err = srt_f_set_GrfnParams(iNumGrfnParams,
								   pszGrfnParamNames,
								   pszGrfnParamValues,
								   &sGrfnParams))
		return err;

	/* Sets the Receiver/Payer type */
	if (err = interp_rec_pay(szCapFloor,&eCapFloor))
		return err;

	/* Gets the Reference Rate information : basis, frequency */
	if(err = swp_f_get_ref_rate_details(szRefRateCode, &eBasis, &eFrequency))
		return err;

	/* Transforms the basis back into a string */
	translate_basis( &szCashBasis, eBasis);

	/* Create a I/O list for the prices */
	if (err = srt_f_IOstructcreate(&psIOList,""))
		return(err);

	/* Remove everything that is before Today */
	while (	plFixingDates[0] < lToday + sGrfnParams.end_of_day_flg)
	{
		lNumDates --;
		plFixingDates ++;
		plStartDates ++;
		plPayDates ++;
		pdCoverages ++;	
	}

	/* Construct a Grfn Tableau that describes the Chooser  */
	(*pppsGrfnTableau)	= NULL;

	if (err = GrfnBarrierMakeTableau(lNumDates,
									 dStrike,
									 dBarrier,
									 eCapFloor,
									 bIsInOrOut,
									 bIsUpOrDown,
									 szCashBasis,
									 szUndName,
									 szRefRateCode,
									 plNumRows,
									 plNumCols,
									 pppsGrfnTableau))
	{
		srt_f_IOstructfree(&psIOList);
		return err;
	}

	/* Sets the Auxiliary ranges ; A[0] = start, A[1] = end, A[2] = coverage */
	(*plNumAuxColumns) = 3;
	(*pplAuxRangesLength) = (long *) srt_malloc ( (*plNumAuxColumns) * sizeof( long));
	(*pplAuxRangesLength)[0] = lNumDates;
	(*pplAuxRangesLength)[1] = lNumDates;
	(*pplAuxRangesLength)[2] = lNumDates;
	(*pppdAuxRanges) = dmatrix(0, (*plNumAuxColumns) - 1, 0, lNumDates -1 );
	
	for ( i =0; i < lNumDates; i++)
	{
		(*pppdAuxRanges)[0][i] = (double) plStartDates[i];
		(*pppdAuxRanges)[1][i] = (double) plPayDates[i];
		(*pppdAuxRanges)[2][i] = (double) pdCoverages[i];
	}


	/* Allocate space for the Grfn Cells (results of the last path) */
	(*pppdLastPath) = dmatrix(0, (*plNumRows) -1 , 0, (*plNumCols) -1 );

	/* Call GRFN to value the deal */
	if(err = srt_f_grfn(sUndPtr,
						&sGrfnParams,
						lNumDates,
						&plFixingDates,
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
						0))
	{
		FREE_GRFN_EXOTIC_CAP_MEMORY;
		return err;
	}

	/* Gets all the columns PV from the I/O list */
	if (err = srt_f_IOstructgetcolpvs(*psIOList, &pdColumnPvs, &lNumPvs))
	{
		FREE_GRFN_EXOTIC_CAP_MEMORY;
		return err;
	}

	/* The Flexi Cap PV is in the last column  */
	*pdBarrierPrice = pdColumnPvs[4];/*(*plNumCols) -1];*/
	
	/* The Full Cap PV is in the second column */
	*pdGrfnFullCapFloorPrice =  pdColumnPvs[2];
	/* The Full Cap PV at strike = barrier is in the third column */
	*pdGrfnCapFloorAtBarrier =  pdColumnPvs[3];

	/* Returns the success message */
	FREE_GRFN_EXOTIC_CAP_MEMORY;
	return err;

} /* End of SrtGrfnBarrierPrice */


void SrtGRFNQuickExoticCapGrfnOptions (String	**grfn_param,
									   String	**grfn_value,
									   long		*lNumParams,
									   SRT_Boolean	UseTree)
{
	/*	No memory error assumed */
	if(UseTree)
		*lNumParams = 2;
	else
		*lNumParams = 6;

	*grfn_param = svector_size (0, *lNumParams, 32);
	*grfn_value = svector_size (0, *lNumParams, 32);

	strcpy ((*grfn_param)[0], "FORCEMC");

	if(UseTree)
		strcpy ((*grfn_value)[0], "NO");
	else
		strcpy ((*grfn_value)[0], "YES");

	strcpy ((*grfn_param)[1], "MAXTIME");
	strcpy ((*grfn_value)[1], "0.08");

	if(UseTree)
	{
		strcpy ((*grfn_param)[2], "MINNODE");
		strcpy ((*grfn_value)[2], "400");
	}
	else
	{
		strcpy ((*grfn_param)[2], "SAMPLETYPE");
		strcpy ((*grfn_value)[2], "ABS");
		strcpy ((*grfn_param)[3], "NUMPATH");
		strcpy ((*grfn_value)[3], "5000");
		strcpy ((*grfn_param)[4], "RENORM");
		strcpy ((*grfn_value)[4], "YES");
		strcpy ((*grfn_param)[5], "RANDSEED");
		strcpy ((*grfn_value)[5], "-987564321");
		strcpy ((*grfn_param)[6], "JUMPING");
		strcpy ((*grfn_value)[6], "YES");
	}
} /* End of SrtGRFNQuickExoticCapGrfnOptions */



/* Quick Price of the Flexi/Chooser using MAD addins */

Err SrtGRFNQuickExoticCapPrice(long		StartDate,
							   long		EndDate,
							   long		lNumExercices,
							   String	szRefRateCode,
							   double	dStrike,
							   double	dEps,
							   String	szCapFloor,
							   String	szYieldCurveName,
							   Err		(*pfGetVol)(double dStart, double dEnd, double dStrike, 
													double dForward, double dSpread, double *pdBsVol),
							   String	szVolType,
							   double	dTau,
							   String	szReturnUnd,
							   double	dAlpha,
							   double	dGamma,
							   double	dRho,
							   SRT_Boolean	bIsChooser,
							   /* Output of the function */
							   String	szUndName,
							   double	*pdGrfnFlexiCapPrice,
							   double	*pdGrfnFullCapFloorPrice,
							   SRT_Boolean	*pbReturnUnd)
{
Err					err=NULL;
SrtBasisCode		SrtFloatBasis;
SrtCompounding		SrtFloatFrequency;
SrtCurvePtr			Crv;
SwapDP				*FloatLegDP;
SrtCcyParam			*ccy_param;
GrfnCell 			**ppsGrfnTableau;
SrtUndListPtr		pSrtUndList;
int					SpotLag,
					NumDates,
					NumPayDates;

Date				*FixingDates, 
					*StartDates, 
					*EndDates, 
					*PayDates,
					Today, 
					SpotDate;

String				FloatBasis,
					FloatFrequency,
					*pszGrfnParamNames,
					*pszGrfnParamValues;

double				*Coverage,
					*pdCapletRealPrices,
					**ppdLastPath,
					dFullCapRealPrice=0,
					**ppdAuxRanges;

long					lNumCaplets,
					lNumRows,
					lNumCols,
					lNumAuxColumns,
					*plAuxRangesLength,
					dNumberOfGrfnParams = 0;	

SRT_Boolean				bUseTwoFactor = SRT_NO, /* One factor by default */
					bReturnUnd = SRT_NO, /* return prices by default */
					UndInList = SRT_NO, /* No existing underlying */
					bNotFullCap = SRT_YES; /* By default the exotic */
										/*	product is not a full cap */

	if (err = swp_f_get_ref_rate_details(szRefRateCode, &SrtFloatBasis, &SrtFloatFrequency))
		return err;

	if (err = translate_basis(&FloatBasis,
							  SrtFloatBasis))
		return err;

	if (err = translate_compounding(&FloatFrequency,
									SrtFloatFrequency))
		return err;

	/* Allocations Sorry No memory error assumed */

	FloatLegDP = (SwapDP*) malloc(sizeof(SwapDP));
	
	/* Get Classics informations from the Yield Curve Name */
	
	Crv = lookup_curve(szYieldCurveName);
	Today = get_clcndate_from_yldcrv(Crv);
	SpotDate = get_spotdate_from_yldcrv(Crv);
	ccy_param = get_ccyparam_from_yldcrv(Crv);
		
	SpotLag = ccy_param->spot_lag;

	if(StartDate < SpotDate)  
		StartDate = SpotDate;

	if (err = swp_f_setSwapDP(StartDate, 
							  EndDate, 
							  SrtFloatFrequency, 
							  SrtFloatBasis, 
							  FloatLegDP))
		return err;
	
	FloatLegDP -> spot_lag = SpotLag;


	if (err = swp_f_make_FloatLegDatesAndCoverages(FloatLegDP, 
												   Today, 
												   &PayDates, 
												   &NumPayDates,
												   &FixingDates, 
												   &StartDates,
												   &EndDates, 
												   &Coverage, 
												   &NumDates))
		return err;


	/* Setup of GRFN options */


	/*	First Automatically detect if Two Factors case : if alpha is different from 0 */

	if(dAlpha != 0.0)
		bUseTwoFactor = SRT_YES;

	strupper(szReturnUnd);

	/* Get rid of the ticker for the underlying name */
	rem_tick_string(szReturnUnd, szReturnUnd);

	/* From string to SRT_Boolean for return underlying: */

	pSrtUndList = get_underlying_list();
	
	UndInList = srt_f_isundinlist(pSrtUndList, szReturnUnd);
		
	if ((!UndInList) && (!strncmp(szReturnUnd, "Y", 1)))
	{
		/* NO underlying in parameter and no ask for returning price 
			-> return the calibrated underlying */
		bReturnUnd = SRT_YES;
	}

	/* Flexi -> Monte Carlo and Chooser -> Tree */

	if (bIsChooser)
		SrtGRFNQuickExoticCapGrfnOptions(&pszGrfnParamNames, 
										 &pszGrfnParamValues,
										 &dNumberOfGrfnParams,
										 SRT_YES);
	else
		SrtGRFNQuickExoticCapGrfnOptions(&pszGrfnParamNames, 
										 &pszGrfnParamValues,
										 &dNumberOfGrfnParams,
										 SRT_NO);


	/* to compute FullCapRealPrice... */

	if (err = swp_f_CapFloor(StartDate, 
							 EndDate,
							 dStrike, 
							 pfGetVol,
							 szCapFloor,
							 szRefRateCode, 
							 szYieldCurveName, 
							 "PREMIUM",
							 szVolType, 
							 &dFullCapRealPrice))
		return err;

	/* Calibration of the underlying */

	if (lNumExercices < NumDates)
	{
		if (!UndInList)
		{
			/* let's give a name to the underlying */
			if (bIsChooser)
				strcpy(szUndName, "Chooser");
			else
				strcpy(szUndName, "Flexi");

			if (err=SrtGrfnFlexiCapCalibrate(FixingDates,
											 StartDates,
											 NumDates,
											 szRefRateCode,
											 dStrike,
											 szCapFloor,
											 "LGM",
											 dTau,
											 bUseTwoFactor,
											 dAlpha,
											 dGamma,
											 dRho,
											 pszGrfnParamNames,
											 pszGrfnParamValues,
											 dNumberOfGrfnParams, 
											 szYieldCurveName,
											 pfGetVol,
											 szVolType,
											 szUndName,
											 &pdCapletRealPrices,
											 &lNumCaplets))
				return err;
		}
		else
		{
			strcpy(szUndName, szReturnUnd);
		}
	}
	else
	/* switch case if the number of exercice of the flexi or chooser 
									>= 
		number of caplets then just return the market capfloor value */
	{
		bReturnUnd = SRT_YES;
		bNotFullCap = SRT_NO;
		*pdGrfnFlexiCapPrice = dFullCapRealPrice;
		*pdGrfnFullCapFloorPrice = dFullCapRealPrice;
	}	
	
	
	/* Compute the Flexi/Chooser price using the Calibrated underlying only 
	   If the user ask for it... */

	if (!bReturnUnd)
	{
		if (bIsChooser)
		{
			if(err= SrtGrfnChooserPrice(FixingDates,
										StartDates,
										EndDates,
										Coverage,
										NumDates,
										szRefRateCode,
										dStrike,
										lNumExercices,
										0,
										szCapFloor,
										szUndName,
										pszGrfnParamNames,
										pszGrfnParamValues,
										dNumberOfGrfnParams, 
										dFullCapRealPrice,
										pdGrfnFlexiCapPrice,
										pdGrfnFullCapFloorPrice,
										&lNumRows,
										&lNumCols,
										&ppsGrfnTableau,
										&ppdLastPath,
										&lNumAuxColumns,
										&plAuxRangesLength,
										&ppdAuxRanges))
			return err;

		}
		else
		{
			if(err = SrtGrfnFlexiCapPrice(FixingDates,
										StartDates,
										EndDates,
										Coverage,
										NumDates,
										szRefRateCode,
										dStrike,
										dEps,
										lNumExercices,
										0,
										szCapFloor,
										szUndName,
										pszGrfnParamNames,
										pszGrfnParamValues,
										dNumberOfGrfnParams, 
										dFullCapRealPrice,
										pdGrfnFlexiCapPrice,
										pdGrfnFullCapFloorPrice,
										&lNumRows,
										&lNumCols,
										&ppsGrfnTableau,
										&ppdLastPath,
										&lNumAuxColumns,
										&plAuxRangesLength,
										&ppdAuxRanges))
			return err;
		}
	
	/* Free all of them */

		if(ppsGrfnTableau)
			grfn_free_GrfnCellmatrix(ppsGrfnTableau, lNumRows, lNumCols);

		ppsGrfnTableau = NULL;
		
		if (plAuxRangesLength)
			srt_free(plAuxRangesLength);

		plAuxRangesLength = NULL;

		if (ppdAuxRanges)
			free_dmatrix(ppdAuxRanges ,0, lNumAuxColumns - 1, 0, NumDates - 1 );
	
		ppdAuxRanges = NULL;
		
		if (ppdLastPath)
			free_dmatrix(ppdLastPath , 0, lNumRows -1 , 0, lNumCols -1 );
		ppdLastPath = NULL;
	}

	*pbReturnUnd = bReturnUnd && bNotFullCap;

	free_svector_size(pszGrfnParamNames, 0, dNumberOfGrfnParams, 32);
	free_svector_size(pszGrfnParamValues, 0, dNumberOfGrfnParams, 32);

	free(FixingDates);
	free(StartDates);
	free(EndDates);
	free(PayDates); 
	free(Coverage); 

	free(FloatLegDP);

	return err;

} /* SrtGRFNQuickExoticCapPrice */

Err SrtGRFNQuickBarrierCapPrice(long		StartDate,
								long		EndDate,
								String		szRefRateCode,
								double		dStrike,
								double		dBarrier,
								String		szCapFloor,
								String		szUpDown,
								String		szInOut,
								String		szYieldCurveName,
								Err			(*pfGetVol)(double dStart, double dEnd, double dStrike, 
														double dForward, double dSpread, double *pdBsVol),
								String		szVolType,
								double		dTau,
								String		szReturnUnd,
								double		dAlpha,
								double		dGamma,
								double		dRho,
								/* Output of the function */
								String		szUndName,
								double		*pdGrfnBarrierPrice,
								double		*pdGrfnFullCapFloorPrice,
								double		*pdGrfnCapFloorAtBarrier,
								SRT_Boolean		*pbReturnUnd)
{
Err					err = NULL;
SrtBasisCode		SrtFloatBasis;
SrtCompounding		SrtFloatFrequency;
SrtCurvePtr			Crv;
SwapDP				*FloatLegDP;
SrtCcyParam			*ccy_param;
GrfnCell 			**ppsGrfnTableau;
SrtUndListPtr			pSrtUndList;
int				NumDates,
				NumPayDates, 
				SpotLag;

Date				*FixingDates, 
					*StartDates, 
					*EndDates, 
					*PayDates,
					Today, 
					SpotDate;

String				*pszGrfnParamNames,
					*pszGrfnParamValues,
					FloatBasis,
					FloatFrequency;

double				*pdCapletRealPrices,
					**ppdLastPath,
					**ppdAuxRanges,
					*Coverage;

long					lNumCaplets,
					lNumRows,
					lNumCols,
					lNumAuxColumns,
					*plAuxRangesLength,
					dNumberOfGrfnParams = 0;	

SRT_Boolean				bUseTwoFactor = SRT_NO, /* One factor by default */
					bReturnUnd = SRT_NO, /* Return prices by default */
					UndInList = SRT_NO, /* No existing underlying */
					bIsUpOrDown = SRT_NO, /* Down by default */
					bIsInOrOut = SRT_NO; /* Out by default */


	if (err = swp_f_get_ref_rate_details(szRefRateCode, &SrtFloatBasis, &SrtFloatFrequency))
		return err;

	if (err = translate_basis(&FloatBasis, SrtFloatBasis))
		return err;

	if (err = translate_compounding(&FloatFrequency, SrtFloatFrequency))
		return err;

	/* Allocations Sorry No memory error assumed */

	FloatLegDP = (SwapDP*) malloc(sizeof(SwapDP));
	
	/* Get Classics informations from the Yield Curve Name */
	
	Crv = lookup_curve(szYieldCurveName);
	Today = get_clcndate_from_yldcrv(Crv);
	SpotDate = get_spotdate_from_yldcrv(Crv);
	ccy_param = get_ccyparam_from_yldcrv(Crv);
		
	SpotLag = ccy_param->spot_lag;

	if(StartDate < SpotDate)  
		StartDate = SpotDate;

	if (err = swp_f_setSwapDP(StartDate, 
							  EndDate, 
							  SrtFloatFrequency, 
							  SrtFloatBasis, 
							  FloatLegDP))
		return err;
	
	FloatLegDP -> spot_lag = SpotLag;


	if (err = swp_f_make_FloatLegDatesAndCoverages(FloatLegDP, 
												   Today, 
												   &PayDates, 
												   &NumPayDates,
												   &FixingDates, 
												   &StartDates,
												   &EndDates, 
												   &Coverage, 
												   &NumDates))
		return err;

	/*	First Automatically detect if Two Factors case : if alpha is different from 0 */

	if(dAlpha != 0.0)
		bUseTwoFactor = SRT_YES;

	strupper(szReturnUnd);

	/* Get rid of the ticker for the underlying name */
	rem_tick_string(szReturnUnd, szReturnUnd);

	/* From string to SRT_Boolean for return underlying: */

	pSrtUndList = get_underlying_list();
	
	UndInList = srt_f_isundinlist(pSrtUndList, szReturnUnd);
		
	if ((!UndInList) && (!strncmp(szReturnUnd, "Y", 1)))
	{
		/* NO underlying in parameter and no ask for returning price 
		-> return the calibrated underlying */
		bReturnUnd = SRT_YES;
	}

	/* Barrier -> Monte Carlo */
	
	SrtGRFNQuickExoticCapGrfnOptions(&pszGrfnParamNames, 
									&pszGrfnParamValues,
									&dNumberOfGrfnParams,
									SRT_NO);
	/*
		SrtGRFNQuickExoticCapGrfnOptions(&pszGrfnParamNames, 
										 &pszGrfnParamValues,
										 &dNumberOfGrfnParams,
										 SRT_NO);

		here follows the transformation from strings to booleans for in/out/up/down ... 
	*/

	strupper(szUpDown);
	strupper(szInOut);

	if(!strncmp(szUpDown, "U", 1))
		bIsUpOrDown = SRT_YES;

	if(!strncmp(szInOut, "I", 1))
		bIsInOrOut = SRT_YES;


	/* Calibration of the underlying */

	if (!UndInList)
	{
		/* let's give a name to the underlying: 
			If 'In then it's Light,
			if 'Out then Extinguish.*/
		if (bIsInOrOut)
			strcpy(szUndName, "Light");
		else
			strcpy(szUndName, "Extinguish");

		if (err=SrtGrfnFlexiCapCalibrate(FixingDates,
										 StartDates,
										 (long) NumDates,
										 szRefRateCode,
										 dStrike,
										 szCapFloor,
										 "LGM",
										 dTau,
										 bUseTwoFactor,
										 dAlpha,
										 dGamma,
										 dRho,
										 pszGrfnParamNames,
										 pszGrfnParamValues,
										 dNumberOfGrfnParams, 
										 szYieldCurveName,
										 pfGetVol,
										 szVolType,
										 szUndName,
										 &pdCapletRealPrices,
										 &lNumCaplets))
			return err;
	}
	else
	{
		strcpy(szUndName, szReturnUnd);
	}
	
	
	/* Compute the Barrier price using the Calibrated underlying only 
	   If the user ask for it... */

	if (!bReturnUnd)
	{
		if(err= SrtGrfnBarrierPrice(FixingDates,
									StartDates,
									EndDates,
									Coverage,
									(long) NumDates,
									szRefRateCode,
									dStrike,
									dBarrier,
									szCapFloor,
									bIsInOrOut,
									bIsUpOrDown,
									szUndName,
									pszGrfnParamNames,
									pszGrfnParamValues,
									dNumberOfGrfnParams, 
									pdGrfnBarrierPrice,
									pdGrfnFullCapFloorPrice,
									pdGrfnCapFloorAtBarrier,
									&lNumRows,
									&lNumCols,
									&ppsGrfnTableau,
									&ppdLastPath,
									&lNumAuxColumns,
									&plAuxRangesLength,
									&ppdAuxRanges))
		return err;
	
	/* Free them all */

		if(ppsGrfnTableau)
			grfn_free_GrfnCellmatrix(ppsGrfnTableau, lNumRows, lNumCols);

		ppsGrfnTableau = NULL;
		
		if (plAuxRangesLength)
			srt_free(plAuxRangesLength);

		plAuxRangesLength = NULL;

		if (ppdAuxRanges)
			free_dmatrix(ppdAuxRanges ,0, lNumAuxColumns - 1, 0, NumDates - 1 );
	
		ppdAuxRanges = NULL;
		
		if (ppdLastPath)
			free_dmatrix(ppdLastPath , 0, lNumRows -1 , 0, lNumCols -1 );
		ppdLastPath = NULL;
	}

	*pbReturnUnd = bReturnUnd;

	free_svector_size(pszGrfnParamNames, 0, dNumberOfGrfnParams, 32);
	free_svector_size(pszGrfnParamValues, 0, dNumberOfGrfnParams, 32);

	free(FixingDates);
	free(StartDates);
	free(EndDates);
	free(PayDates); 
	free(Coverage); 

	free(FloatLegDP);

	return err;

} /* SrtGRFNQuickBarrierCapPrice */


Err GrfnQuickDisplayTableau(String		szProductName,
							long		StartDate,
							long		EndDate,
							long		lNumExercises,
							String		szRefRateCode,
							double		dStrike,
							double		dBarrier,
							String		szCapFloor,
							String		szUpDown,
							String		szInOut,
							String		szYieldCurveName,
							long		*plNumRows,
							long		*plNumCols,
							GrfnCell 	***pppsGrfnTableau)
{
Err					err = NULL;
SrtBasisCode		SrtFloatBasis;
SrtCompounding		SrtFloatFrequency;
SrtCurvePtr			Crv;
SwapDP				*FloatLegDP;
SrtCcyParam			*ccy_param;
SrtReceiverType		eCapFloor;

int				NumPayDates,
				NumDates,
				SpotLag;

Date				*FixingDates, 
					*StartDates, 
					*EndDates, 
					*PayDates,
					Today, 
					SpotDate;

String				FloatBasis;
	
double				*Coverage;

SRT_Boolean				bIsUpOrDown = SRT_NO, /* Down by default */
					bIsInOrOut = SRT_NO; /* Out by default */

char				szUndName[128];

	if (err = swp_f_get_ref_rate_details(szRefRateCode, &SrtFloatBasis, &SrtFloatFrequency))
		return err;

	if (err = translate_basis(&FloatBasis, SrtFloatBasis))
		return err;

	if (err = interp_rec_pay(szCapFloor,&eCapFloor))
		return err;

	/* Allocations Sorry No memory error assumed */

	FloatLegDP = (SwapDP*) malloc(sizeof(SwapDP));
	
	/* Get Classics informations from the Yield Curve Name */
	
	Crv = lookup_curve(szYieldCurveName);
	Today = get_clcndate_from_yldcrv(Crv);
	SpotDate = get_spotdate_from_yldcrv(Crv);
	ccy_param = get_ccyparam_from_yldcrv(Crv);
		
	SpotLag = ccy_param->spot_lag;

	if(StartDate < SpotDate)  
		StartDate = SpotDate;

	if (err = swp_f_setSwapDP(StartDate, 
							  EndDate, 
							  SrtFloatFrequency, 
							  SrtFloatBasis, 
							  FloatLegDP))
		return err;
	
	FloatLegDP -> spot_lag = SpotLag;


	if (err = swp_f_make_FloatLegDatesAndCoverages(FloatLegDP, 
												   Today, 
												   &PayDates, 
												   &NumPayDates,
												   &FixingDates, 
												   &StartDates,
												   &EndDates, 
												   &Coverage, 
												   &NumDates))
		return err;

	/* free everything 'cause we only need the number of dates */
	free(FixingDates);
	free(StartDates);
	free(EndDates);
	free(PayDates); 
	free(Coverage); 
	free(FloatLegDP);


	strupper(szProductName);

	if(!strncmp(szProductName, "B", 1))
	{
		strcpy(szUndName, "BarrierUnd");

		strupper(szUpDown);
		strupper(szInOut);

		if(!strncmp(szUpDown, "U", 1))
			bIsUpOrDown = SRT_YES;

		if(!strncmp(szInOut, "I", 1))
			bIsInOrOut = SRT_YES;

		err = GrfnBarrierMakeTableau(			   (long)NumDates,
								   dStrike,
								   dBarrier,
								   eCapFloor,
								   bIsInOrOut,
								   bIsUpOrDown,
								   FloatBasis,
								   szUndName,
								   szRefRateCode,
								   plNumRows,
								   plNumCols,
								   pppsGrfnTableau);
	}

	if(!strncmp(szProductName, "C", 1))
	{
		strcpy(szUndName, "ChooserUnd");

		err = GrfnChooserMakeTableau(NumDates,
									dStrike,
									lNumExercises,
									eCapFloor,
									FloatBasis,
									szUndName,
									szRefRateCode,
									plNumRows,
									plNumCols,
									pppsGrfnTableau);
	}

	if(!strncmp(szProductName, "F", 1))
	{

		err = GrfnFlexiCapMakeTableau(NumDates,
									dStrike,
									lNumExercises,
									0,
									eCapFloor,
									FloatBasis,
									szUndName,
									szRefRateCode,
									plNumRows,
									plNumCols,
									pppsGrfnTableau);

	}
	return err;
}	
