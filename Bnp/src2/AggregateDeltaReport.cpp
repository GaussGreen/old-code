#include "srt_h_all.h"
#include "AggregateDeltaReport.h"
#include "UTALLHDR.h"

#define SWAP(a,b)  {long tempr; tempr = (a); (a) = (b); (b) = tempr;} 


void   srt_f_UpdateGlobalHedgeUndFraTermStr(Date		*GlobalHedgeUndFraEndDates,
											Date		*LocalHedgeUndFraEndDates,
											long		NumLocalHedgeUndFra,
											long		NumGlobalHedgeUndFra,
											double		*LocalHedgeUndFraPEQ,
											double		**GlobalHedgeUndFraPEQ)
{
long	k,j;

	for(k=0;k<NumLocalHedgeUndFra;k++)
	{
		j=0;
		while((j<NumGlobalHedgeUndFra)&&(GlobalHedgeUndFraEndDates[j]!=LocalHedgeUndFraEndDates[k])) j++;
		(*GlobalHedgeUndFraPEQ)[j]+=LocalHedgeUndFraPEQ[k];
	}

}



void 	srt_f_UpDateHedgeTermStr(Date		 *HedgeTenorDates,
								 long		 NumHedgeTenors,
								 Date		 *SwapHedgeTenorDates,
								 long		 NumSwpHedgeTenors,
								 double		 *SwapHedgePEQ,
								 long		 PrevTenorIndex,
								 long		 NextTenorIndex,

								 double		 **HedgePEQ)
{
long	j;
	
	j=0;
	while(HedgeTenorDates[j]!=SwapHedgeTenorDates[PrevTenorIndex]) j++;
	(*HedgePEQ)[j]+=SwapHedgePEQ[PrevTenorIndex];

	if(PrevTenorIndex!=NextTenorIndex)
	{
		j=0;
		while(HedgeTenorDates[j]!=SwapHedgeTenorDates[NextTenorIndex]) j++;
		(*HedgePEQ)[j]+=SwapHedgePEQ[NextTenorIndex];
	}
}


void srt_f_UpdateLocHedgeUndFraTermStr(Date				*UndFraHedgeEndDates,
									   double			*tmp_HedgeUndFraPEQ,
									   long				PrevTenorIndex,
									   long				NextTenorIndex,
									   SrtUndFRAType    UndFraType,

									   /*OUTPUT */
									   long				*NumLocalHedgeUndFra,
									   double			**LocalHedgeUndFraPEQ,
									   Date				**LocalHedgeUndFraEndDates,
									   SrtUndFRAType	**LocalHedgeUndFraTypes)
{

	/*INITIALISATION */

	(*NumLocalHedgeUndFra)=1;
	
	(*LocalHedgeUndFraPEQ)[0] = tmp_HedgeUndFraPEQ[PrevTenorIndex];
	(*LocalHedgeUndFraEndDates)[0]= UndFraHedgeEndDates[PrevTenorIndex];
	(*LocalHedgeUndFraTypes)[0]= UndFraType;

	if(PrevTenorIndex!=NextTenorIndex)
	{
		(*NumLocalHedgeUndFra)=2;
		(*LocalHedgeUndFraPEQ) = realloc((*LocalHedgeUndFraPEQ),(*NumLocalHedgeUndFra)*sizeof(double));
		(*LocalHedgeUndFraEndDates) = realloc((*LocalHedgeUndFraEndDates),(*NumLocalHedgeUndFra)*sizeof(Date));
		(*LocalHedgeUndFraTypes) = realloc((*LocalHedgeUndFraTypes),(*NumLocalHedgeUndFra)*sizeof(SrtUndFRAType));
		
		(*LocalHedgeUndFraPEQ)[1] = tmp_HedgeUndFraPEQ[NextTenorIndex];
		(*LocalHedgeUndFraEndDates)[1]= UndFraHedgeEndDates[NextTenorIndex];
		(*LocalHedgeUndFraTypes)[1]= UndFraType;
	}

}	



Err	 UpdateUndFRATermStruct(Date				DateToInsert,
							Date				**SetOfRefDates,
							String				UndFRATypeStr,
							SrtUndFRAType		**SetOfRefUndFRATypes,
							long				*NumUndFra)
{

Err		err=NULL;
long	i,j;
long    jmax;
SrtUndFRAType		UndFraType;

		jmax = (*NumUndFra);
		j = 0;
		while((jmax>0)&&(DateToInsert!=(*SetOfRefDates)[j])&&(j<jmax)) j++;

		/*INSERT THE DATES */
		if(j==jmax)
		{
			err = interp_UndFRAType(UndFRATypeStr,&UndFraType);
			if(err)
				return serror("Fatal Error in interp_UndFRAType");

			(*NumUndFra)++;
			jmax++;
			(*SetOfRefDates)		= (Date*)realloc((*SetOfRefDates),(*NumUndFra)*sizeof(Date));
			(*SetOfRefUndFRATypes)	= (SrtUndFRAType*)realloc((*SetOfRefUndFRATypes),(*NumUndFra)*sizeof(SrtUndFRAType));

			(*SetOfRefDates)[jmax-1]	     = DateToInsert; 
			(*SetOfRefUndFRATypes)[jmax-1] = UndFraType;

			/*SORT BY DATES */
			for(j=0;j<jmax;j++)
			{
				for(i=j+1;i<jmax;i++)
				{
					if((*SetOfRefDates)[i]<(*SetOfRefDates)[j])
					{
						SWAP((*SetOfRefDates)[i],(*SetOfRefDates)[j]);
						SWAP((*SetOfRefUndFRATypes)[i],(*SetOfRefUndFRATypes)[j]);
					}
				}
			}

		}

		return err;
}


Err GetTenorsPEQ(String		RefTenors,
				 long		NumMarketTenors,
				 String		*MarketTenors,
				 double		*MarketTenorsPEQ,
				 
				 /*OUTPUT*/
				 double		*RefTenorPEQ)
{
Err			err=NULL;
long		i;

	/*CHECK OF THE INPUT */
	if(!MarketTenors)
		return serror("Fatal:No Market Tenors in GetTenorsPEQ");

	if(!MarketTenorsPEQ)
		return serror("Fatal:No Market Tenors in GetTenorsPEQ");

	i = 0;
	while((i<NumMarketTenors)&&(strcmp(MarketTenors[i],RefTenors)!=0)) i++;

	/*TRUE = 0  FOR strcmp */
	*RefTenorPEQ = MarketTenorsPEQ[i];

	return err;
}


SRT_Boolean ExitsHedgeTenors(Date		MarketTenorDate, /*THE SWAP TO HEDGE */
						long		NumHedgeTenors,
						Date		*HedgeTenorDates)

{

SRT_Boolean ExitsNext = SRT_NO;
SRT_Boolean ExitsPrevious = SRT_NO;
SRT_Boolean Exits = SRT_NO;
long	i;


	for(i=0;i<NumHedgeTenors;i++)
	{
		if(HedgeTenorDates[i]>=MarketTenorDate) ExitsNext = SRT_YES;
		if(HedgeTenorDates[i]<=MarketTenorDate) ExitsPrevious = SRT_YES;
	}

	if((ExitsNext==SRT_YES)&&(ExitsPrevious==SRT_YES)) Exits = SRT_YES;
	
	return Exits;

}
	
Err FindHedges(Date			MarketTenorDate,
			   long			NumHedgeTenors,
			   Date			*HedgeTenorDates,
			   long			*PrevTenorIndex,
			   long			*NextTenorIndex)
{

Err		err = NULL;
long	i;
	
	/*GET THE NextSwapIndex */
	i=0;
	while((i<NumHedgeTenors)&&(HedgeTenorDates[i]<MarketTenorDate))	i++;
	*NextTenorIndex = min(i,NumHedgeTenors-1);

	/*GET THE PrevSwapIndex */
	while((i>0)&&(HedgeTenorDates[i]>MarketTenorDate))	i--;
	*PrevTenorIndex = min(i,NumHedgeTenors-1);
	
	return err;
}

/*LOCAL AGGREGATION OF FRA TENORS */	

Err	srt_f_AggregateFraTenor(Date			FraEndDate,
							double			FraPEQ,
							long			NumHedgeSwpUndFra,
							long			NumHedgeFutUndFra,
							Date			*SwpUndFraHedgeEndDates,
							Date			*FutUndFraHedgeEndDates,
							SrtCompounding	FloatFreq,
							SrtBasisCode	FloatBasis,
							String			FloatRefRate,
							String			YieldCurveName,
	
							/*OUTPUTS */
							long			*NumLocalHedgeUndFra,
							double			**LocalHedgeUndFraPEQ,
							Date			**LocalHedgeUndFraEndDates,
							SrtUndFRAType	**LocalHedgeUndFraTypes)
							
{
Err			err=NULL;
long		PrevTenorIndex,NextTenorIndex;
double	    *tmp_LocalHedgeUndFraPEQ=NULL;
SrtCrvPtr	YieldCrvPtr;
long		SpotLag;

	/*CHECK OF THE INPUTS */
	if((NumHedgeSwpUndFra==0)&&(NumHedgeFutUndFra==0))
		return serror("Hedge Tenors Underlying FRA  are missing ");

	YieldCrvPtr = lookup_curve(YieldCurveName);	
	if(!YieldCrvPtr)
		return serror("Can not get the yield curve");
	SpotLag = get_spotlag_from_curve(YieldCrvPtr);

	/*END OF THE CHECK OF THE INPUTS */
	
	if((NumHedgeSwpUndFra==0)
		||((NumHedgeSwpUndFra>0)&&(NumHedgeFutUndFra>0)
			&&(FraEndDate<=FutUndFraHedgeEndDates[NumHedgeFutUndFra-1])) )
		/*THEN HEDGE WITH THE FUTURE UND FRA */
	{
		FindHedges(FraEndDate,
				   NumHedgeFutUndFra,
				   FutUndFraHedgeEndDates,
				   &PrevTenorIndex,
				   &NextTenorIndex);
		
			tmp_LocalHedgeUndFraPEQ = dvector(0,NumHedgeFutUndFra-1);
			
			err = srt_f_GetHedgePEQ(FraEndDate,
							"FRA",
							FutUndFraHedgeEndDates,
							PrevTenorIndex,
							NextTenorIndex,
							0,
							0,
							FloatFreq,
							FloatBasis,
							FloatRefRate,
							YieldCurveName,
							FraPEQ,
							&tmp_LocalHedgeUndFraPEQ);

			
			srt_f_UpdateLocHedgeUndFraTermStr(FutUndFraHedgeEndDates,
												tmp_LocalHedgeUndFraPEQ,
												PrevTenorIndex,
												NextTenorIndex,
												SRT_UNDFRA_FUT,
											 
												NumLocalHedgeUndFra,
												LocalHedgeUndFraPEQ,
												LocalHedgeUndFraEndDates,
												LocalHedgeUndFraTypes);
			if(err) 
				return err;
		
			free_dvector(tmp_LocalHedgeUndFraPEQ,0,NumHedgeFutUndFra-1);

		
	}
	else if((NumHedgeFutUndFra==0)||
			(NumHedgeFutUndFra>0)&&(NumHedgeSwpUndFra>0)
			&&((FraEndDate>FutUndFraHedgeEndDates[0])))
		/*THEN HEDGE WITH THE SWAP UNDERLYING FRA */
	{

		FindHedges(FraEndDate,
				   NumHedgeSwpUndFra,
				   SwpUndFraHedgeEndDates,
				   &PrevTenorIndex,
				   &NextTenorIndex);

		tmp_LocalHedgeUndFraPEQ = dvector(0,NumHedgeSwpUndFra-1);

		err = srt_f_GetHedgePEQ(FraEndDate,
							"FRA",
							SwpUndFraHedgeEndDates,
							PrevTenorIndex,
							NextTenorIndex,
							0,
							0,
							FloatFreq,
							FloatBasis,
							FloatRefRate,
							YieldCurveName,
							FraPEQ,
							&tmp_LocalHedgeUndFraPEQ);
		if(err) 
			return err;

		srt_f_UpdateLocHedgeUndFraTermStr(SwpUndFraHedgeEndDates,
										tmp_LocalHedgeUndFraPEQ,
										PrevTenorIndex,
										NextTenorIndex,
										SRT_UNDFRA_SWAP,
									 
										NumLocalHedgeUndFra,
										LocalHedgeUndFraPEQ,
										LocalHedgeUndFraEndDates,
										LocalHedgeUndFraTypes);
			if(err) 
				return err;
			free_dvector(tmp_LocalHedgeUndFraPEQ,0,NumHedgeSwpUndFra-1);
	}

	return err;

}

/*HedgePEQ MUST BE ALLOCATED BEFORE */

Err		srt_f_GetHedgePEQ(Date				MKTTenorDate, /*SWAP MAT OR FRA THEO DATE */
						  String			HedgeTypeStr, /* "SWAP" OR "FRA"*/
						  Date				*HedgeTenorDates,
						  long				PrevTenorIndex,
						  long				NextTenorIndex,
						  SrtCompounding	FixedFreq,
						  SrtBasisCode		FixedBasis,
						  SrtCompounding	FloatFreq,
						  SrtBasisCode		FloatBasis,
						  String			FloatRefRate,
						  String			YieldCurveName,
						  double			MarketPEQ,
						  double			**HedgePEQ)
								  
{
Err				err=NULL;
String			FixedFreqStr=NULL,FixedBasisStr=NULL;
String			FloatFreqStr=NULL,FloatBasisStr=NULL;
double			MKTSensitivity,MKTConvexity;
double			PrevHedgeTenorSensitivity,NextHedgeTenorSensitivity;
double			PrevHedgeTenorConvexity,NextHedgeTenorConvexity;
SrtCrvPtr		YieldCrvPtr=NULL;
Date			SpotDate;
long			SpotLag;
Date			today;
Date			FraStartDate;
SrtHedgeType	HedgeType;

	/*GET THE YIELD CURVE PTR, THE DATE TODAY,SPOTLAG,SPOTDATE ...*/
	YieldCrvPtr = lookup_curve(YieldCurveName);	
	if(!YieldCrvPtr)
		return serror("Can not get the yield curve");
	today = get_clcndate_from_yldcrv(YieldCrvPtr);

	SpotDate = get_spotdate_from_yldcrv(YieldCrvPtr);
	SpotLag = get_spotlag_from_curve(YieldCrvPtr);

	err = translate_compounding(&FixedFreqStr,FixedFreq);
	if(err)
		return err;
	
	err = translate_basis(&FixedBasisStr,FixedBasis);
	if(err)
		return err;

	err = translate_compounding(&FloatFreqStr,FloatFreq);
	if(err)
		return err;
	
	err = translate_basis(&FloatBasisStr,FloatBasis);
	if(err)
		return err;

	err = interp_HedgeType(HedgeTypeStr,&HedgeType);
	
	if(HedgeType==SRT_SWAP)
	{
		/*GET THE MARKET SENSITIVITY AND (CONVEXITY) */
		err = swp_f_LevelPayment (SpotDate,
								  MKTTenorDate,
								  FixedFreqStr,
								  FixedBasisStr,
								  YieldCurveName,
								  FloatRefRate,
								  &MKTSensitivity);
		if(err)
			return err;

		/*IF AGGREGATION ON TWO REFERENCES */
		if(PrevTenorIndex!=NextTenorIndex)
		{
			err = srt_f_SwapConvexity(SpotDate,
									  MKTTenorDate,
									  FixedFreq,
									  FixedBasis,
									  FloatFreq,
					   				  FloatBasis,
									  FloatRefRate,
					 				  YieldCurveName,
									 &MKTConvexity);
			if(err)
				return err;
		}

		/*GET THE HEDGE'S SENSTIVITIES  AND (CONVEXITIES) */			
		err = swp_f_LevelPayment (SpotDate,
								  HedgeTenorDates[PrevTenorIndex],
								  FixedFreqStr,
								  FixedBasisStr,
								  YieldCurveName,
								  FloatRefRate,
								  &PrevHedgeTenorSensitivity);
		if(err)
			return err;
		
		/*IF AGGREGATION ON TWO REFERENCES */
		if(PrevTenorIndex!=NextTenorIndex)
		{

			err = swp_f_LevelPayment (SpotDate,
									  HedgeTenorDates[NextTenorIndex],
									  FixedFreqStr,
									  FixedBasisStr,
									  YieldCurveName,
									  FloatRefRate,
									  &NextHedgeTenorSensitivity);
			if(err)
				return err;

			err = srt_f_SwapConvexity(SpotDate,
									  HedgeTenorDates[PrevTenorIndex],
									  FixedFreq,
									  FixedBasis,
									  FloatFreq,
					   				  FloatBasis,
									  FloatRefRate,
					 				  YieldCurveName,
									 &PrevHedgeTenorConvexity);
			if(err)
				return err;

			err = srt_f_SwapConvexity(SpotDate,
									  HedgeTenorDates[NextTenorIndex],
									  FixedFreq,
									  FixedBasis,
									  FloatFreq,
					   				  FloatBasis,
									  FloatRefRate,
					 				  YieldCurveName,
									 &NextHedgeTenorConvexity);
		}
	}

	else if(HedgeType==SRT_FRA)
	{

		/*GET THE MARKET SENSITIVITY AND CONVEXITY */

		FraStartDate = add_unit (MKTTenorDate, 
								  -12/(int)FloatFreq, SRT_MONTH, 								  
								  NO_BUSDAY_CONVENTION);
		
		err = swp_f_LevelPayment (FraStartDate,
								  MKTTenorDate,
								  FloatFreqStr,
								  FloatBasisStr,
								  YieldCurveName,
								  FloatRefRate,
								  &MKTSensitivity);


		if(PrevTenorIndex!=NextTenorIndex)
		{
			err = srt_f_FRAConvexity(FraStartDate,
									 MKTTenorDate,
									 FloatBasis,
									 YieldCurveName,
									 &MKTConvexity);
			if(err)
				return err;
			
		}
		
		/*GET THE HEDGE SENSITIVITY AND CONVEXITY*/

		FraStartDate = add_unit (HedgeTenorDates[PrevTenorIndex], 
								  -12/(int)FloatFreq, SRT_MONTH, 								  
								  NO_BUSDAY_CONVENTION);	
		

		err = swp_f_LevelPayment (FraStartDate,
								  HedgeTenorDates[PrevTenorIndex],
								  FloatFreqStr,
								  FloatBasisStr,
								  YieldCurveName,
								  FloatRefRate,
								  &PrevHedgeTenorSensitivity);


		if(PrevTenorIndex!=NextTenorIndex)
		{
		 	/* FraStartDate IS COMPUTED ABOVE */
			err = srt_f_FRAConvexity(FraStartDate,
									 HedgeTenorDates[PrevTenorIndex],
									 FloatBasis,
									 YieldCurveName,
									 &PrevHedgeTenorConvexity);
			if(err)
				return err;

			FraStartDate = add_unit (HedgeTenorDates[NextTenorIndex], 
								  -12/(int)FloatFreq, SRT_MONTH, 								  
								  NO_BUSDAY_CONVENTION);

			err = swp_f_LevelPayment (FraStartDate,
									  HedgeTenorDates[NextTenorIndex],
									  FloatFreqStr,
									  FloatBasisStr,
									  YieldCurveName,
									  FloatRefRate,
									  &NextHedgeTenorSensitivity);
		
			err = srt_f_FRAConvexity(FraStartDate,
									 HedgeTenorDates[NextTenorIndex],
									 FloatBasis,
									 YieldCurveName,
									 &NextHedgeTenorConvexity);
			if(err)
				return err;
		}

	}

	/*COMPUTATION OF THE HEDGE PEQ */
	if(PrevTenorIndex!=NextTenorIndex)
	{
		(*HedgePEQ)[PrevTenorIndex]=
			MarketPEQ*(NextHedgeTenorConvexity*MKTSensitivity-NextHedgeTenorSensitivity*MKTConvexity)
		 /(PrevHedgeTenorSensitivity*NextHedgeTenorConvexity-PrevHedgeTenorConvexity*NextHedgeTenorSensitivity);

		(*HedgePEQ)[NextTenorIndex]=
		MarketPEQ*(PrevHedgeTenorSensitivity*MKTConvexity-MKTSensitivity*PrevHedgeTenorConvexity)
		/(PrevHedgeTenorSensitivity*NextHedgeTenorConvexity-PrevHedgeTenorConvexity*NextHedgeTenorSensitivity);

	} 

	else if(PrevTenorIndex==NextTenorIndex)
	{

		(*HedgePEQ)[PrevTenorIndex]= MarketPEQ*MKTSensitivity/PrevHedgeTenorSensitivity;
	}

	return err;

}

/* DECOMPOSITION OF A SWAP POSITION ON THE UNDERLYING FRA'S POSITION */	

Err		srt_f_GetSwapUndFRAPrincEq(Date				    SwapStartDate,
								   Date				    SwapEndDate,
								   SrtCompounding		FixedFreq,
								   SrtBasisCode			FixedBasis,
								   SrtCompounding		FloatFreq,
								   SrtBasisCode			FloatBasis,
								   String				FloatRefRate,
								   String				YieldCurveName,
								   double				QSwapPEq,
								   
								   /*OUTPUT	*/
								   long					*NumSwapUndFra,
								   double				**SwapUndFraPEq,
								   long					**SwapUndFraEndDates)
		

								 /*QSwapPEq IS THE QUATERLY PRINCIPAL EQUIVALENT */
{

Err			err = NULL;
Date		*FixingDates=NULL, *StartDates=NULL, *EndDates=NULL;
double		*Coverage=NULL;
long		NumPayDates;
SrtCrvPtr	CrvPtr=NULL; 
long		today;
long		i;								 
SwapDP		*FloatLegDP=NULL;
long		*PayDates;
			       
	/* GET THE SWAP FLOAT LEG DATES AND COVERAGES */
	FloatLegDP = (SwapDP*) malloc(sizeof(SwapDP));

	err = swp_f_setSwapDP(SwapStartDate, 
						  SwapEndDate, 
						  FloatFreq, 
						  FloatBasis, 
						  FloatLegDP); /* INITIALISATION OF THE FLOATING LEG */
	if(err)
		return err;

	CrvPtr = lookup_curve(YieldCurveName);
	today = get_clcndate_from_yldcrv(CrvPtr); /*GET THE DATE TODAY */

	err = swp_f_make_FloatLegDatesAndCoverages(FloatLegDP, 
												   today, 
												   &PayDates, 
												   &NumPayDates,
												   &FixingDates, 
												   &StartDates,
												   &EndDates, 
												   &Coverage, 
												   NumSwapUndFra);

	(*SwapUndFraPEq) = dvector(0,(*NumSwapUndFra)-1);
	(*SwapUndFraEndDates) = lngvector(0,(*NumSwapUndFra)-1);

	if(!SwapUndFraPEq)
		return serror("Allocation Memory Failure in srt_f_GetSwapUndFRAPrincEq");

	for(i=0;i<(*NumSwapUndFra);i++) 
	{
		(*SwapUndFraPEq)[i] = QSwapPEq;
		(*SwapUndFraEndDates)[i] =add_unit (DTOL (SwapStartDate), 
										(i+1)*(12/(int)FloatFreq), SRT_MONTH, 								  
										NO_BUSDAY_CONVENTION);
	}

	if(FloatLegDP) free(FloatLegDP); FloatLegDP=NULL;
	return err;


}

/*COMPUTATION OF THE CONVEXITY OF THE SWAP */
Err		srt_f_SwapConvexity(Date				SwapStartDate,  /*SwapStartDate = SpotDate */
							Date				SwapEndDate,
							SrtCompounding		FixedFreq,
							SrtBasisCode		FixedBasis,
						    SrtCompounding		FloatFreq,
					   	    SrtBasisCode		FloatBasis,
							String				FloatRefRate,
					 	    String				YieldCurveName,
							double				*SwapConvexity)

{

Err			err=NULL;
long		today,Date;
SrtCrvPtr	CrvPtr;
double		**FraSwapWeightsMatrix=NULL;
String		FixedFreqStr=NULL,FixedBasisStr=NULL;
long		FraSwapWeightsSize;
double		*DiscountFactor=NULL;
double		Coverage;

	err = srt_FraSwapWeights(SwapStartDate,
						   SwapEndDate,
						   FixedFreq,
						   FixedBasis,
						   FloatFreq,
						   FloatBasis,
						   YieldCurveName,
						   /*OUTPUT*/
						   &FraSwapWeightsMatrix,
						   &FraSwapWeightsSize);

	if(err)
		return err;


	err = translate_compounding(&FixedFreqStr,FixedFreq);
	if(err)
		return err;
	
	err= translate_basis(&FixedBasisStr,FixedBasis);
	if(err)
		return err;

	CrvPtr = lookup_curve(YieldCurveName);
	if(!CrvPtr)
		return serror("Error in srt_f_SwapConvexity: can not get the yield curve ");

	/* GET THE DATE TODAY */
	today = get_clcndate_from_yldcrv(CrvPtr); 

	(*SwapConvexity) = 0.0;

	Date  = add_unit (SwapStartDate,
						(FraSwapWeightsSize-1)*12/(int)FloatFreq,
						SRT_MONTH, 								  
						NO_BUSDAY_CONVENTION);

	DiscountFactor = dvector(0,1);
	DiscountFactor[0] = swp_f_df((Ddate)today,(Ddate)Date,YieldCurveName);
	DiscountFactor[1] = swp_f_df((Ddate)today,(Ddate)SwapEndDate,YieldCurveName);
	Coverage = coverage(Date,SwapEndDate,FloatBasis);

	(*SwapConvexity) = -2*pow(Coverage*DiscountFactor[1]/DiscountFactor[0],2)
						*FraSwapWeightsMatrix[FraSwapWeightsSize][FraSwapWeightsSize];
	return err;

	if(DiscountFactor) free_dvector(DiscountFactor,0,1);
	

}

/*COMPUTATION OF THE CONVEXITY OF THE FRA */

Err		srt_f_FRAConvexity(long				FraStartDate,
						   long				FraEndDate,
						   SrtBasisCode		FloatBasis,
						   String			YieldCurveName,
						 /*OUTPUT*/
						   double			*FraConvexity)
{

Err			err=NULL;
double		df_Start,df_End;
SrtCrvPtr	CrvPtr;
double		today;

	
	CrvPtr = lookup_curve(YieldCurveName);
	if(!CrvPtr)
		return serror("Error in srt_f_FRAConvexity: can not get the yield curve ");

	/* GET THE DATE TODAY */
	today = get_clcndate_from_yldcrv(CrvPtr); 
	
	/* GET THE ASSOCIATED FRA DISCOUNT FACTORS */
	df_Start = swp_f_df((Ddate)today,(Ddate)FraStartDate,YieldCurveName);
	df_End   = swp_f_df((Ddate)today,(Ddate)FraEndDate,YieldCurveName);

	/* COMPUTE THE CONVEXITY OF THE FRA */

	*FraConvexity = -coverage(FraStartDate,FraEndDate,FloatBasis)*pow(df_End,2)/df_Start;

	return err;
}


/* DECOMPOSITION OF A FRA POSITION ON THE UNDERLYING SWAP'S POSITION */	

Err		srt_f_GetFraUndSWAPPrincEq(Date					SpotDate,
								   Date				    FraEndDate, /*THEO END DATE ? */

								   /*SOME SWAP DETAILS */
								   SrtCompounding		SwpFixedFreq,
								   SrtBasisCode			SwpFixedBasis,
								   SrtCompounding		SwpFloatFreq,
								   SrtBasisCode			SwpFloatBasis,
								   String				SwpFloatRefRate,
								   /*FIXED AND FLOATING FEATURES HAVE TO BE EQUAL */
								   String				YieldCurveName,
								   /*FRA PRINCIPAL EQUIVALENT */
								   double				FraPEQ,
								   
								   /*OUTPUT	*/
								   long					*NumFraUndSwap,
								   double				**FraUndSwapPEq,
								   long					**FraUndSWAPEndDates)
		

{

Err				err=NULL;
SwapDP			*SwpFloatLegDP=NULL;
SrtCrvPtr		CrvPtr;
long			today;
long			NumPayDates;
long			*FixingDates=NULL,*StartDates=NULL,*EndDates=NULL;
double			*Coverage=NULL;
double			**FraSwapWeightsMatrix=NULL;
String			SwpFixedFreqStr=NULL,SwpFixedBasisStr=NULL;
String			SwpFloatFreqStr=NULL,SwpFloatBasisStr=NULL;
int				i;
double			FixedRefSwapLevel,FloatRefSwapLevel;
Date			SwapEndDate ;
long			*PayDates;

	/*GET THE UNDERLYINGS SWAP FEATURES */

	SwpFloatLegDP = (SwapDP*) malloc(sizeof(SwapDP));	
	
	err = swp_f_setSwapDP(SpotDate,	  /* ALL THE SWAP START AT SPOT DATE */
						  FraEndDate, /*THE LAST SWAP MATURES AT THE FRA MATURITY*/
						  SwpFloatFreq, 
						  SwpFloatBasis,
						  SwpFloatLegDP); 
	if(err)
		return err;

	CrvPtr = lookup_curve(YieldCurveName);
	today = get_clcndate_from_yldcrv(CrvPtr);

	/*GET THE DIFFERENT MATURITIES OF THE SWAPS */
	err = swp_f_make_FloatLegDatesAndCoverages(SwpFloatLegDP, 
												   today, 
												   &PayDates,/* THE SWAPS MATURITIES  (PayDates)*/ 
												   &NumPayDates,/* THE NUMBER OF SWAPS */
												   &FixingDates, 
												   &StartDates,
												   &EndDates, 
												   &Coverage, 
												   NumFraUndSwap);

	if(err)
		return err;
	
	/*GET THE SWAP PRINCIPAL EQUIVALENT */
	(*FraUndSwapPEq) = dvector(0,(*NumFraUndSwap));
	(*FraUndSWAPEndDates) =lngvector(0,(*NumFraUndSwap));
	
	for(i=0;i<(*NumFraUndSwap);i++) 	(*FraUndSwapPEq)[i] = 0;
	for(i=0;i<(*NumFraUndSwap)+1;i++)
	{
		(*FraUndSWAPEndDates)[i] = add_unit (DTOL (SpotDate), 
											i*12/(int)SwpFloatFreq, SRT_MONTH, 								  
											NO_BUSDAY_CONVENTION);
	}


	err = translate_compounding(&SwpFixedFreqStr,SwpFixedFreq);
	if(err)
		return err;
	
	err = translate_basis(&SwpFixedBasisStr,SwpFixedBasis);
	if(err)
		return err;

	err = translate_compounding(&SwpFloatFreqStr,SwpFloatFreq);
	if(err)
		return err;
	
	err = translate_basis(&SwpFloatBasisStr,SwpFloatBasis);
	if(err)
		return err;

	
	err = swp_f_LevelPayment(SpotDate,
							 FraEndDate,
							 SwpFixedFreqStr,
							 SwpFixedBasisStr,
							 YieldCurveName,
							 SwpFloatRefRate,
							&FixedRefSwapLevel);
			
	err = swp_f_LevelPayment(SpotDate,
							 FraEndDate,
							 SwpFloatFreqStr,
							 SwpFloatBasisStr,
							 YieldCurveName,
							 SwpFloatRefRate,
							&FloatRefSwapLevel);
	

	(*FraUndSwapPEq)[(*NumFraUndSwap)-1]=FraPEQ*FloatRefSwapLevel/FixedRefSwapLevel;

	if((*NumFraUndSwap)>1)
	{
		SwapEndDate = add_unit (DTOL (SpotDate), 
							  ((*NumFraUndSwap)-1)*12/(int)SwpFloatFreq, SRT_MONTH, 								  
							  NO_BUSDAY_CONVENTION);

		err = swp_f_LevelPayment(SpotDate,
							 SwapEndDate,
							 SwpFixedFreqStr,
							 SwpFixedBasisStr,
							 YieldCurveName,
							 SwpFloatRefRate,
							&FixedRefSwapLevel);
			
		err = swp_f_LevelPayment(SpotDate,
							 SwapEndDate,
							 SwpFloatFreqStr,
							 SwpFloatBasisStr,
							 YieldCurveName,
							 SwpFloatRefRate,
							&FloatRefSwapLevel);

		(*FraUndSwapPEq)[(*NumFraUndSwap)-2]=-FraPEQ*FloatRefSwapLevel/FixedRefSwapLevel;
	}




	if(SwpFloatLegDP) free(SwpFloatLegDP); SwpFloatLegDP=NULL;

}

	


	





