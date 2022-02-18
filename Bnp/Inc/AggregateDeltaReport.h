//
//#ifndef AGGREGATEDELTA_H_
//#define AGGREGATEDELTA_H_
//
// void   srt_f_UpdateGlobalHedgeUndFraTermStr(Date		*GlobalHedgeUndFraEndDates,
//											Date
//*LocalHedgeUndFraEndDates, 											long		 NumLocalHedgeUndFra, 											long
//NumGlobalHedgeUndFra, 											double		*LocalHedgeUndFraPEQ, 											double
//**GlobalHedgeUndFraPEQ);
//
// void 	srt_f_UpDateHedgeTermStr(Date		 *HedgeTenorDates,
//								 long		 NumHedgeTenors,
//								 Date *SwapHedgeTenorDates, 								 long
//NumSwpHedgeTenors, 								 double		 *SwapHedgePEQ, 								 long		 PrevTenorIndex, 								 long
//NextTenorIndex,
//
//								 double		 **HedgePEQ);
//
//
//
// void srt_f_UpdateLocHedgeUndFraTermStr(Date				*UndFraHedgeEndDates,
//									   double
//*tmp_HedgeUndFraPEQ, 									   long				PrevTenorIndex, 									   long
//NextTenorIndex, 									   SrtUndFRAType    UndFraType,
//
//									   /*OUTPUT */
//									   long
//*NumLocalHedgeUndFra, 									   double			**LocalHedgeUndFraPEQ, 									   Date
//**LocalHedgeUndFraEndDates, 									   SrtUndFRAType	**LocalHedgeUndFraTypes);
//
//
//
// Err		srt_f_SwapConvexity(Date				SwapStartDate,
//							Date				SwapEndDate,
//							SrtCompounding		FixedFreq,
//							SrtBasisCode		FixedBasis,
//						    SrtCompounding		FloatFreq,
//					   	    SrtBasisCode		FloatBasis,
//							String FloatRefRate, 					 	    String
//YieldCurveName, 							double				*SwapConvexity);
//
// Err		srt_f_FRAConvexity(long				FraStartDate,
//						   long				FraEndDate,
//						   SrtBasisCode		FloatBasis,
//						   String			YieldCurveName,
//						 /*OUTPUT*/
//						   double			*FraConvexity);
//
//
//
// Err	 UpdateUndFRATermStruct(Date				DateToInsert,
//							Date
//**SetOfRefDates, 							String				UndFRATypeStr, 							SrtUndFRAType
//**SetOfRefUndFRATypes, 							long				*NumUndFra);
//
//
// Err GetTenorsPEQ(String		RefTenors,
//				 long		NumMarketTenors,
//				 String		*MarketTenors,
//				 double		*MarketTenorsPEQ,
//
//				 /*OUTPUT*/
//				 double		*RefTenorPEQ);
//
// SRT_Boolean ExitsHedgeTenors(Date		MarketTenorDate, /*THE SWAP TO HEDGE */
//						long		NumHedgeTenors,
//						Date		*HedgeTenorDates);
//
//
// Err FindHedges(Date			MarketTenorDate,
//				   long			NumHedgeTenors,
//				   Date			*HedgeTenorDates,
//				   long			*PrevTenorIndex,
//				   long			*NextTenorIndex);
//
// Err		srt_f_GetHedgePEQ(Date				MKTTenorDate, /*SWAP MAT OR FRA PAYT DATE
// */ 						  String			HedgeTypeStr, /* "SWAP" OR "FRA"*/ 						  Date
//*HedgeTenorDates, 						  long				PrevTenorIndex, 						  long
//NextTenorIndex, 						  SrtCompounding	FixedFreq, 						  SrtBasisCode		FixedBasis, 						  SrtCompounding
//FloatFreq, 						  SrtBasisCode		FloatBasis, 						  String			FloatRefRate, 						  String
//YieldCurveName, 						  double			MarketPEQ, 						  double			**HedgePEQ);
//
// Err		srt_f_GetSwapUndFRAPrincEq(Date				    SwapStartDate,
//								   Date
//SwapEndDate, 								   SrtCompounding		FixedFreq, 								   SrtBasisCode			FixedBasis,
//								   SrtCompounding
//FloatFreq, 								   SrtBasisCode			FloatBasis, 								   String
//FloatRefRate, 								   String				YieldCurveName, 								   double
//QSwapPEq,
//
//								   /*OUTPUT	*/
//								   long
//*NumSwapUndFra, 								   double				**SwapUndFraPEq, 								   long
//**PayDates);
//
//
// Err		srt_f_GetFraUndSWAPPrincEq(Date					SpotDate,
//								   Date
//FraEndDate, 								   SrtCompounding		SwpFixedFreq, 								   SrtBasisCode SwpFixedBasis, 								   SrtCompounding
//SwpFloatFreq, 								   SrtBasisCode			SwpFloatBasis, 								   String
//SwpFloatRefRate, 								   String				YieldCurveName, 								   double
//FraPEq, 								   long					*NumFraUndSwap, 								   double
//**FraUndSwapPEq, 								   long					**PayDates);
//
// Err	srt_f_AggregateFraTenor(Date			FraTenorDate,
//							double			FraPEQ,
//							long
//NumSwpUndFraHedgeTenor, 							long			NumFutUndFraHedgeTenor, 							Date
//*SwpUndFraHedgeTenorDates, 							Date			*FutUndFraHedgeTenorDates, 							SrtCompounding
//FloatFreq, 							SrtBasisCode	FloatBasis, 							String			FloatRefRate, 							String
//YieldCurveName, 							long			*NumHedgeUndFra, 							double			**HedgeUndFraPEQ,
//							Date
//**HedgeUndFraTenorDates, 							SrtUndFRAType	**HedgeUndFraTypes);
//
//
//
//#endif;