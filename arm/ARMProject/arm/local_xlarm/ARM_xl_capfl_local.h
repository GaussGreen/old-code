#ifndef ARM_XL_CAPFL_LOCAL_H
#define ARM_XL_CAPFL_LOCAL_H


__declspec(dllexport) LPXLOPER WINAPI Local_LIBORCF (LPXLOPER XL_startDate,
													 LPXLOPER XL_endDate,
													 LPXLOPER XL_isItCapOrFloor,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_liborType,
													 LPXLOPER XL_spread,
													 LPXLOPER XL_resetFreq,
													 LPXLOPER XL_payFreq,
													 LPXLOPER XL_ccy);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORCF (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_liborType,
														 LPXLOPER XL_spread,
														 LPXLOPER XL_resetFreq,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_LIBORFLEXCF (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_nbEx,
														 LPXLOPER XL_exerciseType,
														 LPXLOPER XL_liborType,
														 LPXLOPER XL_spread,
														 LPXLOPER XL_resetFreq,
														 LPXLOPER XL_payFreq,
														 LPXLOPER XL_currency);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_LIBORFLEXCF (LPXLOPER XL_startDate,
															 LPXLOPER XL_endDate,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_nbEx,
															 LPXLOPER XL_exerciseType,
															 LPXLOPER XL_liborType,
															 LPXLOPER XL_spread,
															 LPXLOPER XL_resetFreq,
															 LPXLOPER XL_payFreq,
															 LPXLOPER XL_ccy);

__declspec(dllexport) LPXLOPER WINAPI Local_CAPFLOOR (LPXLOPER XL_swapLeg,
													  LPXLOPER XL_isItCapOrFloor,
													  LPXLOPER XL_strike);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CAPFLOOR (LPXLOPER XL_swapLeg,
														  LPXLOPER XL_isItCapOrFloor,
														  LPXLOPER XL_strike);

__declspec(dllexport) LPXLOPER WINAPI Local_FLEXCF (LPXLOPER XL_swapLegId,
													LPXLOPER XL_isItCapOrFloor,
													LPXLOPER XL_strike,
													LPXLOPER XL_nbEx,
													LPXLOPER XL_exerciseType);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FLEXCF (LPXLOPER XL_swapLegId,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_strike,
														LPXLOPER XL_nbEx,
														LPXLOPER XL_exerciseType);

__declspec(dllexport) LPXLOPER WINAPI Local_MATCAPFLOOR(LPXLOPER XL_swapLegId,
														LPXLOPER XL_annuity,
														LPXLOPER XL_initNominal,
														LPXLOPER XL_isTRI,
														LPXLOPER XL_capOrFloor,
														LPXLOPER XL_coeff,
														LPXLOPER XL_firstTRIstrike,
														LPXLOPER XL_minStrikes,
														LPXLOPER XL_isDigitalPayoff,
														LPXLOPER XL_increasingCoef,                                                        
                                                        LPXLOPER XL_maxMaturityDate);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MATCAPFLOOR(LPXLOPER XL_swapLegId,
															LPXLOPER XL_annuity,
															LPXLOPER XL_initNominal,
															LPXLOPER XL_isTRI,
															LPXLOPER XL_capOrFloor,
															LPXLOPER XL_coeff,
															LPXLOPER XL_firstTRIstrike,
															LPXLOPER XL_minStrikes,
															LPXLOPER XL_isDigitalPayoff,
															LPXLOPER XL_increasingCoef);

__declspec(dllexport) LPXLOPER WINAPI Local_STICKY (LPXLOPER XL_swapLeg,
													LPXLOPER XL_isItCapOrFloor,
													LPXLOPER XL_strike,
													LPXLOPER XL_spreadDates,
													LPXLOPER XL_spreadValues,
													LPXLOPER XL_kRefVal);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_STICKY (LPXLOPER XL_swapLeg,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_strike,
														LPXLOPER XL_spreadDates,
														LPXLOPER XL_spreadValues,
														LPXLOPER XL_kRefVal);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADOPTION (LPXLOPER XL_startDate,
															  LPXLOPER XL_endDate,
															  LPXLOPER XL_capOrFloor,
															  LPXLOPER XL_strike,
															  LPXLOPER XL_liborType1,
															  LPXLOPER XL_liborType2,
															  LPXLOPER XL_weight1,
															  LPXLOPER XL_weight2,
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_resetFreq,
															  LPXLOPER XL_payFreq,
															  LPXLOPER XL_resetTiming,
															  LPXLOPER XL_payTiming,
															  LPXLOPER XL_currency,
                                                              LPXLOPER XL_resetGap,
															  LPXLOPER XL_intRule,
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_fixing1,
															  LPXLOPER XL_fixing2,
															  LPXLOPER XL_cptStrikeMethod);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADOPTION (LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_liborType1,
																  LPXLOPER XL_liborType2,
																  LPXLOPER XL_weight1,
																  LPXLOPER XL_weight2,
																  LPXLOPER XL_daycount,
																  LPXLOPER XL_resetFreq,
																  LPXLOPER XL_payFreq,
																  LPXLOPER XL_resetTiming,
																  LPXLOPER XL_payTiming,
																  LPXLOPER XL_currency,
                                                                  LPXLOPER XL_resetGap);
																  LPXLOPER XL_intRule,
																  LPXLOPER XL_stubRule,
																  LPXLOPER XL_fixing1,
															      LPXLOPER XL_fixing2,
																  LPXLOPER XL_cptStrikeMethod);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADOPTION(		LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_leg1Weights,
																		LPXLOPER XL_leg2Weights,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADOPTIONWithLegs (LPXLOPER XL_firstLeg,
																	  LPXLOPER XL_secondLeg,
																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_weight1,
																	  LPXLOPER XL_weight2,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_cptStrikeMethod);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADOPTIONWithLegs (LPXLOPER XL_firstLeg,
																		  LPXLOPER XL_secondLeg,
																		  LPXLOPER XL_capOrFloor,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_weight1,
																		  LPXLOPER XL_weight2,
																		  LPXLOPER XL_fixing1,
																		  LPXLOPER XL_fixing2,
																		  LPXLOPER XL_cptStrikeMethod);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITAL(LPXLOPER XL_startDate,
															  LPXLOPER XL_endDate,
															  LPXLOPER XL_capOrFloor,
															  LPXLOPER XL_strike,
															  LPXLOPER XL_spreadVect,
															  LPXLOPER XL_payoff,
															  LPXLOPER XL_liborType1,
															  LPXLOPER XL_liborType2,
															  LPXLOPER XL_weightVect,
															  LPXLOPER XL_daycount,
															  LPXLOPER XL_resetFreq,
															  LPXLOPER XL_payFreq,
															  LPXLOPER XL_resetTiming,
															  LPXLOPER XL_payTiming,
															  LPXLOPER XL_currency,
                                                              LPXLOPER XL_resetGap,
															  LPXLOPER XL_intRule,
															  LPXLOPER XL_stubRule,
															  LPXLOPER XL_fixing1,
															  LPXLOPER XL_fixing2);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITAL(LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_spreadVect,
																  LPXLOPER XL_payoff,
																  LPXLOPER XL_liborType1,
																  LPXLOPER XL_liborType2,
																  LPXLOPER XL_weightVect,
																  LPXLOPER XL_daycount,
																  LPXLOPER XL_resetFreq,
																  LPXLOPER XL_payFreq,
																  LPXLOPER XL_resetTiming,
																  LPXLOPER XL_payTiming,
																  LPXLOPER XL_currency,
																  LPXLOPER XL_resetGap,
																  LPXLOPER XL_intRule,
																  LPXLOPER XL_stubRule,
																  LPXLOPER XL_fixing1,
																  LPXLOPER XL_fixing2);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADDIGITAL(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_payoff,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_leg1Weights,
																		LPXLOPER XL_leg2Weights,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALWithLegs(LPXLOPER XL_firstLeg,
																	  LPXLOPER XL_secondLeg,
																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_spread1,
																	  LPXLOPER XL_spread2,
																	  LPXLOPER XL_payoff,
																	  LPXLOPER XL_weight1,
																	  LPXLOPER XL_weight2,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_slopeFlag,
																	  LPXLOPER XL_cptStrikeMethod);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALWithLegs(LPXLOPER XL_firstLeg,
																		  LPXLOPER XL_secondLeg,
																		  LPXLOPER XL_capOrFloor,
																		  LPXLOPER XL_strike,
																		  LPXLOPER XL_spread1,
																		  LPXLOPER XL_spread2,
																		  LPXLOPER XL_payoff,
																		  LPXLOPER XL_weight1,
																		  LPXLOPER XL_weight2,
																		  LPXLOPER XL_fixing1,
																		  LPXLOPER XL_fixing2,
																		  LPXLOPER XL_slopeFlag,
																		  LPXLOPER XL_cptStrikeMethod);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_CORRIDORDBLCONDITION(LPXLOPER XL_startDate,
																	 LPXLOPER XL_endDate,
																	 LPXLOPER XL_DigitalcapOrFloor,
																	 LPXLOPER XL_SpreadcapOrFloor,
																	 LPXLOPER XL_digitalBarrier, 
																	 LPXLOPER XL_spreadBarrier,
																	 LPXLOPER XL_spreads,
																	 LPXLOPER XL_payIndexParams,
																	 LPXLOPER XL_spreadIdxTypes,
																	 LPXLOPER XL_Weights,
																	 LPXLOPER XL_prodParamDatas,
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_fixing1,
																	 LPXLOPER XL_fixing2,
																	 LPXLOPER XL_fixing3,
																	 LPXLOPER XL_fixingPay,
																	 LPXLOPER XL_freezeFixing,
																	 LPXLOPER XL_DigitalSpreadCorrel,
																	 LPXLOPER XL_ShiftVolCorrel,
																	 LPXLOPER XL_fixing4);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_CORRIDORDBLCONDITION(LPXLOPER XL_startDate,
																		 LPXLOPER XL_endDate,
																		 LPXLOPER XL_DigitalcapOrFloor,
																		 LPXLOPER XL_SpreadcapOrFloor,
																		 LPXLOPER XL_digitalBarrier, 
																		 LPXLOPER XL_spreadBarrier,
																		 LPXLOPER XL_spreads,  
																		 LPXLOPER XL_payIndexParams,
																		 LPXLOPER XL_spreadIdxTypes,
																		 LPXLOPER XL_Weights,
																		 LPXLOPER XL_prodParamDatas,
																		 LPXLOPER XL_currency,
																		 LPXLOPER XL_fixing1,
																		 LPXLOPER XL_fixing2,
																		 LPXLOPER XL_fixing3,
																		 LPXLOPER XL_fixingPay,
																		 LPXLOPER XL_freezeFixing,
																		 LPXLOPER XL_DigitalSpreadCorrel,
																		 LPXLOPER XL_ShiftVolCorrel,
																		 LPXLOPER XL_fixing4);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALFLT(LPXLOPER XL_startDate,
															     LPXLOPER XL_endDate,
															     LPXLOPER XL_capOrFloor,
															     LPXLOPER XL_strike,
															     LPXLOPER XL_spread,//Vector
															     LPXLOPER XL_payIdxType,
															     LPXLOPER XL_liborType1,
															     LPXLOPER XL_liborType2,
																 LPXLOPER XL_weight,//Vector Weights, slopeflag and cptStrikeMethod
															     LPXLOPER XL_daycount,
															     LPXLOPER XL_resetFreq,
															     LPXLOPER XL_payFreq,
															     LPXLOPER XL_resetTiming,
															     LPXLOPER XL_payTiming,
															     LPXLOPER XL_currency,
                                                                 LPXLOPER XL_resetGap,
																 LPXLOPER XL_intRule,
																 LPXLOPER XL_stubRule,
																 LPXLOPER XL_fixing1,
																 LPXLOPER XL_fixing2);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALFLT(LPXLOPER XL_startDate,
																	 LPXLOPER XL_endDate,
																	 LPXLOPER XL_capOrFloor,
																	 LPXLOPER XL_strike,
																	 LPXLOPER XL_spread,//Vector
																	 LPXLOPER XL_payIdxType,
																	 LPXLOPER XL_liborType1,
																	 LPXLOPER XL_liborType2,
																	 LPXLOPER XL_weight,//Vector Weights, slopeflag and cptStrikeMethod
																	 LPXLOPER XL_daycount,
																	 LPXLOPER XL_resetFreq,
																	 LPXLOPER XL_payFreq,
																	 LPXLOPER XL_resetTiming,
																	 LPXLOPER XL_payTiming,
																	 LPXLOPER XL_currency,
																	 LPXLOPER XL_resetGap,
																	 LPXLOPER XL_intRule,
																	 LPXLOPER XL_stubRule,
																	 LPXLOPER XL_fixing1,
																	 LPXLOPER XL_fixing2);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADDIGITALFLT(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_liborIdxP,
																		LPXLOPER XL_leg1Weights,
																		LPXLOPER XL_leg2Weights,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADCORRIDOR( LPXLOPER XL_startDate,
															    LPXLOPER XL_endDate,
 															    LPXLOPER XL_capOrFloor,
															    LPXLOPER XL_strike,
															    LPXLOPER XL_spreads,
																LPXLOPER XL_payIdxType,
															    LPXLOPER XL_spreadIdxTypes,
																LPXLOPER XL_sprdWeights,
																LPXLOPER XL_cptStrikeMethod,
																LPXLOPER XL_prodParamDatas,
															    LPXLOPER XL_currency,
																LPXLOPER XL_fixing1,
																LPXLOPER XL_fixing2,
																LPXLOPER XL_fixing3,
																LPXLOPER XL_freezeFixing);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADCORRIDOR( LPXLOPER XL_startDate,
																	LPXLOPER XL_endDate,
 																	LPXLOPER XL_capOrFloor,
																	LPXLOPER XL_strike,
																	LPXLOPER XL_spreads,
																	LPXLOPER XL_payIdxType,
																	LPXLOPER XL_spreadIdxTypes,
																	LPXLOPER XL_sprdWeights,
																	LPXLOPER XL_cptStrikeMethod,
																	LPXLOPER XL_prodParamDatas,
																	LPXLOPER XL_currency,
																	LPXLOPER XL_fixing1,
																	LPXLOPER XL_fixing2,
																	LPXLOPER XL_fixing3,
																	LPXLOPER XL_freezeFixing);


//------------------------------------------------------------------------------------------//
// Interface functions for VMS spread options corridor										//
//------------------------------------------------------------------------------------------//
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADCORRIDORVMS(LPXLOPER XL_startDate,
																  LPXLOPER XL_endDate,
 																  LPXLOPER XL_capOrFloor,
																  LPXLOPER XL_strike,
																  LPXLOPER XL_spreads,
																  LPXLOPER XL_payIndexParams,
															      LPXLOPER XL_CMSIndexes1,
																  LPXLOPER XL_CMSIndexes2,
															      LPXLOPER XL_sprdWeights_slope,
															      LPXLOPER XL_prodParamDatas,
															      LPXLOPER XL_currency,
															      LPXLOPER XL_fixing1,
															      LPXLOPER XL_fixing2,
															      LPXLOPER XL_fixingPay,
															      LPXLOPER XL_freezeFixing);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADCORRIDORVMS(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
 																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_spreads,
																	  LPXLOPER XL_payIndexParams,
																	  LPXLOPER XL_CMSIndexes1,
																	  LPXLOPER XL_CMSIndexes2,
																	  LPXLOPER XL_sprdWeights_slope,
																	  LPXLOPER XL_prodParamDatas,
																	  LPXLOPER XL_currency,
																	  LPXLOPER XL_fixing1,
																	  LPXLOPER XL_fixing2,
																	  LPXLOPER XL_fixingPay,
																	  LPXLOPER XL_freezeFixing);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADCORRIDORVMS(LPXLOPER XL_startDate,
																	  LPXLOPER XL_endDate,
 																	  LPXLOPER XL_capOrFloor,
																	  LPXLOPER XL_strike,
																	  LPXLOPER XL_spreads,
																	  LPXLOPER XL_payIndexParams,
																	  LPXLOPER XL_FirstIdxTypes,
																	  LPXLOPER XL_SecondIdxTypes,
																      LPXLOPER XL_sprdWeights_slope,
																      LPXLOPER XL_prodParamDatas,
																      LPXLOPER XL_currency,
																      LPXLOPER XL_fixing1,
																      LPXLOPER XL_fixing2,
																      LPXLOPER XL_fixing3,
																      LPXLOPER XL_freezeFixing);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_SPREADDIGITALFLTWithLegs(LPXLOPER XL_firstLeg,
																		 LPXLOPER XL_secondLeg,
																		 LPXLOPER XL_capOrFloor,
																		 LPXLOPER XL_strike,
																		 LPXLOPER XL_spread1,
																		 LPXLOPER XL_spread2,
																		 LPXLOPER XL_payIdxType,
																		 LPXLOPER XL_weightVec,
																		 LPXLOPER XL_fixing1,
																		 LPXLOPER XL_fixing2,
																		 LPXLOPER XL_slopeFlag,
																		 LPXLOPER XL_cptStrikeMethod);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_SPREADDIGITALFLTWithLegs(LPXLOPER XL_firstLeg,
																			  LPXLOPER XL_secondLeg,
																			  LPXLOPER XL_capOrFloor,
																			  LPXLOPER XL_strike,
																			  LPXLOPER XL_spread1,
																			  LPXLOPER XL_spread2,
																			  LPXLOPER XL_payIdxType,
																			  LPXLOPER XL_weightVec,
																			  LPXLOPER XL_fixing1,
																			  LPXLOPER XL_fixing2,
																			  LPXLOPER XL_slopeFlag,
																			  LPXLOPER XL_cptStrikeMethod);
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_QUANTOSPREADCORRIDOR(	LPXLOPER XL_startDate,
																		LPXLOPER XL_endDate,
 																		LPXLOPER XL_capOrFloor,
																		LPXLOPER XL_strikes,
																		LPXLOPER XL_liborIdx1,
																		LPXLOPER XL_liborIdx2,
																		LPXLOPER XL_liborIdxP,
																		LPXLOPER XL_fixedRateP,
																		LPXLOPER XL_leg1Weights,
																		LPXLOPER XL_leg2Weights,
																		LPXLOPER XL_legPWeights,
																		LPXLOPER XL_leg1Fixings,
																		LPXLOPER XL_leg2Fixings,
																		LPXLOPER XL_legPFixings,
																		LPXLOPER XL_leg1Spread,
																		LPXLOPER XL_leg2Spread,
																		LPXLOPER XL_modelParams,
																		LPXLOPER XL_prodParams,
																		LPXLOPER XL_notional);

__declspec(dllexport) LPXLOPER WINAPI Local_CapLetPrice (LPXLOPER XL_secId,
														 LPXLOPER XL_modId,
														 LPXLOPER XL_numEx);

__declspec(dllexport) LPXLOPER WINAPI Local_RATCHET (LPXLOPER XL_swapLeg,
													 LPXLOPER XL_isItCapOrFloor,
													 LPXLOPER XL_strike,
													 LPXLOPER XL_spreadDates,
													 LPXLOPER XL_spreadValues,
													 LPXLOPER XL_correlDates,
													 LPXLOPER XL_correlValues,
													 LPXLOPER XL_fwdVolsDates,
													 LPXLOPER XL_fwdVolsValues);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_RATCHET (LPXLOPER XL_swapLeg,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_spreadDates,
														 LPXLOPER XL_spreadValues,
														 LPXLOPER XL_correlDates,
														 LPXLOPER XL_correlValues,
														 LPXLOPER XL_fwdVolsDates,
														 LPXLOPER XL_fwdVolsValues);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DIGITAL (LPXLOPER XL_swapLeg,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_spread1,
														 LPXLOPER XL_spread2,
														 LPXLOPER XL_payoff);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_DIGITAL (LPXLOPER XL_swapLeg,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_spread1,
															 LPXLOPER XL_spread2,
															 LPXLOPER XL_payoff);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_DUALCAP (LPXLOPER XL_startDate,
														 LPXLOPER XL_endDate,
														 LPXLOPER XL_isItCapOrFloor,
														 LPXLOPER XL_strike,
														 LPXLOPER XL_liborType1,
														 LPXLOPER XL_liborType2,
														 LPXLOPER XL_daycount,
														 LPXLOPER XL_Freq,
														 LPXLOPER XL_resetTiming1,
														 LPXLOPER XL_resetTiming2,
														 LPXLOPER XL_resetGap1,
														 LPXLOPER XL_resetGap2,
														 LPXLOPER XL_indexCcyId1,
														 LPXLOPER XL_indexCcyId2,
														 LPXLOPER XL_discountCcyId,
														 LPXLOPER XL_dresetCal,
														 LPXLOPER XL_fresetCal);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_DUALCAP (LPXLOPER XL_startDate,
															 LPXLOPER XL_endDate,
															 LPXLOPER XL_isItCapOrFloor,
															 LPXLOPER XL_strike,
															 LPXLOPER XL_liborType1,
															 LPXLOPER XL_liborType2,
															 LPXLOPER XL_daycount,
															 LPXLOPER XL_Freq,
															 LPXLOPER XL_resetTiming1,
															 LPXLOPER XL_resetTiming2,
															 LPXLOPER XL_resetGap1,
															 LPXLOPER XL_resetGap2,
															 LPXLOPER XL_indexCcyId1,
															 LPXLOPER XL_indexCcyId2,
															 LPXLOPER XL_discountCcyId,
															 LPXLOPER XL_dresetCal,
															 LPXLOPER XL_fresetCal);

__declspec(dllexport) LPXLOPER WINAPI Local_GlobalCap ( LPXLOPER XL_swapLeg,
														LPXLOPER XL_isItCapOrFloor,
														LPXLOPER XL_globalCapValue,
														LPXLOPER XL_globalCapSpreads,
														LPXLOPER XL_globalCapFixedRates,
														LPXLOPER XL_globalCapBarriers,
														LPXLOPER XL_finalRatio,
														LPXLOPER XL_MCNbIter,
														LPXLOPER XL_globalCapPastFixings);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GlobalCap ( LPXLOPER XL_swapLeg,
															LPXLOPER XL_isItCapOrFloor,
															LPXLOPER XL_globalCapValue,
															LPXLOPER XL_globalCapSpreads,
															LPXLOPER XL_globalCapFixedRates,
															LPXLOPER XL_globalCapBarriers,
															LPXLOPER XL_finalRatio,
															LPXLOPER XL_MCNbIter,
															LPXLOPER XL_globalCapPastFixings);

#endif /* ARM_XL_CAPFL_LOCAL_H */