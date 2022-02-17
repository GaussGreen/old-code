/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_inflation_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_inflation_local.h
 *
 *  \brief file function of the inflation
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARM_XL_GP_INFLATION_LOCAL_H
#define ARM_XL_GP_INFLATION_LOCAL_H

#include <ARM\libarm_local\firstToBeIncluded.h>

#include <ARM\libarm_local\ARM_local_class.h>
#include <ARM\libarm_local\ARM_local_glob.h>

#include "ARM_local_interglob.h"
#include "ARM_local_interface.h"


//////////////////////////////////
/// 
/// In order to factorise code
/// many function are calling
/// common function
/// these functions are not
/// declared here 
/// 
/// only exported functions are 
/// included here to facilitate 
/// the reading of the header
/// 
//////////////////////////////////

///////////////////////////////////
/// Creates an inflation curve
/// Handles the case of 
///	previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CreateInfCurv(
	 LPXLOPER XL_asOfDate,
	 LPXLOPER XL_indexName,
	 LPXLOPER XL_CPIIndexValue,
	 LPXLOPER XL_CPIIndexDate,
	 LPXLOPER XL_maturities,
	 LPXLOPER XL_values,
	 LPXLOPER XL_MonthlyInterpType,
	 LPXLOPER XL_DailyInterpType,
	 LPXLOPER XL_DCFMonthly,
	 LPXLOPER XL_DCFDaily,
	 LPXLOPER XL_ExtrapolType,
	 LPXLOPER XL_ResetManager,
	 LPXLOPER XL_SeasonManager);

///////////////////////////////////
/// Creates an inflation curve
/// version for non persistent In Xl
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateInfCurv(
	  LPXLOPER XL_asOfDate,
	  LPXLOPER XL_indexName,
	  LPXLOPER XL_CPIIndexValue,
	  LPXLOPER XL_CPIIndexDate,
	  LPXLOPER XL_maturities,
	  LPXLOPER XL_values,
	  LPXLOPER XL_MonthlyInterpType,
	  LPXLOPER XL_DailyInterpType,
	  LPXLOPER XL_DCFMonthly,
	  LPXLOPER XL_DCFDaily,
	  LPXLOPER XL_ExtrapolType,
	  LPXLOPER XL_ResetManager,
	  LPXLOPER XL_SeasonManager);

////////////////////////////
/// interpolation of the CPI
////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_CPIInterp(	
	LPXLOPER XL_curv,
	LPXLOPER XL_CPIDate,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_CPILag,
	LPXLOPER XL_weight );

////////////////////////////////
/// interpolation of the ZC Rate
////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_ZCRateInterp(	
	LPXLOPER XL_curv,
	LPXLOPER XL_CPIDate,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_DailyInterpType,
	LPXLOPER XL_CPILag,
	LPXLOPER XL_weight );


//////////////////////////////////
/// Creates ian inflation index
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfIdxCreate(
	LPXLOPER XL_indexName,
	LPXLOPER XL_resetLag,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_ccyId );

	
///////////////////////////////////
/// Creates ian inflation index
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfIdxCreate(
	LPXLOPER XL_indexName,
	LPXLOPER XL_resetLag,
	LPXLOPER XL_DCFLag,
	LPXLOPER XL_ccyId);


//////////////////////////////////
/// Creates an inflation leg very 
/// very generically using datestrip
/// objects
///
/// version that takes into account 
/// previous creation of object
//////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegwDateStripCreate(
    LPXLOPER XL_infIdx, 
	LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
    LPXLOPER XL_constant,
	LPXLOPER XL_notionalExchange,
	LPXLOPER XL_notionalExchangeType,
	LPXLOPER XL_numStripDateId,
	LPXLOPER XL_denomStripDateId );

///////////////////////////////////
/// Creates an inflation leg very 
/// very generically using datestrip
/// objects
///
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegwDateStripCreate(
    LPXLOPER XL_infIdx, 
	LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
    LPXLOPER XL_constant,
	LPXLOPER XL_notionalExchange,
	LPXLOPER XL_notionalExchangeType,
	LPXLOPER XL_numStripDateId,
	LPXLOPER XL_denomStripDateId );



///////////////////////////////////
/// creates an inflation leg 
/// year to year
///
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegYtYCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_Coleverage);



///////////////////////////////////
/// creates an inflation leg 
/// year to year
///
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegYtYCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate ,
	LPXLOPER XL_Coleverage);

///////////////////////////////////
/// creates an inflation model
/// takes an inflation and a yield curve model
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfYCMOD(
	LPXLOPER XL_zeroCurve,
	LPXLOPER XL_discCurve );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfYCMOD(
	 LPXLOPER XL_zeroCurve,
	 LPXLOPER XL_discCurve);

///////////////////////////////////
/// ZeroCoupon inflation leg
///////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_InfLegZCCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_firstReset );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegZCCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_adjFirstDate,
	LPXLOPER XL_firstReset );


///////////////////////////////////
/// creates an OAT type inflation leg 
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfLegOATCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_firstReset,
	LPXLOPER XL_test );



///////////////////////////////////
/// creates an OAT type inflation leg 
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfLegOATCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_finalNotionalType,
	LPXLOPER XL_firstReset );


///////////////////////////////////
/// creates a fixed zero Coupon leg
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_FixZC(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_fixRate, 
    LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_ccy);

///////////////////////////////////
/// creates a fixed zero Coupon leg
/// PXL version
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FixZC(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_fixRate, 
    LPXLOPER XL_rcvOrPay,
	LPXLOPER XL_dayCount,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar,
	LPXLOPER XL_ccy);


///////////////////////////////////
/// creates a reset manager object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ResetManager(
    LPXLOPER XL_data );

///////////////////////////////////
/// creates a reset manager object
/// PXL version
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ResetManager(
    LPXLOPER XL_data );


///////////////////////////////////
/// creates and fills a sparse vol cube
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_CreateNFill(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_lastKnonwnDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	LPXLOPER XL_strikeType,
	LPXLOPER XL_volType );

///////////////////////////////////
/// creates a sparse vol cube
/// PXL version
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SparseVolCube_CreateNFill(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_lastKnonwnDate,
	LPXLOPER XL_indexName,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols,
	LPXLOPER XL_strikeType,
	LPXLOPER XL_volType );

///////////////////////////////////
/// fills a sparse vol cube already 
/// created
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_SparseVolCube_Fill(
	LPXLOPER XL_sparseVolCubeId,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols );

///////////////////////////////////
/// fills a sparse vol cube already 
/// created. PXL version
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SparseVolCube_Fill(
	LPXLOPER XL_sparseVolCubeId,
	LPXLOPER XL_Dim1Type,
	LPXLOPER XL_Dim1Value,
	LPXLOPER XL_Dim2Type,
	LPXLOPER XL_Dim2Value,
	LPXLOPER XL_strikes,
	LPXLOPER XL_vols );




__declspec(dllexport) LPXLOPER WINAPI Local_VolCubeFromSparseVolCube(
	LPXLOPER XL_sparseVolCubeId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_VolCubeFromSparseVolCube(
	LPXLOPER XL_sparseVolCubeId );


__declspec(dllexport) LPXLOPER WINAPI Local_InfBSMOD(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolCurvId,
	LPXLOPER XL_CorrelManagerId,
	LPXLOPER XL_IRBSModelId,
	LPXLOPER XL_InfSwoptCurveId,
	LPXLOPER XL_IRSwoptCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfBSMOD(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolCurvId,
	LPXLOPER XL_CorrelManagerId,
	LPXLOPER XL_IRBSModelId,
	LPXLOPER XL_InfSwoptCurveId,
	LPXLOPER XL_IRSwoptCurveId );


//////////////////////////////////////////////////////
/// Addin to create an inflation cap floor
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCapFloorCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_swapType,
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar );


//////////////////////////////////////////////////////
/// Version for XL exportation of the fucntion to create
/// and inflation cap and floor in Visual Basic
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCapFloorCreate(
    LPXLOPER XL_startDate, 
	LPXLOPER XL_endDate,
	LPXLOPER XL_infIdx, 
	LPXLOPER XL_capOrFloor,
	LPXLOPER XL_strike,
	LPXLOPER XL_leverage,
	LPXLOPER XL_spread,
    LPXLOPER XL_swapType,
    LPXLOPER XL_rcvOrPay,
    LPXLOPER XL_interpType,
	LPXLOPER XL_multiple,
	LPXLOPER XL_constant,
    LPXLOPER XL_resetFreq,
	LPXLOPER XL_dayCount,
    LPXLOPER XL_resetCalendar,
	LPXLOPER XL_fwdRule,
	LPXLOPER XL_intRule,
	LPXLOPER XL_stubRule,
	LPXLOPER XL_resetNumGap,
	LPXLOPER XL_resetDenomGap,
	LPXLOPER XL_payFreq,
	LPXLOPER XL_payGap,
	LPXLOPER XL_payCalendar );



//////////////////////////////////////////////////////
/// Addin to create a correlation matrix
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelMat_Create(
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z );

//////////////////////////////////////////////////////
/// Addin to create a correlation matrix
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelMat_Create(
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z );


//////////////////////////////////////////////////////
/// Addin to create a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManager_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z );


//////////////////////////////////////////////////////
/// Addin to create a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManager_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_asOfDate,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z );


//////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl mat
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFromMat_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId );


//////////////////////////////////////////////////////
/// Addin to create a correlation manager from a correl mat
/// Version for VBA 
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManagerFromMat_Create(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId );


//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManager_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
    LPXLOPER XL_correlManagerId );



//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManager_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_Z,
    LPXLOPER XL_correlManagerId );


//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CorrelManagerFromMat_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
    LPXLOPER XL_correlManagerId );


//////////////////////////////////////////////////////
/// Addin to fill a correlation manager
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CorrelManagerFromMat_Fill(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_correlMatId,
    LPXLOPER XL_correlManagerId );


//////////////////////////////////////////////////////
/// Addin to compute a correlation from a correlation manager
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ComputeCorrelFromManager(
    LPXLOPER XL_mktTag,
    LPXLOPER XL_intraMktTag,
    LPXLOPER XL_X,
    LPXLOPER XL_Y,
    LPXLOPER XL_correlManagerId );

//Mingzhi
////////////////////////////////////////////////////////////////////////////////////
/// Addin to create and fill in a correlation manager from a correl with the arrays
////////////////////////////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_CreateGenCorrelManager( LPXLOPER XL_mktTags,
																	LPXLOPER XL_intraMktTags,
																	LPXLOPER XL_correlCurveIds);



__declspec(dllexport) LPXLOPER WINAPI Local_PXL_CreateGenCorrelManager( LPXLOPER XL_mktTags,
																	    LPXLOPER XL_intraMktTags,
																		LPXLOPER XL_correlCurveIds);

//////////////////////////////////////////////////////
/// Addin to compute the correlation from volatilities
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Vol_to_Cor(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_MaturityV );

//////////////////////////////////////////////////////
/// Addin to compute the zero coupon volatilities 
/// from correlation and year to year volatilities
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_YtYCor_to_ZC(
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV, 
		LPXLOPER XL_MaturityV );


//////////////////////////////////////////////////////
/// Addin to compute the year to year volatilities 
/// from correlation and zero coupon volatilities
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ZCCor_to_YtY(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV );

//////////////////////////////////////////////////////
/// Addin to provide confidence interval for zero coupon
/// and year to year volatilities and correlations
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_Bounds(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_Choice );

//////////////////////////////////////////////////////
/// Addin to compute the correlation from volatilities
/// in the homogeneous cases
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_HmgVol_to_Cor(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length );

//////////////////////////////////////////////////////
/// Addin to compute the year to year volatilities 
/// from correlation and zero coupon volatilities
/// in the homogeneous cases
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_HmgZCCor_to_YtY(
		LPXLOPER XL_ZCVolV, 
		LPXLOPER XL_CorV,
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length );

//////////////////////////////////////////////////////
/// Addin to compute the zero coupon volatilities 
/// from correlation and year to year volatilities
/// in the homogeneous case
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_HmgYtYCor_to_ZC(
		LPXLOPER XL_YtYVolV, 
		LPXLOPER XL_CorV, 
		LPXLOPER XL_MaturityV,
		LPXLOPER XL_length);


/////////////////////////////////////////////////////
/// Function to compute the inflation swap's 
/// implied volatility from the year-on-year ones
/////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_VolYoY_to_VolSwp(
		LPXLOPER XL_DFactor, 
		LPXLOPER XL_FwdCPI, 
		LPXLOPER XL_Vol_DF,
		LPXLOPER XL_Vol_YoY, 
		LPXLOPER XL_AvgCor, 
		LPXLOPER XL_Dates, 
		LPXLOPER XL_Tenors, 
		LPXLOPER XL_SwpRate );

//////////////////////////////////////////////////////
/// Addin to create a seasonality manager
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_SeasonalityManager_Create(
    LPXLOPER XL_monthList,
    LPXLOPER XL_seasonSpreadList,
	LPXLOPER XL_seasonHorizonList,
    LPXLOPER XL_seasonAdjMode );

//////////////////////////////////////////////////////
/// Addin to create a seasonality manager
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_SeasonalityManager_Create(
    LPXLOPER XL_monthList,
    LPXLOPER XL_seasonSpreadList,
	LPXLOPER XL_seasonHorizonList,
    LPXLOPER XL_seasonAdjMode );


//////////////////////////////////////////////////////
/// Addin to create an inflation curve from Summit
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_GetInfZcFromSummit_Create(
    LPXLOPER XL_index,
    LPXLOPER XL_ccy,
    LPXLOPER XL_cvname,
    LPXLOPER XL_date,
	LPXLOPER XL_seasonAdj,
	LPXLOPER XL_seasonAdjMode);

//////////////////////////////////////////////////////
/// Addin to create an inflation curve from Summit
/// Version for VBA
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetInfZcFromSummit_Create(
    LPXLOPER XL_index,
    LPXLOPER XL_ccy,
    LPXLOPER XL_cvname,
    LPXLOPER XL_date,
	LPXLOPER XL_seasonAdj,
	LPXLOPER XL_seasonAdjMode);



//////////////////////////////////////////////////////
/// Addin for the livreatA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LIVRETACURVE (LPXLOPER XL_asOfDate,
															  LPXLOPER XL_infCurvId,
															  LPXLOPER XL_euribCurvId,
															  LPXLOPER XL_flagRounding,
															  LPXLOPER XL_resetManagerId,
															  LPXLOPER XL_fixingLivretAId,
															  LPXLOPER XL_fixingEuriborId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_LivreACurveGetRateDate(LPXLOPER XL_livretACurvId,
																	   LPXLOPER XL_dateIn);

//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCurveFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries );

//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfSwoVolCurveFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries );

//YK
//////////////////////////////////////////////////////
/// Addin to create an inflation volatility swaption from a model
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfSwoVolCubeFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries, 
	LPXLOPER XL_smiledTenors, 
	LPXLOPER XL_strikes);

//////////////////////////////////////////////////////
/// Addin to create an inflation volatility Cube swaption from a model
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfSwoVolCubeFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries,
	LPXLOPER XL_smiledTenors, 
	LPXLOPER XL_strikes);
//YK
//////////////////////////////////////////////////////
/// Addin to create an OAT inflation volatility swaption from a model
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfOATSwoVolCurveFromModel_Create(
	LPXLOPER XL_asOfDate,
	LPXLOPER XL_InfIRModelId,
	LPXLOPER XL_tenors,
	LPXLOPER XL_expiries, 
	LPXLOPER XL_coupon,
	LPXLOPER XL_choice );

//////////////////////////////////////////////////////
/// Addin to create an inflation multi bs model
/// Version for VBA
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfMultiBSMOD(
	LPXLOPER XL_infMultiBSModIdVec );


////////////////////////////////////////////////
/// Version for XL exportation 
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfMultiBSMOD(
	LPXLOPER XL_infMultiBSModIdVec );



//////////////////////////////////////////////////////
/// Addin to set a reset manager on an inflation curve
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_SetResetManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_ResetManagerId );


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCurv_SetResetManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_ResetManagerId );

//////////////////////////////////////////////////////
/// Addin to set a seasonality manager on an inflation curve
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_InfCurv_SetSeasonalityManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_SeasonalityManagerId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfCurv_SetSeasonalityManager(
	LPXLOPER XL_infCurvId,
	LPXLOPER XL_SeasonalityManagerId );

//////////////////////////////////////////////////////
/// Addin to get a seasonality manager from summit
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GetSeasonMgrFromSummit(LPXLOPER XL_index,
																	   LPXLOPER XL_ccy,
																	   LPXLOPER XL_cvname,
																	   LPXLOPER XL_asof,
																	   LPXLOPER XL_mode);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GetSeasonMgrFromSummit(LPXLOPER XL_index,
																		   LPXLOPER XL_ccy,
																		   LPXLOPER XL_cvname,
																		   LPXLOPER XL_asof,
																		   LPXLOPER XL_mode);

//////////////////////////////////////////////////////
/// Addin to create an INFCAPFLOOR on two legs
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFCAPFLOOR(LPXLOPER XL_swapId,
															   LPXLOPER XL_CAPOrFLOOR,
															   LPXLOPER XL_strike);
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFCAPFLOOR(LPXLOPER XL_swapId,
																   LPXLOPER XL_CAPOrFLOOR,
																   LPXLOPER XL_strike);

//////////////////////////////////////////////////////
/// Addin to create an INFDIGITAL Call Spread on two legs
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFCALLSPREAD(	LPXLOPER XL_swapId,
																	LPXLOPER XL_CAPOrFLOOR,
																	LPXLOPER XL_strike);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFCALLSPREAD(	LPXLOPER XL_swapId,
																		LPXLOPER XL_CAPOrFLOOR,
																		LPXLOPER XL_strike);

//////////////////////////////////////////////////////
/// Addin to create an INFCAPFLOOR on two legs
//////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ARM_GP_INFDIGITAL(LPXLOPER XL_PaySwapId,
															  LPXLOPER XL_DigitSwapId,
															  LPXLOPER XL_PayOffType,
															  LPXLOPER XL_Barrier,
															  LPXLOPER XL_CapOrFloor,
															  LPXLOPER XL_PayOrRec);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_GP_INFDIGITAL(LPXLOPER XL_PaySwapId,
																  LPXLOPER XL_DigitSwapId,
																  LPXLOPER XL_PayOffType,
																  LPXLOPER XL_Barrier,
																  LPXLOPER XL_CapOrFloor,
																  LPXLOPER XL_PayOrRec);

//////////////////////////////////////////////////////
/// Addin to create Hybrid Inf Model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIrMkt_Create		(	LPXLOPER XL_Asof,	
																					LPXLOPER XL_Keys,
																					LPXLOPER XL_ModId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIrMkt_Create	(	LPXLOPER XL_Asof,		
																					LPXLOPER XL_Keys,
																					LPXLOPER XL_ModId);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIr_Load(				LPXLOPER XL_Ins,		
																					LPXLOPER XL_Mkt,
																					LPXLOPER XL_Mod,
																					LPXLOPER XL_Pay);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIr_Load(			LPXLOPER XL_Ins,		
																					LPXLOPER XL_Mkt,
																					LPXLOPER XL_Mod,
																					LPXLOPER XL_Pay);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_HybridInfIr_PayOff(				LPXLOPER XL_CstCpnCoef,
																					LPXLOPER XL_CstOptCoef,

																					LPXLOPER XL_MainCpnName,
																					LPXLOPER XL_MainCpnCoef,
																					LPXLOPER XL_MainOptName,
																					LPXLOPER XL_MainOptCoef,

																					LPXLOPER XL_SubCpnName,
																					LPXLOPER XL_SubCpnCoef,
																					LPXLOPER XL_SubOptName,
																					LPXLOPER XL_SubOptCoef,

																					LPXLOPER XL_SupCpnName,
																					LPXLOPER XL_SupCpnCoef,
																					LPXLOPER XL_SupOptName,
																					LPXLOPER XL_SupOptCoef);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_HybridInfIr_PayOff(			LPXLOPER XL_CstCpnCoef,
																					LPXLOPER XL_CstOptCoef,

																					LPXLOPER XL_MainCpnName,
																					LPXLOPER XL_MainCpnCoef,
																					LPXLOPER XL_MainOptName,
																					LPXLOPER XL_MainOptCoef,

																					LPXLOPER XL_SubCpnName,
																					LPXLOPER XL_SubCpnCoef,
																					LPXLOPER XL_SubOptName,
																					LPXLOPER XL_SubOptCoef,

																					LPXLOPER XL_SupCpnName,
																					LPXLOPER XL_SupCpnCoef,
																					LPXLOPER XL_SupOptName,
																					LPXLOPER XL_SupOptCoef);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_SpreadCap_Create(				LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_Strike,
																					LPXLOPER XL_Notional);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_SpreadCap_Create(			LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_Strike,
																					LPXLOPER XL_Notional);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_Corridor_Create(				LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_StrikeInf,
																					LPXLOPER XL_StrikeSup,
																					LPXLOPER XL_Notional);


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_Corridor_Create(			LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_StrikeInf,
																					LPXLOPER XL_StrikeSup,
																					LPXLOPER XL_Notional);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_SpreadDigital_Create(			LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_Strike,
																					LPXLOPER XL_Notional);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_SpreadDigital_Create(		LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_Strike,
																					LPXLOPER XL_Notional);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_DoubleDigital_Create(			LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_MainStrike,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_SubStrike,
																					LPXLOPER XL_Notional);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_DoubleDigital_Create(		LPXLOPER XL_MainIndex,
																					LPXLOPER XL_MainType,
																					LPXLOPER XL_MainLeverage,
																					LPXLOPER XL_MainStrike,
																					LPXLOPER XL_SubIndex,
																					LPXLOPER XL_SubType,
																					LPXLOPER XL_SubLeverage,
																					LPXLOPER XL_SubStrike,
																					LPXLOPER XL_Spread,
																					LPXLOPER XL_Notional);


__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetPrice(						LPXLOPER XL_InfPricer,
															 						LPXLOPER XL_Key);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetPrice(					LPXLOPER XL_InfPricer,
																					LPXLOPER XL_Key);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetSchedule(					LPXLOPER XL_InfLeg,
															 						LPXLOPER XL_Key);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetSchedule(				LPXLOPER XL_InfLeg,
																					LPXLOPER XL_Key);

//////////////////////////////////////////////////////
/// Addin to create  Inf BS SMILED Model
//////////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_InfBSSmiledModel(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolSigmaId,
	LPXLOPER XL_VolNuId,
	LPXLOPER XL_VolRhoId,
	LPXLOPER XL_VolBetaId,
	LPXLOPER XL_VolAtmIrId,
	LPXLOPER XL_CorrelId,
	LPXLOPER XL_CorrelAdjId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_InfBSSmiledModel(
	LPXLOPER XL_AsOfDate,
	LPXLOPER XL_DiscountCurvId,
	LPXLOPER XL_InfFwdCurvId,
	LPXLOPER XL_VolSigmaId,
	LPXLOPER XL_VolNuId,
	LPXLOPER XL_VolRhoId,
	LPXLOPER XL_VolBetaId,
	LPXLOPER XL_VolAtmIrId,
	LPXLOPER XL_CorrelId,
	LPXLOPER XL_CorrelAdjId);

__declspec(dllexport) LPXLOPER Local_ARM_INF_EQHWSV_Laplace(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_evalTime,
	LPXLOPER XL_startTime,
	LPXLOPER XL_endTime,
	LPXLOPER XL_xt,
	LPXLOPER XL_vt,
	LPXLOPER XL_k_real,
	LPXLOPER XL_k_imag,
	LPXLOPER XL_isReal);

__declspec(dllexport) LPXLOPER Local_ARM_INF_EQHWSV_Density(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_evalTime,
	LPXLOPER XL_startTime,
	LPXLOPER XL_endTime,
	LPXLOPER XL_xt,
	LPXLOPER XL_vt,
	LPXLOPER XL_x,
	LPXLOPER XL_period,	
	LPXLOPER XL_frequency);

__declspec(dllexport) LPXLOPER WINAPI Local_ARM_INF_GetAdjCorrel(
		LPXLOPER XL_InfBsSmiledModel);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ARM_INF_GetAdjCorrel(
		LPXLOPER XL_InfBsSmiledModel);

#endif