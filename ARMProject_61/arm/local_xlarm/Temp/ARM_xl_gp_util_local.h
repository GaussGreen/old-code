/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_xl_gp_util_local.h
 *
 *  \brief file for the various utilities in the generic pricer
 *
 *	\author  J M Prie
 *	\version 1.0
 *	\date January 2004
 */

#ifndef ARM_XL_GP_UTIL_LOCAL_H
#define ARM_XL_GP_UTIL_LOCAL_H

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
/// Trigo Matrix Function
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_TrigoMatrix(
	LPXLOPER XL_Alpha,
    LPXLOPER XL_N);


///////////////////////////////////
/// Local Covariance Computation
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LocalCovariance(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromTime,
    LPXLOPER XL_ToTime,
    LPXLOPER XL_StartTime1,
    LPXLOPER XL_EndTime1,
    LPXLOPER XL_StartTime2,
    LPXLOPER XL_EndTime2,
    LPXLOPER XL_StartTime3,
    LPXLOPER XL_EndTime3,
    LPXLOPER XL_StartTime4,
    LPXLOPER XL_EndTime4);

///////////////////////////////////
/// Local Correlation Computation
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_LocalCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_UnderlyingType,
    LPXLOPER XL_FromTime,
    LPXLOPER XL_ToTime,
    LPXLOPER XL_StartTime1,
    LPXLOPER XL_EndTime1,
    LPXLOPER XL_StartTime2,
    LPXLOPER XL_EndTime2,
    LPXLOPER XL_StartTime3,
    LPXLOPER XL_EndTime3,
    LPXLOPER XL_StartTime4,
    LPXLOPER XL_EndTime4);




////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_Create(
	);

////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_EventViewer_Create(
	);

////////////////////////////////////////////
//// function to reset message in the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_ResetMssg(
	LPXLOPER XL_EventViewerId );

////////////////////////////////////////////
//// function to set verbose mode in the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_EventViewer_SetVerboseMode(
	LPXLOPER XL_verboseMode );


////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_ErrViewer_Create();

////////////////////////////////////////////
//// function to create a representant of the event viewer
////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ErrViewer_Create();


///----------------------------------------------
///----------------------------------------------
///             GP Matrix creation
/// Inputs : Values
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GPMatrix_Create( LPXLOPER XL_values );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GPMatrix_Create( LPXLOPER XL_values );

///----------------------------------------------
///----------------------------------------------
///             GP Vector creation
/// Inputs : Values
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GPVector_Create( LPXLOPER XL_values,LPXLOPER Xl_type );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GPVector_Create( LPXLOPER XL_values,LPXLOPER Xl_type );

///----------------------------------------------
///----------------------------------------------
///             Regression Calculation
/// Inputs : X, Y
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GP_LeastSquareRegression( 
	LPXLOPER XL_X,
	LPXLOPER XL_Y);


///----------------------------------------------
///----------------------------------------------
///             ACP Calculation
/// Inputs : Matrix
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_GP_ACP( 
	LPXLOPER XL_Matrix);

///----------------------------------------------
///----------------------------------------------
///             Local Covariance & Correlation Functions
/// Inputs : 
///----------------------------------------------
///----------------------------------------------

__declspec(dllexport) LPXLOPER WINAPI Local_IntegratedCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Tenor1,
	LPXLOPER XL_Tenors,
	LPXLOPER XL_Expiries );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_IntegratedCorrelation(
	LPXLOPER XL_ModelId,
	LPXLOPER XL_Tenor1,
	LPXLOPER XL_Tenors,
	LPXLOPER XL_Expiries );

///----------------------------------------------
///----------------------------------------------
/// To convert Fx Mkt Data to Totem Format
/// Inputs :  look at .cpp
///----------------------------------------------
///----------------------------------------------
__declspec(dllexport) LPXLOPER WINAPI Local_FxMktToTotemCalibrate(
	LPXLOPER XL_portfolioId,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_DeltaCalls,
	LPXLOPER XL_DeltaPuts,
	LPXLOPER XL_Model,
	LPXLOPER XL_InitPoint,
	LPXLOPER XL_LowerBound,
	LPXLOPER XL_UpperBound,
	LPXLOPER XL_MaxIter,
	LPXLOPER XL_AlgoType,
	LPXLOPER XL_Tolerance,
	LPXLOPER XL_OutPutId);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_FxMktToTotemCalibrate(
	LPXLOPER XL_portfolioId,
	LPXLOPER XL_ATMVol,
	LPXLOPER XL_DeltaCalls,
	LPXLOPER XL_DeltaPuts,
	LPXLOPER XL_Model,
	LPXLOPER XL_InitPoint,
	LPXLOPER XL_LowerBound,
	LPXLOPER XL_UpperBound,
	LPXLOPER XL_MaxIter,
	LPXLOPER XL_AlgoType,
	LPXLOPER XL_Tolerance,
	LPXLOPER XL_OutPutId);

#endif