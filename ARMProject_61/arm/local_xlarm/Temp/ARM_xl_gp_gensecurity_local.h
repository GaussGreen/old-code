/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_gensecurity_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_gensecurity_local.h
 *
 *  \brief file for the generic security part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARM_XL_GP_GENSECURITY_LOCAL_H
#define ARM_XL_GP_GENSECURITY_LOCAL_H

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
/// Creates a deal description object
/// Handles the case of 
///	previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DealDes_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PricedColumns);

///////////////////////////////////
/// Creates a deal description object
/// version for non persistent In Xl
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DealDes_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PricedColumns);


///////////////////////////////////
/// Creates a generic security object
/// Handles the case of 
///	previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PayCurveName,
	LPXLOPER XL_ExercBoundaryResetFlag,
    LPXLOPER XL_OtherPayoffsFlag,
	LPXLOPER XL_PricedColumns);

///////////////////////////////////
/// Creates a generic security object
/// version for non persistent In Xl
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_Create(
	LPXLOPER XL_Data1,
	LPXLOPER XL_Data2,
	LPXLOPER XL_PayCurveName,
	LPXLOPER XL_ExercBoundaryResetFlag,
    LPXLOPER XL_OtherPayoffsFlag,
	LPXLOPER XL_PricedColumns);

///////////////////////////////////
/// Creates a generic security object
/// from Another Generic security
///	(changes a bermudan option into a trigger
/// using its exercise boundary)
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_ChangeAmericanIntoTrigger(
	LPXLOPER XL_GenSecurityId, LPXLOPER XL_PricingModelId );

///////////////////////////////////
/// Creates a generic security object
/// from another generic security
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_ChangeAmericanIntoTrigger(
	LPXLOPER XL_GenSecurityId, LPXLOPER XL_PricingModelId );

///////////////////////////////////
/// Creates a generic security object
/// from Another Generic security
///	(changes a MAX(PV) keywords into exercise
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_ChangeMAXPVIntoExercise(
	LPXLOPER XL_GenSecurityId );

///////////////////////////////////
/// Creates a generic security object
/// from another generic security
///	used for call in VBA
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_ChangeMAXPVIntoExercise(
	LPXLOPER XL_GenSecurityId );

///////////////////////////////////
/// Set the parse Flag on and off
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_SetPTFlag(
	LPXLOPER XL_GenSec,
	LPXLOPER XL_PTFlag );

						 
///////////////////////////////////
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GramHelper_Create(
	LPXLOPER XL_FuncName );


///////////////////////////////////
/// version for non persistentInXl
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GramHelper_Create(
	LPXLOPER XL_FuncName );

///////////////////////////////////
/// Set Look for Circular ref on and off
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_SetCRFlag(
	LPXLOPER XL_CRFlag );

///////////////////////////////////
/// Extract sub deal des from a gen sec
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_DealDes_ExtractSubDealDes(
	LPXLOPER XL_genSecId,
	LPXLOPER XL_cutoffColName,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_discounting );

///////////////////////////////////
/// Extract sub deal des from a gen sec (non persistent version)
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_DealDes_ExtractSubDealDes(
	LPXLOPER XL_genSecId,
	LPXLOPER XL_cutoffColName,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_discounting );


///////////////////////////////////
/// Get the cst manager of a gen sec
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenSec_GetCstManager(
	LPXLOPER XL_GenSec);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenSec_GetCstManager(
	LPXLOPER XL_GenSec);

#endif