/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_pricer_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_pricer_local.h
 *
 *  \brief file for the pricer part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */


#ifndef ARM_XL_GP_PRICER_LOCAL_H
#define ARM_XL_GP_PRICER_LOCAL_H

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
/// Creates an IR Fwd Model
/// version that takes into account 
/// previous creation of object
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_GenPricer_Create(
	LPXLOPER XL_GenSecurityId,
	LPXLOPER XL_PricingModelId,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_columnPrices,
	LPXLOPER XL_refColumn,
	LPXLOPER XL_Beta);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GenPricer_Create(
	LPXLOPER XL_GenSecurityId,
	LPXLOPER XL_PricingModelId,
	LPXLOPER XL_columnNames,
	LPXLOPER XL_columnPrices,
	LPXLOPER XL_refColumn,
	LPXLOPER XL_Beta);

__declspec(dllexport) LPXLOPER WINAPI Local_GetData_FromGenPricer(
	LPXLOPER XL_PricerId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_GetData_FromGenPricer(
	LPXLOPER XL_PricerId,
    LPXLOPER XL_KeyId );

__declspec(dllexport) LPXLOPER WINAPI Local_GenPricer_SetDetailMode(
	LPXLOPER XL_PricerId,
	LPXLOPER XL_DetailMode );

__declspec(dllexport) LPXLOPER WINAPI Local_ChangeSecurityIntoGenSec(
	LPXLOPER XL_SecurityId, LPXLOPER XL_AsOfDate, LPXLOPER XL_ModelName);

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_ChangeSecurityIntoGenSec(
	LPXLOPER XL_SecurityId, LPXLOPER XL_AsOfDate, LPXLOPER XL_ModelName);


#endif

