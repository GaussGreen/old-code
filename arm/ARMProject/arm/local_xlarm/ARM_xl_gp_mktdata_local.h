/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_xl_gp_mktdata_local.h,v $
 * Revision 1.1  2004/01/13 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_xl_gp_mktdata_local.h
 *
 *  \brief file for the mkt data part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARM_XL_GP_MKTDATA_LOCAL_H
#define ARM_XL_GP_MKTDATA_LOCAL_H


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

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_Create(	
	LPXLOPER XL_asOf );
__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_Create(
	LPXLOPER XL_asOf );


///////////////////////////////////
/// ZC Curve addins
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ZCCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_zcCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ZCCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_zcCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ZCCurveGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source );

///////////////////////////////////
/// Vol Curve addins
///////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolCurveGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_volCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolCurveSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType,
	LPXLOPER XL_volCurveId );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolCurveGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volType );

////////////////////////////////////////////////
/// Vol Mkt Model
////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolMktModelGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_source );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volMktModelId );


__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_VolMktModelSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType,
	LPXLOPER XL_volMktModelId );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_VolMktModelGetKey(
	LPXLOPER XL_indexName,
	LPXLOPER XL_ccy,
	LPXLOPER XL_cvName,
	LPXLOPER XL_asOf,
	LPXLOPER XL_source,
	LPXLOPER XL_volMktType );


////////////////////////////////////////////////
/// Mkt Data general
////////////////////////////////////////////////

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_MktDataGet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_MktDataSet(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId );

__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_MktDataFill(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_objectKey,
	LPXLOPER XL_objectId );

/////////////////////////////////////////////////////////////
/// Set Detail mode
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_SetDetailMode(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_detailMode );

/////////////////////////////////////////////////////////////
/// function to give ability to change date and data
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ResetAllData(
	LPXLOPER XL_mktDataManagerId );

/////////////////////////////////////////////////////////////
/// function to give ability to change date and data
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ResetMyData(
	LPXLOPER XL_mktDataManagerId );

/////////////////////////////////////////////////////////////
/// function to give ability to change date and data
/////////////////////////////////////////////////////////////
__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_PXL_ResetMyData(
	LPXLOPER XL_mktDataManagerId );


__declspec(dllexport) LPXLOPER WINAPI Local_MktDataManager_ChangeAsOf(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_asOfData );

__declspec(dllexport) LPXLOPER WINAPI Local_PXL_MktDataManager_ChangeAsOf(
	LPXLOPER XL_mktDataManagerId,
	LPXLOPER XL_asOfData );


#endif