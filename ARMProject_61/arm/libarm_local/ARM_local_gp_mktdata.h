/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_mktdata.h,v $
 * Revision 1.1  2004/03/10 15:08:43  ebenhamou
 * Initial version
 *
 */

/*! \file ARM_local_gp_mktdata.h
 *
 *  \brief file for the various market data object in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date March 2004
 */

#ifndef ARMLOCAL_GP_MKTDATA_H
#define ARMLOCAL_GP_MKTDATA_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"

/// forward declaration
class ARM_Date;

extern long ARMLOCAL_MktDataManager_Create(
    const ARM_Date& asOf, 
    ARM_result&	result, 
    long			objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_RegisterZCCurve(
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	long zcCurveId,
	ARM_result&	result,
    long objId = ARM_NULL_OBJECT_ID);


extern long ARMLOCAL_MktDataManager_GetZCCurve(
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);


extern long ARMLOCAL_MktDataManager_GetZCCurveKey(
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& source, 
	ARM_result&	result ); 


extern long ARMLOCAL_MktDataManager_RegisterVolCurve( 
	long mktDataManagerId,
	const string& indexName, 
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType, 
	const string& volType, 
	const string& source, 
	long volCurveId,
	ARM_result&	result,
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_GetVolCurve( 
	long mktDataManagerId,
	const string& indexName,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOf,
	const string& volMktType,
	const string& volType,
	const string& source,
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_GetVolCurveKey( 
	const string& indexName,
	const string& ccy,
	const string& cvName,
	const ARM_Date& asOf,
	const string& volMktType,
	const string& volType,
	const string& source,
	ARM_result&	result );


extern long ARMLOCAL_MktDataManager_RegisterVolMktModel( 
	long mktDataManagerId,
	const string& indexName,
	const string& ccy, 
	const string& cvName, 
	const ARM_Date& asOf, 
	const string& volMktType,
	const string& source, 
	long volMktModelId,
	ARM_result&	result,
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_GetVolMktModel(
	long mktDataManagerId,
	const string& indexName,
	const string& ccy,
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,
	const string& source,
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_GetVolMktModelKey(
	const string& indexName,
	const string& ccy,
	const string& cvName, 
	const ARM_Date& asOf,
	const string& volMktType,
	const string& source,
	ARM_result&	result );

extern long ARMLOCAL_MktDataManager_RegisterData( 
	long mktDataManagerId,
	const string& objectKey, 
	long objectId,
    bool isFill,
	ARM_result&	result,
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_GetData(
	long mktDataManagerId,
	const string& objectKey,
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);

extern long ARMLOCAL_MktDataManager_SetDetailMode(
	long mktDataManagerId,
	bool detailModeFlag,
	ARM_result&	result );


////////////////////////////////////////////
//// Function to change the date on a mkt data manager
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ChangeDate(
	long mktDataManagerId,
	const ARM_Date& newDate,
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to reset all data
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ResetAllData(
	long mktDataManagerId,
	ARM_result&	result );


////////////////////////////////////////////
//// Function to reset only the current mkt data manager rep
////////////////////////////////////////////
extern long ARMLOCAL_MktDataManager_ResetMyData(
	long mktDataManagerId,
	ARM_result&	result, 
    long objId = ARM_NULL_OBJECT_ID);


#endif

