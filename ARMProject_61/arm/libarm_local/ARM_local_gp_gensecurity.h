/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_gensecurity.h,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


/*! \file ARM_local_gp_gensecurity.h
 *
 *  \brief file for the generic security part in the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARMLOCAL_GP_GENSECURITY_H
#define ARMLOCAL_GP_GENSECURITY_H

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_result.h"
#include <GP_Base\gpbase\valuetype.h>

using ARM::ARM_GP_VALUE_TYPE;

////////////////////////////////////////////
//// Function to create a deal description
////////////////////////////////////////////
extern long ARMLOCAL_DealDes_Create(
	const VECTOR<string>&	vecCCString,
	const VECTOR<ARM_GP_VALUE_TYPE>&	types,
	long rowsNb,
	long colsNb,
    const string& payCurveName,
	long cstManagerId,
	bool ExercBoundaryResetFlag, 
	bool otherPayoffsFlag, 
	const VECTOR<string>&	pricedColumns,
	bool ivFlag,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to create a generic security
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_Create(
	const VECTOR<string>&	vecCCString,
	const VECTOR<ARM_GP_VALUE_TYPE>& types,
	long rowsNb,
	long colsNb,
    const string& payCurveName,
	long cstManagerId,
	bool ExercBoundaryResetFlag,
    bool otherPayoffsFlag,
	const VECTOR<string>&	pricedColumns,
	bool ivFlag,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a new generic security
//// from an old one
////////////////////////////////////////////

extern long ARMLOCAL_GenSec_ChangeAmericanIntoTrigger_Common(
	const long& GenSecId, const long& PricingModelId, 
	ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to create a new generic security
//// from an old one
////////////////////////////////////////////

extern long ARMLOCAL_GenSec_ChangeMAXPVIntoExercise_Common(
	const long& GenSecId, 
	ARM_result&	result, 
	long objId = ARM_NULL_OBJECT_ID);

////////////////////////////////////////////
//// Function to set parse tree on and off
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_SetParseTreeFlag(
	long GenSecId,
	bool ParseTreeFlag,
	ARM_result& result );


////////////////////////////////////////////
//// Function to get help on generic pricer
////////////////////////////////////////////
extern long ARMLOCAL_GramHelper(
	const string& FuncName,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to set look for circular reference on and off
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_SetCircularRefFlag(
	bool CircularRefFlag,
	ARM_result&	result );


////////////////////////////////////////////
//// Function to get the deal description table
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_GetDealDesTable(
	    long GenSecId,
		VECTOR< string >& dealDesText,
		VECTOR< ARM_GP_VALUE_TYPE >& dealDesFormat,
        long& nbRows,
        long& nbCols,
		ARM_result&	result );


////////////////////////////////////////////
//// Function to create a cst manager
////////////////////////////////////////////
extern long ARMLOCAL_CstManager_Create(
	const VECTOR<string>&	cstNames,
	const VECTOR<double>&	values,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to get a cst manager from a GenSec
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_GetCstManager(
	    const long& GenSecId,
	    ARM_result& result, 
	    long objId = ARM_NULL_OBJECT_ID );

////////////////////////////////////////////
//// Function to create an obj manager (cst manager)
////////////////////////////////////////////
extern long ARMLOCAL_ObjManager_Create(
	const VECTOR<string>&	cstNames,
	const VECTOR<long>&    objIds,
	ARM_result& result, 
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to extract a sub deal description in the gen security
////////////////////////////////////////////
extern long ARMLOCAL_GenSec_ExtractSubDealDes(
	const long& GenSecId,
	const string& ColName,
	const vector<string>& otherCols,
	const string& payCurveName,
	ARM_result&	result,
	long objId = ARM_NULL_OBJECT_ID);

#endif