/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_local_gp_pricer.h
 *
 *  \brief file for the pricer part of the generic pricer
 *
 *	\author  E Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#ifndef ARMLOCAL_GP_PRICER_H
#define ARMLOCAL_GP_PRICER_H

#include "firstToBeIncluded.h"
#include "ARM_result.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <GP_Infra\gpinfra\gramfunctorarg.h>
using ARM::ARM_GramFctorArg;

////////////////////////////////////////////
//// Function to create a backward induction
///		numerical method
////////////////////////////////////////////
extern long ARMLOCAL_GenPricer_Create(
	long genSecurityId,
	long pricingModelId,
	const vector<string>& columnNames,
	const vector<double>& columnPrices,
	const string& refColumn,
	const vector<double>& betas,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID);


////////////////////////////////////////////
//// Function to Get the data from the generic pricer
////////////////////////////////////////////
long ARMLOCAL_GenPricer_GetData(
	long genPricerId,
	const string& key,
	const string& columnName,
	ARM_GramFctorArg& argResult,
	ARM_result& result );


////////////////////////////////////////////
//// Function to set detail mode on pricer
////////////////////////////////////////////
extern long ARMLOCAL_GenPricer_SetDetailMode(
	long genPricerId,
	bool detailMode,
	ARM_result& result );

////////////////////////////////////////////
//// Function to change a security into a GenSec
////////////////////////////////////////////
extern long ARMLOCAL_ChangeSecurityIntoGenSec(
	const long& securityId,
	const double& AsofDate,
	const string& modelName,
	ARM_result&	result, 
	long objId )

#endif