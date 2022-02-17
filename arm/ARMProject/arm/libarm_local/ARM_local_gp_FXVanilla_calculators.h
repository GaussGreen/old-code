/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_FXVANILLA_CALCULATORS_H
#define ARMLOCAL_GP_FXVANILLA_CALCULATORS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;


class ARM_FXVanillaCalculator_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_FXVanillaCalculator_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_FXVanillaCalculator_InitFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_FXVanillaCalculator_InitFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


class ARM_FXRACalculator_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_FXRACalculator_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_FXRACalculator_InitFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_FXRACalculator_InitFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


/****************************************************************
	QUANTO
*****************************************************************/
long ARMLOCAL_CRAQuantoCalculator_Create(
		const ARM_Currency& ccyDom,
		const ARM_Currency& ccyFor,
		const double& startDate,
		const double& endDate,
		const int& payReceive,
		const int& callFreq,
		const int& callNotice,
		const string& callCal,
		const int& fundFreq,
		const int& fundDayCount,
		const int& cpnDayCount,
		const int& cpnPayFreq,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const int& cpnResetFreq,
		const int& cpnResetTiming,
		const int& refIndex,
		const int& payIndex,
		const int& payIndexResetTiming,
		const double& notional,
		const long& notionalId,
		const long& callFeesId,
		const double& fundSpread,
		const long& fundSpreadId,
		const double& boostedFix,
		const long& boostedFixId,
		const double& bDown,
		const long& bDownId,
		const double& bUp,
		const long& bUpId,
		const vector<string>& pricingFlags,
		ARM_result&	result,
        long objId);

//===========================================================//
// FXVanillaCalculator										 //
//===========================================================//
extern long ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(
		long FxSpreadId,
        string basketTypeStr,
		string digitTypeStr,
		string vanillaTypeStr,
        ARM_result&	result, 
        long objId = -1);


#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
