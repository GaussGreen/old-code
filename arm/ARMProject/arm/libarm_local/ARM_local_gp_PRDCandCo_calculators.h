/*
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculator.h,v $
 * Revision 1.1  2004/03/02 15:08:43  ebenhamou
 * Initial version
 *
 */

#ifndef ARMLOCAL_GP_PRDCANDCO_CALCULATORS_H
#define ARMLOCAL_GP_PRDCANDCO_CALCULATORS_H

#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_gp_genericaddin.h"
#include "ARM_result.h"
#include <GP_Infra\gpinfra\gramfunctorarg.h>
#include <gpbase/gpvector.h>
#include <gpbase/gpmatrix.h>

using ARM::ARM_GramFctorArg;

class ARM_ReferenceValue;


/****************************************************************************
						PRDC CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_PRDCCalculator_Create(
        const long& prdcId,
        const long& modelId,
        const vector< long >& otherMktDataIds,
        const vector< double >& schedulerDatas,
        const vector< double >& truncatorDatas,
        const vector< string >& columnsToPrice,
        const vector< bool >& prdcFlags,
        const string& calibType,
        const vector< double >& calibDatas,
		const string& marginConvertType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_PRDC_Get(
        const long& prdcId,
        const string& getType,
        ARM_result&	result, 
        long        objId );

extern long ARMLOCAL_PRDC_Set(
        long prdcId,
        long dataId,
        const string& setType,
		const vector< string >& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId );

class ARM_PRDCCalculator_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_PRDCCalculator_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_PRDCCalculator_InitFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_PRDCCalculator_InitFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

/****************************************************************************
						PRCS CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_PRCSCalculator_Create (
        const double&			startDate,															
		const double&			fixEndDate,
        const double&			endDate,
		const string&			domCcy,
		const string&			forCcy,
		const string&			fundCcy,
        const string&			cpnFreqStr,
        const string&			cpnDaycountStr,
		const long&             fxResetGap,  
		const string&			stubRuleStr,
		const string&           resetTimingStr,
        const string&			cpnResetCalStr,
        const string&			cpnPayCalStr,
		const long&				cpnnotionalId,
		const long&				domesticCpnId,
		const long&				foreignCpnId,
		const long&				initialFxId,
        const long&				minCpnId,
        const long&				maxCpnId,
        const string&			fundFreqStr,
        const string&			fundDaycountStr,
		const long&				fundingnotionalId,
		const long&				fundMarginId,
		const long&				fundLvgeId,
        const string&			exerFreqStr,
        const long&				NotifGap,
        const string&			payRecStr,
		const long&				nbNCall,
        const long&				feesId,
		const string&			redempTypeStr,
		const long&				redemptionGap,
		const double&			redemptionStrike, 
		const vector <string >& calibTypesStr,
        const vector< double >&	calibDatasDble,
        const vector< string >&	productsStr,
        const vector< double >&	schedulerDatasDble,
		const vector< double >&	truncatorDatasDble,
        const long &			mktDataManagerId,
        ARM_result&				result, 
        long					objId );

/****************************************************************************
						PRKO CALCULATOR
*****************************************************************************/

extern long ARMLOCAL_PRDKOCalculator_Create (
        const double&			startDate,															
		const double&			fixEndDate,
		const double&			switchDate,
        const double&			endDate,
		const string&			domCcy,
		const string&			forCcy,
		const string&			fundCcy,
        const string&			cpnFreqStr,
        const string&			cpnDaycountStr,
		const long&             fxResetGap,  
		const string&			stubRuleStr,
		const string&           resetTimingStr,
        const string&			cpnResetCalStr,
        const string&			cpnPayCalStr,
		const long&				cpnnotionalId,
		const long&				domesticCpnId,
		const long&				foreignCpnId,
		const long&				initialFxId,
        const long&				minCpnId,
        const long&				maxCpnId,
		const long&				barrierId,
        const string&			fundFreqStr,
        const string&			fundDaycountStr,
		const long&				fundingnotionalId,
		const long&				fundMarginId,
		const long&				fundLvgeId,
        const string&			exerFreqStr,
        const long&				NotifGap,
        const string&			payRecStr,
		const long&				nbNCall,
        const long&				feesId,
		const string&			redempTypeStr,
		const long&				redemptionGap,
		const double&			redemptionStrike, 
		const vector <string >& calibTypesStr,
        const vector< double >&	calibDatasDble,
        const vector< string >&	productsStr,
        const vector< double >&	schedulerDatasDble,
		const vector< double >&	truncatorDatasDble,
        const long &			mktDataManagerId,
        ARM_result&				result, 
        long					objId );


/****************************************************************************
						CCS CALCULATOR
*****************************************************************************/

class ARM_CCSCalculator_CreateFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_CCSCalculator_CreateFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};

class ARM_CCSCalculator_InitFunctor : public ARM_GenericAddinFunctor
{
public:
	ARM_CCSCalculator_InitFunctor() {}

	virtual	long operator()( ARM_result& result, long objId );
};


#endif

//-----------------------------------------------------------------------------
/*---- End of file ----*/
