/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_local_gp_pricer.cpp
 *
 *  \brief file for the pricer part of the generic pricer local addins functions
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2003
 */

#include "firstToBeIncluded.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"

#include <GP_Base\gpbase\singleton.h>
#include <GP_Base\gpbase\autocleaner.h>
#include <GP_Base\gpbase\gpvector.h>

#include <GP_Infra\gpinfra\typedef.h>
#include <GP_Infra\gpinfra\genpricer.h>
#include <GP_Infra\gpinfra\gensecurity.h>
#include <GP_Infra\gpinfra\pricingmodel.h>
#include <GP_Infra\gpinfra\pricingadviser.h>
#include <GP_Infra\gpinfra\modelnrefcall.h>
#include <GP_Infra\gpinfra\gramnode.h>
#include <GP_Infra\gpinfra\pricerinfo.h>
#include <GP_Infra\gpinfra\dealdescription.h>
/// kernel
#include <inst/swaption.h>
#include "gpcalib/vanillaarg.h"
#include "gpcalib/vanillaswaption.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/kerneltogp.h"

/// Objects are in namespace, hence the using directive!
using ARM::ARM_GenSecurity;
using ARM::ARM_GenSecurityPtr;
using ARM::ARM_GenPricer;
using ARM::ARM_PricingModel;
using ARM::ARM_AutoCleaner;
using ARM::ARM_PricerInfo;
using ARM::ARM_GramFctorArg;
using ARM::std::vector<double>;
using ARM::ARM_VanillaArg;
using ARM::ARM_VanillaSwaptionArg;
using ARM::ARM_VanillaSpreadOptionArg;

////////////////////////////////////////////
//// Function to create a generic pricer
////////////////////////////////////////////
extern long ARMLOCAL_GenPricer_Create(
	long genSecurityId,
	long pricingModelId,
	const vector<string>& columnNames,
	const vector<double>& columnPrices,
	const string& refColumn,
	const vector<double>& betas,
	ARM_result& result,
	long objId = ARM_NULL_OBJECT_ID )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenPricer* genpricer = NULL;

	try
	{
		ARM_AutoCleaner<ARM_GenPricer> Hold(genpricer);

		ARM_GenSecurity* genSec = NULL;
		if( !GetObjectFromId( &genSec, genSecurityId, ARM_GENSECURITY ) )
		{
			result.setMsg ("ARM_ERR: generic Security is not of a good type");
			return ARM_KO;
		};

		ARM_PricingModel* pricingModel = NULL;
		if( !GetObjectFromId( &pricingModel, pricingModelId, ARM_PRICINGMODEL ) )
		{
			result.setMsg ("ARM_ERR: pricing model is not of a good type");
			return ARM_KO;
		};
		/// creates locally a genpricer
		/// and uses it to price
		genpricer = new ARM_GenPricer( genSec, pricingModel, columnNames, std::vector<double>(columnPrices), refColumn, std::vector<double>(betas) );

		/// price!
		genpricer->Price();
		
		/// assign object
		if( !assignObject( genpricer, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete genpricer;
		genpricer = NULL;
		x.DebugPrint();
		ARM_RESULT();
	}


	/// catch the rest
	catch (...)
	{
		delete genpricer;
		genpricer = NULL;
		result.setMsg ("ARM_ERR: unrecognized failure in creating a generic pricer" );
		return ARM_KO;
	}
}



////////////////////////////////////////////
//// Function to create a generic pricer
////////////////////////////////////////////
long ARMLOCAL_GenPricer_GetData(
	long genPricerId,
	const string& key,
	const string& columnName,
	ARM_GramFctorArg& argResult,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_GenPricer* genPricer = NULL;
		if( !GetObjectFromId( &genPricer, genPricerId, ARM_GENPRICER ) )
		{
			result.setMsg ("ARM_ERR: generic pricer is not of a good type");
			return ARM_KO;
		};
		
		argResult = genPricer->GetPricerInfo()->GetContents(columnName).GetData(key);

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to set detail flag on/off
////////////////////////////////////////////
extern long ARMLOCAL_GenPricer_SetDetailMode(
	long genPricerId,
	bool detailMode,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{   
		ARM_GenPricer* genPricer = NULL;
		if( !GetObjectFromId( &genPricer, genPricerId, ARM_GENPRICER ) )
		{
			result.setMsg ("ARM_ERR: generic pricer is not of a good type");
			return ARM_KO;
		};

		genPricer->SetDetailMode( detailMode );
		string txt( "Detail Mode:" );
		txt += detailMode ? "On" : "Off";
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a deal description
////////////////////////////////////////////
extern long ARMLOCAL_ChangeSecurityIntoGenSec(
	const long& securityId,
	const double& asOfDate,
	const string& modelName,
	ARM_result&	result, 
	long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GenSecurity* newGenSec = NULL;
	try
	{

		ARM_Security* swaption = NULL;
		ARM_VanillaArg*  arg;
		if( !GetObjectFromId( &swaption, securityId,  ARM_SECURITY ) )
		{
			result.setMsg ("ARM_ERR: Security not of a good type");
			return ARM_KO;
		};

		arg = ARM::ARM_ConverterFromKernel::ConvertSecuritytoArgObject( swaption, XLDateToJulian(asOfDate), modelName);
		ARM_VanillaSwaptionArg* swaptionArg = dynamic_cast<ARM_VanillaSwaptionArg*> (arg);
		ARM_VanillaSpreadOptionArg* spreadArg = dynamic_cast<ARM_VanillaSpreadOptionArg*> (arg);
		ARM_GenSecurityPtr GenSecPtr; 
		if(swaptionArg) 
			GenSecPtr = swaptionArg->VanillaSwaptionToGenSec();
		else if(spreadArg) 
			GenSecPtr = spreadArg->VanillaSpreadOptionToGenSec();
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, " Option is not a swaption nor a spreadoption");

		newGenSec = (ARM_GenSecurity*) GenSecPtr->Clone(); 

		/// assign object
		if( !assignObject( newGenSec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}
	
	catch(Exception& x)
	{
		delete newGenSec;
		x.DebugPrint();
		ARM_RESULT();
	}
}


