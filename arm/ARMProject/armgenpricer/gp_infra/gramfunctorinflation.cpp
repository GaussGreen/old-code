/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorir.cpp
 *
 *  \brief gramfunctorir is the file for interest rate related grammar functor
 *	the functor design allows to have a context and a unified command like
 *	interface
 *	\author  N. Belgrade, A. Schauly
 *	\version 1.0
 *	\date March 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/gramfunctorir.h"
#include "gpinfra/gramfunctorinflation.h"

/// gpbase
#include "gpbase/checkinputs.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestrip.h"
#include "gpbase/gpvector.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctorarghelper.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingfunctionir.h"

/// kernel
#include <glob/expt.h>
//#include <inst/irindex.h>

#include <algorithm>


/// gpinflation
//#include "gpinflation/infidx.h"

/// for easy debugging of shared nodes!
#if defined(__GP_SHOW_SHARED_NODE_COORDINATES)
	#include "gpbase/pair.h"
#endif


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_GramFctor
///	Member : itsFuncName is more for debugging output
////////////////////////////////////////////////////

/// the corresponding static string
string ARM_GP_CPI::itsFuncName    = "CPI";
string ARM_GP_InfSwapCommon::itsFuncName = " InfSwapCommon";
string ARM_GP_InfSwapRate::itsFuncName    = "InfSwapRate";
string ARM_GP_InfSwap::itsFuncName    = "InfSwap";

////////////////////////////////////////////////////
/////////////////////  CPI  ////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CPI
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_CPI::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CPI
///	Routine: ConvCPInterpMethod
///	Returns: long
///	Action : converts string to long
////////////////////////////////////////////////////
long ARM_GP_CPI::ConvCPInterpMethod( const string& InfMethod )
{
	string tmp = InfMethod;

	if( tmp == "STEPWISESTART" || tmp == "START" )
			return K_CPISTEPWISESTART;
	if( tmp == "STEPWISEMID" || tmp == "MID" )
		return K_CPISTEPWISEMIDDLE;
	if( tmp == "STEPWISEEND" || tmp == "END")
		return K_CPISTEPWISEEND;
	if( tmp == "STEPWISE")
		return K_CPISTEPWISE;
	if( tmp == "CPILINEAR")
		return K_CPILINEAR;
	if( tmp == "ZCLINEAR")
		return K_ZCLINEAR;
	if( tmp == "ZCCTFWD")
		return K_ZCCTFWD;
	if( tmp == "DEFAULT_PER_CURRENCY")
		return -1;

	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : You specified an unknown Interpolation Method");
	return -1;
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_CPI
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_CPI::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 5, ARM_GP_CPI::itsFuncName  );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		    ARM_GP_CPI::itsFuncName  );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,			    ARM_GP_CPI::itsFuncName  );
		GPAF_CheckArgType( arg[2], GFAT_MATURITY_TYPE,		    ARM_GP_CPI::itsFuncName  );
		GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,		    ARM_GP_CPI::itsFuncName  );
		GPAF_CheckArgType( arg[4], GFAT_MATURITY_TYPE,		    ARM_GP_CPI::itsFuncName  );


		/// grab inputs
		itsInfCurveName			= arg[0].GetString();
		ARM_Date CPIDate		= arg[1].GetDate();
		itsDCFLag				= arg[2].GetString();
		itsDailyInterp			= ConvCPInterpMethod( arg[3].GetString() );
		itsResetLag				= arg[4].GetString();

		/// Checks that the requested CPI Index has been released at evalDate. 
		ARM_Date NewCPIDate = CPIDate;
		NewCPIDate.AddPeriod( itsDCFLag );

		/// compute inputs
		itsModelInflation				= dynamic_cast< ARM_PricingFuncInflation* >(mod);

        if(!itsModelInflation)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only Inflation model can be used with CPI keyword");
		
		itsCPITime			    = mod->GetTimeFromDate( CPIDate );
		itsEvalTime				= evalDate-mod->GetAsOfDate().GetJulian();

		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_CPI
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the CPI function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_CPI::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	ARM_VectorPtr ReturnedCPI =  itsModelInflation->CPISpot( itsInfCurveName, itsEvalTime, itsCPITime, itsDCFLag, itsDailyInterp, itsResetLag, states );
	/// return result
	return ReturnedCPI;
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CPI
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_CPI::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	CPITimeInfo* cpiTimeInfo = new CPITimeInfo( itsInfCurveName, itsDCFLag, itsResetLag, itsDailyInterp, itsCPITime, 0, 0,0,0,0 );
	return ARM_NodeInfo( ARM_GP_VectorPtr(new std::vector<double>(1,evalDate-mod->GetAsOfDate().GetJulian() )),ARM_AdditionalTimeInfoPtr(cpiTimeInfo) );
}



////////////////////////////////////////////////////
/////////////////InfSwapCommon//////////////////////
////////////////////////////////////////////////////
void ARM_GP_InfSwapCommon::ComputeDateStrips( ARM_PricingModel* mod )
{

		ARM_Currency* ccy = mod->GetCurrency( itsZcCurveName );

		itsModelInflation				= dynamic_cast< ARM_PricingFuncInflation* >(mod);
        if(!itsModelInflation)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only Inflation model can be used with InfSwapRate/InfSwap keyword");

		/// To get the pay calendar
		ARM_InfIdx* tmpInfIdx = new ARM_InfIdx( itsInfCurveName ); 
		const char* ourPayCalendar = tmpInfIdx->GetCurrencyUnit()->GetCcyName();
		ARM_INDEX_TYPE defaultIndex  = GetDefaultIndexFromCurrency( ccy->GetCcyName());
		char* resetCalendar = ccy->GetResetCalName(defaultIndex); 

		/// Creates the datestrips that will be used to price the YoY swap
		ARM_DateStrip * denomDateStrip = new ARM_DateStrip( itsStartDate, itsEndDate, itsFloatFreq, itsFloatDayCount, "INF", K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, itsFloatResetGap,itsFloatFreq, itsFloatPaymGap,ourPayCalendar, K_ADVANCE,K_ARREARS);
		ARM_DateStrip * numDateStrip = new ARM_DateStrip( itsStartDate, itsEndDate, itsFloatFreq, itsFloatDayCount, "INF", K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, itsFloatResetGap,itsFloatFreq, itsFloatPaymGap,ourPayCalendar, K_ARREARS,K_ARREARS);

		if( denomDateStrip->GetResetDates()->Elt(0) - itsFloatResetGap < itsEvalTime + mod->GetAsOfDate().GetJulian() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : You are trying to price a CPIFwd in the past!");


		/// We change numDateStrip's FlowStartDates so that they contain the CPIs FixingDates
		std::vector<double> * numStartDates = numDateStrip->GetFlowStartDates();
		std::vector<double> * numResetDates = numDateStrip->GetResetDates();
		size_t i,vecSize = numStartDates->size();
		for(i=0;i<vecSize;i++)
			(*numStartDates)[i] = (*numResetDates)[i] - itsFloatResetGap;

		/// Stores the corresponding ptrs
		itsInfNumDateStrip = ARM_DateStripPtr( numDateStrip );
		itsInfDenomDateStrip = ARM_DateStripPtr( denomDateStrip );

		/// We now build the fixed leg
        char fixCalendar[100];
        ccy->CalcFixPayCal(fixCalendar);

		ARM_DateStrip * FixedLegDateStrip = new ARM_DateStrip( itsStartDate, itsEndDate, itsFixedFreq, itsFixedDayCount, fixCalendar, K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, itsFixedFreq, GETDEFAULTVALUE,fixCalendar );
		itsFixedLegDateStrip = ARM_DateStripPtr( FixedLegDateStrip );

}

////////////////////////////////////////////////////
/////////////////InfSwapRate ///////////////////////
////////////////////////////////////////////////////


void ARM_GP_InfSwapRate::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 13, ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[0],	GFAT_STRING_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[1],	GFAT_STRING_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[2],	GFAT_STRING_TYPE,		    ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[3],	GFAT_DATE_TYPE,			    ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[4],	GFAT_DATE_TYPE,			    ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[5],	GFAT_DOUBLE_TYPE,		    ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[6],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[7],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[8],	GFAT_STRING_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[9],	GFAT_MULTITOKENSTRING_TYPE,	ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[10], GFAT_STRING_TYPE,		    ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[11], GFAT_MULTITOKENSTRING_TYPE, ARM_GP_InfSwapRate::itsFuncName  );
		GPAF_CheckArgType( arg[12], GFAT_DOUBLE_TYPE,			ARM_GP_InfSwapRate::itsFuncName  );

		/// the reset Time is used to check that the eval time is before resetTime

		/// grabs inputs
		itsZcCurveName					= arg[0].GetString();
		itsInfCurveName					= arg[1].GetString();
		itsSwapType						= arg[2].GetString();
		itsStartDate					= arg[3].GetDate();
		itsEndDate						= arg[4].GetDate();
		itsSpread						= arg[5].GetDouble();
		itsFloatResetGap				= (int)arg[6].GetDouble();
		itsFloatPaymGap					= (int)arg[7].GetDouble();
		
		ARM_Currency* ccy = mod->GetCurrency( itsZcCurveName );
		/// We won't have to delete ccy since it is not cloned

		itsFloatFreq					= GetFrequency(arg[8],	mod, itsZcCurveName, ARM_GP_InfSwapRate::itsFuncName, GetFloatFrequencyFtor(ccy));
		itsFloatDayCount				= GetDayCount( arg[9],	mod, itsZcCurveName, GetFloatDayCountFtor(ccy));
		itsFixedFreq					= GetFrequency(arg[10], mod, itsZcCurveName, ARM_GP_InfSwapRate::itsFuncName, GetFloatFrequencyFtor(ccy));
		itsFixedDayCount				= GetDayCount( arg[11], mod, itsZcCurveName, GetFloatDayCountFtor(ccy));
		itsCoupon						= arg[12].GetDouble();
		itsEvalTime						= evalDate-mod->GetAsOfDate().GetJulian();

		ComputeDateStrips( mod );
		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_InfSwapRate
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_InfSwapRate::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new std::vector<double>( 1, evalDate-mod->GetAsOfDate().GetJulian() )),ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_InfSwapRate
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the InfSwapRate function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_InfSwapRate::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	ARM_VectorPtr ReturnedInfSwapRate;

	if( stringGetUpper( itsSwapType ) == "OAT")
		ReturnedInfSwapRate =  itsModelInflation->OATSwapRate( itsZcCurveName, itsInfCurveName, itsEvalTime, itsInfNumDateStrip, itsInfDenomDateStrip, itsFixedLegDateStrip, itsSpread, states );
	else if ( stringGetUpper( itsSwapType ) == "YOY")
		ReturnedInfSwapRate =  itsModelInflation->YoYSwapRate( itsZcCurveName, itsInfCurveName, itsEvalTime, itsInfNumDateStrip, itsInfDenomDateStrip, itsFixedLegDateStrip, itsSpread, states );
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Unknown swap type :"+itsSwapType );		

	return ARM_GramFctorArg( ReturnedInfSwapRate );
}


////////////////////////////////////////////////////
///////////////////InfSwap /////////////////////////
////////////////////////////////////////////////////


void ARM_GP_InfSwap::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 14, ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[0],	GFAT_STRING_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[1],	GFAT_STRING_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[2],	GFAT_STRING_TYPE,		    ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[3],	GFAT_DATE_TYPE,			    ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[4],	GFAT_DATE_TYPE,			    ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[5],	GFAT_DOUBLE_TYPE,		    ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[6],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[7],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[8],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[9],	GFAT_STRING_TYPE,			ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[10], GFAT_MULTITOKENSTRING_TYPE,	ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[11], GFAT_STRING_TYPE,		    ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[12], GFAT_MULTITOKENSTRING_TYPE, ARM_GP_InfSwap::itsFuncName  );
		GPAF_CheckArgType( arg[13],	GFAT_DOUBLE_TYPE,			ARM_GP_InfSwap::itsFuncName  );

		/// the reset Time is used to check that the eval time is before resetTime

		/// grabs inputs
		itsZcCurveName					= arg[0].GetString();
		itsInfCurveName					= arg[1].GetString();
		itsSwapType						= arg[2].GetString();
		itsStartDate					= arg[3].GetDate();
		itsEndDate						= arg[4].GetDate();
		itsStrike						= arg[5].GetDouble();
		itsSpread						= arg[6].GetDouble();
		itsFloatResetGap				= (int)arg[7].GetDouble();
		itsFloatPaymGap					= (int)arg[8].GetDouble();

		ARM_Currency* ccy = mod->GetCurrency( itsZcCurveName );
		/// We won't have to delete ccy since it is not cloned
		ARM_INDEX_TYPE defaultIndex  = GetDefaultIndexFromCurrency( ccy->GetCcyName());

		itsFloatFreq					= GetFrequency(arg[9],	mod, itsZcCurveName, ARM_GP_InfSwap::itsFuncName, GetFloatFrequencyFtor(ccy));
		itsFloatDayCount				= GetDayCount( arg[10], mod, itsZcCurveName, GetFloatDayCountFtor(ccy));
		itsFixedFreq					= GetFrequency(arg[11], mod, itsZcCurveName, ARM_GP_InfSwap::itsFuncName, GetFloatFrequencyFtor(ccy));
		itsFixedDayCount				= GetDayCount( arg[12], mod, itsZcCurveName, GetFloatDayCountFtor(ccy));
		itsCoupon						= arg[13].GetDouble();
		itsEvalTime						= evalDate-mod->GetAsOfDate().GetJulian();

		ComputeDateStrips( mod );
		SetAlreadyComputed(true);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_InfSwapRate
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : nobody knows
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_InfSwap::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new std::vector<double>( 1, evalDate-mod->GetAsOfDate().GetJulian() )),ARM_AdditionalTimeInfoPtr(NULL) );
}


ARM_GramFctorArg ARM_GP_InfSwap::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	ARM_VectorPtr ReturnedInfSwapRate;

	if( stringGetUpper( itsSwapType ) == "OAT")
		ReturnedInfSwapRate =  itsModelInflation->OATSwap( itsZcCurveName, itsInfCurveName, itsEvalTime, itsStrike, itsSpread, itsInfNumDateStrip, itsInfDenomDateStrip, itsFixedLegDateStrip, itsCoupon, states );
	else if ( stringGetUpper( itsSwapType ) == "YOY")
		ReturnedInfSwapRate =  itsModelInflation->YoYSwap( itsZcCurveName, itsInfCurveName, itsEvalTime, itsStrike, itsSpread, itsInfNumDateStrip, itsInfDenomDateStrip, itsFixedLegDateStrip, itsSpread, states );
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Unknown swap type :"+itsSwapType );

	return ARM_GramFctorArg( ReturnedInfSwapRate );
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

