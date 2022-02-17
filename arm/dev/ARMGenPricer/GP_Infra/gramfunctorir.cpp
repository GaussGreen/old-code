/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorir.cpp
 *
 *  \brief gramfunctorir is the file for interest rate related grammar functor
 *	the functor design allows to have a context and a unified command like
 *	interface
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpinfra/gramfunctorir.h"

/// gpbase
#include "gpbase/gplinalgconvert.h"
#include "gpbase/checkinputs.h"
#include "gpbase/env.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestrip.h"
#include "gpbase/gpvector.h"
#include "gpbase/stringconvert.h"
#include "gpbase/argconvdefault.h"
#include "gpbase/stringmanip.h"
#include "gpbase/numericconstant.h"

/// gpinfra
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gramfunctorarghelper.h"
#include "gpinfra/gramnode.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingmodelir.h"
#include "gpinfra/pricingfunctionir.h"
#include "gpinfra/irrate.h"

/// gpmodel
#include "gpmodels/multiassets.h"
#include "gpmodels/forwardforex.h"

/// kernel
#include <glob/expt.h>
#include <inst/irindex.h>
#include <inst/corridorleg.h>
//#include <inst/corridordblcondition.h>
//#include <inst/spreadoption.h>
#include <util/refvalue.h>
#include <util/fromto.h>

#include <algorithm>
#include <string>

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
string ARM_GP_DF::itsFuncName			= "DF (discount factor)";
string ARM_GP_Libor::itsFuncName		= "Libor";
string ARM_GP_Annuity::itsFuncName		= "Annuity";
string ARM_GP_SwapRate::itsFuncName		= "SwapRate";
string ARM_GP_Swap::itsFuncName			= "Swap";
string ARM_GP_BasisSwap::itsFuncName	= "BasisSwap";
string ARM_GP_Swaption::itsFuncName		= "Swaption";
string ARM_GP_Cap::itsFuncName			= "Cap";
string ARM_GP_SumOption::itsFuncName    = "SumOption";
string ARM_GP_SpreadOption::itsFuncName	= "SpreadOption";
string ARM_GP_Corridor::itsFuncName		= "Corridor";
string ARM_GP_Spread::itsFuncName		= "SpreadCMS";
string ARM_GP_MaxRate::itsFuncName		= "MaxRate";
string ARM_GP_ImpliedVol::itsFuncName	= "ImpliedVol";
string ARM_GP_DoubleDigital::itsFuncName= "DoubleDigital";

////////////////////////////////////////////////////
///												////
///			FINANCIAL FUNCTORS					////
///												////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
////////// Discount Factor functor /////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_DF
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_DF::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 2, ARM_GP_DF::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_DF::itsFuncName );
	itsCurveName = arg[0].GetString();

	/// Modification of nodes if necessary 
	/// set the end date taking the eval date as the start date
	/// a little special compared to other set defaults methods!
	ARM_Currency* ccy			= mod->GetCurrency( itsCurveName );
	char* ccyName				= ccy->GetCcyName();
	ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
	char* payCalendar		    = ccy->GetPayCalName(defaultIndex);

	if( arg[1].GetType() == GFAT_STRING_TYPE )
        nodes[1] = ARM_ExpNodePtr( new ARM_ExpNodeDate( ComputeEndDate(evalDate, arg[1], payCalendar ) ) );

    delete payCalendar;
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DF
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_DF::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, evalDate, nodes );

		/// checking of the size and type
		GPAF_CheckArgSize( arg, 2, ARM_GP_DF::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_DF::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATEORMATU_TYPE,ARM_GP_DF::itsFuncName );

		itsCurveName			= arg[0].GetString();
		itsModel				= mod;
		const char* calendar	= mod->GetCurrency( itsCurveName )->GetCcyName();
		
		/// the second argument can either be a date or a maturity
		ARM_Date endDate= ComputeEndDate( evalDate, arg[1], calendar );
		itsMaturityTime = mod->GetTimeFromDate( endDate );
		itsEvalTime		= GetAndCheckEvalTime( evalDate, mod, ARM_GP_DF::itsFuncName );

		CheckNbSmaller( itsEvalTime, itsMaturityTime, "EvalTime", "MaturityTime", ARM_GP_DF::itsFuncName,__LINE__,__FILE__  );
		
		/// do not do it again!
		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DF
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the DF function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_DF::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	ARM_VectorPtr df =  itsModel->DiscountFactor( itsCurveName, itsEvalTime , itsMaturityTime, states );
	return ARM_GramFctorArg(df);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DF
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_DF::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 1, itsMaturityTime )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////  Forward Rate (LIBOR) //////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_Libor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Libor::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 7, ARM_GP_Libor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Libor::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDateResetCal(		nodes, arg, 2, mod, itsCurveName );
	bool isDefault = SetDayCount(    nodes, arg, 3, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	SetResetDaysGap(nodes, arg, 4, mod, itsCurveName );
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Libor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Libor::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// checking of the size and type
		GPAF_CheckArgSize( arg, 7, ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		    ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,			    ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,	    ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,		    ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[4], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[5], GFAT_DATE_OR_DOUBLE_TYPE,	ARM_GP_Libor::itsFuncName  );
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,		    ARM_GP_Libor::itsFuncName  );

		/// the reset Time is used to check that the eval time is before resetTime
		

		/// grab inputs
		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();

		/// compute inputs
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for LIBOR keyword");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex 
								= GetDefaultIndexFromCurrency( ccyName );

		char* payCalendar = ccy->GetPayCalName(defaultIndex);
		char* resetCalendar = ccy->GetResetCalName(defaultIndex); 

		/// The second argument can either be a date or a maturity
        /// then in case of a maturity, first of all add without
        /// startDate adjustment to get endDate
        /// Becareful : it must be consistent with ARM_GP_Libor::SetDefaults()
        /// ***** Is it the RIGHT default convention ? *****
		ARM_Date endDate	= ComputeEndDate( startDate, arg[2], resetCalendar );

		int dayCount	    = GetDayCount( arg[3], mod, itsCurveName, GetFloatDayCountFtor(ccy));

        /// Compute reset date according to timing & reset gap
        ARM_Date resetDate( ComputeResetOrPayDate(startDate,endDate,K_ADVANCE,resetCalendar,arg[4],
                                                  mod,itsCurveName,ARM_GP_Libor::itsFuncName) );

		/// Compute payment date according to timing & payment gap
		string advArrString = arg[6].GetString();
		int payTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
        ARM_Date payDate( ComputeResetOrPayDate(startDate,endDate,payTiming,payCalendar,arg[5],
                                                mod,itsCurveName,ARM_GP_Libor::itsFuncName) );


		itsPeriod			    = CountYearsWithoutException( dayCount, startDate, endDate );
		itsEvalTime			    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_Libor::itsFuncName );
		itsFwdStartTime		    = mod->GetTimeFromDate( startDate );
		itsFwdEndTime			= mod->GetTimeFromDate( endDate );
		itsResetTime        	= mod->GetTimeFromDate( resetDate );
		itsPayTime		        = mod->GetTimeFromDate( payDate );

		/// some validation
		CheckNbSmaller( itsEvalTime, itsResetTime,	    "EvalTime",     "ResetTime",  ARM_GP_Libor::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTime,itsFwdEndTime,	"FwdStartTime", "FwdEndTime",	  ARM_GP_Libor::itsFuncName,__LINE__,__FILE__ );
		CheckPositiveNb( itsPeriod, "Period", ARM_GP_Libor::itsFuncName,__LINE__,__FILE__ );

		/// remove char for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Libor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the LIBOR function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Libor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	ARM_VectorPtr Libor =  itsModelIR->Libor( itsCurveName, itsEvalTime, itsFwdStartTime, itsFwdEndTime, itsPeriod, itsResetTime, itsPayTime, states );
	/// return result
	return ARM_GramFctorArg(Libor);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Libor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Libor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
    
	vector<double> tmpResult;
	tmpResult.reserve(2);
	
	tmpResult.push_back( itsResetTime );
	if( itsResetTime != itsFwdStartTime )
		tmpResult.push_back(  itsFwdStartTime );
	
	tmpResult.push_back( itsFwdEndTime );
    if( itsPayTime != itsFwdEndTime )
		tmpResult.push_back( itsPayTime );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( tmpResult.size(), tmpResult.begin() )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////		 Annuity	      //////////////////
////////////////////////////////////////////////////



////////////////////////////////////////////////////
///	Class  : ARM_GP_Annuity
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Annuity::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 8, ARM_GP_Annuity::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Annuity::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(		nodes, arg, 2, mod, itsCurveName );
    bool isDefault;
	isDefault=SetFrequency(   nodes, arg, 3, mod, itsCurveName, GetFixedFrequencyFtor(ccy));
	isDefault=SetDayCount(	nodes, arg, 4, mod, itsCurveName, GetFixedDayCountFtor(ccy));
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Annuity
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Annuity::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes  )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// checking of the size and type
		GPAF_CheckArgSize( arg, 8, ARM_GP_Annuity::itsFuncName);
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Annuity::itsFuncName);
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_Annuity::itsFuncName);
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_Annuity::itsFuncName);
		GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,	ARM_GP_Annuity::itsFuncName);
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,	ARM_GP_Annuity::itsFuncName);

		/// grab inputs
		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();

		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingModelIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for ANNUITY keyword");

        size_t dateStripOffset = static_cast<size_t>(arg[7].GetDouble());

        ARM_DateStripPtr sched(NULL);
        if(arg[6].GetType() == GFAT_DATESTRIP_TYPE)
            sched = arg[6].GetDateStrip();

		itsEvalTime			= GetAndCheckEvalTime( evalDate, mod, ARM_GP_Annuity::itsFuncName );

        if(sched == ARM_DateStripPtr(NULL))
        {
		    ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		    char* ccyName			= ccy->GetCcyName();

            char fixCalendar[100];
            ccy->CalcFixPayCal(fixCalendar);

		    /// the second argument can either be a date or a maturity
		    ARM_Date endDate	    = ComputeEndDate( startDate, arg[2], fixCalendar );


            int frequency   = GetFrequency(arg[3], mod, itsCurveName, ARM_GP_Annuity::itsFuncName, GetFixedFrequencyFtor(ccy));
		    int fixDayCount = GetDayCount( arg[4], mod, itsCurveName, GetFixedDayCountFtor(ccy));

		    ARM_DateStrip dateStrip( startDate, endDate, frequency, fixDayCount, fixCalendar,
			    K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, frequency, GETDEFAULTVALUE,
			    fixCalendar );

		    
		    /// not cloned hence no need to delete this!
		    ARM_GP_Vector* pfixPayTimes   = dateStrip.GetFlowEndDates();
		    ARM_GP_Vector* pfixPayPeriods = dateStrip.GetInterestTerms();

		    /// copy constructor
		    itsFixPayTimes	= ARM_GP_Vector( *pfixPayTimes );
		    itsFixPayPeriods= ARM_GP_Vector( *pfixPayPeriods );

		    /// get time from date! 
		    for(int i=0; i<pfixPayTimes->size(); ++i ){
			    itsFixPayTimes[i] = mod->GetTimeFromDate( itsFixPayTimes[i] );
		    }

		    /// Case of Notional Vector
		    size_t sizefixflows = itsFixPayTimes.size();
		    if( arg[5].GetType() ==  GFAT_DOUBLE_TYPE )
		    {
			    double notional = arg[5].GetDouble();
				itsIsVariableNotional =  (notional != 1) ? true :  false;
			    itsVariableNotional = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,notional ));
		    }
		    else if( arg[5].GetType() ==  GFAT_VECTOR_TYPE )
		    {
			   // itsVariableNotional	= arg[5].GetVector();
				itsVariableNotional = GetProfileValues(arg[5],dateStripOffset,itsFixPayTimes,"Notional");
			    itsIsVariableNotional = true;
		    }
		    else if( arg[5].GetType() ==  GFAT_CURVE_TYPE )
		    {
			    ARM_GP_CurvePtr notioCurve	= arg[5].GetCurve();
			    itsVariableNotional = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,1.0)); 
			    /// WARNING We interpolate with respect to FixPayTimes
			    for(i=0; i<sizefixflows; ++i)			
					    (*itsVariableNotional)[i] = notioCurve->Interpolate(itsFixPayTimes[i]);		
			    itsIsVariableNotional = true;
		    }
		    else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "Notional should be double, vector or ARM_Curve" );
        }
        else
        {
            /// Date strip is input
            size_t i,nbFixFlows = sched->GetPaymentDates()->size();
            if(dateStripOffset >= nbFixFlows)
                dateStripOffset = nbFixFlows-1;
            int nbActFixFlows=nbFixFlows-dateStripOffset;

            /// Fill arguments from the dateStrip
		    itsFixPayTimes      = ARM_GP_Vector(nbActFixFlows);
		    itsFixPayPeriods    = ARM_GP_Vector(nbActFixFlows);
		    for(i=0;i+dateStripOffset<nbFixFlows;++i)
            {
			    itsFixPayTimes[i]    = mod->GetTimeFromDate((*(sched->GetPaymentDates()))[i+dateStripOffset]);
			    itsFixPayPeriods[i]  = (*(sched->GetInterestTerms()))[i+dateStripOffset];
            }

            itsVariableNotional = GetProfileValues(arg[5],dateStripOffset,itsFixPayTimes,"Notional");
            itsIsVariableNotional = true;
        }
		
		/// Validation
		CheckNbSmaller( itsEvalTime, itsFixPayTimes[0],	"EvalTime",	"1st Fix PayTimes",	ARM_GP_Annuity::itsFuncName,__LINE__,__FILE__ );
		CheckVectorIncreasing( itsFixPayTimes,	"FixPayTimes",	ARM_GP_Annuity::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFixPayPeriods,"FixPayPeriods",ARM_GP_Annuity::itsFuncName ,__LINE__,__FILE__ );


		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Annuity
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Annuity function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Annuity::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	ARM_VectorPtr Annuity; 
	if (itsIsVariableNotional)
		Annuity =  itsModelIR->AnnuityWithNominal( itsCurveName, itsEvalTime, itsFixPayTimes, itsFixPayPeriods, *itsVariableNotional, states );
	else
		Annuity =  itsModelIR->Annuity( itsCurveName, itsEvalTime, itsFixPayTimes, itsFixPayPeriods, states );

	/// return result
	return ARM_GramFctorArg(Annuity);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Annuity
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Annuity::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( itsFixPayTimes )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////		Swap Rate		  //////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_SwapRate
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_SwapRate::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 16, ARM_GP_SwapRate::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_SwapRate::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(		nodes, arg, 2, mod, itsCurveName );
    bool isDefault,isDbleNotional;
	isDefault=SetFrequency(   nodes, arg, 3, mod, itsCurveName, GetFixedFrequencyFtor(ccy));
	isDefault=SetDayCount(	nodes, arg, 4, mod, itsCurveName, GetFixedDayCountFtor(ccy));

    /// If default values are used set the double notional pricing method flag
	isDbleNotional=SetFrequency(   nodes, arg, 5, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	isDefault=SetDayCount(	nodes, arg, 6, mod, itsCurveName, GetFloatDayCountFtor(ccy));
    isDbleNotional &= isDefault;

    /// Set of the internal argument for double notional approximation
    if(isDbleNotional)
    {
        nodes[9] = ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble(1) );
        arg[9].SetDouble(1);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SwapRate
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_SwapRate::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 16, ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[8], GFAT_STRING_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[9], GFAT_DOUBLE_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[10],GFAT_DATESTRIP_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_DATESTRIP_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[12],GFAT_DOUBLE_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[13],GFAT_DOUBLE_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[14],GFAT_DOUBLE_TYPE,	ARM_GP_SwapRate::itsFuncName );
		GPAF_CheckArgType( arg[15],GFAT_DOUBLE_TYPE,	ARM_GP_SwapRate::itsFuncName );

		/// grab inputs
		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();
        itsDbleNotional         = (arg[9].GetDouble() == 1);

		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for SWAPRATE keyword");
		
		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_SwapRate::itsFuncName );

		itsFloatDateStripOffset = static_cast<size_t>(arg[12].GetDouble());
		itsFixDateStripOffset = static_cast<size_t>(arg[13].GetDouble());
		int floatDateStripWidth = static_cast<int>(arg[14].GetDouble());
		int fixDateStripWidth = static_cast<int>(arg[15].GetDouble());

		ARM_Currency* ccy		    = mod->GetCurrency( itsCurveName );
		char* ccyName			    = ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex =  GetDefaultIndexFromCurrency( ccyName );
		
		char fixCalendar[100];
		ccy->CalcFixPayCal(fixCalendar);
		
		int stubType = ARM_ArgConv_StubRules.GetNumber(arg[8].GetString());

		/// the second argument can either be a date or a maturity
		ARM_Date endDate = ComputeEndDate( startDate, arg[2], fixCalendar );

		ARM_DateStripPtr dateStrip(NULL);
        if(arg[11].GetType() == GFAT_DATESTRIP_TYPE)
            dateStrip = arg[11].GetDateStrip();

		ARM_DateStrip fixDateStrip;

		if(dateStrip == ARM_DateStripPtr(NULL)){

			/// Fixed leg data extraction
			/// get the fixed leg convention
			int fixFreq	    = GetFrequency(arg[3], mod, itsCurveName, ARM_GP_SwapRate::itsFuncName, GetFixedFrequencyFtor(ccy));
			int fixDayCount = GetDayCount( arg[4], mod, itsCurveName, GetFixedDayCountFtor(ccy));

			/// fix leg datestrip
			fixDateStrip = ARM_DateStrip( startDate, endDate, fixFreq, fixDayCount, fixCalendar, K_MOD_FOLLOWING,
				K_ADJUSTED,	stubType, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE, fixCalendar );

		}
		else{
			fixDateStrip = *dateStrip;

			// By default the datestrip until the end
			if (fixDateStripWidth == -1)
				fixDateStripWidth = fixDateStrip.size() - itsFixDateStripOffset;

			fixDateStrip.ResizeAndBuilt(itsFixDateStripOffset, itsFixDateStripOffset+fixDateStripWidth);
		}
		
		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pfixResetTimes  = fixDateStrip.GetResetDates();
		ARM_GP_Vector* pfixPayTimes    = fixDateStrip.GetPaymentDates();
		ARM_GP_Vector* pfixPayPeriods  = fixDateStrip.GetInterestTerms();
		
		/// copy constructor
		itsFixResetTimes	= ARM_GP_Vector( *pfixResetTimes );  
		itsFixPayTimes	    = ARM_GP_Vector( *pfixPayTimes );
		itsFixPayPeriods   = ARM_GP_Vector( *pfixPayPeriods );

		/// get time from date!
		for( int i(0); i<pfixPayTimes->size(); ++i )
		{
			itsFixResetTimes[i] = mod->GetTimeFromDate( itsFixResetTimes[i] );
			itsFixPayTimes[i] = mod->GetTimeFromDate( itsFixPayTimes[i] );
		}

		int fwdDayCount	= GetFloatDayCountFtor(ccy)();
		
		bool withFixDateStrip = false;

		dateStrip = ARM_DateStripPtr(NULL);
		if(arg[10].GetType() == GFAT_DATESTRIP_TYPE)
		{
			dateStrip = arg[10].GetDateStrip();
			withFixDateStrip = true;
		}

		ARM_DateStrip floatDateStrip;

		if(dateStrip == ARM_DateStripPtr(NULL)){
			
			/// the logic here is to get the default calendar for the float leg
			/// these are from the default index
			/// the fix calendar is the currency itself!
			char* floatPayCalendar	= ccy->GetPayCalName(defaultIndex);
			char* floatResetCalendar= ccy->GetResetCalName(defaultIndex);
			CC_NS(std,auto_ptr)<char> holdfloatPayCalendar(floatPayCalendar);
			CC_NS(std,auto_ptr)<char> holdfloatResetCalendar(floatResetCalendar);


			/// Floating leg data extraction
			/// get the floating leg convention
			int floatFreq	    = GetFrequency(arg[5], mod, itsCurveName, ARM_GP_SwapRate::itsFuncName, GetFloatFrequencyFtor(ccy));
			int floatDayCount	= GetDayCount( arg[6], mod, itsCurveName, GetFloatDayCountFtor(ccy));

			/// floating leg datestrip
			floatDateStrip = ARM_DateStrip( startDate, endDate, floatFreq, floatDayCount, floatResetCalendar,
					K_MOD_FOLLOWING, K_ADJUSTED, stubType, GETDEFAULTVALUE, floatFreq, GETDEFAULTVALUE,
					floatPayCalendar );
		}
		else{
			floatDateStrip = *dateStrip;

			// By default the datestrip until the end
			if (floatDateStripWidth == -1)
				floatDateStripWidth = floatDateStrip.size() - itsFloatDateStripOffset;

			floatDateStrip.ResizeAndBuilt(itsFloatDateStripOffset, itsFloatDateStripOffset + floatDateStripWidth);
		}

		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pFwdResetTimes      = floatDateStrip.GetResetDates();
		ARM_GP_Vector* pFwdStartTimes      = floatDateStrip.GetFwdRateStartDates();
		ARM_GP_Vector* pFwdEndTimes        = floatDateStrip.GetFwdRateEndDates();
		ARM_GP_Vector* pFloatPayTimes      = floatDateStrip.GetPaymentDates();
		ARM_GP_Vector* pFloatPayPeriods    = floatDateStrip.GetInterestTerms();

		/// copy constructor
		itsFwdResetTimes	= ARM_GP_Vector( *pFwdResetTimes );
		itsFwdStartTimes	= ARM_GP_Vector( *pFwdStartTimes );
		itsFwdEndTimes      = ARM_GP_Vector( *pFwdEndTimes );
		itsFwdPayPeriods    = ARM_GP_Vector(pFwdStartTimes->size());
		itsFloatPayTimes	= ARM_GP_Vector( *pFloatPayTimes );
		itsFloatPayPeriods  = ARM_GP_Vector( *pFloatPayPeriods );

		/// get time from date! 
		for(i=0; i<pFwdStartTimes->size(); ++i )
		{
		    itsFwdPayPeriods[i] = CountYearsWithoutException( fwdDayCount, ARM_Date(itsFwdStartTimes[i]), ARM_Date(itsFwdEndTimes[i]) );
			itsFwdResetTimes[i] = mod->GetTimeFromDate( itsFwdResetTimes[i] );
			itsFwdStartTimes[i] = mod->GetTimeFromDate( itsFwdStartTimes[i] );
			itsFwdEndTimes[i]   = mod->GetTimeFromDate( itsFwdEndTimes[i] );
			itsFloatPayTimes[i] = mod->GetTimeFromDate( itsFloatPayTimes[i] );
		}

		ARM_Date floatStartDate( startDate );
		itsFloatStartTime   = mod->GetTimeFromDate( floatStartDate );

		ARM_Date floatEndDate( endDate );
		itsFloatEndTime		= mod->GetTimeFromDate( floatEndDate );

		/// Case of Margin Vector
		size_t sizefloatflows = itsFloatPayTimes.size();        
		if( arg[7].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double margin = arg[7].GetDouble();
			itsMarginVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,margin));
		}
		else if( arg[7].GetType() ==  GFAT_VECTOR_TYPE )
		{
			ARM_GP_VectorPtr marginVector = arg[7].GetVector();
			itsMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,0.0)); 
			for(i=0; i<sizefloatflows; ++i)
				(*itsMarginVector)[i] = (*marginVector)[i+itsFloatDateStripOffset];
		}
		else if( arg[7].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr marginCurve	= arg[7].GetCurve();
			itsMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,0.0)); 
			for(i=0; i<sizefloatflows; ++i)			
					(*itsMarginVector)[i] = marginCurve->Interpolate(itsFwdResetTimes[i]);										
		}

		/// Validation
		CheckNbSmaller( itsEvalTime,	itsFloatStartTime,		"EvalTime",			"FloatStartTime",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime, itsFloatEndTime,		"FloatStartTime",	"FloatEndTime",		ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime, itsFixPayTimes[0] ,	"FloatStartTime",   "1st Fix PayTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__ );

        CheckNbSmaller( itsEvalTime,	itsFwdStartTimes[0],		"EvalTime",			    "1st Fwd StartTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTimes[0], itsFixPayTimes[0] ,	"1st Fwd StartTimes", "1st Fix PayTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__ );

		CheckVectorIncreasing( itsFixPayTimes,	"FixPayTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFixPayPeriods,"FixPayPeriods",ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFwdStartTimes,	"FwdStartTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes,	    "FwdEndTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFloatPayTimes,	"FloatPayTimes",	ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFloatPayPeriods,  "FixPayPeriods",ARM_GP_SwapRate::itsFuncName,__LINE__,__FILE__  );

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SwapRate
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Swap Rate function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_SwapRate::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	ARM_VectorPtr SwapRate =  itsModelIR->SwapRate( itsCurveName, itsEvalTime,
        itsFloatStartTime, itsFloatEndTime,
		itsFixPayTimes, itsFixPayPeriods,
        itsFwdStartTimes, itsFwdEndTimes, itsFwdPayPeriods,
        itsFloatPayTimes, itsFloatPayPeriods, *itsMarginVector, itsDbleNotional, states );

	return ARM_GramFctorArg(SwapRate);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SwapRate
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_SwapRate::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

    CC_STL_VECTOR( double ) result;
    result.push_back( itsFloatStartTime );
    result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
    result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
    result.insert( result.end(), itsFloatPayTimes.begin(), itsFloatPayTimes.end() );
    result.push_back( itsFloatEndTime );
    result.insert( result.end(), itsFixPayTimes.begin(), itsFixPayTimes.end() );

	/// sort and remove duplicates!
    CC_NS( std, sort )( result.begin(), result.end() );
    CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
    result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Swap
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Swap::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 19, ARM_GP_Swap::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_Swap::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(		nodes, arg, 2, mod, itsCurveName );
    bool isDefault,isDbleNotional;
	isDefault=SetFrequency(   nodes, arg, 5, mod, itsCurveName, GetFixedFrequencyFtor(ccy));
	isDefault=SetDayCount(	nodes, arg, 6, mod, itsCurveName, GetFixedDayCountFtor(ccy));

    /// If default values are used set the double notional pricing method flag
	isDbleNotional=SetFrequency(   nodes, arg, 7, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	isDefault=SetDayCount(	nodes, arg, 8, mod, itsCurveName, GetFloatDayCountFtor(ccy));
    isDbleNotional &= isDefault;

    /// Set of the internal argument for double notional approximation
    if(isDbleNotional)
    {
        nodes[12] = ARM_ExpNodePtr( new ARM_ExpNodeDateOrDouble(1) );
        arg[12].SetDouble(1);
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Swap
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Swap::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Swap::itsFuncName );
	

	/// if it is not already computed
	if(!GetAlreadyComputed())
	{			
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 19, ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[7], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[8], GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_STRING_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[12],GFAT_DOUBLE_TYPE,	ARM_GP_Swap::itsFuncName );

		GPAF_CheckArgType( arg[13],GFAT_DATESTRIP_TYPE, ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[14],GFAT_DATESTRIP_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[15],GFAT_DOUBLE_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[16],GFAT_DOUBLE_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[17],GFAT_DOUBLE_TYPE,	ARM_GP_Swap::itsFuncName );
		GPAF_CheckArgType( arg[18],GFAT_DOUBLE_TYPE,	ARM_GP_Swap::itsFuncName );

		int argSwapRateSize = 16;
		ARM_GramFctorArgVector argSwapRate(argSwapRateSize);

		vector<ARM_ExpNodePtr> nodesSwapRate(argSwapRateSize); 
		
		for(int i=0; i<3; ++i)
			argSwapRate[i] = arg[i];
	
		for(i=3; i<8; ++i)
			argSwapRate[i] = arg[i+2];

		argSwapRate[8] = arg[11];

		for(i=9; i<16; ++i)
			argSwapRate[i] = arg[i+3];

		ARM_GP_SwapRate::GrabInputs( argSwapRate, mod, evalDate, nodesSwapRate );

		/// handle vectorial strike
		size_t sizefixflows = itsFixPayTimes.size();    
		size_t sizefloatflows = itsFloatPayTimes.size(); 

		bool strikeIsNull (false);

		if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double strike = arg[3].GetDouble();
			itsStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizefixflows,strike));
			if (fabs(strike)<1e-13)
				strikeIsNull = true;
		}
		else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
		{
			ARM_GP_VectorPtr strikePerState	= arg[3].GetVector();
			itsStrikeMatrix = ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,strikePerState->size(),0.0));
			for(size_t i=0; i<strikePerState->size(); ++i)
				(*itsStrikeMatrix)(0,i) = (*strikePerState)[i+itsFixDateStripOffset];
		}
		else if( arg[3].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr strikeCurve	= arg[3].GetCurve();
			itsStrikeMatrix = ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,sizefixflows,0.0));
			for(i=0; i<sizefixflows; ++i)			
					(*itsStrikeMatrix)(0,i) = strikeCurve->Interpolate(itsFixResetTimes[i]);										
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be double, vector or ARM_Curve" );
		

		//Set default for node
		SetDefaults( arg, mod, nodes );
	
		string payRecString	= arg[4].GetString();
		itsPayRec		    = ARM_ArgConv_PayRec.GetNumber( payRecString );

		/// Case of Notional Vector
		if( arg[10].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double notional = arg[10].GetDouble();
			itsFixNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,notional));
			itsFloatNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,notional));

		}
		else if( arg[10].GetType() ==  GFAT_VECTOR_TYPE )
		{
			///
			/// performance issue
			/// case of strike = 0 ---> set fixed leg notional to 0 and don't check size
			/// in later swap PV, the computation of the fixed leg will be avoided 
			if(strikeIsNull)
			{
				itsFixNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows, 0.0));
				ARM_GP_VectorPtr floatNotionalVector = arg[10].GetVector();
				itsFloatNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,0.0)); 
				for(i=0; i<sizefloatflows; ++i)
					(*itsFloatNotionalVector)[i] = (*floatNotionalVector)[i+itsFloatDateStripOffset];
			}
			else 
			{	/// check that fix and float legs are of same size
				if (sizefixflows != sizefloatflows)
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " SWAP keyword : for notional vector, both legs are required to be of same frequency" );


				ARM_GP_VectorPtr notionalVector = arg[10].GetVector();
				itsFixNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,0.0));
				itsFloatNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,0.0)); 
				for(i=0; i<sizefixflows; ++i)
				{
					(*itsFixNotionalVector)[i] = (*notionalVector)[i+itsFixDateStripOffset];
					(*itsFloatNotionalVector)[i] = (*notionalVector)[i+itsFloatDateStripOffset];
				}	
			}
		}
		else if( arg[10].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr notionalCurve	= arg[10].GetCurve();
			itsFixNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,1.0)); 
			itsFloatNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,1.0)); 
			for(i=0; i<sizefixflows; ++i)			
					(*itsFixNotionalVector)[i] = notionalCurve->Interpolate(itsFixPayTimes[i]);
			for(i=0; i<sizefloatflows; ++i)			
					(*itsFloatNotionalVector)[i] = notionalCurve->Interpolate(itsFloatPayTimes[i]);

								
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "Notional should be double, vector or ARM_Curve" );

		
		
		SetAlreadyComputed(true);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Swap
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Cap Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Swap::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
		
	/// Conversion to matriciel strikes
	size_t nbRows = itsStrikeMatrix->GetRowsNb();
    size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
	
	if( nbRows != statesSize && nbRows !=1 )
	{
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " strikeMatrix rows size is failed" );
    }
	else if(nbRows == 1 &&  nbRows != statesSize)
	{
		ARM_GP_Vector* vectStrike = itsStrikeMatrix->GetRow(0);
		itsStrikeMatrix->reserve(statesSize, vectStrike->size() );
		for(size_t i=1; i<statesSize; ++i)		
			itsStrikeMatrix->push_backRow(*vectStrike);
		delete vectStrike;
	}
	   
		ARM_VectorPtr NPVSwap =  itsModelIR->NPVSwap( 
		itsCurveName,
		itsEvalTime,
        itsFloatStartTime,
		itsFloatEndTime,
		itsFixPayTimes,
		itsFixPayPeriods,
        itsFwdStartTimes,
		itsFwdEndTimes,
		itsFwdPayPeriods,
        itsFloatPayTimes,
		itsFloatPayPeriods,		
		*itsMarginVector,
		itsDbleNotional,
		*itsFixNotionalVector,
		*itsFloatNotionalVector,
	    *itsStrikeMatrix,
		itsPayRec,
		states );
	return ARM_GramFctorArg(NPVSwap);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Swap
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Swap::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

    CC_STL_VECTOR( double ) result;
    result.push_back( itsFloatStartTime );
    result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
    result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
    result.insert( result.end(), itsFloatPayTimes.begin(), itsFloatPayTimes.end() );
    result.push_back( itsFloatEndTime );
    result.insert( result.end(), itsFixPayTimes.begin(), itsFixPayTimes.end() );

	/// sort and remove duplicates!
    CC_NS( std, sort )( result.begin(), result.end() );
    CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
    result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////  (SPREAD) //////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_Spread
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Spread::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Spread
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Spread::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// checking of the size and type
		GPAF_CheckArgSize( arg, 6, ARM_GP_Spread::itsFuncName  );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Spread::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_Spread::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_Spread::itsFuncName );
		GPAF_CheckArgType( arg[3], GFAT_DATEORMATU_TYPE,ARM_GP_Spread::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_DOUBLE_TYPE,	ARM_GP_Spread::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_DOUBLE_TYPE,	ARM_GP_Spread::itsFuncName );

		/// grab inputs
		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();
        
		itsCoeff1				= arg[4].GetDouble();
		itsCoeff2				= arg[5].GetDouble();

		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for SPREAD keyword");
		
		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_Spread::itsFuncName );

		ARM_Currency* ccy		    = mod->GetCurrency( itsCurveName );
		char* ccyName			    = ccy->GetCcyName();

		ARM_INDEX_TYPE defaultIndex =  GetDefaultIndexFromCurrency( ccyName );
		
		char fixCalendar[100];
		ccy->CalcFixPayCal(fixCalendar);

		char* floatPayCalendar	= ccy->GetPayCalName(defaultIndex);
		char* floatResetCalendar= ccy->GetResetCalName(defaultIndex);
		CC_NS(std,auto_ptr)<char> holdfloatPayCalendar(floatPayCalendar);
		CC_NS(std,auto_ptr)<char> holdfloatResetCalendar(floatResetCalendar);
		

/**/	int fixFreq	    = GetFixedFrequencyFtor(ccy)();
/**/	int fixDayCount = GetFixedDayCountFtor(ccy)();
/**/	int floatFreq = GetFloatFrequencyFtor(ccy)();
/**/	int floatDayCount = GetFloatDayCountFtor(ccy)();
/**/	int fwdDayCount	= GetFloatDayCountFtor(ccy)();

/// First Swap Rate
		ARM_Date endDate1 = ComputeEndDate( startDate, arg[2], fixCalendar );

		//FIXLEG
		ARM_DateStrip fixDateStrip1 = ARM_DateStrip( startDate, endDate1, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );

		/// copy constructor
		itsFixResetTimes1	= ARM_GP_Vector( *(fixDateStrip1.GetResetDates()) );  
		itsFixPayTimes1	    = ARM_GP_Vector( *(fixDateStrip1.GetPaymentDates()) );
		itsFixPayPeriods1   = ARM_GP_Vector( *(fixDateStrip1.GetInterestTerms()) );

		/// get time from date!
		for( int i(0); i<itsFixPayTimes1.size(); ++i )
		{
			itsFixResetTimes1[i] = mod->GetTimeFromDate( itsFixResetTimes1[i] );
			itsFixPayTimes1[i] = mod->GetTimeFromDate( itsFixPayTimes1[i] );
		}
	
		//FLOATLEG
		ARM_DateStrip floatDateStrip1 = ARM_DateStrip( startDate, endDate1, floatFreq, floatDayCount, floatResetCalendar,
											K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, floatFreq, GETDEFAULTVALUE,
											floatPayCalendar );
		
		/// copy constructor
		itsFwdResetTimes1	= ARM_GP_Vector( *(floatDateStrip1.GetResetDates()) );
		itsFwdStartTimes1	= ARM_GP_Vector( *(floatDateStrip1.GetFwdRateStartDates()) );
		itsFwdEndTimes1     = ARM_GP_Vector( *(floatDateStrip1.GetFwdRateEndDates()) );
		itsFwdPayPeriods1   = ARM_GP_Vector( itsFwdStartTimes1.size());
		itsFloatPayTimes1	= ARM_GP_Vector( *(floatDateStrip1.GetPaymentDates()) );
		itsFloatPayPeriods1 = ARM_GP_Vector( *(floatDateStrip1.GetInterestTerms()) );

		/// get time from date! 
		for(i=0; i<itsFwdStartTimes1.size(); ++i )
		{
		    itsFwdPayPeriods1[i] = CountYearsWithoutException( fwdDayCount, ARM_Date(itsFwdStartTimes1[i]), ARM_Date(itsFwdEndTimes1[i]) );
			itsFwdResetTimes1[i] = mod->GetTimeFromDate( itsFwdResetTimes1[i] );
			itsFwdStartTimes1[i] = mod->GetTimeFromDate( itsFwdStartTimes1[i] );
			itsFwdEndTimes1[i]   = mod->GetTimeFromDate( itsFwdEndTimes1[i] );
			itsFloatPayTimes1[i] = mod->GetTimeFromDate( itsFloatPayTimes1[i] );
		}

		ARM_Date floatStartDate1( startDate );
		itsFloatStartTime1   = mod->GetTimeFromDate( floatStartDate1 );
		ARM_Date floatEndDate1( endDate1 );
		itsFloatEndTime1	= mod->GetTimeFromDate( floatEndDate1 );

		/// Case of Margin Vector
		size_t sizefloatflows1 = itsFloatPayTimes1.size();        

		double margin = 0.;
		itsMarginVector1 = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows1,margin));
		
		/// Validation
		CheckNbSmaller( itsEvalTime,	itsFloatStartTime1,			"EvalTime",			"FloatStartTime1",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime1, itsFloatEndTime1,		"FloatStartTime1",	"FloatEndTime1",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime1, itsFixPayTimes1[0] ,	"FloatStartTime1",  "1st Fix PayTimes1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );

        CheckNbSmaller( itsEvalTime,	itsFwdStartTimes1[0],		"EvalTime",				"1st Fwd StartTimes1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTimes1[0], itsFixPayTimes1[0] ,	"1st Fwd StartTimes1",	"1st Fix PayTimes1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );

		CheckVectorIncreasing( itsFixPayTimes1,		"FixPayTimes1",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFixPayPeriods1,	"FixPayPeriods1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFwdStartTimes1,	"FwdStartTimes1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes1,		"FwdEndTimes1",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFloatPayTimes1,	"FloatPayTimes1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFloatPayPeriods1,	"FixPayPeriods1",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );

/// Second Swap Rate
		ARM_Date endDate2 = ComputeEndDate( startDate, arg[3], fixCalendar );

		//FIXLEG
		ARM_DateStrip fixDateStrip2 = ARM_DateStrip( startDate, endDate2, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );

		/// copy constructor
		itsFixResetTimes2	= ARM_GP_Vector( *(fixDateStrip2.GetResetDates()) );  
		itsFixPayTimes2	    = ARM_GP_Vector( *(fixDateStrip2.GetPaymentDates()) );
		itsFixPayPeriods2   = ARM_GP_Vector( *(fixDateStrip2.GetInterestTerms()) );

		/// get time from date!
		for(i=0; i<itsFixPayTimes2.size(); ++i )
		{
			itsFixResetTimes2[i] = mod->GetTimeFromDate( itsFixResetTimes2[i] );
			itsFixPayTimes2[i] = mod->GetTimeFromDate( itsFixPayTimes2[i] );
		}
	
		//FLOATLEG
		ARM_DateStrip floatDateStrip2 = ARM_DateStrip( startDate, endDate2, floatFreq, floatDayCount, floatResetCalendar,
											K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, floatFreq, GETDEFAULTVALUE,
											floatPayCalendar );
		
		/// copy constructor
		itsFwdResetTimes2	= ARM_GP_Vector( *(floatDateStrip2.GetResetDates()) );
		itsFwdStartTimes2	= ARM_GP_Vector( *(floatDateStrip2.GetFwdRateStartDates()) );
		itsFwdEndTimes2     = ARM_GP_Vector( *(floatDateStrip2.GetFwdRateEndDates()) );
		itsFwdPayPeriods2   = ARM_GP_Vector(itsFwdStartTimes2.size());
		itsFloatPayTimes2	= ARM_GP_Vector( *(floatDateStrip2.GetPaymentDates()) );
		itsFloatPayPeriods2 = ARM_GP_Vector( *(floatDateStrip2.GetInterestTerms()) );

		/// get time from date! 
		for(i=0; i<itsFwdStartTimes2.size(); ++i )
		{
		    itsFwdPayPeriods2[i] = CountYearsWithoutException( fwdDayCount, ARM_Date(itsFwdStartTimes2[i]), ARM_Date(itsFwdEndTimes2[i]) );
			itsFwdResetTimes2[i] = mod->GetTimeFromDate( itsFwdResetTimes2[i] );
			itsFwdStartTimes2[i] = mod->GetTimeFromDate( itsFwdStartTimes2[i] );
			itsFwdEndTimes2[i]   = mod->GetTimeFromDate( itsFwdEndTimes2[i] );
			itsFloatPayTimes2[i] = mod->GetTimeFromDate( itsFloatPayTimes2[i] );
		}

		ARM_Date floatStartDate2( startDate );
		itsFloatStartTime2   = mod->GetTimeFromDate( floatStartDate2 );
		ARM_Date floatEndDate2( endDate2 );
		itsFloatEndTime2	= mod->GetTimeFromDate( floatEndDate2 );

		/// Case of Margin Vector
		size_t sizefloatflows2 = itsFloatPayTimes2.size();        

		itsMarginVector2 = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows2,margin));
		
		/// Validation
		CheckNbSmaller( itsEvalTime,	itsFloatStartTime2,			"EvalTime",			"FloatStartTime2",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime2, itsFloatEndTime2,		"FloatStartTime2",	"FloatEndTime2",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime2, itsFixPayTimes2[0] ,	"FloatStartTime2",  "1st Fix PayTimes2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );

        CheckNbSmaller( itsEvalTime,	itsFwdStartTimes2[0],		"EvalTime",				"1st Fwd StartTimes2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTimes2[0], itsFixPayTimes2[0] ,	"1st Fwd StartTimes2",	"1st Fix PayTimes2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__ );

		CheckVectorIncreasing( itsFixPayTimes2,		"FixPayTimes2",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFixPayPeriods2,	"FixPayPeriods2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFwdStartTimes2,	"FwdStartTimes2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes2,		"FwdEndTimes2",		ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFloatPayTimes2,	"FloatPayTimes2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFloatPayPeriods2,	"FixPayPeriods2",	ARM_GP_Spread::itsFuncName,__LINE__,__FILE__  );


		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Spread
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Spread function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Spread::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	ARM_VectorPtr Spread =  itsModelIR->Spread( itsCurveName, itsEvalTime,
		itsCoeff1, 
        itsFloatStartTime1, itsFloatEndTime1,
		itsFixPayTimes1, itsFixPayPeriods1,
        itsFwdStartTimes1, itsFwdEndTimes1, itsFwdPayPeriods1,
        itsFloatPayTimes1, itsFloatPayPeriods1, *itsMarginVector1,
		itsCoeff2,
		itsFloatStartTime2, itsFloatEndTime2,
		itsFixPayTimes2, itsFixPayPeriods2,
        itsFwdStartTimes2, itsFwdEndTimes2, itsFwdPayPeriods2,
        itsFloatPayTimes2, itsFloatPayPeriods2, *itsMarginVector2,
		states );

	return ARM_GramFctorArg(Spread);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Spread
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Spread::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

    CC_STL_VECTOR( double ) result;

    result.push_back( itsFloatStartTime1 );
    result.insert( result.end(), itsFwdStartTimes1.begin(), itsFwdStartTimes1.end() );
    result.insert( result.end(), itsFwdEndTimes1.begin(), itsFwdEndTimes1.end() );
    result.insert( result.end(), itsFloatPayTimes1.begin(), itsFloatPayTimes1.end() );
    result.push_back( itsFloatEndTime1 );
    result.insert( result.end(), itsFixPayTimes1.begin(), itsFixPayTimes1.end() );

	result.push_back( itsFloatStartTime2 );
    result.insert( result.end(), itsFwdStartTimes2.begin(), itsFwdStartTimes2.end() );
    result.insert( result.end(), itsFwdEndTimes2.begin(), itsFwdEndTimes2.end() );
    result.insert( result.end(), itsFloatPayTimes2.begin(), itsFloatPayTimes2.end() );
    result.push_back( itsFloatEndTime2 );
    result.insert( result.end(), itsFixPayTimes2.begin(), itsFixPayTimes2.end() );

	/// sort and remove duplicates!
    CC_NS( std, sort )( result.begin(), result.end() );
    CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
    result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_BasisSwap
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////
long ARM_GP_BasisSwap::itsNbArg = 27;

void ARM_GP_BasisSwap::SetDefaults( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, ARM_GP_BasisSwap::itsNbArg, ARM_GP_BasisSwap::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_BasisSwap::itsFuncName );
	GPAF_CheckArgType( arg[1], GFAT_STRING_TYPE, ARM_GP_BasisSwap::itsFuncName );
	GPAF_CheckArgType( arg[2], GFAT_STRING_TYPE, ARM_GP_BasisSwap::itsFuncName );

	itsDomCurveName  = arg[0].GetString();
	itsForCurveName  = arg[1].GetString();
	itsFxCurveName	 = arg[2].GetString();

    ARM_Currency* forccy = mod->GetCurrency( itsForCurveName );
	ARM_Currency* domccy = mod->GetCurrency( itsDomCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(nodes, arg, 4, mod, itsDomCurveName );
    bool isDefault,isDbleNotional;
	isDefault=SetFrequency(   nodes, arg, 6, mod, itsDomCurveName, GetFixedFrequencyFtor(domccy));
	isDefault=SetDayCount(	nodes, arg, 7, mod, itsDomCurveName, GetFixedDayCountFtor(forccy));

    /// If default values are used set the double notional pricing method flag
	isDbleNotional=SetFrequency(   nodes, arg, 8, mod, itsForCurveName, GetFloatFrequencyFtor(forccy));
	isDefault=SetDayCount(	nodes, arg, 9, mod, itsForCurveName, GetFloatDayCountFtor(forccy));
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_BasisSwap
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_BasisSwap::GrabInputs( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[15], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_BasisSwap::itsFuncName );
	GPAF_CheckArgType( arg[17], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_BasisSwap::itsFuncName );
	
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_BasisSwap::itsFuncName );

		double asOfDate = mod->GetAsOfDate().GetJulian();
		/// checking of the size and type
		GPAF_CheckArgSize( arg, ARM_GP_BasisSwap::itsNbArg,			ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[0],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[1],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[2],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[3],  GFAT_DATE_TYPE,					ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[4],  GFAT_DATEORMATU_TYPE,			ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[5],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[6],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[7],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[8],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[9],  GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[10], GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[11], GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[12], GFAT_DATE_OR_DOUBLE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[13], GFAT_DATE_OR_DOUBLE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[14], GFAT_VECTOR_OR_CURVE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[16], GFAT_VECTOR_OR_CURVE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[18], GFAT_VECTOR_OR_CURVE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[19], GFAT_VECTOR_OR_CURVE_TYPE,		ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[20], GFAT_STRING_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[21], GFAT_DATESTRIP_TYPE,			ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[22], GFAT_DATESTRIP_TYPE,			ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[23], GFAT_DOUBLE_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[24], GFAT_DOUBLE_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[25], GFAT_DOUBLE_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		GPAF_CheckArgType( arg[26], GFAT_DOUBLE_TYPE,				ARM_GP_BasisSwap::itsFuncName );
		

		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR	= dynamic_cast< ARM_MultiAssetsModel* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only multi asset model can price BASISSWAP keyword");

		/// gram inputs
		itsDomCurveName			= arg[0].GetString();
		itsForCurveName			= arg[1].GetString();
		itsFxCurveName			= arg[2].GetString();				

		///First consistency validation between models
		ARM_Currency* domCcy	= mod->GetCurrency( itsDomCurveName );
		ARM_Currency* forCcy	= mod->GetCurrency( itsForCurveName );
		string domCcyName(domCcy->GetCcyName());
		string forCcyName(forCcy->GetCcyName());

		itsIsDomFlottant =  arg[10].GetString() == "FIXED" ? false : true;
		itsIsForFlottant =  arg[11].GetString() == "FIXED" ? false : true;

		itsExNotionalType = arg[20].GetString();
		stringToUpper(itsExNotionalType);

		ARM_Date startDate		= arg[3].GetDate();
		/// the second argument can either be a date or a maturity
		char* calendar	= domCcy->GetCcyName();
		//CC_NS(std,auto_ptr)<char> holdCalendar(calendar);
		ARM_Date endDate = ComputeEndDate( startDate, arg[4], calendar );

		size_t domDateOffSet = arg[23].GetDouble();
		int domDateWidth = static_cast<int>(arg[24].GetDouble());
		size_t forDateOffSet = arg[25].GetDouble();
		int forDateWidth = static_cast<int>(arg[26].GetDouble());

        /// Domestic leg data extraction, get the fixed leg convention
		ARM_DateStrip dateStrip;
		ARM_DateStripPtr domDateStrip = arg[21].GetType() == GFAT_DATESTRIP_TYPE ? arg[21].GetDateStrip() : ARM_DateStripPtr(NULL);
		if(domDateStrip.IsNull()){
			
			char domPayCalendar[100];
			itsIsDomFlottant ? domCcy->CalcFloatPayCal(domPayCalendar):domCcy->CalcFixPayCal(domPayCalendar);
			char domResetCalendar[100];
			domCcy->CalcResetCal(domResetCalendar);

			GetFrequencyFtor*  freqFtor  = &GetFloatFrequencyFtor(domCcy);
			if(!itsIsDomFlottant) freqFtor = &GetFixedFrequencyFtor(domCcy);

			int domFreq	= GetFrequency(arg[6], mod, itsDomCurveName, ARM_GP_BasisSwap::itsFuncName, *freqFtor);

			GetDayCountFtor*  DcFtor  = &GetFloatDayCountFtor(domCcy);
			if(!itsIsDomFlottant) DcFtor = &GetFixedDayCountFtor(domCcy);

			int domDayCount = GetDayCount( arg[7], mod, itsDomCurveName,*DcFtor);

			dateStrip = ARM_DateStrip( startDate, endDate, domFreq, domDayCount, domResetCalendar,K_MOD_FOLLOWING, 
				K_ADJUSTED, K_SHORTSTART,GETDEFAULTVALUE, domFreq, GETDEFAULTVALUE, domPayCalendar );
			domDateWidth = dateStrip.size();
		}
		else{
			dateStrip = *domDateStrip;
		
			// Resize and store the new dateStrip
			if(domDateWidth == -1)
				domDateWidth = dateStrip.size() - domDateOffSet;			
			dateStrip.ResizeAndBuilt(domDateOffSet, domDateOffSet + domDateWidth);
		}

		/// copy constructor
		itsDomResetTimes		= ARM_GP_Vector( *dateStrip.GetResetDates() )- asOfDate; 
		itsDomFwdStartTimes		= ARM_GP_Vector( *dateStrip.GetFwdRateStartDates() )- asOfDate; 
		itsDomFwdEndTimes		= ARM_GP_Vector( *dateStrip.GetFwdRateEndDates() )- asOfDate; 
		itsDomFlowStartTimes	= ARM_GP_Vector( *dateStrip.GetFlowStartDates() )- asOfDate;
		itsDomFlowEndTimes		= ARM_GP_Vector( *dateStrip.GetFlowEndDates() )- asOfDate;
		itsDomPayTimes			= ARM_GP_Vector( *dateStrip.GetPaymentDates() )- asOfDate;
		itsDomPayPeriods		= ARM_GP_Vector( *dateStrip.GetInterestTerms() );
								
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( domCcy->GetCcyName() );
		int fwdDayCount	= GetFloatDayCountFtor(domCcy)();
		itsDomFwdPeriods	= ARM_GP_Vector(itsDomFwdStartTimes.size()); 
		for(size_t i(0); i<itsDomFwdStartTimes.size(); ++i )
			itsDomFwdPeriods[i]	= CountYearsWithoutException( fwdDayCount, ARM_Date((*dateStrip.GetFwdRateStartDates())[i]), ARM_Date((*dateStrip.GetFwdRateEndDates())[i]) );

		/// Foreign leg data extraction, get the foreign leg convention
		ARM_DateStripPtr forDateStrip = arg[22].GetType() == GFAT_DATESTRIP_TYPE ? arg[22].GetDateStrip() : ARM_DateStripPtr(NULL);

		if(forDateStrip.IsNull()){

			GetFrequencyFtor*  freqFtor  = &GetFloatFrequencyFtor(forCcy);
			if(!itsIsForFlottant) freqFtor = &GetFixedFrequencyFtor(forCcy);

			int forFreq	= GetFrequency(arg[8], mod, itsForCurveName, ARM_GP_BasisSwap::itsFuncName, *freqFtor);

			GetDayCountFtor*  DcFtor  = &GetFloatDayCountFtor(forCcy);
			if(!itsIsDomFlottant) DcFtor = &GetFixedDayCountFtor(forCcy);

			int forDayCount = GetDayCount( arg[9], mod, itsForCurveName,*DcFtor);

			/// the logic here is to get the default calendar for the float leg
			/// these are from the default index
			/// the fix calendar is the currency itself!
			char forPayCalendar[100];
			itsIsForFlottant ? forCcy->CalcFloatPayCal(forPayCalendar):forCcy->CalcFixPayCal(forPayCalendar);
			char forResetCalendar[100];
			forCcy->CalcResetCal(forResetCalendar);

		/// foreign leg datestrip
			dateStrip = ARM_DateStrip( startDate, endDate, forFreq, forDayCount, forResetCalendar,	K_MOD_FOLLOWING,
				K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, forFreq, GETDEFAULTVALUE, forPayCalendar );

			forDateWidth = dateStrip.size();		
		}
		else{
			dateStrip = *forDateStrip;
		
			// Resize and store the new dateStrip
			if(forDateWidth == -1)
				forDateWidth = dateStrip.size() - forDateOffSet;			
			dateStrip.ResizeAndBuilt(forDateOffSet, forDateOffSet + forDateWidth);
		}
		/// Fx reset gap
		int fxResetGap = GetResetDaysGap(arg[12], mod, itsFxCurveName );
		int fxSettlGap = GetResetDaysGap(arg[13], mod, itsFxCurveName );
		char* calendarStr = const_cast<char*> (itsModelIR->GetSettlementCalendar(itsFxCurveName).c_str() );
		//CC_NS(std,auto_ptr)<char> holdCalendarStr(calendarStr);

		/// copy constructor and get time
		itsForResetTimes		= ARM_GP_Vector( *dateStrip.GetResetDates() )- asOfDate; 
		itsForFwdStartTimes		= ARM_GP_Vector( *dateStrip.GetFwdRateStartDates() )- asOfDate; 
		itsForFwdEndTimes		= ARM_GP_Vector( *dateStrip.GetFwdRateEndDates() )- asOfDate; 
		itsForFlowStartTimes	= ARM_GP_Vector( *dateStrip.GetFlowStartDates() )- asOfDate; 
		itsForFlowEndTimes		= ARM_GP_Vector( *dateStrip.GetFlowEndDates() )- asOfDate; 
		itsForPayTimes			= ARM_GP_Vector( *dateStrip.GetPaymentDates() )- asOfDate; 
		itsForPayPeriods		= ARM_GP_Vector( *dateStrip.GetInterestTerms() ); 
								
		defaultIndex = GetDefaultIndexFromCurrency( forCcy->GetCcyName() );
		fwdDayCount	= GetFloatDayCountFtor(forCcy)();
		size_t size = itsForFwdStartTimes.size();
		itsFxResetTimes.resize(size);
		itsFxSettlTimes.resize(size);
		itsForFwdPeriods.resize(size);
		for(i =0; i<size; ++i ){
			itsForFwdPeriods[i]	= CountYearsWithoutException( fwdDayCount, ARM_Date((*dateStrip.GetFwdRateStartDates())[i]), ARM_Date((*dateStrip.GetFwdRateEndDates())[i]) );
			ARM_Date tempDate = ARM_Date((*dateStrip.GetFwdRateStartDates())[i]);
			tempDate.GapBusinessDay(-abs(fxResetGap), calendarStr ); /// UGLY but forced to const cast
			itsFxResetTimes[i] = tempDate.GetJulian()-asOfDate;
			tempDate.GapBusinessDay(abs(fxSettlGap), calendarStr); /// UGLY but forced to const cast
			itsFxSettlTimes[i] = tempDate.GetJulian()- asOfDate;
		}
		itsFxPayTimes = itsForPayTimes;

		itsStartTime =  startDate.GetJulian() - asOfDate;
		itsEndTime = endDate.GetJulian() - asOfDate;
		
		ARM_GP_Vector Vect;
		/// Case of domestic margin Vector
		size = itsDomPayTimes.size();        
		if( arg[14].GetType() ==  GFAT_DOUBLE_TYPE )
			itsDomMarginVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(size,arg[14].GetDouble()));
		else if( arg[14].GetType() ==  GFAT_VECTOR_TYPE ){
			Vect	= *arg[14].GetVector();
			ARM_GP_Vector::iterator begin = Vect.begin(), iter;
			iter = Vect.erase(begin,begin + domDateOffSet);
			iter = Vect.erase(iter + domDateWidth, Vect.end());
			itsDomMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(Vect));
		}
		else if( arg[14].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr marginCurve	= arg[14].GetCurve();
			itsDomMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(size,0.0)); 
			for(int i=0; i<size; ++i)
				(*itsDomMarginVector)[i] = marginCurve->Interpolate(itsDomResetTimes[i]);										
		}
		/// Case of domestic  Notional Vector
		if( arg[18].GetType() ==  GFAT_DOUBLE_TYPE )
			itsDomNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(size,arg[18].GetDouble()));
		else if( arg[18].GetType() ==  GFAT_VECTOR_TYPE )
		{
			Vect	= *arg[18].GetVector();
			ARM_GP_Vector::iterator begin = Vect.begin(), iter;
			iter = Vect.erase(begin,begin + domDateOffSet);
			iter = Vect.erase(iter + domDateWidth, Vect.end());
			itsDomNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(Vect));
		}
		else if( arg[18].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr notionalCurve	= arg[18].GetCurve();
			itsDomNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(size,1.0)); 
			for(i=0; i<size; ++i)	
				(*itsDomNotionalVector)[i] = notionalCurve->Interpolate(itsDomPayTimes[i]);
		}

		/// Case of domestic strike matri if necessairy
		if(!itsIsDomFlottant){
		
			if( arg[15].GetType() ==  GFAT_DOUBLE_TYPE )
				itsDomStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,size,arg[15].GetDouble()));
			else if( arg[15].GetType() ==  GFAT_VECTOR_TYPE )
			{
				Vect = *arg[15].GetVector();
				ARM_GP_Vector::iterator begin = Vect.begin(), iter;
				iter = Vect.erase(begin,begin + domDateOffSet);
				iter = Vect.erase(iter + domDateWidth, Vect.end());
				(*itsDomStrikeMatrix).push_backRow(Vect);
			}
			else if( arg[15].GetType() ==  GFAT_CURVE_TYPE )
			{
				ARM_GP_CurvePtr strikeCurve	= arg[15].GetCurve();
				itsDomStrikeMatrix = ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,size,0.0));
				for(i=0; i<size; ++i)			
					(*itsDomStrikeMatrix)(0,i) = strikeCurve->Interpolate(itsDomResetTimes[i]);										
			}
		}
		/// Case of foreign margin Vector
		 size = itsForPayTimes.size();        
		if( arg[16].GetType() ==  GFAT_DOUBLE_TYPE )
			itsForMarginVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(size,arg[16].GetDouble()));
		else if( arg[16].GetType() ==  GFAT_VECTOR_TYPE ){
			Vect	= *arg[16].GetVector();
			ARM_GP_Vector::iterator begin = Vect.begin(), iter;
			iter = Vect.erase(begin,begin + forDateOffSet);
			iter = Vect.erase(iter + forDateWidth, Vect.end());
			itsForMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(Vect));
		}
		else if( arg[16].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr marginCurve	= arg[16].GetCurve();
			itsForMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(size,0.0)); 
			for(int i=0; i<size; ++i)
				(*itsForMarginVector)[i] = marginCurve->Interpolate(itsForResetTimes[i]);										
		}
		/// Case of foreign Notional Vector
		if( arg[19].GetType() ==  GFAT_DOUBLE_TYPE )
			itsForNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(size,arg[19].GetDouble()));
		else if( arg[19].GetType() ==  GFAT_VECTOR_TYPE )
		{
			Vect	= *arg[19].GetVector();
			ARM_GP_Vector::iterator begin = Vect.begin(), iter;
			iter = Vect.erase(begin,begin + forDateOffSet);
			iter = Vect.erase(iter + forDateWidth, Vect.end());
			itsForNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(Vect));
		}
		else if( arg[19].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr notionalCurve	= arg[19].GetCurve();
			itsForNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(size,1.0)); 
			for(i=0; i<size; ++i)	
				(*itsForNotionalVector)[i] = notionalCurve->Interpolate(itsForPayTimes[i]);
		}
		/// Case of domestic strike matri if necessairy
		if(!itsIsForFlottant){
		
			if( arg[17].GetType() ==  GFAT_DOUBLE_TYPE )
				itsForStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,size,arg[17].GetDouble()));
			else if( arg[17].GetType() ==  GFAT_VECTOR_TYPE )
			{
				Vect = *arg[17].GetVector();
				ARM_GP_Vector::iterator begin = Vect.begin(), iter;
				iter = Vect.erase(begin,begin + forDateOffSet);
				iter = Vect.erase(iter + forDateWidth, Vect.end());
				(*itsForStrikeMatrix).push_backRow(Vect);
			}
			else if( arg[17].GetType() ==  GFAT_CURVE_TYPE )
			{
				ARM_GP_CurvePtr strikeCurve	= arg[17].GetCurve();
				itsForStrikeMatrix = ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,size,0.0));
				for(i=0; i<size; ++i)			
					(*itsForStrikeMatrix)(0,i) = strikeCurve->Interpolate(itsForResetTimes[i]);										
			}
		}

		/// Validation
		CheckNbSmaller( itsEvalTime,  itsStartTime,		"EvalTime",			"StartTime",	ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsStartTime, itsEndTime,		"StartTime",		"EndTime",		ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsStartTime, itsDomPayTimes[0],"FloatStartTime",   "1st Fix PayTimes",	ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__ );

        CheckNbSmaller( itsEvalTime,	itsForFwdStartTimes[0],		"EvalTime",			    "1st Fwd StartTimes",	ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsForFwdStartTimes[0], itsDomPayTimes[0] ,	"1st Fwd StartTimes",	"1st Fix PayTimes",		ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__ );

		CheckVectorPositiveNb( itsForPayPeriods,	"FixPayPeriods",		ARM_GP_BasisSwap::itsFuncName,__LINE__,__FILE__  );

		//Set default for node
		SetDefaults( arg, mod, nodes );
	
		string payRecString	= arg[5].GetString();
		itsPayRec = ARM_ArgConv_PayRec.GetNumber( payRecString );

		SetAlreadyComputed(true);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_BasisSwap
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Cap Function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_BasisSwap::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
		
	/// Conversion to matriciel strikes
	if(!itsIsDomFlottant){
		size_t nbRows = itsDomStrikeMatrix->GetRowsNb();
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		
		if( nbRows != statesSize && nbRows !=1 )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " dom strikeMatrix rows size is failed" );
		else if(nbRows == 1 &&  nbRows != statesSize)
		{
			ARM_GP_Vector* vectStrike = itsDomStrikeMatrix->GetRow(0);
			itsDomStrikeMatrix->reserve(statesSize, vectStrike->size() );
			for(size_t i=1; i<statesSize; ++i)		
				itsDomStrikeMatrix->push_backRow(*vectStrike);
			delete vectStrike;
		}
	}

	if(!itsIsForFlottant){
		size_t nbRows = itsForStrikeMatrix->GetRowsNb();
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		
		if( nbRows != statesSize && nbRows !=1 )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " for strikeMatrix rows size is failed" );
		else if(nbRows == 1 &&  nbRows != statesSize)
		{
			ARM_GP_Vector* vectStrike = itsForStrikeMatrix->GetRow(0);
			itsForStrikeMatrix->reserve(statesSize, vectStrike->size() );
			for(size_t i=1; i<statesSize; ++i)		
				itsForStrikeMatrix->push_backRow(*vectStrike);
			delete vectStrike;
		}
	}
	
	ARM_VectorPtr NPVFloatLeg =  itsModelIR->NPVBasisSwap( 
		itsDomCurveName,
		itsForCurveName,
		itsFxCurveName,
		itsEvalTime,
		itsStartTime,
		itsEndTime,
		itsPayRec,
		itsDomResetTimes,	    
		itsDomFwdStartTimes,
		itsDomFwdEndTimes,
		itsDomFlowStartTimes,			
		itsDomFlowEndTimes,	
		itsDomFwdPeriods,	
		itsDomPayTimes,
		itsDomPayPeriods,
		*itsDomMarginVector,
		*itsDomNotionalVector,
		itsIsDomFlottant,
		itsForResetTimes,       
		itsForFwdStartTimes,
		itsForFwdEndTimes,
		itsForFlowStartTimes,   		
		itsForFlowEndTimes,	    
		itsForFwdPeriods,	
		itsForPayTimes,
		itsForPayPeriods,
		*itsForMarginVector,
		*itsForNotionalVector,
		itsIsForFlottant,
		itsExNotionalType,
		itsFxResetTimes,
		itsFxSettlTimes,  
		itsFxPayTimes,
		*itsDomStrikeMatrix,
		*itsForStrikeMatrix,
		states );

	return ARM_GramFctorArg(NPVFloatLeg);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_BasisSwap
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_BasisSwap::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	ARM_GP_Vector result;
    result.push_back( itsStartTime );

	if(itsIsForFlottant){
		result.insert( result.end(), itsForFwdStartTimes.begin(), itsForFwdStartTimes.end() );
		result.insert( result.end(), itsForFwdEndTimes.begin(), itsForFwdEndTimes.end() );
		result.insert( result.end(), itsForPayTimes.begin(), itsForPayTimes.end() );
	}
	if(itsIsDomFlottant){
		result.insert( result.end(), itsDomFwdStartTimes.begin(), itsDomFwdStartTimes.end() );
		result.insert( result.end(), itsDomFwdEndTimes.begin(), itsDomFwdEndTimes.end() );
	}

	result.push_back( itsEndTime );
    result.insert( result.end(), itsDomPayTimes.begin(), itsDomPayTimes.end() );

	/// sort and remove duplicates!
    CC_NS( std, sort )( result.begin(), result.end() );
    ARM_GP_Vector::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
    result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}

/////////////////////////////////////////////////////////
///////// Caplet	and digital caplet //////////////////
/////////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////
ARM_GP_CapDigital_Common::ARM_GP_CapDigital_Common( 
	const ARM_CapNDigital_PricingFunc& pricingFunc, 
	const string& funcName )
:	ARM_GramFctor(), 
	itsPricingFunc( pricingFunc ),
	itsFuncName( funcName )
{}

////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_CapDigital_Common::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 12, ARM_GP_CapDigital_Common::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_CapDigital_Common::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(		 nodes, arg, 2, mod, itsCurveName );
	bool isDefault=SetDayCount(	 nodes, arg, 5, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	SetResetDaysGap( nodes, arg, 6, mod, itsCurveName );

    SetIndexTerm( nodes, arg, 10, 1, 2, mod, itsCurveName, ARM_GP_CapDigital_Common::itsFuncName);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_CapDigital_Common::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	
	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,	ARM_GP_CapDigital_Common::itsFuncName );

	/// handle vectorial strike
	if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
	{
		itsStrikeDouble = arg[3].GetDouble();
		itsStrikeVector	= ARM_VectorPtr(NULL);
	}
	else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
	{
		const double UNASSIGNED = -1111;
		itsStrikeDouble =  UNASSIGNED;
		itsStrikeVector	= arg[3].GetVector();
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );

	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 12, ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,			    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,	    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_DATE_OR_DOUBLE_TYPE,    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[7], GFAT_DATE_OR_DOUBLE_TYPE,    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[10],GFAT_MATURITY_TYPE,		    ARM_GP_CapDigital_Common::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_STRING_TYPE,	        ARM_GP_CapDigital_Common::itsFuncName );

		/// grab inputs
		itsCurveName			    = arg[0].GetString();
		ARM_Date capletStartDate    = arg[1].GetDate();
		string capFloorString	    = arg[4].GetString();
		itsNotional				    = arg[8].GetDouble();
		string advArrString	        = arg[9].GetString();

		/// compute inputs
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for CAPLET/DIGITAL keywords");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex 
								= GetDefaultIndexFromCurrency( ccyName );

		char* payCalendar		= ccy->GetPayCalName(defaultIndex);
		char* resetCalendar		= ccy->GetResetCalName(defaultIndex); 


		/// The second argument can either be a date or a maturity
        /// then in case of a maturity, first of all add without
        /// startDate adjustment to get endDate
        /// Becareful : it must be consistent with ARM_GP_CapDigital_Common::SetDefaults()
        /// ***** Is it the RIGHT default convention ? *****
		ARM_Date capletEndDate	= ComputeEndDate( capletStartDate, arg[2], payCalendar );

        /// Compute reset date according to timing & reset gap
		int resetTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
        ARM_Date capletResetDate( ComputeResetOrPayDate(capletStartDate,capletEndDate,resetTiming,
                        resetCalendar,arg[6],mod,itsCurveName,ARM_GP_CapDigital_Common::itsFuncName) );


        /// Compute fwd start & end dates according to reset date and standard reset gap
        /// (no caplet on forward/forward rate allowed !)
        int stdSpotGap = ccy->GetSpotDays();
        ARM_Date fwdStartDate(capletResetDate);
        fwdStartDate.GapBusinessDay(stdSpotGap,resetCalendar);
        string indexTerm( ComputeIndexTerm(arg[10],arg[1],arg[2],ARM_GP_CapDigital_Common::itsFuncName) );
        ARM_Date fwdEndDate( ComputeEndDate( fwdStartDate, indexTerm, resetCalendar ) );

		/// Compute payment date according to timing & payment gap
		int payTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
        ARM_Date capletPayDate( ComputeResetOrPayDate(capletStartDate,capletEndDate,K_ARREARS,payCalendar,
                        arg[7],mod,itsCurveName,ARM_GP_CapDigital_Common::itsFuncName) );


		itsCapFloor		    = ARM_ArgConv_CapFloor.GetNumber( capFloorString );
		int cpnDayCount     = GetDayCount( arg[5], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		itsPeriod		    = CountYearsWithoutException( cpnDayCount, capletStartDate, capletEndDate );
		itsEvalTime		    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_CapDigital_Common::itsFuncName );
		itsFwdStartTime	    = mod->GetTimeFromDate( fwdStartDate );
		itsFwdEndTime	    = mod->GetTimeFromDate( fwdEndDate );
		int indexDayCount   = GetDayCount( arg[11], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		itsFwdPeriod        = CountYearsWithoutException( indexDayCount, fwdStartDate, fwdEndDate );
		itsPayTime		    = mod->GetTimeFromDate( capletPayDate );
		itsResetTime	    = mod->GetTimeFromDate( capletResetDate);

		/// some validation
		CheckNbSmaller( itsEvalTime, itsResetTime,      "EvalTime",     "ResetTime", ARM_GP_CapDigital_Common::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTime,itsFwdEndTime,  "FwdStartTime", "FwdEndTime",	ARM_GP_CapDigital_Common::itsFuncName,__LINE__,__FILE__  );
		CheckPositiveNb( itsPeriod,     "Cpn Period", ARM_GP_CapDigital_Common::itsFuncName,__LINE__,__FILE__  );
		CheckPositiveNb( itsFwdPeriod,  "Index Period", ARM_GP_CapDigital_Common::itsFuncName,__LINE__,__FILE__  );

		/// delete char* for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}




////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Swap Rate function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_CapDigital_Common::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	
	/// Conversion to vectorial strikes
	if( ARM_VectorPtr(NULL) == itsStrikeVector)
	{
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		itsStrikeVector		= ARM_VectorPtr( new ARM_GP_Vector( statesSize, itsStrikeDouble ) );
	}

	/// call function
	ARM_VectorPtr CapFloorLet =  (itsModelIR->*itsPricingFunc)( itsCurveName, itsEvalTime, itsPayTime, itsPeriod,
		itsNotional, itsResetTime, itsFwdStartTime, itsFwdEndTime, itsFwdPeriod, *itsStrikeVector, itsCapFloor, states );

	/// return result
	return ARM_GramFctorArg(CapFloorLet);
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_CapDigital_Common::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	
	/// return result
	vector<double> tmpResult;
	tmpResult.reserve(2);
	
	tmpResult.push_back( itsResetTime );
	if( itsResetTime != itsFwdStartTime )
		tmpResult.push_back(  itsFwdStartTime );
	
	tmpResult.push_back( itsFwdEndTime );
    if( itsPayTime != itsFwdEndTime )
		tmpResult.push_back( itsPayTime );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( tmpResult.size(), tmpResult.begin() )), ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////		Swaption		  //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_Swaption
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Swaption::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 9, ARM_GP_Swaption::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Swaption::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(			nodes, arg, 2, mod, itsCurveName );
    bool isDefault;
	isDefault=SetFrequency(       nodes, arg, 5, mod, itsCurveName, GetFixedFrequencyFtor(ccy));
	isDefault=SetDayCount(	    nodes, arg, 6, mod, itsCurveName, GetFixedDayCountFtor(ccy));
	SetResetDaysGap(	nodes, arg, 7, mod, itsCurveName );
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Swaption
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Swaption::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,	ARM_GP_Swaption::itsFuncName );

		/// this is always reasserted!
	GPAF_CheckArgType( arg[8], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Swaption::itsFuncName );
	
	if(!GetAlreadyComputed())
	{
		/// first setDefaults
		SetDefaults( arg, mod, nodes );
		
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 9, ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,	        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,	        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,	        ARM_GP_Swaption::itsFuncName );
		GPAF_CheckArgType( arg[7], GFAT_DATE_OR_DOUBLE_TYPE,    ARM_GP_Swaption::itsFuncName );

		/// grab inputs
		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();
		string payRecString		= arg[4].GetString();
		string dayCountString	= arg[6].GetString();


		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for SWAPTION keyword");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex 
								= GetDefaultIndexFromCurrency( ccyName );

		/// the logic here is to get the default calendar for the float leg
		/// these are from the default index
		/// the fix calendar is the currency itself!
		char* floatPayCalendar	= ccy->GetPayCalName(defaultIndex);
		char* floatResetCalendar= ccy->GetResetCalName(defaultIndex);

		char fixCalendar[100];
        ccy->CalcFixPayCal(fixCalendar);

		/// the second argument can either be a date or a maturity
		ARM_Date endDate = ComputeEndDate( startDate, arg[2], fixCalendar );
		
		int fixFreq		 = GetFrequency(arg[5], mod, itsCurveName, ARM_GP_Swaption::itsFuncName, GetFixedFrequencyFtor(ccy));
		int fixDayCount  = GetDayCount( arg[6], mod, itsCurveName, GetFixedDayCountFtor(ccy));

        /// the convention is to set the option type on the rate :
        ///     payer swaption      =>  call on the swap rate
        ///     receiver swaption   =>  put on the swap rate
		if(ARM_ArgConv_PayRec.GetNumber( payRecString ) == K_PAY)
            itsCallPut = K_CALL;
        else
            itsCallPut = K_PUT;
        

		/// fix leg datestrip
		ARM_DateStrip fixDateStrip( startDate, endDate, fixFreq, fixDayCount, fixCalendar,
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
			fixCalendar );

		/// adjust according to modified following the dates
		ARM_Date floatStartDate( startDate );
		ARM_Date floatEndDate( endDate );

		int floatDayCount = GetFloatDayCountFtor(ccy)();
		int floatFreq = GetFloatFrequencyFtor(ccy)();

		// The swap reset date is the start date ajusted on the payment calendar minus
		// the reset days gap on the reset calendar.
        ARM_Date swapResetDate( ComputeResetOrPayDate(floatStartDate,floatEndDate,K_ADVANCE,
                        floatResetCalendar,arg[7],mod,itsCurveName,ARM_GP_Swaption::itsFuncName) );

		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_Swaption::itsFuncName );

		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pfixPayTimes   = fixDateStrip.GetFlowEndDates();
		ARM_GP_Vector* pfixPayPeriods = fixDateStrip.GetInterestTerms();

		/// copy constructor
		itsFixPayTimes	= ARM_GP_Vector( *pfixPayTimes );
		itsFixPayPeriods= ARM_GP_Vector( *pfixPayPeriods );
		
		/// float leg datestrip  !!! to be changed to take into account the float freq, notional
		ARM_DateStrip floatDateStrip( startDate, endDate, floatFreq, floatDayCount, floatResetCalendar,
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, floatFreq, GETDEFAULTVALUE,
			floatPayCalendar );


		ARM_GP_Vector* pfloatResetTimes	= floatDateStrip.GetResetDates();
		ARM_GP_Vector* pfloatStartTimes	= floatDateStrip.GetFlowStartDates();
		ARM_GP_Vector* pfloatEndTimes	= floatDateStrip.GetFlowEndDates();
		ARM_GP_Vector* pfloatIntTerms   = floatDateStrip.GetInterestTerms();

		itsFloatResetTimes	= ARM_GP_Vector( *pfloatResetTimes );
		itsFloatStartTimes= ARM_GP_Vector( *pfloatStartTimes );
		itsFloatPayTimes	= ARM_GP_Vector( *pfloatEndTimes );
		itsFloatIntTerms= ARM_GP_Vector( *pfloatIntTerms );

		/// get time from date! 
		int i;

		for(i=0; i<pfixPayTimes->size(); ++i )
		{
			itsFixPayTimes[i] = mod->GetTimeFromDate( itsFixPayTimes[i] );
		}

		for(i=0; i<pfloatResetTimes->size(); ++i )
		{
			itsFloatResetTimes[i] = mod->GetTimeFromDate( itsFloatResetTimes[i] );
			itsFloatStartTimes[i] = mod->GetTimeFromDate( itsFloatStartTimes[i] );
			itsFloatPayTimes[i] = mod->GetTimeFromDate( itsFloatPayTimes[i] );
		}

		/// handle reset Days!
		itsSwapResetTime    = mod->GetTimeFromDate( swapResetDate );
		itsFloatStartTime   = mod->GetTimeFromDate( floatStartDate );
		itsFloatEndTime		= mod->GetTimeFromDate( floatEndDate );

		/// Validation
		CheckNbSmaller( itsEvalTime,		itsSwapResetTime,	"EvalTime",			"SwapResetTime",	ARM_GP_Swaption::itsFuncName,__LINE__,__FILE__  );
		CheckNbSmaller( itsFloatStartTime, itsFloatEndTime,		"FloatStartTime",	"FloatEndTime",		ARM_GP_Swaption::itsFuncName,__LINE__,__FILE__  );
		CheckNbSmaller( itsFloatStartTime, itsFixPayTimes[0] ,	"FloatStartTime",   "1st Fix PayTimes",	ARM_GP_Swaption::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFixPayTimes,	"FixPayTimes",	ARM_GP_Swaption::itsFuncName,__LINE__,__FILE__);
		CheckVectorPositiveNb( itsFixPayPeriods,"FixPayPeriods",ARM_GP_Swaption::itsFuncName );

        /// handle vectorial strike
        size_t sizefixflows = itsFixPayTimes.size();        
	    if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
	    {
		    double strike = arg[3].GetDouble();
		    itsStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizefixflows,strike));
	    }
	    else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
	    {
		    ARM_GP_VectorPtr strikePerState	= arg[3].GetVector();
            vector<double> strikes;
            int k =0;
            for(i=0; i<sizefixflows; ++i)
                for(int j=0; j<strikePerState->size();++k, ++j)
                    strikes[k] = (*strikePerState)[i];
            itsStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(strikePerState->size(),sizefixflows,strikes));

	    }
        else if( arg[3].GetType() ==  GFAT_CURVE_TYPE )
	    {
		    ARM_GP_CurvePtr strikeCurve	= arg[3].GetCurve();
            vector<double> strikeStepup(sizefixflows);
            for(i=0; i<sizefixflows; ++i)
                strikeStepup[i] = strikeCurve->Interpolate(itsFixPayTimes[i]);

            itsStrikeMatrix	= ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizefixflows,strikeStepup));
	    }


		/// Case of Notional Vector
		/// Only the fix Frequency is specified

		if( arg[8].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double notional = arg[8].GetDouble();
			itsFixNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,notional));
			itsFloatNotionalVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,notional));
			itsIsConstantNotional = true;

		}
		else if( arg[8].GetType() ==  GFAT_VECTOR_TYPE )
		{
			itsFixNotionalVector	= arg[8].GetVector();
			itsFloatNotionalVector	= arg[8].GetVector();
			itsIsConstantNotional = false;
		}
		else if( arg[8].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr notionalCurve	= arg[8].GetCurve();
			itsFixNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,1.0)); 
			itsFloatNotionalVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefixflows,1.0)); 
			for(i=0; i<sizefixflows; ++i)			
					(*itsFixNotionalVector)[i] = notionalCurve->Interpolate(itsFixPayTimes[i]);
			for(i=0; i<sizefixflows; ++i)			
					(*itsFloatNotionalVector)[i] = notionalCurve->Interpolate(itsFloatPayTimes[i]);
			itsIsConstantNotional = false;
								
		}

	    else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be double, vector or ARM_Curve" );




		/// delete char* for memory leak
		delete floatPayCalendar;
		delete floatResetCalendar;

		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Swaption
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Swaption function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Swaption::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	
	/// Conversion to matriciel strikes
	size_t nbRows = itsStrikeMatrix->GetRowsNb();
    size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
	
	if( nbRows != statesSize && nbRows !=1 )
	{
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + " strikeMatrix rows size is failed" );
    }
	else if(nbRows == 1 &&  nbRows != statesSize)
	{
		ARM_GP_Vector* vectStrike = itsStrikeMatrix->GetRow(0);
		itsStrikeMatrix->reserve(statesSize, vectStrike->size() );
		for(size_t i=1; i<statesSize; ++i)		
			itsStrikeMatrix->push_backRow(*vectStrike);
		delete vectStrike;
	}
	
	ARM_VectorPtr Swaption =  itsModelIR->VanillaSwaption( itsCurveName, itsEvalTime, itsSwapResetTime, *itsFixNotionalVector, *itsFloatNotionalVector,
		itsFloatStartTime, itsFloatEndTime,itsFloatResetTimes,itsFloatStartTimes,itsFloatPayTimes,itsFloatIntTerms, itsFixPayTimes, itsFixPayPeriods, *itsStrikeMatrix, itsCallPut, states, itsIsConstantNotional);

	return ARM_GramFctorArg(Swaption);
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_Swaption
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Swaption::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	/// return result
	vector<double> tmpResult;
	size_t sizeFixPayTimes = itsFixPayTimes.size();
	tmpResult.reserve( sizeFixPayTimes + 1 );
	tmpResult.push_back( itsSwapResetTime );
	if( itsSwapResetTime != itsFloatStartTime )
		tmpResult.push_back( itsFloatStartTime );

	int i;
	for(i=0; i< sizeFixPayTimes; ++i ){
		tmpResult.push_back( itsFixPayTimes[i] );
	}

	if( itsFloatEndTime != itsFixPayTimes[sizeFixPayTimes-1] )
		tmpResult.push_back( itsFloatEndTime  );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( tmpResult.size(), tmpResult.begin() )), ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///////////		Implied Vol		  //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_Swaption
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_ImpliedVol::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 12, ARM_GP_ImpliedVol::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_ImpliedVol::itsFuncName );
	itsCurveName = arg[0].GetString();
    ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(		 nodes, arg, 2, mod, itsCurveName );
	bool isDefault=SetDayCount(	 nodes, arg, 5, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	SetResetDaysGap( nodes, arg, 6, mod, itsCurveName );

    SetIndexTerm( nodes, arg, 10, 1, 2, mod, itsCurveName, ARM_GP_ImpliedVol::itsFuncName);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_ImpliedVol::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	
	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,	ARM_GP_ImpliedVol::itsFuncName );

	/// handle vectorial strike
	if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
	{
		itsStrikeDouble = arg[3].GetDouble();
		itsStrikeVector	= ARM_VectorPtr(NULL);
	}
	else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
	{
		const double UNASSIGNED = -1111;
		itsStrikeDouble =  UNASSIGNED;
		itsStrikeVector	= arg[3].GetVector();
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );

	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 12, ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,			    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,	    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_DATE_OR_DOUBLE_TYPE,    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[7], GFAT_DATE_OR_DOUBLE_TYPE,    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[10],GFAT_MATURITY_TYPE,		    ARM_GP_ImpliedVol::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_STRING_TYPE,	        ARM_GP_ImpliedVol::itsFuncName );

		/// grab inputs
		itsCurveName			    = arg[0].GetString();
		ARM_Date capletStartDate    = arg[1].GetDate();
		string capFloorString	    = arg[4].GetString();
		itsNotional				    = arg[8].GetDouble();
		string advArrString	        = arg[9].GetString();

		/// compute inputs
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for CAPLET/DIGITAL keywords");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex 
								= GetDefaultIndexFromCurrency( ccyName );

		char* payCalendar		= ccy->GetPayCalName(defaultIndex);
		char* resetCalendar		= ccy->GetResetCalName(defaultIndex); 


		/// The second argument can either be a date or a maturity
        /// then in case of a maturity, first of all add without
        /// startDate adjustment to get endDate
        /// Becareful : it must be consistent with ARM_GP_CapDigital_Common::SetDefaults()
        /// ***** Is it the RIGHT default convention ? *****
		ARM_Date capletEndDate	= ComputeEndDate( capletStartDate, arg[2], payCalendar );

        /// Compute reset date according to timing & reset gap
		int resetTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
        ARM_Date capletResetDate( ComputeResetOrPayDate(capletStartDate,capletEndDate,resetTiming,
                        resetCalendar,arg[6],mod,itsCurveName,ARM_GP_ImpliedVol::itsFuncName) );


        /// Compute fwd start & end dates according to reset date and standard reset gap
        /// (no caplet on forward/forward rate allowed !)
        int stdSpotGap = ccy->GetSpotDays();
        ARM_Date fwdStartDate(capletResetDate);
        fwdStartDate.GapBusinessDay(stdSpotGap,resetCalendar);
        string indexTerm( ComputeIndexTerm(arg[10],arg[1],arg[2],ARM_GP_ImpliedVol::itsFuncName) );
        ARM_Date fwdEndDate( ComputeEndDate( fwdStartDate, indexTerm, resetCalendar ) );

		/// Compute payment date according to timing & payment gap
		int payTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
        ARM_Date capletPayDate( ComputeResetOrPayDate(capletStartDate,capletEndDate,K_ARREARS,payCalendar,
                        arg[7],mod,itsCurveName,ARM_GP_ImpliedVol::itsFuncName) );


		itsCapFloor		    = ARM_ArgConv_CapFloor.GetNumber( capFloorString );
		int cpnDayCount     = GetDayCount( arg[5], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		itsPeriod		    = CountYearsWithoutException( cpnDayCount, capletStartDate, capletEndDate );
		itsEvalTime		    = GetAndCheckEvalTime( evalDate, mod, ARM_GP_ImpliedVol::itsFuncName );
		itsFwdStartTime	    = mod->GetTimeFromDate( fwdStartDate );
		itsFwdEndTime	    = mod->GetTimeFromDate( fwdEndDate );
		int indexDayCount   = GetDayCount( arg[11], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		itsFwdPeriod        = CountYearsWithoutException( indexDayCount, fwdStartDate, fwdEndDate );
		itsPayTime		    = mod->GetTimeFromDate( capletPayDate );
		itsResetTime	    = mod->GetTimeFromDate( capletResetDate);

		/// some validation
		CheckNbSmaller( itsEvalTime, itsResetTime,      "EvalTime",     "ResetTime", ARM_GP_ImpliedVol::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTime,itsFwdEndTime,  "FwdStartTime", "FwdEndTime",	ARM_GP_ImpliedVol::itsFuncName,__LINE__,__FILE__  );
		CheckPositiveNb( itsPeriod,     "Cpn Period", ARM_GP_ImpliedVol::itsFuncName,__LINE__,__FILE__  );
		CheckPositiveNb( itsFwdPeriod,  "Index Period", ARM_GP_ImpliedVol::itsFuncName,__LINE__,__FILE__  );

		/// delete char* for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}




////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Swap Rate function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_ImpliedVol::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );

	
	/// Conversion to vectorial strikes
	if( ARM_VectorPtr(NULL) == itsStrikeVector)
	{
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		itsStrikeVector		= ARM_VectorPtr( new ARM_GP_Vector( statesSize, itsStrikeDouble ) );
	}

	/// call function
	ARM_VectorPtr CapFloorLet =  itsModelIR->ImpliedVol( itsCurveName, itsEvalTime, itsPayTime, itsPeriod,
		itsNotional, itsResetTime, itsFwdStartTime, itsFwdEndTime, itsFwdPeriod, *itsStrikeVector, itsCapFloor, states );

	/// return result
	return ARM_GramFctorArg(CapFloorLet);
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_CapDigital_Common
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_ImpliedVol::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	
	
	/// return result
	vector<double> tmpResult;
	tmpResult.reserve(2);
	
	tmpResult.push_back( itsResetTime );
	if( itsResetTime != itsFwdStartTime )
		tmpResult.push_back(  itsFwdStartTime );
	
	tmpResult.push_back( itsFwdEndTime );
    if( itsPayTime != itsFwdEndTime )
		tmpResult.push_back( itsPayTime );

	/// return result
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( tmpResult.size(), tmpResult.begin() )), ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///////////		Cap			      //////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_GP_Cap
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Cap::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 13, ARM_GP_Cap::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_Cap::itsFuncName );
	itsCurveName = arg[0].GetString();
	ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(nodes, arg, 2, mod, itsCurveName );
	bool isDefault;
	isDefault=SetFrequency(   nodes, arg, 5, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	isDefault=SetDayCount(	 nodes, arg, 6, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	SetResetDaysGap( nodes, arg, 7, mod, itsCurveName );																	
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_CAp
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Cap::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,	ARM_GP_Cap::itsFuncName );

	/// handle vectorial strike
	if( arg[3].GetType() ==  GFAT_DOUBLE_TYPE )
	{
		itsStrikeDouble = arg[3].GetDouble();
		itsStrikeVector	= ARM_VectorPtr(NULL);
	}
	else if( arg[3].GetType() ==  GFAT_VECTOR_TYPE )
	{
		const double UNASSIGNED = -1111;
		itsStrikeDouble =  UNASSIGNED;
		itsStrikeVector	= arg[3].GetVector();
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );
	
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 13, ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[7], GFAT_DOUBLE_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[9], GFAT_DOUBLE_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[10],GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_MATURITY_TYPE,	ARM_GP_Cap::itsFuncName );
		GPAF_CheckArgType( arg[12],GFAT_STRING_TYPE,	ARM_GP_Cap::itsFuncName );


		/// grab inputs
		itsCurveName			    = arg[0].GetString();
		ARM_Date capStartDate       = arg[1].GetDate(); 
		string capFloorString	    = arg[4].GetString();
		double payDaysGap		    = arg[8].GetDouble();											
		string advArrString	        = arg[10].GetString();

																			
		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
		if(!itsModelIR)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for Cap keyword");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );

		
		/// the logic here is to get the default calendar for the float leg
		/// these are from the default index
		/// the fix calendar is the currency itself!
		char* payCalendar	= ccy->GetPayCalName(defaultIndex);
		char* resetCalendar= ccy->GetResetCalName(defaultIndex);

		
		ARM_Date capEndDate = ComputeEndDate( capStartDate, arg[2], payCalendar );
		
		/// Floating leg data extraction
		/// get the floating leg convention
		int resetTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
		int resetDaysGap    = GetResetDaysGap(arg[7], mod, itsCurveName );
		int cpnFreq	    = GetFrequency(arg[5], mod, itsCurveName, ARM_GP_Cap::itsFuncName, GetFloatFrequencyFtor(ccy));
		int cpnDayCount	= GetDayCount( arg[6], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		int indexDayCount   = GetDayCount( arg[12], mod, itsCurveName, GetFloatDayCountFtor(ccy));

        int stdSpotGap = ccy->GetSpotDays();

		ARM_DateStrip cpnDateStrip( capStartDate, capEndDate, cpnFreq, cpnDayCount, resetCalendar,
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, resetDaysGap, cpnFreq, payDaysGap,
			payCalendar,resetTiming,K_ARREARS,true,GETDEFAULTVALUESTR,stdSpotGap);

		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pFwdStartTimes      = cpnDateStrip.GetFwdRateStartDates();
		ARM_GP_Vector* pFwdEndTimes        = cpnDateStrip.GetFwdRateEndDates();
		ARM_GP_Vector* pPayTimes           = cpnDateStrip.GetPaymentDates();
		ARM_GP_Vector* pPeriods            = cpnDateStrip.GetInterestTerms();
		ARM_GP_Vector* pFwdResetTimes      = cpnDateStrip.GetResetDates();

        /// Compute correctly fwd end dates w.r.t. index term
        string indexTerm( ComputeIndexTerm(arg[11],arg[1],arg[2],ARM_GP_Cap::itsFuncName) );
        int i;
		for(i=0; i<pFwdStartTimes->size(); ++i )
        {
            ARM_Date fwdEndDate( ComputeEndDate( (*pFwdStartTimes)[i], indexTerm, resetCalendar ) );
            (*pFwdEndTimes)[i] = fwdEndDate.GetJulian();
        }
		

		/// copy constructor
		itsCapFloor		    = ARM_ArgConv_CapFloor.GetNumber( capFloorString );
		itsNotional		    = arg[9].GetDouble();
		itsEvalTime	        = GetAndCheckEvalTime( evalDate, mod, ARM_GP_Cap::itsFuncName );
		itsPayTimes	        = ARM_GP_Vector( *pPayTimes );
		itsPeriods          = ARM_GP_Vector( *pPeriods );
		itsFwdResetTimes    = ARM_GP_Vector(*pFwdResetTimes);
		itsFwdStartTimes	= ARM_GP_Vector( *pFwdStartTimes );
		itsFwdEndTimes      = ARM_GP_Vector( *pFwdEndTimes );
		itsFwdPeriods       = ARM_GP_Vector(pFwdStartTimes->size());										

		/// get time from date! 
		for(i=0; i<pFwdStartTimes->size(); ++i )
		{
			itsFwdPeriods[i] = CountYearsWithoutException( indexDayCount, ARM_Date(itsFwdStartTimes[i]), ARM_Date(itsFwdEndTimes[i]) );
			itsFwdStartTimes[i] = mod->GetTimeFromDate( itsFwdStartTimes[i] );
			itsFwdEndTimes[i]   = mod->GetTimeFromDate( itsFwdEndTimes[i] );
			itsPayTimes[i] = mod->GetTimeFromDate( itsPayTimes[i] );
			itsFwdResetTimes[i] = mod->GetTimeFromDate( itsFwdResetTimes[i] );
		}

		
		/// Validation
		CheckNbSmaller( itsEvalTime, itsFwdResetTimes[0], "EvalTime","1st FwdResetTime", ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );																
		CheckVectorPositiveNb( itsFwdPeriods,"FwdPeriods",ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsPeriods,  "interest Periods",ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFwdResetTimes,	"FwdResetTimes",   ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdStartTimes,	"FwdStartTimes",   ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes,	    "FwdEndTimes",	   ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsPayTimes,	        "PayTimes",	       ARM_GP_Cap::itsFuncName,__LINE__,__FILE__  );
	
		/// delete char* for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Cap
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Cap Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Cap::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	
	/// Conversion to vectorial strikes
	if( ARM_VectorPtr(NULL) == itsStrikeVector)
	{
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		itsStrikeVector		= ARM_VectorPtr( new ARM_GP_Vector( statesSize, itsStrikeDouble ) );
	}

	ARM_VectorPtr Cap =  itsModelIR->VanillaCap( itsCurveName,
									itsEvalTime,
									itsPayTimes,
									itsPeriods,
									itsNotional,
									itsFwdResetTimes,
									itsFwdStartTimes,
									itsFwdEndTimes,
									itsFwdPeriods,
									*itsStrikeVector,
									itsCapFloor,    
									states );

	return ARM_GramFctorArg(Cap);	
}





////////////////////////////////////////////////////
///	Class  : ARM_GP_Cap
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Cap::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	CC_STL_VECTOR( double ) result;
	result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
	result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
	result.insert( result.end(), itsPayTimes.begin(), itsPayTimes.end() );

	/// sort and remove duplicates!
	CC_NS( std, sort )( result.begin(), result.end() );
	CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
	result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///////////		Sum Option		  //////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_SumOption
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_SumOption::SetDefaults(ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 14, ARM_GP_SumOption::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_SumOption::itsFuncName );
	itsCurveName = arg[0].GetString();
	ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(nodes, arg, 2, mod, itsCurveName );
	bool isDefault;
	isDefault=SetFrequency(   nodes, arg, 6, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	SetResetDaysGap( nodes, arg, 7, mod, itsCurveName );
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SumOption
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_SumOption::GrabInputs(ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// this is always reasserted!
	GPAF_CheckArgType( arg[4], GFAT_VECTOR_TYPE,	ARM_GP_SumOption::itsFuncName );

	/// handle vectorial strike
	if( arg[4].GetType() ==  GFAT_DOUBLE_TYPE )
	{
		itsStrikeDouble = arg[4].GetDouble();
		itsStrikeVector	= ARM_VectorPtr(NULL);
	}
	else if( arg[4].GetType() ==  GFAT_VECTOR_TYPE )
	{
		const double UNASSIGNED = -1111;
		itsStrikeDouble =  UNASSIGNED;
		itsStrikeVector	= arg[4].GetVector();
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be vector or double" );
	
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 14, ARM_GP_SumOption::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_SumOption::itsFuncName ); // Curve Name
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_SumOption::itsFuncName ); // Start Date
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_SumOption::itsFuncName ); // End Date
		GPAF_CheckArgType( arg[3], GFAT_DATE_TYPE,ARM_GP_SumOption::itsFuncName );	   // Pay Date
		GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE,	ARM_GP_SumOption::itsFuncName ); // Cap Floor
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,	ARM_GP_SumOption::itsFuncName ); // Freq
		GPAF_CheckArgType( arg[7], GFAT_DOUBLE_TYPE,	ARM_GP_SumOption::itsFuncName ); // Reset Gap
		GPAF_CheckArgType( arg[8],GFAT_STRING_TYPE,	ARM_GP_SumOption::itsFuncName );	   // Reset Timing
		GPAF_CheckArgType( arg[9],GFAT_MATURITY_TYPE,	ARM_GP_SumOption::itsFuncName ); // Index Term
		GPAF_CheckArgType( arg[10],GFAT_STRING_TYPE,	ARM_GP_SumOption::itsFuncName ); // Index Day Count
		GPAF_CheckArgType( arg[12],GFAT_DOUBLE_TYPE,	ARM_GP_SumOption::itsFuncName ); // First idx
		GPAF_CheckArgType( arg[13],GFAT_DOUBLE_TYPE,	ARM_GP_SumOption::itsFuncName ); // Last idx


		/// grab inputs
		itsCurveName			    = arg[0].GetString();
		ARM_Date startDate       = arg[1].GetDate();
		string capFloorString	    = arg[5].GetString();
		string advArrString	        = arg[8].GetString();
																			
		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
		if(!itsModelIR)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for Cap keyword");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );

		
		/// the logic here is to get the default calendar for the float leg
		/// these are from the default index
		/// the fix calendar is the currency itself!
		char* payCalendar	= ccy->GetPayCalName(defaultIndex);
		char* resetCalendar= ccy->GetResetCalName(defaultIndex);
		
		ARM_Date endDate = ComputeEndDate( startDate, arg[2], payCalendar );
		
		/// Floating leg data extraction
		/// get the floating leg convention
		int resetTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
		int resetDaysGap    = GetResetDaysGap(arg[7], mod, itsCurveName );
		int cpnFreq	    = GetFrequency(arg[6], mod, itsCurveName, ARM_GP_SumOption::itsFuncName, GetFloatFrequencyFtor(ccy));
		int indexDayCount   = GetDayCount( arg[10], mod, itsCurveName, GetFloatDayCountFtor(ccy));

        int stdSpotGap = ccy->GetSpotDays();

		ARM_DateStrip cpnDateStrip( 
			startDate, 
			endDate, 
			cpnFreq, 
			K30_360, 
			resetCalendar,
			K_MOD_FOLLOWING, 
			K_ADJUSTED, 
			K_SHORTSTART, 
			resetDaysGap,
			cpnFreq, 
			0,
			payCalendar,
			resetTiming,
			K_ARREARS,
			true,
			GETDEFAULTVALUESTR,
			stdSpotGap);

		ARM_GP_Vector* tmpFwdStartDates = cpnDateStrip.GetFwdRateStartDates();
		ARM_GP_Vector* tmpFwdEndDates = cpnDateStrip.GetFwdRateEndDates();
		ARM_GP_Vector* tmpFwdResetDates = cpnDateStrip.GetResetDates();

		size_t nbPeriods = tmpFwdStartDates->size();

		int firstIdx, lastIdx;

		firstIdx = arg[12].GetDouble();
		if (firstIdx == -1)
			firstIdx = 0;
		lastIdx = arg[13].GetDouble();
		if (lastIdx == -1)
			lastIdx = nbPeriods-1;

		if (firstIdx < 0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : First idx should be greater than 0");

		if (firstIdx > lastIdx )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Last idx should be greater or equal than first idx.");

		if (lastIdx >= nbPeriods )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Last idx should be less than the number of periods of the sum option.");

		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pFwdStartTimes      = new ARM_GP_Vector(lastIdx-firstIdx+1);
		ARM_GP_Vector* pFwdEndTimes        = new ARM_GP_Vector(lastIdx-firstIdx+1);
		itsPayTime = arg[3].GetDate().GetJulian();
		ARM_GP_Vector* pFwdResetTimes      =new ARM_GP_Vector(lastIdx-firstIdx+1);

		int i;
		if(arg[9].GetString() != "")
		{
			string indexTerm( ComputeIndexTerm(arg[9],arg[1],arg[2],ARM_GP_SumOption::itsFuncName) );
			for(i=firstIdx; i<=lastIdx; ++i )
			{
				ARM_Date fwdEndDate( ComputeEndDate( (*tmpFwdStartDates)[i], indexTerm, resetCalendar ) );
				(*pFwdEndTimes)[i-firstIdx] = fwdEndDate.GetJulian();
			}
		}
		else
		{
			for(i=firstIdx; i<=lastIdx; ++i )
			{
				(*pFwdEndTimes)[i-firstIdx] = (*tmpFwdEndDates)[i];
			}
		}
		

		/// copy constructor
		itsCapFloor		    = ARM_ArgConv_CapFloor.GetNumber( capFloorString );
		itsEvalTime	        = GetAndCheckEvalTime( evalDate, mod, ARM_GP_SumOption::itsFuncName );
		itsFwdResetTimes    = ARM_GP_Vector( *pFwdResetTimes );
		itsFwdStartTimes	= ARM_GP_Vector( *pFwdStartTimes );
		itsFwdEndTimes      = ARM_GP_Vector( *pFwdEndTimes );
		itsFwdPeriods       = ARM_GP_Vector(lastIdx-firstIdx+1);

		delete pFwdResetTimes;
		delete pFwdStartTimes;
		delete pFwdEndTimes;

		/// get time from date! 
		for(i=firstIdx; i<=lastIdx; ++i )
		{
			itsFwdPeriods[i-firstIdx] = CountYearsWithoutException( indexDayCount, ARM_Date((*tmpFwdStartDates)[i]), ARM_Date(itsFwdEndTimes[i-firstIdx]) );
			itsFwdStartTimes[i-firstIdx] = mod->GetTimeFromDate( (*tmpFwdStartDates)[i] );
			itsFwdEndTimes[i-firstIdx]   = mod->GetTimeFromDate( itsFwdEndTimes[i-firstIdx] );
			itsFwdResetTimes[i-firstIdx] = mod->GetTimeFromDate( (*tmpFwdResetDates)[i] );
		}

		itsPayTime = mod->GetTimeFromDate( itsPayTime );

		if( arg[11].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double coeff = arg[11].GetDouble();
			itsCoeffVector	= ARM_GP_VectorPtr(new ARM_GP_Vector(lastIdx-firstIdx+1,coeff));
		}
		else if( arg[11].GetType() ==  GFAT_VECTOR_TYPE )
		{
			ARM_GP_Vector* pCoeff      = new ARM_GP_Vector(lastIdx-firstIdx+1);
			ARM_GP_VectorPtr tmpVector	= arg[11].GetVector();
			for(i=firstIdx; i<=lastIdx; ++i )
			{
				(*pCoeff)[i-firstIdx] = (*tmpVector)[i];
			}
			itsCoeffVector = ARM_GP_VectorPtr(pCoeff);
		}
		else if( arg[11].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr coeffCurve	= arg[11].GetCurve();
			itsCoeffVector = ARM_GP_VectorPtr(new ARM_GP_Vector(lastIdx-firstIdx+1,0.0));
			for(i=firstIdx; i<=lastIdx; ++i )
				(*itsCoeffVector)[i-firstIdx] = coeffCurve->Interpolate((*tmpFwdStartDates)[i]);
		}

		/// Validation
		CheckNbSmaller( itsEvalTime, itsFwdResetTimes[0], "EvalTime","1st FwdResetTime", ARM_GP_SumOption::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFwdPeriods,"FwdPeriods",ARM_GP_SumOption::itsFuncName,__LINE__,__FILE__  );		
		CheckVectorIncreasing( itsFwdResetTimes,	"FwdResetTimes",   ARM_GP_SumOption::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdStartTimes,	"FwdStartTimes",   ARM_GP_SumOption::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes,	    "FwdEndTimes",	   ARM_GP_SumOption::itsFuncName,__LINE__,__FILE__  );
	
		/// delete char* for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SumOption
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the SunOption Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_SumOption::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );
	
	/// Conversion to vectorial strikes
	if( ARM_VectorPtr(NULL) == itsStrikeVector)
	{
		size_t statesSize	= states != ARM_PricingStatesPtr( NULL )? states->size() : 1;
		itsStrikeVector		= ARM_VectorPtr( new ARM_GP_Vector( statesSize, itsStrikeDouble ) );
	}

	ARM_VectorPtr SumOption =  itsModelIR->VanillaSumOption( 
									itsCurveName,
									itsEvalTime,
									itsCapFloor,
									*itsCoeffVector,
									itsFwdResetTimes,
									itsFwdStartTimes,
									itsFwdEndTimes,
									itsPayTime,
									itsFwdPeriods,
									*itsStrikeVector,
									1.0,
									NULL,
									NULL,
									states );

	return ARM_GramFctorArg(SumOption);	
}





////////////////////////////////////////////////////
///	Class  : ARM_GP_SumOption
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_SumOption::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	CC_STL_VECTOR( double ) result;
	result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
	result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
	result.push_back( itsPayTime );

	/// sort and remove duplicates!
	CC_NS( std, sort )( result.begin(), result.end() );
	CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
	result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}


//////////////////////////////////////
//// bond analytics
//////////////////////////////////////

//////////////////////////////////////
/// Price To Yield and Yield To Price operator
//////////////////////////////////////


ARM_GP_PTYAndYTPFctor::ARM_GP_PTYAndYTPFctor( 
	const ARM_BondAnalytic_PricingFunc& pricingFunc, 
	const string& funcName )
:	itsPricingFunc( pricingFunc ),
	itsFuncName( funcName ),
	itsModelIR(NULL)
{}


////////////////////////////////////////////////////
///	Class  : ARM_GP_PTYAndYTPFctor
///	Routine: operator() (Price To Yield)
///	Returns: ARM_GramFctorArg
///	Action : computes either 
///				-from the yield the price of a bond
///					with given coupon and redemption value
///					using the given yield!
///				-or from the price the yield of a bond
///					with given coupon and redemption value
///					using a given price! Need to solve the non
///					linear equation!
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_PTYAndYTPFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector<ARM_ExpNodePtr>& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 8, itsFuncName );

	GPAF_CheckArgType(  arg[0], GFAT_DATE_TYPE,			itsFuncName );
	GPAF_CheckArgType(  arg[1], GFAT_DATE_TYPE,			itsFuncName );
	GPAF_CheckArgType(  arg[2], GFAT_DATEORMATU_TYPE,	itsFuncName );
	GPAF_CheckArgType(  arg[3], GFAT_VECTOR_TYPE,		itsFuncName );
	GPAF_CheckArgType(  arg[4], GFAT_VECTOR_TYPE,		itsFuncName );
	GPAF_CheckArgType(  arg[5], GFAT_VECTOR_TYPE,		itsFuncName );
	GPAF_CheckArgType(  arg[6], GFAT_STRING_TYPE,		itsFuncName );
	GPAF_CheckArgType(  arg[7], GFAT_STRING_TYPE,		itsFuncName );

	/// grab the arguments
	ARM_Date settlementDate			= arg[0].GetDate();
	ARM_Date firstCouponDate		= arg[1].GetDate();

	/// compute inputs
	/// currency and calendar are in ARM world similar ...
	itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
    if(!itsModelIR)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for PRICETORATE/RATETOPRICE keywords");

	ARM_Currency* ccy		= mod->GetCurrency( "" );
	char* ccyName			= ccy->GetCcyName();

	char fixCalendar[100];
    ccy->CalcFixPayCal(fixCalendar);

	/// the second argument can either be a date or a maturity
	ARM_Date endDate	    = ComputeEndDate( firstCouponDate, arg[2], fixCalendar );

    int frequency = GetFrequency(arg[6], mod, "", ARM_GP_PTYAndYTPFctor::itsFuncName, GetFixedFrequencyFtor(ccy));
	int dayCount = GetDayCount( arg[7], mod, "", GetFixedDayCountFtor(ccy));

	/// year fraction computation
    ARM_Date issueDate(firstCouponDate);
	issueDate.AddPeriod(-frequency);

	if (settlementDate < issueDate)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": settlement date is before the issue date." );

	if (settlementDate > endDate)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": settlement date is after the maturity date." );
	

	ARM_Date previousCouponDate, nextCouponDate;

	if (settlementDate < firstCouponDate)
	{
		previousCouponDate = issueDate;
		nextCouponDate = firstCouponDate;
	}
	else
	{
		PreviousRegularDate(frequency,issueDate,settlementDate,previousCouponDate);
		nextCouponDate = previousCouponDate;
		nextCouponDate.AddPeriod(frequency);
	}

	/// validation


	/// daycount computation
	double settlToNext		= CountYears(dayCount, settlementDate, nextCouponDate );
	double settlToMaturity	= CountYears(dayCount, settlementDate, endDate );
	double accruingTime		= CountYears(dayCount, previousCouponDate, settlementDate ); 

	/// vectorial formula?
	size_t position[]={3,4,5};
	size_t* positionBegin	= position;
	size_t* positionEnd		= position+3;
	int firstVectorIndex = ARM_GramFctorArgVectorHelper::FindFirstVectorIndex( arg, positionBegin, positionEnd );
	if(	firstVectorIndex != ARM_GramFctorArgVectorHelper::IndexNotFound )
	{
		vector<ARM_VectorPtr> argVec(3);
		size_t refSize = arg[firstVectorIndex].GetVector()->size();
		for( size_t* iter=positionBegin; iter<positionEnd; ++iter )
		{
			if(arg[*iter].GetType() == GFAT_DOUBLE_TYPE )
				argVec[iter-positionBegin] = ARM_VectorPtr( new ARM_GP_Vector(refSize,arg[*iter].GetDouble() ) );
			else
			{
#if defined(__GP_STRICT_VALIDATION)
				if( arg[*iter].GetVector()->size() != refSize )
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": vectorial formula with arguments of different size!");
#endif
				argVec[iter-positionBegin] = arg[*iter].GetVector();
			}
		}

		ARM_VectorPtr result = arg[firstVectorIndex].GetVector();
		/// use in place result and store everything in the vector with firstVectorIndex
		for( size_t i=0;i<refSize; ++i )
		{
			(*result)[i]= (*itsPricingFunc)( (*argVec[0])[i], (*argVec[1])[i], (*argVec[2])[i],
							frequency, settlToNext, settlToMaturity, accruingTime );
		}
		return ARM_GramFctorArg( result );
	}
	else
	{
		double coupon			= arg[3].GetDouble();
		double redemptionValue	= arg[4].GetDouble();
		double yieldOrPrice		= arg[5].GetDouble();

		double result = (*itsPricingFunc)(coupon, redemptionValue, yieldOrPrice,
							frequency, settlToNext, settlToMaturity, accruingTime );

		return ARM_GramFctorArg( result );
	}
}



////////////////////////////////////////////////////
///////////		Spread Option     //////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpreadOption
///	Routine: Desctructor
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
ARM_GP_SpreadOption::~ARM_GP_SpreadOption()
{
	size_t size = itsSwapLongFixPayTimes.size();

	// test
	if (	size!= itsSwapLongFixPayPeriods.size()
		||	size!= itsSwapShortFixPayTimes.size()
		||	size!= itsSwapShortFixPayPeriods.size()  )
		ARM_THROW( ERR_INVALID_ARGUMENT, "unresolved problem in ARM_GP_SpreadOption destructor");


	DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayPeriods);
	DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayPeriods);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpreadOption
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_SpreadOption::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 19, ARM_GP_SpreadOption::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_SpreadOption::itsFuncName );
	itsCurveName = arg[0].GetString();
		ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(nodes, arg, 2, mod, itsCurveName );
	bool isDefault;
	isDefault=SetFrequency(   nodes, arg, 9, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	isDefault=SetDayCount(	 nodes, arg, 10, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	SetResetDaysGap( nodes, arg, 11, mod, itsCurveName );																	
	
	/// set payment frequency (default = cpn frequency)
	GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE, ARM_GP_SpreadOption::itsFuncName );
	SetFrequencyToGramFctorArg(nodes, (ARM_GramFctorArgVector&)arg, 17, mod, itsCurveName, arg[9].GetString());
	
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpreadOption
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_SpreadOption::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	GPAF_CheckArgType( arg[5], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_SpreadOption::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_SpreadOption::itsFuncName );
	GPAF_CheckArgType( arg[7], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_SpreadOption::itsFuncName );
	GPAF_CheckArgType( arg[13], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_SpreadOption::itsFuncName );


	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 19, ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,				ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,		ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[3], GFAT_MATURITY_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[4], GFAT_MATURITY_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		
		GPAF_CheckArgType( arg[8], GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[10],GFAT_MULTITOKENSTRING_TYPE,	ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[11],GFAT_DOUBLE_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[12],GFAT_DOUBLE_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[14],GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[15],GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[16],GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[17],GFAT_STRING_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		GPAF_CheckArgType( arg[18],GFAT_DOUBLE_TYPE,			ARM_GP_SpreadOption::itsFuncName );
		
		/// grab inputs
		itsCurveName				= arg[0].GetString();
		ARM_Date soStartDate		= arg[1].GetDate(); 
		string capFloorString		= arg[8].GetString();
		double payDaysGap			= arg[12].GetDouble();											
		string advArrString			= arg[14].GetString();

		int longIndexType		 = ARM_ArgConv_IndexClass.GetNumber(arg[15].GetString());
		int shortIndexType		 = ARM_ArgConv_IndexClass.GetNumber(arg[16].GetString());

																		
		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
		if(!itsModelIR)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for Cap keyword");

		ARM_Currency* ccy		= mod->GetCurrency( itsCurveName );
		char* ccyName			= ccy->GetCcyName();
		ARM_INDEX_TYPE defaultIndex = GetDefaultIndexFromCurrency( ccyName );
		
		
		/// the logic here is to get the default calendar for the float leg
		/// these are from the default index
		/// the fix calendar is the currency itself!
		char* payCalendar	= ccy->GetPayCalName(defaultIndex);
		char* resetCalendar= ccy->GetResetCalName(defaultIndex);

		ARM_Date soEndDate = ComputeEndDate( soStartDate, arg[2], NULL/*payCalendar*/);
		
		/// Floating leg data extraction
		/// get the floating leg convention
		int resetTiming		= ARM_ArgConv_Timing.GetNumber( advArrString );
		int resetDaysGap    = GetResetDaysGap(arg[11], mod, itsCurveName );
		int cpnFreq	    = GetFrequency(arg[9],  mod, itsCurveName, ARM_GP_SpreadOption::itsFuncName, GetFloatFrequencyFtor(ccy));
		int cpnDayCount	= GetDayCount (arg[10], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		int payFreq	    = GetFrequency(arg[17], mod, itsCurveName, ARM_GP_SpreadOption::itsFuncName, GetFloatFrequencyFtor(ccy));
		
        int stdSpotGap = ccy->GetSpotDays();

		// un peu au pif: à checker...
		ARM_DateStrip cpnDateStrip( soStartDate, soEndDate, cpnFreq, cpnDayCount, resetCalendar,
			K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, resetDaysGap, payFreq, payDaysGap,
			payCalendar,resetTiming,K_ARREARS,true,GETDEFAULTVALUESTR,stdSpotGap);


		/// not cloned hence no need to delete this!
		ARM_GP_Vector* pPayTimes        = cpnDateStrip.GetPaymentDates();
		ARM_GP_Vector* pPayPeriods      = cpnDateStrip.GetInterestTerms();
		ARM_GP_Vector* pResetTimes      = cpnDateStrip.GetResetDates();

		/// copy constructor
		itsEvalTime	        = GetAndCheckEvalTime( evalDate, mod, ARM_GP_SpreadOption::itsFuncName );
		itsCallPut		    = ARM_ArgConv_CapFloor.GetNumber( capFloorString );
		itsResetTimes		= ARM_GP_Vector( *pResetTimes );
		itsPayTimes			= ARM_GP_Vector( *pPayTimes );
		itsPayPeriods       = ARM_GP_Vector( *pPayPeriods );
		
		//
		/// Create underlying CMS schedules
		//

		DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSwapLongFixPayPeriods);
		DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSwapShortFixPayPeriods);

		itsSwapLongFloatStartTime.resize(itsResetTimes.size());
		itsSwapLongFloatEndTime.resize(itsResetTimes.size());
		itsSwapLongFixPayTimes.resize(itsResetTimes.size());
		itsSwapLongFixPayPeriods.resize(itsResetTimes.size());
		
		itsSwapShortFloatStartTime.resize(itsResetTimes.size());
		itsSwapShortFloatEndTime.resize(itsResetTimes.size());
		itsSwapShortFixPayTimes.resize(itsResetTimes.size());
		itsSwapShortFixPayPeriods.resize(itsResetTimes.size());
	
		/// in case of CMS index:
		char fixCalendar[100];
		ccy->CalcFixPayCal(fixCalendar);
		int  fixFreq	 = ccy->GetFixedPayFreq();
		int  fixDayCount = ccy->GetFixedDayCount();
		
		/// in case of LIBOR index:
		char floatCalendar[100];
		ccy->CalcFloatPayCal(floatCalendar);
		int  floatDayCount = ccy->GetLiborIndexDayCount();
		
		ARM_Date endDate;

		for (size_t i(0); i<itsResetTimes.size(); i++)
		{			
			ARM_Date startDate(itsResetTimes[i]);
			startDate.GapBusinessDay(-resetDaysGap, resetCalendar);

			if (longIndexType == K_CMS)
			{
				// CMS #1 (Long)
				endDate = ComputeEndDate(startDate, arg[3], NULL/*fixCalendar*/);

				ARM_DateStrip longDateStrip(startDate, endDate, fixFreq, fixDayCount, fixCalendar,
											K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
											fixCalendar );

				itsSwapLongFixPayTimes[i]    = new ARM_GP_Vector ( *longDateStrip.GetPaymentDates() );
				itsSwapLongFixPayPeriods[i]	 = new ARM_GP_Vector ( *longDateStrip.GetInterestTerms() );
				itsSwapLongFloatStartTime[i] = (*longDateStrip.GetFlowStartDates())[0];
				itsSwapLongFloatEndTime[i]	 = (*longDateStrip.GetFlowEndDates())[longDateStrip.GetFlowEndDates()->size()-1];
			}
			else if (longIndexType == K_LIBOR)
			{
				endDate = ComputeEndDate(startDate, arg[3], floatCalendar);
				double payPeriod = CountYears(floatDayCount, startDate.GetJulian(), endDate.GetJulian());
				itsSwapLongFixPayTimes[i]    = new ARM_GP_Vector (1, endDate.GetJulian());
				itsSwapLongFixPayPeriods[i]  = new ARM_GP_Vector (1, payPeriod);
				itsSwapLongFloatStartTime[i] = startDate.GetJulian();
				itsSwapLongFloatEndTime[i]   = endDate.GetJulian();
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" SpreadOption keyword : bad index type for index #1");
			

			// CMS #2 (Short)
			if (shortIndexType == K_CMS)
			{
				endDate = ComputeEndDate(startDate, arg[4], NULL/*fixCalendar*/);

				ARM_DateStrip shortDateStrip(startDate, endDate, fixFreq, fixDayCount, fixCalendar,
											K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
											fixCalendar );

				itsSwapShortFixPayTimes[i]    = new ARM_GP_Vector ( *shortDateStrip.GetPaymentDates() );
				itsSwapShortFixPayPeriods[i]  = new ARM_GP_Vector ( *shortDateStrip.GetInterestTerms() );
				itsSwapShortFloatStartTime[i] = (*shortDateStrip.GetFlowStartDates())[0];
				itsSwapShortFloatEndTime[i]	  = (*shortDateStrip.GetFlowEndDates())[shortDateStrip.GetFlowEndDates()->size()-1];
			}
			else if (shortIndexType == K_LIBOR)
			{
				endDate = ComputeEndDate(startDate, arg[4], floatCalendar);
				double payPeriod = CountYears(floatDayCount, startDate.GetJulian(), endDate.GetJulian());
				itsSwapShortFixPayTimes[i]    = new ARM_GP_Vector (1, endDate.GetJulian());
				itsSwapShortFixPayPeriods[i]  = new ARM_GP_Vector (1, payPeriod);
				itsSwapShortFloatStartTime[i] = startDate.GetJulian();
				itsSwapShortFloatEndTime[i]   = endDate.GetJulian();
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" SpreadOption keyword : bad index type for index #2");
			
		}
		
		/// get time from date! 
		for(i=0; i<itsPayTimes.size(); ++i )
		{			
			itsPayTimes[i]	 = mod->GetTimeFromDate( itsPayTimes[i] );
			itsResetTimes[i] = mod->GetTimeFromDate( itsResetTimes[i] );
			
			itsSwapLongFloatStartTime[i] = mod->GetTimeFromDate( itsSwapLongFloatStartTime[i] );
			itsSwapLongFloatEndTime[i]   = mod->GetTimeFromDate( itsSwapLongFloatEndTime[i] );
			
			
			itsSwapShortFloatStartTime[i] = mod->GetTimeFromDate( itsSwapShortFloatStartTime[i] );
			itsSwapShortFloatEndTime[i]   = mod->GetTimeFromDate( itsSwapShortFloatEndTime[i] );

			for (size_t j(0); j<itsSwapLongFixPayTimes[i]->size(); j++)
				(*itsSwapLongFixPayTimes[i])[j]  = mod->GetTimeFromDate( (*itsSwapLongFixPayTimes[i])[j] );

			for (j=0; j<itsSwapShortFixPayTimes[i]->size(); j++)
				(*itsSwapShortFixPayTimes[i])[j]  = mod->GetTimeFromDate( (*itsSwapShortFixPayTimes[i])[j] );
			
		}
		
		/// Coeff for CMS #1	
		/// handle vectorial coeffs
		size_t sizeFlows = itsPayTimes.size();    
		
		if( arg[5].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double coeff = arg[5].GetDouble();
			itsCoeffLong = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,coeff));
		}
		else if( arg[5].GetType() ==  GFAT_VECTOR_TYPE )
		{			
			itsCoeffLong= arg[5].GetVector();
			if (itsCoeffLong->size() != sizeFlows)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "SpreadOption keyword: invalid vector size" );

		}
		else if( arg[5].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr coeffCurve	= arg[5].GetCurve();
			itsCoeffLong = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,0.0));
			for(i=0; i<sizeFlows; ++i)			
				(*itsCoeffLong)[i] = coeffCurve->Interpolate(itsResetTimes[i]);										
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "coeff should be double, vector or ARM_Curve" );

		/// Coeff for CMS #2	
		/// handle vectorial coeffs
		
		if( arg[6].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double coeff = arg[6].GetDouble();
			itsCoeffShort = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,coeff));
		}
		else if( arg[6].GetType() ==  GFAT_VECTOR_TYPE )
		{
			itsCoeffShort = arg[6].GetVector();
			if (itsCoeffShort->size() != sizeFlows)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "SpreadOption keyword: invalid vector size" );
		}
		else if( arg[6].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr coeffCurve	= arg[6].GetCurve();
			itsCoeffShort = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,0.0));
			for(i=0; i<sizeFlows; ++i)			
				(*itsCoeffShort)[i] = coeffCurve->Interpolate(itsResetTimes[i]);										
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "coeff should be double, vector or ARM_Curve" );


		/// Spread option strike
		/// handle vectorial strikes
				
		if( arg[7].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double strike = arg[7].GetDouble();
			itsStrikes = ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizeFlows,strike));
		}
		else if( arg[7].GetType() ==  GFAT_VECTOR_TYPE )
		{
			/// Step-up & state dependent strikes not supported
			ARM_VectorPtr inputStrikes(arg[7].GetVector());
			size_t inputSize = inputStrikes->size();
			if(inputSize == sizeFlows)
			{
				/// Vector is assumed to describe step-up deterministic strikes
				itsStrikes = ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizeFlows,(*inputStrikes)));
			}
			else
			{
				/// Vector is assumed to described state dependent strikes but identical for all flows
				itsStrikes = ARM_GP_MatrixPtr(new ARM_GP_Matrix(inputSize,sizeFlows));
				double strike;
				for(size_t j=0;j<inputSize;++j)
				{
					strike = (*inputStrikes)[j];
					for(size_t k=0;k<sizeFlows;++k)
						(*itsStrikes)(j,k) = strike;
				}
			}
		}
		else if( arg[7].GetType() ==  GFAT_CURVE_TYPE )
		{
			/// Curve is interpolated to describe step-up deterministic strikes
			ARM_GP_CurvePtr strikeCurve	= arg[7].GetCurve();
			itsStrikes = ARM_GP_MatrixPtr(new ARM_GP_Matrix(1,sizeFlows,0.0));
			for(i=0; i<sizeFlows; ++i)			
				(*itsStrikes)(0,i) = strikeCurve->Interpolate(itsResetTimes[i]);										
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "strike should be double, vector or ARM_Curve" );

		/// Notional
		/// handle vectorial notionals
		if( arg[13].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			double notio = arg[13].GetDouble();
			itsNotionals = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,notio));
		}
		else if( arg[13].GetType() ==  GFAT_VECTOR_TYPE )
		{
			itsNotionals = arg[13].GetVector();
			if (itsNotionals->size() != sizeFlows)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "SpreadOption keyword: invalid vector size" );
		}
		else if( arg[13].GetType() ==  GFAT_CURVE_TYPE )
		{
			ARM_GP_CurvePtr notioCurve	= arg[13].GetCurve();
			itsNotionals = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeFlows,0.0));
			for(i=0; i<sizeFlows; ++i)			
				(*itsNotionals)[i] = notioCurve->Interpolate(itsResetTimes[i]);										
		}
		else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + "notio should be double, vector or ARM_Curve" );


		itsLeveragePrev = arg[18].GetDouble();


		/// Validation
		CheckNbSmaller( itsEvalTime, itsResetTimes[0], "EvalTime","1st ResetTime", ARM_GP_SpreadOption::itsFuncName,__LINE__,__FILE__  );																
		CheckVectorPositiveNb( itsPayPeriods,  "interest Periods",ARM_GP_SpreadOption::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsResetTimes,	"ResetTimes",   ARM_GP_SpreadOption::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsPayTimes,	    "PayTimes",	       ARM_GP_SpreadOption::itsFuncName,__LINE__,__FILE__  );
	
		
		/// delete char* for memory leak
		delete payCalendar;
		delete resetCalendar;

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_SpreadOption
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Cap Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_SpreadOption::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{	
	GrabInputs( arg, mod, evalDate, nodes );
			
	ARM_VectorPtr SpreadOption =  itsModelIR->VanillaSpreadOption ( itsCurveName, 
																	itsEvalTime, 
																	itsCallPut, 
																	0, 0, // not used (a priori)...
																	itsResetTimes, 
																	itsPayTimes, 
																	itsPayPeriods, 
																	*itsNotionals, 
																	*itsCoeffLong, 
																	*itsCoeffShort, 
																	*itsStrikes, 
																	itsSwapLongFloatStartTime, 
																	itsSwapLongFloatEndTime, 
																	itsSwapLongFixPayTimes, 
																	itsSwapLongFixPayPeriods, 
																	itsSwapLongFloatStartTime, 
																	itsSwapShortFloatEndTime, 
																	itsSwapShortFixPayTimes, 
																	itsSwapShortFixPayPeriods,
																	itsLeveragePrev,
																	states);

	return ARM_GramFctorArg(SpreadOption);	
}





////////////////////////////////////////////////////
///	Class  : ARM_GP_SpreadOption
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_SpreadOption::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	CC_STL_VECTOR( double ) result;
	//result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
	//result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
	// as only used for analytical (local) models, this will not be used
	result.insert( result.end(), itsPayTimes.begin(), itsPayTimes.end() );

	/// sort and remove duplicates!
	/// not
	CC_NS( std, sort )( result.begin(), result.end() );
	CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
	result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}

////////////////////////////////////////////////////
///////////		Corridor     //////////////////
////////////////////////////////////////////////////


// Parameters
// 0 Model						string
// 1 StartDate					date
// 2 EndDate					date
// 3 RcvPay						int
// 4 PaymentIndexType			int
// 5 FixValue					curve
// 6 PayIndexMult				curve
// 7 PaymentFreq				int
// 8 PaymentIndexTerm			string
// 9 PaymentIndexTiming			int
// 10 PaymentDayCount			int
// 11 IntRule					int
// 12 Spread					curve
// 13 Fixing Frequency			int
// 14 Fixing Timing				int
// 15 Reset Gap					int
// 16 FixingIndexType1			int
// 17 FixingIndexTerm1			string
// 18 Coeff1					curve
// 19 BarriersDown				curve
// 20 BarriersUp				curve
// 21 Notional					curve
// 22 FixingIndexType2			int
// 23 FixingIndexTerm2			string
// 24 Coeff2					curve
// 25 FixingIndexType3			int
// 26 FixingIndexTerm3			string
// 27 BarriersDown3				curve
// 28 BarriersUp3				curve
// 29 DateStripStruct			datestrip
// 30 DateStripPay				datestrip
// 31 DateStripFix				datestrip
// 32 Index						vector 
// 33 OffsetIndex				double
// 34 nbflows					double
//
//
// Elementary flow is :
// (PayIndexMult*PayIndex + Spread) or FixValue * Notional * PeriodRatio *
//				Indic{BarrierDown3 <= Index3 <= BarrierUp3} *
//				Indic{BarrierDown <= Coeff1*Index1-Coeff2*Index2 <= BarrierUp}
//
//
////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: Desctructor
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
ARM_GP_Corridor::~ARM_GP_Corridor()
{
	DeletePointorVector<ARM_GP_Vector>(itsFixingTimes);
	DeletePointorVector<ARM_GP_Vector>(itsFixingWeights);
	DeletePointorVector<ARM_GP_Vector>(itsFirstIndexStartTimes);
	DeletePointorVector<ARM_GP_Vector>(itsFirstIndexEndTimes);
	DeletePointorVector<ARM_GP_Vector>(itsFirstIndexTerms);
	DeletePointorVector<ARM_GP_Vector>(itsSecondIndexStartTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSecondIndexEndTimes);
	DeletePointorVector<ARM_GP_Vector>(itsSecondIndexTerms);
	DeletePointorVector<ARM_GP_Vector>(itsBarriersDown);
	DeletePointorVector<ARM_GP_Vector>(itsBarriersUp);
	DeletePointorVector<ARM_GP_Vector>(itsCoeffs1);
	DeletePointorVector<ARM_GP_Vector>(itsCoeffs2);

	DeletePointorVector<ARM_GP_Vector>(itsBarriersDown3);
	DeletePointorVector<ARM_GP_Vector>(itsBarriersUp3);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: SetDefaults
///	Returns: void
///	Action : set the current defaults to the expression 
///				node tree.
////////////////////////////////////////////////////

void ARM_GP_Corridor::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the minimum required size and type
	GPAF_CheckArgSize( arg, 35, ARM_GP_Corridor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE, ARM_GP_Corridor::itsFuncName );
	itsCurveName = arg[0].GetString();
	ARM_Currency* ccy = mod->GetCurrency( itsCurveName );

	/// Modification of nodes if necessary
	SetEndDatePayCal(nodes, arg, 2, mod, itsCurveName );
	bool isDefault;

	// Payment Index Term correspond to payment Frequency
	int payFreq = ARM_ArgConv_StdFrequency.GetNumber(arg[7].GetString());
	isDefault=SetFrequency( nodes, arg, 8, mod, itsCurveName, (payFreq?payFreq:ccy->GetLiborTerm()));

	// Payment DayCount
	isDefault=SetDayCount(	 nodes, arg, 10, mod, itsCurveName, GetFloatDayCountFtor(ccy));
	// First Index Term
	if (arg[17].GetType() != GFAT_CURVE_TYPE)
		isDefault=SetFrequency(	 nodes, arg, 17, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	// Second Index Term
	isDefault=SetFrequency(	 nodes, arg, 23, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
	// Thrid Index Term
	isDefault=SetFrequency(	 nodes, arg, 26, mod, itsCurveName, GetFloatFrequencyFtor(ccy));
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_Corridor::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	// Fix Value
	GPAF_CheckArgType( arg[5], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// PayIndexMult
	GPAF_CheckArgType( arg[6], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Spread
	GPAF_CheckArgType( arg[13], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Coeff1
	GPAF_CheckArgType( arg[18], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Barrier Down
	GPAF_CheckArgType( arg[19], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Barrier Up
	GPAF_CheckArgType( arg[20], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Notional
	GPAF_CheckArgType( arg[21], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Coeff2
	GPAF_CheckArgType( arg[24], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Barrier Down for Index3
	GPAF_CheckArgType( arg[27], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
	// Barrier Up for Index 3
	GPAF_CheckArgType( arg[28], GFAT_VECTOR_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );


	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		/// set defaults first
		SetDefaults( arg, mod, nodes );
	
		/// checking of the size and type
		GPAF_CheckArgSize( arg, 35, ARM_GP_Corridor::itsFuncName );
		// Model Name
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// StartDate
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,				ARM_GP_Corridor::itsFuncName );
		// EndDate
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,		ARM_GP_Corridor::itsFuncName );
		// RcvPay
		GPAF_CheckArgType( arg[3], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Payment Index Type
		GPAF_CheckArgType( arg[4], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		
		// Payment Freq
		GPAF_CheckArgType( arg[7], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Payment Index Term
		GPAF_CheckArgType( arg[8], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Payment Index Timing
		GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Payment Day Count
		GPAF_CheckArgType( arg[10],GFAT_MULTITOKENSTRING_TYPE,	ARM_GP_Corridor::itsFuncName );
		// Reset Gap
		GPAF_CheckArgType( arg[11],GFAT_DOUBLE_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Int Rule
		GPAF_CheckArgType( arg[12],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Frequency
		GPAF_CheckArgType( arg[14],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Timing
		GPAF_CheckArgType( arg[15],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Index Type1
		GPAF_CheckArgType( arg[16],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Index Term1
		GPAF_CheckArgType( arg[17],GFAT_STRING_OR_CURVE_TYPE,	ARM_GP_Corridor::itsFuncName );
		// Fixing Index Type2
		GPAF_CheckArgType( arg[22],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Index Term2
		GPAF_CheckArgType( arg[23],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Index Type3
		GPAF_CheckArgType( arg[25],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		// Fixing Index Term3
		GPAF_CheckArgType( arg[26],GFAT_STRING_TYPE,			ARM_GP_Corridor::itsFuncName );
		
		
		
		ARM_Currency* ccy			= mod->GetCurrency( itsCurveName );
		char* ccyName				= ccy->GetCcyName();
		int stdSpotGap				= ccy->GetSpotDays();

		char fixCalendar[100];
        ccy->CalcFixPayCal(fixCalendar);
		int  fixFreq				= ccy->GetFixedPayFreq();
		int  fixDayCount			= ccy->GetFixedDayCount();

		char floatCalendar[100];
		ccy->CalcFloatPayCal(floatCalendar);
		int  floatDayCount			= ccy->GetLiborIndexDayCount();

		itsEvalTime					= GetAndCheckEvalTime( evalDate, mod, ARM_GP_Corridor::itsFuncName );

/// grab inputs
		itsCurveName				= arg[0].GetString();
		ARM_Date startDate			= arg[1].GetDate();
		ARM_Date endDate			= arg[2].GetDate(); //end date is not adjusted

		itsRcvPay					= ARM_ArgConv_PayRec.GetNumber(arg[3].GetString());


// PAY INDEX

		itsPayIndexType				= ARM_ArgConv_IndexClass.GetNumber(arg[4].GetString());
		string payIndexTerm;

		if (itsPayIndexType == K_LIBOR || itsPayIndexType == K_FIXED)
		{
			int payIndexFreq		= GetFrequency(arg[8], mod, itsCurveName, ARM_GP_Corridor::itsFuncName, GetFloatFrequencyFtor(ccy));
		//  payIndexTerm			= YearTermToStringMaturity(1.0/payIndexFreq);
			payIndexTerm			= ConvertYearTermToStringMatu(1.0/payIndexFreq);
		}	
		else /// K_CMS case
		{
			double payIndexTermDbl	= GetTerm(arg[8], mod, itsCurveName, ARM_GP_Corridor::itsFuncName, GetFloatFrequencyFtor(ccy));
			//payIndexTerm			= YearTermToStringMaturity(payIndexTermDbl);
			payIndexTerm			= ConvertYearTermToStringMatu(payIndexTermDbl);
		}
		int payIndexTiming			= ARM_ArgConv_Timing.GetNumber(arg[9].GetString());

		ARM_IRIndex payIndex;
		if (itsPayIndexType == K_LIBOR)
			payIndex = ARM_IRIndex((ARM_INDEX_TYPE)FromIndexAndTermToIndexType(payIndexTerm, ccyName));
		else /// we don't care, (dates regenerated later)
			payIndex = ARM_IRIndex(EURIBOR1Y);


// FIRST INDEX

		itsFirstIndexType			= ARM_ArgConv_IndexClass.GetNumber(arg[16].GetString());

		double fixIndexTerm1;
		string fixIndexTerm1Str;
		ARM_GP_CurvePtr tenorCurve;

		bool isVms = false;
		if (arg[17].GetType() == GFAT_STRING_TYPE)
		{
			fixIndexTerm1		= GetTerm(arg[17], mod, itsCurveName, ARM_GP_Corridor::itsFuncName, GetFloatFrequencyFtor(ccy));
			fixIndexTerm1Str	= ConvertYearTermToStringMatu(fixIndexTerm1);
		}
		else
		{
			tenorCurve	= GetCurve(arg[17]);
			isVms		= true;
		}

		ARM_INDEX_TYPE fixIndexType1;
		if (itsFirstIndexType == K_LIBOR)
			fixIndexType1 = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(fixIndexTerm1Str, ccyName);
		else
			fixIndexType1 = EURIBOR1Y;
		ARM_IRIndex fixIndex1(fixIndexType1);


// SECOND INDEX

		itsSecondIndexType			= ARM_ArgConv_IndexClass.GetNumber(arg[22].GetString());
		double fixIndexTerm2 		= GetTerm(arg[23], mod, itsCurveName, ARM_GP_Corridor::itsFuncName, GetFloatFrequencyFtor(ccy));
		//string fixIndexTerm2Str		= YearTermToStringMaturity(fixIndexTerm2);
		string fixIndexTerm2Str		= ConvertYearTermToStringMatu(fixIndexTerm2);

		ARM_INDEX_TYPE fixIndexType2;
		if (itsSecondIndexType == K_LIBOR)
			fixIndexType2 = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(fixIndexTerm2Str, ccyName);
		else
			fixIndexType2 = EURIBOR1Y;
		ARM_IRIndex fixIndex2(fixIndexType2);

// THIRD INDEX (for double corridor)

		itsThirdIndexType			= ARM_ArgConv_IndexClass.GetNumber(arg[25].GetString());
		double fixIndexTerm3		= GetTerm(arg[26], mod, itsCurveName, ARM_GP_Corridor::itsFuncName, GetFloatFrequencyFtor(ccy));
		string fixIndexTerm3Str		= ConvertYearTermToStringMatu(fixIndexTerm3);

		ARM_INDEX_TYPE fixIndexType3;
		if (itsThirdIndexType == K_LIBOR)
			fixIndexType3 = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(fixIndexTerm3Str, ccyName);
		else
			fixIndexType3 = EURIBOR1Y;
		ARM_IRIndex fixIndex3(fixIndexType3);

// identify type of corridor

		bool isCorridor = ( (itsPayIndexType == K_FIXED) || (itsPayIndexType == K_LIBOR) ) && (itsFirstIndexType == K_LIBOR) && (itsSecondIndexType == K_LIBOR);
		
		/// variable corridor CMS spread --> itsPayIndexType can be K_LIBOR or K_CMS ...
		bool isCMSSpreadCorridor = (itsFirstIndexType == K_CMS) || (itsSecondIndexType == K_CMS);

		bool isDbleCorridor = (itsThirdIndexType == K_LIBOR || itsThirdIndexType == K_CMS) && isCMSSpreadCorridor;


		// This validation restricts the usage of the keyword 
		if (!isCorridor  && !isCMSSpreadCorridor)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : At this moment the keyword corridor supports:\n_ Fix or Libor paid with condition on a single Libor\n_ LIBOR/CMS paid with condition on a CMS Spread");
		}

		if (isCorridor  && isCMSSpreadCorridor)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Corridor keyword: could not decide between standard LIBOR corridor and CMS spread corridor...");
		}

		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
		if(!itsModelIR)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for Cap keyword");



		int payFreq					= ARM_ArgConv_StdFrequency.GetNumber(arg[7].GetString());
		int payDayCount				= GetDayCount( arg[10], mod, itsCurveName, GetFloatDayCountFtor(ccy));
		int resetGap				= arg[11].GetDouble();
		int intRule					= ARM_ArgConv_IntRule.GetNumber(arg[12].GetString());
		int fixingFreq				= ARM_ArgConv_StdFrequency.GetNumber(arg[14].GetString());
		int fixingTiming			= ARM_ArgConv_Timing.GetNumber(arg[15].GetString());

		payIndex.SetDayCount(payDayCount);

        
		ARM_Vector* refIndexFixingDates = NULL;
		ARM_Vector* refIndexTimes		= NULL;
		ARM_Vector* tmpRefIndexTimes		= NULL;
		ARM_Vector* refIndexStartDates	= NULL;
		ARM_Vector* refIndexEndDates	= NULL;
		ARM_Vector* refIndexTerms		= NULL;
		ARM_Vector* refIndexWeights		= NULL;

		size_t i,j;

		ARM_ReferenceValue barrierDown(0.0,0.0);
		ARM_ReferenceValue barrierUp(0.0,0.0);

		ARM_DateStripPtr schedStruct(NULL);
        if(arg[29].GetType() == GFAT_DATESTRIP_TYPE)
            schedStruct = arg[29].GetDateStrip();
		ARM_DateStripPtr schedPay(NULL);
        if(arg[30].GetType() == GFAT_DATESTRIP_TYPE)
            schedPay = arg[30].GetDateStrip();
		ARM_DateStripPtr schedFix(NULL);
        if(arg[31].GetType() == GFAT_DATESTRIP_TYPE)
            schedFix = arg[31].GetDateStrip();
		ARM_GP_VectorPtr idx(NULL);
		if(arg[32].GetType() == GFAT_VECTOR_TYPE)
            idx = arg[32].GetVector();

		size_t offset = static_cast<size_t>(arg[33].GetDouble());
		size_t nbFlows = static_cast<size_t>(arg[34].GetDouble());

		ARM_CorridorLeg corridorLeg(
				(ARM_Date) startDate,
				(ARM_Date) endDate,
				itsRcvPay,
				&payIndex,
				payFreq,
				0.0,
				&fixIndex1,
				fixingFreq,
				payIndexTiming,
				fixingTiming,
				K_SHORTSTART,
				&barrierDown,
				K_STD,
				&barrierUp,
				K_STD,
				ccy,
				K_LINEAR,
				K_DIGITALE,
				0);

		bool isStandardKeyWord = ( schedStruct == ARM_DateStripPtr(NULL) || schedPay == ARM_DateStripPtr(NULL) || schedFix == ARM_DateStripPtr(NULL) || idx==ARM_GP_VectorPtr(NULL) );

		if ( isStandardKeyWord )
		{
			refIndexTimes			= corridorLeg.GetRefIndexTimes();
			nbFlows					= refIndexTimes->GetSize();
			

			itsStartTimes			= To_ARM_GP_Vector(corridorLeg.GetFlowStartDates());
			itsEndTimes				= To_ARM_GP_Vector(corridorLeg.GetFlowEndDates());
			itsPayTimes				= To_ARM_GP_Vector(corridorLeg.GetPaymentDates());
			//itsPayPeriods			= To_ARM_GP_Vector(corridorLeg.GetPaymentInterestterms());

			itsPayIndexResetTimes	= To_ARM_GP_Vector(corridorLeg.GetPaymentFixingDates());
			itsPayIndexStartTimes	= To_ARM_GP_Vector(corridorLeg.GetPaymentStartDates());
			//itsPayIndexEndTimes		= To_ARM_GP_Vector(corridorLeg.GetPaymentEndDates());
			itsPayIndexTerms		= To_ARM_GP_Vector(corridorLeg.GetPaymentInterestterms());
		
			refIndexFixingDates		= corridorLeg.GetRefIndexFixingDates();
			refIndexStartDates		= corridorLeg.GetRefIndexStartDates();
			if (itsFirstIndexType ==K_LIBOR)
			{
				refIndexEndDates	= corridorLeg.GetRefIndexEndDates();
				refIndexTerms		= corridorLeg.GetRefIndexInterestTerms();
			}
			
			size_t	nbFixings		= refIndexFixingDates->size();
			double* pw				= corridorLeg.GetRefIndexWeight();
			refIndexWeights			= new ARM_Vector(nbFixings);
			for (i = 0; i <nbFixings; ++i)
				(*refIndexWeights)[i]  = pw[i];


			
			// Take care !
			// This is adapted to the strange behaviour of the corridor leg
			for (i = 0; i < nbFlows; ++i)
			{
				ARM_Date date((*refIndexFixingDates)[i]);

				date.GapBusinessDay(resetGap, fixCalendar);

				(*refIndexFixingDates)[i] = date.GetJulian();
			}

		}
		else
		{
			refIndexTimes			= To_pARM_Vector(&*idx);
			tmpRefIndexTimes		= refIndexTimes;
			if (nbFlows==0)
				nbFlows = refIndexTimes->GetSize()-offset;

			size_t nbAux = 0;
			for (i=0;i<nbFlows;i++)
				nbAux += (*refIndexTimes)[offset+i];

			size_t dec=0;
			for (i=0;i<offset;i++)
				dec+=(*refIndexTimes)[i];

			ARM_GP_Vector* aux;

			aux = SubVector(offset,*schedStruct->GetFlowStartDates(),nbFlows);
			itsStartTimes			= ARM_GP_Vector(*aux);	delete aux;
			aux = SubVector(offset,*schedStruct->GetFlowEndDates(),nbFlows);
			itsEndTimes				= ARM_GP_Vector(*aux);	delete aux;
			aux = SubVector(offset,*schedStruct->GetPaymentDates(),nbFlows);
			itsPayTimes				= ARM_GP_Vector(*aux);	delete aux;
			/*aux = SubVector(offset,*schedStruct->GetInterestTerms());
			itsPayPeriods			= ARM_GP_Vector(*aux);	delete aux;*/

			aux = SubVector(offset,*schedPay->GetResetDates(),nbFlows);
			itsPayIndexResetTimes	= ARM_GP_Vector(*aux);	delete aux;
			aux = SubVector(offset,*schedPay->GetFlowStartDates(),nbFlows);
			itsPayIndexStartTimes	= ARM_GP_Vector(*aux);	delete aux;
			/*aux = SubVector(offset,*schedPay->GetFlowEndDates());
			itsPayIndexEndTimes		= ARM_GP_Vector(*aux);	delete aux;*/
			aux = SubVector(offset,*schedPay->GetInterestTerms(),nbFlows);
			itsPayIndexTerms		= ARM_GP_Vector(*aux);	delete aux;

			aux = SubVector(dec,*schedFix->GetResetDates(),nbAux);
			refIndexFixingDates		= To_pARM_Vector(aux);	delete aux;
			//aux = SubVector(dec,*schedFix->GetFlowStartDates());
			aux = SubVector(dec,*schedFix->GetFwdRateStartDates(),nbAux);
			refIndexStartDates		= To_pARM_Vector(aux);	delete aux;

			if (itsFirstIndexType ==K_LIBOR)
			{
				//aux = SubVector(dec,*schedFix->GetFlowEndDates());
				aux = SubVector(dec,*schedFix->GetFwdRateEndDates(),nbAux);
				refIndexEndDates	= To_pARM_Vector(aux);	delete aux;
				aux = SubVector(dec,*schedFix->GetInterestTerms(),nbAux);
				refIndexTerms		= To_pARM_Vector(aux);	delete aux;
			}
			aux = SubVector(dec,*schedFix->GetInterestDays(),nbAux);	//ugly : container for cpn fixing
			refIndexWeights			= To_pARM_Vector(aux);	delete aux;

		}

		
		ARM_GP_CurvePtr fixValCurve			= GetCurve(arg[5]);
		ARM_GP_CurvePtr payIndexMultCurve	= GetCurve(arg[6]);
		ARM_GP_CurvePtr spreadCurve			= GetCurve(arg[13]);
		ARM_GP_CurvePtr coeff1Curve			= GetCurve(arg[18]);
		ARM_GP_CurvePtr barrierDownCurve	= GetCurve(arg[19]);
		ARM_GP_CurvePtr barrierUpCurve		= GetCurve(arg[20]);
		ARM_GP_CurvePtr notionalCurve		= GetCurve(arg[21]);
		ARM_GP_CurvePtr coeff2Curve			= GetCurve(arg[24]);
		/// Double condition case
		ARM_GP_CurvePtr barrierDown3Curve	= GetCurve(arg[27]);
		ARM_GP_CurvePtr barrierUp3Curve		= GetCurve(arg[28]);
		
		DeletePointorVector<ARM_GP_Vector>(itsFixingTimes);
		DeletePointorVector<ARM_GP_Vector>(itsFixingWeights);
		DeletePointorVector<ARM_GP_Vector>(itsFirstIndexStartTimes);
		DeletePointorVector<ARM_GP_Vector>(itsFirstIndexEndTimes);
		DeletePointorVector<ARM_GP_Vector>(itsFirstIndexTerms);
		DeletePointorVector<ARM_GP_Vector>(itsSecondIndexStartTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSecondIndexEndTimes);
		DeletePointorVector<ARM_GP_Vector>(itsSecondIndexTerms);
		DeletePointorVector<ARM_GP_Vector>(itsBarriersDown);
		DeletePointorVector<ARM_GP_Vector>(itsBarriersUp);
		DeletePointorVector<ARM_GP_Vector>(itsCoeffs1);
		DeletePointorVector<ARM_GP_Vector>(itsCoeffs2);
		/// Double condition case
		DeletePointorVector<ARM_GP_Vector>(itsBarriersDown3);
		DeletePointorVector<ARM_GP_Vector>(itsBarriersUp3);

		

		itsFixingTimes.resize(nbFlows);
		itsFixingWeights.resize(nbFlows);
		if (isCorridor)
		{
			itsFirstIndexStartTimes.resize(nbFlows);
			itsFirstIndexEndTimes.resize(nbFlows);
			itsFirstIndexTerms.resize(nbFlows);
// FIXMEFRED: mig.vc8 (30/05/2007 16:40:43):pointer function
			itsCorridorFunc = &ARM::ARM_GP_Corridor::EvalCorridor;
		}
		else if (isCMSSpreadCorridor)
		{
			itsFirstIndexSwapRates.resize(nbFlows);
			itsSecondIndexSwapRates.resize(nbFlows);
// FIXMEFRED: mig.vc8 (30/05/2007 16:40:43):pointer function
			itsCorridorFunc = &ARM::ARM_GP_Corridor::EvalCMSSpreadCorridor;
			itsPayRates.resize(nbFlows);
		}

		/// Vectors are initialized but will contain nothing if double condition is not activated
		itsThirdIndexSwapRates.resize(nbFlows); /// Each element is an empty vector
		itsBarriersDown3.resize(nbFlows);
		itsBarriersUp3.resize(nbFlows);

		itsFixValues.resize(nbFlows);
		itsPayIndexMults.resize(nbFlows);
		itsSpreads.resize(nbFlows);
		itsNotionals.resize(nbFlows);
		itsBarriersDown.resize(nbFlows);
		itsBarriersUp.resize(nbFlows);
		itsCoeffs1.resize(nbFlows);
		itsCoeffs2.resize(nbFlows);
		
		int fixingPos = 0;
		
		for (i = 0; i < nbFlows; ++i)
		{
			int size_i = (*refIndexTimes)[offset+i];
			itsFixingTimes[i] = new ARM_GP_Vector(size_i);
			itsFixingWeights[i] = new ARM_GP_Vector(size_i);
			if (isCorridor)
			{
				itsFirstIndexStartTimes[i] = new ARM_GP_Vector(size_i);
				itsFirstIndexEndTimes[i] = new ARM_GP_Vector(size_i);
				itsFirstIndexTerms[i] = new ARM_GP_Vector(size_i);
			}
			else if (isCMSSpreadCorridor)
			{
				itsFirstIndexSwapRates[i].resize(size_i);
				itsSecondIndexSwapRates[i].resize(size_i);
			}

			if (isDbleCorridor)
			{
				itsThirdIndexSwapRates[i].resize(size_i);
				itsBarriersDown3[i] = new ARM_GP_Vector(size_i);
				itsBarriersUp3[i] = new ARM_GP_Vector(size_i);
			}
			else /// Just for compatibility
			{
				itsBarriersDown3[i] = new ARM_GP_Vector(1,0.0);
				itsBarriersUp3[i] = new ARM_GP_Vector(1,0.0);
			}

			itsBarriersDown[i] = new ARM_GP_Vector(size_i);
			itsBarriersUp[i] = new ARM_GP_Vector(size_i);
			itsCoeffs1[i] = new ARM_GP_Vector(size_i);
			itsCoeffs2[i] = new ARM_GP_Vector(size_i);

			itsStartTimes[i] = mod->GetTimeFromDate(itsStartTimes[i]);

			for (j = 0; j < size_i; ++j)
			{
				(*itsFixingTimes[i])[j] = mod->GetTimeFromDate((*refIndexFixingDates)[fixingPos+j]);
				(*itsFixingWeights[i])[j] = (*refIndexWeights)[fixingPos+j];


				if (isCorridor)
				{
					(*itsFirstIndexStartTimes[i])[j] = mod->GetTimeFromDate((*refIndexStartDates)[fixingPos+j]);
					(*itsFirstIndexEndTimes[i])[j] = mod->GetTimeFromDate((*refIndexEndDates)[fixingPos+j]);
					(*itsFirstIndexTerms[i])[j] = (*refIndexTerms)[fixingPos+j];

					/// specific interp to be consistent with ARM_CorridorLeg
					(*itsBarriersDown[i])[j] = barrierDownCurve->Interpolate((*itsFixingTimes[i])[j]);
					(*itsBarriersUp[i])[j] = barrierUpCurve->Interpolate((*itsFixingTimes[i])[j]);
				
				}
				else if (isCMSSpreadCorridor)
				{	
					/// specific interp to be consistent with ARM_SpreadOption
					(*itsBarriersDown[i])[j] = barrierDownCurve->Interpolate(itsStartTimes[i]);
					(*itsBarriersUp[i])[j] = barrierUpCurve->Interpolate(itsStartTimes[i]);

					ARM_Date startDate ((*refIndexStartDates)[fixingPos+j]);
					
					ARM_Date endDate1(startDate);
					if (isVms)
					{
						int vmsIndex = (int) tenorCurve->Interpolate(itsStartTimes[i] );
						int intVmsType;
						string vmsTerm = FromIndexTypeToTermAndType(vmsIndex, intVmsType);
						endDate1.AddPeriod(vmsTerm);
					}
					else
						endDate1.AddPeriod(fixIndexTerm1Str);

					if (itsFirstIndexType == K_CMS)
					{						
						itsFirstIndexSwapRates[i][j] = ARM_SwapRate::CreateSwapRate(
							mod->GetAsOfDate().GetJulian(),
							(*refIndexStartDates)[fixingPos+j],
							endDate1.GetJulian(),
							fixDayCount,
							fixFreq,
							fixCalendar);
					}
					else if (itsFirstIndexType == K_LIBOR)
					{							
						endDate1.GoodBusinessDay( K_MOD_FOLLOWING, floatCalendar );
						double payPeriod = CountYears(floatDayCount, startDate.GetJulian(), endDate1.GetJulian());
						ARM_SwapRate* sr = new ARM_SwapRate(mod->GetAsOfDate().GetJulian(),
															(*refIndexFixingDates)[fixingPos+j],
															startDate.GetJulian(), 
															endDate1.GetJulian(),
															ARM_GP_Vector(1, endDate1.GetJulian()),
															ARM_GP_Vector(1, payPeriod));
						itsFirstIndexSwapRates[i][j] = ARM_SwapRatePtr(sr);
					}


					ARM_Date endDate2(startDate);
					endDate2.AddPeriod(fixIndexTerm2Str);

					if (itsSecondIndexType == K_CMS)
					{							
						itsSecondIndexSwapRates[i][j] = ARM_SwapRate::CreateSwapRate(
							mod->GetAsOfDate().GetJulian(),
							(*refIndexStartDates)[fixingPos+j],
							endDate2.GetJulian(),
							fixDayCount,
							fixFreq,
							fixCalendar);		
					}
					else if (itsSecondIndexType == K_LIBOR)
					{
						endDate2.GoodBusinessDay( K_MOD_FOLLOWING, floatCalendar );
						double payPeriod = CountYears(floatDayCount, startDate.GetJulian(), endDate2.GetJulian());
						ARM_SwapRate* sr = new ARM_SwapRate(mod->GetAsOfDate().GetJulian(),
															(*refIndexFixingDates)[fixingPos+j],
															startDate.GetJulian(), 
															endDate2.GetJulian(),
															ARM_GP_Vector(1, endDate2.GetJulian()),
															ARM_GP_Vector(1, payPeriod));
						itsSecondIndexSwapRates[i][j] = ARM_SwapRatePtr(sr);
					}
				}
				(*itsCoeffs1[i])[j] = coeff1Curve->Interpolate(itsStartTimes[i]);
				(*itsCoeffs2[i])[j] = coeff2Curve->Interpolate(itsStartTimes[i]);

				if (isDbleCorridor)
				{
					ARM_Date resetDate ((*refIndexFixingDates)[fixingPos+j]);
					ARM_Date startDate ((*refIndexStartDates)[fixingPos+j]);
					
					ARM_Date endDate(startDate);
					endDate.AddPeriod(fixIndexTerm3Str);

					if (itsThirdIndexType == K_CMS)
					{						
						itsThirdIndexSwapRates[i][j] = ARM_SwapRate::CreateSwapRate(
							mod->GetAsOfDate().GetJulian(),
							startDate.GetJulian(),
							endDate.GetJulian(),
							fixDayCount,
							fixFreq,
							fixCalendar);
					}
					else if (itsThirdIndexType == K_LIBOR)
					{							
						endDate.GoodBusinessDay( K_MOD_FOLLOWING, floatCalendar );
						double payPeriod = CountYears(floatDayCount, startDate.GetJulian(), endDate.GetJulian());
						ARM_SwapRate* sr = new ARM_SwapRate(mod->GetAsOfDate().GetJulian(),
															resetDate.GetJulian(),
															startDate.GetJulian(), 
															endDate.GetJulian(),
															ARM_GP_Vector(1, endDate.GetJulian()),
															ARM_GP_Vector(1, payPeriod));
						itsThirdIndexSwapRates[i][j] = ARM_SwapRatePtr(sr);
					}
					(*itsBarriersDown3[i])[j] = barrierDown3Curve->Interpolate(itsStartTimes[i]);
					(*itsBarriersUp3[i])[j] = barrierUp3Curve->Interpolate(itsStartTimes[i]);
				}	
			}
			
			/// case of a CMS spread or double corridor with variable payments (CMS or LIBOR)
			/// ---> build corresp ARM_SwapRate objects
			if (isCMSSpreadCorridor)
			{
				if (itsPayIndexType == K_FIXED)
				{
					itsPayRates[i] = ARM_SwapRatePtr(NULL);
				}
				else if (itsPayIndexType == K_LIBOR)
				{				
					/// at this step, all xxxTimes below are in fact dates...
					ARM_Date endDate (itsPayIndexStartTimes[i]/*in fact it's a date at this step*/);
					endDate.AddPeriod(payIndexTerm);
					endDate.GoodBusinessDay( K_MOD_FOLLOWING, floatCalendar );
					double payPeriod = CountYears(floatDayCount, itsPayIndexStartTimes[i], endDate.GetJulian());

					ARM_SwapRate* swapRate	= new ARM_SwapRate (mod->GetAsOfDate().GetJulian(),
																itsPayIndexResetTimes[i],
																itsPayIndexStartTimes[i],
																endDate.GetJulian(),
																ARM_GP_Vector(1, endDate.GetJulian()),
																ARM_GP_Vector(1, payPeriod) );
				
					itsPayRates[i] = ARM_SwapRatePtr (swapRate);
				}
				else if (itsPayIndexType == K_CMS)
				{
					ARM_Date endDate (itsPayIndexStartTimes[i]/*in fact it's a date at this step*/);
					endDate.AddPeriod(payIndexTerm);
					itsPayRates[i] =  ARM_SwapRate::CreateSwapRate(
															mod->GetAsOfDate().GetJulian(),
															itsPayIndexStartTimes[i],
															endDate.GetJulian(),
															fixDayCount,
															fixFreq,
															fixCalendar);
				}
			}

			itsEndTimes[i]				= mod->GetTimeFromDate(itsEndTimes[i]);
			itsPayIndexResetTimes[i]	= mod->GetTimeFromDate(itsPayIndexResetTimes[i]);
			itsPayIndexStartTimes[i]	= mod->GetTimeFromDate(itsPayIndexStartTimes[i]);
			//itsPayIndexEndTimes[i]		= mod->GetTimeFromDate(itsPayIndexEndTimes[i]);
			itsPayTimes[i]				= mod->GetTimeFromDate(itsPayTimes[i]);

			itsFixValues[i]				= fixValCurve->Interpolate(itsStartTimes[i]);
			itsPayIndexMults[i]			= payIndexMultCurve->Interpolate(itsStartTimes[i]);
			itsSpreads[i]				= spreadCurve->Interpolate(itsStartTimes[i]);
			itsNotionals[i]				= notionalCurve->Interpolate(itsPayTimes[i]);
			
			fixingPos += size_i;
		}

		SetAlreadyComputed(true);
		delete refIndexWeights;
		if (!isStandardKeyWord)
		{
			delete refIndexFixingDates;
			delete refIndexStartDates;
			if (itsFirstIndexType ==K_LIBOR)
			{
				delete refIndexEndDates;
				delete refIndexTerms;
			}
		}

		if (tmpRefIndexTimes)
			delete tmpRefIndexTimes;
		tmpRefIndexTimes = NULL;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: EvalCorridor
///	Returns: ARM_GramFctorArg
///	Action : computes the Corridor Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Corridor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	return (this->*itsCorridorFunc)(arg, mod, evalDate, states, nodes);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: EvalCorridor
///	Returns: ARM_GramFctorArg
///	Action : computes the Corridor Function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Corridor::EvalCorridor( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{	
	GrabInputs( arg, mod, evalDate, nodes );
		
	ARM_VectorPtr result;

	result = itsModelIR->VanillaCorridor(
			itsCurveName,
			itsEvalTime,
			itsPayTimes,
			itsPayIndexResetTimes,
			itsStartTimes,
			itsEndTimes,
			itsPayIndexType,
			itsPayIndexTerms,
			itsFixingTimes,
			itsFirstIndexStartTimes,
			itsFirstIndexEndTimes,
			itsFirstIndexTerms,
			itsFixingWeights,
			(itsPayIndexType == K_FIXED ? itsFixValues : itsSpreads),
			itsBarriersDown,
			itsBarriersUp,
			itsNotionals,
			itsRcvPay,
			states);

	return ARM_GramFctorArg(result);	
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: EvalCMSSpreadCorridor
///	Returns: ARM_GramFctorArg
///	Action : computes the CMS Spread Corridor Function 
/// and return the corresponding values.
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_Corridor::EvalCMSSpreadCorridor( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{	
	GrabInputs( arg, mod, evalDate, nodes );
		
	ARM_VectorPtr result;
	
	/// create payIndexType (vectorial)
	ARM_IntVector payIndexType (itsStartTimes.size(), itsPayIndexType);
		
	result = itsModelIR->VanillaCMSCorridor(
			itsCurveName,
			itsEvalTime,
			itsPayTimes,
			itsPayIndexResetTimes,
			itsStartTimes,
			itsEndTimes,
			itsFixingTimes,
			itsFixingWeights,
			itsCoeffs1,
			itsFirstIndexSwapRates,
			itsCoeffs2,
			itsSecondIndexSwapRates,
			payIndexType,
			itsFixValues,
			itsPayRates,
			itsPayIndexMults,
			itsBarriersDown,
			itsBarriersUp,
			itsNotionals,
			itsRcvPay,
			itsThirdIndexSwapRates,
			itsBarriersDown3,
			itsBarriersUp3,
			states);

	return ARM_GramFctorArg(result);	
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_Corridor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_Corridor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	CC_STL_VECTOR( double ) result;
	result.insert( result.end(), itsStartTimes.begin(), itsStartTimes.end() );
	result.insert( result.end(), itsEndTimes.begin(), itsEndTimes.end() );
	result.insert( result.end(), itsPayTimes.begin(), itsPayTimes.end() );

	for (size_t i = 0; i < itsFirstIndexEndTimes.size(); ++i)
	{
		result.push_back((*itsFirstIndexEndTimes[i])[(*itsFirstIndexEndTimes[i]).size()-1]);
	}

	/// sort and remove duplicates!
	CC_NS( std, sort )( result.begin(), result.end() );
	CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
	result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
///////////		MaxRate     //////////////////
////////////////////////////////////////////////////


void ARM_GP_MaxRate::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
}

////////////////////////////////////////////////////
///	Class  : ARM_GP_MaxRate
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////

void ARM_GP_MaxRate::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// checking of the size and type
		GPAF_CheckArgSize( arg, 12, ARM_GP_MaxRate::itsFuncName  );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,	ARM_GP_MaxRate::itsFuncName ); // modelName
		GPAF_CheckArgType( arg[1], GFAT_DATE_TYPE,		ARM_GP_MaxRate::itsFuncName ); // Swap StartDate
		GPAF_CheckArgType( arg[2], GFAT_DATEORMATU_TYPE,ARM_GP_MaxRate::itsFuncName ); // Swap EndDate/Maturity
		GPAF_CheckArgType( arg[3], GFAT_DATEORMATU_TYPE,ARM_GP_MaxRate::itsFuncName ); // First ResetDate
		GPAF_CheckArgType( arg[4], GFAT_DATEORMATU_TYPE,ARM_GP_MaxRate::itsFuncName ); // First StartDate
		GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE,ARM_GP_MaxRate::itsFuncName );	   // Type MIN,MAX,MAX-MIN option
		GPAF_CheckArgType( arg[7], GFAT_STRING_TYPE,ARM_GP_MaxRate::itsFuncName);	   // Reset frequency D,W,M etc...
		GPAF_CheckArgType( arg[8], GFAT_VECTOR_TYPE,ARM_GP_MaxRate::itsFuncName);	   // Strike (option case)
		GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE,ARM_GP_MaxRate::itsFuncName );	   // Cap/Floor (option case)
		GPAF_CheckArgType( arg[10], GFAT_DOUBLE_TYPE,ARM_GP_MaxRate::itsFuncName );	   // Min vs Max correl (option case)
		GPAF_CheckArgType( arg[11], GFAT_VECTOR_TYPE,ARM_GP_MaxRate::itsFuncName );	   // Min/Max accrued (if activated)

		/// grab inputs
		if( arg[5].GetType() ==  GFAT_VECTOR_TYPE )
			itsFirstRate = arg[5].GetVector();
		else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : previous rate should be vector or double");

		string maxormin			= arg[6].GetString();
		maxormin = stringGetUpper(maxormin);

		if(maxormin != "MAX" && maxormin != "MIN" && maxormin != "MAXMIN")
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " type should be 'max', 'min' or 'maxmin'");
		itsMaxOrMin	= maxormin == "MAX" ? 1 : (maxormin == "MIN" ?  -1 : 0);

		itsResetFreq = ARM_ArgConv_StdFrequency.GetNumber(arg[7].GetString());

		if( arg[8].GetType() ==  GFAT_DOUBLE_TYPE )
			itsStrikes = ARM_VectorPtr(new ARM_GP_Vector(itsFirstRate->size(),arg[8].GetDouble()));
		else
		{
			itsStrikes = arg[8].GetVector();
			if(itsStrikes->size() != itsFirstRate->size())
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between strike and first rate states");
		}

		string capOrFloor		= arg[9].GetString();
		capOrFloor = stringGetUpper(capOrFloor);

		itsRhoMinMax = arg[10].GetDouble();

		if( arg[11].GetType() ==  GFAT_DOUBLE_TYPE )
		{
			itsIsAccrued	= false;
			itsMinAccrued	= ARM_NumericConstants::ARM_INFINITY;
			itsMaxAccrued	= -ARM_NumericConstants::ARM_INFINITY;
		}
		else
		{
			ARM_VectorPtr accruedDatas = arg[11].GetVector();
			itsIsAccrued = true;
			if(accruedDatas->size()<2)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : min and max accrued values are expected");
			itsMinAccrued = (*accruedDatas)[0];
			itsMaxAccrued = (*accruedDatas)[1];
		}

		if(capOrFloor != "CAP" && capOrFloor != "FLOOR")
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "option type should be 'cap' or 'floor'");
		itsCapOrFloor	= capOrFloor == "CAP" ? K_CAP : K_FLOOR;


		itsCurveName			= arg[0].GetString();
		ARM_Date startDate		= arg[1].GetDate();
        
		/// compute inputs
		/// currency and calendar are in ARM world similar ...
		itsModelIR				= dynamic_cast< ARM_PricingFunctionIR* >(mod);
        if(!itsModelIR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only IR model for SPREAD keyword");
		
		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_MaxRate::itsFuncName );

		ARM_Currency* ccy		    = mod->GetCurrency( itsCurveName );
		char* ccyName			    = ccy->GetCcyName();

		ARM_INDEX_TYPE defaultIndex =  GetDefaultIndexFromCurrency( ccyName );
		
		char fixCalendar[100];
		ccy->CalcFixPayCal(fixCalendar);

		char* floatPayCalendar	= ccy->GetPayCalName(defaultIndex);
		char* floatResetCalendar= ccy->GetResetCalName(defaultIndex);
		CC_NS(std,auto_ptr)<char> holdfloatPayCalendar(floatPayCalendar);
		CC_NS(std,auto_ptr)<char> holdfloatResetCalendar(floatResetCalendar);
		

/**/	int fixFreq	    = GetFixedFrequencyFtor(ccy)();
/**/	int fixDayCount = GetFixedDayCountFtor(ccy)();
/**/	int floatFreq = GetFloatFrequencyFtor(ccy)();
/**/	int floatDayCount = GetFloatDayCountFtor(ccy)();
/**/	int fwdDayCount	= GetFloatDayCountFtor(ccy)();

/// First Swap Rate
		ARM_Date endDate = ComputeEndDate( startDate, arg[2], fixCalendar );

		//FIXLEG
		ARM_DateStrip fixDateStrip = ARM_DateStrip( startDate, endDate, fixFreq, fixDayCount, fixCalendar,
										K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
										fixCalendar );

		/// copy constructor
		itsFixResetTimes	= ARM_GP_Vector( *(fixDateStrip.GetResetDates()) );  
		itsFixPayTimes	    = ARM_GP_Vector( *(fixDateStrip.GetPaymentDates()) );
		itsFixPayPeriods    = ARM_GP_Vector( *(fixDateStrip.GetInterestTerms()) );

		/// get time from date!
		for( int i(0); i<itsFixPayTimes.size(); ++i )
		{
			itsFixResetTimes[i] = mod->GetTimeFromDate( itsFixResetTimes[i] );
			itsFixPayTimes[i] = mod->GetTimeFromDate( itsFixPayTimes[i] );
		}
	
		//FLOATLEG
		ARM_DateStrip floatDateStrip = ARM_DateStrip( startDate, endDate, floatFreq, floatDayCount, floatResetCalendar,
											K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, floatFreq, GETDEFAULTVALUE,
											floatPayCalendar );
		
		/// copy constructor
		itsFwdResetTimes	= ARM_GP_Vector( *(floatDateStrip.GetResetDates()) );
		itsFwdStartTimes	= ARM_GP_Vector( *(floatDateStrip.GetFwdRateStartDates()) );
		itsFwdEndTimes		= ARM_GP_Vector( *(floatDateStrip.GetFwdRateEndDates()) );
		itsFwdPayPeriods	= ARM_GP_Vector( itsFwdStartTimes.size());
		itsFloatPayTimes	= ARM_GP_Vector( *(floatDateStrip.GetPaymentDates()) );
		itsFloatPayPeriods	= ARM_GP_Vector( *(floatDateStrip.GetInterestTerms()) );

		/// get time from date! 
		for(i=0; i<itsFwdStartTimes.size(); ++i )
		{
		    itsFwdPayPeriods[i] = CountYearsWithoutException( fwdDayCount, ARM_Date(itsFwdStartTimes[i]), ARM_Date(itsFwdEndTimes[i]) );
			itsFwdResetTimes[i] = mod->GetTimeFromDate( itsFwdResetTimes[i] );
			itsFwdStartTimes[i] = mod->GetTimeFromDate( itsFwdStartTimes[i] );
			itsFwdEndTimes[i]   = mod->GetTimeFromDate( itsFwdEndTimes[i] );
			itsFloatPayTimes[i] = mod->GetTimeFromDate( itsFloatPayTimes[i] );
		}

		ARM_Date floatStartDate( startDate );
		itsFloatStartTime   = mod->GetTimeFromDate( floatStartDate );
		ARM_Date floatEndDate( endDate );
		itsFloatEndTime	= mod->GetTimeFromDate( floatEndDate );

		/// Case of Margin Vector
		size_t sizefloatflows = itsFloatPayTimes.size();        

		double margin = 0.;
		itsMarginVector = ARM_GP_VectorPtr(new ARM_GP_Vector(sizefloatflows,margin));
		
		/// Validation
		CheckNbSmaller( itsEvalTime,	itsFloatStartTime,		"EvalTime",			"FloatStartTime",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime, itsFloatEndTime,		"FloatStartTime",	"FloatEndTime",		ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFloatStartTime, itsFixPayTimes[0],	"FloatStartTime",	"1st Fix PayTimes",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__ );

        CheckNbSmaller( itsEvalTime,	itsFwdStartTimes[0],		"EvalTime",				"1st Fwd StartTimes1",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__ );
		CheckNbSmaller( itsFwdStartTimes[0], itsFixPayTimes[0] ,	"1st Fwd StartTimes",	"1st Fix PayTimes1",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__ );

		CheckVectorIncreasing( itsFixPayTimes,		"FixPayTimes",		ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFixPayPeriods,	"FixPayPeriods",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );

		CheckVectorIncreasing( itsFwdStartTimes,	"FwdStartTimes",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFwdEndTimes,		"FwdEndTimes",		ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorIncreasing( itsFloatPayTimes,	"FloatPayTimes",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );
		CheckVectorPositiveNb( itsFloatPayPeriods,	"FixPayPeriods",	ARM_GP_MaxRate::itsFuncName,__LINE__,__FILE__  );

		ARM_Date firstReset = arg[3].GetDate();
		ARM_Date firstStart = arg[4].GetDate();
		
		itsFirstResetTime = mod->GetTimeFromDate( firstReset );
		itsFirstStartTime = mod->GetTimeFromDate( firstStart );

		SetAlreadyComputed(true);
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_GP_MaxRate
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Spread function and return 
///				the corresponding values
////////////////////////////////////////////////////

ARM_GramFctorArg ARM_GP_MaxRate::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	ARM_VectorPtr MaxRate =  itsModelIR->MaxRate( itsCurveName, itsEvalTime,
        itsFloatStartTime, itsFloatEndTime,
		itsFixPayTimes, itsFixPayPeriods,
        itsFwdStartTimes, itsFwdEndTimes, itsFwdPayPeriods,
        itsFloatPayTimes, itsFloatPayPeriods, *itsMarginVector,
		itsFirstResetTime, itsFirstStartTime, itsFirstRate, itsMaxOrMin,
		itsResetFreq,itsStrikes,itsCapOrFloor,itsRhoMinMax,
		itsIsAccrued,itsMinAccrued,itsMaxAccrued,
		states );

	return ARM_GramFctorArg(MaxRate);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_MaxRate
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_MaxRate::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

    CC_STL_VECTOR( double ) result;

    result.push_back( itsFloatStartTime );
    result.insert( result.end(), itsFwdStartTimes.begin(), itsFwdStartTimes.end() );
    result.insert( result.end(), itsFwdEndTimes.begin(), itsFwdEndTimes.end() );
    result.insert( result.end(), itsFloatPayTimes.begin(), itsFloatPayTimes.end() );
    result.push_back( itsFloatEndTime );
    result.insert( result.end(), itsFixPayTimes.begin(), itsFixPayTimes.end() );

	/// sort and remove duplicates!
    CC_NS( std, sort )( result.begin(), result.end() );
    CC_STL_VECTOR( double )::iterator pos = CC_NS( std, unique )( result.begin(), result.end() );
    result.resize( CC_NS( std, CC_DISTANCE )( result.begin(), pos ) );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}


////////////////////////////////////////////////////
///////////		DoubleDigital     //////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_GP_DoubleDigital
///	Routine: SetDefaults
///	Returns: void
///	Action : no default
////////////////////////////////////////////////////
void ARM_GP_DoubleDigital::SetDefaults( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	vector< ARM_ExpNodePtr >& nodes )
{
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DoubleDigital
///	Routine: GrabInputs
///	Returns: void
///	Action : get the inputs from vector of arg
////////////////////////////////////////////////////
void ARM_GP_DoubleDigital::GrabInputs( const ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,	double evalDate, vector< ARM_ExpNodePtr >& nodes )
{
	/// if it is not already computed
	if(!GetAlreadyComputed())
	{
		SetDefaults( arg, mod, nodes );

		/// Checking of the size and type
		GPAF_CheckArgSize( arg, 9, ARM_GP_DoubleDigital::itsFuncName  );
		GPAF_CheckArgType( arg[0], GFAT_STRING_TYPE,ARM_GP_DoubleDigital::itsFuncName );	// modelName

		GPAF_CheckArgType( arg[1], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 1st Rate
		GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 1st Strike Down
		GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 1st Strike Up

		GPAF_CheckArgType( arg[4], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 2nd Rate
		GPAF_CheckArgType( arg[5], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 2nd Strike Down
		GPAF_CheckArgType( arg[6], GFAT_VECTOR_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 2nd Strike Up

		GPAF_CheckArgType( arg[7], GFAT_DOUBLE_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 1st Strike Spread
		GPAF_CheckArgType( arg[8], GFAT_DOUBLE_TYPE,ARM_GP_DoubleDigital::itsFuncName);		// 2nd Strike Spread

		/// Grad standard inputs
		itsModelName = arg[0].GetString();
		ARM_MultiAssetsModel* mamo	= dynamic_cast< ARM_MultiAssetsModel* >(mod);
        if(!mamo)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only MultiAssetModel allowed for DoubleDigital keyword");
		itsModelIR = static_cast< ARM_PricingModelIR* >(mamo);
		
		itsEvalTime	= GetAndCheckEvalTime( evalDate, mod, ARM_GP_DoubleDigital::itsFuncName );

		/// Grab 1st rate inputs
		itsFirstRate = arg[1].GetVector();
		if(arg[2].GetType() ==  GFAT_DOUBLE_TYPE)
			itsFirstStrikeDown = ARM_VectorPtr(new ARM_GP_Vector(itsFirstRate->size(),arg[2].GetDouble()));
		else
		{
			itsFirstStrikeDown = arg[2].GetVector();
			if(itsFirstStrikeDown->size() != itsFirstRate->size())
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between strike and first rate states");
		}
		if(arg[3].GetType() ==  GFAT_DOUBLE_TYPE)
			itsFirstStrikeUp = ARM_VectorPtr(new ARM_GP_Vector(itsFirstRate->size(),arg[3].GetDouble()));
		else
		{
			itsFirstStrikeUp = arg[3].GetVector();
			if(itsFirstStrikeUp->size() != itsFirstRate->size())
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between strike and first rate states");
		}

		/// Grab 2nd rate inputs
		itsSecondRate = arg[4].GetVector();
		if(arg[5].GetType() ==  GFAT_DOUBLE_TYPE)
			itsSecondStrikeDown = ARM_VectorPtr(new ARM_GP_Vector(itsSecondRate->size(),arg[5].GetDouble()));
		else
		{
			itsSecondStrikeDown = arg[5].GetVector();
			if(itsSecondStrikeDown->size() != itsSecondRate->size())
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between strike and second rate states");
		}
		if(arg[6].GetType() ==  GFAT_DOUBLE_TYPE)
			itsSecondStrikeUp = ARM_VectorPtr(new ARM_GP_Vector(itsSecondRate->size(),arg[6].GetDouble()));
		else
		{
			itsSecondStrikeUp = arg[6].GetVector();
			if(itsSecondStrikeUp->size() != itsSecondRate->size())
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : inconsistency between strike and second rate states");
		}

		itsFirstStrikeSpread	= arg[7].GetDouble();
		itsSecondStrikeSpread	= arg[8].GetDouble();
		if(itsFirstStrikeSpread<0 || itsSecondStrikeSpread<0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike spread must be >= 0");

		SetAlreadyComputed(true);
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DoubleDigital
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the Spread function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_GP_DoubleDigital::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

	ARM_VectorPtr DoubleDigital = itsModelIR->DoubleDigital(itsModelName,itsEvalTime,
        itsFirstRate,*itsFirstStrikeDown,*itsFirstStrikeUp,itsFirstStrikeSpread,
        itsSecondRate,*itsSecondStrikeDown,*itsSecondStrikeUp,itsSecondStrikeSpread,
		states );

	return ARM_GramFctorArg(DoubleDigital);
}


////////////////////////////////////////////////////
///	Class  : ARM_GP_DoubleDigital
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_GP_DoubleDigital::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	GrabInputs( arg, mod, evalDate, nodes );	

    CC_STL_VECTOR( double ) result;
    result.push_back( evalDate );

    return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( result )),ARM_AdditionalTimeInfoPtr(NULL) );
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

