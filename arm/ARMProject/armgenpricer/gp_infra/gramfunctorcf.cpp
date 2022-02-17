
/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpinfra/gramfunctorcf.h"

/// gpbase
#include "gpbase/checkinputs.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/additionaltimeinfo.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

CC_BEGIN_NAMESPACE( ARM )



/// the corresponding static string
string ARM_CFCallFctor::itsFuncName		= "CFCall";
string ARM_CFGreekFctor::itsFuncName	= "CFGreek";


ARM_VectorPtr ARM_GramFctorArg_ConvertToVector( const ARM_GramFctorArg& arg, size_t vecSize )
{
	if( GFAT_VECTOR_TYPE == arg.GetType() )
		return arg.GetVector();
	else if( GFAT_DOUBLE_TYPE == arg.GetType() )
		return ARM_VectorPtr( new std::vector<double>( vecSize, arg.GetDouble() ) );
	else
       ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" : can convert only double to vector!");
}


////////////////////////////////////////////////////
///	Class  : ARM_CFCallFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CFCallFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 7, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_VECTOR_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[1], GFAT_VECTOR_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[3], GFAT_DOUBLE_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[4], GFAT_VECTOR_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE, ARM_CFCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE, ARM_CFCallFctor::itsFuncName );

	ARM_CFDispatcher::DistributionType distType = (ARM_CFDispatcher::DistributionType) ARM_ArgConv_cfDistType.GetNumber( arg[6].GetString() );

	ARM_VectorPtr result;
	bool vectorialFormula = false;
	size_t vecSize = 0;

	for( size_t i=0; i<4; ++i )
	{
		if(	GFAT_VECTOR_TYPE == arg[i].GetType() )
		{
			vectorialFormula = true;
			vecSize = arg[i].GetVector()->size();
		}
	}

	if ( vectorialFormula )
	{
		ARM_VectorPtr fwd		= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike	= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		ARM_VectorPtr vol		= ARM_GramFctorArg_ConvertToVector( arg[2], vecSize );
		ARM_VectorPtr df		= ARM_GramFctorArg_ConvertToVector( arg[4], vecSize );

		double maturity			= arg[3].GetDouble();
		double sqrtMaturity		= sqrt(maturity);

		double callOrPut		= ARM_ArgConv_CallPut.GetNumber( arg[5].GetString() );

		size_t resultSize		= fwd->size();
		result = ARM_VectorPtr( new std::vector<double>( resultSize ) );
		size_t i;
		
		switch( distType )
		{
			case ARM_CFDispatcher::K_LN_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = BlackSholes_Formula( (*fwd)[i], (*vol)[i]*sqrtMaturity, (*df)[i], (*strike)[i], callOrPut );
				break;
			case ARM_CFDispatcher::K_Normal_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = (*df)[i] * VanillaOption_N( (*fwd)[i], (*vol)[i], (*strike)[i], maturity, callOrPut );
				break;
			default:
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
		}

		return ARM_GramFctorArg(result );//merd
	}
	else
	{
		double fwd			= arg[0].GetDouble();
		double strike		= arg[1].GetDouble();
		double vol			= arg[2].GetDouble();
		double maturity		= arg[3].GetDouble();
		double df			= arg[4].GetDouble();
		double sqrtMaturity	= sqrt(maturity);
		double callOrPut	= ARM_ArgConv_CallPut.GetNumber( arg[5].GetString() );
		double resultDouble = 0.0;

		switch( distType )
		{
			case ARM_CFDispatcher::K_LN_DIST:
				resultDouble = BlackSholes_Formula( fwd, vol * sqrtMaturity, df, strike, callOrPut );
				break;

			case ARM_CFDispatcher::K_Normal_DIST:
				resultDouble = df * VanillaOption_N( fwd, vol, strike, maturity, callOrPut );
				break;

			default:
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
		}

		return ARM_GramFctorArg(resultDouble );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CFCallFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CFCallFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo(ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
////////////////////// ARM_CFGreekFctor
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CFGreekFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CFGreekFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 8, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_VECTOR_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[1], GFAT_VECTOR_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[2], GFAT_VECTOR_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[3], GFAT_DOUBLE_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[4], GFAT_VECTOR_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_STRING_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_STRING_TYPE, ARM_CFGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[7], GFAT_STRING_TYPE, ARM_CFGreekFctor::itsFuncName );

	ARM_VectorPtr result;
	bool vectorialFormula = false;
	size_t vecSize = 0;

	for( size_t i=0; i<4; ++i )
	{
		if(	GFAT_VECTOR_TYPE == arg[i].GetType() )
		{
			vectorialFormula = true;
			vecSize = arg[i].GetVector()->size();
		}
	}


	ARM_CFDispatcher::DistributionType distType = (ARM_CFDispatcher::DistributionType) ARM_ArgConv_cfDistType.GetNumber( arg[6].GetString() );
	ARM_CFDispatcher::GreekType greekType		= (ARM_CFDispatcher::GreekType) ARM_ArgConv_cfGreekType.GetNumber( arg[7].GetString() );

	if ( vectorialFormula )
	{
		ARM_VectorPtr fwd		= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike	= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		ARM_VectorPtr vol		= ARM_GramFctorArg_ConvertToVector( arg[2], vecSize );
		ARM_VectorPtr df		= ARM_GramFctorArg_ConvertToVector( arg[4], vecSize );

		double maturity			= arg[3].GetDouble();
		double sqrtMaturity		= sqrt(maturity);
		double callOrPut		= ARM_ArgConv_CallPut.GetNumber( arg[5].GetString() );

		size_t resultSize		= fwd->size();
		result = ARM_VectorPtr( new std::vector<double>( resultSize ) );
		size_t i;

		switch( distType )
		{
		case ARM_CFDispatcher::K_LN_DIST:
			switch( greekType )
			{
			case ARM_CFDispatcher::K_Delta:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = BlackSholes_Derivative_1( (*fwd)[i], (*vol)[i]*sqrtMaturity, (*df)[i], (*strike)[i], callOrPut );
				break;

			case ARM_CFDispatcher::K_Vega:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = BlackSholes_Derivative_2( (*fwd)[i], (*vol)[i]*sqrtMaturity, (*df)[i], (*strike)[i], callOrPut );
				break;

			case ARM_CFDispatcher::K_Gamma:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = BlackSholes_Derivative_1_1( (*fwd)[i], (*vol)[i]*sqrtMaturity, (*df)[i], (*strike)[i], callOrPut );
				break;

			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
			break;
		
		
		case ARM_CFDispatcher::K_Normal_DIST:
			switch( greekType )
			{
			case ARM_CFDispatcher::K_Delta:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = (*df)[i] * DeltaVanillaOption_N( (*fwd)[i], (*vol)[i], (*strike)[i], maturity, callOrPut );
				break;

			case ARM_CFDispatcher::K_Vega:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = (*df)[i] * VegaVanillaOption_N( (*fwd)[i], (*vol)[i], (*strike)[i], maturity, callOrPut );
				break;

			case ARM_CFDispatcher::K_Gamma:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = (*df)[i] * GammaVanillaOption_N( (*fwd)[i], (*vol)[i], (*strike)[i], maturity, callOrPut );
				break;

			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
			break;
	
		default:
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
		}

		return ARM_GramFctorArg(result );
	}
	else
	{
		double fwd			= arg[0].GetDouble();
		double strike		= arg[1].GetDouble();
		double vol			= arg[2].GetDouble();
		double maturity		= arg[3].GetDouble();
		double df			= arg[4].GetDouble();
		double sqrtMaturity	= sqrt(maturity);
		double callOrPut	= ARM_ArgConv_CallPut.GetNumber( arg[5].GetString() );

		double resultDouble = 0.0;


		switch( distType )
		{
		case ARM_CFDispatcher::K_LN_DIST:
			switch( greekType )
			{
			case ARM_CFDispatcher::K_Delta:
				resultDouble = BlackSholes_Derivative_1( fwd, vol *sqrtMaturity, df, strike, callOrPut );
				break;

			case ARM_CFDispatcher::K_Vega:
				resultDouble = BlackSholes_Derivative_2( fwd, vol *sqrtMaturity, df, strike, callOrPut );
				break;

			case ARM_CFDispatcher::K_Gamma:
				resultDouble = BlackSholes_Derivative_1_1( fwd, vol *sqrtMaturity, df, strike, callOrPut );
				break;

			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
			break;
		
		
		case ARM_CFDispatcher::K_Normal_DIST:
			switch( greekType )
			{
			case ARM_CFDispatcher::K_Delta:
				resultDouble= df * DeltaVanillaOption_N( fwd, vol, strike, maturity, callOrPut );
				break;

			case ARM_CFDispatcher::K_Vega:
				resultDouble= df * VegaVanillaOption_N( fwd, vol, strike, maturity, callOrPut );
				break;

			case ARM_CFDispatcher::K_Gamma:
				resultDouble= df * GammaVanillaOption_N( fwd, vol, strike, maturity, callOrPut );;
				break;

			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
		}

		return ARM_GramFctorArg(resultDouble );
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_CFGreekFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CFGreekFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo(ARM_GP_VectorPtr(new ARM_GP_T_Vector<double>( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

