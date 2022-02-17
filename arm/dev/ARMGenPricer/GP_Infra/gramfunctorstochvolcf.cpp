/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file gramfunctorstochvolcf.cpp
 *
 *  \brief gramfunctorcf is the file for closed form related grammar functor
 *	the functor design allows to have a context and a unified command like
 *	interface
 *	\author  O. Croissant
 *	\version 1.0
 *	\date June 2005
 */

#include "gpbase/removeidentifiedwarning.h"

#include "gpinfra/gramfunctorStochVolcf.h"

/// gpbase
#include "gpbase/checkinputs.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/gramfunctorconv.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/gramfunctorargcheck.h"
#include "gpinfra/argconvdefault.h"

/// gpclosedforms
#include "gpclosedforms/stochasticvol_ln.h"
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/sabrimpliedvol.h"
#include "gpclosedforms/sabrvanilla.h"
#include "gpclosedforms/extended_sabr_interface.h"
#include "gpclosedforms/stochasticvol_ln_interface.h"

CC_BEGIN_NAMESPACE( ARM )



/// the corresponding static string
string ARM_CF_StochVolCallFctor::itsFuncName		= "CFStochVolCall";
string ARM_CF_StochVolGreekFctor::itsFuncName		= "CFStochVolCallGreek";

////////////////////////////////////////////////////
///	Class  : ARM_CF_StochVolCallFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CF_StochVolCallFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 12, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[0],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[1],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[2],  GFAT_DOUBLE_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[3],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[4],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[5],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[6],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[7],  GFAT_DOUBLE_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[8],  GFAT_VECTOR_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[9],  GFAT_STRING_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	GPAF_CheckArgType( arg[10],  GFAT_DOUBLE_TYPE, ARM_CF_StochVolCallFctor::itsFuncName );
	

	ARM_CF_StochVolDispatcher::DistributionType distType = (ARM_CF_StochVolDispatcher::DistributionType) ARM_ArgConv_cfMepiDistType.GetNumber( arg[11].GetString() );

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
		ARM_VectorPtr fwd			= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike		= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		double maturity				= arg[2].GetDouble();
		ARM_VectorPtr zero_rate		= ARM_GramFctorArg_ConvertToVector( arg[3], vecSize );
		ARM_VectorPtr Vol_Initial	= ARM_GramFctorArg_ConvertToVector( arg[4], vecSize );
		ARM_VectorPtr VolDrift		= ARM_GramFctorArg_ConvertToVector( arg[5], vecSize );
		ARM_VectorPtr VolOfVol		= ARM_GramFctorArg_ConvertToVector( arg[6], vecSize );
		double avg_time				= arg[7].GetDouble();	
		ARM_VectorPtr reset				= ARM_GramFctorArg_ConvertToVector( arg[8], vecSize );	
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[9].GetString() );
		double hermite_size			= arg[10].GetDouble();		/// in fact integer coded as a double
 

		size_t resultSize			= fwd->size();
		result = ARM_VectorPtr( new ARM_GP_Vector( resultSize ) );
		size_t i;
		
		switch( distType )
		{
			case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset((*fwd)[i],(*strike)[i],maturity,(*zero_rate)[i], 
					(*Vol_Initial)[i], (*VolDrift)[i], (*VolOfVol)[i],
					avg_time,(*reset)[i],(int) callOrPut,(int) hermite_size);
				
				break;
			case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset((*fwd)[i],(*strike)[i],maturity,(*zero_rate)[i], 
					(*Vol_Initial)[i], (*VolDrift)[i], (*VolOfVol)[i],
					avg_time,(*reset)[i],(int) callOrPut,(int) hermite_size);
				
				break;
			case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset((*fwd)[i],(*strike)[i],maturity,(*zero_rate)[i], 
					(*Vol_Initial)[i], (*VolDrift)[i], (*VolOfVol)[i],
					avg_time,(*reset)[i],(int) callOrPut,(int) hermite_size);
				
				break;
			case ARM_CF_StochVolDispatcher::K_SABR_DIST:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = 
					Export_SABR_VanillaOption(
							   (*fwd)[i],			///FORWARD
							   (*strike)[i],		///STRIKE
							   maturity,			///MATURITY
							   (*Vol_Initial)[i],	///ALPHA
							   1.0,					///BETA 
							   0.0,					///RHO 
							   (*VolOfVol)[i],		///NU
							   (int) callOrPut,		///CALLORPUT
							   ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0,  /// Flag
							   (int) hermite_size	///NBSTEPS
							   )*exp(-(*zero_rate)[i]*maturity);
				break;
			default:
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
		}

		return ARM_GramFctorArg(result );
	}
	else
	{
		double fwd					= arg[0].GetDouble();
		double strike				= arg[1].GetDouble();
		double maturity				= arg[2].GetDouble();
		double zero_rate			= arg[3].GetDouble();
		double Vol_Initial			= arg[4].GetDouble();
		double VolDrift				= arg[5].GetDouble();
		double VolOfVol				= arg[6].GetDouble();
		double avg_time				= arg[7].GetDouble();
		double reset				= arg[8].GetDouble();
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[9].GetString() );
		double hermite_size			= arg[10].GetDouble();		/// in fact integer coded as a double
		

 
		double resultDouble = 0.0;
		
		switch( distType )
		{
		case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
			
			resultDouble = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(fwd,strike,maturity,zero_rate, 
				Vol_Initial,VolDrift,VolOfVol,
				avg_time,reset,(int) callOrPut
				,(int) hermite_size);
			break;
		case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST:
			
			resultDouble = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(fwd,strike,maturity,zero_rate, 
				Vol_Initial,VolDrift,VolOfVol,
				avg_time,reset,(int) callOrPut
				,(int) hermite_size);
			break;
		case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST:
			
			resultDouble = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(fwd,strike,maturity,zero_rate, 
				Vol_Initial,VolDrift,VolOfVol,
				avg_time,reset,(int) callOrPut
				,(int) hermite_size);
			break;
		case ARM_CF_StochVolDispatcher::K_SABR_DIST:
			
			resultDouble =
					Export_SABR_VanillaOption(
							   fwd,			///FORWARD
							   strike,		///STRIKE
							   maturity,			///MATURITY
							   Vol_Initial,	///ALPHA
							   1.0,					///BETA 
							   0.0,					///RHO 
							   VolOfVol,		///NU
							   (int) callOrPut,		///CALLORPUT
							   ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0,  /// Flag
							   (int) hermite_size	///NBSTEPS
							   )*exp(-zero_rate*maturity);
				break;


			default:
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
		}

		return ARM_GramFctorArg(resultDouble );
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_CF_StochVolCallFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CF_StochVolCallFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}



////////////////////////////////////////////////////
////////////////////// ARM_CF_StochVolGreekFctor
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_CF_StochVolGreekFctor
///	Routine: operator()
///	Returns: ARM_GramFctorArg
///	Action : computes the closed form function and return 
///				the corresponding values
////////////////////////////////////////////////////
ARM_GramFctorArg ARM_CF_StochVolGreekFctor::operator()( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod,
	double evalDate, const ARM_PricingStatesPtr& states, vector< ARM_ExpNodePtr >& nodes )
{
	/// checking of the size and type
	GPAF_CheckArgSize( arg, 13, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[0], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[1], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[2], GFAT_DOUBLE_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[3], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[4], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[5], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[6], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[7], GFAT_DOUBLE_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[8], GFAT_VECTOR_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[9], GFAT_STRING_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	GPAF_CheckArgType( arg[10], GFAT_DOUBLE_TYPE, ARM_CF_StochVolGreekFctor::itsFuncName );
	

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


	ARM_CF_StochVolDispatcher::DistributionType distType = (ARM_CF_StochVolDispatcher::DistributionType) ARM_ArgConv_cfMepiDistType.GetNumber( arg[11].GetString() );
	ARM_CF_StochVolDispatcher::GreekType greekType		= (ARM_CF_StochVolDispatcher::GreekType) ARM_ArgConv_cfGreekType.GetNumber( arg[12].GetString() );

	if ( vectorialFormula )
	{
		ARM_VectorPtr fwd			= ARM_GramFctorArg_ConvertToVector( arg[0], vecSize );
		ARM_VectorPtr strike		= ARM_GramFctorArg_ConvertToVector( arg[1], vecSize );
		double maturity				= arg[2].GetDouble();
		ARM_VectorPtr zero_rate		= ARM_GramFctorArg_ConvertToVector( arg[3], vecSize );
		ARM_VectorPtr Vol_Initial			= ARM_GramFctorArg_ConvertToVector( arg[4], vecSize );
		ARM_VectorPtr VolDrift		= ARM_GramFctorArg_ConvertToVector( arg[5], vecSize );
		ARM_VectorPtr VolOfVol		= ARM_GramFctorArg_ConvertToVector( arg[6], vecSize );
		double avg_time				= arg[7].GetDouble();
		ARM_VectorPtr reset				= ARM_GramFctorArg_ConvertToVector( arg[8], vecSize );	
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[9].GetString() );
		double hermite_size			= arg[10].GetDouble();		/// in fact integer coded as a double
		size_t resultSize		= fwd->size();
		result = ARM_VectorPtr( new ARM_GP_Vector( resultSize ) );
		size_t i;

		switch( distType )
		{
		case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
			switch( greekType )
			{
			case ARM_CF_StochVolDispatcher::K_Delta:
				for( i=0; i<resultSize; ++i )
					(*result)[i] = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(0,(*fwd)[i],
					(*strike)[i],
					maturity,
					(*zero_rate)[i],
					(*Vol_Initial)[i],
					(*VolDrift)[i],
					(*VolOfVol)[i],
					avg_time,
					(*reset)[i],
					(int) callOrPut,
					(int) hermite_size
					);
				break;
				
			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
			break;
			case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST:
				switch( greekType )
				{
				case ARM_CF_StochVolDispatcher::K_Delta:
					for( i=0; i<resultSize; ++i )
						(*result)[i] = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(0,(*fwd)[i],
						(*strike)[i],
						maturity,
						(*zero_rate)[i],
						(*Vol_Initial)[i],
						(*VolDrift)[i],
						(*VolOfVol)[i],
						avg_time,
						(*reset)[i],
						(int) callOrPut,
						(int) hermite_size
						);
					break;
					
				default:
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
				break;
				case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST:
					switch( greekType )
					{
					case ARM_CF_StochVolDispatcher::K_Delta:
						for( i=0; i<resultSize; ++i )
							(*result)[i] = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(0,(*fwd)[i],
							(*strike)[i],
							maturity,
							(*zero_rate)[i],
							(*Vol_Initial)[i],
							(*VolDrift)[i],
							(*VolOfVol)[i],
							avg_time,
							(*reset)[i],
							(int) callOrPut,
							(int) hermite_size
							);
						break;
						
					default:
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
					}
					break;
					
					case ARM_CF_StochVolDispatcher::K_SABR_DIST:
						switch( greekType )
						{
						case ARM_CF_StochVolDispatcher::K_Delta:
							for( i=0; i<resultSize; ++i )
								
								(*result)[i] = 
								Export_SABR_VanillaOption(0,
								(*fwd)[i],			///FORWARD
								(*strike)[i],		///STRIKE
								maturity,			///MATURITY
								(*Vol_Initial)[i],	///ALPHA
								1.0,					///BETA 
								0.0,					///RHO 
								(*VolOfVol)[i],		///NU
								(int) callOrPut,		///CALLORPUT
								ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0,  /// Flag
								(int) hermite_size	///NBSTEPS
								)*exp(-(*zero_rate)[i]*maturity);
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
		double fwd					= arg[0].GetDouble();
		double strike				= arg[1].GetDouble();
		double maturity				= arg[2].GetDouble();
		double zero_rate			= arg[3].GetDouble();
		double Vol_Initial			= arg[4].GetDouble();
		double VolDrift				= arg[5].GetDouble();
		double VolOfVol				= arg[6].GetDouble();
		double avg_time				= arg[7].GetDouble();	
		double reset				= arg[8].GetDouble();		
		double callOrPut			= ARM_ArgConv_CallPut.GetNumber( arg[9].GetString() );
		double hermite_size			= arg[10].GetDouble();		/// in fact integer coded as a double
		
		double resultDouble = 0.0;
		
		
		switch( distType )
		{
		case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_DIST:
			switch( greekType )
			{
			case ARM_CF_StochVolDispatcher::K_Delta:
				resultDouble  = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(0,fwd,
					strike,
					maturity,
					zero_rate,
					Vol_Initial,
					VolDrift,
					VolOfVol,
					avg_time,
					reset,
					(int) callOrPut,
					(int) hermite_size
					);
				break;
				
			default:
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
			}
			break;
			case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_GEOMETRIC_DIST:
				switch( greekType )
				{
				case ARM_CF_StochVolDispatcher::K_Delta:
					resultDouble  = Export_StochasticVol_LN_Geometric_VanillaOption_with_Reset(0,fwd,
						strike,
						maturity,
						zero_rate,
						Vol_Initial,
						VolDrift,
						VolOfVol,
						avg_time,
						reset,
						(int) callOrPut,
						(int) hermite_size
						);
					break;
					
				default:
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
				}
				break;
				case ARM_CF_StochVolDispatcher::K_STOCHASTIC_BLACKSCHOLES_ARITHMETIC_DIST:
					switch( greekType )
					{
					case ARM_CF_StochVolDispatcher::K_Delta:
						resultDouble  = Export_StochasticVol_LN_Arithmetic_VanillaOption_with_Reset(0,fwd,
							strike,
							maturity,
							zero_rate,
							Vol_Initial,
							VolDrift,
							VolOfVol,
							avg_time,
							reset,
							(int) callOrPut,
							(int) hermite_size
							);
						break;
						
					default:
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
					}
					break;
					case ARM_CF_StochVolDispatcher::K_SABR_DIST:
						switch( greekType )
						{
						case ARM_CF_StochVolDispatcher::K_Delta:
							resultDouble = 
								Export_SABR_VanillaOption(0,
								fwd,			///FORWARD
								strike,		///STRIKE
								maturity,			///MATURITY
								Vol_Initial,	///ALPHA
								1.0,					///BETA 
								0.0,					///RHO 
								VolOfVol,		///NU
								(int) callOrPut,		///CALLORPUT
								ARM_CF_SABR_ImplicitVol_Formula::ANALYTICZP0,  /// Flag
								(int) hermite_size	///NBSTEPS
								)*exp(-zero_rate*maturity);
							break;
							
						default:
							ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": greek type unknown");
						}
						break;
						default:
							ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : distribution type not supported");
							
		}
		
		return ARM_GramFctorArg(resultDouble );
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_CF_StochVolGreekFctor
///	Routine: GetUsedTimeLags
///	Returns: ARM_VectorPtr of time lags
///	Action : computes time lag given the input!
////////////////////////////////////////////////////

ARM_NodeInfo ARM_CF_StochVolGreekFctor::GetUsedTimeLags( ARM_GramFctorArgVector& arg, ARM_PricingModel* mod, 
	double evalDate, vector<ARM_ExpNodePtr>& nodes )
{
	return ARM_NodeInfo( ARM_GP_VectorPtr(new ARM_GP_Vector( 0 )),ARM_AdditionalTimeInfoPtr(NULL) );
}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

