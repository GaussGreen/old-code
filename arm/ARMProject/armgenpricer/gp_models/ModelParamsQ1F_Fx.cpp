/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file ModelParamsQ1F_Fx.cpp
 *
 *  \brief Q model 1 factor FX version
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpvector.h"

/// gpmodels
#include "gpmodels/ModelParamsQ1F_Fx.h"


/// gpinfra
#include "gpinfra/curvemodelparam.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsQ1F_Fx::ARM_ModelParamsQ1F_Fx( const ARM_ModelParamVector& params, ARM_ZeroCurvePtr domCurve, ARM_ZeroCurvePtr fgnCurve, double spot  )
: ARM_ModelParams_Fx(domCurve, fgnCurve, spot),
ARM_ModelParamsQ1F(params)
{
	SetVolatilityType(ARM_ModelParamType::QVol);

	if( params.size() != 3 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected 2 model parameters: Q Vol, Q param and mean reversion!" );

	ValidateModelParams();
}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQ1F_Fx
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_ModelParamsQ1F_Fx::ValidateModelParams() const
{	
    /// checks that the cev model contains a q volatility
	if( !DoesModelParamExist(ARM_ModelParamType::QVol) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: requires a Q parameter!");

	/// checks that the cev model contains a q parameter
	if( !DoesModelParamExist(ARM_ModelParamType::QParameter) )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": QModel1FParam: requires a Q vol!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F_Fx
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsQ1F_Fx::ARM_ModelParamsQ1F_Fx( const ARM_ModelParamsQ1F_Fx& rhs )
: ARM_ModelParams_Fx(rhs),
ARM_ModelParamsQ1F(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F_Fx
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsQ1F_Fx::~ARM_ModelParamsQ1F_Fx()
{}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQ1F_Fx
///	Routine: toString
///	Returns: string
///	Action : Display the contents
////////////////////////////////////////////////////
string ARM_ModelParamsQ1F_Fx::toString(const string& indent,const string& nextIndent) const
{
	string str;

	str += ARM_ModelParams_Fx::toString(indent,nextIndent);
	str += ARM_ModelParamsQ1F::toString(indent,nextIndent);

	return str;
}


ARM_ModelParamVector ARM_Q1FModelParams_FxBuilder::CreateAndValidateModelParams( const ARM_ModelParamVector& params )
{
	/// validates
	if( params.size() != 2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, " expected 2 model parameters: Q Vol and Q param but NO Mean Reversion!" );

	/// copies and clone the first two model parameters
	ARM_ModelParamVector result(3);
	result[0]= params[0]? static_cast<ARM_ModelParam*>( params[0]->Clone() ) : NULL;
	result[1]= params[1]? static_cast<ARM_ModelParam*>( params[1]->Clone() ) : NULL;
	
	/// creates a zero mean reversion
	std::vector<double>& values = new std::vector<double>(1,0.0);
	std::vector<double>& times  = new std::vector<double>(1,0.0);
	result[3]= new ARM_CurveModelParam( ARM_ModelParamType::MeanReversion, values, times );
	delete values;
	delete times;

	return result;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

