/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file TargetFuncHW.cpp
 *
 *  \brief 
 *
 *	\author  E.M Ezzine; E.Benhamou
 *	\version 1.0
 *	\date March 2004
 */



/// this header should come first as it includes
/// some preprocessor constants!

/// gpmodels
#include "gpmodels/TargetFuncHW.h"
#include "gpmodels/modelparamshw1f.h"

///gpbase
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "gpcalib/vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
//////////      VarDiffFunc       //////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : VarDiffFunc
///	Routine: Constructor
///	Returns: 
///	Action : Builds the object!
////////////////////////////////////////////////////
VarDiffFunc::VarDiffFunc(ARM_ModelParams* params,
    ARM_PricingModel* model,
    ARM_ModelParamType::ParamNb paramType,
    DbleBinaryFunctor* Func )
:	itsParams(params),
    itsParamType( paramType), 
    itsBinaryFunc(Func),
    itsPricingModel(model),
    itsDerivative(NULL),
    itsTarget(0), 
    itsIndex(0),
    itsTime(0.0)
{
    itsDerivative = new ARM_NumDerivativeDbleToDbleFunctor(*this);
}


////////////////////////////////////////////////////
///	Class  : FunctionToSolve
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
VarDiffFunc::~VarDiffFunc()
{
    delete itsDerivative;
	itsDerivative = NULL;
}

////////////////////////////////////////////////////
///	Class  : VarDiffFunc
///	Routine: operator() (const double& x)
///	Returns: 
///	Action : compute the value of the target function!
////////////////////////////////////////////////////

double VarDiffFunc::operator() ( double x ) const
{
	/// compute the binary func with the first argument being the variance and the second one being
    /// the target variance

    /// set the value x on at the break point time itsIndex
    /// on the model param with param type = itsParamType
    ARM_ModelParams::iterator foundVol = itsParams->SetModelParamValue( itsParamType, itsIndex, x , itsTime);
	ARM_CurveModelParam* volModelParam  = dynamic_cast<ARM_CurveModelParam*>(*foundVol);
#ifdef __GP_STRICT_VALIDATION
	if( !volModelParam )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
			"Could not cast to an ARM_CurveModelParam*" );
#endif

    double lag = volModelParam->GetCurve()->GetAbscisse(itsIndex);
    double var = ((ARM_CalibParamsHW1FExt*)itsParams)->SquaredIntegral(0,lag,
                        &volModelParam->GetCurve()->GetAbscisses(),
						&volModelParam->GetCurve()->GetOrdinates());
    return (*itsBinaryFunc)(var,itsTarget);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

