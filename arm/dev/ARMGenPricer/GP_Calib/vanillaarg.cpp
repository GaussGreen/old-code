/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillaarg.cpp
 *
 *  \brief vanilla args are object to do the conversion from
 *		kernel object to gp
 *	\author  E.M Ezzine E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalib/vanillaarg.h"

/// gpbase
#include "gpbase/utilityport.h"
#include "gpbase/surface.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/surfacemodelparam.h"
#include "gpinfra/surfacelistmodelparam.h"

/// ARM Kernel
#include <glob/expt.h>


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: ThrowErrorOnEmptyVector
///	Returns: 
///	Action : only throw an exception in strict validation mode
////////////////////////////////////////////////////
void ARM_VanillaArg::ThrowErrorOnEmptyVector( const string& vectorName, const ARM_VectorVector& vec ) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(!vec.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + vectorName + " is a vector of object with null size!" );
#endif
}

////////////////////////////////////////////////////
///	Struct : VanillaArgument
///	Routine: Virtual destructor
///	Returns: 
///	Action : does nothing but necessary for linking
////////////////////////////////////////////////////

ARM_VanillaArg::~ARM_VanillaArg()
{}

////////////////////////////////////////////////////
///	Struct : VanillaArgument
///	Routine: ThrowErrorOnNullObject
///	Returns: 
///	Action : only throw an exception in strict validation mode
////////////////////////////////////////////////////
void ARM_VanillaArg::ThrowErrorOnNullObject( const string& objName, ARM_Object* obj ) const
{
#if defined(__GP_STRICT_VALIDATION)
	if( NULL == obj )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + objName + " is NULL!" );
#endif
}


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: default constructor, copy constructor,
///          assigment, destructor, clone
///	Returns: 
///	Action :
////////////////////////////////////////////////////
// FIXMEFRED: mig.vc8 (22/05/2007 18:04:50):missing return type
void ARM_VanillaArg::CopyNoCleanUp(const ARM_VanillaArg& rhs)
{
    itsCallPut      =   rhs.itsCallPut;	
	itsCurveName    =   rhs.itsCurveName; 
	itsEvalTime     =   rhs.itsEvalTime;
    itsMktPrice     =   rhs.itsMktPrice; 
    itsIndex        =   rhs.itsIndex;   
	itsExpiry		=	rhs.itsExpiry;
}   
 
////////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: operator=
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_VanillaArg& ARM_VanillaArg::operator=(const ARM_VanillaArg& rhs)
{
	if( this != & rhs )
	{
		ARM_RootObject::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}  

////////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: Copy constuctor
///	Returns: 
///	Action :
////////////////////////////////////////////////////
ARM_VanillaArg::ARM_VanillaArg(const ARM_VanillaArg& rhs)
:	
	ARM_RootObject(rhs),
    itsCallPut(K_CALL),
    itsCurveName(string()),
    itsEvalTime(0.0),
    itsMktPrice(0.0),
    itsIndex(0.0)
{
	CopyNoCleanUp(rhs);
}             


////////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: Compute the derivative with respect to a model param
///	Returns: 
///	Action :
////////////////////////////////////////////////////

double ARM_VanillaArg::Derivative(ARM_PricingModel* model, 
                     const ARM_ModelParam& modelParam, 
                     size_t xIndex,
					 size_t yIndex,
					 size_t zIndex, 
                     int factorNb,
                     ARM_MktTargetType targetType) const
{
	if( model->HasClosedFormsDerivatives( modelParam.GetType(), factorNb ) )
	{
		return model->PartialDerivative( modelParam, xIndex, factorNb, *this,targetType );
	}
	else
	{
		return NumericDerivative( model, &const_cast<ARM_ModelParam&>(modelParam), xIndex,yIndex,zIndex,factorNb ,targetType);
	}
}


///////////////////////////////////////////////////
///	Struct : ARM_VanillaArg
///	Routine: Compute the numerical derivative with respect to a model param
///	Returns: 
///	Action :
////////////////////////////////////////////////////

double ARM_VanillaArg::NumericDerivative(ARM_PricingModel* model, 
                                         ARM_ModelParam* modelParam,
                                         size_t xIndex,
										 size_t yIndex,
										 size_t zIndex,
                                         int factorNb,
                                         ARM_MktTargetType targetType) const
{
	double delta = 0.0;
	if ( dynamic_cast<ARM_CurveModelParam*>(modelParam))
	{
		double initParamValue= modelParam->GetValueAtPoint( xIndex );
		double epsilon		 = 0.001 * CC_Max<double>( initParamValue, 1.0 );

		/// Find the equivalnet time lag and index in modelParam of the PricingModel
		double timelag = modelParam->GetTimeAtPoint( xIndex );

		/// centered scheme [f(x+e)-f(x-e)]/2e
		model->GetModelParams()->GetModelParam( modelParam->GetType(), factorNb ).SetValue( timelag, initParamValue+epsilon );
		double priceUp = Price(model);
		model->GetModelParams()->GetModelParam( modelParam->GetType(), factorNb ).SetValue( timelag, initParamValue-epsilon );
		double priceDown = Price(model);
		delta =  (priceUp-priceDown)/(2.0*epsilon);
	}
	else if (ARM_SurfaceModelParam* surfaceCalibParam = dynamic_cast<ARM_SurfaceModelParam*>(modelParam))
	{
		double initParamValue= modelParam->GetValueAtPoint( xIndex, yIndex);
		double epsilon		 = 0.001 * CC_Max<double>( initParamValue, 1.0 );

		/// Find the equivalnet time lag and index in modelParam of the PricingModel
		double timelag = modelParam->GetTimeAtPoint( xIndex);
		double tenor   = modelParam->GetTenorAtPoint( yIndex );

		/// centered scheme [f(x+e)-f(x-e)]/2e
		model->GetModelParams()->GetModelParam( modelParam->GetType(), factorNb ).SetValue( timelag, tenor,initParamValue+epsilon );
		double priceUp = Price(model);
		model->GetModelParams()->GetModelParam( modelParam->GetType(), factorNb ).SetValue( timelag, tenor,initParamValue-epsilon );
		double priceDown = Price(model);
		delta =  (priceUp-priceDown)/(2.0*epsilon);
	
	}
	else if( ARM_SurfaceListModelParam* surfaceListCalibParam = dynamic_cast<ARM_SurfaceListModelParam*>(modelParam))
		throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT, "NumericDerivative : Unimplemented method for ARM_SurfaceListModelParam !" );


	return delta;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

