/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQGM1F.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

#include "gpmodels/modelparamsqgm1f.h"
#include "gpmodels/ouprocess.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"

/// gpcalib
#include "gpcalib/bootstrap1d.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/gpmatrix.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_ModelParamsQGM1F::ARM_ModelParamsQGM1F( const ARM_ModelParamVector& params )
: ARM_ModelParams(params)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsQGM1F::ARM_ModelParamsQGM1F( const ARM_ModelParamsQGM1F& rhs )
: ARM_ModelParams(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsQGM1F::~ARM_ModelParamsQGM1F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ModelParamsQGM1F& ARM_ModelParamsQGM1F::operator=(const ARM_ModelParamsQGM1F& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		/// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsQGM1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsQGM1F::Clone() const
{
	return new ARM_ModelParamsQGM1F(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: PreProcessing
///	Returns: void
///	Action : 
////////////////////////////////////////////////////
void ARM_ModelParamsQGM1F::PreProcessing(ARM_ModelFitter& modelFitter,int factorNb)
{
    /// we put a typeid to implement in a sense a double dispatcher...
    /// we did not use a dynamic to avoid throwing exception of type std::bad_cast
    if( typeid(modelFitter) == typeid(ARM_Bootstrap1D) && 
        modelFitter.GetPortfolio()->GetAsset(0)->GetName() == ARM_SWAPTION &&
		(modelFitter.GetCalibParam())->GetType() == ARM_ModelParamType::Skew) 
    {
        modelFitter.SetCalibDirection(CalibDirection_Backward);
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: StateLocalDrift
///	Returns: value of the drift
///	Action : Relative drift of the state variable
///          from a to b>=a
////////////////////////////////////////////////////
double ARM_ModelParamsQGM1F::StateLocalDrift(double a,double b) const
{
	double MRSValue=((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	return ARM_OUProcess::Drift(a,b,MRSValue);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
///          (alias phi(a,b))
////////////////////////////////////////////////////
double ARM_ModelParamsQGM1F::StateLocalVariance(double a,double b) const
{
    if(b - K_NEW_DOUBLE_TOL <= a)
        return 0.0;

	double MRSValue=((ARM_CurveModelParam&)GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
    const ARM_Curve& sigma = *(((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve());
    return ARM_OUProcess::Variance(a,b,sigma,MRSValue);
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: FwdZcLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of Zc(.,T1)/Zc(.,T2)
///          <=> FwdZcLocalCovariance(a,b,T1,T2,T1,T2)
///           but specialised to save a function call
////////////////////////////////////////////////////
double ARM_ModelParamsQGM1F::FwdZcLocalVariance(double a,double b,double T1,double T2) const
{
    return 0.0;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsQGM1F
///	Routine: FwdZcLocalcovariance
///	Returns: value of the covariance
///	Action : Covariance in [a,b] of Zc(.,T1)/Zc(.,U1)
///          and Zc(.,T2)/Zc(.,U2)
////////////////////////////////////////////////////
double ARM_ModelParamsQGM1F::FwdZcLocalCovariance(double a,double b,double T1,double U1,double T2,double U2,std::vector<double>& vars) const
{
    return 0.0;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

