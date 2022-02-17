/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsMF.cpp
 *
 *  \brief Markov Functional Model Params
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date August 2005
 */

/// gpmodels
#include "gpmodels/ModelParamsMF.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsMF::ARM_ModelParamsMF( const ARM_ModelParamsMF& rhs )
: ARM_ModelParamsHW1FStd(rhs)
{
	
	ARM_Curve* curve = ( (ARM_CurveModelParam&)GetModelParam(GetVolatilityType()) ).GetCurve() ;
	
	if ( !dynamic_cast<ARM_StepUpRightOpenCstExtrapol<double,double>*>(curve->GetInterpolator()) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsMF : vol curve is supposed to be step up right !" );

}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsMF::ARM_ModelParamsMF( const ARM_ModelParamVector& params )
: ARM_ModelParamsHW1FStd(params)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ModelParamsMF::~ARM_ModelParamsMF()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: void
///	Action : SetVolUpToT1AndFreezeGlobVarUpToT2
////////////////////////////////////////////////////
void ARM_ModelParamsMF::SetVolUpToT1AndFreezeGlobVarUpToT2(double T1, double T2, double volUpToT1, bool& varSqueeze)
{
	if ( T1 >= T2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsMF::SetVolUpToT1AndFreezeGlobVarUpToT2 : T1>T2 !" );

	/// Compute initial variance up to T2
	/// That's what we want to leave unchanged
	double globalVar = StateLocalVariance(0.0, T2, T2) ;

	/// Get references to modify sigmaValues
	ARM_CurveModelParam& modelParam = (ARM_CurveModelParam&)GetModelParam(GetVolatilityType());
	ARM_Curve* curve = modelParam.GetCurve() ;
	std::vector<double>& sigmaTimes = curve->GetAbscisses();
	std::vector<double>& sigmaValue = curve->GetOrdinates();


	/// Find index for T1 and T2
	/// could be optimized
	size_t size = sigmaTimes.size();
	
	size_t i;
	int idx1 = -1;
	for (i=0; i<size; i++)
	{
		if( fabs(sigmaTimes[i] - T1) < 0.5 )
		{
			idx1 = i;
			break;
		}
	}

	if (idx1 == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "SetVolUpToT1AndFreezeGlobVarUpToT2 : T1 not found" );

	int idx2 = -1;
	for (i=idx1; i<size; i++)
	{
		if( fabs(sigmaTimes[i] - T2) < 0.5 )
		{
			idx2 = i;
			break;
		}
	}

	if (idx2 == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "SetVolUpToT1AndFreezeGlobVarUpToT2 : T2 not found" );

	
	/// Get Lower Bound
	double lowerBound = modelParam.GetLowerBound()->Elt(idx2);

	/// Get Mean Reversion
	double lambda = GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);

	///
	double t1 = T1 / K_YEAR_LEN;
	double t2 = T2 / K_YEAR_LEN;

	
	/// compute vol from T1 to T2
	double volFromT1toT2, sqrVol, varFromT1toT2Factor, varUpToT1, varUpToT1Factor;
	
	if (fabs(lambda)>1e-6)
	{
		double expT1minusT2 = exp( 2.*lambda*(t1-t2));
		double expminusT2   = exp(-2.*lambda*t2);

		varUpToT1Factor		= (expT1minusT2 - expminusT2) / (2.0 * lambda);
		varUpToT1			=  volUpToT1 * volUpToT1 * varUpToT1Factor;
		varFromT1toT2Factor = (1.0 - expT1minusT2) / (2.0 * lambda);
		
	}
	else
	{
		varUpToT1Factor		= t1;
		varUpToT1			= volUpToT1 * volUpToT1 * t1;
		varFromT1toT2Factor	= t2 - t1 ;
	}

	sqrVol = (globalVar - varUpToT1) / varFromT1toT2Factor;

	
	if ( sqrVol < lowerBound * lowerBound)
	{
		/// set vol T1->T2 to lower bound
		volFromT1toT2 = lowerBound;

		/// set vol 0->T1 so that global var is preserved ( --> input T1 only will be mispriced )
		varUpToT1 = globalVar - varFromT1toT2Factor * volFromT1toT2 * volFromT1toT2  ;
		volUpToT1 = sqrt(varUpToT1/varUpToT1Factor);
		
		varSqueeze = true ;
	}
	else
	{
		volFromT1toT2 = sqrt(sqrVol);
		varSqueeze = false;
	}
	

	/// set vol up to T1 (curve is assumed to be step up right)
	for (i=0; i<=idx1; i++)
		sigmaValue[i] = volUpToT1;

	/// set vol from T1 to T2
	for (i=idx1+1; i<=idx2; i++)
		sigmaValue[i] = volFromT1toT2;

}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: void
///	Action : SetVolFromT1toT2AndFreezeGlobVarUpToT2
////////////////////////////////////////////////////
void ARM_ModelParamsMF::SetVolFromT1toT2AndFreezeGlobVarUpToT2(double T1, double T2, double volFromT1toT2, bool& varSqueeze)
{
	if ( T1 >= T2 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_ModelParamsMF::SetVolUpToT1AndFreezeGlobVarUpToT2 : T1>T2 !" );

	/// Compute initial variance up to T2
	/// That's what we want to leave unchanged
	double globalVar = StateLocalVariance(0.0, T2, T2) ;

	/// Get references to modify sigmaValues
	ARM_CurveModelParam& modelParam = (ARM_CurveModelParam&)GetModelParam(GetVolatilityType());
	ARM_Curve* curve = modelParam.GetCurve() ;
	std::vector<double>& sigmaTimes = curve->GetAbscisses();
	std::vector<double>& sigmaValue = curve->GetOrdinates();


	/// Find index for T1 and T2
	/// could be optimized
	size_t size = sigmaTimes.size();
	
	size_t i;
	int idx1 = -1;
	for (i=0; i<size; i++)
	{
		if( fabs(sigmaTimes[i] - T1) < 0.5 )
		{
			idx1 = i;
			break;
		}
	}

	if (idx1 == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "SetVolFromT1toT2AndFreezeGlobVarUpToT2 : T1 not found" );

	int idx2 = -1;
	for (i=idx1; i<size; i++)
	{
		if( fabs(sigmaTimes[i] - T2) < 0.5 )
		{
			idx2 = i;
			break;
		}
	}

	if (idx2 == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "SetVolFromT1toT2AndFreezeGlobVarUpToT2 : T2 not found" );

	
	/// Get Lower Bound
	double lowerBound = modelParam.GetLowerBound()->Elt(idx1);

	/// Get Mean Reversion
	double lambda = GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);

	///
	double t1 = T1 / K_YEAR_LEN;
	double t2 = T2 / K_YEAR_LEN;

	
	/// compute vol from T1 to T2
	double volUpToT1, sqrVol, varFromT1toT2Factor, varFromT1toT2, varUpToT1Factor;
	
	if (fabs(lambda)>1e-6)
	{
		double expT1minusT2 = exp( 2.*lambda*(t1-t2));
		double expminusT2   = exp(-2.*lambda*t2);
		varUpToT1Factor		= (expT1minusT2 - expminusT2) / (2.0 * lambda);
		varFromT1toT2Factor = (1.0 - expT1minusT2) / (2.0 * lambda);
		varFromT1toT2		= volFromT1toT2 * volFromT1toT2 * varFromT1toT2Factor;
		
	}
	else
	{	varUpToT1Factor		= t1;
		varFromT1toT2Factor	= t2 - t1 ;
		varFromT1toT2		= volFromT1toT2 * volFromT1toT2 * (t2 - t1);
	}

	sqrVol = (globalVar - varFromT1toT2) / varUpToT1Factor;

	
	if ( sqrVol < lowerBound * lowerBound)
	{
		/// set vol 0->T1 to lower bound
		volUpToT1 = lowerBound;
		varSqueeze = true ;
	}
	else
	{
		volUpToT1 = sqrt(sqrVol);
		varSqueeze = false;
	}
	

	/// set vol up to T1 (curve is assumed to be step up right)
	for (i=0; i<=idx1; i++)
		sigmaValue[i] = volUpToT1;

	/// set vol from T1 to T2
	for (i=idx1+1; i<=idx2; i++)
		sigmaValue[i] = volFromT1toT2;


}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: toString
///	Returns: string 
///	Action : computes integrated underlying process variance
////////////////////////////////////////////////////
string ARM_ModelParamsMF::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "ARM_ModelParamsMF\n";
    os << "----------------------\n\n";
    os << ARM_ModelParams::toString();
    return os.str();
}


CC_END_NAMESPACE()