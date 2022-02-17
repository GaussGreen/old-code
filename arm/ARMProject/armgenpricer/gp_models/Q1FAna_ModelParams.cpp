/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1FAna_ModelParams.cpp
 *
 *  \brief
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/Q1FAna_ModelParams.h"
#include "gpmodels/QModelAnalytics.h"

/// gpbase
#include "gpbase/numericconstant.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"

/// gpnumlib
#include "gpnumlib/normalinvcum.h"
#include "gpnumlib/gaussiananalytics.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_ModelParams
///	Routine: Constructor
///	Returns: 
///	Action : Constructor initialising parameters vector
////////////////////////////////////////////////////
ARM_Q1FAna_ModelParams::ARM_Q1FAna_ModelParams( const ARM_ModelParamVector& params )
: ARM_ModelParams(params)
{
	ValidateModelParams();
}



////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_ModelParams
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Q1F model
////////////////////////////////////////////////////
void ARM_Q1FAna_ModelParams::ValidateModelParams() const
{	
    /// checks that the q model is of size 1 since the current model is with cst q!
	if(DoesModelParamExist(ARM_ModelParamType::QParameter) && GetModelParam(ARM_ModelParamType::QParameter).size() != 1)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": currently, only QModel with one q supported!");

	/// check positivity of the q
	if( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0) < -K_NEW_DOUBLE_TOL )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": currently, q should be positive!");

	/// check the vol type
	if( DoesModelParamExist(ARM_ModelParamType::QVol) )
	{
		if( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QVol)).GetCurve()->GetOrdinate(0) < -K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": q volatility should be positive!");
		itsVolType = ARM_ModelParamType::QVol;
	}
	else if( DoesModelParamExist(ARM_ModelParamType::NVol) )
	{
		if( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::NVol)).GetCurve()->GetOrdinate(0) < -K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": normal vol should be positive!");
		itsVolType = ARM_ModelParamType::NVol;
	}
	else if( DoesModelParamExist(ARM_ModelParamType::LNVol) )
	{
		if( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LNVol)).GetCurve()->GetOrdinate(0) < -K_NEW_DOUBLE_TOL )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": lognormal volatility should be positive!");
		itsVolType = ARM_ModelParamType::LNVol;
	}
	else throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unknown vol type, supported is qVol, NVol, LNVol!");
}




////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_ModelParams
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Q1FAna_ModelParams::ARM_Q1FAna_ModelParams( const ARM_Q1FAna_ModelParams& rhs )
: ARM_ModelParams(rhs), itsVolType( rhs.itsVolType)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_ModelParams
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_Q1FAna_ModelParams::~ARM_Q1FAna_ModelParams()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Q1FAna_ModelParams
///	Routine: operator =
///	Returns: itself
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_Q1FAna_ModelParams& ARM_Q1FAna_ModelParams::operator=(const ARM_Q1FAna_ModelParams& rhs)
{
	if(this != &rhs)
	{
		ARM_ModelParams::operator=(rhs);
		itsVolType = rhs.itsVolType;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_ModelParams
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_Q1FAna_ModelParams::Clone() const
{
	return new ARM_Q1FAna_ModelParams(*this);
}



////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_ModelParams
///	Routines: CalibrateModelParams
///	Returns : void
///	Action  : calibrate the Q vol given a Q, N or LN vol
////////////////////////////////////////////////////
void ARM_Q1FAna_ModelParams::CalibrateModelParams( double forward, double maturity, int callPut )
{
	/// the calibration of the qvol is first done closed forms and refined with a newton rhapson procedure
	/// if the newton fails, then takes the closed forms solution!
	CalibrateQVol(forward,maturity,callPut);
}




////////////////////////////////////////////////////
///	Class   : ARM_Q1FAna_ModelParams
///	Routines: CalibrateQVol
///	Returns : void
///	Action  : calibrate the Q vol given a N or LN vol
////////////////////////////////////////////////////
void ARM_Q1FAna_ModelParams::CalibrateQVol( double forward, double maturity, int callPut )
{
	if( itsVolType != ARM_ModelParamType::QVol )
	{
		double qParameter	= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		double qVolValue;
		
		/// eliminitate trivial case:
		if( itsVolType == ARM_ModelParamType::LNVol && qParameter>1-K_NEW_DOUBLE_TOL)
			qVolValue = ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LNVol)).GetCurve()->GetOrdinate(0);
		else if( itsVolType == ARM_ModelParamType::NVol && qParameter<K_NEW_DOUBLE_TOL)
			qVolValue = ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::NVol)).GetCurve()->GetOrdinate(0)/forward;
		else
		{
			double sqrtMaturity = sqrt(maturity);
			double ATMPrice,QVolValueInitial;
			
			switch(itsVolType)
			{
			case ARM_ModelParamType::LNVol:
				ATMPrice	= forward*(2*ARM_GaussianAnalytics::cdfNormal(((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::LNVol)).GetCurve()->GetOrdinate(0)*sqrtMaturity*0.5)-1.0);
				break;
			case ARM_ModelParamType::NVol:
				ATMPrice	= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::NVol)).GetCurve()->GetOrdinate(0) *sqrtMaturity*ARM_NumericConstants::ARM_INVSQRT2PI;
				break;
			default:
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": unknonw vol type!!");
			}
			
			/// 1) get closed forms initial value 
			if( fabs(qParameter) > K_NEW_DOUBLE_TOL )
				QVolValueInitial = ARM_NormalInvCum::Inverse_erf_Moro(0.5*(ATMPrice*qParameter/forward+1.0))
				*2.0/(qParameter*sqrtMaturity);
			else
				QVolValueInitial = ATMPrice*ARM_NumericConstants::ARM_SQRT_2_PI/sqrtMaturity/forward;
			
			size_t maxIter		= 25;
			double epsilon		= 0.00001;
			double precision	= 1.0e-8;
			double invPrecision = 1.0/precision;
			
			/// 2) manual Newton Rhapson
			double targetFunc,targetFuncUp,dTargetFunc,step;
			bool foundSolution	= false;
			qVolValue			= QVolValueInitial;

			for( size_t i=0; i<maxIter; ++i ) 
			{
				targetFunc	= QModelAnalytics::BSQFunction( forward, forward, qVolValue, maturity, qParameter, 1.0, callPut, forward)-ATMPrice;
		
				if( fabs(targetFunc) < precision)
				{
					foundSolution = true;
					break;
				}
				else
				{
					targetFuncUp= QModelAnalytics::BSQFunction( forward, forward, qVolValue+epsilon, maturity, qParameter, 1.0, callPut, forward )-ATMPrice;
					dTargetFunc	= (targetFuncUp-targetFunc)*invPrecision;

					/// test the division by zero
					if( fabs(dTargetFunc) < K_NEW_DOUBLE_TOL )
					{
						break;
					}
					else
					{
						step = epsilon*targetFunc/(targetFuncUp-targetFunc);
						qVolValue -= step;

					}
				}
			}

			/// take the initial value in case of no solution!
			if( !foundSolution )
				qVolValue = QVolValueInitial;
		}

		/// set the corresponding model parameter
		std::vector<double>& values = new std::vector<double>(1,qVolValue);
		std::vector<double>& breakPointTimes = new std::vector<double>(1,0.0);
		ARM_ModelParam* qVol= new ARM_CurveModelParam( ARM_ModelParamType::QVol,
			values, breakPointTimes );
		SetModelParam( qVol );
	}
}





CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

