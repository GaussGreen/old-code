/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file ModelParamsEqHybrid.cpp
 *  \brief
 *	\author  A. Schauly
 *	\version 1.0
 *	\date March 2005
 *
 */


/// this header comes first as it includes some preprocessor constants!

#include "gpmodels/InflationEquityModelParams.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/curve.h"
#include "gpbase/curveutils.h"
#include "gpbase/gpmatrixtriangular.h"

/// gpinfra
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"
#include "gpbase/comparisonfunctor.h"

#include <cmath>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ModelParamsInflationEquity::ARM_ModelParamsInflationEquity( const ARM_ModelParamVector& params )
:	ARM_ModelParams(params)
{
	ValidateModelParams();
	InitSquaredIntegratedVol();
	PreComputeSquaredIntegratedVol();
}


////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsInflationEquity
///	Routines: Validate
///	Returns :
///	Action  : validate the model params to check that this is compatible with the Equity Hybrid model
////////////////////////////////////////////////////
void ARM_ModelParamsInflationEquity::ValidateModelParams() const
{	
	
	if(!DoesModelParamExist(ARM_ModelParamType::Multiplier))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an InflationEquity Model should have a multiplier!");
	if( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Multiplier)).GetCurve()->size() != 1 )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": multiplier should be of size 1 (constant parameter)!");
	
	if(!DoesModelParamExist(ARM_ModelParamType::Volatility))
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": an InflationEquity Model should have a volatility!");
	if(!((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->size() )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the volatility curve should have at least 1 element!");
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: Copy constructor
///	Returns: 
///	Action : copy the object
////////////////////////////////////////////////////
ARM_ModelParamsInflationEquity::ARM_ModelParamsInflationEquity( const ARM_ModelParamsInflationEquity& rhs )
:	ARM_ModelParams(rhs), itsSquaredIntegratedVol( (ARM_Curve*)rhs.itsSquaredIntegratedVol->Clone() )
{
	InitSquaredIntegratedVol();
	PreComputeSquaredIntegratedVol();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsSFRMRow
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_ModelParams::iterator ARM_ModelParamsInflationEquity::SetModelParamValue( int paramType, 
	size_t i,
	double value, 
	double time, 
	double tenor)
{
    ARM_ModelParams::iterator found = ARM_ModelParams::SetModelParamValue(paramType,i,value, time );

    /// should update the vol values
    if( paramType == ARM_ModelParamType::Volatility )
        PreComputeSquaredIntegratedVol();
    return found;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: oeprator=
///	Returns: 
///	Action : assignment operator 
////////////////////////////////////////////////////
ARM_ModelParamsInflationEquity& ARM_ModelParamsInflationEquity::operator=( const ARM_ModelParamsInflationEquity& rhs )
{
	if( this != &rhs )
	{
		ARM_ModelParams::operator=(rhs);
		itsSquaredIntegratedVol = ARM_CurvePtr( (ARM_Curve* )rhs.itsSquaredIntegratedVol->Clone() );
		PreComputeSquaredIntegratedVol();
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: destructor
///	Returns: 
///	Action : destroy the object
////////////////////////////////////////////////////
ARM_ModelParamsInflationEquity::~ARM_ModelParamsInflationEquity()
{}

////////////////////////////////////////////////////
///	Class   : ARM_ModelParamsInflationEquity
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ModelParamsInflationEquity::Clone() const
{
	return new ARM_ModelParamsInflationEquity(*this);
}

///////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: toString
///	Returns: 
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_ModelParamsInflationEquity::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    os << "\n\n =======> Model Params Inflation Equity <====== \n";
    os << "---------------------------------------------\n\n";
    os << ARM_ModelParams::toString();
	return os.str();
}

///////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: toString
///	Returns: 
///	Action : computes the integrated inf volatility basically
///				\latexonly
///				$\exp \left[ \int_{a}^{b}sgm_{u}^2du\right] $
///				\endlatexonly
////////////////////////////////////////////////////
double ARM_ModelParamsInflationEquity::IntegratedVolatility( double a, double b ) const
{
	/// take real year term maturities!
	a = a/K_YEAR_LEN;
	b = b/K_YEAR_LEN;

	int VolSize = ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->size();
	
	for (size_t i=0; i<VolSize; i++);
		((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinate(i)
			*= ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinate(i);

	double integral = IntegrateStepWise( *((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve(), a, b );
	return integral;
}


////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: StateLocalVariance
///	Returns: value of the variance
///	Action : Variance in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_ModelParamsInflationEquity::StateLocalVariance(double a,double b) const
{
    return IntegratedLocalVarianceFromZero(b)-IntegratedLocalVarianceFromZero(a);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: StateLocalDrift
///	Returns: value of the drift
///	Action : Variance in [a,b] of the state variable
////////////////////////////////////////////////////
double ARM_ModelParamsInflationEquity::StateLocalDrift(double a,double b) const
{
    return -0.5*StateLocalVariance( a, b );

}

double ARM_ModelParamsInflationEquity::getMultiplier() const 
{ 
	return GetModelParam( ARM_ModelParamType::Multiplier).GetValue( 0 ); 
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: InitSquaredIntegratedVol
///	Returns: void
///	Action : initializes the squared integrated vol (size)
////////////////////////////////////////////////////
void ARM_ModelParamsInflationEquity::InitSquaredIntegratedVol()
{
	/// Initializes the squaredIntegratedVol
	ARM_GP_Vector times(( (ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
    ARM_GP_Vector values(times.size());
	itsSquaredIntegratedVol = ARM_CurvePtr(new ARM_Curve(times,values, new ARM_LinInterpCstExtrapolDble ));
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: IntegratedLocalVariance
///	Returns: double
///	Action : integrates sum(0,t) of sigma(u)^2*exp(2*MeanRev*u)du
////////////////////////////////////////////////////
void ARM_ModelParamsInflationEquity::PreComputeSquaredIntegratedVol()
{
	/// Computes the squared integrated Vol
	ARM_GP_Vector timesValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses());
	ARM_GP_Vector volValues( ((ARM_CurveModelParam&) GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates());

	/// Initializes the squaredIntegratedVol
    ARM_GP_Vector values(timesValues.size());
	itsSquaredIntegratedVol = ARM_CurvePtr(new ARM_Curve(timesValues,values, new ARM_LinInterpCstExtrapolDble ));

#if defined(__GP_STRICT_VALIDATION)
	CC_NS(ARM_Check,CheckSameArgSize)(volValues,timesValues,"volValues","timesValues");
#endif

	size_t volSize=volValues.size();
	ARM_GP_Vector squaredVol(volSize);
	double value;

	double time = 0.0;
	for(size_t i=0;i<volSize;++i)
	{
        double nextTime = timesValues[i]/K_YEAR_LEN;
		value = volValues[i]*volValues[i]*(nextTime-time);
		time = nextTime;
#if defined(__GP_STRICT_VALIDATION)
		if(value<0.0)
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "negative vol squared!");
#endif
		if(i>0)
			value+=squaredVol[i-1];
		squaredVol[i]=value;
	}

	itsSquaredIntegratedVol->SetOrdinates(squaredVol);
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParamsInflationEquity
///	Routine: IntegratedLocalVariance
///	Returns: double
///	Action : integrates sum(0,t) of sigma(u)^2*exp(2*MeanRev*u)du
////////////////////////////////////////////////////
double ARM_ModelParamsInflationEquity::IntegratedLocalVarianceFromZero(double t) const
{
    if (t<K_NEW_DOUBLE_TOL)
            return  0.0;

	double vol,Ti,sqrtt;
	size_t Size = GetVolCurve()->size();

    ARM_GP_Vector vect(GetVolCurve()->GetAbscisses());
	
	/// to get exactly the index lower than the index, we need to substract one!
	int index = lower_boundPosWithPrecision(vect,t,K_FRM_TOL)-1;
	
	/// t<t0?
	/// if so Sum(0,t) = 
	///		-if( MeanRev==0 ) Sigma(0)^2*t
	///		-else 1.0/(2.0*MeanRev) * Sigma(0)^2 * (exp(2.0*MeanRev*t)-1.0)
    if(index == -1)
	{
	    vol = (GetVolCurve()->GetOrdinates())[0];
	    sqrtt = vol*vol*t/K_YEAR_LEN;
	}
	/// t>t0 : Use SquaredIntegratedVol from 0 to Ti and adds the sum from Ti to t
	else
	{
        double sqrTi = (itsSquaredIntegratedVol->GetOrdinates())[index];
		Ti			 = (GetVolCurve()->GetAbscisses())[index];

		/// t==Ti up to K_NEW_DOUBLE_TOL precision?
		if( fabs(t-Ti)<K_NEW_DOUBLE_TOL )
			sqrtt = sqrTi;	    

		/// otherwise
		/// Get the next one...
		/// and compute the remaining portion with same formula as above!
		/// sqrtt = sqrTi + sqrTitot (sqrtTi)
		else
		{
            vol		= GetVolCurve()->Interpolate(t);
			sqrtt	= vol*vol*(t-Ti)/K_YEAR_LEN + sqrTi;
		}
	}
	return sqrtt;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
