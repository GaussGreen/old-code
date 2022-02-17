/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file AnalyticIRModel.cpp
 *
 *  \brief base class for all the analytic interest rates model
 *	\author  E. Benhamou, O. Croissant
 *	\version 1.0
 *	\date October 2004
 */

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpmodels/AnalyticIRModel.h"

/// gpbase
#include "gpbase/surface.h"
#include "gpbase/surfacetypedef.h"

///gpcalib
#include "gpcalib/calibmethod.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/surfacemodelparam.h"

/// kernel
#include <inst/portfolio.h>
#include <inst/swaption.h>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: CalibSecurityTypeString
///	Returns: static table to convert a CalibSecurityType to its
///				corresponding string
////////////////////////////////////////////////////
string ARM_AnalyticIRModel::CalibSecurityTypeString[] =
{
	"ARM_SWAPTION", 
	"ARM_CAP"
};


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ThrowError
///	Returns: 
///	Action : function to factorise the throw of error
////////////////////////////////////////////////////

void ARM_AnalyticIRModel::ThrowError( const string& msg )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": " + msg );
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ValidatePricingStates
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_AnalyticIRModel::ValidatePricingStates( const ARM_PricingStatesPtr& states ) const
{
#ifdef __GP_STRICT_VALIDATION
	if( states != ARM_PricingStatesPtr(NULL) && 1 != states->size() )
	{
		CC_Ostringstream os;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": states is neither null or of size 1!" );
	}
#endif	
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ARM_VectorPtr 
///	Returns: 
///	Action : function to compute the discount factor from the interest rate curve
///				there is no stochasticity here!
////////////////////////////////////////////////////

ARM_VectorPtr ARM_AnalyticIRModel::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
#ifdef __GP_STRICT_VALIDATION
	ValidatePricingStates(states);

	if( evalTime > maturityTime )
	{
		CC_Ostringstream os;
		os << "Trying to price a discountFactor in the past.\n" 
			<< ARM_USERNAME  << ": please advise\n";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
#endif

	/// does not use modelStates as this is useless for this model
	/// does not use curveName!
	double fwdDf = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
		/ GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
	return new ARM_GP_Vector(1,fwdDf);
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ARM_AnalyticIRModel
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_AnalyticIRModel::ARM_AnalyticIRModel(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, ARM_DensityFunctor* densityFct)
:	ARM_PricingModelIR(zc, &params, densityFct )
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ARM_AnalyticIRModel
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_AnalyticIRModel::ARM_AnalyticIRModel(const ARM_AnalyticIRModel& rhs )
:	ARM_PricingModelIR( rhs) 
{}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ARM_AnalyticIRModel
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_AnalyticIRModel& ARM_AnalyticIRModel::operator =(const ARM_AnalyticIRModel& rhs )
{
	if( this!= &rhs )
		ARM_PricingModelIR::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: Init
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_AnalyticIRModel::Init(
	const string& payModelName, 
	const ARM_TimeInfoPtrVector& timeInfos )
{	
	return ARM_PricingStatesPtr( new ARM_PricingStates(1,0,1) ); 
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: GetMatchingTenor
///	Returns: double
///	Action : computes  the matching tenor
////////////////////////////////////////////////////
double ARM_AnalyticIRModel::GetMatchingTenor( double tenorInit )
{
	double nbYear	= int( tenorInit);
	double tenor	= tenorInit - nbYear;

	const double SevenDaysLag	= 7.0/365.0;
	const double FifteenDaysLag = 15.0/365.0;
	const double OneYear		= 1.0;
	const double SixMonths		= 0.5;
	const double ThreeMonths	= 0.25;
	const double TwoMonths		= 2.0/12.0;
	const double OneMonth		= 1.0/12.0;
	
	/// test 1Y, 6m, 3m ,2m, 1m
	if( fabs( tenor-OneYear) < SevenDaysLag+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
		return nbYear+OneYear;
	else if( fabs( tenor-SixMonths) < SevenDaysLag+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
		return nbYear+SixMonths;
	else if( fabs( tenor-ThreeMonths) < SevenDaysLag +ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE)
		return nbYear+ThreeMonths;
	else if( fabs( tenor-TwoMonths) < SevenDaysLag+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
		return nbYear+TwoMonths;
	else if( fabs( tenor-OneMonth) < SevenDaysLag+ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
		return nbYear+OneMonth;
	else if( fabs( tenor ) < FifteenDaysLag + ARM_NumericConstants::ARM_STD_DOUBLE_TOLERANCE )
		return nbYear;
	else
        return tenorInit;
    ////Fix Fix
	   //throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " unexpected tenor!");
}

////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: GetInstType
///	Returns: CalibSecurityType
///	Action : get the type of an instrument
////////////////////////////////////////////////////

ARM_AnalyticIRModel::CalibSecurityType ARM_AnalyticIRModel::GetInstType( ARM_Security* security )
{
	if( ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(security) )
		return ARM_SWAPTION_TYPE;
	else if(ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(security) )
		return ARM_CAP_TYPE;
	else
	   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " unexpected security: permitted is swaption and capfloor!");
}



////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: ComputeTenorAndExpiry
///	Returns: void
///	Action : from a security, computes the corresponding expriy and tenor
////////////////////////////////////////////////////
void ARM_AnalyticIRModel::ComputeTenorAndExpiry( ARM_Security* security, double& expiry, double& tenor, CalibSecurityType& type, double asOfDate )
{
	double rawTenor;
	if( ARM_Swaption* swaption = dynamic_cast<ARM_Swaption*>(security) )
	{
		expiry		= swaption->GetExpiryDate().GetJulian() - asOfDate;
		rawTenor	= (swaption->GetFloatLeg()->GetEndDateNA().GetJulian()-swaption->GetFloatLeg()->GetStartDateNA().GetJulian())/K_YEAR_LEN;
		tenor		= ARM_AnalyticIRModel::GetMatchingTenor( rawTenor );
		type		= ARM_SWAPTION_TYPE;
	}
	else if(ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>(security) )
	{
		expiry		= capFloor->GetResetDates()->Elt(capFloor->GetResetDates()->size()-1)-asOfDate;
		rawTenor	= capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();
		tenor		= ARM_AnalyticIRModel::GetMatchingTenor(rawTenor);
		type		= ARM_CAP_TYPE;
	}
	else
	   throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, " unexpected security: permitted is swaption and capfloor!");
}


////////////////////////////////////////////////////
///	Class  : ARM_AnalyticIRModel
///	Routine: AdviseBreakPointTimes
///	Returns: void
///	Action  : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_AnalyticIRModel::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_SurfaceModelParam* modelParam = dynamic_cast<ARM_SurfaceModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_SurfaceModelParam!");
	
	size_t portSize = portfolio->GetSize();
	ARM_GP_Vector expiries;
	expiries.reserve(portSize);
	ARM_GP_Vector tenors;
	tenors.reserve(portSize);
    double asOfDate = GetAsOfDate().GetJulian();
	double expiry, tenor;
	ARM_Security* security	= NULL;

	if( portSize )
	{
		CalibSecurityType type,
			firstType = ARM_AnalyticIRModel::GetInstType(portfolio->GetAsset(0));

		for( size_t i=0; i<portSize; ++i )
		{
			ARM_AnalyticIRModel::ComputeTenorAndExpiry( portfolio->GetAsset(i), expiry, tenor, type, asOfDate );
			if( type != firstType )
		       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "only the same type of instrument allowed!, first type =" +
			   ARM_AnalyticIRModel::CalibSecurityTypeString[firstType] + " while found " + ARM_AnalyticIRModel::CalibSecurityTypeString[type]);
			expiries.push_back(expiry);
			tenors.push_back(tenor);
		}
	}

	tenors.sort();
	tenors.unique();
	expiries.sort();
	expiries.unique();
	ARM_GP_T_Matrix<double> values(expiries.size(),tenors.size(), modelParam->GetValue(expiries[0],tenors[0]));

	ARM_InterpolType type = ARM_InterpolationType::linear_column_extrapoleCst;
	ARM_SurfaceWithInterpol* initialSurface;
	if( initialSurface = dynamic_cast<ARM_SurfaceWithInterpol*>(modelParam->GetSurface()) )
		type = initialSurface->GetInterpolType();

	ARM_SurfaceWithInterpol* surface = new ARM_SurfaceWithInterpol( expiries, tenors, values, type);
	modelParam->SetSurface( surface );
}

////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_AnalyticIRModel::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}




////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : cDuration
///	Returns : double
///	Action  : compute actuarial duration like Summit
////////////////////////////////////////////////////
double ARM_AnalyticIRModel::cDuration(double Rate, double Term, int F, double t) const
{
    double duration =1/Rate * ( 1 - 1/pow(1+Rate*t, Term*F) );

    return(duration);
}


/////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : cConvexity
///	Returns : double
///	Action  : compute actuarial convexity like Summit
/////////////////////////////////////////////////////
double ARM_AnalyticIRModel::cConvexity(double Rate, double Term, int F, double t) const
{
    double convexity = 1 + Rate*t * (1 + Term*F);
    convexity = convexity / pow(1+Rate*t, 1+Term*F) - 1;
    convexity *= -2.0 / (Rate*Rate);

    return(convexity);
}


/////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : DecompRate
///	Returns : double
///	Action  : compute actuarial decoupounded Rate
/////////////////////////////////////////////////////
double ARM_AnalyticIRModel::DecompRate(double Rate, int Freq) const
{
    double decomprate = Freq *(pow(1.+Rate,1./Freq)-1.);

    return(decomprate);
}


/////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : DecompRate
///	Returns : double
///	Action  : compute actuarial decoupounded Rate
/////////////////////////////////////////////////////
double ARM_AnalyticIRModel::InverseDecompRate(double Rate, int Freq) const
{
    double compRate = pow((1.+Rate/Freq),Freq)-1.;

    return(compRate);
}


/////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : ForwardDecap
///	Returns : double
///	Action  : compute vol adjusted Decap Adjusted Fwd
/////////////////////////////////////////////////////
double ARM_AnalyticIRModel::ForwardDecap( double ForwardF, double Spread, 
			       double VolF, 
			       double Freq, 
			       double Time ) const
{
	double forwardDecap;

    if (Freq != 1.0)
       forwardDecap = Freq*(pow((1+ForwardF+Spread),(1/Freq))-1) * 
                      exp(.5*pow(VolF*ForwardF,2)*
			  ((1-Freq)/Freq*pow(1+ForwardF+Spread,(1-2*Freq)/Freq))*Time/
			  (Freq*(pow(1+ForwardF+Spread,(1-Freq)/Freq)-1)));
    else
       forwardDecap = ForwardF+Spread;
   
    return(forwardDecap);
}


/////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routine : VolDecap
///	Returns : double
///	Action  : compute adjusted Vol for Decap Fwd
/////////////////////////////////////////////////////
double ARM_AnalyticIRModel::VolDecap( double ForwardF, 
			   double VolF, 
			   double Freq ) const
{
	double volDecap; 
   
    volDecap = VolF* ForwardF* 
               pow(1+ForwardF, (1-Freq)/Freq)/ (Freq*(pow(1+ForwardF, 1/Freq)-1));

    return(volDecap);
}


////////////////////////////////////////////////////
///	Class   : ARM_AnalyticIRModel
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_AnalyticIRModel::VanillaSpreadOptionLet(const string& curveName,
														double evalTime,
														int callPut,
														double startTime,
														double endTime,
														double resetTime,
														double payTime,
														double payPeriod,
														double notional,
														double coeffLong,
														double coeffShort,
														const ARM_GP_Vector& strikes,
														double swapLongFloatStartTime,
														double swapLongFloatEndTime,
														const ARM_GP_Vector& swapLongFixPayTimes,
														const ARM_GP_Vector& swapLongFixPayPeriods,
														double swapShortFloatStartTime,
														double swapShortFloatEndTime,
														const ARM_GP_Vector& swapShortFixPayTimes,
														const ARM_GP_Vector& swapShortFixPayPeriods,
														const ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSpreadOptionLet : unimplemented function for ARM_AnalyticIRModel Model!");
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

