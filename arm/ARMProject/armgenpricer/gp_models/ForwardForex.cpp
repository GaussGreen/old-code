/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ForwardForex.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date July 2004
 */


#include "gpmodels/forwardforex.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/utilityport.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/pricingcontext.h"

/// gpcalib
#include "gpcalib/calibmethod.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ForwardForex::ARM_ForwardForex( const ARM_Forex& forex, 
	const ARM_ZeroCurvePtr& domCurve,
    const ARM_ZeroCurvePtr&foreignCurve )
:   
    ARM_AnalyticPricingModel(domCurve),
	ARM_PricingFunctionEquity(), 
	itsForex(forex),
	itsForeignCurve( foreignCurve )
{
	CheckCurves();
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_ForwardForex::ARM_ForwardForex(const ARM_ForwardForex& rhs)
:
	ARM_AnalyticPricingModel(rhs), 
	ARM_PricingFunctionEquity(rhs), 
	itsForex( rhs.itsForex),
	itsForeignCurve( rhs.itsForeignCurve )
{}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_ForwardForex::~ARM_ForwardForex()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_ForwardForex& ARM_ForwardForex::operator=(const ARM_ForwardForex& rhs)
{
	if(this != &rhs)
	{
		ARM_AnalyticPricingModel::operator=(rhs), 
		ARM_PricingFunctionEquity::operator=(rhs), 
		itsForex		= rhs.itsForex;
        itsForeignCurve = rhs.itsForeignCurve;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: CheckCurves
///	Returns: 
///	Action : Check if curves are consistent with
///          forex object. If not may try to curves
////////////////////////////////////////////////////

void ARM_ForwardForex::CheckCurves()
{
    string domesticForexCcy(const_cast< ARM_Forex& >(itsForex).GetMoneyCurrency()->GetCcyName());
    string domesticCurveCcy(GetZeroCurve()->GetCurrencyUnit()->GetCcyName());

    string foreignForexCcy(const_cast< ARM_Forex& >(itsForex).GetMainCurrency()->GetCcyName());
    string foreignCurveCcy(itsForeignCurve->GetCurrencyUnit()->GetCcyName());

    if(domesticForexCcy != domesticCurveCcy || foreignForexCcy != foreignCurveCcy)
    {
        if(domesticForexCcy == foreignCurveCcy && foreignForexCcy == domesticCurveCcy)
        {
            /// Invert curves
            ARM_ZeroCurvePtr domesticCurve = GetZeroCurve();
            SetZeroCurve(itsForeignCurve);
            itsForeignCurve = domesticCurve;
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : currencies inconstency between Forex and Yield curves");
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: SetCurves
///	Returns: 
///	Action : Set zero curves
////////////////////////////////////////////////////
void ARM_ForwardForex::SetCurves(const ARM_ZeroCurvePtr& domesticCurve, const ARM_ZeroCurvePtr& foreignCurve)
{
    SetZeroCurve(domesticCurve);
    itsForeignCurve = foreignCurve;
    CheckCurves();
}



////////////////////////////////////////////////////
///	Class   : ARM_ForwardForex
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ForwardForex::Clone() const
{
    return new ARM_ForwardForex(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardForex
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_ForwardForex::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_ForwardMargin
///	Routines: toString
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
string ARM_ForwardForex::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Forward Forex Model\n";
    os << indent << "-------------------\n\n";

    /// Pas de toString() pour ARM_Forex !
    os << indent << "Forex :\n";
    os << indent << "Foreign currency :" << string(const_cast< ARM_Forex& >(itsForex).GetMainCurrency()->GetCcyName());
    os << indent << "Domestic currency  :" << string(const_cast< ARM_Forex& >(itsForex).GetMoneyCurrency()->GetCcyName());
    os << indent << "Market Spot       :" << itsForex.GetMarketPrice();

    return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_FXSpot
///	Routine: DiscountFactor
///	Returns: 
///	Action : DiscountFactor
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardForex::DiscountFactor( 
	const string& curveName,
    double evalTime, 
	double maturityTime,
    const ARM_PricingStatesPtr& states) const
{
	/// deterministic discount factor on the zero coupon curve
    double domesticDF   = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
	size_t stateSize;
	if( states != ARM_PricingStatesPtr(NULL) )
		stateSize = states->size();
	else 
		stateSize = 1;
    return ARM_VectorPtr( new std::vector<double>( stateSize, domesticDF ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Forward
///	Returns: double
///	Action : Compute the deterministic forward forex
///          equals to its value at asOfDate
////////////////////////////////////////////////////

double ARM_ForwardForex::Forward(
    double evalTime,
	double maturityTime ) const
{
    double domesticDF   = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);
    double foreignDF    = itsForeignCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
    return itsForex.GetMarketPrice()*foreignDF/domesticDF;
}



////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Forward
///	Returns: ARM_VectorPtr
///	Action : Compute the deterministic forward forex
///          equals to its value at asOfDate
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardForex::Forward(
	const string& curveName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
	size_t stateSize;
	if( states != ARM_PricingStatesPtr(NULL) && states->size())
		stateSize = states->size();
	else 
		stateSize = 1;
	double forward = Forward(evalTime,expiryTime);
    return ARM_VectorPtr( new std::vector<double>( stateSize, forward ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: Forward
///	Returns: ARM_VectorPtr
///	Action : Compute the deterministic forward forex
///          equals to its value at asOfDate
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ForwardForex::CallVectorial(
	const string& curveName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikePerState,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	ARM_VectorPtr call = Forward( curveName, evalTime, expiryTime, settlementTime, payTime, states );
    double domesticDF  = GetZeroCurve()->DiscountPrice(expiryTime/K_YEAR_LEN);

	if( callPut == K_CALL )
		for( size_t i=0; i<call->size(); ++i )
			(*call)[i] = domesticDF * CC_Max( (*call)[i]-strikePerState[i], 0.0 );
	else
		for( size_t i=0; i<call->size(); ++i )
			(*call)[i] = domesticDF * CC_Max( strikePerState[i]-(*call)[i], 0.0 );
	return call;
}



////////////////////////////////////////////////////
///	Class   : ARM_ForwardForex
///	Routines: GetSettlementCalendar
///	Returns : string
///	Action  : get the settlement calendar
////////////////////////////////////////////////////
string ARM_ForwardForex::GetSettlementCalendar(const string& modelName) const
{
	char FXCal[7];
	strcpy(FXCal, GetZeroCurve()->GetCurrencyUnit()->GetCcyName());
	strcat(FXCal, itsForeignCurve->GetCurrencyUnit()->GetCcyName());
	return string(FXCal);
}


////////////////////////////////////////////////////
///	Class   : ARM_ForwardForex
///	Routines: GetSettlementGap
///	Returns : double
///	Action  : get the settlement gap
////////////////////////////////////////////////////
double ARM_ForwardForex::GetSettlementGap(const string& modelName) const
{
	double domGap = GetZeroCurve()->GetCurrencyUnit()->GetSpotDays();
	double forGap = itsForeignCurve->GetCurrencyUnit()->GetSpotDays();
	return CC_Max(domGap,forGap);
}


////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////

int ARM_ForwardForex::GetType() const
{
	return MT_NON_STOCHASTIC_MODEL;
}

////////////////////////////////////////////////////
///	Class  : ARM_ForwardForex
///	Routine: ImpliedVol
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
double ARM_ForwardForex::ImpliedVol(const ARM_VanillaArg& arg) const
{
    return DefaultImpliedVol(arg);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

