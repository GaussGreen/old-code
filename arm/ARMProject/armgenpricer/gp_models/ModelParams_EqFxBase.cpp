/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParams_EqFxBase.cpp
 *
 *  \brief
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2006
 */


#include "gpmodels/ModelParams_EqFxBase.h"

#include "gpinfra/pricingfunctionequity.h"

//#include "crv/zerocurv.h"
#include "ccy/currency.h"

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelParams_EqFxBase::ARM_ModelParams_EqFxBase( const ARM_ModelParams_EqFxBase& rhs )
:	itsSpot(rhs.itsSpot)
{
}
////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
string ARM_ModelParams_EqFxBase::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Eq Fx Modeln";
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelParams_Eq::ARM_ModelParams_Eq( const ARM_ZeroCurvePtr& zeroCurve, double spot )
:	ARM_ModelParams_EqFxBase(spot), 
itsZeroCurve(zeroCurve)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParams_Eq::Forward( double MaturityTime ) const
{
	double df			= itsZeroCurve->DiscountPrice(MaturityTime/K_YEAR_LEN);
	return GetSpot()/df;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
string ARM_ModelParams_Eq::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Eq Model with\n";
	os << indent << " -spot               : " << GetSpot()  << "\n";
	os << indent << " -curve ccy : " << itsZeroCurve->GetCurrencyUnit()->GetCcyName() << "\n";
	return os.str();
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ModelParams_Fx::ARM_ModelParams_Fx( 
	const ARM_ZeroCurvePtr& domCurve, 
	const ARM_ZeroCurvePtr& forCurve, 
	double spot)
	:
	ARM_ModelParams_EqFxBase(spot), 
	itsDomCurve(domCurve), 
	itsForCurve(forCurve)
{
	itsFXCal = ARM_PricingFunctionEquity::GetSettlementCalendar(&*domCurve,&*forCurve);
	itsSpotDays = ARM_PricingFunctionEquity::GetSettlementGap(&*domCurve,&*forCurve);
	ARM_Date spotDate = itsDomCurve->GetAsOfDate();
	spotDate.NextBusinessDay(itsSpotDays, const_cast<char*>(itsFXCal.c_str()));
	itsDfRatioAtSpotDate = itsDomCurve->DiscountPrice(spotDate)/itsForCurve->DiscountPrice(spotDate);
}

	////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: double
///	Action : 
////////////////////////////////////////////////////
double ARM_ModelParams_Fx::Forward( double MaturityTime ) const
{
	/// compute the fwdDate
	double domDfT = itsDomCurve->DiscountPrice(MaturityTime/K_YEAR_LEN);
	double forDfT = itsForCurve->DiscountPrice(MaturityTime/K_YEAR_LEN);
	return GetSpot()*forDfT/domDfT*itsDfRatioAtSpotDate;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: 
///	Returns: double
///	Action : static version of forward() to avoid building an object
////////////////////////////////////////////////////
double ARM_ModelParams_Fx::Forward(ARM_ZeroCurve* domCurve,
   ARM_ZeroCurve* forCurve,
   double spot,
   double MaturityTime)
{
	double domDfT = domCurve->DiscountPrice(MaturityTime/K_YEAR_LEN);
	double forDfT = forCurve->DiscountPrice(MaturityTime/K_YEAR_LEN);
	ARM_Date spotDate = ComputeSettlementDate(domCurve,forCurve,domCurve->GetAsOfDate());
	double domDf0 = domCurve->DiscountPrice(spotDate);
	double forDf0 = forCurve->DiscountPrice(spotDate);
	return spot*forDfT/domDfT*domDf0/forDf0;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: toString
///	Returns: string
///	Action : 
////////////////////////////////////////////////////
ARM_Date ARM_ModelParams_Fx::ComputeSettlementDate(ARM_ZeroCurve* domCurve,
	ARM_ZeroCurve* forCurve,
	ARM_Date& date)
{
	string setCal = ARM_PricingFunctionEquity::GetSettlementCalendar(domCurve,forCurve);
	double setGap = ARM_PricingFunctionEquity::GetSettlementGap(domCurve,forCurve);
	ARM_Date spotDate(date);
	spotDate.NextBusinessDay(setGap, const_cast<char*>(setCal.c_str()));
	return spotDate;
}

////////////////////////////////////////////////////
///	Class  : ARM_ModelParams_Fx
///	Routine: toString
///	Returns: string
///	Action : 
////////////////////////////////////////////////////
string ARM_ModelParams_Fx::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << " Fx Model with\n";
	os << indent << " -spot               : " << GetSpot()  << "\n";
	os << indent << " -fx calendar        : " << itsFXCal << "\n";
	os << indent << " -fx spot days       : " << itsSpotDays << "\n";
	os << indent << " -domDFSpot/forDFSpot: " << itsDfRatioAtSpotDate << "\n";
	os << indent << " -domestic curve ccy : " << itsDomCurve->GetCurrencyUnit()->GetCcyName() << "\n";
	os << indent << " -foreign curve ccy  : " << itsForCurve->GetCurrencyUnit()->GetCcyName() << "\n\n";
	return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
