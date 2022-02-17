/*
 * Copyright (c) CDC IXIS CM July 2005 Paris
 *
 *	\file stripper.cpp
 *
 *  \brief
 *
 *	\author  Y. KHLIF
 *	\version 1.0
 *	\date July 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpvector.h"

#include "gpcalib/stripper.h"
#include "gpbase/curve.h"

#include <iomanip>	

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_STRIPPER::ARM_STRIPPER(const ARM_STRIPPER& rhs):
	itsAsofDate ( rhs.itsAsofDate ),
	itsDatestrip( ARM_DateStripPtr(static_cast<ARM_DateStrip*>(rhs.itsDatestrip->Clone())) ),
	itsSwapRates ( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsSwapRates->Clone())) ),
	itsMaturities( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsMaturities->Clone())) ),  
	itsZcRates( ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsZcRates->Clone())) )	
{
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_STRIPPER& ARM_STRIPPER::operator=(const ARM_STRIPPER& rhs)
{
	if( this != & rhs )
	{
		itsAsofDate		= rhs.itsAsofDate;
		itsDatestrip	= ARM_DateStripPtr(static_cast<ARM_DateStrip*>(rhs.itsDatestrip->Clone()));
		itsSwapRates	= ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsSwapRates->Clone()));
		itsMaturities	= ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsMaturities->Clone()));  
		itsZcRates		= ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(rhs.itsZcRates->Clone()));
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_STRIPPER::~ARM_STRIPPER()
{
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Object* ARM_STRIPPER::Clone() const
{
	return new ARM_STRIPPER(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_STRIPPER::ARM_STRIPPER (double asOfDate, ARM_DateStripPtr datestrip)
:	itsAsofDate(asOfDate),	
	itsDatestrip (datestrip), 
	itsSwapRates( ARM_GP_VectorPtr(new ARM_GP_Vector(datestrip->GetResetDates()->size()))),
	itsMaturities (ARM_GP_VectorPtr(new ARM_GP_Vector(datestrip->GetResetDates()->size()+1))),
	itsZcRates(ARM_GP_VectorPtr(new ARM_GP_Vector(datestrip->GetResetDates()->size()+1)))
{
	std::vector<double>* startDates = datestrip->GetFlowStartDates();
	std::vector<double>* endDates = datestrip->GetFlowEndDates();
	int size = startDates->size(); 
	(*itsMaturities)[0] = (*startDates)[0]-itsAsofDate;
	for (int i=0; i<size; i++)
	{
		(*itsMaturities)[i+1] = (*endDates)[i]-itsAsofDate;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_STRIPPER::setStartDf (double df)
{
	(*itsZcRates)[0] = - 365. * log (df) / double((*itsMaturities)[0]) ;
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_STRIPPER::setSwapRate (int i, double swaprate)
{
	std::vector<double>* startDates = itsDatestrip->GetFlowStartDates();
	int size = startDates->size(); 

	if (i<0 && i>size)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": the maturity Index of the stripper is negative or bigger than the max index");
	}
	(*itsSwapRates)[i] = swaprate;
}


////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_STRIPPER::strip ()
{	
	double Sratio, Df, lastDf;
	double Df0 = exp(-(*itsZcRates)[0]* (*itsMaturities)[0] / 365.0);

	std::vector<double>* interestTerm = itsDatestrip->GetInterestTerms();
	std::vector<double>* endDates = itsDatestrip->GetFlowEndDates();
	int size = interestTerm->size(); 

	for (int i=0; i<size; i++)
	{
		if (i==0)
		{
			Df =  Df0 / ( 1.0 + (*interestTerm)[i] * (*itsSwapRates)[i] );
			lastDf = Df;
		}
		else
		{
			Sratio = (*itsSwapRates)[i] / (*itsSwapRates)[i-1];
			Df  = Sratio * lastDf + (1-Sratio) * Df0;
			Df /= 1.0 + (*interestTerm)[i] * (*itsSwapRates)[i] ;
			lastDf = Df;
		}
		(*itsZcRates)[i+1] = - 365. * log (Df) / double((*itsMaturities)[i+1]) ;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_STRIPPER::df (double maturity) const
{
	double zcrate;

	// maturity out of range
	if (maturity<=(*itsMaturities)[0])
		zcrate = (*itsZcRates)[0];
	else if (maturity>=(*itsMaturities)[(*itsMaturities).size()-1])
		zcrate = (*itsZcRates)[(*itsZcRates).size()-1];
	else // linear interpolation
	{	
		int i(0);
		while (maturity>(*itsMaturities)[i]) i++;

		//maturity is in [Ti-1, Ti]
		double lambda = double (maturity - (*itsMaturities)[i-1]) / ((*itsMaturities)[i] - (*itsMaturities)[i-1]);

		zcrate = lambda * (*itsZcRates)[i] + (1.0-lambda) * (*itsZcRates)[i-1];
	}

	return exp(-zcrate * maturity / 365.0);
}


////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_STRIPPER::swapRate (double startFloatTime,double endFloatTime, const std::vector<double>& FixPayTimes, const std::vector<double>& FixPayPeriods) const
{
	int sizePay = FixPayTimes.size();
	double O1PV = 0.0;
	for(int i=0; i<sizePay;i++)
	{
		double payTime_i = FixPayTimes[i];
		double df_i = df(payTime_i);
		double fixPeriod_i = FixPayPeriods[i];
		O1PV+=df_i*fixPeriod_i;
	}

	double rate = (df(startFloatTime)-df(endFloatTime))/O1PV;

	return rate;
}


////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
double ARM_STRIPPER::stdFundingLeg (const std::vector<double>& startFlowTimes, const std::vector<double>& endFlowTimes,  const std::vector<double>& periods, const std::vector<double>& Notional, const std::vector<double>& Margin, std::vector<double>& LiborLeverage ) const
{
	int sizeFund = startFlowTimes.size();
	double floatPV = 0.0;
	for(int i=0; i<sizeFund;i++)
	{	
		double mult = LiborLeverage.size() > i ? LiborLeverage[i] : 1.0;
		double dfStart   = df(startFlowTimes[i]);
		double dfEnd     = df(endFlowTimes[i]);
		floatPV += mult * (dfStart - dfEnd)  * Notional[i];
		floatPV += periods[i] * dfEnd * Margin[i] * Notional[i];
	}

	return floatPV;
}

////////////////////////////////////////////////////
///	Class  : ARM_STRIPPER
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
string ARM_STRIPPER::toString(const string& indent, const string& nextIndent) const
{ 
	return ""; 
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

