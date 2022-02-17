/*!
 *
 * Copyright (c) IXIS CIB CM January 2005 Paris
 *
 *	\file vanillasecuritydensity.cpp
 *
 *  \brief vanilla security descr + density for HK calib
 *
 *	\author  A. Chaix
 *	\version 1.0
 *	\date November 2005
 */

/// gpcalib
#include "gpbase/utilityport.h"


/// gpcalib
#include "gpcalib/vanillasecuritydensity.h"

/// gpbase
#include "gpbase/datestrip.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/argconvdefault.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// kernel
#include <crv/zerocurv.h>
#include <ccy/currency.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: constructor (default + contextual)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensity::ARM_VanillaSecurityDensity( double resetDate, 
														double startDate,
														double endDate,
														ARM_DensityFunctorPtr densityFunctor, 
														int frequency,
														int daycount,
														int stubRule,
														double weight,
														double adjfwdadd,
														double adjfwdmult,
														ARM_ZeroCurvePtr ZeroCurve
														)
:	itsResetDate		(resetDate),
	itsStartDate		(startDate),
	itsEndDate			(endDate),
	itsDensityFunctor	(densityFunctor),
	itsFrequency		(frequency),
	itsDayCount			(daycount),
	itsStubRule			(stubRule),
	itsZeroCurve		(ZeroCurve),
	itsWeight			(weight),
	itsAdjFwdAdd		(adjfwdadd),
	itsAdjFwdMult		(adjfwdmult),
	itsCalibSchedule	(NULL),
	itsInterestTerms		(0),
	itsPayDatesRelIndexes	(0),
	itsPayDatesAbsIndexes	(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensity::ARM_VanillaSecurityDensity( const ARM_VanillaSecurityDensity& rhs)
:	itsResetDate		(rhs.itsResetDate),
	itsStartDate		(rhs.itsStartDate),
	itsEndDate			(rhs.itsEndDate),
	itsDensityFunctor	(CreateClonedPtr(&*rhs.itsDensityFunctor)),
	itsZeroCurve		(rhs.itsZeroCurve),		 /// shared
	itsWeight			(rhs.itsWeight),
	itsAdjFwdAdd		(rhs.itsAdjFwdAdd),
	itsAdjFwdMult		(rhs.itsAdjFwdMult),
	itsFrequency		(rhs.itsFrequency),
	itsDayCount			(rhs.itsDayCount),
	itsStubRule			(rhs.itsStubRule),
	itsCalibSchedule	(rhs.itsCalibSchedule),	/// ASSOC
	itsInterestTerms	(rhs.itsInterestTerms),
	itsPayDatesRelIndexes (rhs.itsPayDatesRelIndexes),
	itsPayDatesAbsIndexes (rhs.itsPayDatesAbsIndexes)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensity::~ARM_VanillaSecurityDensity()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: AssociateCalibSchedule
///	Returns: 
///	Action : associate calib schedule 
////////////////////////////////////////////////////

void ARM_VanillaSecurityDensity::AssociateCalibSchedule(ARM_DateStripPtr calibSchedule)
{
	if(fabs(itsWeight) > K_DOUBLE_TOL) AssociateCalibSchedule(calibSchedule, itsStartDate, itsEndDate, itsFrequency, itsDayCount, itsInterestTerms, itsPayDatesRelIndexes, itsPayDatesAbsIndexes);
}

void ARM_VanillaSecurityDensity::AssociateCalibSchedule(ARM_DateStripPtr calibSchedule, 
														double StartDate, double EndDate, int Frequency, int DayCount,
														ARM_GP_Vector& InterestTerms, ARM_IntVector& PayDatesRelIndexes, ARM_IntVector& PayDatesAbsIndexes)
{
	if(fabs(itsWeight) < K_DOUBLE_TOL) return;

	/// set calib schedule
	itsCalibSchedule = calibSchedule;

	/// find indexes
	size_t startIndex, endIndex;
	IndexesInCalibSchedule(StartDate, EndDate, startIndex, endIndex);
	/// nb periods in terms of calibSchedule
	/// if a different freq is provided, there will be less periods
	size_t nbPeriods = endIndex - startIndex + 1;
	size_t i;

	/// default case : rate based upon calibSchedule
	if (Frequency == GETDEFAULTVALUE)
	{
		PayDatesRelIndexes.resize(nbPeriods);
		for (i=0; i<nbPeriods; i++)
			PayDatesRelIndexes[i] = i;
	}	
	/// consider only a subset of dates (according to frequency and stub rule)
	else
	{	/// compute freq of calib schedule
		int calibSchedFreq = itsCalibSchedule->GetResetFreq();

		if (calibSchedFreq>14)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::AssociateCalibSchedule : maximum supported frequency is MONTHLY" )
		/// calibSchedFreq = 13, can happen if freq = MONTHLY for the month of february
		else if (calibSchedFreq>12)
			calibSchedFreq = 12;

		
		/// a litte check
		if (Frequency > calibSchedFreq)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::AssociateCalibSchedule : security frequency is required to be < calib freq");
		
		if (calibSchedFreq % Frequency  != 0)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::AssociateCalibSchedule : oups this should not happen");

		/// nice case : no broken period
		if (nbPeriods % (calibSchedFreq / Frequency) == 0)
		{
			int increment = calibSchedFreq / Frequency;
			PayDatesRelIndexes.resize( nbPeriods / increment);
			PayDatesRelIndexes[0] = increment - 1;
			for (i = 1; i<PayDatesRelIndexes.size(); i++)
				PayDatesRelIndexes[i] = PayDatesRelIndexes[i-1] + increment;
		}
		else
		{		
			int increment    = calibSchedFreq / Frequency;
			int residu		 = nbPeriods % increment; /// in terms of nb of periods in calib schedule
			int nbStdPeriods = ( nbPeriods -  residu ) / increment;
			
			PayDatesRelIndexes.resize(1 + nbStdPeriods);

			if (itsStubRule == GETDEFAULTVALUE || itsStubRule == K_SHORTSTART)
			{						
				PayDatesRelIndexes[0] = residu - 1;

				for (i = 1; i<PayDatesRelIndexes.size(); i++)
					PayDatesRelIndexes[i] = PayDatesRelIndexes[i-1] + increment;

			}
			else if (itsStubRule == K_SHORTEND)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::AssociateCalibSchedule : stub rule SHORT END does not work yet");

				if (nbStdPeriods)
					PayDatesRelIndexes[0] = increment - 1;
				else
					PayDatesRelIndexes[0] = 0;
								
				for (i = 1; i<PayDatesRelIndexes.size()-1; i++)
					PayDatesRelIndexes[i] = PayDatesRelIndexes[i-1] + increment;

				if (nbStdPeriods)
					PayDatesRelIndexes[PayDatesRelIndexes.size()-1] 
						= PayDatesRelIndexes[PayDatesRelIndexes.size()-2] + residu;
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::AssociateCalibSchedule : stub rule not supported (use SHORT START or SHORT END)");
		}
	}

	/// compute abs indexes from relative indexes
	PayDatesAbsIndexes = PayDatesRelIndexes + ARM_IntVector(PayDatesRelIndexes.size(), startIndex);

	/// *** compute interest terms ****
	/// case where we use calib schedule
	if (DayCount == GETDEFAULTVALUE && Frequency == GETDEFAULTVALUE)
	{				
		InterestTerms.resize(nbPeriods);
		for (i=0; i<nbPeriods; i++)
			InterestTerms[i] = itsCalibSchedule->GetInterestTerms()->Elt(startIndex + i);
	}
	///  recompute interest terms (daycount or/and freq are not the one of calib schedule)
	else
	{	
		int dayCount;
		if (DayCount == GETDEFAULTVALUE)
			dayCount = itsCalibSchedule->GetDayCount();
		else
			dayCount = DayCount;

		InterestTerms.resize(PayDatesAbsIndexes.size());
		
		double startDate = itsCalibSchedule->GetFlowStartDates()->Elt(startIndex);

		for (i=0; i<InterestTerms.size(); i++)
		{
			double endDate		= itsCalibSchedule->GetFlowEndDates()->Elt  (PayDatesAbsIndexes[i]);
			InterestTerms[i] = CountYears( dayCount, startDate, endDate );
			startDate = endDate;
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: IndexesInCalibSchedule
///	Returns: void
///	Action : find index for reset, start and end date
///			 in itsCalibSchedule
////////////////////////////////////////////////////
void
ARM_VanillaSecurityDensity::IndexesInCalibSchedule(double StartDate, double EndDate, size_t& resetAndStartIndex, size_t& endIndex) const
{ 
	if (itsCalibSchedule.IsNull())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : calib schedule has't been setted");

	if(fabs(itsWeight) < K_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : calib schedule has't been setted (nul weight)");

	/// consistency checks
	ARM_GP_Vector* resetDates	= itsCalibSchedule->GetResetDates();
	ARM_GP_Vector* startDates	= itsCalibSchedule->GetFlowStartDates();
	ARM_GP_Vector* endDates		= itsCalibSchedule->GetFlowEndDates();
	size_t N = resetDates->size();
	
	int resetIndex (-1);
	for (size_t i(0); i<N; i++)
		if (resetDates->Elt(i) == itsResetDate)
		{
			resetIndex = i;
			break;
		}
	
	if (resetIndex == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : reset date not found in calib schedule");

	if (StartDate != startDates->Elt(resetIndex))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : start date does not match");

	resetAndStartIndex = resetIndex;

	int _endIndex (-1);
	for (i = resetIndex; i<N; i++)
		if (endDates->Elt(i) == EndDate)
		{
			_endIndex = i;
			break;
		}

	if (_endIndex == -1)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : end date not found in calib schedule");

	endIndex = _endIndex;	
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: ComputeForwardRate
///	Returns: int
///	Action : Forward (libor or swap rate)
///			 needed to compute rate quantile
////////////////////////////////////////////////////
double ARM_VanillaSecurityDensity::ComputeForwardRate (double StartDate, double EndDate, const ARM_IntVector& PayDatesAbsIndexes, const ARM_GP_Vector& InterestTerms) const
{	
	if (itsCalibSchedule.IsNull())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::ComputeForwardRate : calib schedule has not been set");

	if(fabs(itsWeight) < K_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : calib schedule has't been setted (nul weight)");

	double asOf			= itsZeroCurve->GetAsOfDate().GetJulian();
	ARM_GP_Vector* ends = itsCalibSchedule->GetFlowEndDates();
		
	double fwd =	getZeroCurve()->DiscountPrice((StartDate-asOf)/K_YEAR_LEN)
				-	getZeroCurve()->DiscountPrice((EndDate-asOf)/K_YEAR_LEN);

	double level (0.0);

	size_t idx;
	for (size_t i = 0; i<PayDatesAbsIndexes.size(); i++)
	{
		idx    = PayDatesAbsIndexes[i];
		level += InterestTerms[i] * getZeroCurve()->DiscountPrice(((*ends)[idx]-asOf)/K_YEAR_LEN);
	}

	fwd /= level;
		
	return itsAdjFwdMult * fwd + itsAdjFwdAdd;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: ComputeLevel
///	Returns: int
///	Action : sensi
////////////////////////////////////////////////////
double ARM_VanillaSecurityDensity::ComputeLevel (double startdate, double enddate, const ARM_IntVector& PayDatesAbsIndexes, const ARM_GP_Vector& InterestTerms) const
{	
	if (itsCalibSchedule.IsNull())
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::ComputeLevel : calib schedule has not been set");

	if(fabs(itsWeight) < K_DOUBLE_TOL)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " ARM_VanillaSecurityDensity::indexesInCalibSchedule : calib schedule has't been setted (nul weight)");

	double asOf			= itsZeroCurve->GetAsOfDate().GetJulian();
	ARM_GP_Vector* ends = itsCalibSchedule->GetFlowEndDates();
		
	double level (0.0);

	size_t idx;
	for (size_t i = 0; i<PayDatesAbsIndexes.size(); i++)
	{
		idx    = PayDatesAbsIndexes[i];
		level += InterestTerms[i] * getZeroCurve()->DiscountPrice(((*ends)[idx]-asOf)/K_YEAR_LEN);
	}

	return level;
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: UpdateStates
///	Returns: void
///	Action : computes 1/Numeraire at date lineIdx
////////////////////////////////////////////////////
void ARM_VanillaSecurityDensity::UpdateStates( const ARM_PricingStatesPtr& states, const ARM_GP_VectorPtr& probabilities, size_t lineIdx ) const
{		
	double asOf = itsZeroCurve->GetAsOfDate().GetJulian();

	/// Compute forward rate
	double fwd = ComputeForwardRate();
		
	/// Apply Quantile function
	/// vec below contains Rate( ResetDate, StartDate, EndDate );
	ARM_GP_VectorPtr vec = getDensityFunctor()->Quantile( probabilities, fwd, (itsResetDate - asOf) / K_YEAR_LEN );
	

	/// --------------------------------------
	/// last period : caplet case
	/// --------------------------------------
	if( states->GetPayoffsSize() == 1 )
	{
		//double delta = itsCalibSchedule->GetInterestTerms()->Elt(startIndex);
		double delta = itsInterestTerms[0];

		/// In case we are at the last date. _we basically put one/Zc
		ARM_MemPool_Matrix::iterator iter = states->payoffsBeginIterator(0);
		ARM_GP_Vector::iterator vecIter = vec->begin(), vecEnd = vec->end();

		for( ; vecIter != vecEnd ; ++vecIter,++iter )
			(*iter) = (1.+delta*(*vecIter));
	}
	
	/// --------------------------------------
	/// intermediate period : libor or swaprate
	/// --------------------------------------
	else
	{
		size_t N = states->GetPayoffsSize();
		
		/// iterResult = iterator for one/Numeraire
		double delta;
		ARM_MemPool_Matrix::iterator iterResult, iterDf;
		ARM_GP_Vector::iterator vecIter, vecEnd;
		
		vecEnd	= vec->end();
		size_t nbFlows = itsPayDatesRelIndexes.size();

		for (size_t i = 0; i<nbFlows; i++)
		{
			delta		= itsInterestTerms[i];
			vecIter		= vec->begin(); 
			iterResult	= states->payoffsBeginIterator(N-1);
						
			bool isLastRateFlow  = ( i == itsPayDatesRelIndexes.size() - 1 );
			bool isLastSchedFlow = ( itsPayDatesAbsIndexes[i] == itsCalibSchedule->GetResetDates()->size() - 1 );

			if (!isLastSchedFlow)
				iterDf		= states->payoffsBeginIterator(N - itsPayDatesRelIndexes[i] - 2);
						
			double mult;

			for(; vecIter != vecEnd; ++vecIter, ++iterResult, ++iterDf )
			{
				mult = delta * (*vecIter);
				
				/// Last rate flow : (1+cvg*df) instead of cvg*df
				if (isLastRateFlow) 
					mult += 1.0;
				
				/// Last schedule flow : df(.,T[n+1]) / df(.,T[n+1]) = 1  !!
				if (!isLastSchedFlow)
					(*iterResult) += mult * (*iterDf) ;
				else
					(*iterResult) += mult ;
			}
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: toString
///	Returns: string
///	Action : toString
////////////////////////////////////////////////////
string ARM_VanillaSecurityDensity::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << indent << "\n";
    os << indent << "ARM_VanillaSecurityDensity\n";
    os << indent << "---------------------------\n";

	os << indent << CC_NS(std,endl);
	
	ARM_Date resetDate = itsResetDate;
	ARM_Date startDate = itsStartDate;
	ARM_Date endDate   = itsEndDate;
	
	os << indent << "Reset Date =\t" << resetDate.toString() << CC_NS(std,endl);
	os << indent << "Start Date =\t" << startDate.toString() << CC_NS(std,endl);
	os << indent << "End Date   =\t" << endDate.toString()   << CC_NS(std,endl);
	
	if (itsFrequency != GETDEFAULTVALUE)
		os << indent << "Frequency  =\t" << ARM_ArgConvReverse_LgNameFrequency.GetString( itsFrequency ) << CC_NS(std,endl);
	else
		os << indent << "Frequency  =\t" << "DEFAULT" << CC_NS(std,endl);
	
	if (itsDayCount != GETDEFAULTVALUE)
		os << indent << "Day Count  =\t" << ARM_ArgConvReverse_LgNameDayCount.GetString( itsDayCount )  << CC_NS(std,endl);
	else
		os << indent << "Day Count  =\t" << "DEFAULT" << CC_NS(std,endl);

	if (itsStubRule != GETDEFAULTVALUE)
		os << indent << "Stub Rule  =\t" << ARM_ArgConvReverse_StubRules.GetString( itsStubRule) << CC_NS(std,endl);
	else
		os << indent << "Stub Rule  =\t" << "DEFAULT" << CC_NS(std,endl);

	os << indent << CC_NS(std,endl);

	if (!itsCalibSchedule.IsNull())
	{
		os << indent << "Start Dates\t\t End Dates\t Interest Terms\n";
		
		size_t nFlows = itsInterestTerms.size();
		ARM_GP_Vector* startDates = itsCalibSchedule->GetFlowStartDates();
		ARM_GP_Vector* endDates   = itsCalibSchedule->GetFlowEndDates();

		double startDate  = (*startDates)[ itsPayDatesAbsIndexes[0] - itsPayDatesRelIndexes[0] ];

		for (size_t i = 0; i < nFlows; ++i )
		{			
			ARM_Date tmpDate1 = startDate;
			ARM_Date tmpDate2 = (*endDates)[itsPayDatesAbsIndexes[i]];
			startDate = (*endDates)[itsPayDatesAbsIndexes[i]];
						
			os  << indent << tmpDate1.toString() << "\t " 
				<< tmpDate2.toString() << "\t " 
				<< itsInterestTerms[i] << "\n";
		}
	}
	
	os << indent << "\n" ;
	os << indent << "\n" ;
	os << indent <<  "Associated Density Functor :"<< CC_NS(std,endl);
	os << itsDensityFunctor->toString(indent,nextIndent) << CC_NS(std,endl);
			
	return os.str();	

}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////

ARM_VanillaSecurityDensitySpread::ARM_VanillaSecurityDensitySpread( double resetDate, 
														double startDate,
														double endDate,
														double startDate2,
														double endDate2,
														ARM_DensityFunctorPtr densityFunctor, 
														int frequency,
														int daycount,
														int frequency2,
														int daycount2,
														int stubRule,
														ARM_ZeroCurvePtr ZeroCurve)
:	ARM_VanillaSecurityDensity(resetDate, startDate, endDate, densityFunctor, frequency, daycount, stubRule, 1., 0., 1., ZeroCurve),
	itsStartDate2		(startDate2),
	itsEndDate2			(endDate2),
	itsFrequency2		(frequency2),
	itsDayCount2		(daycount2),
	itsInterestTerms2		(0),
	itsPayDatesRelIndexes2	(0),
	itsPayDatesAbsIndexes2	(0)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensitySpread::ARM_VanillaSecurityDensitySpread( const ARM_VanillaSecurityDensitySpread& rhs)
:	ARM_VanillaSecurityDensity(rhs),
	itsStartDate2		(rhs.itsStartDate2),
	itsEndDate2			(rhs.itsEndDate2),
	itsFrequency2		(rhs.itsFrequency2),
	itsDayCount2		(rhs.itsDayCount2),
	itsInterestTerms2	(rhs.itsInterestTerms2),
	itsPayDatesRelIndexes2 (rhs.itsPayDatesRelIndexes2),
	itsPayDatesAbsIndexes2 (rhs.itsPayDatesAbsIndexes2)
{
}

ARM_VanillaSecurityDensitySpread::~ARM_VanillaSecurityDensitySpread()
{
}

void ARM_VanillaSecurityDensitySpread::AssociateCalibSchedule(ARM_DateStripPtr calibSchedule)
{
	ARM_VanillaSecurityDensity::AssociateCalibSchedule(calibSchedule, itsStartDate, itsEndDate, itsFrequency, itsDayCount, itsInterestTerms, itsPayDatesRelIndexes, itsPayDatesAbsIndexes);
	ARM_VanillaSecurityDensity::AssociateCalibSchedule(calibSchedule, itsStartDate2, itsEndDate2, itsFrequency2, itsDayCount2, itsInterestTerms2, itsPayDatesRelIndexes2, itsPayDatesAbsIndexes2);

	if(itsDensityFunctor != ARM_DensityFunctorPtr(NULL))
	{
		if(dynamic_cast<ARM_BiSABRDensityFunctor*>(&*itsDensityFunctor) != NULL)
		{
			dynamic_cast<ARM_BiSABRDensityFunctor*>(&*itsDensityFunctor)->SetForward1(ComputeForwardRate1());
			dynamic_cast<ARM_BiSABRDensityFunctor*>(&*itsDensityFunctor)->SetForward2(ComputeForwardRate2());
		}
	}
}

string ARM_VanillaSecurityDensitySpread::toString(const string& indent,const string& nextIndent) const
{
	return ARM_VanillaSecurityDensity::toString(indent, nextIndent);
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensityFX
///	Routine: constructor (default + contextual)
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_VanillaSecurityDensityFX::ARM_VanillaSecurityDensityFX(double resetDate,
				ARM_DensityFunctorPtr densityFunctor, 
				ARM_ZeroCurvePtr DomZeroCurve,
				ARM_ZeroCurvePtr FgnZeroCurve,
				double FXSpot)
		:
 ARM_VanillaSecurityDensity(resetDate,0,0,densityFunctor,GETDEFAULTVALUE,GETDEFAULTVALUE,GETDEFAULTVALUE,1.,0.,1.,DomZeroCurve),
	 itsFgnZeroCurve(FgnZeroCurve),
	 itsFXSpot(FXSpot)
{
	string setCal = ARM_PricingFunctionEquity::GetSettlementCalendar(&*DomZeroCurve,&*FgnZeroCurve);
	double setGap = ARM_PricingFunctionEquity::GetSettlementGap(&*DomZeroCurve,&*FgnZeroCurve);

	ARM_Date asOfDate = DomZeroCurve->GetAsOfDate();
	ARM_Date spotDate = asOfDate;
	spotDate.NextBusinessDay(setGap,setCal.c_str());

	itsFXSpot /= FgnZeroCurve->DiscountPrice((spotDate.GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN);
	itsFXSpot *= DomZeroCurve->DiscountPrice((spotDate.GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN);


	ARM_Date settlementDate(resetDate);
	settlementDate.NextBusinessDay(setGap, const_cast<char*>(setCal.c_str()));
	itsSettlementDate = settlementDate.GetJulian();
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensityFX
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_VanillaSecurityDensityFX::ARM_VanillaSecurityDensityFX( const ARM_VanillaSecurityDensityFX& rhs)
:
ARM_VanillaSecurityDensity(rhs),
itsFgnZeroCurve(rhs.itsFgnZeroCurve),
itsFXSpot(rhs.itsFXSpot),
itsSettlementDate(rhs.itsSettlementDate)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensityFX
///	Routine: destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensityFX::~ARM_VanillaSecurityDensityFX()
{
}


/// forward FX Rate
double ARM_VanillaSecurityDensityFX::ComputeForwardRate() const
{

	double asOf			= getZeroCurve()->GetAsOfDate().GetJulian();
	int spotDays		= getZeroCurve()->GetCurrencyUnit()->GetSpotDays();
	string ccyName		= getZeroCurve()->GetStrCurrency();
	ccyName				+=	getFgnZeroCurve()->GetStrCurrency();
	/*ARM_Date fixingDate = ARM_Date(itsResetDate);

	fixingDate.NextBusinessDay(spotDays, ccyName.c_str() ).GetJulian();
	double itsFixingDate = fixingDate.GetJulian();

	double dfRatioAtSpotDate = getFgnZeroCurve()->DiscountPrice((itsFixingDate-asOf)/K_YEAR_LEN)/
								getZeroCurve()->DiscountPrice((itsFixingDate-asOf)/K_YEAR_LEN);*/
	
	double dfRatioAtSpotDate = getFgnZeroCurve()->DiscountPrice((itsResetDate-asOf)/K_YEAR_LEN)/getZeroCurve()->DiscountPrice((itsResetDate-asOf)/K_YEAR_LEN);
	
	double fwdFX =	dfRatioAtSpotDate*itsFXSpot;

	return fwdFX;
}


////////////////////////////////////////////////////
///	Class  : ARM_VanillaSecurityDensity
///	Routine: toString
///	Returns: string
///	Action : toString
////////////////////////////////////////////////////
string ARM_VanillaSecurityDensityFX::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;
    os << indent << "\n";
    os << indent << "ARM_VanillaSecurityDensityFX\n";
    os << indent << "---------------------------\n";

	os << indent << CC_NS(std,endl);
	
	ARM_Date resetDate = getResetDate();
	ARM_Date settlementDate = getSettlementDate();

	
	os << indent << "Reset Date =\t" << resetDate.toString() << CC_NS(std,endl);
	os << indent << "SettlementDate Date =\t" << settlementDate.toString() << CC_NS(std,endl);

	os << indent << CC_NS(std,endl);
	
	os << indent << "\n" ;
	os << indent << "\n" ;
	os << indent <<  "Associated Density Functor :"<< CC_NS(std,endl);
	os << getDensityFunctor()->toString(indent,nextIndent) << CC_NS(std,endl);
			
	return os.str();	

}

CC_END_NAMESPACE()
