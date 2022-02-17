/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  \file ARM_SmiledFRMfactory.cpp
 *  \brief
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 *
 */


/// this header comes first as it include some preprocessor constants
#include "gpmodels/SmiledFRMfactory.h"
#include "gpmodels/SmiledLMM.h"
#include "gpmodels/SmiledSMMcol.h"
#include "gpmodels/SmiledSMMdiag.h"

#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include <util\fromto.h>

CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRMfactory
///	Routine: CreateSmiledMarketModel 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_SmiledMM* ARM_SmiledFRMfactoryImp::CreateSmiledMarketModel(
	  CalibPattern pattern,const ARM_ZeroCurvePtr& zc, const ARM_ModelParams* params, size_t timeStepsNb,size_t gridSize,double stdDevNb,bool skipPDE, bool allowInterpol, ARM_ModelParamsSmiled::CalibProxy calibProxy) const
{
	switch( pattern )
	{
	case LIBOR:
		return new ARM_SmiledLMM(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy);
	case VMS_OLD:
		return new ARM_SmiledSMMdiag(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy,false);
	case VMS:
		return new ARM_SmiledSMMdiag(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy,true);
	case CMS_OLD:
		return new ARM_SmiledSMMcol(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy,false);
	case CMS:
		return new ARM_SmiledSMMcol(zc,params,timeStepsNb,gridSize,stdDevNb,skipPDE,allowInterpol,calibProxy,true);
	default:
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Unknown market model pattern, permitted is LIBOR, SWAPCOL and SWAPDIAG");
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRMfactory
///	Routine: CreateSmiledMarketModelDateStrip 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_DateStrip* ARM_SmiledFRMfactoryImp::CreateSmiledMarketModelDateStrip(
		CalibPattern pattern,
		const ARM_Date& startDate,			/// startDate of the strip
		const ARM_Date& endDate,			/// end Date
		long resetFreq,			
		long indexFreq,			
		int indexType,			
		long resetTiming,
		long dayCount,
		const char* resetCalendar,
		long fwdRule,
		long intRule,
		long stubRule,
		long resetGap) const
{
	bool isArrearsCpn	= (resetTiming == K_ARREARS);
	int maxFreq			= indexFreq<resetFreq?resetFreq:indexFreq;
	
	ARM_Date Date1(startDate);
	ARM_Date Date2(endDate);
	ARM_Date LastFixing;
	ARM_DateStrip* pDS;

	size_t nbmonth ;
	if		(indexFreq == K_ANNUAL)		nbmonth	= 12;
	else if (indexFreq == K_SEMIANNUAL)	nbmonth	= 6;
	else if (indexFreq == K_QUARTERLY)	nbmonth	= 3;
	else if (indexFreq == K_MONTHLY)		nbmonth	= 1;
	else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::GetCalibSchedule : index type not supported" );

	if (isArrearsCpn)
	{
		LastFixing = Date2;
		Date2.AddMonths(nbmonth);
	}
	else
	{
		size_t aux ;
		if (resetFreq == K_ANNUAL)			aux	= 12;
		else if (resetFreq == K_SEMIANNUAL)	aux	= 6;
		else if (resetFreq == K_QUARTERLY)	aux	= 3;
		else if (resetFreq == K_MONTHLY)	aux	= 1;
		else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::GetCalibSchedule : cpn freq not supported" );
		if (indexFreq < resetFreq)
		{
			Date2.AddMonths(-aux);
			LastFixing = Date2;
			Date2.AddMonths(nbmonth);
		}
		else
		{
			LastFixing=Date2;
			LastFixing.AddMonths(-aux);
		}
	}
	
	size_t nper = maxFreq / indexFreq;

	if (pattern == LIBOR || pattern == VMS_OLD || pattern == VMS)
	{
		pDS = new ARM_DateStrip(
			Date1,
			Date2,
			maxFreq,
			dayCount,
			resetCalendar,
			fwdRule,
			intRule,
			stubRule,
			resetGap,
			maxFreq,
			resetGap,
			resetCalendar);
	
		std::vector<double>& fwdEndDates    = pDS->GetFwdRateEndDates();
		std::vector<double>& fwdStartDates  = pDS->GetFwdRateStartDates();
		std::vector<double>& flowStartDates = pDS->GetFlowStartDates();
		std::vector<double>& flowEndDates   = pDS->GetFlowEndDates();

		if (pattern == LIBOR)
		{
			for (size_t i (0); i<fwdEndDates->size(); i++)
			{
				fwdStartDates->Elt(i)	= flowStartDates->Elt(i);
				fwdEndDates->Elt(i)		= flowEndDates->Elt(CC_Min(i + nper - 1, flowEndDates->size() - 1));
			}
		}
		if (pattern == VMS_OLD || pattern == VMS)
		{
			size_t ntot = fwdEndDates->size();
			size_t nbFamily = maxFreq / indexFreq;
			for (size_t i (0); i<ntot; i++)
			{
				fwdStartDates->Elt(i)	= flowStartDates->Elt(i);
				if (i<ntot-nbFamily)
					fwdEndDates->Elt(i)	= flowEndDates->Elt( (ntot-1) - (ntot-i) % nbFamily );
				else
					fwdEndDates->Elt(i)	= flowEndDates->Elt( (ntot-1) ) ;
			}
		}
	}
	else
	{
		ARM_Date Date3(LastFixing);
		Date3.AddMonths(indexType);

		pDS = new ARM_DateStrip(
			Date1,
			Date3,
			maxFreq,
			dayCount,
			resetCalendar,
			fwdRule,
			intRule,
			stubRule,
			resetGap,
			maxFreq,
			resetGap,
			resetCalendar);
	
		std::vector<double>& fwdEndDates    = pDS->GetFwdRateEndDates();
		std::vector<double>& fwdStartDates  = pDS->GetFwdRateStartDates();
		std::vector<double>& flowStartDates = pDS->GetFlowStartDates();
		std::vector<double>& flowEndDates   = pDS->GetFlowEndDates();

		size_t ntot = fwdEndDates->size();
		size_t nbFamily = maxFreq / indexFreq;
		size_t nbCMS = CountYears(K30_360,Date1,LastFixing)*maxFreq+1;
		size_t index = indexType*maxFreq/12;

		for (size_t i (0); i<ntot; i++)
		{
			fwdStartDates->Elt(i)	= flowStartDates->Elt(i);
			if (i<nbCMS)
				fwdEndDates->Elt(i)	= flowEndDates->Elt( i + index-1);
			else if (i<ntot-nbFamily)
				fwdEndDates->Elt(i)	= flowEndDates->Elt( (ntot-1) - (ntot-i) % nbFamily );
			else
				fwdEndDates->Elt(i)	= flowEndDates->Elt( (ntot-1) ) ;
		}
	}

	return pDS;
}


ARM_SingletonHolder<ARM_SmiledFRMfactoryImp> ARM_SmiledFRMfactory;

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
