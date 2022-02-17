/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file infCoupon.h
 *	\author  Francois Poitou
 */

#ifndef _INFCOUPON_H
#define _INFCOUPON_H

/// gpbase
#include <gpbase/countedptr.h>				//ARM_CountedPtr
/// kernel
#include <glob/dates.h>
/// gpinflation
#include "gpinflation/infPerformance.h"
#include "gpinflation/floatingCoupon.h"


CC_BEGIN_NAMESPACE( ARM )

	template<class payOffType>
	class InfFloatingCoupon : public FloatingCoupon {

		public :
			InfFloatingCoupon (){}
		InfFloatingCoupon (	double nominal,
							ARM_CountedPtr<payOffType> payoff,
							ARM_CountedPtr<InfPerformance> inflPerf,
							const ARM_Date& start,
							const ARM_Date& end,
							const Period& numGap,
							const Period& demGap,
							int dayCount,
							int resetRollingConvention,
							char* resetCalendar,
							int payRollingConvention,
							char* payCalendar,
							bool monthBegin)
		:	FloatingCoupon(	nominal, inflPerf,start,end, numGap, dayCount, resetRollingConvention, 
							resetCalendar, payRollingConvention, payCalendar),
			itsPayOff(payoff),itsNumGap(numGap), itsDemGap(demGap), itsMonthBegin(monthBegin)
		{
			itsResetDates.clear() ;
			ARM_Date e(itsEndDate),s(itsStartDate) ;
			itsDemDate 	= s.AddPeriodNoAdj(itsDemGap.GetUnit(), itsDemGap.GetLength());
			itsNumDate 	= e.AddPeriodNoAdj(itsNumGap.GetUnit(), itsNumGap.GetLength());
			if (monthBegin) {
				itsDemDate = ARM_Date(1, itsDemDate.GetMonth(), itsDemDate.GetYear()) ;
				itsNumDate = ARM_Date(1, itsNumDate.GetMonth(), itsNumDate.GetYear()) ;
			}
			itsPaymentDate = itsEndDate ;
			itsPaymentDate.AdjustToBusDate(itsPayCalendar, itsPayRollingConvention);
			itsResetDates.push_back(itsNumDate) ;
		}

		InfFloatingCoupon (const InfFloatingCoupon& inflCF)
		: 	FloatingCoupon(inflCF.itsNominal ,inflCF.itsInflPerf,inflCF.itsStartDate,inflCF.itsEndDate ),
			itsPayOff(inflCF.itsPayOff),itsYieldPeriodlag(inflCF.itsYieldPeriodlag), itsMonthBegin(inflCF.itsMonthBegin),
			itsResetDates(inflCF.itsResetDates)
		{}

		virtual ~InfFloatingCoupon (){}

		double intrinsicValue() const
		{return itsPayOff->operator()(*this) ;}

		double floater() const {
			/*
			double floaterValue =itsIndex->fixing(itsNumDate)/itsIndex->fixing(itsDemDate)-1;
			return floaterValue ;
			*/
			ARMTHROW(ERR_INVALID_ARGUMENT," Not Implemented");

		}


		ARM_Date GetDemDate() const {return itsDemDate;}
		ARM_Date GetNumDate() const {return itsNumDate;}
		ARM_Date GetPayDate() const {return itsPaymentDate ;}
		ARM_GP_T_Vector<ARM_Date> GetResDates()const {return itsResetDates;}

		protected :
			ARM_CountedPtr<payOffType> itsPayOff ;
			ARM_Date itsDemDate, itsNumDate, itsPaymentDate ;
			Period itsDemGap, itsNumGap ;
			bool itsMonthBegin ;
			ARM_GP_T_Vector<ARM_Date> itsResetDates ;


	};
/*

  //PFFFFF......

// Forward Rules 
#define S_FORWARD_RULES      "Forward Rules"
#define K_PREVIOUS         -1
#define K_MOD_PREVIOUS     -2
#define K_FOLLOWING         1
#define K_MOD_FOLLOWING     2
 
// Interest Rules 
#define S_INTEREST_RULES "Interest Rules"
#define K_ADJUSTED          1
#define K_UNADJUSTED        0
#define K_MATUNADJUSTED     2
 
// Stub Rules
#define S_STUB_RULES         "Stub Rules"
#define K_SHORTSTART        1
#define K_LONGSTART         2
#define K_SHORTEND          3
#define K_LONGEND           4
 
*/
template<class payOffType>
class InfComposedCoupon : public InfFloatingCoupon<payOffType> {

		public :

		InfComposedCoupon (	double nominal,
							ARM_CountedPtr<payOffType> payoff,
							ARM_CountedPtr<InfPerformance> inflPerf,
							const ARM_Date& start,
							const ARM_Date& end,
							int resetFreq,
							int dayCount,
							const Period& numGap,
							const Period& demGap,
							int resetRollingConvention,
							char* resetCalendar,
							int payRollingConvention,
							char* payCalendar,
							bool monthBegin)
		:	InfFloatingCoupon<payOffType> (	nominal, payoff, inflPerf,start,end, numGap, demGap,dayCount, resetRollingConvention, 
								resetCalendar, payRollingConvention, payCalendar, monthBegin),
			itsResetFreq(resetFreq)
		{
			itsResetDates.clear() ;
			ARM_DateStrip stripTease(	itsDemDate, itsNumDate, itsResetFreq,
										dayCount, itsResetCalendar,
										K_UNADJUSTED, K_UNADJUSTED, K_SHORTSTART,
										0., itsResetFreq, 0., itsPayCalendar,
										K_ARREARS, K_ARREARS, false );
		
			//PFFFFF
			ARM_GP_Vector v(*stripTease.GetResetDates()) ;
			ARM_GP_Vector::const_iterator v_it = v.begin() ;
			for (;v_it!= v.end(); ++v_it) 
				itsResetDates.push_back(ARM_Date(*v_it));   ;
		}

		InfComposedCoupon (const InfComposedCoupon& inflCF)
		: 	InfFloatingCoupon(	inflCF.itsNominal, inflCF.itsPayOff, inflCF.itsInflPerf,
								inflCF.itsStartDate,inflCF.itsEndDate, inflCF.itsNumGap, 
								inflCF.itsDemGap,inflCF.itsResetRollingConvention,inflCF.itsResetCalendar,
								inflCF.itsPayRollingConvention, inflCF.itsPayCalendar,inflCF.itsMonthBegin), 
								itsResetFreq(inflCF.itsResetFreq),
								itsResetDates(inflCF.itsResetDates)
		{}

		virtual ~InfComposedCoupon (){}

		double intrinsicValue() const
		{return itsPayOff->operator()(*this) ;}

		double floater() const {
			/*
			double floaterValue =itsIndex->fixing(itsNumDate);
			return floaterValue ;
			*/
			ARMTHROW(ERR_INVALID_ARGUMENT," Not Implemented");
		}
		protected :
			int itsResetFreq ;


	};

CC_END_NAMESPACE()

#endif
