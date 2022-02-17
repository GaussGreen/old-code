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
							const std::string& resetCalendar,
							int payRollingConvention,
							const std::string& payCalendar,
							bool monthBegin)
		:	FloatingCoupon(	nominal, inflPerf,start,end, numGap, dayCount, resetRollingConvention, 
							resetCalendar, payRollingConvention, payCalendar),
			itsPayOff(payoff),itsNumGap(numGap), itsDemGap(demGap), itsMonthBegin(monthBegin)
		{
			ARM_Date e(itsEndDate),s(itsStartDate) ;
			itsDemDate 	= s.AddPeriodNoAdj(itsDemGap.GetUnit(), itsDemGap.GetLength());
			itsNumDate 	= e.AddPeriodNoAdj(itsNumGap.GetUnit(), itsNumGap.GetLength());
			if (monthBegin) {
				itsDemDate = ARM_Date(1, itsDemDate.GetMonth(), itsDemDate.GetYear()) ;
				itsNumDate = ARM_Date(1, itsNumDate.GetMonth(), itsNumDate.GetYear()) ;
			}
			itsPaymentDate = itsEndDate ;
			itsPaymentDate.AdjustToBusDate(const_cast<char*>(itsPayCalendar.c_str()), itsPayRollingConvention);
		}

		InfFloatingCoupon (const InfFloatingCoupon& inflCF)
		: 	FloatingCoupon(	inflCF.itsNominal ,inflCF.itsIndex,inflCF.itsStartDate,inflCF.itsEndDate, 
							inflCF.itsResetGap, inflCF.itsDayCount, inflCF.itsResetRollingConvention, 
							inflCF.itsResetCalendar, inflCF.itsPayRollingConvention, inflCF.itsPayCalendar),
			itsDemDate(inflCF.itsDemDate),itsNumDate(inflCF.itsNumDate),itsPaymentDate(inflCF.itsPaymentDate),
			itsPayOff(inflCF.itsPayOff),itsDemGap(inflCF.itsDemGap),itsNumGap(inflCF.itsNumGap), itsMonthBegin(inflCF.itsMonthBegin)
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

		protected :
			ARM_CountedPtr<payOffType> itsPayOff ;
			ARM_Date itsDemDate, itsNumDate, itsPaymentDate ;
			Period itsDemGap, itsNumGap ;
			bool itsMonthBegin ;


	};
/*


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
							const std::string& resetCalendar,
							int payRollingConvention,
							const std::string& payCalendar,
							bool monthBegin)
		:	InfFloatingCoupon<payOffType> (	nominal, payoff, inflPerf,start,end, numGap, demGap,dayCount, resetRollingConvention, 
								resetCalendar, payRollingConvention, payCalendar, monthBegin),
			itsResetFreq(resetFreq)
		{
			ARM_DateStrip stripTease(	itsDemDate, itsNumDate, itsResetFreq,
										dayCount, itsResetCalendar.c_str(),
										K_UNADJUSTED, K_UNADJUSTED, K_SHORTSTART,
										0., itsResetFreq, 0., itsPayCalendar.c_str(),
										K_ARREARS, K_ARREARS, false );
		
			ARM_GP_Vector v(*stripTease.GetResetDates()) ;
			ARM_GP_Vector::const_iterator v_it	= v.begin() ;
			int n			= v.end()-v.begin(),i=0 ; 
			itsResetDates	= ARM_GP_T_Vector<ARM_Date>(n) ;
			for (;v_it!= v.end(); ++v_it, ++i) 
				itsResetDates[i] = ARM_Date(*v_it);   ;
		}

		InfComposedCoupon (const InfComposedCoupon& inflCF)
		: 	InfFloatingCoupon<payOffType>(	inflCF.itsNominal, inflCF.itsPayOff, /*(InfPerformance*)*/inflCF.itsIndex,
								inflCF.itsStartDate,inflCF.itsEndDate, inflCF.itsNumGap, 
								inflCF.itsDemGap,inflCF.itsDayCount,inflCF.itsResetRollingConvention,inflCF.itsResetCalendar,
								inflCF.itsPayRollingConvention, inflCF.itsPayCalendar,inflCF.itsMonthBegin), 
								itsResetFreq(inflCF.itsResetFreq),
								itsResetDates(inflCF.itsResetDates)
		{}

		virtual ~InfComposedCoupon (){}

		ARM_GP_T_Vector<ARM_Date> GetResDates()const {return itsResetDates;}

		double intrinsicValue() const
		{return itsPayOff->operator()(*this) ;}
		double GetAccrualTimeBetweenResets() const 
		{
			if (itsResetDates.size() > 1) {
				ARM_Date firstReset = *(itsResetDates.begin());
				ARM_Date secondReset = *(itsResetDates.begin()+1) ; 
				return CountYearsWithoutException(itsDayCount,firstReset,secondReset);
			}else{
				return GetAccrualTime();
			}
		}
		int GetAccrualDaysBetweenResets() const 
		{
			if (itsResetDates.size() > 1) {
				ARM_Date firstReset = *(itsResetDates.begin());
				ARM_Date secondReset = *(itsResetDates.begin()+1) ; 
				return DaysBetweenDates(itsDayCount,firstReset.GetJulian(),secondReset.GetJulian());
			}else{
				return GetAccrualDays();
			}
		}

		double floater() const {
			/*
			double floaterValue =itsIndex->fixing(itsNumDate);
			return floaterValue ;
			*/
			ARMTHROW(ERR_INVALID_ARGUMENT," Not Implemented");
		}
		protected :
			int itsResetFreq ;
			ARM_GP_T_Vector<ARM_Date> itsResetDates ;



	};

CC_END_NAMESPACE()

#endif
