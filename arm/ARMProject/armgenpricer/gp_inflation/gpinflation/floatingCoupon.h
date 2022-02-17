/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file floatingCoupon.h
 *	\author  Francois Poitou
 */

#ifndef _FLOATINGCOUPON_H
#define _FLOATINGCOUPON_H

/// gpbase
#include <gpbase/countedptr.h>				//ARM_CountedPtr
/// kernel
#include <inst/irindex.h>			//ARM_IRIndex
#include <glob/dates.h>			//ARM_Date

#include "glob/period.h"


CC_BEGIN_NAMESPACE( ARM )

	class CashFlow  : public ARM_Object
	{
		public:
			CashFlow(	double nominal,
						const ARM_Date& start,
						const ARM_Date& end,
						int dayCount):
						itsStartDate(start), itsEndDate(end), itsNominal(nominal), itsDayCount(dayCount)
			{}

			CashFlow(const CashFlow& cf):
				itsStartDate(cf.itsStartDate),
				itsEndDate(cf.itsEndDate),
				itsNominal(cf.itsNominal),
				itsDayCount(cf.itsDayCount)
			{}
			virtual ~CashFlow()
			{}
			CashFlow& operator= (const CashFlow& cf) 
			{
				itsStartDate	= cf.itsStartDate; 
				itsEndDate		= cf.itsEndDate ;
				itsNominal		= cf.itsNominal; 
				itsDayCount		= cf.itsDayCount ;
				return *this;
			}
			double GetNominal() const {return itsNominal;}
			ARM_Date GetStartDate()const {return itsStartDate ;}
			ARM_Date GetEndDate()const {return itsEndDate ;}
			double GetAccrualTime() const 
			{
				return CountYearsWithoutException(itsDayCount,itsStartDate,itsEndDate);
			}
			int GetAccrualDays() const 
			{
				return DaysBetweenDates(itsDayCount,itsStartDate.GetJulian(),itsEndDate.GetJulian());
			}

			
		protected :
			ARM_Date itsStartDate, itsEndDate ;
			double itsNominal;
			int itsDayCount ;

	};

	class FloatingCoupon : public CashFlow 
	{
	public:
			FloatingCoupon(	double nominal,
							ARM_CountedPtr<ARM_IRIndex> index,
							const ARM_Date& start,
							const ARM_Date& end, 
							const Period& resetGap,
							int dayCount,
							int resetRollingConvention,
							const std::string& resetCalendar,
							int payRollingConvention,
							const std::string& payCalendar):
			
							itsResetGap(resetGap), itsResetCalendar(resetCalendar), itsIndex(index),
							itsResetRollingConvention(resetRollingConvention),itsPayCalendar(payCalendar),
							itsPayRollingConvention(payRollingConvention), CashFlow(nominal,start,end, dayCount)
			{
				//itsIndex = (index ? (ARM_IRIndex *) index->Clone() : NULL); 
			}

			FloatingCoupon(const FloatingCoupon& cf):
					
				itsResetGap(cf.itsResetGap),itsIndex(cf.itsIndex),
				itsResetCalendar(cf.itsResetCalendar),
				itsResetRollingConvention(cf.itsResetRollingConvention),
				itsPayCalendar(cf.itsPayCalendar),
				itsPayRollingConvention(cf.itsPayRollingConvention), 
				CashFlow(cf.itsNominal, cf.itsStartDate,cf.itsEndDate, cf.itsDayCount )
			{
				//itsIndex = cf.itsIndex? (ARM_IRIndex *) cf.itsIndex->Clone(): NULL;
			}
			virtual ~FloatingCoupon()
			{
				//if (itsIndex)
				  // delete itsIndex;
			}

			ARM_CountedPtr<ARM_IRIndex> GetIndex() const {return itsIndex ;}
			virtual double floater() const = 0 ;
			virtual	double intrinsicValue() const = 0 ;


		protected :
			ARM_CountedPtr<ARM_IRIndex> itsIndex ;
			Period itsResetGap;
			std::string itsResetCalendar;				
			int itsResetRollingConvention;		
			std::string itsPayCalendar;
			int itsPayRollingConvention;	

	};


std::ostream& operator<<(std::ostream& os, const CashFlow& cf );

CC_END_NAMESPACE()

#endif