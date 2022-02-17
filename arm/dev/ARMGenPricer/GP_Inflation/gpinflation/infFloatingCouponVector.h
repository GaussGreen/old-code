/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file infFloatinCouponVector.h
 *	\author  Francois Poitou
 */

#ifndef _INFCOUPONVECT_H

#define _INFCOUPONVECT_H

/// gpinflation
#include "gpinflation/infFloatingCoupon.h"
///gpbase
#include "gpbase/datestrip.h"			/// for ARM_DateStrip

CC_BEGIN_NAMESPACE( ARM )

template<class payOffType>
ARM_GP_T_Vector<ARM_CountedPtr < CashFlow > >
	InfFloatingCouponVector(	const ARM_GP_Vector& nominals,
								ARM_CountedPtr<payOffType> payoff,
								ARM_CountedPtr<InfPerformance> cpiPerf,
							    const ARM_Date& startDate,
							    const ARM_Date& endDate,
							    int resetFreq,
							    int payFreq,
							    int dayCount,
							    int fwdRule,
							    int stubRule,
							    int intRule,
							    int resetRollingConv,
								char* resetCal,
								int payRollingConv,
								char* payCal,
								const bool& monthBegin,
							    const Period& demGap,
								const Period& numGap)
{
	ARM_DateStrip stripTease( startDate, endDate, payFreq,
						dayCount, resetCal,
						fwdRule, intRule, stubRule,
						0., payFreq, 0., payCal,
						K_ADVANCE, K_ARREARS, false );
	ARM_GP_T_Vector<ARM_CountedPtr < CashFlow > > v ;

	const ARM_GP_Vector&  vEnd		= * stripTease.GetFlowEndDates() ;
	const ARM_GP_Vector&  vStart	= * stripTease.GetFlowStartDates() ;

	int nFlows = stripTease.GetPaymentDates()->size() ;
	ARM_Date sDate, eDate ;
	ARM_GP_Vector  vNominals ;
	if (nominals.size() == 1)
		vNominals = ARM_GP_Vector(nFlows, nominals[0]);
	else if (nominals.size() != nFlows)
		ARM_THROW( ERR_INVALID_ARGUMENT, "wrong size of nominal amount vector");
	double nominal ;


	//ZC CASE
	if (payFreq == K_ZEROCOUPON)
	{
		if (resetFreq != K_ZEROCOUPON)  ARM_THROW( ERR_INVALID_ARGUMENT, "");
		sDate 	= ARM_Date(vStart[0]);
		eDate	= vEnd[0];
		nominal = vNominals[0] ;
		v.push_back(	ARM_CountedPtr < InfFloatingCoupon<payOffType> > (
						new InfFloatingCoupon<payOffType> (	nominal, payoff, cpiPerf, sDate, eDate, 
															demGap, numGap,	dayCount, resetRollingConv, resetCal, 
															payRollingConv, payCal, monthBegin 	)));

	}
	else
	{
		//OAT CASE
		if (resetFreq == K_ZEROCOUPON)
		{
			for( int i=0; i<nFlows; ++i )
			{
				sDate 		= startDate;
				eDate		= ARM_Date( vEnd[i]);
				nominal 	= vNominals[i] ;
				v.push_back(	ARM_CountedPtr < InfFloatingCoupon<payOffType> > (
								new InfFloatingCoupon<payOffType> (	nominal, payoff, cpiPerf, sDate, eDate, 
																	demGap, numGap, dayCount, resetRollingConv, resetCal, 
																	payRollingConv, payCal, monthBegin 	)));

			}
		}
		else  // YOY and OTHER CASES
		{
			// YOY CASE : useless
			/*
			if (resetFreq == payFreq)
			{
				for( int i=0; i<nFlows; ++i )
				{
					sDate 		= vStart[i];
					eDate		= ARM_Date( vEnd[i]);
					nominal 	= vNominals[i] ;
					v.push_back(	ARM_CountedPtr < InfFloatingCoupon<payOffType> > (
									new InfFloatingCoupon<payOffType> (	nominal, payoff, cpiPerf, sDate, eDate, 
																		demGap, numGap,	dayCount, resetRollingConv, resetCal, 
																		payRollingConv, payCal, monthBegin 	)));

				}
			}
			
			else // CORRIDOR or ASIAN CASE
			{
			*/
				for( int i=0; i<nFlows; ++i )
				{
					sDate 		= vStart[i];
					eDate		= ARM_Date( vEnd[i]);
					nominal 	= vNominals[i] ;
					v.push_back(	ARM_CountedPtr < InfFloatingCoupon<payOffType> > (
									new InfComposedCoupon<payOffType> (	nominal, payoff, cpiPerf, sDate, eDate, resetFreq,
															dayCount, demGap, numGap, resetRollingConv, resetCal, 
															payRollingConv, payCal, monthBegin  )));
				}
			/*}*/
		}
	}
	return v ;


}

CC_END_NAMESPACE()

#endif
