/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: bondanalytics.h,v $
 * Revision 1.1  2004/03/01 07:52:17  ebenhamou
 * Initial revision
 *
 */


/*! \file bondanalytics.h
 *
 *  \brief
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date May 2004
 */


#ifndef _INGPCLOSEDFORMS_BONDANALYTIC_H
#define _INGPCLOSEDFORMS_BONDANALYTIC_H

#include "gpbase/port.h"
#include "gpbase/env.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_BondAnalytics
{
	static double PriceToYield( 
		double coupon,
		double redemptionValue,
		double price,
		int frequency,
		double settlToNext,
		double settlToMaturity,
		double accruingTime );

	static double YieldToPrice(
		double coupon,
		double redemptionValue,
		double yield,
		int frequency,
		double settlToNext,
		double settlToMaturity,
		double accruingTime );
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

