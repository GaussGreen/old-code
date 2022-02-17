/*!
 * Copyright (c) NATIXIS  April 2007 Paris
 *
 *	\file infCoupon.h
 *	\author  Francois Poitou
 */

#ifndef _INFPAYOFF_H
#define _INFPAYOFF_H

/// gpbase
#include <gpbase/countedptr.h>				//ARM_CountedPtr
#include "gpinflation/floatingCoupon.h"
/// kernel
/// gpinflation


CC_BEGIN_NAMESPACE( ARM )

/*question : où met-on le strike ?*/

class YOYPayOff {
public :
	YOYPayOff(){}
	virtual ~YOYPayOff(){}
	double operator()(const CashFlow& cf)
	{
		double value = 0.;
		return value ;
	}
	
};
class YOYCallPayOff {
public :
	YOYCallPayOff(){}
	virtual ~YOYCallPayOff(){}
	double operator()(const CashFlow& cf)
	{
		double value = 0.;
		return value;
	}
	
};
class YOYPutPayOff {
public :
	YOYPutPayOff(){}
	virtual ~YOYPutPayOff(){}
	double operator()(const CashFlow& cf)
	{
		double value = 0.;
		return value;
	}
	
};
class YOYDigitalCallPayOff {
public :
	YOYDigitalCallPayOff(){}
	virtual ~YOYDigitalCallPayOff(){}
	double operator()(const CashFlow& cf)
	{
		double value = 0.;
		return value ;
	}
	
};
class YOYDigitalPutPayOff {
public :
	YOYDigitalPutPayOff(){}
	virtual ~YOYDigitalPutPayOff(){}
	double operator()(const CashFlow& cf)
	{
		double value = 0.;
		return value ;
	}
	
};

CC_END_NAMESPACE()

#endif
