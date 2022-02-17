/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file fixzc.h
 *  \brief implements the Fix ZC class which is a very basic leg for inflation swap
 *		Inherit from ARM_SwapLeg
 *
 *	\author  Eric Benhamou
 *	\version 1.0
 *	\date August 2003
 */


#ifndef _INGPINFLATION_FIXZC_H
#define _INGPINFLATION_FIXZC_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/gplinalgtypedef.h"

#include "swapleg.h"


/// ARM namespace
CC_BEGIN_NAMESPACE( ARM )


/*!
 * \class ARM_FixZC
 * \author  Eric Benhamou
 * \version 1.0
 * \date August 2003
 *
 * \brief class for a fix Zero Coupon leg
 * declared as a thin wrapper
 * around ARM_SWAPLEG
 * because there is nothing more than the ARM_SwapLeg
 * the standard function like bitwisecopy and so on
 * are unecessary!
 */
 
class ARM_FixZC: public ARM_SwapLeg
{
private:
	ARM_GP_Vector* itsDiscountFactors;

public:
	ARM_FixZC(
		const ARM_Date& startDate,
		const ARM_Date& endDate,
		double fixRate,
		int rcvOrPay,
		int dayCount,
		int intRule,
		int stubRule,
		int payGap,					/// pay gap default is 0
		const char* payCalendar,	/// calendar used for payment
		const ARM_Currency*	discountCcy = NULL );

	virtual ~ARM_FixZC();
	ARM_FixZC( const ARM_FixZC& rhs);
	ARM_FixZC& operator=( const ARM_FixZC& rhs);


	/// ARM pricing support
	virtual double ComputePrice(int mode = 0);		
	virtual void PropagateModel(ARM_Model *model);
	virtual void CptCashFlowValues(void);

	/// ARM_Object standard support
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
	virtual ARM_Object* Clone();

};

CC_END_NAMESPACE()


#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/


