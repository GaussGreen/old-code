/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *! \file datestripconvert.cpp
 *
 *  \brief file to convert curve into a refvalue!
 *	\author  E.Ezzine
 *	\version 1.0
 *	\date June 2006
 */

#include "gpbase/datestripconvert.h"
#include "gpbase/gplinalgconvert.h"
#include "gpbase/datestrip.h"

CC_BEGIN_NAMESPACE( ARM )

// Function to convert a curve into a reference value
ARM_DateStrip* SwapLegToDateStrip( const ARM_SwapLeg& swapleg)
{
	ARM_GP_Vector resetDates	(To_ARM_GP_Vector(const_cast<ARM_SwapLeg&>(swapleg).GetResetDates()));
	ARM_GP_Vector flowStartDates (To_ARM_GP_Vector(swapleg.GetFlowStartDates()));
	ARM_GP_Vector flowEndDates   ( To_ARM_GP_Vector(swapleg.GetFlowEndDates()));
	ARM_GP_Vector fwdStartDates	( To_ARM_GP_Vector(swapleg.GetFwdRateStartDates()));
	ARM_GP_Vector fwdEndDates		( To_ARM_GP_Vector(swapleg.GetFwdRateEndDates()));
	ARM_GP_Vector payDates			( To_ARM_GP_Vector(swapleg.GetFwdRateEndDates()));
	ARM_GP_Vector interestDays		( To_ARM_GP_Vector(swapleg.GetInterestDays()));
	ARM_GP_Vector interestTerms	( To_ARM_GP_Vector(swapleg.GetInterestTerms()));

ARM_DateStrip*  datestrip = new ARM_DateStrip( &flowStartDates,	
											  &flowEndDates,		
											  &fwdStartDates,	
											  &fwdEndDates,		
											  &resetDates,			
											  &payDates,			
											  &interestDays,			
											  &interestTerms );			



	 return datestrip;
	
}



CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

