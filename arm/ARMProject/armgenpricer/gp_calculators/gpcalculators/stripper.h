/*
 * Copyright (c) CDC IXIS CM July 2005 Paris
 *
 *	\file stripper.h
 *
 *  \brief
 *
 *	\author  Y. KHLIF
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPCALCULATORS_STRIPPER_H
#define _INGPCALCULATORS_STRIPPER_H

#include "gpbase/port.h"
#include "gpcalib/typedef.h"

#include "gpbase/rootobject.h"
#include "gpbase/datestripcombiner.h"

CC_BEGIN_NAMESPACE( ARM )



/////////////////////////////////////////////////////////////////
/// \class ARM_STRIPPER 
/// this class constucts a zero coupon curve  
///from the first zero coupon and the forward swap rates  for a scheduler of dates
/// it helps the CSO calculator to caculte the sensitivities of the underlying
/////////////////////////////////////////////////////////////////

class ARM_STRIPPER : public ARM_RootObject
{
public:
	/// constructor, copy constructor, assignment operator and destructor
	ARM_STRIPPER():itsDatestrip(0),itsSwapRates(0),itsMaturities(0),itsZcRates(0){};
	ARM_STRIPPER(double asOfDate, ARM_DateStripPtr datestrip);
	ARM_STRIPPER( const ARM_STRIPPER& rhs );
	ARM_STRIPPER& operator=( const ARM_STRIPPER& rhs );
	virtual ~ARM_STRIPPER();

	virtual string toString(const string& indent="",const string& nextIndent="") const;
	//virtual ARM_CLASS_NAME GetRootName() { return ARM_STRIPPER; }
	virtual ARM_Object* Clone() const;
	

	void setStartDf (double df);
	void setSwapRate   (int i, double swaprate);
	void strip ();
	double df (double maturity) const;
	double swapRate (double startFloatTime,double endFloatTime, const ARM_GP_Vector& FixPayTimes, const ARM_GP_Vector& FixPayPeriods) const;
	double stdFundingLeg  (const ARM_GP_Vector& startFlowTimes, const ARM_GP_Vector& endFlowTimes,  const ARM_GP_Vector& periods, const ARM_GP_Vector& Notional, const ARM_GP_Vector& Margin, ARM_GP_Vector* LiborLeverage = NULL)  const;


private:

	double					itsAsofDate;
	ARM_DateStripPtr		itsDatestrip;
	ARM_GP_VectorPtr		itsSwapRates;
    ARM_GP_VectorPtr		itsMaturities;
	ARM_GP_VectorPtr		itsZcRates;

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

