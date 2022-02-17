/*
 *
 * Copyright (c) CDC IXIS CM July 2004 Paris
 *
 * $Log: infswaption.h,v $
 * Revision 1.1  2004/09/01 18:53:24  ebenhamou
 *  filtering of blank
 *
 *
 */


/*! \file infswaption.h
 *
 *  \brief the inflation swaption object is a European inflation swaption
 *			that is based on a swap that contains an inflation leg at least!
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPINFLATION_INFSWAPTION_H
#define _INGPINFLATION_INFSWAPTION_H

#include "gpbase/port.h"
#include <inst/swaption.h>
#include <string>
#include "gpbase/env.h"

class ARM_Swap;
class ARM_Date;


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class ARM_TwoAssetsInfo;

/// typedef for better naming
typedef ARM_TwoAssetsInfo ARM_InfSwaptionContext;



/// the ARM_InfSwaption handles only European case
class ARM_InfSwaption : public ARM_Swaption
{
private:
	/// pricing element for view
	ARM_InfSwaptionContext* itsAssetsInfo;
	ARM_SwapLeg* itsFirstInfLeg;		/// simple pointor to the leg... no creation no destruction!
	ARM_SwapLeg* itsSecondLeg;			/// simple pointor to the leg... no creation no destruction!

public:
	ARM_InfSwaption( ARM_Swap* swap, int rcvPay, double strike = K_MARKET_RATE, const ARM_Date& maturity= (ARM_Date) ARM_DEFAULT_DATE );
	ARM_InfSwaption( const ARM_InfSwaption& rhs );
	ARM_InfSwaption& operator=( const ARM_InfSwaption& rhs );
	virtual ~ARM_InfSwaption();

	/// for easy debugging
	virtual CC_NS( std, string ) toString() const;

	/// ARM_Pricing support
	virtual double ComputePrice(int mode = 0);	
	virtual void CptCashFlowValues();

	
	/// ARM_Object support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	/// static function to determine inflation swaption!
	static size_t NbOfInflationLeg( ARM_Swap* swap );
	inline ARM_InfSwaptionContext* GetAssetsInfo() { return itsAssetsInfo; }
};

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
