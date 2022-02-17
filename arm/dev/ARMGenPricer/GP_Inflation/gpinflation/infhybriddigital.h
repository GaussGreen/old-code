/*
 *
 * Copyright (c) CDC IXIS CIB February 2006 Paris
 *
 */


/*! \file infcapfloorrielyield.h
 *
 *  \brief the inflation capfloor on riel yield object is a inflation capfloor
 *			that is based on a riel yield forward flows it contains an inflation leg at and a libor leg! 
 *			each leg fixes on its fixing which gives not necessary a spread option but a path dependent option
 *	\author  Y. Khlif
 *	\version 1.0
 *	\date February 2006
 */


#ifndef _INGPINFLATION_INFHYBRIDDIGITAL_H
#define _INGPINFLATION_INFHYBRIDDIGITAL_H

#include "gpbase/port.h"
#include <inst/swaption.h>
#include <string>
#include "gpbase/env.h"
#include "gpbase/typedef.h"

class ARM_Swap;
class ARM_Date;


CC_BEGIN_NAMESPACE( ARM )

/// the ARM_InfSwaption handles only European case
class ARM_InfHybridDigital : public ARM_Swaption
{
private:

	//ARM_TwoAssetsInfo* itsAssetsInfo;
	/// pricing element for view
	ARM_SwapLeg* itsPayLeg;		/// simple pointor to the leg... no creation no destruction!
	ARM_SwapLeg* itsDigitLeg;			/// simple pointor to the leg... no creation no destruction!

	int itsPayoffType;
	int itsCapOrFloor;

	ARM_VectorPtr itsPricingStrikesForView;
	ARM_VectorPtr itsInflationVolsForView;
	ARM_VectorPtr itsDigitLegVolsForView;
	ARM_VectorPtr itsCorrelationsForView;

public:
	
	ARM_InfHybridDigital( ARM_SwapLeg* payLeg, ARM_SwapLeg* digitLeg, int payOffType, double barrier = 0, int CF= K_CAP, int RecOrPay = K_RCV);
	ARM_InfHybridDigital( ARM_SwapLeg* payLeg, ARM_SwapLeg* digitLeg, int payOffType, ARM_ReferenceValue *barrierProfileint, int CF, int RecOrPay);
	ARM_InfHybridDigital( const ARM_InfHybridDigital& rhs );
	ARM_InfHybridDigital& operator=( const ARM_InfHybridDigital& rhs );
	virtual ~ARM_InfHybridDigital();

	/// for easy debugging
	virtual CC_NS( std, string ) toString() const;

	/// ARM_Pricing support
	virtual double ComputePrice(int mode = 0);	
	virtual void CptCashFlowValues();

	inline int GetPayoffType() const { return itsPayoffType; }
	inline void SetPayoffType(int payoffType) { itsPayoffType = payoffType; }

	inline int GetCapOrFloor() const { return itsCapOrFloor; }
	inline void SetCapOrFloor(int capOrFloor) { itsCapOrFloor = capOrFloor; }


	/// ARM_Object support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);
};

CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
