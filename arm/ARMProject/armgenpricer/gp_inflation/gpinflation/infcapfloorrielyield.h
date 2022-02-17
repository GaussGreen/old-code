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


#ifndef _INGPINFLATION_INFCAPFLOORRIELYIELD_H
#define _INGPINFLATION_INFCAPFLOORRIELYIELD_H

#define  LA_MID_DATE	6
#define  LA_END_DATE	12
#define	 EPSILON		1e-6


#include "gpbase/port.h"
#include <inst/swaption.h>
#include <string>
#include "gpbase/env.h"
#include "gpbase/typedef.h"

class ARM_Swap;
class ARM_Date;
class ARM_SwapLeg;


CC_BEGIN_NAMESPACE( ARM )

class ARM_InfBSModel;
class ARM_InfLeg;
class ARM_TwoAssetsInfo;


/// the ARM_InfSwaption handles only European case
class ARM_InfCapFloorRielYield : public ARM_Swaption
{
private:

	ARM_TwoAssetsInfo* itsAssetsInfo;
	/// pricing element for view
	ARM_SwapLeg* itsFirstInfLeg;		/// simple pointor to the leg... no creation no destruction!
	ARM_SwapLeg* itsSecondLeg;			/// simple pointor to the leg... no creation no destruction!

	ARM_VectorPtr itsPricingStrikesForView;
	ARM_VectorPtr itsInflationVolsForView;
	ARM_VectorPtr itsSecondLegVolsForView;
	ARM_VectorPtr itsCorrelationsForView;

public:
	ARM_InfCapFloorRielYield( ARM_Swap* swap, int CoF, double strike = 0);
	ARM_InfCapFloorRielYield( ARM_Swap* swap, int CoF, ARM_ReferenceValue *strike);
	ARM_InfCapFloorRielYield( const ARM_InfCapFloorRielYield& rhs );
	ARM_InfCapFloorRielYield& operator=( const ARM_InfCapFloorRielYield& rhs );
	
	virtual ~ARM_InfCapFloorRielYield();

	/// for easy debugging
	virtual CC_NS( std, string ) toString() const;

	/// ARM_Pricing support
	virtual double ComputePrice(int mode = 0);	
	virtual void CptCashFlowValues();

	/// ARM_Object support
	virtual ARM_Object* Clone();
	virtual void View(char* id = NULL, FILE* ficOut = NULL);

	inline ARM_TwoAssetsInfo* GetAssetsInfo() { return itsAssetsInfo; }


private:
	
void	ComputeDataFromSwapLeg(	ARM_InfBSModel*			pricingMod,
								ARM_SwapLeg*			irLeg,
								ARM_InfLeg*				infLeg,
								ARM_GP_VectorPtr		InfIntTerm,
								ARM_GP_VectorPtr		InfFwdRate,
								ARM_GP_VectorPtr		IrFwdRate,
								ARM_GP_VectorPtr		InfSpread,
								ARM_GP_VectorPtr		IrResTerm,
								ARM_GP_VectorPtr		IrIntTerm,
								ARM_GP_VectorPtr		IrVolTerm,
								ARM_GP_VectorPtr		IrLevVolTerm,
								ARM_GP_VectorPtr		PriceStrike,
								ARM_GP_VectorPtr		DiscFact	);

void	ComputeDataFromLivretA(	ARM_InfBSModel*			pricingMod,
								ARM_SwapLeg*			irLeg,
								ARM_InfLeg*				infLeg,									   
								ARM_GP_VectorPtr		InfIntTerm,
								ARM_GP_VectorPtr		InfFwdRate,
								ARM_GP_VectorPtr		IrFwdRate,
								ARM_GP_VectorPtr		InfSpread,
								ARM_GP_VectorPtr		IrResTerm,
								ARM_GP_VectorPtr		IrIntTerm,
								ARM_GP_VectorPtr		IrVolTerm,
								ARM_GP_VectorPtr		IrLevVolTerm,
								ARM_GP_VectorPtr		PriceStrike,
								ARM_GP_VectorPtr		DiscFact	);

};



class ARM_InfCallSpreadYield : public ARM_InfCapFloorRielYield
{

public:
	ARM_InfCallSpreadYield	(	 ARM_Swap*	swap, int	CoF, double					strike = 0)
		:ARM_InfCapFloorRielYield(			swap,		CoF,						strike){}

	ARM_InfCallSpreadYield	( ARM_Swap*		swap, int	CoF, ARM_ReferenceValue*	strike)
		:ARM_InfCapFloorRielYield(			swap,		CoF,						strike){}

	ARM_InfCallSpreadYield	( const ARM_InfCallSpreadYield& rhs ):ARM_InfCapFloorRielYield( rhs ){ }
	
	virtual ~ARM_InfCallSpreadYield(){}

	/// ARM_Pricing support
	virtual double ComputePrice(int mode = 0);	

	/// ARM_Object support
	virtual ARM_Object* Clone(){ return new ARM_InfCallSpreadYield(*this); }
};


CC_END_NAMESPACE()

#endif


/*----------------------------------------------------------------------------*/
/*---- End of file ----*/
