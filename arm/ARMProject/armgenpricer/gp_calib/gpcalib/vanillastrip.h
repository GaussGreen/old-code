/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillastriparg.h
 *
 *  \brief vanilla strip arg
 *
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date November 2006
 */


#ifndef _INGPCALIB_VANILLASTRIPARG_H
#define _INGPCALIB_VANILLASTRIPARG_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/typedef.h"
#include "gpbase/assignop.h"

/// gpinfra
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "vanillaarg.h"

CC_BEGIN_NAMESPACE( ARM )

///////////////////////////////////////////////////////////////
/// \struct ARM_VanillaStripArg
/// \brief base class for vanilla instrument
///////////////////////////////////////////////////////////////
struct ARM_VanillaStripArg: public ARM_VanillaArg
{
	ARM_VanillaStripArg( 
		const ARM_GP_T_Vector<ARM_VanillaArg*>& vanillas,
		const std::vector<double>& coeffs)
	:
		ARM_VanillaArg(),
		itsVanillas(vanillas),
		itsCoeffs(coeffs)
	{
        itsAccountingPricingFlag = ARM_DISC_PRICING_METH;
    }

    ARM_VanillaStripArg(const ARM_VanillaStripArg& arg);
    ASSIGN_OPERATOR(ARM_VanillaStripArg)
		virtual ~ARM_VanillaStripArg() {};
    
	virtual double Price(ARM_PricingModel* model) const;
    virtual double ImpliedVol(ARM_PricingModel* model) const;

	virtual VanillaArgType GetType() const { return ARM_VanillaArg::VANILLA_STRIP; }

    void SetDiscPricingMode(int discPricingMode)
    {
        itsAccountingPricingFlag = discPricingMode;
    }

    int GetDiscPricingMode(void) const
    {
        return(itsAccountingPricingFlag);
    }

	virtual ARM_Object* Clone() const { return new ARM_VanillaStripArg(*this); }
    virtual string toString(const string& indent="", const string& nextIndent="") const;

private :

    ARM_GP_T_Vector<ARM_VanillaArg*> itsVanillas;
	std::vector<double> itsCoeffs;

    int itsAccountingPricingFlag;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
