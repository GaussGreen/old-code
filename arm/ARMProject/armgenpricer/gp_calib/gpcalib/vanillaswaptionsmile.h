/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file vanillasecuritydensity.h
 *
 *  \brief vanilla security density
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPCALIB_VANILLASWAPTIONSMILE_H
#define _INGPCALIB_VANILLASWAPTIONSMILE_H

/// this header has to come first

#include "gpbase/env.h"
#include "gpbase/rootobject.h"
#include "gpinfra/typedef.h"
#include "gpcalib/typedef.h"
#include "vanillaswaption.h"

CC_BEGIN_NAMESPACE( ARM )

struct ARM_VanillaSmiledSwaption : public ARM_VanillaSwaptionArg
{
public:
	/// Constructors and Destructors
	ARM_VanillaSmiledSwaption(const ARM_VanillaSwaptionArg& refSwaption,
							const std::vector<double>& target) ;

	ARM_VanillaSmiledSwaption(const ARM_VanillaSmiledSwaption& rhs);
	virtual ~ARM_VanillaSmiledSwaption();
	
	/// Standard ARM Object Support
	virtual ARM_Object* Clone() const {return new ARM_VanillaSmiledSwaption(*this);};
	virtual double Price(ARM_PricingModel* model) const;
	virtual string ExportShortName() const { return "LVSWS";}
	virtual string toString(const string& indent="",const string& nextIndent="") const ;

private: 
	std::vector<double> itsTargetParams;
};


CC_END_NAMESPACE()
#endif