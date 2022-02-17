/*!
 *
 * Copyright (c) IXIS CIB January 2005 Paris
 *
 *	\file SFRMPricingStates.h
 *
 *  \brief SFRMPricingStates
 *	\author  R.Guillemot, A.Schauly
 *	\version 1.0
 *	\date January 2005
 */

#ifndef _INGPMODELS_PRINCINGSTATESSFRM_H
#define _INGPMODELS_PRINCINGSTATESSFRM_H

#include "gpinfra/pricingstates.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_SFRMPricingStatesContext : public ARM_PricingStatesContext
{
private: 
	int itsFwdMinIndex;
	int itsNumeraireTimeIndex;
public: 
	ARM_SFRMPricingStatesContext();
	ARM_SFRMPricingStatesContext( const ARM_SFRMPricingStatesContext& rhs );
	ARM_SFRMPricingStatesContext& operator=( const ARM_SFRMPricingStatesContext& rhs );

	virtual ARM_SFRMPricingStatesContext * ToSFRMPricingStatesContext();

	virtual ARM_Object* Clone() const;
	virtual ~ARM_SFRMPricingStatesContext();

	inline int GetNumeraireTimeIndex() const { return itsNumeraireTimeIndex; }
	inline void SetNumeraireTimeIndex( int NumeraireTimeIndex ) { itsNumeraireTimeIndex = NumeraireTimeIndex; }
	inline int GetFwdMinIndex() const { return itsFwdMinIndex; }
	inline void SetFwdMinIndex( int FwdMinIndex ) { itsFwdMinIndex = FwdMinIndex; }

	virtual string toString(const string& indent="", const string& nextIndent="") const { return "ARM_SFRMPricingStatesContext"; }

};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/