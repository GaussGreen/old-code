/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1FAna_ModelParams.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPMODELS_Q1FANA_MODELPARAMS_H
#define _INGPMODELS_Q1FANA_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Q1FAna_ModelParams
// \brief Class for model parameters of 1F QModel
//-----------------------------------------------------------------------------
class ARM_Q1FAna_ModelParams : public ARM_ModelParams 
{
private:
	void ValidateModelParams() const;
	CC_IS_MUTABLE ARM_ModelParamType::ParamNb itsVolType;
	void CalibrateQVol( double forward, double maturity, int CallPut );

public:
	ARM_Q1FAna_ModelParams( const ARM_Q1FAna_ModelParams& rhs );
	ARM_Q1FAna_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_Q1FAna_ModelParams();
    ARM_Q1FAna_ModelParams& operator = (const ARM_Q1FAna_ModelParams& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Q model specific
	void CalibrateModelParams( double forward, double maturity, int CallPut );

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

