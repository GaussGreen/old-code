/*!
 *
 * Copyright (c) IXIS CIB July 2005 Paris
 *
 *	\file Local_SLN_ModelParams.h
 *
 *  \brief model params for local normal model
 *	\author  J-M Prié
 *	\version 1.0
 *	\date July 2005
 */



#ifndef _INGPMODELS_LOCAL_SLN_MODELPARAMS_H
#define _INGPMODELS_LOCAL_SLN_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Local_SLN_ModelParams
// \brief Class for the model param of the Shifted LogNormal Model
//-----------------------------------------------------------------------------
class ARM_Local_SLN_ModelParams : public ARM_ModelParams
{
private:
	void ValidateModelParams() const;

public:
	ARM_Local_SLN_ModelParams( const ARM_Local_SLN_ModelParams& rhs );
	ARM_Local_SLN_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_Local_SLN_ModelParams();
    ARM_Local_SLN_ModelParams& operator = (const ARM_Local_SLN_ModelParams& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 0; };

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb=0) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
