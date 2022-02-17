/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Local_Normal_ModelParams.h
 *
 *  \brief model params for local normal model
 *	\author  A. Chaix
 *	\version 1.0
 *	\date June 2005
 */



#ifndef _INGPMODELS_LOCAL_NORMAL_MODELPARAMS_H
#define _INGPMODELS_LOCAL_NORMAL_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpbase/gplinalgtypedef.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Local_Normal_ModelParams
// \brief Class for the model param of the Normal Model
//-----------------------------------------------------------------------------
class ARM_Local_Normal_ModelParams : public ARM_ModelParams
{
private:
	void ValidateModelParams(const ARM_IntVector& paramTypes) const;

public:
	ARM_Local_Normal_ModelParams( const ARM_Local_Normal_ModelParams& rhs );
	ARM_Local_Normal_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector(), const ARM_IntVector& paramTypes=ARM_IntVector());
	virtual ~ARM_Local_Normal_ModelParams();
    ARM_Local_Normal_ModelParams& operator = (const ARM_Local_Normal_ModelParams& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 0; };

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
