/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsHeston_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_MODELPARAMSHeston_FX_H
#define _INGPMODELS_MODELPARAMSHeston_FX_H

/// gpbase
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"

/// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/Heston_ModelParams.h"

#include "crv/zerocurv.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_ModelParamsHeston_Fx : public ARM_ModelParams_Fx, public ARM_Heston_ModelParams
{
protected:
	void ValidateModelParams() const;

public:
	ARM_ModelParamsHeston_Fx( const ARM_ModelParamVector& params=ARM_ModelParamVector(), ARM_ZeroCurvePtr domCurve=ARM_ZeroCurvePtr(NULL), ARM_ZeroCurvePtr fgnCurve=ARM_ZeroCurvePtr(NULL), double spot=1.0 );
	ARM_ModelParamsHeston_Fx( const ARM_ModelParamsHeston_Fx& rhs );
	ASSIGN_OPERATOR(ARM_ModelParamsHeston_Fx)
	virtual ~ARM_ModelParamsHeston_Fx();

	/// How many factors?
    virtual size_t FactorCount() const { return 3; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsHeston_Fx(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
