/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsCEV_Fx.h
 *
 *  \brief 
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2006
 */


#ifndef _INGPMODELS_CEVMODELPARAMS_FX_H
#define _INGPMODELS_CEVMODELPARAMS_FX_H

// gpbase
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/assignop.h"

// gpinfra
#include "gpinfra/modelparamtype.h"

// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsHW1F.h"

#include "glob/dates.h"
#include "crv/zerocurv.h"
#include "ccy/currency.h"


CC_BEGIN_NAMESPACE( ARM )
 

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsCEV_Fx
// \brief Class for model parameters of CEV Fx
//-----------------------------------------------------------------------------

class ARM_ModelParamsCEV_Fx : public ARM_ModelParams_Fx, public ARM_ModelParamsHW1FStd
{
private:
	void ValidateModelParams() const;

public:
	ARM_ModelParamsCEV_Fx( const ARM_ModelParamVector& params=ARM_ModelParamVector(), ARM_ZeroCurvePtr domCurve=ARM_ZeroCurvePtr(NULL), ARM_ZeroCurvePtr fgnCurve=ARM_ZeroCurvePtr(NULL), double spot=1.0 );
	ARM_ModelParamsCEV_Fx( const ARM_ModelParamsCEV_Fx& rhs );
	ASSIGN_OPERATOR(ARM_ModelParamsCEV_Fx)
	virtual ~ARM_ModelParamsCEV_Fx();
    

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsCEV_Fx(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif