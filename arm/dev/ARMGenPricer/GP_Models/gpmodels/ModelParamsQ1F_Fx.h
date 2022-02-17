/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQ1F_Fx.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */


#ifndef _INGPMODELS_MODELPARAMSQ1F_FX_H
#define _INGPMODELS_MODELPARAMSQ1F_FX_H

/// gpbase
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"

/// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/ModelParamsQ1F.h"

#include "crv/zerocurv.h"

CC_BEGIN_NAMESPACE( ARM )


class ARM_ModelParamsQ1F_Fx : public ARM_ModelParams_Fx, public ARM_ModelParamsQ1F
{
private:
	void ValidateModelParams() const;

public:
	ARM_ModelParamsQ1F_Fx( const ARM_ModelParamVector& params=ARM_ModelParamVector(), ARM_ZeroCurvePtr domCurve=ARM_ZeroCurvePtr(NULL), ARM_ZeroCurvePtr fgnCurve=ARM_ZeroCurvePtr(NULL), double spot=1.0 );
	ARM_ModelParamsQ1F_Fx( const ARM_ModelParamsQ1F_Fx& rhs );
	ASSIGN_OPERATOR(ARM_ModelParamsQ1F_Fx)
	virtual ~ARM_ModelParamsQ1F_Fx();
    

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_ModelParamsQ1F_Fx(*this); };
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};

struct ARM_Q1FModelParams_FxBuilder
{
	static ARM_ModelParamVector CreateAndValidateModelParams( const ARM_ModelParamVector& params );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
