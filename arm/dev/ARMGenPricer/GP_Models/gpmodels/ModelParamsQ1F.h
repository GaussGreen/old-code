/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQModel.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPMODELS_MODELPARAMSQ1F_H
#define _INGPMODELS_MODELPARAMSQ1F_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparamtype.h"
#include "ModelParamsHW1F.h"


CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsQ1F
// \brief Class for model parameters of 1F QModel
//-----------------------------------------------------------------------------
class ARM_ModelParamsQ1F : public ARM_ModelParamsHW1FStd
{
private:
	void ValidateModelParams() const;

public:
	ARM_ModelParamsQ1F( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	ARM_ModelParamsQ1F( const ARM_ModelParamsQ1F& rhs );
	virtual ~ARM_ModelParamsQ1F();
    ARM_ModelParamsQ1F& operator = (const ARM_ModelParamsQ1F& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;

    size_t GetModelCurves(ARM_GP_Vector& times, ARM_GP_Vector& sigmas, ARM_GP_Vector& qs, ARM_GP_Vector& dqs) const;

	static double Q1FStateCovariance( const ARM_ModelParamsQ1F* lhs, const ARM_ModelParamsQ1F* rhs, double a, double b, double T, bool isrhsQ=true );
    static double Q1FStateZcCovariance( const ARM_ModelParamsQ1F* lhs, const ARM_ModelParamsHW1F* rhs, double a, double b, double T, double U );

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

