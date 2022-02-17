/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ModelParamsQNF.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#ifndef _INGPMODELS_MODELPARAMSQNF_H
#define _INGPMODELS_MODELPARAMSQNF_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "ModelParamsHW1F.h"
#include <vector>

CC_USING_NS(std, vector)

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsQNF
// \brief Class for model parameters of NF QModel
//-----------------------------------------------------------------------------
class ARM_ModelParamsQNF : public ARM_ModelParams
{
private:
	void ValidateModelParams() const;
	vector<ARM_ModelParamsHW1F*> itsFactorParams;
	ARM_CurveModelParam itsQParam;
	ARM_GP_Matrix itsCorrelMatrix;

public:
	/// constructor
	ARM_ModelParamsQNF( const ARM_CurveModelParam& QParam, const vector<ARM_ModelParamsHW1F*>& inputs, const ARM_GP_Matrix& correlMatrix);
	ARM_ModelParamsQNF( const ARM_ModelParamsQNF& rhs );
	virtual ~ARM_ModelParamsQNF();
    ARM_ModelParamsQNF& operator = (const ARM_ModelParamsQNF& rhs);

	/// accessors (rference to avoid copy, but beware that once the object is destroyed, these go away as well!
	const ARM_ModelParamsHW1F* GetModelParamsPerDim( size_t dim ) const;
	inline const ARM_CurveModelParam& GetQParam() const { return itsQParam; }
	inline const ARM_GP_Matrix& GetCorrelMatrix() const { return itsCorrelMatrix; }

	/// How many factors?
    virtual size_t FactorCount() const { return itsFactorParams.size(); };
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {}
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

