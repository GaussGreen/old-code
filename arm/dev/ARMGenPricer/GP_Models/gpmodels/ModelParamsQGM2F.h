/*!
 *
 * Copyright (c) CDC IXIS CM May 2006 Paris
 *
 *	\file ModelParamsQGM2F.cpp
 *
 *  \brief 
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date May 2006
 */

#ifndef _INGPMODELS_MODELPARAMSQGM2F_H
#define _INGPMODELS_MODELPARAMSQGM2F_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamsvec.h"


CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsQGM2FOneDim
// \brief Class for model parameters of each factor of the 2Factor Quadratic 
// \Gaussian Model
//-----------------------------------------------------------------------------
class ARM_ModelParamsQGM2FOneDim : public ARM_ModelParams 
{
public:
	ARM_ModelParamsQGM2FOneDim( const ARM_ModelParamsQGM2FOneDim& rhs );
	ARM_ModelParamsQGM2FOneDim( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ModelParamsQGM2FOneDim();
    ARM_ModelParamsQGM2FOneDim& operator = (const ARM_ModelParamsQGM2FOneDim& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 1; };
   
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0);
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb=0) {};

	// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


//-----------------------------------------------------------------------------
// \class ARM_ModelParamsQGM2F
// \brief Class for model parameters of 1F Quadratic Gaussian Model
//-----------------------------------------------------------------------------
class ARM_ModelParamsQGM2F : public ARM_ModelParamsVec 
{
public:
	ARM_ModelParamsQGM2F( const ARM_ModelParamsQGM2F& rhs );
	ARM_ModelParamsQGM2F( const vector<ARM_ModelParams*>& paramsVec  );
	virtual ~ARM_ModelParamsQGM2F();
    ARM_ModelParamsQGM2F& operator = (const ARM_ModelParamsQGM2F& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 2; };

   	// Standard ARM object support
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

