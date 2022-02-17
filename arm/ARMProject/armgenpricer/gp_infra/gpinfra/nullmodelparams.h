/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file ARM_NullModelParams.h
 *
 *  \brief 
 *	\author  A. Chaix
 *	\version 1.0
 *	\date July 2005
 */


#ifndef _INGPINFRA_NULLMODELPARAMS_H
#define _INGPINFRA_NULLMODELPARAMS_H


#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_NullModelParams
// \brief This is a default instantiable class for ARM_ModelParams
//-----------------------------------------------------------------------------
class ARM_NullModelParams : public ARM_ModelParams
{
public:
	ARM_NullModelParams() {};
	ARM_NullModelParams( const ARM_NullModelParams& rhs ) : ARM_ModelParams (rhs) {};
	virtual ~ARM_NullModelParams(){};
    ARM_NullModelParams& operator = (const ARM_NullModelParams& rhs) {ARM_ModelParams::operator =(rhs);}
	virtual ARM_Object* Clone() const {return new ARM_NullModelParams(*this);}
    
	/// calibration part
	/// mandatory because pure virtual in ARM_ModelParams 
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb = 0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model ,int factorNb = 0) {};
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
