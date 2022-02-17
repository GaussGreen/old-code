/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: Merton_ModelParams.h,v $
 * Revision 1.1  2004/07/30 09:52:19  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file Merton_ModelParams.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPMODELS_MERTON_MODELPARAMS_H
#define _INGPMODELS_MERTON_MODELPARAMS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Merton_ModelParams
// \brief Class for the model param of the Merton Model
//-----------------------------------------------------------------------------
class ARM_Merton_ModelParams : public ARM_AnalyticModelParams
{
private:
	void ValidateModelParams() const;

public:
	ARM_Merton_ModelParams( const ARM_Merton_ModelParams& rhs );
	ARM_Merton_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_Merton_ModelParams();
    ARM_Merton_ModelParams& operator = (const ARM_Merton_ModelParams& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 2; };

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
