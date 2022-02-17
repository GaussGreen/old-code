/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: Heston_ModelParams.h,v $
 */

/*! \file ShiftedHeston_ModelParams.h
 *
 *  \brief 
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_SHIFTEDHESTON_MODELPARALS_H
#define _INGPMODELS_SHIFTEDHESTON_MODELPARALS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "AnalyticModelParams.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Heston_ModelParams
// \brief Class for the model param of the Heston Model
//-----------------------------------------------------------------------------
class ARM_ShiftedHeston_ModelParams : public ARM_AnalyticModelParams
{
private:
	void ValidateModelParams() const;

public:
	ARM_ShiftedHeston_ModelParams ( const ARM_ShiftedHeston_ModelParams & rhs );
	ARM_ShiftedHeston_ModelParams ( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_ShiftedHeston_ModelParams ();
    ARM_ShiftedHeston_ModelParams & operator = (const ARM_ShiftedHeston_ModelParams & rhs);

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

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

