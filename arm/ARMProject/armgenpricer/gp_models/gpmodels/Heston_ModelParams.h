/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: Heston_ModelParams.h,v $
 * Revision 1.1  2004/09/22 14:45:48  ebenhamou, ocroissant
 * Initial revision
 *
 *
 */

/*! \file Heston_ModelParams.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou, O Croissant
 *	\version 1.0
 *	\date September 2004
 */


#ifndef _INGPMODELS_HESTON_MODELPARALS_H
#define _INGPMODELS_HESTON_MODELPARALS_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpinfra/modelparams.h"
#include "ModelParamsHW1F.h"

CC_BEGIN_NAMESPACE( ARM )

//-----------------------------------------------------------------------------
// \class ARM_Heston_ModelParams
// \brief Class for the model param of the Heston Model
//-----------------------------------------------------------------------------

class ARM_Heston_ModelParams : public ARM_ModelParamsHW1FStd
{
protected:
	void ValidateModelParams() const;

public:
	ARM_Heston_ModelParams( const ARM_Heston_ModelParams& rhs );
	ARM_Heston_ModelParams( const ARM_ModelParamVector& params=ARM_ModelParamVector() );
	virtual ~ARM_Heston_ModelParams();
    ARM_Heston_ModelParams& operator = (const ARM_Heston_ModelParams& rhs);

	/// How many factors?
    virtual size_t FactorCount() const { return 4; };

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

	/// Standard ARM object support
	virtual ARM_Object* Clone() const { return new ARM_Heston_ModelParams(*this);  }

	ARM_CurveModelParam GetVolAlpha() const;

	double Scaling(double t) const;
	double InitialVol() const;
	double MC_scheme() const;
	double LongTermVol(double t) const;
	double VolOfVol(double t) const;
	double VolMeanReversion(double t) const;
	double Correlation(double t) const;
	double Beta(double t) const;
	double Alpha() const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

