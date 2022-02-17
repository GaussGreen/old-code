/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *


/*! \file ModelParamsHWSV.h
 *
 *  \brief 
 *
 *	\author J-M Prié
 *	\version 1.0
 *	\date September 2006
 */


#ifndef _INGPMODELS_MODELPARAMSHWSV_H
#define _INGPMODELS_MODELPARAMSHWSV_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/gpvector.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"		/// for ARM_ModelParamType
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/modelparamshw2f.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;
class ARM_ModelFitter;

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsHWSV
// \brief Interface class for model parameters of the Hull & White SV models
//-----------------------------------------------------------------------------
class ARM_ModelParamsHWSV : public ARM_ModelParams
{
protected :
	mutable std::vector<double> itsSchedule;

public:
	ARM_ModelParamsHWSV( const ARM_ModelParamsHWSV& rhs );
	ARM_ModelParamsHWSV( const ARM_ModelParamVector& params=ARM_ModelParamVector());
	virtual ~ARM_ModelParamsHWSV();
    ARM_ModelParamsHWSV& operator = (const ARM_ModelParamsHWSV& rhs);

	/// Standard ARM object support
	virtual ARM_Object* Clone() const;

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const;
		
	const std::vector<double>& GetSchedule() const { return itsSchedule;}


	std::vector<double>& BetastT(double t,double T) const
	{ return ARM_ModelParamsHW2F::BetatT( GetModelParam(ARM_ModelParamType::MeanReversion),GetModelParam(ARM_ModelParamType::MeanReversionSpread),t,T); }

	virtual void GenerateSchedule();

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

