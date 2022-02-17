/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *


/*! \file ModelParamsHWSV1F.h
 *
 *  \brief 
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date November 2005
 */


#ifndef _INGPMODELS_MODELPARAMSHWSV1F_H
#define _INGPMODELS_MODELPARAMSHWSV1F_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"
#include "modelparamsmsv.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"			/// for calibParamVector
#include "gpinfra/modelparamtype.h"		/// for ARM_ModelParamType
#include "gpmodels/modelparamshw1f.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMSV1F
// \brief Interface class for model parameters of the MSV 1F model
//-----------------------------------------------------------------------------
class ARM_ModelParamsHWSV1F : public ARM_ModelParams
{
protected :
	mutable ARM_GP_Vector itsSchedule;

public:
	ARM_ModelParamsHWSV1F( const ARM_ModelParamsHWSV1F& rhs );
	ARM_ModelParamsHWSV1F( const ARM_ModelParamVector& params=ARM_ModelParamVector());
	virtual ~ARM_ModelParamsHWSV1F();
    ARM_ModelParamsHWSV1F& operator = (const ARM_ModelParamsHWSV1F& rhs);

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 2; }
		
	const ARM_GP_Vector& GetSchedule() const { return itsSchedule;}


	/// Standard ARM object support
	virtual ARM_Object* Clone() const;

    /// Coefficient of the state variable (for Zc closed form formula)
    virtual double BetatT(double t,double T)
	{ return ARM_ModelParamsHW1F::BetatT( GetModelParam(ARM_ModelParamType::MeanReversion),t,T); }

	void GenerateSchedule();

	/// calibration part
    virtual void PreProcessing(ARM_ModelFitter& modelFitter,int factorNb=0) {};
    virtual void PostProcessing(const ARM_ModelFitter& modelFitter,ARM_PricingModel* model,int factorNb=0 ) {};

};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

