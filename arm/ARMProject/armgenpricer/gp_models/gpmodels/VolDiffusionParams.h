/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *


/*! \file ModelParamsMSV1F.h
 *
 *  \brief 
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */


#ifndef _INGPMODELS_VOLDIFFUSIONPARAMS_H
#define _INGPMODELS_VOLDIFFUSIONPARAMS_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "modelparamsmsv.h"
#include "gpbase/port.h"
#include "gpcalib/typedef.h"			/// for calibParamVector
#include "gpinfra/modelparamtype.h"		/// for ARM_ModelParamType


/// forward declaration
class ARM_IRIndex;

CC_BEGIN_NAMESPACE( ARM )

class ARM_CurveModelParam;

//-----------------------------------------------------------------------------
// \class ARM_ModelParamsMSV1F
// \brief Interface class for model parameters of the MSV 1F model
//-----------------------------------------------------------------------------
class ARM_VolParams : public ARM_ModelParams
{
public:
	ARM_VolParams( const ARM_VolParams& rhs );
	ARM_VolParams( const ARM_ModelParamVector& params = ARM_ModelParamVector() );
	virtual ~ARM_VolParams();
    ARM_VolParams& operator = (const ARM_VolParams& rhs);

	/// Dates at which modelparams values change
	virtual ARM_GP_VectorPtr ModelParamsTimeSteps() const;

	/// How many factors?
	virtual size_t FactorCount() const { return 1; }

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

