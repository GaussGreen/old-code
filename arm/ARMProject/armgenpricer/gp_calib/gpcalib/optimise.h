/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file optimise.h
 *
 *  \brief modelfitter that optimise with no knowledge 
 *		of the derivative function
 *	\author  E.M Ezzine, E. Benhamou
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_OPTIMISE_H
#define _INGPCALIB_OPTIMISE_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "optimisebase.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class WeightedSquareFunc;

class ARM_Optimise : public ARM_OptimiseBase
{
private:
	enum Algorithm
	{
		NAG_OPT_BOUNDS_NO_DERIV,
		NAG_OPT_NLIN_LSQ,
		NAG_OPT_LSQ_DERIV,
		NAG_OPT_LSQ_CHECK_DERIV
	};

	Algorithm itsAlgorithm;
	
	/// specific to ND optimizer
	void OptimizerND(ARM_GP_Vector* Var,ARM_GP_Vector* boundLower,ARM_GP_Vector* boundUpper);
	void OptimizerNDWithDeriv(ARM_GP_Vector* Var,ARM_GP_Vector* boundLower,ARM_GP_Vector* boundUpper);
	void SetTargetFunction();

public:
    ARM_Optimise( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& modelParams,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter = ARM_ModelFitterPtr(NULL ),
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER,
		bool getDetails							= false,
		double tolerance						= 1.0e-3,
		double stepMax							= 1.0e+002,
		bool localSearch						= false,
		size_t factorNb							= 0,
		size_t NbCalibIter						= 1,
		int  algorithm					        = NAG_OPT_BOUNDS_NO_DERIV ,
        ARM_MktTargetType  typeTarget           = ARM_CalibrationTarget::PriceTarget );

    ARM_Optimise(const ARM_Optimise& rhs);
	ARM_Optimise& operator = (const ARM_Optimise& rhs);
    virtual ~ARM_Optimise();

	/// calibration loop
	virtual void Calibrate();

    /// standard ARM Object support
    virtual ARM_Object* Clone() const;
	virtual string ModelFitterName() const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
