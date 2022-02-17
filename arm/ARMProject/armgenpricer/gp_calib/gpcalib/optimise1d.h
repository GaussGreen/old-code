/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 *
 *	\file optimise1d.h
 *
 *  \brief optimisation in one dimension
 *	\author  E.M Ezzine
 *	\version 1.0
 *	\date November 2003
 */


#ifndef _INGPCALIB_OPTIMISE1D_H
#define _INGPCALIB_OPTIMISE1D_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for std::vector<double> and ARM_Matrix
#include "gpbase/env.h"
#include "gpbase/port.h"
#include "optimisebase.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class WeightedSquareFunc;

class ARM_Optimise1D : public ARM_OptimiseBase
{
    /// function that calls the NAG optimizer
    void Optimizer1D(double Var,double boundLower,double boundUpper);
	double itsRelativeTolerance;
	double itsAbsoluteTolerance;
    
	/// function for the constructor
	void Validate();

public:
    ARM_Optimise1D( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& modelParam ,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter = ARM_ModelFitterPtr(NULL ),
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER ,
		bool getDetails							= false,
		size_t factorNb							= 0,
		double relativeTolerance				= 0.01,
		double absoluteTolerance				= 0.0001,
		size_t NbCalibIter						= 1,
        ARM_MktTargetType       = ARM_CalibrationTarget::PriceTarget
        );

    ARM_Optimise1D(const ARM_Optimise1D& rhs);
	ARM_Optimise1D& operator = (const ARM_Optimise1D& rhs);
    virtual ~ARM_Optimise1D();

	/// calibration
    virtual void Calibrate();

	/// standard ARM Object support
    virtual ARM_Object* Clone() const;
	virtual string ModelFitterName() const;
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
