/*!
 *
 * Copyright (c) CDC IXIS CM November 2004 Paris
 *
 *	\file optimisewithbrent.h
 *
 *  \Brent optimizer 
 *	\author  A. Triki
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPCALIB_OPTIMISEWITHBRENT_H
#define _INGPCALIB_OPTIMISEWITHBRENT_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "optimisebase.h"
#include "targetfunc.h"


CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
class WeightedSquareFunc;

class ARM_OptimiseWithBrent : public ARM_OptimiseBase
{
	double itsPrecision;
	/// specific to 1D optimizer
    void Optimizer1DWithBrent(double Var,double boundLower,double boundUpper,double Max_Iter,double Precision);

public:
    ARM_OptimiseWithBrent( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& modelParam ,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter = ARM_ModelFitterPtr(NULL ),
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER,
		bool getDetails							= false,
		double tolerance						= 1.0e-3,
		double stepMax							= 1.0e+002,
		size_t factorNb							= 0,
		size_t NbCalibIter						= 1,
        ARM_MktTargetType                       = ARM_CalibrationTarget::PriceTarget);

    ARM_OptimiseWithBrent(const ARM_OptimiseWithBrent& rhs);
	ARM_OptimiseWithBrent& operator = (const ARM_OptimiseWithBrent& rhs);
    virtual ~ARM_OptimiseWithBrent();//OK

    virtual void Calibrate();

    /// standard ARM Object support
    virtual ARM_Object* Clone() const;
	virtual string ModelFitterName() const;
};



////////////////////////////////////////////////////
///	Class  : WeightedSquareCalculate
///	Routine: Call Brent Optimizer
///	Returns: 
///	Action : optimisation routine
////////////////////////////////////////////////////
class WeightedSquareCalculate
{
public: 
	MultiDimFunc* itsFunction;

public: 
	WeightedSquareCalculate(const MultiDimFunc* Func=NULL);
	WeightedSquareCalculate(const WeightedSquareCalculate& rhs );
	WeightedSquareCalculate& operator=( const WeightedSquareCalculate& rhs );
	~WeightedSquareCalculate();
	double operator () ( double x ) const;
};


CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
