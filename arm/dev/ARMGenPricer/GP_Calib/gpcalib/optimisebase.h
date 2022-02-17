/*!
 *
 * Copyright (c) CDC IXIS CM November 2003 Paris
 * 
 * \file optimisebase.h
 *  \brief nagfunction provides some simple function for
 *		nag optimisation
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2004
 */


#ifndef _INGPCALIB_OPTIMISEBASE_H
#define _INGPCALIB_OPTIMISEBASE_H

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "modelfitter.h"

/// nag headers
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nage04.h"


CC_BEGIN_NAMESPACE( ARM )

class MultiDimFunc;

class ARM_OptimiseBase : public ARM_ModelFitter
{
private:
	double itsTolerance;
	double itsStepMax;
	bool itsLocalSearch;
	ARM_GP_Vector* itsInitialValuesVec;
	ARM_GP_Matrix* itsInitialValuesMat;
	MultiDimFunc* itsFunction;

	/// method for the constructor
	void Validate();
    void Initialise();
	void SetCurrentAsInitialValue();

public:
	static const bool DefaultGetDetails;
	static const double DefaultTolerance;
	static const double DefaultStepMax;
	static const bool LocalSearch;
	static const double DefaultRelativeTolerance;
	static const double DefaultAbsoluteTolerance;

    ARM_OptimiseBase( ARM_PricingModel* model,
        const ARM_StdPortfolioPtr portfolio, 
        const ARM_ModelParamVector& modelParams ,
        ARM_ModelFitterPtr& linkedModelFitter	= ARM_ModelFitterPtr(NULL ),
        ARM_ModelFitterPtr& previousModelFitter = ARM_ModelFitterPtr(NULL ),
        const size_t max_iter					= ARM_NumericConstants::ARM_GP_MAX_ITER ,
		bool getDetails							= ARM_OptimiseBase::DefaultGetDetails,
		double tolerance						= ARM_OptimiseBase::DefaultTolerance ,
		double stepMax							= ARM_OptimiseBase::DefaultStepMax,
		bool localSearch						= ARM_OptimiseBase::LocalSearch,
		size_t factorNb							= 0,
		size_t NbCalibIter						= 1,
        ARM_MktTargetType                       = ARM_CalibrationTarget::PriceTarget );

    ARM_OptimiseBase(const ARM_OptimiseBase& rhs);
	ARM_OptimiseBase& operator = (const ARM_OptimiseBase& rhs);
    virtual ~ARM_OptimiseBase();

	/// the core of the calibration
    virtual void Calibrate() = 0;

	/// ----------- general optimizer part
    /// Merge all parameters to input its in Nag algorithm
    ARM_VectorVector  MergeVariables();
	/// Split all parameters afer optimisation
    void SplitVariables(const ARM_GP_Vector& variables);
	/// Set the fileName;
	void StoreDetailsInFileIfRequired( Nag_E04_Opt& options, NagError& fail );
	
	/// accessors
	inline double GetTolerance() const {return itsTolerance; }
	inline double GetStepMax() const {return itsStepMax; }
	inline bool GetLocalSearch() const { return itsLocalSearch; }
	inline const MultiDimFunc* const GetFunction() const { return itsFunction; }
	void SetFunctionNoClone( MultiDimFunc* function );

    /// standard ARM Object support
	virtual string ModelFitterName() const = 0;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	void ProcessNagWarning( NagError fail );
};

CC_END_NAMESPACE()

#endif
//-----------------------------------------------------------------------------
/*---- End of file ----*/
