/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file mcmethod.h
 *	\author  E Benhamou
 *	\version 1.0
 *	\date January 2004
 */


#ifndef _INGPINFRA_MCMETHOD_H
#define _INGPINFRA_MCMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/typedef.h"
#include "gpnumlib/typedef.h"
#include "typedef.h"
#include "impsampler.h"
#include "pathscheme.h"

CC_BEGIN_NAMESPACE( ARM )

/// forward declaration
struct ARM_DiscretisationScheme;
struct ARM_SamplerBase;
struct ARM_SchedulerBase;


//////////////////////////////////////////////
/// \class ARM_MCMethod
/// \brief class for Monte Carlo
/// the design of the generic Monte Carlo is 
///		1) to be model free...
///		2) to support super bucket Monte Carlo. 
///		Default is super bucket Monte Carlo with all paths
///		to do a Monte Carlo in one go...
///		3) fast to run
///
///	1) Model free is achieve by communicating to the model
///		requiring the model to compute drift and volatility
///		part. If a model supports the various function defined
///		in pricingmodel.h, then it supports the generic Monte Carlo
///
///	2) super bucket Monte Carlo is monitored by the vector of buckets
///		defined in itsBuckets with itsCurrentBucket being the current
///		bucket number
///
///	3) fast to run mean that we cache lots of values... and the caching of stateLocalVars and
///		local stdDev and localDrift
///
///	Like any numerical method, it allows to induct from one time to the next time
//////////////////////////////////////////////
class ARM_MCMethod : public ARM_NumMethod
{
private:
	/// random number generator ...
	ARM_RandomGeneratorPtrVector itsRandGenVector;
	// Process
	ARM_MatrixPtrVector itsProcessStates;

    ///... indexed by times, matrix for each state variable (row var) of its local variance between 
	/// current time and previous time for each factor (col var)
    //ARM_SimpleMatrixVector itsModelStateLocalVars;
	double itsFirstInductTime;

	size_t GetNbOfPaths() const;
	bool itsOtherPayoffsFlag;
	ARM_ImpSamplerPtr itsImpSampler;
	ARM_PathSchemePtr itsPathScheme;

protected: 
	ARM_VectorPtr itsBuckets;	/// bucket specification for super bucket monte carlo
	int itsBucketIndex;			/// index of the current bucket (int to allow -1 index!)

public:
	ARM_MCMethod( 
		size_t itersNb,
		const ARM_RandomGeneratorPtrVector& randGenVector,
		ARM_SamplerBase* sampler = NULL,
		size_t MaxBucketSize = 100000,
		ARM_ImpSamplerPtr impSampler = ARM_ImpSamplerPtr(NULL),
		ARM_PathSchemePtr pathScheme = ARM_PathSchemePtr(NULL));
	ARM_MCMethod( 
		size_t itersNb, 
		const ARM_RandomGeneratorPtr& randGen, 
		ARM_SamplerBase* sampler = NULL,
		size_t MaxBucketSize = 100000,
		ARM_ImpSamplerPtr impSampler = ARM_ImpSamplerPtr(NULL),
		ARM_PathSchemePtr pathScheme = ARM_PathSchemePtr(NULL));

	ARM_MCMethod(const ARM_MCMethod& rhs);
	ARM_MCMethod& operator=(const ARM_MCMethod& rhs);
	void Initialize(
		size_t itersNb,
		const ARM_RandomGeneratorPtrVector& randGenVector,
		size_t MaxBucketSize );

	virtual ~ARM_MCMethod();

	/// forward pricing
    virtual GP_PricingDirection GetPricingDirection() const {return ARM_NumMethod::GP_FWDLOOKING;}
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return ARM_NumMethod::GP_FWDLOOKING;}
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[GetTimeSteps()->size()-1]; }
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const;

    /// Initialisation of the numerical method
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
	virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model );
	// Build the process states which will be set in the num method states
	void BuildProcessStates( const ARM_PricingModel& model );

	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {}

    /// Numerical induct for a method
	virtual ARM_PricingStatesPtr Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, double toTime);

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// Does we do not compute otherpayoffs for too large simluations
	virtual bool GetOtherPayoffsFlag() const { return itsOtherPayoffsFlag; }

	/// accessors to the buckets
	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual const ARM_MatrixVector& GetNumMethodStateGlobalVars() const;
	virtual ARM_GP_MatrixPtr GetSpotProbabilities(const std::vector<double>& eventTimes) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const;

	// accessor for the random generator
	ARM_RandomGeneratorPtr GetRandGen(int i=0) const;
	ARM_RandomGeneratorPtrVector GetRandGenVector() const {return itsRandGenVector;}

	const ARM_ImpSamplerPtr& GetImpSampler() const { return itsImpSampler; };
	const ARM_PathSchemePtr& GetPathScheme() const { return itsPathScheme; };

	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;

	// Those functions are used to apply Importance Sampling
	virtual void ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;
	virtual void ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

