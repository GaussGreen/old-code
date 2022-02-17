/*
 *
 * Copyright (c) IXIS CIB November 2004 Paris
 *
 *	\file amcmethod.h
 *	\author  A Schauly
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPINFRA_AMCMETHOD_H
#define _INGPINFRA_AMCMETHOD_H

#include "typedef.h"
#include "mcmethod.h"
#include "amc_exercboundcalc.h"
#include "gpnumlib/random.h"
#include "gpbase/assignop.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_AMCMethod : public ARM_MCMethod
{
private:
	ARM_ExerciseBoundaryCalcPtr itsExerBoundCalc;
	GP_PricingDirection itsPricingDirCurrLoop;
	ARM_MatrixPtrVector modelStatesVector;
	std::vector<ARM_PricingStatesContextPtrVectorPtr> pricingStatesContextVector;

	ARM_MatrixPtrVector::reverse_iterator iter;
	std::vector<ARM_PricingStatesContextPtrVectorPtr>::reverse_iterator iter2;
	double itsPrevLastInductTime, itsLastInductTime;

	size_t itsLoopNb;

public: 
	/// Constructors/Destructor
	ARM_AMCMethod( 
		size_t itersNb,
		const ARM_RandomGeneratorPtrVector& randGenVector,
		ARM_SamplerBase* sampler,
		ARM_ExerciseBoundaryCalc* exerBoundCal, 
		size_t MaxBucketSize = 100000,
		const ARM_ImpSamplerPtr& impSampler = ARM_ImpSamplerPtr(NULL),
		const ARM_PathSchemePtr& pathScheme = ARM_PathSchemePtr(NULL));

	ARM_AMCMethod( 
		size_t itersNb,
		const ARM_RandomGeneratorPtr& randGen, 
		ARM_SamplerBase* sampler,
		ARM_ExerciseBoundaryCalc* exerBoundCal, 
		size_t MaxBucketSize = 100000,
		const ARM_ImpSamplerPtr& impSampler = ARM_ImpSamplerPtr(NULL),
		const ARM_PathSchemePtr& pathScheme = ARM_PathSchemePtr(NULL));

	virtual ~ARM_AMCMethod();

	//ASSIGN_OPERATOR(ARM_AMCMethod);

	ARM_AMCMethod(const ARM_AMCMethod& rhs);

	// Pricing Direction
	virtual GP_PricingDirection GetPricingDirection() const { return ARM_NumMethod::GP_FWDBCKWDLOOKING; }
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return itsPricingDirCurrLoop; }
	virtual inline size_t GetLoopNb() const { return itsLoopNb; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);

	// Init before Processing ; Reinit ; NodeReset
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model);
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const;
	virtual bool NeedToCreateDefaultArgument() const;

	// Free memory at end of pricing
	virtual void Finalize();

	/// Standard ARM support
	virtual ARM_GP_MatrixPtr GetSpotProbabilities(const std::vector<double>& u) const
	{
		return ARM_GP_MatrixPtr(NULL);;
	};
	virtual ARM_Object* Clone() const { return new ARM_AMCMethod(*this); };
    virtual string toString(const string& indent="", const string& nextIndent="") const; 

	// Induct
	virtual ARM_PricingStatesPtr Induct(const ARM_PricingModel& model,ARM_PricingStatesPtr& states,double toTime );
	virtual ARM_PricingStatesPtr InductBackward(const ARM_PricingModel& model,ARM_PricingStatesPtr& states,double toTime );

	// Compute the exercise Boundary of the exercise node
	ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, 
		const ARM_GP_MatrixPtr& StatesVector );

};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

