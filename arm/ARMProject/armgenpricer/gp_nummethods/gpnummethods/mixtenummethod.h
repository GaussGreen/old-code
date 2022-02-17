/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: mxitenummethod.h,v $
 * Revision 1.1  2003/12/30 16:45:13  ebenhamou
 * Initial revision
 *
 *
 */

/*! \file mixtenummethod.h
 *
 *  \brief 
 *	\author  R. Guillemot A. Schauly
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPINFRA_MIXTENUMMETHOD_H
#define _INGPINFRA_MIXTENUMMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"

CC_BEGIN_NAMESPACE(ARM)

class ARM_PricingModel;

/// \class ARM_MixteNumMethod 
/// \brief Forward/Backward induction numerical method object

class ARM_MixteNumMethod : public ARM_NumMethod
{
public:
	ARM_MixteNumMethod();
	virtual ~ARM_MixteNumMethod();
	ARM_MixteNumMethod( const ARM_MixteNumMethod& rhs);
	ARM_MixteNumMethod& operator=( const ARM_MixteNumMethod& rhs);

    virtual GP_PricingDirection GetPricingDirection() const { return ARM_NumMethod::GP_FWDBCKWDLOOKING; }
	virtual GP_PricingDirection GetPricingDirCurrLoop() const { return itsPricingDirCurrLoop; }
	virtual inline size_t GetLoopNb() const { return 2; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[GetTimeSteps()->size()-1]; }

	/// init before processing
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model);
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {};
	virtual void Finalize( const ARM_PricingModel& model, ARM_PricingStatesPtr& states ) const {};
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const {};
	
	// Do we have to reset the precomputed nodes between each iterations
	virtual void ComputeTimeSteps(ARM_PricingModel& model);
	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const std::vector<double>& eventTimes) const;

	/// forward induct
	virtual ARM_PricingStatesPtr Induct(
        const ARM_PricingModel& model,
		ARM_PricingStatesPtr& states,
		double toTime );

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary* ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;

private:
	GP_PricingDirection itsPricingDirCurrLoop;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

