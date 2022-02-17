/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file finummethod.h
 *
 *  \brief 
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2003
 */


#ifndef _INGPINFRA_FINUMMETHOD_H
#define _INGPINFRA_FINUMMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"

CC_BEGIN_NAMESPACE(ARM)

class ARM_PricingModel;

/// \class ARM_FINumMethod 
/// \brief Backward induction numerical method object

class ARM_FINumMethod : public ARM_NumMethod
{
public:
	ARM_FINumMethod();
	virtual ~ARM_FINumMethod();
	ARM_FINumMethod( const ARM_FINumMethod& rhs);
	ARM_FINumMethod& operator=( const ARM_FINumMethod& rhs);

    virtual GP_PricingDirection GetPricingDirection() const { return ARM_NumMethod::GP_FWDLOOKING; }
	virtual GP_PricingDirection GetPricingDirCurrLoop() const { return ARM_NumMethod::GP_FWDLOOKING; }
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[GetTimeSteps()->size()-1]; }
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const {};

	/// init before processing
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model);
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {};
	virtual void Finalize( const ARM_PricingModel& model, ARM_PricingStatesPtr& states ) const {};
	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const;
	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary* ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// forward induct
	virtual ARM_PricingStatesPtr Induct(
        const ARM_PricingModel& model,
		ARM_PricingStatesPtr& states,
		double toTime );

	/// standard ARM Object support
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="", const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

