/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file binummethod.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */



#ifndef _INGPINFRA_BINUMMETHOD_H
#define _INGPINFRA_BINUMMETHOD_H

#include "gpbase/port.h"
#include "gpinfra/nummethod.h"


CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;

/// \class ARM_BINumMethod 
/// \brief Backward induction numerical method object

class ARM_BINumMethod : public ARM_NumMethod
{
public:
	ARM_BINumMethod();
	virtual ~ARM_BINumMethod();

	ARM_BINumMethod( const ARM_BINumMethod& rhs);
	ARM_BINumMethod& operator=( const ARM_BINumMethod& rhs);


    virtual GP_PricingDirection GetPricingDirection() const { return ARM_NumMethod::GP_BCKWDLOOKING; }
	virtual GP_PricingDirection GetPricingDirCurrLoop() const { return ARM_NumMethod::GP_BCKWDLOOKING; };
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[0]; }

	/// init before processing
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model);
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {};
	virtual void Finalize( const ARM_PricingModel& model, ARM_PricingStatesPtr& states ) const {};
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const {};

	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const std::vector<double>& eventTimes) const;

	/// backward induct
	virtual ARM_PricingStatesPtr Induct(
        const ARM_PricingModel& model,
		ARM_PricingStatesPtr& states,
		double toTime );

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	/// standard ARM Object support
	virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_Object* Clone() const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

