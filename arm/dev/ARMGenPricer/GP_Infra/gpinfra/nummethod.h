/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file nummethod.h
 *
 *  \brief 
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_NUMMETHOD_H
#define _INGPINFRA_NUMMETHOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix

#include "gpbase/env.h"
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/cloneutilityfunc.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/gpmatrixtriangular.h"
#include "typedef.h"

#include "samplerbase.h"
#include "schedulerbase.h"

#include <string>
CC_USING_NS( std, string )

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

struct ARM_PricerInfo;
class ARM_PricingAdviser;

/// To allow using of new tree method
extern const bool isTree2G;

//////////////////////////////////////////////
/// \class ARM_NumMethod 
/// \brief numerical method object
//////////////////////////////////////////////

class ARM_NumMethod: public ARM_RootObject
{
public:
	/// enum for direction
	/// might be useful to have rather derived class for avoiding 
	///	ugly switch!
    enum GP_PricingDirection
    {
        GP_AMBIGUOUS = 0,
        GP_FWDLOOKING,
        GP_BCKWDLOOKING,
        GP_FWDBCKWDLOOKING,
        GP_UNKNOWN,
    };

	/// text corresponding to the direction enum
	static const string GP_PricingDirectionTxt[];

private:
    /// Number of steps for schedule building
    int itsNbSteps;

    /// Schedule (tp) with t0=0
    int itsLastTimeIdx;
	ARM_GP_Vector* itsTimeSteps;

	ARM_SamplerBase* itsSampler;

    void CopyNoCleanUp(const ARM_NumMethod& rhs);
    void CleanUp();

public:
	ARM_NumMethod(const ARM_SamplerBase* sampler = NULL);
	ARM_NumMethod( const ARM_NumMethod& rhs);
	virtual ~ARM_NumMethod();
	ARM_NumMethod& operator=( const ARM_NumMethod& rhs);

    /// ----------------- Accessors
    inline int GetNbSteps() const {return itsNbSteps;}
    inline void SetNbSteps(int nbSteps) {itsNbSteps=nbSteps;}
    inline int GetLastTimeIdx() const {return itsLastTimeIdx;};
    inline void SetLastTimeIdx(int lastTimeIdx) {itsLastTimeIdx=lastTimeIdx;}

    inline const ARM_GP_Vector* const GetTimeSteps() const {return itsTimeSteps;}
    inline ARM_GP_Vector* GetTimeSteps() { return itsTimeSteps;}
    inline double GetTimeStep(int timeIdx) const { return (*itsTimeSteps)[timeIdx]; }
    inline double GetLastTimeStep() const { return (*itsTimeSteps)[itsLastTimeIdx];}
    void SetTimeSteps( const ARM_GP_Vector& timeSteps );

	const ARM_SamplerBase* const GetSampler() const { return itsSampler; }
    ARM_SamplerBase* GetSampler() { return itsSampler; }
    void SetSampler(const ARM_SamplerBase* const sampler);

	/// ----------- Init methods
	/// To be redefined by every method to
    /// set its internal datas
	/// model is not const ARM_PricingModel& but rather ARM_PricingModel& 
	/// to allow to call back PostInit to do some caching!
    virtual ARM_PricingStatesPtr Init(ARM_PricingModel& model, double firstInductTime) = 0;
	virtual ARM_PricingStatesPtr ReInit(const ARM_PricingModel& model) = 0;
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const = 0;
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const = 0;
	virtual void Finalize() {};

    virtual bool SupportIntegratedSampling() const { return true; }
	virtual bool SupportFrontierComputation() const { return false; }

	/// ---------------------- bucket part of the numerical method!
	/// accessors to the buckets
	/// returns the vector of buckets
	virtual ARM_VectorPtr GetBuckets() const = 0; 
	/// returns the current bucket index
	virtual size_t GetBucketIndex() const = 0;
	/// returns the corresponding pricer info
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const = 0;
	
    /// ---------------- Induction 
	/// induction style implemented
    virtual GP_PricingDirection GetPricingDirection() const = 0;
	virtual GP_PricingDirection GetPricingDirCurrLoop() const = 0;
	/// Either backward or forward induct
	virtual ARM_PricingStatesPtr Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states, double toTime ) = 0;

	/// ---------------- Loop management
	virtual size_t GetLoopNb() const = 0;
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model) = 0;
    virtual double GetPricingFinalTimeStep() const = 0;
	virtual bool NeedToCreateDefaultArgument() const { return false; }

	/// ---------------- Information management
	// Global Variance
	virtual const ARM_MatrixVector& GetNumMethodStateGlobalVars() const;
    /// Probabilities to access futur states
    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const = 0;
	virtual ARM_GP_VectorPtr GetSpotProbabilities(size_t timeIdx) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const;

	/// Does model advise to use other payoffs ?
	virtual bool GetOtherPayoffsFlag() const { return true; }

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary * ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector ) = 0;

	// Used for the closed form num method
	virtual double ConvertEvalDate(double evalDate, double asOfDate) const;

	/// Standard ARM Object support 
    virtual string toString(const string& indent="", const string& nextIndent="") const;
	virtual ARM_CLASS_NAME GetRootName() { return ARM_NUMMETHOD; }

	// Those functions are used to apply Importance Sampling
	virtual void ProcessPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;
	virtual void ProcessUnPaidPayoffs(ARM_VectorPtr& payoffs, double evalTime) const;
};



CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

