/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file treebase.h
 *
 *  \brief
 *
 *	\author  JM Prie, E. Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#ifndef _INGPNUMMETHODS_TREEBASE_H
#define _INGPNUMMETHODS_TREEBASE_H

/// gpbase
#include "gpbase/port.h"

/// gpinfra
#include "gpinfra/nummethod.h"

#include "typedef.h"

CC_BEGIN_NAMESPACE( ARM )

class ARM_PricingModel;
struct ARM_SliceBase;
struct ARM_SamplerBase;
struct ARM_TruncatorBase;
struct ARM_ReconnectorBase;
struct ARM_SmootherBase;


class ARM_TreeBase : public ARM_NumMethod
{
	ARM_SliceVectorPtr itsSlices;
    ARM_TruncatorBase* itsTruncator;
    ARM_ReconnectorBase* itsReconnector;
    ARM_SmootherBase* itsSmoother;

    bool itsComputeSpotProbas;

    ARM_PricingStatesPtr itsStates;

    void CopyNoCleanUp(const ARM_TreeBase& rhs);
    void CleanUp();

protected:
    void ExerciseSmoothing1D(const ARM_GP_Vector& exerFct, double coef, ARM_GP_Vector& smoothValues, ARM_GP_Vector& exerStates) const;

    /// Default smoothing implentation is a 1D one
    virtual void ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const { ExerciseSmoothing1D(*exerFct,1.0,*smoothValues,*exerStates); } 

public:

    /// To indicate the arbitrary number of days to add to the last slice
    /// to calibrate short rate on it fulfilling
    /// the DF(0,sliceTime+LastSliceCalibrationTermInDays)
    static const double LastSliceCalibrationTermInDays;

	ARM_TreeBase( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas=false);
	ARM_TreeBase(const ARM_TreeBase& rhs);
	ARM_TreeBase& operator = (const ARM_TreeBase& rhs);
	virtual ~ARM_TreeBase();

    /// Accessors
    ARM_SliceVectorPtr GetSlices() { return itsSlices; }
    const ARM_SliceVectorPtr GetSlices() const { return itsSlices; }
    void SetSlices( const ARM_SliceVectorPtr sliceVector );

    const ARM_ReconnectorBase* const GetReconnector() const { return itsReconnector; }
    ARM_ReconnectorBase* GetReconnector() { return itsReconnector; }

    const ARM_TruncatorBase* const GetTruncator() const { return itsTruncator; }
    ARM_TruncatorBase* GetTruncator() { return itsTruncator; }
    
    const ARM_SmootherBase* const GetSmoother() const { return itsSmoother; }
    ARM_SmootherBase* GetSmoother() { return itsSmoother; }
    
    void SetComputeSpotProbas(bool probasFlag) { itsComputeSpotProbas=probasFlag; }
    bool GetComputeSpotProbas() const { return itsComputeSpotProbas; }

    /// Compute a schedule
    ARM_GP_VectorPtr ComputeTimeSteps(const ARM_PricingModel& model) const;
	virtual void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const {};

    /// Maximum number of nodes connected by a transition
    virtual size_t GetTransitionSize() const = 0;

    /// Generation of numerical states at a time step
    virtual void ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const = 0;

    virtual size_t dim() const = 0;
	virtual bool SupportFrontierComputation() const { return dim()==1; }

    virtual const ARM_MatrixVector& GetNumMethodStateGlobalVars() const;

    /// Initialisation of the tree
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime)
    { return Init(model,firstInductTime,ARM_VectorPtrVector(0),0,ARM_GP_Vector(0),ARM_GP_VectorPtr(NULL)); }

    ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime, const ARM_VectorPtrVector& prevDriftCorrections, size_t driftOffset, const ARM_GP_Vector& targetLocalPayoffs, ARM_GP_VectorPtr& timeSteps);

    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model); // not supported
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {};
	virtual void Finalize( const ARM_PricingModel& model, ARM_PricingStatesPtr& states ) const {};
	
	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const;

    // Backward pricing
	virtual ARM_PricingStatesPtr Induct( const ARM_PricingModel& model, ARM_PricingStatesPtr& states,  double toTime);
    virtual GP_PricingDirection GetPricingDirection() const {return ARM_NumMethod::GP_BCKWDLOOKING;}
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return ARM_NumMethod::GP_BCKWDLOOKING;}
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);
    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[0]; }

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const;
	virtual ARM_GP_VectorPtr GetSpotProbabilities(size_t timeIdx) const;

	/// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary* ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

    /// Collect drift correction
    void GetDriftCorrections(ARM_VectorPtrVector& driftCorrections ) const;

	/// Standard ARM support
    virtual string toString(const string& indent="",const string& nextIndent="") const;
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

