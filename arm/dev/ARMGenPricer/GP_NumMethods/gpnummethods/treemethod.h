/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: treemethod.h,v $
 * Revision 1.1  2003/10/13 07:52:03  jmprie
 * Initial revision
 *
 *
 */



/*! \file treemethod.h
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


#ifndef _INGPINFRA_TREEMETHOD_H
#define _INGPINFRA_TREEMETHOD_H

/// this header has to come first
/// as env.h defines pre-processor constants for checking and validation
/// namely for ARM_GP_Vector and ARM_Matrix
#include "gpbase/env.h"

/// gpbase
#include "gpbase/port.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/mempoolmatrix.h"

/// gpinfra
#include "gpinfra/typedef.h"
#include "gpinfra/nummethod.h"

#include "gridindex.h"

#include <vector>
CC_USING_NS(std,vector)

CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////
/// \class ARM_TreeMethod
/// \brief
/// ARM_TreeMethod class computes conditional
/// expectations using a recombining trinomial
/// tree with a constant space sampling in
/// each slice and each independant dimension
//////////////////////////////////////////////
class ARM_TreeMethod : public ARM_NumMethod
{
private:

    /// Space sampling value in each dimension
    /// for each time step
    /// Used with the index list it defines the
    /// discrete mapping to state variables
    ARM_TriangularMatrixVector itsSpaceSteps;

    /// For each time step, maximum value of the tree index
    ARM_GridIndexVector itsMaxIndex; 

   
   
    ///itsNumMethodStateLocalStdDevs (tp-1 -> tp);
    //..itsNumMethodStateGlobalStdDevs mapping (0->tp)
    /// (for truncation purpose)
    ARM_TriangularMatrixVector itsNumMethodStateGlobalStdDevs;
	ARM_MatrixVector itsNumMethodStateGlobalVars;

	ARM_VectorVector itsProbaSpotForward;

	ARM_TreeTransitionsVector itsTransitionStates;

    /// Constant truncation factor in standard deviation
    /// of state variables of the pricing model
    double itsStdDevRatio;

    /// Minimum value of annualized local standard deviation
    /// Below this value, space step resizing occurs
    double itsMinStdDev;

    /// Minimum number of time steps required before first
    /// date in the schedule given by the model at initialisation stage
    int itsNbMinSteps;


    /// Flag to force probability computation
    bool itsIsProbabilityComputation;

    void CopyNoCleanUp(const ARM_TreeMethod& rhs);
    void CleanUp();

    // Schedule uniform in time
    void ComputeCstTimeSchedule();

    // Schedule tends to be uniform in variance
    void ComputeCstVarSchedule(const ARM_PricingModel& model);

    // Schedule mixing time and local variance
    void ComputeMixTimeVarSchedule(const ARM_PricingModel& model);
    double YtoX(ARM_GP_Vector& X,ARM_GP_Vector& Y,double y,int& lastIdx,int nextIdx);

    // For each time step, compute the continuous mapping i.e.
    // cholesky transformation of local and global
    // variances/covariances matrix
    void InitContinuousMapping( const ARM_PricingModel& model );

    /// For each time step, compute space steps i.e. discrete mapping
    void InitDiscreteMapping();

    /// For each time step, get the maximum value of the tree index
    void ComputeMaxIndex();

    /// Degenerate a transition state to a binomial one
    void ComputeBinomialTransition( int dirIdx,
                                    double driftErr,
                                    ARM_GP_Vector& proba,
                                    ARM_GP_Vector& nextStateVariance) const;

    /// Build all transition states from the current slice
    /// indexed by timeIdx to the next one
    ARM_TreeTransitions* ComputeTransitionStates(int timeIdx) const;

    /// Compute the probability of a transition
    double ComputeStateProba(
        ARM_GridIndex& relIndex,
        const ARM_VectorVector& elemProbas) const;

    double ComputeStateProba(
        ARM_GridIndex& relIndex,
        const ARM_VectorVector& elemProbas, 
        ARM_GridIndex& lastIndex,
        ARM_GP_Vector& lastProbas) const;

    /// Compute the state variables from a position in the grid
    void ComputeStateVariables(
        int timeIdx,
        const ARM_GridIndex& maxIndex,
        ARM_PricingStatesPtr& states) const;

	void ComputeAndSetTimeSteps( const ARM_PricingModel& model) const;

public:
	ARM_TreeMethod(double stdDevRatio=5.0,double minStdDev=0.001,int nbMinStep=0,bool isProbabilityComputation=false);
	ARM_TreeMethod(const ARM_TreeMethod& rhs);
	ARM_TreeMethod& operator = (const ARM_TreeMethod& rhs);
	virtual ~ARM_TreeMethod();

    /// Accessors
    ARM_GridIndex& GetMaxIndex(int i) const {return *(itsMaxIndex[i]);}
    int GetNbDir() const {(itsMaxIndex.size()>0 && itsMaxIndex[0] ?  itsMaxIndex[0]->size() : 0);}
    const ARM_GP_TriangularMatrix& GetSpaceStep(int i) {return *(itsSpaceSteps[i]);}
  	const ARM_GP_Vector& GetProbaSpotForward(int i) {return *(itsProbaSpotForward[i]);}
    const ARM_GP_TriangularMatrix& GetNumMethodStateGlobalStdDevs(int i) {return *(itsNumMethodStateGlobalStdDevs[i]);}
	virtual const ARM_MatrixVector& GetNumMethodStateGlobalVars() const {return itsNumMethodStateGlobalVars;}
	const ARM_TreeTransitions& GetTransitionStates(int i) {return *(itsTransitionStates[i]);}
    void SetProbilityComputation(bool isProbabilityComputation) {itsIsProbabilityComputation=isProbabilityComputation;}

    // Backward pricing
    virtual GP_PricingDirection GetPricingDirection() const {return ARM_NumMethod::GP_BCKWDLOOKING;}
	virtual GP_PricingDirection GetPricingDirCurrLoop() const {return ARM_NumMethod::GP_BCKWDLOOKING;}

	// Loop management
	virtual inline size_t GetLoopNb() const { return 1; }
	virtual ARM_PricingStatesPtr ReInitLoop(const ARM_PricingModel& model);

    virtual double GetPricingFinalTimeStep() const { return (*(GetTimeSteps()))[0]; }

    /// Initialisation of the tree
    virtual ARM_PricingStatesPtr Init( ARM_PricingModel& model, double firstInductTime);
    virtual ARM_PricingStatesPtr ReInit( const ARM_PricingModel& model);
	virtual void ModifyParseTree( ARM_PricingAdviser * ) const {};
	virtual void Finalize( const ARM_PricingModel& model, ARM_PricingStatesPtr& states ) const {};

	virtual ARM_VectorPtr GetBuckets() const; 
	virtual size_t GetBucketIndex() const;
	virtual ARM_PricerInfo* CreatePricerInfo( const ARM_PricingModel& model ) const;

    virtual ARM_GP_MatrixPtr GetSpotProbabilities(const ARM_GP_Vector& eventTimes) const;
	virtual ARM_GP_VectorPtr GetArrowDebreuPrices(size_t timeIdx, const ARM_PricingModel& model ) const
    {return ARM_GP_VectorPtr( timeIdx > 0 ? new ARM_GP_Vector(*(itsProbaSpotForward[timeIdx-1])) : NULL );}

    /// Numerical expectation of discounted payoffs
    /// from lastTimeIdx to endTimeIdx
    /// (no event date allowed between)
	virtual ARM_PricingStatesPtr Induct(
        const ARM_PricingModel& model,
		ARM_PricingStatesPtr& states,
		double toTime);

	// Compute the exercise Boundary of the exercise node
	virtual ARM_ExerciseBoundary* ComputeExercise(const ARM_GP_VectorPtr& payoff, const ARM_GP_VectorPtr& contOpt, const ARM_GP_MatrixPtr& StatesVector );

	void InductOneStep(
        ARM_GP_Matrix& payoffs,
		ARM_TreeTransitions* transitionStates,
		int timeIndex);

	void InductForwardOneStep(
        ARM_GP_Vector& payoffs,
		ARM_TreeTransitions* transitionStates,
		int timeIndex);

	void InductForwardFromToNextStep(
        ARM_GP_Vector& payoffs,
		int timeIdx_1,
		int timeIdx_2);

	//Construct the tree (at each time steps  itsTransitionStates and the forward probabilities itsProbaSpotForward) 
	void TreeConstruction(const ARM_GP_Vector& timeSteps, vector< ARM_GP_Vector* >& probaSoptFwd, 
										vector< ARM_TreeTransitions* >& transitionStates, bool needTreeProbaForward);
	
	/// Standard ARM support
	virtual ARM_Object* Clone() const;
    virtual string toString(const string& indent="", const string& nextIndent="") const;
};

CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

