/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree2d.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/tree2d.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/slice.h"

#define FIRST_STATE_VARIABLE    0
#define SECOND_STATE_VARIABLE   1

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree2D::ARM_Tree2D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas )
:	ARM_TreeBase( sampler,truncator,reconnector,smoother,computeSpotProbas )
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree2D::ARM_Tree2D(const ARM_Tree2D& rhs)
:   ARM_TreeBase(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: Assignement operator
///	Returns: ARM_Tree2D&
///	Action : 
////////////////////////////////////////////////////
ARM_Tree2D& ARM_Tree2D::operator= (const ARM_Tree2D& rhs )
{
    if( this != &rhs )
        ARM_TreeBase::operator =(rhs);
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: ComputeStateVariables
///	Returns: List of state variables
///	Action : Compute the state variables for a
///          given position in the tree (defined
///          by its timeIdx)
////////////////////////////////////////////////////
void ARM_Tree2D::ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const
{
    ARM_SliceNDBase* slice = (*GetSlices())[timeIdx]->ToSliceNDBase();
    size_t nbStates = slice->size();

    /// Restore X states from the slice
    ARM_GP_MatrixPtr xStates = slice->GetXStates();
    if(xStates == ARM_GP_MatrixPtr(NULL))
    {
        /// Not previously computed then compute Z states...
        /// Compute Z states
        ARM_GP_MatrixPtr zStates = slice->GetZStates();

        /// ... and restore X states through the sampler
        const ARM_SamplerNDBase* sampler = GetSampler()->ToSamplerNDBase();
        xStates = sampler->ComputeZtoXStates(timeIdx,zStates);
    }

    /// Set numerical method states = X states
    if(states == ARM_PricingStatesPtr(NULL))
        states = ARM_PricingStatesPtr( new ARM_PricingStates(nbStates,2,0,2) );
    else if(states->size() != nbStates || states->NumMethodStatesSize() != 2)
        states->resizeNumMethodStates(2,nbStates);
	
	/// Add the drift correction to the first variable
    ARM_GP_VectorPtr driftCorrection = slice->GetDriftCorrectionVect();
    for(size_t stateIdx=0; stateIdx<nbStates; ++stateIdx )
    {
        states->SetNumMethodState(stateIdx,FIRST_STATE_VARIABLE,(*xStates)(FIRST_STATE_VARIABLE,stateIdx)+(*driftCorrection)[FIRST_STATE_VARIABLE]);
        states->SetNumMethodState(stateIdx,SECOND_STATE_VARIABLE,(*xStates)(SECOND_STATE_VARIABLE,stateIdx)+(*driftCorrection)[SECOND_STATE_VARIABLE]);
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: ExerciseSmoothing
///	Returns: void
///	Action : Compute a smoothed exercise function
////////////////////////////////////////////////////
void ARM_Tree2D::ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const
{
    int timeIdx = GetLastTimeIdx();
    const ARM_SliceNDBase* slice = (*(GetSlices()))[timeIdx]->ToSliceND();
    const ARM_IntVector& minIdx = slice->GetMin();
    const ARM_IntVector& maxIdx = slice->GetMax();
    ARM_TreeIndex index0(minIdx,maxIdx);
    ARM_TreeIndex index1;

    size_t i0,range0=maxIdx[0]-minIdx[0]+1;
    size_t i1,range1=maxIdx[1]-minIdx[1]+1;
    size_t maxRange = (range0 > range1 ? range0 : range1);

    /// 2D slice is smoothed by a 1D smoothing for each column then for each row
    /// Global correction is divided by 2
    double coef = 0.5;

#ifdef __GP_STRICT_VALIDATION
    if(exerFct->size() != range0*range1)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : slice size mismatch in Tree2D" );
#endif

    std::vector<double> exerFct1D,smooth1D,unusedExer1D;
    ARM_IntVector indexPos(maxRange);

    exerFct1D.reserve(maxRange);
    smooth1D.reserve(maxRange);

    /// First 1D smoothing in the 2nd direction
    exerFct1D.resize(range1);
    smooth1D.resize(range1);
    int pos;
    for(i0=0,index0.Reset();i0<range0;++i0,index0.NthIncr(0))
    {
        /// Get the 1D function to smooth
        for(i1=0,index1=index0,pos=index1.GetPosition();i1<range1;++i1,pos=index1.NthIncr(1))
        {
            indexPos[i1] = pos;
            exerFct1D[i1]=(*exerFct)[pos];
        }

        /// 1D smoothing
        ExerciseSmoothing1D(exerFct1D,coef,smooth1D,unusedExer1D);

        /// Save correction in 2D
        for(i1=0;i1<range1;++i1)
            (*smoothValues)[indexPos[i1]] += smooth1D[i1];
    }

    /// Second 1D smoothing in the 1st direction
    exerFct1D.resize(range0);
    smooth1D.resize(range0);
    for(i1=0,index1.Reset();i1<range1;++i1,index1.NthIncr(1))
    {
        /// Get the 1D function to smooth
        exerFct1D.resize(range0);
        smooth1D.resize(range0);
        for(i0=0,index0=index1,pos=index0.GetPosition();i0<range0;++i0,pos=index0.NthIncr(0))
        {
            indexPos[i0] = pos;
            exerFct1D[i0]=(*exerFct)[pos];
        }

        /// 1D smoothing
        ExerciseSmoothing1D(exerFct1D,coef,smooth1D,unusedExer1D);

        /// Save correction in 2D
        for(i0=0;i0<range0;++i0)
            (*smoothValues)[indexPos[i0]] += smooth1D[i0];
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree2D
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Tree2D::toString(const string& indent, const string& nextIndent ) const
{
    CC_Ostringstream os;
    os << indent << "\n\n";
    os << indent << "2D Tree\n";
    os << indent << "-------\n\n";

    os << ARM_TreeBase::toString(indent,nextIndent);

    return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/