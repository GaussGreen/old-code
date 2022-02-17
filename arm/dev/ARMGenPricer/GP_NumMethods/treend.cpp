/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree3d.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#include "gpnummethods/treeNd.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/slice.h"

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeND::ARM_TreeND( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, size_t dim, const ARM_SmootherBase* smoother, bool computeSpotProbas  )
:	ARM_TreeBase( sampler,truncator,reconnector,smoother,computeSpotProbas ), itsDim( dim )
{}



////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_TreeND::ARM_TreeND(const ARM_TreeND& rhs)
:   ARM_TreeBase(rhs), itsDim( rhs.itsDim )
{}


////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: Assignement operator
///	Returns: ARM_TreeND&
///	Action : 
////////////////////////////////////////////////////
ARM_TreeND& ARM_TreeND::operator= (const ARM_TreeND& rhs )
{
    if( this != &rhs )
	{
        ARM_TreeBase::operator =(rhs);
		itsDim = rhs.itsDim;
	}
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: ComputeStateVariables
///	Returns: List of state variables
///	Action : Compute the state variables for a
///          given position in the tree (defined
///          by its timeIdx)
////////////////////////////////////////////////////
void ARM_TreeND::ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const
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
        states = ARM_PricingStatesPtr( new ARM_PricingStates(nbStates,itsDim,0,itsDim) );
    else if(states->size() != nbStates || states->NumMethodStatesSize() != itsDim)
        states->resizeNumMethodStates(itsDim,nbStates);
	
	/// Add the drift correction to the first variable
    ARM_GP_VectorPtr driftCorrection = slice->GetDriftCorrectionVect();
    for(size_t stateIdx=0; stateIdx<nbStates; ++stateIdx )
    {
		for (size_t i=0; i<itsDim; ++i )
		    states->SetNumMethodState(stateIdx,i, (*xStates) (i,stateIdx)+ (*driftCorrection)[i]);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: ExerciseSmoothing
///	Returns: void
///	Action : Compute a smoothed exercise function
////////////////////////////////////////////////////
void ARM_TreeND::ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const
{
    /// No smoothing at the moment !
}


////////////////////////////////////////////////////
///	Class  : ARM_TreeND
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_TreeND::toString(const string& indent, const string& nextIndent ) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << itsDim << "D Tree\n";
    os << "-------\n\n";

    os << ARM_TreeBase::toString(indent,nextIndent);

    return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/