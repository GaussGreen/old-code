/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tree1d.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/tree1d.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/slice.h"

#define FIRST_STATE_VARIABLE    0

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_Tree1D
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree1D::ARM_Tree1D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas )
:	ARM_TreeBase( sampler,truncator,reconnector,smoother,computeSpotProbas )
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree1D
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree1D::ARM_Tree1D(const ARM_Tree1D& rhs)
:   ARM_TreeBase(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree1D
///	Routine: Assignement operator
///	Returns: ARM_Tree1D&
///	Action : 
////////////////////////////////////////////////////
ARM_Tree1D& ARM_Tree1D::operator= (const ARM_Tree1D& rhs )
{
    if( this != &rhs )
        ARM_TreeBase::operator =(rhs);
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree1D
///	Routine: ComputeStateVariables
///	Returns: List of state variables
///	Action : Compute the state variables for a
///          given position in the tree (defined
///          by its timeIdx)
////////////////////////////////////////////////////
void ARM_Tree1D::ComputeStateVariables(int timeIdx, ARM_PricingStatesPtr& states) const
{
    ARM_Slice1DBase* slice = (*GetSlices())[timeIdx]->ToSlice1DBase();
    size_t nbStates = slice->size();

    /// Restore X states from the slice
    ARM_GP_VectorPtr xStates = slice->GetXStates();
    if(xStates == ARM_GP_VectorPtr(NULL))
    {
        /// Not previously computed then compute Z states...
        ARM_GP_VectorPtr zStates = slice->GetZStates();

        /// ... and restore X states through the sampler
        const ARM_Sampler1DBase* sampler = GetSampler()->ToSampler1DBase();
        xStates = sampler->ComputeZtoXStates(timeIdx,zStates);
    }

    /// Set numerical method states = X states
    if(states == ARM_PricingStatesPtr(NULL))
        states = ARM_PricingStatesPtr( new ARM_PricingStates(nbStates,1,0,1) );
    else if(states->size() != nbStates || states->NumMethodStatesSize() != 1)
        states->resizeNumMethodStates(1,nbStates);
	
	/// add the drift correction
    for(size_t i=0; i<nbStates; ++i )
        states->SetNumMethodState(i,FIRST_STATE_VARIABLE,(*xStates)[i]+slice->GetDriftCorrection());
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree1D
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Tree1D::toString(const string& indent, const string& nextIndent ) const
{
    CC_Ostringstream os;
    os << indent << "\n\n";
    os << indent << "1D Tree\n";
    os << indent << "-------\n\n";

    os << ARM_TreeBase::toString(indent,nextIndent);

    return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/