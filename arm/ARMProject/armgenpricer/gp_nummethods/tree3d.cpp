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
 *	\date November 2004
 */


#include "gpnummethods/tree3d.h"

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
#define THIRD_STATE_VARIABLE    2

CC_BEGIN_NAMESPACE( ARM )

int ARM_Tree3D::SmoothingDirection[3][3]={{0,1,2},{2,0,1},{1,2,0}};

////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: Default constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree3D::ARM_Tree3D( const ARM_SamplerBase* sampler, const ARM_TruncatorBase* truncator, const ARM_ReconnectorBase* reconnector, const ARM_SmootherBase* smoother, bool computeSpotProbas )
:	ARM_TreeBase( sampler,truncator,reconnector,smoother,computeSpotProbas )
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: Copy constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_Tree3D::ARM_Tree3D(const ARM_Tree3D& rhs)
:   ARM_TreeBase(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: Assignement operator
///	Returns: ARM_Tree3D&
///	Action : 
////////////////////////////////////////////////////
ARM_Tree3D& ARM_Tree3D::operator= (const ARM_Tree3D& rhs )
{
    if( this != &rhs )
        ARM_TreeBase::operator =(rhs);
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: ComputeStateVariables
///	Returns: List of state variables
///	Action : Compute the state variables for a
///          given position in the tree (defined
///          by its timeIdx)
////////////////////////////////////////////////////
void ARM_Tree3D::ComputeStateVariables( int timeIdx, ARM_PricingStatesPtr& states ) const
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
        states = ARM_PricingStatesPtr( new ARM_PricingStates(nbStates,3,0,3) );
    else if(states->size() != nbStates || states->NumMethodStatesSize() != 3)
        states->resizeNumMethodStates(3,nbStates);
	
	/// Add the drift correction to the first variable
    ARM_GP_VectorPtr driftCorrection = slice->GetDriftCorrectionVect();
    for(size_t stateIdx=0; stateIdx<nbStates; ++stateIdx )
    {
        states->SetNumMethodState(stateIdx,FIRST_STATE_VARIABLE,(*xStates)(FIRST_STATE_VARIABLE,stateIdx)+(*driftCorrection)[FIRST_STATE_VARIABLE]);
        states->SetNumMethodState(stateIdx,SECOND_STATE_VARIABLE,(*xStates)(SECOND_STATE_VARIABLE,stateIdx)+(*driftCorrection)[SECOND_STATE_VARIABLE]);
        states->SetNumMethodState(stateIdx,THIRD_STATE_VARIABLE,(*xStates)(THIRD_STATE_VARIABLE,stateIdx)+(*driftCorrection)[THIRD_STATE_VARIABLE]);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: ExerciseSmoothing
///	Returns: void
///	Action : Compute a smoothed exercise function
////////////////////////////////////////////////////
void ARM_Tree3D::ExerciseSmoothing(const ARM_GP_VectorPtr& exerFct, ARM_GP_VectorPtr& smoothValues, ARM_GP_VectorPtr& exerStates) const
{
    int timeIdx = GetLastTimeIdx();
    const ARM_SliceNDBase* slice = (*(GetSlices()))[timeIdx]->ToSliceND();
    const ARM_IntVector& minIdx = slice->GetMin();
    const ARM_IntVector& maxIdx = slice->GetMax();

    vector< ARM_TreeIndex > index(3);
    index[0]=ARM_TreeIndex(minIdx,maxIdx);

    ARM_IntVector range(3);
    int maxRange=0;
    for(size_t i=0;i<3;++i)
    {
        range[i] = maxIdx[i]-minIdx[i]+1;
        if(range[i] > maxRange)
            maxRange=range[i];
    }

    /// For a given the 3rd choice (#1 or #2 or #3 direction),
    /// we chose only 1 among 2 possible permutations.
    /// 3D slice is smoothed by a 1D smoothing for each point of space (#1,#2)
    /// then (#3,#1) then (#2,#3) (refer to SmoothingDirection static array)
    /// Global correction is divided by 3
    size_t nbPerm=sizeof(SmoothingDirection)/sizeof(SmoothingDirection[0]);
    double coef = 1.0/nbPerm;

#ifdef __GP_STRICT_VALIDATION
    if(exerFct->size() != range[0]*range[1]*range[2])
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : slice size mismatch in Tree3D" );
#endif

    std::vector<double> exerFct1D,smooth1D,unusedExer1D;
    ARM_IntVector indexPos(maxRange);

    exerFct1D.reserve(maxRange);
    smooth1D.reserve(maxRange);

    size_t permIdx,i0,i1,i2;
    int dir0,dir1,dir2,pos;
    for(permIdx=0;permIdx<nbPerm;++permIdx)
    {
        dir0 = SmoothingDirection[permIdx][0];
        dir1 = SmoothingDirection[permIdx][1];
        dir2 = SmoothingDirection[permIdx][2];
        exerFct1D.resize(range[dir2]);
        smooth1D.resize(range[dir2]);
        for(i0=0,index[0].Reset();i0<range[dir0];++i0,index[0].NthIncr(dir0))
        {
            for(i1=0,index[1]=index[0];i1<range[dir1];++i1,index[1].NthIncr(dir1))
            {
                /// Get the 1D function to smooth in the dir2 th direction
                for(i2=0,index[2]=index[1],pos=index[2].GetPosition();i2<range[dir2];++i2,pos=index[2].NthIncr(dir2))
                {
                    indexPos[i2] = pos;
                    exerFct1D[i2]=(*exerFct)[pos];
                }

                /// 1D smoothing
                ExerciseSmoothing1D(exerFct1D,coef,smooth1D,unusedExer1D);

                /// Save correction in 3D
                for(i2=0;i2<range[dir2];++i2)
                    (*smoothValues)[indexPos[i2]] += smooth1D[i2];

            } // for dir1 states
        } // for dir2 states
    } // for permutation nb
}


////////////////////////////////////////////////////
///	Class  : ARM_Tree3D
///	Routine: toString
///	Returns: string
///	Action : stringify the object
////////////////////////////////////////////////////
string ARM_Tree3D::toString(const string& indent, const string& nextIndent ) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << "3D Tree\n";
    os << "-------\n\n";

    os << ARM_TreeBase::toString(indent,nextIndent);

    return os.str();
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/