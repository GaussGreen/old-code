/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file transitor.cpp
 *
 *  \brief
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#include "gpnummethods/transitor.h"

/// gpinfra
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/slice.h"

CC_BEGIN_NAMESPACE( ARM )

#define DOWN    -1
#define MIDDLE  0
#define UP		1


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Transition2D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Transitor2D
///	Routine: ComputeTransitions
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Transitor2D::ComputeTransitions (const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* tmpNextSlice, ARM_TransitionStates& transStates) const
{
    const ARM_SliceNDBase* slice = tmpSlice->ToSliceNDBase();
    const ARM_SliceNDBase* nextSlice = tmpNextSlice->ToSliceNDBase();

    /// Get the tree index of the transition node in the next slice
    ARM_IntVector nextNodeIndex(2);
    nextNodeIndex[0] = slice->GetNextNodeIndex(stateIdx,0);
    nextNodeIndex[1] = slice->GetNextNodeIndex(stateIdx,1);

    ARM_TreeIndex nextIndex(nextSlice->GetMin(),nextSlice->GetMax());
    nextIndex.SetIndex(nextNodeIndex);

	int start0  = slice->GetProbaDown(stateIdx,0) > 0 ? DOWN    : MIDDLE;
	int end0    = slice->GetProbaUp(stateIdx,0) >   0 ?   UP    : MIDDLE;
	int start1  = slice->GetProbaDown(stateIdx,1) > 0 ? DOWN    : MIDDLE;
	int end1    = slice->GetProbaUp(stateIdx,1) >   0 ?   UP    : MIDDLE;

	size_t nbNodes = (end0 - start0 + 1)*(end1 - start1 + 1);

    transStates.SetUsedSize(nbNodes);

	int i0,i1;
    size_t refPos=nextIndex.GetPosition(), refPos0;
    size_t transIdx=0;
    size_t offset0 = 1;
    size_t offset1 = nextSlice->GetMax(0) - nextSlice->GetMin(0) + 1;
    double transProba, proba0;
    size_t transStateIdx;
	for(i0=start0; i0 <= end0; ++i0)
	{
        /// Get position and probas for 1st direction
		switch(i0)
        {
		case DOWN:
			refPos0 = refPos - offset0;
			proba0 = slice->GetProbaDown(stateIdx,0);
            break;

		case MIDDLE:
            refPos0 = refPos;
			proba0 = 1.0 - slice->GetProbaDown(stateIdx,0) - slice->GetProbaUp(stateIdx,0);
            break;

		case UP:
			refPos0 = refPos + offset0;
			proba0 = slice->GetProbaUp(stateIdx,0);
		    break;
        }

	    for(i1=start1; i1 <= end1; ++i1)
	    {
            /// Get position and probas for 2nd direction
		    switch(i1)
            {
		    case DOWN:
			    transStateIdx = refPos0 - offset1;
			    transProba = proba0 * slice->GetProbaDown(stateIdx,1);
                break;

		    case MIDDLE:
			    transStateIdx = refPos0;
			    transProba = proba0 * (1.0 - slice->GetProbaDown(stateIdx,1) - slice->GetProbaUp(stateIdx,1));
                break;

		    case UP:
			    transStateIdx = refPos0 + offset1;
			    transProba = proba0 * slice->GetProbaUp(stateIdx,1);
                break;
            }

            transStates.SetTransition(transIdx,transStateIdx,transProba);
            ++transIdx;

        } // 2nd direction
    } // 1st direction
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Transitor3D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Transitor3D
///	Routine: ComputeTransitions
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_Transitor3D::ComputeTransitions (const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* tmpNextSlice, ARM_TransitionStates& transStates) const
{
    const ARM_SliceNDBase* slice = tmpSlice->ToSliceNDBase();
    const ARM_SliceNDBase* nextSlice = tmpNextSlice->ToSliceNDBase();

    /// Get the tree index of the transition node in the next slice
    ARM_IntVector nextNodeIndex(3);
    nextNodeIndex[0] = slice->GetNextNodeIndex(stateIdx,0);
    nextNodeIndex[1] = slice->GetNextNodeIndex(stateIdx,1);
    nextNodeIndex[2] = slice->GetNextNodeIndex(stateIdx,2);

    ARM_TreeIndex nextIndex(nextSlice->GetMin(),nextSlice->GetMax());
    nextIndex.SetIndex(nextNodeIndex);

	int start0 = slice->GetProbaDown(stateIdx,0) > 0 ? DOWN : MIDDLE;
	int end0 = slice->GetProbaUp(stateIdx,0) > 0 ? UP : MIDDLE;
	int start1 = slice->GetProbaDown(stateIdx,1) > 0 ? DOWN : MIDDLE;
	int end1 = slice->GetProbaUp(stateIdx,1) > 0 ? UP : MIDDLE;
	int start2 = slice->GetProbaDown(stateIdx,2) > 0 ? DOWN : MIDDLE;
	int end2 = slice->GetProbaUp(stateIdx,2) > 0 ? UP : MIDDLE;

	size_t nbNodes = (end0 - start0 + 1)*(end1 - start1 + 1)*(end2 - start2 + 1);

    transStates.SetUsedSize(nbNodes);

	int i0,i1,i2;
    size_t refPos=nextIndex.GetPosition(),refPos0,refPos1;
    size_t transIdx=0;
    size_t offset0 = 1;
    size_t offset1 = nextSlice->GetMax(0) - nextSlice->GetMin(0) + 1;
    size_t offset2 = offset1 * (nextSlice->GetMax(1) - nextSlice->GetMin(1) + 1);
    double transProba,proba0,proba1;
    size_t transStateIdx;
	for(i0=start0; i0 <= end0; ++i0)
	{
        /// Get position and probas for 1st direction
		switch(i0)
        {
		case DOWN:
			refPos0 = refPos - offset0;
			proba0 = slice->GetProbaDown(stateIdx,0);
            break;

		case MIDDLE:
            refPos0 = refPos;
			proba0 = 1.0 - slice->GetProbaDown(stateIdx,0) - slice->GetProbaUp(stateIdx,0);
            break;

		case UP:
			refPos0 = refPos + offset0;
			proba0 = slice->GetProbaUp(stateIdx,0);
            break;
		}

	    for(i1=start1; i1 <= end1; ++i1)
	    {
            /// Get position and probas for 2nd direction
		    switch(i1)
            {
		    case DOWN:
			    refPos1 = refPos0 - offset1;
			    proba1 = proba0 * slice->GetProbaDown(stateIdx,1);
                break;

		    case MIDDLE:
                refPos1 = refPos0;
			    proba1 = proba0 * (1.0 - slice->GetProbaDown(stateIdx,1) - slice->GetProbaUp(stateIdx,1));
                break;

		    case UP:
			    refPos1 = refPos0 + offset1;
			    proba1 = proba0 * slice->GetProbaUp(stateIdx,1);
                break;
		    }

	        for(i2=start2; i2 <= end2; ++i2)
	        {
                /// Get position and probas for 3rd direction
		        switch(i2)
                {
		        case DOWN:
			        transStateIdx = refPos1 - offset2;
			        transProba = proba1 * slice->GetProbaDown(stateIdx,2);
                    break;

		        case MIDDLE:
			        transStateIdx = refPos1;
			        transProba = proba1 * (1.0 - slice->GetProbaDown(stateIdx,2) - slice->GetProbaUp(stateIdx,2));
                    break;

		        case UP:
			        transStateIdx = refPos1 + offset2;
			        transProba = proba1 * slice->GetProbaUp(stateIdx,2);
                    break;
		        }

                transStates.SetTransition(transIdx,transStateIdx,transProba);
                ++transIdx;

            } // 3rd direction
        } // 2nd direction
    } // 1st direction
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_TransitorND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_TransitorND
///	Routine: ComputeTransitions
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_TransitorND::ComputeTransitions (const ARM_SliceBase* tmpSlice, size_t stateIdx, const ARM_SliceBase* tmpNextSlice, ARM_TransitionStates& transStates) const
{
    const ARM_SliceNDBase* slice = tmpSlice->ToSliceNDBase();
    const ARM_SliceNDBase* nextSlice = tmpNextSlice->ToSliceNDBase();

    /// Get the tree index of the transition node in the next slice
    ARM_IntVector nextNodeIndex(itsDim);
	size_t i,j;
	int k;
	for( i=0; i<itsDim; ++i )
		nextNodeIndex[i] = slice->GetNextNodeIndex(stateIdx,i);
	
    ARM_TreeIndex nextIndex(nextSlice->GetMin(),nextSlice->GetMax());
    nextIndex.SetIndex(nextNodeIndex);
	
	ARM_IntVector down(itsDim,0), up(itsDim,0);
	ARM_IntVector OffsetVec(itsDim, 1);
	ARM_IntVector iVec( itsDim );
	size_t nbNodes = 1;
    size_t transIdx=0;
	
	for( i=0; i<itsDim; ++i )
	{
		down[i]  = slice->GetProbaDown(stateIdx,i)> 0 ? -1 : 0;
		up[i]	 = slice->GetProbaUp(stateIdx,i)  > 0 ?  1 : 0;
		nbNodes *= (up[i] - down[i] + 1);
		iVec[i]  = down[i];
		
		if( 0==i)
			OffsetVec[i] = 1;
		else
			OffsetVec[i] = OffsetVec[i-1]*(nextSlice->GetMax(i-1) - nextSlice->GetMin(i-1) + 1);
	}
	
    transStates.SetUsedSize(nbNodes);
	
	ARM_IntVector refPosVec(itsDim, 1);
	ARM_GP_Vector proba( itsDim, 0.0);
	
	for( i=0; i<nbNodes; ++i )
	{
		switch( iVec[0] )
		{
		case DOWN:
			refPosVec[0]= nextIndex.GetPosition() - OffsetVec[0];
			proba[0]	= slice->GetProbaDown(stateIdx,0);
			break;
		case MIDDLE:
			refPosVec[0]= nextIndex.GetPosition();
			proba[0]	= 1.0 - slice->GetProbaDown(stateIdx,0) - slice->GetProbaUp(stateIdx,0);
			break;
		case UP:
			refPosVec[0]= nextIndex.GetPosition() + OffsetVec[0];
			proba[0]	= slice->GetProbaUp(stateIdx,0);
			break;
		}	
		
		for( j=1; j<itsDim; ++j )
		{
			switch( iVec[j] )
			{
			case DOWN:
				refPosVec[j] = refPosVec[j-1] - OffsetVec[j];
				proba[j]	 = proba[j-1] * slice->GetProbaDown(stateIdx,j);
				break;
			case MIDDLE:
				refPosVec[j] = refPosVec[j-1];
				proba[j]	 = proba[j-1] * ( 1.0 - slice->GetProbaDown(stateIdx,j) - slice->GetProbaUp(stateIdx,j) );
				break;
			case UP:
				refPosVec[j] = refPosVec[j-1] + OffsetVec[j];
				proba[j]	 = proba[j-1] * slice->GetProbaUp(stateIdx,j);
				break;
			}	
		}
		
        transStates.SetTransition(transIdx,refPosVec[itsDim-1],proba[itsDim-1]);
		++transIdx;
		
		k = itsDim-1;
		bool addOneTo_iVec = true;
		while( addOneTo_iVec )
		{
			++iVec[k];
			if( iVec[k] > up[k] )
			{
				iVec[k] = down[k];
				--k;
				if( k<0)
					addOneTo_iVec = false;
			}
			else
				addOneTo_iVec = false;	
		}
	}
}

#undef DOWN
#undef MIDDLE
#undef UP

CC_END_NAMESPACE()


#undef FIRST_STATE_VARIABLE

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

