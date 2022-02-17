/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file truncator.cpp
 *
 *  \brief
 *
 *	\author  JM Prie E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/truncator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/utilityport.h"  /// for CC_Round

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodel.h"

/// gpnummethods
#include "gpnummethods/slice.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/treeindex.h"


#include <iomanip>
CC_USING_NS(std,dec)
CC_USING_NS(std,fixed)
CC_USING_NS(std,setw)
CC_USING_NS(std,setprecision)
CC_USING_NS(std,scientific)

#include <cmath>

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_TruncatorBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

const size_t ARM_TruncatorBase::MAX_NBSTEP_EXPANSION [] = {100,10,3,2};
const int ARM_TruncatorBase::NO_MAXMAX_INDEX            = 10000;

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorBase
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_TruncatorBase::ARM_TruncatorBase( const ARM_TruncatorBase& rhs )
:   ARM_RootObject(rhs),
    itsReconnector( rhs.itsReconnector ? static_cast< ARM_ReconnectorBase* >(rhs.itsReconnector->Clone()) : NULL )
{}


ARM_TruncatorBase& ARM_TruncatorBase::operator=(const ARM_TruncatorBase& rhs )
{
    if( this != &rhs )
    {
        ARM_RootObject::operator =( rhs );
        delete itsReconnector;
        itsReconnector = rhs.itsReconnector ? static_cast< ARM_ReconnectorBase* >(rhs.itsReconnector->Clone()) : NULL;
    }
    return *this;
}

ARM_TruncatorBase::~ARM_TruncatorBase()
{
    delete itsReconnector;
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorBase
///	Routine: ToTruncator???
///	Returns: ARM_Truncator1D*,ARM_TruncatorND*,
///          ToTruncator1DArrowDebreu*,ToTruncatorNDArrowDebreu*
///	Action : downcast
////////////////////////////////////////////////////

ARM_Truncator1D* ARM_TruncatorBase::ToTruncator1D()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Truncator1D" );
}

const ARM_Truncator1D* ARM_TruncatorBase::ToTruncator1D() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Truncator1D" );
}

ARM_TruncatorND* ARM_TruncatorBase::ToTruncatorND()
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to TruncatorND" );
}

const ARM_TruncatorND* ARM_TruncatorBase::ToTruncatorND() const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to TruncatorND" );
}

ARM_Truncator1DArrowDebreu* ARM_TruncatorBase::ToTruncator1DArrowDebreu() 
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Truncator1DArrowDebreu" );
}

const ARM_Truncator1DArrowDebreu* ARM_TruncatorBase::ToTruncator1DArrowDebreu()  const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to Truncator1DArrowDebreu" );
}

ARM_TruncatorNDArrowDebreu* ARM_TruncatorBase::ToTruncatorNDArrowDebreu() 
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to TruncatorNDArrowDebreu" );
}

const ARM_TruncatorNDArrowDebreu* ARM_TruncatorBase::ToTruncatorNDArrowDebreu()  const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot cast to TruncatorNDArrowDebreu" );
}

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorBase
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
void ARM_TruncatorBase::SetReconnector(const ARM_ReconnectorBase* reconnector)
{
	delete itsReconnector;
    itsReconnector = static_cast< ARM_ReconnectorBase* >(reconnector->Clone());
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorBase
///	Routine: CorrespondingTruncator1D
///	Returns: 
///	Action : creates the corresponding truncator 1D
////////////////////////////////////////////////////
ARM_Truncator1D* ARM_TruncatorBase::CorrespondingTruncator1D(size_t dim) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": cannot create the corresponding truncator 1D" );
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Truncator1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_Truncator1D::ARM_Truncator1D( const ARM_Truncator1D& rhs )
:   ARM_TruncatorBase( rhs ), 
    itsStdDevRatio(rhs.itsStdDevRatio),
    itsMinIndex(rhs.itsMinIndex),
    itsMaxIndex(rhs.itsMaxIndex)
{}


ARM_Truncator1D& ARM_Truncator1D::operator=(const ARM_Truncator1D& rhs )
{
    if( this != &rhs )
    {
        ARM_TruncatorBase::operator =( rhs );
        itsStdDevRatio=rhs.itsStdDevRatio;
        itsMinIndex=rhs.itsMinIndex;
        itsMaxIndex=rhs.itsMaxIndex;
    }
    return *this;
}

ARM_Truncator1D::~ARM_Truncator1D()
{}

////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine: ComputeCenter
///	Returns: ARM_IntVectorPtr
///	Action : Compute center indexes of the diffusion
////////////////////////////////////////////////////
ARM_IntMatrixPtr ARM_Truncator1D::ComputeCenter(const ARM_SamplerBase* sampler, const ARM_SliceVectorPtr& slices) const
{
    const ARM_PricingModel& model = * sampler->GetModel();
    size_t nbSlices = slices->size();
    ARM_IntMatrixPtr center( new ARM_IntMatrix(nbSlices,1,0) );

    ARM_GP_MatrixPtr globalDrifts;
	model.IntegratedGlobalDrifts(model.GetNumMethod()->GetTimeSteps()->GetValues(),globalDrifts);

    if(globalDrifts != ARM_GP_MatrixPtr(NULL))
    {
        /// Convert drift from X to Z space
        sampler->ComputeXtoZCenters(globalDrifts);

        /// Rounded center to indexes in the discrete Z space
        for(size_t i=0;i<nbSlices;++i)
        {
            const ARM_Slice1DBase* slice1D = (*slices)[i]->ToSlice1DBase();
            (*center)(i,0) = CC_Round((*globalDrifts)(i,0)/slice1D->GetSpaceStep());
        }
    }

    return center;
}

////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine: Init
///	Returns: void
///	Action : Initialisation of truncation decision datas
////////////////////////////////////////////////////
void ARM_Truncator1D::Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center)
{
    const ARM_Sampler1DBase* sampler1D = sampler->ToSampler1DBase();
    double stdDev;
    int max,lastMax=1;
    size_t nbSlices = slices->size();
    itsMinIndex.resize(nbSlices);
    itsMinIndex[0]=0;
    itsMaxIndex.resize(nbSlices);
    itsMaxIndex[0]=0;

    /// Compute diffusion center if necessary
    ARM_IntMatrixPtr localCenter(center);
    if(localCenter == ARM_IntMatrixPtr(NULL))
        localCenter = ComputeCenter(sampler,slices);

    double theoMax;
    int maxMaxIndex = GetMaxMaxIndex();

    ARM_Slice1DBase* slice1D;
    for(size_t i=1;i<nbSlices;++i)
    {
        slice1D = (*slices)[i]->ToSlice1DBase();
        stdDev = sqrt(sampler1D->GetGlobalVar(i));
        max  = static_cast<int>(floor(stdDev*itsStdDevRatio/slice1D->GetSpaceStep())+1);
        if(max > maxMaxIndex && max > (theoMax=lastMax+ARM_TruncatorBase::MAX_NBSTEP_EXPANSION[0]))
        {
            /// Increase space step to avoid uncontrolled expansion
            /// Connections inside tree may be degenerated in binomial
            slice1D->SetSpaceStep(stdDev*itsStdDevRatio/theoMax);
        }
        itsMinIndex[i] = (*localCenter)(i,0) - max;
        itsMaxIndex[i] = (*localCenter)(i,0) + max;
        lastMax = max;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine: NeedTruncation
///	Returns: void
///	Action : detect if truncation is needed
////////////////////////////////////////////////////
bool ARM_Truncator1D::NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, int nextNodeIdx) const
{
    size_t nextSliceIdx = tmpNextSlice->GetIndex();

    return nextNodeIdx <= itsMinIndex[nextSliceIdx] || nextNodeIdx >= itsMaxIndex[nextSliceIdx];
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine: UpdateNode
///	Returns: void
///	Action : truncates connection if necessary
///          Slices are necessary 1D slices with explicit
///          nodes (ARM_Slice1D)
////////////////////////////////////////////////////
void ARM_Truncator1D::UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, double meanRel, double varRel ) const
{
    ARM_Slice1D* slice      = tmpSlice->ToSlice1D();
    ARM_Slice1D* nextSlice  = tmpNextSlice->ToSlice1D();

    size_t nextSliceIdx = nextSlice->GetIndex();
    int nextNodeIdx     = slice->GetNextNodeIndex(stateIdx);
    double probaUp      = slice->GetProbaUp(stateIdx);
    double probaDown    = slice->GetProbaDown(stateIdx);
    int center          = (itsMaxIndex[nextSliceIdx]+itsMinIndex[nextSliceIdx])/2;
    bool isReconnected  = false;
    if(nextNodeIdx <= center)
    {
        if( nextNodeIdx < itsMinIndex[nextSliceIdx] )
        {
            /// Two nodes away => out of lattice lower bound
            GetReconnector()->ReconnectOutOfBoundDown( nextNodeIdx, itsMinIndex[nextSliceIdx], meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }

        else if( nextNodeIdx == itsMinIndex[nextSliceIdx] )
        {
                /// Central connection on the lattice lower edge
            GetReconnector()->ReconnectOnTheHedgeDown( nextNodeIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }
    }
    else if(nextNodeIdx > center)
    {
        if( nextNodeIdx > itsMaxIndex[nextSliceIdx] )
        {
            /// Two nodes away => out of lattice upper bound
            GetReconnector()->ReconnectOutOfBoundUp( nextNodeIdx, itsMaxIndex[nextSliceIdx], meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }

        else if( nextNodeIdx == itsMaxIndex[nextSliceIdx] )
        {
                /// Central connection on the lattice upper edge
            GetReconnector()->ReconnectOnTheHedgeUp( nextNodeIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }
    }

    if(!isReconnected && (probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0))
        /// Last chance : try a binomial reconnection inside the tree
        GetReconnector()->ReconnectInside( nextNodeIdx, meanRel, varRel, probaUp, probaDown );

    if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
    {
        /// Nothing else to do !!
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in the tree" );
    }
    else
    {
        slice->SetNextNodeIndex(stateIdx,nextNodeIdx);
        slice->SetProbaUp(stateIdx,probaUp);
        slice->SetProbaDown(stateIdx,probaDown);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1D
///	Routine: toString
///	Returns: 
///	Action : object dump
////////////////////////////////////////////////////
string ARM_Truncator1D::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << "1D StdDev Truncator :  MaxStdDev = " << fixed << itsStdDevRatio << "\n";
	return os.str();
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_Truncator1DArrowDebreu
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_Truncator1DArrowDebreu
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_Truncator1DArrowDebreu::ARM_Truncator1DArrowDebreu( const ARM_Truncator1DArrowDebreu& rhs )
:   ARM_Truncator1D( rhs ), 
    itsMaxMaxIndex(rhs.itsMaxMaxIndex),
    itsArrowDebreuThreshold(rhs.itsArrowDebreuThreshold),
    itsTimeThreshold(rhs.itsTimeThreshold),
    itsCenter(ARM_IntMatrixPtr(rhs.itsCenter != ARM_IntMatrixPtr(NULL) ? (ARM_IntMatrix*) rhs.itsCenter->Clone() : NULL))
{}


ARM_Truncator1DArrowDebreu& ARM_Truncator1DArrowDebreu::operator=(const ARM_Truncator1DArrowDebreu& rhs )
{
    if( this != &rhs )
    {
        ARM_Truncator1D::operator =( rhs );
        itsMaxMaxIndex=rhs.itsMaxMaxIndex;
        itsArrowDebreuThreshold=rhs.itsArrowDebreuThreshold;
        itsTimeThreshold=rhs.itsTimeThreshold;
        itsCenter = ARM_IntMatrixPtr(rhs.itsCenter != ARM_IntMatrixPtr(NULL) ? (ARM_IntMatrix*) rhs.itsCenter->Clone() : NULL) ;
    }
    return *this;
}

ARM_Truncator1DArrowDebreu::~ARM_Truncator1DArrowDebreu()
{}

////////////////////////////////////////////////////
///	Class  : ARM_Truncator1DArrowDebreu
///	Routine: Init
///	Returns: void
///	Action : initialisation of truncation decision datas
////////////////////////////////////////////////////
void ARM_Truncator1DArrowDebreu::Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center)
{
    /// Compute diffusion center if necessary
    if(center == ARM_IntMatrixPtr(NULL))
        itsCenter = ComputeCenter(sampler,slices);
    else
        itsCenter = center;

    /// Initialise the StdDev truncator to get min and max indexes
    ARM_Truncator1D::Init(slices,sampler,itsCenter);
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1DArrowDebreu
///	Routine: NeedTruncation
///	Returns: void
///	Action : detect if truncation is needed
////////////////////////////////////////////////////
bool ARM_Truncator1DArrowDebreu::NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, int nextNodeIdx) const
{
    size_t nextSliceIdx = tmpNextSlice->GetIndex();

    bool isADTrunc = sliceTime > itsTimeThreshold ||
                     ( tmpSlice->GetArrowDebreuPrices() != ARM_GP_VectorPtr(NULL) &&
                      tmpSlice->GetArrowDebreuPrices(stateIdx) < itsArrowDebreuThreshold );

    int nextMinIdx = GetMinIndex(nextSliceIdx);
    int center = (*itsCenter)(nextSliceIdx,0);
    if(isADTrunc && nextMinIdx < center-itsMaxMaxIndex)
        nextMinIdx = center-itsMaxMaxIndex;

    if(nextNodeIdx <= nextMinIdx)
        return true;

    else
    {
        int nextMaxIdx = GetMaxIndex(nextSliceIdx);
        if(isADTrunc && nextMaxIdx > center+itsMaxMaxIndex)
            nextMaxIdx = center+itsMaxMaxIndex;

        return nextNodeIdx >= nextMaxIdx;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1DArrowDebreu
///	Routine: UpdateNode
///	Returns: void
///	Action : truncates if necessary the node
///          Slices are necessary 1D slices with explicit
///          nodes (ARM_Slice1D)
////////////////////////////////////////////////////
void ARM_Truncator1DArrowDebreu::UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, double meanRel, double varRel ) const
{
    ARM_Slice1D* slice      = tmpSlice->ToSlice1D();
    ARM_Slice1D* nextSlice  = tmpNextSlice->ToSlice1D();

    size_t nextSliceIdx = nextSlice->GetIndex();
    int nextNodeIdx     = slice->GetNextNodeIndex(stateIdx);

    /// Activation of AD truncation
    bool isADTrunc = sliceTime > itsTimeThreshold ||
                    (  slice->GetArrowDebreuPrices() != ARM_GP_VectorPtr(NULL) &&
                      slice->GetArrowDebreuPrices(stateIdx) < itsArrowDebreuThreshold );

    int nextMinIdx,nextMaxIdx;
    double probaUp      = slice->GetProbaUp(stateIdx);
    double probaDown    = slice->GetProbaDown(stateIdx);
    int center          = (*itsCenter)(nextSliceIdx,0);
    bool isReconnected  = false;

    if(nextNodeIdx <= center)
    {
        nextMinIdx = GetMinIndex(nextSliceIdx);
        if(isADTrunc && nextMinIdx < center-itsMaxMaxIndex)
            nextMinIdx = center-itsMaxMaxIndex;

        if( nextNodeIdx < nextMinIdx)
        {
            /// Two nodes away => out of lattice lower bound
            GetReconnector()->ReconnectOutOfBoundDown( nextNodeIdx, nextMinIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }

        else if( nextNodeIdx == nextMinIdx )
        {
            /// Central connection on the lattice lower edge
            GetReconnector()->ReconnectOnTheHedgeDown( nextNodeIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }
    }
    else if(nextNodeIdx > center)
    {
        nextMaxIdx = GetMaxIndex(nextSliceIdx);
        if(isADTrunc && nextMaxIdx > center+itsMaxMaxIndex)
            nextMaxIdx = center+itsMaxMaxIndex;

        if( nextNodeIdx > nextMaxIdx )
        {
            /// Two nodes away => out of lattice upper bound
            GetReconnector()->ReconnectOutOfBoundUp( nextNodeIdx, nextMaxIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }

        else if( nextNodeIdx == nextMaxIdx )
        {
            /// Central connection on the lattice upper edge
            GetReconnector()->ReconnectOnTheHedgeUp( nextNodeIdx, meanRel, varRel, probaUp, probaDown );
            isReconnected = true;
        }
    }


    if(!isReconnected && (probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0))
        /// Last chance : try a binomial reconnection inside the tree
        GetReconnector()->ReconnectInside( nextNodeIdx, meanRel, varRel, probaUp, probaDown );


    if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
    {
        /// Nothing else to do !!
	    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in the tree" );
    }
    else
    {
        slice->SetNextNodeIndex(stateIdx,nextNodeIdx);
        slice->SetProbaUp(stateIdx,probaUp);
        slice->SetProbaDown(stateIdx,probaDown);
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_Truncator1DArrowDebreu
///	Routine: toString
///	Returns: 
///	Action : object dump
////////////////////////////////////////////////////
string ARM_Truncator1DArrowDebreu::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
	os << indent << "1D Arrow-Debreu Truncator :\n";
    os << indent << "    MaxStdDev   = " << fixed << setw(5) << setprecision(2) << GetStdDevRatio() << "\n";
    os << indent << "    MaxMaxIndex = " << dec   << setw(5) << itsMaxMaxIndex << "\n";
    os << indent << "    ADLimit     = " << scientific << itsArrowDebreuThreshold << "\n";
    os << indent << "    TimeLimit   = " << fixed << setw(5) << setprecision(2) << itsTimeThreshold << "\n";
	return os.str();
}





////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_TruncatorND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: ComputeCenter
///	Returns: ARM_IntVectorPtr
///	Action : Compute center indexes of the diffusion
////////////////////////////////////////////////////
ARM_IntMatrixPtr ARM_TruncatorND::ComputeCenter(const ARM_SamplerBase* sampler, const ARM_SliceVectorPtr& slices) const
{
    const ARM_PricingModel& model = * sampler->GetModel();
    size_t nbSlices = slices->size(), nbDims=dim();
    ARM_IntMatrixPtr center( new ARM_IntMatrix(nbSlices,nbDims,0) );

    ARM_GP_MatrixPtr globalDrifts;
	model.IntegratedGlobalDrifts(model.GetNumMethod()->GetTimeSteps()->GetValues(),globalDrifts);

    if(globalDrifts != ARM_GP_MatrixPtr(NULL))
    {
        /// Convert drift from X to Z space
        sampler->ComputeXtoZCenters(globalDrifts);

        /// Rounded center to indexes in the discrete Z space
        for(size_t i=0;i<nbSlices;++i)
        {
            const ARM_SliceNDBase* sliceND = (*slices)[i]->ToSliceNDBase();
            for(size_t j=0;j<nbDims;++j)
                (*center)(i,j) = CC_Round((*globalDrifts)(i,j)/sliceND->GetSpaceStep(j));
        }
    }

    return center;
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: Init
///	Returns: void
///	Action : initialisation of truncation decision datas
////////////////////////////////////////////////////
void ARM_TruncatorND::Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center)
{
    const ARM_SamplerNDBase* samplerND = sampler->ToSamplerNDBase();
	size_t nbDims = dim();
    size_t nbSlices = slices->size();
	size_t i,j;
	double stdDev;
    int max;
    ARM_IntVector lastMax(nbDims,1);

    itsMinIndex.resize(nbSlices,nbDims);
    itsMaxIndex.resize(nbSlices,nbDims);
	for( j=0; j<nbDims; ++j )
    {
        itsMinIndex(0,j)=0;
        itsMaxIndex(0,j)=0;
    }

    /// Compute diffusion center if necessary
    ARM_IntMatrixPtr localCenter(center);
    if(localCenter == ARM_IntMatrixPtr(NULL))
        localCenter = ComputeCenter(sampler,slices);

    double theoMax;
    ARM_IntVector maxMaxIndex;
    GetMaxMaxIndex(maxMaxIndex);

    ARM_SliceNDBase* sliceND;
    int maxExpansion = ARM_TruncatorBase::MAX_NBSTEP_EXPANSION[nbDims<=4 ? nbDims-1 : 3];
    for(i=1; i<nbSlices; ++i)
    {
		sliceND = (*slices)[i]->ToSliceNDBase();
		const std::vector<double>& globalVarVec = samplerND->GetGlobalVar(i);
#if defined(__GP_STRICT_VALIDATION)
		if( nbDims != globalVarVec.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": nbDims != globalVarVec.size()" );
		if( nbDims != itsStdDevRatio.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": nbDims != itsStdDevRatio.size()" );
#endif

		/// loop over the dimension
		for( j=0; j<nbDims; ++j )
		{
			stdDev	= sqrt(globalVarVec[j]);
			max	= static_cast<int>(floor(stdDev*itsStdDevRatio[j]/sliceND->GetSpaceStep(j))+1);
            if(max > maxMaxIndex[j] && max > (theoMax=lastMax[j]+maxExpansion))
            {
                /// Increase space step to avoid uncontrolled expansion
                /// Connections inside tree may be degenerated in binomial
                sliceND->SetSpaceStep(j,stdDev*itsStdDevRatio[j]/theoMax);
            }
			itsMinIndex(i,j) = (*localCenter)(i,j) - max;
			itsMaxIndex(i,j) = (*localCenter)(i,j) + max;
            lastMax[j] = max;
		}
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: NeedTruncation
///	Returns: void
///	Action : detect if truncation is needed
////////////////////////////////////////////////////
ARM_BoolVector ARM_TruncatorND::NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const ARM_TreeIndex& nextNodeIdx) const
{
    size_t nbDims = nextNodeIdx.size();

    size_t nextSliceIdx = tmpNextSlice->GetIndex();

    ARM_BoolVector status(nbDims);
	for(size_t i=0; i<nbDims; ++i )
        status[i] = nextNodeIdx[i] <= itsMinIndex(nextSliceIdx,i) || nextNodeIdx[i] >= itsMaxIndex(nextSliceIdx,i);

    return status;
}

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: UpdateNode
///	Returns: void
///	Action : truncates if necessary the node
///          Slices are necessary ND slices with explicit
///          nodes (ARM_SliceND)
////////////////////////////////////////////////////
void ARM_TruncatorND::UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const std::vector<double>& meanRel, const std::vector<double>& varRel ) const
{
    ARM_SliceND* slice      = tmpSlice->ToSliceND();
	size_t nbDims		    = slice->dim();
    ARM_SliceND* nextSlice  = tmpNextSlice->ToSliceND();
    size_t nextSliceIdx		= nextSlice->GetIndex();
	
    double probaUp,probaDown;
    int nextNodeIdx,center;
    bool isReconnected;

	for(size_t i=0; i<nbDims; ++i )
	{
        nextNodeIdx     = slice->GetNextNodeIndex(stateIdx,i);
        probaUp         = slice->GetProbaUp(stateIdx,i);
        probaDown       = slice->GetProbaDown(stateIdx,i);
        center          = (itsMaxIndex(nextSliceIdx,i)+itsMinIndex(nextSliceIdx,i))/2;
        isReconnected   = false;
		if( nextNodeIdx <= center)
		{
			if( nextNodeIdx < itsMinIndex(nextSliceIdx,i) )
            {
                /// Two nodes away => out of lattice lower bound
				GetReconnector()->ReconnectOutOfBoundDown( nextNodeIdx, itsMinIndex(nextSliceIdx,i), meanRel[i], varRel[i], probaUp, probaDown);
                isReconnected = true;
            }

			else if( nextNodeIdx == itsMinIndex(nextSliceIdx,i) )
            {
                /// Central connection on the lattice lower edge
				GetReconnector()->ReconnectOnTheHedgeDown( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }
		}
		else if( nextNodeIdx > center)
		{
			if( nextNodeIdx > itsMaxIndex(nextSliceIdx,i) )
            {
                /// Two nodes away => out of lattice upper bound
				GetReconnector()->ReconnectOutOfBoundUp( nextNodeIdx, itsMaxIndex(nextSliceIdx,i), meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }

			else if( nextNodeIdx  == itsMaxIndex(nextSliceIdx,i) )
            {
                /// Central connection on the lattice upper edge
				GetReconnector()->ReconnectOnTheHedgeUp( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }
		}
		
        if(!isReconnected && (probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0))
            /// Last chance : try a binomial reconnection inside the tree
            GetReconnector()->ReconnectInside( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );

        if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
        {
            /// Nothing else to do !!
	        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in the tree" );
        }
        else
        {
            slice->SetNextNodeIndex(stateIdx,i,nextNodeIdx);
            slice->SetProbaUp(stateIdx,i,probaUp);
            slice->SetProbaDown(stateIdx,i,probaDown);
        }
	} // for nbDims
}



////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: toString
///	Returns: 
///	Action : object dump
////////////////////////////////////////////////////
string ARM_TruncatorND::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    size_t nbDims = itsStdDevRatio.size();
	os << indent << dec << setw(1) << nbDims << "D StdDev Truncator : MaxStdDev = (" << fixed;
    for(size_t i=0;i<nbDims-1;++i)
        os << itsStdDevRatio[i] << ",";
    os << itsStdDevRatio[nbDims-1] << ")\n";
	return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: toString
///	Returns: 
///	Action : object dump
////////////////////////////////////////////////////
ARM_Truncator1D* ARM_TruncatorND::CorrespondingTruncator1D(size_t dim) const
{
	ARM_Truncator1D* result = new ARM_Truncator1D( itsStdDevRatio[dim] );
	result->SetReconnector( GetReconnector() );
    return result;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
///	ARM_TruncatorNDArrowDebreu
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_TruncatorNDArrowDebreu
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_TruncatorNDArrowDebreu::ARM_TruncatorNDArrowDebreu( const ARM_TruncatorNDArrowDebreu& rhs )
:   ARM_TruncatorND( rhs ), 
    itsMaxMaxIndex(rhs.itsMaxMaxIndex),
    itsArrowDebreuThreshold(rhs.itsArrowDebreuThreshold),
    itsTimeThreshold(rhs.itsTimeThreshold),
    itsCenter(ARM_IntMatrixPtr(rhs.itsCenter != ARM_IntMatrixPtr(NULL) ? (ARM_IntMatrix*) rhs.itsCenter->Clone() : NULL))
{}


ARM_TruncatorNDArrowDebreu& ARM_TruncatorNDArrowDebreu::operator=(const ARM_TruncatorNDArrowDebreu& rhs )
{
    if( this != &rhs )
    {
        ARM_TruncatorND::operator =( rhs );
        itsMaxMaxIndex=rhs.itsMaxMaxIndex;
        itsArrowDebreuThreshold=rhs.itsArrowDebreuThreshold;
        itsTimeThreshold=rhs.itsTimeThreshold;
        itsCenter = ARM_IntMatrixPtr(rhs.itsCenter != ARM_IntMatrixPtr(NULL) ? (ARM_IntMatrix*) rhs.itsCenter->Clone() : NULL);
    }
    return *this;
}

ARM_TruncatorNDArrowDebreu::~ARM_TruncatorNDArrowDebreu()
{}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorNDArrowDebreu
///	Routine: Init
///	Returns: void
///	Action : initialisation of truncation decision datas
////////////////////////////////////////////////////
void ARM_TruncatorNDArrowDebreu::Init(const ARM_SliceVectorPtr& slices, const ARM_SamplerBase* sampler, const ARM_IntMatrixPtr& center)
{
    /// Compute diffusion center if necessary
    if(center == ARM_IntMatrixPtr(NULL))
        itsCenter = ComputeCenter(sampler,slices);
    else
        itsCenter = center;

    /// Initialise the StdDev truncator to get min and max indexes
    ARM_TruncatorND::Init(slices,sampler,itsCenter);
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorNDArrowDebreu
///	Routine: NeedTruncation
///	Returns: void
///	Action : detect if truncation is needed
////////////////////////////////////////////////////
ARM_BoolVector ARM_TruncatorNDArrowDebreu::NeedTruncation(double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const ARM_TreeIndex& nextNodeIdx) const
{
    size_t nbDims = nextNodeIdx.size();

    bool isADTrunc = sliceTime > itsTimeThreshold ||
                     ( tmpSlice->GetArrowDebreuPrices() != ARM_GP_VectorPtr(NULL) &&
                      tmpSlice->GetArrowDebreuPrices(stateIdx) < itsArrowDebreuThreshold );

    size_t nextSliceIdx = tmpNextSlice->GetIndex();

    int nextMinIdx,nextMaxIdx,center;
    ARM_BoolVector status(nbDims);
	for(size_t i=0; i<nbDims; ++i )
    {
        nextMinIdx = GetMinIndex(nextSliceIdx,i);
        center = (*itsCenter)(nextSliceIdx,i);
        if(isADTrunc && nextMinIdx < center-itsMaxMaxIndex[i])
            nextMinIdx = center-itsMaxMaxIndex[i];
        if(nextNodeIdx[i] <= nextMinIdx)
            status[i] = true;
        
        else
        {
            nextMaxIdx = GetMaxIndex(nextSliceIdx,i);
            if(isADTrunc && nextMaxIdx > center+itsMaxMaxIndex[i])
                nextMaxIdx = center+itsMaxMaxIndex[i];

            status[i] = (nextNodeIdx[i] >= nextMaxIdx);
        }
    }

    return status;
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorNDArrowDebreu
///	Routine: UpdateNode
///	Returns: void
///	Action : truncates if necessary the node
///          Slices are necessary ND slices with explicit
///          nodes (ARM_SliceND)
////////////////////////////////////////////////////
void ARM_TruncatorNDArrowDebreu::UpdateNode( double sliceTime, ARM_SliceBase* tmpSlice, size_t stateIdx, ARM_SliceBase* tmpNextSlice, const std::vector<double>& meanRel, const std::vector<double>& varRel ) const
{
    ARM_SliceND* slice      = tmpSlice->ToSliceND();
    ARM_SliceND* nextSlice  = tmpNextSlice->ToSliceND();

	size_t nbDims		    = slice->dim();
    size_t nextSliceIdx     = nextSlice->GetIndex();

    /// Activation of AD truncation
    bool isADTrunc = sliceTime > itsTimeThreshold ||
                     ( slice->GetArrowDebreuPrices() != ARM_GP_VectorPtr(NULL) &&
                      slice->GetArrowDebreuPrices(stateIdx) < itsArrowDebreuThreshold );

    int nextMinIdx,nextMaxIdx,nextNodeIdx,center;
    double probaUp,probaDown;
    bool isReconnected;
	for(size_t i=0; i<nbDims; ++i )
	{
        nextNodeIdx     = slice->GetNextNodeIndex(stateIdx,i);
        center          = (*itsCenter)(nextSliceIdx,i);
        probaUp         = slice->GetProbaUp(stateIdx,i);
        probaDown       = slice->GetProbaDown(stateIdx,i);
        isReconnected   = false;
        if(nextNodeIdx <= center)
        {
            nextMinIdx = GetMinIndex(nextSliceIdx,i);
            if(isADTrunc && nextMinIdx < center-itsMaxMaxIndex[i])
                nextMinIdx = center-itsMaxMaxIndex[i];

            if( nextNodeIdx < nextMinIdx)
            {
                /// Two nodes away => out of lattice lower bound
                GetReconnector()->ReconnectOutOfBoundDown( nextNodeIdx, nextMinIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }

            else if( nextNodeIdx == nextMinIdx )
            {
                /// Central connection on the lattice lower edge
                GetReconnector()->ReconnectOnTheHedgeDown( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }
        }
        else if(nextNodeIdx > center)
        {
            nextMaxIdx = GetMaxIndex(nextSliceIdx,i);
            if(isADTrunc && nextMaxIdx > center+itsMaxMaxIndex[i])
                nextMaxIdx = center+itsMaxMaxIndex[i];

            if( nextNodeIdx > nextMaxIdx )
            {
                /// Two nodes away => out of lattice upper bound
                GetReconnector()->ReconnectOutOfBoundUp( nextNodeIdx, nextMaxIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }

            else if( nextNodeIdx == nextMaxIdx )
            {
                /// Central connection on the lattice upper edge
                GetReconnector()->ReconnectOnTheHedgeUp( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );
                isReconnected = true;
            }
        }

        if(!isReconnected && (probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0))
            /// Last chance : try a binomial reconnection inside the tree
            GetReconnector()->ReconnectInside( nextNodeIdx, meanRel[i], varRel[i], probaUp, probaDown );

        if(probaUp < 0.0 || probaDown < 0.0 || probaUp + probaDown > 1.0)
        {
            /// Nothing else to do !!
	        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in the tree" );
        }
        else
        {
            slice->SetNextNodeIndex(stateIdx,i,nextNodeIdx);
            slice->SetProbaUp(stateIdx,i,probaUp);
            slice->SetProbaDown(stateIdx,i,probaDown);
        }
    } // for nbDims
}


////////////////////////////////////////////////////
///	Class  : ARM_TruncatorND
///	Routine: toString
///	Returns: 
///	Action : object dump
////////////////////////////////////////////////////
string ARM_TruncatorNDArrowDebreu::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;
    size_t nbDims = GetStdDevRatio().size();
	os << indent << dec << setw(1) << nbDims << "D Arrow-Debreu Truncator :\n";
    os << indent <<  "    MaxStdDev    = (" << fixed;
    for(size_t i=0;i<nbDims-1;++i)
        os << setw(5) << setprecision(2) << GetStdDevRatio(i) << ",";
    os << setw(5) << setprecision(2) << GetStdDevRatio(nbDims-1) << ")\n";
    os << indent <<  "    MaxMaxIndex  = (" << dec;
    for(size_t i=0;i<nbDims-1;++i)
        os << setw(4) << itsMaxMaxIndex[i] << ",";
    os << setw(4) << itsMaxMaxIndex[nbDims-1] << ")\n";
    os << indent <<  "    ADLimit      = " << scientific << itsArrowDebreuThreshold << "\n";
    os << indent <<  "    TimeLimit(y) = " << fixed << setw(5) << setprecision(2) << itsTimeThreshold << "\n";
	return os.str();
}



////////////////////////////////////////////////////
///	Class  : ARM_TruncatorNDArrowDebreu
///	Routine: CorrespondingTruncator1D
///	Returns: ARM_Truncator1D*
///	Action : creates the corresponding truncator 1D
////////////////////////////////////////////////////
ARM_Truncator1D* ARM_TruncatorNDArrowDebreu::CorrespondingTruncator1D(size_t dim) const
{
	ARM_Truncator1D* result = new ARM_Truncator1DArrowDebreu( GetStdDevRatio(dim), itsMaxMaxIndex[dim], itsArrowDebreuThreshold, itsTimeThreshold );
	result->SetReconnector( GetReconnector() );
	return result;
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/