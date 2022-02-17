/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file normalcentredsampler.cpp
 *
 *  \brief sampler object are used to discretise in time according to 
 *      various hypotheses
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/normalcentredsampler.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/numericconstant.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/slice.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_NormalCentredSamplerBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerBase
///	Routine: ARM_SamplerBase
///	Returns: ARM_Sampler1DBase*
///	Action : convert to the 1D equivalent for hybrid model
////////////////////////////////////////////////////

ARM_SamplerBase* ARM_NormalCentredSamplerBase::CorrespondingSampler1D(size_t dim) const
{
	return new ARM_NormalCentredSampler1D( &*GetScheduler());
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_NormalCentredSampler1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSampler1D
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NormalCentredSampler1D::ARM_NormalCentredSampler1D(const ARM_NormalCentredSampler1D& rhs)
: ARM_NormalCentredSamplerBase(rhs),
itsLocalVarX(rhs.itsLocalVarX),
itsGlobalVarX(rhs.itsGlobalVarX)
{
	
}

////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSampler1D
///	Routine: copy constructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_NormalCentredSampler1D::~ARM_NormalCentredSampler1D()
{}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSampler1D
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_NormalCentredSampler1D::Clone() const
{   return new ARM_NormalCentredSampler1D(*this); }



////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSampler1D
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices
///	Action : Init the sampler
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_NormalCentredSampler1D::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t nbSteps = timeSteps->size();

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get VCVs from model according to the schedule
    ARM_MatrixVector localVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,localVCV,globalVCV);
    SetLocalVarX(localVCV);
    SetGlobalVarX(globalVCV);

    if(isInitSlices)
    {
        /// Compute constant space step
        double maxLocalVar = 0.0;
        size_t i;
        for( i=0; i<localVCV.size(); ++i )
        {
            if( (*localVCV[i])(0,0) > maxLocalVar )
                maxLocalVar = (*localVCV[i])(0,0);
        }
        maxLocalVar *= 3.0; // for 4th moment matching

        double spaceStep = sqrt(maxLocalVar);

        /// Build slices
        ARM_SliceVector* slices = new ARM_SliceVector(timeSteps->size());
        double proba;
        for( i=0; i<nbSteps-1; ++i )
        {
            proba    = (*localVCV[i])(0,0) / (2.0*maxLocalVar);

#ifdef __GP_STRICT_VALIDATION
            if(proba < 0.0 || proba > 0.5)
            {
	            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in slice transitions" );
            }
#endif

            (*slices)[i]= new ARM_Slice1DCstSymProba( i,-i,i,spaceStep,proba);
        }

        /// Last slice (proba is set arbitrary to 0)
        i = nbSteps-1;
        (*slices)[i]= new ARM_Slice1DCstSymProba( i,-i,i,spaceStep,0.0);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    /// Keep record of initial global VCV and delete local one
    SetGlobalVCV(globalVCV);
	DeletePointorVector<ARM_GP_Matrix>(localVCV);

    return result;
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_NormalCentredSamplerND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_NormalCentredSamplerND::CopyNoCleanUp(const ARM_NormalCentredSamplerND& rhs)
{
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsZtoX, itsZtoX );
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsXtoZ, itsXtoZ );

    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsLocalVarZ, itsLocalVarZ );
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsGlobalVarZ, itsGlobalVarZ );
}

////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_NormalCentredSamplerND::CleanUp()
{
    DeletePointorVector<ARM_GP_Matrix>( itsZtoX );
    DeletePointorVector<ARM_GP_Matrix>( itsXtoZ );

    DeletePointorVector<ARM_GP_Vector>( itsLocalVarZ );
    DeletePointorVector<ARM_GP_Vector>( itsGlobalVarZ );
}

////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: copy constructor
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_NormalCentredSamplerND::ARM_NormalCentredSamplerND( const ARM_NormalCentredSamplerND& rhs )
:   ARM_SamplerNDBase(rhs), ARM_NormalCentredSamplerBase( rhs ),
    itsZtoX(0), itsXtoZ(0),
    itsLocalVarZ(0), itsGlobalVarZ(0),
	itsNbRank(rhs.itsNbRank)
{
    CopyNoCleanUp(rhs);
}


ARM_NormalCentredSamplerND::~ARM_NormalCentredSamplerND()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_NormalCentredSamplerND::Clone() const
{   return new ARM_NormalCentredSamplerND(*this); }


////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices
///	Action : Init the sampler
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_NormalCentredSamplerND::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t nbSteps = timeSteps->size();

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get VCVs from model according to the schedule
    ARM_MatrixVector localVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,localVCV,globalVCV);

    size_t nbDims = localVCV[0]->rows();

    /// Compute rotation matrixes and variances in Z space
    size_t i,j,k;

    itsZtoX.resize(nbSteps);
    itsXtoZ.resize(nbSteps);

	// The dimension simulated by the process cannot be less than the 
	// rank
	size_t nbEffectiveDims;

	if (itsNbRank == -1)
	{
		nbEffectiveDims = nbDims;
	}
	else
	{
		nbEffectiveDims = (nbDims<itsNbRank?nbDims:itsNbRank);
	}

	itsZtoX[0] = new ARM_GP_Matrix(nbEffectiveDims,nbDims,0.0);
	itsXtoZ[0] = new ARM_GP_Matrix(nbEffectiveDims,nbDims,0.0);

    itsLocalVarZ.resize(nbSteps-1);
    itsGlobalVarZ.resize(nbSteps);

//  ARM_GP_Vector eigenValues(nbDims);
//  ARM_GP_Matrix* eigenMatrix;
//  ARM_GP_Matrix locVCV(nbDims,nbDims),gloVCV(nbDims,nbDims);
    
	ARM_GP_Vector eigenValues;
	ARM_GP_Matrix * eigenMatrix;
	ARM_GP_Matrix locVCV, gloVCV;

	for(i=0;i<nbSteps-1;++i)
    {
		nbDims = localVCV[i]->rows();
		nbEffectiveDims = itsNbRank == -1 ? nbDims : nbDims < itsNbRank ? nbDims : itsNbRank;

		eigenValues.resize(nbDims);
		locVCV.resize(nbDims,nbDims);
		gloVCV.resize(nbDims,nbDims);

        /// Convert triangular to symetric matrix
        for(j=0;j<nbDims;++j)
        {
            locVCV(j,j) = (*(localVCV[i]))(j,j);
            gloVCV(j,j) = (*(globalVCV[i]))(j,j);
            for(k=0;k<j;++k)
            {
                locVCV(j,k)=(*(localVCV[i]))(j,k);
                locVCV(k,j)=locVCV(j,k);
                gloVCV(j,k)=(*(globalVCV[i]))(j,k);
                gloVCV(k,j)=gloVCV(j,k);
            }
        }

        /// Diagonalize VCV matrix
        itsZtoX[i+1] = ACPTransformation(&locVCV,eigenValues,nbEffectiveDims);
		ARM_GP_Matrix* tmpMatrix;
		
		tmpMatrix = static_cast< ARM_GP_Matrix* >(itsZtoX[i+1]->Clone());
		tmpMatrix->transpose();
        itsXtoZ[i+1] = tmpMatrix;
        itsLocalVarZ[i] = new ARM_GP_Vector(eigenValues);

        eigenMatrix = ACPTransformation(&gloVCV,eigenValues);
        itsGlobalVarZ[i] = new ARM_GP_Vector(eigenValues);
        delete eigenMatrix;
    }

    /// Last step for GlobalVar only
    for(j=0;j<nbDims;++j)
    {
        gloVCV(j,j) = (*(globalVCV[nbSteps-1]))(j,j);
        for(k=0;k<j;++k)
        {
            gloVCV(j,k)=(*(globalVCV[nbSteps-1]))(j,k);
            gloVCV(k,j)=gloVCV(j,k);
        }
    }

    eigenMatrix = ACPTransformation(&gloVCV,eigenValues);
    itsGlobalVarZ[nbSteps-1] = new ARM_GP_Vector(eigenValues);
    delete eigenMatrix;

    if(isInitSlices)
    {
        /// Compute constant space steps in Z space
        ARM_GP_Vector maxLocalVarZ(nbDims,0.0);
        for( i=0; i<itsLocalVarZ.size(); ++i )
        {
            for( j=0; j<nbDims; ++j )
                if( (*itsLocalVarZ[i])[j] > maxLocalVarZ[j] )
                    maxLocalVarZ[j] = (*itsLocalVarZ[i])[j];
        }
        maxLocalVarZ *= 3.0; // for 4th moment matching

        ARM_GP_Vector spaceStep(maxLocalVarZ);
        for( j=0; j<nbDims; ++j )
            spaceStep[j] = sqrt(spaceStep[j]);

        /// Build slices
        ARM_SliceVector* slices = new ARM_SliceVector(timeSteps->size());

        /// Last slice (probas are set arbitrary to 0)
        i = nbSteps-1;
        ARM_GP_Vector proba(nbDims,0.0);
        ARM_IntVector min(nbDims,-i);
        ARM_IntVector max(nbDims,i);
		(*slices)[i]= new ARM_SliceNDCstSymProba(nbDims,i,min,max,spaceStep.GetValues(),proba.GetValues());

        for( i=0; i<nbSteps-1; ++i )
        {
            for( j=0; j<nbDims; ++j )
            {
                min[j] = -i;
                max[j] = i;
                proba[j] = (*itsLocalVarZ[i])[j] / (2.0*maxLocalVarZ[j]);

#ifdef __GP_STRICT_VALIDATION
                if(proba[j] < 0.0 || proba[j] > 0.5)
                {
	                ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": negative probabilities in slice transitions" );
                }
#endif

            }
			(*slices)[i] = new ARM_SliceNDCstSymProba(nbDims,i,min,max,spaceStep.GetValues(),proba.GetValues());
        }

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    /// Keep record of initial global VCV and delete local one
    SetGlobalVCV(globalVCV);
	DeletePointorVector<ARM_GP_Matrix>(localVCV);

	itsFakeRelativeDrift = ARM_GP_Matrix(nbSteps, nbDims, 1.0);

	// Construction des dimensions
	itsNbFactors.resize(itsLocalVarZ.size());
	for(k = 0; k < itsLocalVarZ.size(); k++) itsNbFactors[k] = itsLocalVarZ[k]->size();

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_NormalCentredSamplerND
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_VectorPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_NormalCentredSamplerND::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates) const
{
    size_t i;

    size_t stateIdx,nbStates = zStates->cols();
    size_t k;

    /// Get rotation matrix from Z to X
    ARM_GP_Matrix *zToX = itsZtoX[sliceIdx];

	size_t nbDims = zToX->rows();
	size_t nbEffectiveDims = zToX->cols();

	ARM_GP_MatrixPtr xStates( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        for(i=0;i<nbDims;++i)
            for(k=0;k<nbEffectiveDims;++k)
                (*xStates)(i,stateIdx) += (*zToX)(i,k) * (*zStates)(k,stateIdx);

    /// Y = X because no volatility change
    yStates = xStates;

    return xStates;
}


CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

