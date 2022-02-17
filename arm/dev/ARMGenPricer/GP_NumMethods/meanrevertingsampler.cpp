/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file meanrevertingsampler.cpp
 *
 *  \brief sampler object are used to discretise in time according to 
 *      various hypotheses
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpbase/removeidentifiedwarning.h"

#include "gpnummethods/meanrevertingsampler.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/numericconstant.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/slice.h"


CC_BEGIN_NAMESPACE( ARM )

const size_t FIRST_STATE_VARIABLE   = 0;

/// Space step factors to keep transisiton probas between [0,1]
const double KURTOSIS_SPACE_STEP_FACTOR = sqrt(3.0);  // fulfills the 4th normal moment 
const double MAX_SPACE_STEP_FACTOR      = 2.0;
const double MIN_SPACE_STEP_FACTOR      = 2.0/sqrt(3.0);


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MeanRevertingSamplerBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerBase
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MeanRevertingSamplerBase::ARM_MeanRevertingSamplerBase( const ARM_MeanRevertingSamplerBase& rhs )
:   ARM_SamplerBase( rhs ),
    itsMinStdDev( rhs.itsMinStdDev )
{}


ARM_MeanRevertingSamplerBase& ARM_MeanRevertingSamplerBase::operator=(const ARM_MeanRevertingSamplerBase& rhs )
{
    if( this != &rhs )
    {
        ARM_SamplerBase::operator =( rhs );
        itsMinStdDev = rhs.itsMinStdDev;
    }
    return *this;
}

ARM_MeanRevertingSamplerBase::~ARM_MeanRevertingSamplerBase()
{}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerBase
///	Routine: CorrespondingSampler1D
///	Returns: ARM_SamplerBase*
///	Action : convert to the 1D equivalent for hybrid model
////////////////////////////////////////////////////

ARM_SamplerBase* ARM_MeanRevertingSamplerBase::CorrespondingSampler1D(size_t dim) const
{
	return new ARM_MeanRevertingSampler1D( &*GetScheduler(), GetMinStdDev() );
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MeanRevertingSampler1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_MeanRevertingSampler1D::CopyNoCleanUp(const ARM_MeanRevertingSampler1D& rhs)
{
    itsRelDriftX    = rhs.itsRelDriftX;
    itsAbsDriftX    = rhs.itsAbsDriftX;
    itsLocalVarX    = rhs.itsLocalVarX;
    itsGlobalVarX   = rhs.itsGlobalVarX;
}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MeanRevertingSampler1D::ARM_MeanRevertingSampler1D( const ARM_MeanRevertingSampler1D& rhs )
:   ARM_MeanRevertingSamplerBase( rhs ),
    itsRelDriftX(0),itsAbsDriftX(0),itsLocalVarX(0),itsGlobalVarX(0)
{
    CopyNoCleanUp(rhs);
}


ARM_MeanRevertingSampler1D& ARM_MeanRevertingSampler1D::operator=(const ARM_MeanRevertingSampler1D& rhs )
{
    if( this != &rhs )
    {
        ARM_MeanRevertingSamplerBase::operator =( rhs );
        CopyNoCleanUp(rhs);
    }
    return *this;
}

ARM_MeanRevertingSampler1D::~ARM_MeanRevertingSampler1D()
{}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_MeanRevertingSampler1D::Clone() const
{   return new ARM_MeanRevertingSampler1D(*this); }


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Init the sampler
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MeanRevertingSampler1D::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t nbSteps = timeSteps->size();

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get local relative & absolute drifts from model according to the schedule
	ARM_GP_MatrixPtr relativeDrifts;
	ARM_GP_MatrixPtr absoluteDrifts;
	model.IntegratedLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);
	SetRelDriftX(relativeDrifts);
	SetAbsDriftX(absoluteDrifts);

    /// Get VCVs from model according to the schedule
    ARM_MatrixVector localVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,localVCV,globalVCV);
    SetLocalVarX(localVCV);
    SetGlobalVarX(globalVCV);

    if(isInitSlices)
    {
        /// Build slices
        double minVar = GetMinStdDev()*GetMinStdDev();

        double minTimeVar,localVar,spaceStep = 0;
        ARM_SliceVector* slices = new ARM_SliceVector(timeSteps->size());
        (*slices)[0] = new ARM_Slice1D(0,0,0,spaceStep);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        size_t i;
        for( i=1; i<nbSteps; ++i )
        {
            minTimeVar = minVar * ( (*timeSteps)[i]-(*timeSteps)[i-1] ) / K_YEAR_LEN;
            localVar = itsLocalVarX[i-1];
            if(localVar < minTimeVar)
            {
                if(localVar < K_DOUBLE_TOL)
                {
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                        "Can't resize space steps, local variance too low");
                }
                spaceStep = KURTOSIS_SPACE_STEP_FACTOR * sqrt(minTimeVar);
            }
            else
                spaceStep = KURTOSIS_SPACE_STEP_FACTOR * sqrt(localVar);

            (*slices)[i]= new ARM_Slice1D(i,-i,i,spaceStep);
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    /// Keep record of initial global VCV and delete local one
    SetGlobalVCV(globalVCV);
	DeletePointorVector<ARM_GP_Matrix>(localVCV);

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: ApplyRelDriftToIntXStates
///	Returns: void
///	Action : Apply drift on the integrated X 
/// process
////////////////////////////////////////////////////
void ARM_MeanRevertingSampler1D::ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_VectorPtr& XStates) const
{
	double absDrift = (*itsAbsDriftX)(sliceIdx,0);

	int nbStates = XStates->size();

    for(size_t i=0;i<nbStates;++i)
        (*XStates)[i] = absDrift + (*XStates)[i];

}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: ApplyAbsDriftToIntXStates
///	Returns: void
///	Action : Apply drift on X 
/// process
////////////////////////////////////////////////////
void ARM_MeanRevertingSampler1D::ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_VectorPtr& integXStates) const
{
	double relDrift = (*itsRelDriftX)(sliceIdx,0);

	int nbStates = integXStates->size();

    for(size_t i=0;i<nbStates;++i)
        (*integXStates)[i] = relDrift * (*integXStates)[i];

}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSampler1D
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_VectorPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MeanRevertingSampler1D::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const
{
    size_t i,nbStates = zStates->size();
    ARM_GP_Vector* mean = new ARM_GP_Vector(nbStates);

    /// Z = X because 1D
    double relDrift = (*itsRelDriftX)(sliceIdx,0);
    double absDrift;
    if( itsAbsDriftX != ARM_GP_MatrixPtr(NULL) )
    {
        absDrift = (*itsAbsDriftX)(sliceIdx,0);
        for(i=0;i<nbStates;++i)
            (*mean)[i] = relDrift * (*zStates)[i] + absDrift;
    }
    else
    {
        for(i=0;i<nbStates;++i)
            (*mean)[i] = relDrift * (*zStates)[i];
    }

    xStates = ComputeZtoXStates(sliceIdx,zStates);

    return ARM_GP_VectorPtr(mean);
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_DriftedMeanRevertingSampler1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////


////////////////////////////////////////////////////
///	Class  : ARM_DriftedMeanRevertingSampler1D
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_VectorPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_DriftedMeanRevertingSampler1D::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const
{
    size_t i,nbStates = zStates->size();
    ARM_GP_Vector* mean = new ARM_GP_Vector(nbStates);

    /// Get markovian drift part
    xStates = ComputeZtoXStates(sliceIdx,zStates);
    ARM_GP_MatrixPtr markovDriftX = GetModel()->MarkovianDrift( sliceIdx, ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,nbStates,*xStates) ) );

    /// Z = X because 1D
    double relDrift = (*GetRelDriftX())(sliceIdx,FIRST_STATE_VARIABLE);
    double absDrift;
    if( GetAbsDriftX() != ARM_GP_MatrixPtr(NULL) )
    {
        absDrift = (*GetAbsDriftX())(sliceIdx,FIRST_STATE_VARIABLE);
        for(i=0;i<nbStates;++i)
            (*mean)[i] = relDrift * (*zStates)[i] + absDrift + (*markovDriftX)(FIRST_STATE_VARIABLE,i);
    }
    else
    {
        for(i=0;i<nbStates;++i)
            (*mean)[i] = relDrift * (*zStates)[i]  + (*markovDriftX)(FIRST_STATE_VARIABLE,i);
    }

    return ARM_GP_VectorPtr(mean);
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MeanRevertingSamplerND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_MeanRevertingSamplerND::CopyNoCleanUp(const ARM_MeanRevertingSamplerND& rhs)
{
	/// clone to avoid side effect!
	itsRelDriftX = rhs.itsRelDriftX != ARM_GP_MatrixPtr(NULL)? ARM_GP_MatrixPtr( static_cast<ARM_GP_Matrix*>( rhs.itsRelDriftX->Clone() ) ) : ARM_GP_MatrixPtr(NULL);
	itsAbsDriftX = rhs.itsAbsDriftX != ARM_GP_MatrixPtr(NULL)? ARM_GP_MatrixPtr( static_cast<ARM_GP_Matrix*>( rhs.itsAbsDriftX->Clone() ) ) : ARM_GP_MatrixPtr(NULL);

	itsNbRank = rhs.itsNbRank;

    DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsZtoX, itsZtoX );
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsXtoZ, itsXtoZ );

    DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsGlobalXtoZ, itsGlobalXtoZ );

    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsLocalVarZ, itsLocalVarZ );
    DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsGlobalVarZ, itsGlobalVarZ );
}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_MeanRevertingSamplerND::CleanUp()
{
    DeletePointorVector<ARM_GP_Matrix>( itsZtoX );
    DeletePointorVector<ARM_GP_Matrix>( itsXtoZ );

    DeletePointorVector<ARM_GP_Matrix>( itsGlobalXtoZ );

    DeletePointorVector<ARM_GP_Vector>( itsLocalVarZ );
    DeletePointorVector<ARM_GP_Vector>( itsGlobalVarZ );
}

////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MeanRevertingSamplerND::ARM_MeanRevertingSamplerND( const ARM_MeanRevertingSamplerND& rhs )
:   ARM_SamplerNDBase(rhs), ARM_MeanRevertingSamplerBase( rhs ),
    itsRelDriftX(0), itsAbsDriftX(0), itsZtoX(0), itsXtoZ(0), itsGlobalXtoZ(0), itsLocalVarZ(0), itsGlobalVarZ(0)
{
    CopyNoCleanUp(rhs);
}


ARM_MeanRevertingSamplerND& ARM_MeanRevertingSamplerND::operator=(const ARM_MeanRevertingSamplerND& rhs )
{
    if( this != &rhs )
    {
        ARM_MeanRevertingSamplerBase::operator =( rhs );
        CleanUp();
        CopyNoCleanUp(rhs);
    }
    return *this;
}

ARM_MeanRevertingSamplerND::~ARM_MeanRevertingSamplerND()
{
    CleanUp();
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: Clone
///	Returns: ARM_Object*
///	Action : Clones the object
////////////////////////////////////////////////////
ARM_Object* ARM_MeanRevertingSamplerND::Clone() const
{   return new ARM_MeanRevertingSamplerND(*this); }


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Init the sampler
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MeanRevertingSamplerND::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t nbSteps = timeSteps->size();

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get local relative & absolute drifts from model according to the schedule
	ARM_GP_MatrixPtr relativeDrifts;
	ARM_GP_MatrixPtr absoluteDrifts;
	model.IntegratedLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);
	SetRelDriftX(relativeDrifts);
	SetAbsDriftX(absoluteDrifts);

    /// Get VCVs from model according to the schedule
    ARM_MatrixVector localVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,localVCV,globalVCV);

    size_t nbDims = localVCV[0]->rows();
	size_t nbEffectiveDims;

	if (itsNbRank == -1)
	{
		nbEffectiveDims = nbDims;
	}
	else
	{
		nbEffectiveDims = (nbDims<itsNbRank?nbDims:itsNbRank);
	}

    /// Compute rotation matrixes and variances in Z space
    size_t i,j,k;

    if(itsZtoX.size() > 0)
	    DeletePointorVector<ARM_GP_Matrix>(itsZtoX);
    if(itsXtoZ.size() > 0)
	    DeletePointorVector<ARM_GP_Matrix>(itsXtoZ);
    itsZtoX.resize(nbSteps);
    itsXtoZ.resize(nbSteps);
    itsZtoX[0] = new ARM_GP_Matrix(nbEffectiveDims,nbDims,0.0);
    itsXtoZ[0] = new ARM_GP_Matrix(nbEffectiveDims,nbDims,0.0);

    if(itsGlobalXtoZ.size() > 0)
	    DeletePointorVector<ARM_GP_Matrix>(itsGlobalXtoZ);
    itsGlobalXtoZ.resize(nbSteps);

    if(itsLocalVarZ.size() > 0)
	    DeletePointorVector<ARM_GP_Vector>(itsLocalVarZ);
    if(itsGlobalVarZ.size() > 0)
	    DeletePointorVector<ARM_GP_Vector>(itsGlobalVarZ);
    itsLocalVarZ.resize(nbSteps-1);
    itsGlobalVarZ.resize(nbSteps);

    ARM_GP_Vector eigenValues(nbDims);
    ARM_GP_Matrix* eigenMatrix;
    ARM_GP_Matrix locVCV(nbDims,nbDims),gloVCV(nbDims,nbDims);
    bool isGlobalVCV;
    for(i=0;i<nbSteps-1;++i)
    {
        /// Convert triangular to symetric matrix
        isGlobalVCV = globalVCV[i]!=NULL;
        if(isGlobalVCV)
        {
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
        }
        else
        {
            for(j=0;j<nbDims;++j)
            {
                locVCV(j,j) = (*(localVCV[i]))(j,j);
                for(k=0;k<j;++k)
                {
                    locVCV(j,k)=(*(localVCV[i]))(j,k);
                    locVCV(k,j)=locVCV(j,k);
                }
            }
        }

        /// Diagonalize VCV matrix
        /// itsZToX[i+1] is used to convert Z states at time i+1 to X states
        /// It depends on localVCV between time i & i+1
        itsZtoX[i+1] = ACPTransformation(&locVCV,eigenValues,nbEffectiveDims);
        for(j=0;j<nbEffectiveDims;++j)
        {
            if(eigenValues[j] < K_NEW_DOUBLE_TOL)
                eigenValues[j] = K_NEW_DOUBLE_TOL; // numerical noise
        }
		ARM_GP_Matrix* tmpMatrix = static_cast< ARM_GP_Matrix* >(itsZtoX[i+1]->Clone());
		tmpMatrix->transpose();
        itsXtoZ[i+1] = tmpMatrix;
        itsLocalVarZ[i] = static_cast< ARM_GP_Vector* >(eigenValues.Clone());

        if(isGlobalVCV)
        {
            /// Global 
            eigenMatrix = ACPTransformation(&gloVCV,eigenValues);
            for(j=0;j<nbDims;++j)
            {
                if(eigenValues[j] < K_NEW_DOUBLE_TOL)
                    eigenValues[j] = K_NEW_DOUBLE_TOL; // numerical noise
            }
            itsGlobalVarZ[i] = static_cast< ARM_GP_Vector* >(eigenValues.Clone());

            eigenMatrix->transpose();
            itsGlobalXtoZ[i] = eigenMatrix;
        }
    }

    if(globalVCV[nbSteps-1])
    {
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
        for(j=0;j<nbDims;++j)
        {
            if(eigenValues[j] < K_NEW_DOUBLE_TOL)
                eigenValues[j] = K_NEW_DOUBLE_TOL; // numerical noise
        }
        itsGlobalVarZ[nbSteps-1] = static_cast< ARM_GP_Vector* >(eigenValues.Clone());
        eigenMatrix->transpose();
        itsGlobalXtoZ[nbSteps-1] = eigenMatrix;
    }

    if(isInitSlices)
    {
        /// Build slices
        double minVar = GetMinStdDev()*GetMinStdDev();

        double minTimeVar,localVar;
        ARM_GP_Vector spaceStep(nbDims,0.0);
        ARM_IntVector min(nbDims,0);
        ARM_IntVector max(nbDims,0);
        ARM_SliceVector* slices = new ARM_SliceVector(timeSteps->size());
        (*slices)[0]= new ARM_SliceND(nbDims,0,min,max,spaceStep);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        double yfStep,lastYfStep = (*timeSteps)[0] / K_YEAR_LEN;
        for( i=1; i<nbSteps; ++i )
        {
            yfStep = (*timeSteps)[i] / K_YEAR_LEN;
            minTimeVar = minVar * ( yfStep-lastYfStep );
            for( j=0; j<nbDims; ++j )
            {
                min[j]=-i;
                max[j]=i;
                localVar = (*(itsLocalVarZ[i-1]))[j];
                if(localVar < minTimeVar)
                {
                    if(localVar < K_DOUBLE_TOL)
                    {
                        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+
                            "Can't resize space steps, local variance too low");
                    }
                    spaceStep[j] = MAX_SPACE_STEP_FACTOR * sqrt(minTimeVar);
                }
                else
                    spaceStep[j] = KURTOSIS_SPACE_STEP_FACTOR * sqrt(localVar);
            }

            (*slices)[i]= new ARM_SliceND(nbDims,i,min,max,spaceStep);

            lastYfStep = yfStep;
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    /// Keep record of initial global VCV and delete local one
    SetGlobalVCV(globalVCV);
	DeletePointorVector<ARM_GP_Matrix>(localVCV);

	// Construction des dimensions
	itsNbFactors.resize(itsLocalVarZ.size());
	for(k = 0; k < itsLocalVarZ.size(); k++) itsNbFactors[k] = itsLocalVarZ[k]->size();

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MeanRevertingSamplerND::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
{
    /// No use of driftCorrection because no markovian drift

    size_t i,nbDims = zStates->rows();
    size_t stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr zMean( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    /// Rotation to get correlated X variables
    xStates = ComputeZtoXStates(sliceIdx,zStates);

    /// Get conditional mean in X space then rotate to Z one at next slice time
    ARM_GP_Matrix *xToZ = itsXtoZ[sliceIdx+1];
    double xMean,absDrift,relDrift;
    size_t k;
    if( itsAbsDriftX != ARM_GP_MatrixPtr(NULL) )
    {
        /// Use of an absolute deterministic drift on X variables
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        {
            for(i=0;i<nbDims;++i)
            {
                absDrift = (*itsAbsDriftX)(sliceIdx,i);
                relDrift = (*itsRelDriftX)(sliceIdx,i);
                xMean = absDrift + relDrift * (*xStates)(i,stateIdx);
                for(k=0;k<nbDims;++k)
                    (*zMean)(k,stateIdx) += (*xToZ)(k,i) * xMean;
            }
            
        }
    }
    else
    {
        /// No absolute deterministic drift on X variables
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        {
            for(i=0;i<nbDims;++i)
            {
                relDrift = (*itsRelDriftX)(sliceIdx,i);
                xMean = relDrift * (*xStates)(i,stateIdx);
                for(k=0;k<nbDims;++k)
                    (*zMean)(k,stateIdx) += (*xToZ)(k,i) * xMean;
            }
            
        }
    }

    return zMean;
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: ApplyRelDriftToIntXStates
///	Returns: void
///	Action : Apply drift on the integrated X 
/// process
////////////////////////////////////////////////////
void ARM_MeanRevertingSamplerND::ApplyRelDriftToIntXStates(size_t sliceIdx, ARM_GP_MatrixPtr& integXStates) const
{
	size_t i,nbDims = integXStates->rows();
    size_t stateIdx,nbStates = integXStates->cols();
    double relDrift;
 
	if (!itsRelDriftX.IsNull())
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx)
		{
			for(i=0;i<nbDims;++i)
			{
				relDrift = (*itsRelDriftX)(sliceIdx,i);
				(*integXStates)(i,stateIdx) = relDrift * (*integXStates)(i,stateIdx);
			}   
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: ApplyAbsDriftToXStates
///	Returns: void
///	Action : Apply drift on X  states
/// process
////////////////////////////////////////////////////
void ARM_MeanRevertingSamplerND::ApplyAbsDriftToXStates(size_t sliceIdx, ARM_GP_MatrixPtr& XStates) const
{
	size_t i,nbDims = XStates->rows();
    size_t stateIdx,nbStates = XStates->cols();
    double absDrift;

	if (!itsAbsDriftX.IsNull())
	{
		for(stateIdx=0;stateIdx<nbStates;++stateIdx)
		{
			for(i=0;i<nbDims;++i)
			{
				absDrift = (*itsAbsDriftX)(sliceIdx,i);
				(*XStates)(i,stateIdx) = absDrift + (*XStates)(i,stateIdx);
			}   
		}
	}
}




////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_MatrixPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MeanRevertingSamplerND::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates) const
{
    size_t i;
    size_t stateIdx,nbStates = zStates->cols();
    size_t k;

	ARM_GP_Matrix *zToX = itsZtoX[sliceIdx];
	size_t nbDims = zToX->rows();
	size_t nbEffectiveDims = zToX->cols();

	ARM_GP_MatrixPtr xStates( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    /// Get rotation matrix from Z to X
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        for(i=0;i<nbDims;++i)
            for(k=0;k<nbEffectiveDims;++k)
                (*xStates)(i,stateIdx) += (*zToX)(i,k) * (*zStates)(k,stateIdx);

    /// Y = X because no volatility change
    yStates = xStates;

    return xStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_MeanRevertingSamplerND
///	Routine: ComputeXtoZCenters
///	Returns: void
///	Action : Convert centers in X space to centers in Z space
////////////////////////////////////////////////////
void ARM_MeanRevertingSamplerND::ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const
{
    if(globalDrifts == ARM_GP_MatrixPtr(NULL))
        return;

    size_t i,nbDims = globalDrifts->cols();
    size_t sliceIdx,nbSlices = globalDrifts->rows();
    size_t k;

    /// Rotate centers using rotation matrix from X to Z
    ARM_GP_Matrix *xToZ;
    ARM_GP_Vector centers(nbDims);
    for(sliceIdx=0;sliceIdx<nbSlices;++sliceIdx)
    {
        for(i=0;i<nbDims;++i)
            centers[i] = (*globalDrifts)(sliceIdx,i);

        xToZ = itsGlobalXtoZ[sliceIdx];

        for(i=0;i<nbDims;++i)
        {
            (*globalDrifts)(sliceIdx,i) = 0.0;
            for(k=0;k<nbDims;++k)
                (*globalDrifts)(sliceIdx,i) += (*xToZ)(i,k) * centers[k];
        }
    }

}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_DriftedMeanRevertingSamplerND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_DriftedMeanRevertingSamplerND
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_DriftedMeanRevertingSamplerND::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
{
    size_t i,j,nbDims = zStates->rows(),stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr zMean( new ARM_GP_Matrix(nbDims,nbStates,0.0) );

    /// Rotation to get correlated X variables
    xStates = ComputeZtoXStates(sliceIdx,zStates);

    /// Compute markovian drift then restore initial states
    ARM_GP_MatrixPtr markovDriftX = GetModel()->IntegratedMarkovianDrift(sliceIdx,xStates,driftCorrection);

    /// Get conditional mean in X space then rotate to Z one at next slice time
    const ARM_GP_Matrix *xToZ = GetXtoZ()[sliceIdx+1];
    double xMean,absDrift,relDrift;
    if( GetAbsDriftX() != ARM_GP_MatrixPtr(NULL) )
    {
        /// Use of an absolute deterministic drift on X variables
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        {
            for(i=0;i<nbDims;++i)
            {
                absDrift = (*GetAbsDriftX())(sliceIdx,i);
                relDrift = (*GetRelDriftX())(sliceIdx,i);
                xMean = absDrift + relDrift * (*xStates)(i,stateIdx) + (*markovDriftX)(i,stateIdx);
                for(j=0;j<nbDims;++j)
                    (*zMean)(j,stateIdx) += (*xToZ)(j,i) * xMean;
            }
            
        }
    }
    else
    {
        /// No absolute deterministic drift on X variables
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
        {
            for(i=0;i<nbDims;++i)
            {
                relDrift = (*GetRelDriftX())(sliceIdx,i);
                xMean = relDrift * (*xStates)(i,stateIdx) + (*markovDriftX)(i,stateIdx);
                for(j=0;j<nbDims;++j)
                    (*zMean)(j,stateIdx) += (*xToZ)(j,i) * xMean;
            }
            
        }
    }

    return zMean;
}



CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

