/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file markoviandriftsampler.cpp
 *
 *  \brief sampler object are used to discretise in time according to 
 *      various hypotheses
 *
 *	\author  JM Prie, E Benhamou
 *	\version 1.0
 *	\date November 2004
 */


#include "gpnummethods/markoviandriftsampler.h"

/// gpbase
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/numericconstant.h"
#include "gpbase/env.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingstates.h"

/// gpnummethods
#include "gpnummethods/slice.h"

#define MDS_GLOBAL_ROTATION


CC_BEGIN_NAMESPACE( ARM )

const size_t FIRST_STATE_VARIABLE   = 0;
const size_t SECOND_STATE_VARIABLE  = 1;
const size_t THIRD_STATE_VARIABLE   = 2;


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSamplerBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerBase
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MarkovianDriftSamplerBase::ARM_MarkovianDriftSamplerBase( const ARM_MarkovianDriftSamplerBase& rhs )
:   ARM_SamplerBase( rhs ),
	itsMinStdDev( rhs.itsMinStdDev)
{}


ARM_MarkovianDriftSamplerBase& ARM_MarkovianDriftSamplerBase::operator=(const ARM_MarkovianDriftSamplerBase& rhs )
{
    if( this != &rhs )
    {
        ARM_SamplerBase::operator =( rhs );
        itsMinStdDev = rhs.itsMinStdDev;
    }
    return *this;
}

ARM_MarkovianDriftSamplerBase::~ARM_MarkovianDriftSamplerBase()
{}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerBase
///	Routine: CorrespondingSampler1D
///	Returns: ARM_SamplerBase*
///	Action : convert to the 1D equivalent for hybrid model
////////////////////////////////////////////////////

ARM_SamplerBase* ARM_MarkovianDriftSamplerBase::CorrespondingSampler1D(size_t dim) const
{
	return new ARM_MarkovianDriftSampler1D( &*GetScheduler(), GetMinStdDev() );
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSampler1D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_MarkovianDriftSampler1D::CopyNoCleanUp(const ARM_MarkovianDriftSampler1D& rhs)
{
    itsVolX			= rhs.itsVolX == ARM_GP_VectorPtr(NULL)? ARM_GP_VectorPtr(NULL) : ARM_GP_VectorPtr( (ARM_GP_Vector*) rhs.itsVolX->Clone() );
    itsRelDriftZ	= rhs.itsRelDriftZ;
    itsAbsDriftX    = rhs.itsAbsDriftX;
    itsVolZ			= rhs.itsVolZ;
    itsVolZ2		= rhs.itsVolZ2;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine:
///	Returns: 
///	Action : copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MarkovianDriftSampler1D::ARM_MarkovianDriftSampler1D( const ARM_MarkovianDriftSampler1D& rhs )
:   ARM_MarkovianDriftSamplerBase( rhs ), 
    itsVolX(ARM_GP_VectorPtr(NULL)),
    itsRelDriftZ(0),itsAbsDriftX(0)
{
    CopyNoCleanUp(rhs);
}


ARM_MarkovianDriftSampler1D& ARM_MarkovianDriftSampler1D::operator=(const ARM_MarkovianDriftSampler1D& rhs )
{
    if( this != &rhs )
    {
        ARM_MarkovianDriftSamplerBase::operator =( rhs );
        CopyNoCleanUp(rhs);
    }
    return *this;
}

ARM_MarkovianDriftSampler1D::~ARM_MarkovianDriftSampler1D()
{}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Initialise internal structures :
///          - volatilities for variable change
///          - relative drift contributions
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MarkovianDriftSampler1D::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t i,nbSteps=timeSteps->size();
    ARM_GP_Vector yfSteps(*timeSteps);
    yfSteps /= K_YEAR_LEN;

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get instantaneous volatilities and derivates from model
    /// sampled by the schedule
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls);

    itsVolX = ARM_GP_VectorPtr( vols->GetRow(FIRST_STATE_VARIABLE) );

/***
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"1D MarkovianDrift Sampled Volatilities\n");
fprintf(f,"--------------------------------------\n");
fprintf(f,"    \t       t\t                vol\t                d1Vol\n");
for(i=0;i<timeSteps->size();++i)
    fprintf(f,"#%3d\t  %6.2lf\t    %15.10lf\t      %15.10lf\n",i,(*timeSteps)[i],(*vols)(0,i),(*d1Vols)(0,i));
fclose(f);
***/

    /// Set constant volatility for variable change
    itsVolZ    = (*itsVolX)[vols->cols()-1];
    if(itsVolZ < GetMinStdDev())
        itsVolZ = GetMinStdDev();

    itsVolZ2   = itsVolZ*itsVolZ;

    /// Get instantaneous drifts (deterministic part)
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);

    /// Compute drift contributions of Y to Z and absolute drift
    itsRelDriftZ.resize(nbSteps-1);
    itsAbsDriftX.resize(nbSteps-1);
    double dt;
    for(i=0;i<nbSteps-1;++i)
    {
        dt = yfSteps[i+1]-yfSteps[i];

        itsRelDriftZ[i] = (*relativeDrifts)(i,FIRST_STATE_VARIABLE)
            -dt * (*d1Vols)(FIRST_STATE_VARIABLE,i) / (*vols)(FIRST_STATE_VARIABLE,i);

        itsAbsDriftX[i] = (*absoluteDrifts)(i,FIRST_STATE_VARIABLE);
    }

    if(isInitSlices)
    {
        double spaceStep = 0.0;
        ARM_SliceVector* slices = new ARM_SliceVector(nbSteps);
        (*slices)[0] = new ARM_Slice1D(0,0,0,spaceStep);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );


        /// Slice space step at ti is computed from standard deviation
        /// of Y process from ti-1 to ti
        for(i=1;i<nbSteps;++i)
        {
            dt = yfSteps[i]-yfSteps[i-1];
            spaceStep       = ARM_NumericConstants::ARM_SQRT_3 * itsVolZ * sqrt(dt);
            (*slices)[i]    = new ARM_Slice1D(i,-i,i,spaceStep);
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }
	
    return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_VectorPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MarkovianDriftSampler1D::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_VectorPtr& zStates, ARM_GP_VectorPtr& xStates) const
{
    /// Get markovian drift part
    xStates = ComputeZtoXStates(sliceIdx,zStates);
    ARM_GP_MatrixPtr markovDriftX = GetModel()->MarkovianDrift(sliceIdx,ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,xStates->size(),*xStates) ));

    size_t i,nbStates = zStates->size();
    ARM_GP_Vector* zMean = new ARM_GP_Vector(nbStates);

    /// no rotation because 1D
    double volRatio = itsVolZ/(*itsVolX)[sliceIdx];
    for(i=0;i<nbStates;++i)
        (*zMean)[i] = (1+itsRelDriftZ[sliceIdx]) * (*zStates)[i] + 
            + ( itsAbsDriftX[sliceIdx] + (*markovDriftX)(FIRST_STATE_VARIABLE,i) ) * volRatio;

    return ARM_GP_VectorPtr(zMean);
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: GetLocalVar
///	Returns: double
///	Action : Get local variance between current & next slice
////////////////////////////////////////////////////
double ARM_MarkovianDriftSampler1D::GetLocalVar(size_t sliceIdx) const
{
    return itsVolZ2 * (GetModel()->GetNumMethod()->GetTimeStep(sliceIdx+1)-GetModel()->GetNumMethod()->GetTimeStep(sliceIdx))/K_YEAR_LEN;
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: GetGlobalVar
///	Returns: double
///	Action : Get global variance up to current slice
////////////////////////////////////////////////////
double ARM_MarkovianDriftSampler1D::GetGlobalVar(size_t sliceIdx) const
{
    return itsVolZ2 * GetModel()->GetNumMethod()->GetTimeStep(sliceIdx) / K_YEAR_LEN;
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_VectorPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MarkovianDriftSampler1D::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_VectorPtr& zStates) const
{
    size_t nbStates = zStates->size();
    ARM_GP_Vector* xStates = new ARM_GP_Vector(nbStates,0.0);

    /// Back from constant volatility variable Yto non constant volatility X variables
    double volRatio = (*itsVolX)[sliceIdx] / itsVolZ;
    for( size_t i=0; i<nbStates; ++i )
        (*xStates)[i] = (*zStates)[i] * volRatio;

// FIXMEFRED: mig.vc8 (30/05/2007 16:39:59):cast
	return static_cast<ARM_GP_VectorPtr>(xStates); 
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler1D
///	Routine: ComputeXtoZCenters
///	Returns: void
///	Action : Convert centers in X space to centers in Z space
////////////////////////////////////////////////////
void ARM_MarkovianDriftSampler1D::ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const
{
    if(globalDrifts == ARM_GP_MatrixPtr(NULL) )
        return;

    size_t sliceIdx,nbSlices = globalDrifts->rows();
    double volRatio;

    /// Change X to Y due to volatility change then return because Z=Y
    for(sliceIdx=0;sliceIdx<nbSlices;++sliceIdx)
    {
        volRatio = itsVolZ / (*itsVolX)[sliceIdx];
        (*globalDrifts)(sliceIdx,0) *= volRatio;
    }
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSamplerNDBase
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerNDBase
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_MarkovianDriftSamplerNDBase::CopyNoCleanUp(const ARM_MarkovianDriftSamplerNDBase& rhs)
{
    itsVolX = (rhs.itsVolX == ARM_GP_MatrixPtr(NULL)
                    ? ARM_GP_MatrixPtr(NULL)
                    : ARM_GP_MatrixPtr( static_cast< ARM_GP_Matrix* >(rhs.itsVolX->Clone()) ));
    itsVolY = rhs.itsVolY;

    DuplicateCloneablePtrAndNullVectorInPlace<ARM_GP_Matrix>(rhs.itsYtoZ,itsYtoZ);
    DuplicateCloneablePtrAndNullVectorInPlace<ARM_GP_Matrix>(rhs.itsZtoY,itsZtoY);

    DuplicateCloneablePtrAndNullVectorInPlace<ARM_GP_Matrix>(rhs.itsGlobalYtoZ,itsGlobalYtoZ);

    DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Vector>( rhs.itsLocalVarZ, itsLocalVarZ );
    DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Vector>( rhs.itsGlobalVarZ, itsGlobalVarZ );

    DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Vector>(rhs.itsRelDriftY,itsRelDriftY);

    DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Vector>(rhs.itsAbsDriftX,itsAbsDriftX);
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerNDBase
///	Routine: CleanUp
///	Returns: 
///	Action : Arguments destruction
////////////////////////////////////////////////////
void ARM_MarkovianDriftSamplerNDBase::CleanUp()
{
    DeletePointorVector<ARM_GP_Vector>( itsLocalVarZ );
    DeletePointorVector<ARM_GP_Vector>( itsGlobalVarZ );

    DeletePointorVector<ARM_GP_Vector>( itsRelDriftY );

    DeletePointorVector<ARM_GP_Vector>( itsAbsDriftX );
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerNDBase
///	Routine:
///	Returns: 
///	Action : constructor, copy constructor, operator=, destructor
////////////////////////////////////////////////////
ARM_MarkovianDriftSamplerNDBase::ARM_MarkovianDriftSamplerNDBase( const ARM_SchedulerBase* scheduler,double minStdDev )
:   ARM_MarkovianDriftSamplerBase(scheduler,minStdDev),
    itsVolX(ARM_GP_MatrixPtr(NULL)),
    itsVolY(0),itsYtoZ(0),itsZtoY(0),itsGlobalYtoZ(0),
    itsLocalVarZ(0),itsGlobalVarZ(0),
    itsRelDriftY(0),itsAbsDriftX(0)
{}

ARM_MarkovianDriftSamplerNDBase::ARM_MarkovianDriftSamplerNDBase( const ARM_MarkovianDriftSamplerNDBase& rhs )
:   ARM_MarkovianDriftSamplerBase(rhs),
    itsVolX(ARM_GP_MatrixPtr(NULL)),
    itsVolY(0),itsYtoZ(0),itsZtoY(0),itsGlobalYtoZ(0),
    itsLocalVarZ(0),itsGlobalVarZ(0),
    itsRelDriftY(0),itsAbsDriftX(0)
{
    CopyNoCleanUp(rhs);
}


ARM_MarkovianDriftSamplerNDBase& ARM_MarkovianDriftSamplerNDBase::operator=(const ARM_MarkovianDriftSamplerNDBase& rhs )
{
    if( this != &rhs )
    {
        ARM_MarkovianDriftSamplerBase::operator =( rhs );
        CleanUp();
        CopyNoCleanUp(rhs);
    }
    return *this;
}

ARM_MarkovianDriftSamplerNDBase::~ARM_MarkovianDriftSamplerNDBase()
{
    CleanUp();
}



////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSampler2D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler2D
///	Routine:
///	Returns: 
///	Action : operator=
////////////////////////////////////////////////////
ARM_MarkovianDriftSampler2D& ARM_MarkovianDriftSampler2D::operator=(const ARM_MarkovianDriftSampler2D& rhs )
{
    if( this != &rhs )
    {
        ARM_MarkovianDriftSamplerNDBase::operator =( rhs );
    }
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler2D
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Initialise internal structures :
///          - volatilities for variable change
///          - PCA
///          - relative drift contributions
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MarkovianDriftSampler2D::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t i,j,k,nbSteps=timeSteps->size();
    ARM_GP_Vector yfSteps(*timeSteps);
    yfSteps /= K_YEAR_LEN;

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get instantaneous volatilities and derivates from model
    /// sampled by the schedule
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls);

    size_t nbDims = vols->rows();
    if(nbDims != 2)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : two factors model is required to use MD 2D sampler");

    itsVolX = vols;

#ifdef MDS_GLOBAL_ROTATION
    /// Get VCVs from model according to the schedule
    ARM_MatrixVector unsedLocalVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,unsedLocalVCV,globalVCV);
    ARM_GP_Matrix gloVCV(nbDims,nbDims);
#endif

    /// Set constant volatility for variable change
    itsVolY.resize(nbDims);
    double minVolY = GetMinStdDev();
    for(i=0;i<nbDims;++i)
    {
        itsVolY[i] = (*itsVolX)(i,vols->cols()-1);
        if(itsVolY[i] < minVolY)
            itsVolY[i] = minVolY;
    }

    /// Compute rotation matrixes Y (correlated cst vol) to Z (independant) spaces
    double rho,lastRho = -2.0;
    ARM_GP_MatrixPtr lastZtoY,lastYtoZ,lastGlobalYtoZ;
    ARM_GP_Matrix vcv(2,2);
    ARM_GP_Vector eigenValues(2),globalEigenValues(2);
    ARM_GP_Vector* lastEigenValues;
    ARM_GP_Matrix* eigenVectors;
    ARM_GP_Matrix* globalEigenVectors;

    itsZtoY.resize(nbSteps);
    itsYtoZ.resize(nbSteps);
    itsZtoY[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,2,0.0) );
    itsYtoZ[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,2,0.0) );

    itsGlobalYtoZ.resize(nbSteps);

    ARM_GP_Vector volRatio(2);

    itsLocalVarZ.resize(nbSteps-1);
    itsGlobalVarZ.resize(nbSteps);
    double dt;
    for(i=0;i<nbSteps;++i)
    {
        if( lastRho != (rho = (*correls)(0,i)) )
        {
            vcv(0,0) = itsVolY[0] * itsVolY[0];
            vcv(0,1) = itsVolY[0] * itsVolY[1] * rho;
            vcv(1,0) = vcv(0,1);
            vcv(1,1) = itsVolY[1] * itsVolY[1];

            eigenVectors = ACPTransformation(&vcv,eigenValues);

            if(i<nbSteps-1)
            {
                /// itsZToY[i+1] is used to convert Z states at slice i+1 to Y states
                /// It depends on localVCV = ((ti+1-ti)*vols(ti)) between slice ti & ti+1
                itsZtoY[i+1]    = ARM_GP_MatrixPtr( eigenVectors );
                lastZtoY        = itsZtoY[i+1];

		        ARM_GP_Matrix* tmpMatrix = static_cast< ARM_GP_Matrix* >(eigenVectors->Clone());
		        tmpMatrix->transpose();
                itsYtoZ[i+1]    = ARM_GP_MatrixPtr( tmpMatrix );
                lastYtoZ        = itsYtoZ[i+1];
            }

            lastEigenValues = &eigenValues;
            lastRho = rho;
        }
        else if(i<nbSteps-1)
        {
            /// Keep previous rotations
            itsZtoY[i+1]  = lastZtoY;
            itsYtoZ[i+1]  = lastYtoZ;
        }

        /// Local variances in Z space
        if(i<nbSteps-1)
        {
            itsLocalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
            dt = yfSteps[i+1] - yfSteps[i];
            (*(itsLocalVarZ[i]))[0] *= dt;
            (*(itsLocalVarZ[i]))[1] *= dt;
        }

#ifdef MDS_GLOBAL_ROTATION

        /// Correct computation using the diagonalized global VCV on Y processes
        for(j=0;j<nbDims;++j)
        {
            volRatio[j] = itsVolY[j] / (*itsVolX)(j,i);
            gloVCV(j,j) = volRatio[j] * volRatio[j] * (*(globalVCV[i]))(j,j);
            for(k=0;k<j;++k)
            {
                gloVCV(j,k)= volRatio[j] * volRatio[k] * (*(globalVCV[i]))(j,k);
                gloVCV(k,j)=gloVCV(j,k);
            }
        }
        globalEigenVectors = ACPTransformation(&gloVCV,globalEigenValues);
        itsGlobalVarZ[i] = static_cast< ARM_GP_Vector* >(globalEigenValues.Clone());

        globalEigenVectors->transpose();
        itsGlobalYtoZ[i] = ARM_GP_MatrixPtr( globalEigenVectors );

#else

        /// Global variances in Z space (not correct because correlation, mean reversion...
        /// but something near Dimitri's code)
        itsGlobalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
        (*(itsGlobalVarZ[i]))[0] *= yfSteps[i];
        (*(itsGlobalVarZ[i]))[1] *= yfSteps[i];

        itsGlobalYtoZ[i] = lastYtoZ;

#endif

    }

#ifdef MDS_GLOBAL_ROTATION
    /// Free temporary VCV matrixes returned by the model
    DeletePointorVector<ARM_GP_Matrix>(unsedLocalVCV);
    DeletePointorVector<ARM_GP_Matrix>(globalVCV);
#endif

    /// Get instantaneous drifts (deterministic part)
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);

    /// Compute absolute drift of X and relative drift of Y
    itsAbsDriftX.resize(nbSteps-1);
    itsRelDriftY.resize(nbSteps-1);
    for(i=0;i<nbSteps-1;++i)
    {
        dt = yfSteps[i+1] - yfSteps[i];

        itsAbsDriftX[i] = new ARM_GP_Vector(2);
		itsRelDriftY[i] = new ARM_GP_Vector(2);

        (*(itsRelDriftY[i]))[0] = (*relativeDrifts)(i,0) - dt*(*d1Vols)(0,i)/(*vols)(0,i);
        (*(itsRelDriftY[i]))[1] = (*relativeDrifts)(i,1) - dt*(*d1Vols)(1,i)/(*vols)(1,i);

        (*(itsAbsDriftX[i]))[0] = (*absoluteDrifts)(i,0);
        (*(itsAbsDriftX[i]))[1] = (*absoluteDrifts)(i,1);
    }

    if(isInitSlices)
    {
        ARM_GP_Vector spaceSteps(2,0.0);
        ARM_SliceVector* slices = new ARM_SliceVector(nbSteps);
        (*slices)[0] = new ARM_SliceND(2,0,ARM_IntVector(2,0),ARM_IntVector(2,0),spaceSteps);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        for(i=1;i<nbSteps;++i)
        {
            spaceSteps[0] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[0]);
            spaceSteps[1] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[1]);

            (*slices)[i]    = new ARM_SliceND(2,i,ARM_IntVector(2,-i),ARM_IntVector(2,i),spaceSteps);
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler2D
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSampler2D::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
{
    size_t stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr zMean( new ARM_GP_Matrix(2,nbStates,0.0) );

    /// Get original X variable states to compute absolute markovian drift from model
    ARM_GP_MatrixPtr yStates;
    xStates = ComputeZtoXStates(sliceIdx,zStates,yStates);

    /// Add drift corrections, compute markovian drift then restore initial states
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*xStates)(FIRST_STATE_VARIABLE,stateIdx)   += (*driftCorrection)[FIRST_STATE_VARIABLE];
        (*xStates)(SECOND_STATE_VARIABLE,stateIdx)  += (*driftCorrection)[SECOND_STATE_VARIABLE];
    }
    ARM_GP_MatrixPtr markovDriftX = GetModel()->MarkovianDrift(sliceIdx,xStates);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*xStates)(FIRST_STATE_VARIABLE,stateIdx)   -= (*driftCorrection)[FIRST_STATE_VARIABLE];
        (*xStates)(SECOND_STATE_VARIABLE,stateIdx)  -= (*driftCorrection)[SECOND_STATE_VARIABLE];
    }

    /// Get relative drift of Y
    ARM_GP_Vector *relDriftY = itsRelDriftY[sliceIdx];

    /// Get absolute drift of X
    ARM_GP_Vector* absDriftX = itsAbsDriftX[sliceIdx];

    /// Compute vol ratio for conversion from X to Y
    ARM_GP_Vector volRatio(2);
    volRatio[FIRST_STATE_VARIABLE] = itsVolY[FIRST_STATE_VARIABLE] / (*itsVolX)(FIRST_STATE_VARIABLE,sliceIdx);
    volRatio[SECOND_STATE_VARIABLE] = itsVolY[SECOND_STATE_VARIABLE] / (*itsVolX)(SECOND_STATE_VARIABLE,sliceIdx);

    /// Get rotation from Y to Z at next slice time
    ARM_GP_MatrixPtr yToZ = itsYtoZ[sliceIdx+1];

    /// Compute conditionnal mean
    ARM_GP_Vector driftY(2);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        /// Absolute part
        driftY[FIRST_STATE_VARIABLE] = volRatio[FIRST_STATE_VARIABLE] * 
            ( (*absDriftX)[FIRST_STATE_VARIABLE] + (*markovDriftX)(FIRST_STATE_VARIABLE,stateIdx) );
        driftY[SECOND_STATE_VARIABLE] = volRatio[SECOND_STATE_VARIABLE] * 
            ( (*absDriftX)[SECOND_STATE_VARIABLE] + (*markovDriftX)(SECOND_STATE_VARIABLE,stateIdx) );

        /// Relative part
        driftY[FIRST_STATE_VARIABLE] +=
            (*relDriftY)[FIRST_STATE_VARIABLE] * (*yStates)(FIRST_STATE_VARIABLE,stateIdx);
        driftY[SECOND_STATE_VARIABLE] +=
            (*relDriftY)[SECOND_STATE_VARIABLE] * (*yStates)(SECOND_STATE_VARIABLE,stateIdx);

        /// Contributions in Z space
        (*zMean)(FIRST_STATE_VARIABLE,stateIdx) = (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*yToZ)(0,0) * driftY[FIRST_STATE_VARIABLE] +
            (*yToZ)(0,1) * driftY[SECOND_STATE_VARIABLE];

        (*zMean)(SECOND_STATE_VARIABLE,stateIdx) = (*zStates)(SECOND_STATE_VARIABLE,stateIdx) +
            (*yToZ)(1,0) * driftY[FIRST_STATE_VARIABLE] +
            (*yToZ)(1,1) * driftY[SECOND_STATE_VARIABLE];
    }

    return zMean;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler2D
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_VectorPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSampler2D::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates) const
{
    size_t stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr xStates( new ARM_GP_Matrix(2,nbStates) );
    yStates = ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,nbStates) );

    /// Get rotation matrix from Z to Y
    ARM_GP_MatrixPtr zToY = itsZtoY[sliceIdx];

    /// Compute vol ratio for conversion from Y to X
    ARM_GP_Vector volRatio(2);
    volRatio[FIRST_STATE_VARIABLE] = (*itsVolX)(FIRST_STATE_VARIABLE,sliceIdx) / itsVolY[FIRST_STATE_VARIABLE];
    volRatio[SECOND_STATE_VARIABLE] = (*itsVolX)(SECOND_STATE_VARIABLE,sliceIdx) / itsVolY[SECOND_STATE_VARIABLE];
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*yStates)(FIRST_STATE_VARIABLE,stateIdx) =
            (*zToY)(0,0) * (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*zToY)(0,1) * (*zStates)(SECOND_STATE_VARIABLE,stateIdx);

        (*xStates)(FIRST_STATE_VARIABLE,stateIdx) = volRatio[FIRST_STATE_VARIABLE] * 
            (*yStates)(FIRST_STATE_VARIABLE,stateIdx);

        (*yStates)(SECOND_STATE_VARIABLE,stateIdx) = 
            (*zToY)(1,0) * (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*zToY)(1,1) * (*zStates)(SECOND_STATE_VARIABLE,stateIdx);

        (*xStates)(SECOND_STATE_VARIABLE,stateIdx) = volRatio[SECOND_STATE_VARIABLE] *
            (*yStates)(SECOND_STATE_VARIABLE,stateIdx);
    }

    return xStates; 
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler2D
///	Routine: ComputeXtoZCenters
///	Returns: void
///	Action : Convert centers in X space to centers in Z space
////////////////////////////////////////////////////
void ARM_MarkovianDriftSampler2D::ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const
{
    if(globalDrifts == ARM_GP_MatrixPtr(NULL) )
        return;

    size_t sliceIdx,nbSlices = globalDrifts->rows();
    ARM_GP_MatrixPtr yToZ;
    double cy0,cy1;

    /// Change X to Y due to volatility change then rotate Y to Z
    for(sliceIdx=0;sliceIdx<nbSlices;++sliceIdx)
    {
        cy0 = itsVolY[0] / (*itsVolX)(0,sliceIdx) * (*globalDrifts)(sliceIdx,0);
        cy1 = itsVolY[1] / (*itsVolX)(1,sliceIdx) * (*globalDrifts)(sliceIdx,1);
        yToZ = itsGlobalYtoZ[sliceIdx];
        (*globalDrifts)(sliceIdx,0) = (*yToZ)(0,0) * cy0 + (*yToZ)(0,1) * cy1;
        (*globalDrifts)(sliceIdx,1) = (*yToZ)(1,0) * cy0 + (*yToZ)(1,1) * cy1;
    }
}

////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSampler3D
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler3D
///	Routine:
///	Returns: 
///	Action : coperator=
////////////////////////////////////////////////////
ARM_MarkovianDriftSampler3D& ARM_MarkovianDriftSampler3D::operator=(const ARM_MarkovianDriftSampler3D& rhs )
{
    if( this != &rhs )
    {
        ARM_MarkovianDriftSamplerNDBase::operator =( rhs );
    }
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler3D
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Initialise internal structures :
///          - volatilities for variable change
///          - PCA
///          - relative drift contributions
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MarkovianDriftSampler3D::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t i,j,k,nbSteps=timeSteps->size();
    ARM_GP_Vector yfSteps(*timeSteps);
    yfSteps /= K_YEAR_LEN;

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

    /// Get instantaneous volatilities and derivates from model
    /// sampled by the schedule
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls);

    size_t nbDims = vols->rows();
    if(nbDims != 3)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : three factors model is required to use MD 3D sampler");

/****
FILE* f=fopen("c:\\temp\\dumpTree2G.txt","a");
fprintf(f,"3D MarkovianDrift Sampled Volatilities\n");
fprintf(f,"--------------------------------------\n");
for( j=0; j<nbSteps; ++j )
{
    fprintf(f,"#%3d\tt=\t%7.1lf\t",j,(*timeSteps)[j]);
    fprintf(f,"domVol=\t%9.7lf\t%9.7lf\t",(*vols)(0,j),(*d1Vols)(0,j));
    fprintf(f,"forVol=\t%9.7lf\t%9.7lf\t",(*vols)(1,j),(*d1Vols)(1,j));
    fprintf(f,"fxVol =\t%9.7lf\t%9.7lf\t",(*vols)(2,j),(*d1Vols)(2,j));
    fprintf(f,"\n");
}
fprintf(f,"\n");
fclose(f);
****/

    itsVolX = vols;

#ifdef MDS_GLOBAL_ROTATION
    /// Get VCVs from model according to the schedule
    ARM_MatrixVector unsedLocalVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,unsedLocalVCV,globalVCV);
    ARM_GP_Matrix gloVCV(nbDims,nbDims);
#endif

    /// Set constant volatility for variable change
    itsVolY.resize(3);
    double minVolY = GetMinStdDev();
    for(i=0;i<nbDims;++i)
    {
        itsVolY[i] = (*itsVolX)(i,vols->cols()-1);
        if(itsVolY[i] < minVolY)
            itsVolY[i] = minVolY;
    }

    /// Compute rotation matrixes Y (correlated cst vol) to Z (independant) spaces
    double rho12,lastRho12 = -2.0;
    double rho13,lastRho13 = -2.0;
    double rho23,lastRho23 = -2.0;
    ARM_GP_MatrixPtr lastZtoY,lastYtoZ;
    ARM_GP_Matrix vcv(3,3);
    ARM_GP_Vector eigenValues(3),globalEigenValues(3);
    ARM_GP_Vector* lastEigenValues;
    ARM_GP_Matrix* eigenVectors;
    ARM_GP_Matrix* globalEigenVectors;

    itsZtoY.resize(nbSteps);
    itsYtoZ.resize(nbSteps);
    itsZtoY[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(3,3,0.0) );
    itsYtoZ[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(3,3,0.0) );

    itsGlobalYtoZ.resize(nbSteps);

    ARM_GP_Vector volRatio(3);

    itsLocalVarZ.resize(nbSteps-1);
    itsGlobalVarZ.resize(nbSteps);
    double dt;
    for(i=0;i<nbSteps;++i)
    {
        rho12 = (*correls)(0,i);
        rho13 = (*correls)(1,i);
        rho23 = (*correls)(2,i);
        if( lastRho12 != rho12 || lastRho13 != rho13 || lastRho23 != rho23 )
        {
            vcv(0,0) = itsVolY[0] * itsVolY[0];
            vcv(0,1) = itsVolY[0] * itsVolY[1] * rho12;
            vcv(0,2) = itsVolY[0] * itsVolY[2] * rho13;
            vcv(1,0) = vcv(0,1);
            vcv(2,0) = vcv(0,2);
            vcv(1,1) = itsVolY[1] * itsVolY[1];
            vcv(1,2) = itsVolY[1] * itsVolY[2] * rho23;
            vcv(2,1) = vcv(1,2);
            vcv(2,2) = itsVolY[2] * itsVolY[2];

            eigenVectors = ACPTransformation(&vcv,eigenValues);

            if(i<nbSteps-1)
            {
                /// itsZToY[i+1] is used to convert Z states at slice i+1 to Y states
                /// It depends on localVCV = ((ti+1-ti)*vols(ti)) between slice ti & ti+1
                itsZtoY[i+1]  = ARM_GP_MatrixPtr( eigenVectors );
                lastZtoY    = itsZtoY[i+1];

		        ARM_GP_Matrix* tmpMatrix = static_cast< ARM_GP_Matrix* >(eigenVectors->Clone());
		        tmpMatrix->transpose();
                itsYtoZ[i+1]  = ARM_GP_MatrixPtr( tmpMatrix );
                lastYtoZ    = itsYtoZ[i+1];
            }

            lastEigenValues = &eigenValues;
            lastRho12 = rho12;
            lastRho13 = rho13;
            lastRho23 = rho23;
        }
        else if(i<nbSteps-1)
        {
            /// Keep previous rotations
            itsZtoY[i+1]  = lastZtoY;
            itsYtoZ[i+1]  = lastYtoZ;
        }

        /// Local & global variances in Z space
        if(i<nbSteps-1)
        {
            itsLocalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
            dt = yfSteps[i+1] - yfSteps[i];
            (*(itsLocalVarZ[i]))[0] *= dt;
            (*(itsLocalVarZ[i]))[1] *= dt;
            (*(itsLocalVarZ[i]))[2] *= dt;
        }

#ifdef MDS_GLOBAL_ROTATION

        /// Correct computation using diagonalized global VCV
        for(j=0;j<nbDims;++j)
        {
            volRatio[j] = itsVolY[j] / (*itsVolX)(j,i);
            gloVCV(j,j) = volRatio[j] * volRatio[j] * (*(globalVCV[i]))(j,j);
            for(k=0;k<j;++k)
            {
                gloVCV(j,k)= volRatio[j] * volRatio[k] * (*(globalVCV[i]))(j,k);
                gloVCV(k,j)=gloVCV(j,k);
            }
        }
        globalEigenVectors = ACPTransformation(&gloVCV,globalEigenValues);
        itsGlobalVarZ[i] = static_cast< ARM_GP_Vector* >(globalEigenValues.Clone());

        globalEigenVectors->transpose();
        itsGlobalYtoZ[i] = ARM_GP_MatrixPtr( globalEigenVectors );

#else

        /// Global variances in Z space (not correct because correlation, mean reversion...
        /// but something near Dimitri's code)
        itsGlobalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
        (*(itsGlobalVarZ[i]))[0] *= yfSteps[i];
        (*(itsGlobalVarZ[i]))[1] *= yfSteps[i];
        (*(itsGlobalVarZ[i]))[2] *= yfSteps[i];

        itsGlobalYtoZ[i] = lastYtoZ;

#endif

    }

#ifdef MDS_GLOBAL_ROTATION
    /// Free temporary VCV matrixes returned by the model
    DeletePointorVector<ARM_GP_Matrix>(unsedLocalVCV);
    DeletePointorVector<ARM_GP_Matrix>(globalVCV);
#endif

    /// Get instantaneous drifts (deterministic part)
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);

    /// Compute absolute drift of X and relative drift of Y
    itsAbsDriftX.resize(nbSteps-1);
    itsRelDriftY.resize(nbSteps-1);
    for(i=0;i<nbSteps-1;++i)
    {
        dt = yfSteps[i+1] - yfSteps[i];

        itsAbsDriftX[i] = new ARM_GP_Vector(3);
		itsRelDriftY[i] = new ARM_GP_Vector(3);

        (*(itsRelDriftY[i]))[0] = (*relativeDrifts)(i,0) - dt*(*d1Vols)(0,i)/(*vols)(0,i);
        (*(itsRelDriftY[i]))[1] = (*relativeDrifts)(i,1) - dt*(*d1Vols)(1,i)/(*vols)(1,i);
        (*(itsRelDriftY[i]))[2] = (*relativeDrifts)(i,2) - dt*(*d1Vols)(2,i)/(*vols)(2,i);

        (*(itsAbsDriftX[i]))[0] = (*absoluteDrifts)(i,0);
        (*(itsAbsDriftX[i]))[1] = (*absoluteDrifts)(i,1);
        (*(itsAbsDriftX[i]))[2] = (*absoluteDrifts)(i,2);
    }

    if(isInitSlices)
    {
        ARM_GP_Vector spaceSteps(3,0.0);
        ARM_SliceVector* slices = new ARM_SliceVector(nbSteps);
        (*slices)[0] = new ARM_SliceND(3,0,ARM_IntVector(3,0),ARM_IntVector(3,0),spaceSteps);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        for(i=1;i<nbSteps;++i)
        {
            spaceSteps[0] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[0]);
            spaceSteps[1] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[1]);
            spaceSteps[2] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[2]);

            (*slices)[i]    = new ARM_SliceND(3,i,ARM_IntVector(3,-i),ARM_IntVector(3,i),spaceSteps);
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler3D
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSampler3D::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
{
    size_t stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr zMean( new ARM_GP_Matrix(3,nbStates,0.0) );

    /// Get original X variable states to compute absolute markovian drift from model
    ARM_GP_MatrixPtr yStates;
    xStates = ComputeZtoXStates(sliceIdx,zStates,yStates);

    /// Add drift corrections, compute markovian drift then restore initial states
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*xStates)(FIRST_STATE_VARIABLE,stateIdx)   += (*driftCorrection)[FIRST_STATE_VARIABLE];
        (*xStates)(SECOND_STATE_VARIABLE,stateIdx)  += (*driftCorrection)[SECOND_STATE_VARIABLE];
        (*xStates)(THIRD_STATE_VARIABLE,stateIdx)   += (*driftCorrection)[THIRD_STATE_VARIABLE];
    }
    ARM_GP_MatrixPtr markovDriftX = GetModel()->MarkovianDrift(sliceIdx,xStates);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*xStates)(FIRST_STATE_VARIABLE,stateIdx)   -= (*driftCorrection)[FIRST_STATE_VARIABLE];
        (*xStates)(SECOND_STATE_VARIABLE,stateIdx)  -= (*driftCorrection)[SECOND_STATE_VARIABLE];
        (*xStates)(THIRD_STATE_VARIABLE,stateIdx)   -= (*driftCorrection)[THIRD_STATE_VARIABLE];
    }

    /// Get relative drift of Y
    ARM_GP_Vector *relDriftY = itsRelDriftY[sliceIdx];

    /// Get absolute drift of X
    ARM_GP_Vector* absDriftX = itsAbsDriftX[sliceIdx];

    /// Compute vol ratio for conversion from X to Y
    ARM_GP_Vector volRatio(3);
    volRatio[FIRST_STATE_VARIABLE] = itsVolY[FIRST_STATE_VARIABLE] / (*itsVolX)(FIRST_STATE_VARIABLE,sliceIdx);
    volRatio[SECOND_STATE_VARIABLE] = itsVolY[SECOND_STATE_VARIABLE] / (*itsVolX)(SECOND_STATE_VARIABLE,sliceIdx);
    volRatio[THIRD_STATE_VARIABLE] = itsVolY[THIRD_STATE_VARIABLE] / (*itsVolX)(THIRD_STATE_VARIABLE,sliceIdx);

    /// Get rotation from Y to Z at next slice time
    ARM_GP_MatrixPtr yToZ = itsYtoZ[sliceIdx+1];

    /// Compute conditionnal mean
    ARM_GP_Vector driftY(3);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        /// Absolute part
        driftY[FIRST_STATE_VARIABLE] = volRatio[FIRST_STATE_VARIABLE] * 
            ( (*absDriftX)[FIRST_STATE_VARIABLE] + (*markovDriftX)(FIRST_STATE_VARIABLE,stateIdx) );
        driftY[SECOND_STATE_VARIABLE] = volRatio[SECOND_STATE_VARIABLE] * 
            ( (*absDriftX)[SECOND_STATE_VARIABLE] + (*markovDriftX)(SECOND_STATE_VARIABLE,stateIdx) );
        driftY[THIRD_STATE_VARIABLE] = volRatio[THIRD_STATE_VARIABLE] * 
            ( (*absDriftX)[THIRD_STATE_VARIABLE] + (*markovDriftX)(THIRD_STATE_VARIABLE,stateIdx) );

        /// Relative part
        driftY[FIRST_STATE_VARIABLE] +=
            (*relDriftY)[FIRST_STATE_VARIABLE] * (*yStates)(FIRST_STATE_VARIABLE,stateIdx);
        driftY[SECOND_STATE_VARIABLE] +=
            (*relDriftY)[SECOND_STATE_VARIABLE] * (*yStates)(SECOND_STATE_VARIABLE,stateIdx);
        driftY[THIRD_STATE_VARIABLE] +=
            (*relDriftY)[THIRD_STATE_VARIABLE] * (*yStates)(THIRD_STATE_VARIABLE,stateIdx);

        /// Contributions in Z space
        (*zMean)(FIRST_STATE_VARIABLE,stateIdx) = (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*yToZ)(0,0) * driftY[FIRST_STATE_VARIABLE] +
            (*yToZ)(0,1) * driftY[SECOND_STATE_VARIABLE] +
            (*yToZ)(0,2) * driftY[THIRD_STATE_VARIABLE];

        (*zMean)(SECOND_STATE_VARIABLE,stateIdx) = (*zStates)(SECOND_STATE_VARIABLE,stateIdx) +
            (*yToZ)(1,0) * driftY[FIRST_STATE_VARIABLE] +
            (*yToZ)(1,1) * driftY[SECOND_STATE_VARIABLE] +
            (*yToZ)(1,2) * driftY[THIRD_STATE_VARIABLE];

        (*zMean)(THIRD_STATE_VARIABLE,stateIdx) = (*zStates)(THIRD_STATE_VARIABLE,stateIdx) +
            (*yToZ)(2,0) * driftY[FIRST_STATE_VARIABLE] +
            (*yToZ)(2,1) * driftY[SECOND_STATE_VARIABLE] +
            (*yToZ)(2,2) * driftY[THIRD_STATE_VARIABLE];
    }

    return zMean;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler3D
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_VectorPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSampler3D::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates) const
{
    size_t stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr xStates( new ARM_GP_Matrix(3,nbStates) );
    yStates = ARM_GP_MatrixPtr( new ARM_GP_Matrix(3,nbStates) );

    /// Get rotation matrix from Z to Y
    ARM_GP_MatrixPtr zToY = itsZtoY[sliceIdx];

    /// Compute vol ratio for conversion from Y to X
    ARM_GP_Vector volRatio(3);
    volRatio[FIRST_STATE_VARIABLE] = (*itsVolX)(FIRST_STATE_VARIABLE,sliceIdx) / itsVolY[FIRST_STATE_VARIABLE];
    volRatio[SECOND_STATE_VARIABLE] = (*itsVolX)(SECOND_STATE_VARIABLE,sliceIdx) / itsVolY[SECOND_STATE_VARIABLE];
    volRatio[THIRD_STATE_VARIABLE] = (*itsVolX)(THIRD_STATE_VARIABLE,sliceIdx) / itsVolY[THIRD_STATE_VARIABLE];
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        (*yStates)(FIRST_STATE_VARIABLE,stateIdx) =
            (*zToY)(0,0) * (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*zToY)(0,1) * (*zStates)(SECOND_STATE_VARIABLE,stateIdx) +
            (*zToY)(0,2) * (*zStates)(THIRD_STATE_VARIABLE,stateIdx);

        (*xStates)(FIRST_STATE_VARIABLE,stateIdx) = volRatio[FIRST_STATE_VARIABLE] * 
            (*yStates)(FIRST_STATE_VARIABLE,stateIdx);

        (*yStates)(SECOND_STATE_VARIABLE,stateIdx) = 
            (*zToY)(1,0) * (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*zToY)(1,1) * (*zStates)(SECOND_STATE_VARIABLE,stateIdx) +
            (*zToY)(1,2) * (*zStates)(THIRD_STATE_VARIABLE,stateIdx);

        (*xStates)(SECOND_STATE_VARIABLE,stateIdx) = volRatio[SECOND_STATE_VARIABLE] *
            (*yStates)(SECOND_STATE_VARIABLE,stateIdx);

        (*yStates)(THIRD_STATE_VARIABLE,stateIdx) = 
            (*zToY)(2,0) * (*zStates)(FIRST_STATE_VARIABLE,stateIdx) +
            (*zToY)(2,1) * (*zStates)(SECOND_STATE_VARIABLE,stateIdx) +
            (*zToY)(2,2) * (*zStates)(THIRD_STATE_VARIABLE,stateIdx);

        (*xStates)(THIRD_STATE_VARIABLE,stateIdx) = volRatio[THIRD_STATE_VARIABLE] *
            (*yStates)(THIRD_STATE_VARIABLE,stateIdx);
    }

    return xStates; 
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSampler3D
///	Routine: ComputeXtoZCenters
///	Returns: void
///	Action : Convert centers in X space to centers in Z space
////////////////////////////////////////////////////
void ARM_MarkovianDriftSampler3D::ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const
{
    if(globalDrifts == ARM_GP_MatrixPtr(NULL) )
        return;

    size_t sliceIdx,nbSlices = globalDrifts->rows();
    ARM_GP_MatrixPtr yToZ;
    double cy0,cy1,cy2;

    /// Change X to Y due to volatility change then rotate Y to Z
    for(sliceIdx=0;sliceIdx<nbSlices;++sliceIdx)
    {
        cy0 = itsVolY[FIRST_STATE_VARIABLE] / (*itsVolX)(FIRST_STATE_VARIABLE,sliceIdx) * (*globalDrifts)(sliceIdx,FIRST_STATE_VARIABLE);
        cy1 = itsVolY[SECOND_STATE_VARIABLE] / (*itsVolX)(SECOND_STATE_VARIABLE,sliceIdx) * (*globalDrifts)(sliceIdx,SECOND_STATE_VARIABLE);
        cy2 = itsVolY[THIRD_STATE_VARIABLE] / (*itsVolX)(THIRD_STATE_VARIABLE,sliceIdx) * (*globalDrifts)(sliceIdx,THIRD_STATE_VARIABLE);
        yToZ = itsGlobalYtoZ[sliceIdx];
        (*globalDrifts)(sliceIdx,FIRST_STATE_VARIABLE)  = (*yToZ)(0,0) * cy0 + (*yToZ)(0,1) * cy1 + (*yToZ)(0,2) * cy2;
        (*globalDrifts)(sliceIdx,SECOND_STATE_VARIABLE) = (*yToZ)(1,0) * cy0 + (*yToZ)(1,1) * cy1 + (*yToZ)(1,2) * cy2;
        (*globalDrifts)(sliceIdx,THIRD_STATE_VARIABLE)  = (*yToZ)(2,0) * cy0 + (*yToZ)(2,1) * cy1 + (*yToZ)(2,2) * cy2;
    }
}


////////////////////////////////////////////////////
////////////////////////////////////////////////////
/// ARM_MarkovianDriftSamplerND
////////////////////////////////////////////////////
////////////////////////////////////////////////////

////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerND
///	Routine:
///	Returns: 
///	Action : coperator=
////////////////////////////////////////////////////
ARM_MarkovianDriftSamplerND& ARM_MarkovianDriftSamplerND::operator=(const ARM_MarkovianDriftSamplerND& rhs )
{
    if( this != &rhs )
	{
        ARM_MarkovianDriftSamplerNDBase::operator =( rhs );
		itsDim = rhs.itsDim;
	}
    return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerND
///	Routine: Init
///	Returns: ARM_TimeStepsAndSlices*
///	Action : Initialise internal structures :
///          - volatilities for variable change
///          - PCA
///          - relative drift contributions
////////////////////////////////////////////////////
ARM_TimeStepsAndSlices* ARM_MarkovianDriftSamplerND::Init(const ARM_PricingModel& model, ARM_GP_VectorPtr& timeSteps, bool isInitSlices)
{
	/// initialise the dimension!
	itsDim = model.FactorCount();

    /// Set model
    SetModel(&model);

    /// Compute numerical method schedule
    if(timeSteps == ARM_GP_VectorPtr(NULL))
        timeSteps = GetScheduler()->ComputeTimeSteps( model );
    size_t i,j,k,nbSteps=timeSteps->size();
    ARM_GP_Vector yfSteps(*timeSteps);
    yfSteps /= K_YEAR_LEN;

    ARM_TimeStepsAndSlices* result = new ARM_TimeStepsAndSlices;
    result->itsTimeSteps = timeSteps;

#ifdef MDS_GLOBAL_ROTATION
    /// Get VCVs from model according to the schedule
    ARM_MatrixVector unsedLocalVCV;
    ARM_MatrixVector globalVCV;
    model.NumMethodStateLocalGlobalVariances(*timeSteps,unsedLocalVCV,globalVCV);
    ARM_GP_Matrix gloVCV(itsDim,itsDim);
#endif

    /// Get instantaneous volatilities and derivates from model
    /// sampled by the schedule
    ARM_GP_MatrixPtr vols,d1Vols,correls;
    model.VolatilitiesAndCorrelations(*timeSteps,vols,d1Vols,correls);

    size_t nbDims = vols->rows();
    if(nbDims != itsDim )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
        " : inconsistency between factors number of the model and the MD sampler dimension");

    itsVolX = vols;

    /// Set constant volatility for variable change
    itsVolY.resize(itsDim);
    double minVolY = GetMinStdDev();
	for( i=0; i<itsDim; ++i )
    {
		itsVolY[i] = (*itsVolX)(i,vols->cols()-1);
        if(itsVolY[i] < minVolY)
            itsVolY[i] = minVolY;
    }

    /// Compute rotation matrixes Y (correlated cst vol) to Z (independant) spaces
	size_t correlationSize = itsDim *(itsDim-1)/2;
	ARM_GP_Vector lastRho( correlationSize, -2.0 );
	ARM_GP_Vector Rho( correlationSize, -2.0 );
    ARM_GP_MatrixPtr lastZtoY,lastYtoZ;
    ARM_GP_Matrix vcv(itsDim,itsDim);
    ARM_GP_Vector eigenValues(itsDim),globalEigenValues(itsDim);
    ARM_GP_Vector* lastEigenValues;
    ARM_GP_Matrix* eigenVectors;
    ARM_GP_Matrix* globalEigenVectors;

    itsZtoY.resize(nbSteps);
    itsYtoZ.resize(nbSteps);
    itsZtoY[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(itsDim,itsDim,0.0) );
    itsYtoZ[0] = ARM_GP_MatrixPtr( new ARM_GP_Matrix(itsDim,itsDim,0.0) );

    itsGlobalYtoZ.resize(nbSteps);

    ARM_GP_Vector volRatio(itsDim);

    itsLocalVarZ.resize(nbSteps-1);
    itsGlobalVarZ.resize(nbSteps);
    double dt;
	size_t offseti,offset;
	bool SameCorrelationMatrix;
    for(i=0;i<nbSteps;++i)
    {
		SameCorrelationMatrix = true;
		for( j=0; j<correlationSize; ++j )
		{
            Rho[j] = (*correls)(j,i);
			if( SameCorrelationMatrix && fabs( lastRho[j]-Rho[j]) > K_DOUBLE_TOL )
				SameCorrelationMatrix = false;
		}

        if( !SameCorrelationMatrix )
        {
            /// Generate new rotations
			offseti = 0;
			for( j=0; j<itsDim; ++j )
			{
				for( k=0; k<j; ++k)
					vcv(j,k) = vcv(k,j);

				vcv(j,j) = itsVolY[j] * itsVolY[j];

				for( k=j+1, offset=offseti; k<itsDim; ++k, ++offset )
		            vcv(j,k) = itsVolY[j] * itsVolY[k] * Rho[ offset ];
				offseti += itsDim-j-1;
			}


            eigenVectors = ACPTransformation(&vcv,eigenValues);

            if(i<nbSteps-1)
            {
                /// itsZToY[i+1] is used to convert Z states at slice i+1 to Y states
                /// It depends on localVCV = ((ti+1-ti)*vols(ti)) between slice ti & ti+1
                itsZtoY[i+1]  = ARM_GP_MatrixPtr( eigenVectors );
                lastZtoY    = itsZtoY[i+1];

		        ARM_GP_Matrix* tmpMatrix = static_cast< ARM_GP_Matrix* >(eigenVectors->Clone());
		        tmpMatrix->transpose();
                itsYtoZ[i+1]  = ARM_GP_MatrixPtr( tmpMatrix );
                lastYtoZ    = itsYtoZ[i+1];
            }

            lastEigenValues = &eigenValues;
            lastRho = Rho;
        }
        else if(i<nbSteps-1)
        {
            /// Keep previous rotations
            itsZtoY[i+1]  = lastZtoY;
            itsYtoZ[i+1]  = lastYtoZ;
        }

        /// Local & global variances in Z space
        if(i<nbSteps-1)
        {
            itsLocalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
            dt = yfSteps[i+1] - yfSteps[i];

			for( j=0;j<itsDim; ++j )
	            (*(itsLocalVarZ[i]))[j] *= dt;
        }

#ifdef MDS_GLOBAL_ROTATION

        /// Correct computation using diagonalized global VCV
        for(j=0;j<itsDim;++j)
        {
            volRatio[j] = itsVolY[j] / (*itsVolX)(j,i);
            gloVCV(j,j) = volRatio[j] * volRatio[j] * (*(globalVCV[i]))(j,j);
            for(k=0;k<j;++k)
            {
                gloVCV(j,k)= volRatio[j] * volRatio[k] * (*(globalVCV[i]))(j,k);
                gloVCV(k,j)=gloVCV(j,k);
            }
        }
        globalEigenVectors = ACPTransformation(&gloVCV,globalEigenValues);
        itsGlobalVarZ[i] = static_cast< ARM_GP_Vector* >(globalEigenValues.Clone());

        globalEigenVectors->transpose();
        itsGlobalYtoZ[i] = ARM_GP_MatrixPtr( globalEigenVectors );

#else

        /// Global variances in Z space (not correct because correlation, mean reversion...
        /// but something near Dimitri's code)
        itsGlobalVarZ[i]  = static_cast< ARM_GP_Vector* >(lastEigenValues->Clone());
		for( j=0;j<itsDim; ++j )
			(*(itsGlobalVarZ[i]))[j] *= yfSteps[i];

        itsGlobalYtoZ[i] = lastYtoZ;

#endif

    }

#ifdef MDS_GLOBAL_ROTATION
    /// Free temporary VCV matrixes returned by the model
    DeletePointorVector<ARM_GP_Matrix>(unsedLocalVCV);
    DeletePointorVector<ARM_GP_Matrix>(globalVCV);
#endif

    /// Get instantaneous drifts (deterministic part)
    ARM_GP_MatrixPtr relativeDrifts;
    ARM_GP_MatrixPtr absoluteDrifts;
    model.EulerLocalDrifts(*timeSteps,relativeDrifts,absoluteDrifts);

    /// Compute absolute drift of X and relative drift of Y
    itsAbsDriftX.resize(nbSteps-1);
    itsRelDriftY.resize(nbSteps-1);
    for(i=0;i<nbSteps-1;++i)
    {
        dt = yfSteps[i+1] - yfSteps[i];

        itsAbsDriftX[i] = new ARM_GP_Vector(itsDim);
		itsRelDriftY[i] = new ARM_GP_Vector(itsDim);

		for( j=0;j<itsDim; ++j )
		{
			(*(itsRelDriftY[i]))[j] = (*relativeDrifts)(i,j) - dt*(*d1Vols)(j,i)/(*vols)(j,i);
	        (*(itsAbsDriftX[i]))[j] = (*absoluteDrifts)(i,j);
		}
    }

    if(isInitSlices)
    {
        ARM_GP_Vector spaceSteps(itsDim,0.0);
        ARM_SliceVector* slices = new ARM_SliceVector(nbSteps);
        (*slices)[0] = new ARM_SliceND(itsDim,0,ARM_IntVector(itsDim,0),ARM_IntVector(itsDim,0),spaceSteps);

	    /// First probas and Arrow Debreu price are trivial and equal to 1.0
	    (*slices)[0]->SetSpotProbas( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );
	    (*slices)[0]->SetArrowDebreuPrices( ARM_GP_VectorPtr( new ARM_GP_Vector(1,1.0) ) );

        for(i=1;i<nbSteps;++i)
        {
		    for( j=0;j<itsDim; ++j )
	            spaceSteps[j] = ARM_NumericConstants::ARM_SQRT_3 * sqrt((*(itsLocalVarZ[i-1]))[j]);
            (*slices)[i]    = new ARM_SliceND(itsDim,i,ARM_IntVector(itsDim,-i),ARM_IntVector(itsDim,i),spaceSteps);
        }

        result->itsSlices    = ARM_SliceVectorPtr(slices);
    }

    return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerND
///	Routine: ComputeConditionalMean
///	Returns: ARM_GP_MatrixPtr
///	Action : Compute conditional mean in Z space from current slice
///          (referred by sliceIdx) to next one
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSamplerND::ComputeConditionalMean(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, const ARM_GP_VectorPtr& driftCorrection, ARM_GP_MatrixPtr& xStates) const
{
    size_t i,j,stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr zMean( new ARM_GP_Matrix(itsDim,nbStates,0.0) );

    /// Get original X variable states to compute absolute markovian drift from model
    ARM_GP_MatrixPtr yStates;
    xStates = ComputeZtoXStates(sliceIdx,zStates,yStates);

    /// Add drift corrections, compute markovian drift then restore initial states
	for( i=0; i<itsDim; ++i )
    {
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
            (*xStates)(i,stateIdx)   += (*driftCorrection)[i];
    }
    ARM_GP_MatrixPtr markovDriftX = GetModel()->MarkovianDrift(sliceIdx,xStates);
	for( i=0; i<itsDim; ++i )
    {
        for(stateIdx=0;stateIdx<nbStates;++stateIdx)
            (*xStates)(i,stateIdx)   -= (*driftCorrection)[i];
    }

    /// Get relative drift of Y
    ARM_GP_Vector *relDriftY = itsRelDriftY[sliceIdx];

    /// Get absolute drift of X
    ARM_GP_Vector* absDriftX = itsAbsDriftX[sliceIdx];

    /// Compute vol ratio for conversion from X to Y
    ARM_GP_Vector volRatio(itsDim);
	
	for( i=0; i<itsDim; ++i )
		volRatio[i] = itsVolY[i] / (*itsVolX)(i,sliceIdx);

    /// Get rotation from Y to Z at next slice time
    ARM_GP_MatrixPtr yToZ = itsYtoZ[sliceIdx+1];

    /// Compute conditionnal mean
    ARM_GP_Vector driftY(itsDim);
    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
        /// Absolute & relative parts...
		for( i=0; i<itsDim; ++i )
            driftY[i] = volRatio[i] * ((*absDriftX)[i] + (*markovDriftX)(i,stateIdx)) +
                (*relDriftY)[i] * (*yStates)(i,stateIdx);

        /// ... then do the rotation to get contribution in Z space
		for( i=0; i<itsDim; ++i )
		{
			(*zMean)(i,stateIdx) = (*zStates)(i,stateIdx);
			for( j=0; j<itsDim; ++j )
				(*zMean)(i,stateIdx) += (*yToZ)(i,j) * driftY[j];
		}
    }

    return zMean;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerND
///	Routine: ComputeZtoXStates
///	Returns: ARM_GP_VectorPtr
///	Action : Conversion from the numerical state variables
///          (Z space) to the process states variables (X space)
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovianDriftSamplerND::ComputeZtoXStates(size_t sliceIdx, const ARM_GP_MatrixPtr& zStates, ARM_GP_MatrixPtr& yStates) const
{
    size_t i,j,stateIdx,nbStates = zStates->cols();
    ARM_GP_MatrixPtr xStates(   new ARM_GP_Matrix(itsDim,nbStates) );
    yStates = ARM_GP_MatrixPtr( new ARM_GP_Matrix(itsDim,nbStates) );

    /// Get rotation matrix from Z to Y
    ARM_GP_MatrixPtr zToY = itsZtoY[sliceIdx];

    /// Compute vol ratio for conversion from Y to X
    ARM_GP_Vector volRatio(itsDim);
	
	for( i=0; i<itsDim; ++i )
		volRatio[i] = (*itsVolX)(i,sliceIdx) / itsVolY[i];

    for(stateIdx=0;stateIdx<nbStates;++stateIdx)
    {
		for( i=0; i<itsDim; ++i )
		{
			(*yStates)(i,stateIdx) = 0;
			for(j=0; j<itsDim; ++j )
				(*yStates)(i,stateIdx) += (*zToY)(i,j) * (*zStates)(j,stateIdx);
	        (*xStates)(i,stateIdx) = volRatio[i] * (*yStates)(i,stateIdx);
		}
	}
    return xStates; 
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovianDriftSamplerND
///	Routine: ComputeXtoZCenters
///	Returns: void
///	Action : Convert centers in X space to centers in Z space
////////////////////////////////////////////////////
void ARM_MarkovianDriftSamplerND::ComputeXtoZCenters(ARM_GP_MatrixPtr& globalDrifts) const
{
    if(globalDrifts == ARM_GP_MatrixPtr(NULL) )
        return;

    size_t sliceIdx,nbSlices = globalDrifts->rows();
    ARM_GP_MatrixPtr yToZ;
    size_t i,j;
    ARM_GP_Vector cy(itsDim);

    /// Change X to Y due to volatility change then rotate Y to Z
    for(sliceIdx=0;sliceIdx<nbSlices;++sliceIdx)
    {
        for(i=0;i<itsDim;++i)
            cy[i] = itsVolY[i] / (*itsVolX)(i,sliceIdx) * (*globalDrifts)(sliceIdx,i);
        yToZ = itsGlobalYtoZ[sliceIdx];
        for(i=0;i<itsDim;++i)
        {
            (*globalDrifts)(sliceIdx,i) = 0.0;
            for(j=0;j<itsDim;++j)
                (*globalDrifts)(sliceIdx,i)  += (*yToZ)(i,j) * cy[j];
        }
    }
}


CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

