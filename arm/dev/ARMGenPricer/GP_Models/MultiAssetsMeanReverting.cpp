/*!
 *
 * Copyright (c) IXIS CIB Paris 2005 Paris
 *
 *	\file MultiAssetsReverting.cpp
 *
 *  \brief multi assets mean reverting model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date January 2005
 */

#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/MultiAssetsMeanReverting.h"


/// gpinfra
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingstates.h"

/// gpmodels
#include "gpmodels/HW.h"
#include "gpmodels/ModelParamsHW1F.h"

#include "gpbase/cloneutilityfunc.h"


CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsMeanReverting
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
////////////////////////////////////////////////////
ARM_MultiAssetsMeanReverting::ARM_MultiAssetsMeanReverting(
		const ARM_ModelNameMap*	modelNameMap,
		const ARM_CurveMatrix* correlMatrix)
:	
	ARM_MultiAssetsModel( modelNameMap, correlMatrix )
{}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsMeanReverting
///	Routines: Constructor
///	Returns :
///	Action  : builds the object
////////////////////////////////////////////////////
ARM_MultiAssetsMeanReverting::ARM_MultiAssetsMeanReverting(
		const ARM_ModelNameMap*	modelNameMap,
		const ARM_GP_Matrix* correlMatrix)
:	
	ARM_MultiAssetsModel( modelNameMap, correlMatrix )
{}
	
////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsMeanReverting
///	Routines: Copy Constructor
///	Returns :
///	Action  : builds the object
////////////////////////////////////////////////////
ARM_MultiAssetsMeanReverting::ARM_MultiAssetsMeanReverting(const ARM_MultiAssetsMeanReverting& rhs)
:	ARM_MultiAssetsModel(rhs)
{
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routines: ~ARM_MultiAssetsModel
///	Returns :
///	Action  : destructor
////////////////////////////////////////////////////	
ARM_MultiAssetsMeanReverting::~ARM_MultiAssetsMeanReverting()	
{
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsModel
///	Routine : MeanRevertingCompatible
///	Returns : void
///	Action  : tells whether all the base models are mean reverting compatible
////////////////////////////////////////////////////
bool ARM_MultiAssetsMeanReverting::MeanRevertingCompatible(const ARM_ModelNameMap& modelMap)
{
	if (modelMap.UsedModelsSize()==1)
		return false;

	for( ARM_ModelNameMap::const_iterator iter=modelMap.begin(); iter!=modelMap.end(); modelMap.getNextUsedIter(iter) )
	{
		if (!((*iter).Model()->IsMeanRevertingCompatible()) )
			return false;
	}

	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsMeanReverting
///	Routine : NumMethodStateLocalVariances
///	Returns : void
///	Action  : computes the integrated local covariances 
/// of mean reverting processes
////////////////////////////////////////////////////
void ARM_MultiAssetsMeanReverting::NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t factorNb	= FactorCount();
	const ARM_CurveMatrix* const correlCurveMatrix  = GetCorrelMatrix();
	ARM_GP_Matrix correlMatrix;
	const ARM_ModelNameMap* const modelMap = GetModelMap();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size() != 0" );
#endif

	localVariances.resize(nbSteps-1);

	size_t i,j,k;
	ARM_ModelNameMap::const_iterator iter,iter2;
	double elem;
	
	/// Basic MeanReverting VarCovar Matrix
	for(i=0;i<nbSteps-1;++i)
	{
		nextStep			= timeSteps[i+1];
		correlMatrix = correlCurveMatrix->Interpolate(step);
		localVariances[i]	= new ARM_GP_TriangularMatrix(factorNb,0.0);

		j = 0;
		for( iter=modelMap->begin(); iter!=modelMap->end(); modelMap->getNextUsedIter(iter), ++j )
		{
			const ARM_ModelParams* const modelParams = (*iter).Model()->GetModelParams();
			const ARM_ModelParamsHW1FStd* const modelParamsHW1FStd = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( modelParams );

			k=0;
			if( modelParamsHW1FStd )
			{
				(*(localVariances[i]))(j,j) = modelParamsHW1FStd->StateLocalVariance(step,nextStep,nextStep);

				for( iter2=modelMap->begin(); k<j; modelMap->getNextUsedIter(iter2), ++k )
				{
					const ARM_ModelParams* const modelParams_2 = (*iter2).Model()->GetModelParams();
					const ARM_ModelParamsHW1FStd* const modelParamsHW1FStd_2 = dynamic_cast<const ARM_ModelParamsHW1FStd* const>( modelParams_2 );

					elem = (modelParamsHW1FStd_2 ? (correlMatrix)(j,k)* ARM_ModelParamsHW1F::HW1FStateCovariance(modelParamsHW1FStd, modelParamsHW1FStd_2,step,nextStep,nextStep) : 0.0 );
					(*(localVariances[i]))(j,k)	= elem;
					(*(localVariances[i]))(k,j)	= elem;
				}
			}
			else
			{
				for( iter2=modelMap->begin(); k<j; modelMap->getNextUsedIter(iter2), ++k )
				{
					(*(localVariances[i]))(j,k)	= 0.0;
					(*(localVariances[i]))(k,j)	= 0.0;
				}
				(*(localVariances[i]))(j,j)	= 1.;
			}
		}
		step=nextStep;
	}

	/// Models bring their corrections
	for( iter=modelMap->begin(); iter!= modelMap->end(); GetModelMap()->getNextUsedIter(iter) )
		(*iter).Model()->NumMethodStateLocalCovariances( timeSteps, localVariances, (*this) );
}

////////////////////////////////////////////////////
///	Class   : ARM_MultiAssetsMeanReverting
///	Routine : NumMethodStateGlobalVariances
///	Returns : void
///	Action  : computes the integrated global covariances of 
/// mean reverting processes
////////////////////////////////////////////////////
void ARM_MultiAssetsMeanReverting::NumMethodStateGlobalVariances( 
		const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	size_t factorNb	= FactorCount();
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= 0 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size() != 0" );
#endif

	globalVariances.resize(nbSteps);

	// There is no variance in the numerical method

	ARM_GP_Matrix identity(factorNb,factorNb);

	size_t i,j;

	for (i = 0; i < factorNb; ++i)
	{
		for(j = 0; j < i; ++j)
		{
			identity(i,j) = identity(j,i) = 0.0;
		}
		identity(i,i) = 1.0;
	}
	
	for(i=0;i<nbSteps;++i)
	{
		nextStep			= timeSteps[i];
		globalVariances[i]	= static_cast<ARM_GP_Matrix*>(identity.Clone());
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : call MCModelStatesFromToNextTime on the sub models
////////////////////////////////////////////////////

void ARM_MultiAssetsMeanReverting::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#ifdef __GP_STRICT_VALIDATION

	ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	ARM_GP_MatrixPtr& modelStates = states->GetModelStates();
	size_t bucketSize =  states->size();

	if( numMethodStates->rows() != FactorCount() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "numMethodStates->rows() != FactorCount()!" );
	if( modelStates->rows() != FactorCount() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelStates->rows() != FactorCount()!" );
	if( modelStates->cols() != bucketSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelStates->cols() != bucketSize!" );
	if( numMethodStates->cols() != bucketSize )
		ARM_THROW( ERR_INVALID_ARGUMENT, "numMethodStates->cols() != bucketSize!" );
#endif

	for( ARM_ModelNameMap::const_iterator iter=GetModelMap()->sortedBegin(); iter!=GetModelMap()->sortedEnd(); GetModelMap()->getSortedNextUsedIter(iter) )
			(*iter).Model()->MCModelStatesFromToNextTime( states, timeIndex);
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_MultiAssetsMeanReverting::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

