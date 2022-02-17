/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 */

/// this header comes first as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/datestrip.h"
#include "gpbase/vectormanip.h"

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/SmiledFRM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpmodels/ProcessBuilderSmiledFRM.h"
#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillasecuritydensity.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <inst/swaption.h>




CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SmiledFRM::ARM_SmiledFRM( const ARM_ZeroCurvePtr& zc,
		 const ARM_ModelParams* params, 
		 size_t timeStepsNb,
		 size_t gridSize,
		 double stdDevNb,
		 bool skipPDE,
		 bool allowInterpol,
		 ARM_ModelParamsSmiled::CalibProxy swaptionApprox,
		 bool withRescalling)
:	ARM_PricingModelIR(zc,params),
	itsNumericalModelFitter(NULL), 
	itsResetDates(NULL),
	itsStartDates(NULL),
	itsEndDates(NULL),
	itsProcess(0),
	itsEigenValues(0),
	itsNumeraireTimeIndex(0),
	itsCalibrationStatus(true),
	itsCalibrated(false),
	itsSkipPDE(skipPDE),
	itsAllowInterpol(allowInterpol),
	itsTimeStepsNb(timeStepsNb),
	itsGridSize(gridSize),
	itsStdDevNb(stdDevNb),
	itsCurrentArg(NULL),
	itsNewArgNeeded(false),
	itsSwaptionApprox(swaptionApprox),
	itsWeight(0),
	itsWithRescalling(withRescalling),
	itsCachedModelStatesLocalVars(0)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SmiledFRM::ARM_SmiledFRM(const ARM_SmiledFRM& rhs)
:	ARM_PricingModelIR(rhs),
	itsNumericalModelFitter(CreateClone(rhs.itsNumericalModelFitter)), 
	itsResetDates(rhs.itsResetDates),
	itsStartDates(rhs.itsStartDates),
	itsEndDates(rhs.itsEndDates),
	itsNumeraireTimeIndex(rhs.itsNumeraireTimeIndex),
	itsCalibrationStatus(rhs.itsCalibrationStatus),
	itsCalibrated(rhs.itsCalibrated),
	itsSkipPDE(rhs.itsSkipPDE),
	itsAllowInterpol(rhs.itsAllowInterpol),
	itsTimeStepsNb(rhs.itsTimeStepsNb),
	itsGridSize(rhs.itsGridSize),
	itsStdDevNb(rhs.itsStdDevNb),
	itsCurrentArg(rhs.itsCurrentArg),
	itsNewArgNeeded(rhs.itsNewArgNeeded),
	itsSwaptionApprox(rhs.itsSwaptionApprox),
	itsWeight(rhs.itsWeight),
	itsWithRescalling(rhs.itsWithRescalling),
	itsCachedModelStatesLocalVars(rhs.itsCachedModelStatesLocalVars)
{
	DuplicateCloneablePointorVectorInPlace<ARM_ProcessBuilder>( rhs.itsProcess, itsProcess );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SmiledFRM::~ARM_SmiledFRM()
{
	DeletePointorVector<ARM_ProcessBuilder>( itsProcess );
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: ComputeNumeraireTimeIndex
///	Returns : int
///	Action  : Compute the numeraire time index based on the numeraire
/// current date. 
////////////////////////////////////////////////////
void ARM_SmiledFRM::ComputeNumeraireTimeIndex()
{
	double numeraireTime = GetNumeraire()->GetMaturity();	
	itsNumeraireTimeIndex = lower_boundPosWithPrecision(itsEndDates,numeraireTime)-1;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: ComputeModelTimes
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_SmiledFRM::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetDates.Clone());
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::DiscountFactorNoInterpol( const string& curveName, 
	  double evalTime,
	  double maturityTime, 
	  const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( maturityTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,DiscountFactor: MaturityTime after schedule's last date");


	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );
	
	if( DoesResetDateExist( evalTime ) )
	{
		ARM_GP_Vector::iterator iter;
		ARM_GP_Vector::iterator iter1;
		ARM_GP_Vector::iterator iter2;

		size_t ResetIdx = IdxFromValueWithTol( itsResetDates, evalTime, 1. );
		size_t StartIdx = ResetIdx;

		if ( maturityTime >= itsStartDates[itsStartDates.size()-1] )
		{
			StartIdx = itsStartDates.size()-1;
			while (!IsOnSamePath( ResetIdx, StartIdx )) {StartIdx--;}			
		}
		else
		{
			if ( itsStartDates[StartIdx] > maturityTime )
			{
			}
			else {
				while (itsStartDates[StartIdx] <= maturityTime) { StartIdx++;};
				StartIdx -= 1;
			}
		}

		if ( IsOnSamePath( ResetIdx, StartIdx ) && FastCompute() )
		{
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, StartIdx, modelNb, states);

			for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
				(*iter) *= (*iter1);
		}
		else
		{
			ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, StartIdx, itsResetDates.size(), modelNb, states);
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, itsResetDates.size(), modelNb, states);

			for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
				(*iter) *= (*iter2)/(*iter1);
		}

		if (maturityTime != itsStartDates[StartIdx])
		{
			double theta = (maturityTime - itsStartDates[StartIdx])/360.;
			double delta = itsProcess[StartIdx]->GetDelta();
			double mu =		GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
						/	GetZeroCurve()->DiscountPrice(itsStartDates[StartIdx]/K_YEAR_LEN)
						*	( 1 + theta * itsProcess[StartIdx]->GetFwdRate());
			
			ARM_VectorPtr dfratio = itsProcess[StartIdx]->DFRatio(evalTime,states,modelNb,StartIdx);

			for (iter = result->begin(),iter2 = dfratio->begin() ; iter != result->end() ; ++iter,++iter2)
				(*iter) *= mu / ( 1. + theta * ((*iter2)-1.)/delta );
		}
		
		return static_cast<ARM_VectorPtr>(result);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,DiscountFactor: evalTime not in schedule!!");		
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if (itsAllowInterpol)
		return ForwardDiscountFactor(curveName,evalTime,evalTime,maturityTime,states);
	else
		return DiscountFactorNoInterpol(curveName,evalTime,maturityTime,states);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::ForwardDiscountFactor( const string& curveName, double evalTime, double startTime, double endTime, const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN) ));
	}

	if( endTime > itsEndDates[itsEndDates.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,DiscountFactor: endTime after schedule's last date");

	if( evalTime < itsResetDates[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,DiscountFactor: evalTime before first reset");

	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );

	ARM_GP_Vector::iterator iter;
	ARM_GP_Vector::iterator iter1;
	ARM_GP_Vector::iterator iter2;

	size_t FirstIdx = IndexOfFirstHigherEqInVector_DefaultLast(startTime, itsStartDates); //cas où start dans dernière période!
	if (startTime > itsStartDates[FirstIdx])
		FirstIdx++;

	size_t LastIdx  = IndexOfLastLowerEqInVector(endTime,itsStartDates);
	
	if ( (LastIdx<FirstIdx) || ( IsOnSamePath( FirstIdx, LastIdx ) ) && FastCompute())
	{
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, LastIdx, modelNb, states);

		for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
			(*iter) *= (*iter1);
	}
	else
	{
		ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, LastIdx, itsResetDates.size(), modelNb,states);
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, itsResetDates.size(), modelNb, states);

		for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
			(*iter) *= (*iter2)/(*iter1);
	}

	if ( (LastIdx < FirstIdx) || (startTime < itsStartDates[FirstIdx]))
	{
		double mat				= LastIdx < FirstIdx?endTime:itsStartDates[FirstIdx];
		double theta			= (mat - startTime)/360.;
		size_t idx				= FirstIdx>0?FirstIdx-1:FirstIdx;
		double delta1			= itsProcess[idx]->GetDelta();
		double fwd1				= itsProcess[idx]->GetFwdRate();
		double fwd				= (GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)/GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)-1.)/theta;//( 1. - weight2 ) * fwd1 + weight2 * fwd2;

		ARM_VectorPtr dfratio1	= itsProcess[idx]->DFRatio(evalTime,states,modelNb,idx,true);
		
		double mu				=		GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)
									/	GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)
									*	( 1 + theta * fwd);

		for (iter = result->begin(),iter1 = dfratio1->begin(); iter != result->end() ; ++iter,++iter1)
				(*iter) *= mu / ( 1. + theta * fwd/fwd1 * ((*iter1)-1.)/delta1 );
	}

	if ( (endTime > itsStartDates[LastIdx]) && (LastIdx >= FirstIdx))
	{
			double theta = (endTime - itsStartDates[LastIdx])/360.;
			double delta = itsProcess[LastIdx]->GetDelta();
				
			double mu =		GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)
						/	GetZeroCurve()->DiscountPrice(itsStartDates[LastIdx]/K_YEAR_LEN)
						*	( 1 + theta * itsProcess[LastIdx]->GetFwdRate());
					
			ARM_VectorPtr dfratio = itsProcess[LastIdx]->DFRatio(evalTime,states,modelNb,LastIdx,true);

			for (iter = result->begin(),iter2 = dfratio->begin() ; iter != result->end() ; ++iter,++iter2)
				(*iter) *= mu / ( 1. + theta * ((*iter2)-1.)/delta );
	}

	return static_cast<ARM_VectorPtr>(result);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t IdxFrom, size_t IdxTo , size_t modelNb, const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr dfratio;
	ARM_GP_Vector* result	= new ARM_GP_Vector( states->size(),1.0 );
		
	for (size_t k = IdxFrom; k < IdxTo ; k++ )
	{
		if ( itsWeight(IdxFrom,k) == 1. )
		{
			dfratio = itsProcess[k]->DFRatio(evalTime,states,modelNb,k);

			ARM_GP_Vector::iterator iter2=dfratio->begin();

			for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter, ++iter2)
				(*iter) /= (*iter2);
		}
	}
	return static_cast<ARM_VectorPtr>(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PostInit
///	Returns: nothing
///	Action : Function to init after the numerical method
///			 has been initialised.. non const
////////////////////////////////////////////////////
void ARM_SmiledFRM::PostInit()
{
	ComputeNumeraireTimeIndex();
}



////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledFRM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double currTime  = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime  = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb = FactorCount();
	size_t statesNb  = states->size();
	size_t modelNb	 = GetModelNb();

	size_t j,k;
	size_t fwdsNb  = itsResetDates.size();
	size_t iLast   = fwdsNb-1;

	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetDates[iFirst]<nextTime)
	{
		++iFirst;
	}
	if (itsAllowInterpol && iFirst>0) iFirst--;

	//parait débile mais indispensable pour les échanciers non triviaux!!!
	size_t iFirstMartingale = iFirst;
	while( itsEndDates[iFirstMartingale] < itsEndDates[iLast]) iFirstMartingale++;
	
	const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();
	ARM_GP_Vector* eigen = itsEigenValues[timeIndex];

	int i;
	double old,next,modStd;

	ARM_VectorPtrVector dxLog(iLast-iFirst+1);

	k = iFirst;
	for (ARM_VectorPtrVector::iterator iter = dxLog.begin() ; iter != dxLog.end() ; ++iter , ++k){
		(*iter) = itsProcess[k]->DxLogDFRatio(currTime,states,modelNb,k,true);
	}

	ARM_GP_Matrix E(fwdsNb,factorsNb,0.);
	
	for(k = 0 ; k < statesNb ; ++k )
	{
		for ( i=iFirstMartingale ; i<fwdsNb ; i++ )
			for (j = 0 ; j < factorsNb ; j++ )
				E(i,j) = 0.;

		i = itsNumeraireTimeIndex;
		for (ARM_VectorPtrVector::reverse_iterator iterDx = dxLog.rbegin(); iterDx != dxLog.rend() ; ++iterDx,--i)
		{
			old			= states->GetModelState(k,i+modelNb);
			next		= old;
			j			= 0;
			size_t prec	= i+1;

			if (itsEndDates[i] == itsEndDates[iLast])
				prec = i;
			else while (!IsOnSamePath(i,prec)) prec++;
			
			for (ARM_GP_Vector::iterator iter = eigen->begin(); iter != eigen->end() ; ++iter,++j )
			{
				modStd	=	(*modelLocalVar[timeIndex])(i,j);
				next	+=	modStd * ( states->GetNumMethodState(k,j+modelNb) + (*iter) * E(prec,j));
				E(i,j)	=	E(prec,j) - modStd * (*iterDx)->Elt(k);
			}
			
			states->SetModelState(k,i+modelNb,next);
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_SmiledFRM::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(!params.DoesModelParamExist(ARM_ModelParamType::Hump) || 
       !params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ValidateModelParams: hump and beta correlation structures are needed for Smiled FRM" );
	}
	
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routine : SetNumeraire
///	Returns : bool
///	Action  : 
////////////////////////////////////////////////////
void ARM_SmiledFRM::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	if(	numerairePtr->GetType() != ARM_Numeraire::TerminalZc )
    {		
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,SetNumeraire: only numeraire supported is Terminal ZC" );
    }
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routine : ValidateCalibMethod
///	Returns : bool
///	Action  : validate calib method
////////////////////////////////////////////////////

void ARM_SmiledFRM::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);
	if (!itsCalibrated && !isNumerical)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ValidateCalibMethod: should calib functionals first!!!" );

	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_SmiledFRM::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
	/// nothing as we need the full path to compute localDrift
	/// and compute drift on the fly
	size_t factorsNb = FactorCount();
	size_t nbSteps	 = timeSteps.size();
	size_t modelNb	 = GetModelNb();

    /// Initialise a Relative Drift matrix to 1 (no drift)
    relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( (nbSteps-1)*(modelNb+1), factorsNb, 1.0 ) );
    absoluteDrifts	= ARM_GP_MatrixPtr( NULL);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_SmiledFRM::NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t totalFwds     = itsResetDates.size(),
		   factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   startTimePos  = 0;
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	if(itsEigenValues.empty())
		const_cast< ARM_SmiledFRM* >(this)->ModelStateLocalVariancesAndStdDev(timeSteps);

	for(size_t i=0;i<timeStepsSize-1 ;++i)
	{
		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);
		for (size_t k=0; k<factorsNb;k++)
			(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_SmiledFRM::NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	size_t totalFwds     = itsResetDates.size(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (size_t i=0;i<timeStepsSize;i++){
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds,1.0);
	}
	
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_SmiledFRM::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
//	ARM_PricingModel::ModelStateLocalVariancesAndStdDev(timeSteps);	
	
	ARM_MatrixVector auxLocalCov;

	ModelStateLocalVariances( timeSteps, auxLocalCov);

	ModelStateLocalStdDev(timeSteps, auxLocalCov);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: NumMethodStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_SmiledFRM::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{
	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_SmiledFRM::ModelStateLocalStdDev(const ARM_GP_Vector& timeSteps, const ARM_MatrixVector& localVariances)
{
	ARM_MatrixVector modStateLocalVars;

	if (!itsCachedTimeSteps.IsNull() && (timeSteps == *(itsCachedTimeSteps)))
	{
		DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Matrix>(itsCachedModelStatesLocalVars,modStateLocalVars);
	}
	else
	{
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>(localVariances, modStateLocalVars);

		size_t totalFwds     = itsResetDates.size(),
		   factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb,
		   startTimePos  = 0,
		   i,j,k, ii;
		double	fromTime= timeSteps[0];

		DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
		itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

		ARM_GP_Vector eigenValues(totalFwds);		
		for(i=0, ii = offsetIndex;i<timeStepsSize-1 ;++i, ++ii)
		{
			/// initialize everything
			itsEigenValues[i] = new ARM_GP_Vector(factorsNb,0.);
			
			ARM_GP_Matrix* ACPMatrix=ACPTransformation(modStateLocalVars[ii],eigenValues);

			int maxRank = -1;
			for(k=0;k<factorsNb;++k)
				if(fabs(eigenValues[k]) < K_NEW_DOUBLE_TOL)
				{
					eigenValues[k] = 0.;
					if (maxRank == -1)
						maxRank = k;
				}
				else
				{
					if (eigenValues[k] < 0.)
					{
	#if defined(__GP_STRICT_VALIDATION)
							ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariancesAndStdDev: negative eigenvalue" );
	#endif
							eigenValues[k] = 0.;
					}
				}
			
			if (maxRank == -1) maxRank = factorsNb;
			size_t effectiveRank = ( maxRank < factorsNb ? maxRank : factorsNb);
			//rescaling
			for (j=0;j<totalFwds;j++)
			{
				double sum=0.;
				for(k=0;k<effectiveRank-1;++k){
					sum+=(*ACPMatrix)(j,k)*(*ACPMatrix)(j,k)*eigenValues[k];
				}
				double sgn = ((*ACPMatrix)(j,k)>0.?1.:-1.);
				double res = (*modStateLocalVars[ii])(j,j)-sum;
				if (res>K_NEW_DOUBLE_TOL){
					(*ACPMatrix)(j,k)=sgn*sqrt(res/eigenValues[k]);
				}
			}

			for(k=0;k<factorsNb;++k)
			{
				(*itsEigenValues[i])(k) =	eigenValues[k];
				for(j=0; j<totalFwds; ++j)	{	
					(*modStateLocalVars[ii])(j,k) =	(*ACPMatrix)(j,k);
				}
			}

			delete ACPMatrix;
		}

		itsCachedTimeSteps = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(timeSteps.Clone()));
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>(modStateLocalVars, itsCachedModelStatesLocalVars);
	}

	SetModelStateLocalVars(modStateLocalVars);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: ModelStateLocalVariances
///	Returns : void
///	Action  : ModelStateLocalVariancess
////////////////////////////////////////////////////
void ARM_SmiledFRM::ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	/// The model local vars is now cached very usefull
	/// for the importance sampling optimization
	if (!itsCachedTimeSteps.IsNull() && (timeSteps == *(itsCachedTimeSteps)))
	{
		 //DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Matrix>(itsCachedModelStatesLocalVars,localVariances);
		return;
	}
	else
	{
		size_t totalFwds     = itsResetDates.size(),
		   factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb,
		   startTimePos  = 0,
		   i,j,k, ii;
		double	fromTime= timeSteps[0],
			toTime;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariances: localDrifts.size() != offsetIndex" );
#endif

    localVariances.resize((timeStepsSize-1)*(modelNb+1));

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariances: first time step != 0" );
#endif

		for(i=0, ii = offsetIndex;i<timeStepsSize-1 ;++i, ++ii)
		{
			/// initalize the toTime
			toTime  = timeSteps[i+1];

			/// initialize everything
			localVariances[ii] = new ARM_GP_Matrix(totalFwds,totalFwds);
			
			/// get the first bigger time
			while(startTimePos< itsResetDates.size()
				&& itsResetDates[startTimePos] < toTime )
				++startTimePos;

			if (itsAllowInterpol && startTimePos>0) startTimePos--;	

			for(j=0; j<startTimePos; ++j)
			{
				for(k=j;k<totalFwds;++k)
					(*localVariances[ii])(j,k)=(*localVariances[ii])(k,j)=0.;
			}
			for(j=startTimePos; j<totalFwds; ++j)
			{
				for(k=j;k<totalFwds;++k)
					(*localVariances[ii])(j,k)=(*localVariances[ii])(k,j)=((ARM_ModelParamsSmiled*) GetModelParams())->IntegratedCovariance(fromTime,toTime,j,k);
			}

			fromTime= toTime;
		}

		/*
		DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
		itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

		ARM_GP_Vector eigenValues(totalFwds);		
		for(i=0, ii = offsetIndex;i<timeStepsSize-1 ;++i, ++ii)
		{
			/// initialize everything
			itsEigenValues[i] = new ARM_GP_Vector(factorsNb,0.);
			
			ARM_GP_Matrix* ACPMatrix=ACPTransformation(localVariances[ii],eigenValues);

			int maxRank = -1;
			for(k=0;k<factorsNb;++k)
				if(fabs(eigenValues[k]) < K_NEW_DOUBLE_TOL)
				{
					eigenValues[k] = 0.;
					if (maxRank == -1)
						maxRank = k;
				}
				else
				{
					if (eigenValues[k] < 0.)
					{
	#if defined(__GP_STRICT_VALIDATION)
							ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariancesAndStdDev: negative eigenvalue" );
	#endif
							eigenValues[k] = 0.;
					}
				}
			
			if (maxRank == -1) maxRank = factorsNb;
			size_t effectiveRank = ( maxRank < factorsNb ? maxRank : factorsNb);
			//rescaling
			for (j=0;j<totalFwds;j++)
			{
				double sum=0.;
				for(k=0;k<effectiveRank-1;++k){
					sum+=(*ACPMatrix)(j,k)*(*ACPMatrix)(j,k)*eigenValues[k];
				}
				double sgn = ((*ACPMatrix)(j,k)>0.?1.:-1.);
				double res = (*localVariances[ii])(j,j)-sum;
				if (res>K_NEW_DOUBLE_TOL){
					(*ACPMatrix)(j,k)=sgn*sqrt(res/eigenValues[k]);
				}
			}

			for(k=0;k<factorsNb;++k)
			{
				(*itsEigenValues[i])(k) =	eigenValues[k];
				for(j=0; j<totalFwds; ++j)	{	
					(*localVariances[ii])(j,k) =	(*ACPMatrix)(j,k);
				}
			}

			delete ACPMatrix;
		}

		itsCachedTimeSteps = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(timeSteps.Clone()));
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>(localVariances, itsCachedModelStatesLocalVars);
		*/
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : AdviseBreakPointTimes
////////////////////////////////////////////////////
void ARM_SmiledFRM::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    ARM_GP_Vector  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
	case ARM_ModelParamType::BetaCorrelation:
	    {
			double date = portfolio->GetAsset(0)->GetResetDates()->Elt(0) - asOfDate;
			tmpdates.push_back(date);
			for(i=1; i<size1; i++) 
			{
				double resetlag = portfolio->GetAsset(i)->GetResetDates()->Elt(0) - asOfDate;
				if(fabs (date - resetlag) > FRMVOL_LAG_THRESHOLD)
				{
					tmpdates.push_back(resetlag);
					date = resetlag;
				}
				else
				{
					portfolio->SetWeight(0.0,i);
				}
			}
			modelParam->UpdateValues(&tmpdates);
        }
        break;
		
		/// just return NULL
	case ARM_ModelParamType::Hump:
		break;
	default:
        ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,AdviseBreakPointTimes: only hump or beta correlation" );
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledFRM::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SmiledFRM,TreeStatesToModelStates :  not implemented ARM_SmiledFRM Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_SmiledFRM::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SmiledFRM,VarianceToTime not implemented ARM_SmiledFRM Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: GetEquivalentShift
///	Returns : double
///	Action  : checks if schedule is coherent with model
////////////////////////////////////////////////////

double ARM_SmiledFRM::GetEquivalentShift(const ARM_GP_Vector& weight) const
{
	double shift = 0;
	size_t k = 0;

	for (ARM_GP_Vector::const_iterator iter = weight.begin() ; iter != weight.end() ;++iter,++k)
		shift += (*iter) * itsProcess[k]->GetEqShift();

	return shift;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: ValidateFixSchedule
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_SmiledFRM::GetDiffusionCoeff() const
{
	ARM_GP_Vector* res=new ARM_GP_Vector(itsResetDates.size());

	size_t k = 0;
	
	for (ARM_GP_Vector::iterator iter = res->begin() ; iter != res->end() ;++iter,++k)
		(*iter) = itsProcess[k]->GetDiffusion();

	return ARM_GP_VectorPtr(res);
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SmiledFRM::PreProcessing(ARM_ModelFitter& modelFitter)
{
	itsCurrentArg = NULL;
	((ARM_ModelParamsSmiled* ) GetModelParams())->AdviseSwaptionApprox(itsSwaptionApprox);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SmiledFRM::PostProcessing(const ARM_ModelFitter& modelFitter)
{
    itsCurrentArg	= NULL;
	itsNewArgNeeded = false;
	((ARM_ModelParamsSmiled* ) GetModelParams())->PreComputeIntegratedCovariance();
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: AdviseCurrentCalibSecIndex
///	Returns : void
///	Action  : Advises the model that the current index of the calibration security
///             is the index given ... The model advises just the model params
////////////////////////////////////////////////////
void ARM_SmiledFRM::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	if (itsCurrentArg)
	{
		delete itsCurrentArg;
		itsCurrentArg = NULL;
	}
    itsNewArgNeeded = true;
	ARM_GP_Vector date = ((ARM_ModelParamsSmiled* ) GetModelParams())->GetCorrelCurve()->GetAbscisses();

	if (date.size()>1)
		((ARM_ModelParamsSmiled* ) GetModelParams())->AdviseCurrentCalib(date[index],date[index+1]);
	else
		((ARM_ModelParamsSmiled* ) GetModelParams())->ForOptimize();

}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: GetVanillaSwaptionArg
///	Returns: ARM_VanillaSwaptionArgSmiledFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM* ARM_SmiledFRM::GetVanillaSwaptionArg( 
	const string& curveName,
	double evalTime,
	double swapResetTime,
	double swapNotional,
    const ARM_GP_Vector& fixNotional,
	const ARM_GP_Vector& floatNotional,
	double floatStartTime,
	double floatEndTime,
	const ARM_GP_Vector& floatResetTimes,
	const ARM_GP_Vector& floatStartTimes,
	const ARM_GP_Vector& floatEndTimes,
	const ARM_GP_Vector& floatIntTerms,
	const ARM_GP_Vector& fixPayTimes,
	const ARM_GP_Vector& fixPayPeriods,
    const ARM_GP_Matrix& strikesPerState,
    int callPut,
	const ARM_PricingStatesPtr& states,
	bool isConstantNotional	) const
{
	if( abs(evalTime) > K_DOUBLE_TOL )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,GetVanillaSwaptionArg: evalTime != 0" );

	if( fixNotional.size() != fixPayTimes.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,GetVanillaSwaptionArg: fixNotional.size() != fixPayTimes()" );

	if( fixPayTimes.size() != fixPayPeriods.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,GetVanillaSwaptionArg: fixPayTimes.size() != fixPayPeriods()" );

	ARM_VectorPtr level		= Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);

	ARM_VectorPtr fwd		= SwapRate(	curveName, evalTime,floatStartTime,floatEndTime,
										fixPayTimes,fixPayPeriods, 
										floatStartTimes, floatEndTimes, floatIntTerms, floatEndTimes, floatIntTerms, 
										ARM_GP_Vector(1,0.), false, states);

	ARM_VectorPtr weight	= ComputeDSwapDLib(floatStartTimes, floatEndTimes, fixPayTimes,fixPayPeriods);

	double shift			= GetEquivalentShift( *weight );

	ARM_VectorPtr mu		= ComputeDLogShiftSwapDLib(shift, floatStartTimes, floatEndTimes, fixPayTimes,fixPayPeriods);

	ARM_VanillaSwaptionArgSmiledFRM* arg = new ARM_VanillaSwaptionArgSmiledFRM(shift, ARM_VectorPtr(0), level, fwd);	

	arg->SetMu(mu);
	
	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : computes price approx for standard swaption (using "vol swap/vol fra"-like relationship)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledFRM::VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
	if( !(isConstantNotional && isConstantSpread && isConstantStrike) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,VanillaSwaption: notional, spread, strike must be const" );

	double swapNotional = fixNotional[0];

	size_t nbStates = states->size();
	
	ARM_VanillaSwaptionArgSmiledFRM* arg;
	bool isArgLocal = false;
	if( NULL == itsCurrentArg )
	{
		if (itsNewArgNeeded)
		{
			itsCurrentArg=GetVanillaSwaptionArg( curveName, evalTime, swapResetTime,swapNotional,
				fixNotional, floatNotional, floatStartTime, floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
				fixPayTimes,fixPayPeriods, strikesPerState, callPut, states,isConstantNotional );
			arg = (ARM_VanillaSwaptionArgSmiledFRM*) &*itsCurrentArg;
			itsNewArgNeeded = false;
		}
		else
		{
			arg = GetVanillaSwaptionArg( curveName, evalTime, swapResetTime,swapNotional,
				fixNotional, floatNotional, floatStartTime, floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
				fixPayTimes,fixPayPeriods, strikesPerState, callPut, states,isConstantNotional );
			isArgLocal = true;
		}
	}
	else
	{
		arg = (ARM_VanillaSwaptionArgSmiledFRM*) &*itsCurrentArg;
	}
	
	int approxSwapResetIdx = IdxFromValueWithTol( itsResetDates, swapResetTime, 7. );
	double approxSwapResetDate = itsResetDates[approxSwapResetIdx];

	ARM_VectorPtr coeff		= GetDiffusionCoeff();

	double vol = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(evalTime,approxSwapResetDate,0,0,*arg->GetMu(),*coeff,itsSwaptionApprox);
	vol *= sqrt((swapResetTime-evalTime)/K_YEAR_LEN);

	ARM_VectorPtr level  = arg->GetFixAnnuity();
	ARM_VectorPtr fwd	 = arg->GetSwapFwd();
	double shift		 = arg->GetAverageShift();

	ARM_VectorPtr result(new ARM_GP_Vector(nbStates));
	ARM_GP_Vector::iterator iterFwd = fwd->begin();
	ARM_GP_Vector::iterator iterLevel = level->begin();

	size_t k=0;
	if( vol <= K_NEW_DOUBLE_TOL )
	{
		for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter,++iterFwd,++iterLevel,++k)
			(*result) = swapNotional*(*iterLevel)*CC_Max<double>(callPut*( (*iterFwd) - strikesPerState(k,0)),0.);
	}
	else
	{
		for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter,++iterFwd,++iterLevel,++k)
			(*result) =		swapNotional*BlackSholes_Formula( (*iterFwd) + shift ,vol, (*iterLevel),strikesPerState(k,0) + shift,callPut);
	}

	if (isArgLocal)
		delete arg;
	
	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: ComputeDZCDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

double ARM_SmiledFRM::DZCDLib( double maturity, size_t i) const
{
	if ( maturity < itsStartDates[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,DZCDLib : maturity < first startDate" );

	size_t StartIdx = 0;

	if ( maturity >= itsStartDates[itsStartDates.size()-1] )
		StartIdx = itsStartDates.size()-1;
	else
	{
		if ( itsStartDates[StartIdx] > maturity)
		{
		}
		else 
		{
			while (itsStartDates[StartIdx] <= maturity) { StartIdx++;};
			StartIdx -= 1;
		}
	}

	double res = 0.;

	if ( IsOnSamePath(0,StartIdx) && IsOnSamePath(0,i) )
	{
		if ( i < StartIdx )
			res =	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsProcess[i]->GetDelta() / (1 + itsProcess[i]->GetDelta() * itsProcess[i]->GetFwdRate() ) ;
		else
			if ( StartIdx == i )
			{
				double theta = (maturity - itsStartDates[i])/360.;
				res =	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
						* theta / (1 + theta * itsProcess[i]->GetFwdRate() );
			}
	}
	else
	{
		if ( IsOnSamePath(0,i))
			res = - GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsProcess[i]->GetDelta() / (1 + itsProcess[i]->GetDelta() * itsProcess[i]->GetFwdRate() ) ;
		else 
			if ( IsOnSamePath(StartIdx,i))
				res = GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
					* itsProcess[i]->GetDelta() / (1 + itsProcess[i]->GetDelta() * itsProcess[i]->GetFwdRate() ) ;

		if ( StartIdx == i )
			{
				double theta = (maturity - itsStartDates[i])/360.;
				res +=	- GetZeroCurve()->DiscountPrice(maturity/K_YEAR_LEN)
						* theta / (1 + theta * itsProcess[i]->GetFwdRate() );
			}
	}

	return res;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: ComputeDSwapDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::ComputeDSwapDLib(	const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
					const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const
{
	size_t	i,j,k,
			fwdsNb		= itsResetDates.size(),
			size		= startDates.size(),
			nbFixPay	= fixPayTimes.size();

	double	level		= 0,
			rate		= GetZeroCurve()->DiscountPrice(startDates[0]/K_YEAR_LEN)
						- GetZeroCurve()->DiscountPrice(endDates[size-1]/K_YEAR_LEN);

	for ( k = 0 ; k<nbFixPay ; ++k)
		level += GetZeroCurve()->DiscountPrice(fixPayTimes[k] /K_YEAR_LEN) * fixPayPeriods[k];
	rate /= level;

	ARM_GP_Vector dswapdzc(nbFixPay+2);
	for (i = 0 ; i<nbFixPay ; i++)
		dswapdzc[i]	= - fixPayPeriods[i] * rate / level;

	dswapdzc[i++] = 1. / level;
	dswapdzc[i++] = - 1. / level;

	ARM_GP_Vector* result	=	new ARM_GP_Vector(fwdsNb,0.);
	for (j = 0; j < fwdsNb ; j++)
	{
		for (i = 0 ; i < nbFixPay; i++)
			(*result)[j] += dswapdzc[i]*DZCDLib(fixPayTimes[i],j);

		(*result)[j] += dswapdzc[i++]*DZCDLib(startDates[0],j);
		(*result)[j] += dswapdzc[i++]*DZCDLib(endDates[size-1],j);
	}
	
	
	return ARM_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: ComputeDSwapDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledFRM::ComputeDLogShiftSwapDLib(double shift, const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
					const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const
{
	size_t	k,
			size		= startDates.size(),
			nbFixPay	= fixPayTimes.size();

	double	level		= 0,
			rate		= GetZeroCurve()->DiscountPrice(startDates[0]/K_YEAR_LEN)
						- GetZeroCurve()->DiscountPrice(endDates[size-1]/K_YEAR_LEN);

	for ( k = 0 ; k<nbFixPay ; ++k)
		level += GetZeroCurve()->DiscountPrice(fixPayTimes[k] /K_YEAR_LEN) * fixPayPeriods[k];
	rate /= level;

	ARM_VectorPtr result = ComputeDSwapDLib(startDates,endDates,fixPayTimes,fixPayPeriods);
	for (ARM_GP_Vector::iterator iter = result->begin(); iter != result->end() ; ++iter)
		(*iter) /= ( shift + rate);
	return result;
}
////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledFRM::VanillaDigital(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
        double fwdPeriod,
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,VanillaDigital not implemented" );
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_SmiledFRM::VanillaSpreadOptionLet(const string& curveName,
													double evalTime,
													int callPut,
													double startTime,
													double endTime,
													double resetTime,
													double payTime,
													double payPeriod,
													double notional,
													double coeffLong,
													double coeffShort,
													const ARM_GP_Vector& strikes,
													double swapLongFloatStartTime,
													double swapLongFloatEndTime,
													const ARM_GP_Vector& swapLongFixPayTimes,
													const ARM_GP_Vector& swapLongFixPayPeriods,
													double swapShortFloatStartTime,
													double swapShortFloatEndTime,
													const ARM_GP_Vector& swapShortFixPayTimes,
													const ARM_GP_Vector& swapShortFixPayPeriods,
													const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,VanillaSpreadOptionLet not implemented" );
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Semi-Closed form formula for caplet/floorlet
///	  using gauss-legender
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledFRM::VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const ARM_GP_Vector& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	if ( abs(evalTime-fwdResetTime) < K_DOUBLE_TOL )
    {
		ARM_VectorPtr libor = Libor(curveName, evalTime, fwdStartTime, fwdEndTime, fwdPeriod, fwdResetTime, payTime, states);

        size_t stateSize    = states->size();

        ARM_VectorPtr dfPay = DiscountFactor(curveName, evalTime, payTime, states);

        ARM_GP_VectorPtr result( new ARM_GP_Vector( stateSize ) );

        double payoff;

        for (size_t i(0); i<stateSize; i++)
        {
			payoff = capFloor * ( (*libor)[i] - strikesPerState[i] );

            (*result)[i] = (payoff > 0) ? period * payNotional * (*dfPay)[i] * payoff : 0.0;

		}

		return result;

	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,VanillaCaplet not implemented" );
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
///  Warning : Enhancements are needed!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledFRM::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
	if (fabs(fwdEndTime-payTime)>7 && fabs(evalTime-fwdResetTime)>7)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,Libor : case payTime != fwdEndTime not handled." );

	return DefaultLibor(curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : ARM_SmiledFRM --> num method is always under terminal measuer
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledFRM::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,getTerminalTime()));
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: FirstPricingStates
///	Returns : void
///	Action  : FirstPricingStates
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledFRM::FirstPricingStates( size_t bucketSize ) const
{
	const size_t nbPayoffs			= 0;
	size_t nbModelStates			= ModelStatesSize();
	size_t factorsNb				= FactorCount();

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	//Each model state is initialized to 0.0
	for(size_t i=0; i<bucketSize; ++i )
	{
		for(size_t j=0;j<nbModelStates; ++j)
		{
			initStates->SetModelState(i,j,0.);
		}
	}
	return initStates;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledFRM::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	if (itsCalibrationStatus){
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}else{
		// Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,Init: numerical method not set in SmiledFRM model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,Init: numeraire not set in the SmiledFRM model!");

		/// creates the model schedule (smart pointor for exception safety!)
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps(PricingTimeSteps(timeInfos));

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,ComputeNumeraireTimes( timeInfos ));

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: Induct
///	Returns : ARM_PricingStatesPtr
///	Action  : induct...
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledFRM::Induct(ARM_PricingStatesPtr& states,double toTime)
{

	if (itsCalibrationStatus){
		if (!itsCalibrated){
			ARM_ModelParamsSmiled* params	=	dynamic_cast<ARM_ModelParamsSmiled*>(GetModelParams());
			if( !params )
			   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,Induct : ARM_ModelParamsSmiled* needed" );

			params->setResetDates(ARM_GP_VectorPtr((ARM_GP_Vector*) itsResetDates.Clone()));
			
			size_t nb			=	itsResetDates.size();
			size_t nbDates		=	itsEndDates.size();
			size_t nbTimeSteps	=	itsTimeStepsNb;
			
			ARM_GP_Vector timeSteps(nbTimeSteps);
			double dt			=	itsEndDates[nbDates-1]/K_YEAR_LEN/(nbTimeSteps-1);

			timeSteps[0]		=	0.;
			for	(	size_t i=1 ; i< nbTimeSteps ; ++i)
				timeSteps[i]	=	timeSteps[i-1]+dt;

			ARM_GP_Vector auxDates(itsResetDates.size());
			for (i=0;i<itsResetDates.size();i++)
				auxDates[i]=itsResetDates[i]/K_YEAR_LEN;

			ARM_GP_Vector * newTimeSteps = MergeSortedVectorNoDuplicates( timeSteps, auxDates);
			ARM_VanillaSecDensityPtrVector densityVector=itsNumericalModelFitter->getCalibSecDensities();

			DeletePointorVector<ARM_ProcessBuilder>( itsProcess );
			itsProcess.resize(0);
	
			ARM_ShiftedLNDensityFunctor* pdensity;
			for (size_t k=0;k<nb;k++)
			{
				pdensity = dynamic_cast <ARM_ShiftedLNDensityFunctor*> (&*densityVector[k]->getDensityFunctor());
				if (!(pdensity && itsSkipPDE)){
					ARM_ProcessBuilderPDE* process=new ARM_ProcessBuilderPDE(params->meanReversion(k),(*params->GetVolCurve(k)), itsWithRescalling);
					process->SetPDEParams(itsGridSize,itsStdDevNb,0.5);
					itsProcess.push_back(process);
				}
				else
				{
					itsProcess.push_back(new ARM_ProcessBuilderSLN(params->meanReversion(k),(*params->GetVolCurve(k))));
				}
			}
			
			for (k=0;k<nb;k++)
			{
				itsProcess[k]->Calibrate(itsResetDates[k],(*densityVector[k]),itsResetDates,*newTimeSteps);
				itsProcess[k]->ComputeProxyForCalib(itsSwaptionApprox,(*densityVector[k]),params->GetVolCurveCalib(k));
			}

			itsCalibrated=true;
			delete newTimeSteps;
		}
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}
	
	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );
	return newStates;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: PricingTimeSteps
///	Returns :
///	Action  : returns PricingTimeSteps
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_SmiledFRM::PricingTimeSteps(const ARM_TimeInfoPtrVector& timeInfos)
{
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventAndModelTime();
	ARM_GP_Vector* ptimeSteps = discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this );
	return ptimeSteps;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: computeWeights
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
	
void ARM_SmiledFRM::computeWeights() 
{
	size_t size		= itsResetDates.size();
	double terminal	= itsEndDates[size-1];
	double term;

	itsWeight.resize( size, size);

	for (size_t j=0;j<size;j++)
	{
		size_t k = 0;
		term = itsEndDates[j];
		while (term<terminal)
		{
			while (itsStartDates[k]<term)
			{
				itsWeight(j,k)=0.;
				k++;
			}
			itsWeight(j,k)=1.;
			term = itsEndDates[k];
			k++;
		}
		itsWeight(j,j)=1.;
	}


}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledFRM
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter + reset/start/end dates to DF map
////////////////////////////////////////////////////
	
void ARM_SmiledFRM::setNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{

	if (numericalModelFitter)
	{
		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		size_t size = densityVector.size();

		ARM_GP_Vector resetDates(size); 
		ARM_GP_Vector startDates(size); 
		ARM_GP_Vector endDates(size); 

		for (size_t k = 0 ; k < size ; k++)
		{
			resetDates[k] = densityVector[k]->getResetDate();
			startDates[k] = densityVector[k]->getStartDate();
			endDates[k] = densityVector[k]->getEndDate();
			if (densityVector[k]->getInterestTerms().size()!=1)
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::setNumericalModelFitter : only LIBOR allowed!" );
			
		}
		
		double asof = GetAsOfDate().GetJulian();

		setResetTimes( resetDates - asof);
		setStartTimes( startDates - asof);
		setEndTimes  ( endDates   - asof);

		computeWeights();
		itsCalibrated = false;
	}
	itsNumericalModelFitter = numericalModelFitter;

}




////////////////////////////////////////////////////
///	Class   : ARM_SmiledFRM
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SmiledFRM::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Smiled FRM Model\n";
    os << indent << "----------------------------\n";
	os << ARM_PricingModel::toString(indent);
	os << "\n\n";

	os << indent << "Martingale Processes\n";
    os << indent << "----------------------------\n";
	size_t size=itsProcess.size();
	if (size>0)
		for (size_t k=0;k<size;k++)
			os << itsProcess[k]->toString(indent);

	return os.str();	
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

