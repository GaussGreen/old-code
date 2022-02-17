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
#include "gpbase/stringmanip.h"

/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/SmiledMM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpmodels/ProcessBuilderSmiledFRM.h"
#include "gpmodels/VanillaSwaptionArgSmiledFRM.h"
#include "gpmodels/VanillaSpreadOptionArgSmiledFRM.h"

/// gpcalib
#include "gpcalib/calibmethod.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

/// kernel
#include <inst/swaption.h>




CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SmiledMM::ARM_SmiledMM( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, size_t timeStepsNb,size_t gridSize,double stdDevNb,bool skipPDE, bool allowInterpol, ARM_ModelParamsSmiled::CalibProxy calibProxy,bool cache)
:	ARM_PricingModelIR(zc,params),
	itsNumericalModelFitter(NULL), 
	itsCalibSecDensities(0),
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
	itsCalibProxy(calibProxy),
	itsCache(cache),
	itsRateWeight(0)
{

}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SmiledMM::ARM_SmiledMM(const ARM_SmiledMM& rhs)
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
	itsCalibProxy(rhs.itsCalibProxy),
	itsCache(rhs.itsCache),
	itsRateWeight(rhs.itsRateWeight)
{
	DuplicateCloneablePointorVectorInPlace<ARM_ProcessBuilder>( rhs.itsProcess, itsProcess );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsCalibSecDensities, itsCalibSecDensities);
}	



////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SmiledMM::~ARM_SmiledMM()
{
	DeletePointorVector<ARM_ProcessBuilder>( itsProcess );
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ComputeNumeraireTimeIndex
///	Returns : int
///	Action  : Compute the numeraire time index based on the numeraire
/// current date. 
////////////////////////////////////////////////////
void ARM_SmiledMM::ComputeNumeraireTimeIndex()
{
	double numeraireTime = GetNumeraire()->GetMaturity();	
	itsNumeraireTimeIndex = lower_boundPosWithPrecision(itsEndDates,numeraireTime)-1;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : ARM_SmiledMM --> num method is always under terminal measuer
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledMM::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
// FIXMEFRED: mig.vc8 (30/05/2007 16:33:45):caast
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,getTerminalTime()));
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ComputeModelTimes
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_SmiledMM::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetDates.Clone());
}




////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: PostInit
///	Returns: nothing
///	Action : Function to init after the numerical method
///			 has been initialised.. non const
////////////////////////////////////////////////////
void ARM_SmiledMM::PostInit()
{
	ComputeNumeraireTimeIndex();
}




////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_SmiledMM::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(!params.DoesModelParamExist(ARM_ModelParamType::Hump) || 
       !params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::ValidateModelParams: hump and beta correlation structures are needed for Smiled FRM" );
	}
	
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routine : SetNumeraire
///	Returns : bool
///	Action  : 
////////////////////////////////////////////////////
void ARM_SmiledMM::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	if(	numerairePtr->GetType() != ARM_Numeraire::TerminalZc )
    {		
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::SetNumeraire: only numeraire supported is Terminal ZC" );
    }
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routine : ValidateCalibMethod
///	Returns : bool
///	Action  : validate calib method
////////////////////////////////////////////////////

void ARM_SmiledMM::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);
	if (!itsCalibrated && !isNumerical)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::ValidateCalibMethod: should calib functionals first!!!" );

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
void ARM_SmiledMM::IntegratedLocalDrifts(
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
///	Class   : ARM_SmiledMM
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_SmiledMM::NumMethodStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t totalFwds     = itsResetDates.size(),
		   factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   startTimePos  = 0;
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	for(size_t i=0;i<timeStepsSize-1 ;++i)
	{
		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);
		for (size_t k=0; k<factorsNb;k++)
			(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_SmiledMM::NumMethodStateGlobalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	size_t totalFwds     = itsResetDates.size(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (size_t i=0;i<timeStepsSize;i++){
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds,1.0);
	}
	
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ModelStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_SmiledMM::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
	ARM_MatrixVector modStateLocalVars;
	// The model local vars is now cached
	// very usefull for the importance sampling
	// optimization
	if (!itsCachedTimeSteps.IsNull() && (timeSteps == *(itsCachedTimeSteps)))
	{
		 DuplicateCloneablePointorAndNullVectorInPlace<ARM_GP_Matrix>(itsCachedModelStatesLocalVars,modStateLocalVars);
	}
	else
	{
		ARM_MatrixVector auxLocalCov;
		ModelStateLocalVariances( timeSteps, auxLocalCov);

		size_t totalFwds     = itsResetDates.size(),
			   timeStepsSize = timeSteps.size(),
			   factorsNb     = FactorCount(),
			   modelNb		 = GetModelNb(),
			   startTimePos  = 0,
			   i,j,k;

		double	fromTime= timeSteps[0],
				toTime;

		modStateLocalVars.resize((timeStepsSize-1)*(modelNb+1));

		DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
		itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

		ARM_GP_Vector eigenValues(totalFwds);
		
		for(i=0;i<timeStepsSize-1 ;++i)
		{
			/// initalize the toTime
			toTime  = timeSteps[i+1];

			/// initialize everything
			modStateLocalVars[i] = new ARM_GP_Matrix(totalFwds,factorsNb);
			itsEigenValues[i] = new ARM_GP_Vector(factorsNb);
			
			/// get the first bigger time
			while(startTimePos< itsResetDates.size()
				&& itsResetDates[startTimePos] < toTime )
				++startTimePos;

			ARM_GP_Matrix* ACPMatrix=ACPTransformation(auxLocalCov[i],eigenValues);

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
							ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::ModelStateLocalVariancesAndStdDev: negative eigenvalue" );
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
				double res = (*auxLocalCov[i])(j,j)-sum;
				if (res>K_NEW_DOUBLE_TOL)
					(*ACPMatrix)(j,k)=sgn*sqrt(res/eigenValues[k]);
			}


			for(k=0;k<factorsNb;++k)
			{
				(*itsEigenValues[i])(k) =	eigenValues[k];
				for(j=0; j<totalFwds; ++j)
				{	
					(*modStateLocalVars[i])(j,k) =	(*ACPMatrix)(j,k);
				}
			}

			delete ACPMatrix;
			fromTime= toTime;
		}

		itsCachedTimeSteps = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(timeSteps.Clone()));
		DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>(modStateLocalVars, itsCachedModelStatesLocalVars);
	}

	/// set the result
	SetModelStateLocalVars(modStateLocalVars);
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: NumMethodStateLocalVariancesAndStdDev
///	Returns : void
///	Action  : computes the local variance and std dev 
////////////////////////////////////////////////////
void ARM_SmiledMM::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ModelStateLocalVariances
///	Returns : void
///	Action  : ModelStateLocalVariancess
////////////////////////////////////////////////////
void ARM_SmiledMM::ModelStateLocalVariances( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t totalFwds     = itsResetDates.size(),
		   factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb,
		   startTimePos  = 0;
	double	fromTime= timeSteps[0],
			toTime;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::ModelStateLocalVariances: localDrifts.size() != offsetIndex" );
#endif

    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	size_t i,j,k;


#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::ModelStateLocalVariances: first time step != 0" );
#endif

	for(i=0;i<timeStepsSize-1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(totalFwds,totalFwds);
		
		/// get the first bigger time
		while(startTimePos< itsResetDates.size()
			&& itsResetDates[startTimePos] < toTime )
			++startTimePos;

		if (itsAllowInterpol && startTimePos>0) startTimePos--;	

		for(j=0; j<startTimePos; ++j)
		{
			for(k=j;k<totalFwds;++k)
				(*localVariances[i])(j,k)=(*localVariances[i])(k,j)=0.;
		}
		for(j=startTimePos; j<totalFwds; ++j)
		{
			for(k=j;k<totalFwds;++k)
				(*localVariances[i])(j,k)=(*localVariances[i])(k,j)=((ARM_ModelParamsSmiled*) GetModelParams())->IntegratedCovariance(fromTime,toTime,j,k);
		}

		fromTime= toTime;
	}
}



////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : AdviseBreakPointTimes
////////////////////////////////////////////////////
void ARM_SmiledMM::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
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
        ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::AdviseBreakPointTimes: only hump or beta correlation" );
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledMM::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SmiledMM::TreeStatesToModelStates :  not implemented ARM_SmiledMM Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_SmiledMM::VarianceToTime(double var,double minTime,double maxTime) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT,"ARM_SmiledMM::VarianceToTime not implemented ARM_SmiledMM Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: UnderlyingCovariance
///	Returns : void
///	Action  :  compute the covariance between two ZC forwards
////////////////////////////////////////////////////

double ARM_SmiledMM::UnderlyingCovariance( string	underlyingType, double fromTime, double toTime,
	double startTime1, double endTime1, double startTime2, double endTime2,
	double startTime3, double endTime3, double startTime4, double endTime4) const
{
	if ( stringGetUpper(underlyingType)!="FWD" && stringGetUpper(underlyingType)!= "CMS" )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only FWD or CMS type supported!");

	ARM_GP_Vector VarCovar(3);
	ARM_GP_Vector swapfwd(4);
	ComputeUnderlyingVarCovar( fromTime,toTime,startTime1,endTime1,startTime2,endTime2,VarCovar,
		swapfwd);

	return VarCovar[2];

}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: GetEquivalentShift
///	Returns : double
///	Action  : checks if schedule is coherent with model
////////////////////////////////////////////////////

double ARM_SmiledMM::GetEquivalentShift(const ARM_GP_Vector& weight) const
{
	double shift = 0;
	size_t k = 0;

	for (ARM_GP_Vector::const_iterator iter = weight.begin() ; iter != weight.end() ;++iter,++k)
		shift += (*iter) * itsProcess[k]->GetEqShift();

	return shift;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ValidateFixSchedule
///	Returns : void
///	Action  : 
////////////////////////////////////////////////////

ARM_GP_VectorPtr ARM_SmiledMM::GetDiffusionCoeff() const
{
	ARM_GP_Vector* res=new ARM_GP_Vector(itsResetDates.size());

	size_t k = 0;
	
	for (ARM_GP_Vector::iterator iter = res->begin() ; iter != res->end() ;++iter,++k)
		(*iter) = itsProcess[k]->GetDiffusion();
		//(*iter) = (itsProcess[k]->GetEqShift() + itsProcess[k]->GetFwdRate()) * itsProcess[k]->GetScaling();

	return ARM_GP_VectorPtr(res);
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SmiledMM::PreProcessing(ARM_ModelFitter& modelFitter)
{
	itsCurrentArg = NULL;
	((ARM_ModelParamsSmiled* ) GetModelParams())->AdviseSwaptionApprox(itsCalibProxy);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_SmiledMM::PostProcessing(const ARM_ModelFitter& modelFitter)
{
    itsCurrentArg	= NULL;
	itsNewArgNeeded = false;
	((ARM_ModelParamsSmiled* ) GetModelParams())->PreComputeIntegratedCovariance();
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: AdviseCurrentCalibSecIndex
///	Returns : void
///	Action  : Advises the model that the current index of the calibration security
///             is the index given ... The model advises just the model params
////////////////////////////////////////////////////
void ARM_SmiledMM::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
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
///	Class  : ARM_SmiledMM
///	Routine: ComputeDSwapDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledMM::ComputeDSwapDRate(	const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
					const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const
{
	size_t	i,j,k,
			fwdsNb		= itsResetDates.size(),
			size		= startDates.size(),
			nbFixPay	= fixPayTimes.size();

	ARM_GP_Vector* result	=	new ARM_GP_Vector(fwdsNb,0.);

	if (nbFixPay>0)
	{

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

	
	
		for (j = 0; j < fwdsNb ; j++)
		{
			for (i = 0 ; i < nbFixPay; i++)
				(*result)[j] += dswapdzc[i]*DZCDRate(fixPayTimes[i],j);

			(*result)[j] += dswapdzc[i++]*DZCDRate(startDates[0],j);
			(*result)[j] += dswapdzc[i++]*DZCDRate(endDates[size-1],j);
		}
	}		
	
	return ARM_VectorPtr(result);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: ComputeDSwapDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledMM::ComputeDLogShiftSwapDRate(double shift, const ARM_GP_Vector& startDates,const ARM_GP_Vector& endDates,
					const ARM_GP_Vector& fixPayTimes,const ARM_GP_Vector& fixPayPeriods) const
{
	size_t	k,
			size		= startDates.size(),
			nbFixPay	= fixPayTimes.size();

	ARM_VectorPtr result = ComputeDSwapDRate(startDates,endDates,fixPayTimes,fixPayPeriods);
	
	if (nbFixPay>0)
	{
		double	level		= 0,
				rate		= GetZeroCurve()->DiscountPrice(startDates[0]/K_YEAR_LEN)
							- GetZeroCurve()->DiscountPrice(endDates[size-1]/K_YEAR_LEN);

		for ( k = 0 ; k<nbFixPay ; ++k)
			level += GetZeroCurve()->DiscountPrice(fixPayTimes[k] /K_YEAR_LEN) * fixPayPeriods[k];
		rate /= level;

		for (ARM_GP_Vector::iterator iter = result->begin(); iter != result->end() ; ++iter)
			(*iter) /= ( shift + rate);
	}

	return result;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: ComputeDZCDLib
///	Returns: double
///	Action : 
////////////////////////////////////////////////////

double ARM_SmiledMM::DZCDRate( double maturity, size_t i) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::DZCDRate not implemented" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: GetVanillaSwaptionArg
///	Returns: ARM_VanillaSwaptionArgSmiledFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSmiledFRM* ARM_SmiledMM::GetVanillaSwaptionArg( 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::GetVanillaSwaptionArg: evalTime != 0" );

/*	if( fixNotional.size() != fixPayTimes.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::GetVanillaSwaptionArg: fixNotional.size() != fixPayTimes()" );*/

	if( fixPayTimes.size() != fixPayPeriods.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::GetVanillaSwaptionArg: fixPayTimes.size() != fixPayPeriods()" );

	ARM_VectorPtr level		= Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);

	ARM_VectorPtr fwd		= SwapRate(	curveName, evalTime,floatStartTime,floatEndTime,
										fixPayTimes,fixPayPeriods, 
										floatStartTimes, floatEndTimes, floatIntTerms, floatEndTimes, floatIntTerms, 
										ARM_GP_Vector(1,0.), false, states);

	ARM_VectorPtr weight	= ComputeDSwapDRate(floatStartTimes, floatEndTimes, fixPayTimes,fixPayPeriods);

	itsRateWeight = weight;

	double shift			= GetEquivalentShift( *weight );

	ARM_VectorPtr mu		= ComputeDLogShiftSwapDRate(shift, floatStartTimes, floatEndTimes, fixPayTimes,fixPayPeriods);

	ARM_VanillaSwaptionArgSmiledFRM* arg = new ARM_VanillaSwaptionArgSmiledFRM(shift, ARM_VectorPtr(0), level, fwd);	

	arg->SetMu(mu);
	
	return arg;
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : computes price approx for standard swaption (using "vol swap/vol fra"-like relationship)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledMM::VanillaSwaption(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::VanillaSwaption: notional, spread, strike must be const" );

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
	
	
	ARM_VectorPtr level  = arg->GetFixAnnuity();
	ARM_VectorPtr fwd	 = arg->GetSwapFwd();
	double shift		 = arg->GetAverageShift();

	int approxSwapResetIdx = IdxFromValueWithTol( itsResetDates, swapResetTime, 7. );
	double approxSwapResetDate = itsResetDates[approxSwapResetIdx];

	ARM_VectorPtr result(new ARM_GP_Vector(nbStates));
	ARM_GP_Vector::iterator iterFwd = fwd->begin();
	ARM_GP_Vector::iterator iterLevel = level->begin();
	
	if (itsCalibProxy==ARM_ModelParamsSmiled::GaussBasketMoneyness || itsCalibProxy==ARM_ModelParamsSmiled::GaussBasketAtm)
	{
		for (size_t i=0;i<itsResetDates.size();i++)
			itsProcess[i]->ComputeProxyForCalib(itsCalibProxy,(*itsCalibSecDensities[i]),0);

		ARM_VectorPtr coeff		= GetDiffusionCoeff();
		double vol = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(evalTime,approxSwapResetDate,0,0,*arg->GetMu(),*coeff,itsCalibProxy);
		vol *= sqrt((swapResetTime-evalTime)/K_YEAR_LEN);
		vol *= (*fwd)[0];

		if (itsCalibProxy==ARM_ModelParamsSmiled::GaussBasketMoneyness)
		{
			double moneyness = (strikesPerState(0,0)-(*fwd)[0])/vol;
			for (size_t i=0;i<itsResetDates.size();i++)
			{
				itsProcess[i]->ComputeProxyForCalib(itsCalibProxy,(*itsCalibSecDensities[i]),0,moneyness);
			}
			coeff = GetDiffusionCoeff();
			vol = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(evalTime,approxSwapResetDate,0,0,*arg->GetMu(),*coeff,itsCalibProxy);
			vol *= sqrt((swapResetTime-evalTime)/K_YEAR_LEN);
			vol *= (*fwd)[0];
		}

		size_t k=0;
		if( vol <= K_NEW_DOUBLE_TOL )
		{
			for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter,++iterFwd,++iterLevel,++k)
				(*result) = swapNotional*(*iterLevel)*CC_Max<double>(callPut*( (*iterFwd) - strikesPerState(k,0)),0.);
		}
		else
		{
			for (ARM_GP_Vector::iterator iter = result->begin() ; iter != result->end() ; ++iter,++iterFwd,++iterLevel,++k)
				(*result) =		swapNotional*(*iterLevel)*VanillaOption_N((*iterFwd),vol,strikesPerState(k,0),1.,callPut);
		}

	}
	else
	{
		ARM_VectorPtr coeff		= GetDiffusionCoeff();
		double vol = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(evalTime,approxSwapResetDate,0,0,*arg->GetMu(),*coeff,itsCalibProxy);
		vol *= sqrt((swapResetTime-evalTime)/K_YEAR_LEN);

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
	}

	if (isArgLocal)
		delete arg;
	
	return result;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: UnderlyingCorrelation
///	Returns : void
///	Action  :  compute the correlation between two underlying
////////////////////////////////////////////////////

double ARM_SmiledMM::UnderlyingCorrelation(	string	underlyingType,
	double fromTime, double toTime, double	startTime1, double endTime1,
	double startTime2, double endTime2, double startTime3, double endTime3, double startTime4, double endTime4) const
{

	if ( stringGetUpper(underlyingType)!="FWD" && stringGetUpper(underlyingType)!= "CMS" )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only FWD or CMS type supported!");

	ARM_GP_Vector VarCovar(3);
	ARM_GP_Vector swapfwd(4);
	
	ComputeUnderlyingVarCovar( fromTime,toTime,startTime1,endTime1,startTime2,endTime2,VarCovar,
		swapfwd);
	double stdev1 = VarCovar[0];
	double stdev2 = VarCovar[1];
	double covariance_1_2 = VarCovar[2]; 
	double correlation_1_2 = covariance_1_2/(stdev1*stdev2);
	return  correlation_1_2 ;	
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: ComputeUnderlyingVarCovar
///	Returns : void
///	Action  :  compute the correlation between two underlying
////////////////////////////////////////////////////

void ARM_SmiledMM::ComputeUnderlyingVarCovar( double fromTime, double toTime, 
	double startTime1, double endTime1, double startTime2,
	double endTime2, ARM_GP_Vector& VarCovar, ARM_GP_Vector& swapfwd) const
{
	//ARM_ModelParams* modParamsSFRM = ((ARM_ModelParamsSmiled*)GetModelParams() );

	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	string curveCcy				= Ccy->GetCcyName();
	int fixFreq					= Ccy->GetFixedPayFreq();
	long fixDayCount			= Ccy->GetFixedDayCount();
	char fixCalendar[100];
	Ccy->CalcFixPayCal(fixCalendar);
	int fwdDayCount				= Ccy->GetLiborIndexDayCount();

	double evalTime				= 0.0;
	double swapNotional			= 100.0;
	int callPut					= 1;
	ARM_GP_Matrix strikesPerState(1,0.0);
	size_t factorsNb			= FactorCount();
	ARM_PricingStatesPtr states = ARM_PricingStatesPtr( new ARM_PricingStates(1,1,1,factorsNb) );

	//CMS1
	ARM_Date startDate1(asOfDate+startTime1);
	ARM_Date endDate1(asOfDate+endTime1);
	ARM_DateStrip fixDateStrip1( startDate1, endDate1, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );

	ARM_GP_Vector* fixPayTimes1		= fixDateStrip1.GetPaymentDates();
	ARM_GP_Vector* fixPayPeriods1	= fixDateStrip1.GetInterestTerms();
	ARM_GP_Vector* floatStartTimes1	= fixDateStrip1.GetFlowStartDates();
	ARM_GP_Vector* floatPayTimes1	= fixDateStrip1.GetFlowEndDates();
	ARM_GP_Vector* floatIntTerms1	= fixDateStrip1.GetInterestTerms();
	ARM_GP_Vector* floatResetDate1	= fixDateStrip1.GetResetDates();

	double swapResetTime1	= (*floatResetDate1)[0]-asOfDate;	
	int size1				= floatResetDate1->size();
	int i,nbFixFlows1		= fixPayTimes1->size();    
	
	for(i=0;i<nbFixFlows1;++i)
		(*fixPayTimes1)[i] = (*fixPayTimes1)[i]-asOfDate;

	ARM_GP_Vector fixNotional1 (nbFixFlows1,1.0);
	ARM_GP_Vector floatNotional1 (nbFixFlows1,1.0);
	
	for(i=0;i<size1;++i)
	{
		(*floatResetDate1)[i] = (*floatResetDate1)[i]-asOfDate;
		(*floatStartTimes1)[i] = (*floatStartTimes1)[i]-asOfDate;
		(*floatPayTimes1)[i] = (*floatPayTimes1)[i]-asOfDate;
	}
	
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArg1;
	swaptionArg1 = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime1,
		swapNotional,fixNotional1,floatNotional1, startTime1, endTime1,*floatResetDate1,*floatStartTimes1,*floatPayTimes1,*floatIntTerms1, *fixPayTimes1, *fixPayPeriods1, strikesPerState,
		callPut, states);

	ARM_VectorPtr coeff		= GetDiffusionCoeff();
	ARM_GP_VectorPtr mu1 = swaptionArg1->GetMu();
	double var1 = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(fromTime,toTime,0,0,*mu1,*coeff,itsCalibProxy);
	var1 = var1 * sqrt((toTime-fromTime)/K_YEAR_LEN);

	ARM_VectorPtr swapfwd1	= swaptionArg1->GetSwapFwd();
	double shift1 = swaptionArg1->GetAverageShift();

	//CMS2		
	ARM_Date startDate2(asOfDate+startTime2);
	ARM_Date endDate2(asOfDate+endTime2);
	ARM_DateStrip fixDateStrip2( startDate2, endDate2, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );

	ARM_GP_Vector* fixPayTimes2		= fixDateStrip2.GetPaymentDates();
	ARM_GP_Vector* fixPayPeriods2	= fixDateStrip2.GetInterestTerms();
	ARM_GP_Vector* floatStartTimes2	= fixDateStrip2.GetFlowStartDates();
	ARM_GP_Vector* floatPayTimes2	= fixDateStrip2.GetFlowEndDates();
	ARM_GP_Vector* floatIntTerms2	= fixDateStrip2.GetInterestTerms();
	ARM_GP_Vector* floatResetDate2	= fixDateStrip2.GetResetDates();

	double swapResetTime2 = (*floatResetDate2)[0]-asOfDate;
	int size2 = floatResetDate2->size();
	int nbFixFlows2=fixPayTimes2->size();    
	
	for(i=0;i<nbFixFlows2;++i)
		(*fixPayTimes2)[i] = (*fixPayTimes2)[i]-asOfDate;
	
	for(i=0;i<size2;++i)
	{
		(*floatResetDate2)[i] = (*floatResetDate2)[i]-asOfDate;
		(*floatStartTimes2)[i] = (*floatStartTimes2)[i]-asOfDate;
		(*floatPayTimes2)[i] = (*floatPayTimes2)[i]-asOfDate;
	}
	ARM_GP_Vector fixNotional2 (nbFixFlows2,1.0);
	ARM_GP_Vector floatNotional2 (nbFixFlows2,1.0);
	
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArg2;
	swaptionArg2 = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime2,swapNotional,fixNotional2,floatNotional2,
		startTime2,endTime2,*floatResetDate2,*floatStartTimes2,*floatPayTimes2,*floatIntTerms2,*fixPayTimes2,*fixPayPeriods2,strikesPerState,
		callPut, states);

	ARM_GP_VectorPtr mu2 = swaptionArg2->GetMu();
	double var2 = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(fromTime,toTime,0,0,*mu2,*coeff,itsCalibProxy);
	var2 = var2 * sqrt((toTime-fromTime)/K_YEAR_LEN);
	
	ARM_VectorPtr swapfwd2	= swaptionArg2->GetSwapFwd();
	double shift2 = swaptionArg2->GetAverageShift();
	
	double covariance_1_2 =((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentCov(fromTime,toTime,*mu1,*mu2,*coeff,itsCalibProxy);
	
	VarCovar[0] = var1;
	VarCovar[1] = var2;
	VarCovar[2] = covariance_1_2;

	swapfwd.resize(4);
	swapfwd[0]=(*swapfwd1)[0];
	swapfwd[1]=(*swapfwd2)[0];
	swapfwd[2]=shift1;
	swapfwd[3]=shift2;

	delete swaptionArg1;
	delete swaptionArg2;
}
////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledMM::VanillaDigital(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::VanillaDigital not implemented" );
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_SmiledMM::VanillaSpreadOptionLet(const string& curveName,
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
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (30/05/2007 16:34:04):cast
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));
	
	size_t nbStates = states->size();
	
	ARM_VanillaSpreadOptionArgSmiledFRM* arg;
	bool isArgLocal = false;
	if( NULL == itsCurrentArg )
	{
		if (itsNewArgNeeded)
		{
			itsCurrentArg=GetVanillaSpreadOptionArg(	curveName,evalTime,callPut,startTime,endTime,resetTime,payTime,
														payPeriod,notional,coeffLong,coeffShort,strikes,swapLongFloatStartTime,
														swapLongFloatEndTime,swapLongFixPayTimes,swapLongFixPayPeriods,
														swapShortFloatStartTime,swapShortFloatEndTime,swapShortFixPayTimes,swapShortFixPayPeriods);
			arg = (ARM_VanillaSpreadOptionArgSmiledFRM*) &*itsCurrentArg;
			itsNewArgNeeded = false;
		}
		else
		{
			arg = GetVanillaSpreadOptionArg(	curveName,evalTime,callPut,startTime,endTime,resetTime,payTime,
												payPeriod,notional,coeffLong,coeffShort,strikes,swapLongFloatStartTime,
												swapLongFloatEndTime,swapLongFixPayTimes,swapLongFixPayPeriods,
												swapShortFloatStartTime,swapShortFloatEndTime,swapShortFixPayTimes,swapShortFixPayPeriods);
			isArgLocal = true;
		}
	}
	else
	{
		arg = (ARM_VanillaSpreadOptionArgSmiledFRM*) &*itsCurrentArg;
	}

	
    ARM_VanillaSwaptionArgSmiledFRM* swaptionArg_Long				= ((ARM_VanillaSpreadOptionArgSmiledFRM*)arg)->GetSwaptionArg_Long();
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArg_Short				= ((ARM_VanillaSpreadOptionArgSmiledFRM*)arg)->GetSwaptionArg_Short();
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArg_StartTimePayTime	= ((ARM_VanillaSpreadOptionArgSmiledFRM*)arg)->GetSwaptionArg_StartTimePayTime();

	//Pricing
	double asOfDate		= GetAsOfDate().GetJulian();
	double fromTime		= evalTime;
	double toTime		= resetTime;
	double startTime1	= swapShortFloatStartTime;
	double endTime1		= swapShortFloatEndTime;
	double startTime2	= swapLongFloatStartTime;
	double endTime2		= swapLongFloatEndTime;
	
	double   coeff1  = coeffShort;
	double   coeff2  = coeffLong;

	ARM_VectorPtr coeff		= GetDiffusionCoeff();

	//long Swaption	
	ARM_GP_VectorPtr mu_Long = swaptionArg_Long->GetMu();
	double vol_Long = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(fromTime,toTime,0,0,*swaptionArg_Long->GetMu(),*coeff,itsCalibProxy);
	vol_Long = vol_Long*sqrt((toTime-fromTime)/K_YEAR_LEN);
	ARM_VectorPtr swapfwd_Long	= swaptionArg_Long->GetSwapFwd();
	double shift_Long = swaptionArg_Long->GetAverageShift();
	

	//shortSwaption
	ARM_GP_VectorPtr mu_Short = swaptionArg_Short->GetMu();
	double vol_Short = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(fromTime,toTime,0,0,*swaptionArg_Short->GetMu(),*coeff,itsCalibProxy);
	vol_Short = vol_Short*sqrt((toTime-fromTime)/K_YEAR_LEN);
	ARM_VectorPtr swapfwd_Short	= swaptionArg_Short->GetSwapFwd();
	double shift_Short = swaptionArg_Short->GetAverageShift();

	//StartTimePayTime Swaption
	ARM_GP_VectorPtr mu_StartTimePayTime = swaptionArg_StartTimePayTime->GetMu();
	double vol_StartTimePayTime = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentVol(fromTime,toTime,0,0,*swaptionArg_StartTimePayTime->GetMu(),*coeff,itsCalibProxy);
	vol_StartTimePayTime = vol_StartTimePayTime*sqrt((toTime-fromTime)/K_YEAR_LEN);
	ARM_VectorPtr swapfwd_StartTimePayTime	= swaptionArg_StartTimePayTime->GetSwapFwd();
	double shift_StartTimePayTime = swaptionArg_StartTimePayTime->GetAverageShift();

	
	//covariance_long_Short	&& covariance_StartTimePayTime_Long
	double covariance_Long_Short			 = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentCov(fromTime,toTime,*mu_Long,*mu_Short,*coeff,itsCalibProxy);
	double covariance_StartTimePayTime_Long  = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentCov(fromTime,toTime,*mu_Long,*mu_StartTimePayTime,*coeff,itsCalibProxy);
	double covariance_StartTimePayTime_Short = ((ARM_ModelParamsSmiled*) GetModelParams())->GetEquivalentCov(fromTime,toTime,*mu_StartTimePayTime,*mu_Short,*coeff,itsCalibProxy);

	double theta_R = (toTime-evalTime)/K_YEAR_LEN;
	double vol1 = vol_Short/sqrt(theta_R);
	double vol2 = vol_Long/sqrt(theta_R);
	double fwd1 = (*swapfwd_Short)[0];
	double fwd2 = (*swapfwd_Long)[0];
	double shift1 = shift_Short;
	double shift2 = shift_Long;
	double vol1ForP = vol1*(fwd1+shift1)/fwd1;
	double vol2ForP = vol2*(fwd2+shift2)/fwd2;

	//BS convexity Adjustment (As ARM&Summit)
	
	double payLagConv_1 = 1.;
	double payLagConv_2 = 1.;

	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	
	int decompfreq = Ccy->GetFixedPayFreq();
	int Tenor_1 = static_cast<int>(floor((endTime1-startTime1)/K_YEAR_LEN));
	int Tenor_2 = static_cast<int>(floor((endTime2-startTime2)/K_YEAR_LEN));

	double natLagConv_1 = 1-Tenor_1*fwd1/((1+fwd1/decompfreq)*(pow(1+fwd1/decompfreq,Tenor_1*decompfreq)-1));
	double natLagConv_2 = 1-Tenor_2*fwd2/((1+fwd2/decompfreq)*(pow(1+fwd2/decompfreq,Tenor_2*decompfreq)-1));
	natLagConv_1 *=(fwd1+shift1)/fwd1*(fwd1+shift1)/fwd1*(exp(vol1*vol1*theta_R)-1);
	natLagConv_2 *=(fwd2+shift2)/fwd2*(fwd2+shift2)/fwd2*(exp(vol2*vol2*theta_R)-1);
	natLagConv_1 +=1;
	natLagConv_2 +=1;

	double AdjFwd1 = (fwd1)*natLagConv_1;
	double AdjFwd2 = (fwd2)*natLagConv_2;		

	if (fabs(covariance_StartTimePayTime_Long)>K_DOUBLE_TOL)
	{
		double theta_P = (payTime-startTime)/360.0;
		double Vol_P = vol_StartTimePayTime/sqrt(theta_R);
		double swap_P =(*swapfwd_StartTimePayTime)[0];
		double shift_P = shift_StartTimePayTime;
		Vol_P=Vol_P*(shift_P+swap_P)/swap_P;

		double Correl_1 = covariance_StartTimePayTime_Short/(vol_StartTimePayTime*vol_Short);
		double Correl_2 = covariance_StartTimePayTime_Long/(vol_StartTimePayTime*vol_Long);

		double payLagConv = sqrt(exp(pow(theta_P*Vol_P*swap_P*sqrt(theta_R)/(1+theta_P*swap_P),2))-1);
		payLagConv_1 =1-Correl_1*payLagConv*sqrt(exp(vol1ForP*vol1ForP*theta_R)-1);
		payLagConv_2 =1-Correl_2*payLagConv*sqrt(exp(vol2ForP*vol2ForP*theta_R)-1);
		/*payLagConv_1 = exp(-Correl_1*vol1ForP*Vol_P*theta_P*swap_P/(1+theta_P*swap_P)*theta_R);
		payLagConv_2 = exp(-Correl_2*vol2ForP*Vol_P*theta_P*swap_P/(1+theta_P*swap_P)*theta_R);*/
		/*payLagConv_1 = 1.+theta_P*swap_P/(1+theta_P*swap_P)*(1.-(1+theta_P*swap_P/(1+theta_P*swap_P)*(exp(Vol_P*Vol_P*theta_R)-1))*exp(Correl_1*Vol_P*vol1ForP*theta_R));
		payLagConv_2 = 1.+theta_P*swap_P/(1+theta_P*swap_P)*(1.-(1+theta_P*swap_P/(1+theta_P*swap_P)*(exp(Vol_P*Vol_P*theta_R)-1))*exp(Correl_2*Vol_P*vol2ForP*theta_R));*/
	}

	
	AdjFwd1 *= payLagConv_1;
	AdjFwd2 *= payLagConv_2;		
	
	double Correl = covariance_Long_Short/(vol_Short*vol_Long);
	double optMat = toTime/K_YEAR_LEN;
	ARM_VectorPtr spreadOption(new ARM_GP_Vector(nbStates));
	size_t i;

	double zcpay	=  GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);

	for(i=0;i<nbStates;i++)
		(*spreadOption)[i]=notional*zcpay*payPeriod*SpreadOption(coeff1*(AdjFwd1+shift1), coeff2*(AdjFwd2+shift2), 
			vol1, vol2,0.0, 0.0, Correl, 0.0,strikes[i]+coeff2*shift2-coeff1*shift1, optMat, callPut);
	
	if (isArgLocal)
		delete arg;

	return spreadOption;
}

////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: GetVanillaSpreadOptionletArg
///	Returns: ARM_VanillaSwaptionArgSFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArgSmiledFRM* ARM_SmiledMM::GetVanillaSpreadOptionArg( 
		const string& curveName,
		double evalTime,
		int	callPut,
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
		const ARM_GP_Vector& swapShortFixPayPeriods) const
{
	if( abs(evalTime) > K_DOUBLE_TOL )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::GetVanillaSpreadOptionArg: evalTime != 0" );

	double asOfDate				= GetAsOfDate().GetJulian();
	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	string curveCcy				= Ccy->GetCcyName();
	
	int fixFreq					= Ccy->GetFixedPayFreq();
	long fixDayCount			= Ccy->GetFixedDayCount();
	
	char fixCalendar[100];
	Ccy->CalcFixPayCal(fixCalendar);
	int fwdDayCount				= Ccy->GetLiborIndexDayCount();

//CMSLong (Fix Leg)
	ARM_Date startDate1(asOfDate+swapLongFloatStartTime);
	ARM_Date endDate1(asOfDate+swapLongFloatEndTime);
	ARM_DateStrip fixDateStrip1( startDate1, endDate1, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,fixCalendar );

	ARM_GP_Vector* longFixPayTimes	 = fixDateStrip1.GetPaymentDates();
	ARM_GP_Vector* longFixPayPeriods = fixDateStrip1.GetInterestTerms(); 
	ARM_GP_Vector* floatStartTimes1	= fixDateStrip1.GetFlowStartDates();
	ARM_GP_Vector* floatPayTimes1	= fixDateStrip1.GetFlowEndDates();
	ARM_GP_Vector* floatIntTerms1	= fixDateStrip1.GetInterestTerms();
	ARM_GP_Vector* floatResetDate1	= fixDateStrip1.GetResetDates();

	int i,nbFixFlowsLong=longFixPayTimes->size();    
	
	for(i=0;i<nbFixFlowsLong;++i)
	{
		(*longFixPayTimes)[i]	= (*longFixPayTimes)[i]-asOfDate;
		(*floatResetDate1)[i] = (*floatResetDate1)[i]-asOfDate;
		(*floatStartTimes1)[i] = (*floatStartTimes1)[i]-asOfDate;
		(*floatPayTimes1)[i] = (*floatPayTimes1)[i]-asOfDate;
	}

//CMSShort (Fix leg)
	ARM_Date startDate2(asOfDate+swapShortFloatStartTime);
	ARM_Date endDate2(asOfDate+swapShortFloatEndTime);
	ARM_DateStrip fixDateStrip2( startDate2, endDate2, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,fixCalendar );

	ARM_GP_Vector* shortFixPayTimes   = fixDateStrip2.GetPaymentDates();
	ARM_GP_Vector* shortFixPayPeriods = fixDateStrip2.GetInterestTerms();
	ARM_GP_Vector* floatStartTimes2	= fixDateStrip2.GetFlowStartDates();
	ARM_GP_Vector* floatPayTimes2	= fixDateStrip2.GetFlowEndDates();
	ARM_GP_Vector* floatIntTerms2	= fixDateStrip2.GetInterestTerms();
	ARM_GP_Vector* floatResetDate2	= fixDateStrip2.GetResetDates();

	int nbFixFlowsShort=shortFixPayTimes->size();    
	
	for(i=0;i<nbFixFlowsShort;++i)
	{
		(*shortFixPayTimes)[i]	= (*shortFixPayTimes)[i]-asOfDate;
		(*floatResetDate2)[i] = (*floatResetDate2)[i]-asOfDate;
		(*floatStartTimes2)[i] = (*floatStartTimes2)[i]-asOfDate;
		(*floatPayTimes2)[i] = (*floatPayTimes2)[i]-asOfDate;
	}

//Create ARM_VanillaSpreadOptionArg
	ARM_GP_Vector* NullStrikes = new ARM_GP_Vector(strikes);
	ARM_GP_Vector* resetTimesVector					= new ARM_GP_Vector(1,resetTime);
	ARM_GP_Vector* payTimesVector					= new ARM_GP_Vector(1,payTime);
	ARM_GP_Vector* payPeriodsVector					= new ARM_GP_Vector(1,payPeriod);
	ARM_GP_Vector* notionalVector					= new ARM_GP_Vector(1,notional);
	ARM_GP_Vector* coeffLongVector					= new ARM_GP_Vector(1,coeffLong);
	ARM_GP_Vector* coeffShortVector					= new ARM_GP_Vector(1,coeffShort);
	ARM_GP_Vector* swapLongFloatStartTimeVector		= new ARM_GP_Vector(1,swapLongFloatStartTime);
	ARM_GP_Vector* swapLongFloatEndTimeVector		= new ARM_GP_Vector(1,swapLongFloatEndTime);
	ARM_GP_Vector* swapShortFloatStartTimeVector	= new ARM_GP_Vector(1,swapShortFloatStartTime);
	ARM_GP_Vector* swapShortFloatEndTimeVector		= new ARM_GP_Vector(1,swapShortFloatEndTime);	

	ARM_VectorVector swapLongFixPayTimesVector;
	ARM_VectorVector swapLongFixPayPeriodsVector;
	ARM_VectorVector swapShortFixPayTimesVector;
	ARM_VectorVector swapShortFixPayPeriodsVector;
	
	swapLongFixPayTimesVector.push_back(longFixPayTimes);
	swapLongFixPayPeriodsVector.push_back(longFixPayPeriods);
	swapShortFixPayTimesVector.push_back(shortFixPayTimes);
	swapShortFixPayPeriodsVector.push_back(shortFixPayPeriods);
	
	
	ARM_VanillaSpreadOptionArgSmiledFRM* arg = new ARM_VanillaSpreadOptionArgSmiledFRM();
	
	//create the swaptions arg	
	double swapNotional			= 100.0;
	ARM_GP_Vector fixNotional (1,0.0);
	ARM_GP_Vector floatNotional (1,0.0);
	ARM_GP_Matrix strikesPerState(1,0.0);
	size_t factorsNb			= FactorCount();
	ARM_PricingStatesPtr states = ARM_PricingStatesPtr( new ARM_PricingStates(1,1,1,factorsNb) );

	double swapResetTime1 = (*floatResetDate1)[0];
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArgLong;
	swaptionArgLong = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime1,
		swapNotional,fixNotional,floatNotional, swapLongFloatStartTime, swapLongFloatEndTime,*floatResetDate1,*floatStartTimes1,*floatPayTimes1,*floatIntTerms1, *longFixPayTimes, *longFixPayPeriods, strikesPerState,
		callPut, states);
	arg->SetSwaptionArg_Long(*swaptionArgLong);

	
	double swapResetTime2 = (*floatResetDate2)[0];	
		
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArgShort;
	swaptionArgShort = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTime2,
		swapNotional,fixNotional,floatNotional, swapShortFloatStartTime, swapShortFloatEndTime, *floatResetDate2,*floatStartTimes2,*floatPayTimes2,*floatIntTerms2,*shortFixPayTimes, *shortFixPayPeriods, strikesPerState,
		callPut, states);
	arg->SetSwaptionArg_Short(*swaptionArgShort);


	//itsSwaptionArg_ToTimePayTime
	ARM_Date startDateToTime(asOfDate+resetTime);
	ARM_Date endDateToTime(asOfDate+payTime);
	ARM_DateStrip fixDateStripToTime( startDateToTime, endDateToTime, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );
	ARM_GP_Vector* fixPayTimesToTime	= fixDateStripToTime.GetPaymentDates();
	ARM_GP_Vector* fixPayPeriodsToTime	= fixDateStripToTime.GetInterestTerms();
	ARM_GP_Vector* floatResetDateToTime	= fixDateStripToTime.GetResetDates();
	ARM_GP_Vector* floatStartDateToTime	= fixDateStripToTime.GetFlowStartDates();
	ARM_GP_Vector* floatEndDateToTime	= fixDateStripToTime.GetFlowEndDates();
	ARM_GP_Vector* floatIntTermsToTime	= fixDateStripToTime.GetInterestTerms();

	int nbFixFlowsToTime=fixPayTimesToTime->size();    
	
	for(i=0;i<nbFixFlowsToTime;++i)
	{
		(*fixPayTimesToTime)[i] = (*fixPayTimesToTime)[i]-asOfDate;
		(*floatResetDateToTime)[i] = (*floatResetDateToTime)[i]-asOfDate;
		(*floatStartDateToTime)[i] = CC_Max<double>((*floatStartDateToTime)[i]-asOfDate,itsStartDates[0]);
		(*floatEndDateToTime)[i] = (*floatEndDateToTime)[i]-asOfDate;
	}

	
	double swapResetTimeToTime = (*floatResetDateToTime)[0];	
		
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArgToTime;
	swaptionArgToTime = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTimeToTime,
		swapNotional,fixNotional,floatNotional, resetTime, payTime, *floatResetDateToTime,*floatStartDateToTime,*floatEndDateToTime,*floatIntTermsToTime,*fixPayTimesToTime, *fixPayPeriodsToTime, strikesPerState,
		callPut, states);

	arg->SetSwaptionArg_ToTimePayTime(*swaptionArgToTime);

	//itsSwaptionArg_StartTimePayTime
	ARM_Date startDateStartTime(asOfDate+swapLongFloatStartTime);
	ARM_Date endDateStartTime(asOfDate+payTime);
	ARM_DateStrip fixDateStripStartTime( startDateStartTime, endDateStartTime, 1, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, 1, GETDEFAULTVALUE,
		fixCalendar );
	ARM_GP_Vector* fixPayTimesStartTime		= fixDateStripStartTime.GetPaymentDates();
	ARM_GP_Vector* fixPayPeriodsStartTime	= fixDateStripStartTime.GetInterestTerms();
	ARM_GP_Vector* floatResetDateStartTime	= fixDateStripStartTime.GetResetDates();
	ARM_GP_Vector* floatStartDateStartTime	= fixDateStripStartTime.GetFlowStartDates();
	ARM_GP_Vector* floatEndDateStartTime	= fixDateStripStartTime.GetFlowEndDates();
	ARM_GP_Vector* floatIntTermsStartToTime	= fixDateStripStartTime.GetInterestTerms();


	int nbFixFlowsStartTime=fixPayTimesStartTime->size();    
	
	for(i=0;i<nbFixFlowsStartTime;++i)
	{
		(*fixPayTimesStartTime)[i] = (*fixPayTimesStartTime)[i]-asOfDate;
		(*floatResetDateStartTime)[i] = (*floatResetDateStartTime)[i]-asOfDate;
		(*floatStartDateStartTime)[i] = CC_Max<double>((*floatStartDateStartTime)[i]-asOfDate,itsStartDates[0]);
		(*floatEndDateStartTime)[i] = (*floatEndDateStartTime)[i]-asOfDate;
	}

	double swapResetTimeStartTime = (*floatResetDateStartTime)[0];	
	
	ARM_VanillaSwaptionArgSmiledFRM* swaptionArgStartTime;
	swaptionArgStartTime = GetVanillaSwaptionArg(curveCcy, evalTime, swapResetTimeStartTime,
		swapNotional,fixNotional,floatNotional, swapLongFloatStartTime, payTime, *floatResetDateStartTime,*floatStartDateStartTime,*floatEndDateStartTime,*floatIntTermsStartToTime,*fixPayTimesStartTime, *fixPayPeriodsStartTime, strikesPerState,
		callPut, states);

	arg->SetSwaptionArg_StartTimePayTime(*swaptionArgStartTime);

	return arg;
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: DiscountFactor
///	Returns : ARM_VectorPtr
///	Action  : 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_SmiledMM::DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::DiscountFactor not implemented" );
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
///  Warning : Enhancements are needed!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledMM::Libor( 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Libor : case payTime != fwdEndTime not handled." );

	return DefaultLibor(curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Semi-Closed form formula for caplet/floorlet
///	  using gauss-legender
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledMM::VanillaCaplet(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::VanillaCaplet not implemented" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SmiledMM
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledMM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::MCModelStatesFromToNextTime not implemented" );
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: FirstPricingStates
///	Returns : void
///	Action  : FirstPricingStates
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledMM::FirstPricingStates( size_t bucketSize ) const
{
	const size_t nbPayoffs			= 0;
	size_t nbModelStates			= (itsCache?2.*itsResetDates.size():itsResetDates.size());
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
///	Class  : ARM_SmiledMM
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledMM::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	if (itsCalibrationStatus){
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}else{
		// Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Init: numerical method not set in SmiledFRM model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Init: numeraire not set in the SmiledFRM model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventAndModelTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,ARM_GP_VectorPtr( new ARM_GP_Vector( 1, getTerminalTime() )));

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: Induct
///	Returns : ARM_PricingStatesPtr
///	Action  : induct...
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SmiledMM::Induct(ARM_PricingStatesPtr& states,double toTime)
{

	if (itsCalibrationStatus){
		if (!itsCalibrated){
			ARM_ModelParamsSmiled* params	=	dynamic_cast<ARM_ModelParamsSmiled*>(GetModelParams());
			if( !params )
			   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Induct : ARM_ModelParamsSmiled* needed" );

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
					ARM_ProcessBuilderPDE* process=new ARM_ProcessBuilderPDE(params->meanReversion(k),(*params->GetVolCurve(k)));
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
				itsProcess[k]->ComputeProxyForCalib(itsCalibProxy,(*densityVector[k]),params->GetVolCurveCalib(k));
			}

			computeDWeights();
			itsCalibrated=true;
			delete newTimeSteps;
		}
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}
	
	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );
	return newStates;
}









////////////////////////////////////////////////////
///	Class   : ARM_SmiledMM
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SmiledMM::toString(const string& indent, const string& nextIndent) const
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

	os << indent << "Rate weight\n";
    os << indent << "----------------------------\n";
	if (itsRateWeight != ARM_VectorPtr(NULL))
	{
		for (size_t i=0;i<itsRateWeight->size();i++)
		{
			os << itsRateWeight->Elt(i) << "\n";
		}
	}

	return os.str();	
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

