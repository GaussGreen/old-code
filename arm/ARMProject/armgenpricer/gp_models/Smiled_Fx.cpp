/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Smiled_Fx.cpp
 *
 *  \brief Smiled FX model
 *
 *	\author  R. Guillemot
 *	\version 1.0
 *	\date June 2006
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"

/// gpmodel
#include "gpmodels/Smiled_Fx.h"

// gpbase
#include "gpbase/vectormanip.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/interpolatorvector.h"

/// gpinfra
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelnamemap.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillasecuritydensity.h"

/// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/2IRFXModel.h"

#include "gpnummethods/treebase.h"

/// to look at the aleas only
#include <iostream>
#include <fstream>
#include <ios>
#include <iomanip>
using namespace std;

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SmiledModel_Fx::ARM_SmiledModel_Fx(const ARM_ZeroCurvePtr& zc, 
					 ARM_ModelParamsSmiled_Fx* modelParam,
					 size_t timeStepsNb,
					 size_t gridSize,
					 double stdDevNb,
					 bool skipPDE,
					 bool withRescalling,
					 const ARM_2IRFXModel* Model2IRFX)
:	ARM_EqFxBase(zc,modelParam),
	itsNumericalModelFitter(NULL), 
	itsEigenValues(NULL),
	itsResetDates(NULL),
	itsSettlementDates(NULL),
	itsProcess(0),
	itsCalibrationStatus(true),
	itsIsCalibrated(false),
	itsSkipPDE(skipPDE),
	itsTimeStepsNb(timeStepsNb),
	itsGridSize(gridSize),
	itsStdDevNb(stdDevNb),
	itsWithRescalling(withRescalling),
	itsIRModel(NULL),
	its2IRFXModel(CreateClonedPtr(const_cast<ARM_2IRFXModel*>(Model2IRFX))),
	itsTotalVar(0)
{
	if (!zc.IsNull() && modelParam)
		ARM_EqFxBase::Init();

	if(	!dynamic_cast<const ARM_ModelParamsSmiled_Fx*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only Smiled Fx model param" );
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: Copy Constructor
///	Returns: 
///	Action : Copy Constructor
////////////////////////////////////////////////////


ARM_SmiledModel_Fx::ARM_SmiledModel_Fx( const ARM_SmiledModel_Fx& rhs )
:	ARM_EqFxBase(rhs),
itsNumericalModelFitter(rhs.itsNumericalModelFitter),
itsResetDates(rhs.itsResetDates),
itsSettlementDates(rhs.itsSettlementDates),
itsProcess(rhs.itsProcess),
itsCalibrationStatus(rhs.itsCalibrationStatus),
itsIsCalibrated(rhs.itsIsCalibrated),
itsSkipPDE(rhs.itsSkipPDE),
itsTimeStepsNb(rhs.itsTimeStepsNb),
itsGridSize(rhs.itsGridSize),
itsStdDevNb(rhs.itsStdDevNb),
itsWithRescalling(rhs.itsWithRescalling),
itsIRModel(rhs.itsIRModel),
its2IRFXModel(rhs.its2IRFXModel),
itsTotalVar(rhs.itsTotalVar)
{
	DuplicateCloneablePointorVectorInPlace<ARM_ProcessBuilder>( rhs.itsProcess, itsProcess );
	DuplicateCloneablePointorVectorInPlace<std::vector<double>>( rhs.itsEigenValues, itsEigenValues );
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////

ARM_SmiledModel_Fx::~ARM_SmiledModel_Fx()
{
	DeletePointorVector<ARM_ProcessBuilder>( itsProcess );
	DeletePointorVector<std::vector<double>>( itsEigenValues );
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx

///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledModel_Fx::Forward(
	const string& curveName, 
	double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	ARM_GP_MatrixPtr modelStates = states->GetModelStates();
	ARM_VectorPtr fwdStates;

	const ARM_ModelParamsSmiled_Fx* modelParams = dynamic_cast<const ARM_ModelParamsSmiled_Fx*>(&*GetModelParams());

	if (!modelParams)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Smiled FX Model Params is not of good type" );

	double forwardValuefwdTime	= ComputeFwdAtTime( settlementTime );
	
	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,forwardValuefwdTime) );
	}
	else
	{
		if( DoesResetDateExist( evalTime ) )
		{
			size_t stateSize = states->size();
			ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new std::vector<double>(stateSize, 0.0 ));

			size_t ResetIdx = IdxFromValueWithTol( itsResetDates, expiryTime, 1. );

			if (modelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
			{
				int k;

				if( (GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDLOOKING)
				||(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_FWDBCKWDLOOKING) )
				{
					double fwdValue = ComputeFwdAtTime(evalTime);

					for (k = 0; k < stateSize; ++k)
						(*modelStates)(modelNb+1, k) = log((*modelStates)(modelNb+1, k)/fwdValue)+0.5*itsTotalVar[ResetIdx];

					fwdVector = itsProcess[ResetIdx]->Rate(evalTime,states,modelNb+1,0,false);
			
					for (k = 0; k < stateSize; ++k)
						(*modelStates)(modelNb+1, k) = fwdValue*exp(-0.5*itsTotalVar[ResetIdx]+(*modelStates)(modelNb+1, k));				
				} 
				else 
				{
					double lagdrift = 0.0;
					if (GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
					{
						double numTime = GetNumeraire()->GetMaturity();
						lagdrift = its2IRFXModel->ComputeAbsTimeLag(expiryTime,numTime);
					}

					for (k = 0; k < stateSize; ++k)
						(*modelStates)(modelNb+1, k) = (*modelStates)(modelNb+1, k)+lagdrift;

					fwdVector = itsProcess[ResetIdx]->Rate(evalTime,states,modelNb+1,0,false);

					for (k = 0; k < stateSize; ++k)
						(*modelStates)(modelNb+1, k) = (*modelStates)(modelNb+1, k)-lagdrift;
				}
			}
			else
				fwdVector = itsProcess[ResetIdx]->Rate(evalTime,states,modelNb,ResetIdx,false);

			return fwdVector;
		}
		else
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << "the reset time " << evalTime <<" does'nt exist in the Smile FX Mode";
			ARM_THROW( ERR_INVALID_ARGUMENT, os.str() );
		}
	}
	
};

////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledModel_Fx::CallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikePerState,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_SmiledModel_Fx::CallVectorial is not implemented");

	ARM_GP_VectorPtr result(NULL);

	return result;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SmiledModel_Fx::DigitalVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikePerState,
	double notional,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_SmiledModel_Fx::VolatilitiesAndCorrelations is not implemented");
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: FirstPricingStates
///	Returns : void
///	Action  : FirstPricingStates
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_SmiledModel_Fx::FirstPricingStates( size_t bucketSize ) const
{
	const size_t nbPayoffs			= 0;
	
	const ARM_ModelParamsSmiled_Fx* modelParams = dynamic_cast<const ARM_ModelParamsSmiled_Fx*>(&*GetModelParams());
	ARM_PricingStatesPtr initStates;

	if (modelParams->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
	{
		size_t nbModelStates			= 2;
		size_t factorsNb				= 2;
		initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );
		double initValue = ComputeFwdAtTime(0.0);
		for(size_t i=0; i<bucketSize; ++i )
		{
			initStates->SetModelState(i,0,0.0);
			initStates->SetModelState(i,1,initValue);
		}

	}
	else
	{
		size_t nbModelStates			= itsResetDates.size();
		size_t factorsNb				= FactorCount();

		initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

		//Each model state is initialized to 0.0
		for(size_t i=0; i<bucketSize; ++i )
		{
			for(size_t j=0;j<nbModelStates; ++j)
			{
				initStates->SetModelState(i,j,0.);
			}
		}
	}

	return initStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_SmiledModel_Fx::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
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
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventAndModelTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,ARM_GP_VectorPtr( new std::vector<double>( 1, getTerminalTime() )));

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: Induct
///	Returns : ARM_PricingStatesPtr
///	Action  : induct...
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_SmiledModel_Fx::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	if (itsCalibrationStatus){
		if (!itsIsCalibrated){
			ARM_ModelParamsSmiled* params	=	dynamic_cast<ARM_ModelParamsSmiled*>(GetModelParams());
			if( !params )
			   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::Induct : ARM_ModelParamsSmiled* needed" );

			ARM_GP_Matrix volatilities;
			ARM_MatrixVector correlMatrix;

			if (params->GetCorrelType() == ARM_ModelParamsSmiled::CorrelMatrix)
			{
				if (its2IRFXModel.IsNull())
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::Induct : The CorrelType=CorrelMatrix needs a 2IRFX Model to calculate correlation." );

				std::vector<double> totalVolatilites;
				its2IRFXModel->ComputeVolatilitiesAndCorrelMatrix(itsResetDates,itsSettlementDates,volatilities,correlMatrix,totalVolatilites,false);
				params->setVolatilitiesAndCorrelMatrix(volatilities,correlMatrix);
			}
			else if (params->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
			{ 
				if (its2IRFXModel.IsNull())
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::Induct : The CorrelType=CorrelMatrix needs a 2IRFX Model to calculate correlation." );

				its2IRFXModel->ComputeVolatilities(itsResetDates,itsSettlementDates,volatilities,itsTotalVar);
				params->setVolatilities(volatilities); //.transpose());
			}
				
			params->setResetDates(ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*> itsResetDates.Clone()));
			
			size_t nb			=	itsResetDates.size();
			size_t nbDates		=	itsResetDates.size();
			size_t nbTimeSteps	=	itsTimeStepsNb;
			
			std::vector<double> timeSteps(nbTimeSteps);
			double dt			=	itsResetDates[nb-1]/K_YEAR_LEN/(nbTimeSteps-1);

			timeSteps[0]		=	0.;
			for	(	size_t i=1 ; i< nbTimeSteps ; ++i)
				timeSteps[i]	=	timeSteps[i-1]+dt;

			std::vector<double> auxDates(itsResetDates.size());
			for (i=0;i<itsResetDates.size();i++)
				auxDates[i]=itsResetDates[i]/K_YEAR_LEN;

			std::vector<double> * newTimeSteps = MergeSortedVectorNoDuplicates( timeSteps, auxDates);
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
				itsProcess[k]->Calibrate(itsResetDates[k],(*densityVector[k]),itsResetDates,*newTimeSteps,!itsSkipPDE);
			}

			delete newTimeSteps;
			itsIsCalibrated = true;
		}

		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}
	
	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );
	return newStates;
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter
////////////////////////////////////////////////////
	
void ARM_SmiledModel_Fx::setNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{
	if (numericalModelFitter)
	{
		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		ARM_VanillaSecurityDensityFX* vsec;
		size_t size = densityVector.size();

		std::vector<double> resetDates(size); 
		std::vector<double> settlementDates(size);

		for (size_t k = 0 ; k < size ; k++)
		{
			vsec =  dynamic_cast<ARM_VanillaSecurityDensityFX*>(&*densityVector[k]);
			if (vsec)
			{
				resetDates[k] = vsec->getResetDate();
				settlementDates[k] = vsec->getSettlementDate();
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_Smiled_Fx::setNumericalModelFitter: calibration is just possible with vanilla security density." );
		}
		
		double asof = GetAsOfDate().GetJulian();

		setResetTimes( resetDates - asof);
		setSettlementTimes( settlementDates - asof);
	}
	itsNumericalModelFitter = numericalModelFitter;

}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t factorsNb     = FactorCount(),
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
///	Class   : ARM_SmiledModel_Fx
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	size_t factorsNb     = FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (size_t i=0;i<timeStepsSize;i++){
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(factorsNb,1.0);
	}
	
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: ModelStateLocalVariances
///	Returns : void
///	Action  : ModelStateLocalVariancess
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::ModelStateLocalVariances( const std::vector<double>& timeSteps,
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	DeletePointorVector<std::vector<double>>( itsEigenValues );
	itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

	size_t i,j,k;

	// Variance/Covariance matrix before ACP
	ARM_GP_Matrix auxLocalVariances(totalFwds,totalFwds);

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,ModelStateLocalVariances: first time step != 0" );
#endif

	for(i=0;i<timeStepsSize-1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(totalFwds,factorsNb);
		itsEigenValues[i] = new std::vector<double>(factorsNb,0.);
		
		/// get the first bigger time
		while(startTimePos< itsResetDates.size()
			&& itsResetDates[startTimePos] < toTime )
			++startTimePos;

		for(j=0; j<startTimePos; ++j)
		{
			for(k=j;k<totalFwds;++k)
				auxLocalVariances(j,k)=auxLocalVariances(k,j)=0.;
		}
		for(j=startTimePos; j<totalFwds; ++j)
		{
			for(k=j;k<totalFwds;++k)
				auxLocalVariances(j,k)=auxLocalVariances(k,j)=((ARM_ModelParamsSmiled*) GetModelParams())->IntegratedCovariance(fromTime,toTime,j,k);
		}

		ACPTransformationWithRescalling(&auxLocalVariances,*itsEigenValues[i],*localVariances[i]);

		fromTime= toTime;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
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
///	Class  : ARM_SmiledModel_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SmiledModel_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double currTime  = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime  = GetNumMethod()->GetTimeStep(timeIndex+1);

	const ARM_ModelParamsSmiled* params	=	dynamic_cast< const ARM_ModelParamsSmiled*>(GetModelParams());
	if( !params )
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM::Induct : ARM_ModelParamsSmiled* needed" );

	if (params->GetCorrelType() != ARM_ModelParamsSmiled::Fwd)
	{
		size_t factorsNb = FactorCount();
		size_t statesNb  = states->size();
		size_t modelNb	 = GetModelNb();

		size_t i,j,k;
		size_t fwdsNb  = itsResetDates.size();
		size_t iLast   = fwdsNb-1;

		int iFirst  = 0;
		while( iFirst < fwdsNb 	&& itsResetDates[iFirst]<nextTime)
		{
			++iFirst;
		}
		
		const ARM_MatrixVector& modelLocalVar	= GetModelStateLocalVars();
		double next,modStd;

		std::vector<double>& eigen = itsEigenValues[timeIndex];
		ARM_GP_Matrix E(fwdsNb,factorsNb,0.);
		
		for(k = 0 ; k < statesNb ; ++k )
		{
			for (i=iFirst;i<=iLast;++i)
			{
				next		= states->GetModelState(k,i+modelNb);
				for (j=0; j<factorsNb; ++j )
				{
					modStd	=	(*modelLocalVar[timeIndex])(i,j);
					next	+=	modStd * states->GetNumMethodState(k,j+modelNb);
				}
				states->SetModelState(k,i+modelNb,next);
			}
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: ComputeModelTimes
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
std::vector<double>& ARM_SmiledModel_Fx::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetDates.Clone());
}


////////////////////////////////////////////////////
///	Class  : ARM_SmiledModel_Fx
///	Routine: UpdateLinks
///	Returns: void
///	Action : maintains links with itsIRModel
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();

	string IRModelName, ModelName = GetModelName();

	/// Find the name of its IRModel
	if( itsIRModel.IsNull() )
	{
		/// If there is no IRModel, its name is found through its Other Models in the ModelNameMap
		ARM_IntVector otherModels = (*modelMap)[ModelName]->OtherModelRefNb();

		if( otherModels.size() != 1 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Wrong number of otherModels" );

		IRModelName = (*modelMap)[otherModels[0]]->ModelName();
	}
	else
		IRModelName = itsIRModel->GetModelName();

	/// The updated IRModel is set. 
	itsIRModel = (*modelMap)[IRModelName]->Model();
}


////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: ModelStatesSize
///	Returns :
///	Action  : Return the size of the model
////////////////////////////////////////////////////

size_t ARM_SmiledModel_Fx::ModelStatesSize() const
{
	const ARM_ModelParamsSmiled* params	=	dynamic_cast< const ARM_ModelParamsSmiled*>(GetModelParams());

	if (params->GetCorrelType() == ARM_ModelParamsSmiled::Fwd)
		return FactorCount();
	else
		return itsResetDates.size();
}

////////////////////////////////////////////////////
///	Class   : ARM_SmiledModel_Fx
///	Routines: ModelStatesSize
///	Returns :
///	Action  : Volatilities & correlation in the trre
/// Nabyl
////////////////////////////////////////////////////
void ARM_SmiledModel_Fx::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
														ARM_GP_MatrixPtr& vols,
														ARM_GP_MatrixPtr& d1Vols, 
														ARM_GP_MatrixPtr& correls,
														bool linearVol) const
{
	//ARM_THROW( ERR_INVALID_ARGUMENT, " unimplemented function for ARM_SmiledModel_Fx model!");
	
	if (!linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_SmiledModel_Fx::VolatilitiesAndCorrelations : only linearVol case is implemented");
	
	/// to speed up, compute in one go the vol and vol derivatives
	std::vector<double> times		= GetModelParams()->GetModelParam(ARM_ModelParamType::Hump).ToCurveModelParam().GetCurve()->GetAbscisses();
	std::vector<double> values	= GetModelParams()->GetModelParam(ARM_ModelParamType::Hump).ToCurveModelParam().GetCurve()->GetOrdinates();

	/// First step : shift values to be equivalent as a constant left interpolation
	size_t nbT ( times.size() );
	std::vector<double> newTimes;
	std::vector<double> newValues;
	if (times[0]>0.)
	{
		newTimes.push_back(0.);
		newValues.push_back(values[0]);
	}

	size_t i;
	for(i=0; i<nbT; ++i)
	{
		newTimes.push_back( times[i] );
		newValues.push_back( values[i+1 < nbT ? i+1 : nbT-1] );
	}

	/// Second step : use special interpolation to compute vol values and 1st derivatives
	nbT = timeSteps.size();
	std::vector<double> volsVec(nbT),d1VolsVec(nbT);

	for(i=0; i<nbT; ++i)
		volsVec[i]	 = FunctionSpecialInterpolation(timeSteps[i], newTimes, newValues);
	
	for(i=0; i<nbT; ++i)
		d1VolsVec[i] = DerivativeSpecialInterpolation(timeSteps[i], timeSteps, volsVec)*K_YEAR_LEN;

	/// Factor by line and zero correl because in one factor!
	vols = ARM_GP_MatrixPtr		( new ARM_GP_Matrix(1, volsVec.size(),	 &volsVec[0]) );
	d1Vols = ARM_GP_MatrixPtr	( new ARM_GP_Matrix(1, d1VolsVec.size(), &d1VolsVec[0]) );
	correls	= ARM_GP_MatrixPtr( NULL );
}

void ARM_SmiledModel_Fx::EulerLocalDrifts(const std::vector<double>& timeSteps, ARM_GP_MatrixPtr& relativeDrifts,
										  ARM_GP_MatrixPtr& absoluteDrifts) const
{
	// Nabyl
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(), 1, -GetModelParams()->GetModelParam( 
											ARM_ModelParamType::BetaCorrelation).GetValueAtPoint(0) ) );

	for(size_t i=0; i<timeSteps.size()-1; ++i)
		(*relativeDrifts)(i,0) = (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;

	absoluteDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0. ) );
}

////////////////////////////////////////////////////	
///	Class   : ARM_SmiledModel_Fx
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_SmiledModel_Fx::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Smiled Fx Model\n";
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

