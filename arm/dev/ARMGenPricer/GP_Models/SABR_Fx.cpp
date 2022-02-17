/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file SABR_Fx.cpp
 *
 *  \brief SABR FX version
 *
 *	\author  E. Ezzine
 *	\version 1.0
 *	\date May 2006
 */

/// to remove identified warning
#include "gpbase/removeidentifiedwarning.h"

/// gpmodel
#include "gpmodels/SABR_Fx.h"


/// gpinfra
#include "gpinfra/pricingmodeltype.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricingcontext.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillafxoption.h"

/// gpmodels
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/SABR_ModelParams.h"
#include "gpmodels/multiassets.h"

/// gpnummethods
#include "gpnummethods/tree1D.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/meanrevertingsampler.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/sabrvanilla.h"


/// to look at the aleas only
#include <iostream>
#include <fstream>
#include <ios>
#include <iomanip>
using namespace std;

CC_BEGIN_NAMESPACE( ARM )


///  define for code clarity
#if defined(FIRST_STATE_VARIABLE)
	#undef	FIRST_STATE_VARIABLE
#endif
#define FIRST_STATE_VARIABLE 0

////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SABRModel_Fx::ARM_SABRModel_Fx(const ARM_ZeroCurvePtr& zc, 
					 ARM_ModelParamsSABR_Fx* modelParam)
:	ARM_EqFxBase(zc,modelParam,ARM_GP_Matrix(NULL),new ARM_SABRDensityFunctor()),
	itsIntegrationStep(140),
	itsTypeOfModel(ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT)
{
	if (!zc.IsNull() && modelParam)
		ARM_EqFxBase::Init();

	if(	!dynamic_cast<const ARM_ModelParamsSABR_Fx*>( modelParam ) )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": accepted model params are only SABR Fx model param" );
}

////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: Forward
///	Returns: a vector of forward (t,T)
///	Action : computes the forward of the equity or fx model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SABRModel_Fx::Forward(
	const string& curveName, 
	double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
	const ARM_PricingStatesPtr& states) const
{
	double forward	= ComputeFwdAtTime( settlementTime );
	
	if(		evalTime<=K_NEW_DOUBLE_TOL
		 || states == ARM_PricingStatesPtr(NULL) )
	{
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,forward) );
	}
	else
	{
		size_t stateSize = states->size();
		ARM_GP_VectorPtr fwdVector = ARM_GP_VectorPtr( new ARM_GP_Vector(stateSize, 0.0 ));
		double stateLocalVar = static_cast< const ARM_SABR_ModelParams* const>(GetModelParams())->StateLocalVariance(0.0,evalTime);
		double C = 0.5*stateLocalVar;

		for( size_t i=0; i<stateSize; ++i )
			(*fwdVector )[i] = forward * exp( states->GetModelState(i,FIRST_STATE_VARIABLE)-C );
		return fwdVector;
	}
};

////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: CallVectorial
///	Returns: a vector of call vectorial
///	Action : computes equity or fx call
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SABRModel_Fx::CallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikePerState,
	int callPut,
	double payTime,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
    /// FIX FIX : review formula with payTime !!
	ARM_GP_VectorPtr FwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, states );
	ARM_VectorPtr DF		= DiscountFactor( modelName, evalTime, payTime, states );

#if defined(__GP_STRICT_VALIDATION)
	if( strikePerState.size() != FwdFx->size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": strikePerState.size() != forward.size()" );
#endif

	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	size_t nbStates = states!=ARM_PricingStatesPtr(NULL)? states->size(): 1;

	double time		= expiryTime-evalTime;
	double tenor	= 0.0;
	int TypeOfModel = ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT;

	size_t size = FwdFx->size();
	ARM_GP_VectorPtr result( new ARM_GP_Vector(size));
	for( size_t i=0; i<size; ++i )
	{
		double impliedVol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol((*FwdFx)[i],
			strikePerState[i],time,tenor,TypeOfModel,140.0/*itsIntegrationStep */);
		
		(*result)[i] = BlackSholes_Formula((*FwdFx)[i],impliedVol*sqrt(time/K_YEAR_LEN),(*DF)[i],
			strikePerState[i],callPut);
	}


	return result;
}

///////////////////////////////////////////////////
///	Class  : ARM_SABR_Model
///	Routine: ImpliedVol
///	Returns: double
///	Action : To Calculate the Implied Volatility
///  By defaut using BS Formula
////////////////////////////////////////////////////
double ARM_SABRModel_Fx::ImpliedVol(const ARM_VanillaArg& arg) const
{
    double fxFwd,expiryTime,strike;
	const ARM_VanillaFxOption* fxOption = dynamic_cast< const ARM_VanillaFxOption*>(&arg);
    if(fxOption)
	{
		string modelName = arg.GetCurveName();
		double evalTime    = arg.GetEvalTime();
		expiryTime  = arg.GetExpiry();
		int callPut = arg.GetCallPut();

		double settlementTime = ((ARM_VanillaEqFxOption&)arg).GetFwdTime();
		strike = ((ARM_VanillaEqFxOption&)arg).GetStrike();
		double payTime = ((ARM_VanillaEqFxOption&)arg).GetPayTime();
		
		ARM_GP_VectorPtr FwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, ARM_PricingStatesPtr(NULL) );
		
		fxFwd =(*FwdFx)[0];
	}
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, "ImpliedVol: Only Fx Option is supported" );

	double tenor = 0.0;
    double value = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol(fxFwd,
                    strike,expiryTime,tenor,itsTypeOfModel,itsIntegrationStep);

    return value;
}

////////////////////////////////////////////////////
///	Class   : ARM_SABRModel_Fx
///	Routine : ComputeDerivatives
///	Returns : double
///	Action  : computes the derivatives of the SABR model
////////////////////////////////////////////////////
double ARM_SABRModel_Fx::PartialDerivative( const ARM_ModelParam& modelParam, 
       size_t number, 
       size_t factorNb,
       const ARM_VanillaArg& arg ,
       ARM_MktTargetType targetFuncType)
{
    double value;	
	const ARM_VanillaFxOption* fxOption = dynamic_cast< const ARM_VanillaFxOption*>(&arg);
	if(fxOption)
	{
		string modelName = arg.GetCurveName();
		double evalTime    = arg.GetEvalTime();
		double expiryTime  = arg.GetExpiry();
		int callPut = arg.GetCallPut();

		double settlementTime = ((ARM_VanillaEqFxOption&)arg).GetFwdTime();
		double strike = ((ARM_VanillaEqFxOption&)arg).GetStrike();
		double payTime = ((ARM_VanillaEqFxOption&)arg).GetPayTime();
		
		ARM_GP_VectorPtr FwdFx = Forward(modelName, evalTime, expiryTime, settlementTime, payTime, ARM_PricingStatesPtr(NULL) );
		
		double fxFwd =(*FwdFx)[0];  
		
		double tenor = 0.0;
		value = ((ARM_SABR_ModelParams*)GetModelParams())->PartialDerivative(fxFwd,
                        strike,evalTime,tenor,itsTypeOfModel,itsIntegrationStep, modelParam.GetType());

        if(targetFuncType == ARM_CalibrationTarget::PriceTarget)
		{
			ARM_VectorPtr DF = DiscountFactor( modelName,
                    evalTime, 
                    payTime,
                    ARM_PricingStatesPtr(NULL)  );
	
                double impliedvol = ((ARM_SABR_ModelParams*)GetModelParams())->ImpliedVol(fxFwd,
                                       strike,expiryTime,tenor,itsTypeOfModel,itsIntegrationStep);

                double vega = BlackSholes_Derivative_2(fxFwd,
                    impliedvol*sqrt(expiryTime/K_YEAR_LEN),
                    (*DF)[0],
                    strike,
                    callPut);
                value*=vega;
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "PartialDerivative: Only Fx Option is supported"  );

    return value;
}

////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: DigitalVectorial
///	Returns: a vector of digital vectorial
///	Action : computes equity or fx digital
////////////////////////////////////////////////////
ARM_VectorPtr ARM_SABRModel_Fx::DigitalVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const ARM_GP_Vector& strikePerState,
	double notional,
	int callPut,
	double payTime,
	ARM_DigitType digitType,
	double epsilon,
    const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_SABRModel_Fx::VolatilitiesAndCorrelations is not implemented");
}

////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: UpdateDensityFunctor
///	Returns: void
///	Action : Update the density functor at expiryTime		
////////////////////////////////////////////////////

void ARM_SABRModel_Fx::UpdateDensityFunctor(double forward, double expiryTime) 
{
	ARM_SABRDensityFunctor* densityfunctor = dynamic_cast<ARM_SABRDensityFunctor*>(&*GetDensityFunctor());
	if(!densityfunctor)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : the density functor must be of the type ARM_SABRDensityFunctor");
	double tenor = 0;
	double alpha = GetModelParams()->GetModelParam(ARM_ModelParamType::Alpha).GetValue(expiryTime,tenor);
	double rho = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(expiryTime,tenor);
	double nu = GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol).GetValue(expiryTime,tenor);
	double beta = GetModelParams()->GetModelParam(ARM_ModelParamType::Beta).GetValue(expiryTime,tenor);
	int sabrType = ARM_CF_SABR_ImplicitVol_Formula::DIRECTEXACT;

	densityfunctor->SetAlpha( alpha );
	densityfunctor->SetBeta( beta );
	densityfunctor->SetRho( rho );
	densityfunctor->SetNu( nu );
	densityfunctor->SetSabrType( sabrType );
}


////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_SABRModel_Fx::ModelStateLocalVariances(
	const ARM_GP_Vector& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelRank= GetModelRank();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelRank;
	localVariances.resize((nbSteps-1)*(modelRank+1));
	size_t i,j;

	int factorNb = FactorCount();

	ARM_GP_Matrix identity(factorNb,factorNb);

	for (i = 0; i < factorNb; ++i)
	{
		for(j = 0; j < i; ++j)
		{
			identity(i,j) = identity(j,i) = 0.0;
		}
		identity(i,i) = 1.0;
	}

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = static_cast<ARM_GP_Matrix*>(identity.Clone());
		
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_SABRModel_Fx::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	 /// test the numeraire and its type!
	ARM_NumerairePtr numeraire=GetNumeraire();
    if( numeraire == ARM_NumerairePtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the FxEQuity model!");

	//// Initialise the numeraire
	numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);
	size_t modelNb = GetModelNb();

    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in FxEQuity model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<ARM_GP_Vector> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );
    
		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		/// ...initialise it
		return numMethod->Init(*this,firstInductTime);
    }
    else
    {
        // Compute a single model states set to (0.0,...,0.0)
        int nbDir = FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_SABRModel_Fx::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	/// FIX FIX: Olivier compléter par l'algorithme d'Euler
	/// il faut corréler les BM
	/// les BM sont rangés dans l'order (i,nb) pour trajectoire i et BM nb

	/// les states comprennent les model states et les num method states
	/// il faut récupérer les BM des num method states et peupler les model states
	/// les valeurs model précédentes sont donc dans les model states
	

#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	double nextTime		= GetNumMethod()->GetTimeStep(timeIndex+1);
	double currentTime	= GetNumMethod()->GetTimeStep(timeIndex);
	double dt			= (nextTime-currentTime)/K_YEAR_LEN;
	double sqrDt		= sqrt(dt);

	size_t factorsNb	= FactorCount();
	size_t modelNb		= GetModelNb();
	size_t statesNb		= states->size();
	double rho			= GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation ).GetValue(currentTime);
	double beta			= GetModelParams()->GetModelParam( ARM_ModelParamType::Beta ).GetValue(currentTime);
	double voldrift		= GetModelParams()->GetModelParam( ARM_ModelParamType::VolDrift ).GetValue(currentTime);
	double volOfVol		= GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).GetValue(currentTime);
	double divPart		= IntegrateStepWise( *GetModelParams()->GetModelParam(ARM_ModelParamType::Dividend).ToCurveModelParam().GetCurve(), currentTime/K_YEAR_LEN, nextTime/K_YEAR_LEN );
	double expDiv		= exp(-divPart);

	ARM_ZeroCurvePtr ZcCurve 
						= GetZeroCurve();
	double zcT			= ZcCurve->DiscountPrice(currentTime/K_YEAR_LEN);
	double zcTPlus		= ZcCurve->DiscountPrice(nextTime/K_YEAR_LEN);
	double irTerm		= zcT/zcTPlus;
	double rhoOrt		= sqrt(1.0-rho*rho);
	double volOfVolSq	= volOfVol*volOfVol*dt;

	double f,sigma,W1,W2,fBetaMinus1,volTerm,drift;
	size_t i;

//// Checks if the model works with an IRModel
	const ARM_PricingModelIR* itsRefModel = getIRModel();
	bool worksWithIRModel = (itsRefModel ? true : false );

	if( worksWithIRModel )  /// the correlation between brownians is managed externally (we do not use the rho parameters)
	{
		/// We are supposed to be in fwd neutral probability there (checked in firstpricing states)... 
		ARM_GP_VectorPtr DFs = DiscountFactor( GetModelName(), currentTime,nextTime,states);
		double probaDate = GetNumeraire()->GetMaturity();
		double FwdDrift = itsRefModel->IntegratedBondCovariance( currentTime, nextTime, nextTime, probaDate );
		FwdDrift -= itsRefModel->IntegratedBondSquaredVol( currentTime, nextTime, nextTime );
		FwdDrift = exp( FwdDrift );

		ARM_GP_Vector times(1), values(1);
		values.Elt(0) = 1; times.Elt(0) = 0;
		ARM_CurveModelParam DummyModelParam( ARM_ModelParamType::Volatility, &times, &values);
		double IntegratedBondVol = itsRefModel->VolatilityScalarProduct( currentTime, nextTime, probaDate, DummyModelParam );

		for( i=0;i<statesNb; ++i )
		{
			f			= states->GetModelState(i,modelNb);
			sigma		= states->GetModelState(i,modelNb+1);
			W1			= states->GetNumMethodState(i,modelNb);
			W2			= states->GetNumMethodState(i,modelNb+1);
			if(f<=K_NEW_DOUBLE_TOL)
			{
				f=K_NEW_DOUBLE_TOL;
				sigma	   *= exp( -0.5*volOfVolSq+voldrift*dt+volOfVol*sqrDt*W2) ;
			}
			else
			{
				if(beta==1) fBetaMinus1=1.; else fBetaMinus1	= pow(f,beta-1.0);
				volTerm		= sigma*fBetaMinus1*sqrDt;
				drift		= -divPart-0.5*volTerm*volTerm;
				f		   *= exp(drift+volTerm*W1);
				f          *= FwdDrift / DFs->Elt(i);
				sigma	   *= exp( -0.5*volOfVolSq+voldrift*dt+volOfVol*sqrDt*W2 );
			}

			/// store it in the model state
			states->SetModelState(i,modelNb,f);
			states->SetModelState(i,modelNb+1,sigma);
		}
	}
	else  ////  we use here the parameter rho to create the correlation between the brownians
	{
		for( i=0;i<statesNb; ++i )
		{
			f			= states->GetModelState(i,modelNb);
			sigma		= states->GetModelState(i,modelNb+1);
			W1			= states->GetNumMethodState(i,modelNb);
			W2			= states->GetNumMethodState(i,modelNb+1);
			if(f<=K_NEW_DOUBLE_TOL)
			{
				f=K_NEW_DOUBLE_TOL;
				sigma	   *= exp( -0.5*volOfVolSq+voldrift*dt+volOfVol*sqrDt*(rho*W1+rhoOrt*W2) );
			}
			else
			{
				if(beta==1) fBetaMinus1=1.; else fBetaMinus1	= pow(f,beta-1.0);
				volTerm		= sigma*fBetaMinus1*sqrDt;
				drift		= -divPart-0.5*volTerm*volTerm;
				f		   *= irTerm * exp(drift+volTerm*W1); 
				sigma	   *= exp( -0.5*volOfVolSq+voldrift*dt+volOfVol*sqrDt*(rho*W1+rhoOrt*W2) );
			}

			/// store it in the model state
			states->SetModelState(i,modelNb,f);
			states->SetModelState(i,modelNb+1,sigma);
		}
	}
}

/*
////////////////////////////////////////////////////
///	Class  : ARM_SABRModel_Fx
///	Routine: UpdateLinks
///	Returns: void
///	Action : maintains links with itsIRModel
////////////////////////////////////////////////////
void ARM_SABRModel_Fx::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{
	const ARM_PricingModelIR* itsRefModel = getIRModel();
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();

	string IRModelName, 
		itsModelName = GetModelName();

	/// Find the name of its IRModel
	if( !itsRefModel )
	{
		/// If there is no IRModel, its name is found through its Other Models in the ModelNameMap
		ARM_IntVector itsOtherModels = (*modelMap)[itsModelName]->OtherModelRefNb();

		if( itsOtherModels.size() != 1 )
			ARM_THROW( ERR_INVALID_ARGUMENT, "Wrong number of otherModels" );

		IRModelName = (*modelMap)[itsOtherModels[0]]->ModelName();
	}
	else
		IRModelName = itsRefModel->GetModelName();

	/// The updated IRModel is set. 
	setIRModel( (*modelMap)[IRModelName]->Model() );

	/// Correlation between EquityInflation and IRModel is updated. 
	ARM_GP_MatrixPtr itsCorrelMatrix = multiAssetsModel.GetCorrelSubMatrix( IRModelName, itsModelName );

	/// Check that the IR/Inflation CorrelMatrix has the right size
	if( itsCorrelMatrix->cols() != 3 || itsCorrelMatrix->rows() != 3 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Inflation/IR Correl Matrix has the wrong size!" );

	size_t irModelNb = getIRModel()->GetModelNb();
	size_t SABRModelNb = GetModelNb();

	if( irModelNb < SABRModelNb )
	{
		itsIRSpotCorrel = itsCorrelMatrix->Elt(0,1);
		itsIRVolCorrel = itsCorrelMatrix->Elt(0,2);
		itsSpotVolCorrel = itsCorrelMatrix->Elt(1,2);
	}
	else
	{
		itsIRSpotCorrel = itsCorrelMatrix->Elt(2,0);
		itsIRVolCorrel = itsCorrelMatrix->Elt(2,1);
		itsSpotVolCorrel = itsCorrelMatrix->Elt(0,1);
	}
	/// Check that the correl provided as modelparam is the same as the correl provided in model map
		double rho			= GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation ).GetValue(0);

	if (fabs(itsSpotVolCorrel-rho)>1e-10)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "rho is different from the correl given in the model map!" );

	}
}
*/

#undef	FIRST_STATE_VARIABLE

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

