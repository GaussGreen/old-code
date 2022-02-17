#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/hwsvbgmqtomodel.h"
#include "gpmodels/bgmsv1f.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/LN_Fx.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/datestrip.h"
#include "gpbase/vectormanip.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"

#include "gpclosedforms/normal.h"
#include "gpnumlib/gaussiananalytics.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_HWSVBGMQtoModel::ARM_HWSVBGMQtoModel(const ARM_ModelNameMap& modelNameMap, 
									   const ARM_CurveMatrix& correlCurveMatrix) : 
ARM_MultiAssetsModel(&modelNameMap, &correlCurveMatrix), itsEigenValues(0), 
itsFakeParam(ARM_ModelParamType::Volatility, 0.)
{
	if(modelNameMap.size() != 3)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : must have three models to build hw x sbgm quanto");
	}

	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(modelNameMap)[0]->Model());
	itsSVMMModel	= dynamic_cast<ARM_BGMSV1F*>(&*(modelNameMap)[1]->Model());
	itsFXModel		= dynamic_cast<ARM_LN_Fx*>(&*(modelNameMap)[2]->Model());

	itsSVMMModel->SetProxyStatus(true);

	if(itsHWModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" first model in map should be hull white model");
	}
	if(itsSVMMModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" second model in map should be bgm sv model");
	}
	if(itsFXModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" third model in map should be fx model");
	}

	itsResetTimes = itsSVMMModel->GetResetTimes();
	itsNbEffectiveReset = itsSVMMModel->GetNbEffectiveReset();

	(*GetModelMap())[2]->SetUsedInPricing( false );
	itsFXModel->SetFactorCount(0);

	ARM_GP_Vector times(1,0.), values(1, 1.);
	itsFakeParam.SetValuesAndTimes(&times, &values);
}

ARM_HWSVBGMQtoModel::ARM_HWSVBGMQtoModel(const ARM_HWSVBGMQtoModel& rhs) : ARM_MultiAssetsModel(rhs),
itsNbEffectiveReset(rhs.itsNbEffectiveReset), itsResetTimes(rhs.itsResetTimes), itsFakeParam(rhs.itsFakeParam)
{
	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(*GetModelMap())[0]->Model());
	itsSVMMModel	= dynamic_cast<ARM_BGMSV1F*>(&*(*GetModelMap())[1]->Model());
	itsFXModel		= dynamic_cast<ARM_LN_Fx*>(&*(*GetModelMap())[2]->Model());

	itsSVMMModel->SetProxyStatus(true);

	(*GetModelMap())[2]->SetUsedInPricing( false );
	itsFXModel->SetFactorCount(0);

	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
}

ARM_HWSVBGMQtoModel::~ARM_HWSVBGMQtoModel()
{
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );	
}

string ARM_HWSVBGMQtoModel::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "HW x SVBGM Quanto Model\n";
    os << indent << "----------------------------------------------\n\n";

	if( GetRefModel() )
	{
		os << "Payment model : " << GetRefModel()->GetModelName() << "\n\n";
	}

	if( GetCorrelMatrix() )
	{
		os << "Correlation matrix\n";
		os << GetCorrelMatrix()->toString(indent,nextIndent);
	}

    os << indent << "\n\n------> HW Model <------\n";
    os << modelMap[0]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------> SVBGM Model  <------\n";
    os << modelMap[1]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------> FX Model  <------\n";
    os << modelMap[2]->Model()->toString(indent,nextIndent);

    return os.str();
}

ARM_BoolVector ARM_HWSVBGMQtoModel::NeedMCIntegProcess() const
{
	ARM_BoolVector ret(FactorCount(), false);

	return ret;
}

void ARM_HWSVBGMQtoModel::NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
{
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));

	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		/// initialize everything
		// localVariances[i] = new ARM_GP_Matrix(1,1,1.); continue;

		int factorsNb = itsEigenValues[i]->size() + 1 + 1;
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		(*localVariances[i])(0,0) = 1.;

		for (int k = 0; k < factorsNb-2; k++)
			(*localVariances[i])(k+1,k+1) = (*itsEigenValues[i])[k];

		(*localVariances[i])(factorsNb-1,factorsNb-1) = 1.;
	}
}

void ARM_HWSVBGMQtoModel::NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (int i=0;i<timeStepsSize;i++)
	{
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds+1+1,1.0);
	}
}

void ARM_HWSVBGMQtoModel::ModelStateLocalVariancesAndStdDev(const ARM_GP_Vector& timeSteps)
{
	itsSVMMModel->CalcNbEffectiveReset(timeSteps);

	itsNbEffectiveReset = itsSVMMModel->GetNbEffectiveReset();

	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances(timeSteps, auxLocalCov);

	itsSVMMModel->ModelStateLocalStdDev(timeSteps, auxLocalCov);

	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(itsSVMMModel->GetEigenValues(), itsEigenValues);
}

void ARM_HWSVBGMQtoModel::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_HWSVBGMQtoModel::ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= 0.;//GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;


#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ModelStateLocalVariances: first time step != 0" );
#endif

    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	
	int i,j,k, idx = 0;

	int nbFwds;

	for(i = 0;i < timeStepsSize - 1 ; ++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		
		/// get the first bigger time
		while(idx< itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		if(nbFwds <= 0)
		{
			localVariances[i] = new ARM_GP_Matrix(1,1);
			(*localVariances[i])(0,0) = 0.;
			continue;
		}

		localVariances[i] = new ARM_GP_Matrix(nbFwds, nbFwds);	

		double rhoj, rhok, rhojk, rhojhw, rhokhw, rhoj_, rhok_, rhojk_;

		ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(toTime);

		rhojhw = rhokhw = correlMatrix(1,0);

		for(j = 0; j < nbFwds; j++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV1F*) itsSVMMModel->GetModelParams())->GetRho(j + idx);
			rhoj_	= rhoj / sqrt(1. - rhojhw*rhojhw);

			(*localVariances[i])(j,j) = 1.;

			for(k = j+1; k < nbFwds; k++)
			{
				rhok	= ((ARM_ModelParamsBGMSV1F*) itsSVMMModel->GetModelParams())->GetRho(k + idx);
				rhok_	= rhok / sqrt(1. - rhokhw*rhokhw);

				rhojk	= ((ARM_ModelParamsBGMSV1F*) itsSVMMModel->GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);
				rhojk_	= (rhojk - rhojhw*rhokhw)/sqrt(1.- rhojhw*rhojhw)/sqrt(1.- rhokhw*rhokhw);

				(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk_ - rhoj_*rhok_)/sqrt(1.- rhoj_*rhoj_)/sqrt(1.- rhok_*rhok_);
			}
		}

		fromTime = toTime;
	}
}

void ARM_HWSVBGMQtoModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	int statesNb	= states->size();
	int fwdsNb		= itsNbEffectiveReset;
	int iLast		= fwdsNb-1;
	
	int i, j, n;
	
	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetTimes[iFirst]<nextTime)
	{
		++iFirst;
	}

	int aliveFwds	= iLast - iFirst + 1;

	if(aliveFwds <= 0) return;

	const ARM_MatrixVector& modelLocalStdev	= itsSVMMModel->GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols(); 

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(nextTime);

	// Corrélation des gaussiennes
	ARM_GP_Matrix x(statesNb, aliveFwds + 1);

	ARM_GP_Vector fachw(aliveFwds), facx(aliveFwds);
	double rhoi, rhoi_;
	double rhoihw = correlMatrix(1,0);

	for(i = 0; i < aliveFwds; i++)
	{
		rhoi		= ((ARM_ModelParamsBGMSV1F*) itsSVMMModel->GetModelParams())->GetRho(i + iFirst);
		rhoi_		= rhoi / sqrt(1. - rhoihw*rhoihw);

		fachw[i]	= rhoihw / sqrt(1. - rhoi*rhoi);
		facx[i]		= sqrt(1.- rhoi_*rhoi_)*sqrt(1.- rhoihw*rhoihw)/sqrt(1.- rhoi*rhoi);
	}

	for(n = 0; n < statesNb; n++)
	{
		for(i = 0; i < aliveFwds; i++)
		{
			x(n,i) = 0.;
			for(j = 0; j < eigensNb; j++) x(n,i) += states->GetNumMethodState(n, j+1) * (*modelLocalStdev[timeIndex])(i,j);
			
			x(n,i) = fachw[i] * states->GetNumMethodState(n, 0) + facx[i] * x(n,i);
		}
		x(n, aliveFwds) = states->GetNumMethodState(n, eigensNb+1);
	}

	// Simulation HW
	HWMCFromToNextTime(states, timeIndex);
	
	// Simulation SV BGM
	computeConvexAdj(states, GetNumMethod()->GetTimeStep(timeIndex), GetNumMethod()->GetTimeStep(timeIndex+1), iFirst, iLast);

	itsSVMMModel->MCFromToNextTime(states, timeIndex, x);
}

void ARM_HWSVBGMQtoModel::HWMCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
{
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt		= (nextTime - currTime) / K_YEAR_LEN;
	double racdt	= sqrt(dt);
	int statesNb	= states->size();
	int modelNb		= GetModelNb();

	double mr		= ((ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	double adj		= exp(-mr*dt);
	double var		= ((const ARM_ModelParamsHW1F* const) itsHWModel->GetModelParams())->StateLocalVariance(currTime, nextTime, nextTime);
	double stddev	= sqrt(var);

	double x;

	for(int n = 0; n < statesNb; n++)
	{
		x = states->GetModelState(n, 0);
		x = adj * x + stddev * states->GetNumMethodState(n, 0);

		states->SetModelState(n, 0, x);
	}
}

void ARM_HWSVBGMQtoModel::computeConvexAdj(ARM_PricingStatesPtr& states, double t, double T, int iFirst, int iLast) const
{
	double volfx = ((ARM_CurveModelParam&) itsFXModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->Interpolate(t);
	
	double voli, libi;
	double sigi, adji, Ti, Tn = GetNumeraire()->GetMaturity();

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(T);

	double corrDomFor = correlMatrix(0,1);
	double corrForFX = correlMatrix(1,2);

	const ARM_ModelParamsHW1F * const hwModelParams = dynamic_cast<const ARM_ModelParamsHW1F * const>(itsHWModel->GetModelParams());

	double convZCTn = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Tn);

	int i, n , statesNb = states->size();

	for(i = iFirst; i <= iLast; i++)
	{
		Ti = itsSVMMModel->GetEndTimes()[i];
		
		adji = corrForFX * volfx * (T - t)/365. + corrDomFor * (convZCTn - ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Ti));

		sigi =  ((ARM_ModelParamsBGMSV1F*) itsSVMMModel->GetModelParams())->GetLevel(i);

		adji = adji * sigi;

		for(n = 0; n < statesNb; n++)
		{
			voli = sqrt(states->GetModelState(n, itsNbEffectiveReset + 1));
			libi = states->GetModelState(n, i+1);
			states->SetModelState(n,i+1,libi*exp(adji*voli));
		}
	}
}


CC_END_NAMESPACE()