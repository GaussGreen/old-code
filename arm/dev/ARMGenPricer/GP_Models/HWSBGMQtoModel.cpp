
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/hwsbgmqtomodel.h"
#include "gpmodels/svmmspread.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/LN_Fx.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/SmiledFRM.h"
#include "gpmodels/ModelParamsSmiledFRM.h"
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

ARM_HWSBGMQtoModel::ARM_HWSBGMQtoModel(const ARM_ModelNameMap& modelNameMap, 
									   const ARM_CurveMatrix& correlCurveMatrix) : 
ARM_MultiAssetsModel(&modelNameMap, &correlCurveMatrix), itsFakeParam(ARM_ModelParamType::Volatility, 0.),
itsSBGMoffsetIdx(0), itsEigenValues(0)
{
	if(modelNameMap.size() != 3)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : must have three models to build hw x sbgm quanto");
	}

	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(modelNameMap)[0]->Model());
	itsSBGMModel	= dynamic_cast<ARM_SmiledFRM*>(&*(modelNameMap)[1]->Model());
	itsFXModel		= dynamic_cast<ARM_LN_Fx*>(&*(modelNameMap)[2]->Model());

	if(itsHWModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" first model in map should be hull white model");
	}
	if(itsSBGMModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" second model in map should be sbgm model");
	}
	if(itsFXModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" third model in map should be fx model");
	}

	ARM_GP_Vector times(1,0.), values(1, 1.);
	itsFakeParam.SetValuesAndTimes(&times, &values);

	// pour ne pas le simuler
	(*GetModelMap())[2]->SetUsedInPricing( false );
	itsFXModel->SetFactorCount(0);
}


ARM_HWSBGMQtoModel::ARM_HWSBGMQtoModel(const ARM_HWSBGMQtoModel& rhs) : ARM_MultiAssetsModel(rhs),
itsFakeParam(rhs.itsFakeParam), itsSBGMoffsetIdx(rhs.itsSBGMoffsetIdx)
{
	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(*GetModelMap())[0]->Model());
	itsSBGMModel	= dynamic_cast<ARM_SmiledFRM*>(&*(*GetModelMap())[1]->Model());
	itsFXModel		= dynamic_cast<ARM_LN_Fx*>(&*(*GetModelMap())[2]->Model());

	(*GetModelMap())[2]->SetUsedInPricing( false );
	itsFXModel->SetFactorCount(0);

	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
}

ARM_HWSBGMQtoModel::~ARM_HWSBGMQtoModel()
{
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );	
}

string ARM_HWSBGMQtoModel::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "HW x SBGM Quanto Model\n";
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

    os << indent << "\n\n------> SBGM Model  <------\n";
    os << modelMap[1]->Model()->toString(indent,nextIndent);

    os << indent << "\n\n------> FX Model  <------\n";
    os << modelMap[2]->Model()->toString(indent,nextIndent);

    return os.str();
}

ARM_BoolVector ARM_HWSBGMQtoModel::NeedMCIntegProcess() const
{
	ARM_BoolVector ret(FactorCount(), false);

	return ret;
}

void ARM_HWSBGMQtoModel::NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
{
	size_t totalFwds     = itsSBGMModel->GetResetDates().size(),
		   factorsNb     = itsSBGMModel->FactorCount(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   startTimePos  = 0;
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));

	for(size_t i=0;i<timeStepsSize-1 ;++i)
	{
		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb+1,factorsNb+1,0.);
		(*localVariances[i])(0,0) = 1.;
		for (size_t k=0; k<factorsNb;k++)
			(*localVariances[i])(k+1,k+1) = (*itsEigenValues[i])[k];
	}
}

void ARM_HWSBGMQtoModel::NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances) const
{
	size_t totalFwds     = itsSBGMModel->GetResetDates().size(),
		   timeStepsSize = timeSteps.size(),
		   modelNb		 = GetModelNb(),
		   offsetIndex	 = (timeStepsSize-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledFRM,NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (size_t i=0;i<timeStepsSize;i++){
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds+1,1.0);
	}
}

void ARM_HWSBGMQtoModel::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);

	itsSBGMoffsetIdx = timeSteps.size()-1;
}

void ARM_HWSBGMQtoModel::ModelStateLocalVariancesAndStdDev(const ARM_GP_Vector& timeSteps)
{
	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances(timeSteps, auxLocalCov);

	itsSBGMModel->ModelStateLocalStdDev(timeSteps, auxLocalCov);

	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(itsSBGMModel->GetEigenValues(), itsEigenValues);
}

void ARM_HWSBGMQtoModel::ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
{
	size_t nbSteps = timeSteps.size();
	localVariances.reserve( (nbSteps-1)*2);	

	itsHWModel->ModelStateLocalVariances(timeSteps, localVariances);

	// la partie SBGM
	size_t totalFwds = itsSBGMModel->GetResetDates().size(),
	   factorsNb     = itsSBGMModel->FactorCount(),
	   timeStepsSize = timeSteps.size(),
	   modelNb		 = 1,
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

	ARM_ModelParamsSmiled * SBGMparams = dynamic_cast<ARM_ModelParamsSmiled*>(itsSBGMModel->GetModelParams());

	for(i=0, ii = offsetIndex;i<timeStepsSize-1 ;++i, ++ii)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		/// initialize everything
		localVariances[ii] = new ARM_GP_Matrix(totalFwds,totalFwds);
		
		/// get the first bigger time
		while(startTimePos< itsSBGMModel->GetResetDates().size()
			&& itsSBGMModel->GetResetDates()[startTimePos] < toTime )
			++startTimePos;

		if (itsSBGMModel->GetAllowInterpol() && startTimePos>0) startTimePos--;	

		for(j=0; j<startTimePos; ++j)
		{
			for(k=j;k<totalFwds;++k)
				(*localVariances[ii])(j,k)=(*localVariances[ii])(k,j)=0.;
		}

		double rhojk, rhoDomFor;

		ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(toTime);

		rhoDomFor = correlMatrix(1,0);

		for(j=startTimePos; j<totalFwds; ++j)
		{
			(*localVariances[ii])(j,j) = 1.;

			for(k=j+1;k<totalFwds;++k)
			{
				rhojk = SBGMparams->InstantCorrel(toTime,j,k);
				(*localVariances[ii])(j,k)=(*localVariances[ii])(k,j)=(rhojk - rhoDomFor*rhoDomFor)/(1. - rhoDomFor*rhoDomFor);
			}
		}

		fromTime= toTime;
	}

}

void ARM_HWSBGMQtoModel::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
{
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt		= (nextTime - currTime) / K_YEAR_LEN;
	double racdt	= sqrt(dt);
	int statesNb	= states->size();
	int modelNb		= GetModelNb();

	// les params HW
	double mr		= ((ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	double adj		= exp(-mr*dt);
	double var		= ((const ARM_ModelParamsHW1F* const) itsHWModel->GetModelParams())->StateLocalVariance(currTime, nextTime, nextTime);
	double hstddev	= sqrt(var);

	// les params SBGM
	size_t fwdsNb  = itsSBGMModel->GetResetDates().size();
	size_t iLast   = fwdsNb-1;

	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsSBGMModel->GetResetDates()[iFirst]<nextTime)
	{
		++iFirst;
	}
	if (itsSBGMModel->GetAllowInterpol() && iFirst>0) iFirst--;

	//parait débile mais indispensable pour les échanciers non triviaux!!!
	size_t iFirstMartingale = iFirst;
	while( itsSBGMModel->GetEndDates()[iFirstMartingale] < itsSBGMModel->GetEndDates()[iLast]) iFirstMartingale++;
	
	const ARM_MatrixVector& modelLocalVar	= itsSBGMModel->GetModelStateLocalVars();
	ARM_GP_Vector* eigen = itsEigenValues[timeIndex];

	// les ajustements de convexité quanto
	ARM_GP_Vector qtoAdj = GetQuantoConvexAdj(currTime, nextTime, iFirst, fwdsNb-1);

	int i,j;
	double old,next,modStd;
	double x;

	ARM_GP_Vector stddev(fwdsNb,0.);

	for(i = iFirst; i <= itsSBGMModel->GetNumeraireTimeIndex(); i++)
	{
		stddev[i] = sqrt( dynamic_cast<ARM_ModelParamsSmiled*>(itsSBGMModel->GetModelParams())->IntegratedVariance(currTime, nextTime, i) );
	}

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(nextTime);

	double corrDomFor = correlMatrix(0,1);
	double UnMoinsCorr2 = sqrt(1. - corrDomFor*corrDomFor);

	for(int n = 0; n < statesNb; n++)
	{
		// la simulation Hull White
		x = states->GetModelState(n, 0);
		x = adj * x + hstddev * states->GetNumMethodState(n, 0);

		states->SetModelState(n, 0, x);

		// la simulation SBGM
		for (i = itsSBGMModel->GetNumeraireTimeIndex(); i >= iFirst ;--i)
		{
			old			= states->GetModelState(n,i+1+modelNb);
			next		= old;
			j			= 0;
			size_t prec	= i+1;

			if (itsSBGMModel->GetEndDates()[i] == itsSBGMModel->GetEndDates()[iLast])
				prec = i;
			else while (!itsSBGMModel->IsOnSamePath(i,prec)) prec++;
			
			x = 0.;
			for (ARM_GP_Vector::iterator iter = eigen->begin(); iter != eigen->end() ; ++iter,++j )
			{
				modStd	=	(*modelLocalVar[timeIndex+itsSBGMoffsetIdx])(i,j);
				x		+=	modStd * ( states->GetNumMethodState(n,j+1+modelNb));
			}
			
			next += qtoAdj[i] + stddev[i] * (corrDomFor * states->GetNumMethodState(n,0) + UnMoinsCorr2 * x);

			states->SetModelState(n,i+1+modelNb,next);
		}
		
	}
}

ARM_GP_Vector ARM_HWSBGMQtoModel::GetQuantoConvexAdj(double t, double T, int iFirst, int iLast) const
{
	ARM_GP_Vector adj(iLast+1);

	double volfx = ((ARM_CurveModelParam&) itsFXModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).GetCurve()->Interpolate(t);
	
	double voli;
	double Ti, Tn = GetNumeraire()->GetMaturity();

	ARM_GP_Matrix correlMatrix = GetCorrelMatrix()->Interpolate(T);

	double corrDomFor = correlMatrix(0,1);
	double corrForFX = correlMatrix(1,2);

	const ARM_ModelParamsHW1F * const hwModelParams = dynamic_cast<const ARM_ModelParamsHW1F * const>(itsHWModel->GetModelParams());

	double convZCTn = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Tn);

	for(int i = iFirst; i <= iLast; i++)
	{
		Ti = itsSBGMModel->GetEndDates()[i];

		voli = ((ARM_ModelParamsSmiled*) itsSBGMModel->GetModelParams())->GetVolCurve(i)->Interpolate(T);

		// ajustement ti forward neutre foreign -> ti forward neutre domestique
		adj[i] = corrForFX * volfx * voli * (T - t) / 365.;

		// ajustement ti forward neutre domestique -> terminal forward neutre domestique
		adj[i] += corrDomFor * voli * (convZCTn - ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Ti));
	}

	return adj;
}

CC_END_NAMESPACE()