
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/hwxsvmmspread.h"
#include "gpmodels/svmmspread.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/hw2f.h"
#include "gpmodels/modelparamshw1f.h"
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

ARM_HWxSVMMSpread::ARM_HWxSVMMSpread(const ARM_ModelNameMap& modelNameMap, ARM_HullWhite2F * correlModel,
									 const ARM_GP_Vector& corrEndTimes, double ConstantCrossCorrel) : 
ARM_MultiAssetsModel(&modelNameMap), itsNbEffectiveReset(0), itsEigenValues(0), itsxCorrel(0), itsCorrEndTimes(corrEndTimes),
itsConvexAdjDone(false), itsFakeParam(ARM_ModelParamType::Volatility, 0.)
{
	if(modelNameMap.size() != 2)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : must have two models to build hw x svmm spread");
	}

	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(modelNameMap)[0]->Model());
	itsSpreadModel	= dynamic_cast<ARM_SVMMSpread*>(&*(modelNameMap)[1]->Model());

	if(itsHWModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" first model in map should be hull white model");
	}
	if(itsSpreadModel == NULL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" second model in map should be spread model");
	}

	itsResetTimes = itsSpreadModel->GetResetTimes();
	itsNbEffectiveReset = itsSpreadModel->GetNbEffectiveReset();

	if(itsCorrEndTimes.size() > 0 && itsCorrEndTimes.size() != itsResetTimes.size())
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" End Times Index for Cross Correlation should have same size as spread reset times");
	}

	computeCrossCorrel(itsCorrEndTimes.size() == 0 ? NULL : correlModel, ConstantCrossCorrel);

	itsLastTHW = ( (ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )).GetCurve()->GetAbscisses()[itsResetTimes.size()];

	itsSpreadGVol = itsSpreadModel->GetATMGaussianVol();

	ARM_GP_Vector times(1,0.), values(1, 1.);
	itsFakeParam.SetValuesAndTimes(&times, &values);
}

ARM_HWxSVMMSpread::ARM_HWxSVMMSpread(const ARM_HWxSVMMSpread& rhs) : ARM_MultiAssetsModel(rhs),
itsNbEffectiveReset(rhs.itsNbEffectiveReset), itsResetTimes(rhs.itsResetTimes), itsCorrEndTimes(rhs.itsCorrEndTimes),
itsConvexAdjDone(false), itsLastTHW(rhs.itsLastTHW), itsSpreadGVol(rhs.itsSpreadGVol), itsFakeParam(rhs.itsFakeParam)
{
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsxCorrel, itsxCorrel );

	itsHWModel		= dynamic_cast<ARM_HullWhite1F*>(&*(*GetModelMap())[0]->Model() );
	itsSpreadModel	= dynamic_cast<ARM_SVMMSpread*>(&*(*GetModelMap())[1]->Model() );
}

ARM_HWxSVMMSpread::~ARM_HWxSVMMSpread()
{
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );	
	DeletePointorVector<ARM_GP_Vector>( itsxCorrel );
}

void ARM_HWxSVMMSpread::computeCrossCorrel(ARM_HullWhite2F * correlModel, double ConstantCrossCorrel)
{
	itsxCorrel.resize(itsResetTimes.size());

	int i, k, dim, size = (int)itsResetTimes.size();

	double start, startk, end1k, end2k;

	for(i = 0; i < size; i++)
	{
		dim = size - i;

		itsxCorrel[i] = new ARM_GP_Vector(dim, 0.);

		if(correlModel == NULL) 
		{
			for(k = i; k < size; k++)
			{
				(*itsxCorrel[i])[k-i] = ConstantCrossCorrel;
			}
		}
		else
		{
			start = itsResetTimes[i];

			// calcul de la correl entre diag et spread
			for(k = i; k < size; k++)
			{
				startk	= itsSpreadModel->GetStartTimes()[k];
				end1k	= itsSpreadModel->GetEndTimes1()[k];
				end2k	= itsSpreadModel->GetEndTimes2()[k];

				(*itsxCorrel[i])[k-i] = correlModel->UnderlyingCorrelation("SO/CMS",0., start, 
										startk, end1k, startk, end2k, start, itsCorrEndTimes[i], startk, end1k);
			}
		}
	}
}

void ARM_HWxSVMMSpread::computeConvexAdj() const
{
	// on commence par calculer les vols gaussiennes de spread atm

    ARM_GP_Vector sigmaTimes( ( (ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )).GetCurve()->GetAbscisses());  // Ui
    ARM_GP_Vector sigmaValues(( (ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )).GetCurve()->GetOrdinates()); // Sigma(Ui)
	
	double mr = ((ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
	
	int k, size = itsResetTimes.size();

	if(sigmaTimes.size() - 1 != size)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT, ARM_USERNAME+" hw and svmm should have same calibration dates");
	}

	double integ = 0;
	double lastT = sigmaTimes[size]/365.;//GetNumeraire()->GetMaturity() / 365.;
	double t, T, rhokhw;

	for(k = 0; k < size; k++)
	{
		t = sigmaTimes[k]/365.;
		T = sigmaTimes[k+1]/365.;

		if(fabs(mr) < K_NEW_DOUBLE_TOL)
		{
			integ += sigmaValues[k+1] * 0.5 * ((lastT - t)*(lastT - t) - (lastT - T)*(lastT - T));
		}
		else
		{
			integ += sigmaValues[k+1] * ((T - t) / mr - (exp(-mr*(lastT-T)) - exp(-mr*(lastT-t)))/mr/mr);
		}

		rhokhw = (*itsxCorrel[k])[0];

		itsSpreadModel->CorrectIthFwd(k,- itsSpreadGVol[k] * integ * rhokhw);
	}
}

void ARM_HWxSVMMSpread::computeConvexAdj(ARM_PricingStatesPtr& states, double t, double T, int iFirst, int iLast) const
{
	bool proxy = false;

	if(proxy)
	{
		double hwsig = ((ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility )).GetCurve()->Interpolate(T);
		double mr = ((ARM_CurveModelParam&) itsHWModel->GetModelParams()->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

		double adj;

		if(fabs(mr) < K_DOUBLE_TOL)
		{
			adj = hwsig * 0.5 * ((itsLastTHW - t)*(itsLastTHW - t) - (itsLastTHW - T)*(itsLastTHW - T)) / 365. / 365.;
		}
		else
		{
			adj = hwsig * ((T - t) / 365. / mr - (exp(-mr*(itsLastTHW - T) / 365.) - exp(-mr*(itsLastTHW - t)/365.))/mr/mr);
		}

		int k, n, statesNb = states->size();

		double spreadk, adjk;
		int fwdsNb = itsNbEffectiveReset;

		for(k = iFirst; k <= iLast; k++)
		{
			adjk = - itsSpreadGVol[k] * (*itsxCorrel[iFirst])[k-iFirst] * adj;

			for(n = 0; n < statesNb; n++)
			{
				spreadk = states->GetModelState(n, k+1) + adjk;
				states->SetModelState(n, k+1, spreadk);	
			}
		}
	}
	else
	{
		double Ti, Tn = GetNumeraire()->GetMaturity();

		const ARM_ModelParamsHW1F * const hwModelParams = dynamic_cast<const ARM_ModelParamsHW1F * const>(itsHWModel->GetModelParams());

		double convZCTn = ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Tn);

		int i, n , statesNb = states->size();

		double adji, sigi, voli, spreadi;

		for(i = iFirst; i <= iLast; i++)
		{
			Ti = itsSpreadModel->GetStartTimes()[i];
			
			adji = (*itsxCorrel[iFirst])[i-iFirst] * (convZCTn - ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(itsFakeParam, hwModelParams, t, T, Ti));

			sigi =  ((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->GetLevel(i);

			adji = adji * sigi;

			for(n = 0; n < statesNb; n++)
			{
				voli = sqrt(states->GetModelState(n, itsNbEffectiveReset + 1));
				spreadi = states->GetModelState(n, i+1);
				states->SetModelState(n,i+1,spreadi+adji*voli);
			}
		}

	}
}

string ARM_HWxSVMMSpread::toString(const string& indent, const string& nextIndent) const
{
 	const ARM_ModelNameMap& modelMap = *GetModelMap();

    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Hybrid HW x SVMM Spread Model\n";
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

    os << indent << "\n\n------> SVMM Spread Model  <------\n";
    os << modelMap[1]->Model()->toString(indent,nextIndent);

	os << indent << "----------------------\n";
	os << indent << " Cross Correlation \n";

	int size = (int)itsResetTimes.size();

	for(int i = 0; i < size; i++)
	{
		os << indent << "( ResetDate : \t" << CC_NS(std,setprecision)(0) <<itsResetTimes[i] << " ) \t"<<CC_NS(std,setprecision)(6);
		for(int j = 0; j < (int)(*itsxCorrel[i]).size(); j++)
		{
			os	<< CC_NS(std,fixed)   << CC_NS(std,setprecision)(5) 
				<< CC_NS(std,setw)(8) << (*itsxCorrel[i])[j] << "\t";
		}
		os << CC_NS(std,endl);
	}

    return os.str();
}

ARM_BoolVector ARM_HWxSVMMSpread::NeedMCIntegProcess() const
{
	ARM_BoolVector ret(FactorCount(), false);

	return ret;
}

ARM_PricingStatesPtr ARM_HWxSVMMSpread::FirstPricingStates(size_t bucketSize) const
{
	if(itsConvexAdjDone == false)
	{
//		computeConvexAdj();
		itsConvexAdjDone = true;
	}

	return ARM_MultiAssetsModel::FirstPricingStates(bucketSize);
}

void ARM_HWxSVMMSpread::NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
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

void ARM_HWxSVMMSpread::NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances) const
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

void ARM_HWxSVMMSpread::ModelStateLocalVariancesAndStdDev(const ARM_GP_Vector& timeSteps)
{
//	((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->SetFactorCount(0);
//	return;

	itsSpreadModel->CalcNbEffectiveReset(timeSteps);

	itsNbEffectiveReset = itsSpreadModel->GetNbEffectiveReset();

	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances(timeSteps, auxLocalCov);

	itsSpreadModel->ModelStateLocalStdDev(timeSteps, auxLocalCov);

	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>(itsSpreadModel->GetEigenValues(), itsEigenValues);

//	SetModelStateLocalVars(itsSpreadModel->GetModelStateLocalVars());
}

void ARM_HWxSVMMSpread::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_HWxSVMMSpread::ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances) const
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

		for(j = 0; j < nbFwds; j++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->GetRho(j + idx);
			rhojhw	= (*itsxCorrel[idx])[j];
			rhoj_	= rhoj / sqrt(1. - rhojhw*rhojhw);

			(*localVariances[i])(j,j) = 1.;

			for(k = j+1; k < nbFwds; k++)
			{
				rhok	= ((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->GetRho(k + idx);
				rhokhw	= (*itsxCorrel[idx])[k];
				rhok_	= rhok / sqrt(1. - rhokhw*rhokhw);

				rhojk	= ((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);
				rhojk_	= (rhojk - rhojhw*rhokhw)/sqrt(1.- rhojhw*rhojhw)/sqrt(1.- rhokhw*rhokhw);

				(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk_ - rhoj_*rhok_)/sqrt(1.- rhoj_*rhoj_)/sqrt(1.- rhok_*rhok_);
			}
		}

		fromTime = toTime;
	}
}

void ARM_HWxSVMMSpread::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
//	HWMCFromToNextTime(states, timeIndex);
//	return;

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

	const ARM_MatrixVector& modelLocalStdev	= itsSpreadModel->GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols(); 

	// Corrélation des gaussiennes
	ARM_GP_Matrix x(statesNb, aliveFwds + 1);

	ARM_GP_Vector fachw(aliveFwds), facx(aliveFwds);
	double rhoi, rhoihw, rhoi_;

	for(i = 0; i < aliveFwds; i++)
	{
		rhoi		= ((ARM_ModelParamsBGMSV1F*) itsSpreadModel->GetModelParams())->GetRho(i + iFirst);
		rhoihw		= (*itsxCorrel[iFirst])[i];
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
	
	// Simulation Spread Model
	computeConvexAdj(states, GetNumMethod()->GetTimeStep(timeIndex), GetNumMethod()->GetTimeStep(timeIndex+1), iFirst, iLast);

	itsSpreadModel->MCFromToNextTime(states, timeIndex, x);
}

void ARM_HWxSVMMSpread::HWMCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex) const
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

CC_END_NAMESPACE()