/// gpmodels
#include "gpmodels/svmmspread.h"
/// gpbase
#include "gpbase/removeidentifiedwarning.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/comparisonfunctor.h"
#include "gpbase/datestrip.h"
#include "gpbase/vectormanip.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/discretisationscheme.h"

#include "gpmodels/ModelParamsSmiledFRM.h"
#include "gpclosedforms/normal.h"
#include "gpnumlib/gaussiananalytics.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/vanilla_normal.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_SVMMSpread::ARM_SVMMSpread( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params)
:	ARM_PricingModelIR(zc,params), 
	itsResetTimes(NULL),
	itsStartTimes(NULL),
	itsEndTimes1(NULL),
	itsEndTimes2(NULL),
	itsDelta(NULL),
	itsFwdRate(NULL),
	itsEigenValues(0),
	itsModelLocalRealVar(0),
	itsCalibSecDensities(0),
	itsCalibrationStatus(true),
	itsCalibratedStatus(false)
{
	itsSimpleEulerScheme = false;
}


////////////////////////////////////////////////////
///	Class  : ARM_SVMMSpread
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SVMMSpread::ARM_SVMMSpread(const ARM_SVMMSpread& rhs)
:	ARM_PricingModelIR(rhs),
	itsResetTimes(rhs.itsResetTimes),
	itsStartTimes(rhs.itsStartTimes),
	itsEndTimes1(rhs.itsEndTimes1),
	itsEndTimes2(rhs.itsEndTimes2),
	itsDelta(rhs.itsDelta),
	itsFwdRate(rhs.itsFwdRate),
	itsSimpleEulerScheme(rhs.itsSimpleEulerScheme),
	itsCalibrationStatus(rhs.itsCalibrationStatus),
	itsCalibratedStatus(rhs.itsCalibratedStatus)
{
	DuplicateCloneablePointorVectorInPlace<std::vector<double>>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelLocalRealVar, itsModelLocalRealVar);
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsCalibSecDensities, itsCalibSecDensities);
}	


ARM_SVMMSpread& ARM_SVMMSpread::operator = (const ARM_SVMMSpread& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_SVMMSpread();

		new (this) ARM_SVMMSpread (rhs);
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_SVMMSpread
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SVMMSpread::~ARM_SVMMSpread()
{
	DeletePointorVector<std::vector<double>>( itsEigenValues );
	DeletePointorVector<ARM_GP_Matrix>( itsModelLocalRealVar );
}

bool ARM_SVMMSpread::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(params.DoesModelParamExist(ARM_ModelParamType::Alpha) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: alpha structure (level) is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::QParameter) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: Q structure (rho Rate/variance) is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: VolOfVol structure is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::VolMeanReversion) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: VolMeanReversion structure is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::InitialVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: InitialVol structure is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::LongTermVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: LongTermVol structure is needed for SVMMSpread" );

	if(params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateModelParams: beta correlation structure is needed for SVMMSpread" );

	return true;
}

void ARM_SVMMSpread::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	/*
	if(	numerairePtr->GetType() != ARM_Numeraire::TerminalZc && numerairePtr->GetType() != ARM_Numeraire::TerminalEventZc )
    {		
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::SetNumeraire: only numeraire supported is Terminal ZC" );
    }
	*/
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

void ARM_SVMMSpread::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);

//	if (!itsCalibratedStatus && !isNumerical)
//		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::ValidateCalibMethod: should calib functionals first!!!" );

	calibMethod.DefaultValidateWithModel(*this);
}

void ARM_SVMMSpread::PreProcessing(ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

void ARM_SVMMSpread::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	// TO DO 
}


void ARM_SVMMSpread::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	// TO DO 
}


ARM_PricingStatesPtr ARM_SVMMSpread::FirstPricingStates(size_t bucketSize) const
{
	const size_t nbPayoffs			= 0;
	
	size_t nbModelStates			= itsNbEffectiveReset + 1;
	/*
	size_t nbModelStates			= itsNbEffectiveReset + 1 +1 ;
	*/
	size_t factorsNb				= FactorCount();

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	int i,k;

	
	for(i = 0; i < nbModelStates - 1; i++)
	/*
	for(i = 0; i < nbModelStates - 2; i++)
	*/
	{
		for(k = 0; k < bucketSize; k++) 
		{
			initStates->SetModelState(k, i, itsFwdRate[i]);
		}
	}

	double initVarValue		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();

	for(k = 0; k < bucketSize; k++)
	{
		initStates->SetModelState(k, itsNbEffectiveReset, initVarValue);
	}

	/*
	for(k = 0; k < bucketSize; k++) initStates->SetModelState(k,itsNbEffectiveReset+1,0.);
	*/

	return initStates;
}

ARM_PricingStatesPtr ARM_SVMMSpread::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	if(itsCalibrationStatus)
	{
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}

	// Numerical method and numeraire are needed to go on
    ARM_NumMethodPtr numMethod=GetNumMethod();
	if( numMethod == ARM_NumMethodPtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Init: numerical method not set in SmiledFRM model!");

    /// test the numeraire and its type!
	ARM_NumerairePtr numeraire=GetNumeraire();
    if( numeraire == ARM_NumerairePtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledMM::Init: numeraire not set in the SmiledFRM model!");

	/// creates the model schedule (smart pointor for exception safety!)
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime(); //ARM_EventAndModelTime();
	CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

    //// Initialise the numeraire
	numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,ARM_GP_VectorPtr( new std::vector<double>( 1, getTerminalTime() )));

	/// Set the basic schedule in the numerical method and...
	numMethod->SetTimeSteps(*ptimeSteps);

	double firstInductTime = timeInfos[0]->GetEventTime();

	/// ...initialise it
	return numMethod->Init(*this,firstInductTime);
}

ARM_PricingStatesPtr ARM_SVMMSpread::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	if(itsCalibrationStatus)
	{	
		if(itsCalibratedStatus == false)
		{
			itsCalibratedStatus = true;
		}

		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}

	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );

	return newStates;
}

ARM_VectorPtr ARM_SVMMSpread::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
// FIXMEFRED: mig.vc8 (30/05/2007 16:17:49):cast
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,getTerminalTime()));
}

std::vector<double>& ARM_SVMMSpread::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetTimes.Clone());
}


void ARM_SVMMSpread::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);

	if( !modelParam )
    {
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");
	}

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    std::vector<double>  tmpdates;
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
	}
}


string ARM_SVMMSpread::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "\n\n";
    os << indent << "SVBGM1F : Shifted Stochastic Volatility BGM with One process for variance (Heston time dependant) \n";
    os << indent << "----------------------------\n";
	os << ARM_PricingModel::toString(indent);
	os << "\n\n";

	return os.str();
}

void ARM_SVMMSpread::NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{
//	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	
	if(itsSimpleEulerScheme)
	{
		for(int i = 0; i < timeStepsSize-1 ; ++i)
		{
			int factorsNb = itsEigenValues[i]->size();

			/// initialize everything
			localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

			for (int k = 0; k < factorsNb; k++)
				(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
		}
	}
	else
	{
		for(int i = 0; i < timeStepsSize-1 ; ++i)
		{
			
			int factorsNb = itsEigenValues[i]->size() + 1 ;

			/// initialize everything
			localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

			for (int k = 0; k < factorsNb-1; k++)
				(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];

			(*localVariances[i])(factorsNb-1,factorsNb-1) = 1.;
			
			/*
			int factorsNb = itsEigenValues[i]->size() + 1 + 1;

			/// initialize everything
			localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

			(*localVariances[i])(0,0) = 1.;

			for (int k = 0; k < factorsNb-2; k++)
				(*localVariances[i])(k+1,k+1) = (*itsEigenValues[i])[k];

			(*localVariances[i])(factorsNb-1,factorsNb-1) = 1.;
			*/
		}
	}
}

void ARM_SVMMSpread::NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
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
		
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds+1,1.0);
		/*
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds+1+1,1.0);
		*/
	}
}

void ARM_SVMMSpread::CalcNbEffectiveReset(const std::vector<double>& timeSteps)
{
	double lastTimeStep = timeSteps[timeSteps.size() - 1];
	itsNbEffectiveReset = 0;

	while(itsResetTimes[itsNbEffectiveReset] < lastTimeStep)
	{
		itsNbEffectiveReset ++;
		if(itsNbEffectiveReset == itsResetTimes.size()-1) break;
	}

	itsNbEffectiveReset ++;

	/*
	while(!(itsResetTimes[itsNbEffectiveReset] > lastTimeStep))
	{
		itsNbEffectiveReset ++;
		if(itsNbEffectiveReset == itsResetTimes.size()) break;
	}
	*/
	((ARM_ModelParamsBGMSV1F*) GetModelParams())->checkFactorCount(itsNbEffectiveReset);

}

void ARM_SVMMSpread::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	CalcNbEffectiveReset(timeSteps);

	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances( timeSteps, auxLocalCov);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);

	ModelStateLocalStdDev(timeSteps, auxLocalCov);
}

void ARM_SVMMSpread::ModelStateLocalStdDev( const std::vector<double>& timeSteps, const ARM_MatrixVector& auxLocalCov)
{
	ARM_MatrixVector modStateLocalVars;
//	ARM_MatrixVector auxLocalCov;

//	ModelStateLocalVariances( timeSteps, auxLocalCov);
//	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);
	
	int totalFwds		= itsNbEffectiveReset;
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= 0.; // GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;

	int i, j, k, startTimePos = 0;

	modStateLocalVars.resize((timeStepsSize-1)*(modelNb+1));

	DeletePointorVector<std::vector<double>>( itsEigenValues );
	itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

	int nbFwds, idx = 0, currfactorsNb = factorsNb;
	int maxfactorsNb = 0;
	double minratio = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetMinRatio();
	bool resizefactors = minratio > 1. || fabs(minratio - 1.) < K_DOUBLE_TOL ? false : true;

	// on commence par calculer toutes les matrices, et les vecteurs propres
	std::vector<double> * eigenValues = new std::vector<double>[timeStepsSize - 1];
	ARM_GP_Matrix ** ACPMatrix = new ARM_GP_Matrix * [timeStepsSize - 1];
	std::vector<double> factorsize(timeStepsSize-1, factorsNb);

	for(i = 0; i < timeStepsSize - 1; i++)
	{
		toTime  = timeSteps[i+1];

		/// get the first bigger time
		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		int esize = itsSimpleEulerScheme ? nbFwds + 1 : nbFwds;

		eigenValues[i].resize(esize,0.);

		ACPMatrix[i] = ACPTransformation(auxLocalCov[i],eigenValues[i]);

		// calcul de la dimension effective
		if(resizefactors)
		{
			// variance total
			double varTotal = 0.;
			for(k = 0; k < esize; k++) 
			{
				if(eigenValues[i][k] > 0.) varTotal += eigenValues[i][k];
			}

			double varExpliq = 0.;
			for(k = 0; k < esize; k++)
			{
				if(varExpliq > minratio * varTotal) break;
				if(eigenValues[i][k] > 0.)
				{
					varExpliq += eigenValues[i][k];
				}
			}

			factorsize[i]= k;
			maxfactorsNb = maxfactorsNb > k ? maxfactorsNb : k;
		}
		else
		{
			maxfactorsNb = factorsNb;
		}
	}

	
	if(itsSimpleEulerScheme)
		((ARM_ModelParamsBGMSV1F*) GetModelParams())->SetFactorCount(maxfactorsNb);
	else
		((ARM_ModelParamsBGMSV1F*) GetModelParams())->SetFactorCount(maxfactorsNb + 1);


	factorsNb = FactorCount();

	idx = 0;
	for(i = 0; i < timeStepsSize - 1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		int esize = itsSimpleEulerScheme ? nbFwds + 1 : nbFwds;

		currfactorsNb = factorsize[i] > esize ? esize : factorsize[i];

		/// initialize everything
		modStateLocalVars[i] = new ARM_GP_Matrix(esize, currfactorsNb);

		itsEigenValues[i] = new std::vector<double>(currfactorsNb);
		
		int maxRank = -1;
		for(k=0;k<currfactorsNb;++k)
		{
			if(eigenValues[i][k] < K_NEW_DOUBLE_TOL)
			{
				eigenValues[i][k] = 0.;
				if (maxRank == -1)
					maxRank = k;
			}
		}

		if (maxRank == -1) maxRank = currfactorsNb;
		size_t effectiveRank = ( maxRank < currfactorsNb ? maxRank : currfactorsNb);

		//rescaling
		
		for (j = 0;j < esize; j++)
		{
			double sum=0.;
			for(k=0;k<effectiveRank-1;++k){
				sum+=(*ACPMatrix[i])(j,k)*(*ACPMatrix[i])(j,k)*eigenValues[i][k];
			}
			double sgn = ((*ACPMatrix[i])(j,k)>0.?1.:-1.);
			double res = (*auxLocalCov[i])(j,j)-sum;
			if (res>K_NEW_DOUBLE_TOL)
				(*ACPMatrix[i])(j,k)=sgn*sqrt(res/eigenValues[i][k]); 
		}
		
		for(k = 0; k < currfactorsNb; ++k)
		{
			(*itsEigenValues[i])(k) = eigenValues[i][k];
			for(j = 0; j < esize; ++j)
			{	
				(*modStateLocalVars[i])(j,k) =	(*ACPMatrix[i])(j,k);
			}
		}

		fromTime= toTime;
	}

	delete [] eigenValues;
	for(i = 0; i < timeStepsSize - 1; i++) delete [] ACPMatrix[i];
	delete [] ACPMatrix;

	/// set the result
	SetModelStateLocalVars(modStateLocalVars);
}

void ARM_SVMMSpread::NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_SVMMSpread::ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{

	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
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

		if(itsSimpleEulerScheme)
		{
			localVariances[i] = new ARM_GP_Matrix(nbFwds+1, nbFwds+1);	

			// les correls taux taux
			
			for(j = 0; j < nbFwds; ++j)
			{
				for(k = j; k < nbFwds; ++k)
					(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);
			}
			
			// les correls taux vol
			for(j = 0; j < nbFwds; j++)
			{
					(*localVariances[i])(j, nbFwds) = (*localVariances[i])(nbFwds, j) = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(j + idx);
			}
			

			(*localVariances[i])(nbFwds, nbFwds) = 1.;
		}
		else
		{
			localVariances[i] = new ARM_GP_Matrix(nbFwds, nbFwds);	

			double rhoj, rhok, rhojk, fac;
			for(j = 0; j < nbFwds; j++)
			{
				rhoj	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(j + idx);

				(*localVariances[i])(j,j) = 1.;

				for(k = j+1; k < nbFwds; k++)
				{
					rhok	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(k + idx);
					rhojk	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);

					fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
					
					if(fabs(fac) < K_DOUBLE_TOL)
					{
						rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
						fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
					}

					(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk - rhoj*rhok) / fac;
				}
			}
		}

		fromTime = toTime;
	}
}

void ARM_SVMMSpread::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt		= (nextTime - currTime) / K_YEAR_LEN;
	double racdt	= sqrt(dt);
	int factorsNb	= FactorCount();
	int statesNb	= states->size();
	int fwdsNb		= itsNbEffectiveReset;
	int iLast		= fwdsNb-1;
	
	int i, j, n;
	
	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetTimes[iFirst]<nextTime)
	{
		++iFirst;
	}

	//parait débile mais indispensable pour les échéanciers non triviaux!!!
	int iFirstMartingale = iFirst;
	while( itsEndTimes1[iFirstMartingale] < itsEndTimes1[iLast]) iFirstMartingale++;
	
	int aliveFwds	= iLast - iFirst + 1;

	if(aliveFwds <= 0) return;

	const ARM_MatrixVector& modelLocalStdev	= GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols(); //factorsNb > 2 * aliveFwds ? 2 * aliveFwds : factorsNb;
	
	double rdrift = 0., vdrift = 0., adj = 0.;
	double var, spdi;
	double vvol = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(currTime);
	double kappa = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
	double theta = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
	double scale;
	double E_var, E2_var, V_var, wE2_var = 1. + 0.5 * vvol * vvol / (kappa * theta);;
	double adjt = exp(- kappa * dt);

	double rhov, rhov_;

	/************************************
	/*
	/*			old version mc simple euler scheme
	/*
	/************************************/

	if(itsSimpleEulerScheme)
	{
		std::vector<double> x(aliveFwds + 1);
	
		for(n = 0; n < statesNb; n++)
		{
			// corrélation des gaussiennes
			
			for(i = 0; i < aliveFwds + 1; i++)
			{
				x[i] = 0.;
				for(j = 0; j < eigensNb; j++)
				{
					x[i] += states->GetNumMethodState(n, j) * (*modelLocalStdev[timeIndex])(i,j);
				}
			}
			
			var		= states->GetModelState(n, fwdsNb);

			// simulation des taux
			for(i = iFirst; i <= iLast; i++)
			{
				rhov	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(i);
				rhov_	= sqrt(1. - rhov*rhov);

				spdi	= states->GetModelState(n, i);

				scale	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * racdt;
				spdi	+= scale * sqrt(var) * x[i - iFirst];

				states->SetModelState(n, i, spdi);
			}

			// simulation de la variance
			E_var	= var * adjt + theta * (1. - adjt);
			E2_var	= wE2_var * E_var * E_var + (1. - wE2_var) * adjt * adjt * var * var;
			V_var	= log(E2_var / (E_var * E_var));
			var		= V_var > 0. ? E_var * exp(-0.5 * V_var + sqrt(V_var) * x[aliveFwds]) : E_var;

			states->SetModelState(n, fwdsNb, var);
		}
	}
	else
	{

	/***************************************************
	/*
	/*		new version : andreasen scheme (Quadratic Exponential)
	/*
	/***************************************************/

		// corrélation des gaussiennes
		ARM_GP_Matrix x(statesNb, aliveFwds + 1);

		for(n = 0; n < statesNb; n++)
		{
			for(i = 0; i < aliveFwds; i++)
			{
				x(n,i) = 0.;
				for(j = 0; j < eigensNb; j++) x(n,i) += states->GetNumMethodState(n, j) * (*modelLocalStdev[timeIndex])(i,j);
			}
			x(n, aliveFwds) = states->GetNumMethodState(n, eigensNb);
		}
		
		/*
		for(n = 0; n < statesNb; n++)
		{
			for(i = 0; i < aliveFwds; i++)
			{
				x(n,i) = 0.;
				for(j = 0; j < eigensNb; j++) x(n,i) += states->GetNumMethodState(n, j + 1) * (*modelLocalStdev[timeIndex])(i,j);
			}
			x(n, aliveFwds) = states->GetNumMethodState(n, eigensNb + 1);
		}
		*/

		MCFromToNextTime(states, timeIndex, x);
	}
}

void ARM_SVMMSpread::MCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex, const ARM_GP_Matrix& x) const
{
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt		= (nextTime - currTime) / K_YEAR_LEN;
	double racdt	= sqrt(dt);
	int statesNb	= states->size();
	int fwdsNb		= itsNbEffectiveReset;
	int modelNb		= GetModelNb();
	int iLast		= fwdsNb-1;

	int i, n;
	
	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetTimes[iFirst]<nextTime)
	{
		++iFirst;
	}

	//parait débile mais indispensable pour les échéanciers non triviaux!!!
	int iFirstMartingale = iFirst;
	while( itsEndTimes1[iFirstMartingale] < itsEndTimes1[iLast]) iFirstMartingale++;
	
	int aliveFwds = iLast - iFirst + 1;

	double var, spdi, rhoi;
	double vvol = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(nextTime);
	double kappa = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
	double theta = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
	double w1E2_var = 0.5 * theta * vvol * vvol * SQR(1. - exp(-kappa*dt));
	double w2E2_var = vvol * vvol * exp(-kappa*dt)*(1.-exp(-kappa*dt))/kappa;
	double E_var, E2_var, psi, scale, tipmo, a_var, b_var, invbeta, p_var, u_var;
	double varf;
	double K0, K1, K2, K3, K4, w1 = 0.5, w2 = 0.5;
	double adjt = exp(- kappa * dt);

	for(n = 0; n < statesNb; n++)
	{
		var		= states->GetModelState(n, fwdsNb + modelNb);

		// simulation de la variance
		
		E_var		= var * adjt + theta * (1. - adjt);
		w1E2_var	= 0.5 * theta * vvol * vvol * SQR(1. - adjt);
		w2E2_var	= vvol * vvol * adjt * (1.- adjt) / kappa;
		E2_var		= w1E2_var + w2E2_var * var;
		
		psi			= E2_var / (E_var*E_var);

		if(psi < 1.5)
		{
			tipmo	= 2./psi - 1.;
			b_var	= sqrt(tipmo + sqrt(tipmo * tipmo + tipmo));
			a_var	= E_var / (1. + b_var*b_var);
			varf	= a_var * (b_var + x(n,aliveFwds)) * (b_var + x(n,aliveFwds));
		}
		else
		{
			invbeta	= 0.5 * E_var * (psi + 1.);
			p_var	= (psi - 1.)/(psi + 1.);
			u_var	= NormalCDF(x(n,aliveFwds));
			varf	= u_var < p_var ? 0. : invbeta * log((1.- p_var)/(1.- u_var));
		}

		states->SetModelState(n, fwdsNb + modelNb, varf);

		// simulation des taux
		for(i = iFirst; i <= iLast; i++)
		{	
			spdi	= states->GetModelState(n, i + modelNb);
			
			scale	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i);

			rhoi	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(i);

			K0	= - rhoi * kappa * theta * dt * scale / vvol;
			K1	= w1 * dt * (kappa * rhoi * scale / vvol) - rhoi * scale/ vvol;
			K2	= w1 * dt * (kappa * rhoi * scale / vvol) + rhoi * scale/ vvol;
			K3	= w1 * dt * (1. - rhoi*rhoi) * scale * scale;
			K4	= w2 * dt * (1. - rhoi*rhoi) * scale * scale;

			spdi	+= K0 + K1 * var + K2 * varf + sqrt(K3 * var + K4 * varf) * x(n,i - iFirst);

			states->SetModelState(n, i + modelNb, spdi);
		}
	}
}

void ARM_SVMMSpread::setNumericalModelFitter(ARM_NumericalModelFitter * numericalModelFitter)
{
	if(numericalModelFitter)
	{
		double asof = GetAsOfDate().GetJulian();

		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		int k, size = 0;
		for(k = densityVector.size()-1; k >= 0; k--)
		{
			if(fabs(densityVector[k]->getWeight()) > K_DOUBLE_TOL) break;
		}
		size = k+1;
		
		if(size == 0)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : nul weight for each density" );
		}

		itsResetTimes.resize(size); 
		itsStartTimes.resize(size); 
		itsEndTimes1.resize(size); 
		itsEndTimes2.resize(size); 
		itsFwdRate.resize(size);

		double v0, kappa, theta, vvol;
		std::vector<double> rhos(size), levels(size), vvols(size);	

		ARM_VanillaSecurityDensitySpread * vanillasec = dynamic_cast<ARM_VanillaSecurityDensitySpread*>(&*densityVector[0]);

		if(vanillasec == NULL)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : vanilla security must of spread type");
		}

		ARM_NormalHestonDensityFunctor* density = dynamic_cast<ARM_NormalHestonDensityFunctor*>(&*(densityVector[0]->getDensityFunctor()));

		if(density == NULL)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : density must be of normal heston type");
		}

		v0					= density->GetV0();
		kappa				= density->GetKappa();
		theta				= density->GetTheta();
		vvols[0]			= density->GetVVol();
		rhos[0]				= density->GetRho();
		levels[0]			= density->GetLevel();
		itsFwdRate[0]		= density->GetFwd();
		itsResetTimes[0]	= densityVector[0]->getResetDate() - asof;
		itsStartTimes[0]	= densityVector[0]->getStartDate() - asof;
		itsEndTimes1[0]		= vanillasec->getEndDate1() - asof;
		itsEndTimes2[0]		= vanillasec->getEndDate2() - asof;

		for(k = 1; k < size; k++)
		{
			if((vanillasec = dynamic_cast<ARM_VanillaSecurityDensitySpread*>(&*densityVector[k])) == NULL)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : vanilla security must of spread type");
			}
			if((density = dynamic_cast<ARM_NormalHestonDensityFunctor*>(&*(densityVector[k]->getDensityFunctor()))) == NULL)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : density must be of normal heston type");
			}

			if(fabs(density->GetV0() - v0) > K_DOUBLE_TOL)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : all densities must have same initial var");
			}
			if(fabs(density->GetKappa() - kappa) > K_DOUBLE_TOL)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : all densities must have same var mean reversion");
			}
			if(fabs(density->GetTheta() - theta) > K_DOUBLE_TOL)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : all densities must have same long term var");
			}
			//if(fabs(density->GetVVol() - vvol) > K_DOUBLE_TOL)
			//{
			//	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread : all densities must have same vol of var");
			//}

			rhos[k]				= density->GetRho();
			levels[k]			= density->GetLevel();
			vvols[k]			= density->GetVVol();
			itsFwdRate[k]		= density->GetFwd();
			itsResetTimes[k]	= vanillasec->getResetDate() - asof;
			itsStartTimes[k]	= vanillasec->getStartDate() - asof;
			itsEndTimes1[k]		= vanillasec->getEndDate1() - asof;
			itsEndTimes2[k]		= vanillasec->getEndDate2() - asof;
		}

		((ARM_ModelParamsBGMSV1F*) GetModelParams())->SetParams(itsResetTimes, v0, kappa, theta, std::vector<double>(1,0.), levels, vvols, rhos);

//		((ARM_ModelParamsBGMSV1F*) GetModelParams())->UpdateParamValues(&itsResetTimes, NULL, &levels, &rhos, &std::vector<double>(1, vvol),
//				&std::vector<double>(1,theta), v0);

		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);
	}
}

ARM_VectorPtr ARM_SVMMSpread::DiscountFactor(const string& curveName, double evalTime, double maturityTime, 
										const ARM_PricingStatesPtr& states) const
{
	if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
// FIXMEFRED: mig.vc8 (30/05/2007 16:18:15):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
	else
		return static_cast<ARM_VectorPtr>(new std::vector<double>(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
}

ARM_VectorPtr ARM_SVMMSpread::VanillaCaplet(
		const string& curveName, 
		double evalTime,
		double payTime, 
		double period,
        double payNotional,
		double fwdResetTime, 
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
		const std::vector<double>& strikesPerState,
        int capFloor,
		const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::VanillaCaplet not implemented" );
}

ARM_VectorPtr ARM_SVMMSpread::VanillaSwaption(
		const string& curveName,
		double evalTime,
		double swapResetTime,
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& floatResetTimes,
		const std::vector<double>& floatStartTimes,
		const std::vector<double>& floatEndTimes,
		const std::vector<double>& floatIntTerms,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::VanillaSwaption not implemented" );
}

ARM_VectorPtr ARM_SVMMSpread::Spread(
		const string& curveName, 
		double evalTime,
		double coeff1,
		double floatStartTime1,
		double floatEndTime1, 
		const std::vector<double>& fixPayTimes1,
		const std::vector<double>& fixPayPeriods1,
		const std::vector<double>& fwdStartTimes1,
		const std::vector<double>& fwdEndTimes1,
		const std::vector<double>& fwdPayPeriods1, 
		const std::vector<double>& floatPayTimes1,
		const std::vector<double>& floatPayPeriods1,
		const std::vector<double>& margin1,
		double coeff2,
		double floatStartTime2,
		double floatEndTime2, 
		const std::vector<double>& fixPayTimes2,
		const std::vector<double>& fixPayPeriods2,
		const std::vector<double>& fwdStartTimes2,
		const std::vector<double>& fwdEndTimes2,
		const std::vector<double>& fwdPayPeriods2, 
		const std::vector<double>& floatPayTimes2,
		const std::vector<double>& floatPayPeriods2,
		const std::vector<double>& margin2,
		const ARM_PricingStatesPtr& states) const
{
	int modelNb = GetModelNb();

	if(fabs(floatStartTime1 - fwdStartTimes1[0]) < 7
	&& fabs(floatEndTime1 - fwdEndTimes1[fwdEndTimes1.size()-1]) < 7
	&& fabs(floatStartTime2 - fwdStartTimes2[0]) < 7
	&& fabs(floatEndTime2 - fwdEndTimes2[fwdEndTimes2.size()-1]) < 7
	&& fabs(floatStartTime1 - floatStartTime2) < 7
	&& fabs(coeff1 - coeff2) < K_DOUBLE_TOL)
	{
		ARM_VectorPtr result = cmsSpread(curveName, evalTime, floatStartTime1, floatEndTime1, floatEndTime2, states);

		if(result != ARM_VectorPtr(NULL)) return result;
	}

	return ARM_PricingModelIR::Spread(curveName, evalTime, coeff1, floatStartTime1, floatEndTime1, fixPayTimes1, fixPayPeriods1,
		fwdStartTimes1, fwdEndTimes1, fwdPayPeriods1, floatPayTimes1, floatPayPeriods1, margin1,
		coeff2,
		floatStartTime2, floatEndTime2, fixPayTimes2, fixPayPeriods2, fwdStartTimes2, fwdEndTimes2, fwdPayPeriods2, 
		floatPayTimes2, floatPayPeriods2, margin2, states);
}

ARM_VectorPtr ARM_SVMMSpread::VanillaSpreadOptionLet(const string& curveName,
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
	const std::vector<double>& strikes,
	double swapLongFloatStartTime,
	double swapLongFloatEndTime,
	const std::vector<double>& swapLongFixPayTimes,
	const std::vector<double>& swapLongFixPayPeriods,
	double swapShortFloatStartTime,
	double swapShortFloatEndTime,
	const std::vector<double>& swapShortFixPayTimes,
	const std::vector<double>& swapShortFixPayPeriods,
	const ARM_PricingStatesPtr& states) const
{
	if(fabs(evalTime - swapLongFloatStartTime) < 7
	&& fabs(swapLongFloatStartTime - swapShortFloatStartTime) < 7)
	{
		ARM_VectorPtr spread = cmsSpread(curveName, evalTime, swapLongFloatStartTime, swapLongFloatEndTime, swapShortFloatEndTime, states);

        int i, stateSize = states->size();

        ARM_VectorPtr dfPay = DiscountFactor(curveName, evalTime, payTime, states);

        ARM_GP_VectorPtr result( new std::vector<double>( stateSize ) );

        double payoff;

        for (i = 0; i<stateSize; i++)
        {
			payoff = callPut * ( (*spread)[i] - strikes[i] );

            (*result)[i] = (payoff > 0) ? payPeriod * notional * (*dfPay)[i] * payoff : 0.0;

		}

		return result;
	}

	if(evalTime > K_DOUBLE_TOL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::VanillaSpreadOptionLet implemented only for spot date" );
	}

	if(fabs(swapLongFloatStartTime - swapShortFloatStartTime) < 7)
	{
		int idx = IdxFromValue(itsStartTimes, swapLongFloatStartTime,7);

		if(fabs(itsEndTimes1[idx] - itsStartTimes[idx] + swapLongFloatStartTime - swapLongFloatEndTime) < 7
		&& fabs(itsEndTimes2[idx] - itsStartTimes[idx] + swapLongFloatStartTime - swapShortFloatEndTime) < 7)
		{

			double v0		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();
			double kappa	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
			double theta	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
			double vvol		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(idx);
			double rho		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx);
			double level	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx);

			double spreadopt = NormalHeston(rho, kappa, theta * SQR(level), vvol * level, v0 * SQR(level), 
							itsFwdRate[idx], strikes[0], swapLongFloatStartTime / 365., callPut, -1.);

			ARM_VectorPtr result(new std::vector<double>(1, spreadopt * notional * GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN)));

			return result;
		}
		else if(fabs(itsEndTimes1[idx] - itsStartTimes[idx] + swapLongFloatStartTime - swapShortFloatEndTime) < 7
		&& fabs(itsEndTimes2[idx] - itsStartTimes[idx] + swapLongFloatStartTime - swapLongFloatEndTime) < 7)
		{
			double v0		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();
			double kappa	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
			double theta	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
			double vvol		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(idx);
			double rho		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx);
			double level	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx);

			double spreadopt = NormalHeston(rho, kappa, theta * SQR(level), vvol * level, v0 * SQR(level), 
							itsFwdRate[idx], strikes[0], swapLongFloatStartTime / 365., -callPut, -1.);

			ARM_VectorPtr result(new std::vector<double>(1, spreadopt * notional * GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN)));

			return result;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVMMSpread::VanillaSpreadOptionLet implemented only for calibrated Spread" );
}

ARM_VectorPtr ARM_SVMMSpread::cmsSpread(const string& curveName, double evalTime, double startTime, double endTime1, double endTime2,
										const ARM_PricingStatesPtr& states) const
{
	int modelNb = GetModelNb();

	int idx = IdxFromValue(itsStartTimes, startTime,7);

	if(idx > -1)
	{
		if(fabs(itsEndTimes1[idx] - itsStartTimes[idx] + startTime - endTime1) < 7
		&& fabs(itsEndTimes2[idx] - itsStartTimes[idx] + startTime - endTime2) < 7)
		{
			std::vector<double>& result = new std::vector<double>( states->size() );
			
			for(int i = 0; i < states->size(); i++) 
				(*result)[i] = states->GetModelState(i, idx + modelNb);

			return static_cast<ARM_VectorPtr>(result);
		}
		else if(fabs(itsEndTimes1[idx] - itsStartTimes[idx] + startTime - endTime2) < 7
		&& fabs(itsEndTimes2[idx] - itsStartTimes[idx] + startTime - endTime1) < 7)
		{
			std::vector<double>& result = new std::vector<double>( states->size() );

			for(int i = 0; i < states->size(); i++)
				(*result)[i] = - states->GetModelState(i, idx + modelNb);

			return static_cast<ARM_VectorPtr>(result);
		}
	}
	else
	{
		int infidx = IndexOfLastLowerEqInVector_DefaultFirst(startTime, itsStartTimes);
		int supidx = IndexOfFirstHigherEqInVector_DefaultLast(startTime, itsStartTimes);

		double length1 = endTime1 - startTime;
		double length2 = endTime2 - startTime;

		if(fabs(itsEndTimes1[infidx] - itsStartTimes[infidx] - length1) < 7
		&& fabs(itsEndTimes1[supidx] - itsStartTimes[supidx] - length1) < 7
		&& fabs(itsEndTimes2[infidx] - itsStartTimes[infidx] - length2) < 7
		&& fabs(itsEndTimes2[supidx] - itsStartTimes[supidx] - length2) < 7)
		{
			double w = (startTime - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

			std::vector<double>& result = new std::vector<double>( states->size() );

			for(int i = 0; i < states->size(); i++)
			{
				
				(*result)[i] = states->GetModelState(i, infidx + modelNb) 
					+ w * (states->GetModelState(i, supidx + modelNb) - states->GetModelState(i, infidx + modelNb));
			}

			return static_cast<ARM_VectorPtr>(result);
		}
		else if(fabs(itsEndTimes1[infidx] - itsStartTimes[infidx] - length2) < 7
		&& fabs(itsEndTimes1[supidx] - itsStartTimes[supidx] - length2) < 7
		&& fabs(itsEndTimes2[infidx] - itsStartTimes[infidx] - length1) < 7
		&& fabs(itsEndTimes2[supidx] - itsStartTimes[supidx] - length1) < 7)
		{
			double w = (startTime - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

			std::vector<double>& result = new std::vector<double>( states->size() );

			for(int i = 0; i < states->size(); i++)
			{
				
				(*result)[i] = - states->GetModelState(i, infidx + modelNb) 
					- w * (states->GetModelState(i, supidx + modelNb) - states->GetModelState(i, infidx + modelNb));
			}

			return static_cast<ARM_VectorPtr>(result);
		}
	}

	return ARM_VectorPtr(NULL);
}

std::vector<double> ARM_SVMMSpread::GetATMGaussianVol()
{
	int k, size = itsResetTimes.size();

	std::vector<double> gvol(size);

	double v0, kappa, theta, rho, vvol, level, fwd, price;

	for(k = 0; k < size; k++)
	{
		v0		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();
		kappa	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
		theta	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
		vvol	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(k);
		rho		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(k);
		level	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(k);
		fwd		= itsFwdRate[k];

		price	= NormalHeston(rho, kappa, theta * SQR(level), vvol * level, v0 * SQR(level), 
							fwd, fwd, itsResetTimes[k] / 365., -1, -1.);

		gvol[k]	= VanillaImpliedVol_N(fwd, price, fwd, itsResetTimes[k] / 365., -1);
	}

	return gvol;
}

CC_END_NAMESPACE()
