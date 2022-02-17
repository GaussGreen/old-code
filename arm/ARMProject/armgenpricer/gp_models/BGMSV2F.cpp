/// gpmodels
#include "gpmodels/bgmsv2f.h"
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
#include "gpclosedforms/heston_pricer.h"

CC_BEGIN_NAMESPACE( ARM )

ARM_BGMSV2F::ARM_BGMSV2F( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, bool AllowInterpol, bool ComputeDrift, bool Proxy)
:	ARM_PricingModelIR(zc,params), 
	itsAllowInterpol(AllowInterpol),
	itsResetTimes(NULL),
	itsStartTimes(NULL),
	itsEndTimes(NULL),
	itsCalibrationStatus(true),
	itsCalibratedStatus(false),
	itsWeight(NULL),
	itsDelta(NULL),
	itsFwdRate(NULL),
	itsATMVols(NULL),
	itsFromShiftedToRate(NULL),
	itsComputeDriftAdj(ComputeDrift),
	itsEigenValues(0),
	itsModelLocalRealVar(0),
	itsCalibSecDensities(0),
	itsProxyStatus(Proxy),
	itsRandGen(-156)
{
	itsSpreadCalib = false;
	if(itsProxyStatus) itsComputeDriftAdj = false;
	//itsAllowInterpol = true;
	itsInterpolZC = true;
	itsPrevEvalTime = 0.;
}


////////////////////////////////////////////////////
///	Class  : ARM_BGMSV2F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BGMSV2F::ARM_BGMSV2F(const ARM_BGMSV2F& rhs)
:	ARM_PricingModelIR(rhs),
	itsAllowInterpol(rhs.itsAllowInterpol),
	itsResetTimes(rhs.itsResetTimes),
	itsStartTimes(rhs.itsStartTimes),
	itsEndTimes(rhs.itsEndTimes),
	itsCalibrationStatus(rhs.itsCalibrationStatus),
	itsCalibratedStatus(rhs.itsCalibratedStatus),
	itsWeight(rhs.itsWeight),
	itsDelta(rhs.itsDelta),
	itsFwdRate(rhs.itsFwdRate),
	itsATMVols(rhs.itsATMVols),
	itsFromShiftedToRate(rhs.itsFromShiftedToRate),
	itsComputeDriftAdj(rhs.itsComputeDriftAdj),
	itsProxyStatus(rhs.itsProxyStatus),
	itsSpreadCalib(rhs.itsSpreadCalib),
	itsCalibWeight(rhs.itsCalibWeight),
	itsInterpolZC(rhs.itsInterpolZC),
	itsRandGen(rhs.itsRandGen),
	itsPrevEvalTime(rhs.itsPrevEvalTime)
{
	DuplicateCloneablePointorVectorInPlace<std::vector<double>>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelLocalRealVar, itsModelLocalRealVar);
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsCalibSecDensities, itsCalibSecDensities);
}	


ARM_BGMSV2F& ARM_BGMSV2F::operator = (const ARM_BGMSV2F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_BGMSV2F();

		new (this) ARM_BGMSV2F (rhs);
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_BGMSV2F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_BGMSV2F::~ARM_BGMSV2F()
{
	DeletePointorVector<std::vector<double>>( itsEigenValues );
	DeletePointorVector<ARM_GP_Matrix>( itsModelLocalRealVar );
}

bool ARM_BGMSV2F::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::ValidateModelParams: beta correlation structure is needed for BGMSV1F" );

	return true;
}

void ARM_BGMSV2F::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

void ARM_BGMSV2F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);

	calibMethod.DefaultValidateWithModel(*this);
}

void ARM_BGMSV2F::PreProcessing(ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

void ARM_BGMSV2F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	// TO DO 
}


void ARM_BGMSV2F::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

ARM_PricingStatesPtr ARM_BGMSV2F::FirstPricingStates(size_t bucketSize) const
{
	const size_t nbPayoffs			= 0;
	size_t nbModelStates			= itsNbEffectiveReset + 2;
	size_t factorsNb				= FactorCount();

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	for(int i = 0; i < itsNbEffectiveReset; i++)
	{
		// initialisation des taux shiftés
		double initRateValue	= itsFwdRate[i] / ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(i);

		for(int k = 0; k < bucketSize; k++) 
		{
			initStates->SetModelState(k, i, initRateValue);
		}
	}

	double v01 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV01();
	double v02 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV02();

	for(int k = 0; k < bucketSize; k++)
	{
		initStates->SetModelState(k, itsNbEffectiveReset, v01);
		initStates->SetModelState(k, itsNbEffectiveReset+1, v02);
	}

	return initStates;
}

ARM_PricingStatesPtr ARM_BGMSV2F::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	if(itsCalibrationStatus)
	{
		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}
	else
	{
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
}

ARM_PricingStatesPtr ARM_BGMSV2F::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	if(itsCalibrationStatus)
	{	
		if(itsCalibratedStatus == false)
		{
			((ARM_ModelParamsBGMSV2F*) GetModelParams())->Calibrate(itsResetTimes, itsFwdRate, itsCalibSecDensities);

			itsFromShiftedToRate.resize(itsResetTimes.size());

			// pour passer des libors shiftés aux libors : fwd x (shift - 1) / shift
			for(int k = 0; k < itsResetTimes.size(); k++) 
			{
				itsFromShiftedToRate[k] = (((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(k) - 1.) * itsFwdRate[k] / ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(k);
			}

			itsCalibratedStatus = true;
		}

		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}

	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );

	return newStates;
}

ARM_VectorPtr ARM_BGMSV2F::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,getTerminalTime()));
}

std::vector<double>& ARM_BGMSV2F::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetTimes.Clone());
}


void ARM_BGMSV2F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
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


string ARM_BGMSV2F::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "\n\n";
    os << indent << "SVBGM2F : Shifted Stochastic Volatility BGM with two processes for variance (Heston time dependant) \n";
    os << indent << "----------------------------\n";
	os << ARM_PricingModel::toString(indent);
	os << "\n\n";

	return os.str();
}

///////////////////////////
//
//		Monte Carlo
//
///////////////////////////

void ARM_BGMSV2F::CalcNbEffectiveReset(const std::vector<double>& timeSteps)
{
	if(itsProxyStatus)
	{
		double lastTimeStep = timeSteps[timeSteps.size() - 1];
		itsNbEffectiveReset = 0;
		while(!(itsResetTimes[itsNbEffectiveReset] > lastTimeStep))
		{
			itsNbEffectiveReset ++;
			if(itsNbEffectiveReset == itsResetTimes.size()) break;
		}
	}
	else
		itsNbEffectiveReset = itsResetTimes.size();

	((ARM_ModelParamsBGMSV2F*) GetModelParams())->checkFactorCount(itsNbEffectiveReset);
}

void ARM_BGMSV2F::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	CalcNbEffectiveReset(timeSteps);

	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances( timeSteps, auxLocalCov);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);

	ModelStateLocalStdDev(timeSteps, auxLocalCov);
}

void ARM_BGMSV2F::NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{
//	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	
	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		int factorsNb = itsEigenValues[i]->size() + 2 ;

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		for (int k = 0; k < factorsNb-2; k++)
			(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
		

		(*localVariances[i])(factorsNb-2,factorsNb-2) = 1.;
		(*localVariances[i])(factorsNb-1,factorsNb-1) = 1.;
	}
}

void ARM_BGMSV2F::NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (int i=0;i<timeStepsSize;i++)
	{
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(2*totalFwds+2,1.0);
	}
}


void ARM_BGMSV2F::ModelStateLocalStdDev( const std::vector<double>& timeSteps, const ARM_MatrixVector& auxLocalCov)
{
	ARM_MatrixVector modStateLocalVars;
	
	int totalFwds		= itsNbEffectiveReset;
	int factorsNb		= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->EffFactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;

	int i, j, k, startTimePos = 0;

	modStateLocalVars.resize((timeStepsSize-1)*(modelNb+1));

	DeletePointorVector<std::vector<double>>( itsEigenValues );
	itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

	int nbFwds, idx = 0, currfactorsNb = factorsNb;
	int maxfactorsNb = 0;
	double minratio = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetMinRatio();
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

		int esize = 2*nbFwds;

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

	((ARM_ModelParamsBGMSV2F*) GetModelParams())->SetFactorCount(maxfactorsNb);

	idx = 0;
	for(i = 0; i < timeStepsSize - 1 ;++i)
	{
		/// initalize the toTime
		toTime  = timeSteps[i+1];

		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		int esize = 2*nbFwds;

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

void ARM_BGMSV2F::NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_BGMSV2F::ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;


#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::ModelStateLocalVariances: first time step != 0" );
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

		localVariances[i] = new ARM_GP_Matrix(2*nbFwds, 2*nbFwds);	

		double rhoj, rhok, rhojk, fac;
		for(j = 0; j < nbFwds; j++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(j + idx);

			(*localVariances[i])(j,j) = 1.;

			for(k = j+1; k < nbFwds; k++)
			{
				rhok	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(k + idx);
				rhojk	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);

				//fac		= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
				fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				
				if(fabs(fac) < K_DOUBLE_TOL)
				{
					rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
					//fac	= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
					fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				}

				(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk - rhoj*rhok) / fac;
			}
		}
		for(j = 0; j < nbFwds; j++)
		{
			rhoj	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(j + idx);

			(*localVariances[i])(j+nbFwds,j+nbFwds) = 1.;

			for(k = j+1; k < nbFwds; k++)
			{
				rhok	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(k + idx);
				rhojk	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);

				//fac		= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
				fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				
				if(fabs(fac) < K_DOUBLE_TOL)
				{
					rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
					//fac	= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
					fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
				}

				(*localVariances[i])(j+nbFwds,k+nbFwds) = (*localVariances[i])(k+nbFwds,j+nbFwds) = (rhojk - rhoj*rhok) / fac;
			}
		}

		fromTime = toTime;
	}
}

void ARM_BGMSV2F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
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
	if (itsAllowInterpol && iFirst>0) iFirst--;

	//parait débile mais indispensable pour les échéanciers non triviaux!!!
	int iFirstMartingale = iFirst;
	while( itsEndTimes[iFirstMartingale] < itsEndTimes[iLast]) iFirstMartingale++;
	
	int aliveFwds	= iLast - iFirst + 1;

	if(aliveFwds <= 0) return;

	const ARM_MatrixVector& modelLocalStdev	= GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols(); //factorsNb > 2 * aliveFwds ? 2 * aliveFwds : factorsNb;
	

	/***************************************************
	/*
	/*		new version : andreasen scheme (Quadratic Exponential)
	/*
	/***************************************************/

	// corrélation des gaussiennes
	ARM_GP_Matrix x(statesNb, 2*aliveFwds + 2);

	for(n = 0; n < statesNb; n++)
	{
		for(i = 0; i < 2*aliveFwds; i++)
		{
			x(n,i) = 0.;
			for(j = 0; j < eigensNb; j++) x(n,i) += states->GetNumMethodState(n, j) * (*modelLocalStdev[timeIndex])(i,j);
		}

		x(n, 2*aliveFwds) = states->GetNumMethodState(n, eigensNb);
		x(n, 2*aliveFwds+1) = states->GetNumMethodState(n, eigensNb+1);
	}

	MCFromToNextTime(states, timeIndex, x);
}

void ARM_BGMSV2F::MCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex, const ARM_GP_Matrix& x) const
{
	double currTime = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt		= (nextTime - currTime) / K_YEAR_LEN;
	double racdt	= sqrt(dt);
	int statesNb	= states->size();
	int fwdsNb		= itsNbEffectiveReset;
	int modelNb		= GetModelNb();
	int iLast		= fwdsNb-1;

	int i, j, n;
	
	int iFirst  = 0;
	while( iFirst < fwdsNb 	&& itsResetTimes[iFirst]<nextTime)
	{
		++iFirst;
	}
	if (itsAllowInterpol && iFirst>0) iFirst--;

	//parait débile mais indispensable pour les échéanciers non triviaux!!!
	int iFirstMartingale = iFirst;
	while( itsEndTimes[iFirstMartingale] < itsEndTimes[iLast]) iFirstMartingale++;
	
	int aliveFwds = iLast - iFirst + 1;

	double rcov, vcov1, vcov2;
	double rdrift = 0., vdrift1 = 0., vdrift2, adj = 0.;
	double var1, var2, libi, rhoi1, rhoi2;

	double vvol1 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol1();
	double kappa1 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa1();
	double theta1 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta1();
	
	double vvol2 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol2();
	double kappa2 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa2();
	double theta2 = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta2();

	double w1E2_var, w2E2_var;
	double E_var, E2_var, psi, adjt, scale, tipmo, a_var, b_var, invbeta, p_var, u_var;
	double akappa, atheta;
	double varf1, varf2;
	double K01, K11, K21, K31, K41, K02, K12, K22, K32, K42, w1 = 0.5, w2 = 0.5;

	for(n = 0; n < statesNb; n++)
	{
		var1	= states->GetModelState(n, fwdsNb + modelNb);
		var2	= states->GetModelState(n, fwdsNb + 1 + modelNb);

		// simulation de la variance
		vdrift1 = vdrift2 = 1.;

		if(itsComputeDriftAdj)
		{
			for(i = iFirst; i <= iLast; i++)
			{
				if(IsOnSamePath(i,iFirst) == false) continue;

				libi	= states->GetModelState(n, i + modelNb);
				vcov1	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(i);
				vcov2	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(i);
				
				adj		= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(i) * states->GetModelState(n, i + modelNb) * itsDelta[i]
						/ (1. + itsDelta[i] * (states->GetModelState(n, i) + itsFromShiftedToRate[i]));

				vdrift1	+= adj * vcov1 * vvol1 / kappa1;
				vdrift2	+= adj * vcov2 * vvol2 / kappa2;
			}
		}
		
		akappa		= kappa1 * vdrift1;
		atheta		= theta1 / vdrift1;
		adjt		= exp(- akappa * dt);
		E_var		= var1 * adjt + atheta * (1. - adjt);
		w1E2_var	= 0.5 * atheta * vvol1 * vvol1 * SQR(1. - adjt);
		w2E2_var	= fabs(akappa) < K_DOUBLE_TOL ? vvol1*vvol1*dt : vvol1 * vvol1 * adjt * (1.- adjt) / akappa;
		E2_var		= w1E2_var + w2E2_var * var1;
		
		psi			= E2_var / (E_var*E_var);

		if(psi < 1.5)
		{
			tipmo	= fabs(psi) < K_DOUBLE_TOL ? 0. : 2./psi - 1.;
			b_var	= sqrt(tipmo + sqrt(tipmo * tipmo + tipmo));
			a_var	= E_var / (1. + b_var*b_var);
			varf1	= a_var * (b_var + x(n,2*aliveFwds)) * (b_var + x(n,2*aliveFwds));
		}
		else
		{
			invbeta	= 0.5 * E_var * (psi + 1.);
			p_var	= (psi - 1.)/(psi + 1.);
			u_var	= NormalCDF(x(n,2*aliveFwds));
			varf1	= u_var < p_var ? 0. : invbeta * log((1.- p_var)/(1.- u_var));
		}

		states->SetModelState(n, fwdsNb + modelNb, varf1);

		akappa		= kappa2 * vdrift2;
		atheta		= theta2 / vdrift2;
		adjt		= exp(- akappa * dt);
		E_var		= var2 * adjt + atheta * (1. - adjt);
		w1E2_var	= 0.5 * atheta * vvol2 * vvol2 * SQR(1. - adjt);
		w2E2_var	= fabs(akappa) < K_DOUBLE_TOL ? vvol2*vvol2*dt : vvol2 * vvol2 * adjt * (1.- adjt) / akappa;
		E2_var		= w1E2_var + w2E2_var * var2;
		
		psi			= E2_var / (E_var*E_var);

		if(psi < 1.5)
		{
			tipmo	= fabs(psi) < K_DOUBLE_TOL ? 0. : 2./psi - 1.;
			b_var	= sqrt(tipmo + sqrt(tipmo * tipmo + tipmo));
			a_var	= E_var / (1. + b_var*b_var);
			varf2	= a_var * (b_var + x(n,2*aliveFwds+1)) * (b_var + x(n,2*aliveFwds+1));
		}
		else
		{
			invbeta	= 0.5 * E_var * (psi + 1.);
			p_var	= (psi - 1.)/(psi + 1.);
			u_var	= NormalCDF(x(n,2*aliveFwds+1));
			varf2	= u_var < p_var ? 0. : invbeta * log((1.- p_var)/(1.- u_var));
		}

		states->SetModelState(n, fwdsNb + 1 + modelNb, varf2);


		// simulation des taux
		for(i = iFirst; i <= iLast; i++)
		{
			rdrift	= 0.;
			
			libi	= states->GetModelState(n, i + modelNb);

			if(itsComputeDriftAdj)
			{
				for(j = i + 1; j <= iLast; j++)
				{
					if(IsOnSamePath(i,j) == false) continue;

					rcov = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->RateRateCorrel(nextTime, itsResetTimes[i],i,itsResetTimes[j],j) * dt;

					adj		= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(j) * states->GetModelState(n, j + modelNb) * itsDelta[j]
							/ (1. + itsDelta[j] * (states->GetModelState(n, j + modelNb) + itsFromShiftedToRate[j]));

					rdrift	-= adj * rcov;
				}

				rdrift *= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(i) 
						* ((w1*var1 + w2*varf1) + (w1*var2 + w2*varf2));
			}
			
			scale	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(i);

			rhoi1	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(i);
			rhoi2	= ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(i);

			K01	= - rhoi1 * kappa1 * theta1 * dt * scale / vvol1;
			K11	= w1 * dt * (kappa1 * rhoi1 * scale / vvol1 - 0.5 * scale * scale) - rhoi1 * scale/ vvol1;
			K21	= w1 * dt * (kappa1 * rhoi1 * scale / vvol1 - 0.5 * scale * scale) + rhoi1 * scale/ vvol1;
			K31	= w1 * dt * (1. - rhoi1*rhoi1) * scale * scale;
			K41	= w2 * dt * (1. - rhoi1*rhoi1) * scale * scale;

			K02	= - rhoi2 * kappa2 * theta2 * dt * scale / vvol2;
			K12	= w1 * dt * (kappa2 * rhoi2 * scale / vvol2 - 0.5 * scale * scale) - rhoi2 * scale/ vvol2;
			K22	= w1 * dt * (kappa2 * rhoi2 * scale / vvol2 - 0.5 * scale * scale) + rhoi2 * scale/ vvol2;
			K32	= w1 * dt * (1. - rhoi2*rhoi2) * scale * scale;
			K42	= w2 * dt * (1. - rhoi2*rhoi2) * scale * scale;

			libi	*= exp(rdrift + K01 + K11 * var1 + K21 * varf1 + sqrt(K31 * var1 + K41 * varf1) * x(n,i - iFirst)
						+ K02 + K12 * var2 + K22 * varf2 + sqrt(K32 * var2 + K42 * varf2) * x(n, i - iFirst + aliveFwds));

			states->SetModelState(n, i + modelNb, libi);
		}
	}
}

void ARM_BGMSV2F::setNumericalModelFitter(ARM_NumericalModelFitter * numericalModelFitter)
{
	if(numericalModelFitter)
	{
		double asof = GetAsOfDate().GetJulian();

		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		int i, k, size = 0;
		for(k = densityVector.size()-1; k >= 0; k--)
		{
			if(fabs(densityVector[k]->getWeight()) > K_DOUBLE_TOL) break;
		}
		size = k+1;
		
		if(size == 0)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F : nul weight for each density" );
		}

		std::vector<double> resetDates(size); 
		std::vector<double> startDates(size); 
		std::vector<double> endDates(size); 
		std::vector<double> delta(size); 
		itsFwdRate.resize(size);

		for (k = 0, i = 0; i < densityVector.size() ; i++)
		{
			if(fabs(densityVector[i]->getWeight()) < K_DOUBLE_TOL) continue;

			resetDates[k]			= densityVector[i]->getResetDate();
			startDates[k]			= densityVector[i]->getStartDate();
			endDates[k]				= densityVector[i]->getEndDate();

			delta[k]				= densityVector[i]->getInterestTerms()[0];
			itsFwdRate[k]			= densityVector[i]->ComputeForwardRate();

			if(densityVector[i]->getInterestTerms().size() > 1)
			{
				itsComputeDriftAdj	= false;
				itsProxyStatus		= true;
			}

			k++;
		}
		

		itsResetTimes			= resetDates - asof;
		itsStartTimes			= startDates - asof;
		itsEndTimes				= endDates   - asof;
		itsDelta				= delta;

		computeWeights();

		((ARM_ModelParamsBGMSV2F*) GetModelParams())->checkFactorCount(size);

		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);

		itsCalibratedStatus = false;
	}
}

void ARM_BGMSV2F::computeWeights() 
{
	size_t size		= itsResetTimes.size();
	double terminal	= itsEndTimes[size-1];
	double term;

	itsWeight.resize( size, size);

	for (size_t j=0;j<size;j++)
	{
		size_t k = 0;
		term = itsEndTimes[j];
		while (term<terminal)
		{
			while (itsStartTimes[k]<term)
			{
				itsWeight(j,k) = 0.;
				if(k == size-1) break;
				k++;
			}
			itsWeight(j,k) = 1.;
			term = itsEndTimes[k];
			if(k == size-1) break;
			k++;
		}

		itsWeight(j,j) = 1.;
	}
}

ARM_VectorPtr ARM_BGMSV2F::DiscountFactor(const string& curveName, double evalTime, double maturityTime, 
										const ARM_PricingStatesPtr& states) const
{
	// if(itsAllowInterpol)
	if(itsInterpolZC)
		return ForwardDiscountFactor(curveName,evalTime,evalTime,maturityTime,states);
	else
		return DiscountFactorNoInterpol(curveName,evalTime,maturityTime,states);
}

ARM_VectorPtr ARM_BGMSV2F::ForwardDiscountFactor( const string& curveName, double evalTime, double startTime, double endTime, const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	if( evalTime < K_DOUBLE_TOL || itsProxyStatus == true)
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new std::vector<double>(1,GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new std::vector<double>(states->size(), GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN) ));
	}

	if( endTime > itsEndTimes[itsEndTimes.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: endTime after schedule's last date");

	if( evalTime < itsResetTimes[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime before first reset");

	std::vector<double>& result = new std::vector<double>( states->size(),1.0 );

	std::vector<double>::iterator iter;
	std::vector<double>::iterator iter1;
	std::vector<double>::iterator iter2;

	size_t FirstIdx = IndexOfFirstHigherEqInVector_DefaultLast(startTime, itsStartTimes); //cas où start dans dernière période!
	if (startTime > itsStartTimes[FirstIdx])
		FirstIdx++;

	size_t LastIdx  = IndexOfLastLowerEqInVector(endTime,itsStartTimes);
	
	if ( (LastIdx < FirstIdx) || ( IsOnSamePath(FirstIdx, LastIdx) ))
	{
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, LastIdx, modelNb, states);

		for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
			(*iter) *= (*iter1);
	}
	else
	{
		ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, LastIdx, itsResetTimes.size(), modelNb, states);
		ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, FirstIdx, itsResetTimes.size(), modelNb, states);

		for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
			(*iter) *= (*iter2)/(*iter1);
	}

	if ( (LastIdx < FirstIdx) || (startTime < itsStartTimes[FirstIdx]))
	{
		double mat				= LastIdx < FirstIdx ? endTime : itsStartTimes[FirstIdx];
		double theta			= (mat - startTime)/360.;
		size_t idx				= FirstIdx > 0 ? FirstIdx - 1 : FirstIdx;
		double fwd1				= itsFwdRate[idx];
		double fwd				= (GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)/GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)-1.)/theta;//( 1. - weight2 ) * fwd1 + weight2 * fwd2;
		double coeff			= theta * (fwd / fwd1);

		double mu				= GetZeroCurve()->DiscountPrice(mat/K_YEAR_LEN)
								/ GetZeroCurve()->DiscountPrice(startTime/K_YEAR_LEN)
								* ( 1 + theta * fwd);

		int k;

		for (iter = result->begin(), k = 0; iter != result->end() ; ++iter, k++)
				(*iter) *= mu / ( 1. + coeff * ( states->GetModelState(k, idx + modelNb) + itsFromShiftedToRate[idx]) );
	}

	if ( (endTime > itsStartTimes[LastIdx]) && (LastIdx >= FirstIdx))
	{
		double theta	= (endTime - itsStartTimes[LastIdx])/360.;

		double mu		= GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)
						/ GetZeroCurve()->DiscountPrice(itsStartTimes[LastIdx]/K_YEAR_LEN)
						* ( 1 + theta * itsFwdRate[LastIdx]);
				
		int k;

		for(iter = result->begin(), k = 0 ; iter != result->end() ; ++iter, k++)
			(*iter) *= mu / ( 1. + theta * (states->GetModelState(k, LastIdx + modelNb) + itsFromShiftedToRate[LastIdx]));
	}

	return static_cast<ARM_VectorPtr>(result);
}

ARM_VectorPtr ARM_BGMSV2F::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t IdxFrom, size_t IdxTo , size_t modelNb, const ARM_PricingStatesPtr& states) const
{
	std::vector<double>& result	= new std::vector<double>( states->size(),1.0 );
	
	size_t size = states->size();

	for (size_t k = IdxFrom; k < IdxTo ; k++ )
	{
		if ( itsWeight(IdxFrom,k) == 1. )
		{
			for(size_t i = 0; i < size; i++)
			{
				(*result)[i] /= (1. + itsDelta[k] * (states->GetModelState(i, k + modelNb) + itsFromShiftedToRate[k]));
			}
		}
	}

	return static_cast<ARM_VectorPtr>(result);
}

ARM_VectorPtr ARM_BGMSV2F::DiscountFactorNoInterpol(const string& curveName, double evalTime, double maturityTime, 
												  const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();
	
	if( evalTime < K_DOUBLE_TOL || itsProxyStatus == true)
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new std::vector<double>(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new std::vector<double>(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( maturityTime > itsEndTimes[itsEndTimes.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::DiscountFactor: MaturityTime after schedule's last date");


	std::vector<double>& result = new std::vector<double>( states->size(),1.0 );
	
	if( DoesResetTimeExist( evalTime ) )
	{
		std::vector<double>::iterator iter;
		std::vector<double>::iterator iter1;
		std::vector<double>::iterator iter2;

		size_t ResetIdx = IdxFromValueWithTol( itsResetTimes, evalTime, 1. );
		size_t StartIdx = ResetIdx;

		if ( maturityTime >= itsStartTimes[itsStartTimes.size()-1] )
		{
			StartIdx = itsStartTimes.size()-1;
			while (!IsOnSamePath( ResetIdx, StartIdx )) {StartIdx--;}
			
		}
		else
		{
			if ( itsStartTimes[StartIdx] > maturityTime )
			{
			}
			else 
			{
				while (itsStartTimes[StartIdx] <= maturityTime) { StartIdx++;};
				StartIdx -= 1;
			}
		}

		if ( IsOnSamePath( ResetIdx, StartIdx ))
		{
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, StartIdx, modelNb, states);

			for (iter = result->begin(), iter1 = zc->begin(); iter1 != zc->end() ; ++iter, ++iter1)
				(*iter) *= (*iter1);
		}
		else
		{
			ARM_VectorPtr fwd = ForwardDiscountFactorFromIdx(curveName, evalTime, StartIdx, itsResetTimes.size(), modelNb, states);
			ARM_VectorPtr zc  = ForwardDiscountFactorFromIdx(curveName, evalTime, ResetIdx, itsResetTimes.size(), modelNb, states);

			for (iter = result->begin(), iter1 = fwd->begin(), iter2 = zc->begin() ; iter1 != fwd->end() ; ++iter, ++iter1, ++iter2)
				(*iter) *= (*iter2)/(*iter1);
		}

		if (maturityTime != itsStartTimes[StartIdx])
		{
			double theta	= (maturityTime - itsStartTimes[StartIdx])/360.;
			double mu		= GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)
							/ GetZeroCurve()->DiscountPrice(itsStartTimes[StartIdx]/K_YEAR_LEN)
							* ( 1 + theta * itsFwdRate[StartIdx]);
			
			int k;

			for (iter = result->begin(), k = 0 ; iter != result->end() ; ++iter, k++)
				(*iter) *= mu / (1. + theta * (states->GetModelState(k, StartIdx + modelNb) + itsFromShiftedToRate[StartIdx]));
		}
		
		return static_cast<ARM_VectorPtr>(result);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime not in schedule!!");		
	}
}

ARM_VectorPtr ARM_BGMSV2F::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
	if (fabs(fwdEndTime-payTime)>7)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::Libor : case payTime != fwdEndTime not handled." );

	int modelNb = GetModelNb();

	if(ExistsInVector(itsResetTimes, fwdResetTime))
	{
		int idx = IdxFromValue(itsResetTimes, fwdResetTime);

		if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			std::vector<double>& result = new std::vector<double>( states->size() );

			for(int i = 0; i < states->size(); i++) 
				(*result)[i] = states->GetModelState(i, idx + modelNb) + itsFromShiftedToRate[idx];

			return static_cast<ARM_VectorPtr>(result);
		}
	}
	else
	{
		int infidx = IndexOfLastLowerEqInVector_DefaultFirst(fwdResetTime, itsResetTimes);
		int supidx = IndexOfFirstHigherEqInVector_DefaultLast(fwdResetTime, itsResetTimes);

		if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[infidx] - itsEndTimes[infidx]) < 7
		&& fabs(fwdEndTime - fwdStartTime + itsStartTimes[supidx] - itsEndTimes[supidx]) < 7)
		{
			double fwdrate = (GetZeroCurve()->DiscountPrice(fwdStartTime/K_YEAR_LEN) / GetZeroCurve()->DiscountPrice(fwdEndTime/K_YEAR_LEN) - 1.) / period;

			double shiftinf = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(infidx);
			double shiftsup = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(supidx);
			
			double midshift = shiftinf + (shiftsup - shiftinf) * (fwdStartTime - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

			double decal	= (midshift - 1.) * fwdrate / midshift;
			double weight	= (fwdrate/midshift - itsFwdRate[infidx]/shiftinf) / (itsFwdRate[supidx]/shiftsup - itsFwdRate[infidx]/shiftinf);

			std::vector<double>& result = new std::vector<double>( states->size() );

			for(int i = 0; i < states->size(); i++)
			{
				
				(*result)[i] = states->GetModelState(i, infidx + modelNb) 
					+ weight * (states->GetModelState(i, supidx + modelNb) - states->GetModelState(i, infidx + modelNb))
					+ decal;
			}

			return static_cast<ARM_VectorPtr>(result);
		}
	}

	return DefaultLibor(curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
}

ARM_VectorPtr ARM_BGMSV2F::SwapRate(
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const std::vector<double>& fixPayTimes,
        const std::vector<double>& fixPayPeriods,
		const std::vector<double>& fwdStartTimes,
        const std::vector<double>& fwdEndTimes,
        const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& floatPayTimes,
        const std::vector<double>& floatPayPeriods,
        const std::vector<double>& margin,
        bool isDbleNotional,
        const ARM_PricingStatesPtr& states) const
{
	int modelNb = GetModelNb();

	if(fabs(floatStartTime - fwdStartTimes[0]) < 7
	&& fabs(floatEndTime - fwdEndTimes[fwdEndTimes.size()-1]) < 7)
	{
		int idx = IdxFromValueWithTol(itsStartTimes, fwdStartTimes[0],7, true);

		if(idx > -1)
		{
			if(fabs(itsEndTimes[idx] - itsStartTimes[idx] + fwdStartTimes[0] - fwdEndTimes[fwdEndTimes.size()-1]) < 7)
			{
				std::vector<double>& result = new std::vector<double>( states->size() );
				
				bool evalMargin = false;
				for(int i = 0; i < margin.size(); i++)
				{
					if(fabs(margin[i]) > K_DOUBLE_TOL) 
					{
						evalMargin = true;
						break;
					}
				}

				if(evalMargin)
				{
					ARM_GP_VectorPtr num = AnnuityWithNominal(curveName, evalTime, floatPayTimes, floatPayPeriods, margin, states);
					ARM_GP_VectorPtr denom = Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);

					for(i = 0; i < states->size(); i++)
						(*result)[i] = states->GetModelState(i, idx + modelNb) + itsFromShiftedToRate[idx] + (*num)[i]/(*denom)[i];
				}
				else
				{
					for(i = 0; i < states->size(); i++) 
						(*result)[i] = states->GetModelState(i, idx + modelNb) + itsFromShiftedToRate[idx];
				}

				return static_cast<ARM_VectorPtr>(result);
			}
		}
		else
		{
			int infidx = IndexOfLastLowerEqInVector_DefaultFirst(fwdStartTimes[0], itsStartTimes);
			int supidx = IndexOfFirstHigherEqInVector_DefaultLast(fwdStartTimes[0], itsStartTimes);

			double length = fwdEndTimes[fwdEndTimes.size()-1] - fwdStartTimes[0];

			if(fabs(itsEndTimes[infidx] - itsStartTimes[infidx] - length) < 7
			&& fabs(itsEndTimes[supidx] - itsStartTimes[supidx] - length) < 7)
			{
				double flleg = GetZeroCurve()->DiscountPrice(fwdStartTimes[0]/K_YEAR_LEN) - GetZeroCurve()->DiscountPrice(fwdEndTimes[fwdEndTimes.size()-1]/K_YEAR_LEN);

				double annuity = 0.;
				for(int i = 0; i < fixPayPeriods.size(); i++)
					annuity += fixPayPeriods[i] * GetZeroCurve()->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN);

				double fwdrate = flleg / annuity;
								
				// petit check bidouille
				if((fwdrate < itsFwdRate[infidx] && fwdrate < itsFwdRate[supidx]) || (fwdrate > itsFwdRate[infidx] && fwdrate > itsFwdRate[supidx]))
				{
					fwdrate = itsFwdRate[infidx] + (floatStartTime - itsStartTimes[infidx]) * (itsFwdRate[supidx] - itsFwdRate[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);
				}

				bool evalMargin = false;
				for(i = 0; i < margin.size(); i++)
				{
					if(fabs(margin[i]) > K_DOUBLE_TOL) 
					{
						evalMargin = true;
						break;
					}
				}

				double shiftinf = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(infidx);
				double shiftsup = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(supidx);
				
				double midshift = shiftinf + (shiftsup - shiftinf) * (fwdStartTimes[0] - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

				double decal	= (midshift - 1.) * fwdrate / midshift;
				double weight	= (fwdrate/midshift - itsFwdRate[infidx]/shiftinf) / (itsFwdRate[supidx]/shiftsup - itsFwdRate[infidx]/shiftinf);

				std::vector<double>& result = new std::vector<double>( states->size() );

				if(evalMargin)
				{
					ARM_GP_VectorPtr num = AnnuityWithNominal(curveName, evalTime, floatPayTimes, floatPayPeriods, margin, states);
					ARM_GP_VectorPtr denom = Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);

					for(i = 0; i < states->size(); i++)
					{
						
						(*result)[i] = states->GetModelState(i, infidx + modelNb) 
							+ weight * (states->GetModelState(i, supidx + modelNb) - states->GetModelState(i, infidx + modelNb))
							+ decal + (*num)[i]/(*denom)[i];
					}
				}
				else
				{
					for(i = 0; i < states->size(); i++)
					{
						
						(*result)[i] = states->GetModelState(i, infidx + modelNb) 
							+ weight * (states->GetModelState(i, supidx + modelNb) - states->GetModelState(i, infidx + modelNb))
							+ decal;
					}
				}

				return static_cast<ARM_VectorPtr>(result);
			}
		}
	}

	return ARM_PricingModelIR::SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, 
		fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, isDbleNotional, states);
}

ARM_VectorPtr ARM_BGMSV2F::VanillaCaplet(
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
	if ( fabs(evalTime-fwdResetTime) < K_DOUBLE_TOL )
    {
		ARM_VectorPtr libor = Libor(curveName, evalTime, fwdStartTime, fwdEndTime, fwdPeriod, fwdResetTime, payTime, states);

        size_t stateSize    = states->size();

        ARM_VectorPtr dfPay = DiscountFactor(curveName, evalTime, payTime, states);

        ARM_GP_VectorPtr result( new std::vector<double>( stateSize ) );

        double payoff;

        for (size_t i(0); i<stateSize; i++)
        {
			payoff = capFloor * ( (*libor)[i] - strikesPerState[i] );

            (*result)[i] = (payoff > 0) ? period * payNotional * (*dfPay)[i] * payoff : 0.0;

		}

		return result;

	}
	else
	{
		if(ExistsInVector(itsResetTimes, fwdResetTime))
		{
			int idx = IdxFromValueWithTol(itsResetTimes, fwdResetTime, 7);

			if(fabs(payTime - fwdEndTime) < 7 && evalTime < K_DOUBLE_TOL)
			{
				if(itsCalibSecDensities[idx]->getInterestTerms().size() > 1)
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaCaplet called but it is a swaption" );
				}

				if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
				{
					double caplet = 0.;

					ARM_Heston2BOptionPricer pricer(fwdResetTime / 365., itsFwdRate[idx], 
						strikesPerState[0], capFloor,
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV01(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa1(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta1(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(idx),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol1(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV02(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa2(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta2(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(idx),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol2(),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(idx),
						((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idx)
						);

					caplet = pricer.price();

					ARM_VectorPtr result(new std::vector<double>(1, caplet * itsCalibSecDensities[idx]->ComputeLevel() * payNotional));

					return result;
				}
			}
			else if(evalTime > K_DOUBLE_TOL)
			{
				if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
				{
					ARM_VectorPtr result(new std::vector<double>(states->size()));
					ARM_VectorPtr fwdrate = Libor(curveName, evalTime, fwdStartTime, fwdEndTime, fwdPeriod, 
						fwdResetTime, fwdEndTime, states);

					int i, size = (int)states->size();
					int vidx = itsNbEffectiveReset + (int)GetModelNb();

					double opt, K;

					for(i = 0; i < size; i++)
					{
						K = strikesPerState.size() == 1 ? strikesPerState[0] : strikesPerState[i];

						ARM_Heston2BOptionPricer pricer((fwdResetTime - evalTime) / 365., (*fwdrate)[i], 
							K, capFloor,
							states->GetModelState(i, vidx),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa1(),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta1(),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(idx),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol1(),
							states->GetModelState(i, vidx+1),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa2(),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta2(),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(idx),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol2(),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(idx),
							((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idx)
							);

						opt = pricer.price();

						(*result)[i] = opt < 0. ? capFloor * ((*fwdrate)[i] - K) > 0. ? capFloor * ((*fwdrate)[i] - K) : 0. : opt;

					}

					return result;
				}
			}
		}

		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaCaplet implemented only for calibrated Libor" );
	}
}

ARM_VectorPtr ARM_BGMSV2F::VanillaSwaption(
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
	if( !(isConstantSpread && isConstantStrike) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSwaption : notional, spread, strike must be const" );

	if(itsSpreadCalib)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSwaption not implemented while calibrating on spread" );
	}

	if(evalTime > K_DOUBLE_TOL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSwaption implemented only for spot date" );
	}

	int idx;
	if((idx = IdxFromValueWithTol(itsResetTimes, swapResetTime, 7, true)) != -1)
	{
		double swapNotional = fixNotional[0];

		if(fabs(floatEndTime - floatStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			if(!isConstantNotional)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSwaption : notional must be const" );
			}

			double swopt = 0.;
			
			ARM_Heston2BOptionPricer pricer(swapResetTime / 365., itsFwdRate[idx], 
				strikesPerState(0,0), callPut,
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV01(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa1(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta1(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(idx),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol1(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV02(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa2(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta2(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(idx),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol2(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(idx),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idx)
				);

			swopt = pricer.price();

			ARM_VectorPtr result(new std::vector<double>(1, swopt * itsCalibSecDensities[idx]->ComputeLevel() * swapNotional));

			return result;
		}
		else if(itsCalibSecDensities[idx]->getInterestTerms().size() == 1)
		{
			// approx à partir des libors : vol swap vol fra
			std::vector<double> weight;
			ARM_IntVector idxs(floatResetTimes.size());
			
			int idxend = IdxFromValueWithTol(itsEndTimes, floatEndTimes[floatEndTimes.size()-1], 7, true);

			if(idxend == -1)
			{
				ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV2F::VanillaSwaption swaption schedule not in calibration schedule");
			}

			for(int i = 0; i < idxs.size(); i++)
			{
				idxs[i] = IdxFromValueWithTol(itsResetTimes, floatResetTimes[i], 7, true);

				if(idxs[i] == -1)
				{
					ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_BGMSV2F : libor in swaption not in calibration schedule");
				}
			}
			

			double fwdrate = (*SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, floatStartTimes, 
								floatEndTimes, floatIntTerms, floatEndTimes, floatIntTerms, std::vector<double>(1,0.),isConstantNotional, states))[0];

			double annuity = (*AnnuityWithNominal(curveName, evalTime, fixPayTimes, fixPayPeriods, fixNotional, states))[0];

			GetWeightSwaptionApprox(curveName, evalTime, floatNotional, floatStartTimes, floatEndTimes, floatIntTerms, idx, fwdrate, annuity, weight);

			// calcul du shift équivalent et du rho équivalent
			double swshift, swrho1, swrho2, swlevel;

			swshift = GetEquivSwaptionShift(idx, idxs, weight);
			swlevel	= GetEquivSwaptionLevel(idx, idxs, weight);
			swrho1	= GetEquivSwaptionRho(idx, idxs, weight, swlevel, 1);
			swrho2	= GetEquivSwaptionRho(idx, idxs, weight, swlevel, 2);

			double swopt = 0.;
			
			ARM_Heston2BOptionPricer pricer(swapResetTime / 365., itsFwdRate[idx], 
				strikesPerState(0,0), callPut,
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV01(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa1(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta1(),
				swrho1,
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol1(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetV02(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetKappa2(),
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetTheta2(),
				swrho2,
				((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetVVol2(),
				swshift,
				swlevel
				);

			swopt = pricer.price();

			ARM_VectorPtr result(new std::vector<double>(1, swopt * annuity));

			return result;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSwaption implemented only for calibrated CMS" );
}

ARM_VectorPtr ARM_BGMSV2F::VanillaSpreadOptionLet(const string& curveName,
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
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV2F::VanillaSpreadOptionLet not implemented" );
}


void ARM_BGMSV2F::GetWeightSwaptionApprox(const string& curveName, double evalTime, const std::vector<double>& floatNotional, 
										  const std::vector<double>& floatStartTimes, const std::vector<double>& floatEndTimes, const std::vector<double>& floatPayPeriods,
										  int idx, double fwdrate, double annuity, std::vector<double>& weight) const
{
	weight.resize(floatStartTimes.size());

	for(int k = 0; k < floatStartTimes.size(); k++)
	{
		weight[k] = floatNotional[k] * ((*DiscountFactor(curveName, evalTime, floatStartTimes[k], ARM_PricingStatesPtr(NULL)))[0] - (*DiscountFactor(curveName, evalTime, floatEndTimes[k], ARM_PricingStatesPtr(NULL)))[0]) / (annuity * fwdrate);
	}
}

double ARM_BGMSV2F::GetEquivSwaptionShift(int idx, const ARM_IntVector& idxs, std::vector<double>& weight) const
{
	double shift = 0.;
	for(int k = 0; k < weight.size(); k++)
	{
		double shk =  idxs.size() == 0 ? ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(idx+k) : ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetShift(idxs[k]);

		shift += weight[k] * (shk - 1.) / shk;

		// modif des poids pour la suite
		weight[k] /= shk;
	}

	shift = 1. / (1. - shift);

	weight *= shift;

	return shift;
}

double ARM_BGMSV2F::GetEquivSwaptionLevel(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, const std::vector<double> * optweight) const
{
	int i, j, k;
	double dt, corr, var = 0.;
	double voli, volj;

	const std::vector<double> * tmpweight = optweight == NULL ? &weight : optweight;

	for(k = 0; k <= idx; k++)
	{
		dt	= k == 0 ? itsResetTimes[k] / 365. : (itsResetTimes[k] - itsResetTimes[k-1]) / 365.;

		for(i = 0; i < weight.size(); i++)
		{
			if(idxs.size() == 0)
				voli = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(i + idx);
			else
				voli = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idxs[i]);

			for(j = 0; j < tmpweight->size(); j++)
			{
				if(idxs.size() == 0)
					volj = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(j + idx);
				else
					volj = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idxs[j]);

				if(idxs.size() == 0)
					corr = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->RateRateCorrel(itsResetTimes[k], itsResetTimes[i + idx], i + idx, itsResetTimes[j + idx], j + idx);
				else
					corr = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->RateRateCorrel(itsResetTimes[k], itsResetTimes[idxs[i]], idxs[i], itsResetTimes[idxs[j]], idxs[j]);

				var += weight[i] * (*tmpweight)[j] * voli * volj * corr * dt;
			}
		}
	}

	return sqrt(var * 365. / itsResetTimes[idx]);
}

double ARM_BGMSV2F::GetEquivSwaptionRho(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, double swlevel,
										int UnOuDeux) const
{
	int k;
	double rhok, volk, rho = 0.;

	for(k = 0; k < weight.size(); k++)
	{
		if(idxs.size() == 0)
		{
			rhok = UnOuDeux == 1 ? ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(k + idx) : ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(k + idx);
			volk = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(k + idx);
		}
		else
		{
			rhok = UnOuDeux == 1 ? ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho1(idxs[k]) : ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetRho2(idxs[k]);
			volk = ((ARM_ModelParamsBGMSV2F*) GetModelParams())->GetSigma(idxs[k]);
		}

		rho += volk * rhok * weight[k];
	}

	rho /= swlevel;

	return rho;
}

CC_END_NAMESPACE()