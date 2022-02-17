/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// gpmodels
#include "gpmodels/bgmsv1f.h"
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


CC_BEGIN_NAMESPACE( ARM )



////////////////////////////////////////////////////
///	Class  : ARM_BGMSV1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_BGMSV1F::ARM_BGMSV1F( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, bool AllowInterpol, bool ComputeDrift, bool Proxy)
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
	itsSimpleEulerScheme = false;
	//itsAllowInterpol = true;
	itsInterpolZC = true;
	itsPrevEvalTime = 0.;
}


////////////////////////////////////////////////////
///	Class  : ARM_BGMSV1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_BGMSV1F::ARM_BGMSV1F(const ARM_BGMSV1F& rhs)
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
	itsSimpleEulerScheme(rhs.itsSimpleEulerScheme),
	itsCalibWeight(rhs.itsCalibWeight),
	itsInterpolZC(rhs.itsInterpolZC),
	itsRandGen(rhs.itsRandGen),
	itsPrevEvalTime(rhs.itsPrevEvalTime)
{
	DuplicateCloneablePointorVectorInPlace<std::vector<double>>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelLocalRealVar, itsModelLocalRealVar);
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsCalibSecDensities, itsCalibSecDensities);
}	


ARM_BGMSV1F& ARM_BGMSV1F::operator = (const ARM_BGMSV1F& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_BGMSV1F();

		new (this) ARM_BGMSV1F (rhs);
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_BGMSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_BGMSV1F::~ARM_BGMSV1F()
{
	DeletePointorVector<std::vector<double>>( itsEigenValues );
	DeletePointorVector<ARM_GP_Matrix>( itsModelLocalRealVar );
}

bool ARM_BGMSV1F::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(params.DoesModelParamExist(ARM_ModelParamType::Shift) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: shift structure is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::Alpha) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: alpha structure (level) is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::QParameter) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: Q structure (rho Rate/variance) is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: VolOfVol structure is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::VolMeanReversion) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: VolMeanReversion structure is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::InitialVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: InitialVol structure is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::LongTermVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: LongTermVol structure is needed for BGMSV1F" );

	if(params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateModelParams: beta correlation structure is needed for BGMSV1F" );


	return true;
}

void ARM_BGMSV1F::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	/*
	if(	numerairePtr->GetType() != ARM_Numeraire::TerminalZc && numerairePtr->GetType() != ARM_Numeraire::TerminalEventZc )
    {		
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::SetNumeraire: only numeraire supported is Terminal ZC" );
    }
	*/
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

void ARM_BGMSV1F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);

//	if (!itsCalibratedStatus && !isNumerical)
//		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ValidateCalibMethod: should calib functionals first!!!" );

	calibMethod.DefaultValidateWithModel(*this);
}

void ARM_BGMSV1F::PreProcessing(ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

void ARM_BGMSV1F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	// TO DO 
}


void ARM_BGMSV1F::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	// TO DO 
}


ARM_PricingStatesPtr ARM_BGMSV1F::FirstPricingStates(size_t bucketSize) const
{
	const size_t nbPayoffs			= 0;
	size_t nbModelStates			= itsNbEffectiveReset + 1;
	size_t factorsNb				= FactorCount();

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	for(int i = 0; i < nbModelStates - 1; i++)
	{
		// initialisation des taux shiftés
		double initRateValue	= itsFwdRate[i] / ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(i);

		for(int k = 0; k < bucketSize; k++) 
		{
			initStates->SetModelState(k, i, initRateValue);
		}
	}

	double initVarValue		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();

	for(int k = 0; k < bucketSize; k++)
	{
		initStates->SetModelState(k, itsNbEffectiveReset, initVarValue);
	}

	return initStates;
}

ARM_PricingStatesPtr ARM_BGMSV1F::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
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

ARM_PricingStatesPtr ARM_BGMSV1F::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	if(itsCalibrationStatus)
	{	
		if(itsCalibratedStatus == false)
		{
			int k, size = itsResetTimes.size();		

			/*
			std::vector<double> calibTimes(size);
			calibTimes[0] = 0.;
			for(k = 1; k < size; k++) calibTimes[k] = itsResetTimes[k-1];

			ARM_VectorVector strikes, mktvols;
			std::vector<double> resettimes, forwards, atmvols;
			std::vector<double> shifts, levels, weights;

			getCalibrationData(resettimes, forwards, strikes, mktvols, atmvols, weights);

			{
				resettimes /= 365.;

				double v0 = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0();
				double kappa = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
				double rho = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(0);
				double shift = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(0);
				double theta = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
				double nu = 0.;

				bool calibkappa = false;
				bool calibtheta = theta < K_DOUBLE_TOL ? true : false;
				bool calibrho = fabs(rho) > 0.999 ? true : false;
				bool calibshift = fabs(shift) < K_DOUBLE_TOL ? true : false;
				bool calibnu = true;
				bool localrho = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLocalRhoCalib();

				ARM_SmileCalibration_Params_Heston params(v0, kappa, theta, rho, nu, shift, 1., 
					calibkappa, calibtheta, calibrho, calibnu, calibshift, localrho);

				ARM_SmileCalibration_Heston func;
				
				func.Init(resettimes, forwards, mktvols, strikes, true, forwards, atmvols, weights, &params);
				func.Calibrate();

				if(calibshift == false) shifts.resize(levels.size(), ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(0));

				finishCalibration(&calibTimes, params);
			}
			*/

			((ARM_ModelParamsBGMSV1F*) GetModelParams())->Calibrate(itsResetTimes, itsFwdRate, itsCalibSecDensities);
			
			itsFromShiftedToRate.resize(size);

			// pour passer des libors shiftés aux libors : fwd x (shift - 1) / shift
			for(k = 0; k < size; k++) 
			{
				itsFromShiftedToRate[k] = (((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(k) - 1.) * itsFwdRate[k] / ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(k);
			}

			itsCalibratedStatus = true;
		}

		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}

	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );

	return newStates;
}

ARM_VectorPtr ARM_BGMSV1F::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,getTerminalTime()));
}

std::vector<double>& ARM_BGMSV1F::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetTimes.Clone());
}


void ARM_BGMSV1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
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


string ARM_BGMSV1F::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "\n\n";
    os << indent << "SVBGM1F : Shifted Stochastic Volatility BGM with One process for variance (Heston time dependant) \n";
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

void ARM_BGMSV1F::CalcNbEffectiveReset(const std::vector<double>& timeSteps)
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

	((ARM_ModelParamsBGMSV1F*) GetModelParams())->checkFactorCount(itsNbEffectiveReset);
}

void ARM_BGMSV1F::ModelStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps)
{
	CalcNbEffectiveReset(timeSteps);

	ARM_MatrixVector auxLocalCov;
	ModelStateLocalVariances( timeSteps, auxLocalCov);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);

	ModelStateLocalStdDev(timeSteps, auxLocalCov);
}

void ARM_BGMSV1F::NumMethodStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
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
		}
	}
}

void ARM_BGMSV1F::NumMethodStateGlobalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (int i=0;i<timeStepsSize;i++)
	{
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(totalFwds+1,1.0);
	}
}


void ARM_BGMSV1F::ModelStateLocalStdDev( const std::vector<double>& timeSteps, const ARM_MatrixVector& auxLocalCov)
{
	ARM_MatrixVector modStateLocalVars;
	
	int totalFwds		= itsNbEffectiveReset;
	int factorsNb		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->EffFactorCount();
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

	((ARM_ModelParamsBGMSV1F*) GetModelParams())->SetFactorCount(maxfactorsNb);

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

void ARM_BGMSV1F::NumMethodStateLocalVariancesAndStdDev( const std::vector<double>& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_BGMSV1F::ModelStateLocalVariances(const std::vector<double>& timeSteps, ARM_MatrixVector& localVariances ) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;


#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ModelStateLocalVariances: first time step != 0" );
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

					//fac		= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
					fac		= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
					
					if(fabs(fac) < K_DOUBLE_TOL)
					{
						rhok = rhok > 0.9999 - K_DOUBLE_TOL ? rhok - 0.0001 : rhok + 0.0001;
						//fac	= (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));
						fac	= sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok));
					}

					// (*localVariances[i])(j,k) = (*localVariances[i])(k,j) = rhojk / (rhoj * rhok + sqrt((1. - rhoj*rhoj) * (1. - rhok*rhok)));

					(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = (rhojk - rhoj*rhok) / fac;
				}
			}
		}

		fromTime = toTime;
	}
}

void ARM_BGMSV1F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
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
	
	double rcov, vcov;
	double rdrift = 0., vdrift = 0., adj = 0.;
	double var, libi;
	double vvol = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(nextTime);
	double kappa = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
	double theta = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
	double wesp2 = 1. + 0.5 * vvol * vvol / (kappa * theta);
	double w1E2_var = 0.5 * theta * vvol * vvol * SQR(1. - exp(-kappa*dt));
	double w2E2_var = vvol * vvol * exp(-kappa*dt)*(1.-exp(-kappa*dt))/kappa;
	double vesp, vesp2, vvar, adjt, scale;
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

				libi	= states->GetModelState(n, i);

				if(itsComputeDriftAdj)
				{
					rdrift = 0.;

					for(j = i + 1; j <= iLast; j++)
					{
						if(IsOnSamePath(i,j) == false) continue;

						rcov	= (*itsModelLocalRealVar[timeIndex])(i - iFirst, j - iFirst) * dt;				// corrélation taux i, taux j

						adj		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(j) * states->GetModelState(n, j) * itsDelta[j]
								/ (1. + itsDelta[j] * (states->GetModelState(n, j) + itsFromShiftedToRate[j]));

						rdrift	-= adj * rcov;
					}

					rdrift *= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * var;
				}

				scale	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * racdt;
				libi	*= exp(rdrift - 0.5 * scale * scale * var + scale * sqrt(var) * x[i - iFirst]);

				states->SetModelState(n, i, libi);
			}

			// simulation de la variance
			vdrift = 1.;

			if(itsComputeDriftAdj)
			{
				for(i = iFirst; i <= iLast; i++)
				{
					if(IsOnSamePath(i,iFirst) == false) continue;

					libi	= states->GetModelState(n, i);
					vcov	= (*itsModelLocalRealVar[timeIndex])(i - iFirst, aliveFwds);
					
					adj		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * states->GetModelState(n, i) * itsDelta[i]
							/ (1. + itsDelta[i] * (states->GetModelState(n, i) + itsFromShiftedToRate[i]));

					vdrift	+= adj * vcov * vvol / kappa;
				}
			}
		

			adjt	= exp(- kappa * vdrift * dt);
			vesp	= var * adjt + (theta / vdrift) * (1. - adjt);
			vesp2	= wesp2 * vesp * vesp + (1. - wesp2) * adjt * adjt * var * var;
			vvar	= log(vesp2 / (vesp * vesp));
			var		= vvar > 0. ? vesp * exp(-0.5 * vvar + sqrt(vvar) * x[aliveFwds]) : vesp;

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

		MCFromToNextTime(states, timeIndex, x);


	}
}

void ARM_BGMSV1F::MCFromToNextTime(ARM_PricingStatesPtr& states, int timeIndex, const ARM_GP_Matrix& x) const
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

	double rcov, vcov;
	double rdrift = 0., vdrift = 0., adj = 0.;
	double var, libi, rhoi;
	double vvol = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVol(nextTime);
	double kappa = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa();
	double theta = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta();
	double wesp2 = 1. + 0.5 * vvol * vvol / (kappa * theta);
	double w1E2_var = 0.5 * theta * vvol * vvol * SQR(1. - exp(-kappa*dt));
	double w2E2_var = vvol * vvol * exp(-kappa*dt)*(1.-exp(-kappa*dt))/kappa;
	double E_var, E2_var, psi, adjt, scale, tipmo, a_var, b_var, invbeta, p_var, u_var;
	double akappa, atheta;
	double varf;
	double K0, K1, K2, K3, K4, w1 = 0.5, w2 = 0.5;

	for(n = 0; n < statesNb; n++)
	{
		var		= states->GetModelState(n, fwdsNb + modelNb);

		// simulation de la variance
		vdrift = 1.;

		if(itsComputeDriftAdj)
		{
			for(i = iFirst; i <= iLast; i++)
			{
				if(IsOnSamePath(i,iFirst) == false) continue;

				libi	= states->GetModelState(n, i + modelNb);
				vcov	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(i);
				
				adj		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * states->GetModelState(n, i + modelNb) * itsDelta[i]
						/ (1. + itsDelta[i] * (states->GetModelState(n, i) + itsFromShiftedToRate[i]));

				vdrift	+= adj * vcov * vvol / kappa;
			}
		}
		
		akappa		= kappa * vdrift;
		atheta		= theta / vdrift;
		adjt		= exp(- akappa * dt);
		E_var		= var * adjt + atheta * (1. - adjt);
		w1E2_var	= 0.5 * theta * vvol * vvol * SQR(1. - adjt);
		w2E2_var	= vvol * vvol * adjt * (1.- adjt) / akappa;
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
			rdrift	= 0.;
			
			libi	= states->GetModelState(n, i + modelNb);

			if(itsComputeDriftAdj)
			{
				for(j = i + 1; j <= iLast; j++)
				{
					if(IsOnSamePath(i,j) == false) continue;

					rcov = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->RateRateCorrel(nextTime, itsResetTimes[i],i,itsResetTimes[j],j) * dt;

					adj		= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(j) * states->GetModelState(n, j + modelNb) * itsDelta[j]
							/ (1. + itsDelta[j] * (states->GetModelState(n, j + modelNb) + itsFromShiftedToRate[j]));

					rdrift	-= adj * rcov;
				}

				rdrift *= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i) * (w1*var + w2*varf);
			}
			
			scale	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i);

			rhoi	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(i);

			K0	= - rhoi * kappa * theta * dt * scale / vvol;
			K1	= w1 * dt * (kappa * rhoi * scale / vvol - 0.5 * scale * scale) - rhoi * scale/ vvol;
			K2	= w1 * dt * (kappa * rhoi * scale / vvol - 0.5 * scale * scale) + rhoi * scale/ vvol;
			K3	= w1 * dt * (1. - rhoi*rhoi) * scale * scale;
			K4	= w2 * dt * (1. - rhoi*rhoi) * scale * scale;

			libi	*= exp(rdrift + K0 + K1 * var + K2 * varf + sqrt(K3 * var + K4 * varf) * x(n,i - iFirst));

			states->SetModelState(n, i + modelNb, libi);
		}
	}
}

///////////////////////
//
//		les prix
//
///////////////////////

void ARM_BGMSV1F::setNumericalModelFitter(ARM_NumericalModelFitter * numericalModelFitter)
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
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F : nul weight for each density" );
		}

		std::vector<double> resetDates(size); 
		std::vector<double> startDates(size); 
		std::vector<double> endDates(size); 
		std::vector<double> delta(size); 
		itsFwdRate.resize(size);
		itsCalibWeight.resize(size,1);
		itsATMVols.resize(size,0.);

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

		((ARM_ModelParamsBGMSV1F*) GetModelParams())->checkFactorCount(size);

		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);

		itsCalibratedStatus = false;
	}
}

void ARM_BGMSV1F::computeWeights() 
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

ARM_VectorPtr ARM_BGMSV1F::DiscountFactor(const string& curveName, double evalTime, double maturityTime, 
										const ARM_PricingStatesPtr& states) const
{
	// if(itsAllowInterpol)
	if(itsInterpolZC)
		return ForwardDiscountFactor(curveName,evalTime,evalTime,maturityTime,states);
	else
		return DiscountFactorNoInterpol(curveName,evalTime,maturityTime,states);
}

ARM_VectorPtr ARM_BGMSV1F::ForwardDiscountFactor( const string& curveName, double evalTime, double startTime, double endTime, const ARM_PricingStatesPtr& states) const
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

ARM_VectorPtr ARM_BGMSV1F::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t IdxFrom, size_t IdxTo , size_t modelNb, const ARM_PricingStatesPtr& states) const
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

ARM_VectorPtr ARM_BGMSV1F::DiscountFactorNoInterpol(const string& curveName, double evalTime, double maturityTime, 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::DiscountFactor: MaturityTime after schedule's last date");


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

ARM_VectorPtr ARM_BGMSV1F::Libor( 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::Libor : case payTime != fwdEndTime not handled." );

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

			double shiftinf = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(infidx);
			double shiftsup = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(supidx);
			
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

ARM_VectorPtr ARM_BGMSV1F::SwapRate(
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
		int idx = IdxFromValue(itsStartTimes, fwdStartTimes[0],7);

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

				double shiftinf = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(infidx);
				double shiftsup = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(supidx);
				
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

ARM_VectorPtr ARM_BGMSV1F::MaxRate(		
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
        double firstReset,
		double firstStart,
		const ARM_VectorPtr& firstRate,
		int MaxOrMin,
		int ResetFreq,
		const ARM_VectorPtr& strikes,
		int CapOrFloor,
		double RhoMinMax,
		bool IsAccrued,
		double MinAccrued,
		double MaxAccrued,
        const ARM_PricingStatesPtr& states) const
{
	int modelNb = GetModelNb();

	if(fabs(floatStartTime - fwdStartTimes[0]) < 7
	&& fabs(floatEndTime - fwdEndTimes[fwdEndTimes.size()-1]) < 7
	&& firstRate->size() == states->size())
	{
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
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::MaxRate not implemented for non null margin" );
		}

		ARM_VectorPtr rate = SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, 
								fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, true, states);

		if(rate->size() == firstRate->size())
		{
			double dt = (evalTime - firstReset) / 365.;

			double firstshift, firstdecal, firstsigma, shift, decal, sigma, firstfwd, fwd;
			int infidx, supidx;

			// le décalage précédent
			infidx = supidx	= IdxFromValue(itsStartTimes, firstStart,7);
			
			if(infidx == -1)
			{
				infidx = IndexOfLastLowerEqInVector_DefaultFirst(firstStart, itsStartTimes);
				supidx = IndexOfFirstHigherEqInVector_DefaultLast(firstStart, itsStartTimes);
			}
			
			if(supidx > infidx)
			{
				double weight = (firstStart - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

				firstdecal = itsFromShiftedToRate[infidx] + (itsFromShiftedToRate[supidx] - itsFromShiftedToRate[infidx]) * weight;
				firstsigma = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx) + 
					weight * (((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(supidx) - ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx));

			}
			else
			{
				firstdecal	= itsFromShiftedToRate[infidx];
				firstsigma	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx);
				firstfwd	= itsFwdRate[infidx];
				firstshift	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(infidx);
			}

			// le décalage courant
			infidx = supidx	= IdxFromValue(itsStartTimes, fwdStartTimes[0],7);

			if(infidx == -1)
			{
				infidx = IndexOfLastLowerEqInVector_DefaultFirst(fwdStartTimes[0], itsStartTimes);
				supidx = IndexOfFirstHigherEqInVector_DefaultLast(fwdStartTimes[0], itsStartTimes);
			}

			if(supidx > infidx)
			{
				double weight = (fwdStartTimes[0] - itsStartTimes[infidx]) / (itsStartTimes[supidx] - itsStartTimes[infidx]);

				decal = itsFromShiftedToRate[infidx] + (itsFromShiftedToRate[supidx] - itsFromShiftedToRate[infidx]) * weight;
				sigma = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx) + 
					weight * (((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(supidx) - ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx));
			}
			else
			{
				decal	= itsFromShiftedToRate[infidx];
				sigma	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(infidx);
				fwd		= itsFwdRate[infidx];
				shift	= ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(infidx);
			}

			// double maxdecal	= 0.5 * (decal + firstdecal);
			// double maxsigma	= 0.5 * (sigma + firstsigma);
			double maxdecal = MaxOrMin * (fwd - firstfwd) > 0 ? decal : firstdecal;
			double maxsigma = MaxOrMin * (fwd - firstfwd) > 0 ? sigma : firstsigma;
			double maxshift	= MaxOrMin * (fwd - firstfwd) > 0 ? shift : firstshift;
				
			if(evalTime > itsPrevEvalTime)
			{
				itsU.resize(states->size());
				itsRandGen.draw(itsU);
			}

			double x0, x1, x;
			double var;
			int fwdsNb = itsNbEffectiveReset;
			double a, b, c, d, r;
			
			std::vector<double>& result = new std::vector<double>( states->size() );
			
			for(i = 0; i < states->size(); i++) 
			{
				var = states->GetModelState(i, fwdsNb + modelNb) * maxsigma * maxsigma * dt;
				x0	= ((*firstRate)[i] - firstdecal) * firstshift / maxshift;
				x1	= ((*rate)[i] - decal) * shift / maxshift;
				x	= MaxOrMin * (x1 - x0) > 0 ? x1 : x0;

				a	= 1.;
				b	= - log(x0) - log(x1);
				c	= 0.5 * log(1.- itsU[i]) * var + log(x0) * log(x1);
				d	= b * b - 4. * a * c;

				r	= d < 0. ? x : exp(0.5*(-b + MaxOrMin * sqrt(d)));
				
				(*result)[i] = maxdecal + (MaxOrMin * (x - r) > 0. ? x : r);
			}

			itsPrevEvalTime = evalTime;

			return static_cast<ARM_VectorPtr>(result);
		}
	}

	return ARM_PricingModelIR::MaxRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, 
		fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, firstReset, firstStart, firstRate,
		MaxOrMin, ResetFreq, strikes,CapOrFloor,RhoMinMax,IsAccrued,MinAccrued,MaxAccrued,states);
}

ARM_VectorPtr ARM_BGMSV1F::VanillaCaplet(
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
		if(itsSpreadCalib)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaCaplet not implemented while calibrating on spread" );
		}

		if(ExistsInVector(itsResetTimes, fwdResetTime))
		{
			int idx = IdxFromValue(itsResetTimes, fwdResetTime);

			if(fabs(payTime - fwdEndTime) < 7 && evalTime < K_DOUBLE_TOL)
			{
				if(itsCalibSecDensities[idx]->getInterestTerms().size() > 1)
				{
					ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaCaplet called but it is a swaption" );
				}

				if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
				{
					double caplet = 0.;

					ARM_HestonOptionPricerVVolt pricer(fwdResetTime / 365., itsFwdRate[idx], strikesPerState[0], capFloor,
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0(),
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa(),
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta(),
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx),
						*((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVolt(),
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idx),
						((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx)
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

					ARM_Curve * vvolt = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVolt(evalTime);

					for(i = 0; i < size; i++)
					{
						K = strikesPerState.size() == 1 ? strikesPerState[0] : strikesPerState[i];

						ARM_HestonOptionPricerVVolt pricer((fwdResetTime - evalTime) / 365., (*fwdrate)[i], K, capFloor,
							states->GetModelState(i, vidx),
							((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa(),
							((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta(),
							((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx),
							*vvolt,
							((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idx),
							((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx)
							);

						opt = pricer.price();

						(*result)[i] = opt < 0. ? capFloor * ((*fwdrate)[i] - K) > 0. ? capFloor * ((*fwdrate)[i] - K) : 0. : opt;

					}
					
					delete vvolt;

					return result;
				}
			}
		}

		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaCaplet implemented only for calibrated Libor" );
	}
}

ARM_VectorPtr ARM_BGMSV1F::VanillaSwaption(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSwaption : notional, spread, strike must be const" );

	if(itsSpreadCalib)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSwaption not implemented while calibrating on spread" );
	}

	if(evalTime > K_DOUBLE_TOL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSwaption implemented only for spot date" );
	}

	int idx;
	if((idx = IdxFromValue(itsResetTimes, swapResetTime, 7)) != -1)
	{
		double swapNotional = fixNotional[0];

		if(fabs(floatEndTime - floatStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			if(!isConstantNotional)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSwaption : notional must be const" );
			}

			double swopt = 0.;
			
			ARM_HestonOptionPricerVVolt pricer(swapResetTime / 365., itsFwdRate[idx], strikesPerState(0,0), callPut,
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx),
				*((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVolt(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idx),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx)
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
			
			int idxend = IdxFromValue(itsEndTimes, floatEndTimes[floatEndTimes.size()-1], 7);

			if(idxend == -1)
			{
				ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::VanillaSwaption swaption schedule not in calibration schedule");
			}

			for(int i = 0; i < idxs.size(); i++)
			{
				idxs[i] =IdxFromValue(itsResetTimes, floatResetTimes[i], 7);

				if(idxs[i] == -1)
				{
					ARM_THROW(ERR_INVALID_ARGUMENT, "ARM_BGMSV1F : libor in swaption not in calibration schedule");
				}
			}
			

			double fwdrate = (*SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, floatStartTimes, 
								floatEndTimes, floatIntTerms, floatEndTimes, floatIntTerms, std::vector<double>(1,0.),isConstantNotional, states))[0];

			double annuity = (*AnnuityWithNominal(curveName, evalTime, fixPayTimes, fixPayPeriods, fixNotional, states))[0];

			GetWeightSwaptionApprox(curveName, evalTime, floatNotional, floatStartTimes, floatEndTimes, floatIntTerms, idx, fwdrate, annuity, weight);

			// calcul du shift équivalent et du rho équivalent
			double swshift, swrho, swlevel;

			swshift = GetEquivSwaptionShift(idx, idxs, weight);
			swlevel	= GetEquivSwaptionLevel(idx, idxs, weight);
			swrho	= GetEquivSwaptionRho(idx, idxs, weight, swlevel);

			double swopt = 0.;
			
			ARM_HestonOptionPricerVVolt pricer(swapResetTime / 365., fwdrate, strikesPerState(0,0), callPut,
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetV0(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa(),
				((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta(),
				swrho,
				*((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVolt(),
				swshift,
				swlevel
				);

			swopt = pricer.price();

			ARM_VectorPtr result(new std::vector<double>(1, swopt * annuity));

			if(itsCalibrationStatus)
			{
				bool success;

				(*result)[0] = VanillaImpliedVol_BS(fwdrate,strikesPerState(0,0),swapResetTime/365.,(*result)[0] / 100.,callPut,NULL,&success);
			}

			return result;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSwaption implemented only for calibrated CMS" );
}

ARM_VectorPtr ARM_BGMSV1F::ImpliedVol(
		const string& curveName, 
		double evalTime,
		double payTime,
		double period,
        double payNotional,
		double fwdResetTime,	/// used for volatility computation
		double fwdStartTime,
        double fwdEndTime,
		double fwdPeriod,
        const std::vector<double>& strikesPerState,
        int capFloor,
        const ARM_PricingStatesPtr& states) const
{
	if(fabs(evalTime - fwdResetTime) < 7)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F : Implied Vol not computed for spot swaption");

	int idx;
	if((idx = IdxFromValue(itsResetTimes, fwdResetTime, 7)) != -1)
	{
		if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			ARM_VectorPtr fwdrate = Libor(curveName, evalTime, fwdStartTime, fwdEndTime, fwdPeriod, fwdResetTime,
				fwdEndTime, states);

			int k, size = (int)states->size();

			ARM_VectorPtr result(new std::vector<double>(size));

			int vidx = itsNbEffectiveReset + (int)GetModelNb();

			double K, impvol;
			bool success;

			ARM_Curve * vvolt = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetVVolt(evalTime);

			for(k = 0; k < size; k++)
			{
				double swopt = 0.;

				K = strikesPerState.size() == 1 ? strikesPerState[0] : strikesPerState[k];

				ARM_HestonOptionPricerVVolt pricer((fwdResetTime - evalTime) / 365., (*fwdrate)[k], K, capFloor,
					states->GetModelState(k, vidx),
					((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetKappa(),
					((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetTheta(),
					((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idx),
					*vvolt,
					((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idx),
					((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idx)
					);

				swopt = pricer.price();

				ARM_ImpliedVolBS inverse((*fwdrate)[k], K, (fwdResetTime - evalTime)/365., capFloor);

				impvol = inverse.vol(swopt,&success);

				(*result)[k] = success == false ? itsATMVols[idx] : impvol;

			}

			delete vvolt;

			return result;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::ImpliedVol implemented only for calibrated CMS" );
}

double ARM_BGMSV1F::UnderlyingCorrelation(string underlyingType,
										  double fromTime, double toTime,
										  double startTime1, double endTime1, 
										  double startTime2, double endTime2,
										  double startTime3, double endTime3,
										  double startTime4, double endTime4) const
{
	if ( stringGetUpper(underlyingType)!="FWD" && stringGetUpper(underlyingType)!= "CMS" )
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only FWD or CMS type supported!");
	
	std::vector<double> VarCovar(3);
	std::vector<double> swapfwd(4);
	
	ComputeUnderlyingVarCovar( fromTime,toTime,startTime1,endTime1,startTime2,endTime2,VarCovar,
		swapfwd);

	double stdev1 = VarCovar[0];
	double stdev2 = VarCovar[1];
	double covariance_1_2 = VarCovar[2]; 
	double correlation_1_2 = covariance_1_2/(stdev1*stdev2);
	return  correlation_1_2 ;	
}

void ARM_BGMSV1F::ComputeUnderlyingVarCovar(double fromTime, double toTime, 
											double startTime1, double endTime1, 
											double startTime2, double endTime2,
											std::vector<double>& varcovar, std::vector<double>& swapfwd) const
{
	double asOfDate	= GetAsOfDate().GetJulian();

	if(fromTime > K_DOUBLE_TOL)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar defined only from AsOfDate");
	}

	if(fabs(startTime1 - toTime) > 7 || fabs(startTime2 - toTime) > 7)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar defined only for spot integrated covariance");
	}

	int idx = IdxFromValue(itsResetTimes, toTime, 7);

	if(idx == -1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar defined only on calibrated schedule");
	}

	if(itsCalibSecDensities[0]->getInterestTerms().size() > 1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar defined only for libor model");
	}

	ARM_ZeroCurvePtr ZcCurve	= GetZeroCurve();
	ARM_Currency* Ccy			= ZcCurve->GetCurrencyUnit();
	string curveCcy				= Ccy->GetCcyName();
	int fixFreq					= itsCalibSecDensities[0]->getFrequency();
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

	std::vector<double>& fixPayTimes1		= fixDateStrip1.GetPaymentDates();
	std::vector<double>& fixPayPeriods1	= fixDateStrip1.GetInterestTerms();
	std::vector<double>& floatStartTimes1	= fixDateStrip1.GetFlowStartDates();
	std::vector<double>& floatPayTimes1	= fixDateStrip1.GetFlowEndDates();
	std::vector<double>& floatIntTerms1	= fixDateStrip1.GetInterestTerms();
	std::vector<double>& floatResetDate1	= fixDateStrip1.GetResetDates();

	double swapResetTime1	= (*floatResetDate1)[0]-asOfDate;	
	int size1				= floatResetDate1->size();
	int i,nbFixFlows1		= fixPayTimes1->size();    
	
	for(i=0;i<nbFixFlows1;++i)
		(*fixPayTimes1)[i] = (*fixPayTimes1)[i]-asOfDate;

	std::vector<double> fixNotional1 (nbFixFlows1,1.0);
	std::vector<double> floatNotional1 (nbFixFlows1,1.0);
	
	for(i=0;i<size1;++i)
	{
		(*floatResetDate1)[i] = (*floatResetDate1)[i]-asOfDate;
		(*floatStartTimes1)[i] = (*floatStartTimes1)[i]-asOfDate;
		(*floatPayTimes1)[i] = (*floatPayTimes1)[i]-asOfDate;
	}

	// CMS2
	ARM_Date startDate2(asOfDate+startTime2);
	ARM_Date endDate2(asOfDate+endTime2);
	ARM_DateStrip fixDateStrip2( startDate2, endDate2, fixFreq, fixDayCount, fixCalendar,
		K_MOD_FOLLOWING, K_ADJUSTED, K_SHORTSTART, GETDEFAULTVALUE, fixFreq, GETDEFAULTVALUE,
		fixCalendar );

	std::vector<double>& fixPayTimes2		= fixDateStrip2.GetPaymentDates();
	std::vector<double>& fixPayPeriods2	= fixDateStrip2.GetInterestTerms();
	std::vector<double>& floatStartTimes2	= fixDateStrip2.GetFlowStartDates();
	std::vector<double>& floatPayTimes2	= fixDateStrip2.GetFlowEndDates();
	std::vector<double>& floatIntTerms2	= fixDateStrip2.GetInterestTerms();
	std::vector<double>& floatResetDate2	= fixDateStrip2.GetResetDates();

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
	std::vector<double> fixNotional2 (nbFixFlows2,1.0);
	std::vector<double> floatNotional2 (nbFixFlows2,1.0);

	int idxend1 = IdxFromValue(itsEndTimes, (*floatPayTimes1)[floatPayTimes1->size()-1], 7);
	int idxend2 = IdxFromValue(itsEndTimes, (*floatPayTimes2)[floatPayTimes2->size()-1], 7);

	if(idxend1 == -1 || idxend2 == -1)
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar swaption schedule not in calibration schedule");
	}

	if((idxend1 - idx + 1 != floatPayTimes1->size()) || (idxend2 - idx + 1 != floatPayTimes2->size()))
	{
		ARM_THROW(ERR_INVALID_ARGUMENT,"ARM_BGMSV1F::ComputeUnderlyingVarCovar float schedule must have same period as calibration schedule");
	}

	double fwdrate1 = (*ARM_PricingModelIR::SwapRate(curveCcy, evalTime, startTime1, endTime1, 
								*fixPayTimes1, *fixPayPeriods1, 
								*floatStartTimes1, *floatPayTimes1, *floatIntTerms1, *floatPayTimes1, *floatIntTerms1,
								std::vector<double>(1,0.),true, states))[0];


	double annuity1 = (*ARM_PricingModelIR::AnnuityWithNominal(curveCcy, evalTime, *fixPayTimes1, *fixPayPeriods1, fixNotional1, states))[0];

	double fwdrate2	= (*ARM_PricingModelIR::SwapRate(curveCcy, evalTime, startTime2, endTime2,
								*fixPayTimes2, *fixPayPeriods2, 
								*floatStartTimes2, *floatPayTimes2, *floatIntTerms2, *floatPayTimes2, *floatIntTerms2,
								std::vector<double>(1,0.),true, states))[0];

	double annuity2	= (*ARM_PricingModelIR::AnnuityWithNominal(curveCcy, evalTime, *fixPayTimes2, *fixPayPeriods2, fixNotional2, states))[0];

	std::vector<double> weight1, weight2;

	GetWeightSwaptionApprox(curveCcy, evalTime, floatNotional1, *floatStartTimes1, *floatPayTimes1, *floatIntTerms1,
							idx, fwdrate1, annuity1, weight1);

	GetWeightSwaptionApprox(curveCcy, evalTime, floatNotional2, *floatStartTimes2, *floatPayTimes2, *floatIntTerms2,
							idx, fwdrate2, annuity2, weight2);

	ARM_IntVector idxs(0);

	varcovar[0]	= GetEquivSwaptionLevel(idx, idxs, weight1);
	varcovar[1]	= GetEquivSwaptionLevel(idx, idxs, weight2);
	double cov	= GetEquivSwaptionLevel(idx, idxs, weight1, &weight2);
	varcovar[2]	= cov * cov;
}


ARM_VectorPtr ARM_BGMSV1F::VanillaSpreadOptionLet(const string& curveName,
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
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_BGMSV1F::VanillaSpreadOptionLet implemented only for calibrated Spread" );
}

/*
void ARM_BGMSV1F::getCalibrationData(std::vector<double>& resettimes, std::vector<double>& forwards, ARM_VectorVector& strikes, ARM_VectorVector& vols, std::vector<double>& atmvols, std::vector<double>& weights)
{
	double asof = GetAsOfDate().GetJulian();

	int i, k, size = itsResetTimes.size();
	int rsize = 0;
	double delta, prevtime = 0;
	
	itsCalibWeight.resize(size, 1);

	for(k = 0; k < size; k++)
	{
		delta = k == 0 ? itsResetTimes[k] : itsResetTimes[k] - prevtime;

		if(fabs(itsCalibSecDensities[k]->getWeight()) < K_DOUBLE_TOL)
		{
			itsCalibWeight[k] = 0;
		}
		else
		{
			if(itsCalibSecDensities[k]->getInterestTerms().size() > 1 && delta < 360)
			{
				itsCalibWeight[k] = 0;
			}
			else if(itsCalibSecDensities[k]->getInterestTerms().size() == 1 && delta < 150)
			{
				itsCalibWeight[k] = 0;
			}
			else
			{
				itsCalibWeight[k] = 1;
				rsize ++;
				prevtime = itsResetTimes[k];
			}
		}
	}

	if(rsize == 0)
	{
		rsize ++;
		itsCalibWeight[size-1] = 1;
	}

	resettimes.resize(rsize);
	forwards.resize(rsize);
	strikes.resize(rsize);
	vols.resize(rsize);
	atmvols.resize(rsize);
	weights.resize(rsize);

	double price;

	for(i = 0, k = 0; i < size; i++)
	{
		if(itsCalibWeight[i] == 0) 
		{
			price	= itsCalibSecDensities[i]->getDensityFunctor()->Call_Option(itsFwdRate[i], itsFwdRate[i], itsResetTimes[i] / 365.);
			itsATMVols[i] = VanillaImpliedVol_BS(itsFwdRate[i], itsFwdRate[i], itsResetTimes[i] / 365., price, 1);
		}
		else
		{
			strikes[k] = new std::vector<double>;
			vols[k] = new std::vector<double>;

			resettimes[k] = itsResetTimes[i];
			forwards[k] = itsFwdRate[i];
			weights[k] = 1.0;

			getkthStrikesAndVols(i, (*strikes[k]), (*vols[k]), atmvols[k]);

			itsATMVols[i] = atmvols[k];

			k++;
		}
	}
}


void ARM_BGMSV1F::finishCalibration(std::vector<double> * calibTimes, ARM_SmileCalibration_Params_Heston& params)
{
	ARM_SmileCalibration_Heston func;

	int k, size = itsResetTimes.size(), idx = 0;
	std::vector<double> calibShifts(size), calibLevels(size), calibRhos(size);
	std::vector<double> strikes, mktvols;
	double atmvol, locshift, locrho;
	bool calib, locrhoCalib = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLocalRhoCalib();
	int idxinf, idxsup;
	
	double infshift = 0., supshift = 0.;

	for(k = 0; k < size; k++)
	{
		if(itsCalibWeight[k] == 1)
		{
			calibShifts[k]	= params.shifts()[idx];
			calibLevels[k]	= params.levels()[idx];
			calibRhos[k]	= locrhoCalib ? params.rhos()[idx] : params.rho();
			idx ++;
			if(fabs(supshift) < K_DOUBLE_TOL) supshift = calibShifts[k];
			infshift = calibShifts[k]; 
			idxsup = k;
		}
	}

	for(k = idxsup; k < size; k++) calibShifts[k] = infshift;

	for(k = 0; k < size; k++)
	{
		if(itsCalibWeight[k] == 0)
		{
			if(k == 0 || k == size - 1)
				calib = true;
			else 
			{	
				idxinf = -1;
				idxsup = size;

				for(idxinf = k-1; idxinf > -1; idxinf --)
				{
					if(itsCalibWeight[idxinf] == 1)
					{
						break;
					}
				}

				for(idxsup = k+1; idxsup < size; idxsup++)
				{
					if(itsCalibWeight[idxsup] == 1)
					{
						break;
					}
				}

				if(idxinf > -1 && idxsup < size)
					calib = params.calibshift() == false ? true : false;
				else 
					calib = true;
			}
			
			if(calib && (k == 0 || k == size-1))
			{
				// récupération des strikes et des vols
				getkthStrikesAndVols(k, strikes, mktvols, atmvol);

				locshift = k == 0 ? supshift : infshift;
				
				locrho = locrhoCalib ? k == 0 ? params.rhos()[0] : params.rhos()[params.rhos().size()-1] : params.rho();

				// calibration du shift et du level
				ARM_SmileCalibration_Params_Heston locparams(params.v0(), params.kappa(), params.theta(), locrho, params.nu(), 
					locshift, 1., false, false, false, false, params.calibshift());
				
				func.Init(itsResetTimes[k] / 365., itsFwdRate[k], mktvols, strikes, true, itsFwdRate[k], atmvol, &locparams);

				bool success = func.Calibrate();

				if(success == false && locparams.calibshift() == true)
				{
					locparams.calibshift() = false;
					func.Calibrate();

					calibShifts[k]	= locparams.shifts()[0];
					calibLevels[k]	= locparams.levels()[0];
					calibRhos[k]	= locparams.rho();
				}
				else
				{
					calibShifts[k]	= locparams.shifts()[0];
					calibLevels[k]	= locparams.levels()[0];
					calibRhos[k]	= locparams.rho();
				}
			}
			else if(calib)
			{
				if(idxinf == -1) idxinf = k-1;
				if(idxsup == size) idxsup = k+1;

				// récupération des strikes et des vols
				getkthStrikesAndVols(k, strikes, mktvols, atmvol);

				// interpolation linéaire en temps pour le shift
				locshift = calibShifts[k] = calibShifts[idxinf] + (itsResetTimes[k] - itsResetTimes[idxinf]) * (calibShifts[idxsup] - calibShifts[idxinf]) / (itsResetTimes[idxsup] - itsResetTimes[idxinf]);

				// interpolation linéaire de la correl
				locrho = calibRhos[k] = calibRhos[idxinf] + (itsResetTimes[k] - itsResetTimes[idxinf]) * (calibRhos[idxsup] - calibRhos[idxinf]) / (itsResetTimes[idxsup] - itsResetTimes[idxinf]);

				// calibration du shift et du level
				ARM_SmileCalibration_Params_Heston locparams(params.v0(), params.kappa(), params.theta(), locrho, params.nu(), 
					locshift, 1., false, false, false, false, false);

				func.Init(itsResetTimes[k] / 365., itsFwdRate[k], mktvols, strikes, true, itsFwdRate[k], atmvol, &locparams);

				bool success = func.Calibrate();

				if(success == false)
				{
					// interpolation linéaire sur les forwards pour le level
					calibLevels[k] = calibLevels[idxinf] + (itsFwdRate[k] - itsFwdRate[idxinf]) * (calibLevels[idxsup] - calibLevels[idxinf]) / (itsFwdRate[idxsup] - itsFwdRate[idxinf]);
				}
				else
				{
					calibLevels[k] = locparams.levels()[0];
				}
				
			}
			else
			{
				// interpolation linéaire en temps pour le shift
				calibShifts[k] = calibShifts[idxinf] + (itsResetTimes[k] - itsResetTimes[idxinf]) * (calibShifts[idxsup] - calibShifts[idxinf]) / (itsResetTimes[idxsup] - itsResetTimes[idxinf]);

				// interpolation linéaire sur les forwards pour le level
				calibLevels[k] = calibLevels[idxinf] + (itsFwdRate[k] - itsFwdRate[idxinf]) * (calibLevels[idxsup] - calibLevels[idxinf]) / (itsFwdRate[idxsup] - itsFwdRate[idxinf]);

				// interpolation linéaire pour la correl
				calibRhos[k] = calibRhos[idxinf] + (itsResetTimes[k] - itsResetTimes[idxinf]) * (calibRhos[idxsup] - calibRhos[idxinf]) / (itsResetTimes[idxsup] - itsResetTimes[idxinf]);

			}
		}
	}

	std::vector<double> calibTheta(size, params.theta()), calibvvol(size, params.nu());

	((ARM_ModelParamsBGMSV1F*) GetModelParams())->UpdateParamValues(calibTimes, &calibShifts, &calibLevels, &calibRhos, &calibvvol, &calibTheta, params.v0());
}


void ARM_BGMSV1F::getkthStrikesAndVols(int k, std::vector<double>& strikes, std::vector<double>& vols,double& atmvol)
{
	double price;
	// vol monnaie
	price	= itsCalibSecDensities[k]->getDensityFunctor()->Call_Option(itsFwdRate[k], itsFwdRate[k], itsResetTimes[k] / 365.);
	atmvol	= VanillaImpliedVol_BS(itsFwdRate[k], itsFwdRate[k], itsResetTimes[k] / 365., price, 1);

	int i, nbStrikes = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetStddevCalib().size();

	strikes.resize(nbStrikes);
	vols.resize(nbStrikes);

	for(i = 0; i < nbStrikes; i++)
	{
		strikes[i] = itsFwdRate[k] * exp(((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetStddevCalib()[i] * sqrt(itsResetTimes[k] / 365.) * atmvol);
	}

	for(i = 0; i < nbStrikes; i++)
	{
		price	= itsCalibSecDensities[k]->getDensityFunctor()->Call_Option(strikes[i], itsFwdRate[k], itsResetTimes[k] / 365.);
		vols[i]	= VanillaImpliedVol_BS(itsFwdRate[k], strikes[i], itsResetTimes[k] / 365., price, 1);
	}
}
*/

void ARM_BGMSV1F::GetWeightSwaptionApprox(const string& curveName, double evalTime, const std::vector<double>& floatNotional, 
										  const std::vector<double>& floatStartTimes, const std::vector<double>& floatEndTimes, const std::vector<double>& floatPayPeriods,
										  int idx, double fwdrate, double annuity, std::vector<double>& weight) const
{
	weight.resize(floatStartTimes.size());

	for(int k = 0; k < floatStartTimes.size(); k++)
	{
		weight[k] = floatNotional[k] * ((*DiscountFactor(curveName, evalTime, floatStartTimes[k], ARM_PricingStatesPtr(NULL)))[0] - (*DiscountFactor(curveName, evalTime, floatEndTimes[k], ARM_PricingStatesPtr(NULL)))[0]) / (annuity * fwdrate);
	}
}

double ARM_BGMSV1F::GetEquivSwaptionShift(int idx, const ARM_IntVector& idxs, std::vector<double>& weight) const
{
	double shift = 0.;
	for(int k = 0; k < weight.size(); k++)
	{
		double shk =  idxs.size() == 0 ? ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idx+k) : ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetShift(idxs[k]);

		shift += weight[k] * (shk - 1.) / shk;

		// modif des poids pour la suite
		weight[k] /= shk;
	}

	shift = 1. / (1. - shift);

	weight *= shift;

	return shift;
}

double ARM_BGMSV1F::GetEquivSwaptionLevel(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, const std::vector<double> * optweight) const
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
				voli = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(i + idx);
			else
				voli = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idxs[i]);

			for(j = 0; j < tmpweight->size(); j++)
			{
				if(idxs.size() == 0)
					volj = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(j + idx);
				else
					volj = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idxs[j]);

				if(idxs.size() == 0)
					corr = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->RateRateCorrel(itsResetTimes[k], itsResetTimes[i + idx], i + idx, itsResetTimes[j + idx], j + idx);
				else
					corr = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->RateRateCorrel(itsResetTimes[k], itsResetTimes[idxs[i]], idxs[i], itsResetTimes[idxs[j]], idxs[j]);

				var += weight[i] * (*tmpweight)[j] * voli * volj * corr * dt;
			}
		}
	}

	return sqrt(var * 365. / itsResetTimes[idx]);
}

double ARM_BGMSV1F::GetEquivSwaptionRho(int idx, const ARM_IntVector& idxs, const std::vector<double>& weight, double swlevel) const
{
	int k;
	double rhok, volk, rho = 0.;

	for(k = 0; k < weight.size(); k++)
	{
		if(idxs.size() == 0)
		{
			rhok = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(k + idx);
			volk = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(k + idx);
		}
		else
		{
			rhok = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetRho(idxs[k]);
			volk = ((ARM_ModelParamsBGMSV1F*) GetModelParams())->GetLevel(idxs[k]);
		}

		rho += volk * rhok * weight[k];
	}

	rho /= swlevel;

	return rho;
}

CC_END_NAMESPACE()
