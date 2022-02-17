/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
*/

/// gpmodels
#include "gpmodels/SVBGM.h"
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


CC_BEGIN_NAMESPACE( ARM )


void ARM_ShiftedSABR::calibre(double resetTime, double forward, double atmvol,
							  ARM_GP_Vector * mktVols, ARM_GP_Vector * strikes,
							  double shift, double rho, double nu)
{
	bool calibshift = true; //shift == -999.;
	bool calibrho	= true; //rho == -999.;
	bool calibnu	= true; //nu == -999.;

//	build(resetTime, forward, forward, 1, atmvol, DefaultAlpha, calibrho ? 0. : rho, calibnu ? 0.1 : nu, calibshift ? 1. : shift);
	build(resetTime, forward, forward, 1, atmvol, DefaultAlpha, rho, nu, shift);

	double eps = 1e-3;

	ObjectiveFunction func(this, calibshift, calibrho, calibnu, eps, strikes, mktVols);

	ARM_GP_Vector fguess(calibshift + calibrho + calibnu);
	ARM_GP_Vector ubound(fguess.size(), 100.), lbound(fguess.size(), -100.);
	ARM_GP_Vector fatfguess(strikes->size(), 0.);

	int k = 0;
	
	if(calibshift)	fguess[k++] = shift;
	if(calibrho)	fguess[k++] = rho;
	if(calibnu)		fguess[k]	= nu;

	double info[LM_INFO_SZ];

	int status = LEVMARQMinization_WithNumDerivatives(func, fguess, fatfguess, info);
}



////////////////////////////////////////////////////
///	Class  : ARM_SVBGM
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_SVBGM::ARM_SVBGM( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, bool AllowInterpol, bool ComputeDrift, bool Proxy)
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
	itsFromShiftedToRate(NULL),
	itsComputeDriftAdj(ComputeDrift),
	itsEigenValues(0),
	itsModelLocalRealVar(0),
	itsCalibSecDensities(0),
	itsProxyStatus(Proxy)
{
	if(itsProxyStatus) itsComputeDriftAdj = false;
}


////////////////////////////////////////////////////
///	Class  : ARM_SVBGM
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_SVBGM::ARM_SVBGM(const ARM_SVBGM& rhs)
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
	itsFromShiftedToRate(rhs.itsFromShiftedToRate),
	itsComputeDriftAdj(rhs.itsComputeDriftAdj),
	itsProxyStatus(rhs.itsProxyStatus)
{
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Vector>( rhs.itsEigenValues, itsEigenValues );
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( rhs.itsModelLocalRealVar, itsModelLocalRealVar);
	DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (rhs.itsCalibSecDensities, itsCalibSecDensities);
}	


ARM_SVBGM& ARM_SVBGM::operator = (const ARM_SVBGM& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_SVBGM();

		new (this) ARM_SVBGM (rhs);
	}

	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_SVBGM
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_SVBGM::~ARM_SVBGM()
{
	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
	DeletePointorVector<ARM_GP_Matrix>( itsModelLocalRealVar );
}

bool ARM_SVBGM::ValidateModelParams(const ARM_ModelParams& params) const
{
	if(params.DoesModelParamExist(ARM_ModelParamType::Shift) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: shift structure is needed for SVBGM" );

	if(params.DoesModelParamExist(ARM_ModelParamType::Alpha) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: alpha structure is needed for SVBGM" );

	if(params.DoesModelParamExist(ARM_ModelParamType::QParameter) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: Q structure (rho sabr) is needed for SVBGM" );

	if(params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: VolOfVol structure (nu sabr) is needed for SVBGM" );

	if(params.DoesModelParamExist(ARM_ModelParamType::BetaCorrelation) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: beta correlation structure is needed for SVBGM" );

	if(params.DoesModelParamExist(ARM_ModelParamType::CrossFactor) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: CrossFactor structure (rate(i) / vol(k) correlation) is needed for SVBGM" );
	
	if(params.DoesModelParamExist(ARM_ModelParamType::ReCorrelation) == false)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateModelParams: ReCorrelation structure (vol(i) / vol(k) correlation) is needed for SVBGM" );

	return true;
}

void ARM_SVBGM::SetNumeraire(const ARM_NumerairePtr& numerairePtr)
{
	if(	numerairePtr->GetType() != ARM_Numeraire::TerminalZc )
    {		
	   ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::SetNumeraire: only numeraire supported is Terminal ZC" );
    }
   	ARM_PricingModel::SetNumeraire(numerairePtr);
}

void ARM_SVBGM::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	bool isNumerical = (calibMethod.GetMethodType()== ARM_CalibMethodType::Numerical);

	if (!itsCalibratedStatus && !isNumerical)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ValidateCalibMethod: should calib functionals first!!!" );

	calibMethod.DefaultValidateWithModel(*this);
}

void ARM_SVBGM::PreProcessing(ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

void ARM_SVBGM::PostProcessing(const ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

void ARM_SVBGM::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{
	// TO DO 
}

ARM_PricingStatesPtr ARM_SVBGM::FirstPricingStates(size_t bucketSize) const
{
	const size_t nbPayoffs			= 0;
	size_t nbModelStates			= 2 * itsNbEffectiveReset;
	size_t factorsNb				= FactorCount();

	ARM_PricingStatesPtr initStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,nbModelStates,nbPayoffs,factorsNb) );

	for(int i = 0; i < nbModelStates / 2; i++)
	{
		// initialisation des taux shiftés
		double initRateValue	= itsFwdRate[i] / ((ARM_ModelParamsSVBGM*) GetModelParams())->GetShift(i);
		double initVolValue		= ((ARM_ModelParamsSVBGM*) GetModelParams())->GetAlpha(i);

		for(int k = 0; k < bucketSize; k++) 
		{
			initStates->SetModelState(k, i, initRateValue);
			initStates->SetModelState(k, i + nbModelStates / 2, initVolValue);
		}
	}

	return initStates;
}

ARM_PricingStatesPtr ARM_SVBGM::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
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

ARM_PricingStatesPtr ARM_SVBGM::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	if(itsCalibrationStatus)
	{		
		if(itsCalibratedStatus == false)
		{
			int k, size = itsResetTimes.size();
			
			ARM_GP_Vector calibshift(size), calibalpha(size), calibrho(size), calibnu(size);

			for(k = 0; k < size; k++)
			{
				kthLocalCalibration(k, calibshift[k], calibalpha[k], calibrho[k], calibnu[k]);
			}

			((ARM_ModelParamsSVBGM*) GetModelParams())->UpdateParamValues(&itsResetTimes, &calibshift, &calibalpha, &calibrho, &calibnu);
			
			
			itsFromShiftedToRate.resize(size);

			// pour passer des libors shiftés aux libors : fwd x (shift - 1) / shift
			for(k = 0; k < size; k++) 
			{
				itsFromShiftedToRate[k] = (((ARM_ModelParamsSVBGM*) GetModelParams())->GetShift(k) - 1.) * itsFwdRate[k] / ((ARM_ModelParamsSVBGM*) GetModelParams())->GetShift(k);
			}

			itsCalibratedStatus = true;
		}

		return ARM_PricingStatesPtr(new ARM_PricingStates(0));
	}

	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );

	return newStates;
}

ARM_VectorPtr ARM_SVBGM::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
// FIXMEFRED: mig.vc8 (30/05/2007 16:22:22):cvast
	return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,getTerminalTime()));
}

ARM_GP_Vector* ARM_SVBGM::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	return static_cast<ARM_GP_Vector*>(itsResetTimes.Clone());
}

void ARM_SVBGM::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	// TO DO    
}

string ARM_SVBGM::toString(const string& indent,const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "\n\n";
    os << indent << "SVBGM : full Stochastic Volatilty BGM with log normal volatility for each underlying \n";
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

void ARM_SVBGM::NumMethodStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const
{
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	
    localVariances.resize((timeStepsSize-1)*(modelNb+1));
	for(int i = 0; i < timeStepsSize-1 ; ++i)
	{
		int factorsNb = itsEigenValues[i]->size();

		/// initialize everything
		localVariances[i] = new ARM_GP_Matrix(factorsNb,factorsNb,0.);

		for (int k = 0; k < factorsNb; k++)
			(*localVariances[i])(k,k) = (*itsEigenValues[i])[k];
	}
}

void ARM_SVBGM::NumMethodStateGlobalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& globalVariances ) const
{
	// USELESS !!!!! But we have to compute it for the sampler
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::NumMethodStateGlobalVariances: globalVariances.size() != offsetIndex" );
#endif
    globalVariances.resize(timeStepsSize*(modelNb+1));
	for (int i=0;i<timeStepsSize;i++)
	{
		globalVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(2*totalFwds,1.0);
	}
}


void ARM_SVBGM::ModelStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps)
{
	ARM_MatrixVector modStateLocalVars;
	ARM_MatrixVector auxLocalCov;

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

	((ARM_ModelParamsSVBGM*) GetModelParams())->checkFactorCount(itsNbEffectiveReset);

	ModelStateLocalVariances( timeSteps, auxLocalCov);
	DuplicateCloneablePointorVectorInPlace<ARM_GP_Matrix>( auxLocalCov, itsModelLocalRealVar);
	
	int totalFwds		= itsNbEffectiveReset;
	int factorsNb		= FactorCount();
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;

	int i, j, k, startTimePos = 0;

	modStateLocalVars.resize((timeStepsSize-1)*(modelNb+1));

	DeletePointorVector<ARM_GP_Vector>( itsEigenValues );
	itsEigenValues.resize((timeStepsSize-1)*(modelNb+1));

	int nbFwds, idx = 0, currfactorsNb = factorsNb;
	int maxfactorsNb = 0;
	double minratio = ((ARM_ModelParamsSVBGM*) GetModelParams())->GetMinRatio();
	bool resizefactors = minratio > 1. || fabs(minratio - 1.) < K_DOUBLE_TOL ? false : true;

	// on commence par calculer toutes les matrices, et les vecteurs propres
	ARM_GP_Vector * eigenValues = new ARM_GP_Vector[timeStepsSize - 1];
	ARM_GP_Matrix ** ACPMatrix = new ARM_GP_Matrix * [timeStepsSize - 1];
	ARM_GP_Vector factorsize(timeStepsSize-1, factorsNb);

	for(i = 0; i < timeStepsSize - 1; i++)
	{
		toTime  = timeSteps[i+1];

		/// get the first bigger time
		while(idx < itsResetTimes.size()
			&& itsResetTimes[idx] < toTime )
			++idx;

		nbFwds = totalFwds - idx;

		eigenValues[i].resize(2 * nbFwds,0.);
		ACPMatrix[i] = ACPTransformation(auxLocalCov[i],eigenValues[i]);

		// calcul de la dimension effective
		if(resizefactors)
		{
			// variance total
			double varTotal = 0.;
			for(k = 0; k < 2 * nbFwds; k++) 
			{
				if(eigenValues[i][k] > 0.) varTotal += eigenValues[i][k];
			}

			double varExpliq = 0.;
			for(k = 0; k < 2 * nbFwds; k++)
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

	((ARM_ModelParamsSVBGM*) GetModelParams())->SetFactorCount(maxfactorsNb);
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

		currfactorsNb = factorsize[i] > 2 * nbFwds ? 2 * nbFwds : factorsize[i];

		/// initialize everything
		modStateLocalVars[i] = new ARM_GP_Matrix(2 * nbFwds, currfactorsNb);

		itsEigenValues[i] = new ARM_GP_Vector(currfactorsNb);
		
		
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
		
		for (j=0;j<2*nbFwds;j++)
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
		
		for(k=0;k<currfactorsNb;++k)
		{
			(*itsEigenValues[i])(k) = eigenValues[i][k];
			for(j=0; j<2*nbFwds; ++j)
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

void ARM_SVBGM::NumMethodStateLocalVariancesAndStdDev( const ARM_GP_Vector& timeSteps )
{

	ARM_MatrixVector numMethodStateLocalVars;
	
	/// computes the local variance
	NumMethodStateLocalVariances( timeSteps, numMethodStateLocalVars);
	
	/// set the result
	SetNumMethodStateLocalVars(numMethodStateLocalVars);
	
}

void ARM_SVBGM::ModelStateLocalVariances(const ARM_GP_Vector& timeSteps, ARM_MatrixVector& localVariances ) const
{
	int totalFwds		= itsNbEffectiveReset;
	int timeStepsSize	= timeSteps.size();
	int modelNb			= GetModelNb();
	int offsetIndex		= (timeStepsSize - 1) * modelNb;
	
	double fromTime		= timeSteps[0];
	double toTime;


#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ModelStateLocalVariances: localVariances.size() != offsetIndex" );
#endif

#if defined(__GP_STRICT_VALIDATION)
	if( fromTime != 0.0 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::ModelStateLocalVariances: first time step != 0" );
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

		localVariances[i] = new ARM_GP_Matrix(2 * nbFwds, 2 * nbFwds);	

		// les correls taux taux
		for(j = 0; j < nbFwds; ++j)
		{
			for(k = j; k < nbFwds; ++k)
				(*localVariances[i])(j,k) = (*localVariances[i])(k,j) = ((ARM_ModelParamsSVBGM*) GetModelParams())->RateRateCorrel(toTime, itsResetTimes[k + idx], k + idx, itsResetTimes[j + idx], j + idx);
		}
		
		// les correls vols vols
		for(j = 0; j < nbFwds; j++)
		{
			for(k = j; k < nbFwds; k++)
				(*localVariances[i])(j + nbFwds, k + nbFwds) = (*localVariances[i])(k + nbFwds, j + nbFwds) = ((ARM_ModelParamsSVBGM*) GetModelParams())->VolVolCorrel(toTime, k + idx, j + idx);
		}

		// les correls taux vols
		for(j = 0; j < nbFwds; j++)
		{
			for(k = j; k < nbFwds; k++)
			{
				double correl = ((ARM_ModelParamsSVBGM*) GetModelParams())->RateVolCorrel(toTime, k + idx, j + idx);

				(*localVariances[i])(j + nbFwds, k) = (*localVariances[i])(k + nbFwds, j) = correl;
				(*localVariances[i])(j, k + nbFwds) = (*localVariances[i])(k, j + nbFwds) = correl;
			}
		}

		fromTime = toTime;
	}

}

void ARM_SVBGM::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
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

	ARM_GP_Vector x(2 * aliveFwds);

	const ARM_MatrixVector& modelLocalStdev	= GetModelStateLocalVars();

	int eigensNb	= (*modelLocalStdev[timeIndex]).cols(); //factorsNb > 2 * aliveFwds ? 2 * aliveFwds : factorsNb;
	
	double rcov, vcov;
	double rdrift = 0., vdrift = 0., adj = 0.;

	for(n = 0; n < statesNb; n++)
	{
		// corrélation des gaussiennes
		for(i = 0; i < 2 * aliveFwds; i++)
		{
			x[i] = 0.;
			for(j = 0; j < eigensNb; j++)
			{
				x[i] += states->GetNumMethodState(n, j) * (*modelLocalStdev[timeIndex])(i,j);
			}
		}

		for(i = iLast; i >= iFirst; i--)
		{
			double voli	= states->GetModelState(n, i + fwdsNb);
			double libi	= states->GetModelState(n, i);
			double vvol = ((ARM_ModelParamsSVBGM*) GetModelParams())->GetNu(i);
			double vstd	= vvol * racdt;

			if(itsComputeDriftAdj)
			{
				rdrift = vdrift = 0.;

				for(j = i + 1; j <= iLast; j++)
				{
					rcov	= (*itsModelLocalRealVar[timeIndex])(i - iFirst, j - iFirst) * dt;				// corrélation taux i, taux j
					vcov	= (*itsModelLocalRealVar[timeIndex])(i + aliveFwds - iFirst, j - iFirst) * dt;	// corrélation taux i, vol j
					
					adj		= states->GetModelState(n, j) * itsDelta[j] * states->GetModelState(n, j + fwdsNb)
							/ (1. + itsDelta[j] * (states->GetModelState(n, j) + itsFromShiftedToRate[j]));

					rdrift	-= adj * rcov;
					vdrift	-= adj * vcov;
				}

				rdrift *= voli;
				vdrift *= vvol;
			}

			libi	*= exp(rdrift - 0.5 * voli * voli * dt + voli * racdt * x[i - iFirst]);
			voli	*= exp(vdrift - 0.5 * vstd * vstd + vstd * x[i + aliveFwds - iFirst]);

			states->SetModelState(n, i, libi);
			states->SetModelState(n, i + fwdsNb, voli);
		}
	}
}

///////////////////////
//
//		les prix
//
///////////////////////

void ARM_SVBGM::setNumericalModelFitter(ARM_NumericalModelFitter * numericalModelFitter)
{
	if(numericalModelFitter)
	{
		ARM_VanillaSecDensityPtrVector densityVector = numericalModelFitter->getCalibSecDensities();
		size_t size = densityVector.size();

		ARM_GP_Vector resetDates(size); 
		ARM_GP_Vector startDates(size); 
		ARM_GP_Vector endDates(size); 
		ARM_GP_Vector delta(size); 
		ARM_GP_Vector fwdRate(size);

		for (size_t k = 0 ; k < size ; k++)
		{
			resetDates[k]			= densityVector[k]->getResetDate();
			startDates[k]			= densityVector[k]->getStartDate();
			endDates[k]				= densityVector[k]->getEndDate();

			delta[k]				= densityVector[k]->getInterestTerms()[0];
			fwdRate[k]				= densityVector[k]->ComputeForwardRate();

			if(densityVector[k]->getInterestTerms().size() > 1)
			{
				itsComputeDriftAdj	= false;
				itsProxyStatus		= true;
			}
		}
		
		double asof = GetAsOfDate().GetJulian();

		itsResetTimes			= resetDates - asof;
		itsStartTimes			= startDates - asof;
		itsEndTimes				= endDates   - asof;
		itsDelta				= delta;
		itsFwdRate				= fwdRate;

		computeWeights();

		((ARM_ModelParamsSVBGM*) GetModelParams())->checkFactorCount(size);

		DuplicateCloneablePtrVectorInPlace<ARM_VanillaSecurityDensity> (densityVector, itsCalibSecDensities);

		itsCalibratedStatus = false;
	}
}

void ARM_SVBGM::computeWeights() 
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
				k++;
			}
			itsWeight(j,k) = 1.;
			term = itsEndTimes[k];
			k++;
		}

		itsWeight(j,j) = 1.;
	}
}

ARM_VectorPtr ARM_SVBGM::DiscountFactor(const string& curveName, double evalTime, double maturityTime, 
										const ARM_PricingStatesPtr& states) const
{
	if(itsAllowInterpol)
		return ForwardDiscountFactor(curveName,evalTime,evalTime,maturityTime,states);
	else
		return DiscountFactorNoInterpol(curveName,evalTime,maturityTime,states);
}

ARM_VectorPtr ARM_SVBGM::ForwardDiscountFactor( const string& curveName, double evalTime, double startTime, double endTime, const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();

	if( evalTime < K_DOUBLE_TOL || itsProxyStatus == true)
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN) ));
	}

	if( endTime > itsEndTimes[itsEndTimes.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: endTime after schedule's last date");

	if( evalTime < itsResetTimes[0] )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime before first reset");

	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );

	ARM_GP_Vector::iterator iter;
	ARM_GP_Vector::iterator iter1;
	ARM_GP_Vector::iterator iter2;

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
				(*iter) *= mu / ( 1. + coeff * ( states->GetModelState(k, idx) + itsFromShiftedToRate[idx]) );
	}

	if ( (endTime > itsStartTimes[LastIdx]) && (LastIdx >= FirstIdx))
	{
		double theta	= (endTime - itsStartTimes[LastIdx])/360.;

		double mu		= GetZeroCurve()->DiscountPrice(endTime/K_YEAR_LEN)
						/ GetZeroCurve()->DiscountPrice(itsStartTimes[LastIdx]/K_YEAR_LEN)
						* ( 1 + theta * itsFwdRate[LastIdx]);
				
		int k;

		for(iter = result->begin(), k = 0 ; iter != result->end() ; ++iter, k++)
			(*iter) *= mu / ( 1. + theta * (states->GetModelState(k, LastIdx) + itsFromShiftedToRate[LastIdx]));
	}

	return static_cast<ARM_VectorPtr>(result);
}

ARM_VectorPtr ARM_SVBGM::ForwardDiscountFactorFromIdx( const string& curveName, double evalTime, size_t IdxFrom, size_t IdxTo , size_t modelNb, const ARM_PricingStatesPtr& states) const
{
	ARM_GP_Vector* result	= new ARM_GP_Vector( states->size(),1.0 );
	
	size_t size = states->size();

	for (size_t k = IdxFrom; k < IdxTo ; k++ )
	{
		if ( itsWeight(IdxFrom,k) == 1. )
		{
			for(size_t i = 0; i < size; i++)
			{
				(*result)[i] /= (1. + itsDelta[k] * (states->GetModelState(i, k) + itsFromShiftedToRate[k]));
			}
		}
	}

	return static_cast<ARM_VectorPtr>(result);
}

ARM_VectorPtr ARM_SVBGM::DiscountFactorNoInterpol(const string& curveName, double evalTime, double maturityTime, 
												  const ARM_PricingStatesPtr& states) const
{
	size_t modelNb = GetModelNb();
	
	if( evalTime < K_DOUBLE_TOL || itsProxyStatus == true)
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( maturityTime > itsEndTimes[itsEndTimes.size() -1] + 30 )
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::DiscountFactor: MaturityTime after schedule's last date");


	ARM_GP_Vector* result = new ARM_GP_Vector( states->size(),1.0 );
	
	if( DoesResetTimeExist( evalTime ) )
	{
		ARM_GP_Vector::iterator iter;
		ARM_GP_Vector::iterator iter1;
		ARM_GP_Vector::iterator iter2;

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
				(*iter) *= mu / (1. + theta * (states->GetModelState(k, StartIdx) + itsFromShiftedToRate[StartIdx]));
		}
		
		return static_cast<ARM_VectorPtr>(result);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SmiledLMM::DiscountFactor: evalTime not in schedule!!");		
	}
}

ARM_VectorPtr ARM_SVBGM::Libor( 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::Libor : case payTime != fwdEndTime not handled." );

	if(ExistsInVector(itsResetTimes, fwdResetTime))
	{
		int idx = IdxFromValue(itsResetTimes, fwdResetTime);

		if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			ARM_GP_Vector* result = new ARM_GP_Vector( states->size() );
			
			for(int i = 0; i < states->size(); i++) 
				(*result)[i] = states->GetModelState(i, idx) + itsFromShiftedToRate[idx];

			return static_cast<ARM_VectorPtr>(result);
		}
	}

	return DefaultLibor(curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
}

ARM_VectorPtr ARM_SVBGM::SwapRate(
		const string& curveName, 
        double evalTime,
		double floatStartTime,
        double floatEndTime, 
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
		const ARM_GP_Vector& fwdStartTimes,
        const ARM_GP_Vector& fwdEndTimes,
        const ARM_GP_Vector& fwdPayPeriods, 
		const ARM_GP_Vector& floatPayTimes,
        const ARM_GP_Vector& floatPayPeriods,
        const ARM_GP_Vector& margin,
        bool isDbleNotional,
        const ARM_PricingStatesPtr& states) const
{
	if(fabs(floatStartTime - fwdStartTimes[0]) < 7
	&& fabs(floatEndTime - fwdEndTimes[fwdEndTimes.size()-1]) < 7)
	{
		int idx = IdxFromValue(itsStartTimes, fwdStartTimes[0],7);

		if(idx > -1)
		{
			if(fabs(itsEndTimes[idx] - itsStartTimes[idx] + fwdStartTimes[0] - fwdEndTimes[fwdEndTimes.size()-1]) < 7)
			{
				ARM_GP_Vector* result = new ARM_GP_Vector( states->size() );
				
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
						(*result)[i] = states->GetModelState(i, idx) + itsFromShiftedToRate[idx] + (*num)[i]/(*denom)[i];
				}
				else
				{
					for(i = 0; i < states->size(); i++) 
						(*result)[i] = states->GetModelState(i, idx) + itsFromShiftedToRate[idx];
				}

				return static_cast<ARM_VectorPtr>(result);
			}
		}
	}

	return ARM_PricingModelIR::SwapRate(curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes, fixPayPeriods, 
		fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, margin, isDbleNotional, states);
}

ARM_VectorPtr ARM_SVBGM::VanillaCaplet(
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
	if ( fabs(evalTime-fwdResetTime) < K_DOUBLE_TOL )
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
	{
		if(evalTime > K_DOUBLE_TOL)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaCaplet implemented only for spot date" );
		}

		if(ExistsInVector(itsResetTimes, fwdResetTime))
		{
			int idx = IdxFromValue(itsResetTimes, fwdResetTime);

			if(itsCalibSecDensities[idx]->getInterestTerms().size() > 1)
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaCaplet called but it is a swaption" );
			}

			if(fabs(payTime - fwdEndTime) < 7)
			{
				if(fabs(fwdEndTime - fwdStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
				{
					ARM_ShiftedSABR vanillapricer;
					double caplet = vanillapricer.price(fwdResetTime / 365., itsFwdRate[idx], strikesPerState[0], capFloor,0., 
						((ARM_ModelParamsSVBGM*) GetModelParams())->GetAlpha(idx),
						((ARM_ModelParamsSVBGM*) GetModelParams())->GetRho(idx),
						((ARM_ModelParamsSVBGM*) GetModelParams())->GetNu(idx),
						((ARM_ModelParamsSVBGM*) GetModelParams())->GetShift(idx) );
					
					ARM_VectorPtr result(new ARM_GP_Vector(1, caplet * itsCalibSecDensities[idx]->ComputeLevel() * payNotional));

					return result;
				}
			}
		}

		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaCaplet implemented only for calibrated Libor" );
	}
}

ARM_VectorPtr ARM_SVBGM::VanillaSwaption(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaSwaption : notional, spread, strike must be const" );

	if(evalTime > K_DOUBLE_TOL)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaSwaption implemented only for spot date" );
	}

	if(ExistsInVector(itsResetTimes, swapResetTime))
	{
		double swapNotional = fixNotional[0];

		int idx = IdxFromValue(itsResetTimes, swapResetTime);

		if(fabs(floatEndTime - floatStartTime + itsStartTimes[idx] - itsEndTimes[idx]) < 7)
		{
			ARM_ShiftedSABR vanillapricer;
			double swopt = vanillapricer.price(swapResetTime / 365., itsFwdRate[idx], strikesPerState(0,0), callPut,0., 
				((ARM_ModelParamsSVBGM*) GetModelParams())->GetAlpha(idx),
				((ARM_ModelParamsSVBGM*) GetModelParams())->GetRho(idx),
				((ARM_ModelParamsSVBGM*) GetModelParams())->GetNu(idx),
				((ARM_ModelParamsSVBGM*) GetModelParams())->GetShift(idx) );
			
			ARM_VectorPtr result(new ARM_GP_Vector(1, swopt * itsCalibSecDensities[idx]->ComputeLevel() * swapNotional));

			return result;
		}
	}

	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaSwaption implemented only for calibrated CMS" );
}

ARM_VectorPtr ARM_SVBGM::VanillaSpreadOptionLet(const string& curveName,
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
	ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_SVBGM::VanillaSpreadOptionLet not implemented" );
}

void ARM_SVBGM::kthLocalCalibration(int k, double& shift, double& alpha, double& rho, double& nu)
{
	ARM_GP_Vector strikes, vols;
	double atmvol;

	getkthStrikesAndVols(k, strikes, vols, atmvol);

	ARM_ShiftedSABR calibsabr;

	double initshift, initrho, initnu;

	if(dynamic_cast<ARM_SABRDensityFunctor*>(&*itsCalibSecDensities[k]->getDensityFunctor()) != NULL)
	{
		initshift	= ((ARM_SABRDensityFunctor*)(&*itsCalibSecDensities[k]->getDensityFunctor()))->GetBeta();
		initrho		= ((ARM_SABRDensityFunctor*)(&*itsCalibSecDensities[k]->getDensityFunctor()))->GetRho();
		initnu		= ((ARM_SABRDensityFunctor*)(&*itsCalibSecDensities[k]->getDensityFunctor()))->GetNu();
	}
	else
	{
		initshift	= 1.;
		initrho		= 0.;
		initnu		= 0.;
	}

	calibsabr.calibre(itsResetTimes[k] / 365., itsFwdRate[k], atmvol, &vols, &strikes, initshift, initrho, initnu);

	shift	= calibsabr.GetShift();
	alpha	= calibsabr.GetAlpha();
	rho		= calibsabr.GetRho();
	nu		= calibsabr.GetNu();
}

void ARM_SVBGM::getkthStrikesAndVols(int k, ARM_GP_Vector& strikes, ARM_GP_Vector& vols,double& atmvol)
{
	double price;
	// vol monnaie
	price	= itsCalibSecDensities[k]->getDensityFunctor()->Call_Option(itsFwdRate[k], itsFwdRate[k], itsResetTimes[k] / 365.);
	atmvol	= VanillaImpliedVol_BS(itsFwdRate[k], itsFwdRate[k], itsResetTimes[k] / 365., price, 1);

	int i, nbStrikes = 11;
	/*
	double moneyNess = 0.5;
	double step = 2. * (1. - moneyNess) / (nbStrikes - 1);


	for(i = 0; i < nbStrikes; i++, moneyNess += step) strikes[i] = itsFwdRate[k] * moneyNess;
	*/

	strikes.resize(nbStrikes);
	vols.resize(nbStrikes);

	double L = 2.;
	double h = 2 * L / (nbStrikes - 1);
	double stdev = - L;

	for(i = 0; i < nbStrikes; i++, stdev += h)
	{
		strikes[i]		= itsFwdRate[k] * exp(stdev * sqrt(itsResetTimes[k] / 365.) * atmvol);
		price			= itsCalibSecDensities[k]->getDensityFunctor()->Call_Option(strikes[i], itsFwdRate[k], itsResetTimes[k] / 365.);
		vols[i]			= VanillaImpliedVol_BS(itsFwdRate[k], strikes[i], itsResetTimes[k] / 365., price, 1);
	}
}

CC_END_NAMESPACE()