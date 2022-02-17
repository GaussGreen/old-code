/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MarkovFunctional.cpp
 *
 *  \brief Markov Functional Model 1 Factor
 *
 *	\author  A Schauly
 *	\version 1.0
 *	\date August 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/MarkovFunctional.h"
#include "gpmodels/ModelParamsMF.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/curve.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/curvemodelparam.h"

/// gpclosedforms
#include "gpclosedforms/normal.h"
#include "gpclosedforms/gaussian_integrals.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/numerical.h"

/// gpnummethods
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/pdemethod.h"


#define CALIBGLPOINTSNB 8
#define PRICGLPOINTSNB	120
#define INVSQRT2PI		0.398942280401433
#define INVSQRT2		0.707106781186547


CC_BEGIN_NAMESPACE( ARM )

//////////////////////////////////////////////////////////
/// numerical recipes normal cdf, for testing purpose
//////////////////////////////////////////////////////////
double NormalCDF_NR (double x)
{
	x *= INVSQRT2;
	float t,z,ans;
	z=fabs(x);
	t=1.0/(1.0+0.5*z);
	ans=t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+
	t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+
	t*(-0.82215223+t*0.17087277)))))))));
	double erfc = x >= 0.0 ? ans : 2.0-ans;
	return 0.5 * ( 2.0 - erfc );   
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarkovFunctional::ARM_MarkovFunctional( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params )
:	ARM_PricingModelIR(zc,params), 
	itsNumericalModelFitter(NULL), 
	itsCalibrationStatus(false), 
	itsCurrentModelIsCalibrationModel(false), 
	itsCalibrationNbStdev(0.0),
	itsDiscountFactorMap(),
	itsPrevToTime(-1)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarkovFunctional::ARM_MarkovFunctional(const ARM_MarkovFunctional& rhs)
:	ARM_PricingModelIR(rhs), 
	itsNumericalModelFitter(CreateClone(rhs.itsNumericalModelFitter)), 
	itsCalibrationStatus(rhs.itsCalibrationStatus), 
	itsCurrentModelIsCalibrationModel(rhs.itsCurrentModelIsCalibrationModel), 
	itsCalibrationNbStdev(rhs.itsCalibrationNbStdev),
	itsDiscountFactorMap(rhs.itsDiscountFactorMap),
	itsPrevToTime(rhs.itsPrevToTime)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_MarkovFunctional::~ARM_MarkovFunctional()
{
	delete itsNumericalModelFitter;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_MarkovFunctional& ARM_MarkovFunctional::operator=(const ARM_MarkovFunctional& rhs)
{
	if (&rhs != this)
	{ 
		this->~ARM_MarkovFunctional();
		new (this) ARM_MarkovFunctional (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: setNumericalModelFitter
///	Returns: void
///	Action : sets Num model fitter + reset/start/end dates to DF map
////////////////////////////////////////////////////
void ARM_MarkovFunctional::setNumericalModelFitter( ARM_NumericalModelFitter* numericalModelFitter ) 
{
	if (numericalModelFitter)
	{
		/// set reset/start/end dates to DF map
		ARM_DateStripPtr sched = numericalModelFitter->getCalibSchedule();

		std::vector<double>& resetDates = sched->GetResetDates(); 
		std::vector<double>& startDates = sched->GetFlowStartDates(); 
		std::vector<double>& endDates	  = sched->GetFlowEndDates(); 

		double asof = GetAsOfDate().GetJulian();

		itsDiscountFactorMap.setResetTimes( ARM_GP_VectorPtr(new ARM_GP_Vector((*resetDates) - asof)) );
		itsDiscountFactorMap.setStartTimes( ARM_GP_VectorPtr(new ARM_GP_Vector((*startDates) - asof)) );
		itsDiscountFactorMap.setEndTimes  ( ARM_GP_VectorPtr(new ARM_GP_Vector((*endDates)   - asof)) );
		
	}
		
	/// assign model fitter
	itsNumericalModelFitter = numericalModelFitter; 
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MarkovFunctional::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    // Numerical method and numeraire are needed to go on
    ARM_NumMethodPtr numMethod=GetNumMethod();
	if( numMethod == ARM_NumMethodPtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in MF model!");

    /// test the numeraire and its type!
	ARM_NumerairePtr numeraire=GetNumeraire();
    if( numeraire == ARM_NumerairePtr(NULL) )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the MF model!");

	/// check that sampler type is MeanReverting
	ARM_SamplerBase* sampler = numMethod->GetSampler();
	if (sampler)
	{
		ARM_MeanRevertingSampler1D* sampler1D = dynamic_cast<ARM_MeanRevertingSampler1D*>(sampler);
		ARM_MeanRevertingSamplerND* samplerND = dynamic_cast<ARM_MeanRevertingSamplerND*>(sampler);
		if (!sampler1D && !samplerND)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Markov Functional model : a MeanRevertingSampler is required");
	}

	/// check that if we are calibrating, the num method is a PDE
	ARM_PDEMethod* pdeMethod = dynamic_cast<ARM_PDEMethod*>(&*numMethod);
	if (itsCalibrationStatus)
	{
		if (!pdeMethod)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Markov Functional model : calibration requires a PDE method");
		
		/// Store nb of stdev used in PDE for calibration
		itsCalibrationNbStdev = pdeMethod->getStdDevNb();
	}

	/// creates the model schedule (smart pointor for exception safety!)
	ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();

	CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps;

	if( !itsCalibrationStatus )
		ptimeSteps = CC_NS(std,auto_ptr)<std::vector<double>>( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );
	else
		ptimeSteps = CC_NS(std,auto_ptr)<std::vector<double>>( static_cast<ARM_GP_Vector*> (itsDiscountFactorMap.getResetTimes()->Clone()) );
	
	/// Initialise the numeraire
	numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos,ARM_GP_VectorPtr( new std::vector<double>( 1, itsDiscountFactorMap.getTerminalTime() )));
    	
	/// Little magouille so that we use exactly the same PDE state discretisation
	/// when repricing a security with last event date != schedule last reset date
	/// --> NB : the time discretisation will not be the same...
	if (!itsCalibrationStatus && pdeMethod)
	{ 
		const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (GetModelParams());
		double secVar = modelParamsMF->StateLocalVariance(0.0, ptimeSteps->Elt(ptimeSteps->size()-1),   ptimeSteps->Elt(ptimeSteps->size()-1));
		double calVar = modelParamsMF->StateLocalVariance(0.0, itsDiscountFactorMap.getLastResetTime(), itsDiscountFactorMap.getLastResetTime()); 
		pdeMethod->setStdDevNb(itsCalibrationNbStdev * sqrt(calVar/secVar));
	}
			
	/// Set the basic schedule in the numerical method and...
	numMethod->SetTimeSteps(*ptimeSteps);

	/// this flag changes to false when a new nummethod is set only
	if( itsCalibrationStatus )
		itsCurrentModelIsCalibrationModel = true;

	double firstInductTime = timeInfos.size() != 0 ? timeInfos[0]->GetEventTime() : itsDiscountFactorMap.getLastResetTime();

	/// ...initialise it
	return numMethod->Init(*this,firstInductTime);
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Default implementation is no markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovFunctional::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix(numMethodStates->rows(), numMethodStates->cols(), 0.0 ) );
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
///  Warning : Enhancements are needed!
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::Libor( 
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarkovFunctional::Libor : case payTime != fwdEndTime not handled." );

	return DefaultLibor(curveName, evalTime, fwdStartTime, fwdEndTime, period, fwdResetTime, payTime, states);
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Semi-Closed form formula for caplet/floorlet
///	  using gauss-legender
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::VanillaCaplet(
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
	if( abs(evalTime) < K_DOUBLE_TOL )
	{
		ARM_GP_VectorPtr resultVector( new std::vector<double>( MAX(1,states->size() )));
		
		double numStateMin = itsDiscountFactorMap.getNumMethStateMin(fwdResetTime);
		double numStateMax = itsDiscountFactorMap.getNumMethStateMax(fwdResetTime);
		
		GaussLegendre_Coefficients glc( PRICGLPOINTSNB, numStateMin, numStateMax);
		
		double Strike = strikesPerState[0];

		const ARM_ModelParams* modelParams = GetModelParams();
		const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (modelParams);

		if( !modelParamsMF )	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"MF Model params Expected !");

		double result = 0;
		double resultByState=0;
		double state = 0;
		double OneOnVar = 1./modelParamsMF->StateLocalVariance(0.0, fwdResetTime, fwdResetTime);
		double OneOnStdDev = sqrt( OneOnVar );
		
		std::vector<double> states (PRICGLPOINTSNB);
		for (size_t i(0); i<PRICGLPOINTSNB; i++)
			states[i] = glc.get_point(i);


		ARM_GP_VectorPtr dfRatioStart = itsDiscountFactorMap.DiscountFactorInterpolate(fwdResetTime, fwdStartTime, states, true);
		ARM_GP_VectorPtr dfRatioEnd   = itsDiscountFactorMap.DiscountFactorInterpolate(fwdResetTime, fwdEndTime,   states, true);
		
		/// standard caplet case: pay date = libor end date
		if ( fabs( payTime - fwdEndTime ) < 1 )
		{
			for( i=0 ; i<PRICGLPOINTSNB; ++i )
			{
				state = glc.get_point(i);
				resultByState = capFloor * ( dfRatioStart->Elt(i) - (1. + fwdPeriod*Strike) * dfRatioEnd->Elt(i) );
				result += (resultByState > 0) ? glc.get_weight(i) * resultByState * exp(-0.5*state*state*OneOnVar) : 0;
			}

			result *= OneOnStdDev*INVSQRT2PI;

			state = glc.get_point(0);
			resultByState = capFloor * ( dfRatioStart->Elt(0) - (1. + fwdPeriod*Strike) * dfRatioEnd->Elt(0) );
			result += (resultByState > 0) ? NormalCDF(state*OneOnStdDev) * resultByState : 0;
			
			state = glc.get_point(PRICGLPOINTSNB-1);
			resultByState = capFloor * ( dfRatioStart->Elt(PRICGLPOINTSNB-1) - (1. + fwdPeriod*Strike) * dfRatioEnd->Elt(PRICGLPOINTSNB-1) );
			result += (resultByState > 0) ? (1.-NormalCDF(state*OneOnStdDev)) * resultByState : 0;
		
		}
		/// non standard caplet: pay date != libor end date
		/// works like for standard caplet but takes more time
		else
		{
			ARM_GP_VectorPtr dfRatioPay = itsDiscountFactorMap.DiscountFactorInterpolate(fwdResetTime, payTime, states, true);

			for( i=0 ; i<PRICGLPOINTSNB; ++i )
			{
				state = glc.get_point(i);
				resultByState = capFloor * dfRatioPay->Elt(i) * ( dfRatioStart->Elt(i)/dfRatioEnd->Elt(i) - (1. + fwdPeriod*Strike) );
				result += (resultByState > 0) ? glc.get_weight(i) * resultByState * exp(-0.5*state*state*OneOnVar) : 0;
			}

			result *= OneOnStdDev*INVSQRT2PI;

			state = glc.get_point(0);
			resultByState = capFloor * dfRatioPay->Elt(0) * ( dfRatioStart->Elt(0)/dfRatioEnd->Elt(0) - (1. + fwdPeriod*Strike) );
			result += (resultByState > 0) ? NormalCDF(state*OneOnStdDev) * resultByState : 0;
			
			state = glc.get_point(PRICGLPOINTSNB-1);
			resultByState = capFloor * dfRatioPay->Elt(PRICGLPOINTSNB-1) * ( dfRatioStart->Elt(PRICGLPOINTSNB-1)/dfRatioEnd->Elt(PRICGLPOINTSNB-1) - (1. + fwdPeriod*Strike) );
			result += (resultByState > 0) ? (1.-NormalCDF(state*OneOnStdDev)) * resultByState : 0;
		}

		
		result *= period/fwdPeriod;
		result *= GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getTerminalTime()/K_YEAR_LEN);

		for( std::vector<double>::iterator iter = resultVector->begin() ; iter != resultVector->end() ; ++iter )
			(*iter) = result;

		return resultVector;
	}
	/// eval = reset -> return intrinsic value
	else if ( abs(evalTime-fwdResetTime) < K_DOUBLE_TOL )
	{
		ARM_VectorPtr libor = Libor(curveName, evalTime, fwdStartTime, fwdEndTime, fwdPeriod, fwdResetTime, payTime, states);
		size_t stateSize    = states->size();
		ARM_VectorPtr dfPay = DiscountFactor(curveName, evalTime, payTime, states);
		ARM_GP_VectorPtr result( new std::vector<double>( stateSize ) );
		double payoff, strike = strikesPerState[0];
		
		for (size_t i(0); i<stateSize; i++)
		{
			payoff = capFloor * ( (*libor)[i] - strikesPerState[i] );
			(*result)[i] = (payoff > 0) ? period * payNotional * (*dfPay)[i] * payoff : 0.0;
		}

		return result;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Unimplement VanillaCaplet for evalDate!=0. Please advise." );
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::VanillaDigital(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::VanillaSwaption(
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
	if( abs(evalTime) > K_DOUBLE_TOL )
		ARM_THROW( ERR_INVALID_ARGUMENT, "evalTime != 0" );

	if( fixNotional.size() != fixPayTimes.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "fixNotional.size() != fixPayTimes()" );

	if( fixPayTimes.size() != fixPayPeriods.size() )
		ARM_THROW( ERR_INVALID_ARGUMENT, "( fixPayTimes.size() != fixPayPeriods() )" );


	if (isConstantNotional)
	{
		ARM_GP_VectorPtr resultVector( new std::vector<double>( MAX(1,states->GetNumMethodStates()->size() )));
			
		double numStateMin = itsDiscountFactorMap.getNumMethStateMin(swapResetTime);
		double numStateMax = itsDiscountFactorMap.getNumMethStateMax(swapResetTime);
		
		GaussLegendre_Coefficients glc( PRICGLPOINTSNB, numStateMin, numStateMax);
		
		double Strike	= strikesPerState(0,0);
		double Notional = fixNotional[0];

		const ARM_ModelParams* modelParams = GetModelParams();
		const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (modelParams);

		if( !modelParamsMF )	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"MF Model params Expected !");

		double result = 0;
		double resultByState=0;
		double swapPv = 0;
		double correction;
		double terminalDf = GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getTerminalTime()/K_YEAR_LEN);
		double state = 0;
		double OneOnVar = 1./modelParamsMF->StateLocalVariance(0.0, swapResetTime, swapResetTime);
		double OneOnStdDev = sqrt( OneOnVar );
		double cdf, mktSwapPv;

		/// set this to true to perform rescaling
		bool doRescaling = false;
		
		std::vector<double> states (PRICGLPOINTSNB);
		for (size_t i(0); i<PRICGLPOINTSNB; i++)
			states[i] = glc.get_point(i);


		ARM_GP_VectorPtr startDf = itsDiscountFactorMap.DiscountFactorInterpolate(swapResetTime, floatStartTime, states, true);
		
		size_t N = fixPayTimes.size();
		std::vector<ARM_GP_VectorPtr> payDfs(N);

		double lastStateResult, firstStateResult;
			
		for (size_t j(0); j<N; j++)
			payDfs[j] = itsDiscountFactorMap.DiscountFactorInterpolate(swapResetTime, fixPayTimes[j], states, true);
		

		for( i=0 ; i<PRICGLPOINTSNB; ++i )
		{
			state = glc.get_point(i);
						
			/// float leg 
			resultByState = startDf->Elt(i) - payDfs[N-1]->Elt(i);

			/// fixed leg		
			for (j=0; j<N; j++)
				resultByState -=  fixPayPeriods[j] * Strike * payDfs[j]->Elt(i);
						
			resultByState *= callPut * glc.get_weight(i) * exp(-0.5*state*state*OneOnVar) ;
			
			result += (resultByState > 0) ? resultByState : 0;
			swapPv += resultByState ;

			if (i==0)
			{
				firstStateResult = resultByState;
			}
			if (i==PRICGLPOINTSNB-1)
			{
				lastStateResult  = resultByState;
			}
		}
		
		swapPv *= OneOnStdDev*INVSQRT2PI;
		result *= OneOnStdDev*INVSQRT2PI;
	
		state = glc.get_point(0);
		cdf = NormalCDF(state*OneOnStdDev) ;
		result += (firstStateResult > 0) ? cdf * firstStateResult : 0;
		swapPv += cdf * firstStateResult;

		state = glc.get_point(PRICGLPOINTSNB-1);
		cdf = NormalCDF(state*OneOnStdDev) ;
		result  += (lastStateResult > 0) ? (1.- cdf) * lastStateResult : 0;
		swapPv  += (1.- cdf) * lastStateResult;

		result *= terminalDf;
		swapPv *= terminalDf;

		/// rescale swap pv
		if (doRescaling)
		{
			/// float leg
			mktSwapPv =	GetZeroCurve()->DiscountPrice(floatStartTime/K_YEAR_LEN)
					  -	GetZeroCurve()->DiscountPrice(floatEndTime/K_YEAR_LEN)  ;

			/// fixed leg		
			for (j=0; j<N; j++)
				mktSwapPv -=  fixPayPeriods[j] * Strike * GetZeroCurve()->DiscountPrice(fixPayTimes[j]/K_YEAR_LEN);

			correction = callPut * mktSwapPv / swapPv;	
			result *= correction;
		}

		result *= Notional;
				
		for( std::vector<double>::iterator iter = resultVector->begin() ; iter != resultVector->end() ; ++iter )
			(*iter) = result;

		return resultVector;
		
	}	
	else
	{
		ARM_GP_VectorPtr resultVector( new std::vector<double>( MAX(1,states->GetNumMethodStates()->size() )));
			
		double numStateMin = itsDiscountFactorMap.getNumMethStateMin(swapResetTime);
		double numStateMax = itsDiscountFactorMap.getNumMethStateMax(swapResetTime);
		
		GaussLegendre_Coefficients glc( PRICGLPOINTSNB, numStateMin, numStateMax);
		
		double Strike	= strikesPerState(0,0);
		double Notional = fixNotional[0];

		const ARM_ModelParams* modelParams = GetModelParams();
		const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (modelParams);

		if( !modelParamsMF )	
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"MF Model params Expected !");

		double result = 0;
		double resultByState=0;
		double swapPv = 0;
		double terminalDf = GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getTerminalTime()/K_YEAR_LEN);
		double state = 0;
		double OneOnVar = 1./modelParamsMF->StateLocalVariance(0.0, swapResetTime, swapResetTime);
		double OneOnStdDev = sqrt( OneOnVar );
		double cdf;

		/// set this to true to perform rescaling
		bool doRescaling = false;
		
		std::vector<double> states (PRICGLPOINTSNB);
		for (size_t i(0); i<PRICGLPOINTSNB; i++)
			states[i] = glc.get_point(i);

		double lastStateResult, firstStateResult;

		size_t N = fixPayTimes.size();
		std::vector<ARM_GP_VectorPtr> payDfs(N);
			
		for (size_t j(0); j<N; j++)
			payDfs[j] = itsDiscountFactorMap.DiscountFactorInterpolate(swapResetTime, fixPayTimes[j], states, true);
		
		size_t M = floatStartTimes.size();
		std::vector<ARM_GP_VectorPtr> startDfs(M);
		std::vector<ARM_GP_VectorPtr> endDfs(M);

		for (j=0; j<M; j++)
		{
			startDfs[j] = itsDiscountFactorMap.DiscountFactorInterpolate(swapResetTime, floatStartTimes[j], states, true);
			endDfs[j]	= itsDiscountFactorMap.DiscountFactorInterpolate(swapResetTime, floatEndTimes[j], states, true);
		}

		for( i=0 ; i<PRICGLPOINTSNB; ++i )
		{
			state = glc.get_point(i);
				
			resultByState = 0.;
			/// float leg
			for (j=0; j<M; j++)
				resultByState +=  floatNotional[j] * (startDfs[j]->Elt(i)-endDfs[j]->Elt(i));
			
			/// fixed leg		
			for (j=0; j<N; j++)
				resultByState -=  fixNotional[j] *   fixPayPeriods[j] * Strike * payDfs[j]->Elt(i);
						
			resultByState *= callPut * glc.get_weight(i) * exp(-0.5*state*state*OneOnVar) ;
			
			result += (resultByState > 0) ? resultByState : 0;
			swapPv += resultByState ;

			if (i==0)
			{
				firstStateResult = resultByState;
			}
			if (i==PRICGLPOINTSNB-1)
			{
				lastStateResult  = resultByState;
			}
		}
		
		swapPv *= OneOnStdDev*INVSQRT2PI;
		result *= OneOnStdDev*INVSQRT2PI;
	
		state = glc.get_point(0);
		cdf = NormalCDF(state*OneOnStdDev) ;
		result += (firstStateResult > 0) ? cdf * firstStateResult : 0;
		swapPv += cdf * firstStateResult;

		state = glc.get_point(PRICGLPOINTSNB-1);
		cdf = NormalCDF(state*OneOnStdDev) ;
		result  += (lastStateResult > 0) ? (1.- cdf) * lastStateResult : 0;
		swapPv  += (1.- cdf) * lastStateResult;

		result *= terminalDf;
		swapPv *= terminalDf;

		for( std::vector<double>::iterator iter = resultVector->begin() ; iter != resultVector->end() ; ++iter )
			(*iter) = result;

		return resultVector;
	}
	
	

	return ARM_GP_VectorPtr(NULL);

}



////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
std::vector<double>& ARM_MarkovFunctional::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new std::vector<double>(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MarkovFunctional::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MarkovFunctional::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex > GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_GP_MatrixPtr& numMethodStates = states->GetNumMethodStates();
	states->SetModelStates(numMethodStates);
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_MarkovFunctional::VanillaSpreadOptionLet(const string& curveName,
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
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSpreadOption : unimplemented function for ARM_MarkovFunctional Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_MarkovFunctional::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	if (linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_MarkovFunctional::VolatilitiesAndCorrelations : linearVol not implemented" );

	/// Convention : 
	/// vols[i] = instant. vol of process between timeSteps[i] and timeSteps[i+1]
	///
	ARM_CurveModelParam& mp = (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility);
	size_t timeStepsSize = timeSteps.size() ; 
	std::vector<double> volsVec(timeStepsSize);
	double epsToBeSafeWithInterpol = 1.0e-5;

	for (size_t j=0; j<timeStepsSize; j++)
	{
		volsVec[j] = mp.GetValue(timeSteps[j]+epsToBeSafeWithInterpol);
	}

	vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,volsVec.size(),&volsVec[0]) );
	d1Vols	= ARM_GP_MatrixPtr( NULL );
	correls	= ARM_GP_MatrixPtr( NULL );
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: ValidateModelParams
///	Returns : void
///	Action  : ValidateModelParams
////////////////////////////////////////////////////
bool ARM_MarkovFunctional::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*>(&params);

	if( !modelParamsMF )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsMF" );

    if(!params.DoesModelParamExist(ARM_ModelParamType::Volatility) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsMF" );

	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: PreProcessing
///	Returns : void
///	Action  : PreProcessing
////////////////////////////////////////////////////
void ARM_MarkovFunctional::PreProcessing(ARM_ModelFitter& modelFitter)
{
/// Not Implemented Yet
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: PostProcessing
///	Returns : void
///	Action  : PostProcessing
////////////////////////////////////////////////////
void ARM_MarkovFunctional::PostProcessing(const ARM_ModelFitter& modelFitter)
{

}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: ValidateCalibMethod
///	Returns : void
///	Action  : ValidateCalibMethod
////////////////////////////////////////////////////
void ARM_MarkovFunctional::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : AdviseBreakPointTimes
////////////////////////////////////////////////////
void ARM_MarkovFunctional::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, ARM_ModelParam* inputModelParam, size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    std::vector<double>  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
    case ARM_ModelParamType::Volatility:
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
					/// ignore this instrument
					portfolio->SetWeight(0.0,i);
				}
            }
			modelParam->UpdateValues(&tmpdates);
        }
        break;
    case ARM_ModelParamType::MeanReversion:
        {
            double date = portfolio->GetAsset(0)->GetFlowEndDates()->Elt(0) - asOfDate; 
            tmpdates.push_back(date);
            for(i=1; i<size1; i++)
            {
                double startlag = portfolio->GetAsset(i)->GetFlowEndDates()->Elt(0) - asOfDate;
                if(fabs (date - startlag) > FRMVOL_LAG_THRESHOLD)
                {
                    tmpdates.push_back(startlag);
                    date = startlag;
                }
				else
				{
					/// ignore this instrument
					portfolio->SetWeight(0.0,i);
				}
            }
			modelParam->UpdateValues(&tmpdates);
        }
    default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... an HW1F model only supports mean reversion and volatility" );
    }
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: Induct
///	Returns : ARM_PricingStatesPtr
///	Action  : induct...
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MarkovFunctional::Induct(ARM_PricingStatesPtr& states,double toTime)
{
	/// Keep memory of payoffs in the boostrap case
	if (itsCalibrationStatus && itsNumericalModelFitter->BoostrapMode())
	{
		/// Could be optimized avoiding recopies
		if (itsPrevToTime != toTime)
			itsPrevPayoffs = states->GetPayoffs(); /// save payoffs
		else
			states->SetPayoffs(itsPrevPayoffs);    /// set payoffs if Induct on same dates
	}
	
	/// Backward PDE induction
	ARM_PricingStatesPtr newStates( ARM_PricingModel::Induct( states, toTime ) );

	/// Extract 1/Numeraire
	if( itsCalibrationStatus )
		BackwardCalibration( newStates, toTime );

	/// update itsPrevToTime
	itsPrevToTime = toTime;

	/// v'la l'travail
	return newStates;
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: FirstPricingStates
///	Returns : void
///	Action  : FirstPricingStates
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MarkovFunctional::FirstPricingStates( size_t bucketSize ) const
{
	ARM_PricingStatesPtr newPricingStates;

	if( itsCalibrationStatus )
	{
		newPricingStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,1,itsDiscountFactorMap.getResetTimesSize(),1) );
		newPricingStates->resizePayoffs(0);
	}
	else
		newPricingStates = ARM_PricingStatesPtr( new ARM_PricingStates(bucketSize,1,0,1) );
	
	return newPricingStates;
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: DiscountFactor
///	Returns : void
///	Action  : DiscountFactor
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::DiscountFactor( const string& curveName, double evalTime, double maturityTime, const ARM_PricingStatesPtr& states) const
{
	if( evalTime < K_DOUBLE_TOL )
	{
		if( states == ARM_PricingStatesPtr(NULL) || states->size() == 0 ) 
			return static_cast<ARM_VectorPtr>(new std::vector<double>(1,GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN)));
		else
			return static_cast<ARM_VectorPtr>(new std::vector<double>(states->size(), GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN) ));
	}

	if( abs( evalTime - maturityTime ) < K_DOUBLE_TOL )
		return static_cast<ARM_VectorPtr>(new std::vector<double>(states->size(), 1. ));

	ARM_GP_VectorPtr result(0);

	result = itsDiscountFactorMap.DiscountFactor( evalTime, maturityTime, states );

	ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();

	size_t resetIdx = itsDiscountFactorMap.NDaysFromResetTimeIdx( evalTime, 7. );
	double correction = ZcCurve->DiscountPrice(itsDiscountFactorMap.getStartTime(resetIdx)/K_YEAR_LEN) / ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);

	for (std::vector<double>::iterator iter = result->begin() ; iter != result->end() ; ++iter)
		(*iter) *= correction;

	return result;
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: NumMethodStateLocalVariances
///	Returns : void
///	Action  : NumMethodStateLocalVariances
////////////////////////////////////////////////////
void ARM_MarkovFunctional::NumMethodStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelRank= GetModelRank();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelRank;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelRank+1));
	size_t i;

	// All the variance is in the numerical method

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsMF* const) GetModelParams())->StateLocalVariance(step,nextStep,nextStep));
		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: NumMethodStateGlobalVariances
///	Returns : void
///	Action  : NumMethodStateGlobalVariances
////////////////////////////////////////////////////
void ARM_MarkovFunctional::NumMethodStateGlobalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& globalVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( globalVariances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	globalVariances.resize(nbSteps*(modelNb+1));

	/// fills the variance
    globalVariances[offsetIndex2+0]=new ARM_GP_TriangularMatrix(1,0.0);

    for(size_t i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
        
		/// [i+1] => variance from 0 -> ti+1
        /// we can't sum up local variance !
        globalVariances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsMF* const) GetModelParams())->StateLocalVariance(0.0,nextStep,nextStep));
        step=nextStep;
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: EulerLocalDrifts
///	Returns : void
///	Action  : EulerLocalDrifts
////////////////////////////////////////////////////
void ARM_MarkovFunctional::EulerLocalDrifts(const std::vector<double>& timeSteps,ARM_GP_MatrixPtr& relativeDrifts,ARM_GP_MatrixPtr& absoluteDrifts) const
{
	/// FIX FIX : only valid for cst MRS !!!
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1,-GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0) ) );
    for(size_t i=0;i<timeSteps.size()-1;++i)
        (*relativeDrifts)(i,0) *= (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;

	absoluteDrifts= ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0.0 ) );
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: ModelStateLocalVariances
///	Returns : void
///	Action  : ModelStateLocalVariancess
////////////////////////////////////////////////////
void ARM_MarkovFunctional::ModelStateLocalVariances( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localVariances ) const
{
	NumMethodStateLocalVariances( timeSteps, localVariances );
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: ModelStateLocalCorrels
///	Returns : void
///	Action  : ModelStateLocalCorrels
////////////////////////////////////////////////////
void ARM_MarkovFunctional::ModelStateLocalCorrels( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& ma )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"ModelStateLocalCorrels : unimplemented function for ARM_MarkovFunctional Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: VarianceToTime
///	Returns : void
///	Action  : VarianceToTime
////////////////////////////////////////////////////
double ARM_MarkovFunctional::VarianceToTime(double var,double minTime,double maxTime) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VarianceToTime : unimplemented function for ARM_MarkovFunctional Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: toString
///	Returns : void
///	Action  : toString
////////////////////////////////////////////////////
string ARM_MarkovFunctional::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Markov Functional Model 1F\n";
    os << indent << "----------------------------\n";
	os << ARM_PricingModel::toString(indent);
	//os << itsDiscountFactorMap.toString(indent, nextIndent);
	os << itsDiscountFactorMap.toString(nextIndent, nextIndent);
	return os.str();	
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: BackwardCalibration
///	Returns : void
///	Action  : BackwardCalibration
////////////////////////////////////////////////////
void ARM_MarkovFunctional::BackwardCalibration( const ARM_PricingStatesPtr& states, double toTime )
{
	if( itsDiscountFactorMap.DoesResetTimeExist(toTime) ) 
	{
		/// toTime is a ResetDate -> Computation and storage
		states->resizePayoffs( itsDiscountFactorMap.dfsNbAtTime(toTime) );
		BuildDFs( states, toTime, itsDiscountFactorMap.getResetTimeIdx(toTime) );
		itsDiscountFactorMap.StoreDiscountFactors( toTime, states );
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : void
///	Action  : VolatilitiesAndCorrelationTimesSteps for PDE
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MarkovFunctional::VolatilitiesAndCorrelationTimesSteps() const
{
	const ARM_ModelParams* modelParams = GetModelParams();
	const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (modelParams);

	if( !modelParamsMF )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"MF Model params Expected !");

	return modelParamsMF->ModelParamsTimeSteps();

}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: BuildDFs
///	Returns : void
///	Action  : Build discount factors
////////////////////////////////////////////////////
void ARM_MarkovFunctional::BuildDFs( const ARM_PricingStatesPtr& states, double toTime, size_t lineIdx )
{
	/// This vector will contain probabilities
	ARM_GP_VectorPtr probas( new std::vector<double>( states->GetNumMethodStates()->size() ));

	/// Fills vec with P^Q_i+1 ( 1_X<x )		
	buildNormalProbabilities( probas, states, toTime );

	/// Computes new discount factors and fills pricingstates with them
	itsNumericalModelFitter->UpdateStates( states, probas, lineIdx );

	/// Correct probabilites so that E^N+1 [1/B(T_j,T_N+1)] = B(0,T_j)/B(0,T_N+1). Multiplicative correction
	CorrectProbaChanges( states, toTime, lineIdx );
} 

////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: CorrectProbaChanges
///	Returns : void
///	Action  : Build discount factors
////////////////////////////////////////////////////
void ARM_MarkovFunctional::CorrectProbaChanges( const ARM_PricingStatesPtr& states, double toTime, size_t lineIdx, bool correctNumeraireOnly ) const
{
	/// Target
	size_t lastIdx = itsDiscountFactorMap.getResetTimesSize()-1;
	double lastDf  = GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getEndTime(lastIdx)/K_YEAR_LEN) ;
	
	const ARM_GP_MatrixPtr nummethodStates = states->GetNumMethodStates();
	size_t i,N = nummethodStates->cols();

	/// 1/StdDev, 1/Var, etc
	const ARM_ModelParams* modelParams = GetModelParams();
	const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (modelParams);

	if( !modelParamsMF )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"MF Model params Expected !");

	double OneOnStdDev = 1./sqrt(modelParamsMF->StateLocalVariance(0.0, toTime, toTime));
	double OneOnVar = OneOnStdDev*OneOnStdDev;
		
	/// What we want to compute
	size_t j=0;
	int payoffIndex = itsDiscountFactorMap.dfsNbAtTime(toTime)-1;

	int bound = correctNumeraireOnly ? payoffIndex : 0;

	while (payoffIndex>=bound)
	{
		double correction = GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getStartTime(lineIdx)/K_YEAR_LEN) / lastDf;
		double point = 0.;
		double byComputation = 0.;

		for( i = 1 ; i < N ; ++i )
		{
			GaussLegendre_Coefficients glc( CALIBGLPOINTSNB, nummethodStates->Elt(0,i-1), nummethodStates->Elt(0,i) );

			for( j=0 ; j<CALIBGLPOINTSNB ; ++j )
			{
				point = glc.get_point(j);
				byComputation += glc.get_weight(j)*exp(-0.5*point*point*OneOnVar)*itsDiscountFactorMap.DFRatioInterpolateFromStates( states, point, i, payoffIndex);
			}

		}

		byComputation *= OneOnStdDev*INVSQRT2PI;

		byComputation += NormalCDF( nummethodStates->Elt(0,0)*OneOnStdDev ) * states->GetPayoff(0,payoffIndex);
		byComputation += (1.00-NormalCDF(nummethodStates->Elt(0,i-1)*OneOnStdDev))*states->GetPayoff(i-1,payoffIndex);
		
		correction /= byComputation;

		/// states correction
		ARM_MemPool_Matrix::iterator statesIterator = states->payoffsBeginIterator(payoffIndex);

		for( i=0 ; i<N ; ++i, ++statesIterator)
			(*statesIterator) *=correction;

		payoffIndex --;
		lineIdx ++;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: buildNormalProbabilities
///	Returns : void
///	Action  : Build discount factors
////////////////////////////////////////////////////
void ARM_MarkovFunctional::buildNormalProbabilities( ARM_GP_VectorPtr& probas, const ARM_PricingStatesPtr& states, double toResetTime ) const
{
	const ARM_ModelParamsMF* modelParamsMF = dynamic_cast<const ARM_ModelParamsMF*> (GetModelParams());

	if( !modelParamsMF )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"MF Model params Expected !");

	/// General stuff used to build probabilities
	const ARM_GP_MatrixPtr nummethodStates = states->GetNumMethodStates();
	double OneOnStdDev = 1./sqrt(modelParamsMF->StateLocalVariance(0.0, toResetTime, toResetTime));
	double OneOnVar = OneOnStdDev*OneOnStdDev;
	size_t lastIdx  = itsDiscountFactorMap.getResetTimesSize()-1;
	size_t resetIdx = itsDiscountFactorMap.getResetTimeIdx(toResetTime);
	size_t N = probas->size();

	if( lastIdx == resetIdx)
	{
		/// Probabilities for the terminal date. 
		for( size_t i = 0 ; i < N ; ++i )
			probas->Elt(i) = NormalCDF( nummethodStates->Elt(0,i)*OneOnStdDev );
	}
	else
	{		
		/// Global remarks
		/// Here we compute P^Q^i+1 [ 1_{X_i<x} ] = E^{Q_{N+1}} [\frac{B(T_i,T_i+1)*B(0,T_N+1)}{B(0,T_i+1)*B(T_i,T_N+1)} 1_{{X_i} < x}]
		/// 
		/// sigma/sqrt(2pi)
		
		double lastDf		  = GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getEndTime(lastIdx)/K_YEAR_LEN);
		double lastDfOnLevel0 = 0.0;
		const ARM_IntVector&  absIdx = itsNumericalModelFitter->getCalibSecPayDatesAbsIndexes(resetIdx);
		const ARM_IntVector&  relIdx = itsNumericalModelFitter->getCalibSecPayDatesRelIndexes(resetIdx);
		const std::vector<double>&  interestTerms = itsNumericalModelFitter->getCalibSecInterestTerms(resetIdx);

		for (size_t k(0); k<interestTerms.size(); k++)
		{
			lastDfOnLevel0 +=	interestTerms[k]
							*	GetZeroCurve()->DiscountPrice(itsDiscountFactorMap.getEndTime(absIdx[k])/K_YEAR_LEN);

		}

		lastDfOnLevel0 = lastDf / lastDfOnLevel0;
		
		/// mult factor for integrand
		double MultiplicationFactor  = OneOnStdDev * INVSQRT2PI;
		MultiplicationFactor		*= lastDfOnLevel0;

		/// Used for computation
		size_t			 j = 0;
		size_t payoffIndex = itsDiscountFactorMap.dfsNbAtTime(toResetTime) - 2;
		
		double IntegratedStuff	= 0.;
		double point			= 0.;

		/// First probability
				
		double levelOnLastDf = itsDiscountFactorMap.LevelOnLastDfInterpolateFromStates (
																	states, 
																	nummethodStates->Elt(0,0), 
																	0, 
																	payoffIndex, 
																	relIdx, 
																	interestTerms);

		probas->Elt(0) = NormalCDF( nummethodStates->Elt(0,0)*OneOnStdDev ) * levelOnLastDf * lastDfOnLevel0;
		probas->Elt(0) = MAX(1e-15, probas->Elt(0) );
		
		/// Filling the vector
		for( size_t i = 1 ; i < probas->size() ; ++i )
		{
			GaussLegendre_Coefficients glc( CALIBGLPOINTSNB, nummethodStates->Elt(0,i-1), nummethodStates->Elt(0,i) );

			/// IntegratedStuff will contain \int_{x_i-1}^{x_i} exp(-x^2/2sigma^2)*B(T_lineIdx,T_lineIdx+1)/B(T_lineIdx, T_numeraireDate) dx * 1/(sigma*sqrt(2Pi) )* SamedFratio(0)
			IntegratedStuff = 0.0;
			
			for( j=0 ; j<CALIBGLPOINTSNB ; ++j )
			{
				point = glc.get_point(j);
				
				levelOnLastDf = itsDiscountFactorMap.LevelOnLastDfInterpolateFromStates (
														states, 
														point, 
														i, 
														payoffIndex, 
														relIdx, 
														interestTerms);


				IntegratedStuff +=		glc.get_weight(j) 
									*	exp(-0.5*point*point*OneOnVar)
									*	levelOnLastDf;
			}

			/// add small integral to get \int_{-infty}^{x_i} ... dx
			probas->Elt(i) = probas->Elt(i-1) + IntegratedStuff * MultiplicationFactor;
		}

						
		/// normalize
		levelOnLastDf = itsDiscountFactorMap.LevelOnLastDfInterpolateFromStates (
														states, 
														nummethodStates->Elt(0,i-1), 
														i-1, 
														payoffIndex, 
														relIdx, 
														interestTerms);

		MultiplicationFactor =  1./( probas->Elt( probas->size()-1 ) + (1.00-NormalCDF(nummethodStates->Elt(0,i-1)*OneOnStdDev))*levelOnLastDf*lastDfOnLevel0);

		
		std::vector<double>::iterator iter;

		/// Correction to be sure that \int_{-infty}^{+infty} exp(-x^2/2)*lastDfOnLevel(x,T_i,T_i+1,T_N+1)/dFratio(0)/{sigma*sqrt(2pi)} dx = 1
		for(iter = probas->begin(); iter != probas->end() ; ++iter )
			(*iter) *= MultiplicationFactor;

		/// Maj/Min to avoid numerical problems
		/// double normalMin = NormalCDF(-7.), normalMax = NormalCDF(7.);
		double	normalMin = 1.e-7, 
				normalMax = 1.-1.e-7;

		for( iter = probas->begin(); iter != probas->end() ; ++iter )
			(*iter) = MAX( normalMin, MIN(normalMax,(*iter)) );
		
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_MarkovFunctional::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 1, 0.0 ) );
	absoluteDrifts	= ARM_GP_MatrixPtr(NULL);

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		(*relativeDrifts)(i,0) = ((const ARM_ModelParamsMF* const) GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: SetNumMethod
///	Returns: void
///	Action : sets nummethodtomodel, and changes the 
///  currentModelIsCalibratingModel flag
////////////////////////////////////////////////////
void ARM_MarkovFunctional::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	ARM_PricingModel::SetNumMethod( numMethodPtr );
	itsCurrentModelIsCalibrationModel = false;
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovFunctional
///	Routines: AdviseCurrentCalib
///	Returns : void
///	Action  : Advises the model that of the calibration  ... 
///           The model advises just the model params
////////////////////////////////////////////////////
void ARM_MarkovFunctional::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
	GetModelParams()->PostProcessing(modelFitter,this);
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovFunctional
///	Routine: ComputeNumeraireTimes
///	Returns: 
///	Action : HK --> num method is always under terminal measuer
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovFunctional::ComputeNumeraireTimes( const ARM_TimeInfoPtrVector& timeInfos ) const
{
	return static_cast<ARM_VectorPtr>(new std::vector<double>(1,itsDiscountFactorMap.getTerminalTime()));
}


CC_END_NAMESPACE()
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

