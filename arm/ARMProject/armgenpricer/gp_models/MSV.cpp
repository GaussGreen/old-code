/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV.cpp
 *
 *  \brief file for the Markovian Stochastic volatility model
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */



/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/MSV.h"
#include "gpmodels/ModelParamsMSV.h"
#include "gpmodels/riccatimsv.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"


#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/riccati.h"


/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/rungekutta.h"


CC_BEGIN_NAMESPACE( ARM )

#define FIRST_STATE_VARIABLE    0
const double VOL_LIMIT			  = 0.000001;
const double GL_NBPOINTS_FIRST    = 4; // Number of points per Year
const double GL_NBPOINTS_SECOND	  = 4; // Number of points per Year
const double ODE_PRECISION		  = 1e-6; 


const double ARM_MarkovSV::VOL_NORMAL_MAX = 0.15;

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarkovSV::ARM_MarkovSV( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params )
:	ARM_SVModels(zc,params)

{
	int sizeSchedule = ((ARM_ModelParamsMSV*) GetModelParams())->GetSchedule().size();
	itsStoredData = ARM_GeneralSwaptionData(sizeSchedule); 

}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarkovSV::ARM_MarkovSV(const ARM_MarkovSV& rhs)
: ARM_SVModels(rhs)
{
    // Copy class attributes if any
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_MarkovSV::~ARM_MarkovSV()
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_MarkovSV& ARM_MarkovSV::operator=(const ARM_MarkovSV& rhs)
{
	if(this != &rhs)
	{
		ARM_SVModels::operator=(rhs);
		// Copy class attributes if any
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : void
///	Action  : VolatilitiesAndCorrelationTimesSteps for PDE
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_MarkovSV::VolatilitiesAndCorrelationTimesSteps() const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VolatilitiesAndCorrelationTimesSteps : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Default implementation is no markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_MarkovSV::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix(numMethodStates->rows(), numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MarkovSV::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
	
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in MSV1F model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the MSV1F model!");

		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

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
///	Class   : ARM_MarkovSV
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
   	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"Libor : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV::VanillaCaplet(
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
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaCaplet : unimplemented function for ARM_MarkovSV Model!");
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV::VanillaDigital(
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
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaDigital : unimplemented function for ARM_MarkovSV Model!");
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV::VanillaSwaption(
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
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VanillaSwaption : unimplemented function for ARM_MarkovSV Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
std::vector<double>& ARM_MarkovSV::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new std::vector<double>(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MarkovSV::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"MCModelStatesFromToNextTime : unimplemented function for ARM_MarkovSV Model!");
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MarkovSV::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
{   
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_MarkovSV Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_MarkovSV::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_MarkovSV Model!");
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

