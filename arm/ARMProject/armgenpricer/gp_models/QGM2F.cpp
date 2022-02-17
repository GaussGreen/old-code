/*!
 *
 * Copyright (c) CDC IXIS CM May 2006 Paris
 *
 *	\file QGM2F.cpp
 *
 *  \brief 
 *
 *	\author  Y KHLIF
 *	\version 1.0
 *	\date May 2006
 */



/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/QGM2F.h"
#include "gpmodels/ModelParamsQGM2F.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"

/// gpcalib
#include "gpcalib/modelfitter.h"
#include "gpcalib/calibmethod.h"

/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <inst/portfolio.h>
#include <inst/swaption.h>

/// nag
#include "gpbase/removenagwarning.h"
#include "nag.h"
#include "nags.h"

#include <iomanip> /// for setprecision()
#include <algorithm>
CC_USING_NS(std,pair)
#include <set>
CC_USING_NS(std,set)

CC_BEGIN_NAMESPACE( ARM )

const short BCOEF_STORED    = 1;
const short B_STORED        = 2;
const short CACOEF_STORED   = 4;
const short VAR_STORED      = 8;

const double VOL_LIMIT      = 0.000001;
const double MRS_LIMIT      = 1.0e-6;
const double SKEW_LIMIT      = 1.0e-6;
const double DELTA_LIMIT    = 1.0e-10;

const double STDDEV_RATIO           = 5.0;
const double MIN_SLOPE_STEP_ZERO	= 0.001;
const double MIN_PRICE_ZERO	        = 0.0000001;
const double RATIO_NR_ROOT		    = 0.000001;
const double RATIO_NR_MAX_MOVE		= 0.01;
const double STEP_NR_1D		        = 0.000001;
const double MAX_NR_ITER		    = 20;

const int NBPOINT_YEAR_GL   = 2;
const int NBPOINT_MIN_GL    = 2;



////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_QGM2F::CopyNoCleanUp(const ARM_QGM2F& rhs)
{
    itsTime             = rhs.itsTime;
    itsMaturity         = rhs.itsMaturity;
    itsA                = rhs.itsA;
    itsB                = rhs.itsB;
    itsC                = rhs.itsC;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QGM2F::ARM_QGM2F(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQGM2F& params) :
ARM_PricingModelIR(zc,&params)
{
    CC_ARM_SETNAME(ARM_QGM2F_MODEL);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QGM2F::ARM_QGM2F(const ARM_QGM2F& rhs)
: ARM_PricingModelIR(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_QGM2F::~ARM_QGM2F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QGM2F& ARM_QGM2F::operator=(const ARM_QGM2F& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_QGM2F::Clone() const
{
	return new ARM_QGM2F(*this);
}


double ARM_QGM2F::A(double t, double T) const
{
	return 0.0;
}


double ARM_QGM2F::B(double t, double T) const
{
    return 0.0;
}


double ARM_QGM2F::C(double t, double T) const
{
	return 0.0;
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: 
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM2F::InitABC(const std::vector<double>& times, const vector< std::vector<double> >& maturities)
{
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QGM2F::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in QGM model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the QGM model!");

        if(numeraire->GetType() != ARM_Numeraire::TerminalZc )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            ": only TerminalZc numeraire supported by QGM model at the moment!");
		
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
        int nbDir = GetModelParams()->FactorCount();
        ARM_PricingStatesPtr initStates(new ARM_PricingStates(1,nbDir,0));
        for(int i=0;i<nbDir;++i)
            initStates->SetModelState(0,i,0.0);
		
        return initStates;
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Pre-Processing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM2F::PreProcessing(ARM_ModelFitter& modelFitter)
{
///After validate modelfitter, we call this function to manage correctly optimisation.
    GetModelParams()->PreProcessing(modelFitter,modelFitter.GetFactorNb());
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: Post Processing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM2F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: AdviseCurrentCalib
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM2F::AdviseCurrentCalib(ARM_ModelFitter& modelFitter)
{
   
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: AdviseCurrentCalibSecIndex
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_QGM2F::AdviseCurrentCalibSecIndex(size_t index,ARM_ModelFitter& modelFitter)
{}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: AdviseBreakPointTimes
///	Returns: void
///	Action : sets the corresponding suggested break point times to the model param
////////////////////////////////////////////////////
void ARM_QGM2F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio,
				  ARM_ModelParam* inputModelParam,
				  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int pfSize = portfolio->GetSize();
    std::vector<double> sched;
    double prevTime = 0.0,endTime,payTime,curTime;
    ARM_Vector *endDates,*payDates,*resetDates;

    ARM_Security *sec;
    ARM_CapFloor *cap;
    ARM_Swaption *swaption;

    switch( modelParam->GetType() )
    {
	case ARM_ModelParamType::Volatility:
        {
            for(int i=0; i<pfSize; ++i)
            {
                /// Get the last expiry of the product
                sec = portfolio->GetAsset(i);
                if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
                {
					resetDates = cap->GetResetDates();
                    curTime = (*resetDates)[resetDates->size()-1] - asOfDate;
                }
				
                else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
                    curTime = swaption->GetExpiryDate().GetJulian() - asOfDate;
				
                else
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                    ": QGM2F may be calibrated only on cap/floor or swaption");
				
                if(curTime > prevTime)
                {
                    sched.push_back(curTime);
                    prevTime = curTime;
                }
            }
			modelParam->UpdateValues(&sched);
        }
        break;
		
	case ARM_ModelParamType::Skew:
		{
            for(int i=0; i<pfSize; ++i)
            {
                /// Get the very last Zc maturity of the product
                sec = portfolio->GetAsset(i);
                if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
                {
					endDates = cap->GetSwapLeg()->GetFwdRateEndDates();
                    endTime  = (*endDates)[endDates->size()-1] - asOfDate;
					payDates = cap->GetPaymentDates();
					payTime  = (*payDates)[payDates->size()-1] - asOfDate;
					curTime     = (endTime < payTime ? payTime : endTime);
                }
                else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
                {
					endTime  = swaption->GetFloatLeg()->GetEndDateNA().GetJulian() - asOfDate;
					payDates = swaption->GetFixedLeg()->GetFlowEndDates();
					//get the next payment date to avoid taking into account days adjustments in swaption
					if(((*payDates)[0]-swaption->GetExpiryDate().GetJulian())<15) 
						payTime  = (*payDates)[1] - asOfDate;
					else payTime  = (*payDates)[0] - asOfDate;
					
					curTime     = payTime;
                }
                else
                    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
                    ": QGM2F may be calibrated only on cap/floor or swaption");
				
                if(curTime > prevTime)
                {
                    sched.push_back(curTime);
                    prevTime = curTime;
                }
            }
			modelParam->UpdateValues(&sched);
			
        }
        break;
		
	case ARM_ModelParamType::MeanReversion:
        {
			for(int i=0; i<pfSize; ++i)
			{
				/// Get the very last Zc maturity of the product
				sec = portfolio->GetAsset(i);
				if((cap = dynamic_cast< ARM_CapFloor* >(sec)) != NULL)
				{
					endDates = cap->GetSwapLeg()->GetFwdRateEndDates();
					endTime  = (*endDates)[endDates->size()-1] - asOfDate;
					payDates = cap->GetPaymentDates();
				}
				else if((swaption = dynamic_cast< ARM_Swaption* >(sec)) != NULL)
				{
					endTime  = swaption->GetFloatLeg()->GetEndDateNA().GetJulian() - asOfDate;
					payDates = swaption->GetFixedLeg()->GetFlowEndDates();
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					": QGM2F may be calibrated only on cap/floor or swaption");
				
				payTime  = (*payDates)[payDates->size()-1] - asOfDate;
				
				curTime     = (endTime < payTime ? payTime : endTime);
				if(curTime > prevTime)
				{
					sched.push_back(curTime);
					prevTime = curTime;
				}
				modelParam->UpdateValues(&sched);
			}
        }
        break;
        
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
            "Unknown type : QGM2F model only supports mean reversion, volatility and skew" );
    }

}



////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_QGM2F::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM2F::ModelStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F //Tested
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local  variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM2F::NumMethodStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	
}

////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_QGM2F::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
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
///	Class  : ARM_QGM2F
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QGM2F::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances) const
{	
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM2F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");

	int i,nbStates=states->size();
     ARM_VectorPtr values(new std::vector<double>(nbStates));
    for(i=0;i<nbStates;++i)
        (*values)[i]=0.0;

    return values;
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM2F::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");
	
	/// handle libor with no convexity correction
	if (	(	evalTime <= K_NEW_DOUBLE_TOL 
	 	||	states   == ARM_PricingStatesPtr(NULL) ) 
		&&  (-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5) )
    {
		ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		double libor = ((*ZcStart)[0]/(*ZcEnd)[0]-1.0)/period;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,libor) );
    }

    int i,nbStates=states->size();
    ARM_VectorPtr values( new std::vector<double>(nbStates) );
    if(-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5)
    {
		ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
		ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
        /// No convexity
        for(i=0;i<nbStates;++i)
            (*values)[i]=((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
    }
    return values;
}


///////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_QGM2F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
// FIXMEFRED: mig.vc8 (25/05/2007 15:28:30):cast
	return static_cast<ARM_PricingStatesPtr>(new ARM_PricingStates(bucketSize,1,0,1));
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
std::vector<double>& ARM_QGM2F::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new std::vector<double>(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void 
///	Action : 
////////////////////////////////////////////////////

void ARM_QGM2F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
	const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	size_t factorsNb= FactorCount();
	size_t statesNb = states->size();
	double currentState,stdDev;
	size_t modelNb	= GetModelNb();
	
	for( size_t i=0;i<statesNb; ++i )
	{
		for( size_t j=0;  j<factorsNb; ++j )
		{
			currentState = 0.0;
			for( size_t k =0; k<=j; ++k )
			{
				stdDev			= GetLocalMatrixElemWithModelNb(localStdDev,timeIndex,modelNb,j,k);
				double gaussian = states->GetNumMethodState(i,modelNb+k);
				currentState  += stdDev*gaussian;
			}
			states->SetModelState(i,j+modelNb,currentState);
		}
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM2F::VanillaCaplet(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");
	/// Handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:28:40):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

    int nbStates=states->size();

#if defined(__GP_STRICT_VALIDATION)
    if(nbStates != strikesPerState.size())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": inconsistency between model & strikes states");
#endif
	
	ARM_VectorPtr values( new std::vector<double>(nbStates) );

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM2F::VanillaDigital(
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
	ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:28:49):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

    int nbStates=states->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates));


    return values;
}



////////////////////////////////////////////////////
///	Class  : ARM_QGM2F
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QGM2F::VanillaSwaption(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": not implemented yet");

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:29:04):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

    int nbStates=states->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates));

    return values;
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routine : toString
///	Returns : string
///	Action  : object dump into a string
////////////////////////////////////////////////////
string ARM_QGM2F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << indent << "QGM2F Model\n";
    os << indent << "-----------\n";

    os << ARM_PricingModel::toString(indent);

    os << indent << "Function Datas\n";
    os << indent << setw(5)     << "time";
    os << indent << setw(10)    << "mrs";
    os << indent << setw(10)    << "2Vol2";
    os << indent << setw(10)    << "Delta";
    os << indent << setw(10)    << "NbR";
    os << indent << setw(10)    << "RInf";
    os << indent << setw(10)    << "Rsup";
    os << endl;
    

    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_QGM2F::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsQGM2F* modelParamsQGM2F = dynamic_cast<const ARM_ModelParamsQGM2F*>(&params);
	if( !modelParamsQGM2F )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsQGM" );
	return true;
}


////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routine : ValidateCalibMethod
///	Returns : void
///	Action  : validate calibMethod
////////////////////////////////////////////////////

void ARM_QGM2F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	/// checks that if the CalibParam that we are trying 
	/// to calibrate exist into Model
	ARM_ModelParamVector CalibParams = calibMethod.GetCalibParams();
	size_t sizeCalibParams = CalibParams.size();
	for(size_t i=0; i< sizeCalibParams; ++i)
	{
		if(!CalibParams[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "You are trying to validate calibMethod with modelParam NULL, please advice!" );

		if( !GetModelParams()->DoesModelParamExist(ARM_ModelParamType::ParamNb(CalibParams[i]->GetType())) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, ARM_USERNAME + ": a Calib Method should contain a model Param of pricing model!");
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_QGM2F
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_QGM2F::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_QGM2F Model!");
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

