/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW.cpp
 *
 *  \brief file for the HW Model
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */



/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"

/// this header comes firts as it includes some preprocessor constants!
#include "gpmodels/HW.h"
#include "gpmodels/ModelParamsHW.h"
#include "gpmodels/ModelParamsHW1F.h"
#include "gpmodels/ModelParamsHW2F.h"
#include "gpmodels/Local_Normal_Model.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringmanip.h"

/// gpinfra
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/irrate.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/vanilla_normal.h"

/// gpnumlib
#include "gpnumlib/gaussiananalytics.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"

/// gpnummethods
#include "gpnummethods/meanrevertingsampler.h"

/// to look at the aleas only
#include <iostream>
#include <fstream>
#include <ios>
#include <iomanip>
using namespace std;


CC_BEGIN_NAMESPACE( ARM )

#define FIRST_STATE_VARIABLE    0
const double VOL_LIMIT      = 0.000001;

const double ARM_HullWhite::VOL_NORMAL_MAX = 150.0;

////////////////////////////////////////////////////
///	Class  : HWBoundaryD1Function
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
HWBoundaryD1Function::HWBoundaryD1Function(int size)
: itsPolyCoef(size), itsExpCoef(size)
{}

HWBoundaryD1Function::~HWBoundaryD1Function()
{}

void HWBoundaryD1Function::SetPolyCoef(const vector< double >& polyCoef)
{
    if(itsPolyCoef.size() != polyCoef.size())
        itsPolyCoef.resize(polyCoef.size());

    if(itsExpCoef.size() == polyCoef.size())
    {
        for(int i=0;i<itsPolyCoef.size();++i)
            itsPolyCoef[i] = polyCoef[i] * itsExpCoef[i];
    }
    else
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
        "Inconsistency in swaption exercise boundary search");
    }
}

void HWBoundaryD1Function::SetExpCoef(const vector< double >& expCoef)
{
    if(itsExpCoef.size() == expCoef.size())
    {
        for(int i=0;i<itsExpCoef.size();++i)
            itsExpCoef[i] = expCoef[i];
    }
    else
    {
        itsExpCoef.resize(expCoef.size());
        for(int i=0;i<itsExpCoef.size();++i)
            itsExpCoef[i] = expCoef[i];
    }
}

double HWBoundaryD1Function::operator () ( double x ) const
{
    double value=0.0;
    for(int i=0;i<itsPolyCoef.size();++i)
        value += itsPolyCoef[i] * exp( itsExpCoef[i] * x );
    return value;
}


////////////////////////////////////////////////////
///	Class  : HWBoundaryFunction
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

HWBoundaryFunction::HWBoundaryFunction(int size)
: itsPolyCoef(size), itsExpCoef(size)
{
    itsBoundaryD1 = new HWBoundaryD1Function(size);
}

HWBoundaryFunction::~HWBoundaryFunction()
{
    delete itsBoundaryD1;
}

void HWBoundaryFunction::SetPolyCoef(const vector< double >& polyCoef)
{
    if(itsPolyCoef.size() != polyCoef.size())
        itsPolyCoef.resize(polyCoef.size());

    for(int i=0;i<polyCoef.size();++i)
        itsPolyCoef[i] = polyCoef[i];

    itsBoundaryD1->SetPolyCoef(polyCoef); 
}

void HWBoundaryFunction::SetExpCoef(const vector< double >& expCoef)
{
    if(itsExpCoef.size() != expCoef.size())
        itsExpCoef.resize(expCoef.size());

    for(int i=0;i<itsExpCoef.size();++i)
        itsExpCoef[i] = expCoef[i];

    itsBoundaryD1->SetExpCoef(expCoef); 
}

double HWBoundaryFunction::operator () ( double x ) const
{
    double value=0.0;
    for(int i=0;i<itsPolyCoef.size();++i)
        value += itsPolyCoef[i] * exp( itsExpCoef[i] * x );

    return value-1.0;
}

HWBoundaryD1Function* HWBoundaryFunction::Derivative() const
{
    return itsBoundaryD1;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HullWhite::ARM_HullWhite( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, const ARM_BoolVector& soFormulaFlags )
:	ARM_PricingModelIR(zc,params),itsSOFormulaFlags(soFormulaFlags)
{}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HullWhite::ARM_HullWhite(const ARM_HullWhite& rhs)
: ARM_PricingModelIR(rhs),itsSOFormulaFlags(rhs.itsSOFormulaFlags)
{
    // Copy class attributes if any
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_HullWhite::~ARM_HullWhite()
{
	DeletePointorVector<ARM_GP_Matrix>( itsLocalStdDevs );
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HullWhite& ARM_HullWhite::operator=(const ARM_HullWhite& rhs)
{
	if(this != &rhs)
	{
		ARM_PricingModelIR::operator=(rhs);
		itsSOFormulaFlags = rhs.itsSOFormulaFlags;
		// Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: Init
///	Returns: ARM_PricingStatesPtr
///	Action : Default initialisation of the model and the
///          associated numerical method
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HullWhite::Init(const string& payModelName, const ARM_TimeInfoPtrVector& timeInfos)
{
    int nbEvents=timeInfos.size();
    bool isSpotUse = nbEvents == 0 || (nbEvents==1 && timeInfos[0]->GetEventTime() <= K_NEW_DOUBLE_TOL);
	
    if(!isSpotUse)
    {
        // Numerical method and numeraire are needed to go on
        ARM_NumMethodPtr numMethod=GetNumMethod();
		if( numMethod == ARM_NumMethodPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numerical method not set in HW model!");

        /// test the numeraire and its type!
		ARM_NumerairePtr numeraire=GetNumeraire();
        if( numeraire == ARM_NumerairePtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": numeraire not set in the HW model!");

		/// check that sampler type is MeanReverting
		ARM_SamplerBase* sampler = numMethod->GetSampler();
		if (sampler)
		{
			ARM_MeanRevertingSampler1D* sampler1D = dynamic_cast<ARM_MeanRevertingSampler1D*>(sampler);
			ARM_MeanRevertingSamplerND* samplerND = dynamic_cast<ARM_MeanRevertingSamplerND*>(sampler);
			if (!sampler1D && !samplerND)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Hull & White model : a MeanRevertingSampler is required!");
		}
		
		/// creates the model schedule (smart pointor for exception safety!)
		ARM_DiscretisationScheme& discretisationScheme = ARM_EventTime();
		CC_NS(std,auto_ptr)<std::vector<double>> ptimeSteps( discretisationScheme.ModelTimesFromTimesInfo(timeInfos, *this) );

        //// Initialise the numeraire
		numeraire->Init(*(GetDiscountFunctor()),payModelName,timeInfos);

		/// Set the basic schedule in the numerical method and...
		numMethod->SetTimeSteps(*ptimeSteps);

		double firstInductTime = timeInfos[0]->GetEventTime();

		ARM_PricingStatesPtr firstStates = numMethod->Init(*this,firstInductTime);

		std::vector<double>& timeSteps = numMethod->GetTimeSteps();

		ARM_MatrixVector localVariances;

		NumMethodStateLocalVariances(*timeSteps, localVariances);
		ComputestdDevMatrixVector( localVariances, itsLocalStdDevs, NeedsToCholeskyDecomposeFactors() );

		DeletePointorVector<ARM_GP_Matrix>( localVariances );

		ARM_GP_MatrixPtr absoluteDrifts;

		IntegratedLocalDrifts(
				*timeSteps,
				itsRelativeDrifts,
				absoluteDrifts);

		/// ...initialise it
		return firstStates;
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
///	Class   : ARM_PricingModel
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : Default implementation is no markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_HullWhite::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix(numMethodStates->rows(), numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{
    ARM_VectorPtr ZcStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
    ARM_VectorPtr ZcEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

	/// handle libor with no convexity correction
	if (	(	evalTime <= K_NEW_DOUBLE_TOL 
	 	||	states   == ARM_PricingStatesPtr(NULL) ) 
		&&  (-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		double libor = ((*ZcStart)[0]/(*ZcEnd)[0]-1.0)/period;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,libor) );
    }

    int i,nbStates=states->size();
    ARM_VectorPtr values( new std::vector<double>(nbStates) );
    if(-5 <= fwdEndTime - payTime && fwdEndTime - payTime <= 5)
    {
        /// No convexity
        for(i=0;i<nbStates;++i)
            (*values)[i]=((*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
    }
    else
    {
        /// Specialised version to compute payment convexity
        std::vector<double> unusedVars;
        double payConvex = exp( ((const ARM_ModelParamsHW* const)GetModelParams())->FwdZcLocalCovariance(evalTime,fwdResetTime,fwdStartTime,fwdEndTime,payTime,fwdEndTime,unusedVars) );
        for(i=0;i<nbStates;++i)
            (*values)[i]=(payConvex*(*ZcStart)[i]/(*ZcEnd)[i]-1.0)/period;
    }

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for caplet/floorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite::VanillaCaplet(
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
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:25:45):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

    double amount=payNotional;
    if(period != fwdPeriod)
        amount *= period/fwdPeriod;

    // Fwd Zc volatility computation (ARM_ModelParamsHW1F class is pure virtual)
    double zcStdDev=sqrt(((const ARM_ModelParamsHW* const)GetModelParams())->FwdZcLocalVariance(evalTime,fwdResetTime,fwdStartTime,fwdEndTime));

    // Payment convexity
    double payConvex=1.0;
    if(fwdEndTime - payTime < -5 ||  fwdEndTime - payTime > 5)
    {
        std::vector<double> unusedVars;
        payConvex = exp( ((const ARM_ModelParamsHW* const)GetModelParams())->FwdZcLocalCovariance(evalTime,fwdResetTime,fwdStartTime,fwdEndTime,payTime,fwdEndTime,unusedVars) );
    }

    /// Option computation
    ARM_VectorPtr zcFwdStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
    ARM_VectorPtr zcFwdEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

    /// Payment Zc
	ARM_VectorPtr zcPay;
	if( fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        zcPay = zcFwdEnd;
    else
		zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);


    int i,nbStates=zcPay->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates));

    for(i=0;i<nbStates;i++)
        (*values)[i]=amount*BlackSholes_Formula(payConvex*(*zcFwdStart)[i]/(*zcFwdEnd)[i],zcStdDev,(*zcPay)[i], 1.0+fwdPeriod*strikesPerState[i],capFloor);

    return values;
}



////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: VanillaDigital
///	Returns: a vector of Digital(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          digital caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite::VanillaDigital(
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
	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:25:56):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

    double amount=payNotional*period;

    // Fwd Zc volatility computation (ARM_ModelParamsHW1F class is pure virtual)
    double zcStdDev=sqrt(((const ARM_ModelParamsHW* const)GetModelParams())->FwdZcLocalVariance(evalTime,fwdResetTime,fwdStartTime,fwdEndTime));

    // Payment convexity
    double payConvex=1.0;
    if(fwdEndTime - payTime < -5 ||  fwdEndTime - payTime > 5)
    {
        std::vector<double> unusedVars;
        payConvex = exp( ((const ARM_ModelParamsHW* const)GetModelParams())->FwdZcLocalCovariance(evalTime,fwdResetTime,fwdStartTime,fwdEndTime,payTime,fwdEndTime,unusedVars) );
    }

    /// Option computation
    ARM_VectorPtr zcFwdStart= GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
    ARM_VectorPtr zcFwdEnd  = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);

    /// Payment Zc
	ARM_VectorPtr zcPay;
	if( fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        zcPay = zcFwdEnd;
    else
		zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

    int i,nbStates=zcPay->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates));

	double d2;
    for(i=0;i<nbStates;i++)
	{
		d2		= (log(payConvex*(*zcFwdStart)[i]/(*zcFwdEnd)[i]/(1.0+fwdPeriod*strikesPerState[i]))-zcStdDev*zcStdDev*.5)/zcStdDev;
        (*values)[i]=amount*(*zcPay)[i]*ARM_GaussianAnalytics::cdfNormal(capFloor*d2);
	}
    return values;
}



////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite::VanillaSwaption(
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
	/// TO BE UPDATED
	/// Check that the notional is constant
	double swapNotional = fixNotional[0];
	if (!(isConstantNotional&&isConstantSpread&&isConstantStrike))
				ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with variable notional, Spread or Strike!" );

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:26:07):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
#if defined(__GP_STRICT_VALIDATION)
	if( strikesPerState.GetRowsNb() != states->size() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": strike per states size != states size" );
	if( strikesPerState.GetColsNb() != fixPayTimes.size() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": strike nb != fixed flow nb" );
#endif

    // It will be usefull to test if some datas (fwdZcVariance(t,T), beta(t,Ti)...)
    // are already computed and saved... especially in MC use

    if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + 
            "Swaption pricing not implemented for differents discount & fixing curves" );

	/// check that strikePerState is meaningful

    /// Compute Zc and bond datas
    ARM_VectorPtr zcFloatStart;
    vector< ARM_VectorPtr > zcFlowCouponPay;
    vector< double > couponPayPeriod;
    vector< double > flowPayTimes;

    ComputeFlowZc(curveName,evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
		states,zcFloatStart,zcFlowCouponPay,flowPayTimes,couponPayPeriod);

    // Volatility computation (ARM_ModelParamsHW class is pure virtual)
    int i,nbFlows=zcFlowCouponPay.size();
    vector< double > stdDev(nbFlows);
    vector< double > drift(nbFlows);
    for(i=0;i<nbFlows;++i)
    {
        stdDev[i]	= ((const ARM_ModelParamsHW* const) GetModelParams())->FwdZcLocalVariance(evalTime,swapResetTime,floatStartTime,flowPayTimes[i]);
        drift[i]	= exp(-0.5*stdDev[i]);
        stdDev[i]	= sqrt(stdDev[i]);
    }

    /// Compute prices
    return ARM_VectorPtr( VanillaSwaptionPrice(fixNotional,floatNotional,zcFloatStart,zcFlowCouponPay,
        couponPayPeriod,stdDev,drift,strikesPerState,callPut,states,isConstantNotional,isConstantSpread,isConstantStrike) );
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: ComputeFlowZc
///	Returns: floating leg start, payment date Zc,
///          flow times and bond coefficients
///          (by argument references)
///	Action : Compute zero coupond prices at evalTime
///          , flow times and bond coefficients
////////////////////////////////////////////////////
void ARM_HullWhite::ComputeFlowZc(
		const string& curveName, 
		double evalTime,
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
        const ARM_PricingStatesPtr& states,
		ARM_VectorPtr& zcFloatStart,
		vector< ARM_VectorPtr >& zcFlowCouponPay,
        vector< double >& flowPayTimes,
        vector< double >& couponPayPeriod ) const
{
    int i,nbFixFlows=fixPayTimes.size();
    int nbFlows=nbFixFlows;
    if(floatEndTime != fixPayTimes[nbFixFlows-1])
        ++nbFlows;
    bool isSameEnd=(nbFlows == nbFixFlows);

    flowPayTimes.resize(nbFlows);
    couponPayPeriod.resize(nbFlows);
    zcFlowCouponPay.resize(nbFlows);

    zcFloatStart = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);

    for(i=0;i<nbFixFlows;++i)
    {
        flowPayTimes[i]	   = fixPayTimes[i];
        couponPayPeriod[i] = fixPayPeriods[i]; 
	    zcFlowCouponPay[i] = GetFixingFunctor()->DiscountFactor(curveName,evalTime,flowPayTimes[i],states);
	}

	/// no difference between last fixPayTimes and floatEndTime!
	if(isSameEnd)
    {
        flowPayTimes[nbFlows-1]		= fixPayTimes[nbFlows-1];
        couponPayPeriod[nbFlows-1]	= fixPayPeriods[nbFlows-1]; 
	    zcFlowCouponPay[nbFlows-1]	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,flowPayTimes[nbFlows-1],states);
    }
	else
    {
        // Add the last notional flow at the end of flow list
        flowPayTimes[nbFlows-1]		= floatEndTime;
        couponPayPeriod[nbFlows-1]	= 1.0;
	    zcFlowCouponPay[nbFlows-1]	= GetFixingFunctor()->DiscountFactor(curveName,evalTime,flowPayTimes[nbFlows-1],states);
	}

}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: isFwdStart
///	Returns: bool
///	Action : checks notional to identify fwd start swaptions
////////////////////////////////////////////////////
bool ARM_HullWhite::isFwdStart(const std::vector<double>& fixNotional, int& idx) const
{
	bool ok  = false;
	bool cst = true;
	int lastPosIdx = -1;
	for (int i=1;i<fixNotional.size();i++)
	{
		if (fabs(fixNotional[i]-fixNotional[i-1])>K_NEW_DOUBLE_TOL)
		{
			if (fabs(fixNotional[i-1])<K_NEW_DOUBLE_TOL && cst)
			{
				lastPosIdx = i-1;
				ok = true;
			}
			else
				ok = false;

			cst = false;
		}
	}
	if (ok)
		idx = lastPosIdx;
	return ok;
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: VanillaSwaptionPrice
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite::VanillaSwaptionPrice(
		const std::vector<double>& fixNotional,
		const std::vector<double>& floatNotional,
        const ARM_VectorPtr& zcFloatStart,
        const vector< ARM_VectorPtr >& zcFlowCouponPay,
        const vector< double >& couponPayPeriod,
        const vector< double >& stdDev,
        const vector< double >& drift,
		const ARM_GP_Matrix& strikesPerState,
        int callPut,
		const ARM_PricingStatesPtr& states,
		bool isConstantNotional,
		bool isConstantSpread,
		bool isConstantStrike,
		int idx) const
{
	/// Check that the notional is constant
	if (!(isConstantSpread && isConstantStrike))
			ARM_THROW( ERR_INVALID_ARGUMENT, "The Model can not price a swaption with Spread or Strike!" );

	/// Take care of the case where FixNotio != FloatNotio !!
	double swapNotional	= floatNotional[floatNotional.size()-1];
	double notioRatio	= fixNotional[fixNotional.size()-1]/swapNotional;

	bool isFwdStartNotional=(idx>-1);
	
	if (!(isConstantNotional || isFwdStartNotional))
		ARM_THROW( ERR_INVALID_ARGUMENT, "swaption with variable notional only for fwd start" );
	
	/// for parsing of the function using dummyStates!
	if( states == ARM_PricingStatesPtr(NULL) )
// FIXMEFRED: mig.vc8 (25/05/2007 15:26:45):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));


    const double stdDevLimit = 1.0e-10;

    int payRec,callPutBond;
    if(callPut == K_CALL)
    {
        payRec=K_PAY;
        callPutBond=K_PUT;
    }
    else
    {
        payRec=K_RCV;
        callPutBond=K_CALL;
    }

    int nbStates=zcFloatStart->size();
    ARM_VectorPtr values = ARM_VectorPtr(new std::vector<double>(nbStates));

    // Flow discount factors
    int i,nbFlows = zcFlowCouponPay.size();
    vector< double > expCoef(nbFlows);
    vector< double > logNorDrift(nbFlows);
    bool isAllNull=true;
    for(i=idx+1;i<nbFlows;++i)
    {
        expCoef[i]=stdDev[i];
        logNorDrift[i]=exp(0.5*stdDev[i]*stdDev[i]);
        if(isAllNull && (stdDev[i] < -stdDevLimit || stdDevLimit < stdDev[i]))
           isAllNull=false;
    }

    /// Remember that last notional flow and last fixed flow date may differ
    int nbFixedFlows = strikesPerState.cols();
    int flowNotionalIndex = nbFlows-1;

    // Solver init for exercise boundary search
    // (one for all iterations on states)
    HWBoundaryFunction boundaryFct(nbFlows);
    boundaryFct.SetExpCoef(expCoef);

    /// Create a local solver
	const double target = 0.0;
	const double ftol = 1.0e-8;
	const double xtol = 1.0e-8;
	const size_t max_iter = 30;
    T_NewtonRaphsonSolverNoThrow<HWBoundaryFunction> solver(boundaryFct,target,ftol,xtol,max_iter);

    double sum,weightedSum,x0,eqStrike,eqFwd;
    ARM_VectorPtr curzcFlowCouponPay;
    vector< double > polyCoef(nbFlows);
	double polyCoeffNotional;
    try
    {
        for(int stIdx=0;stIdx<nbStates;++stIdx)
        {
            double strikemoy =0.0;
            for(i=idx+1;i<nbFixedFlows;++i)
                strikemoy += strikesPerState(stIdx,i);
            strikemoy/=nbFixedFlows;

            if(isAllNull )
            {
                // Always ITM or OTM
                sum=0.0;
                for(i=idx+1;i<nbFixedFlows;++i)
                {
                    curzcFlowCouponPay=zcFlowCouponPay[i];
                    sum += notioRatio*couponPayPeriod[i]*strikesPerState(stIdx,i)*(*curzcFlowCouponPay)[stIdx];
                }
				/// last bit from notional
                sum += (*(zcFlowCouponPay[flowNotionalIndex]))[stIdx];

                sum = payRec*(sum - (*zcFloatStart)[stIdx]);
                if(sum<0.0)
                    sum=0.0;
            }
			/// Handle special case of high negative strikes (-7%)
			else if( strikemoy < -0.07 )
            {
				if( payRec == K_RCV )
				{
					sum = 0.0;
				}
				else
				{
					sum=0.0;
					for(i=idx+1;i<nbFixedFlows;++i)
					{
						curzcFlowCouponPay=zcFlowCouponPay[i];
						sum += notioRatio*couponPayPeriod[i]*strikesPerState(stIdx,i)*(*curzcFlowCouponPay)[stIdx];
					}
					/// last bit from notional
					sum += (*(zcFlowCouponPay[flowNotionalIndex]))[stIdx];
					sum = payRec*(sum - (*zcFloatStart)[stIdx]);
				}
            }
            else
            {
                // Coefficient for exercise boundary function
                sum=0.0;
                weightedSum=0.0;
                for(i=idx+1;i<nbFixedFlows;++i)
                {
                    curzcFlowCouponPay=zcFlowCouponPay[i];
                    polyCoef[i]=notioRatio*couponPayPeriod[i]*strikesPerState(stIdx,i)*drift[i]*(*curzcFlowCouponPay)[stIdx]/(*zcFloatStart)[stIdx];

                    sum += polyCoef[i];
                    weightedSum += polyCoef[i]*expCoef[i];
                }
				
				/// last bit from notional
                polyCoeffNotional = drift[flowNotionalIndex]* (*(zcFlowCouponPay[flowNotionalIndex]))[stIdx]/(*zcFloatStart)[stIdx];
				if(flowNotionalIndex == nbFixedFlows-1)
                    polyCoef[flowNotionalIndex] += polyCoeffNotional;
                else
                    polyCoef[flowNotionalIndex] = polyCoeffNotional;

                sum += polyCoeffNotional;
                weightedSum += polyCoeffNotional*expCoef[flowNotionalIndex];


                // Update context for solving and find exercise boundary
                boundaryFct.SetPolyCoef(polyCoef);
                if(-stdDevLimit <= weightedSum && weightedSum <= stdDevLimit)
                    x0=0.0;
                else
                {
                    /// Locale a smart initial guess beginning with
                    /// a Taylor expansion near 0
                    x0=(1.0-sum)/weightedSum;
                    if(x0<-0.75 || x0>0.75)
                    {
                         /// Approx not relevant then locate roughly
                        double dx0=0.75;
                        x0 = (x0>0 ? dx0 : -dx0);
                        double fx0=boundaryFct(x0);
                        double signfx0 = (fx0>0 ? 1 : -1);
                        double dfx0=(*(boundaryFct.Derivative()))(x0);
                        double step=((fx0>0 && dfx0>0) || (fx0<0 && dfx0<0) ? -dx0 : dx0);
                        while(-4.5 < x0 && x0 < 4.5)
                        {
                            x0 += dx0;
                            fx0=boundaryFct(x0);
                            if(signfx0 != (fx0>0 ? 1 : -1))
                                break;
                        }
                    }
                }
                
                try
                {
                    solver.setInitialGuess(x0);

                    x0=solver.Solve();

                    // Sum each option on fwd Zc at equivalent strike
                    sum=0.0;
                    for(i=idx+1;i<nbFlows;++i)
                    {
                        curzcFlowCouponPay = zcFlowCouponPay[i];
                        if( i!= flowNotionalIndex )
					    {
						    if( polyCoef[i] < -1.0e-14 )
						    {
								/// Take care of the case polyCoef[i] < 0 !!
							    eqStrike	= exp(expCoef[i]*x0);
							    eqFwd		= logNorDrift[i];
							    sum			+= polyCoef[i]*BlackSholes_Formula(eqFwd,stdDev[i],(*zcFloatStart)[stIdx],eqStrike,callPutBond);
						    }
						    else if( polyCoef[i] > 1.0e-14 )
						    {
								/// Because a.BS(F,K,StdDev) may slightly differ from BS(a.F,a.K,StdDev)...
								/// and cause regression !
							    eqStrike	= polyCoef[i]*exp(expCoef[i]*x0);
							    eqFwd		= polyCoef[i]*logNorDrift[i];
							    sum			+= BlackSholes_Formula(eqFwd,stdDev[i],(*zcFloatStart)[stIdx],eqStrike,callPutBond);
						    }
					    }
					    /// part with notional ... hence condition always met!
					    else
					    {
						    if( polyCoef[i] < -1.0e-14 )
						    {
								/// Take care of the case polyCoef[i] < 0 !!
								eqStrike	= exp(expCoef[i]*x0);
								eqFwd		= logNorDrift[i];
								sum			+= polyCoef[i]*BlackSholes_Formula(eqFwd,stdDev[i],(*zcFloatStart)[stIdx],eqStrike,callPutBond);
							}
						    else if( polyCoef[i] > 1.0e-14 )
						    {
								eqStrike	= polyCoef[i]*exp(expCoef[i]*x0);
								eqFwd		= polyCoef[i]*logNorDrift[i];
								sum			+= BlackSholes_Formula(eqFwd,stdDev[i],(*zcFloatStart)[stIdx],eqStrike,callPutBond);
							}
					    }
                    }
                }
                catch(Exception& err)
                {
                    if(strikemoy < 0.0)
                    {
                        /// Still pb with negative strikes then use intrinic value
				        if( payRec == K_RCV )
				        {
					        sum = 0.0;
				        }
				        else
				        {
					        sum=0.0;
					        for(i=idx+1;i<nbFixedFlows;++i)
					        {
						        curzcFlowCouponPay=zcFlowCouponPay[i];
						        sum += notioRatio*couponPayPeriod[i]*strikesPerState(stIdx,i)*(*curzcFlowCouponPay)[stIdx];
					        }
					        /// last bit from notional
					        sum += (*(zcFlowCouponPay[flowNotionalIndex]))[stIdx];
					        sum = payRec*(sum - (*zcFloatStart)[stIdx]);
				        }
                    }
                    else
                        throw(err);
                }
           }
            (*values)[stIdx]=swapNotional*sum;

        } // for state idx

    } // try

    catch(Exception& err)
	{
        double signPoly = (polyCoef[0] > 0 ? 1 : -1);
        double signExp = (expCoef[0] > 0 ? 1 : -1);
        bool isCstSignPoly=true;
        bool isCstSignExp=true;
        for(int k=1;k<nbFlows;++k)
        {
            if(isCstSignPoly && signPoly != (polyCoef[k] > 0 ? 1 : -1))
                isCstSignPoly=false;
            if(isCstSignExp && signExp != (expCoef[k] > 0 ? 1 : -1))
                isCstSignExp=false;
        }

        CC_Ostringstream os;
        os << "Max iter reached in boundary search (H&W) : "
            << (isCstSignPoly ?  "PolySign=Cst," : "PolySign=NotCst,")
            << (isCstSignExp ?  "ExpSign=Cst" : "ExpSign=NotCst");
        os << ", Previous error msge : " << err.GetMessage() ;
        throw Exception( __LINE__, __FILE__, ERR_INVALID_ARGUMENT,  os.str() );
    }

    return values;
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routine : ComputeModelTimes
///	Returns : an empty vector since in HW there is not
///				such a thing as model times
////////////////////////////////////////////////////
std::vector<double>& ARM_HullWhite::ComputeModelTimes(  const ARM_TimeInfoPtrVector& timeInfos )
{
	/// since there is no concept of model time
	/// returns an empty vector
	return new std::vector<double>(0);
}



////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_HullWhite::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
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
	double currentState;
	size_t modelNb	= GetModelNb();
	
	for( size_t i=0;i<statesNb; ++i )
	{
		for( size_t j=0;  j<factorsNb; ++j )
		{
			currentState = states->GetNumMethodState(i,modelNb+j);
/*
			ofstream testfile("C:\\Documents and Settings\\ocroissant\\My Documents\\toto.txt",ios_base::app);
			testfile.setf(ios::fixed, ios::floatfield);
			testfile<<setprecision(18)<<",{"<<currentState<<"}";
			testfile.close();
*/

			states->SetModelState(i,j+modelNb,currentState);
		}
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : Take care this function is just used to 
/// make the equity model work !!!!
/// It will diseapear quite soon.
////////////////////////////////////////////////////

void ARM_HullWhite::MCModelStatesFromToNextTimeSpecialEquity(ARM_PricingStatesPtr& states,int timeIndex) const
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
///	Class  : ARM_HullWhite
///	Routine: TreeStatesToModelStates
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_HullWhite::TreeStatesToModelStates(ARM_PricingStatesPtr& states, int timeIndex) const
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
///	Class   : ARM_HullWhite
///	Routines: UnderlyingCorrelation
///	Returns : void
///	Action  :  compute the covariance between two ZC Forward
////////////////////////////////////////////////////

double ARM_HullWhite::UnderlyingCorrelation(   string	underlyingType,
											   double	fromTime,
											   double	toTime,
											   double	startTime1,
											   double   endTime1,
											   double	startTime2,
											   double   endTime2,
											   double	startTime3,
											   double   endTime3,
											   double	startTime4,
											   double   endTime4) const
{
	stringGetUpper(underlyingType);

	if (underlyingType == "ZC")
	{
		double varCovarStepTime = 0.0001;
		if(fromTime == toTime)
		{
			/// Instantaneous values
			if(fromTime>=varCovarStepTime) fromTime-=varCovarStepTime;
			else if(toTime<=MIN(startTime1,startTime2)-varCovarStepTime) toTime=fromTime+varCovarStepTime;
		}

		/// Call the local variance & covariance functions
		ARM_ModelParams* modParamsHW = const_cast<ARM_HullWhite* const>(this)->GetModelParams();
		//ARM_ModelParamsHW *modParamsHW;
		double covar;
		std::vector<double> vars(2);
		if((modParamsHW = dynamic_cast<ARM_ModelParamsHW *>(modParamsHW)) != NULL)
		{
			covar= dynamic_cast<ARM_ModelParamsHW *>(modParamsHW)->FwdZcLocalCovariance(fromTime,toTime,startTime1,endTime1,startTime2,endTime2,vars);
			if(vars[0] < K_DOUBLE_TOL || vars[1] < K_DOUBLE_TOL)
			{
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
						"variances too low HW Model!");		   
			}
			else
				return covar/sqrt(vars[0]*vars[1]);

		}
		else
			return 0.0;
	}
	else if (underlyingType == "CMS" || underlyingType == "SO/CMS" || underlyingType == "SO/SO")
	{
		ARM_ModelParamsHW2F* params2F = dynamic_cast<ARM_ModelParamsHW2F*>((ARM_ModelParams*)GetModelParams());
		

		ARM_Currency* ccy	= GetZeroCurve()->GetCurrencyUnit();
		int stdFixFreq		= ccy->GetFixedPayFreq();
		int stdFixDayCount	= ccy->GetFixedDayCount();
		
		ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
		char* resetCalendar		 = ccy->GetResetCalName(indexType);
		char* payCalendar		 = ccy->GetPayCalName(indexType);
		int resetGap			 = - ccy->GetSpotDays(); 
		int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
		int intRule				 = K_ADJUSTED;			// for interest dates
		int stubRule			 = K_SHORTSTART;
		int resetTiming			 = K_ADVANCE;
		int payTiming			 = K_ARREARS;
				
		double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
		double longStart  = asOf + startTime1;
		double longEnd    = asOf + endTime1;
		double shortStart = asOf + startTime2;
		double shortEnd   = asOf + endTime2;
		
		ARM_DateStrip longSched  (longStart,  longEnd,  stdFixFreq, stdFixDayCount, resetCalendar, fwdRule, intRule, stubRule, resetGap, stdFixFreq);
		ARM_DateStrip shortSched (shortStart, shortEnd, stdFixFreq, stdFixDayCount, resetCalendar, fwdRule, intRule, stubRule, resetGap, stdFixFreq);

		std::vector<double> longFixPayDfs    (longSched.GetFlowEndDates()->size()) ;
		std::vector<double> longFixPayTimes  (longSched.GetFlowEndDates()->size()) ;
		std::vector<double> shortFixPayDfs   (shortSched.GetFlowEndDates()->size()) ;
		std::vector<double> shortFixPayTimes (shortSched.GetFlowEndDates()->size()) ;


		/// compute dfs
		for (size_t i(0); i<longFixPayDfs.size(); i++)
		{
			longFixPayDfs[i]   = GetZeroCurve()->DiscountPrice((longSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
			longFixPayTimes[i] = longSched.GetFlowEndDates()->Elt(i) - asOf;
		}
			
		for (i=0; i<shortFixPayDfs.size(); i++)
		{
			shortFixPayDfs[i] = GetZeroCurve()->DiscountPrice((shortSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
			shortFixPayTimes[i] = shortSched.GetFlowEndDates()->Elt(i) - asOf;
		}

		double longFloatStartTime	= longSched.GetFlowStartDates()->Elt(0) - asOf;
		double longFloatEndTime		= longSched.GetFlowEndDates()->Elt(longFixPayDfs.size()-1) - asOf;
		double longDfStart			= GetZeroCurve()->DiscountPrice((longSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
		double longDfEnd			= GetZeroCurve()->DiscountPrice((longSched.GetFlowEndDates()->Elt(longFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);
		double shortFloatStartTime	= shortSched.GetFlowStartDates()->Elt(0) - asOf;
		double shortFloatEndTime	= shortSched.GetFlowEndDates()->Elt(shortFixPayDfs.size()-1) - asOf;
		double shortDfStart			= GetZeroCurve()->DiscountPrice((shortSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
		double shortDfEnd			= GetZeroCurve()->DiscountPrice((shortSched.GetFlowEndDates()->Elt(shortFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);

		double lambda1 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);

		double longVolFactor1, longVolFactor2;
		double shortVolFactor1, shortVolFactor2;
		
		longVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																toTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longDfStart, 
																longDfEnd, 
																longFixPayTimes, 
																*longSched.GetInterestTerms(), 
																longFixPayDfs);
		
		shortVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																toTime, 
																shortFloatStartTime, 
																shortFloatEndTime, 
																shortDfStart, 
																shortDfEnd, 
																shortFixPayTimes, 
																*shortSched.GetInterestTerms(), 
																shortFixPayDfs);

		std::vector<double> rateLongFixPayDfs,rateLongFixPayTimes,rateLongInterestTerms;
		double rateLongFloatStartTime,rateLongFloatEndTime,rateLongDfStart,rateLongDfEnd;
		if(underlyingType == "SO/CMS" || underlyingType == "SO/SO")
		{
			/// CMS Spread = LongCMS - ShortCMS
			longVolFactor1 -= shortVolFactor1;

			double rateLongStart = asOf + startTime3;
			double rateLongEnd   = asOf + endTime3;
			ARM_DateStrip rateLongSched(rateLongStart,rateLongEnd,stdFixFreq,stdFixDayCount,resetCalendar,fwdRule,intRule,stubRule,resetGap,stdFixFreq);

			rateLongFixPayDfs.resize(rateLongSched.GetFlowEndDates()->size());
			rateLongFixPayTimes.resize(rateLongSched.GetFlowEndDates()->size());
			rateLongInterestTerms = *rateLongSched.GetInterestTerms();
			for(i=0;i<rateLongFixPayDfs.size();i++)
			{
				rateLongFixPayDfs[i] = GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
				rateLongFixPayTimes[i] = rateLongSched.GetFlowEndDates()->Elt(i) - asOf;
			}
			rateLongFloatStartTime	= rateLongSched.GetFlowStartDates()->Elt(0) - asOf;
			rateLongFloatEndTime	= rateLongSched.GetFlowEndDates()->Elt(rateLongFixPayDfs.size()-1) - asOf;
			rateLongDfStart			= GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
			rateLongDfEnd			= GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowEndDates()->Elt(rateLongFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);

			shortVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																	toTime, 
																	rateLongFloatStartTime, 
																	rateLongFloatEndTime, 
																	rateLongDfStart, 
																	rateLongDfEnd, 
																	rateLongFixPayTimes, 
																	rateLongInterestTerms, 
																	rateLongFixPayDfs);
		}
		std::vector<double> rateShortFixPayDfs,rateShortFixPayTimes,rateShortInterestTerms;
		double rateShortFloatStartTime,rateShortFloatEndTime,rateShortDfStart,rateShortDfEnd;
		if(underlyingType == "SO/SO")
		{
			double rateShortStart = asOf + startTime4;
			double rateShortEnd   = asOf + endTime4;
			ARM_DateStrip rateShortSched(rateShortStart,rateShortEnd,stdFixFreq,stdFixDayCount,resetCalendar,fwdRule,intRule,stubRule,resetGap,stdFixFreq);

			rateShortFixPayDfs.resize(rateShortSched.GetFlowEndDates()->size());
			rateShortFixPayTimes.resize(rateShortSched.GetFlowEndDates()->size());
			rateShortInterestTerms = *rateShortSched.GetInterestTerms();
			for(i=0;i<rateShortFixPayDfs.size();i++)
			{
				rateShortFixPayDfs[i] = GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
				rateShortFixPayTimes[i] = rateShortSched.GetFlowEndDates()->Elt(i) - asOf;
			}
			rateShortFloatStartTime	= rateShortSched.GetFlowStartDates()->Elt(0) - asOf;
			rateShortFloatEndTime	= rateShortSched.GetFlowEndDates()->Elt(rateShortFixPayDfs.size()-1) - asOf;
			rateShortDfStart			= GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
			rateShortDfEnd			= GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowEndDates()->Elt(rateShortFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);

			shortVolFactor1 -= ARM_HullWhite::SwapRateVolFactor (lambda1, 
																	toTime, 
																	rateShortFloatStartTime, 
																	rateShortFloatEndTime, 
																	rateShortDfStart, 
																	rateShortDfEnd, 
																	rateShortFixPayTimes, 
																	rateShortInterestTerms, 
																	rateShortFixPayDfs);
		}

		/// HW 1F case
		if (params2F == NULL)
			return (longVolFactor1*shortVolFactor1 > 0.0 ? 1.0 : -1.0);

		/// HW 2F case
		double lambda2 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversionSpread).GetValueAtPoint(0);
		lambda2 += lambda1;

		longVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																toTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longDfStart, 
																longDfEnd, 
																longFixPayTimes, 
																*longSched.GetInterestTerms(), 
																longFixPayDfs);

		shortVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																toTime, 
																shortFloatStartTime, 
																shortFloatEndTime, 
																shortDfStart, 
																shortDfEnd, 
																shortFixPayTimes, 
																*shortSched.GetInterestTerms(), 
																shortFixPayDfs);
		

		ARM_GP_TriangularMatrix* phi = params2F->StateLocalVariance(fromTime, toTime);
		double phi11 = (*phi)(0, 0);
		double phi22 = (*phi)(1, 1);
		double phi12 = (*phi)(1, 0);
		delete phi;

		if(underlyingType == "SO/CMS" || underlyingType == "SO/SO")
		{
			/// CMS Spread = LongCMS - ShortCMS
			longVolFactor2 -= shortVolFactor2;

			shortVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																	toTime, 
																	rateLongFloatStartTime, 
																	rateLongFloatEndTime, 
																	rateLongDfStart, 
																	rateLongDfEnd, 
																	rateLongFixPayTimes, 
																	rateLongInterestTerms, 
																	rateLongFixPayDfs);
		}

		if(underlyingType == "SO/SO")
		{
			shortVolFactor2 -= ARM_HullWhite::SwapRateVolFactor (lambda2, 
																	toTime, 
																	rateShortFloatStartTime, 
																	rateShortFloatEndTime, 
																	rateShortDfStart, 
																	rateShortDfEnd, 
																	rateShortFixPayTimes, 
																	rateShortInterestTerms, 
																	rateShortFixPayDfs);
		}

		double correl =		 longVolFactor1 * shortVolFactor1 * phi11
						+	 longVolFactor2 * shortVolFactor2 * phi22
						+	(longVolFactor1 * shortVolFactor2 + longVolFactor2 * shortVolFactor1) * phi12 ;

		correl /= sqrt(		longVolFactor1 * longVolFactor1 * phi11
						+	longVolFactor2 * longVolFactor2 * phi22
						+	2.0 * longVolFactor1 * longVolFactor2 * phi12 );

		correl /= sqrt(		shortVolFactor1 * shortVolFactor1 * phi11
						+	shortVolFactor2 * shortVolFactor2 * phi22
							+	2.0 * shortVolFactor1 * shortVolFactor2 * phi12 );

		return correl;
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "UnderlyingCorrelation : Only ZC, CMS or CMS Spread underlyingType are available for  HW Model");
		return 0.0;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routines: UnderlyingCovariance
///	Returns : void
///	Action  :  compute the covariance between two ZC forwards
////////////////////////////////////////////////////

double ARM_HullWhite::UnderlyingCovariance(    string	underlyingType,
											   double	fromTime,
											   double	toTime,
											   double	startTime1,
											   double   endTime1,
											   double	startTime2,
											   double   endTime2,
											   double	startTime3,
											   double   endTime3,
											   double	startTime4,
											   double   endTime4) const
{	
	stringGetUpper(underlyingType);

	if (underlyingType == "ZC")
	{
		double varCovarStepTime = 0.0001;
		bool isInst;
		if(isInst=(fromTime == toTime))
		{
			/// Instantaneous values
			if(fromTime>=varCovarStepTime) fromTime-=varCovarStepTime;
			else if(toTime<=MIN(startTime1,startTime2)-varCovarStepTime) toTime=fromTime+varCovarStepTime;
		}

		/// Call the local variance & covariance functions
		ARM_ModelParams* modParamsHW = const_cast<ARM_HullWhite* const>(this)->GetModelParams();
		//ARM_ModelParamsHW *modParamsHW;
		double var;
		std::vector<double> vars(2);

		if((modParamsHW = dynamic_cast<ARM_ModelParamsHW *>(modParamsHW)) != NULL)
		{
			var= dynamic_cast<ARM_ModelParamsHW *>(modParamsHW)->FwdZcLocalCovariance(fromTime,toTime,startTime1,endTime1,startTime2,endTime2,vars);
			if(isInst)
			{
				var = var*K_YEAR_LEN/varCovarStepTime;		   
			}
			else
				return var;

		}
		
		return 0.0;
	}
	else if (underlyingType == "CMS" || underlyingType == "SO/CMS" || underlyingType == "SO/SO")
	{
		ARM_Currency* ccy	= GetZeroCurve()->GetCurrencyUnit();
		int stdFixFreq		= ccy->GetFixedPayFreq();
		int stdFixDayCount	= ccy->GetFixedDayCount();
		
		ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
		char* resetCalendar		 = ccy->GetResetCalName(indexType);
		char* payCalendar		 = ccy->GetPayCalName(indexType);
		int resetGap			 = - ccy->GetSpotDays(); 
		int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
		int intRule				 = K_ADJUSTED;			// for interest dates
		int stubRule			 = K_SHORTSTART;
		int resetTiming			 = K_ADVANCE;
		int payTiming			 = K_ARREARS;
				
		double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
		double longStart  = asOf + startTime1;
		double longEnd    = asOf + endTime1;
		double shortStart = asOf + startTime2;
		double shortEnd   = asOf + endTime2;
		
		ARM_DateStrip longSched  (longStart,  longEnd,  stdFixFreq, stdFixDayCount, resetCalendar, fwdRule, intRule, stubRule, resetGap, stdFixFreq);
		ARM_DateStrip shortSched (shortStart, shortEnd, stdFixFreq, stdFixDayCount, resetCalendar, fwdRule, intRule, stubRule, resetGap, stdFixFreq);

		std::vector<double> longFixPayDfs    (longSched.GetFlowEndDates()->size()) ;
		std::vector<double> longFixPayTimes  (longSched.GetFlowEndDates()->size()) ;
		std::vector<double> shortFixPayDfs   (shortSched.GetFlowEndDates()->size()) ;
		std::vector<double> shortFixPayTimes (shortSched.GetFlowEndDates()->size()) ;

		/// compute dfs
		for (size_t i(0); i<longFixPayDfs.size(); i++)
		{
			longFixPayDfs[i]   = GetZeroCurve()->DiscountPrice((longSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
			longFixPayTimes[i] = longSched.GetFlowEndDates()->Elt(i) - asOf;
		}
			
		for (i=0; i<shortFixPayDfs.size(); i++)
		{
			shortFixPayDfs[i] = GetZeroCurve()->DiscountPrice((shortSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
			shortFixPayTimes[i] = shortSched.GetFlowEndDates()->Elt(i) - asOf;
		}

		double longFloatStartTime	= longSched.GetFlowStartDates()->Elt(0) - asOf;
		double longFloatEndTime		= longSched.GetFlowEndDates()->Elt(longFixPayDfs.size()-1) - asOf;
		double longDfStart			= GetZeroCurve()->DiscountPrice((longSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
		double longDfEnd			= GetZeroCurve()->DiscountPrice((longSched.GetFlowEndDates()->Elt(longFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);
		double shortFloatStartTime	= shortSched.GetFlowStartDates()->Elt(0) - asOf;
		double shortFloatEndTime	= shortSched.GetFlowEndDates()->Elt(shortFixPayDfs.size()-1) - asOf;
		double shortDfStart			= GetZeroCurve()->DiscountPrice((shortSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
		double shortDfEnd			= GetZeroCurve()->DiscountPrice((shortSched.GetFlowEndDates()->Elt(shortFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);
		
		std::vector<double> rateLongFixPayDfs,rateLongFixPayTimes,rateLongInterestTerms;
		double rateLongFloatStartTime,rateLongFloatEndTime,rateLongDfStart,rateLongDfEnd;
		if(underlyingType == "SO/CMS" || underlyingType == "SO/SO" )
		{
			double rateLongStart = asOf + startTime3;
			double rateLongEnd   = asOf + endTime3;
			ARM_DateStrip rateLongSched(rateLongStart,rateLongEnd,stdFixFreq,stdFixDayCount,resetCalendar,fwdRule,intRule,stubRule,resetGap,stdFixFreq);

			rateLongFixPayDfs.resize(rateLongSched.GetFlowEndDates()->size());
			rateLongFixPayTimes.resize(rateLongSched.GetFlowEndDates()->size());
			rateLongInterestTerms = *rateLongSched.GetInterestTerms();
			for(i=0;i<rateLongFixPayDfs.size();i++)
			{
				rateLongFixPayDfs[i] = GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
				rateLongFixPayTimes[i] = rateLongSched.GetFlowEndDates()->Elt(i) - asOf;
			}
			rateLongFloatStartTime	= rateLongSched.GetFlowStartDates()->Elt(0) - asOf;
			rateLongFloatEndTime	= rateLongSched.GetFlowEndDates()->Elt(rateLongFixPayDfs.size()-1) - asOf;
			rateLongDfStart			= GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
			rateLongDfEnd			= GetZeroCurve()->DiscountPrice((rateLongSched.GetFlowEndDates()->Elt(rateLongFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);
		}

		std::vector<double> rateShortFixPayDfs,rateShortFixPayTimes,rateShortInterestTerms;
		double rateShortFloatStartTime,rateShortFloatEndTime,rateShortDfStart,rateShortDfEnd;
		if(underlyingType == "SO/SO")
		{
			double rateShortStart = asOf + startTime4;
			double rateShortEnd   = asOf + endTime4;
			ARM_DateStrip rateShortSched(rateShortStart,rateShortEnd,stdFixFreq,stdFixDayCount,resetCalendar,fwdRule,intRule,stubRule,resetGap,stdFixFreq);

			rateShortFixPayDfs.resize(rateShortSched.GetFlowEndDates()->size());
			rateShortFixPayTimes.resize(rateShortSched.GetFlowEndDates()->size());
			rateShortInterestTerms = *rateShortSched.GetInterestTerms();
			for(i=0;i<rateShortFixPayDfs.size();i++)
			{
				rateShortFixPayDfs[i] = GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowEndDates()->Elt(i)-asOf)/K_YEAR_LEN);
				rateShortFixPayTimes[i] = rateShortSched.GetFlowEndDates()->Elt(i) - asOf;
			}
			rateShortFloatStartTime	= rateShortSched.GetFlowStartDates()->Elt(0) - asOf;
			rateShortFloatEndTime	= rateShortSched.GetFlowEndDates()->Elt(rateShortFixPayDfs.size()-1) - asOf;
			rateShortDfStart			= GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowStartDates()->Elt(0)-asOf)/K_YEAR_LEN);
			rateShortDfEnd			= GetZeroCurve()->DiscountPrice((rateShortSched.GetFlowEndDates()->Elt(rateShortFixPayDfs.size()-1)-asOf)/K_YEAR_LEN);
		}
		
		ARM_ModelParamsHW1F* params1F = dynamic_cast<ARM_ModelParamsHW1F*>((ARM_ModelParams*)GetModelParams());
		ARM_ModelParamsHW2F* params2F = dynamic_cast<ARM_ModelParamsHW2F*>((ARM_ModelParams*)GetModelParams());
		
		/// little check
		if ( (!params1F && !params2F) || (params1F && params2F) )
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "UnderlyingCovariance : invalid model param");


		/// HW 1F case
		if (params2F == NULL)
		{
			double lambda = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);

			double longVolFactor, shortVolFactor;
		
			longVolFactor   = ARM_HullWhite::SwapRateVolFactor (lambda, 
																toTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longDfStart, 
																longDfEnd, 
																longFixPayTimes, 
																*longSched.GetInterestTerms(), 
																longFixPayDfs);

			shortVolFactor  = ARM_HullWhite::SwapRateVolFactor (lambda, 
																toTime, 
																shortFloatStartTime, 
																shortFloatEndTime, 
																shortDfStart, 
																shortDfEnd, 
																shortFixPayTimes, 
																*shortSched.GetInterestTerms(), 
																shortFixPayDfs);

			if(underlyingType == "SO/CMS" || underlyingType == "SO/SO" )
			{
				/// CMS Spread = LongCMS - ShortCMS
				longVolFactor -= shortVolFactor;

				shortVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
																		toTime, 
																		rateLongFloatStartTime, 
																		rateLongFloatEndTime, 
																		rateLongDfStart, 
																		rateLongDfEnd, 
																		rateLongFixPayTimes, 
																		rateLongInterestTerms, 
																		rateLongFixPayDfs);
			}
			if(underlyingType == "SO/SO" )
			{
				/// CMS Spread = LongCMS - ShortCMS
				shortVolFactor -= ARM_HullWhite::SwapRateVolFactor (lambda, 
																		toTime, 
																		rateShortFloatStartTime, 
																		rateShortFloatEndTime, 
																		rateShortDfStart, 
																		rateShortDfEnd, 
																		rateShortFixPayTimes, 
																		rateShortInterestTerms, 
																		rateShortFixPayDfs);
			}


			double phi = params1F->StateLocalVariance(fromTime, toTime, toTime);

			double cov = longVolFactor * shortVolFactor * phi;

			return cov;
		}
		/// HW 2F case
		else
		{
			double lambda1 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
			double lambda2 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversionSpread).GetValueAtPoint(0);
			lambda2 += lambda1;

			double longVolFactor1, longVolFactor2;
			double shortVolFactor1, shortVolFactor2;

			longVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																toTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longDfStart, 
																longDfEnd, 
																longFixPayTimes, 
																*longSched.GetInterestTerms(), 
																longFixPayDfs);
		
			longVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																	toTime, 
																	longFloatStartTime, 
																	longFloatEndTime, 
																	longDfStart, 
																	longDfEnd, 
																	longFixPayTimes, 
																	*longSched.GetInterestTerms(), 
																	longFixPayDfs);

			shortVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																	toTime, 
																	shortFloatStartTime, 
																	shortFloatEndTime, 
																	shortDfStart, 
																	shortDfEnd, 
																	shortFixPayTimes, 
																	*shortSched.GetInterestTerms(), 
																	shortFixPayDfs);

			shortVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																	toTime, 
																	shortFloatStartTime, 
																	shortFloatEndTime, 
																	shortDfStart, 
																	shortDfEnd, 
																	shortFixPayTimes, 
																	*shortSched.GetInterestTerms(), 
																	shortFixPayDfs);
			
			if(underlyingType == "SO/CMS" || underlyingType == "SO/SO")
			{
				/// CMS Spread = LongCMS - ShortCMS
				longVolFactor1 -= shortVolFactor1;
				longVolFactor2 -= shortVolFactor2;

				shortVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
																		toTime, 
																		rateLongFloatStartTime, 
																		rateLongFloatEndTime, 
																		rateLongDfStart, 
																		rateLongDfEnd, 
																		rateLongFixPayTimes, 
																		rateLongInterestTerms, 
																		rateLongFixPayDfs);

				shortVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
																		toTime, 
																		rateLongFloatStartTime, 
																		rateLongFloatEndTime, 
																		rateLongDfStart, 
																		rateLongDfEnd, 
																		rateLongFixPayTimes, 
																		rateLongInterestTerms, 
																		rateLongFixPayDfs);
			}
			if(underlyingType == "SO/SO")
			{
				shortVolFactor1 -= ARM_HullWhite::SwapRateVolFactor (lambda1, 
																		toTime, 
																		rateShortFloatStartTime, 
																		rateShortFloatEndTime, 
																		rateShortDfStart, 
																		rateShortDfEnd, 
																		rateShortFixPayTimes, 
																		rateShortInterestTerms, 
																		rateShortFixPayDfs);

				shortVolFactor2 -= ARM_HullWhite::SwapRateVolFactor (lambda2, 
																		toTime, 
																		rateShortFloatStartTime, 
																		rateShortFloatEndTime, 
																		rateShortDfStart, 
																		rateShortDfEnd, 
																		rateShortFixPayTimes, 
																		rateShortInterestTerms, 
																		rateShortFixPayDfs);
			}

			ARM_GP_TriangularMatrix* phi = params2F->StateLocalVariance(fromTime, toTime);
			double phi11 = (*phi)(0, 0);
			double phi22 = (*phi)(1, 1);
			double phi12 = (*phi)(1, 0);
			delete phi;

			double cov	=		 longVolFactor1 * shortVolFactor1 * phi11
							+	 longVolFactor2 * shortVolFactor2 * phi22
							+	(longVolFactor1 * shortVolFactor2 + longVolFactor2 * shortVolFactor1) * phi12 ;

		
			return cov;

		}
		
	}
	else
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "UnderlyingCovariance : Only ZC and CMS underlyingType are available for  HW Model");
		return 0.0;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_HullWhite::VanillaSpreadOptionLet(const string& curveName,
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
	/// Need to extend all Local_Normal_Model_Calibration_Helper functions
	/// to deal with state dependency
	if ( evalTime > K_NEW_DOUBLE_TOL )
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                    "VanillaSpreadOptionLet only works for evalTime = 0");
	}

	/// compute swap schedule since they have not been computed yet !
	ARM_Currency* ccy = GetCurrency(GetModelName());
	double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
	char fixCalendar[100];
    ccy->CalcFixPayCal(fixCalendar);
	int  fixFreq	 = ccy->GetFixedPayFreq();
	int  fixDayCount = ccy->GetFixedDayCount();
	
	ARM_SwapRatePtr longSwapRate = ARM_SwapRate::CreateSwapRate(asOf, 
																asOf+swapLongFloatStartTime, 
																asOf+swapLongFloatEndTime, 
																fixDayCount, 
																fixFreq, 
																fixCalendar);
	
	ARM_SwapRatePtr shortSwapRate = ARM_SwapRate::CreateSwapRate(asOf, 
																asOf+swapShortFloatStartTime, 
																asOf+swapShortFloatEndTime, 
																fixDayCount, 
																fixFreq, 
																fixCalendar);

	
	Local_Normal_Model_Calibration_Helper helper ((ARM_PricingModel*)this);

	double strike = strikes[0];

		
	double cmsLong = helper.CmsRateFromNumericalModel (	resetTime, 
														payTime, 
														longSwapRate->floatStartTime, 
														longSwapRate->floatEndTime, 
														longSwapRate->fixPayTimes, 
														longSwapRate->fixPayPeriods);

	double cmsShort = helper.CmsRateFromNumericalModel (resetTime, 
														payTime, 
														shortSwapRate->floatStartTime, 
														shortSwapRate->floatEndTime, 
														shortSwapRate->fixPayTimes, 
														shortSwapRate->fixPayPeriods);


	double spreadVol = helper.SpreadNormalVolatilityFromNumericalModel(	resetTime, 
																		payTime, 
																		coeffLong, 
																		coeffShort, 
																		strike, 
																		longSwapRate->floatStartTime, 
																		longSwapRate->floatEndTime, 
																		longSwapRate->fixPayTimes, 
																		longSwapRate->fixPayPeriods, 
																		shortSwapRate->floatStartTime, 
																		shortSwapRate->floatEndTime, 
																		shortSwapRate->fixPayTimes, 
																		shortSwapRate->fixPayPeriods);

	
	
	double mat = resetTime / K_YEAR_LEN;
	double fwd = coeffLong * cmsLong - coeffShort * cmsShort;

	double optionLet  = VanillaOption_N(fwd, spreadVol, strike, mat, callPut);
	double payDf	  = GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	double price      = notional * payPeriod * payDf * optionLet;
	
	size_t payoffSize = (states != ARM_PricingStatesPtr(NULL)) ? states->size() : 1 ;
	
	return ARM_VectorPtr( new std::vector<double>(payoffSize, price) );
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routines: VolatilitiesAndCorrelationTimesSteps
///	Returns : void
///	Action  : VolatilitiesAndCorrelationTimesSteps for PDE
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_HullWhite::VolatilitiesAndCorrelationTimesSteps() const
{
	const ARM_ModelParams* modelParams = GetModelParams();
	const ARM_ModelParamsHW* modelParamsHW = dynamic_cast<const ARM_ModelParamsHW*> (modelParams);

	if( !modelParamsHW )
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"HW Model params Expected !");

	return modelParamsHW->ModelParamsTimeSteps();

}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routine : SwapRateVolFactor (static method)
///	Returns : double
///	Action  : factor in front of sigma(t)*exp(-lambda*(Te-t)) in swap rate normal dynamics
////////////////////////////////////////////////////
double ARM_HullWhite::SwapRateVolFactor (   double lambda,
											double expiryTime,
											double swapFloatStartTime,	// adjusted ...
											double swapFloatEndTime,	// adjusted ...
											double dfStart,
											double dfEnd,
											const std::vector<double>& swapFixPayTimes,
											const std::vector<double>& swapFixPayPeriods,
											const std::vector<double>& swapFixPayDfs,
											// optional 
											double* swapRate, 
											double payTime,
											double* numeraireVolFactor) 

{
	size_t i, size = swapFixPayTimes.size();
	double Tstart  = swapFloatStartTime / K_YEAR_LEN;
	double Tend    = swapFloatEndTime / K_YEAR_LEN;
	double Te	   = expiryTime / K_YEAR_LEN;
		
	std::vector<double> T   (size);

	for (i=0; i<size; i++)
		T[i]   = swapFixPayTimes[i]/K_YEAR_LEN ;
	
	double level (0.0);
	for (i=0; i<size; i++)
		level += swapFixPayPeriods[i] * swapFixPayDfs[i] ;

	double _swapRate = (dfStart - dfEnd) / level;
	if (swapRate)
		*swapRate = _swapRate;

	double swapVolFactor (0);
		
	if ( fabs(lambda)<=K_NEW_DOUBLE_TOL )
		lambda = 1e-7; // pas très fair play mais bon
	
	for (i=0; i<size; i++)
		swapVolFactor -= swapFixPayPeriods[i] * swapFixPayDfs[i] * exp(-lambda * T[i]);

	swapVolFactor /= level;

	// Optional result:
	// computes factor in front of sigma(t)*exp(-lambda*(Te-t)) 
	// in the lognormal dynamics of B(t,Tp) / LVL(t)
	if (numeraireVolFactor)
	{	
		double Tp = payTime / K_YEAR_LEN;
		*numeraireVolFactor = swapVolFactor;
		*numeraireVolFactor += exp(-lambda * Tp);
		*numeraireVolFactor *= exp(lambda * Te) / lambda; // not multiplied by swaprate because lognormal vol
	}
	
	swapVolFactor += (exp(-lambda*Tstart) * dfStart - exp(-lambda*Tend) * dfEnd) / (dfStart - dfEnd);
	swapVolFactor *= _swapRate * exp(lambda * Te) / lambda;
	
	return swapVolFactor;
}
							



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

