/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Q1F.cpp
 *
 *  \brief Q model 1 factor
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date September 2004
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/Q1F.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/interpolatorvector.h"

/// gpinfra
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/nummethod.h"

/// gpmodels
#include "gpmodels/QModelAnalytics.h"
#include "gpmodels/ModelParamsQ1F.h"

/// gpnummethods
#include "gpnummethods/tree1D.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/meanrevertingsampler.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"


///  define for code clarity
#if defined(FIRST_STATE_VARIABLE)
	#undef	FIRST_STATE_VARIABLE
#endif
#define FIRST_STATE_VARIABLE 0

CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModel1F::ARM_QModel1F(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQ1F& params, bool degenerateInHW) 
:	ARM_QModelBaseIR(zc,params,degenerateInHW)
{}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModel1F::ARM_QModel1F(const ARM_QModel1F& rhs)
:	ARM_QModelBaseIR(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QModel1F& ARM_QModel1F::operator=(const ARM_QModel1F& rhs)
{
	if(this != &rhs)
		ARM_QModelBaseIR::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_QModel1F::Clone() const
{
	return new ARM_QModel1F(*this);
}



////////////////////////////////////////////////////
///	Class   : ~ARM_QModel1F
///	Routines: void
///	Returns :
///	Action  : Destructor
////////////////////////////////////////////////////
ARM_QModel1F::~ARM_QModel1F()
{}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Approximated closed form formula for DF
///          with recalibration using Arrow-Debreu prices
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    // CurveName is not used (only for multi-currencies model)

    ARM_GP_VectorPtr null(NULL);
    bool asOfPrice = evalTime <= K_NEW_DOUBLE_TOL || states == ARM_PricingStatesPtr(NULL);

    const ARM_QModelBaseIR::ARM_QDfTarget* qDfTarget=GetTargetPayoff(maturityTime);
    bool isValidPayoff = qDfTarget != NULL;

    double targetPayoff;
    ARM_GP_VectorPtr payoffFactor;
    if(isValidPayoff && !asOfPrice)
    {
        /// Use a previously set target payoff
        targetPayoff = qDfTarget->GetTargetPayoff();
        payoffFactor = qDfTarget->GetPayoffFactor();
    }
    else
        /// No special target payoff then use Zc on current yield curve
        targetPayoff = GetZeroCurve()->DiscountPrice(maturityTime/K_YEAR_LEN);

	size_t nbStates = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
    size_t i,modelNb = GetModelNb();

	if(maturityTime - K_NEW_DOUBLE_TOL <= evalTime)
	    /// Maturity overcome then return 1 by default
		return ARM_VectorPtr( new std::vector<double>(nbStates,1.0) );

	else if(asOfPrice)
	    /// AsOf Df
		return ARM_VectorPtr( new std::vector<double>(nbStates,targetPayoff) );

    else if(!isValidPayoff)
    {
        /// Closed form formula only valid for a H&W degenerated Q Model...
        ///... and Cash numeraire at the moment

        if(!IsDegenerateInHW())
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Closed form for DF only available for H&W like Q model" );

		if( GetNumeraire()->GetType() != ARM_Numeraire::Cash &&
            GetNumeraire()->GetType() != ARM_Numeraire::RollingEvent &&
            GetNumeraire()->GetType() != ARM_Numeraire::RollingCash && 
			GetNumeraire()->GetType() != ARM_Numeraire::TerminalZc)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Closed form for DF only available for Cash or RollingEvent or TerminalZc numeraire" );

        ARM_VectorPtr values;
        if(evalTime < maturityTime)
        {
            double zcT = targetPayoff;
            double zct = GetZeroCurve()->DiscountPrice(evalTime/K_YEAR_LEN);
            double stateFact = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(evalTime,maturityTime);
            double zcVar;
				
			if ((GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc) || 
				(GetNumeraire()->GetType() == ARM_Numeraire::RollingEvent) ||
				(GetNumeraire()->GetType() == ARM_Numeraire::RollingPayment))
			{
				double numMatTime = GetNumeraire()->GetMaturity();
				double drift	  = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(evalTime,numMatTime);
				zcVar= (stateFact*stateFact-2.0*drift*stateFact)*((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(0.0,evalTime,evalTime);
			}
			else
			{
				zcVar=((const ARM_ModelParamsHW1F* const) GetModelParams())->ZcVarianceSpread(evalTime,evalTime,maturityTime,evalTime);
			}

            values = ARM_VectorPtr(new std::vector<double>(nbStates));
            for(i=0;i<nbStates;++i)
                (*values)[i]=zcT/zct*exp(-0.5*zcVar-stateFact * states->GetModelState(i,modelNb));
        }
        else
            values = ARM_VectorPtr(new std::vector<double>(nbStates,1.0));

        return values;
    }
    else
    {
	    /// Approximated closed form formula & Arrow-Debreu calibration
		int timeIdx		= GetNumMethod()->GetLastTimeIdx(); /// Set by the numerical method <=> evalTime
		double f0		= (*GetPrecomputedFwds())[timeIdx];

		ARM_GP_VectorPtr arrowDebreuPrices = GetNumMethod()->GetArrowDebreuPrices( timeIdx, *this );

        if(arrowDebreuPrices == null)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Arrow-Debreu prices are missing for DF recalibration" );

		double expTerm, sum = 0.0;
		ARM_VectorPtr values(new std::vector<double>(nbStates));

/****
        /// Simple mapping version
		double betatT   = ((const ARM_ModelParamsQ1F* const) GetModelParams())->BetatT(evalTime,maturityTime);
	   	double q0		= GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).GetValueAtPoint(0);

		for(i=0;i<nbStates;++i)
		{
			expTerm = exp(-betatT*MappingFunction(states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE),f0,q0));
			sum	   += (*ArrowDebreuPrices)[i]*expTerm;
			(*values)[i] = expTerm;
        }
****/

/****/
        /// Zero-coupon vol approximation : forward rate in Zc vol is
        /// set to its spot value
        evalTime /= K_YEAR_LEN;
        maturityTime /= K_YEAR_LEN;
        if(isValidPayoff)
        {
            size_t nbPayoffFactors = payoffFactor->size();

#if defined( __GP_STRICT_VALIDATION )
            if(nbPayoffFactors != nbStates && nbPayoffFactors != 1)
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " payoff factor not correctly set" );
#endif
            /// Use the external factor to compute calibration payoff
		    for(i=0;i<nbStates;++i)
		    {
			    expTerm = exp(-DecayFactor(evalTime,maturityTime,evalTime,GetModelParams()) 
                    * states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE));
                sum	   += (*arrowDebreuPrices)[i] * (*payoffFactor)[i<nbPayoffFactors ? i : 0] * expTerm;
			    (*values)[i] = expTerm;
		    }
        }
        else
        {
            /// Standard calibration
		    for(i=0;i<nbStates;++i)
		    {
			    expTerm = exp(-DecayFactor(evalTime,maturityTime,evalTime,GetModelParams()) 
                    * states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE));
			    sum	   += (*arrowDebreuPrices)[i] * expTerm;
			    (*values)[i] = expTerm;
		    }
        }
/****/

/****
        /// Enhanced zero-coupon vol approximation : forward rate in
        /// Zc vol is set to its spot expectation
        evalTime /= K_YEAR_LEN;
        maturityTime /= K_YEAR_LEN;
		for(i=0;i<nbStates;++i)
		{
			expTerm = exp(-SuperDecayFactor(evalTime,maturityTime,evalTime,GetModelParams()) 
                * states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE));
			sum	   += (*ArrowDebreuPrices)[i]*expTerm;
			(*values)[i] = expTerm;
		}
****/

		for( i=0;i<nbStates;++i)
			(*values)[i] *= targetPayoff/sum;

		return values;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F::VanillaCaplet(
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
// FIXMEFRED: mig.vc8 (25/05/2007 15:48:07):cast
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));
		
	size_t i,nbStates = states->size();

#if defined( __GP_STRICT_VALIDATION )
	if( strikesPerState.size() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
#endif

    /// Conversion to year fraction
    double fwdResetYf = fwdResetTime/K_YEAR_LEN;
    double fwdStartYf = fwdStartTime/K_YEAR_LEN;
    double fwdEndYf = fwdEndTime/K_YEAR_LEN;

    /// Compute the flat libor rate
	ARM_VectorPtr zcEnd,libor(new std::vector<double>(nbStates));
	ARM_VectorPtr zcPay	= GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);
	ARM_VectorPtr zcStart = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdStartTime,states);
	if( fwdEndTime == payTime && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
		zcEnd = zcPay;
    else
		zcEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fwdEndTime,states);
	for(i=0;i<nbStates; ++i)
		(*libor)[i] = ((*zcStart)[i]/(*zcEnd)[i]-1.0)/fwdPeriod;

	ARM_VectorPtr caplet(new std::vector<double>(nbStates));

	if( fwdResetTime <= K_NEW_DOUBLE_TOL || evalTime >= fwdResetTime-K_NEW_DOUBLE_TOL )
	{
        /// Just compute the intrinsic value
		for(i=0;i<nbStates;++i)
			(*caplet)[i]=payNotional*period*CC_Max<double>(capFloor*(*libor)[i]-capFloor*strikesPerState[i],0) *(*zcPay)[i];
	}
	else
	{
        /// Compute vol FRA / vol Zc relation
        /// The Zc vol uses the "Mean field theory approximation"

        double decay = DecayFactor(fwdStartYf,fwdEndYf,fwdStartYf,GetModelParams());
        double mrs = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinate(0);
        double mrStdDev = static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(evalTime,fwdResetTime,fwdStartTime);
        mrStdDev = sqrt(mrStdDev);

        double payConvex=1.0;
        if(fwdEndTime - payTime < -5 ||  fwdEndTime - payTime > 5)
        {
            /// Compute a convexity adjustment using the forward Q volatility and
            /// the payment Zc one... to do !
        }

        if(IsDegenerateInHW())
        {
            /// ATM normal vol <=> Zc vol relation + normal price
            double norStdDev;
		    for(i=0;i<nbStates;++i)
            {
                norStdDev = ( (*libor)[i] + 1/fwdPeriod ) * decay * mrStdDev;
			    (*caplet)[i] = payNotional * period * (*zcPay)[i] *
                    VanillaOption_N( (*libor)[i], norStdDev, strikesPerState[i], 1.0, capFloor );
            }
/****
            /// Other method : price by an H&W like method with fwd Zc approx
            double amount=payNotional;
            if(period != fwdPeriod)
                amount *= period/fwdPeriod;
            double zcStdDev = decay * mrStdDev;
		    for(i=0;i<nbStates;++i)
	            (*caplet)[i] = amount * BlackSholes_Formula((*zcStart)[i]/(*zcEnd)[i],zcStdDev,(*zcPay)[i], 1.0+fwdPeriod*strikesPerState[i],capFloor);
****/
        }
        else
        {
            /// Q Vol <=> Zc vol relation
            double qParameter	= ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
            double qStdDev;
		    for(i=0;i<nbStates;++i)
            {
                qStdDev = ( 1+1/(fwdPeriod*(*libor)[i]) ) * decay * mrStdDev;
			    (*caplet)[i]=payNotional*period*QModelAnalytics::BSQFunction((*libor)[i], strikesPerState[i], qStdDev, 1.0, qParameter, (*zcPay)[i], capFloor, (*libor)[i] );
            }

/****
            /// Other method : ATM normal vol <=> Zc vol relation + conversion to ATM Q Vol + Q price
            double qParameter	= ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
            double norStdDev;
            double qStdDev;
		    for(i=0;i<nbStates;++i)
            {
                norStdDev = ( (*libor)[i] + 1/fwdPeriod ) * decay * mrStdDev;
                qStdDev = QModelAnalytics::NormalVolToQVol((*libor)[i],norStdDev,1.0,qParameter);
			    (*caplet)[i]=payNotional*period*QModelAnalytics::BSQFunction((*libor)[i], strikesPerState[i], qStdDev, 1.0, qParameter, (*zcPay)[i], capFloor, (*libor)[i] );
            }
****/
        }
    }

	return caplet;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////

ARM_VectorPtr ARM_QModel1F::VanillaSwaption(
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
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

	size_t i,nbStates = states->size();

#if defined( __GP_STRICT_VALIDATION )
	if( strikesPerState.rows() != nbStates )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " strike vector size != states size" );
	if( strikesPerState.cols() != fixPayTimes.size() )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": strike nb != fixed flow nb" );
#endif

    /// Use basis curve effect ?
    bool isBasis = !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor());

    /// Conversion to year fraction
    double swapResetYf	= swapResetTime/K_YEAR_LEN;
    double floatStartYf = floatStartTime/K_YEAR_LEN;
    double floatEndYf	= floatEndTime/K_YEAR_LEN;

    /// Compute swap rate, fixed Zc, annuity and decay factors
    std::vector<double> fixAnnuity(nbStates,0.0);
    size_t iFix,nbFixFlows=fixPayTimes.size();
    std::vector<double> decayFixLeg(nbStates,0.0);
    ARM_GP_VectorPtr zcStart( GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states) );
    vector< ARM_VectorPtr > zcFixFlow(nbFixFlows);
    ARM_VectorPtr zcEnd;
    double decay,flow;

    /// Create a single state to compute deterministic basis margins if necessary
    ARM_PricingStatesPtr singleState( new ARM_PricingStates(1,1) );
    ARM_GP_VectorPtr zcBasisFix,zcFix;
    double basisMargin,nextBasisMargin,nextFlowTime;
    std::vector<double> decayBasisLeg(nbStates,0.0);
    std::vector<double> basisFlow(nbFixFlows,0.0);

    /// Decay factor : fixed and pure basis leg part
    if(isBasis)
    {
        zcBasisFix      = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,singleState);
        zcFix           = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,singleState);
        basisMargin     = (*zcBasisFix)[0] / (*zcFix)[0];

        zcBasisFix      = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[0],singleState);
        zcFix           = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[0],singleState);
        nextBasisMargin = (*zcBasisFix)[0] / (*zcFix)[0];

        basisFlow[0]    = (nextBasisMargin/basisMargin - 1.0);

        basisMargin     = nextBasisMargin;
    }
    for(iFix=0;iFix<nbFixFlows;++iFix)
    {
        zcFixFlow[iFix] = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
        decay = DecayFactor(floatStartYf,fixPayTimes[iFix]/K_YEAR_LEN,floatStartYf,GetModelParams());

        if(isBasis && iFix+1 < nbFixFlows)
        {
            if(iFix+1 == nbFixFlows-1 && floatEndTime != fixPayTimes[nbFixFlows-1])
                nextFlowTime = floatEndTime;
            else
                nextFlowTime = fixPayTimes[iFix+1];

            zcBasisFix          = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,nextFlowTime,singleState);
            zcFix               = GetFixingFunctor()->DiscountFactor(curveName,evalTime,nextFlowTime,singleState);
            nextBasisMargin     = (*zcBasisFix)[0] / (*zcFix)[0];
            basisFlow[iFix+1]   = (nextBasisMargin/basisMargin - 1.0);
            basisMargin         = nextBasisMargin;
        }

        for(i=0;i<nbStates;++i)
        {
            flow = fixPayPeriods[iFix] * (*(zcFixFlow[iFix]))[i];
            fixAnnuity[i]       += flow;
            decayFixLeg[i]      += flow * decay;
            if(isBasis && iFix+1 < nbFixFlows)
                decayBasisLeg[i]    += basisFlow[iFix+1] * decay * (*(zcFixFlow[iFix]))[i];
        }
    }

    /// Decay factor : last nominal part
    double endDecay;
    if(floatEndTime != fixPayTimes[nbFixFlows-1])
    {
        zcEnd       = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatEndTime,states);
        endDecay    = DecayFactor(floatStartYf,floatEndYf,floatStartYf,GetModelParams());
    }
    else
    {
        zcEnd       = zcFixFlow[nbFixFlows-1];
        endDecay    = decay;
    }



    /// Approximation of the floating leg computation (no spread input neither frequency and basis)
    std::vector<double> floatLeg(nbStates,0.0);
    std::vector<double> basisAnnuity(nbStates,0.0);
    if(isBasis)
    {
        /// Use basis effect on a pseudo floating leg synchronised to the fixed one
        /// That begins as the double notional but basis margin flows must be added
        /// at each forward starting date (this sum is called the "basisAnnuity")
        for(i=0;i<nbStates;++i)
        {
            basisAnnuity[i] = basisFlow[0] * (*zcStart)[i];
            for(iFix=0;iFix+1<nbFixFlows;++iFix)
                basisAnnuity[i] += basisFlow[iFix+1] * (*(zcFixFlow[iFix]))[i];
            floatLeg[i] = (*zcStart)[i] - (*zcEnd)[i] + basisAnnuity[i];
        }
    }
    else
    {
        /// Standard double notional
        for(i=0;i<nbStates;++i)
            floatLeg[i] = (*zcStart)[i] - (*zcEnd)[i];
    }


    /// Compute option price
	ARM_VectorPtr swaption(new std::vector<double>(nbStates));
    std::vector<double> swapRate(nbStates,0.0);
    std::vector<double> swapDecay(nbStates,0.0);
    if(IsDegenerateInHW())
    {
        /// ATM normal vol <=> Zc vol relation + normal price
        for(i=0;i<nbStates;++i)
        {
            swapRate[i]     = floatLeg[i] / fixAnnuity[i];
            swapDecay[i]    = (endDecay * (*zcEnd)[i] + decayFixLeg[i]*swapRate[i] - decayBasisLeg[i]) / fixAnnuity[i];
        }

	    if( swapResetTime <= K_NEW_DOUBLE_TOL || evalTime >= swapResetTime-K_NEW_DOUBLE_TOL )
	    {
            /// Just compute the intrinsic value
		    for(i=0;i<nbStates;i++)
			    (*swaption)[i]=swapNotional*fixAnnuity[i]*CC_Max<double>(callPut*(swapRate[i]- strikesPerState(i,0)),0.0);
	    }
	    else
	    {
            double mrs = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinate(0);
            double mrStdDev = static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(evalTime,swapResetTime,floatStartTime);
            mrStdDev = sqrt(mrStdDev);

            double norStdDev;
		    for(i=0;i<nbStates;i++)
            {
                norStdDev = mrStdDev * swapDecay[i];
			    (*swaption)[i] = swapNotional * fixAnnuity[i] *
                    VanillaOption_N( swapRate[i], norStdDev, strikesPerState(i,0), 1.0, callPut );
            }
        }
    }
    else
    {
        /// Q Vol <=> Zc vol relation
        for(i=0;i<nbStates;++i)
        {
            swapRate[i]     = floatLeg[i] / fixAnnuity[i];
            swapDecay[i]    = ( (*zcEnd)[i] * endDecay - decayBasisLeg[i] ) / ((*zcStart)[i]-(*zcEnd)[i]+basisAnnuity[i])
                              + decayFixLeg[i]/fixAnnuity[i];
        }


	    if( swapResetTime <= K_NEW_DOUBLE_TOL || evalTime >= swapResetTime-K_NEW_DOUBLE_TOL )
	    {
            /// Just compute the intrinsic value
		    for(i=0;i<nbStates;i++)
			    (*swaption)[i]=swapNotional*fixAnnuity[i]*CC_Max<double>(callPut*(swapRate[i]- strikesPerState(i,0)),0.0);
	    }
	    else
	    {
            double mrs = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinate(0);
            double mrStdDev = static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(evalTime,swapResetTime,floatStartTime);
            mrStdDev = sqrt(mrStdDev);

            double qParameter = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);

            double qStdDev;
		    for(i=0;i<nbStates;i++)
            {
                qStdDev = mrStdDev * swapDecay[i];
			    (*swaption)[i] = swapNotional* 
				    QModelAnalytics::BSQFunction( swapRate[i], strikesPerState(i,0), qStdDev, 1.0, qParameter, fixAnnuity[i], callPut, swapRate[i] );
            }
        }
    }

	return swaption;    
}





////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_QModel1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    if(IsDegenerateInHW() )
        os << indent << "1F Q Model degenerated in H&W1F\n";
    else
        os << indent << "1F Q Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_QModel1F::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	if (!linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_QModel1F::VolatilitiesAndCorrelations : only linearVol case is implemented");

	/// to speed up, compute in one go the vol and vol derivatives
	std::vector<double> times = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetAbscisses();
	std::vector<double> values= ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetOrdinates();

/****
	std::vector<double> volsVec,d1VolsVec;
	VectorValuesAndDerivativesLinearMidPoints(times,values,timeSteps,volsVec,d1VolsVec);
	for(size_t i=0; i<d1VolsVec.size(); ++i )
		d1VolsVec[i] *= K_YEAR_LEN;
****/

    /// First step : shift values to be equivalent as a constant left interpolation
    size_t nbT=times.size();
    std::vector<double> newTimes;
    std::vector<double> newValues;
    if(times[0]>0.0)
    {
        newTimes.push_back(0.0);
        newValues.push_back(values[0]);
    }
    size_t i;
    for(i=0;i<nbT;++i)
    {
        newTimes.push_back(times[i]);
        newValues.push_back(values[i+1 < nbT ? i+1 : nbT-1]);
    }

    /// Second step : use special interpolation to compute vol values and 1st derivatives
    nbT = timeSteps.size();
	std::vector<double> volsVec(nbT),d1VolsVec(nbT);
    for(i=0;i<nbT;++i)
        volsVec[i] = FunctionSpecialInterpolation(timeSteps[i],newTimes,newValues);
    for(i=0;i<nbT;++i)
        d1VolsVec[i] = DerivativeSpecialInterpolation(timeSteps[i],timeSteps,volsVec) * K_YEAR_LEN;

	/// factor by line and zero correl because in one factor!
	vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,volsVec.size(),&volsVec[0]) );
	d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,d1VolsVec.size(),&d1VolsVec[0]));
	correls	= ARM_GP_MatrixPtr( NULL );
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_QModel1F::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1,-GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0) ) );
    for(size_t i=0;i<timeSteps.size()-1;++i)
        (*relativeDrifts)(i,FIRST_STATE_VARIABLE) *= (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;
	absoluteDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(), 1, 0.0 ) );
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_QModel1F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	ARM_PricingStatesPtr initStates( new ARM_PricingStates(bucketSize,1,0,1) );
    for(size_t i=0;i<bucketSize;++i)
        initStates->SetModelState(i,0,0.0);

    return initStates;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result;
	if( IsDegenerateInHW() )
	{
		size_t statesSize = states->size();
		const std::vector<double>& const timeSteps = GetNumMethod()->GetTimeSteps();

		if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING)
		{
			double f0			= (*GetPrecomputedFwds())[timeIdx];
   			double q0			= GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter).GetValueAtPoint(0);
			result				= ARM_GP_VectorPtr( new std::vector<double>(statesSize,0.0) );
			double startTime	= (*timeSteps)[timeIdx];
			size_t modelNb		= GetModelNb();

			double r;
			dt /= K_YEAR_LEN;
			for( size_t i=0; i<statesSize; ++i )
			{
				r = MappingFunction( states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE), f0, q0 );
				(*result)[i] = exp(-r*dt);
			}
		}
		else 
		{
			double t = (*timeSteps)[timeIdx];

			/// Compute an accurate value of f(0,t) with through the discount functor
			ARM_VectorPtr zct	= GetDiscountFunctor()->DiscountFactor("",0.0,t,states);
			double bm = (*zct)[0];

			ARM_VectorPtr zctdt	= GetDiscountFunctor()->DiscountFactor("",0.0,t+dt,states);
			double bp = (*zctdt)[0];

			double f0t = log(bm/bp);

			bool isDriftAdded = true; // to get the drift part related to volatility
			result = RiskNeutralDrift(timeIdx,states->GetModelStates(),isDriftAdded);
			dt /= K_YEAR_LEN;
			for( size_t i=0; i<statesSize; ++i )
				(*result)[i] = exp(-((*result)[i]) * dt - f0t);
		}
	}
	else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": LocalDiscounts() only holds for Q model degenerated to HW");

	return result;
}



////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_QModel1F::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 1, 0.0 ) );
	absoluteDrifts	= ARM_GP_MatrixPtr(NULL);

	for(size_t i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		(*relativeDrifts)(i,0) = ((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModel1F::ModelStateLocalVariances(
	const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	
	localVariances.resize((nbSteps-1)*(modelNb+1));
	size_t i;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalVariance(step,nextStep,nextStep));
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModel1F::NumMethodStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	ModelStateLocalVariances( timeSteps, localVariances );
}



////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModel1F::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex1	= (nbSteps-1)*modelNb;
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	variances.resize(nbSteps*(modelNb+1));

	/// fills the variance
    variances[offsetIndex2+0]=new ARM_GP_TriangularMatrix(1,0.0);
    for(size_t i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
		
		/// [i+1] => variance from 0 -> ti+1
        /// we can't sum up local variance !
        variances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsQ1F* const) GetModelParams())->StateLocalVariance(0.0,nextStep,nextStep));
        step=nextStep;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_QModel1F
///	Routine: SetNumMethod
///	Returns: void 
///	Action : Set the num method
////////////////////////////////////////////////////
void ARM_QModel1F::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	if( !GetFromMultiFactor() )
	{
		const ARM_TreeBase* treeBase = dynamic_cast<const ARM_TreeBase*>(&*numMethodPtr);
		if( treeBase)
		{
			if( dynamic_cast<const ARM_Tree1D*>(treeBase) )
			{
				if( !dynamic_cast<const ARM_MarkovianDriftSampler1D*>(treeBase->GetSampler()) &&
					!dynamic_cast<const ARM_MeanRevertingSampler1D*>(treeBase->GetSampler()))
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
					" Q model requires Markovian drift or Mean revering samplers !" );
			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 1D compatible" );
		}
	}
	ARM_PricingModel::SetNumMethod( numMethodPtr);
}



////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_QModel1F::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsQ1F* ModelParamsQ1F = dynamic_cast<const ARM_ModelParamsQ1F*>(&params);
	if( !ModelParamsQ1F )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsQ1F" );
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_QModel1F::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_QModel1F Model!");
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: IntegratedRiskNeutralDrift
///	Returns: a vector
///	Action : computes the integrated zero-coupon
///          risk neutral drift (only part which
///          depends on volatility parameters i.e.
///          Integ{s=evalTime->nextTime, [r(s)-f(0,s)]ds}
///          - diffusion part
///          Deterministic drift may be computed and
///          added if required (lattice may not need
///          this computation because of its internal calibration)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModel1F::IntegratedRiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded) const
{
	if( IsDegenerateInHW() )
	{
        double evalTime = GetNumMethod()->GetTimeStep(timeIdx);
        double nextTime = GetNumMethod()->GetTimeStep(timeIdx+1);
        double beta = ((const ARM_ModelParamsQ1F* const) GetModelParams())->BetatT(evalTime,nextTime);

        double varianceSpread=0.0;
        if(isDriftAdded)
            /// Compute the drift part depending of volatility data
            varianceSpread = 0.5 * ( ((const ARM_ModelParamsQ1F* const) GetModelParams())->ZcVarianceSpread(nextTime,evalTime,nextTime,evalTime) );

        size_t nbStates = numMethodStates->cols();
	    size_t modelNb	= GetModelNb();

	    ARM_VectorPtr result( new std::vector<double>(nbStates,0.0) );
	    for( size_t i=0; i<nbStates; ++i )
		    (*result)[i] = varianceSpread + beta * (*numMethodStates)(modelNb,i);

        return result;
    }
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": IntegratedRiskNeutralDrift() only holds for Q model degenerated to HW");
}


////////////////////////////////////////////////////
///	Class  : ARM_PricingModel
///	Routine: RiskNeutralDrift
///	Returns: a vector
///	Action : Compute the curve independant "short rate"
///          i.e. r(t)-f(0,t)
///          Deterministic drift may be computed and
///          added if required (lattice may not need
///          this computation because of its internal calibration)
////////////////////////////////////////////////////
ARM_GP_VectorPtr ARM_QModel1F::RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded ) const
{
	if( IsDegenerateInHW() )
	{
	    size_t nbStates = numMethodStates->cols();
	    size_t modelNb	= GetModelNb();

        double integPhi = 0.0;
        if(isDriftAdded)
        {
            /// Compute the drift part depending of volatility data
            double t = GetNumMethod()->GetTimeStep(timeIdx);
            const ARM_ModelParamsQ1F* const modelParams = dynamic_cast< const ARM_ModelParamsQ1F* const>(GetModelParams());
            integPhi = - ARM_ModelParamsHW1F::HW1FStateZcCovariance( modelParams, modelParams, 0.0, t, t, t );
        }

	    ARM_VectorPtr result( new std::vector<double>(nbStates,0.0) );
	    for( size_t i=0; i<nbStates; ++i )
		    (*result)[i] = integPhi + (*numMethodStates)(modelNb,i);

	    return result;
    }
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": RiskNeutralDrift() only holds for Q model degenerated to HW");
}


////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: IntegratedVolZcBond
///	Returns : ARM_VectorPtr
///	Action  : computes the risk neutral drift
////////////////////////////////////////////////////
double ARM_QModel1F::IntegratedVolZcBond( double a, double b, double c) const
{
    double decay	= DecayFactor(a/K_YEAR_LEN,b/K_YEAR_LEN,b/K_YEAR_LEN,GetModelParams());
    double mrStdDev = static_cast< const ARM_ModelParamsQ1F* const>(GetModelParams())->StateLocalVariance(a,b,c);
	mrStdDev = sqrt(mrStdDev);
	return decay * mrStdDev;
}

#undef FIRST_STATE_VARIABLE


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

