/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QNF.cpp
 *
 *  \brief Q model 1 factor
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/QNF.h"

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
#include "gpmodels/ModelParamsQNF.h"

/// gpnummethods
#include "gpnummethods/tree1d.h"
#include "gpnummethods/tree2d.h"
#include "gpnummethods/tree3d.h"
#include "gpnummethods/treend.h"
#include "gpnummethods/markoviandriftsampler.h"
#include "gpnummethods/meanrevertingsampler.h"


CC_BEGIN_NAMESPACE( ARM )

////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelNF::ARM_QModelNF(const ARM_ZeroCurvePtr& zc, const ARM_ModelParamsQNF& params, bool degenerateInHW ) 
:	ARM_QModelBaseIR(zc,params,degenerateInHW)
{}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelNF::ARM_QModelNF(const ARM_QModelNF& rhs)
:	ARM_QModelBaseIR(rhs)
{}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QModelNF& ARM_QModelNF::operator=(const ARM_QModelNF& rhs)
{
	if(this != &rhs)
		ARM_QModelBaseIR::operator=(rhs);
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_QModelNF::Clone() const
{
	return new ARM_QModelNF(*this);
}


////////////////////////////////////////////////////
///	Class   : ~ARM_QModelNF
///	Routines: void
///	Returns :
///	Action  : Destructor
////////////////////////////////////////////////////
ARM_QModelNF::~ARM_QModelNF()
{}




////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routines: toString
///	Returns :
///	Action  : stringify the object
////////////////////////////////////////////////////
string ARM_QModelNF::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << FactorCount() << "F Q Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: LocalDrifts
///	Returns: void
///	Action : local drifts
////////////////////////////////////////////////////
void ARM_QModelNF::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
    double step		= timeSteps[0],nextStep;
	size_t i,j;
	size_t factorNb = FactorCount();

	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(nbSteps-1,factorNb,0.0) );
	absoluteDrifts	= ARM_GP_MatrixPtr( NULL );

	for( i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		for( j=0; j<factorNb; ++j )
			(*relativeDrifts)(i,j) = ((const ARM_ModelParamsQNF* const) GetModelParams())->GetModelParamsPerDim(j)->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelNF::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	const std::vector<double>& const timeSteps = GetNumMethod()->GetTimeSteps();
	double f0		        = (*GetPrecomputedFwds())[timeIdx];
   	double q0               = ((const ARM_ModelParamsQNF* const) GetModelParams())->GetQParam().GetValueAtPoint(0);
	size_t statesSize		= states->size();
	std::vector<double>& result	= new std::vector<double>(statesSize,0.0);
	double startTime		= (*timeSteps)[timeIdx];
	size_t modelNb			= GetModelNb();
 
    /// Simple mapping version
	size_t factorNb = FactorCount();
    double factor,r;
    dt /= K_YEAR_LEN;
	size_t i,j;

	for( i=0; i<statesSize; ++i )
	{
		factor = 0.0;
		for( j=0; j<factorNb; ++j )
			factor += states->GetModelState(i,modelNb+j);

		r= MappingFunction( factor, f0, q0 );
		(*result)[i] = exp(-r*dt);
    }
	return ARM_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Approximated closed form formula for DF
///          with recalibration using Arrow-Debreu prices
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelNF::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    /// CurveName is not used (only for multi-currencies model)

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
	size_t modelNb	= GetModelNb();

	if(maturityTime - K_NEW_DOUBLE_TOL <= evalTime)
	    /// Maturity overcome then return 1 by default
		return ARM_VectorPtr( new std::vector<double>(nbStates,1.0) );

	else if(asOfPrice)
	    /// AsOf Df
		return ARM_VectorPtr( new std::vector<double>(nbStates,targetPayoff) );

	else
	{
	    /// Approximated closed form formula
		size_t factorNb	= FactorCount();
		ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );

		int timeIdx		= GetNumMethod()->GetLastTimeIdx(); /// Set by the numerical method <=> evalTime
		double f0		= (*GetPrecomputedFwds())[timeIdx];
		size_t i,j;
		ARM_GP_VectorPtr arrowDebreuPrices = GetNumMethod()->GetArrowDebreuPrices( timeIdx, *this );

        if(arrowDebreuPrices == null)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Arrow-Debreu prices are missing for DF recalibration" );

		double expTerm, sum = 0.0;
		ARM_VectorPtr values(new std::vector<double>(nbStates));

        /// Zero-coupon vol approximation : forward rate in Zc vol is
        /// set to its spot value
        evalTime     /= K_YEAR_LEN;
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
			    expTerm = 0.0;
			    for(j=0;j<factorNb; ++j )
				    expTerm += -DecayFactor(evalTime,maturityTime,evalTime,modelParams->GetModelParamsPerDim(j) )*states->GetModelState(i,modelNb+j);
			    expTerm = exp(expTerm);
			    sum	   += (*arrowDebreuPrices)[i] * (*payoffFactor)[i<nbPayoffFactors ? i : 0] * expTerm;
			    (*values)[i] = expTerm;
		    }
        }
        else
        {
            /// Standard calibration
		    for(i=0;i<nbStates;++i)
		    {
			    expTerm = 0.0;
			    for(j=0;j<factorNb; ++j )
				    expTerm += -DecayFactor(evalTime,maturityTime,evalTime,modelParams->GetModelParamsPerDim(j) )*states->GetModelState(i,modelNb+j);
			    expTerm = exp(expTerm);
			    sum	   += (*arrowDebreuPrices)[i] * expTerm;
			    (*values)[i] = expTerm;
		    }
        }

		for( i=0;i<nbStates;++i)
			(*values)[i] *= targetPayoff/sum;

		return values;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: VanillaCaplet
///	Returns: a vector of Caplet(t,L(R,S),K,S-E)
///	Action : Closed form formula for standard
///          caplet/floorlet (i.e. on libor resetting
///          in advance, paying in arrears)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelNF::VanillaCaplet(
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
		size_t iFactor,jFactor,nbFactors = FactorCount();
		std::vector<double> decay( nbFactors, 0.0 );
		std::vector<double> stdDevs( nbFactors, 0.0 );
		double stdDev = 0.0,var;
		ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );
		const ARM_GP_Matrix& correlMatrix = modelParams->GetCorrelMatrix();

		for( iFactor=0; iFactor<nbFactors; ++iFactor )
		{
			decay[iFactor]	= DecayFactor(fwdStartYf,fwdEndYf,fwdStartYf, modelParams->GetModelParamsPerDim(iFactor));
			var			    = modelParams->GetModelParamsPerDim(iFactor)->StateLocalVariance(evalTime,fwdResetTime,fwdStartTime);
			stdDev		    += decay[iFactor]*decay[iFactor]*var;
			for( jFactor=0; jFactor<iFactor; ++jFactor )
				stdDev += 2.0* correlMatrix(iFactor,jFactor)*decay[iFactor]*decay[jFactor]*ARM_ModelParamsHW1F::HW1FStateCovariance(modelParams->GetModelParamsPerDim(iFactor), modelParams->GetModelParamsPerDim(jFactor),evalTime, fwdResetTime, fwdStartTime );
		}
		stdDev = sqrt(stdDev);
		double payConvex= 1.0;
        if(fwdEndTime - payTime < -5 ||  fwdEndTime - payTime > 5)
        {
            /// Compute a convexity adjustment using the forward Q volatility and
            /// the payment Zc one... to do !
        }

        /// Method: Q Vol <=> Zc vol relation
        double qParameter = modelParams->GetQParam().GetCurve()->GetOrdinate(0);
        double qStdDev;
		for(i=0;i<nbStates;++i)
        {
            qStdDev = ( 1+1/(fwdPeriod*(*libor)[i]) )*stdDev;
			(*caplet)[i]=payNotional*period*QModelAnalytics::BSQFunction((*libor)[i], strikesPerState[i], qStdDev, 1.0, qParameter, (*zcPay)[i], capFloor, (*libor)[i] );
        }

    }

	return caplet;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : 1 factor like closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////

ARM_VectorPtr ARM_QModelNF::VanillaSwaption(
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

    /// Conversion to year fraction
    double swapResetYf	= swapResetTime/K_YEAR_LEN;
    double floatStartYf = floatStartTime/K_YEAR_LEN;
    double floatEndYf	= floatEndTime/K_YEAR_LEN;

    ARM_GP_VectorPtr swapRate( new std::vector<double>(nbStates,0.0) );
    ARM_GP_VectorPtr fixAnnuity( new std::vector<double>(nbStates,0.0) );
    size_t iFix,iLastFix,nbFixFlows=fixPayTimes.size();
    ARM_GP_VectorPtr zcStart( GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states) );
    vector< ARM_VectorPtr > zcFixFlow(nbFixFlows);
    ARM_VectorPtr zcEnd;
    double flow;
	size_t factorNb	= FactorCount(),iFactor,jFactor;
	ARM_GP_Matrix decay( factorNb, nbFixFlows+1, 0.0 );
	std::vector<double> var( factorNb, 0.0 );
	ARM_GP_Matrix covar( factorNb, factorNb, 0.0 );
	ARM_GP_Matrix stdDev( factorNb, nbStates, 0.0 );

	ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );
	const ARM_GP_Matrix& correlMatrix = modelParams->GetCorrelMatrix();

    /// State variances/covariances and decay factors
	for( iFactor=0; iFactor<factorNb; ++iFactor )
    {
		var[iFactor]	= modelParams->GetModelParamsPerDim(iFactor)->StateLocalVariance(evalTime,swapResetTime,floatStartTime);

		for( jFactor=0; jFactor<iFactor; ++jFactor )
            covar(iFactor,jFactor) = 2.0*correlMatrix(iFactor,jFactor)*ARM_ModelParamsHW1F::HW1FStateCovariance(modelParams->GetModelParamsPerDim(iFactor), modelParams->GetModelParamsPerDim(jFactor), evalTime, swapResetTime, floatStartTime );

        for(iFix=0;iFix<nbFixFlows;++iFix)
 			decay(iFactor,iFix) = DecayFactor(floatStartYf,fixPayTimes[iFix]/K_YEAR_LEN,floatStartYf,modelParams->GetModelParamsPerDim(iFactor));
    }

    /// Decay factors : fixed leg part
    iLastFix=nbFixFlows-1;
    for(iFix=0;iFix<nbFixFlows;++iFix)
    {
        zcFixFlow[iFix] = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);

        for(i=0;i<nbStates;++i)
        {
            flow = fixPayPeriods[iFix] * (*(zcFixFlow[iFix]))[i];
            (*fixAnnuity)[i] += flow;
	        for( iFactor=0; iFactor<factorNb; ++iFactor )
                stdDev(iFactor,i) += flow * decay(iFactor,iFix);
        }
    }

    /// Decay factors : last nominal part
    if(floatEndTime != fixPayTimes[nbFixFlows-1])
    {
        iLastFix = nbFixFlows;
        zcEnd = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,floatEndTime,states);
	    for( iFactor=0; iFactor<factorNb; ++iFactor )
 			decay(iFactor,nbFixFlows) = DecayFactor(floatStartYf,floatEndYf,floatStartYf,modelParams->GetModelParamsPerDim(iFactor));
    }
    else
        zcEnd = zcFixFlow[nbFixFlows-1];


    /// Approximation of the floating leg computation (no spread input neither frequency and basis)
    std::vector<double> floatLeg(nbStates,0.0);
    if(GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()))
    {
        /// Standard double notional
        for(i=0;i<nbStates;++i)
            floatLeg[i] = (*zcStart)[i] - (*zcEnd)[i];
    }
    else
    {
        /// Use basis effect on a pseudo floating leg synchronised to the fixed one
        ARM_VectorPtr zcFwdStart,zcFwdEnd;
        zcFwdStart = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
        for(iFix=0;iFix<nbFixFlows-1;++iFix)
        {
            zcFwdEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
            for(i=0;i<nbStates;++i)
                floatLeg[i] += ( (*zcFwdStart)[i]/(*zcFwdEnd)[i] - 1.0 ) * (*(zcFixFlow[iFix]))[i];
            zcFwdStart = zcFwdEnd;
        }
        if(floatEndTime != fixPayTimes[nbFixFlows-1])
        {
            zcFwdEnd = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatEndTime,states);
            for(i=0;i<nbStates;++i)
                floatLeg[i] += ( (*zcFwdStart)[i]/(*zcFwdEnd)[i] - 1.0 ) * (*zcEnd)[i];
        }
    }


    /// Compute option price
	ARM_VectorPtr swaption(new std::vector<double>(nbStates));


    /// Method #1 : Q Vol <=> Zc vol relation
    for(i=0;i<nbStates;++i)
    {
        flow = 1.0 / ( (*zcStart)[i]/(*zcEnd)[i] - 1.0 );
        (*swapRate)[i] = floatLeg[i] / (*fixAnnuity)[i];
	    for( iFactor=0; iFactor<factorNb; ++iFactor )
            stdDev(iFactor,i) = flow * decay(iFactor,iLastFix) + stdDev(iFactor,i)/(*fixAnnuity)[i];
    }


	if( swapResetTime <= K_NEW_DOUBLE_TOL || evalTime >= swapResetTime-K_NEW_DOUBLE_TOL )
	{
        /// Just compute the intrinsic value
		for(i=0;i<nbStates;i++)
			(*swaption)[i]=swapNotional*(*fixAnnuity)[i]*CC_Max<double>(callPut*((*swapRate)[i]- strikesPerState(i,0)),0.0);
	}
	else
	{
        double qParameter = modelParams->GetQParam().GetCurve()->GetOrdinate(0);
        double qStdDev,x;
		for(i=0;i<nbStates;i++)
        {
            qStdDev = 0.0;
	        for( iFactor=0; iFactor<factorNb; ++iFactor )
            {
                x=0.0;
	            for( jFactor=0; jFactor<iFactor; ++jFactor )
                    x += stdDev(jFactor,i)*covar(iFactor,jFactor);
                qStdDev += stdDev(iFactor,i) * (x + stdDev(iFactor,i)*var[iFactor]);
            }
            qStdDev = sqrt(qStdDev);

			(*swaption)[i] = swapNotional *
				QModelAnalytics::BSQFunction( (*swapRate)[i], strikesPerState(i,0), qStdDev, 1.0, qParameter, (*fixAnnuity)[i], callPut, (*swapRate)[i] );
        }
    }

	return swaption;    
}





////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_QModelNF::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	if (!linearVol)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_QModelNF::VolatilitiesAndCorrelations : only linearVol case is implemented");

	/// nb of factors
	size_t factorNb	= FactorCount();
	vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb,timeSteps.size(),0.0) );
	d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb,timeSteps.size(),0.0) );
	correls	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(factorNb*(factorNb-1)/2,timeSteps.size(),0.0) );
	ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );

	/// get correlation (the matrix is at this stage constant per factor!)
	const ARM_GP_Matrix& correlMatrix = modelParams->GetCorrelMatrix();

	/// to speed up, compute in one go the vol and vol derivatives
	size_t i,j,k,offseti=0,offset;
	for( i=0; i<factorNb; ++i )
	{
		std::vector<double> times = ((ARM_CurveModelParam&) modelParams->GetModelParamsPerDim(i)->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetAbscisses();
		std::vector<double> values= ((ARM_CurveModelParam&) modelParams->GetModelParamsPerDim(i)->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetOrdinates();
		std::vector<double> volsVec,d1VolsVec;
		VectorValuesAndDerivativesLinearMidPoints(times,values,timeSteps,volsVec,d1VolsVec);

#if defined( __GP_STRICT_VALIDATION )
		if( timeSteps.size()!= volsVec.size() )
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " timeSteps.size()!= volsVec.size()" );
#endif
		for( j=0; j<timeSteps.size(); ++j )
		{
			d1VolsVec[j]   *= K_YEAR_LEN;
			(*vols)(i,j)	= volsVec[j];
			(*d1Vols)(i,j)	= d1VolsVec[j];
			
			for( k=i+1,offset=offseti; k<factorNb; ++k,++offset )
				(*correls)(offset,j) = correlMatrix(i,k);
		}
        offseti += factorNb-1-i;
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_QModelNF::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t factorNb	= FactorCount();
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),factorNb,0.0) );
	absoluteDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),factorNb,0.0) );
	ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );
	size_t i,j;

    std::vector<double> yfSteps(timeSteps);
    yfSteps /= K_YEAR_LEN;

	for( i=0; i<timeSteps.size()-1;++i)
		for( j=0; j<factorNb; ++j )
			(*relativeDrifts)(i,j) = -modelParams->GetModelParamsPerDim(j)->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0) * (yfSteps[i+1]-yfSteps[i]);
}



////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: FirstPricingStates
///	Returns: ARM_PricingStatesPtr
///	Action  : create the first pricing state
////////////////////////////////////////////////////

ARM_PricingStatesPtr ARM_QModelNF::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = NF , nbPayoffs = 0)
	size_t factorNb	= FactorCount();
// FIXMEFRED: mig.vc8 (25/05/2007 15:44:23):cast
	return static_cast<ARM_PricingStatesPtr>(new ARM_PricingStates(bucketSize,factorNb,0,1));
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModelNF::ModelStateLocalVariances(
	const std::vector<double>& timeSteps,
	ARM_MatrixVector& localVariances ) const
{
	ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );
	const ARM_GP_Matrix& correlMatrix= modelParams->GetCorrelMatrix();

	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0], nextStep;
	size_t offsetIndex	= (nbSteps-1)*modelNb;
	size_t factorNb = FactorCount();

#if defined(__GP_STRICT_VALIDATION)
	if( localVariances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelNb+1));
	size_t i,j,k;
	double elem;
	
	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(factorNb,1.0);

		for( j=0; j<factorNb; ++j )
		{
			localVariances[offsetIndex+i]->Elt(j,j) = modelParams->GetModelParamsPerDim(j)->StateLocalVariance(step,nextStep,nextStep);
			for( k=0; k<j; ++k )
			{
				elem = 	correlMatrix(j,k)* ARM_ModelParamsHW1F::HW1FStateCovariance(modelParams->GetModelParamsPerDim(j), modelParams->GetModelParamsPerDim(k),step,nextStep,nextStep);
				localVariances[offsetIndex+i]->Elt(j,k) = elem;
				localVariances[offsetIndex+i]->Elt(k,j) = elem;
			}
		}
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModelNF::NumMethodStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
	ModelStateLocalVariances( timeSteps, localVariances );
}



////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_QModelNF::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances) const
{
	size_t factorNb	= FactorCount();
	ARM_ModelParamsQNF* modelParams = dynamic_cast<ARM_ModelParamsQNF*>( const_cast<ARM_ModelParams *> ( GetModelParams()) );
	const ARM_GP_Matrix& correlMatrix = modelParams->GetCorrelMatrix();

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
	variances[offsetIndex2+0]=new ARM_GP_TriangularMatrix(factorNb,0.0);
	size_t i,j,k;
	double elem;
	
	for(i=0;i<nbSteps-1;++i)
	{
		nextStep						= timeSteps[i+1];
		variances[offsetIndex2+i+1]		= new ARM_GP_TriangularMatrix(factorNb,0.0);

		for( j=0; j<factorNb; ++j )
		{
			variances[offsetIndex2+i+1]->Elt(j,j)    = modelParams->GetModelParamsPerDim(j)->StateLocalVariance(0,nextStep,nextStep);
			for( k=0; k<j; ++k )
			{
				elem= 	correlMatrix(j,k)*ARM_ModelParamsHW1F::HW1FStateCovariance(modelParams->GetModelParamsPerDim(j), modelParams->GetModelParamsPerDim(k),0,nextStep,nextStep);
	
				variances[offsetIndex2+i+1]->Elt(j,k)		= elem;
				variances[offsetIndex2+i+1]->Elt(k,j)		= elem;
			}
		}
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelNF
///	Routine: SetNumMethod
///	Returns: void 
///	Action : Set the num method
////////////////////////////////////////////////////
void ARM_QModelNF::SetNumMethod(const ARM_NumMethodPtr& numMethodPtr)
{
	const ARM_TreeBase* treeBase = dynamic_cast<const ARM_TreeBase*>(&*numMethodPtr);

	if(treeBase)
	{
		size_t factorNb	= FactorCount();
		switch( factorNb)
		{
		case 1:
			{
				if( dynamic_cast<const ARM_Tree2D*>(treeBase) )
				{
					if( !dynamic_cast<const ARM_MarkovianDriftSampler1D*>(treeBase->GetSampler()) &&
						!dynamic_cast<const ARM_MeanRevertingSampler1D*>(treeBase->GetSampler()))
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
						" Q model requires Markovian drift or Mean revering samplers !" );
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 1D compatible" );
			}
			break;
		case 2:
			{
				if( dynamic_cast<const ARM_Tree2D*>(treeBase) )
				{
					if( !dynamic_cast<const ARM_MarkovianDriftSampler2D*>(treeBase->GetSampler()) &&
					    !dynamic_cast<const ARM_MarkovianDriftSamplerND*>(treeBase->GetSampler()) &&
						!dynamic_cast<const ARM_MeanRevertingSamplerND*>(treeBase->GetSampler()))
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
						" Q model requires Markovian drift or Mean revering samplers !" );
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 2D compatible" );
			}
			break;
		case 3:
			{
				if( dynamic_cast<const ARM_Tree3D*>(treeBase) )
				{
					if( !dynamic_cast<const ARM_MarkovianDriftSampler3D*>(treeBase->GetSampler()) &&
					    !dynamic_cast<const ARM_MarkovianDriftSamplerND*>(treeBase->GetSampler()) &&
						!dynamic_cast<const ARM_MeanRevertingSamplerND*>(treeBase->GetSampler()))
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
						" Q model requires Markovian drift or Mean revering samplers !" );
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not 3D compatible" );
			}
			break;
		default:
			{
				if( dynamic_cast<const ARM_TreeND*>(treeBase) )
				{
					if( !dynamic_cast<const ARM_MarkovianDriftSamplerND*>(treeBase->GetSampler()) &&
						!dynamic_cast<const ARM_MeanRevertingSamplerND*>(treeBase->GetSampler()))
						ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
						" Q model requires Markovian drift or Mean revering samplers !" );
				}
				else
					ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " tree is not ND compatible" );
			}
			break;
		}
	}

	ARM_PricingModel::SetNumMethod( numMethodPtr);
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelNF
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_QModelNF::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ModelParamsQNF* ModelParamsQNF = dynamic_cast<const ARM_ModelParamsQNF*>(&params);
	if( !ModelParamsQNF )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ModelParamsQNF" );
	return true;
}

////////////////////////////////////////////////////
///	Class   : ARM_QModel1F
///	Routines: VanillaSpreadOption
///	Returns : void
///	Action  : No default implementation
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_QModelNF::VanillaSpreadOptionLet(const string& curveName,
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
		"VanillaSpreadOption : unimplemented function for ARM_QModelNF Model!");
}




CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/


