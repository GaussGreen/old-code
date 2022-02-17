/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW1F.cpp
 *  \brief Hull and white 1 factor model
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date September 2003
 */

/// this header comes first as it include some preprocessor constants
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/HW1F.h"
#include "gpmodels/ModelParamsHW1f.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/datestrip.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/irrate.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/bootstrapnd.h"
#include "gpcalib/vanillaarg.h"


/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"

/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )

const double quadMinx=-5.0;     /// -5.stdDev
const double quadMaxx=5.0;      ///  5.stdDev

////////////////////////////////////////////////////
///	Class  : HW1FVarianceFunction
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
HW1FVarianceFunction::HW1FVarianceFunction(const ARM_ModelParams* modelParams)
: itsModelParams(modelParams)
{}

HW1FVarianceFunction::~HW1FVarianceFunction()
{}

double HW1FVarianceFunction::operator() ( double x ) const
{
    return ((const ARM_ModelParamsHW1F*)itsModelParams)->StateLocalVariance(0.0,x,x);
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_HullWhite1F::ARM_HullWhite1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params,const ARM_BoolVector& soFormulaFlags )
:
	ARM_HullWhite(zc,params,soFormulaFlags)
{	
	CC_ARM_SETNAME(ARM_HW1FMODEL);
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HullWhite1F::ARM_HullWhite1F(const ARM_HullWhite1F& rhs)
: ARM_HullWhite(rhs) 
{
    // Copy class attributes if any
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_HullWhite1F::~ARM_HullWhite1F( )
{}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HullWhite1F& ARM_HullWhite1F::operator=(const ARM_HullWhite1F& rhs)
{
	if(this != &rhs)
	{
		ARM_HullWhite::operator=(rhs);		
        // Copy class attributes if any
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_HullWhite1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F Hull & White Model\n";
    os << indent << "---------------------\n";
	if(IsApproxSOFormula())
		os << indent << "\nApprox CMS Spread Option Formula\n\n";
	else
	{
		os << indent << "\nExact CMS Spread Option Formula, ";
		if(IsDeepITMSOFormula())
			os << indent << "Speed integration\n\n";
		else
			os << indent << "Accurate integration\n\n";
	}

    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_HullWhite1F::Clone() const
{
	return new ARM_HullWhite1F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HullWhite1F::PreProcessing(ARM_ModelFitter& modelFitter)
{ 
    /// we put a typeid to implement in a sense a double dispatcher...
    /// we did not use a dynamic to avoid throwing exception of type std::bad_cast
    if( typeid(modelFitter) == typeid(ARM_BootstrapND) )
    {
        ARM_CalibParamsHW1FExt CalibParamsHW1FExt(GetModelParams()->GetModelParams());
        CalibParamsHW1FExt.PreProcessing(modelFitter);
        SetModelParams(CalibParamsHW1FExt);
    }
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_HullWhite1F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
    /// we put a typeid to implement in a sense a double dispatcher...
    /// we did not use a dynamic to avoid throwing exception of type std::bad_cast
    if( typeid(modelFitter) == typeid(ARM_BootstrapND) )
    {
        ((ARM_CalibParamsHW1FExt*)GetModelParams())->PostProcessing(modelFitter,this);
        ARM_ModelParamsHW1FExt* ModelParamsHW1FExt = new ARM_ModelParamsHW1FExt(GetModelParams()->GetModelParams());
        SetModelParams(*(const ARM_ModelParamsHW1FExt *)(ModelParamsHW1FExt->Clone()));
    
        delete ModelParamsHW1FExt;
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_HullWhite1F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_HullWhite1F::ValidateModelParams(const ARM_ModelParams& params) const
{
    if(!params.DoesModelParamExist(ARM_ModelParamType::Volatility) || 
       !params.DoesModelParamExist(ARM_ModelParamType::MeanReversion))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "A vol and a MRS model parameters are required for 1F Hull & White model");
    }
	
	const ARM_CurveModelParam& volCurve = params.GetModelParam( ARM_ModelParamType::Volatility ).ToCurveModelParam();
    if( volCurve.GetCurve()->GetOrdinates().size() > 0 && 
        *volCurve.GetCurve() > ARM_HullWhite::VOL_NORMAL_MAX )
	{
		CC_Ostringstream os;
		os << "A normal volatility cannot be more than " << ARM_HullWhite::VOL_NORMAL_MAX;
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_HullWhite1F::IntegratedLocalDrifts(
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
		(*relativeDrifts)(i,0) = ((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalDrift(step,nextStep);
		step=nextStep;
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_HullWhite1F::ModelStateLocalVariances(
    const std::vector<double>& timeSteps,
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

	// There is no model local variance

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,1.0);
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_HullWhite1F::NumMethodStateLocalVariances(
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

	// All the variance is in the numerical method

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(step,nextStep,nextStep));
		step=nextStep;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_HullWhite1F::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& globalVariances) const
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
        globalVariances[offsetIndex2+i+1] = new ARM_GP_TriangularMatrix(1,((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(0.0,nextStep,nextStep));
        step=nextStep;
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_HullWhite1F::VarianceToTime(double var,double minTime,double maxTime) const
{
    HW1FVarianceFunction f(GetModelParams());
    T_BrentSolver<HW1FVarianceFunction> solver(f,var,0.0001);
	solver.setInitialGuess(0.0,minTime,maxTime);
    return solver.Solve();
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF. With
///          the Cash numeraire could be optimised
///          using phi(t) (see code). With Zc
///          numeraires, it uses a HJM like version.
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite1F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    // Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,zcT) );
    }

    double zct,zcVar,stateFact;
    int i,nbStates=states->size();
    if(    GetNumeraire() == ARM_NumerairePtr(NULL) 
		|| GetNumeraire()->GetType() == ARM_Numeraire::Cash
		|| GetNumeraire()->GetType() == ARM_Numeraire::RollingCash)
    {
        if(evalTime < maturityTime)
        {
            /// Could be optimized using phi(t) of dX(t)=[phi(t)-MRS.X(t)]dt + sigma(t).dWt
            stateFact = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(evalTime,maturityTime);
            zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
            zcVar=((const ARM_ModelParamsHW1F* const) GetModelParams())->ZcVarianceSpread(evalTime,evalTime,maturityTime,evalTime);
        }
        else
            return ARM_VectorPtr(new std::vector<double>(nbStates,1.0));
    }
    else
    {
        /// TerminalZc cases
        if(evalTime < maturityTime)
            stateFact = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(evalTime,maturityTime);
        else
            stateFact=0.0;
        zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);

		/// convexity change..
        double numMatTime = GetNumeraire()->GetMaturity();
        double drift	  = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(evalTime,numMatTime);
        zcVar= (stateFact*stateFact-2.0*drift*stateFact)*((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(0.0,evalTime,evalTime);
    }

	size_t modelNb = GetModelNb();
    ARM_VectorPtr values(new std::vector<double>(nbStates));
    for(i=0;i<nbStates;++i)
        (*values)[i]=zcT/zct*exp(-0.5*zcVar-stateFact * states->GetModelState(i,modelNb));

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite1F
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
///			 If the swaption is standard, this method
///          calls ARM_HullWhite::VanillaSwaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite1F::VanillaSwaption(
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
	if (isConstantNotional)
		return ARM_HullWhite::VanillaSwaption (	curveName, 
												evalTime, 
												swapResetTime, 
												fixNotional, 
												floatNotional, 
												floatStartTime, 
												floatEndTime, 
												floatResetTimes,
												floatStartTimes,
												floatEndTimes,
												floatIntTerms,
												fixPayTimes, 
												fixPayPeriods, 
												strikesPerState, 
												callPut, 
												states, 
												isConstantNotional, 
												isConstantSpread, 
												isConstantStrike);


	/// some validations...
	if (!isConstantSpread)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F: variable spread not supported)" );

	if (!isConstantStrike)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F: variable strike not supported)" );

	if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Swaption pricing not implemented for differents discount & fixing curves" );

	if (states->size() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  multi state not supported" );

	if (strikesPerState.rows() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  strike per state not supported" );

	if (strikesPerState.cols() != fixPayTimes.size()) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption:  strikesPerState: bad size" );
	
	if ( fabs(floatEndTime - fixPayTimes[fixPayTimes.size()-1])>0.001)
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption: float and fixed end dates are not matching" );

	if (floatNotional.size() != fixNotional.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW1F / Variable Notional Swaption: float an fix legs are required to be of same frequency." );
	
	
	/// compute coupons
	size_t i, j, size_cpn   = fixPayTimes.size();
	vector<double> coupons (size_cpn);

	for (i=0; i<size_cpn-1; i++)
		coupons[i] = floatNotional[i] + strikesPerState(0,i) * fixPayPeriods[i] * fixNotional[i] - floatNotional[i+1] ;
	
	coupons[size_cpn-1] = floatNotional[size_cpn-1] + fixPayPeriods[size_cpn-1] * strikesPerState(0,size_cpn-1) * fixNotional[size_cpn-1];

	/// bond option strike
	double strike = floatNotional[0];
			
	/// precomputations for numerical integral
	vector<double> dfs (size_cpn);
	vector<double> factors (size_cpn);
	double betai;
	double beta0 = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(swapResetTime, floatStartTime);
	double phi	 = ((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(0.0, swapResetTime, swapResetTime);
	double startdf = GetZeroCurve()->DiscountPrice(floatStartTime/K_YEAR_LEN);
	
	double EX   = -phi * beta0; // E^Q(T0)[X(Te)] = - phi(Te) * beta(Te,T0)
	double varX = phi;
	double invVarX = 1.0 / varX;
	double factor =  1.0 / sqrt(2.0 * ARM_NumericConstants::ARM_PI * varX) ;
		
	double	nstdev	= 6.0 ;
	int		nx		= 201;
	double  dX      = 2.0 * nstdev * sqrt(varX) / (nx - 1.0);
	double  Xmin    = EX -  nstdev * sqrt(varX) ;
	
	for (i=0; i<size_cpn; i++)
	{
		betai = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(swapResetTime, fixPayTimes[i]);
		dfs[i] = GetZeroCurve()->DiscountPrice (fixPayTimes[i]/K_YEAR_LEN) / startdf;
		dfs[i] *= exp (-0.5 * phi * (betai*betai - beta0*beta0) - (betai - beta0) * Xmin);
		factors[i] = exp(-(betai - beta0) * dX);
	}

	double expect (0.0);
	double X (Xmin);
	double argexp, coeff, payoff;
	double sgn = -(double)callPut; // payer swaption = put bond option

	for (j=0; j<nx; j++)
	{
		X += dX;
		coeff = (j==0||j==nx-1) ? 0.5 * dX : dX ;
		argexp = -0.5 * (X-EX) * (X-EX) * invVarX ;
	
		payoff = 0.0;
		for (i=0; i<size_cpn; i++)
		{
			dfs[i] *= factors[i];
			payoff += coupons[i] * dfs[i] ; 
		}

		payoff -= strike;
		payoff *= sgn;
		if (payoff<0.0) payoff = 0.0;
		
		expect += coeff * payoff * factor * exp (argexp);
	}
	
	double price = startdf * expect ;
	
	return ARM_VectorPtr (new std::vector<double>(1,price));
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: VanillaSpreadOption
///	Returns : ARM_VectorPtr
///	Action  : Price an option on a.LongCMS - b.ShortCMS where
///			  LongCMS and ShortCMS are two swap rates paid at time Te
///			  Useful if degenerated in CMS caplet with b=0
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite1F::VanillaSpreadOptionLet(
			const string& curveName,
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
    double var=((const ARM_ModelParamsHW1F* const) GetModelParams())->StateLocalVariance(evalTime,resetTime,resetTime);
	double stdDev = sqrt(var);

	double stdDevLimit = 1.0e-4*sqrt((resetTime-evalTime)/K_YEAR_LEN); /// volLimit=1bp

    if( evalTime < K_NEW_DOUBLE_TOL && (stdDev < stdDevLimit || IsApproxSOFormula()) )
    {
		/// Only as of evaluation is allowed (approx pricing is not upgraded
		/// to support state dependency)
        /// Standard normal spread approximation based evaluation
		return ARM_HullWhite::VanillaSpreadOptionLet(
			curveName,
			evalTime,
			callPut,
			startTime,
			endTime,
			resetTime,
			payTime,
			payPeriod,
			notional,
			coeffLong,
			coeffShort,
			strikes,
			swapLongFloatStartTime,
			swapLongFloatEndTime,
			swapLongFixPayTimes,
			swapLongFixPayPeriods,
			swapShortFloatStartTime,
			swapShortFloatEndTime,
			swapShortFixPayTimes,
			swapShortFixPayPeriods,
			states);
    }

	/// Compute an exact price by 1D numerical integation
	double longFloatStartTime = swapLongFloatStartTime;
	double longFloatEndTime = swapLongFloatEndTime;
	double shortFloatStartTime = swapShortFloatStartTime;
	double shortFloatEndTime = swapShortFloatEndTime;
	std::vector<double> *longFixPayTimes=NULL;
	std::vector<double> *longFixPayPeriods=NULL;
	std::vector<double> *shortFixPayTimes=NULL;
	std::vector<double> *shortFixPayPeriods=NULL;

	bool isLong		= coeffLong!=0.0;
	bool isShort	= coeffShort!=0.0;

	/// compute swap schedule since they have not been computed yet !
	ARM_Currency* ccy = GetCurrency(GetModelName());
	double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
	char fixCalendar[100];
	ccy->CalcFixPayCal(fixCalendar);
	int  fixFreq	 = ccy->GetFixedPayFreq();
	int  fixDayCount = ccy->GetFixedDayCount();

	if(isLong)
	{
		if(&swapLongFixPayTimes != NULL)
		{
			longFixPayTimes		= const_cast<ARM_GP_Vector*>(&swapLongFixPayTimes);
			longFixPayPeriods	= const_cast<ARM_GP_Vector*>(&swapLongFixPayPeriods);
		}
		else
		{
			ARM_SwapRatePtr longSwapRate(ARM_SwapRate::CreateSwapRate(asOf,
														asOf+swapLongFloatStartTime, 
														asOf+swapLongFloatEndTime, 
														fixDayCount, 
														fixFreq, 
														fixCalendar));

			longFloatStartTime	= longSwapRate->floatStartTime;
			longFloatEndTime	= longSwapRate->floatEndTime;
			longFixPayTimes		= & longSwapRate->fixPayTimes;
			longFixPayPeriods	= & longSwapRate->fixPayPeriods;
		}
	}
	
	if(isShort)
	{
		if(&swapShortFixPayTimes != NULL)
		{
			shortFixPayTimes	= const_cast<ARM_GP_Vector*>(&swapShortFixPayTimes);
			shortFixPayPeriods	= const_cast<ARM_GP_Vector*>(&swapShortFixPayPeriods);
		}
		else
		{
			ARM_SwapRatePtr shortSwapRate(ARM_SwapRate::CreateSwapRate(asOf,
														asOf+swapShortFloatStartTime, 
														asOf+swapShortFloatEndTime, 
														fixDayCount, 
														fixFreq, 
														fixCalendar));

			shortFloatStartTime	= shortSwapRate->floatStartTime;
			shortFloatEndTime	= shortSwapRate->floatEndTime;
			shortFixPayTimes	= & shortSwapRate->fixPayTimes;
			shortFixPayPeriods	= & shortSwapRate->fixPayPeriods;
		}
	}

	/// Beta(Te,Tpay)
	double betatTeTp = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(resetTime,payTime);

	size_t i;

	/// Get long standard swap datas
	size_t longNbFlows=0;
	double longBetatTeT0=0.0,longLogDriftT0=0.0,longStdDevT0=0.0;
	vector< double > longFlowPeriod;
	vector< double > longLogDrift;
	vector< double > longStdDev;
	ARM_VectorPtr longZcFloatStart(NULL);
	vector< ARM_VectorPtr > longZcFlowPay;
	if(isLong)
	{
		vector< double > longFlowTimes;
		ComputeFlowZc(curveName,evalTime,longFloatStartTime,longFloatEndTime,*longFixPayTimes,
			*longFixPayPeriods,states,longZcFloatStart,longZcFlowPay,longFlowTimes,longFlowPeriod);
		longNbFlows=longZcFlowPay.size();
		bool isLongSameEnd = (longNbFlows==swapLongFixPayTimes.size());

		/// Drift and diffusion factor of B(Te,TlongStart)
		longBetatTeT0 = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(resetTime,swapLongFloatStartTime);
		longBetatTeT0	-= betatTeTp;

		longStdDevT0 = longBetatTeT0 * stdDev;
		longLogDriftT0 = -0.5*longStdDevT0*longStdDevT0;

		/// Drift and diffusion factor of B(Te,TlongFlowi)
		double longBetatTeTi;
		longLogDrift.resize(longNbFlows);
		longStdDev.resize(longNbFlows);
		for(i=0;i<longNbFlows;++i)
		{
			longBetatTeTi = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(resetTime,longFlowTimes[i]);
			longBetatTeTi	-= betatTeTp;

			longStdDev[i] = longBetatTeTi * stdDev;
			longLogDrift[i] = -0.5*longStdDev[i]*longStdDev[i];
		}

		if(longNbFlows == longFixPayTimes->size())
		{
			/// Add nominal payback flow to simplify further loops
			longZcFlowPay.push_back(longZcFlowPay[longNbFlows-1]);
			longFlowPeriod.push_back(1.0);
			longStdDev.push_back(longStdDev[longNbFlows-1]);
			longLogDrift.push_back(longLogDrift[longNbFlows-1]);
		}
	}


	/// Get short standard swap datas
	size_t shortNbFlows=0;
	double shortLogDriftT0=0.0,shortStdDevT0=0.0;
	vector< double > shortFlowPeriod;
	vector< double > shortLogDrift;
	vector< double > shortStdDev;
	ARM_VectorPtr shortZcFloatStart(NULL);
	vector< ARM_VectorPtr > shortZcFlowPay;
	if(isShort)
	{
		vector< double > shortFlowTimes;
		ComputeFlowZc(curveName,evalTime,shortFloatStartTime,shortFloatEndTime,*shortFixPayTimes,
			*shortFixPayPeriods,states,shortZcFloatStart,shortZcFlowPay,shortFlowTimes,shortFlowPeriod);
		shortNbFlows=shortZcFlowPay.size();

		/// Drift and diffusion factor of B(Te,TshortStart)
		double shortBetatTeT0=longBetatTeT0;
		shortStdDevT0=longStdDevT0;
		shortLogDriftT0=longLogDriftT0;
		if(swapShortFloatStartTime != swapLongFloatStartTime)
		{
			shortBetatTeT0 = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(resetTime,swapShortFloatStartTime);
			shortBetatTeT0	-= betatTeTp;

			shortStdDevT0 = shortBetatTeT0 * stdDev;
			shortLogDriftT0 = -0.5*shortStdDevT0*shortStdDevT0;
		}

		/// Drift and diffusion factor of B(Te,TshortFlowi)
		double shortBetatTeTi;
		shortLogDrift.resize(shortNbFlows);
		shortStdDev.resize(shortNbFlows);
		for(i=0;i<shortNbFlows;++i)
		{
			shortBetatTeTi = ((const ARM_ModelParamsHW1F* const) GetModelParams())->BetatT(resetTime,shortFlowTimes[i]);
			shortBetatTeTi -= betatTeTp;

			shortStdDev[i] = shortBetatTeTi * stdDev;
			shortLogDrift[i] = -0.5*shortStdDev[i]*shortStdDev[i];
		}

		if(shortNbFlows == shortFixPayTimes->size())
		{
			/// Add nominal payback flow to simplify further loops
			shortZcFlowPay.push_back(shortZcFlowPay[shortNbFlows-1]);
			shortFlowPeriod.push_back(1.0);
			shortStdDev.push_back(shortStdDev[shortNbFlows-1]);
			shortLogDrift.push_back(shortLogDrift[shortNbFlows-1]);
		}
	}


	int nbStates=(isLong ? longZcFloatStart->size() : shortZcFloatStart->size());
	ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));

	double ratio = payPeriod*notional;
	ARM_VectorPtr zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

	size_t j,k;
	double vi,sLong,sShort;

	if(stdDev<=ARM_NumericConstants::ARM_TOLERENCE)
	{
		/// Case of intrinsic value
		for(k=0;k<nbStates;++k)
		{
			sLong=0;
			if(isLong)
			{
				for(i=0;i<longNbFlows;++i)
					sLong += longFlowPeriod[i] * (*(longZcFlowPay[i]))[k];
				sLong = ( (*longZcFloatStart)[k] - (*(longZcFlowPay[longNbFlows]))[k])/sLong;
			}

			sShort=0;
			if(isShort)
			{
				for(i=0;i<shortNbFlows;++i)
					sShort += shortFlowPeriod[i] * (*(shortZcFlowPay[i]))[k];
				sShort = ( (*shortZcFloatStart)[k] - (*(shortZcFlowPay[shortNbFlows]))[k])/sShort;
			}

			vi = callPut*(coeffLong*sLong - coeffShort*sShort - strikes[k]);
			(*values)[k] = (vi<0 ? 0.0 : vi) * (*zcPay)[k] * ratio;
		}
		return values;
	}

	/// Compute an exact price of the payoff = MAX(callPut*(a1.S1-a2.K),0)
	/// using a 1D numerical integration for NC(0,1) between [quadMinx,quadMaxx]
	/// ------------------------------------------------------------------------

	int nbRefPts=8,nbStrats=5;
	if(IsDeepITMSOFormula())
	{
		/// Rough 1D integration because the payoff is no more conditionnal
		/// Only convexity adjustment is really computed
		nbRefPts=7;
		nbStrats=3;
	}

	GaussLegendre_Coefficients GL(nbRefPts); /// becareful : never add both last arguments quadMinx && quadMaxx here !!
	GaussStratifiedHermiteLegendre_Coefficients GsHL(&GL,nbStrats,quadMinx,quadMaxx);

	double z,w;
	std::vector<double> longFlowStoc(longNbFlows+1,0.0),shortFlowStoc(shortNbFlows+1,0.0);
	double longStartStoc=0.0,shortStartStoc=0.0;
	std::vector<double> condValues(nbStates);
	for(i=0;i<GsHL.get_order();++i)
	{
		z = GsHL.get_point(i);
		w = GsHL.get_weight(i);

		if(isLong)
		{
			longStartStoc = exp(longLogDriftT0 + longStdDevT0*z);
			for(j=0;j<longNbFlows+1;++j)
				longFlowStoc[j] = longFlowPeriod[j] * exp(longLogDrift[j] + longStdDev[j]*z);
		}

		if(isShort)
		{
			shortStartStoc = exp(shortLogDriftT0 + shortStdDevT0*z);
			for(j=0;j<shortNbFlows+1;++j)
				shortFlowStoc[j] = shortFlowPeriod[j] * exp(shortLogDrift[j] + shortStdDev[j]*z);
		}

		for(k=0;k<nbStates;++k)
		{
			/// Only Zc ratios appears in swap rate so it
			/// not needed to divide all B(Te,Ti) by B(Te,Tp)
			sLong=0;
			if(isLong)
			{
				for(j=0;j<longNbFlows;++j)
					sLong += (*(longZcFlowPay[j]))[k] * longFlowStoc[j];
				sLong = ( (*longZcFloatStart)[k] * longStartStoc - (*(longZcFlowPay[longNbFlows]))[k] * longFlowStoc[longNbFlows] )/sLong;
			}

			sShort=0;
			if(isShort)
			{
				for(j=0;j<shortNbFlows;++j)
					sShort += (*(shortZcFlowPay[j]))[k] * shortFlowStoc[j];
				sShort = ( (*shortZcFloatStart)[k] * shortStartStoc - (*(shortZcFlowPay[shortNbFlows]))[k] * shortFlowStoc[shortNbFlows] )/sShort;
			}

			vi = callPut*(coeffLong*sLong - coeffShort*sShort - strikes[k]);
			(*values)[k] += (vi<0 ? 0.0 : vi) * w;
		}
	}

	for(k=0;k<nbStates;++k)
		(*values)[k] *= ((*zcPay)[k] * ratio);

    return values;
}


///////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HullWhite1F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	return new ARM_PricingStates(bucketSize,1,0,1);
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_HullWhite1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
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
///	Class   : ARM_HullWhite1F
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_HullWhite1F::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
    /// FIX FIX : only valid for cst MRS !!!
	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1,-GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0) ) );
    for(size_t i=0;i<timeSteps.size()-1;++i)
        (*relativeDrifts)(i,0) *= (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;

	absoluteDrifts= ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),1, 0.0 ) );
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_HullWhite1F::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
	/// used for markovian drift sampler
	/// vol is made piecewise linear
	/// derivatives are computed
	if (linearVol)
	{
		/// to speed up, compute in one go the vol and vol derivatives
		std::vector<double> times = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility)).GetCurve()->GetAbscisses();
		std::vector<double> values= ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility)).GetCurve()->GetOrdinates();

		/// Since volatilies are stepUpRight!!
		for( std::vector<double>::iterator iter = values.begin() ; iter != values.end()-1 ; ++iter)
			(*iter) = *(iter+1);

		std::vector<double> volsVec,d1VolsVec;
		VectorValuesAndDerivativesLinearMidPoints(times,values,timeSteps,volsVec,d1VolsVec);
		for(size_t i=0; i<d1VolsVec.size(); ++i )
			d1VolsVec[i] *= K_YEAR_LEN;

		/// factor by line and zero correl because in one factor!
		vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,volsVec.size(),&volsVec[0]) );
		d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1,d1VolsVec.size(),&d1VolsVec[0]));
		correls	= ARM_GP_MatrixPtr( NULL );
	}
	/// Instant. vols for PDE schemes
	else
	{
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
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite1F::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	ARM_GP_VectorPtr result;

	size_t statesSize = states->size();
	const std::vector<double>& const timeSteps = GetNumMethod()->GetTimeSteps();

	if(GetNumMethod()->GetPricingDirection() == ARM_NumMethod::GP_BCKWDLOOKING)
	{
		result				= ARM_GP_VectorPtr( new std::vector<double>(statesSize,0.0) );
		double startTime	= (*timeSteps)[timeIdx];
		size_t modelNb		= GetModelNb();

		double r;
		dt /= K_YEAR_LEN;
		for( size_t i=0; i<statesSize; ++i )
		{
			r = states->GetModelState(i,modelNb+0);
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
			(*result)[i] = exp(-(*result)[i] * dt);
	}
	
	return result;
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
ARM_GP_VectorPtr ARM_HullWhite1F::RiskNeutralDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates, bool isDriftAdded ) const
{
    size_t nbStates = numMethodStates->cols();
    size_t modelNb	= GetModelNb();

    double integPhi = 0.0;

	//dates needed
	double ti = GetNumMethod()->GetTimeStep(timeIdx);
	double tk = GetNumMethod()->GetTimeStep(timeIdx+1);
	double dt = (tk-ti)/K_YEAR_LEN;

	//Calculation of forward instaneous rate
	ARM_ZeroCurvePtr& ZcCurve=GetZeroCurve();
    double Bdom0ti = ZcCurve->DiscountPrice(ti/K_YEAR_LEN);
    double Bdom0tk = ZcCurve->DiscountPrice(tk/K_YEAR_LEN);

    double domFwdRate = log(Bdom0ti/Bdom0tk)/dt;

    if(isDriftAdded)
    {
    /// Compute the drift part depending of volatility data
        double t = GetNumMethod()->GetTimeStep(timeIdx);
        const ARM_ModelParamsHW1F* const modelParams = dynamic_cast< const ARM_ModelParamsHW1F* const>(GetModelParams());
        integPhi = - ARM_ModelParamsHW1F::HW1FStateZcCovariance( modelParams, modelParams, 0.0, t, t, t );
    }

	ARM_VectorPtr result( new std::vector<double>(nbStates,0.0) );
	for( size_t i=0; i<nbStates; ++i )
		(*result)[i] = domFwdRate + integPhi + (*numMethodStates)(modelNb,i);

	return result;
}

////////////////////////////////////////////////////
///	Class   : ARM_PricingModel
///	Routines: UpdatePDE3DCoeffs
///	Returns : void
///	Action  : Update the PDE 3D coefficients
////////////////////////////////////////////////////
void ARM_HullWhite1F::UpdatePDE3DCoeffs(
		size_t timeIdx,
		const ARM_GP_MatrixPtr& numMethodStates,
		ARM_GP_VectorPtr& qxx,
		ARM_GP_VectorPtr& qyy,
		ARM_GP_VectorPtr& qzz,
		ARM_GP_VectorPtr& qxy,
		ARM_GP_VectorPtr& qyz,
		ARM_GP_VectorPtr& qzx,
		ARM_GP_VectorPtr& px,
		ARM_GP_VectorPtr& py,
		ARM_GP_VectorPtr& pz,
		ARM_GP_VectorPtr& o
		) const
{
	size_t N = numMethodStates->cols();
	double time = GetNumMethod()->GetTimeStep(timeIdx);
	
	ARM_CurveModelParam& curve_vol = (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility);
	ARM_CurveModelParam& curve_mr = (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion);
	double minus_square_vol_t=-pow ( curve_vol.GetValue(time), 2);
	double mr_t=curve_mr.GetValue(time);

	int K;
	if(GetNumeraire()->GetType() == ARM_Numeraire::TerminalZc)
	{
		for (K = 0; K < N; ++K)
		{
			(*qxx)[K] = minus_square_vol_t;
			(*qyy)[K] = 0.0;
			(*qzz)[K] = 0.0;
			(*qxy)[K] = 0.0;
			(*qyz)[K] = 0.0;
			(*qzx)[K] = 0.0;
			(*px)[K] = mr_t*(*numMethodStates)(0,K);
			(*py)[K] = 0.0;
			(*pz)[K] = 0.0;
			(*o)[K] = 0.0;
		}
	}
	if(GetNumeraire()->GetType() == ARM_Numeraire::Cash)
	{
		for (K = 0; K < N; ++K)
		{
			(*qxx)[K] = minus_square_vol_t;
			(*qyy)[K] = 0.0;
			(*qzz)[K] = 0.0;
			(*qxy)[K] = 0.0;
			(*qyz)[K] = 0.0;
			(*qzx)[K] = 0.0;
			(*px)[K] = mr_t*(*numMethodStates)(0,K);
			(*py)[K] = 0.0;
			(*pz)[K] = 0.0;
		}
			*o = *(RiskNeutralDrift(timeIdx,numMethodStates,true));

	}
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: IntegratedBondSquaredVol
///	Returns : double
///	Action  : Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_HullWhite1F::IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const
{
	const ARM_ModelParamsHW1F * modelParams = static_cast<const ARM_ModelParamsHW1F*> (GetModelParams());
	return ARM_ModelParamsHW1F::HW1FZcCovariance( modelParams, modelParams, startTime, endTime, bondMaturity );
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: IntegratedBondCovariance
///	Returns : double
///	Action  : Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_HullWhite1F::IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const
{
	double result = 0;
	const ARM_ModelParamsHW1F * modelParams = static_cast<const ARM_ModelParamsHW1F*> (GetModelParams());
	result += ARM_ModelParamsHW1F::HW1FZcCovariance( modelParams, modelParams, startTime, endTime, bondMaturity1 );
	result += ARM_ModelParamsHW1F::HW1FZcCovariance( modelParams, modelParams, startTime, endTime, bondMaturity2 );
	result -= modelParams->FwdZcLocalVariance( startTime, endTime, bondMaturity1, bondMaturity2 );
	
	return 0.5*result;
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite1F
///	Routines: VolatilityScalarProduct
///	Returns : double (qui l'eut cru)
///	Action  :  Int_startTime^endTime Gamma(s,bondMaturity) * dW_s 
////////////////////////////////////////////////////
double ARM_HullWhite1F::VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const
{
	const ARM_CurveModelParam * cmp = dynamic_cast<const ARM_CurveModelParam *> (&otherModelVolatility);

	if ( !cmp )
		ARM_THROW( ERR_INVALID_ARGUMENT, "VolatilityScalarProduct: wrong ModelParam Type" );

	const ARM_ModelParamsHW1F* modp = static_cast<const ARM_ModelParamsHW1F*> (GetModelParams());
	return modp->HW1FEqFxZcCovariance( *cmp, modp, startTime, endTime, bondMaturity );
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

