/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file HW2F.cpp
 *
 *  \brief
 *
 *	\author  JM Prie
 *	\version 1.0
 *	\date October 2003
 */


/// this header comes first as it includes some preprocessor constants!
#include "gpmodels/HW2F.h"
#include "gpmodels/ModelParamsHW2F.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/irrate.h"


/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"

/// closed forms
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/spreadoption_lognormal_interface.h"
#include "gpclosedforms/spreadoption_shiftedlognormal_interface.h"
#include "gpclosedforms/vanilla_normal.h"

/// kernel
#include <inst/portfolio.h>
#include <util/foncstd.h>

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillaswaption.h"


CC_BEGIN_NAMESPACE( ARM )

#define FIRST_STATE_VARIABLE    0
#define SECOND_STATE_VARIABLE   1

const int nbQuadPoints=64;      /// Nb quadrature pts
const double quadMinx=-5.0;     /// -5.stdDev
const double quadMaxx=5.0;      ///  5.stdDev


//GaussLegendre_Coefficients Quadrature(nbQuadPoints);


////////////////////////////////////////////////////
///	Class  : HW2FVarianceFunction
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
HW2FVarianceFunction::HW2FVarianceFunction(const ARM_ModelParams* modelParams)
: itsModelParams(modelParams)
{}

HW2FVarianceFunction::~HW2FVarianceFunction()
{}

double HW2FVarianceFunction::operator () ( double x ) const
{
    ARM_GP_TriangularMatrix* varCovar = ((const ARM_ModelParamsHW2F* )itsModelParams)->StateLocalVariance(0.0,x);
    double value = varCovar->trace();

    delete varCovar;

    return value;
}




////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HullWhite2F::ARM_HullWhite2F( const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, const ARM_BoolVector& soFormulaFlags )
:	
	ARM_HullWhite(zc,params,soFormulaFlags),
	itsCalibParam(0),
	itsPF(0),
	itsSize(0)
{
	CC_ARM_SETNAME(ARM_HW2FMODEL);
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_HullWhite2F::ARM_HullWhite2F(const ARM_HullWhite2F& rhs)
: ARM_HullWhite(rhs),
	itsCalibParam(0),
	itsPF(0),
	itsSize(0)
{
    // Copy class attributes if any
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_HullWhite2F::~ARM_HullWhite2F()
{
	FreeFastCalib();
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_HullWhite2F& ARM_HullWhite2F::operator=(const ARM_HullWhite2F& rhs)
{
	if(this != &rhs)
	{
		ARM_HullWhite::operator=(rhs);
        // Copy class attributes if any
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_HullWhite2F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "2F Hull & White Model\n";
    os << indent << "---------------------\n";
	if(IsApproxSOFormula())
		os << indent << "\nApprox CMS Spread Option Formula\n\n";
	else
	{
		os << indent << "\nExact CMS Spread Option Formula, ";
		if(IsDeepITMSOFormula())
			os << indent << "Speed 2D integration\n\n";
		else
			os << indent << "Accurate 2D integration\n\n";
	}
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_HullWhite2F::Clone() const
{
	return new ARM_HullWhite2F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_HullWhite2F::ValidateModelParams(const ARM_ModelParams& params) const
{
    if( !params.DoesModelParamExist(ARM_ModelParamType::Volatility)				|| 
        !params.DoesModelParamExist(ARM_ModelParamType::MeanReversion)			||
        !params.DoesModelParamExist(ARM_ModelParamType::VolatilityRatio)		||
        !params.DoesModelParamExist(ARM_ModelParamType::MeanReversionSpread)	||
        !params.DoesModelParamExist(ARM_ModelParamType::Correlation))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "A vol, a MRS, a vol ratio, a MRS spread and a correl model parameters are required for 2F Hull & White model");
    }
    else if(
        params.GetModelParam(ARM_ModelParamType::MeanReversion).size() != 1			||
        params.GetModelParam(ARM_ModelParamType::MeanReversionSpread).size() != 1)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "This version of 2F Hull & White model needs constant MRS and MRS spread");
    }
    else
        return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: LocalDrifts
///	Returns: For each time step, a vector with the
///          local drift of state variables
///	Action : Compute local relative drifts of state
///          variables between time steps
////////////////////////////////////////////////////
void ARM_HullWhite2F::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts) const
{
	size_t nbSteps	= timeSteps.size();
	relativeDrifts	= ARM_GP_MatrixPtr( new ARM_GP_Matrix( nbSteps-1, 2, 0.0 ) );
	relativeDrifts  = ( (const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalDriftVec(timeSteps);
	absoluteDrifts	= ARM_GP_MatrixPtr(NULL);
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 2D triangular matrix
///	Action : Compute local variances & 
///          covariances (from 0) of state variables
///          between each time step
////////////////////////////////////////////////////
void ARM_HullWhite2F::ModelStateLocalVariances(
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
		ARM_THROW( ERR_INVALID_ARGUMENT, "localVariances.size() != offsetIndex" );
#endif
	localVariances.resize((nbSteps-1)*(modelRank+1));
	size_t i;

	// There is no model state

	ARM_GP_TriangularMatrix identity(2,1.0);
	identity(0,1)=identity(1,0)=0.0;

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = static_cast<ARM_GP_TriangularMatrix*>(identity.Clone());
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of 2D triangular matrix
///	Action : Compute local variances & 
///          covariances (from 0) of state variables
///          between each time step
////////////////////////////////////////////////////
void ARM_HullWhite2F::NumMethodStateLocalVariances(
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
		localVariances[offsetIndex+i] = ((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(step,nextStep)
			;

		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: NumMethodStateGlobalVariances
///	Returns: A vector of 2D triangular matrix
///	Action : global variances & 
///          covariances (from 0) of state variables
///          between each time step
////////////////////////////////////////////////////
void ARM_HullWhite2F::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances ) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],nextStep;
	size_t offsetIndex	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	variances.resize(nbSteps*(modelNb+1));
	size_t i;

    variances[offsetIndex]=new ARM_GP_TriangularMatrix(2,0.0);
	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];

        /// [i+1] => variances/covariances from 0 -> ti+1
        variances[offsetIndex+i+1] = ((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(0.0,nextStep);

		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_HullWhite2F::VarianceToTime(double var,double minTime,double maxTime) const
{
    HW2FVarianceFunction f(GetModelParams());
    T_BrentSolver<HW2FVarianceFunction> solver(f,var,0.0001);
	solver.setInitialGuess(0.0,minTime,maxTime);
    return solver.Solve();
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF. With
///          the Cash numeraire could be optimised
///          using phi(t) (see code). With Zc
///          numeraires, it uses a HJM like version.
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
    // Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
    double zct,numMatTime;
    double zcVar;

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new std::vector<double>(payoffSize,zcT) );
    }

    // Volatility computation (ARM_ModelParamsHW1F class is pure virtual)
    std::vector<double>& stateFact;
    int i,nbStates=states->size();
    if(    GetNumeraire() == ARM_NumerairePtr(NULL)
		|| GetNumeraire()->GetType() == ARM_Numeraire::Cash
		|| GetNumeraire()->GetType() == ARM_Numeraire::RollingCash
	  )
    {
        if(evalTime < maturityTime)
        {
            // Could be optimized using phis(t)...?
            stateFact = ((const ARM_ModelParamsHW2F*) GetModelParams())->BetatT(evalTime,maturityTime);
            zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
            zcVar=((const ARM_ModelParamsHW2F* const) GetModelParams())->ZcVarianceSpread(evalTime,evalTime,maturityTime,evalTime);
        }
        else
            return ARM_VectorPtr(new std::vector<double>(nbStates,1.0));
    }
    else
    {
        /// TerminalZc cases
        if(evalTime < maturityTime)
            stateFact = ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(evalTime,maturityTime);
        else
            stateFact = new std::vector<double>(2,0.0);
        zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
        numMatTime=GetNumeraire()->GetMaturity();
        zcVar = ((const ARM_ModelParamsHW2F* const) GetModelParams())->FwdZcLocalVariance(0.0,evalTime,maturityTime,numMatTime)
			-   ((const ARM_ModelParamsHW2F* const) GetModelParams())->FwdZcLocalVariance(0.0,evalTime,evalTime,numMatTime);
    }

	size_t modelNb = GetModelNb();
    ARM_VectorPtr values(new std::vector<double>(nbStates));
    for(i=0;i<nbStates;++i)
        (*values)[i]=zcT/zct*exp( -0.5*zcVar
                        - (*stateFact)[FIRST_STATE_VARIABLE] * states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE)
                        - (*stateFact)[SECOND_STATE_VARIABLE] * states->GetModelState(i,modelNb+SECOND_STATE_VARIABLE) );

    delete stateFact;
    return values;
}



////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: LogNorBondVanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : Log normal bond assuption to
///          get a closed form formula for
///          vanillia swaption price
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::LogNorBondVanillaSwaption(
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
    /// Save times and factors of fwd Zc
    /// Compute Zc and bond datas
    ARM_VectorPtr zcFloatStart;
    vector< ARM_VectorPtr > zcFlowCouponPay;
    vector< double > couponPayPeriod;
    vector< double > flowPayTimes;

    ComputeFlowZc(curveName,evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
        states,zcFloatStart,zcFlowCouponPay,flowPayTimes,couponPayPeriod);

    int i,nbFlows=zcFlowCouponPay.size();
    int nbStates=zcFloatStart->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates));


    /// The forward bond (sum(i=0,nbFlows) couponPayPeriod[i].zcFlowCouponPay[i]/zcFloatStart)
    /// is assumed to have a lognormal distribution. We fit its 1st
    /// and 2nd moments to get its distribution parameters
    std::vector<double> unusedVars;
    int j;
    double ti,ci,invZc0,zci,fwdZcCovar,fwdBond,fwdBondTerm;
    double stdDev;
    double price;
    for(int stIdx=0;stIdx<nbStates;++stIdx)
    {
        stdDev=0.0;
        fwdBond=0.0;
        invZc0 = 1.0/(*zcFloatStart)[stIdx];
        for(i=0;i<nbFlows;++i)
        {
            ti=flowPayTimes[i];
            ci=couponPayPeriod[i];
            zci = (*(zcFlowCouponPay[i]))[stIdx];
            fwdBondTerm = ci*zci*invZc0;
            fwdBond += fwdBondTerm;
            for(j=0;j<i;++j)
            {
                fwdZcCovar=((const ARM_ModelParamsHW* const) GetModelParams())->FwdZcLocalCovariance(evalTime,swapResetTime,floatStartTime,ti,floatStartTime,flowPayTimes[j],unusedVars);
                stdDev += 2*fwdBondTerm * exp(fwdZcCovar) * couponPayPeriod[j] * (*(zcFlowCouponPay[j]))[stIdx] * invZc0;
            }
            fwdZcCovar=((const ARM_ModelParamsHW* const) GetModelParams())->FwdZcLocalCovariance(evalTime,swapResetTime,floatStartTime,ti,floatStartTime,ti,unusedVars);
            stdDev += exp(fwdZcCovar) * fwdBondTerm * fwdBondTerm;
        }
        stdDev = log(stdDev/(fwdBond*fwdBond));

        /// Compute Black price
        if(stdDev > stdDevLimit)
        {
            stdDev = sqrt(stdDev);
            double vol=stdDev/sqrt((swapResetTime-evalTime)/K_YEAR_LEN); // for debug
            price = BlackSholes_Formula(fwdBond,stdDev,(*zcFloatStart)[stIdx],1.0,callPutBond);
        }
        else if(payRec==K_RCV)
            price = (fwdBond > 1 ? (*zcFloatStart)[stIdx] * (fwdBond-1) : 0.0);
        else
            price = (fwdBond < 1 ? (*zcFloatStart)[stIdx] * (1-fwdBond) : 0.0);

        (*values)[stIdx] = swapNotional * price;
    }

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: IntegrationVanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : Log normal bond assuption to
///          get a closed form formula for
///          vanillia swaption price
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::IntegrationVanillaSwaption(
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
    const double stdDevLimit = 1.0e-10;
    const double invSqr2Pi = 0.398942280402;
    const double invSqrPi = 0.564189583546;

    int i,j,nbFixFlows=fixPayTimes.size();

    /// Check if we can degenerate in the 1F case
    double varLimit=stdDevLimit*stdDevLimit;

    ARM_GP_TriangularMatrix* StateLocVar=((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(evalTime,swapResetTime);

    if( (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE) < varLimit ||
        (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE) < varLimit )
    {
        /// The 1F case will give a very good price
        return ARM_VectorPtr( ARM_HullWhite::VanillaSwaption(
            curveName,
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
			isConstantStrike) );
    }
	/// now it is computed only once 
	GaussLegendre_Coefficients Quadrature(nbQuadPoints);
    /// Compute Zc and bond datas
    ARM_VectorPtr zcFloatStart;
    vector< ARM_VectorPtr > zcFlowCouponPay;
    vector< double > couponPayPeriod;
    vector< double > flowPayTimes;

	int idx=-1;
	double effFloatStartTime = floatStartTime;
	bool isFwdStartNotional = isFwdStart(fixNotional,idx);
	if (isFwdStartNotional)
		effFloatStartTime = fixPayTimes[idx];

    ComputeFlowZc(curveName,evalTime,effFloatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
        states,zcFloatStart,zcFlowCouponPay,flowPayTimes,couponPayPeriod);

    int nbFlows=zcFlowCouponPay.size();
    int nbStates=zcFloatStart->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));


    /// Orthogonalisation
    (*StateLocVar)(FIRST_STATE_VARIABLE, FIRST_STATE_VARIABLE)  = sqrt((*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE));
    (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE) /= (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE);
    (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE)-= (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE) *
                                                                  (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE);
    (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE) = sqrt((*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE));


    /// Compute fwd Zc coef in a uncorrelated world (i.e. ACP)
    std::vector<double> *betatT0=((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime,effFloatStartTime);
    ARM_VectorVector betatTiT0(nbFlows);
    vector< double > logDrift(nbFlows);
    vector< double > stdDev1(nbFlows);
    vector< double > stdDev2(nbFlows);
    for(i=0;i<nbFlows;++i)
    {
        betatTiT0[i] = ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime,flowPayTimes[i]);
        (*(betatTiT0[i]))[FIRST_STATE_VARIABLE] -= (*betatT0)[FIRST_STATE_VARIABLE];
        (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] -= (*betatT0)[SECOND_STATE_VARIABLE];

        stdDev1[i] = (*(betatTiT0[i]))[FIRST_STATE_VARIABLE] * (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE)
                     + (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] * (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE);

        stdDev2[i] = (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] * (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE);

        logDrift[i] = -0.5*(stdDev1[i]*stdDev1[i] + stdDev2[i]*stdDev2[i]);
    }

    /// Gauss Legendre numerical integration loop
    double scalex = 0.5*(quadMaxx-quadMinx);
    double scalew = scalex*invSqr2Pi;
    double x,w;


    vector< double > conditionalDrift(nbFlows);
    ARM_VectorPtr prices;
    for(i=0;i<nbQuadPoints;++i)
    {
        /// Rescaling in [quadMinx,quadMaxx]
        x = quadMinx + (Quadrature.get_point(i)+1.0)*scalex;
        w = Quadrature.get_weight(i) * scalew;

        /// Update drift term
        for(j=0;j<nbFlows;++j)
            conditionalDrift[j] = exp(logDrift[j] + stdDev1[j] * x);

        /// Compute conditionnal price
        prices = VanillaSwaptionPrice(fixNotional,floatNotional,zcFloatStart,zcFlowCouponPay,
			couponPayPeriod,stdDev2,conditionalDrift,strikesPerState,callPut,states,isConstantNotional,
			isConstantSpread,
			isConstantStrike,idx);

        /// Update GL sum
        for(j=0;j<nbStates;++j)
            (*values)[j] += (*prices)[j] * w * exp(-0.5*x*x);
    }

    /// Free memory
    delete betatT0;
    for(i=0;i<nbFlows;++i)
        delete betatTiT0[i];

    delete StateLocVar;

    return values;
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite
///	Routines: VanillaSpreadOption
///	Returns : ARM_VectorPtr
///	Action  : Price spread option by a 2D GL numerical
///			  integration
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_HullWhite2F::VanillaSpreadOptionLet(const string& curveName,
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
    ARM_GP_TriangularMatrix* StateLocVar=((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(evalTime,resetTime);

	double varLimit = 1.0e-4*sqrt((resetTime-evalTime)/K_YEAR_LEN); /// volLimit=1bp
	varLimit *= varLimit;

    if( evalTime < K_NEW_DOUBLE_TOL &&
		( ((*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE) < varLimit &&
           (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE) < varLimit) ||
		  IsApproxSOFormula() ) )
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

	/// Compute an exact price by 2D numerical integation
	double longFloatStartTime = swapLongFloatStartTime;
	double longFloatEndTime = swapLongFloatEndTime;
	double shortFloatStartTime = swapShortFloatStartTime;
	double shortFloatEndTime = swapShortFloatEndTime;
	std::vector<double> *longFixPayTimes=NULL;
	std::vector<double> *longFixPayPeriods=NULL;
	std::vector<double> *shortFixPayTimes=NULL;
	std::vector<double> *shortFixPayPeriods=NULL;

	ARM_SwapRatePtr longSwapRate(NULL),shortSwapRate(NULL);
	if(&swapLongFixPayTimes && &swapShortFixPayTimes)
	{
		longFixPayTimes = const_cast<ARM_GP_Vector*>(&swapLongFixPayTimes);
		longFixPayPeriods = const_cast<ARM_GP_Vector*>(&swapLongFixPayPeriods);
		shortFixPayTimes = const_cast<ARM_GP_Vector*>(&swapShortFixPayTimes);
		shortFixPayPeriods = const_cast<ARM_GP_Vector*>(&swapShortFixPayPeriods);
	}
	else
	{
		/// compute swap schedule since they have not been computed yet !
		ARM_Currency* ccy = GetCurrency(GetModelName());
		double asOf = GetZeroCurve()->GetAsOfDate().GetJulian();
		char fixCalendar[100];
		ccy->CalcFixPayCal(fixCalendar);
		int  fixFreq	 = ccy->GetFixedPayFreq();
		int  fixDayCount = ccy->GetFixedDayCount();
		
		longSwapRate = ARM_SwapRate::CreateSwapRate(asOf,
													asOf+swapLongFloatStartTime, 
													asOf+swapLongFloatEndTime, 
													fixDayCount, 
													fixFreq, 
													fixCalendar);

		longFloatStartTime	= longSwapRate->floatStartTime;
		longFloatEndTime	= longSwapRate->floatEndTime;
		longFixPayTimes		= & longSwapRate->fixPayTimes;
		longFixPayPeriods	= & longSwapRate->fixPayPeriods;
		
		shortSwapRate = ARM_SwapRate::CreateSwapRate(asOf,
													asOf+swapShortFloatStartTime, 
													asOf+swapShortFloatEndTime, 
													fixDayCount, 
													fixFreq, 
													fixCalendar);

		shortFloatStartTime	= shortSwapRate->floatStartTime;
		shortFloatEndTime	= shortSwapRate->floatEndTime;
		shortFixPayTimes	= & shortSwapRate->fixPayTimes;
		shortFixPayPeriods	= & shortSwapRate->fixPayPeriods;
	}

	/// Compute an exact price using a 2D numerical integration of
	/// the payoff = MAX(callPut*(a1.S1-a2.K),0)


	/// Orthogonalisation
	double v11 = (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE);
	double v12 = (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE);
	double v22 = (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE);

	double stdDev11=0.0,stdDev12=0.0,stdDev21=0.0,stdDev22=0.0;
	if(v11>ARM_NumericConstants::ARM_TOLERENCE || v22>ARM_NumericConstants::ARM_TOLERENCE)
	{
		bool isPCA=false;
		if(isPCA)
		{
			/// 2D PCA
			double vTemp,e11,e21,e12,e22,stdDev1,stdDev2;
			if(v12 < -ARM_NumericConstants::ARM_TOLERENCE || v12 > ARM_NumericConstants::ARM_TOLERENCE)
			{
				vTemp = v11+v22;
				double delta = vTemp*vTemp - 4*(v11*v22-v12*v12);
				double v1 = 0.5*(vTemp+sqrt(delta));
				double v2 = 0.5*(vTemp-sqrt(delta));
				if(v1<v2)
				{
					/// v1 is always the greater eigen value
					vTemp=v2;
					v2=v1;
					v1=vTemp;
				}
				e21 = (v1-v11)/v12;
				double invNorm=1/sqrt(1+e21*e21);
				e11 = invNorm;
				e21 *= invNorm;
				e22 = (v2-v11)/v12;
				invNorm=1/sqrt(1+e22*e22);
				e12 = invNorm;
				e22 *= invNorm;
				stdDev1 = sqrt(v1);
				stdDev2 = sqrt(v2);
				stdDev11 = e11*stdDev1;
				stdDev12 = e12*stdDev2;
				stdDev21 = e21*stdDev1;
				stdDev22 = e22*stdDev2;
			}
			else
			{
				/// Already diagonalized
				if(v11<v22)
				{
					/// v1 is always the greater eigen value
					vTemp	= v22;
					stdDev2	= sqrt(v11);
					stdDev1	= sqrt(vTemp);
					stdDev11= 0.0;
					stdDev12= stdDev2;
					stdDev21= stdDev1;
					stdDev22= 0.0;
				}
				else
				{
					stdDev1	= sqrt(v11);
					stdDev2 = sqrt(v22);
					stdDev11= stdDev1;
					stdDev12= 0.0;
					stdDev21= 0.0;
					stdDev22= stdDev2;
				}
			}
		}
		else
		{
			/// Cholesky
			if(v11>ARM_NumericConstants::ARM_TOLERENCE)
			{
				stdDev11 = sqrt(v11);
				stdDev12 = 0.0;
				stdDev21 = v12 / stdDev11;
				stdDev22 = sqrt(v22 - stdDev21 * stdDev21);
			}
			else
			{
				stdDev22 = sqrt(v22);
				stdDev21 = 0.0;
				stdDev12 = 0.0; // to avoid pb
				stdDev11 = 0.0; // to avoid pb
			}
		}
	}



	/// Beta(Te,Tpay)
	ARM_VectorPtr betatTeTp(((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(resetTime,payTime));


	/// Long swap calculation
	/// ---------------------

	/// Get long standard swap datas
	ARM_VectorPtr longZcFloatStart;
	vector< ARM_VectorPtr > longZcFlowPay;
	vector< double > longFlowTimes;
	vector< double > longFlowPeriod;
	ComputeFlowZc(curveName,evalTime,longFloatStartTime,longFloatEndTime,*longFixPayTimes,
		*longFixPayPeriods,states,longZcFloatStart,longZcFlowPay,longFlowTimes,longFlowPeriod);
	size_t longNbFlows=longZcFlowPay.size();
	bool isLongSameEnd = (longNbFlows==swapLongFixPayTimes.size());

	/// Drift and diffusion factor of B(Te,TlongStart)
	ARM_VectorPtr longBetatTeT0(((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(resetTime,swapLongFloatStartTime));
	(*longBetatTeT0)[FIRST_STATE_VARIABLE]	-= (*betatTeTp)[FIRST_STATE_VARIABLE];
	(*longBetatTeT0)[SECOND_STATE_VARIABLE]	-= (*betatTeTp)[SECOND_STATE_VARIABLE];

	double longStdDev1T0 = (*longBetatTeT0)[FIRST_STATE_VARIABLE] * stdDev11 +
						   (*longBetatTeT0)[SECOND_STATE_VARIABLE] * stdDev21;

	double longStdDev2T0 = (*longBetatTeT0)[FIRST_STATE_VARIABLE] * stdDev12 +
						   (*longBetatTeT0)[SECOND_STATE_VARIABLE] * stdDev22;

	double longLogDriftT0 = -0.5*(longStdDev1T0*longStdDev1T0 + longStdDev2T0*longStdDev2T0);

	/// Drift and diffusion factor of B(Te,TlongFlowi)
	ARM_VectorPtr longBetatTeTi(NULL);
	vector< double > longLogDrift(longNbFlows);
	vector< double > longStdDev1(longNbFlows);
	vector< double > longStdDev2(longNbFlows);
	size_t i;
	for(i=0;i<longNbFlows;++i)
	{
		longBetatTeTi = ARM_VectorPtr( ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(resetTime,longFlowTimes[i]) );
		(*longBetatTeTi)[FIRST_STATE_VARIABLE]	-= (*betatTeTp)[FIRST_STATE_VARIABLE];
		(*longBetatTeTi)[SECOND_STATE_VARIABLE]	-= (*betatTeTp)[SECOND_STATE_VARIABLE];

		longStdDev1[i] = (*longBetatTeTi)[FIRST_STATE_VARIABLE] * stdDev11 +
						 (*longBetatTeTi)[SECOND_STATE_VARIABLE] * stdDev21;

		longStdDev2[i] = (*longBetatTeTi)[FIRST_STATE_VARIABLE] * stdDev12 +
						 (*longBetatTeTi)[SECOND_STATE_VARIABLE] * stdDev22;

		longLogDrift[i] = -0.5*(longStdDev1[i]*longStdDev1[i] + longStdDev2[i]*longStdDev2[i]);
	}

	if(longNbFlows == longFixPayTimes->size())
	{
		/// Add nominal payback flow to simplify further loops
		longZcFlowPay.push_back(longZcFlowPay[longNbFlows-1]);
		longFlowPeriod.push_back(1.0);
		longStdDev1.push_back(longStdDev1[longNbFlows-1]);
		longStdDev2.push_back(longStdDev2[longNbFlows-1]);
		longLogDrift.push_back(longLogDrift[longNbFlows-1]);
	}


	/// Short swap calculation
	/// ----------------------

	/// Get short standard swap datas
	ARM_VectorPtr shortZcFloatStart;
	vector< ARM_VectorPtr > shortZcFlowPay;
	vector< double > shortFlowTimes;
	vector< double > shortFlowPeriod;
	ComputeFlowZc(curveName,evalTime,shortFloatStartTime,shortFloatEndTime,*shortFixPayTimes,
		*shortFixPayPeriods,states,shortZcFloatStart,shortZcFlowPay,shortFlowTimes,shortFlowPeriod);
	size_t shortNbFlows=shortZcFlowPay.size();

	/// Drift and diffusion factor of B(Te,TshortStart)
	ARM_VectorPtr shortBetatTeT0(longBetatTeT0);
	double shortStdDev1T0=longStdDev1T0;
	double shortStdDev2T0=longStdDev2T0;
	double shortLogDriftT0=longLogDriftT0;
	if(swapShortFloatStartTime != swapLongFloatStartTime)
	{
		shortBetatTeT0 = ARM_VectorPtr( ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(resetTime,swapShortFloatStartTime) );
		(*shortBetatTeT0)[FIRST_STATE_VARIABLE]	-= (*betatTeTp)[FIRST_STATE_VARIABLE];
		(*shortBetatTeT0)[SECOND_STATE_VARIABLE]-= (*betatTeTp)[SECOND_STATE_VARIABLE];

		shortStdDev1T0 = (*shortBetatTeT0)[FIRST_STATE_VARIABLE] * stdDev11 +
						 (*shortBetatTeT0)[SECOND_STATE_VARIABLE] * stdDev21;

		shortStdDev2T0 = (*shortBetatTeT0)[FIRST_STATE_VARIABLE] * stdDev12 +
						 (*shortBetatTeT0)[SECOND_STATE_VARIABLE] * stdDev22;

		shortLogDriftT0 = -0.5*(shortStdDev1T0*shortStdDev1T0 + shortStdDev2T0*shortStdDev2T0);
	}

	/// Drift and diffusion factor of B(Te,TshortFlowi)
	ARM_VectorPtr shortBetatTeTi(NULL);
	vector< double > shortLogDrift(shortNbFlows);
	vector< double > shortStdDev1(shortNbFlows);
	vector< double > shortStdDev2(shortNbFlows);
	for(i=0;i<shortNbFlows;++i)
	{
		shortBetatTeTi = ARM_VectorPtr( ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(resetTime,shortFlowTimes[i]) );
		(*shortBetatTeTi)[FIRST_STATE_VARIABLE] -= (*betatTeTp)[FIRST_STATE_VARIABLE];
		(*shortBetatTeTi)[SECOND_STATE_VARIABLE]-= (*betatTeTp)[SECOND_STATE_VARIABLE];

		shortStdDev1[i] = (*shortBetatTeTi)[FIRST_STATE_VARIABLE] * stdDev11 +
						  (*shortBetatTeTi)[SECOND_STATE_VARIABLE] * stdDev21;

		shortStdDev2[i] = (*shortBetatTeTi)[FIRST_STATE_VARIABLE] * stdDev12 +
						  (*shortBetatTeTi)[SECOND_STATE_VARIABLE] * stdDev22;

		shortLogDrift[i] = -0.5*(shortStdDev1[i]*shortStdDev1[i] + shortStdDev2[i]*shortStdDev2[i]);
	}

	if(shortNbFlows == shortFixPayTimes->size())
	{
		/// Add nominal payback flow to simplify further loops
		shortZcFlowPay.push_back(shortZcFlowPay[shortNbFlows-1]);
		shortFlowPeriod.push_back(1.0);
		shortStdDev1.push_back(shortStdDev1[shortNbFlows-1]);
		shortStdDev2.push_back(shortStdDev2[shortNbFlows-1]);
		shortLogDrift.push_back(shortLogDrift[shortNbFlows-1]);
	}


	int nbStates=longZcFloatStart->size();
	ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));

    /// Free memory
    delete StateLocVar;

	double ratio = payPeriod*notional;
	ARM_VectorPtr zcPay = GetDiscountFunctor()->DiscountFactor(curveName,evalTime,payTime,states);

	size_t i1,i2,k;
	double vi,sLong,sShort;

	if(v11<=ARM_NumericConstants::ARM_TOLERENCE && v22<=ARM_NumericConstants::ARM_TOLERENCE)
	{
		/// Case of intrinsic value
		for(k=0;k<nbStates;++k)
		{
			sLong=0;
			for(i=0;i<longNbFlows;++i)
				sLong += longFlowPeriod[i] * (*(longZcFlowPay[i]))[k];
			sLong = ( (*longZcFloatStart)[k] - (*(longZcFlowPay[longNbFlows]))[k])/sLong;

			sShort=0;
			for(i=0;i<shortNbFlows;++i)
				sShort += shortFlowPeriod[i] * (*(shortZcFlowPay[i]))[k];
			sShort = ( (*shortZcFloatStart)[k] - (*(shortZcFlowPay[shortNbFlows]))[k])/sShort;

			vi = callPut*(coeffLong*sLong - coeffShort*sShort - strikes[k]);
			(*values)[k] = (vi<0 ? 0.0 : vi) * (*zcPay)[k] * ratio;
		}
		return values;
	}

	/// 2D numerical integration for NC(0,1) between [quadMinx,quadMaxx]
	/// ----------------------------------------------------------------

	double z1,z2,w1,w2;
	int nb1RefPts=6,nb2RefPts=8,nb1Strats=5,nb2Strats=5;
	if(IsDeepITMSOFormula())
	{
		/// Rough 2D integration because the payoff is no more conditionnal
		/// Only convexity adjustment is really computed
		nb1RefPts=7;
		nb2RefPts=7;
		nb1Strats=3;
		nb2Strats=3;
	}

	GaussLegendre_Coefficients GL1(nb1RefPts); /// becareful : never add both last arguments quadMinx && quadMaxx here !!
	GaussStratifiedHermiteLegendre_Coefficients GsHL1(&GL1,nb1Strats,quadMinx,quadMaxx);
	int n1 = GsHL1.get_order();

	GaussLegendre_Coefficients GL2(nb2RefPts);
	GaussStratifiedHermiteLegendre_Coefficients GsHL2(&GL2,nb2Strats,quadMinx,quadMaxx);
	int n2 = GsHL2.get_order();

	std::vector<double> longFlowStoc(longNbFlows+1),shortFlowStoc(shortNbFlows+1);
	double longStartStoc,shortStartStoc;
	std::vector<double> condValues(nbStates);
	for(i1=0;i1<n1;++i1)
	{
		z1 = GsHL1.get_point(i1);
		w1 = GsHL1.get_weight(i1);
		for(i2=0;i2<n2;++i2)
		{
			z2 = GsHL2.get_point(i2);
			w2 = GsHL2.get_weight(i2);

			longStartStoc = exp(longLogDriftT0 + longStdDev1T0*z1 + longStdDev2T0*z2);
			for(i=0;i<longNbFlows+1;++i)
				longFlowStoc[i] = longFlowPeriod[i] * exp(longLogDrift[i] + longStdDev1[i]*z1 + longStdDev2[i]*z2);

			shortStartStoc = exp(shortLogDriftT0 + shortStdDev1T0*z1 + shortStdDev2T0*z2);
			for(i=0;i<shortNbFlows+1;++i)
				shortFlowStoc[i] = shortFlowPeriod[i] * exp(shortLogDrift[i] + shortStdDev1[i]*z1 + shortStdDev2[i]*z2);

			for(k=0;k<nbStates;++k)
			{
				/// Only Zc ratios appears in swap rate so it
				/// not needed to divide all B(Te,Ti) by B(Te,Tp)
				sLong=0;
				for(i=0;i<longNbFlows;++i)
					sLong += (*(longZcFlowPay[i]))[k] * longFlowStoc[i];
				sLong = ( (*longZcFloatStart)[k] * longStartStoc - (*(longZcFlowPay[longNbFlows]))[k] * longFlowStoc[longNbFlows] )/sLong;

				sShort=0;
				for(i=0;i<shortNbFlows;++i)
					sShort += (*(shortZcFlowPay[i]))[k] * shortFlowStoc[i];
				sShort = ( (*shortZcFloatStart)[k] * shortStartStoc - (*(shortZcFlowPay[shortNbFlows]))[k] * shortFlowStoc[shortNbFlows] )/sShort;

				vi = callPut*(coeffLong*sLong - coeffShort*sShort - strikes[k]);
				condValues[k] += (vi<0 ? 0.0 : vi) * w2;
			}
		}
		for(k=0;k<nbStates;++k)
		{
			(*values)[k] += condValues[k] * w1;
			condValues[k]=0.0;
		}
	}

	for(k=0;k<nbStates;++k)
		(*values)[k] *= ((*zcPay)[k] * ratio);

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: IntegrationVanillaSwaptionNew
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : integration of HW1F-like swaptions price
///          (achaix's method)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::IntegrationVanillaSwaptionNew(
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
    
	/// NOT IMPLEMENTED YET
	const double stdDevLimit = 1.0e-10;
    const double invSqr2Pi = 0.398942280402;
    const double invSqrPi = 0.564189583546;

    int i,j,nbFixFlows=fixPayTimes.size();

    /// Check if we can degenerate in the 1F case
    double varLimit=stdDevLimit*stdDevLimit;
	/// now it is computed oly once
	GaussLegendre_Coefficients Quadrature(nbQuadPoints);

    ARM_GP_TriangularMatrix* StateLocVar=((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(evalTime,swapResetTime);

    if( (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE) < varLimit ||
        (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE) < varLimit )
    {
        /// The 1F case will give a very good price
        return ARM_VectorPtr( ARM_HullWhite::VanillaSwaption(
            curveName,
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
			isConstantStrike) );
    }


    /// Compute Zc and bond datas
    ARM_VectorPtr zcFloatStart;
    vector< ARM_VectorPtr > zcFlowCouponPay;
    vector< double > couponPayPeriod;
    vector< double > flowPayTimes;

    ComputeFlowZc(curveName,evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,
        states,zcFloatStart,zcFlowCouponPay,flowPayTimes,couponPayPeriod);

    int nbFlows=zcFlowCouponPay.size();
    int nbStates=zcFloatStart->size();
    ARM_VectorPtr values(new std::vector<double>(nbStates,0.0));


    /// Orthogonalisation
    (*StateLocVar)(FIRST_STATE_VARIABLE, FIRST_STATE_VARIABLE)  = sqrt((*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE));
    (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE) /= (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE);
    (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE)-= (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE) *
                                                                  (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE);
    (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE) = sqrt((*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE));


    /// Compute fwd Zc coef in a uncorrelated world (i.e. ACP)
    std::vector<double> *betatT0=((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime,floatStartTime);
    ARM_VectorVector betatTiT0(nbFlows);
    vector< double > logDrift(nbFlows);
    vector< double > stdDev1(nbFlows);
    vector< double > stdDev2(nbFlows);
    for(i=0;i<nbFlows;++i)
    {
        betatTiT0[i] = ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime,flowPayTimes[i]);
        (*(betatTiT0[i]))[FIRST_STATE_VARIABLE] -= (*betatT0)[FIRST_STATE_VARIABLE];
        (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] -= (*betatT0)[SECOND_STATE_VARIABLE];

        stdDev1[i] = (*(betatTiT0[i]))[FIRST_STATE_VARIABLE] * (*StateLocVar)(FIRST_STATE_VARIABLE,FIRST_STATE_VARIABLE)
                     + (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] * (*StateLocVar)(SECOND_STATE_VARIABLE,FIRST_STATE_VARIABLE);

        stdDev2[i] = (*(betatTiT0[i]))[SECOND_STATE_VARIABLE] * (*StateLocVar)(SECOND_STATE_VARIABLE,SECOND_STATE_VARIABLE);

        logDrift[i] = -0.5*(stdDev1[i]*stdDev1[i] + stdDev2[i]*stdDev2[i]);
    }

    /// Gauss Legendre numerical integration loop
    double scalex = 0.5*(quadMaxx-quadMinx);
    double scalew = scalex*invSqr2Pi;
    double x,w;


    vector< double > conditionalDrift(nbFlows);
    ARM_VectorPtr prices;
    for(i=0;i<nbQuadPoints;++i)
    {
        /// Rescaling in [quadMinx,quadMaxx]
        x = quadMinx + (Quadrature.get_point(i)+1.0)*scalex;
        w = Quadrature.get_weight(i) * scalew;

        /// Update drift term
        for(j=0;j<nbFlows;++j)
            conditionalDrift[j] = exp(logDrift[j] + stdDev1[j] * x);

        /// Compute conditionnal price
        prices = VanillaSwaptionPrice(fixNotional,floatNotional,zcFloatStart,zcFlowCouponPay,
			couponPayPeriod,stdDev2,conditionalDrift,strikesPerState,callPut,states,isConstantNotional,
			isConstantSpread,
			isConstantStrike);

        /// Update GL sum
        for(j=0;j<nbStates;++j)
            (*values)[j] += (*prices)[j] * w * exp(-0.5*x*x);
    }

    /// Free memory
    delete betatT0;
    for(i=0;i<nbFlows;++i)
        delete betatTiT0[i];

    delete StateLocVar;

    return values;
}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: VanillaSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : Approximated closed form formula for standard
///          swaption (i.e. on standard swap with
///          a "double notional" evaluation of its
///          floating leg)
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::VanillaSwaption(
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
	int idx = -1;
    bool isFwdStartNotional = isFwdStart(fixNotional,idx);
	if (!(isConstantNotional || isFwdStartNotional))
		return VariableNotionalSwaption (curveName, 
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

	
	if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
    {
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "Swaption pricing not implemented for different discount & fixing curves" );
    }

    // Flags for testing purpose
    bool is1FApprox=false; 
    bool isLNBondApprox=false;

    if(is1FApprox)
        /// 1 factor pricing approximation i.e. variance/covariance
        /// matrix of fwd Zc has only one eigen value > 0...
        /// ...then call the default implementation
        return ARM_VectorPtr( ARM_HullWhite::VanillaSwaption(
            curveName,
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
			isConstantStrike) );

    else if(isLNBondApprox)
        /// The fwd bond lognormal distribution is assumed LN
        return ARM_VectorPtr( LogNorBondVanillaSwaption(
            curveName,
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
			isConstantStrike) );

    /// 2D integration using Gauss-Legendre and conditionnal
    /// analytic prices
    return ARM_VectorPtr( IntegrationVanillaSwaption(
        curveName,
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
		isConstantStrike) );

}


////////////////////////////////////////////////////
///	Class  : ARM_HullWhite2F
///	Routine: VariableNotionalSwaption
///	Returns: a vector of Swaption(t,S(R,Ti),K)
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_HullWhite2F::VariableNotionalSwaption(
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
	size_t i, j, k,  size_cpn   = fixPayTimes.size();
	vector<double>  E(size_cpn), beta1(size_cpn), beta2(size_cpn);
	ARM_GP_Matrix Exp(size_cpn, size_cpn);
	double startdf = GetZeroCurve()->DiscountPrice(floatStartTime/K_YEAR_LEN);
	double flow;
	
	for (i=0; i<size_cpn-1; i++)
	{
		flow  = floatNotional[i] + strikesPerState(0,i) * fixPayPeriods[i] * fixNotional[i] - floatNotional[i+1] ;
		flow *= GetZeroCurve()->DiscountPrice (fixPayTimes[i]/K_YEAR_LEN) / startdf;
		E[i]  = flow;
	}
	
	flow = floatNotional[size_cpn-1] + fixPayPeriods[size_cpn-1] * strikesPerState(0,size_cpn-1) * fixNotional[size_cpn-1];
	flow *= GetZeroCurve()->DiscountPrice (fixPayTimes[size_cpn-1]/K_YEAR_LEN) / startdf;
	E[size_cpn-1]  = flow;

	/// bond option strike
	double strike = floatNotional[0];

	/// phi
	ARM_GP_TriangularMatrix* phi = ((const ARM_ModelParamsHW2F* const) GetModelParams())->StateLocalVariance(0.0, swapResetTime);
	double phi11 = (*phi)(0, 0);
	double phi22 = (*phi)(1, 1);
	double phi12 = (*phi)(1, 0);
	delete phi;
		
	/// beta
	std::vector<double>& beta; 
	beta = ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime, floatStartTime);
	double beta1start = (*beta)[0];
	double beta2start = (*beta)[1];
	delete beta;

	for (i=0; i<size_cpn; i++)
	{
		beta = ((const ARM_ModelParamsHW2F* const) GetModelParams())->BetatT(swapResetTime, fixPayTimes[i]);
		beta1[i] = (*beta)[0];
		beta2[i] = (*beta)[1];
		delete beta;
	}

	for (i=0; i<size_cpn; i++)
	{
		for (j=i; j<size_cpn; j++)
		{
				
			Exp(i,j) = exp(	(beta1[i] - beta1start) * (beta1[j] - beta1start) * phi11
						  +	(beta2[i] - beta2start) * (beta2[j] - beta2start) * phi22
						  +	(beta2[i] - beta2start) * (beta1[j] - beta1start) * phi12
						  +	(beta1[i] - beta1start) * (beta2[j] - beta2start) * phi12  );
		}
	}

	for (i=0; i<size_cpn; i++)
		for (j=0; j<i; j++)
			Exp(i,j) = Exp(j,i);

		
	/// decl
	double Tmp1, Tmp2;
	double Eplus, Eminus;
	double Vplus, Vminus;
	double Splus, Sminus;
	double CoefPlus, CoefMinus;
	double ShiftPlus, ShiftMinus;
	double SigPlus, SigMinus;
	double CovPlusMinus;
	double Correl;
	

	/// Moments of order 1, 2, 3
	double	M1plus	(0.0), 
			M2plus	(0.0), 
			M3plus	(0.0),
			M1minus (0.0),
			M2minus (0.0),
			M3minus (0.0);


	// -- M1 -
	for (i=0; i<size_cpn; i++)
	{
		if (E[i]>0) 
			M1plus  += E[i];
		else
			M1minus -= E[i];
	}
	
	// -- M2 --
	
	// non diag terms
	for (i=0; i<size_cpn; i++)
	{
		if (E[i]>0)
		{
			for (j=0; j<i; j++)
			{
				if (E[j]>0)
					M2plus += E[i] * E[j] * Exp(i,j);
			}
		}
		else if (E[i]<0)
		{
			for (j=0; j<i; j++)
			{
				if (E[j]<0)
					M2minus += E[i] * E[j] * Exp(i,j);
			}
		}
	}

	M2plus  *= 2.0;
	M2minus *= 2.0;
	
	// diag terms
	for (i=0; i<size_cpn; i++)
	{
		if (E[i]>0)
			M2plus  += E[i] * E[i] * Exp(i,i);
		else
			M2minus += E[i] * E[i] * Exp(i,i);
	}


	
	// -- M3 --	
	// diag terms
	for (i=0; i<size_cpn; i++)
	{
		Tmp1 = E[i];
		Tmp2 = Exp(i,i);

		if (E[i]>0)
			M3plus  +=  Tmp1 * Tmp1 * Tmp1 * Tmp2 * Tmp2 * Tmp2;
		else
			M3minus -=  Tmp1 * Tmp1 * Tmp1 * Tmp2 * Tmp2 * Tmp2;
	}

	// a^2 * b
	for (i=0; i<size_cpn; i++)
	{
		if (E[i]>0)
		{
			for (j=0; j<size_cpn; j++)
			{
				if (j!=i)
				{
					if (E[j]>0)
					{
						Tmp1 = E[i];
						Tmp2 = Exp(i,j);
						M3plus += 3.0 * Tmp1 * Tmp1 * E[j] * Exp(i,i) * Tmp2 * Tmp2;
					}
				}
			}
		}
		else if (E[i]<0)
		{
			for (j=0; j<size_cpn; j++)
			{
				if (j!=i)
				{					
					if (E[j]<0)
					{
						Tmp1 = E[i];
						Tmp2 = Exp(i,j);
						M3minus -= 3.0 * Tmp1 * Tmp1 * E[j] * Exp(i,i) * Tmp2 * Tmp2;
					}
				}
			}
		}
	}

	// a * b * c
	for (i=0; i<size_cpn; i++)
	{
		if (E[i]>0)
		{
			for (j=0; j<i; j++)
			{
				if (E[j]>0)
				{
					for (k=0; k<j; k++)
					{
						if (E[k]>0)
							M3plus += 6.0 * E[i] * E[j] * E[k] * Exp(i,j) * Exp(i,k) * Exp(j,k);
					}
				}
			}
		}
		else if (E[i]<0)
		{
			for (j=0; j<i; j++)
			{
				if (E[j]<0)
				{
					for (k=0; k<j; k++)
					{
						if (E[k]<0)
							M3minus -= 6.0 * E[i] * E[j] * E[k] * Exp(i,j) * Exp(i,k) * Exp(j,k);
					}
				}
			}
		}
	}


	Eplus	= M1plus;
	Eminus	= M1minus;
	Vplus	= M2plus  - Eplus  * Eplus;
	Vminus	= M2minus - Eminus * Eminus;
	Splus   = M3plus  - Eplus  * Eplus  * Eplus  - 3 * Eplus  * Vplus;
	Sminus  = M3minus - Eminus * Eminus * Eminus - 3 * Eminus * Vminus;
		
	CoefPlus  = 4.0 * Splus ;
	CoefPlus += 4.0 * sqrt(4.0 * Vplus * Vplus * Vplus + Splus * Splus);
	CoefPlus  = pow(CoefPlus, 1./3.);

	CoefMinus  = 4.0 * Sminus ;
	CoefMinus += 4.0 * sqrt(4.0 * Vminus * Vminus * Vminus + Sminus * Sminus);
	CoefMinus  = pow(CoefMinus, 1./3.);

	if (Eplus)
	{
		ShiftPlus  = Vplus * CoefPlus;
		ShiftPlus /= 0.5 * CoefPlus * CoefPlus - 2.0 * Vplus;
		ShiftPlus -= Eplus;
		SigPlus  = sqrt( log(1.0 + Vplus  / pow(Eplus  + ShiftPlus,  2) ) );
	}
	else
	{
		ShiftPlus = 0.0;
		SigPlus   = 0.0;
	}

	if (Eminus)
	{
		ShiftMinus  = Vminus * CoefMinus;
		ShiftMinus /= 0.5 * CoefMinus * CoefMinus - 2.0 * Vminus;
		ShiftMinus -= Eminus;
		SigMinus = sqrt( log(1.0 + Vminus / pow(Eminus + ShiftMinus, 2) ) );
	}
	else
	{
		ShiftMinus = 0.0;
		SigMinus   = 0.0;
	}

		
	CovPlusMinus = 0.0;

	if (Eminus && Eplus)
	{		
		for (i=0; i<size_cpn; i++)
			for (j=0; j<i; j++)
			{
				if (E[i] * E[j]<0) 
				{					
					CovPlusMinus -= E[i] * E[j] * Exp(i,j);
				}
			}
	
		CovPlusMinus -= Eplus * Eminus;
		Correl = log( 1.0 + CovPlusMinus / ((Eplus + ShiftPlus)*(Eminus + ShiftMinus)) ) / (SigPlus * SigMinus) ;
	}
	else
		Correl = 0.0;

	
	// il se peut que ca dépasse légèrement....
	if (Correl>1)
	{
		if (Correl>1.10)
			ARM_THROW( ERR_INVALID_ARGUMENT, "HW2F / Variable Notional Swaption: unresolved problem" )
		else 
			Correl = 1.0;
	}

	if (Correl<-1)
	{
		if (Correl<-1.10)
			ARM_THROW( ERR_INVALID_ARGUMENT, "HW2F / Variable Notional Swaption: unresolved problem" )
		else 
			Correl = -1.0;
	}
	
		
	double Te = swapResetTime / K_YEAR_LEN;
	double sqrtTe = sqrt(Te);
	double VolPlus  = SigPlus  / sqrtTe;
	double VolMinus = SigMinus / sqrtTe;

	double BsSpreadStrike = strike - ShiftMinus + ShiftPlus;
	double FwdPlus  = Eplus  + ShiftPlus;
	double FwdMinus = Eminus + ShiftMinus;
	double BsResult;

	if (Eminus && Eplus)
	{	
		BsResult = SpreadOption(FwdMinus, FwdPlus, VolMinus, VolPlus, 0.0, 0.0, Correl, 0.0, BsSpreadStrike, Te, - callPut);
	}
	else if (Eplus)
	{
		BsResult = BS (FwdPlus, BsSpreadStrike, Te, VolPlus, -callPut);
	}
	else if (Eminus)
	{
		BsResult = BS (FwdPlus, -BsSpreadStrike, Te, VolMinus, callPut);
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "HW2F / Variable Notional Swaption: unresolved problem" )
		
	
	double price = startdf * BsResult;

	return ARM_VectorPtr(new std::vector<double>(1,price));
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_HullWhite2F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_HullWhite2F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 2F , nbPayoffs = 0)
	return new ARM_PricingStates(bucketSize,2,0,2);
}



////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routines: AdviseBreakPointTimes
///	Returns : void
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_HullWhite2F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->size();  
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
        /// just return NULL
	case ARM_ModelParamType::MeanReversion:
    case ARM_ModelParamType::VolatilityRatio:
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
    case ARM_ModelParamType::MeanReversionSpread:
    case ARM_ModelParamType::Correlation:
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
    default:
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, 
            "Unknown type... an HW2F only supports volatility, meanreversion, volatilityRatio, mean reversion spread and correlation" );
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_HullWhite2F::EulerLocalDrifts( const std::vector<double>& timeSteps,
		ARM_GP_MatrixPtr& relativeDrifts,
		ARM_GP_MatrixPtr& absoluteDrifts) const
{
    /// FIX FIX : only valid for cst MRS !!!
	double MR1 = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
	double MR2 = -(GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversionSpread).GetValueAtPoint(0) + MR1);
	double timeStep;

	relativeDrifts = ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),2,-MR1 ) );
    for(size_t i=0;i<timeSteps.size()-1;++i)
	{
		timeStep = (timeSteps[i+1]-timeSteps[i])/K_YEAR_LEN;
        (*relativeDrifts)(i,0) *= timeStep;
		(*relativeDrifts)(i,1) = MR2*timeStep;
	}

	absoluteDrifts= ARM_GP_MatrixPtr( new ARM_GP_Matrix(timeSteps.size(),2, 0.0 ) );
}

////////////////////////////////////////////////////
///	Class   : ARM_HullWhite2F
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_HullWhite2F::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
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
		double VolRatio = (GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio)).GetValueAtPoint(0);

		/// Since volatilies are stepUpRight!!
		for( std::vector<double>::iterator iter = values.begin() ; iter != values.end()-1 ; ++iter)
			(*iter) = *(iter+1);

		std::vector<double> volsVec,d1VolsVec;
		VectorValuesAndDerivativesLinearMidPoints(times,values,timeSteps,volsVec,d1VolsVec);

		/// factor by line and zero correl because in one factor!
		vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,volsVec.size()) );
		d1Vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,d1VolsVec.size() ) );
		correls	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1, volsVec.size(), (GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation)).GetValueAtPoint(0) ));

		for(size_t i=0; i<d1VolsVec.size(); ++i )
		{
			d1VolsVec[i] *= K_YEAR_LEN;
			vols->Elt(0,i) = volsVec[i];
			vols->Elt(1,i) = VolRatio * volsVec[i];
			d1Vols->Elt(0,i) = d1VolsVec[i];
			d1Vols->Elt(1,i) = VolRatio*d1VolsVec[i];
		}
	}
	/// Instant. vols for PDE schemes
	else
	{
		/// Convention : 
		/// vols[i] = instant. vol of process between timeSteps[i] and timeSteps[i+1]
		///
		double VolRatio = (GetModelParams()->GetModelParam( ARM_ModelParamType::VolatilityRatio)).GetValueAtPoint(0);
		ARM_CurveModelParam& mp = (ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::Volatility);
		size_t timeStepsSize = timeSteps.size() ; 
		std::vector<double> volsVec(timeStepsSize);
		double epsToBeSafeWithInterpol = 1.0e-5;

		vols	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(2,volsVec.size()) );
		d1Vols	= ARM_GP_MatrixPtr( NULL );
		correls	= ARM_GP_MatrixPtr( new ARM_GP_Matrix(1, volsVec.size(), (GetModelParams()->GetModelParam( ARM_ModelParamType::Correlation)).GetValueAtPoint(0) ));

		for (size_t j=0; j<timeStepsSize; j++)
		{
			vols->Elt(0,j) = mp.GetValue(timeSteps[j]+epsToBeSafeWithInterpol);
			vols->Elt(1,j) = VolRatio * vols->Elt(0,j);
		}
	}
}

////////////////// FOR FAST CALIBRATION
////////////////// FOR FAST CALIBRATION
////////////////// FOR FAST CALIBRATION
void ARM_HullWhite2F::FreeFastCalib()
{
	itsCalibParam.resize(0);
	for (int i=0;i<itsPF.size();i++)
		delete itsPF[i];
	itsPF.resize(0);
	itsSize.resize(0);
}

std::vector<double>& csecurity::getRateWeight(double t_start,double t_end,const std::vector<double>& payTimes,const std::vector<double>& periods,double mrs1,double mrs2)
{
	ARM_ZeroCurvePtr ZcCurve = m_pmodel->GetZeroCurve();

	double tStart		= t_start/K_YEAR_LEN;
	double tEnd			= t_end/K_YEAR_LEN;
	double zcStart		= ZcCurve->DiscountPrice(tStart);
	double zcEnd		= ZcCurve->DiscountPrice(tEnd);
	double eStart1		= exp(-mrs1*tStart);
	double eEnd1		= exp(-mrs1*tEnd);
	double eStart2		= exp(-mrs2*tStart);
	double eEnd2		= exp(-mrs2*tEnd);
					
	double t,zc;
	double num1=0,num2=0,level=0,numbis=0.;
	for (size_t k=0;k<payTimes.size();k++)
	{
		t		= payTimes[k]/K_YEAR_LEN;
		zc		= ZcCurve->DiscountPrice(t);
		num1	+= exp(-mrs1*t)*periods[k]*zc;
		num2	+= exp(-mrs2*t)*periods[k]*zc;
		numbis	+= periods[k]*zc*t;
		level	+= periods[k]*zc;
	}
	double fwd = (zcStart-zcEnd)/level;
	double w1  = fabs(mrs1)<K_NEW_DOUBLE_TOL?
					-fwd*((tStart*zcStart-tEnd*zcEnd)/(zcStart-zcEnd)-numbis/level)
				:	fwd/mrs1*((eStart1*zcStart-eEnd1*zcEnd)/(zcStart-zcEnd)-num1/level);
	double w2  = fabs(mrs2)<K_NEW_DOUBLE_TOL?
					-fwd*((tStart*zcStart-tEnd*zcEnd)/(zcStart-zcEnd)-numbis/level)
				:	fwd/mrs2*((eStart2*zcStart-eEnd2*zcEnd)/(zcStart-zcEnd)-num2/level);

	std::vector<double>& pRes = new std::vector<double>(4,0.);
	(*pRes)[0]=fwd;
	(*pRes)[1]=w1;
	(*pRes)[2]=w2;
	(*pRes)[3]=level;
	return pRes;
}

std::vector<double>& csecurity::getRateWeight(double t_start,const std::vector<double>& floatPayTimes,const std::vector<double>& floatNotionals,
		const std::vector<double>& fixPayTimes,const std::vector<double>& fixNotionals,const std::vector<double>& fixPeriods,
		double mrs1,double mrs2)
{
	ARM_ZeroCurvePtr ZcCurve = m_pmodel->GetZeroCurve();

	size_t k;

	/// MRS low levels
	if(fabs(mrs1)<K_NEW_DOUBLE_TOL)
		mrs1 = (mrs1 < 0 ? -K_NEW_DOUBLE_TOL : K_NEW_DOUBLE_TOL);
	if(fabs(mrs2)<K_NEW_DOUBLE_TOL)
		mrs2 = (mrs2 < 0 ? -K_NEW_DOUBLE_TOL : K_NEW_DOUBLE_TOL);

	/// Variable notional float leg
	size_t nbFloatFlows = floatPayTimes.size();
	double t		= t_start/K_YEAR_LEN;
	double zcNtl	= ZcCurve->DiscountPrice(t) * floatNotionals[0];
					
	double floatNum1	= zcNtl * exp(-mrs1*t);
	double floatNum2	= zcNtl * exp(-mrs2*t);
	double floatLeg		= zcNtl;
	for(k=0;k+1<nbFloatFlows;++k)
	{
		t			= floatPayTimes[k]/K_YEAR_LEN;
		zcNtl		= ZcCurve->DiscountPrice(t) * (floatNotionals[k+1]-floatNotionals[k]);
		floatNum1	+= zcNtl * exp(-mrs1*t);
		floatNum2	+= zcNtl * exp(-mrs2*t);
		floatLeg	+= zcNtl;
	}
	t			= floatPayTimes[nbFloatFlows-1]/K_YEAR_LEN;
	zcNtl		= ZcCurve->DiscountPrice(t)	 * floatNotionals[nbFloatFlows-1];
	floatNum1	-= zcNtl * exp(-mrs1*t);
	floatNum2	-= zcNtl * exp(-mrs2*t);
	floatLeg	-= zcNtl;

	/// Variable notional fixed leg
	double fixNum1=0,fixNum2=0,level=0;
	size_t nbFixFlows = fixPayTimes.size();
	for(k=0;k<nbFixFlows;++k)
	{
		t		= fixPayTimes[k]/K_YEAR_LEN;
		zcNtl	= ZcCurve->DiscountPrice(t) * fixNotionals[k] * fixPeriods[k];
		fixNum1	+= zcNtl * exp(-mrs1*t);
		fixNum2	+= zcNtl * exp(-mrs2*t);
		level	+= zcNtl;
	}
	double fwd = floatLeg/level;
	double w1  = fwd/mrs1*(floatNum1/floatLeg-fixNum1/level);
	double w2  = fwd/mrs2*(floatNum2/floatLeg-fixNum2/level);

	std::vector<double>& pRes = new std::vector<double>(4,0.);
	(*pRes)[0]=fwd;
	(*pRes)[1]=w1;
	(*pRes)[2]=w2;
	(*pRes)[3]=0.01*level; /// to be consistent with the constant notional getRateWeight()
	return pRes;
}


std::vector<double>& csecurity::getChangeNumWeight(double t_pay,double t_start,double t_end,const std::vector<double>& payTimes,const std::vector<double>& periods,double mrs1,double mrs2)
{
	ARM_ZeroCurvePtr ZcCurve = m_pmodel->GetZeroCurve();

	double tPay			= t_pay/K_YEAR_LEN;
	double tStart		= t_start/K_YEAR_LEN;
	double tEnd			= t_end/K_YEAR_LEN;
	double zcStart		= ZcCurve->DiscountPrice(tStart);
	double zcEnd		= ZcCurve->DiscountPrice(tEnd);
	double zcPay		= ZcCurve->DiscountPrice(tPay);
	double eStart1		= exp(-mrs1*tStart);
	double eEnd1		= exp(-mrs1*tEnd);
	double eStart2		= exp(-mrs2*tStart);
	double eEnd2		= exp(-mrs2*tEnd);

	double ePay1		= exp(-mrs1*tPay);
	double ePay2		= exp(-mrs2*tPay);
					
	double t,zc;
	double num1=0,num2=0,level=0,numbis=0.;
	for (size_t k=0;k<payTimes.size();k++)
	{
		t		= payTimes[k]/K_YEAR_LEN;
		zc		= ZcCurve->DiscountPrice(t);
		num1	+= exp(-mrs1*t)*periods[k]*zc;
		num2	+= exp(-mrs2*t)*periods[k]*zc;
		numbis	+= periods[k]*zc*t;
		level	+= periods[k]*zc;
	}
	double w1  = fabs(mrs1)<K_NEW_DOUBLE_TOL?
					-(tPay-numbis/level)
				:	1./mrs1*(ePay1-num1/level);
	double w2  = fabs(mrs2)<K_NEW_DOUBLE_TOL?
					-(tPay-numbis/level)
				:	1./mrs2*(ePay2-num2/level);

	std::vector<double>& pRes = new std::vector<double>(2,0.);
	(*pRes)[0]=w1;
	(*pRes)[1]=w2;
	return pRes;
}

cswaption::cswaption(ARM_PricingModel* model, double price, ARM_VanillaSwaptionArg* argSW, double mrs1, double mrs2)
	:csecurity(model,price)
{
	m_reset		= argSW->GetResetTime();
	m_callPut	= argSW->GetCallPut();

	std::vector<double>& rateWeight;
	if(argSW->GetIsConstantNotional())
		rateWeight = getRateWeight(	argSW->GetStartTime(),
									argSW->GetEndTime(),
									*argSW->GetFixPayTimes(),
									*argSW->GetFixPayPeriods(),
									mrs1,mrs2);
	else
		rateWeight = getRateWeight(	argSW->GetStartTime(),
									*argSW->GetFloatEndTimes(),
									*argSW->GetFloatNotional(),
									*argSW->GetFixPayTimes(),
									*argSW->GetFixNotional(),
									*argSW->GetFixPayPeriods(),
									mrs1,mrs2);

	m_fwd = rateWeight->Elt(0);
	m_w1 = rateWeight->Elt(1);
	m_w2 = rateWeight->Elt(2);
	m_strike = argSW->GetStrikes()->Elt(0);
	m_level = rateWeight->Elt(3);
	delete rateWeight;
}

void cswaption::BuildCoeffs(const std::vector<double>& lv, const std::vector<double>& lv0, const std::vector<double>& lv1)
{
	m_A  = lv[0]*m_w1*m_w1;
	m_B  = lv[1]*m_w2*m_w2;
	m_C  = 2.*lv[2]*m_w1*m_w2;

	double cov	= lv0[0]*m_w1*m_w1 + lv0[1]*m_w2*m_w2 + 2*lv0[2]*m_w1*m_w2;

	bool success(true);
	double mktvol = VanillaImpliedVol_N(m_fwd,m_price/100/m_level,m_strike,m_reset/K_YEAR_LEN,m_callPut,NULL, &success);
	if(success)
		m_D = mktvol*mktvol*m_reset/K_YEAR_LEN-cov;
	else
		m_D = 0.0;
}

cspreadoption::cspreadoption(ARM_PricingModel* pmodel, double price, ARM_VanillaSpreadOptionArg* argSO, double mrs1, double mrs2,ARM_Currency* ccy)
	:csecurity(pmodel,price)
{
/*****
	m_reset		= argSO->GetResetTimes()->Elt(0);
	m_levLong	= argSO->GetCoeffLong()->Elt(0);
	m_levShort	= argSO->GetCoeffShort()->Elt(0);
	m_callPut	= argSO->GetCallPut();

	double payTime					= argSO->GetPayTimes()->Elt(0);
	double swapLongFloatStartTime	= argSO->GetSwapLongFloatStartTime()->Elt(0);
	double swapLongFloatEndTime		= argSO->GetSwapLongFloatEndTime()->Elt(0);
	double swapShortFloatStartTime	= argSO->GetSwapShortFloatStartTime()->Elt(0);
	double swapShortFloatEndTime	= argSO->GetSwapShortFloatEndTime()->Elt(0);

	const std::vector<double>& swapLongFixPayTimes	= *(argSO->GetSwapLongFixPayTimes()[0]);
	const std::vector<double>& swapLongFixPayPeriods	= *(argSO->GetSwapLongFixPayPeriods()[0]);
	const std::vector<double>& swapShortFixPayTimes	= *(argSO->GetSwapShortFixPayTimes()[0]);
	const std::vector<double>& swapShortFixPayPeriods	= *(argSO->GetSwapShortFixPayPeriods()[0]);
	

	std::vector<double>& rateShortWeight = getRateWeight(	swapShortFloatStartTime, 
													swapShortFloatEndTime,
													swapShortFixPayTimes,
													swapShortFixPayPeriods,
													mrs1,mrs2);
			
	std::vector<double>& chNumShortWeight = getChangeNumWeight(	payTime,
															swapShortFloatStartTime, 
															swapShortFloatEndTime,
															swapShortFixPayTimes,
															swapShortFixPayPeriods,
															mrs1,mrs2);


	std::vector<double>& rateLongWeight = getRateWeight(	swapLongFloatStartTime, 
													swapLongFloatEndTime,
													swapLongFixPayTimes, 
													swapLongFixPayPeriods,
													mrs1,mrs2);
			
	std::vector<double>& chNumLongWeight = getChangeNumWeight(	payTime,
															swapLongFloatStartTime, 
															swapLongFloatEndTime,
															swapLongFixPayTimes, 
															swapLongFixPayPeriods,
															mrs1,mrs2);

	m_fwdLong = rateLongWeight->Elt(0);
	m_fwdShort = rateShortWeight->Elt(0);

	m_wLong1 = rateLongWeight->Elt(1);
	m_wLong2 = rateLongWeight->Elt(2);
	m_wShort1 = rateShortWeight->Elt(1);
	m_wShort2 = rateShortWeight->Elt(2);

	m_cLong1 = chNumLongWeight->Elt(0);
	m_cLong2 = chNumLongWeight->Elt(1);
	m_cShort1 = chNumShortWeight->Elt(0);
	m_cShort2 = chNumShortWeight->Elt(1);

	m_strike = argSO->GetStrikes()->Elt(0);
	m_level  = m_pmodel->GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN)*argSO->GetPayPeriods()->Elt(0);

	delete rateShortWeight;
	delete chNumShortWeight;
	delete rateLongWeight;
	delete chNumLongWeight;
*****/
	m_reset		= argSO->GetResetTimes()->Elt(0);
	m_levLong	= argSO->GetCoeffLong()->Elt(0);
	m_levShort	= argSO->GetCoeffShort()->Elt(0);
	m_callPut	= argSO->GetCallPut();

	double payTime					= argSO->GetPayTimes()->Elt(0);
	double swapLongFloatStartTime	= argSO->GetSwapLongFloatStartTime()->Elt(0);
	double swapLongFloatEndTime		= argSO->GetSwapLongFloatEndTime()->Elt(0);
	double swapShortFloatStartTime	= argSO->GetSwapShortFloatStartTime()->Elt(0);
	double swapShortFloatEndTime	= argSO->GetSwapShortFloatEndTime()->Elt(0);
	
			
	double asOf = m_pmodel->GetZeroCurve()->GetAsOfDate().GetJulian();
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

	std::vector<double>& rateShortWeight = getRateWeight(	shortSwapRate->floatStartTime, 
													shortSwapRate->floatEndTime,
													shortSwapRate->fixPayTimes,
													shortSwapRate->fixPayPeriods,
													mrs1,mrs2);
			
	std::vector<double>& chNumShortWeight = getChangeNumWeight(	payTime,
															shortSwapRate->floatStartTime, 
															shortSwapRate->floatEndTime,
															shortSwapRate->fixPayTimes,
															shortSwapRate->fixPayPeriods,
															mrs1,mrs2);


	std::vector<double>& rateLongWeight = getRateWeight(	longSwapRate->floatStartTime, 
													longSwapRate->floatEndTime,
													longSwapRate->fixPayTimes, 
													longSwapRate->fixPayPeriods,
													mrs1,mrs2);
			
	std::vector<double>& chNumLongWeight = getChangeNumWeight(	payTime,
															longSwapRate->floatStartTime, 
															longSwapRate->floatEndTime,
															longSwapRate->fixPayTimes, 
															longSwapRate->fixPayPeriods,
															mrs1,mrs2);

	m_fwdLong = rateLongWeight->Elt(0);
	m_fwdShort = rateShortWeight->Elt(0);

	m_wLong1 = rateLongWeight->Elt(1);
	m_wLong2 = rateLongWeight->Elt(2);
	m_wShort1 = rateShortWeight->Elt(1);
	m_wShort2 = rateShortWeight->Elt(2);

	m_cLong1 = chNumLongWeight->Elt(0);
	m_cLong2 = chNumLongWeight->Elt(1);
	m_cShort1 = chNumShortWeight->Elt(0);
	m_cShort2 = chNumShortWeight->Elt(1);

	m_strike = argSO->GetStrikes()->Elt(0);
	m_level  = m_pmodel->GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN)*argSO->GetPayPeriods()->Elt(0);

	delete rateShortWeight;
	delete chNumShortWeight;
	delete rateLongWeight;
	delete chNumLongWeight;

}

void cspreadoption::BuildCoeffs(const std::vector<double>& lv, const std::vector<double>& lv0, const std::vector<double>& lv1)
{
	double conv1 = lv1[0]*m_cLong1*m_wLong1
				 + lv1[1]*m_cLong2*m_wLong2
				 + lv1[2]*m_cLong1*m_wLong2
				 + lv1[2]*m_cLong2*m_wLong1;
					
	double conv2 = lv1[0]*m_cShort1*m_wShort1
				 + lv1[1]*m_cShort2*m_wShort2
				 + lv1[2]*m_cShort1*m_wShort2
				 + lv1[2]*m_cShort2*m_wShort1;

	double fwd	 = m_levLong*(m_fwdLong+conv1)-m_levShort*(m_fwdShort+conv2);

	double cov   =		lv0[0]*(m_levLong*m_wLong1-m_levShort*m_wShort1)*(m_levLong*m_wLong1-m_levShort*m_wShort1)
				 +		lv0[1]*(m_levLong*m_wLong2-m_levShort*m_wShort2)*(m_levLong*m_wLong2-m_levShort*m_wShort2)
				 + 2*	lv0[2]*(m_levLong*m_wLong1-m_levShort*m_wShort1)*(m_levLong*m_wLong2-m_levShort*m_wShort2);

	m_A  =		lv[0]*(m_levLong*m_wLong1-m_levShort*m_wShort1)*(m_levLong*m_wLong1-m_levShort*m_wShort1);
	m_B  =		lv[1]*(m_levLong*m_wLong2-m_levShort*m_wShort2)*(m_levLong*m_wLong2-m_levShort*m_wShort2);
	m_C  = 2.*	lv[2]*(m_levLong*m_wLong1-m_levShort*m_wShort1)*(m_levLong*m_wLong2-m_levShort*m_wShort2);

	bool success(true);
	double targetPrice = m_price/100/m_level;
	double vi = m_callPut * (fwd-m_strike);
	vi = vi < 0 ? 0 : vi;

	/// Take care if model convexity is inconsistent with market price
	m_D = 0.0;
	if(targetPrice > vi)
	{
		double mktvol = VanillaImpliedVol_N(fwd,targetPrice,m_strike,m_reset/K_YEAR_LEN,m_callPut,NULL, &success);
		if(success)
			m_D = mktvol*mktvol*m_reset/K_YEAR_LEN-cov;
	}
}

#undef FIRST_STATE_VARIABLE
#undef SECOND_STATE_VARIABLE

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
