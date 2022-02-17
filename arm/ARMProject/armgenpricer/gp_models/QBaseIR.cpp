/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file QBaseIR.cpp
 *
 *  \brief base class for the q model
 *
 *	\author  E. Benhamou
 *	\version 1.0
 *	\date December 2004
 */

#include "gpmodels/QBaseIR.h"

/// gpbase
#include "gpbase/curve.h"

/// gpinfra
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingmodeltype.h"

/// gpclosedforms
#include "gpclosedforms/gaussian_integrals.h"

/// kernel
#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

const size_t GL_PY_NBPOINTS = 4;

const double ARM_QModelBaseIR::RatePeriod=0.25;           /// 3M
const double ARM_QModelBaseIR::RateUpperBound=5.0;        /// 500%
const double ARM_QModelBaseIR::Rate0LowerBound=0.0020;    /// 20bp



////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: CopyNoCleanUp
///	Returns: 
///	Action : Arguments copy
////////////////////////////////////////////////////
void ARM_QModelBaseIR::CopyNoCleanUp(const ARM_QModelBaseIR& rhs)
{
    itsDegenerateInHW       = rhs.itsDegenerateInHW;
    itsQDfTargetMap         = rhs.itsQDfTargetMap;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelBaseIR::ARM_QModelBaseIR(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params, bool degenerateInHW ) 
:	ARM_QModelBase(zc,params), itsDegenerateInHW( degenerateInHW )
{
	if( itsDegenerateInHW )
	{
		double qParameter = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		if( fabs(qParameter) > K_NEW_DOUBLE_TOL)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": qPArameter != 0 && flag to degenerate in HW On!");
	}
}



////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_QModelBaseIR::ARM_QModelBaseIR(const ARM_QModelBaseIR& rhs)
:	ARM_QModelBase(rhs)
{
    CopyNoCleanUp(rhs);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_QModelBaseIR& ARM_QModelBaseIR::operator=(const ARM_QModelBaseIR& rhs)
{
	if(this != &rhs)
	{
		ARM_QModelBase::operator=(rhs);
        CopyNoCleanUp(rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_QModelBaseIR::~ARM_QModelBaseIR()
{}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: MappingFunction
///	Returns :
///	Action  : mapping function of the Q model
////////////////////////////////////////////////////

double ARM_QModelBaseIR::MappingFunction( double x, double x0, double q0 ) const
{
	double r;

	if( itsDegenerateInHW )
	{
#if defined(__GP_STRICT_VALIDATION)
		double qParameter = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		if( fabs(qParameter) > K_NEW_DOUBLE_TOL)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": qPArameter != 0 && flag to degenerate in HW On!");
#endif
		r=x;
	}
	else
	{
		if(fabs(q0)<K_NEW_DOUBLE_TOL)
			r=x0*(1+x);
		else if( x0 > ARM_QModelBaseIR::Rate0LowerBound)
			r=x0*(1.0+(exp(q0*x)-1.0)/q0);
		else
			r=ARM_QModelBaseIR::Rate0LowerBound*(1.0+(exp(q0*x)-1.0)/q0);

		if(r > ARM_QModelBaseIR::RateUpperBound)
			r = ARM_QModelBaseIR::RateUpperBound;
	}
	return r;
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: ComputeFwdAtTime
///	Returns: double 
///	Action : computes a given short rate at the time evalTime
////////////////////////////////////////////////////
double ARM_QModelBaseIR::ComputeFwdAtTime( double evalTime ) const
{
	if( itsDegenerateInHW )
	{
#if defined(__GP_STRICT_VALIDATION)
		double qParameter = ((ARM_CurveModelParam&) GetModelParams()->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);
		if( fabs(qParameter) > K_NEW_DOUBLE_TOL)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": qPArameter != 0 && flag to degenerate in HW On!");
#endif
		return 1.0;
	}
	else
	{
		ARM_ZeroCurvePtr zc	= GetZeroCurve();
		double df1			= zc->DiscountPrice(evalTime/K_YEAR_LEN);
		double df2			= zc->DiscountPrice(evalTime/K_YEAR_LEN+ARM_QModelBaseIR::RatePeriod);
		return (df1/df2-1.0)/ARM_QModelBaseIR::RatePeriod;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: DecayFactor
///	Returns : double
///	Action  : Compute the decay factor between [a,b] :
///           integral s=a->b of f0(s)*exp(-mrs*(s-T)) 
////////////////////////////////////////////////////
double ARM_QModelBaseIR::DecayFactor(double a, double b, double T, const ARM_ModelParams* modelParams ) const
{
    double mrs = modelParams->GetModelParam( ARM_ModelParamType::MeanReversion).ToCurveModelParam().GetCurve()->GetOrdinate(0);

    double decay;
	if( itsDegenerateInHW )
	{
        if(fabs(mrs)>K_NEW_DOUBLE_TOL)
            decay = (exp(-mrs*(a-T))-exp(-mrs*(b-T)))/mrs;
        else
            decay = b-a;
    }
    else
    {
        /// Compute a Gauss-Legendre discretization
        size_t i,nbPoints=static_cast<int>(ceil(b-a)*GL_PY_NBPOINTS);
        GaussLegendre_Coefficients GLPoints(nbPoints,a,b);

        /// Gauss-Legendre summation
        double t,f0;
        decay=0.0;
        for(i=0;i<nbPoints;++i)
        {
            t       = GLPoints.get_point(i);
            f0      = ComputeFwdAtTime( t*K_YEAR_LEN );
            decay  += f0 * exp(-mrs*(t-T)) * GLPoints.get_weight(i);
        }
    }

    return decay;
}



////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: PartialIntegratedDrift
///	Returns : double
///	Action  : Compute the integral between [a,b]
///         of sigma(s)^2 * exp(-lambda * (T-s) ) * DecayFactor(s,T,s)
///         volatility is cst in the interval
////////////////////////////////////////////////////

double ARM_QModelBaseIR::PartialIntegratedDrift(int i, double mrs, double a, double b, double T, const ARM_ModelParams* modelParams ) const
{
    /// Compute a Gauss-Legendre discretization
    size_t nbPoints=static_cast<int>(ceil(b-a)*GL_PY_NBPOINTS);
    GaussLegendre_Coefficients GLPoints(nbPoints,a,b);

    double GLPoint;
    double sum = 0.0;
    for(size_t j=0;j<nbPoints;++j)
    {
        GLPoint = GLPoints.get_point(j);
        sum    += exp(-mrs*(T-GLPoint) ) * DecayFactor(GLPoint,T,GLPoint, modelParams);
    }
    double vol = ((ARM_CurveModelParam&) modelParams->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetOrdinate(i);
    sum *= vol*vol;
    return sum;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: IntegratedDrift
///	Returns : double
///	Action  : Compute the integral between [0,t]
///         of sigma(s)^2 * exp(-mrs * (T-s) ) * DecayFactor(s,T,s)
///         use breakpoint times of the volatility
////////////////////////////////////////////////////

double ARM_QModelBaseIR::IntegratedDrift(double t,double T, const ARM_ModelParams* modelParams ) const
{
    t /= K_YEAR_LEN;
    T /= K_YEAR_LEN;

	const std::vector<double>& BreakPointTimes = ((ARM_CurveModelParam&) modelParams->GetModelParam( ARM_ModelParamType::QVol)).GetCurve()->GetAbscisses();
    std::vector<double> MaturityTimes( BreakPointTimes );
    MaturityTimes /= K_YEAR_LEN;

    int position = CC_NS(std,lower_bound)(MaturityTimes.begin(), MaturityTimes.end(), t) - MaturityTimes.begin()-1;
    double mrs = ((ARM_CurveModelParam&) modelParams->GetModelParam( ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinate(0);

    double sumTotal = 0.0;
    if( -1 == position )
    {
        sumTotal = PartialIntegratedDrift(0,mrs,0,t,T, modelParams);
    }
    else
    {
        sumTotal += PartialIntegratedDrift(0,mrs,0,MaturityTimes[0],T,modelParams);
        for(int i=1; i<position; ++i )
            sumTotal += PartialIntegratedDrift(i,mrs,MaturityTimes[i],MaturityTimes[i+1],T,modelParams);
        sumTotal += PartialIntegratedDrift(position,mrs,MaturityTimes[i],t,T,modelParams);
    }

    return sumTotal;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: SuperDecayFactor
///	Returns : double
///	Action  : Compute the integral between [a,b] :
///           integral s=a->b of f0(s)*exp(-mrs*(s-T))
///           *q*Integral of drift 0 -> s
////////////////////////////////////////////////////
double ARM_QModelBaseIR::SuperDecayFactor(double a, double b, double T, const ARM_ModelParams* modelParams ) const
{
    double qParameter	= ((ARM_CurveModelParam&) modelParams->GetModelParam( ARM_ModelParamType::QParameter)).GetCurve()->GetOrdinate(0);

    if( fabs(qParameter) < K_NEW_DOUBLE_TOL )
        return DecayFactor(a,b,T,modelParams);
    else
    {
        double mrs = ((ARM_CurveModelParam&) modelParams->GetModelParam( ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinate(0);

        /// Compute a Gauss-Legendre discretization
        size_t i,nbPoints=static_cast<int>(ceil(b-a)*GL_PY_NBPOINTS);
        GaussLegendre_Coefficients GLPoints(nbPoints,a,b);

        /// Gauss-Legendre summation
        double t,decay = 0.0,f0;
        for(i=0;i<nbPoints;++i)
        {
            t       = GLPoints.get_point(i);
            f0      = ComputeFwdAtTime( t*K_YEAR_LEN );
            decay  += f0 * exp(-mrs*(t-T)) * exp( qParameter * IntegratedDrift(a/2,t,modelParams) ) * GLPoints.get_weight(i);
        }
        return decay;
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: VolZc
///	Returns : double
///	Action  : Compute the Zc instantaneous volatility  :
///           volZc(t,T)=-sigma(t)*decay(t,T)
////////////////////////////////////////////////////
double ARM_QModelBaseIR::VolZc(double t, double T) const
{
    const ARM_ModelParams* modelParams = GetModelParams();
    double sigma    = modelParams->GetModelParam(ARM_ModelParamType::QVol).ToCurveModelParam().GetCurve()->Interpolate(t);
    double yft=t/K_YEAR_LEN;
    double decay    = DecayFactor(yft,T/K_YEAR_LEN,yft,modelParams);

    return -sigma*decay;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: MarkovianDrift
///	Returns : void
///	Action  : No Markovian drift
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_QModelBaseIR::MarkovianDrift( size_t timeIdx, const ARM_GP_MatrixPtr& numMethodStates ) const
{
	return ARM_GP_MatrixPtr( new ARM_GP_Matrix( numMethodStates->rows() , numMethodStates->cols(), 0.0 ) );
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: Libor
///	Returns : a vector of libor values
///	Action  : Default Libor computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelBaseIR::Libor( 
		const string& curveName, 
		double evalTime,
		double fwdStartTime,
		double fwdEndTime,
		double period,
		double fwdResetTime,
		double payTime,
		const ARM_PricingStatesPtr& states) const
{    
    /// Get Libor Zc through the fixing functor
    return DefaultLibor(curveName,evalTime,fwdStartTime,fwdEndTime,period,fwdResetTime,payTime,states);
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: Annuity
///	Returns : a vector of annuity
///	Action  : Default Annuity Computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelBaseIR::Annuity(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const ARM_PricingStatesPtr& states) const
{
     return DefaultAnnuity(curveName,evalTime,fixPayTimes,fixPayPeriods,states);
}



////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: SwapRateInPlaceWithComputedAnnuity
///	Returns : a vector of swap rate values
///	Action  : Default Swap Rate computation
///           using double notional method
///				WARNING: need to clone the annuity as the computation
///				is done in place... I REPEAT
///				the annuity argument will be modified by the computation as it is in place
///				so to avoid side effect on the annuity (or if you want to keep the value of the annuity)
///				one needs to first clone the value of the annuity
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelBaseIR::SwapRateInPlaceWithComputedAnnuity(
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
		const ARM_VectorPtr& FixedComputedAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
		const ARM_PricingStatesPtr& states ) const
{
	ARM_GP_VectorPtr FloatComputedAnnuity = ARM_GP_VectorPtr( NULL );
    if( isDbleNotional && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
		FloatComputedAnnuity = Annuity( curveName, evalTime, floatPayTimes, floatPayPeriods, states);

	return DefaultSwapRateInPlaceWithComputedAnnuity( curveName, evalTime, floatStartTime, floatEndTime, 
		fixPayTimes, fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods,
		margin, isDbleNotional, FixedComputedAnnuity, FloatComputedAnnuity, states);
}


////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: NPVSwap
///	Returns: a vector of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwap computation
///           using double notional method
/// 
////////////////////////////////////////////////////
ARM_VectorPtr ARM_QModelBaseIR::NPVSwap(
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
		const std::vector<double>& FixNotional,
		const std::vector<double>& FloatNotional,
		const ARM_GP_Matrix& strikesPerState,
		int payRec,
		const ARM_PricingStatesPtr& states) const
{
	ARM_VectorPtr floatAnnuity = ARM_VectorPtr (NULL);
    if( isDbleNotional && GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) && margin[0] != 0.0)
		floatAnnuity = Annuity( curveName, evalTime, floatPayTimes, floatPayPeriods, states);
    
	return DefaultNPVSwapWithComputedAnnuity( curveName, evalTime, floatStartTime, floatEndTime, fixPayTimes,
		fixPayPeriods, fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods, 
		margin, isDbleNotional, FixNotional,FloatNotional, strikesPerState, payRec, floatAnnuity,states);
}

////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: DefaultNPVSwapLeg
///	Returns: a matrix of NPVSwap(t,F(R,Ti),K)
///	Action : 
/// Default: Default NPVSwalLeg computation
////////////////////////////////////////////////////
ARM_GP_MatrixPtr ARM_QModelBaseIR::NPVSwapLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fwdStartTimes, 
		const std::vector<double>& fwdEndTimes, 
		const std::vector<double>& fwdPayPeriods, 
		const std::vector<double>& payTimes, 
		const std::vector<double>& payPeriods, 
		const std::vector<double>& margin, 
		const std::vector<double>& notional, 
		const ARM_PricingStatesPtr& states) const
{
	ARM_GP_MatrixPtr result= DefaultNPVSwapLeg(curveName, evalTime,fwdStartTimes,fwdEndTimes, fwdPayPeriods, 
		 payTimes, payPeriods, margin, notional,  states); 

	return result;
}
////////////////////////////////////////////////////
///	Class   : ARM_ForwardMarginBasis
///	Routines: NPVFixLeg
///	Returns : matrix ptr
///	Action  : Fix leg computing
ARM_GP_MatrixPtr ARM_QModelBaseIR::NPVFixLeg(
		const string& curveName, 
		double evalTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		const std::vector<double>& FixNotional,
		const ARM_GP_Matrix& strikesPerState,
		int   payRec,
		const ARM_PricingStatesPtr& states) const
{
	ARM_GP_MatrixPtr result= DefaultNPVFixLeg(curveName,evalTime,fixPayTimes,fixPayPeriods,
		FixNotional, strikesPerState,payRec,states);
	
	return result;
}
///////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: ImpliedVol
///	Returns: double
///	Action : To Calculate the Implied Volatility
///  By defaut using BS Formula
////////////////////////////////////////////////////
double ARM_QModelBaseIR::ImpliedVol(const ARM_VanillaArg& arg) const
{
    CC_Ostringstream os;
	os << ARM_USERNAME << " : Function is not implimented yet";
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
}
////////////////////////////////////////////////////
///	Class  : ARM_QModelBaseIR
///	Routine: GetType
///	Returns: int
///	Action : tells the type of the model
////////////////////////////////////////////////////
int ARM_QModelBaseIR::GetType() const
{
	return MT_EQUITY_MODEL;
}


////////////////////////////////////////////////////
///	Class   : ARM_QModelBaseIR
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_QModelBaseIR::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
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
    case ARM_ModelParamType::QVol:
        {
            double date = portfolio->GetAsset(0)->GetExpiryDate().GetJulian() - asOfDate;
            tmpdates.push_back(date);
            for(i=1; i<size1; i++) 
            {
                double resetlag = portfolio->GetAsset(i)->GetExpiryDate().GetJulian() - asOfDate;
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


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

