/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file MSV1F.cpp
 *  \brief Markov Stochastic Volatility 1 factor model
 *
 *	\author  A Triki
 *	\version 1.0
 *	\date October 2005
 */

/// this header comes first as it include some preprocessor constants
#include "gpbase/removeidentifiedwarning.h"
#include "gpmodels/FRMSV1F.h"
#include "gpmodels/riccatimsv.h"
#include "gpmodels/ModelParamsSFRM.h"
#include "gpmodels/VanillaSwaptionArgSFRM.h"
#include "gpmodels/ModelParamsFRMSV1F.h"
#include "gpmodels/VolDiffusionParams.h"

/// gpbase
#include "gpbase/ostringstream.h"
#include "gpbase/gpmatrixtriangular.h"
#include "gpbase/curve.h"
#include "gpbase/interpolatorvector.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/datestrip.h"


/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/modelparams.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/nummethod.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/discretisationscheme.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/bootstrapnd.h"
#include "gpcalib/vanillaarg.h"


/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/heston_interface.h"
#include "gpclosedforms/riccati.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/rungekutta.h"
#include "gpnumlib/solver.h"
#include "gpnumlib/newtonraphson.h"
#include "gpnumlib/numfunction.h"


//gpmodel
#include "gpmodels/AnalyticIRModel.h"

/// kernel
#include <inst/portfolio.h>


CC_BEGIN_NAMESPACE( ARM )


#define FIRST_STATE_VARIABLE    0 // for Xt
#define SECOND_STATE_VARIABLE   1 // for Yt
#define THIRD_STATE_VARIABLE    2 // for Zt (Volatility Variable)

///Constants for Heston Closed Form Pricng Formula
const double DEFAULT_NBOSCILLATIONS		    = 3; 
const double DEFAULT_PRECISION				= 1.e-6;
const double DEFAULT_NBSTAGE				= -1;
const double DEFAULT_NBPOINTS_STAGE1		= 120;
const double DEFAULT_NBPOINTS_PER_STAGE		= 30;

///For Numerical Derivatives Calculation
const double DERIVATIVE_EPS = 1.e-6;
///For MC simulation
const double FWD_MC_MAX = 1;

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_FRMSV1F::ARM_FRMSV1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params)
:ARM_SVModels(zc,params)
{	
}

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_FRMSV1F::ARM_FRMSV1F(const ARM_FRMSV1F& rhs)
: ARM_SVModels(rhs) 
{
}


/////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_FRMSV1F::~ARM_FRMSV1F( )
{
}


////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_FRMSV1F& ARM_FRMSV1F::operator=(const ARM_FRMSV1F& rhs)
{
	if(this != &rhs)
	{
		ARM_SVModels::operator=(rhs);	
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_FRMSV1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "FRM Stochastic Volatility Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_FRMSV1F::Clone() const
{
	return new ARM_FRMSV1F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_FRMSV1F::PreProcessing(ARM_ModelFitter& modelFitter)
{ 
}

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_FRMSV1F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
}

////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_FRMSV1F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_FRMSV1F::ValidateModelParams(const ARM_ModelParams& params) const
{
	return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_FRMSV1F::IntegratedLocalDrifts(
	const std::vector<double>& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
}



////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_FRMSV1F::ModelStateLocalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
}


////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of ND matrix
///	Action : Compute local variances of the state variable 
/// between each time step
////////////////////////////////////////////////////
void ARM_FRMSV1F::NumMethodStateLocalVariances(
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
	size_t i,j;

	size_t factorNb = FactorCount();

	ARM_GP_Matrix identity(factorNb,factorNb);

	for (i = 0; i < factorNb; ++i)
	{
		for(j = 0; j < i; ++j)
		{
			identity(i,j) = identity(j,i) = 0.0;
		}
		identity(i,i) = 1.0;
	}

	// All the variance is in the numerical method

	for(i=0;i<nbSteps-1;++i)
	{
		nextStep=timeSteps[i+1];
		/// [i] => local variance from ti->ti+1
		localVariances[offsetIndex+i] = static_cast<ARM_GP_Matrix*>(identity.Clone());
		step=nextStep;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: NumMethodStatesGlobalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute global of the state 
/// variable between each time step
////////////////////////////////////////////////////
void ARM_FRMSV1F::NumMethodStateGlobalVariances(
    const std::vector<double>& timeSteps,
    ARM_MatrixVector& variances) const
{
	size_t nbSteps	= timeSteps.size();
	size_t modelNb	= GetModelNb();
    double step		= timeSteps[0],
		   nextStep;
	size_t offsetIndex2	= nbSteps*modelNb;

#if defined(__GP_STRICT_VALIDATION)
	if( variances.size()!= offsetIndex2 ) 
		ARM_THROW( ERR_INVALID_ARGUMENT, "localDrifts.size() != offsetIndex" );
#endif

	variances.resize(nbSteps*(modelNb+1));
	size_t i,j;

	size_t factorNb = FactorCount();

	ARM_GP_Matrix identity(factorNb,factorNb);

	for (i = 0; i < factorNb; ++i)
	{
		for(j = 0; j < i; ++j)
		{
			identity(i,j) = identity(j,i) = 0.0;
		}
		identity(i,i) = 1.0;
	}

	/// fills the variance
    variances[offsetIndex2+0]=static_cast<ARM_GP_Matrix*>(identity.Clone());

    for(i=0;i<nbSteps-1;++i)
    {
        nextStep=timeSteps[i+1];
        
		/// [i+1] => variance from 0 -> ti+1
        /// we can't sum up local variance !
        variances[offsetIndex2+i+1] = static_cast<ARM_GP_Matrix*>(identity.Clone());
        step=nextStep;
    }

}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: VarianceToTime
///	Returns: a time
///	Action : Compute the time such that
///          var(t)=var
////////////////////////////////////////////////////
double ARM_FRMSV1F::VarianceToTime(double var,double minTime,double maxTime) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_MarkovSV Model!");

}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: DiscountFactor
///	Returns: a vector of Zc(t,T)
///	Action : Closed form formula for DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_FRMSV1F::DiscountFactor( 
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

	size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	return ARM_VectorPtr( new std::vector<double>(payoffSize,zcT) );

}

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_FRMSV1F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
{
#if defined(__GP_STRICT_VALIDATION)
	if(timeIndex<0)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index is negative!");
	if( timeIndex >= GetNumMethod()->GetTimeSteps()->size()-1 )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "time index bigger than max size!");
#endif

	const ARM_MatrixVector& localVar	= GetModelStateLocalVars();
	const ARM_MatrixVector& localStdDev = GetModelStateLocalStdDevs();

	double time = GetNumMethod()->GetTimeStep(timeIndex);
	double nextTime = GetNumMethod()->GetTimeStep(timeIndex+1);
	double dt = (nextTime - time)/K_YEAR_LEN;
	size_t factorsNb= FactorCount();
	size_t statesNb = states->size();
	double X_State,Y_State,Z_State;
	size_t modelNb	= GetModelNb();

	const ARM_VectorPtr EtaVector = ARM_VectorPtr(NULL);//ComputeEtaStates(GetNumMethod()->GetTimeStep(timeIndex), states);

	////// We have constant Mean Reversion
	double Lambda_t = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
	double VolMeanReversion_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).GetValueAtPoint(0);
	double LongTermVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).GetValueAtPoint(0);


	/// We have to Interpolate
	double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(time);
	double VolOfVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(time);
	

	// We have to do the maximum of calculus outside the following routine
	double VolOfVol_dt = VolOfVol_t*sqrt(dt);
	double VolMeanReversion_dt = VolMeanReversion_t * dt;

	///// Varaibles for the calculation of the discretisation of Zt using the Andreasen method

	///// var_1 = exp(-Beta * dt) 
	///// var_2 = 1-exp(-Beta * dt) 
	///// var_3 = Epsilon²/Beta
	///// var_4 = [exp(-Beta * dt) - exp(-2 Beta * dt) ] * Epsilon²/Beta
	///// var_5 = [1-exp(-Beta * dt)]* Epsilon²/2Beta
	///// VolOfVol_2_t = Epsilon²

	double var_1 = exp(-VolMeanReversion_dt);
	double var_2 = 1 - var_1;
	double VolOfVol_2_t = VolOfVol_t*VolOfVol_t;
	double var_3 = VolOfVol_2_t/VolMeanReversion_t;
	double var_4 = (var_1 - var_1 * var_1) * var_3;
	double var_5 = var_2 * var_3 * 0.5;
	double var_6 = 1 - var_1 * var_1;
	double var_7 = 0.5 * var_3 * var_6;


	for( size_t i=0;i<statesNb; ++i )
	{
		X_State = states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE);
		Y_State = states->GetModelState(i,modelNb+SECOND_STATE_VARIABLE);
		Z_State = states->GetModelState(i,modelNb+THIRD_STATE_VARIABLE);

		double Eta = (*EtaVector)[i];
		/// First Variable Xt
		double gaussian_1 = states->GetNumMethodState(i,modelNb+FIRST_STATE_VARIABLE);
		X_State = X_State + (Y_State - Lambda_t * X_State)*dt + sqrt(dt)*gaussian_1*Eta ;
		states->SetModelState(i,modelNb+FIRST_STATE_VARIABLE,X_State);

		/// Second Variable Yt
		Y_State = Y_State + (Eta*Eta - 2*Lambda_t*Y_State)*dt;
		states->SetModelState(i,modelNb+SECOND_STATE_VARIABLE,Y_State);

		/// Third Variable Zt
		////Standard Euler Discretisation approach
		
		double gaussian_2 = states->GetNumMethodState(i,modelNb+SECOND_STATE_VARIABLE);
		
		/// Third Variable Zt
		//// Approach 3 : Old Andreasen
		double Z_State_bis = Z_State*var_1 + var_2;	
		double var_i = 1 + 1/(Z_State_bis*Z_State_bis)* var_7 * Z_State;
		double Nu_2 = log(var_i);
		double Nu = sqrt(Nu_2);

		Z_State = Z_State_bis*exp(-0.5*Nu_2+ Nu*gaussian_2);
		states->SetModelState(i,modelNb+THIRD_STATE_VARIABLE,Z_State);


	}
}

////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
///			 If the swaption is standard, this method
///          calls ARM_HullWhite::VanillaSwaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_FRMSV1F::VanillaSwaption(
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
	if (!isConstantNotional)
		ARM_THROW( ERR_INVALID_ARGUMENT, "MSV1F: variable spread not supported)" );

	if (!isConstantSpread)
		ARM_THROW( ERR_INVALID_ARGUMENT, "MSV1F: variable spread not supported)" );

	if (!isConstantStrike)
		ARM_THROW( ERR_INVALID_ARGUMENT, "MSV1F: variable strike not supported)" );

	if( !GetFixingFunctor()->IsSameModel(*GetDiscountFunctor()) )
        /// We need to compute the floating leg by forward method and no more by double notional
        /// but we have not at the moment all floating leg datas => throw an error
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "Swaption pricing not implemented for differents discount & fixing curves" );

	/// handle the case of dummy states!
	if( states == ARM_PricingStatesPtr(NULL) )
		return static_cast<ARM_VectorPtr>(new std::vector<double>(1,0.0));

	/////// Compute the SwapRate at 0; the Annuity.
	std::vector<double> dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;
	std::vector<double> margin = std::vector<double>(1,0.0);

	ARM_VectorPtr swapRate = SwapRate(
		curveName, 
		evalTime,
		floatStartTime, 
		floatEndTime, 
		fixPayTimes,
		fixPayPeriods,
		dummyFwdStartTimes,
        dummyFwdEndTimes,
        dummyFwdPayTimes,
        dummyFloatPayTimes,
        dummyFloatPayPeriods,
		margin,	/// margin
		true,	/// isDbleNotional to avoid computing the float cash flows piece by piece
		ARM_PricingStatesPtr(NULL) );

	ARM_VectorPtr annuity = Annuity(
		curveName, 
        evalTime,
		fixPayTimes,
        fixPayPeriods,
		ARM_PricingStatesPtr(NULL) );
	
	//////////////////////////////////
	itsCurrentArg = NULL;
	itsCurrentArg = GetVanillaSwaptionArg( curveName, evalTime, swapResetTime,fixNotional[0],
			fixNotional, floatNotional, floatStartTime, floatEndTime,floatResetTimes,floatStartTimes,floatEndTimes,floatIntTerms,
			fixPayTimes,fixPayPeriods, strikesPerState, callPut, states,isConstantNotional );

	////// Given the constant params, we price the swaption using the closed Heston Formula	
	double swapNotional = fixNotional[0];
	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double speed		= ((ARM_CurveModelParam&) (((ARM_ModelParamsFRMSV1F*) (GetModelParams()))->GetVolParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);
	double volvol		= ((ARM_CurveModelParam&) (((ARM_ModelParamsFRMSV1F*) (GetModelParams()))->GetVolParams()->GetModelParam(ARM_ModelParamType::VolOfVol))).GetValueAtPoint(0);// Equivalent Vol Of Vol
	double rho			= 0;
	double initegratedVol = 0.;
	
/*	(const_cast<ARM_FRMSV1F*> (this))->PrecomputeDatas( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;
	double cstShiftValue = ComputeEquivalentCstShift( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods,&initegratedVol) ;
	double cstVolOfVolValue = ComputeEquivalentCstVolOfVol( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;
	double Mu = ComputeMu(cstShiftValue, initegratedVol);


	double V0			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol))).GetValueAtPoint(0);// Equivalent Vol Of Vol	
	ARM_RiccatiMSV* test = new ARM_RiccatiMSV((*GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
										const_cast<ARM_FRMSV1F*> (this),curveName,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,states,cstShiftValue,Mu);

	test->InitFunctionDatas();

	double target = ComputePhi( test,cstShiftValue,swapResetTime);
	ConstantVolatilityFinder func(test,0.0,swapResetTime/K_YEAR_LEN,target,1.0);

	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
	
	///To initialize departure point
	solver.setInitialGuess(0.2);

	double result = solver.Solve();

	double Var = result * result;

	double optionValue = swapNotional * (*annuity)[0] * Export_Shifted_Heston_VanillaOption(
		(*swapRate)[0],strikesPerState(0,0), Var,
		time/K_YEAR_LEN, Var, speed, volvol * result, rho, cstShiftValue, callPut, DEFAULT_NBPOINTS_STAGE1,DEFAULT_NBPOINTS_PER_STAGE, 
		DEFAULT_NBSTAGE,DEFAULT_NBOSCILLATIONS,DEFAULT_PRECISION);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));*/
	
	double volatility_t = ComputeEquivalentVolatilityt(curveName,evalTime,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods);
	return ARM_VectorPtr( new std::vector<double>(1,volatility_t));
}


////////////////////////////////////////////////////
///	Class  : ARM_SFRM
///	Routine: GetVanillaSwaptionArg
///	Returns: ARM_VanillaSwaptionArgSFRM
///	Action : computes all the required data for the pricing of a swaption
///				except the vol
////////////////////////////////////////////////////
ARM_VanillaSwaptionArgSFRM* ARM_FRMSV1F::GetVanillaSwaptionArg( 
	const string& curveName,
	double evalTime,
	double swapResetTime,
	double swapNotional,
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
	bool isConstantNotional	) const
{
	double asOfDate	= GetAsOfDate().GetJulian();
	ARM_Date startDate(asOfDate+floatStartTime);
	ARM_Date endDate(asOfDate+floatEndTime);

	double resetTime = swapResetTime;
	std::vector<double>& pFixPayTimes	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixPayTimes).Clone());
	std::vector<double>& pFixPayPeriods	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixPayPeriods).Clone());
	std::vector<double>& pFloatResetTimes = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatResetTimes).Clone());
	std::vector<double>& pFloatStartTimes = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatStartTimes).Clone());
	std::vector<double>& pFloatEndTimes	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());
	std::vector<double>& pFloatIntTerms	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatIntTerms).Clone());

	std::vector<double>& pfixNotional		 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(fixNotional).Clone());
	std::vector<double>& pfloatNotional	 = static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatNotional).Clone());

    ///Compute the fixed frequency
    int i;
    double periodmax = (*pFixPayPeriods)[0];
    for(i=1; i<pFixPayTimes->size(); ++i)
        periodmax = CC_Max(periodmax,(*pFixPayPeriods)[i]);

    int fixFreq = ROUND(1.0/periodmax);
	/// period of more than every two months, than it should be one month
	if( fixFreq > 8 )
		fixFreq = 12;

	size_t nbFlow = pFloatEndTimes->size();

	/// computes the averageShift
	double averageShift		 = ((ARM_ModelParamsSFRM*) GetModelParams())->AverageShift(*pFloatResetTimes);
    double startTime		 = (*pFloatStartTimes)[0];
	double endTime			 = (*pFloatEndTimes)[nbFlow-1];


	ARM_Currency* pCcy	= const_cast<ARM_FRMSV1F* const>(this)->GetZeroCurve()->GetCurrencyUnit();
	int floatDayCount	= pCcy->GetLiborIndexDayCount();
	std::vector<double> margin = std::vector<double>(1,0.0);
	bool isDbleNotional =false;
	
	std::vector<double>& fwdStartTimes	= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatStartTimes).Clone());
	std::vector<double>& fwdEndTimes		= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());
	std::vector<double>& floatPayTimes	= static_cast<ARM_GP_Vector*>(const_cast<ARM_GP_Vector*>(floatEndTimes).Clone());


	std::vector<double> fwdPayPeriods = std::vector<double>(fwdStartTimes->size());

	for(i=0; i<fwdStartTimes->size(); ++i )
	{
		fwdPayPeriods.Elt(i) = CountYearsWithoutException( floatDayCount, ARM_Date((fwdStartTimes->Elt(i)+asOfDate)), ARM_Date((fwdEndTimes->Elt(i)+asOfDate)) );
	}

    
	ARM_VectorPtr fixAnnuity (NULL);
	ARM_VectorPtr swapFwd (NULL); 
	ARM_VectorPtr fixAnnuityWithNominal (NULL);

	if (isConstantNotional)
	{
		fixAnnuity = Annuity(curveName, evalTime, fixPayTimes, fixPayPeriods, states);
		ARM_VectorPtr clonedfixAnnuity = ARM_CountedPtr<std::vector<double>>( (std::vector<double>&) fixAnnuity->Clone() );
		///Compute a stub to add to fixAnnuity
		ARM_Date tmpDate(asOfDate + pFixPayTimes->Elt(0));
		char* CurrencyName = pCcy->GetCcyName();  
		tmpDate.AddPeriod(-fixFreq,CurrencyName);
		if(fabs(tmpDate-startDate)>FRMVOL_LAG_THRESHOLD)
		{  
			int fixDayCount = pCcy->GetFixedDayCount();
			double delta = CountYears(fixDayCount,tmpDate,startDate);
			ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
			ARM_VectorPtr ZcFwdPay    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,pFixPayTimes->Elt(0),states);
			
			for(int i=0; i<states->size(); ++i)
			{
				(*clonedfixAnnuity)[i] = (*clonedfixAnnuity)[i] + delta*((*ZcFwdPay)[i]-(*ZcFwdStart)[i]);
			}
		} 
		
		swapFwd = SwapRateInPlaceWithComputedAnnuity(	curveName, evalTime,
			startTime, endTime, *pFixPayTimes, *pFixPayPeriods, *fwdStartTimes,
			*fwdEndTimes, fwdPayPeriods, *floatPayTimes, *pFloatIntTerms,
			margin, isDbleNotional, clonedfixAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
			states);
	}
	else
	{
		fixAnnuityWithNominal = AnnuityWithNominal(curveName, evalTime, fixPayTimes, fixPayPeriods, fixNotional,states)		;
		
		ARM_VectorPtr clonedfixAnnuity = ARM_CountedPtr<std::vector<double>>( (std::vector<double>&) fixAnnuityWithNominal->Clone() );
		///Compute a stub to add to fixAnnuity
		ARM_Date tmpDate(asOfDate + pFixPayTimes->Elt(0));
		char* CurrencyName = pCcy->GetCcyName();  
		tmpDate.AddPeriod(-fixFreq,CurrencyName);
		double Notional0 = fixNotional[0];
		if(fabs(tmpDate-startDate)>FRMVOL_LAG_THRESHOLD)
		{  
			int fixDayCount = pCcy->GetFixedDayCount();
			double delta = CountYears(fixDayCount,tmpDate,startDate);
			ARM_VectorPtr ZcFwdStart    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,floatStartTime,states);
			ARM_VectorPtr ZcFwdPay    = GetFixingFunctor()->DiscountFactor(curveName,evalTime,pFixPayTimes->Elt(0),states);
			
			for(int i=0; i<states->size(); ++i)
			{
				(*clonedfixAnnuity)[i] = (*clonedfixAnnuity)[i] + delta*Notional0*((*ZcFwdPay)[i]-(*ZcFwdStart)[i]);
			}
		} 

		swapFwd = SwapRateInPlaceWithComputedAnnuityAndNominal(	curveName, evalTime,
			startTime, endTime, *pFixPayTimes, *pFixPayPeriods, *fwdStartTimes,
			*fwdEndTimes, fwdPayPeriods, *floatPayTimes, *pFloatIntTerms,
			margin, clonedfixAnnuity, //// make sure you clone it if neeeded to keep the value of the annuity
			floatNotional,states);

	}
    /// creates the swaption arg
	ARM_VectorPtr nullstrike(new std::vector<double>(1));
	ARM_VanillaSwaptionArgSFRM* arg = new ARM_VanillaSwaptionArgSFRM(
		resetTime, startTime, endTime, averageShift, pFixPayTimes, pFixPayPeriods, pFloatResetTimes, 
		pFloatStartTimes, pFloatEndTimes, pFloatIntTerms, ARM_VectorPtr(NULL), fixAnnuity, swapFwd,fixFreq, fixAnnuityWithNominal,swapFwd,pfixNotional,pfloatNotional);	

	/// computes the mu
	ARM_VectorPtr mu = ((ARM_ModelParamsSFRM*) GetModelParams())->ComputeVolSwapvolFRA(*arg,*this,isConstantNotional);

	/// AAAAAAAAAAAARRRRRRRRRGGGGGGGGG the ComputeVolSwapvolFRA change the swap fwd!
	arg->SetMu( mu );
	arg->SetSwapFwd( swapFwd );

	return arg;
}


///////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_FRMSV1F::FirstPricingStates( size_t bucketSize ) const
{
	/// ARM_PricingStates(nbStates = bucketSize, nbModelStates = 1F , nbPayoffs = 0)
	ARM_PricingStatesPtr states(new ARM_PricingStates(bucketSize,FactorCount() + 1 ,0,FactorCount()));
	int statesNb = states->size();
	for( size_t i=0;i<statesNb; ++i )
	{
		states->SetModelState(i,THIRD_STATE_VARIABLE,1.0);
	}
	return states;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: ComputeEquivalentLambdat
///	Returns: a double
///	Action : Compute the equivalent Lambdat in the Heston Swap Rate diffusion
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ARM_FRMSV1F::ComputeEquivalentVolatilityt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods) const
{
   /// Global attributes to use the SFRM Code
   bool itsUpdateVolSwapvolFRA = true;
   if(itsUpdateVolSwapvolFRA)
   {
        /// The average shift is set at 0. to keep the same formulae as the SFRM model
	    ((ARM_VanillaSwaptionArgSFRM*)itsCurrentArg)->SetAverageShift(0.0);
        /// Computes  Mu
        ((ARM_VanillaSwaptionArgSFRM*)itsCurrentArg)->SetMu( 
			(((ARM_ModelParamsFRMSV1F*) GetModelParams())->GetSFRMParams())->ComputeVolSwapvolFRA(*(ARM_VanillaSwaptionArgSFRM*)itsCurrentArg,*this,true) ); // true for cst notional Flag
    }
	

	///// the function is valid only at 0
	ARM_VanillaSwaptionArgSFRM* arg = (ARM_VanillaSwaptionArgSFRM*) &*itsCurrentArg;
	double swapVol = (((ARM_ModelParamsFRMSV1F*) GetModelParams())->GetSFRMParams())->LocalVolatity(0.0,evalTime,*arg,arg->GetMu() );
	swapVol *= sqrt((floatStartTime-evalTime)/K_YEAR_LEN);

	return swapVol;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_FRMSV1F
///	Routine: ComputeEquivalentShiftt
///	Returns: a double
///	Action : Compute the equivalent Shiftt in the Heston Swap Rate diffusion
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ARM_FRMSV1F::ComputeEquivalentShiftt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const std::vector<double>& fixPayTimes,
		const std::vector<double>& fixPayPeriods,
		double lambdaValue) const
{
	/// For Test Purpose
	return (GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->Interpolate(evalTime));	
}


////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_FRMSV1F::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"LocalDiscounts : unimplemented function for ARM_FRMSV1F Model!");
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_FRMSV1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_MarkovSV Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_FRMSV1F::EulerLocalDrifts( const std::vector<double>& timeSteps,
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
///	Class   : ARM_FRMSV1F
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_FRMSV1F::VolatilitiesAndCorrelations( const std::vector<double>& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
}


////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: IntegratedBondSquaredVol
///	Returns : double
///	Action  : Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_FRMSV1F::IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_FRMSV1F Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: IntegratedBondCovariance
///	Returns : double
///	Action  : Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_FRMSV1F::IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_FRMSV1F Model!");

}

////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: VolatilityScalarProduct
///	Returns : double (qui l'eut cru)
///	Action  :  Int_startTime^endTime Gamma(s,bondMaturity) * dW_s 
////////////////////////////////////////////////////
double ARM_FRMSV1F::VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_FRMSV1F Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_FRMSV1F
///	Routines: ModelStateLocalCorrels
///	Returns : void
///	Action  :  
////////////////////////////////////////////////////
void ARM_FRMSV1F::ModelStateLocalCorrels( const std::vector<double>& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& multiassets )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_FRMSV1F Model!");


}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

