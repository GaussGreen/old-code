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
#include "gpmodels/MSV1F.h"
#include "gpmodels/ModelParamsMSV1F.h"
#include "gpmodels/riccatimsv.h"

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
const double INTERPOL_EPS = 1.e-5;
///For MC simulation
const double FWD_MC_MAX = 0.5;

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: Default Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////

ARM_MarkovSV1F::ARM_MarkovSV1F(const ARM_ZeroCurvePtr& zc,const ARM_ModelParams* params, double forwardTerm, bool isSwapRate )
:ARM_MarkovSV(zc,params), itsForwardTerm (forwardTerm * K_YEAR_LEN), itsFixDateStrip(NULL), itsSwapRateCurve(NULL),itsSwapRateDerivativeCurve(NULL)
{	
	itsIsSwapRate = isSwapRate;
	if (itsIsSwapRate )
	{
		double asOfDate	= GetAsOfDate().GetJulian();
		double endTime = asOfDate + itsForwardTerm;
		itsFixDateStrip = GetFloatDateStrip(asOfDate,endTime,((const ARM_ModelParamsMSV1F*) GetModelParams())->GetIRIndex());
	
		// Create the SwapRateCurve
		ARM_PricingStatesPtr tempStates(NULL);
		double annuity =0,annuityDerivative=0, swapRate=0,evalTime;
		ARM_DateStrip* pFixDateStrip = GetDateStrip();
		ARM_GP_Vector* fixPayTimes   = pFixDateStrip->GetPaymentDates();
		ARM_GP_Vector* fixPayPeriods = pFixDateStrip->GetInterestTerms();
		ARM_GP_Vector* resetTimes    = pFixDateStrip->GetResetDates();
		ARM_GP_Vector* startTimes	 = pFixDateStrip->GetFlowStartDates();

		size_t nbFixFlows = fixPayTimes->size();
		ARM_GP_Vector rates(nbFixFlows, 0.0);
		ARM_GP_Vector ratesDerivatives(nbFixFlows, 0.0);
		ARM_GP_Vector dates(nbFixFlows, 0.0);

		for (int iFix=0;iFix<nbFixFlows;++iFix)
		{
			annuity = 0;
			annuityDerivative = 0;
			//Compute The swapRate
			evalTime = (*startTimes)[iFix]-asOfDate;
			double floatStartTime = evalTime;
			
			if (iFix >= nbFixFlows)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"The forward Swap End Time is inferior to the Reset Date!");
			
			for(int l=iFix;l<nbFixFlows;++l)
			{
				double payTime = (*fixPayTimes)[l]-asOfDate;
				ARM_VectorPtr ZcFixFlow = DiscountFactor("Toto",evalTime,payTime,tempStates);
				double BtTi = (*ZcFixFlow)[0];
				double BetatTi = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,(*fixPayTimes)[iFix]);
				double flow_i = (*fixPayPeriods)[l]*BtTi;
				annuity+= flow_i;
				annuityDerivative -= BetatTi * flow_i;

			}
			if (annuity < K_NEW_DOUBLE_TOL)
				throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
				"ComputeFwdTerm  : Annuity is Null!");
			
			////// Compute the float leg value and Derivative
			double floatLegValue = 0;
			ARM_VectorPtr ZcVarStrat = DiscountFactor("Toto",evalTime,floatStartTime,tempStates);
			ARM_VectorPtr ZcVarEnd	 = DiscountFactor("Toto",evalTime,itsForwardTerm,tempStates);
			double BtT0 = (*ZcVarStrat)[0];
			double BtTN = (*ZcVarEnd)[0];
			double BetatT0 = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,floatStartTime);
			double BetatTN = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,itsForwardTerm);
			double FloatLegDerivative = BtTN * BetatTN - BtT0 * BetatT0;


			floatLegValue = BtT0 - BtTN;
			swapRate = floatLegValue/annuity;
			double swapRateDerivative = (annuity * FloatLegDerivative  -  floatLegValue * annuityDerivative)/(annuity * annuity);
			rates[iFix] = swapRate;
			dates[iFix] = evalTime;
			ratesDerivatives[iFix] = swapRate;
		}
		ARM_StepUpRightOpenCstExtrapolDble* interpolator = new ARM::ARM_StepUpRightOpenCstExtrapolDble();
		itsSwapRateCurve = new ARM_Curve(dates,rates, interpolator);	
		itsSwapRateDerivativeCurve = new ARM_Curve(dates,rates, interpolator);	

	}
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: Copy Constructor
///	Returns: 
///	Action : Constructor
////////////////////////////////////////////////////
ARM_MarkovSV1F::ARM_MarkovSV1F(const ARM_MarkovSV1F& rhs)
: ARM_MarkovSV(rhs) 
{
	itsForwardTerm  = rhs.itsForwardTerm;
	itsFixDateStrip = rhs.itsFixDateStrip? (ARM_DateStrip*) rhs.itsFixDateStrip->Clone() : NULL;
	itsSwapRateCurve = rhs.itsSwapRateCurve? (ARM_Curve*) rhs.itsSwapRateCurve->Clone() : NULL;
	itsSwapRateDerivativeCurve = rhs.itsSwapRateDerivativeCurve? (ARM_Curve*) rhs.itsSwapRateDerivativeCurve->Clone() : NULL;
	itsIsSwapRate	= rhs.itsIsSwapRate; 
}


/////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: Destructor
///	Returns: 
///	Action : Destructor
////////////////////////////////////////////////////
ARM_MarkovSV1F::~ARM_MarkovSV1F( )
{
	delete itsFixDateStrip;
	delete itsSwapRateCurve;
	delete itsSwapRateDerivativeCurve;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: operator =
///	Returns: this
///	Action : Affectation of a rhs object
////////////////////////////////////////////////////
ARM_MarkovSV1F& ARM_MarkovSV1F::operator=(const ARM_MarkovSV1F& rhs)
{
	if(this != &rhs)
	{
		ARM_MarkovSV::operator=(rhs);	
		itsForwardTerm		= rhs.itsForwardTerm;
		itsFixDateStrip		= rhs.itsFixDateStrip? (ARM_DateStrip*) rhs.itsFixDateStrip->Clone() : NULL;
		itsSwapRateCurve	= rhs.itsSwapRateCurve? (ARM_Curve*) rhs.itsSwapRateCurve->Clone() : NULL;
		itsSwapRateDerivativeCurve = rhs.itsSwapRateDerivativeCurve? (ARM_Curve*) rhs.itsSwapRateDerivativeCurve->Clone() : NULL;
		itsIsSwapRate		= rhs.itsIsSwapRate;
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_MarkovSV1F::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "1F Markov Stochastic Volatility Model\n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_MarkovSV1F::Clone() const
{
	return new ARM_MarkovSV1F(*this);
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_MarkovSV1F::PreProcessing(ARM_ModelFitter& modelFitter)
{ 
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: PostProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_MarkovSV1F::PostProcessing(const ARM_ModelFitter& modelFitter)
{
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_MarkovSV1F::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ValidateModelParams
///	Returns: true/false
///	Action : Check the consistency of the model
///          parameters
////////////////////////////////////////////////////
bool ARM_MarkovSV1F::ValidateModelParams(const ARM_ModelParams& params) const
{
    if(!params.DoesModelParamExist(ARM_ModelParamType::MeanReversion) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::InitialVol) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolOfVol) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::VolMeanReversion) ||
	   !params.DoesModelParamExist(ARM_ModelParamType::Shift)||
	   !params.DoesModelParamExist(ARM_ModelParamType::LongTermVol))
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "At least 1 Model Param is not of a good type!");
    }
	return true;
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: LocalDrifts
///	Returns: A vector saving the local drift of
///          the state variable
///	Action : Compute local relative drifts of the
///          state variable between each time step
////////////////////////////////////////////////////
void ARM_MarkovSV1F::IntegratedLocalDrifts(
	const ARM_GP_Vector& timeSteps,
	ARM_GP_MatrixPtr& relativeDrifts,
	ARM_GP_MatrixPtr& absoluteDrifts ) const
{
}



////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ModelStateLocalVariances
///	Returns: A vector of 1D matrix
///	Action : Compute local and global variances
///          (from 0) of the state variable between
///          each time step
////////////////////////////////////////////////////
void ARM_MarkovSV1F::ModelStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
    ARM_MatrixVector& localVariances ) const
{
}


////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: NumMethodStateLocalVariances
///	Returns: A vector of ND matrix
///	Action : Compute local variances of the state variable 
/// between each time step
////////////////////////////////////////////////////
void ARM_MarkovSV1F::NumMethodStateLocalVariances(
    const ARM_GP_Vector& timeSteps,
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
void ARM_MarkovSV1F::NumMethodStateGlobalVariances(
    const ARM_GP_Vector& timeSteps,
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
double ARM_MarkovSV1F::VarianceToTime(double var,double minTime,double maxTime) const
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
ARM_VectorPtr ARM_MarkovSV1F::DiscountFactor( 
		const string& curveName, 
		double evalTime, 
		double maturityTime, 
		const ARM_PricingStatesPtr& states) const
{
	// Waiting for the access to the yield curve with curveName
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
    double zcT=ZcCurve->DiscountPrice(maturityTime/K_YEAR_LEN);
    double zct;
	double Beta;
	double Beta_2;

	if(		evalTime <= K_NEW_DOUBLE_TOL
		 || states   == ARM_PricingStatesPtr(NULL) )
    {
		size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
		return ARM_VectorPtr( new ARM_GP_Vector(payoffSize,zcT) );
    }

    // Volatility computation (ARM_ModelParamsMSV1F class is pure virtual)
    int i,nbStates=states->size();
    if(    GetNumeraire() == ARM_NumerairePtr(NULL)
		|| GetNumeraire()->GetType() == ARM_Numeraire::Cash)
    {
        if(evalTime < maturityTime)
        {
            // Could be optimized using phis(t)...?
            Beta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,maturityTime);
            zct=ZcCurve->DiscountPrice(evalTime/K_YEAR_LEN);
            Beta_2 = Beta*Beta;
        }
        else
            return ARM_VectorPtr(new ARM_GP_Vector(nbStates,1.0));
    }


	size_t modelNb = GetModelNb();
    ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
	double temp = 0.;
    for(i=0;i<nbStates;++i)
    {
		(*values)[i] = zcT/zct*exp( -Beta * states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE)
                        - 0.5 * Beta_2 * states->GetModelState(i,modelNb+SECOND_STATE_VARIABLE) );
	}
    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: MCModelStatesFromToNextTime
///	Returns: void
///	Action : 
////////////////////////////////////////////////////

void ARM_MarkovSV1F::MCModelStatesFromToNextTime(ARM_PricingStatesPtr& states,int timeIndex) const
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


	////// We have constant Mean Reversion
	double Lambda_t = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
	double VolMeanReversion_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).GetValueAtPoint(0);
	double LongTermVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::LongTermVol).GetValueAtPoint(0);


	/// We have to Interpolate
	double Volatility_t,VolOfVol_t;

	/////// To Avoid Discontinuities in the derivative calculation///////////////////////////////////////
	for(size_t j=0;j<((ARM_ModelParamsMSV*) GetModelParams())->GetSchedule().size();++j)
	{
		double DiscontTime = ((ARM_ModelParamsMSV*) GetModelParams())->GetSchedule().Elt(j);
		if(((DiscontTime + K_NEW_DOUBLE_TOL)<time)&&(time<(DiscontTime + K_NEW_DOUBLE_TOL)))
				time+=INTERPOL_EPS;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////
	
	Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(time);
	VolOfVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()->Interpolate(time);
	
	const ARM_VectorPtr EtaVector = ComputeEtaStates(time , states);

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

		double y1 = Y_State;
		/// Second Variable Yt
		Y_State = Y_State + (Eta*Eta - 2*Lambda_t*Y_State)*dt;
		states->SetModelState(i,modelNb+SECOND_STATE_VARIABLE,Y_State);

		X_State = X_State + (0.5*(Y_State + y1) - Lambda_t * X_State)*dt + sqrt(dt)*gaussian_1*Eta ;
		states->SetModelState(i,modelNb+FIRST_STATE_VARIABLE,X_State);


		/// Third Variable Zt
		////Standard Euler Discretisation approach
		
		double gaussian_2 = states->GetNumMethodState(i,modelNb+SECOND_STATE_VARIABLE);
		/// Third Variable Zt
		/////Standard Euler Discretisation approach  + MilsteinD - 
	//	Z_State = Z_State + VolMeanReversion_dt * (1- Z_State) + VolOfVol_dt * sqrt(Z_State)*gaussian_2 + 0.25 * VolOfVol_2_t * dt*(gaussian_2*gaussian_2 - 1);
	//	if (Z_State <0)
	//		Z_State = - Z_State;

/*		///// Approach 1 to avoid Negative Values of Zt
		Z_State = Z_State + VolMeanReversion_dt * (1- Z_State) + VolOfVol_dt * sqrt(Z_State)*gaussian_2;
		if (Z_State<0)
			Z_State = - Z_State;
		states->SetModelState(i,modelNb+THIRD_STATE_VARIABLE,Z_State);*/
	
/*		double temp =  VolOfVol_dt * sqrt(Z_State)*gaussian_2;
		Z_State = Z_State + VolMeanReversion_dt * (1- Z_State) + temp;
		if (Z_State<0)
			Z_State-= temp;
		states->SetModelState(i,modelNb+THIRD_STATE_VARIABLE,Z_State);*/


		//// New Approach : Andreasen
		//// Intermediate Variables
	/*	double Z_State_bis = Z_State*var_1 + var_2;	
		double var_i = 1 + 1/(Z_State_bis*Z_State_bis)*(var_5+(Z_State-1)*var_4);
		double Nu_2 = log(fabs(var_i-1)+1);
		double Nu = sqrt(Nu_2);

		Z_State = Z_State_bis*exp(-0.5*Nu_2+ Nu*gaussian_2);
		states->SetModelState(i,modelNb+THIRD_STATE_VARIABLE,Z_State);*/

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
///	Class  : ARM_MarkovSV1F
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : Pricing of a variable notional swaption
///          via numerical integration. 
///			 If the swaption is standard, this method
///          calls ARM_HullWhite::VanillaSwaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV1F::VanillaSwaption(
				const string& curveName,
				double evalTime,
				double swapResetTime,
				const ARM_GP_Vector& fixNotional,
				const ARM_GP_Vector& floatNotional,
				double floatStartTime,
				double floatEndTime,
				const ARM_GP_Vector& floatResetTimes,
				const ARM_GP_Vector& floatStartTimes,
				const ARM_GP_Vector& floatEndTimes,
				const ARM_GP_Vector& floatIntTerms,
				const ARM_GP_Vector& fixPayTimes,
				const ARM_GP_Vector& fixPayPeriods,
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
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/////// Compute the SwapRate at 0; the Annuity.
	ARM_GP_Vector dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;
	ARM_GP_Vector margin = ARM_GP_Vector(1,0.0);

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

		
	////// Given the constant params, we price the swaption using the closed Haston Formula	
	double swapNotional = fixNotional[0];
	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double speed		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);
	double rho			= 0;
	double initegratedVol = 0.;
	
	(const_cast<ARM_MarkovSV1F*> (this))->PrecomputeDatas( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;
	double cstVolOfVolValue = ComputeEquivalentCstVolOfVol( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;
	double cstShiftValue = ComputeEquivalentCstShift( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods,cstVolOfVolValue,&initegratedVol) ;
	double Mu = ComputeMu(cstShiftValue, initegratedVol);


	double V0			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol))).GetValueAtPoint(0);	
	ARM_RiccatiMSV test((*GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()),
										(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
										const_cast<ARM_MarkovSV1F*> (this),curveName,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,states,cstShiftValue,Mu);

	test.InitFunctionDatas();
	double at = 0.5 *  cstVolOfVolValue * cstVolOfVolValue ;
	test.SetTempAt(at);

	double target = ComputePhi(test,cstShiftValue,swapResetTime);
	ConstantVolatilityFinder func(&test,0.0,swapResetTime/K_YEAR_LEN,target,1.0);

	UnaryFuncWithNumDerivative<double> funcWithDev(func);

	T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
	
	///To initialize departure point
	solver.setInitialGuess(0.2);

	double result = solver.Solve();
	double Var = result * result;

	double optionValue = swapNotional * (*annuity)[0] * Export_Shifted_Heston_VanillaOption(
		(*swapRate)[0],strikesPerState(0,0), Var,
		time/K_YEAR_LEN, Var, speed, cstVolOfVolValue * result, rho, cstShiftValue, callPut, DEFAULT_NBPOINTS_STAGE1,DEFAULT_NBPOINTS_PER_STAGE, 
		DEFAULT_NBSTAGE,DEFAULT_NBOSCILLATIONS,DEFAULT_PRECISION);

	return ARM_VectorPtr( new ARM_GP_Vector(1,optionValue));
}

////////////////////////////////////////////////////
///	Class  : ARM_PricingFunctionIR
///	Routine: VanillaSmiledSwaption
///	Returns: a vector of sum option values
///	Action : Closed form formula for Smiled Swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_MarkovSV1F::VanillaSmiledSwaption(
		const string& curveName,
        double evalTime,
		double swapResetTime,
		const ARM_GP_Vector& fixNotional,
		const ARM_GP_Vector& floatNotional,
		double floatStartTime,
        double floatEndTime,
		const ARM_GP_Vector& floatResetTimes,
		const ARM_GP_Vector& floatStartTimes,
		const ARM_GP_Vector& floatEndTimes,
		const ARM_GP_Vector& floatIntTerms,
		const ARM_GP_Vector& fixPayTimes,
        const ARM_GP_Vector& fixPayPeriods,
        const ARM_GP_Vector& strikesPerState,
        int callPut,
        const ARM_PricingStatesPtr& states,
		const ARM_GP_Vector& data,
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
		return static_cast<ARM_VectorPtr>(new ARM_GP_Vector(1,0.0));

	/////// Compute the SwapRate at 0; the Annuity.
	ARM_GP_Vector dummyFwdStartTimes, dummyFwdEndTimes, dummyFwdPayPeriods, dummyFwdPayTimes, dummyFloatPayTimes, dummyFloatPayPeriods;
	ARM_GP_Vector margin = ARM_GP_Vector(1,0.0);

		
	////// Given the constant params, we price the swaption using the closed Haston Formula	
	double swapNotional = fixNotional[0];
	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double speed		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetValueAtPoint(0);
	double rho			= 0;
	double initegratedVol = 0.;


	/////////// Targets and Weights of the Smiled Swaption //////////////////////
	if (data.size()<6)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "At Least 1 target Param is missing : Volatility, Shift , VolOfVol !" );
	double volTarget = data[0];
	double varTarget = volTarget * volTarget;
	double shiftTarget = data[1];
	double volOfVolTarget = data[2];

	double weight1 = data[3];
	double weight2 = data[4];
	double weight3 = data[5];
	//////////////////////////////////////////////////////////////////////////////
	
	(const_cast<ARM_MarkovSV1F*> (this))->PrecomputeDatas( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;
	double cstVolOfVolValue = ComputeEquivalentCstVolOfVol( curveName,swapResetTime, floatStartTime,floatEndTime,
		fixPayTimes,fixPayPeriods) ;

	double cstShiftValue = 0.;
	////// Temp Variables
	double Mu,V0,target,effectiveVol = 0.;
	if(!((weight1<K_NEW_DOUBLE_TOL)&&(weight2<K_NEW_DOUBLE_TOL)))
	{
		cstShiftValue = ComputeEquivalentCstShift( curveName,swapResetTime, floatStartTime,floatEndTime,
			fixPayTimes,fixPayPeriods,cstVolOfVolValue,&initegratedVol) ;
		
		if (!(weight1<K_NEW_DOUBLE_TOL))
		{
			
			Mu = ComputeMu(cstShiftValue, initegratedVol);
			V0			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol))).GetValueAtPoint(0);// Equivalent Vol Of Vol	
			ARM_RiccatiMSV test((*GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).ToCurveModelParam().GetCurve()),
				(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
				(*GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()),
				(*GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).ToCurveModelParam().GetCurve()),
				const_cast<ARM_MarkovSV1F*> (this),curveName,floatStartTime,floatEndTime,fixPayTimes,fixPayPeriods,states,cstShiftValue,Mu);
			
			test.InitFunctionDatas();
			double at = 0.5 *  cstVolOfVolValue * cstVolOfVolValue ;
			test.SetTempAt(at);
			
			target = ComputePhi(test,cstShiftValue,swapResetTime);
			ConstantVolatilityFinder func(&test,0.0,swapResetTime/K_YEAR_LEN,target,1.0);
			
			UnaryFuncWithNumDerivative<double> funcWithDev(func);
			T_NewtonRaphsonSolver< UnaryFuncWithNumDerivative<double> > solver(funcWithDev,0,DEFAULT_PRECISION,DEFAULT_PRECISION);
			
			///To initialize departure point
			solver.setInitialGuess(0.2);
			effectiveVol = solver.Solve();
		}
	}
		
	double result = weight1 * (volTarget - effectiveVol)*(volTarget - effectiveVol)+
			weight2 * (shiftTarget - cstShiftValue)*(shiftTarget - cstShiftValue)+
			weight3 * (volOfVolTarget - cstVolOfVolValue)*(volOfVolTarget - cstVolOfVolValue);
	return ARM_VectorPtr( new ARM_GP_Vector(1,result));

}

////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV1F::VanillaCaplet(
	const string& curveName, 
	double evalTime,
	double payTime, 
	double period,
    double payNotional,
	double fwdResetTime, 
	double fwdStartTime,
    double fwdEndTime,
	double fwdPeriod,
	const ARM_GP_Vector& strikesPerState,
    int capFloor,
	const ARM_PricingStatesPtr& states) const
{

	double delta = (fwdEndTime-fwdStartTime)/K_YEAR_LEN;

	ARM_GP_Vector fixNotional(1,payNotional);
	ARM_GP_Vector floatNotional(1,payNotional);
	ARM_GP_Vector floatResetTimes(1,fwdResetTime);
	ARM_GP_Vector floatStartTimes(1,fwdStartTime);
	ARM_GP_Vector floatEndTimes(1,fwdEndTime);
	ARM_GP_Vector floatIntTerms(1,delta);
	ARM_GP_Vector fixPayTimes(1,fwdEndTime);
	ARM_GP_Vector fixPayPeriods(1,delta);
	ARM_GP_Matrix strikesPerStateMatrix; 
	
	return VanillaSwaption(
				curveName,
				evalTime,
				fwdResetTime,
				fixNotional,
				floatNotional,
				fwdStartTime,
				fwdEndTime,
				floatResetTimes,
				floatStartTimes,
				floatEndTimes,
				floatIntTerms,
				fixPayTimes,
				fixPayPeriods,
				strikesPerStateMatrix,
				capFloor,
				states,
				true,
				true,
				true) ;
}


///////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routine : FirstPricingStates,
///	Returns :
///	Action  : create the first pricing state
////////////////////////////////////////////////////
ARM_PricingStatesPtr ARM_MarkovSV1F::FirstPricingStates( size_t bucketSize ) const
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

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeFwdTerm 
///	Returns: F(0,T) or S(0,T)
///	Action : Compute the Forward Term
////////////////////////////////////////////////////
double ARM_MarkovSV1F::ComputeFwdTerm(
					double evalTime) const
{
    /// Test JMP
//    return 1.0;

	if (!itsIsSwapRate)
	{
		/*double evalTerm = evalTime/K_YEAR_LEN;
		double delta =  itsForwardTerm/K_YEAR_LEN;
		double zcStart	= GetZeroCurve()->DiscountPrice(evalTerm);
		double zcEnd	= GetZeroCurve()->DiscountPrice(evalTerm + delta);
		double F0 = (zcStart/zcEnd-1)/delta;*/

		//// For continuus forward rate 
		double endTime = (itsForwardTerm + evalTime)/K_YEAR_LEN;		
		double epsilon = 1e-4; //DERIVATIVE_EPS;
		double Zc =  GetZeroCurve()->DiscountPrice(endTime);
		double Zc_eps =  GetZeroCurve()->DiscountPrice(endTime + epsilon);
		double F0 = (log(Zc) - log(Zc_eps))/epsilon;
		return F0;
	}
	else // SwapRate
	{
		return GetSwapRateCurve()->Interpolate(evalTime);
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeRate_Zero 
///	Returns: F(0,T) or S(0,T)
///	Action : Compute the Forward Term at time t
////////////////////////////////////////////////////
ARM_GP_Vector* ARM_MarkovSV1F::ComputeFwdTerm_t(
					double evalTime,
					double F0,
					const ARM_PricingStatesPtr& states) const
{	
	size_t payoffSize = states != ARM_PricingStatesPtr(NULL) ? states->size(): 1;
	ARM_GP_Vector* values = new ARM_GP_Vector(payoffSize);
	double Ft,Xt,Yt;
	if (!itsIsSwapRate)
	{
		double evalTerm = evalTime/K_YEAR_LEN;
		double delta =  itsForwardTerm/K_YEAR_LEN;
		double endFwdTime = evalTime + itsForwardTerm;
		
		double zcStart	= GetZeroCurve()->DiscountPrice(evalTerm);
		double zcEnd	= GetZeroCurve()->DiscountPrice(evalTerm + delta);

		double Beta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,endFwdTime);
		double DBeta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->Deriv_BetatT(evalTime,endFwdTime);
		double Beta_2 = Beta * Beta;

		for(size_t i=0;i<payoffSize;++i)
		{ 
			Xt = states->GetModelState(i,GetModelNb()+FIRST_STATE_VARIABLE);
			Yt = states->GetModelState(i,GetModelNb()+SECOND_STATE_VARIABLE);
		//	Ft = (zcStart/zcEnd * exp(Beta * Xt + 0.5 * Beta_2 * Yt) -1)/delta;
			//// For instantaneous fwd rate
			Ft = F0 + DBeta * Xt + DBeta * Beta * Yt;
			(*values)[i] = Ft ;
		}
	}
	else // SwapRate 
	{
		////////////////////////////////////////////////
		//// WARNING : FUNCTION TO CHECK
		////////////////////////////////////////////////

		double asOfDate = GetAsOfDate().GetJulian();
		double endTime = itsForwardTerm;
		double startDate = evalTime;
		int i=0;
		while ((i<itsSwapRateDerivativeCurve->size()) && (itsSwapRateDerivativeCurve->GetAbscisse(i)>=startDate))
		{
			startDate = itsSwapRateDerivativeCurve->GetAbscisse(i);
			++i;
		}
		ARM_DateStrip* dateStrip = GetFloatDateStrip(startDate + asOfDate,endTime + asOfDate ,((const ARM_ModelParamsMSV1F*) GetModelParams())->GetIRIndex());
	
		// Create the SwapRateCurve
		ARM_PricingStatesPtr tempStates(NULL);
		double annuity =0,annuityDerivative=0, swapRate=0;
		ARM_DateStrip* pFixDateStrip = GetDateStrip();
		ARM_GP_Vector* fixPayTimes   = pFixDateStrip->GetPaymentDates();
		ARM_GP_Vector* fixPayPeriods = pFixDateStrip->GetInterestTerms();
		ARM_GP_Vector* resetTimes    = pFixDateStrip->GetResetDates();
		ARM_GP_Vector* startTimes	 = pFixDateStrip->GetFlowStartDates();
		ARM_GP_Vector margin(startTimes->size(),0.0);

		ARM_GP_VectorPtr swapRateVector = SwapRate("Toto", evalTime,startDate,endTime,(*fixPayTimes),(*fixPayPeriods),(*startTimes),(*fixPayTimes),(*fixPayPeriods),(*fixPayTimes),(*fixPayPeriods),margin,true,states);
		for(i=0;i<payoffSize;++i)
		{ 
			(*values)[i] = (*swapRateVector)[i] ;
		}
	}
	return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeEta 
///	Returns: Eta(t,T)
///	Action : Compute the volatility function of the process Xt
///			for quick calculation of the shift_t
////////////////////////////////////////////////////
double ARM_MarkovSV1F::ComputeEta( 
		double evalTime, 
		double x) const
{
	double evalTerm = evalTime/K_YEAR_LEN;
	double zcStart	= GetZeroCurve()->DiscountPrice(evalTerm);
	double delta =  itsForwardTerm/K_YEAR_LEN;
    double zcEnd	= GetZeroCurve()->DiscountPrice(evalTerm + delta);
	double F0 = ComputeFwdTerm(evalTime);
	double endFwdTime = itsForwardTerm + evalTime;
	double interpolTime = evalTime;
	if (evalTime <K_NEW_DOUBLE_TOL )
		interpolTime+=INTERPOL_EPS;

	double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(interpolTime);	
	double DBeta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->Deriv_BetatT(evalTime,(endFwdTime));
	double Eta;	

    ///// test JMP
//    return Volatility_t;
    

	if (x > K_NEW_DOUBLE_TOL || x < - K_NEW_DOUBLE_TOL)
	{
		////// We have constant Mean Reversions
		double Lambda_t = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
		
		////// We interpolate
		double Shift_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->Interpolate(interpolTime);
		double var_3 = Volatility_t /DBeta;
		double Beta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,endFwdTime);
	//	double Ft = (zcStart/zcEnd * exp(Beta *x) -1)/delta; 
		///// For instantaneous forward rate 
		double Ft = F0 + DBeta * x ;
		Eta = var_3 * (Shift_t * Ft + (1-Shift_t)*F0);
	}
	else
		Eta = Volatility_t /DBeta * F0;
	return Eta;
}

////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeEtaStates // The forwardTime has to be changed by the chosen forward // parameter of the model
///	Returns: a vector of Eta(t,T)
///	Action : Compute the volatility function of the process Xt
///			Eta is function of Xt , Yt and Zt
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV1F::ComputeEtaStates( 
		double evalTime, 
		const ARM_PricingStatesPtr& states) const
{
	double endFwdTime = itsForwardTerm + evalTime;
	if(	states   == ARM_PricingStatesPtr(NULL) )
    {
		double DBeta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->Deriv_BetatT(evalTime,(endFwdTime));
		double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
		double var_3 = Volatility_t /DBeta;		
		double F0 = ComputeFwdTerm(evalTime);
		double result = var_3 * F0;
		return ARM_VectorPtr( new ARM_GP_Vector(1,result) );
    }

	int i,nbStates=states->size();
	ARM_VectorPtr values(new ARM_GP_Vector(nbStates));
	double sqrt_Zt, Ft, F0;
	size_t modelNb = GetModelNb();

	F0 = ComputeFwdTerm(evalTime);
	////// We have constant Mean Reversions
	double Lambda_t = GetModelParams()->GetModelParam( ARM_ModelParamType::MeanReversion).GetValueAtPoint(0);
	////// We interpolate
	double Shift_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
	double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime);
	
	double DBeta;
	if (!itsIsSwapRate)
	{
		DBeta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->Deriv_BetatT(evalTime,endFwdTime);
	}
	else
		DBeta = GetSwapRateDerivativeCurve()->Interpolate(evalTime);

	double Beta = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,endFwdTime);	
	double var_3 = Volatility_t /DBeta;
	ARM_GP_Vector* FtVector = ComputeFwdTerm_t(evalTime,F0,states);

    for(i=0;i<nbStates;++i)
	{ 
		sqrt_Zt =states->GetModelState(i,modelNb+THIRD_STATE_VARIABLE);
		sqrt_Zt = sqrt(sqrt_Zt);
		Ft = (*FtVector)[i];
		/// We delimit the forward at FWD_MC_MAX to avoid explosions
		Ft = (Ft<FWD_MC_MAX)?Ft:FWD_MC_MAX;
		Ft = (Ft<-FWD_MC_MAX)?-FWD_MC_MAX:Ft;
        (*values)[i]= var_3 * sqrt_Zt* (Shift_t * Ft + (1-Shift_t) * F0);
	}
	delete FtVector;
    return values;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeDerivativeSt
///	Returns: a double
///	Action : Compute the derivative of the swapRate
///          at X=x0; Y=y0
///			This derivative is used for the calculation of time-dependent
///         values of Lambda and Shift in the determination of the swapRate diffusion
///         WARNING : the states contains only the state of X and Y at which we want to calculate the derivative
///						This state could be 0 or better the approximated levels under the annuity measure
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ARM_MarkovSV1F::ComputeDerivativeSt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		double* annuityValue,
		double* swapRate,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		const ARM_PricingStatesPtr& states,
		double x) const
{
	/////// We calculate the Annuity Value inside the function to avoid repetition in ZC Calculation
	/////// We save the ZC values into a vector
	/////// We have ONLY ONE STATE x0, y0

    double Annuity = 0,AnnuityDerivative = 0;

    //Compute The annuity value and serivative
    size_t iFix,nbFixFlows=fixPayTimes.size();
    for(iFix=0;iFix<nbFixFlows;++iFix)
    {
        ARM_VectorPtr ZcFixFlow = DiscountFactor(curveName,evalTime,fixPayTimes[iFix],states);
		double BtTi = (*ZcFixFlow)[0];
		double BetatTi = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,fixPayTimes[iFix]);
		if (x > K_NEW_DOUBLE_TOL || x < -K_NEW_DOUBLE_TOL)
		{
			BtTi = BtTi * exp( - BetatTi * x);
		}
		double flow_i = fixPayPeriods[iFix]*BtTi;
		Annuity += flow_i;
		AnnuityDerivative -= BetatTi * flow_i;
    }
	if (Annuity < K_NEW_DOUBLE_TOL)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Swap Rate Derivative  : Annuity is Null!");
    
	////// Compute the float leg value and Derivative
	double FloatLegValue = 0,FloatLegDerivative = 0.0;
	ARM_VectorPtr ZcVarStrat = DiscountFactor(curveName,evalTime,floatStartTime,states);
	ARM_VectorPtr ZcVarEnd	 = DiscountFactor(curveName,evalTime,floatEndTime,states);
	double BtT0 = (*ZcVarStrat)[0];
	double BtTN = (*ZcVarEnd)[0];

	double BetatT0 = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,floatStartTime);
	double BetatTN = ((const ARM_ModelParamsMSV1F*) GetModelParams())->BetatT(evalTime,floatEndTime);

	
	if (x > K_NEW_DOUBLE_TOL || x < -K_NEW_DOUBLE_TOL)
	{
		BtT0 = BtT0 * exp( - BetatT0 * x);
		BtTN = BtTN * exp( - BetatTN * x);
	}
	FloatLegValue = BtT0 - BtTN;
	FloatLegDerivative = BtTN * BetatTN - BtT0 * BetatT0;

	*annuityValue = Annuity;
	*swapRate = FloatLegValue/Annuity;

	double result = (Annuity * FloatLegDerivative  -  FloatLegValue * AnnuityDerivative)/(Annuity * Annuity);
	return (result);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeEquivalentLambdat
///	Returns: a double
///	Action : Compute the equivalent Lambdat in the Heston Swap Rate diffusion
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double ARM_MarkovSV1F::ComputeEquivalentVolatilityt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods) const
{
	double annuity, swapRate;
	ARM_PricingStatesPtr tempStates(NULL);
	double swapDerivative = ComputeDerivativeSt(curveName,evalTime,floatStartTime,floatEndTime,&annuity,&swapRate,fixPayTimes,fixPayPeriods,tempStates);
	
#if defined(__GP_STRICT_VALIDATION)
	if (annuity < K_NEW_DOUBLE_TOL)
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
			"Annuity is Null at 0!");
# endif

	double Eta = 	ComputeEta( evalTime, 0.) ;
	double result = fabs((swapDerivative * Eta)/swapRate) ;

	return result ;
//	return (GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(evalTime));	
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_MarkovSV1F
///	Routine: ComputeEquivalentShiftt
///	Returns: a double
///	Action : Compute the equivalent Lambdat in the Heston Swap Rate diffusion
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double ARM_MarkovSV1F::ComputeEquivalentShiftt( 
		const string& curveName,
		double evalTime, 
		double floatStartTime,
		double floatEndTime,
		const ARM_GP_Vector& fixPayTimes,
		const ARM_GP_Vector& fixPayPeriods,
		double lambdaValue) const
{
	ARM_PricingStatesPtr tempStates(NULL);
	//// for derivative calculation
	double epsilon = DERIVATIVE_EPS;

	double Eta_0	= ComputeEta( evalTime, 0.) ;
	double Eta_eps1	= ComputeEta( evalTime, epsilon ) ;
	double Eta_eps2	= ComputeEta( evalTime,  - epsilon ) ;

	double annuity,swapRate;

	double SwapRate_Derivative_eps2	= ComputeDerivativeSt(curveName,evalTime,floatStartTime,floatEndTime,&annuity,&swapRate,fixPayTimes,fixPayPeriods,tempStates,- epsilon );
	double SwapRate_Derivative_eps1	= ComputeDerivativeSt(curveName,evalTime,floatStartTime,floatEndTime,&annuity,&swapRate,fixPayTimes,fixPayPeriods,tempStates,epsilon);
	double SwapRate_Derivative_0	= ComputeDerivativeSt(curveName,evalTime,floatStartTime,floatEndTime,&annuity,&swapRate,fixPayTimes,fixPayPeriods,tempStates);

	double factor = swapRate * lambdaValue * lambdaValue * SwapRate_Derivative_0;
	double swapVolDerivative_2 = 	(SwapRate_Derivative_eps1 * Eta_eps1) * (SwapRate_Derivative_eps1 * Eta_eps1) - (SwapRate_Derivative_eps2 * Eta_eps2) * (SwapRate_Derivative_eps2 * Eta_eps2);
	swapVolDerivative_2/=( 2 * epsilon);

#if defined(__GP_STRICT_VALIDATION)
	if( factor < K_NEW_DOUBLE_TOL) 
		ARM_THROW( ERR_INVALID_ARGUMENT, " Volatility, SwapRate or SwapRate Derivative is NULL!!" );
#endif

	double result = swapVolDerivative_2/(2 * factor) ;
	/// For Test Purpose
//	return (GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->Interpolate(evalTime));	
	return result;
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: LocalDiscounts
///	Returns : void
///	Action  : Computes the LocalDiscounts
////////////////////////////////////////////////////
ARM_VectorPtr ARM_MarkovSV1F::LocalDiscounts(
	size_t timeIdx, 
	double dt, 
	const ARM_PricingStatesPtr& states) const
{
	const ARM_GP_Vector* const timeSteps = GetNumMethod()->GetTimeSteps();

	// Compute the deterministic part
    ARM_ZeroCurvePtr ZcCurve = GetZeroCurve();
	double startTime		= (*timeSteps)[timeIdx];
	double endTime			= startTime + dt;

	double zcStart	= ZcCurve->DiscountPrice(startTime/K_YEAR_LEN);
    double zcEnd	= ZcCurve->DiscountPrice(endTime/K_YEAR_LEN);
	double discount = zcEnd/zcStart;

	size_t statesSize		= states->size();
	ARM_GP_Vector* result	= new ARM_GP_Vector(statesSize,0.0);
	size_t modelNb			= GetModelNb();

    /// Simple mapping version
    double Xt;
    dt /= K_YEAR_LEN;
	for( size_t i=0; i<statesSize; ++i )
	{
		Xt = states->GetModelState(i,modelNb+FIRST_STATE_VARIABLE);
		(*result)[i] = discount * exp(-Xt*dt);
    }

	return ARM_VectorPtr(result);
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_MarkovSV1F::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
							  ARM_ModelParam* inputModelParam, 
							  size_t factorNb )
{
	ARM_CurveModelParam* modelParam = dynamic_cast<ARM_CurveModelParam*>(inputModelParam);
	if( !modelParam )
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "expected an ARM_CurveModelParam!");

    double asOfDate = GetAsOfDate().GetJulian();
    int size1       = portfolio->GetSize();  
    ARM_GP_Vector  tmpdates;
    int i;
    
    switch( modelParam->GetType() )
    {
    case ARM_ModelParamType::InitialVol:
	case ARM_ModelParamType::Shift:
	case ARM_ModelParamType::VolOfVol:
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
            "Unknown type... Model Param Not Supported by MSV1F" );
    }
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: EulerLocalDrifts
///	Returns :
///	Action  : computes the relative and absolute drift
////////////////////////////////////////////////////

void ARM_MarkovSV1F::EulerLocalDrifts( const ARM_GP_Vector& timeSteps,
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
///	Class   : ARM_MarkovSV1F
///	Routines: VolatilitiesAndCorrelations
///	Returns :
///	Action  : computes the volatilities its derivatives and the correlation
////////////////////////////////////////////////////
void ARM_MarkovSV1F::VolatilitiesAndCorrelations( const ARM_GP_Vector& timeSteps, 
	ARM_GP_MatrixPtr& vols,
	ARM_GP_MatrixPtr& d1Vols,
	ARM_GP_MatrixPtr& correls,
	bool linearVol) const
{
}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: IntegratedBondSquaredVol
///	Returns : double
///	Action  : Int_startTime^endTime Gamma(s,bondMaturity)^2 ds 
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_MarkovSV1F::IntegratedBondSquaredVol( double startTime, double endTime, double bondMaturity ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"IntegratedBondSquaredVol : unimplemented function for ARM_MarkovSV Model!");
}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: IntegratedBondCovariance
///	Returns : double
///	Action  : Int_startTime,endTime,gamma(s,bondMaturity1)*gamma(s,bondMaturity2)ds
///      Where dB(t,T)/B(t,T) = r dt + Gamma(t,T) dW_t
////////////////////////////////////////////////////
double ARM_MarkovSV1F::IntegratedBondCovariance( double startTime, double endTime, double bondMaturity1, double bondMaturity2 ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"IntegratedBondCovariance : unimplemented function for ARM_MarkovSV Model!");

}

////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: VolatilityScalarProduct
///	Returns : double (qui l'eut cru)
///	Action  :  Int_startTime^endTime Gamma(s,bondMaturity) * dW_s 
////////////////////////////////////////////////////
double ARM_MarkovSV1F::VolatilityScalarProduct( double startTime, double endTime, double bondMaturity, const ARM_ModelParam& otherModelVolatility ) const
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"VolatilityScalarProduct : unimplemented function for ARM_MarkovSV Model!");

}


////////////////////////////////////////////////////
///	Class   : ARM_MarkovSV1F
///	Routines: ModelStateLocalCorrels
///	Returns : void
///	Action  :  
////////////////////////////////////////////////////
void ARM_MarkovSV1F::ModelStateLocalCorrels( const ARM_GP_Vector& timeSteps,
		ARM_MatrixVector& localCorrels, const ARM_MultiAssetsModel& multiassets )
{
	throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
		"TreeStatesToModelStates : unimplemented function for ARM_MarkovSV Model!");


}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

