/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Heston_Model.cpp
 *
 *  \brief base class for Heston Model
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */

#include "gpmodels/ShiftedHeston_Model.h"

/// gpbase
#include "gpbase/eventviewerfwd.h"
#include "gpbase/removeidentifiedwarning.h"

#include "gpbase/curve.h"
#include "gpbase/singleton.h"
#include "gpbase/numericconstant.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelparamtype.h"
#include "gpinfra/modelparam.h"
#include "gpinfra/curvemodelparam.h"

/// gpclosedforms
#include "gpclosedforms/heston_interface.h"
#include "gpnumlib/gaussiananalytics.h"

/// gpmodels
#include "gpmodels/ShiftedHeston_ModelParams.h"

/// gpnumlib
#include "gpnumlib/solver.h"
#include "gpnumlib/brent.h"
#include "gpnumlib/boxmuller.h"
#include "gpnumlib/uniform.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"
#include "gpnumlib/argconvdefault.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/bootstrapnd.h"
#include "gpcalib/vanillaarg.h"

/// kernel
#include <inst/portfolio.h>

CC_BEGIN_NAMESPACE( ARM )

const double DEFAULT_NBOSCILLATIONS		    = 3; 
const double DEFAULT_PRECISION				= 1.e-6;
const double DEFAULT_NBSTAGE				= -1;
const double DEFAULT_NBPOINTS_STAGE1		= 120;
const double DEFAULT_NBPOINTS_PER_STAGE		= 30;

#define FIRST_STATE_VARIABLE    0 // for Xt
#define SECOND_STATE_VARIABLE   1 // for Yt




////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_ShiftedHeston_Model
////////////////////////////////////////////////////

ARM_ShiftedHeston_Model::ARM_ShiftedHeston_Model(const ARM_ZeroCurvePtr& zc, const ARM_ShiftedHeston_ModelParams& params, size_t IntegrationStep , bool isMcPrice, size_t nbSteps, size_t nbSimul )
:	ARM_AnalyticIRModel( zc, params ), itsIntegrationStep(IntegrationStep),itsIsMCPrice(isMcPrice),itsNbStepsMC(nbSteps),itsNbSimulations(nbSimul)
{
}

////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////

ARM_ShiftedHeston_Model::ARM_ShiftedHeston_Model(const ARM_ShiftedHeston_Model& rhs)
:	ARM_AnalyticIRModel( rhs ), 
	itsIntegrationStep( rhs.itsIntegrationStep),
	itsIsMCPrice(rhs.itsIsMCPrice),
	itsNbStepsMC(rhs.itsNbStepsMC),
	itsNbSimulations(rhs.itsNbSimulations)
{}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_ShiftedHeston_Model& ARM_ShiftedHeston_Model::operator=(const ARM_ShiftedHeston_Model& rhs)
{
	if( this != &rhs )
	{
		ARM_AnalyticIRModel::operator =(rhs);
		itsIntegrationStep	= rhs.itsIntegrationStep;
		itsIsMCPrice	= rhs.itsIsMCPrice;
		itsNbStepsMC	= rhs.itsNbStepsMC;
		itsNbSimulations	= rhs.itsNbSimulations;
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_ShiftedHeston_Model::~ARM_ShiftedHeston_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_ShiftedHeston_Model::VanillaCaplet(
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
	ARM_VectorPtr liborRate	= Libor(curveName, evalTime,fwdStartTime, 
		fwdEndTime, period, fwdResetTime, payTime, states );
	ARM_VectorPtr DF		= DiscountFactor( curveName, evalTime, payTime, states );
	
	double time				= fwdResetTime-evalTime;
	double tenor			= ARM_AnalyticIRModel::GetMatchingTenor((fwdEndTime-fwdStartTime)/K_YEAR_LEN);

	double s0				= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol))).GetCurve()->Interpolate(time);
	double V0				= s0*s0;
	double longtermS		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol))).GetCurve()->Interpolate(time);
	double longtermV		= longtermS*longtermS;
	double speed			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetCurve()->Interpolate(time);
	double volvol			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol))).GetCurve()->Interpolate(time);
	double rho				= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation))).GetCurve()->Interpolate(time);
	double shift			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::Shift))).GetCurve()->Interpolate(time);

	double optionValue = payNotional* (*DF)[0] * period * Export_Shifted_Heston_VanillaOption(
		(*liborRate)[0], strikesPerState[0], V0,
		time/K_YEAR_LEN, longtermV, speed, volvol, rho, shift, capFloor, DEFAULT_NBPOINTS_STAGE1,DEFAULT_NBPOINTS_PER_STAGE, 
		DEFAULT_NBSTAGE,DEFAULT_NBOSCILLATIONS,DEFAULT_PRECISION);

	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}



////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routines: VanillaSwaption
///	Returns :
///	Action  : computes a vanilla swaption
////////////////////////////////////////////////////

ARM_VectorPtr ARM_ShiftedHeston_Model::VanillaSwaption(
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


	/// not necesary to use  fwdStartTimes, fwdEndTimes, fwdPayPeriods, floatPayTimes, floatPayPeriods
	/// as we are in the case of a vanilla swap!
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

	double time			= swapResetTime-evalTime;
	double tenor		= ARM_AnalyticIRModel::GetMatchingTenor((floatEndTime-floatStartTime)/K_YEAR_LEN );
	double s0			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::InitialVol))).GetCurve()->Interpolate(time);
	double V0			= s0*s0;
	double longtermS	= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::LongTermVol))).GetCurve()->Interpolate(time);
	/// To be Updated
	/// We fix the longTermVol at InitialVol
	double longtermV	= V0;
	double speed		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolMeanReversion))).GetCurve()->Interpolate(time);
	double volvol		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::VolOfVol))).GetCurve()->Interpolate(time);
	double rho			= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation))).GetCurve()->Interpolate(time);
	double shift		= ((ARM_CurveModelParam&) (GetModelParams()->GetModelParam(ARM_ModelParamType::Shift))).GetCurve()->Interpolate(time);
  
	bool isMCPrice = true;
	ARM_VectorPtr values; 
	double MCValue = 0.;
	double MCValue_2 = 0.;
	double stdDev;

	if(itsIsMCPrice)
	{
		values = MCPrice(swapResetTime,(*swapRate)[0],strikesPerState(0,0),callPut);
		////// Option Value
		int nbPaths = values->size();
		for (size_t i=0; i<nbPaths; ++i)
		{
			MCValue+= (*values)[i];
			MCValue_2+=(*values)[i]*(*values)[i];
		}
		MCValue*=(*annuity)[0] * swapNotional;
		MCValue_2*=(*annuity)[0]*(*annuity)[0];
		MCValue_2/=nbPaths;
		MCValue/= nbPaths;
		////// StdDev Estimation
		stdDev = sqrt((swapNotional * swapNotional * MCValue_2 - MCValue*MCValue)/(nbPaths-1));
	

		CC_Ostringstream os;

		os.flush();
		os << "Option Value  = " << MCValue << std::endl;
		os << "Option StdDev = " << stdDev << std::endl;
		os << "End Of Pricing"<< std::endl;
		ARM_TheEventViewer.Instance()->AddToMessage(os.str());
		return ARM_VectorPtr( new std::vector<double>(1,MCValue));

	}
	volvol*=s0;
	double 	optionValue = swapNotional * (*annuity)[0] * Export_Shifted_Heston_VanillaOption(
			(*swapRate)[0],strikesPerState(0,0), V0,
			time/K_YEAR_LEN, longtermV, speed, volvol, rho, shift, callPut, DEFAULT_NBPOINTS_STAGE1,DEFAULT_NBPOINTS_PER_STAGE, 
			DEFAULT_NBSTAGE,DEFAULT_NBOSCILLATIONS,DEFAULT_PRECISION);
//	}
	return ARM_VectorPtr( new std::vector<double>(1,optionValue));
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: MCPrice
///	Returns: void
///	Action : WARNING!!!!! We suppose the correlation = 1
//////////////////////////////////////////////////////////////////////////////////////////////////////////
ARM_VectorPtr ARM_ShiftedHeston_Model::MCPrice(double swapResetTime, double swapRate_Zero, double strike, int callPut) const //ARM_PricingStatesPtr& states,int timeIndex) const
{	
	// MRGK5 Random Generator
	ARM_RandomGeneratorPtr  pBaseRandomGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
				(ARM_RandGenFactoryImp::BaseGenType)	 ARM_ArgConv_BaseGenAlgoType.GetNumber("NR_RAN2"),
				ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	ARM_RandomGeneratorPtr normRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm, ARM_RandGenFactoryImp::BoxMuller,	pBaseRandomGen ) );

	/// antithetic variates!
	ARM_RandomGeneratorPtr numRandGen = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
			ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm, ARM_RandGenFactoryImp::AntitheticOne,normRandGen ) );



	int nbPaths = itsNbSimulations;
	int nbSteps = itsNbStepsMC;
	nbSteps= floor(nbSteps * swapResetTime/K_YEAR_LEN);
	double dt  = swapResetTime/K_YEAR_LEN/nbSteps;
	////// We do the MonteCarlo

	std::vector<double> S_State (nbPaths,swapRate_Zero);
	std::vector<double> Z_State (nbPaths,1.0);
	double time = 0;

	//// We handl only time-dependent shift and volatility 
	//// Test of the MSV purpose

	////// We have constant Mean Reversion
	double VolMeanReversion_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::VolMeanReversion).GetValueAtPoint(0);
	double VolOfVol_t = GetModelParams()->GetModelParam( ARM_ModelParamType::VolOfVol).GetValueAtPoint(0);
	double LongTermVol_t = 1.0;

	////// We suppose that the long term vol = 1
	////// We develop a MonteCarlo of the form dS = Sigma (m St + (1-m) S0) sqrt(Vt) dWt
	//////									   dVt = a(V0 - Vt) dt + eps sqrt(Vt) dW2t
	
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


	for( size_t j=0;j<nbSteps; ++j )
	{
		std::vector<double>& gaussian_Vector_1 = (*normRandGen)(nbPaths);
		std::vector<double>& gaussian_Vector_2 = (*normRandGen)(nbPaths);
		std::vector<double>& gaussian_Vector_3 = (*normRandGen)(nbPaths);
//		std::vector<double>& uniform_Vector = (*pBaseRandomGen)(nbPaths);

		for( size_t i=0;i<nbPaths; ++i )
		{			
			/// First Variable St
			double gaussian_1 = (*gaussian_Vector_1)[i];
			double gaussian_2 = (*gaussian_Vector_2)[i];

			

			/// We have to Interpolate
			double Volatility_t = GetModelParams()->GetModelParam( ARM_ModelParamType::InitialVol).ToCurveModelParam().GetCurve()->Interpolate(time * K_YEAR_LEN);
			double Shift_t  = GetModelParams()->GetModelParam( ARM_ModelParamType::Shift).ToCurveModelParam().GetCurve()->Interpolate(time * K_YEAR_LEN);
			

		/*	double gaussian_3 = (*gaussian_Vector_3)[i];//gaussian_1 * gaussian_2;//
			double uniform = ARM_GaussianAnalytics::cdfNormal(gaussian_3);//(*uniform_Vector)[i];
			if ((uniform<0.001) || (uniform>0.999))
				uniform = ARM_GaussianAnalytics::cdfNormal(gaussian_2);

			double R_2 = dt * (gaussian_1 *gaussian_1 + gaussian_2 * gaussian_2);
			double X = dt/ARM_NumericConstants::ARM_PI * log(uniform / (1-uniform));
			double Y = sqrt(R_2*dt/3) * gaussian_3;
			*/
			
			/////// Exact Fourier Levy Formula
			S_State[i] = S_State[i] +  (S_State[i] * Shift_t + (1 - Shift_t) * swapRate_Zero) * Volatility_t * ((sqrt(dt)*gaussian_1*sqrt(Z_State[i]) ));// + 0.5 * Z_State[i] * Volatility_t * dt * (gaussian_1*gaussian_1 - 1) + 0.25  * VolOfVol_t * (gaussian_1 * gaussian_2 * dt - X - Y));
			
			/////// Standard Euler Approach
			//S_State[i] = S_State[i] +  (S_State[i] * Shift_t + (1 - Shift_t) * swapRate_Zero) * Volatility_t * sqrt(dt)*gaussian_1*sqrt(Z_State[i]) ; 

			/////// 1D Milstein 
			//S_State[i] = S_State[i] +  (S_State[i] * Shift_t + (1 - Shift_t) * swapRate_Zero) * Volatility_t * ( sqrt(dt) * gaussian_1 * sqrt(Z_State[i])  + 0.5 * Z_State[i] * Volatility_t * dt * (gaussian_1*gaussian_1 - 1) );

			/////// 2D Milstein  Simple Approx of the Integral I12 = 0.5 * dW1 * dW2
			//S_State[i] = S_State[i] +  (S_State[i] * Shift_t + (1 - Shift_t) * swapRate_Zero) * Volatility_t * ((sqrt(dt)*gaussian_1*sqrt(Z_State[i]) ) + 0.5 * Z_State[i] * Volatility_t * dt * (gaussian_1*gaussian_1 - 1) + 0.25  * VolOfVol_t * gaussian_1 * gaussian_2 * dt);



			/// Third Variable Zt
			/////Standard Euler Discretisation approach  + Milstein
		/*	Z_State[i] = Z_State[i] + VolMeanReversion_dt * (1- Z_State[i]) + VolOfVol_dt * sqrt(Z_State[i])*gaussian_2 + 0.25 * VolOfVol_2_t * dt*(gaussian_2*gaussian_2 - 1);
			if (Z_State[i] <0)
				Z_State[i] = - Z_State[i];*/
			
			//// New Approach : Andreasen
			//// Intermediate Variables
/*			double gaussian_2 = (*gaussian_Vector_2)[i+ j*nbPaths];		
			double Z_State_bis = Z_State[i]*var_1 + var_2;	
			double var_i = 1 + 1/(Z_State_bis*Z_State_bis)*(var_5+(Z_State[i]-1)*var_4);
			double Nu_2 = log(fabs(var_i-1)+1);
			double Nu = sqrt(Nu_2);
			
			Z_State[i] = Z_State_bis*exp(-0.5*Nu_2+ Nu*gaussian_2);*/

			double Z_State_bis = Z_State[i] * var_1 + var_2;	
			double var_i = 1 + 1/(Z_State_bis*Z_State_bis)* var_7 * Z_State[i];
			double Nu_2 = log(var_i);
			double Nu = sqrt(Nu_2);

			Z_State[i] = Z_State_bis*exp(-0.5*Nu_2+ Nu*gaussian_2);
		}
		time += dt;
		delete gaussian_Vector_1;
		delete gaussian_Vector_2;
		delete gaussian_Vector_3;
	}
	std::vector<double>& values = new std::vector<double>(nbPaths);
	double optionValue = 0.;
	for(size_t h=0;h<nbPaths; ++h)
	{
		double value = S_State[h] - strike;
		if (callPut == 1)
			(*values)[h] = (value>0.0) ? value : 0.;
		else
			(*values)[h] = (value<0.0) ? -value : 0.;
	}
	optionValue/=nbPaths;
	return ARM_VectorPtr( values );
}


    
////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routines: Clone
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_ShiftedHeston_Model::Clone() const
{
	return new ARM_ShiftedHeston_Model(*this);
}


////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routines: void 
///	Returns :
///	Action  : sets the corresponding suggested break point times to the model param
///////////////////////////////////////
void ARM_ShiftedHeston_Model::AdviseBreakPointTimes( const ARM_StdPortfolioPtr& portfolio, 
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
    case ARM_ModelParamType::InitialVol:
	case ARM_ModelParamType::Shift:
	case ARM_ModelParamType::VolOfVol:
	case ARM_ModelParamType::VolMeanReversion:
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
            "Unknown type... a Shifted Heston Model only supports InitialVol" );
    }
}


////////////////////////////////////////////////////
///	Class  : ARM_ShiftedHeston_Model
///	Routine: PreProcessing
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
void ARM_ShiftedHeston_Model::PreProcessing(ARM_ModelFitter& modelFitter)
{ 
	GetModelParams()->PreProcessing(modelFitter,modelFitter.GetFactorNb());
}


////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routine : ValidateCalibMethod
///	Returns :void
///	Action  : call DefaultValidateWithModel
////////////////////////////////////////////////////
void ARM_ShiftedHeston_Model::ValidateCalibMethod(ARM_CalibMethod& calibMethod)
{
	calibMethod.DefaultValidateWithModel(*this);
}

////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_ShiftedHeston_Model::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Generalized Shifted Heston Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}



////////////////////////////////////////////////////
///	Class   : ARM_ShiftedHeston_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_ShiftedHeston_Model::ValidateModelParams(const ARM_ModelParams& params) const
{
	const ARM_ShiftedHeston_ModelParams* ShiftedHeston_ModelParams = dynamic_cast<const ARM_ShiftedHeston_ModelParams*>(&params);
	if( !ShiftedHeston_ModelParams )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_ShiftedHeston_ModelParams" );
	return true;
}
CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

