/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file Local_Model.cpp
 *
 *  \brief base class for local analytical model
 *	\author  A. Chaix
 *	\version 1.0
 *	\date June 2005
 */

/// remove identified warning
#include "gpbase/removeidentifiedwarning.h"

// nag
#include "nag.h"
#include "nage04.h"

/// gpmodels
#include "gpmodels/local_normal_model.h"
#include "gpmodels/Local_Normal_ModelParams.h"
#include "gpmodels/ForwardMarginBasis.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/zccrvfunctor.h"
#include "gpinfra/surfacelistmodelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/irrate.h"
#include "gpinfra/modelnamemap.h"

/// gpbase
#include "gpbase/curve.h"
#include "gpbase/surface.h"
#include "gpbase/datestrip.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/vectormanip.h"
#include "gpbase/timer.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/normal.h"

// gpmodels (calibration purpose)
#include "gpmodels/hw.h"
#include "gpmodels/hw1f.h"
#include "gpmodels/hw2f.h"
#include "gpmodels/qgm1f.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/modelparamshw2f.h"
#include "gpmodels/multiassets.h"

// gpcalib
#include "gpcalib/densityfunctors.h"

// kernel
#include <inst/spreadoption.h>
#include <inst/corridordblcondition.h>


/// gpinflation
#include "gpinflation/infcapfloor.h"

/// gpnummethods (for dump only)
#include "gpnummethods/treebase.h"

// gpcalculators
#include "gpcalculators/fxracalculator.h"

///----------------------------------------------------------------------
///-- macros to handle fixing curves (CSO with foreign funding)
///----------------------------------------------------------------------
#define CHANGE_MODEL_FUNCTORS() \
    	ARM_ZeroCurveFunctor* oldDiscountFunctor = GetNumericalModel()->GetDiscountFunctor();\
        bool isSameMod = GetFixingFunctor()->IsSameModel(*GetDiscountFunctor());\
    if(!isSameMod){\
		const_cast<ARM_PricingModel*>( GetNumericalModel() )->SetDiscountFunctor(GetDiscountFunctor() ); \
	}
#define CAST_TO_IRMODEL() \
	ARM_PricingModelIR* IRMODEL = dynamic_cast<ARM_PricingModelIR*>( const_cast<ARM_PricingModel *>( GetNumericalModel() ));  \
	if( !IRMODEL ) \
		ARM_THROW(  ERR_INVALID_ARGUMENT, "could not find a numerical model for " + curveName + "!" );

#define RESTORE_MODEL_FUNCTORS() \
	if(!isSameMod){\
		const_cast<ARM_PricingModel*>( GetNumericalModel() )->SetDiscountFunctor(oldDiscountFunctor);\
	}

CC_BEGIN_NAMESPACE( ARM )


///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ
////////////////////////////////////////////////////////
/// ARM_Local_Normal_Model implementation //////////////
////////////////////////////////////////////////////////
///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ

///----------------------------------------------------------------------
///-- nb stev above which we do not try to calibrate forward vol
///----------------------------------------------------------------------
const double ARM_Local_Normal_Model::NSTDEV_NO_CALIB = 8;

///-----------------------------------------------
///--- Indexes in surface lists where
///--- are stored local vols & adjustments
///-----------------------------------------------
///--------------------------------------------
/// Libor & VanillaCaplet
const size_t ARM_Local_Normal_Model::CAPLET_ADJ		 = 0;
const size_t ARM_Local_Normal_Model::CAPLET_ADJ_SIZE = 1;
/// - - - - - - - - - - - 
const size_t ARM_Local_Normal_Model::CAPLET_VOL		 = 0;
const size_t ARM_Local_Normal_Model::CAPLET_VOL_SIZE = 1;
///--------------------------------------------
/// VanillaSpreadOptionlet
const size_t ARM_Local_Normal_Model::SO_ADJ_LONG	= 0;
const size_t ARM_Local_Normal_Model::SO_ADJ_SHORT	= 1;
const size_t ARM_Local_Normal_Model::SO_ADJ_SIZE	= 2;
/// - - - - - - - - - - - 
const size_t ARM_Local_Normal_Model::SO_VOL			= 0;
const size_t ARM_Local_Normal_Model::SO_VOL_SIZE	= 1;
///--------------------------------------------
/// VanillaCMSCorridorlet (single condition)
const size_t ARM_Local_Normal_Model::SOCOR_DOWN		= 0;
const size_t ARM_Local_Normal_Model::SOCOR_UP		= 1;

const size_t ARM_Local_Normal_Model::SOCOR_ADJ_LONG				= 0;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_SHORT			= 1;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_LONG_UP_FLT		= 2;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_SHORT_UP_FLT		= 3;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_LONG_DOWN_FLT	= 4;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_SHORT_DOWN_FLT	= 5;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_PAY_FLT			= 6;
const size_t ARM_Local_Normal_Model::SOCOR_ADJ_SIZE				= 7;
/// - - - - - - - - - - - 
const size_t ARM_Local_Normal_Model::SOCOR_VOL_UPLEFT		= 0;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_UPRIGHT		= 1;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_DOWNLEFT		= 2;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_DOWNRIGHT	= 3;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_UPLEFT_FLT	= 4;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_UPRIGHT_FLT	= 5;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_DOWNLEFT_FLT	= 6;
const size_t ARM_Local_Normal_Model::SOCOR_VOL_DOWNRIGHT_FLT= 7;

const size_t ARM_Local_Normal_Model::SOCOR_VOL_SIZE			= 8;
/// - - - - - - - - - - - 
const double ARM_Local_Normal_Model::SOCOR_STRIKE_SPREAD    = 0.0001; /// 1 BP
/// - - - - - - - - - - - 
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_UPLEFT					= 0;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_UPRIGHT					= 1;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_DOWNLEFT				= 2;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_DOWNRIGHT				= 3;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_UPLEFT_FLT				= 4;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_UPRIGHT_FLT				= 5;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_DOWNLEFT_FLT			= 6;
const size_t ARM_Local_Normal_Model::SOCOR_KADJ_DOWNRIGHT_FLT			= 7;

// One more surface to save market target residual RSO let at each (notice,reset idx) dates
const size_t ARM_Local_Normal_Model::SOCOR_MKT_RESIDUAL_RSO_PRICE		= 8;

// One more surface to save multiplicative PV correction for residual RSO let at at each (notice,reset flow) dates
const size_t ARM_Local_Normal_Model::SOCOR_RESIDUAL_RSO_PV_CORRECTION	= 9;

const size_t ARM_Local_Normal_Model::SOCOR_KADJ_SIZE					= 10;
///--------------------------------------------

///--------------------------------------------
/// VanillaCMSCorridorlet (double condition)

const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN	= 0;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP		= 1;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN		= 2;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP		= 3;

const size_t ARM_Local_Normal_Model::DBLECOR_NB_BARRIERS	= 4;

/// Adj rates
const size_t ARM_Local_Normal_Model::DBLECOR_ADJ_SPREAD							= 0;
const size_t ARM_Local_Normal_Model::DBLECOR_ADJ_RATE							= 1;

const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_ADJ_SPREAD_FLT			= 2;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_ADJ_SPREAD_FLT			= 3;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_ADJ_SPREAD_FLT			= 4;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_ADJ_SPREAD_FLT				= 5;

const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_ADJ_RATE_FLT			= 6;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_ADJ_RATE_FLT				= 7;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_ADJ_RATE_FLT				= 8;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_ADJ_RATE_FLT				= 9;

const size_t ARM_Local_Normal_Model::DBLECOR_ADJ_PAY_FLT						= 10;
														
const size_t ARM_Local_Normal_Model::DBLECOR_ADJ_SIZE							= 11;
														
/// Adj strikes											
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_KADJ_SPREAD			= 0;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_KADJ_SPREAD				= 1;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_KADJ_SPREAD				= 2;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_KADJ_SPREAD				= 3;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_KADJ_RATE				= 4;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_KADJ_RATE				= 5;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_KADJ_RATE				= 6;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_KADJ_RATE					= 7;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_KADJ_SPREAD_FLT		= 8;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_KADJ_SPREAD_FLT			= 9;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_KADJ_SPREAD_FLT			= 10;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_KADJ_SPREAD_FLT			= 11;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_KADJ_RATE_FLT			= 12;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_KADJ_RATE_FLT			= 13;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_KADJ_RATE_FLT			= 14;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_KADJ_RATE_FLT				= 15;
														
const size_t ARM_Local_Normal_Model::DBLECOR_KADJ_SIZE							= 16;
														
/// Adj vols											
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_VOL_SPREAD				= 0;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_VOL_SPREAD				= 1;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_VOL_SPREAD				= 2;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_VOL_SPREAD					= 3;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_VOL_RATE				= 4;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_VOL_RATE					= 5;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_VOL_RATE					= 6;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_VOL_RATE					= 7;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_VOL_SPREAD_FLT			= 8;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_VOL_SPREAD_FLT			= 9;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_VOL_SPREAD_FLT			= 10;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_VOL_SPREAD_FLT				= 11;
														
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_VOL_RATE_FLT			= 12;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_VOL_RATE_FLT				= 13;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_VOL_RATE_FLT				= 14;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_VOL_RATE_FLT				= 15;
														
const size_t ARM_Local_Normal_Model::DBLECOR_VOL_SIZE							= 16;

/// Adj correls
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL		= 0;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL		= 1;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL		= 2;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_SPREAD_RATE_CORREL			= 3;

const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL_FLT	= 4;
const size_t ARM_Local_Normal_Model::DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL_FLT	= 5;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL_FLT	= 6;
const size_t ARM_Local_Normal_Model::DBLECOR_SUP_RUP_SPREAD_RATE_CORREL_FLT		= 7;

// One more surface to save market target residual RA2 let at each (notice,reset idx) dates
const size_t ARM_Local_Normal_Model::DBLECOR_MKT_RESIDUAL_RA2_PRICE				= 8;

// One more surface to save multiplicative PV correction for residual RA2 let at at each (notice,reset flow) dates
const size_t ARM_Local_Normal_Model::DBLECOR_RESIDUAL_RA2_PV_CORRECTION			= 9;


const size_t ARM_Local_Normal_Model::DBLECOR_CORREL_SIZE						= 10;
///--------------------------------------------

/// Gaussian correlations for copula based calibration
const size_t ARM_Local_Normal_Model::COPULA_MARKET_CORREL						= 0;
const size_t ARM_Local_Normal_Model::COPULA_MODEL_CORREL						= 1;
const size_t ARM_Local_Normal_Model::COPULA_CORREL_SIZE							= 2;


const double LOCAL_SPREAD_VOL_MIN	= 0.0002; /// 2bp of normal vol
const double LOCAL_RATE_VOL_MIN		= 0.0005; /// 5bp of normal vol

//const double COR_OPTIM_EVAL_AVE_YEARLEN		= 365.25;
//const double COR_OPTIM_EVAL_DAILYLIMIT		= 2.0/52.0*COR_OPTIM_EVAL_AVE_YEARLEN; // 2W
//const double COR_OPTIM_EVAL_WEEKLYLIMIT		= 1.5/12.0*COR_OPTIM_EVAL_AVE_YEARLEN; // 1,5M

#define PRICGLPOINTSNB	120
#define INVSQRT2PI		0.398942280401433
#define INVSQRT2		0.707106781186547

// Functional extrapolation allowed 6M before 1st or 6M after last calibrated reset
const double MAX_FUNCTIONAL_MISMATCH = 182.0;

/****
double lastEval=0.0;
size_t nbEvals=0;
double computingTime=0.0;
double computingTimeAll=0.0;
****/

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Local_Normal_Model
////////////////////////////////////////////////////
ARM_Local_Normal_Model::ARM_Local_Normal_Model(const ARM_ZeroCurvePtr& zc, const ARM_Local_Normal_ModelParams& params)
:	ARM_Local_Model( zc, params ),
	itsParamTypes(0),itsParamNbSurfaces(0),
	itsLiborAsShiftedLognormal (true)
{
	if( GetModelParams()->DoesModelParamExist(ARM_ModelParamType::ForwardAdjustment) &&
		GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Volatility) &&
		GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Shift) )
	{
		/// Default parameters
		itsParamTypes.resize(3);
		itsParamTypes[0] = ARM_ModelParamType::ForwardAdjustment;
		itsParamTypes[1] = ARM_ModelParamType::Volatility;
		itsParamTypes[2] = ARM_ModelParamType::Shift;

		itsParamNbSurfaces.resize(3);
		itsParamNbSurfaces[0] = SOCOR_ADJ_SIZE;
		itsParamNbSurfaces[1] = SOCOR_VOL_SIZE;
		itsParamNbSurfaces[2] = SOCOR_KADJ_SIZE;
	}
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Local_Normal_Model
////////////////////////////////////////////////////
ARM_Local_Normal_Model::ARM_Local_Normal_Model(const ARM_ZeroCurvePtr& zc, const ARM_Local_Normal_ModelParams& params,
						const ARM_IntVector& paramTypes, const ARM_IntVector& paramNbSurfaces)
:	ARM_Local_Model( zc, params ),
	itsLiborAsShiftedLognormal (true),
	itsParamTypes(paramTypes),
	itsParamNbSurfaces(paramNbSurfaces)
{
	if(paramTypes.size() != paramNbSurfaces.size())
		ARM_THROW(  ERR_INVALID_ARGUMENT, " : param types and param nb surfaces must be equal" );
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_Local_Normal_Model::ARM_Local_Normal_Model(const ARM_Local_Normal_Model& rhs)
:	ARM_Local_Model (rhs),
	itsLiborAsShiftedLognormal (rhs.itsLiborAsShiftedLognormal),
	itsParamTypes(rhs.itsParamTypes),
	itsParamNbSurfaces(rhs.itsParamNbSurfaces)

{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_Local_Normal_Model& ARM_Local_Normal_Model::operator = (const ARM_Local_Normal_Model& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_Local_Normal_Model();
		new (this) ARM_Local_Normal_Model (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Local_Normal_Model::~ARM_Local_Normal_Model()
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: Libor
///	Returns: ARM_VectorPtr
///	Action : computes the Libor rate in local normal model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const 
{
	if (GetFunctionalFlag())
	{
		CHANGE_MODEL_FUNCTORS();
		CAST_TO_IRMODEL();
	
		ARM_VectorPtr stdLibor = IRMODEL->Libor(	curveName, 
												evalTime,
												fwdStartTime, 
												fwdEndTime,
												period,
												resetTime,
												fwdEndTime,
												states );

		ARM_VectorPtr newLibor = Func(evalTime,stdLibor);
	
		RESTORE_MODEL_FUNCTORS();
		return newLibor;
	}
	else
	{
		// downcast to get functor
		ARM_PricingFunctionIR* irFunctor = dynamic_cast < ARM_PricingFunctionIR* > (GetNumericalModel());

		if( !irFunctor )
			ARM_THROW(  ERR_INVALID_ARGUMENT, "IR functions are not supported by numerical model" );

		
		// compute standard libor rate from model (payTime = fwdEndTime)
		ARM_VectorPtr liborRate	= irFunctor->Libor(curveName, evalTime,fwdStartTime, 
															 fwdEndTime, period, resetTime, fwdEndTime, states );

		// adjust forward (additive, same adjustment in all states)
		// libor adjustment supposed to be stored in first surface of the list
		double adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(CAPLET_ADJ, evalTime, resetTime);
		(*liborRate) += adjustment;

		return liborRate;
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::VanillaCaplet(
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

	CHANGE_MODEL_FUNCTORS();
	CAST_TO_IRMODEL();
		
	// compute libor rate with previous method (additive adjustement on standard FRA)
	ARM_VectorPtr liborRate = ARM_Local_Normal_Model::Libor (curveName, evalTime,fwdStartTime, 
															 fwdEndTime, fwdPeriod, fwdResetTime, payTime, states );
	
	// compute discount factor with num model
	ARM_VectorPtr df	= IRMODEL->GetDiscountFunctor()->DiscountFactor( curveName, evalTime, payTime, states );

	// compute option maturity
	double maturity	;
	if (fwdResetTime > evalTime)
		maturity = (fwdResetTime - evalTime) / K_YEAR_LEN;
	else if (fwdResetTime == evalTime)
		maturity = 1./K_YEAR_LEN;
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaSpreadOptionLet: reset date is prior eval date !" );
	
	
	// get volatility (same vol in all states)
	// vols are supposed to be stored in the first surface of the list
	double volatility	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(CAPLET_VOL, evalTime, fwdResetTime);
	
	// compute caplet/floorlet for all states
    int i, nbStates = df->size();
    ARM_VectorPtr values (new ARM_GP_Vector(nbStates));

	if (!LiborAsShiftedLognormal())
		for (i=0; i<nbStates; i++)
			(*values)[i]  =		payNotional 
							 *	period 
							 *  (*df)[i] 
							 *  VanillaOption_N  (	(*liborRate)[i], 
													volatility,
													strikesPerState[i], 
													maturity,  
													capFloor );
	else
		for (i=0; i<nbStates; i++)
			(*values)[i]  =		payNotional 
							 *	period 
							 *  (*df)[i] 
							 *  BS(	(*liborRate)[i]     + 1./fwdPeriod, 
									 strikesPerState[i] + 1./fwdPeriod,
									 maturity,
									 volatility,
									 capFloor );
		
		RESTORE_MODEL_FUNCTORS();
				
    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: VanillaSwaption
///	Returns: ARM_VectorPtr
///	Action : computes the price of a swaption
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::VanillaSwaption(
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
	
	ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaSwaption not implemented for Local_Normal_Model" );
	return ARM_VectorPtr( new ARM_GP_Vector (1, 0.0) );
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routines: Spread
///	Returns : a vector of spread values
///	Action  : Default Spread computation
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::Spread(
		const string& curveName, 
		double evalTime,
		double coeff1,
		double floatStartTime1, 
		double floatEndTime1, 
		const ARM_GP_Vector& fixPayTimes1,
		const ARM_GP_Vector& fixPayPeriods1,
		const ARM_GP_Vector& fwdStartTimes1,
        const ARM_GP_Vector& fwdEndTimes1,
        const ARM_GP_Vector& fwdPayPeriods1,
		const ARM_GP_Vector& floatPayTimes1,
        const ARM_GP_Vector& floatPayPeriods1,
        const ARM_GP_Vector& margin1,
		double coeff2,
		double floatStartTime2, 
		double floatEndTime2, 
		const ARM_GP_Vector& fixPayTimes2,
		const ARM_GP_Vector& fixPayPeriods2,
		const ARM_GP_Vector& fwdStartTimes2,
        const ARM_GP_Vector& fwdEndTimes2,
        const ARM_GP_Vector& fwdPayPeriods2,
		const ARM_GP_Vector& floatPayTimes2,
        const ARM_GP_Vector& floatPayPeriods2,
        const ARM_GP_Vector& margin2,
		const ARM_PricingStatesPtr& states) const
{
	if (GetFunctionalFlag())
	{
		CHANGE_MODEL_FUNCTORS();
		CAST_TO_IRMODEL();
		
		ARM_VectorPtr stdSpread = IRMODEL->Spread(	curveName,
											evalTime,
											coeff1,
											floatStartTime1, 
											floatEndTime1, 
											fixPayTimes1,
											fixPayPeriods1,
											fwdStartTimes1,
											fwdEndTimes1,
											fwdPayPeriods1,
											floatPayTimes1,
											floatPayPeriods1,
											margin1,
											coeff2,
											floatStartTime2, 
											floatEndTime2, 
											fixPayTimes2,
											fixPayPeriods2,
											fwdStartTimes2,
											fwdEndTimes2,
											fwdPayPeriods2,
											floatPayTimes2,
											floatPayPeriods2,
											margin2,
											states);

		ARM_VectorPtr newSpread = Func(evalTime,stdSpread);
		
		RESTORE_MODEL_FUNCTORS();
		return newSpread;
	}
	else
		ARM_THROW(  ERR_INVALID_ARGUMENT, "spread for normal model not ready yet!!" );

}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: VanillaSpreadOptionLet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a spread option
////////////////////////////////////////////////////
ARM_VectorPtr  ARM_Local_Normal_Model::VanillaSpreadOptionLet(
							const string& curveName,
							double evalTime,
							int callPut,
							double startTime, // NOT USED !
							double endTime,	  // NOT USED !
							double resetTime,
							double payTime,
							double payPeriod,
							double notional,
							double coeffLong,
							double coeffShort,
							const ARM_GP_Vector& strikes,
							double swapLongFloatStartTime,
							double swapLongFloatEndTime,
							const ARM_GP_Vector& swapLongFixPayTimes,
							const ARM_GP_Vector& swapLongFixPayPeriods,
							double swapShortFloatStartTime,
							double swapShortFloatEndTime,
							const ARM_GP_Vector& swapShortFixPayTimes,
							const ARM_GP_Vector& swapShortFixPayPeriods,
							const ARM_PricingStatesPtr& states)  const
{

	int i, nbStates = states->size();
	ARM_PricingFunctionIR* irFunctor = dynamic_cast < ARM_PricingFunctionIR* > (GetNumericalModel());
	if( !irFunctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "IR functions are not supported by numerical model" );


	/// compute caplet/floorlet for all states
    ARM_VectorPtr values (new ARM_GP_Vector(nbStates));

	CHANGE_MODEL_FUNCTORS();
	CAST_TO_IRMODEL();
	// compute discount factor with num model 
	ARM_VectorPtr df	= IRMODEL->GetDiscountFunctor()->DiscountFactor( curveName, evalTime, payTime, states );
	RESTORE_MODEL_FUNCTORS();

	/// particular case of 0 leverages
	/// if both leverages = 0, the SO are not calibrated in CSO calculator
	/// --> used fwd vols are obtained by extrapolation and are not, as required, equal to 0
	if (	fabs(coeffLong)<K_NEW_DOUBLE_TOL
		&&	fabs(coeffShort)<K_NEW_DOUBLE_TOL)
	{
		   for (i=0; i<nbStates; i++)
				(*values)[i]  =		notional 
								 *	payPeriod 
								 *  (*df)[i] 
								 * ( (callPut*strikes[i]<0) ? fabs(strikes[i]) : 0.0 );

	}
	else{

		// compute long and short swap rates for all states
		ARM_VectorPtr cmsRateLong (new ARM_GP_Vector(nbStates)), cmsRateShort (new ARM_GP_Vector(nbStates));;
		ARM_VectorPtr dfStart, dfEnd, level;
		
		dfStart = GetNumericalModel()->GetFixingFunctor()->DiscountFactor(curveName, evalTime, swapLongFloatStartTime, states);
		dfEnd   = GetNumericalModel()->GetFixingFunctor()->DiscountFactor(curveName, evalTime, swapLongFloatEndTime, states);
		level   = irFunctor->Annuity(curveName, evalTime, swapLongFixPayTimes, swapLongFixPayPeriods, states);
		for (i=0; i<nbStates; i++)
			(*cmsRateLong)[i] = ( (*dfStart)[i] - (*dfEnd)[i] ) / (*level)[i];


		dfStart = GetNumericalModel()->GetFixingFunctor()->DiscountFactor(curveName, evalTime, swapShortFloatStartTime, states);
		dfEnd   = GetNumericalModel()->GetFixingFunctor()->DiscountFactor(curveName, evalTime, swapShortFloatEndTime, states);
		level   = irFunctor->Annuity(curveName, evalTime, swapShortFixPayTimes, swapShortFixPayPeriods, states);
		for (i=0; i<nbStates; i++)
			(*cmsRateShort)[i] = ( (*dfStart)[i] - (*dfEnd)[i] ) / (*level)[i];

		
		// decl
		double adjustment;

		// maybe we should check here that the ForwardAdjustement surface list has 2 elements

		// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
		// adjustment for cms #1 supposed to be stored in first surface of the list
		adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SO_ADJ_LONG, evalTime, resetTime);
		(*cmsRateLong) += adjustment;

		// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
		// adjustment for cms #2 supposed to be stored in first surface of the list
		adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SO_ADJ_SHORT, evalTime, resetTime);
		(*cmsRateShort) += adjustment;

		// spread computation
		ARM_VectorPtr spread ( new ARM_GP_Vector(nbStates) );
		//(*spread) = coeffLong * (*cmsRateLong) - coeffShort * (*cmsRateShort);

		for (i=0; i<nbStates; i++)
			(*spread)[i] = coeffLong * (*cmsRateLong)[i] - coeffShort * (*cmsRateShort)[i] ;
		
		// compute maturity
		double maturity	;
		if (resetTime > evalTime)
			maturity = (resetTime - evalTime) / K_YEAR_LEN;
		else if (resetTime == evalTime)
			maturity = 1./K_YEAR_LEN;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaSpreadOptionLet: reset date is prior eval date !" );
		
		// get spread volatility (same vol in all states)
		// vols are supposed to be stored in the first surface of the list
		double volatility	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SO_VOL, evalTime, resetTime);
		
		
		for (i=0; i<nbStates; i++)
			(*values)[i]  =		notional 
							 *	payPeriod 
							 *  (*df)[i] 
							 *  VanillaOption_N  (	(*spread)[i], 
													volatility,
													strikes[i],
													maturity,  
													callPut );
	}



    return values;
}


////////////////////////////////////////////////////
///	Routine: ComputeSwapRates
///	Returns: ARM_VectorPtr
///	Action : Compute the swap rate with ARM_SwapRate
/// structure
////////////////////////////////////////////////////
ARM_VectorPtr ComputeSwapRates(	ARM_PricingModel* numericalModel,
								const string& curveName,
								double evalTime,
								double payTime,
								const ARM_SwapRate& swapRate,
								const   ARM_PricingStatesPtr& states)
{
	int i, nbStates = states->size();

	ARM_VectorPtr cmsRate (new ARM_GP_Vector(nbStates));

	ARM_PricingFunctionIR* irFunctor = dynamic_cast < ARM_PricingFunctionIR* > (numericalModel);

	ARM_VectorPtr dfStart, dfEnd, level;

	dfStart = numericalModel->DiscountFactor(curveName, evalTime, swapRate.floatStartTime, states);
	dfEnd   = numericalModel->DiscountFactor(curveName, evalTime, swapRate.floatEndTime, states);
	level   = irFunctor->Annuity(curveName, evalTime, swapRate.fixPayTimes, swapRate.fixPayPeriods, states);
		for (i=0; i<nbStates; i++)
			(*cmsRate)[i] = ( (*dfStart)[i] - (*dfEnd)[i] ) / (*level)[i];

	return cmsRate;
}

void ARM_Local_Normal_Model::ComputeDoubleCorridorDigital(
				ARM_VectorPtr& spreads,
				ARM_VectorPtr& rates,
				ARM_GP_Vector& spreadAdj,
				ARM_GP_Vector& rateAdj,
				ARM_VectorPtr& payedRates,
				double spreadDownBarrier,
				double spreadUpBarrier,
				double rateDownBarrier,
				double rateUpBarrier,
				double evalTime,
				double resetTime,
				double sqrTtoE,
				double coefCoupon,
				size_t KAdjSpreadSdRd,
				size_t KAdjSpreadSdRu,
				size_t KAdjSpreadSuRd,
				size_t KAdjSpreadSuRu,
				size_t KAdjRateSdRd,
				size_t KAdjRateSdRu,
				size_t KAdjRateSuRd,
				size_t KAdjRateSuRu,
				size_t VolSpreadSdRd,
				size_t VolSpreadSdRu,
				size_t VolSpreadSuRd,
				size_t VolSpreadSuRu,
				size_t VolRateSdRd,
				size_t VolRateSdRu,
				size_t VolRateSuRd,
				size_t VolRateSuRu,
				size_t CorrelSpreadRateSdRd,
				size_t CorrelSpreadRateSdRu,
				size_t CorrelSpreadRateSuRd,
				size_t CorrelSpreadRateSuRu,
				ARM_VectorPtr& values) const
{
	/// Constants for NC(x,y,rho) integration
	double nbGL1=5,nbGL2=3,maxX=5.0;

	size_t nbStates = spreads->size();

	double actBarLimit = 0.01 * ARM_CorridorDblCondition::BarrierLimit; // no more in %

	bool isRateBarDownOff	= rateDownBarrier < -actBarLimit;
	bool isRateBarUpOff		= rateUpBarrier > actBarLimit;
	bool isSpreadBarDownOff	= spreadDownBarrier < -actBarLimit;
	bool isSpreadBarUpOff	= spreadUpBarrier > actBarLimit;

	/// Get adjusted strikes for each spread/rate barriers
	double spreadSdRdKAdj = spreadDownBarrier	+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjSpreadSdRd, evalTime, resetTime);
	double spreadSdRuKAdj = spreadDownBarrier	+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjSpreadSdRu, evalTime, resetTime);
	double spreadSuRdKAdj = spreadUpBarrier		+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjSpreadSuRd, evalTime, resetTime);
	double spreadSuRuKAdj = spreadUpBarrier		+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjSpreadSuRu, evalTime, resetTime);

	double rateSdRdKAdj = rateDownBarrier	+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjRateSdRd, evalTime, resetTime);
	double rateSdRuKAdj = rateUpBarrier		+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjRateSdRu, evalTime, resetTime);
	double rateSuRdKAdj = rateDownBarrier	+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjRateSuRd, evalTime, resetTime);
	double rateSuRuKAdj = rateUpBarrier		+	GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KAdjRateSuRu, evalTime, resetTime);

	/// Get adjusted volatilities for each spread/rate barriers
	double spreadSdRdVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolSpreadSdRd, evalTime, resetTime);
	double spreadSdRuVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolSpreadSdRu, evalTime, resetTime);
	double spreadSuRdVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolSpreadSuRd, evalTime, resetTime);
	double spreadSuRuVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolSpreadSuRu, evalTime, resetTime);

	double rateSdRdVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolRateSdRd, evalTime, resetTime);
	double rateSdRuVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolRateSdRu, evalTime, resetTime);
	double rateSuRdVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolRateSuRd, evalTime, resetTime);
	double rateSuRuVolAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VolRateSuRu, evalTime, resetTime);

	/// Get adjusted spread/rate correlations for each spread/rate barriers
	double spreadRateSdRdCorrelAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(CorrelSpreadRateSdRd, evalTime, resetTime);
	double spreadRateSdRuCorrelAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(CorrelSpreadRateSdRu, evalTime, resetTime);
	double spreadRateSuRdCorrelAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(CorrelSpreadRateSuRd, evalTime, resetTime);
	double spreadRateSuRuCorrelAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(CorrelSpreadRateSuRu, evalTime, resetTime);

	double invStdDevMinValue = 1.0e+10;
	double invSpreadSdRdStdDev	= (spreadSdRdVolAdj != 0.0 ? 1/(spreadSdRdVolAdj*sqrTtoE) : invStdDevMinValue);
	double invSpreadSdRuStdDev	= (spreadSdRuVolAdj != 0.0 ? 1/(spreadSdRuVolAdj*sqrTtoE) : invStdDevMinValue);
	double invSpreadSuRdStdDev	= (spreadSuRdVolAdj != 0.0 ? 1/(spreadSuRdVolAdj*sqrTtoE) : invStdDevMinValue);
	double invSpreadSuRuStdDev	= (spreadSuRuVolAdj != 0.0 ? 1/(spreadSuRuVolAdj*sqrTtoE) : invStdDevMinValue);

	double invRateSdRdStdDev	= (rateSdRdVolAdj != 0.0 ? 1/(rateSdRdVolAdj*sqrTtoE) : invStdDevMinValue);
	double invRateSdRuStdDev	= (rateSdRuVolAdj != 0.0 ? 1/(rateSdRuVolAdj*sqrTtoE) : invStdDevMinValue);
	double invRateSuRdStdDev	= (rateSuRdVolAdj != 0.0 ? 1/(rateSuRdVolAdj*sqrTtoE) : invStdDevMinValue);
	double invRateSuRuStdDev	= (rateSuRuVolAdj != 0.0 ? 1/(rateSuRuVolAdj*sqrTtoE) : invStdDevMinValue);
	double priceSdRd,priceSdRu,priceSuRd,priceSuRu;

	bool isFloatPay = payedRates != ARM_VectorPtr(NULL);
	bool isAdjArray = spreadAdj.size()!=1;

	double spread,rate,factor,adjS,adjR;
	if((isRateBarDownOff || isRateBarUpOff) && (isSpreadBarDownOff || isSpreadBarUpOff))
	{
		/// Only one double condition is active
		if(isRateBarDownOff)
		{
			if(isSpreadBarDownOff)
			{
				/// S<=Us & R<=Ur => NC(Us,Ur)
				adjS = spreadAdj[isAdjArray ? DBLECOR_SUP_RUP : 0];
				adjR = rateAdj[isAdjArray ? DBLECOR_SUP_RUP : 0];
				for (size_t j=0; j<nbStates; j++)
				{
					spread	= (*spreads)[j] + adjS;
					rate	= (*rates)[j] + adjR;
					factor	= isFloatPay ? coefCoupon * (*payedRates)[j] : coefCoupon;
					(*values)[j] += factor * NormalCDF((spreadSuRuKAdj - spread)*invSpreadSuRuStdDev, (rateSuRuKAdj - rate)*invRateSuRuStdDev, spreadRateSuRuCorrelAdj,nbGL1,nbGL2,-maxX,maxX);

				}
			}
			else
			{
				/// Ls<=S & R<=Ur => NC(Ur) - NC(Ls,Ur)
				adjS = spreadAdj[isAdjArray ? DBLECOR_SDOWN_RUP : 0];
				adjR = rateAdj[isAdjArray ? DBLECOR_SDOWN_RUP : 0];
				for (size_t j=0; j<nbStates; j++)
				{
					spread	= (*spreads)[j] + adjS;
					rate	= (*rates)[j] + adjR;
					factor	= isFloatPay ? coefCoupon * (*payedRates)[j] : coefCoupon;
					(*values)[j] += factor * (NormalCDF((rateSdRuKAdj - rate)*invRateSdRuStdDev)
												  - NormalCDF((spreadSdRuKAdj - spread)*invSpreadSdRuStdDev, (rateSdRuKAdj - rate)*invRateSdRuStdDev, spreadRateSdRuCorrelAdj,nbGL1,nbGL2,-maxX,maxX));
				}
			}
		}
		else
		{
			if(isSpreadBarDownOff)
			{
				/// S<=Us & Lr<=R => NC(Us) - NC(Us,Lr)
				adjS = spreadAdj[isAdjArray ? DBLECOR_SUP_RDOWN : 0];
				adjR = rateAdj[isAdjArray ? DBLECOR_SUP_RDOWN : 0];
				for (size_t j=0; j<nbStates; j++)
				{
					spread	= (*spreads)[j] + adjS;
					rate	= (*rates)[j] + adjR;
					factor	= isFloatPay ? coefCoupon * (*payedRates)[j] : coefCoupon;
					(*values)[j] += factor * (NormalCDF((spreadSuRdKAdj - spread)*invSpreadSuRdStdDev)
												  - NormalCDF((spreadSuRdKAdj - spread)*invSpreadSuRdStdDev, (rateSuRdKAdj - rate)*invRateSuRdStdDev, spreadRateSuRdCorrelAdj,nbGL1,nbGL2,-maxX,maxX));
				}
			}
			else
			{
				/// Ls<=S & Lr<=R => 1 - NC(Ls) - NC(Lr) + NC(Ls,Lr)
				adjS = spreadAdj[isAdjArray ? DBLECOR_SDOWN_RDOWN : 0];
				adjR = rateAdj[isAdjArray ? DBLECOR_SDOWN_RDOWN : 0];
				for (size_t j=0; j<nbStates; j++)
				{
					spread	= (*spreads)[j] + adjS;
					rate	= (*rates)[j] + adjR;
					factor	= isFloatPay ? coefCoupon * (*payedRates)[j] : coefCoupon;
					(*values)[j] += factor * (1
												  - NormalCDF((rateSdRdKAdj - rate)*invRateSdRdStdDev)
												  - NormalCDF((spreadSdRdKAdj - spread)*invSpreadSdRdStdDev)
												  + NormalCDF((spreadSdRdKAdj - spread)*invSpreadSdRdStdDev, (rateSdRdKAdj - rate)*invRateSdRdStdDev, spreadRateSdRdCorrelAdj,nbGL1,nbGL2,-maxX,maxX));
				}
			}
		}
	}
	else
	{
		/// Full double condition ranges
		double adjSSdRd = spreadAdj[isAdjArray ? DBLECOR_SDOWN_RDOWN : 0];
		double adjSSdRu = spreadAdj[isAdjArray ? DBLECOR_SDOWN_RUP : 0];
		double adjSSuRd = spreadAdj[isAdjArray ? DBLECOR_SUP_RDOWN : 0];
		double adjSSuRu = spreadAdj[isAdjArray ? DBLECOR_SUP_RUP : 0];

		double adjRSdRd = rateAdj[isAdjArray ? DBLECOR_SDOWN_RDOWN : 0];
		double adjRSdRu = rateAdj[isAdjArray ? DBLECOR_SDOWN_RUP : 0];
		double adjRSuRd = rateAdj[isAdjArray ? DBLECOR_SUP_RDOWN : 0];
		double adjRSuRu = rateAdj[isAdjArray ? DBLECOR_SUP_RUP : 0];
		for (size_t j=0; j<nbStates; j++)
		{
			factor	= isFloatPay ? coefCoupon * (*payedRates)[j] : coefCoupon;
			priceSdRd = NormalCDF( (spreadSdRdKAdj - ((*spreads)[j]+adjSSdRd))*invSpreadSdRdStdDev, (rateSdRdKAdj - ((*rates)[j]+adjRSdRd))*invRateSdRdStdDev, spreadRateSdRdCorrelAdj,nbGL1,nbGL2,-maxX,maxX);
			priceSdRu = NormalCDF( (spreadSdRuKAdj - ((*spreads)[j]+adjSSdRu))*invSpreadSdRuStdDev, (rateSdRuKAdj - ((*rates)[j]+adjRSdRu))*invRateSdRuStdDev, spreadRateSdRuCorrelAdj,nbGL1,nbGL2,-maxX,maxX);
			priceSuRd = NormalCDF( (spreadSuRdKAdj - ((*spreads)[j]+adjSSuRd))*invSpreadSuRdStdDev, (rateSuRdKAdj - ((*rates)[j]+adjRSuRd))*invRateSuRdStdDev, spreadRateSuRdCorrelAdj,nbGL1,nbGL2,-maxX,maxX);
			priceSuRu = NormalCDF( (spreadSuRuKAdj - ((*spreads)[j]+adjSSuRu))*invSpreadSuRuStdDev, (rateSuRuKAdj - ((*rates)[j]+adjRSuRu))*invRateSuRuStdDev, spreadRateSuRuCorrelAdj,nbGL1,nbGL2,-maxX,maxX);

			(*values)[j] += factor * (priceSuRu - priceSuRd - priceSdRu + priceSdRd);
		}
	}
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: VanillaCMSCorridorlet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a CMS corridorlet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::VanillaCMSCorridorlet(
	const string& curveName,
	double evalTime,
	double payTime,
	double resetTime,
	double startTime,
	double endTime,
	const ARM_GP_Vector& refIdxResettimes,
	const ARM_GP_Vector& refIndexWeights,
	const ARM_GP_Vector& coeff1,
	const ARM_SwapRatePtrVector& firstIndex,
	const ARM_GP_Vector& coeff2,
	const ARM_SwapRatePtrVector& secondIndex,
	int		payIndexType,			/// K_FIXED, K_LIBOR or K_CMS
	double	coupon,					/// fixed rate (K_FIXED) or spread (K_LIBOR or K_CMS)
	const	ARM_SwapRate& payRate,	/// rate description (K_LIBOR or K_CMS)
	double  payIndexLeverage,
	const ARM_GP_Vector& spreadDownBarriers,
    const ARM_GP_Vector& spreadUpBarriers,
    double  payNotional,
    int     rcvPay,
	const ARM_SwapRatePtrVector& thirdIndex, // Single rate index for double condition
	const ARM_GP_Vector& rateDownBarriers,
	const ARM_GP_Vector& rateUpBarriers,
    const   ARM_PricingStatesPtr& states) const
{
/****
if(lastEval > evalTime)
{
	FILE* f=fopen("c:\\temp\\dumpCRA2Timer.txt","a");
	fprintf(f,"Eval=%6.1lf\tDuration(ms) =\t%10.5lf\t%10.5lf\tNbEvals=%3d\n",evalTime,computingTime,computingTimeAll,nbEvals);
	fclose(f);
	computingTime=0.0;
	computingTimeAll=0.0;
	nbEvals=0;
}
ARM_Timer timerAll;
timerAll.ClockStartTime();
****/
	CHANGE_MODEL_FUNCTORS();
	CAST_TO_IRMODEL();
	// downcast to get functor
	// performance issue : could be done in SetNumericalModel ....
	ARM_PricingFunctionIR* irFunctor = dynamic_cast < ARM_PricingFunctionIR* > (GetNumericalModel());

	if( !irFunctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "IR functions are not supported by numerical model" );
	
	int i, j, nbStates = states->size(), nbFixings = refIdxResettimes.size();

	bool isDbleCorridor = (thirdIndex.size() > 0);
	bool isPVCalib=false;
	ARM_GP_VectorPtr spotProbas(NULL);

	if(	GetPVAdjuster() &&
		GetNumericalModel()->GetNumMethod() != ARM_NumMethodPtr(NULL) &&
		GetNumericalModel()->GetNumeraire() != ARM_NumerairePtr(NULL) &&
		dynamic_cast<const ARM_NumeraireFwd *const>(&*(GetNumericalModel())->GetNumeraire()))
	{
		/// Arrow-Debreu computation can't be forced (the model only decides) then
		/// spot probabilities may be used but only for forward like numeraire
		int timeIdx	= GetNumericalModel()->GetNumMethod()->GetLastTimeIdx();
		spotProbas	= GetNumericalModel()->GetNumMethod()->GetSpotProbabilities(timeIdx);
		isPVCalib	= (spotProbas != ARM_GP_VectorPtr(NULL));
	}


	double longAdjustment, shortAdjustment, payAdjustment;
	double longAdjustmentDown, shortAdjustmentDown;
	double longAdjustmentUp,   shortAdjustmentUp;
	double volDownLeft,volDownRight,volUpLeft,volUpRight;
	double kAdjDownLeft,kAdjDownRight,kAdjUpLeft,kAdjUpRight;

	// compute long and short swap rates for all the fixing and for all states
	ARM_VectorPtr cmsRateLong,cmsRateShort,rates;
	ARM_VectorPtr spreads ( new ARM_GP_Vector(nbStates) );

	ARM_GP_Vector spreadAdjFlt(DBLECOR_NB_BARRIERS),rateAdjFlt(DBLECOR_NB_BARRIERS);
	ARM_GP_Vector spreadAdj(1),rateAdj(1);

	// compute caplet/floorlet for all states
	ARM_VectorPtr values (new ARM_GP_Vector(nbStates,0.0));

	// compute discount factor with num model 
	ARM_VectorPtr df	= GetNumericalModel()->GetDiscountFunctor()->DiscountFactor( curveName, evalTime, payTime, states );


	// We assume that the cap spread is computed with a gap of 1bp 
	const double strikeSpread = SOCOR_STRIKE_SPREAD;

	/// ----------------------------
	/// case of a floating payment
	/// ----------------------------
	ARM_VectorPtr payedRates(NULL);
	bool isVariableCorridor = (payIndexType == K_LIBOR || payIndexType == K_CMS);
	if (isVariableCorridor)
	{
		payedRates = ComputeSwapRates(GetNumericalModel(),curveName,evalTime,payTime,payRate,states);
		// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
		// adjustment for payed rate supposed to be stored in fifth surface of the list
		payAdjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(isDbleCorridor ? DBLECOR_ADJ_PAY_FLT : SOCOR_ADJ_PAY_FLT, evalTime, resetTime);
		*payedRates += payAdjustment;

/***
FILE* f;
if(isDbleCorridor)
	f=fopen("c:\\temp\\dumpCRA2.txt","a");
else
	f=fopen("c:\\temp\\dumpCCSO.txt","a");
fprintf(f,"notice=%6d\tpay=%6d\tpayAdj=%10.7lf\n",(int)(evalTime),(int)(payTime),payAdjustment);
fclose(f);
***/

	}

	double maturity,sqrTtoE,coefCoupon;
	double resetWeight=0.0;
	size_t lastFixIdx=0;
	bool isSkippedFix=false;
	for (i = 0; i < nbFixings; ++i)
	{
		//isSkippedFix = i>0 && i+1<nbFixings && i-lastFixIdx < 5 &&
		//				((refIdxResettimes[i]-refIdxResettimes[i-1] < 5.5 && refIdxResettimes[i] - evalTime > COR_OPTIM_EVAL_DAILYLIMIT) ||
		//				 (refIdxResettimes[i]-refIdxResettimes[i-1] < 15.0 && refIdxResettimes[i] - evalTime > COR_OPTIM_EVAL_WEEKLYLIMIT));

		resetWeight += refIndexWeights[i];

		if(isSkippedFix)
		{
			/// This reset is not evaluated but its weight is added to next one
			continue;
		}
		else
		{
			lastFixIdx = i;
		}

		coefCoupon = resetWeight * coupon;

		cmsRateLong  = ComputeSwapRates(GetNumericalModel(),curveName,evalTime,payTime,*firstIndex[i],states);
		cmsRateShort = ComputeSwapRates(GetNumericalModel(),curveName,evalTime,payTime,*secondIndex[i],states);

		if(isDbleCorridor)
		{
			rates = ComputeSwapRates(GetNumericalModel(),curveName,evalTime,payTime,*thirdIndex[i],states);

			/// Adjustments are saved for first condition spreads and second condition rates
			spreadAdj[0]	= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_ADJ_SPREAD, evalTime, refIdxResettimes[i]);
			rateAdj[0]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_ADJ_RATE, evalTime, refIdxResettimes[i]);
			for (j=0; j<nbStates; j++)
				(*spreads)[j]	= coeff1[i] * (*cmsRateLong)[j] - coeff2[i] * (*cmsRateShort)[j];
		}
		else
		{
			// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
			// adjustment for cms #1 supposed to be stored in first surface of the list
			longAdjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_LONG, evalTime, refIdxResettimes[i]);
			
			// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
			// adjustment for cms #2 supposed to be stored in second surface of the list
			shortAdjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_SHORT, evalTime, refIdxResettimes[i]);
			
			for (j=0; j<nbStates; j++)
				(*spreads)[j] = coeff1[i] * ((*cmsRateLong)[j] + longAdjustment) - coeff2[i] * ((*cmsRateShort)[j] + shortAdjustment);
		}


		// compute maturity
		if (refIdxResettimes[i] > evalTime)
			maturity = (refIdxResettimes[i] - evalTime) / K_YEAR_LEN;
		else if (refIdxResettimes[i] == evalTime)
			maturity = 1./K_YEAR_LEN;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaCMSCorridorlet: reset date is prior eval date !" );
		sqrTtoE = sqrt(maturity);

	
		if(isDbleCorridor)
		{
			if(coefCoupon != 0.0)
			{
//ARM_Timer timer;
//timer.ClockStartTime();
				ComputeDoubleCorridorDigital(spreads,rates,spreadAdj,rateAdj,ARM_VectorPtr(NULL),
					spreadDownBarriers[i],spreadUpBarriers[i],rateDownBarriers[i],rateUpBarriers[i],
					evalTime,refIdxResettimes[i],sqrTtoE,coefCoupon,
					DBLECOR_SDOWN_RDOWN_KADJ_SPREAD,DBLECOR_SDOWN_RUP_KADJ_SPREAD,DBLECOR_SUP_RDOWN_KADJ_SPREAD,DBLECOR_SUP_RUP_KADJ_SPREAD,
					DBLECOR_SDOWN_RDOWN_KADJ_RATE,DBLECOR_SDOWN_RUP_KADJ_RATE,DBLECOR_SUP_RDOWN_KADJ_RATE,DBLECOR_SUP_RUP_KADJ_RATE,
					DBLECOR_SDOWN_RDOWN_VOL_SPREAD,DBLECOR_SDOWN_RUP_VOL_SPREAD,DBLECOR_SUP_RDOWN_VOL_SPREAD,DBLECOR_SUP_RUP_VOL_SPREAD,
					DBLECOR_SDOWN_RDOWN_VOL_RATE,DBLECOR_SDOWN_RUP_VOL_RATE,DBLECOR_SUP_RDOWN_VOL_RATE,DBLECOR_SUP_RUP_VOL_RATE,
					DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL,DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL,DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL,DBLECOR_SUP_RUP_SPREAD_RATE_CORREL,
					values);
/***
timer.ClockEndTime();
computingTime += timer.GetDuration()*1000.0;
++nbEvals;
lastEval=evalTime;
***/
			}
		}
		else
		{
//ARM_Timer timer;
//timer.ClockStartTime();
			// get the 2 spread volatilities for the 2 barriers : 4 volatilites  (same vol in all states)
			// vols are supposed to be stored in the 4 first surfaces of the list
			volDownLeft = 0.0, volDownRight = 0.0, volUpLeft = 0.0, volUpRight = 0.0;
					
			volUpLeft	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_UPLEFT, evalTime, refIdxResettimes[i]);
			volUpRight	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_UPRIGHT, evalTime, refIdxResettimes[i]);
			volDownLeft  = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_DOWNLEFT, evalTime, refIdxResettimes[i]);
			volDownRight = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_DOWNRIGHT, evalTime, refIdxResettimes[i]);

			kAdjUpLeft	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_UPLEFT, evalTime, refIdxResettimes[i]);
			kAdjUpRight	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_UPRIGHT, evalTime, refIdxResettimes[i]);
			kAdjDownLeft  = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_DOWNLEFT, evalTime, refIdxResettimes[i]);
			kAdjDownRight = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_DOWNRIGHT, evalTime, refIdxResettimes[i]);

/*FILE* f=fopen("c:\\temp\\dumpCCSO.txt","a");
fprintf(f,"reset=%6d\tAUL=%10.7lf\tAUR=%10.7lf\tADL=%10.7lf\tADR=%10.7lf\n",
		(int)(refIdxResettimes[i]),
		kAdjUpLeft,
		kAdjUpRight,
		kAdjDownLeft,
		kAdjDownRight);
fclose(f);*/

			for (j=0; j<nbStates; j++)
			{
				(*values)[j]  +=	coefCoupon *
									(VanillaOption_N  (	(*spreads)[j],
														volDownLeft,
														spreadDownBarriers[i]-strikeSpread+kAdjDownLeft,
														maturity, 
														K_CAP )-
									VanillaOption_N  (	(*spreads)[j],
														volDownRight,
														spreadDownBarriers[i]+strikeSpread+kAdjDownRight,
														maturity, 
														K_CAP ))/2/strikeSpread;


				(*values)[j]  -=	coefCoupon *
									((VanillaOption_N  (	(*spreads)[j],
														volUpLeft,
														spreadUpBarriers[i]-strikeSpread+kAdjUpLeft,
														maturity, 
														K_CAP )-
									VanillaOption_N  (	(*spreads)[j],
														volUpRight,
														spreadUpBarriers[i]+strikeSpread+kAdjUpRight,
														maturity, 
														K_CAP ))/2/strikeSpread);
			}
/***
timer.ClockEndTime();
computingTime += timer.GetDuration()*1000.0;
++nbEvals;
lastEval=evalTime;
***/
		}

		/// ----------------------------
		/// case of a floating payment
		/// ----------------------------
		if (isVariableCorridor)
		{
			coefCoupon = resetWeight * payIndexLeverage;

			if(isDbleCorridor)
			{
				if(coefCoupon != 0.0)
				{
					/// Adjustments are saved for each barriers and first condition spreads and second condition rates
					spreadAdjFlt[DBLECOR_SDOWN_RDOWN]	= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SDOWN_RDOWN_ADJ_SPREAD_FLT, evalTime, refIdxResettimes[i]);
					spreadAdjFlt[DBLECOR_SDOWN_RUP]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SDOWN_RUP_ADJ_SPREAD_FLT, evalTime, refIdxResettimes[i]);
					spreadAdjFlt[DBLECOR_SUP_RDOWN]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SUP_RDOWN_ADJ_SPREAD_FLT, evalTime, refIdxResettimes[i]);
					spreadAdjFlt[DBLECOR_SUP_RUP]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SUP_RUP_ADJ_SPREAD_FLT, evalTime, refIdxResettimes[i]);

					rateAdjFlt[DBLECOR_SDOWN_RDOWN]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SDOWN_RDOWN_ADJ_RATE_FLT, evalTime, refIdxResettimes[i]);
					rateAdjFlt[DBLECOR_SDOWN_RUP]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SDOWN_RUP_ADJ_RATE_FLT, evalTime, refIdxResettimes[i]);
					rateAdjFlt[DBLECOR_SUP_RDOWN]		= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SUP_RDOWN_ADJ_RATE_FLT, evalTime, refIdxResettimes[i]);
					rateAdjFlt[DBLECOR_SUP_RUP]			= GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DBLECOR_SUP_RUP_ADJ_RATE_FLT, evalTime, refIdxResettimes[i]);

					ComputeDoubleCorridorDigital(spreads,rates,spreadAdjFlt,rateAdjFlt,payedRates,
						spreadDownBarriers[i],spreadUpBarriers[i],rateDownBarriers[i],rateUpBarriers[i],
						evalTime,refIdxResettimes[i],sqrTtoE,coefCoupon,
						DBLECOR_SDOWN_RDOWN_KADJ_SPREAD_FLT,DBLECOR_SDOWN_RUP_KADJ_SPREAD_FLT,DBLECOR_SUP_RDOWN_KADJ_SPREAD_FLT,DBLECOR_SUP_RUP_KADJ_SPREAD_FLT,
						DBLECOR_SDOWN_RDOWN_KADJ_RATE_FLT,DBLECOR_SDOWN_RUP_KADJ_RATE_FLT,DBLECOR_SUP_RDOWN_KADJ_RATE_FLT,DBLECOR_SUP_RUP_KADJ_RATE_FLT,
						DBLECOR_SDOWN_RDOWN_VOL_SPREAD_FLT,DBLECOR_SDOWN_RUP_VOL_SPREAD_FLT,DBLECOR_SUP_RDOWN_VOL_SPREAD_FLT,DBLECOR_SUP_RUP_VOL_SPREAD_FLT,
						DBLECOR_SDOWN_RDOWN_VOL_RATE_FLT,DBLECOR_SDOWN_RUP_VOL_RATE_FLT,DBLECOR_SUP_RDOWN_VOL_RATE_FLT,DBLECOR_SUP_RUP_VOL_RATE_FLT,
						DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL_FLT,DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL_FLT,DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL_FLT,DBLECOR_SUP_RUP_SPREAD_RATE_CORREL_FLT,
						values);
				}

/***
FILE* f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"reset=%6d\tcoefCpn=%10.7lf\tSdRd=%10.7lf\tSdRu=%10.7lf\tSuRd=%10.7lf\tSuRu=%10.7lf\n",
		(int)(refIdxResettimes[i]),coefCoupon,spreadSdRdAdj,spreadSdRuAdj,spreadSuRdAdj,spreadSuRuAdj);
fclose(f);
***/

			}
			else
			{
				// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
				// adjustment for cms #1 supposed to be stored in third surface of the list
				longAdjustmentDown = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_LONG_DOWN_FLT, evalTime, refIdxResettimes[i]);
				longAdjustmentUp   = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_LONG_UP_FLT,   evalTime, refIdxResettimes[i]);
				
				// adjust forward swap rate to get cms rate (additive, same adjustment in all states)
				// adjustment for cms #2 supposed to be stored in fourth surface of the list
				shortAdjustmentDown = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_SHORT_DOWN_FLT, evalTime, refIdxResettimes[i]);
				shortAdjustmentUp   = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(SOCOR_ADJ_SHORT_UP_FLT,   evalTime, refIdxResettimes[i]);

/***
FILE* f=fopen("c:\\temp\\dumpCCSO.txt","a");
fprintf(f,"reset=%6d\tcoefCpn=%10.7lf\tSd=%10.7lf\tSu=%10.7lf\n",
		(int)(refIdxResettimes[i]),coefCoupon,
		coeff1[i]*longAdjustmentDown-coeff2[i]*shortAdjustmentDown,
		coeff1[i]*longAdjustmentUp-coeff2[i]*shortAdjustmentUp);
fclose(f);
***/

				/// spread (adjusted down
				for (j=0; j<nbStates; j++)
					(*spreads)[j] = coeff1[i] * ((*cmsRateLong)[j] + longAdjustmentDown) - coeff2[i] * ((*cmsRateShort)[j] + shortAdjustmentDown);
					
				
				// vols are supposed to be stored in the 4 following surfaces of the list
				volUpLeft	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_UPLEFT_FLT, evalTime, refIdxResettimes[i]);
				volUpRight	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_UPRIGHT_FLT, evalTime, refIdxResettimes[i]);
				volDownLeft  = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_DOWNLEFT_FLT, evalTime, refIdxResettimes[i]);
				volDownRight = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(SOCOR_VOL_DOWNRIGHT_FLT, evalTime, refIdxResettimes[i]);

				kAdjUpLeft	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_UPLEFT_FLT, evalTime, refIdxResettimes[i]);
				kAdjUpRight	 = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_UPRIGHT_FLT, evalTime, refIdxResettimes[i]);
				kAdjDownLeft  = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_DOWNLEFT_FLT, evalTime, refIdxResettimes[i]);
				kAdjDownRight = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_KADJ_DOWNRIGHT_FLT, evalTime, refIdxResettimes[i]);

/*FILE* f=fopen("c:\\temp\\dumpCCSO.txt","a");
fprintf(f,"reset=%6d\tAUL=%10.7lf\tAUR=%10.7lf\tADL=%10.7lf\tADR=%10.7lf\n",
		(int)(refIdxResettimes[i]),
		kAdjUpLeft,
		kAdjUpRight,
		kAdjDownLeft,
		kAdjDownRight);
fclose(f);*/

				//coefCoupon /= (2*strikeSpread);
				for (j=0; j<nbStates; j++)
				{
					(*values)[j]  +=	coefCoupon * (*payedRates)[j] *
										(VanillaOption_N  (	(*spreads)[j],
															volDownLeft,
															spreadDownBarriers[i]-strikeSpread+kAdjDownLeft,
															maturity, 
															K_CAP )-
										VanillaOption_N  (	(*spreads)[j],
															volDownRight,
															spreadDownBarriers[i]+strikeSpread+kAdjDownRight,
															maturity, 
															K_CAP ))/2/strikeSpread;



				}

				/// spread (adjusted up)
				for (j=0; j<nbStates; j++)
					(*spreads)[j]   = coeff1[i] * ((*cmsRateLong)[j] + longAdjustmentUp)   - coeff2[i] * ((*cmsRateShort)[j] + shortAdjustmentUp);

				for (j=0; j<nbStates; j++)
				{
					(*values)[j]  -=	coefCoupon * (*payedRates)[j] *
										((VanillaOption_N  (	(*spreads)[j],
															volUpLeft,
															spreadUpBarriers[i]-strikeSpread+kAdjUpLeft,
															maturity, 
															K_CAP )-
										VanillaOption_N  (	(*spreads)[j],
															volUpRight,
															spreadUpBarriers[i]+strikeSpread+kAdjUpRight,
															maturity, 
															K_CAP ))/2/strikeSpread);
				}
			}
		}

		if(!isSkippedFix)
			resetWeight = 0.0;
	}

	
	coefCoupon = payNotional *	rcvPay;
	for (i=0; i<nbStates; i++)
		(*values)[i]  *= (coefCoupon * (*df)[i]);


	if(isPVCalib)
	{
		size_t correctionParamIdx = (isDbleCorridor ? ARM_ModelParamType::Correlation : ARM_ModelParamType::Shift);
		size_t targetIdx = (isDbleCorridor ? DBLECOR_MKT_RESIDUAL_RA2_PRICE : SOCOR_MKT_RESIDUAL_RSO_PRICE);
		size_t correctionIdx = (isDbleCorridor ? DBLECOR_RESIDUAL_RA2_PV_CORRECTION : SOCOR_RESIDUAL_RSO_PV_CORRECTION);

		/// Correct state values to exactly fit the residual RA2 through a relative shift
		double targetPrice = 0.0;
		for(i=0;i<nbFixings;++i)
			targetPrice += GetModelParams()->GetModelParam(correctionParamIdx).GetValue(targetIdx,evalTime,refIdxResettimes[i]);
		double modelPrice = 0.0;
		ARM_GP_Vector initialValues(*values);
		GetNumeraire()->ProcessPaidPayoffs(curveName,values,evalTime,states,*this);
		for(i=0;i<nbStates;++i)
			modelPrice += (*values)[i] * (*spotProbas)[i];
		double correction = fabs(modelPrice) > ARM_NumericConstants::ARM_TOLERENCE ? targetPrice/modelPrice : 1.0;

		const_cast<ARM_Local_Normal_Model*>(this)->GetModelParams()->GetModelParam(correctionParamIdx).SetValue(correctionIdx,evalTime,resetTime,correction);

		for(i=0;i<nbStates;++i)
			(*values)[i] = initialValues[i] * correction;
	}
	
	RESTORE_MODEL_FUNCTORS();

//timerAll.ClockEndTime();
//computingTimeAll += timerAll.GetDuration()*1000.0;

    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_MultiAssetsModel
///	Routine: RangeAccrualVectorial
///	Returns: void (go to the local normal model)
///	Action : compute the corridor double condition: one condition
///				on FX, one condition on Libor, both fixing 
///				at the same date
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::RangeAccrualVectorial(
		const string& curveName,
		double evalTime,
		double startTime,
		double endTime,
		double payTime,
		const  ARM_GP_Vector& fixingTimes,
		int    payIndexType, 
        double payIndexTerm,
		const  string& fxModelName,
		int    irIndexType, 
		const  ARM_GP_Vector& irIndexResetTimes,
		const  ARM_GP_Vector& irIndexStartTimes,
		const  ARM_GP_Vector& irIndexEndTimes,
		const  ARM_GP_Vector& irIndexTerms,
		const  ARM_GP_Vector& fxDownBarriers,
		const  ARM_GP_Vector& fxUpBarriers,
		const  ARM_GP_Vector& irDownBarriers,
		const  ARM_GP_Vector& irUpBarriers,
		const  ARM_GP_Vector& notionals,
		const  ARM_PricingStatesPtr& states,
		ARM_Vector* eachFixingPrices,
        ARM_PricingContext* context) const

{
	//ARM_DateStrip dateStrip();
	/*ARM_Date asOfDate	= GetMktDataManager()->GetAsOfDate();
	ARM_Date startDate	= ARM_Date(startTime);
	ARM_Date endDate	= ARM_Date(endTime);

	if (	(fxDownBarriers.size()	!=	1) || 
			(fxUpBarriers.size()	!=	1) ||
			(irUpBarriers.size()	!=	1) ||
			(notionals.size()		!=	1) ||
			(irDownBarriers.size()	!=	1) )

		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Local Normal Model : only one barrier allowed by periods");

	ARM_GP_Vector dates(1,0.0);
	ARM_GP_Vector zeroValues(1,0.0);
	ARM_GP_Vector oneValues(1,1.0);
	
	ARM_Curve fxUpbarrierCv		= ARM_Curve(dates,fxUpBarriers);
	ARM_Curve fxDownbarrierCv	= ARM_Curve(dates,fxDownBarriers);
	ARM_Curve irUpbarrierCv		= ARM_Curve(dates,irUpBarriers);
	ARM_Curve irDownbarrierCv	= ARM_Curve(dates,irDownBarriers);
	ARM_Curve notionalCv		= ARM_Curve(dates,notionals);

	//Only used in fxvanillacalculator, so initialized at default values (not used)
	ARM_Curve strike			= ARM_Curve(dates,zeroValues);
	ARM_Curve alpha				= ARM_Curve(dates,zeroValues);
	ARM_Curve beta				= ARM_Curve(dates,zeroValues);
	ARM_Curve strike2			= ARM_Curve(dates,zeroValues);
	ARM_VanillaType vanillaType = ARM_FXVanillaType::vanilla;
	ARM_BasketType minMax		= ARM_FXBasketType::max;
	ARM_DigitType digitType		= ARM_FXDigitType::centred;
	string fx2Name				= "";
	string payIdx				= "";
	string payIdxIT				= "";
	double payIdxSpread			= 0.0;
	int callPut					= 1;
	int callPut2				= 1;

	//calculator constructor need leverage, but not the key world, only the notional is adjusted
	ARM_Curve leverage			= ARM_Curve(dates,oneValues);

	string fx1Name				= fxModelName;*/
	
	
	return ARM_VectorPtr(0);
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: DoubleDigital
///	Returns: double digital condition values
///	Action : Implicitly option expiry = eval date
///			 Computes the double digital price on
///			 both rates given their values
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::DoubleDigital(
		const string& modelName, 
		double evalTime,
		const ARM_VectorPtr& firstRate,
        const ARM_GP_Vector& firstStrikeDown,
        const ARM_GP_Vector& firstStrikeUp,
		double firstStrikeSpread,
		const ARM_VectorPtr& secondRate,
        const ARM_GP_Vector& secondStrikeDown,
        const ARM_GP_Vector& secondStrikeUp,
		double secondStrikeSpread,
        const ARM_PricingStatesPtr& states) const
{
	/// Check inputs
	if(!GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Correlation))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "Correlation data is missing in Double Digital evaluation");
	}

	ARM_MultiAssetsModel* refModel = dynamic_cast<ARM_MultiAssetsModel*>(GetNumericalModel());
	if(!refModel)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "MultiAssets models to restore marginal models is missing");
	}

	const ARM_GP_Vector& resets = ResetTimes();
	size_t prevIdx,nextIdx,nbResets=resets.size();

	if(evalTime<0 || evalTime < resets[0] - MAX_FUNCTIONAL_MISMATCH)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "Double digital evalTime negative or below first calibrated reset");
	if(evalTime>resets[nbResets-1] + MAX_FUNCTIONAL_MISMATCH)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Double digital valTime beyond horizon");


	/// Find reset dates around eval date
	double ratio=1.0;
	if(ExistsInVector(resets,evalTime))
	{
		nextIdx = IdxFromValue(resets,evalTime);
		prevIdx = nextIdx;
	}
	else
	{
		nextIdx=0;
		while(nextIdx < nbResets && evalTime > resets[nextIdx])
			++nextIdx;
		if(nextIdx>=nbResets)
		{
			prevIdx=nbResets-1;
			nextIdx=prevIdx;
		}
		else if(nextIdx==0)
			prevIdx=0;
		else
		{
			prevIdx = nextIdx - 1;
			ratio = (evalTime - resets[prevIdx])/(resets[nextIdx]-resets[prevIdx]);
		}
	}

	/// Restore from the multi-asset reference model both marginal models.
	/// The diffusion model is not used because to simplify DoubleDigital keyword
	/// both rate values are input and assumed to be generated by this model.
	/// Alternative : input all features of both rates (= spread rates !)
	/// an compute their values through the diffusion model and states

	const ARM_ModelNameMap& modelMap = *(refModel->GetModelMap());
	size_t i,nbModels = modelMap.size();
	if(nbModels < 3)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "Marginal models are missing");
	}
	ARM_Local_Normal_Model *localModel1	= dynamic_cast<ARM_Local_Normal_Model*>(&* (modelMap[1]->Model()) );
	ARM_Local_Normal_Model *localModel2	= dynamic_cast<ARM_Local_Normal_Model*>(&* (modelMap[2]->Model()) );
	if(!localModel1 || !localModel2)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "Marginal models are not local normal models");
	}

	/// Restore gaussian draws of each rates, decorrelate then recorrelate
	/// w.r.t. market normal correlation
	double fwd1		= localModel1->GetFwd(nextIdx);
	double vol1		= localModel1->GetVol(nextIdx);
	double shift1	= localModel1->GetShift(nextIdx);

	double fwd2		= localModel2->GetFwd(nextIdx);
	double vol2		= localModel2->GetVol(nextIdx);
	double shift2	= localModel2->GetShift(nextIdx);

	double modelCorrel	= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(COPULA_MODEL_CORREL,resets[nextIdx],0);
	double mktCorrel	= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(COPULA_MARKET_CORREL,resets[nextIdx],0);

	/// Interpolation if eval date is not a calibrated reset date
	if(ratio != 1.0)
	{
		double prevRatio = 1-ratio;
		fwd1	= ratio*fwd1 + prevRatio*localModel1->GetFwd(prevIdx);
		vol1	= ratio*vol1 + prevRatio*localModel1->GetVol(prevIdx);
		shift1	= ratio*shift1 + prevRatio*localModel1->GetShift(prevIdx);

		fwd2	= ratio*fwd2 + prevRatio*localModel2->GetFwd(prevIdx);
		vol2	= ratio*vol2 + prevRatio*localModel2->GetVol(prevIdx);
		shift2	= ratio*shift2 + prevRatio*localModel2->GetShift(prevIdx);

		double prevModelCorrel	= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(COPULA_MODEL_CORREL,resets[prevIdx],0);
		double prevMktCorrel	= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(COPULA_MARKET_CORREL,resets[prevIdx],0);
		modelCorrel				= ratio*modelCorrel + prevRatio*prevModelCorrel;
		mktCorrel				= ratio*mktCorrel + prevRatio*prevMktCorrel;

	}

	double racVar1=0,lnDrift1=0.0;
	bool isSLN1=fabs(shift1)>=1.e-12;
	if(isSLN1)
	{
		racVar1	= vol1*shift1;
		lnDrift1= 0.5*racVar1*racVar1;
	}
	else
	{
		racVar1 = vol1*fwd1;
	}

	double racVar2=0,lnDrift2=0.0;
	bool isSLN2=fabs(shift2)>=1.e-12;
	if(isSLN2)
	{
		racVar2	= vol2*shift2;
		lnDrift2= 0.5*racVar2*racVar2;
	}
	else
	{
		racVar2 = vol2*fwd2;
	}

	/// Recorrelate w.r.t market correlation
	double x,y,yuc;
	size_t nbStates=states->size();
	if(-0.999 < modelCorrel && modelCorrel < 0.999)
	{
		double modelOrthCoef	= sqrt(1-modelCorrel*modelCorrel);
		double mktOrthCoef		= sqrt(1-mktCorrel*mktCorrel);
		for(i=0;i<nbStates;++i)
		{
			if(isSLN1)
				x = (log((shift1*(*firstRate)[i]+(1.-shift1)*fwd1)/fwd1)+lnDrift1)/racVar1;
			else
				x = ((*firstRate)[i]-fwd1)/racVar1;

			if(isSLN2)
				y = (log((shift2*(*secondRate)[i]+(1.-shift2)*fwd2)/fwd2)+lnDrift2)/racVar2;
			else
				y = ((*secondRate)[i]-fwd2)/racVar2;

			yuc = (y - modelCorrel*x)/modelOrthCoef;
			y	= mktCorrel*x + mktOrthCoef*yuc;

			if(isSLN2)
				(*secondRate)[i] = fwd2/shift2*( exp(-lnDrift2+racVar2*y) - 1 + shift2 );
			else
				(*secondRate)[i] = fwd2 + y*racVar2;
		}
	}

	/// Draw market rates and compute payoff
	ARM_VectorPtr mktRate1 = localModel1->Func(evalTime,firstRate);
	ARM_VectorPtr mktRate2 = localModel2->Func(evalTime,secondRate);

	ARM_GP_Vector* result = new ARM_GP_Vector(nbStates);
	if(firstStrikeSpread==0.0 && secondStrikeSpread==0.0)
	{
		/// Digital on both sides
		bool cond1,cond2;
		for(i=0;i<nbStates;++i)
		{
			cond1 = firstStrikeDown[i] < (*mktRate1)[i] && (*mktRate1)[i] < firstStrikeUp[i];
			cond2 = secondStrikeDown[i] < (*mktRate2)[i] && (*mktRate2)[i] < secondStrikeUp[i];
			(*result)[i] = cond1 && cond2 ? 1.0 : 0.0;
		}
	}
	else
	{
		/// Call spreads
		double vi,k,kD,kU,payoff1,payoff2;
		double range1 = 2*firstStrikeSpread;
		double range2 = 2*secondStrikeSpread;
		for(i=0;i<nbStates;++i)
		{
			kD = firstStrikeDown[i];
			kU = firstStrikeUp[i];
			if(kD>kU)
			{
				k = kD; kD = kU; kU = k;
			}
			if(firstStrikeSpread>0)
			{
				vi = (*mktRate1)[i] - kD;
				if(vi < -firstStrikeSpread)
					payoff1 = 0.0;
				else if(vi < firstStrikeSpread)
					payoff1 = (vi + firstStrikeSpread)/range1;
				else
					payoff1 = 1.0;

				vi = (*mktRate1)[i] - kU;
				if(vi >= -firstStrikeSpread && vi < firstStrikeSpread)
					payoff1 -= (vi + firstStrikeSpread)/range1;
				else if(vi >= firstStrikeSpread)
					payoff1 -= 1.0;
			}
			else
				payoff1 = (kD <= (*mktRate1)[i] && (*mktRate1)[i] <= kU) ? 1 : 0;

			kD = secondStrikeDown[i];
			kU = secondStrikeUp[i];
			if(kD>kU)
			{
				k = kD; kD = kU; kU = k;
			}
			if(secondStrikeSpread>0)
			{
				vi = (*mktRate2)[i] - kD;
				if(vi < -secondStrikeSpread)
					payoff2 = 0.0;
				else if(vi < secondStrikeSpread)
					payoff2 = (vi + secondStrikeSpread)/range2;
				else
					payoff2 = 1.0;

				vi = (*mktRate2)[i] - kU;
				if(vi >= -secondStrikeSpread && vi < secondStrikeSpread)
					payoff2 -= (vi + secondStrikeSpread)/range2;
				else if(vi >= secondStrikeSpread)
					payoff2 -= 1.0;
			}
			else
				payoff2 = (kD <= (*mktRate2)[i] && (*mktRate2)[i] <= kU) ? 1 : 0;

			(*result)[i] = payoff1 * payoff2;
		}
	}

	return ARM_GP_VectorPtr(result);
}



////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Local_Normal_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Local Normal Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Local_Normal_Model::ValidateModelParams(const ARM_ModelParams& params) const
{	
	if( !dynamic_cast<const ARM_Local_Normal_ModelParams*>(&params) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_Local_Normal_ModelParams" );
	return true;
}



////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultModelParams
///	Returns : ARM_Local_Normal_ModelParams*
///	Action  : static method to build a default ModelParams 
///			  for ARM_Local_Normal_Model 
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams* ARM_Local_Normal_Model::CreateDefaultModelParams ()
{
	ARM_ModelParam* adj = CreateDefaultForwardAdjustmentModelParam ();
	ARM_ModelParam* vol = CreateDefaultVolatilityModelParam ();
	ARM_ModelParam* shift = CreateDefaultShiftModelParam ();
	ARM_ModelParamVector params(3);
	params[0] = adj;
	params[1] = vol;
	params[2] = shift;

	ARM_Local_Normal_ModelParams* modelParams =  new ARM_Local_Normal_ModelParams (params);
	
	/// params[0] & params[1] are clone in ARM_ModelParams constructor...
	delete adj;
	delete vol;
	delete shift;

	return modelParams;
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultForwardAdjustmentModelParam
///	Returns : ARM_ModelParam*
///	Action  : static method to build a default 
///			  ForwardAdjustment ModelParam 
////////////////////////////////////////////////////
ARM_ModelParam* ARM_Local_Normal_Model::CreateDefaultForwardAdjustmentModelParam ()
{
	// sizeAdj need to be increased if we need more that 2 adj surfaces
	size_t i ;
	size_t sizeAdj = SOCOR_ADJ_SIZE;
		
	ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
 	ARM_SurfaceWithInterpol surface (ARM_GP_T_Vector<double> (0),
									 ARM_GP_T_Vector<double> (0), 
									 ARM_GP_T_Matrix<double> (0, 0),
									 type);
	
	// build ForwardAdjustment model param
	ARM_SurfacePtrVector surfaceList(0);
	ARM_GP_Vector index(0);

	for (i=0; i<sizeAdj; i++)
	{		
		ARM_SurfacePtr surfacePtr =  ARM_SurfacePtr (static_cast <ARM_SurfaceWithInterpol*> (surface.Clone()));
		surfaceList.push_back (surfacePtr);
		index.push_back(i);
	}
		
	return new ARM_SurfaceListModelParam(   ARM_ModelParamType::ForwardAdjustment,
											index,
											surfaceList);
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultShiftModelParam
///	Returns : ARM_ModelParam*
///	Action  : static method to build a default 
///			  Barrier Shift ModelParam 
////////////////////////////////////////////////////
ARM_ModelParam* ARM_Local_Normal_Model::CreateDefaultShiftModelParam ()
{
	// sizeAdj need to be increased if we need more that 2 adj surfaces
	size_t i ;
	size_t sizeAdj = SOCOR_KADJ_SIZE;
		
	ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
 	ARM_SurfaceWithInterpol surface (ARM_GP_T_Vector<double> (0),
									 ARM_GP_T_Vector<double> (0), 
									 ARM_GP_T_Matrix<double> (0, 0),
									 type);
	
	// build ForwardAdjustment model param
	ARM_SurfacePtrVector surfaceList(0);
	ARM_GP_Vector index(0);

	for (i=0; i<sizeAdj; i++)
	{		
		ARM_SurfacePtr surfacePtr =  ARM_SurfacePtr (static_cast <ARM_SurfaceWithInterpol*> (surface.Clone()));
		surfaceList.push_back (surfacePtr);
		index.push_back(i);
	}
		
	return new ARM_SurfaceListModelParam(   ARM_ModelParamType::Shift,
											index,
											surfaceList);
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultVolatilityModelParam
///	Returns : ARM_ModelParam*
///	Action  : static method to build a default 
///			  Volatility ModelParam 
////////////////////////////////////////////////////
ARM_ModelParam* ARM_Local_Normal_Model::CreateDefaultVolatilityModelParam ()
{
	// sizeAdj need to be increased if we need more that 2 vol surfaces
	size_t i ;
	size_t sizeVol = SOCOR_VOL_SIZE;
		
	ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
 	ARM_SurfaceWithInterpol surface (	ARM_GP_T_Vector<double> (0),
										ARM_GP_T_Vector<double> (0), 
										ARM_GP_T_Matrix<double> (0, 0),
										type);
	
	// build ForwardAdjustment model param
	ARM_SurfacePtrVector surfaceList(0);
	ARM_GP_Vector index(0);

	for (i=0; i<sizeVol; i++)
	{		
		ARM_SurfacePtr surfacePtr =  ARM_SurfacePtr (static_cast <ARM_SurfaceWithInterpol*> (surface.Clone()));
		surfaceList.push_back (surfacePtr);
		index.push_back(i);
	}
		
	return new ARM_SurfaceListModelParam(   ARM_ModelParamType::Volatility,
											index,
											surfaceList);
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultModelParams
///	Returns : ARM_Local_Normal_ModelParams*
///	Action  : static method to build a default ModelParams 
///			  for ARM_Local_Normal_Model 
////////////////////////////////////////////////////
ARM_Local_Normal_ModelParams* ARM_Local_Normal_Model::CreateDefaultModelParams(const ARM_IntVector& paramTypes, const ARM_IntVector& paramNbSurfaces)
{
	size_t i,nbParams = paramTypes.size();
	if(nbParams != paramNbSurfaces.size())
		ARM_THROW( ERR_INVALID_ARGUMENT, " : ParamTypes & ParamNbSurfaces vectors must have the same size" );

	ARM_ModelParamVector params;
	for(i=0;i<nbParams;++i)
		params.push_back(CreateDefaultModelParam(static_cast<ARM_ModelParamType::ParamNb>(paramTypes[i]),paramNbSurfaces[i]));

	ARM_Local_Normal_ModelParams* modelParams =  new ARM_Local_Normal_ModelParams(params,paramTypes);
	
	/// Free memory because cloned by previous constructor
	DeletePointorVector<ARM_ModelParam>(params);

	return modelParams;
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CreateDefaultModelParam
///	Returns : ARM_ModelParam*
///	Action  : static method to build a default ModelParam
///			  with the input type
////////////////////////////////////////////////////
ARM_ModelParam* ARM_Local_Normal_Model::CreateDefaultModelParam(ARM_ModelParamType::ParamNb paramType, int nbSurfaces)
{		
	ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;
 	ARM_SurfaceWithInterpol surface (ARM_GP_T_Vector<double> (0),
									 ARM_GP_T_Vector<double> (0), 
									 ARM_GP_T_Matrix<double> (0, 0),
									 type,0.0);
	
	// Build each model param with its required surfaces
	ARM_SurfacePtrVector surfaceList(0);
	ARM_GP_Vector index(0);

	for(size_t i=0; i<nbSurfaces; i++)
	{		
		ARM_SurfacePtr surfacePtr =  ARM_SurfacePtr (static_cast <ARM_SurfaceWithInterpol*> (surface.Clone()));
		surfaceList.push_back (surfacePtr);
		index.push_back(i);
	}
		
	return new ARM_SurfaceListModelParam(paramType,index,surfaceList);
}



////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration (given a calculator -> computation of target prices)
///				used for FX and Libor dual condition
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateLocalModel(const ARM_GenCalculator& calculator, const ARM_GP_Vector& evalTimes) 
{
	if(ARM_FXRACalculator* fxraCalculator = dynamic_cast<ARM_FXRACalculator*>((ARM_GenCalculator*)&calculator))
	{
		fxraCalculator->Price();
		const ARM_Vector* allFixingPrices = fxraCalculator->GetvAllFixingPrices();

		const size_t EVALTIME=1;
		const size_t FXPRICEBDOWN=2;
		const size_t FXPRICEBUP=3;
		const size_t IRPRICEBDOWN=4;
		const size_t IRPRICEBUP=5;
		const size_t TOTAL=5;

		//Get the barriers from the calculator
		ARM_Curve fxDownBarrierCv	= fxraCalculator->GetFXDownBarrier();
		ARM_Curve fxUpBarrierCv		= fxraCalculator->GetFXUpBarrier();
		ARM_Curve irDownBarrierCv	= fxraCalculator->GetIRDownBarrier();
		ARM_Curve irUpBarrierCv		= fxraCalculator->GetIRUpBarrier();

		size_t nbEval = evalTimes.size();
		double pricefxBDown, pricefxBUp, priceirBDown, priceirBUp;
		double volfxBDownTgt, volfxBUpTgt, volirBDownTgt, volirBUpTgt;
		double fxDownBarrier, fxUpBarrier, irDownBarrier, irUpBarrier;
		size_t size = allFixingPrices->size();

		double fixingTime;
		size_t nbFixIdx=0;
		size_t nbFixings=0; 

		for (size_t i=0; i<nbEval; i++)
		{	
			double evalTime = evalTimes[i];
			fxDownBarrier	= fxDownBarrierCv.Interpolate(evalTime);
			fxUpBarrier		= fxUpBarrierCv.Interpolate(evalTime);
			irDownBarrier	= irDownBarrierCv.Interpolate(evalTime);
			irUpBarrier		= irUpBarrierCv.Interpolate(evalTime);


			nbFixIdx = (nbFixings*TOTAL + 1)*i;
			nbFixings=(*allFixingPrices)[nbFixIdx];

			for (size_t p=0; p<nbFixings; p++)
			{
				
				fixingTime = (*allFixingPrices)[nbFixIdx + p*TOTAL + EVALTIME];

				//target prices 
				pricefxBDown		= (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPRICEBDOWN];
				pricefxBUp			= (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPRICEBUP];
				priceirBDown		= (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPRICEBDOWN];
				priceirBUp			= (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPRICEBUP];

				//Compute target normal volatilities associated with theses prices
				//volfxBDownTgt = VanillaImpliedVol_N(fwdFX, pricefxBDown, fxBDown, maturity, 1);

			}
			
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: calculator not supported" );
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : ResetModelParams
///	Returns : void
///	Action  : resets model params (to be used before
///			  calibration)
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::ResetModelParams()
{
	size_t i,nbParams = itsParamTypes.size();
	for(i=0;i<nbParams;++i)
	{
		if (!GetModelParams()->DoesModelParamExist(itsParamTypes[i]))
		{
			string msge("ARM_Local_Model : Missing ModelParam #" + i);
			ARM_THROW( ERR_INVALID_ARGUMENT, msge );
		}

		ARM_ModelParam& param = GetModelParams()->GetModelParam(itsParamTypes[i]);

		ARM_SurfaceListModelParam* paramSurf = dynamic_cast < ARM_SurfaceListModelParam* > (&param);
		if (!paramSurf)
		{
			string msge("ARM_Local_Model : ModelParam #" + i);
			msge += " is required to be a ARM_SurfaceListModelParam";
			ARM_THROW( ERR_INVALID_ARGUMENT, msge  );
		}

		// Create default params
		ARM_ModelParam* newParam = CreateDefaultModelParam(static_cast<ARM_ModelParamType::ParamNb>(itsParamTypes[i]),itsParamNbSurfaces[i]);
		GetModelParams()->DeleteModelParam(itsParamTypes[i]);
		GetModelParams()->SetModelParam(newParam);

		/// Free new param because cloned by previous method
		delete newParam;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration 
///			  (for 1 given ARM_Security)
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateLocalModel(const ARM_Security& security, 
	double targetPrice, 
	const ARM_GP_Vector& evalTimes,
	size_t secIdx)
{	
	/// Switch on security type
	if(ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>((ARM_Security*)&security))
	{			
		//----------------------------
		//--- SpreadOption case
		//----------------------------

		ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>((ARM_Security*)&security);

		if (spreadOption && !(spreadOption->IsCorridorSpread()) && !(spreadOption->IsDigitalFLT()))
		{
			CalibrateSpreadOption(	spreadOption, 
									targetPrice,
									evalTimes);
		}
		//----------------------------
		//--- Spread Option / Corridor Spread Option / Double Condition Corridor Option case
		//----------------------------
		else if (spreadOption && (spreadOption->IsCorridorSpread()))
		{
			if(spreadOption->IsCorridorDblCondition())
				CalibrateDoubleCorridorOption(static_cast<ARM_CorridorDblCondition*>(spreadOption),targetPrice,evalTimes,secIdx);
			else
				CalibrateCorridorSpreadOption(spreadOption,targetPrice,evalTimes,secIdx);
		}
		//------------------------------
		//--- Inflation: not supported
		//------------------------------
		else if(dynamic_cast<ARM_InfCapFloor*>((ARM_Security*)&security))
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model: calibration to inflation cap/floors not supported" );
		}
		//----------------------------
		//--- VanillaCaplet case
		//----------------------------
		else
		{
			CalibrateCapFloor(	capFloor, 
								targetPrice,
								evalTimes);
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: instrument not supported" );
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : Func
///	Returns : ARM_VectorPtr
///	Action  : evaluation of functional
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Normal_Model::Func(double evalTime,const ARM_GP_VectorPtr& values) const
{
	ARM_GP_Vector reset = ResetTimes();
	if(evalTime<0 || evalTime < reset[0] - MAX_FUNCTIONAL_MISMATCH)
	    ARM_THROW( ERR_INVALID_ARGUMENT, "Func: evalTime negative or below first calibrated reset!");
	size_t nbResets = reset.size();
	if(evalTime>reset[nbResets-1] + MAX_FUNCTIONAL_MISMATCH)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Func: evalTime beyond horizon!");

	size_t newSize = values->size();
	size_t gridSize = GetGridSize();
	double storageXmin = GetGrid(0);
	double storageXmax = GetGrid(gridSize-1);
	double storageDx = (storageXmax - storageXmin) / (gridSize-1);

	if (ExistsInVector(reset,evalTime))
	{
		size_t index = IdxFromValue(reset,evalTime);
	
		double curNumMethState, prevNumMethState, nextNumMethState;
		double fwd = GetFwd(index);
		double vol = GetVol(index);
		double shift = GetShift(index);

		double racVar = vol*shift;
		double lndrift = 0.5*racVar*racVar;

		ARM_GP_Vector* result = new ARM_GP_Vector( newSize );

		size_t j;
		for (size_t i = 0; i<newSize; i++)
		{	
			//NEW : index is shifted lognormal with relative shift
			//		dS = (mS+(1-m)S0)*vol*dW
			if (fabs(shift)<1.e-12)
				curNumMethState = ((*values)[i]-fwd)/vol/fwd;
			else
				curNumMethState = (log((shift*(*values)[i]+(1.-shift)*fwd)/fwd)+lndrift)/racVar;
			
			/// flat extrapol if out of range
			if (curNumMethState<=storageXmin)
				result->Elt(i) =  FuncValue(index,0);
			else if (curNumMethState>=storageXmax)
				result->Elt(i) =  FuncValue(index,gridSize-1);
			else
			{
				//standard interpol
				j = (size_t)ceil((curNumMethState - storageXmin) / storageDx) ;
				nextNumMethState = GetGrid(j);
				prevNumMethState = GetGrid(j-1);
				
				result->Elt(i) = (  (curNumMethState - prevNumMethState) * FuncValue(index,j)
								  + (nextNumMethState - curNumMethState) * FuncValue(index,j-1) ) / (nextNumMethState - prevNumMethState);
			}
		}
// FIXMEFRED: mig.vc8 (30/05/2007 16:15:42):cast
		return static_cast<ARM_VectorPtr>(result);
		return ARM_GP_VectorPtr(result);


	}
	else
	{
		size_t iLeft,iRight=0;
		while(iRight < nbResets && evalTime > reset[iRight] )
			iRight++;
		if(iRight==0)
			iLeft = 0;
		else if(iRight>=nbResets)
		{
			iLeft = nbResets - 1;
			iRight = iLeft;
		}
		else
			iLeft = iRight - 1;

		double ratio = iLeft == iRight ? 1.0 : ( evalTime - reset[iLeft] ) / ( reset[iRight] - reset[iLeft] ) ;
		
		ARM_GP_VectorPtr result = Func(reset[iLeft],values);
		if(ratio!=1.0)
		{
			ARM_GP_VectorPtr resRight = Func(reset[iRight],values);
			for (size_t i = 0; i<newSize; i++)
				result->Elt(i) = ( 1. - ratio ) * result->Elt(i) + ratio * resRight->Elt(i);
		}

		return result;
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration  to functional (for 1 given ARM_Security)
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateLocalModel (const ARM_Security& security, ARM_DensityFunctor& density, bool rescaling)
{
	// will generate an exception if numerical model is not supported for calibration
	Local_Normal_Model_Calibration_Helper helper (GetNumericalModel());

	double fwd,vol,fwdTH,shift;
	double resetTime;

	ARM_SpreadOption* spreadOption = dynamic_cast<ARM_SpreadOption*>((ARM_Security*)&security);
	ARM_CapFloor* capFloor = dynamic_cast<ARM_CapFloor*>((ARM_Security*)&security);
	ARM_CorridorDblCondition* ra2 = dynamic_cast<ARM_CorridorDblCondition*>((ARM_Security*)&security);
	if(ra2)
	{
		/// Very simple mapping at each expiry date : marginals are assumed
		/// already calibrated and market co-distribution is given via
		/// gaussian copula :
		///		1)	restore the gaussian correlation used to price
		///			from the RA2 input security
		///		2)	compute the numerical model (gaussian) correlation between
		///			both indexes of the RA2

		/// Build correlation surface if necessary
		if(!GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Correlation))
		{
			ARM_ModelParam* correlParam = CreateDefaultModelParam(ARM_ModelParamType::Correlation,COPULA_CORREL_SIZE);
			GetModelParams()->SetModelParam(correlParam);
			delete correlParam;
			itsParamTypes.push_back(ARM_ModelParamType::Correlation);
			itsParamNbSurfaces.push_back(COPULA_CORREL_SIZE);
		}

		// Get product datas
		ARM_SpreadOption* liborCondLeg	= ra2->GetSpreadDigital();
		ARM_SpreadOption* soCondLeg		= static_cast<ARM_SpreadOption*>(ra2);

		int callPut;
		double asOfDate = GetAsOfDate().GetJulian();
		double resetTime,payTime,payPeriod,payNotional;

		double soCoeffLong,soCoeffShort,soStrike;
		double soLongFloatStartTime,soLongFloatEndTime,soShortFloatStartTime,soShortFloatEndTime;
		ARM_GP_Vector soLongFixPayTimes,soLongFixPayPeriods,soShortFixPayTimes,soShortFixPayPeriods; 

		double libCoeffLong,libCoeffShort,libStrike;
		double libLongFloatStartTime,libLongFloatEndTime,libShortFloatStartTime,libShortFloatEndTime;
		ARM_GP_Vector libLongFixPayTimes,libLongFixPayPeriods,libShortFixPayTimes,libShortFixPayPeriods; 

		ARM_GP_Vector resets(To_ARM_GP_Vector(*(ra2->GetSwapLeg()->GetResetDates())));
		size_t i,nbResets = resets.size();

		const ARM_Vector& marketNormalCorrels = ra2->GetSpreadDigitalCorrel();
		double soVolModel,libVolModel,soLibCovModel,modelCorrel;
		for(i=0;i<nbResets;++i)
		{
			/// Save market gaussian correlation in a single column correlation surface
			resetTime = resets[i] - asOfDate;
			GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(COPULA_MARKET_CORREL,resetTime,0,marketNormalCorrels[i]);
			setTime(resetTime); // doublon... no matter

			GetSpreadOptionLetData(	soCondLeg,
									i,
									callPut, 
									resetTime, 
									payTime, 
									payPeriod, 
									payNotional, 
									soCoeffLong, 
									soCoeffShort, 
									soStrike, 
									soLongFloatStartTime, 
									soLongFloatEndTime, 
									soLongFixPayTimes, 
									soLongFixPayPeriods, 
									soShortFloatStartTime, 
									soShortFloatEndTime, 
									soShortFixPayTimes, 
									soShortFixPayPeriods);

			soVolModel = helper.SpreadNormalVolatilityFromNumericalModel(
									resetTime,
									payTime,
									soCoeffLong,
									soCoeffShort,
									soStrike, 
									soLongFloatStartTime, 
									soLongFloatEndTime, 
									soLongFixPayTimes, 
									soLongFixPayPeriods, 
									soShortFloatStartTime, 
									soShortFloatEndTime, 
									soShortFixPayTimes, 
									soShortFixPayPeriods);

			GetSpreadOptionLetData(	liborCondLeg,
									i,
									callPut, 
									resetTime, 
									payTime, 
									payPeriod, 
									payNotional, 
									libCoeffLong, 
									libCoeffShort, 
									libStrike, 
									libLongFloatStartTime, 
									libLongFloatEndTime, 
									libLongFixPayTimes, 
									libLongFixPayPeriods, 
									libShortFloatStartTime, 
									libShortFloatEndTime, 
									libShortFixPayTimes, 
									libShortFixPayPeriods);

			libVolModel = helper.SpreadNormalVolatilityFromNumericalModel(
									resetTime,
									payTime,
									libCoeffLong,
									libCoeffShort,
									libStrike, 
									libLongFloatStartTime, 
									libLongFloatEndTime, 
									libLongFixPayTimes, 
									libLongFixPayPeriods, 
									libShortFloatStartTime, 
									libShortFloatEndTime, 
									libShortFixPayTimes, 
									libShortFixPayPeriods);

			soLibCovModel = helper.SpreadRateNormalCovarianceFromNumericalModel(
									resetTime,
									payTime,
									soCoeffLong,
									soCoeffShort,
									soStrike,
									soLongFloatStartTime, 
									soLongFloatEndTime, 
									soLongFixPayTimes, 
									soLongFixPayPeriods, 
									soShortFloatStartTime, 
									soShortFloatEndTime, 
									soShortFixPayTimes,
									soShortFixPayPeriods, 
									libStrike,
									libLongFloatStartTime, 
									libLongFloatEndTime, 
									libLongFixPayTimes, 
									libLongFixPayPeriods, 
									libCoeffLong,
									libCoeffShort,
									libShortFloatStartTime, 
									libShortFloatEndTime, 
									libShortFixPayTimes, 
									libShortFixPayPeriods);

			/// Save model gaussian correlation in another single column correlation surface
			modelCorrel = soLibCovModel/(soVolModel*libVolModel*resetTime/K_YEAR_LEN);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(COPULA_MODEL_CORREL,resetTime,0,modelCorrel);
		}

	}
	else if (spreadOption)
	{
		size_t soSize = spreadOption->GetResetDates()->GetSize();
		if (soSize==1)
		{
			double fwd1 = (*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()))[0];
			double fwd2 = (*(spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()))[0];
			double weight1 = spreadOption->GetWeight1();
			double weight2 = spreadOption->GetWeight2();
			fwdTH = (weight2*fwd2 - weight1*fwd1)/100.;

			// get product data
			int callPut;
			double payTime, coeffLong, coeffShort, strike, payPeriod, payNotional;
			double longFloatStartTime, longFloatEndTime, shortFloatStartTime, shortFloatEndTime;
			ARM_GP_Vector longFixPayTimes, longFixPayPeriods, shortFixPayTimes, shortFixPayPeriods; 

			GetSpreadOptionLetData   (	spreadOption,
										0,
										callPut, 
										resetTime, 
										payTime, 
										payPeriod, 
										payNotional, 
										coeffLong, 
										coeffShort, 
										strike, 
										longFloatStartTime, 
										longFloatEndTime, 
										longFixPayTimes, 
										longFixPayPeriods, 
										shortFloatStartTime, 
										shortFloatEndTime, 
										shortFixPayTimes, 
										shortFixPayPeriods);

			double longCmsModel, shortCmsModel;
				
			longCmsModel = helper.CmsRateFromNumericalModel (   resetTime, 
																payTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longFixPayTimes, 
																longFixPayPeriods );

			shortCmsModel = helper.CmsRateFromNumericalModel (  resetTime, 
																payTime, 
																shortFloatStartTime, 
																shortFloatEndTime, 
																shortFixPayTimes, 
																shortFixPayPeriods );

			fwd = coeffLong * longCmsModel - coeffShort * shortCmsModel;

			double modelSpreadVol = helper.SpreadNormalVolatilityFromNumericalModel (  resetTime,
																					payTime,
																					coeffLong,
																					coeffShort,
																					strike,
																					longFloatStartTime,
																					longFloatEndTime,
																					longFixPayTimes,
																					longFixPayPeriods,
																					shortFloatStartTime,
																					shortFloatEndTime,
																					shortFixPayTimes,
																					shortFixPayPeriods);
			modelSpreadVol /= fwd;
			shift = 0.;

			vol = modelSpreadVol*sqrt(resetTime/K_YEAR_LEN);
			
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: too many flows in spread option" );
	}
	else if (capFloor)
	{
		size_t capSize = capFloor->GetResetDates()->GetSize();
		if (capSize==1)
		{
			
			// get computed libor rate (potentially adjusted)
			fwdTH = capFloor->GetSwapLeg()->GetFwdRates()->Elt(0) * 0.01;
			
			double fwdStartTime, fwdEndTime, fwdPeriod, payTime, payPeriod, payNotional, strike;
			int callPut;
	
			GetCapFloorLetData(	capFloor, 
								0,
								resetTime,
								payTime,
								payPeriod,
								payNotional,
								fwdStartTime,
								fwdEndTime,
								fwdPeriod,
								strike,
								callPut);
	
			/// adjusment on libor fwd
			fwd = helper.LiborRateFromNumericalModel(	resetTime, 
														payTime, 
														fwdStartTime, 
														fwdEndTime, 
														fwdPeriod, 
														strike);

			double modelVol;

			if (!LiborAsShiftedLognormal())
			{
				modelVol = helper.LiborNormalVolatilityFromNumericalModel (	resetTime, 
																			fwdStartTime, 
																			fwdEndTime, 
																			fwdPeriod, 
																			strike);
				modelVol /= fwd;
				shift = 0.;
			}
			else
			{
				modelVol = helper.FwdZcVolatilityFromNumericalModel (	resetTime, 
																		fwdStartTime, 
																		fwdEndTime, 
																		fwdPeriod, 
																		strike);
				shift = 1./fwdPeriod; 
				shift = fwd/(fwd+shift);
				modelVol /= shift;
			}
			vol = modelVol*sqrt(resetTime/K_YEAR_LEN);
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: too many flows in cap" );
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: only SO, caplet or RA2!!" );

	if(!ra2)
	{
		double mat = resetTime/K_YEAR_LEN;

		ARM_GP_VectorPtr proba = ARM_GP_VectorPtr(new ARM_GP_Vector(GetGridSize()));
		ARM_GP_VectorPtr grid = Grid();

		ARM_GP_Vector::iterator iter = proba->begin();
		ARM_GP_Vector::iterator iter2 = grid->begin();
		for (;iter!=proba->end();++iter,++iter2)
			(*iter) = NormalCDF((*iter2));

		ARM_GP_VectorPtr func = density.Quantile(proba,fwdTH,mat);

		setVol(vol);
		setFwd(fwd);
		setShift(shift);
		setFunc(func);
		setTime(resetTime);

		if (rescaling)
		{
			double numStateMin = grid->Elt(0);
			double numStateMax = grid->Elt(grid->size()-1);

			GaussLegendre_Coefficients glc( PRICGLPOINTSNB, numStateMin, numStateMax);
			ARM_GP_Vector states (PRICGLPOINTSNB);
			for (size_t i(0); i<PRICGLPOINTSNB; i++)
				states[i] = glc.get_point(i);

			ARM_GP_VectorPtr rate(new ARM_GP_Vector(PRICGLPOINTSNB));
			if (fabs(shift)<1e-12)
				for (i=0; i<PRICGLPOINTSNB; i++)
					(*rate)[i] = fwd*(1.+vol*states[i]);
			else
				for (i=0; i<PRICGLPOINTSNB; i++)
					(*rate)[i] = fwd/shift*(exp(shift*vol*states[i]-0.5*shift*shift*vol*vol)-(1.-shift));

			double state,result=0.,resultByState;
			ARM_GP_VectorPtr index = Func(resetTime,rate);
			for( i=0 ; i<PRICGLPOINTSNB; ++i )
			{
				state = glc.get_point(i);
				resultByState = index->Elt(i);
				result += glc.get_weight(i) * resultByState * exp(-0.5*state*state);
			}
			result *= INVSQRT2PI;
			double coeff = fwdTH/result;

			int idx = itsFunctionals.size()-1;
			for (iter=itsFunctionals[idx]->begin();iter!=itsFunctionals[idx]->end();++iter)
				(*iter) *= coeff;
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateSpreadOption
///	Returns : void
///	Action  : Local model calibration to a 
///			  spread option 
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateSpreadOption(	ARM_SpreadOption* spreadOption,
													double targetPrice,
													const ARM_GP_Vector& evalTimes)
{
	/// target price is ignored because we need the pv of all spread optionlets
	/// --> use GetCashFlowValues()
	size_t soSize = spreadOption->GetResetDates()->GetSize();
	ARM_Vector* cashFlows = spreadOption->GetCashFlowValues();

	for (size_t i(0); i<soSize; i++)
		CalibrateSpreadOptionLet(spreadOption, i, cashFlows->Elt(i), evalTimes);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateSpreadOptionLet
///	Returns : void
///	Action  : Local model calibration to a 
///			  spread option let
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateSpreadOptionLet(	ARM_SpreadOption* spreadOption,
														size_t periodIdx,
														double targetPrice,
														const ARM_GP_Vector& evalTimes)
{	
	// will generate an exception if numerical model is not supported for calibration
	Local_Normal_Model_Calibration_Helper helper (GetNumericalModel());
			
	size_t i, sizeEvalTimes = evalTimes.size();	
	// get convexified CMS rates from ARM_SpreadOption
	double longCmsTarget  = spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(periodIdx);
	double shortCmsTarget = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(periodIdx);
	longCmsTarget  *= 0.01;
	shortCmsTarget *= 0.01;

	// get product data
	double payTime, resetTime, coeffLong, coeffShort, strike, payPeriod, payNotional;
	double longFloatStartTime, longFloatEndTime, shortFloatStartTime, shortFloatEndTime;
	int callPut;
	ARM_GP_Vector longFixPayTimes, longFixPayPeriods, shortFixPayTimes, shortFixPayPeriods; 

	GetSpreadOptionLetData   (	spreadOption,
								periodIdx,
								callPut, 
								resetTime, 
								payTime, 
								payPeriod, 
								payNotional, 
								coeffLong, 
								coeffShort, 
								strike, 
								longFloatStartTime, 
								longFloatEndTime, 
								longFixPayTimes, 
								longFixPayPeriods, 
								shortFloatStartTime, 
								shortFloatEndTime, 
								shortFixPayTimes, 
								shortFixPayPeriods);
		
	/// compute approx atm spread vol
	double Tfix	 = resetTime / K_YEAR_LEN;
	double atmShortVol  = spreadOption->GetVol1ATM(Tfix, (shortFloatEndTime - shortFloatStartTime)/K_YEAR_LEN);
	double atmLongVol   = spreadOption->GetVol2ATM(Tfix, (longFloatEndTime - longFloatStartTime)/K_YEAR_LEN);
	double correl		= spreadOption->GetSpreadCorrel(Tfix, (shortFloatEndTime - shortFloatStartTime)/K_YEAR_LEN, (longFloatEndTime - longFloatStartTime)/K_YEAR_LEN, 0.0);
	atmShortVol			*= 0.01 * coeffShort * shortCmsTarget ;
	atmLongVol			*= 0.01 * coeffLong  * longCmsTarget ;
	correl				*= 0.01;
	double atmSpreadVol   = sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
	double atmSpreadStdev = atmSpreadVol * sqrt(Tfix);
	
	double forwardSpread	= coeffLong * longCmsTarget - coeffShort * shortCmsTarget;

	bool calibrateVol (true);

	if (fabs(strike - forwardSpread)> NSTDEV_NO_CALIB * atmSpreadStdev)
		calibrateVol = false;

	/// compute target spread vol only if 
	/// not too far from atm
	double targetSpreadVol (0.0);
	
	if (calibrateVol)
	{
		double dfPay = GetNumericalModel()->GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
		double targetPrice_N = targetPrice / (dfPay * payNotional * payPeriod);
		targetSpreadVol		 = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strike, Tfix, callPut, &atmSpreadVol);
	}

	///
	///  compute forward spread vol & forward adjustments
	///
	ARM_GP_Vector VarianceSqueezInfo(3);
	VarianceSqueezInfo[0]=resetTime;
	for (i=0; i<sizeEvalTimes; i++)
	{
		// calibration only if option maturity is after evalTime
		if (resetTime >= evalTimes[i] )
		{
			VarianceSqueezInfo[1]=evalTimes[i];
			/// compute CMS adjustements
			double longCmsModel, shortCmsModel;
			
			longCmsModel = helper.CmsRateFromNumericalModel (   evalTimes[i], 
																payTime, 
																longFloatStartTime, 
																longFloatEndTime, 
																longFixPayTimes, 
																longFixPayPeriods );

			shortCmsModel = helper.CmsRateFromNumericalModel (  evalTimes[i], 
																payTime, 
																shortFloatStartTime, 
																shortFloatEndTime, 
																shortFixPayTimes, 
																shortFixPayPeriods );
				

			double longAdjustment  = longCmsTarget  - longCmsModel;
			double shortAdjustment = shortCmsTarget - shortCmsModel;

			// set adjusments to model params
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SO_ADJ_LONG, evalTimes[i], resetTime, longAdjustment);
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SO_ADJ_SHORT, evalTimes[i], resetTime, shortAdjustment);
			
			///
			/// calibrate forward spread vol adjustment
			/// only if strike is not too far from atm
			/// if too far, a forward vol = 0 is inputed
			///
			double forwardSpreadVol = 0.0;

			if (calibrateVol)
			{
				/// compute spread forward volatility
				double modelSpreadVol ; /// model spread vol up to evalTimes[i] 

				modelSpreadVol = helper.SpreadNormalVolatilityFromNumericalModel (  evalTimes[i],
																					payTime,
																					coeffLong,
																					coeffShort,
																					strike,
																					longFloatStartTime,
																					longFloatEndTime,
																					longFixPayTimes,
																					longFixPayPeriods,
																					shortFloatStartTime,
																					shortFloatEndTime,
																					shortFixPayTimes,
																					shortFixPayPeriods);
				
				double Tcall = evalTimes[i] / K_YEAR_LEN;
				if (Tcall==Tfix)
					Tfix += 1./K_YEAR_LEN;
				
				double var = (	targetSpreadVol * targetSpreadVol * Tfix 
							  -	modelSpreadVol  * modelSpreadVol  * Tcall) ;

				/// everything is fine
				if (var>=0)
					forwardSpreadVol = sqrt( var / (Tfix - Tcall) );
				
				/// oups: variance squeeze
				else 
				{
					VarianceSqueezInfo[2]=1.0;
					SetVarianceSqueezeStatus(true);
					CC_Ostringstream os;
					os << ARM_USERNAME << " : Local_Normal_Model / Spread Option Calibration : variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(Tfix*10.0)/10.0 << " Y]"<< endl;
					ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				}
				push_backVarianceSqueezInfo(VarianceSqueezInfo);
			}
			else
			{
#ifdef _DEBUG
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Local_Normal_Model / Spread Option Calibration : strike seems to be to far from ATM. Null forward vol inputed..."<< endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
#endif
			}

			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(SO_VOL, evalTimes[i], resetTime, forwardSpreadVol);
		}
		else
		{	
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SO_ADJ_LONG, evalTimes[i], resetTime, 0.0);
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SO_ADJ_SHORT, evalTimes[i], resetTime, 0.0);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(SO_VOL, evalTimes[i], resetTime, 0.0);
		}
	}

}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateCorridorSpreadOption
///	Returns : void
///	Action  : Local model calibration to a corridor
///			  spread option wich is in fact a corridor
///			  spread option
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateCorridorSpreadOption(ARM_SpreadOption* spreadOption, 
		double targetPrice, 
		const ARM_GP_Vector& evalTimes,
		size_t oldbarrierType)
{
	size_t barrierType = oldbarrierType%2;
	// will generate an exception if numerical model is not supported for calibration
	Local_Normal_Model_Calibration_Helper helper (GetNumericalModel());
			
	size_t i, j, sizeEvalTimes = evalTimes.size();

	// get product data
	int callPut;
	double spread1, spread2;
	ARM_GP_Vector payTimes, resetTimes, coeffsLong, coeffsShort, strikes, payPeriods, payNotionals, fixValues, payIndexLeverages;
	ARM_GP_Vector longFloatStartTimes, longFloatEndTimes, shortFloatStartTimes, shortFloatEndTimes;
	ARM_GP_Vector payFloatStartTimes, payFloatEndTimes;
	vector<ARM_GP_Vector> longSwapsFixPayTimes, longSwapsFixPayPeriods, shortSwapsFixPayTimes, shortSwapsFixPayPeriods; 
	vector<ARM_GP_Vector> paySwapsFixPayTimes, paySwapsFixPayPeriods; 
	ARM_GP_Vector payIndexResetTimes;
	ARM_IntVector periodIndexes;
	int payIndexType(0);	

	/// Unused double condition datas
	int rateCallPut;
	ARM_GP_Vector rateStrikes,rateFloatStartTimes,rateFloatEndTimes;
	vector<ARM_GP_Vector> rateFixPayTimes,rateFixPayPeriods;
	double rateSpread1,rateSpread2;


	GetCorridorSpreadOptionData   (	
							spreadOption, 
							callPut, 
							resetTimes, 
							payTimes, 
							payPeriods, 
							payNotionals, 
							coeffsLong, 
							coeffsShort, 
							strikes, 
							longFloatStartTimes,
							longFloatEndTimes,
							longSwapsFixPayTimes,
							longSwapsFixPayPeriods,
							shortFloatStartTimes,
							shortFloatEndTimes,
							shortSwapsFixPayTimes,
							shortSwapsFixPayPeriods,
							fixValues,
							spread1,
							spread2,
							payIndexType,
							payIndexLeverages,
							payFloatStartTimes, 
							payFloatEndTimes,
							paySwapsFixPayTimes,
							paySwapsFixPayPeriods,
							payIndexResetTimes,
							periodIndexes,
							rateCallPut,
							rateStrikes,
							rateFloatStartTimes,
							rateFloatEndTimes,
							rateFixPayTimes,
							rateFixPayPeriods,
							rateSpread1,rateSpread2);

	/// spreads are required to be 10 bp for the moment
	if (	(fabs(spread1 - (-SOCOR_STRIKE_SPREAD))>1e-12) 
		||	(fabs(spread2 -  SOCOR_STRIKE_SPREAD)>1e-12)   )
		ARM_THROW( ERR_INVALID_ARGUMENT, "CalibrateCorridorSpreadOption : corridor strike spreads are required to be = 1 bp" );

	
	/// switch index depending on CAP or FLOOR
	size_t VOL_LEFT;
	size_t VOL_RIGHT;
	size_t VOL_LEFT_FLT;
	size_t VOL_RIGHT_FLT;
	size_t ADJ_LONG_FLT;
	size_t ADJ_SHORT_FLT;
	size_t KADJ_LEFT;
	size_t KADJ_RIGHT;
	size_t KADJ_LEFT_FLT;
	size_t KADJ_RIGHT_FLT;

	// When we calibrate the low barrier volatility
	if (barrierType == SOCOR_DOWN)
	{
		VOL_LEFT		= SOCOR_VOL_DOWNLEFT;
		VOL_RIGHT		= SOCOR_VOL_DOWNRIGHT;
		VOL_LEFT_FLT	= SOCOR_VOL_DOWNLEFT_FLT;
		VOL_RIGHT_FLT	= SOCOR_VOL_DOWNRIGHT_FLT;
		ADJ_LONG_FLT	= SOCOR_ADJ_LONG_DOWN_FLT;
		ADJ_SHORT_FLT	= SOCOR_ADJ_SHORT_DOWN_FLT;
		KADJ_LEFT		= SOCOR_KADJ_DOWNLEFT;
		KADJ_RIGHT		= SOCOR_KADJ_DOWNRIGHT;
		KADJ_LEFT_FLT	= SOCOR_KADJ_DOWNLEFT_FLT;
		KADJ_RIGHT_FLT	= SOCOR_KADJ_DOWNRIGHT_FLT;
		

	}
	// When we calibrate the high barrier volatility
	else if (barrierType == SOCOR_UP)
	{
		VOL_LEFT		= SOCOR_VOL_UPLEFT;
		VOL_RIGHT		= SOCOR_VOL_UPRIGHT;
		VOL_LEFT_FLT	= SOCOR_VOL_UPLEFT_FLT;
		VOL_RIGHT_FLT	= SOCOR_VOL_UPRIGHT_FLT;
		ADJ_LONG_FLT	= SOCOR_ADJ_LONG_UP_FLT;
		ADJ_SHORT_FLT	= SOCOR_ADJ_SHORT_UP_FLT;
		KADJ_LEFT		= SOCOR_KADJ_UPLEFT;
		KADJ_RIGHT		= SOCOR_KADJ_UPRIGHT;
		KADJ_LEFT_FLT	= SOCOR_KADJ_UPLEFT_FLT;
		KADJ_RIGHT_FLT	= SOCOR_KADJ_UPRIGHT_FLT;
	}

	size_t nbFlows = spreadOption->GetNumFlows();

	ARM_Vector* Cap1Prices = spreadOption->GetCap1Prices();
	ARM_Vector* Cap2Prices = spreadOption->GetCap2Prices();
	ARM_Vector* Cap1PricesFLT = spreadOption->GetCap1PricesFLT();
	ARM_Vector* Cap2PricesFLT = spreadOption->GetCap2PricesFLT();
	ARM_Vector* longCmsTargetFLTs  = spreadOption->GetSecondFwdRateFLT();
	ARM_Vector* shortCmsTargetFLTs = spreadOption->GetFirstFwdRateFLT();
	ARM_Vector* payFwdTargetFLTs   = spreadOption->GetPayFwdRateFLT();
	
	bool varianceSqueeze = false;
	bool nullVolInputed = false;
	int periodIndex ;

	/// To save mkt target residual CSO prices
	ARM_GP_Matrix residualPrices(sizeEvalTimes,nbFlows,0.0);
	ARM_Vector& csoPrices = * spreadOption->GetCapletValues();

	for (i = 0; i < nbFlows; ++i)
	{
		double dfPay = GetNumericalModel()->GetZeroCurve()->DiscountPrice(payTimes[i]/K_YEAR_LEN);
		double div	(1.0);///  = dfPay * payNotionals[i] * payPeriods[i] * fixValues[i];
		double divFLT(1.0);/// = dfPay * payNotionals[i] * payPeriods[i] * payIndexLeverages[i]; 

		
		bool isNewPeriod = (i==0) || ( (i>0) && (periodIndexes[i]!=periodIndexes[i-1]) ) ;
		periodIndex = periodIndexes[i];
		bool isVariableCorridor = ((payIndexType == K_LIBOR || payIndexType == K_CMS)) && (payIndexLeverages[periodIndex] != 0);
		
		// We don'care of null flows
		if (fabs(payPeriods[i]) > K_NEW_DOUBLE_TOL )
		{
			// get convexified CMS rates from ARM_SpreadOption
			double longCmsTarget  = spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(i);
			double shortCmsTarget = spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(i);
			longCmsTarget  *= 0.01;
			shortCmsTarget *= 0.01;

			double longCmsTargetFLT (0.0);
			double shortCmsTargetFLT(0.0);
			double payFwdTargetFLT	(0.0);
			
			/// Float payments case
			if (isVariableCorridor)
			{				
				longCmsTargetFLT  = (*longCmsTargetFLTs)[i] ;
				shortCmsTargetFLT = (*shortCmsTargetFLTs)[i] ;
				payFwdTargetFLT   = (*payFwdTargetFLTs)[i] ;
				longCmsTargetFLT  *= 0.01;
				shortCmsTargetFLT *= 0.01;	
				payFwdTargetFLT	  *= 0.01;
			}
			
			
			/// compute approx atm spread vol
			double Tfix = resetTimes[i] / K_YEAR_LEN;
			double atmShortVol  = spreadOption->GetVol1ATM(Tfix, (shortFloatEndTimes[i] - shortFloatStartTimes[i])/K_YEAR_LEN);
			double atmLongVol   = spreadOption->GetVol2ATM(Tfix, (longFloatEndTimes[i] - longFloatStartTimes[i])/K_YEAR_LEN);
			double correl		= spreadOption->GetSpreadCorrel(Tfix, (shortFloatEndTimes[i] - shortFloatStartTimes[i])/K_YEAR_LEN, (longFloatEndTimes[i] - longFloatStartTimes[i])/K_YEAR_LEN, 0.0);
			atmShortVol			*= 0.01 * coeffsShort[i] * shortCmsTarget;
			atmLongVol			*= 0.01 * coeffsLong[i]  * longCmsTarget;
			correl				*= 0.01;
			double atmSpreadVol   = sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
			double atmSpreadStdev = atmSpreadVol * sqrt(Tfix);
			
			double forwardSpread	= coeffsLong[i] * longCmsTarget - coeffsShort[i] * shortCmsTarget;

			bool calibrateVol(true);

			if (fabs(strikes[i] - forwardSpread)> NSTDEV_NO_CALIB * atmSpreadStdev)
				calibrateVol = false;
			
			double targetSpreadVol1, targetSpreadVol2;
			double targetSpreadVol1FLT, targetSpreadVol2FLT;

			if (calibrateVol)
			{	
				double targetPrice1_N, targetPrice2_N;
				targetPrice1_N = (*Cap1Prices)[i] / div;
				targetSpreadVol1	 = VanillaImpliedVol_N (forwardSpread, targetPrice1_N, strikes[i]+spread1, Tfix, callPut, &atmSpreadVol);

				targetPrice2_N = (*Cap2Prices)[i] / div;
				targetSpreadVol2	 = VanillaImpliedVol_N (forwardSpread, targetPrice2_N, strikes[i]+spread2, Tfix, callPut, &atmSpreadVol);

				if (isVariableCorridor)
				{
					forwardSpread	= coeffsLong[i] * longCmsTargetFLT - coeffsShort[i] * shortCmsTargetFLT;
					double targetPrice1_N_FLT, targetPrice2_N_FLT;
					targetPrice1_N_FLT = (*Cap1PricesFLT)[i] / divFLT;
					targetSpreadVol1FLT   = VanillaImpliedVol_N (forwardSpread, targetPrice1_N_FLT, strikes[i]+spread1, Tfix, callPut, &atmSpreadVol);

					targetPrice2_N_FLT = (*Cap2PricesFLT)[i] / divFLT;
					targetSpreadVol2FLT   = VanillaImpliedVol_N (forwardSpread, targetPrice2_N_FLT, strikes[i]+spread2, Tfix, callPut, &atmSpreadVol);
				}
			}
		
			for (j=0; j<sizeEvalTimes; j++)
			{
			
				// calibration only if option maturity is after evalTime
				if (resetTimes[i] >= evalTimes[j] )
				{
					/// Save the current corridorlet price
					residualPrices(j,i) = csoPrices[i];

					/// compute CMS adjustements
					double longCmsModel, shortCmsModel;
					
					longCmsModel = helper.CmsRateFromNumericalModel (   evalTimes[j], 
																		payTimes[i], 
																		longFloatStartTimes[i], 
																		longFloatEndTimes[i], 
																		longSwapsFixPayTimes[i], 
																		longSwapsFixPayPeriods[i] );

					shortCmsModel = helper.CmsRateFromNumericalModel (  evalTimes[j], 
																		payTimes[i], 
																		shortFloatStartTimes[i], 
																		shortFloatEndTimes[i], 
																		shortSwapsFixPayTimes[i], 
																		shortSwapsFixPayPeriods[i] );
						

					double longAdjustment     = longCmsTarget  - longCmsModel;
					double shortAdjustment    = shortCmsTarget - shortCmsModel;
					

					// set adjusments to model params
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_LONG,  evalTimes[j], resetTimes[i], longAdjustment);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_SHORT, evalTimes[j], resetTimes[i], shortAdjustment);
				
					//
					// if variable corridor
					// o compute adjustment for the 2 CMS
					// o compute adjustment for pay rate
					if (isVariableCorridor)
					{	
						if (calibrateVol)
						{
							longCmsModel = helper.FLTconvexifiedRateFromNumericalModel(longCmsModel, 
																							evalTimes[j],
																							longFloatStartTimes[i], 
																							longFloatEndTimes[i], 
																							longSwapsFixPayTimes[i], 
																							longSwapsFixPayPeriods[i],
																							payFloatStartTimes[periodIndex], 
																							payFloatEndTimes[periodIndex], 
																							paySwapsFixPayTimes[periodIndex], 
																							paySwapsFixPayPeriods[periodIndex]);
								
							shortCmsModel = helper.FLTconvexifiedRateFromNumericalModel(shortCmsModel, 
																							evalTimes[j],
																							shortFloatStartTimes[i], 
																							shortFloatEndTimes[i], 
																							shortSwapsFixPayTimes[i], 
																							shortSwapsFixPayPeriods[i],
																							payFloatStartTimes[periodIndex], 
																							payFloatEndTimes[periodIndex], 
																							paySwapsFixPayTimes[periodIndex], 
																							paySwapsFixPayPeriods[periodIndex]);

																							
							double longAdjustmentFLT  = longCmsTargetFLT  - longCmsModel;
							double shortAdjustmentFLT = shortCmsTargetFLT - shortCmsModel;
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_LONG_FLT, evalTimes[j], resetTimes[i], longAdjustmentFLT);
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SHORT_FLT, evalTimes[j], resetTimes[i], shortAdjustmentFLT);
						}
						else
						{
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_LONG_FLT, evalTimes[j], resetTimes[i], 0.0);
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SHORT_FLT, evalTimes[j], resetTimes[i], 0.0);
						}
						
						
						/// Pay Index Ajustment : only if new period
						if (isNewPeriod)
						{
							double payFwdModelFLT = helper.CmsRateFromNumericalModel (  evalTimes[j], 
																						payTimes[i], 
																						payFloatStartTimes[periodIndex], 
																						payFloatEndTimes[periodIndex], 
																						paySwapsFixPayTimes[periodIndex], 
																						paySwapsFixPayPeriods[periodIndex] );

							double payAdjustment  = payFwdTargetFLT - payFwdModelFLT;
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_PAY_FLT, evalTimes[j], payIndexResetTimes[periodIndex], payAdjustment);
						}
					}
					
					///
					/// calibrate forward spread vol adjustment
					/// only if strike is not too far from atm
					/// if too far, a forward vol = 0 is inputed
					///
					double forwardSpreadVol1, forwardSpreadVol2;
					double forwardSpreadVol1FLT(0.0), forwardSpreadVol2FLT(0.0);
					double forwardSpreadKadj1, forwardSpreadKadj2;
					double forwardSpreadKadj1FLT(0.0), forwardSpreadKadj2FLT(0.0);
					
					if (calibrateVol)
					{
						/// compute spread forward volatility
						double modelSpreadVol1, modelSpreadVol2 ; /// model spread vol up to evalTimes[i] 

						modelSpreadVol1 = helper.SpreadNormalVolatilityFromNumericalModel (  evalTimes[j],
																							payTimes[i],
																							coeffsLong[i],
																							coeffsShort[i],
																							strikes[i]+spread1, 
																							longFloatStartTimes[i],
																							longFloatEndTimes[i], 
																							longSwapsFixPayTimes[i], 
																							longSwapsFixPayPeriods[i], 
																							shortFloatStartTimes[i], 
																							shortFloatEndTimes[i], 
																							shortSwapsFixPayTimes[i], 
																							shortSwapsFixPayPeriods[i]);

						modelSpreadVol2 = helper.SpreadNormalVolatilityFromNumericalModel ( evalTimes[j], 
																							payTimes[i], 
																							coeffsLong[i], 
																							coeffsShort[i], 
																							strikes[i]+spread2, 
																							longFloatStartTimes[i],
																							longFloatEndTimes[i], 
																							longSwapsFixPayTimes[i], 
																							longSwapsFixPayPeriods[i], 
																							shortFloatStartTimes[i], 
																							shortFloatEndTimes[i], 
																							shortSwapsFixPayTimes[i],
																							shortSwapsFixPayPeriods[i]);
						
						double Tcall = evalTimes[j] / K_YEAR_LEN;
						if (Tcall==Tfix)
							Tfix += 1./K_YEAR_LEN;
						
						double var1 = (	targetSpreadVol1 * targetSpreadVol1 * Tfix 
									  -	modelSpreadVol1  * modelSpreadVol1  * Tcall) ;

						double var2 = (	targetSpreadVol2 * targetSpreadVol2 * Tfix 
									  -	modelSpreadVol2  * modelSpreadVol2  * Tcall) ;

					

						/// everything is fine
						if (var1>=0)
						{
							forwardSpreadVol1 = sqrt( var1 / (Tfix - Tcall) );
							forwardSpreadKadj1 = 0.;
						}
						
						/// oups: variance squeeze
						else 
						{
							varianceSqueeze=true;
							forwardSpreadVol1 = 0.0;
							forwardSpreadKadj1 = 0.;

							if (GetVolUnSqueezer())
							{
								double newAdj=0.;
								if (VolUnSqueezer(coeffsLong[i] * longCmsTarget - coeffsShort[i] * shortCmsTarget, (*Cap1Prices)[i] / div, strikes[i]+spread1, Tfix, callPut,modelSpreadVol1,newAdj))
								{
									forwardSpreadKadj1 = newAdj;
									forwardSpreadVol1 = modelSpreadVol1;

								}
							}
						}

						/// everything is fine
						if (var2>=0)
						{
							forwardSpreadVol2 = sqrt( var2 / (Tfix - Tcall) );
							forwardSpreadKadj2 = 0.;
						}
						
						/// oups: variance squeeze
						else 
						{
							varianceSqueeze=true;
							forwardSpreadVol2 = 0.0;
							forwardSpreadKadj2 = 0.0;

							if (GetVolUnSqueezer())
							{
								double newAdj=0.;
								if (VolUnSqueezer(coeffsLong[i] * longCmsTarget - coeffsShort[i] * shortCmsTarget, (*Cap2Prices)[i] / div, strikes[i]+spread2, Tfix, callPut,modelSpreadVol2,newAdj))
								{
									forwardSpreadKadj2 = newAdj;
									forwardSpreadVol2 = modelSpreadVol2;
								}
							}
						}

						if (isVariableCorridor)
						{
							var1 = (	targetSpreadVol1FLT * targetSpreadVol1FLT * Tfix 
									  -	modelSpreadVol1  * modelSpreadVol1  * Tcall) ;

							var2 = (	targetSpreadVol2FLT * targetSpreadVol2FLT * Tfix 
										  -	modelSpreadVol2  * modelSpreadVol2  * Tcall) ;

						

							/// everything is fine
							if (var1>=0)
							{
								forwardSpreadVol1FLT = sqrt( var1 / (Tfix - Tcall) );
								forwardSpreadKadj1FLT = 0.;
							}
							
							/// oups: variance squeeze
							else 
							{
								varianceSqueeze=true;
								forwardSpreadVol1FLT = 0.0;
								forwardSpreadKadj1FLT = 0.;

								if (GetVolUnSqueezer())
								{
									double newAdj;
									if (VolUnSqueezer(coeffsLong[i] * longCmsTargetFLT - coeffsShort[i] * shortCmsTargetFLT, (*Cap1PricesFLT)[i] / divFLT, strikes[i]+spread1, Tfix, callPut,modelSpreadVol1,newAdj))
									{
										forwardSpreadKadj1FLT = newAdj;
										forwardSpreadVol1FLT = modelSpreadVol1;
									}
								}
							}

							/// everything is fine
							if (var2>=0)
							{
								forwardSpreadVol2FLT = sqrt( var2 / (Tfix - Tcall) );
								forwardSpreadKadj2FLT = 0.;
							}
							/// oups: variance squeeze
							else 
							{
								varianceSqueeze=true;
								forwardSpreadVol2FLT = 0.0;
								forwardSpreadKadj2FLT = 0.;

								if (GetVolUnSqueezer())
								{
									double newAdj;
									if (VolUnSqueezer(coeffsLong[i] * longCmsTargetFLT - coeffsShort[i] * shortCmsTargetFLT, (*Cap2PricesFLT)[i] / divFLT, strikes[i]+spread2, Tfix, callPut,modelSpreadVol2,newAdj))
									{
										forwardSpreadKadj2FLT = newAdj;
										forwardSpreadVol2FLT = modelSpreadVol2;
									}
								}

							}

						}
					}
					else
					{
						nullVolInputed=true;
						forwardSpreadVol1 = 0.0;
						forwardSpreadVol2 = 0.0;
						forwardSpreadVol1FLT = 0.0;
						forwardSpreadVol2FLT = 0.0;
						forwardSpreadKadj1 = 0.0;
						forwardSpreadKadj2 = 0.0;
						forwardSpreadKadj1FLT = 0.0;
						forwardSpreadKadj2FLT = 0.0;
					}

					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_LEFT, evalTimes[j], resetTimes[i], forwardSpreadVol1);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RIGHT, evalTimes[j], resetTimes[i], forwardSpreadVol2);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_LEFT, evalTimes[j], resetTimes[i], forwardSpreadKadj1);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RIGHT, evalTimes[j], resetTimes[i], forwardSpreadKadj2);

					if (isVariableCorridor)
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_LEFT_FLT, evalTimes[j], resetTimes[i], forwardSpreadVol1FLT);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RIGHT_FLT, evalTimes[j], resetTimes[i], forwardSpreadVol2FLT);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_LEFT_FLT, evalTimes[j], resetTimes[i], forwardSpreadKadj1FLT);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RIGHT_FLT, evalTimes[j], resetTimes[i], forwardSpreadKadj2FLT);
					}

/*FILE* f=fopen("c:\\temp\\dumpCCSOcalib.txt","a");
fprintf(f,"reset=%6d\tAUL=%10.7lf\tAUR=%10.7lf\tADL=%10.7lf\tADR=%10.7lf\n",
		(int)(resetTimes[i]),
		forwardSpreadKadj1,
		forwardSpreadKadj2,
		forwardSpreadKadj1FLT,
		forwardSpreadKadj2FLT);
fclose(f);*/

				}
				else
				{	
					
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_LONG, evalTimes[j], resetTimes[i], 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_SHORT, evalTimes[j], resetTimes[i], 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_LEFT, evalTimes[j], resetTimes[i], 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RIGHT, evalTimes[j], resetTimes[i], 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_LEFT, evalTimes[j], resetTimes[i], 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RIGHT, evalTimes[j], resetTimes[i], 0.0);

					if (isVariableCorridor)
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_LONG_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SHORT_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_LEFT_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RIGHT_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_LEFT_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RIGHT_FLT, evalTimes[j], resetTimes[i], 0.0);
						
						if (isNewPeriod)
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(SOCOR_ADJ_PAY_FLT, evalTimes[j], payIndexResetTimes[periodIndex], 0.0);
					}
				
				
				}
			}

		}

		if ( (i == nbFlows-1) || ( (i<nbFlows-1) && (periodIndexes[i+1] !=periodIndexes[i]) ) )
		{
			if (varianceSqueeze)
			{
				SetVarianceSqueezeStatus(true);
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Local_Normal_Model / Corridor Spread Option Calibration : variance squeeze in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
			}

			if (nullVolInputed)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Local_Normal_Model / Corridor Spread Option Calibration : strike seems to be to far from ATM. Null forward vol inputed in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
			}

			varianceSqueeze = false;
			nullVolInputed  = false;
		}
	}

	/// Save residual prices at the right place
	double price;
	size_t nbValues;
	bool isNewEvent;
	for(j=0;j<sizeEvalTimes;++j)
	{
		nbValues = static_cast<const ARM_SurfaceListModelParam&>(GetModelParams()->GetModelParam(ARM_ModelParamType::Shift)).rows(SOCOR_MKT_RESIDUAL_RSO_PRICE);
		isNewEvent = (nbValues<=j);
		for(i=0;i<nbFlows;++i)
		{
			if(resetTimes[i]>=evalTimes[j])
			{
				if(isNewEvent)
					price = 0.0;
				else
					price = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(SOCOR_MKT_RESIDUAL_RSO_PRICE,evalTimes[j],resetTimes[i]);

				if(barrierType == SOCOR_DOWN)
					price += residualPrices(j,i);
				else
					price -= residualPrices(j,i);
			}
			else
				price = 0.0;

			GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(SOCOR_MKT_RESIDUAL_RSO_PRICE,evalTimes[j],resetTimes[i],price);
		}
	}

}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateDoubleCorridorOption
///	Returns : void
///	Action  : Local model calibration to a corridor
///			  spread option wich is in fact a corridor
///			  spread option
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateDoubleCorridorOption(ARM_CorridorDblCondition* dblCorridor, 
		double targetPrice, 
		const ARM_GP_Vector& evalTimes,
		size_t barrierType)
{
	/// Static variable to avoid computing local spread & rate adjustments
	/// for corridors that differs only by the barrier profiles
	static double prevFirstResetTime = -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;

	// will generate an exception if numerical model is not supported for calibration
	Local_Normal_Model_Calibration_Helper helper (GetNumericalModel());
			
	size_t i, j, nbEvalTimes = evalTimes.size();

	// get product data
	int spreadCallPut;
	double spread1, spread2;
	ARM_GP_Vector payTimes, resetTimes, coeffsLong, coeffsShort, spreadStrikes, payPeriods, payNotionals, fixValues, payIndexLeverages;
	ARM_GP_Vector longFloatStartTimes, longFloatEndTimes, shortFloatStartTimes, shortFloatEndTimes;
	ARM_GP_Vector payFloatStartTimes, payFloatEndTimes;
	vector<ARM_GP_Vector> longSwapsFixPayTimes, longSwapsFixPayPeriods, shortSwapsFixPayTimes, shortSwapsFixPayPeriods; 
	vector<ARM_GP_Vector> paySwapsFixPayTimes, paySwapsFixPayPeriods; 
	ARM_GP_Vector payIndexResetTimes;
	ARM_IntVector periodIndexes;
	int payIndexType(0);

	/// Double condition datas
	int rateCallPut;
	ARM_GP_Vector rateStrikes,rateFloatStartTimes,rateFloatEndTimes;
	vector<ARM_GP_Vector> rateFixPayTimes,rateFixPayPeriods;
	double rateSpread1,rateSpread2;


	GetCorridorSpreadOptionData   (	
							dblCorridor, 
							spreadCallPut, /// secIdx give barrier types but may be used if conditions degenerate
							resetTimes, 
							payTimes, 
							payPeriods, 
							payNotionals, 
							coeffsLong, 
							coeffsShort, 
							spreadStrikes, 
							longFloatStartTimes,
							longFloatEndTimes,
							longSwapsFixPayTimes,
							longSwapsFixPayPeriods,
							shortFloatStartTimes,
							shortFloatEndTimes,
							shortSwapsFixPayTimes,
							shortSwapsFixPayPeriods,
							fixValues,
							spread1,
							spread2,
							payIndexType,
							payIndexLeverages,
							payFloatStartTimes, 
							payFloatEndTimes,
							paySwapsFixPayTimes,
							paySwapsFixPayPeriods,
							payIndexResetTimes,
							periodIndexes,
							rateCallPut, /// secIdx give barrier types but may be used if conditions degenerate
							rateStrikes,
							rateFloatStartTimes,
							rateFloatEndTimes,
							rateFixPayTimes,
							rateFixPayPeriods,
							rateSpread1,rateSpread2);


	size_t nbFlows = dblCorridor->GetNumFlows();

	/// Build fictive rate schedules for single rate vol model computation
	/// Only the 1st period is build for computational efficency
	ARM_GP_Vector unsedRateFloatStartTimes(rateFloatStartTimes);
	ARM_GP_Vector unsedRateFloatEndTimes(nbFlows);
	vector<ARM_GP_Vector> unsedRateFixPayTimes(nbFlows),unsedRateFixPayPeriods(nbFlows);
	for(i=0;i<nbFlows;++i)
	{
		unsedRateFloatEndTimes[i]=rateFixPayTimes[i][0];
		unsedRateFixPayTimes[i]=ARM_GP_Vector(1,rateFixPayTimes[i][0]);
		unsedRateFixPayPeriods[i]=ARM_GP_Vector(1,rateFixPayPeriods[i][0]);
	}

	size_t ADJ_SPREAD_FLT;
	size_t ADJ_RATE_FLT;

	size_t VOL_SPREAD;
	size_t VOL_RATE;
	size_t CORREL_SPREAD_RATE;

	size_t VOL_SPREAD_FLT;
	size_t VOL_RATE_FLT;
	size_t CORREL_SPREAD_RATE_FLT;

	size_t KADJ_SPREAD;
	size_t KADJ_RATE;

	size_t KADJ_SPREAD_FLT;
	size_t KADJ_RATE_FLT;

	/// 4 cases depending of lower barrier or upper barrier
	/// for CMS spread or rate indexes
	switch(barrierType)
	{
	case DBLECOR_SDOWN_RDOWN:
		{
			/// LS <= CMS Spread
			ADJ_SPREAD_FLT			= DBLECOR_SDOWN_RDOWN_ADJ_SPREAD_FLT;
			VOL_SPREAD				= DBLECOR_SDOWN_RDOWN_VOL_SPREAD;
			VOL_SPREAD_FLT			= DBLECOR_SDOWN_RDOWN_VOL_SPREAD_FLT;
			KADJ_SPREAD				= DBLECOR_SDOWN_RDOWN_KADJ_SPREAD;
			KADJ_SPREAD_FLT			= DBLECOR_SDOWN_RDOWN_KADJ_SPREAD_FLT;

			/// LR <= Rate 
			ADJ_RATE_FLT			= DBLECOR_SDOWN_RDOWN_ADJ_RATE_FLT;
			VOL_RATE				= DBLECOR_SDOWN_RDOWN_VOL_RATE;
			VOL_RATE_FLT			= DBLECOR_SDOWN_RDOWN_VOL_RATE_FLT;
			KADJ_RATE				= DBLECOR_SDOWN_RDOWN_KADJ_RATE;
			KADJ_RATE_FLT			= DBLECOR_SDOWN_RDOWN_KADJ_RATE_FLT;

			CORREL_SPREAD_RATE		= DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL;
			CORREL_SPREAD_RATE_FLT	= DBLECOR_SDOWN_RDOWN_SPREAD_RATE_CORREL_FLT;

			break;
		}
	case DBLECOR_SDOWN_RUP:
		{
			/// LS <= CMS Spread
			ADJ_SPREAD_FLT			= DBLECOR_SDOWN_RUP_ADJ_SPREAD_FLT;
			VOL_SPREAD				= DBLECOR_SDOWN_RUP_VOL_SPREAD;
			VOL_SPREAD_FLT			= DBLECOR_SDOWN_RUP_VOL_SPREAD_FLT;
			KADJ_SPREAD				= DBLECOR_SDOWN_RUP_KADJ_SPREAD;
			KADJ_SPREAD_FLT			= DBLECOR_SDOWN_RUP_KADJ_SPREAD_FLT;

			/// Rate <= UR 
			ADJ_RATE_FLT			= DBLECOR_SDOWN_RUP_ADJ_RATE_FLT;
			VOL_RATE				= DBLECOR_SDOWN_RUP_VOL_RATE;
			VOL_RATE_FLT			= DBLECOR_SDOWN_RUP_VOL_RATE_FLT;
			KADJ_RATE				= DBLECOR_SDOWN_RUP_KADJ_RATE;
			KADJ_RATE_FLT			= DBLECOR_SDOWN_RUP_KADJ_RATE_FLT;

			CORREL_SPREAD_RATE		= DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL;
			CORREL_SPREAD_RATE_FLT	= DBLECOR_SDOWN_RUP_SPREAD_RATE_CORREL_FLT;

			break;
		}
	case DBLECOR_SUP_RDOWN:
		{
			/// CMS Spread <= US
			ADJ_SPREAD_FLT			= DBLECOR_SUP_RDOWN_ADJ_SPREAD_FLT;
			VOL_SPREAD				= DBLECOR_SUP_RDOWN_VOL_SPREAD;
			VOL_SPREAD_FLT			= DBLECOR_SUP_RDOWN_VOL_SPREAD_FLT;
			KADJ_SPREAD				= DBLECOR_SUP_RDOWN_KADJ_SPREAD;
			KADJ_SPREAD_FLT			= DBLECOR_SUP_RDOWN_KADJ_SPREAD_FLT;

			/// LR <= Rate 
			ADJ_RATE_FLT			= DBLECOR_SUP_RDOWN_ADJ_RATE_FLT;
			VOL_RATE				= DBLECOR_SUP_RDOWN_VOL_RATE;
			VOL_RATE_FLT			= DBLECOR_SUP_RDOWN_VOL_RATE_FLT;
			KADJ_RATE				= DBLECOR_SUP_RDOWN_KADJ_RATE;
			KADJ_RATE_FLT			= DBLECOR_SUP_RDOWN_KADJ_RATE_FLT;

			CORREL_SPREAD_RATE		= DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL;
			CORREL_SPREAD_RATE_FLT	= DBLECOR_SUP_RDOWN_SPREAD_RATE_CORREL_FLT;

			break;
		}
	case DBLECOR_SUP_RUP:
		{
			/// CMS Spread <= US
			ADJ_SPREAD_FLT			= DBLECOR_SUP_RUP_ADJ_SPREAD_FLT;
			VOL_SPREAD				= DBLECOR_SUP_RUP_VOL_SPREAD;
			VOL_SPREAD_FLT			= DBLECOR_SUP_RUP_VOL_SPREAD_FLT;
			KADJ_SPREAD				= DBLECOR_SUP_RUP_KADJ_SPREAD;
			KADJ_SPREAD_FLT			= DBLECOR_SUP_RUP_KADJ_SPREAD_FLT;

			/// Rate <= UR 
			ADJ_RATE_FLT			= DBLECOR_SUP_RUP_ADJ_RATE_FLT;
			VOL_RATE				= DBLECOR_SUP_RUP_VOL_RATE;
			VOL_RATE_FLT			= DBLECOR_SUP_RUP_VOL_RATE_FLT;
			KADJ_RATE				= DBLECOR_SUP_RUP_KADJ_RATE;
			KADJ_RATE_FLT			= DBLECOR_SUP_RUP_KADJ_RATE_FLT;

			CORREL_SPREAD_RATE		= DBLECOR_SUP_RUP_SPREAD_RATE_CORREL;
			CORREL_SPREAD_RATE_FLT	= DBLECOR_SUP_RUP_SPREAD_RATE_CORREL_FLT;

			break;
		}
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, " : unknown spread/rate barrier type");

	}

	/// Restore convexified forward rates of the underlying corridor
	const ARM_Vector& longCmsTargets	= * dblCorridor->GetSpreadLeg()->GetSecondLeg()->GetFwdRates();
	const ARM_Vector& shortCmsTargets	= * dblCorridor->GetSpreadLeg()->GetFirstLeg()->GetFwdRates();
	const ARM_Vector& longCmsTargetsFLT	= * dblCorridor->GetSecondFwdRateFLT();
	const ARM_Vector& shortCmsTargetsFLT= * dblCorridor->GetFirstFwdRateFLT();

	const ARM_Vector& rateTargetsFLT	= * dblCorridor->GetSpreadDigital()->GetSecondFwdRateFLT();
	const ARM_Vector& rateTargets		= * dblCorridor->GetSpreadDigital()->GetSpreadLeg()->GetSecondLeg()->GetFwdRates();

	const ARM_Vector& payFwdTargetsFLT	= * dblCorridor->GetPayFwdRateFLT();

	/// Restore adjusted strikes, vols and correls of the underlying corridor
	const ARM_Vector& spreadKAdj		= dblCorridor->GetSpreadKAdj();
	const ARM_Vector& spreadVolTarget	= dblCorridor->GetSpreadVol();
	const ARM_Vector& rateKAdj			= dblCorridor->GetDigitalKAdj();
	const ARM_Vector& rateVolTarget		= dblCorridor->GetDigitalVol();

	const ARM_Vector& spreadKAdjFLT			= dblCorridor->GetSpreadKAdjFLT();
	const ARM_Vector& spreadVolTargetFLT	= dblCorridor->GetSpreadVolFLT();
	const ARM_Vector& rateKAdjFLT			= dblCorridor->GetDigitalKAdjFLT();
	const ARM_Vector& rateVolTargetFLT		= dblCorridor->GetDigitalVolFLT();

	const ARM_Vector& spreadRateCorrelTarget	= dblCorridor->GetSpreadDigitalCorrel();

	/// To save mkt target residual RA2 prices
	ARM_GP_Matrix residualPrices(nbEvalTimes,nbFlows,0.0);
	ARM_Vector& dblCorridorPrices		= * dblCorridor->GetCapletValues();


	bool isNewPeriod,isVariableCorridor,isCalibrateSpreadVol,isCalibrateRateVol;
	bool isSpreadOut,isRateOut,isFixedPart;
	bool spreadSqueeze=false,rateSqueeze=false,correlSqueeze=false;
	bool spreadLocalSqueeze,rateLocalSqueeze,spreadLocalSqueezeFLT,rateLocalSqueezeFLT,isVolUnsqueezer;
	bool nullSpreadVolInput=false,nullRateVolInput=false,nullSpreadRateCorrelInput=false;
	double maxCorrelSqueeze = -ARM_NumericConstants::ARM_BIGGEST_POSITIVE_NUMBER;

	/// Allow spread/rate local adjustments only if it is a brand new corridor !
	bool isNewCorridor = (prevFirstResetTime != resetTimes[0]);
	bool isCalibrateSpreadRateAdj=true;
	if(GetResetCalib())
		SetResetCalib(false);
	else
		isCalibrateSpreadRateAdj = isNewCorridor;

	prevFirstResetTime = resetTimes[0];

	int periodIndex ;

	double dfPay,Tfix,sTfix,Tcall,sTfixTcall;
	double longCmsTarget,shortCmsTarget,spreadTarget,rateTarget;
	double longCmsTargetFLT,shortCmsTargetFLT,spreadTargetFLT,rateTargetFLT,payFwdTargetFLT;
	double longCmsModel,shortCmsModel,spreadModel,rateModel;
	double spreadVolModel,rateVolModel,spreadVarModel,rateVarModel;
	double localSpreadVar,localSpreadVol,localSpreadVolFLT;
	double localRateVar,localRateVol,localRateVolFLT;
	double spreadRateCovModel,localSpreadRateCov,localSpreadRateCorrel,localSpreadRateCorrelFLT;
	size_t iter;

	double newSpreadKAdj,newRateKAdj,newSpreadStdDev,newRateStdDev;
	double newLocalSpreadStdDev,newLocalRateStdDev;
	double newSpreadKAdjFLT,newRateKAdjFLT,newSpreadStdDevFLT,newRateStdDevFLT;
	double newLocalSpreadStdDevFLT,newLocalRateStdDevFLT;
	double residualValue;

/****
FILE* f2=fopen("c:\\temp\\dumpCRA2Ref.txt","a");
double xSRef0,xRRef0;
for(i = 0; i < nbFlows; ++i)
{
	/// Double condition reset time
	double Tfix = resetTimes[i] / K_YEAR_LEN;
	double sTfix = sqrt(Tfix);

	periodIndex			= periodIndexes[i];
	isVariableCorridor	= ((payIndexType == K_LIBOR || payIndexType == K_CMS)) && (payIndexLeverages[periodIndex] != 0);
	newSpreadStdDev		= spreadVolTarget[i] * sTfix;
	newRateStdDev		= rateVolTarget[i] * sTfix;
	longCmsTarget		= 0.01*longCmsTargets[i];
	shortCmsTarget		= 0.01*shortCmsTargets[i];
	spreadTarget		= coeffsLong[i] * longCmsTarget - coeffsShort[i] * shortCmsTarget;
	rateTarget			= 0.01*rateTargets[i];

	if(newSpreadStdDev!=0.0)
		xSRef0 = (spreadKAdj[i]-spreadTarget)/newSpreadStdDev;
	else
		xSRef0 = spreadKAdj[i] > spreadTarget ? ARM_NumericConstants::ARM_INFINITY : -ARM_NumericConstants::ARM_INFINITY;
	double NSRef0=NormalCDF(xSRef0);
	if(newRateStdDev!=0.0)
		xRRef0 = (rateKAdj[i]-rateTarget)/newRateStdDev;
	else
		xRRef0 = rateKAdj[i] > rateTarget ? ARM_NumericConstants::ARM_INFINITY : -ARM_NumericConstants::ARM_INFINITY;
	double NRRef0=NormalCDF(xRRef0);
	double refPrice0=NormalCDF(xRRef0,xSRef0,spreadRateCorrelTarget[i]);
	refPrice0 += 1 - NRRef0 - NSRef0;
	fprintf(f2,"Tfix=\t%5d\tSK=\t%10.7lf\tRK=\t%10.7lf\tAdjSK=\t%10.7lf\tAdjRK=\t%10.7lf\tSVol=\t%10.7lf\tRVol=\t%10.7lf\trefPrice=\t%15.10lf\n",
			(int)(Tfix*365),spreadStrikes[i]*100,rateStrikes[i]*100,spreadKAdj[i]*100,rateKAdj[i]*100,
			spreadVolTarget[i]*100,rateVolTarget[i]*100,refPrice0);
}
fclose(f2);
****/

	for(i = 0; i < nbFlows; ++i)
	{
		dfPay = GetNumericalModel()->GetZeroCurve()->DiscountPrice(payTimes[i]/K_YEAR_LEN);
		
		periodIndex = periodIndexes[i];

		isNewPeriod = (i==0) || ( (i>0) && (periodIndexes[i]!=periodIndexes[i-1]) ) ;
		isVariableCorridor = ((payIndexType == K_LIBOR || payIndexType == K_CMS)) && (payIndexLeverages[periodIndex] != 0);
		
		// We don't care of null flows
		if (fabs(payPeriods[i]) > K_NEW_DOUBLE_TOL )
		{
			// Get convexified rates 
			longCmsTarget	= 0.01*longCmsTargets[i];
			shortCmsTarget	= 0.01*shortCmsTargets[i];
			spreadTarget	= coeffsLong[i] * longCmsTarget - coeffsShort[i] * shortCmsTarget;

			rateTarget		= 0.01*rateTargets[i];
			
			/// Case of floating payment rates
			if(isVariableCorridor)
			{				
				longCmsTargetFLT	= 0.01*longCmsTargetsFLT[i];
				shortCmsTargetFLT	= 0.01*shortCmsTargetsFLT[i];
				spreadTargetFLT		= coeffsLong[i] * longCmsTargetFLT - coeffsShort[i] * shortCmsTargetFLT;

				payFwdTargetFLT		= 0.01*payFwdTargetsFLT[i];
				rateTargetFLT		= 0.01*rateTargetsFLT[i];
			}
			
			
			/// Check if market vols will to be calibrated
			isSpreadOut = fabs(spreadStrikes[i] - spreadTarget) >= NSTDEV_NO_CALIB * spreadVolTarget[i];
			if( (!isVariableCorridor && isSpreadOut) ||
				( isVariableCorridor && isSpreadOut && fabs(spreadStrikes[i] - spreadTargetFLT) >= NSTDEV_NO_CALIB * spreadVolTargetFLT[i]) )
				isCalibrateSpreadVol = false;
			else
				isCalibrateSpreadVol = true;

			isRateOut = fabs(rateStrikes[i] - rateTarget)>= NSTDEV_NO_CALIB * rateVolTarget[i];
			if( (!isVariableCorridor && isRateOut) ||
				( isVariableCorridor && isRateOut && fabs(rateStrikes[i] - rateTargetFLT)>= NSTDEV_NO_CALIB * rateVolTargetFLT[i]) )
				isCalibrateRateVol = false;
			else
				isCalibrateRateVol = true;

			isFixedPart = fixValues[i] != 0.0; /// fixed rate or margin of floating rate paid
			
			/// Double condition reset time
			Tfix = resetTimes[i] / K_YEAR_LEN;
			sTfix = sqrt(Tfix);

			for(j=0; j<nbEvalTimes; j++)
			{
				Tcall = evalTimes[j] / K_YEAR_LEN;
				if (Tcall==Tfix)
					Tfix += 1./K_YEAR_LEN; // to allow artificial VCV local calibration

				sTfixTcall = sqrt(Tfix-Tcall);

				// Local model calibrations only if option maturity is after evalTime
				if(resetTimes[i] >= evalTimes[j] )
				{
					/// Save the current corridorlet price
					residualPrices(j,i) = dblCorridorPrices[i];

					newSpreadStdDev		= spreadVolTarget[i] * sTfix;
					newSpreadStdDevFLT	= isVariableCorridor ? spreadVolTargetFLT[i] * sTfix : 0.0;

					newRateStdDev		= rateVolTarget[i] * sTfix;
					newRateStdDevFLT	= isVariableCorridor ? rateVolTargetFLT[i] * sTfix : 0.0;

/****
double xSRef,xRRef;
if(newSpreadStdDev!=0.0)
	xSRef = (spreadKAdj[i]-spreadTarget)/newSpreadStdDev;
else
xSRef = spreadKAdj[i] > spreadTarget ? ARM_NumericConstants::ARM_INFINITY : -ARM_NumericConstants::ARM_INFINITY;
double NSRef=NormalCDF(xSRef);
if(newRateStdDev!=0.0)
	xRRef = (rateKAdj[i]-rateTarget)/newRateStdDev;
else
	xRRef = rateKAdj[i] > rateTarget ? ARM_NumericConstants::ARM_INFINITY : -ARM_NumericConstants::ARM_INFINITY;
double NRRef=NormalCDF(xRRef);
double refPrice=NormalCDF(xRRef,xSRef,spreadRateCorrelTarget[i]);
refPrice += 1 - NRRef - NSRef;
FILE* f0=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f0,"Tcall=\t%5d\tTfix=\t%5d\tSK=\t%10.7lf\tRK=\t%10.7lf\tAdjSK=\t%10.7lf\tAdjRK=\t%10.7lf\tSVol=\t%10.7lf\tRVol=\t%10.7lf\trefPrice=\t%15.10lf\tSVar=\t%15.10lf\tRVar=\t%15.10lf\tCov=\t%15.10lf\n",
		(int)(Tcall*365),(int)(Tfix*365),
		spreadStrikes[i]*100,rateStrikes[i]*100,spreadKAdj[i]*100,rateKAdj[i]*100,
		newSpreadStdDev/sTfix*100,newRateStdDev/sTfix*100,refPrice,
		newSpreadStdDev*newSpreadStdDev,newRateStdDev*newRateStdDev,
		newSpreadStdDev*newRateStdDev*spreadRateCorrelTarget[i]);
fclose(f0);
****/

					if( isVariableCorridor || (!isVariableCorridor && isCalibrateSpreadRateAdj) )
					{
						/// Compute forward spread model adjustement					
						longCmsModel = helper.CmsRateFromNumericalModel (   evalTimes[j], 
																			payTimes[i], 
																			longFloatStartTimes[i], 
																			longFloatEndTimes[i], 
																			longSwapsFixPayTimes[i], 
																			longSwapsFixPayPeriods[i] );

						shortCmsModel = helper.CmsRateFromNumericalModel (  evalTimes[j], 
																			payTimes[i], 
																			shortFloatStartTimes[i], 
																			shortFloatEndTimes[i], 
																			shortSwapsFixPayTimes[i], 
																			shortSwapsFixPayPeriods[i] );
							


						/// Compute forward rate model adjustment
						rateModel = helper.CmsRateFromNumericalModel (  evalTimes[j], 
																		payTimes[i], 
																		rateFloatStartTimes[i], 
																		rateFloatEndTimes[i], 
																		rateFixPayTimes[i], 
																		rateFixPayPeriods[i] );
					}
						

					if(isCalibrateSpreadRateAdj)
					{
						/// Set adjusments to model params only for brand new corridor else saved the first time
						spreadModel = coeffsLong[i] * longCmsModel - coeffsShort[i] * shortCmsModel;
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_SPREAD,evalTimes[j],resetTimes[i],spreadTarget-spreadModel);
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_RATE,evalTimes[j],resetTimes[i],rateTarget-rateModel);

						if(isVariableCorridor)
						{
							if(isNewPeriod)
							{
								/// Compute and set forward payment rate model adjustment only for brand new corridor and new paid flow
								double payFwdModel = helper.CmsRateFromNumericalModel(	evalTimes[j], 
																						payTimes[i], 
																						payFloatStartTimes[periodIndex], 
																						payFloatEndTimes[periodIndex], 
																						paySwapsFixPayTimes[periodIndex], 
																						paySwapsFixPayPeriods[periodIndex] );

								GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_PAY_FLT, evalTimes[j], payIndexResetTimes[periodIndex],payFwdTargetFLT - payFwdModel);
							}
						}
						else
							GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_PAY_FLT, evalTimes[j], payIndexResetTimes[periodIndex], 0.0);
					}

					if(isVariableCorridor)
					{
						/// Here spread/rate adjustement depends of the current barriers

						/// Compute forward spread model adjustement					
						longCmsModel = helper.FLTconvexifiedRateFromNumericalModel(longCmsModel, 
																					evalTimes[j],
																					longFloatStartTimes[i], 
																					longFloatEndTimes[i], 
																					longSwapsFixPayTimes[i], 
																					longSwapsFixPayPeriods[i],
																					payFloatStartTimes[periodIndex], 
																					payFloatEndTimes[periodIndex], 
																					paySwapsFixPayTimes[periodIndex], 
																					paySwapsFixPayPeriods[periodIndex]);
							
						shortCmsModel = helper.FLTconvexifiedRateFromNumericalModel(shortCmsModel, 
																					evalTimes[j],
																					shortFloatStartTimes[i], 
																					shortFloatEndTimes[i], 
																					shortSwapsFixPayTimes[i], 
																					shortSwapsFixPayPeriods[i],
																					payFloatStartTimes[periodIndex], 
																					payFloatEndTimes[periodIndex], 
																					paySwapsFixPayTimes[periodIndex], 
																					paySwapsFixPayPeriods[periodIndex]);

						spreadModel = coeffsLong[i] * longCmsModel - coeffsShort[i] * shortCmsModel;

						/// Compute forward rate model adjustment
						rateModel = helper.FLTconvexifiedRateFromNumericalModel(rateModel, 
																				evalTimes[j],
																				rateFloatStartTimes[i], 
																				rateFloatEndTimes[i], 
																				rateFixPayTimes[i], 
																				rateFixPayPeriods[i],
																				payFloatStartTimes[periodIndex], 
																				payFloatEndTimes[periodIndex], 
																				paySwapsFixPayTimes[periodIndex], 
																				paySwapsFixPayPeriods[periodIndex]);

						/// Set convexity adjustments
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],spreadTargetFLT-spreadModel);
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_RATE_FLT,evalTimes[j],resetTimes[i],rateTargetFLT-rateModel);
					}
					else
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SPREAD_FLT, evalTimes[j], resetTimes[i], 0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_RATE_FLT, evalTimes[j], resetTimes[i], 0.0);
					}


					/// Save model strike adjustements for spread & rate digitals
					/// ---------------------------------------------------------

					if(isFixedPart)
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i],spreadKAdj[i]-spreadStrikes[i]);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE,evalTimes[j],resetTimes[i],rateKAdj[i]-rateStrikes[i]);
					}
					else
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i],0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE,evalTimes[j],resetTimes[i],0.0);
					}
					
					if(isVariableCorridor)
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],spreadKAdjFLT[i]-spreadStrikes[i]);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE_FLT,evalTimes[j],resetTimes[i],rateKAdjFLT[i]-rateStrikes[i]);
					}
					else
					{
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],0.0);
						GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE_FLT,evalTimes[j],resetTimes[i],0.0);
					}



					/// Compute local volatility of the spread (fixed or floating
					/// rate payment) to fulfill market variance target
					/// ---------------------------------------------------------
					spreadVarModel			= 0.0;
					localSpreadVol			= 0.0;
					localSpreadVolFLT		= 0.0;
					spreadLocalSqueeze		= false;
					spreadLocalSqueezeFLT	= false;
					if(isCalibrateSpreadVol)
					{
						spreadVolModel = helper.SpreadNormalVolatilityFromNumericalModel(	evalTimes[j],
																							payTimes[i],
																							coeffsLong[i],
																							coeffsShort[i],
																							spreadStrikes[i], 
																							longFloatStartTimes[i],
																							longFloatEndTimes[i], 
																							longSwapsFixPayTimes[i], 
																							longSwapsFixPayPeriods[i], 
																							shortFloatStartTimes[i], 
																							shortFloatEndTimes[i], 
																							shortSwapsFixPayTimes[i], 
																							shortSwapsFixPayPeriods[i]);
						spreadVarModel = spreadVolModel  * spreadVolModel  * Tcall;

						if(isFixedPart)
						{
							localSpreadVar = newSpreadStdDev * newSpreadStdDev - spreadVarModel;

							if(localSpreadVar>0)
								localSpreadVol = sqrt( localSpreadVar / (Tfix - Tcall) );
							else // Variance squeeze
								spreadLocalSqueeze=true;
						}

						if(isVariableCorridor)
						{
							/// Model spread vol not recomputed (no dependency with our implicit H&W model !)
							localSpreadVar = newSpreadStdDevFLT * newSpreadStdDevFLT - spreadVarModel;
							if(localSpreadVar>0.0)
								localSpreadVolFLT = sqrt( localSpreadVar / (Tfix - Tcall) );
							else // Variance squeeze
								spreadLocalSqueezeFLT=true;
						}
					}
					else
						nullSpreadVolInput=true;

					spreadSqueeze = spreadSqueeze|| spreadLocalSqueeze || spreadLocalSqueezeFLT;

					/// Compute local volatility of the rate (fixed or floating
					/// rate payment) to fulfill market variance target
					/// -------------------------------------------------------
					rateVarModel		= 0.0;
					localRateVol		= 0.0;
					localRateVolFLT		= 0.0;
					rateLocalSqueeze	= false;
					rateLocalSqueezeFLT	= false;
					if(isCalibrateRateVol)
					{
						/// Same spread vol function but second rate is unused !
						rateVolModel = helper.SpreadNormalVolatilityFromNumericalModel(	evalTimes[j], 
																						payTimes[i], 
																						1.0, 
																						0.0, 
																						rateStrikes[i], 
																						rateFloatStartTimes[i],
																						rateFloatEndTimes[i], 
																						rateFixPayTimes[i], 
																						rateFixPayPeriods[i], 
																						unsedRateFloatStartTimes[i], 
																						unsedRateFloatEndTimes[i], 
																						unsedRateFixPayTimes[i],
																						unsedRateFixPayPeriods[i]);
						rateVarModel = rateVolModel  * rateVolModel  * Tcall;

						if(isFixedPart)
						{
							localRateVar = newRateStdDev * newRateStdDev - rateVarModel;

							localRateVol=0.0;
							if(localRateVar>0.0)
								localRateVol = sqrt( localRateVar / (Tfix - Tcall) );
							else // Variance squeeze
								rateLocalSqueeze=true;
						}

						if(isVariableCorridor)
						{
							/// Model rate vol not recomputed (no dependency with our implicit H&W model !)
							localRateVar = newRateStdDevFLT * newRateStdDevFLT - rateVarModel;
							if (localRateVar>0)
								localRateVolFLT = sqrt( localRateVar / (Tfix - Tcall) );
							else // Variance squeeze
								rateLocalSqueezeFLT=true;
						}
					}
					else
						nullRateVolInput=true;
;
					rateSqueeze = rateSqueeze || rateLocalSqueeze || rateLocalSqueezeFLT;

					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD,evalTimes[j],resetTimes[i],localSpreadVol);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE,evalTimes[j],resetTimes[i],localRateVol);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD_FLT,evalTimes[j],resetTimes[i],localSpreadVolFLT);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE_FLT,evalTimes[j],resetTimes[i],localRateVolFLT);

/****
f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"SVolMod=%15.10lf\tRVolMod=%15.10lf\t",spreadVolModel,rateVolModel);
fclose(f);
****/

					spreadRateCovModel		= 0.0;
					localSpreadRateCorrel	= 0.0;
					localSpreadRateCorrelFLT= 0.0;
					if(isCalibrateSpreadVol && isCalibrateRateVol)
					{
						/// Compute local correlation between spread and rate (fixed/margin or
						/// floating rate payment) to fulfill market covariance target
						/// ------------------------------------------------------------------
						spreadRateCovModel = helper.SpreadRateNormalCovarianceFromNumericalModel(
							evalTimes[j],
							payTimes[i],
							coeffsLong[i],
							coeffsShort[i],
							spreadStrikes[i],
							longFloatStartTimes[i],
							longFloatEndTimes[i], 
							longSwapsFixPayTimes[i], 
							longSwapsFixPayPeriods[i], 
							shortFloatStartTimes[i], 
							shortFloatEndTimes[i], 
							shortSwapsFixPayTimes[i], 
							shortSwapsFixPayPeriods[i],
							rateStrikes[i],
							rateFloatStartTimes[i], 
							rateFloatEndTimes[i], 
							rateFixPayTimes[i], 
							rateFixPayPeriods[i]);

/****
f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"SRCorMod=%15.10lf\t",spreadRateCovModel/(spreadVolModel*rateVolModel*Tcall));
fclose(f);
****/

						if(isFixedPart)
						{
							localSpreadRateCov = spreadRateCorrelTarget[i]*spreadVolTarget[i]*rateVolTarget[i]*Tfix - spreadRateCovModel;
							isVolUnsqueezer = GetVolUnSqueezer() && (spreadLocalSqueeze || rateLocalSqueeze);
							if((localSpreadVol > 0.0 && localRateVol > 0.0) || isVolUnsqueezer)
							{
								if(localSpreadVol > 0.0 && localRateVol > 0.0)
									localSpreadRateCorrel = localSpreadRateCov / (localSpreadVol*localRateVol*(Tfix - Tcall));
								else
									localSpreadRateCorrel = localSpreadRateCov > 0 ? 101 : -101;

								double localSpreadRateCorrelInit = localSpreadRateCorrel;

								if(fabs(localSpreadRateCorrel) > 1)
								{
									correlSqueeze=true;
									residualValue = (localSpreadRateCorrel>0 ? localSpreadRateCorrel-1 : localSpreadRateCorrel+1);
									maxCorrelSqueeze = (residualValue > maxCorrelSqueeze ? residualValue : maxCorrelSqueeze);

									/// Covariance squeeze
									newSpreadKAdj	= spreadKAdj[i];
									newRateKAdj		= rateKAdj[i];
									if( GetCorrelUnSqueezer() || isVolUnsqueezer )
									{
										/// Try to fit variance/covariance by modifying
										/// reference vols and consequently adjusted strikes
										if( CorrelUnSqueezer(Tcall,Tfix,spreadTarget,rateTarget,spreadRateCorrelTarget[i],
											spreadVarModel,rateVarModel,spreadRateCovModel,sTfixTcall,
											newSpreadKAdj,newRateKAdj,newSpreadStdDev,newRateStdDev,
											newLocalSpreadStdDev,newLocalRateStdDev,localSpreadRateCorrel,iter))
										{

/***
FILE* f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"Call=%6.1lf\t,Fix=%6.1lf\t => locCorInit=%10.7lf\tlocCorUnsq=%10.7lf\titer=%3d\n",
		Tcall*365,Tfix*365,localSpreadRateCorrelInit,localSpreadRateCorrel,iter);
fclose(f);
***/
										}
										else
										{
/***
FILE* f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"Call=%6.1lf\t,Fix=%6.1lf\t => locCorInit=%10.7lf\tUnable to unsqueeze : locCor=%10.7lf\titer=%3d\t\t",
		Tcall*365,Tfix*365,localSpreadRateCorrelInit,localSpreadRateCorrel,iter);
fprintf(f,"SVarMod=%15.10lf\tRVarMod=%15.10lf\tSRCovMod=%15.10lf\tSStdTarg=%15.10lf\tRStdTarg=%15.10lf\tSRCorTarg=%15.10lf\tS=%15.10lf\tR=%15.10lf\tKS=%15.10lf\tKR=%15.10lf\n",
		spreadVarModel,rateVarModel,spreadRateCovModel,
		spreadVolTarget[i] * sTfix,rateVolTarget[i] * sTfix,spreadRateCorrelTarget[i],
		spreadTarget,rateTarget,spreadKAdj[i],rateKAdj[i]);
fclose(f);
***/
										}
										if( (localSpreadRateCorrelInit > 1 && -1 < localSpreadRateCorrel && localSpreadRateCorrel<localSpreadRateCorrelInit) ||
											(localSpreadRateCorrelInit < -1 && localSpreadRateCorrelInit < localSpreadRateCorrel && localSpreadRateCorrel < 1) )
										{
											/// Better solution => change local saved values : KAdj & localVol
											GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i],newSpreadKAdj-spreadStrikes[i]);
											GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE,evalTimes[j],resetTimes[i],newRateKAdj-rateStrikes[i]);

											GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD,evalTimes[j],resetTimes[i],newLocalSpreadStdDev/sTfixTcall);
											GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE,evalTimes[j],resetTimes[i],newLocalRateStdDev/sTfixTcall);
										}
									}
									else
									{
/***
FILE* f=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f,"Call=%6.1lf\t,Fix=%6.1lf\t => locCorInit=%10.7lf\tCorrel unsqueeze not activated : locCor=%10.7lf\t\t",
		Tcall*365,Tfix*365,localSpreadRateCorrelInit,localSpreadRateCorrel);
fprintf(f,"SVarMod=%15.10lf\tRVarMod=%15.10lf\tSRCovMod=%15.10lf\tSStdTarg=%15.10lf\tRStdTarg=%15.10lf\tSRCorTarg=%15.10lf\tS=%15.10lf\tR=%15.10lf\tKS=%15.10lf\tKR=%15.10lf\n",
		spreadVarModel,rateVarModel,spreadRateCovModel,
		spreadVolTarget[i] * sTfix,rateVolTarget[i] * sTfix,spreadRateCorrelTarget[i],
		spreadTarget,rateTarget,spreadKAdj[i],rateKAdj[i]);
fclose(f);
***/
									}
									localSpreadRateCorrel = (localSpreadRateCorrel > 1 ? 0.99 : -0.99);
								}
							}
							else
								localSpreadRateCorrel = (localSpreadRateCov > 0 ? 1 : -1);
						} // if fixed rate or margin payment


						if(isVariableCorridor)
						{
							localSpreadRateCov = spreadRateCorrelTarget[i]*spreadVolTargetFLT[i]*rateVolTargetFLT[i]*Tfix - spreadRateCovModel;
							isVolUnsqueezer = GetVolUnSqueezer() && (spreadLocalSqueezeFLT || rateLocalSqueezeFLT);
							if((localSpreadVolFLT > 0.0 && localRateVolFLT > 0.0) || isVolUnsqueezer)
							{
								if(localSpreadVolFLT > 0.0 && localRateVolFLT > 0.0)
									localSpreadRateCorrelFLT = localSpreadRateCov / (localSpreadVolFLT*localRateVolFLT*(Tfix - Tcall));
								else
									localSpreadRateCorrelFLT = localSpreadRateCov > 0 ? 101 : -101;

								double localSpreadRateCorrelInitFLT = localSpreadRateCorrelFLT;

								if(fabs(localSpreadRateCorrelFLT) > 1)
								{
									/// Covariance squeeze
									correlSqueeze=true;
									residualValue = (localSpreadRateCorrel>0 ? localSpreadRateCorrel-1 : localSpreadRateCorrelFLT+1);
									maxCorrelSqueeze = (residualValue > maxCorrelSqueeze ? residualValue : maxCorrelSqueeze);

									/// Covariance squeeze
									if(GetCorrelUnSqueezer() || isVolUnsqueezer)
									{
										/// Try to fit variance/covariance by modifying
										/// reference vols and consequently adjusted strikes
										newSpreadKAdjFLT	= spreadKAdj[i];
										newRateKAdjFLT		= rateKAdj[i];
										CorrelUnSqueezer(Tcall,Tfix,spreadTargetFLT,rateTargetFLT,spreadRateCorrelTarget[i],
											spreadVarModel,rateVarModel,spreadRateCovModel,sTfixTcall,
											newSpreadKAdjFLT,newRateKAdjFLT,newSpreadStdDevFLT,newRateStdDevFLT,
											newLocalSpreadStdDevFLT,newLocalRateStdDevFLT,localSpreadRateCorrelFLT,iter);

										if( (localSpreadRateCorrelInitFLT > 1 && -1 < localSpreadRateCorrelFLT && localSpreadRateCorrelFLT<localSpreadRateCorrelInitFLT) ||
											(localSpreadRateCorrelInitFLT < -1 && localSpreadRateCorrelInitFLT < localSpreadRateCorrelFLT && localSpreadRateCorrelFLT < 1) )
										{
											/// Change local saved values : KAdj & localVol
											GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],newSpreadKAdjFLT-spreadStrikes[i]);
											GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE_FLT,evalTimes[j],resetTimes[i],newRateKAdjFLT-rateStrikes[i]);

											GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD_FLT,evalTimes[j],resetTimes[i],newLocalSpreadStdDevFLT/sTfixTcall);
											GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE_FLT,evalTimes[j],resetTimes[i],newLocalRateStdDevFLT/sTfixTcall);
										}
									}
									localSpreadRateCorrelFLT = (localSpreadRateCorrelFLT>0 ? 0.99 : -0.99);
								}
							}
							else
								localSpreadRateCorrelFLT = (localSpreadRateCov > 0 ? 1 : -1);
						} // if variable payment
					} // if spread & rate vol calibration
					else
					{
						nullSpreadRateCorrelInput=true;

						if(GetVolUnSqueezer())
						{
							/// Unsqueeze volatility by moving strikes so that
							/// local vols are just set to the minimum value =>
							/// VarTarget(0->Fixing) = VarModel(0->Notice) + VarMin(Notice->Fixing)
							double xS,xR;

							double localSpreadVarMin = LOCAL_SPREAD_VOL_MIN*sTfixTcall;
							localSpreadVarMin *= localSpreadVarMin;
							if(spreadLocalSqueeze && newSpreadStdDev != 0.0)
							{
								xS = (spreadTarget-spreadKAdj[i])/newSpreadStdDev;
								newSpreadStdDev = sqrt(spreadVarModel+localSpreadVarMin);
								newSpreadKAdj = spreadTarget - xS * newSpreadStdDev;
								GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i],newSpreadKAdj-spreadStrikes[i]);
								GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD,evalTimes[j],resetTimes[i],LOCAL_SPREAD_VOL_MIN);
							}
							if(spreadLocalSqueezeFLT && newSpreadStdDevFLT != 0.0)
							{
								xS = (spreadTarget-spreadKAdjFLT[i])/newSpreadStdDevFLT;
								newSpreadStdDevFLT = sqrt(spreadVarModel+localSpreadVarMin);
								newSpreadKAdjFLT = spreadTargetFLT - xS * newSpreadStdDevFLT;
								GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],newSpreadKAdjFLT-spreadStrikes[i]);
								GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD_FLT,evalTimes[j],resetTimes[i],LOCAL_SPREAD_VOL_MIN);
							}

							double localRateVarMin = LOCAL_RATE_VOL_MIN*sTfixTcall;
							localRateVarMin *= localRateVarMin;
							if(rateLocalSqueeze && newRateStdDev != 0.0)
							{
								xR = (rateTarget-rateKAdj[i])/newRateStdDev;
								newRateStdDev = sqrt(rateVarModel+localRateVarMin);
								newRateKAdj = rateTarget - xR * newRateStdDev;
								GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE,evalTimes[j],resetTimes[i],newRateKAdj-rateStrikes[i]);
								GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE,evalTimes[j],resetTimes[i],LOCAL_RATE_VOL_MIN);
							}
							if(rateLocalSqueezeFLT && newRateStdDevFLT != 0.0)
							{
								xR = (rateTargetFLT-rateKAdjFLT[i])/newRateStdDevFLT;
								newRateStdDevFLT = sqrt(rateVarModel+localRateVarMin);
								newRateKAdjFLT = rateTargetFLT - xR * newRateStdDevFLT;
								GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE_FLT,evalTimes[j],resetTimes[i],newRateKAdjFLT-rateStrikes[i]);
								GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE_FLT,evalTimes[j],resetTimes[i],LOCAL_RATE_VOL_MIN);
							}
						}
					}

					GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(CORREL_SPREAD_RATE,evalTimes[j],resetTimes[i],localSpreadRateCorrel);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(CORREL_SPREAD_RATE_FLT,evalTimes[j],resetTimes[i],localSpreadRateCorrelFLT);

/****
double SKAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i]);
double RKAdj = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(KADJ_RATE,evalTimes[j],resetTimes[i]);
double SVolNew	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VOL_SPREAD,evalTimes[j],resetTimes[i]);
double RVolNew	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(VOL_RATE,evalTimes[j],resetTimes[i]);
double SRCorNew	= GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(CORREL_SPREAD_RATE,evalTimes[j],resetTimes[i]);

double xSNew = (spreadStrikes[i]+SKAdj-spreadTarget)/sqrt(spreadVarModel+SVolNew*SVolNew*(Tfix-Tcall));
double NSNew=NormalCDF(xSNew);
double xRNew = (rateStrikes[i]+RKAdj-rateTarget)/sqrt(rateVarModel+RVolNew*RVolNew*(Tfix-Tcall));
double NRNew=NormalCDF(xRNew);

double refPriceNew=NormalCDF(xRNew,xSNew,spreadRateCorrelTarget[i]);
refPriceNew += 1 - NRNew - NSNew;
FILE* f1=fopen("c:\\temp\\dumpCRA2.txt","a");
fprintf(f1,"      \t%5d\t     \t%5d\tSK=\t%10.7lf\tRK=\t%10.7lf\tAdjSK=\t%10.7lf\tAdjRK=\t%10.7lf\tSVol=\t%10.7lf\tRVol=\t%10.7lf\trefPriceNew=\t%15.10lf\tSVar=\t%15.10lf\tRVar=\t%15.10lf\tCov=\t%15.10lf\t%s\t%s\n",
		(int)(Tcall*365),(int)(Tfix*365),
		spreadStrikes[i]*100,rateStrikes[i]*100,(spreadStrikes[i]+SKAdj)*100,(rateStrikes[i]+RKAdj)*100,
		newSpreadStdDev/sTfix*100,newRateStdDev/sTfix*100,refPriceNew,
		spreadVarModel+SVolNew*SVolNew*(Tfix-Tcall),rateVarModel+RVolNew*RVolNew*(Tfix-Tcall),
		spreadRateCovModel+SVolNew*RVolNew*SRCorNew*(Tfix-Tcall),
		(isCalibrateSpreadVol && isCalibrateRateVol) ? "SVol + RVol Calib" : (isCalibrateSpreadVol ? "SVol Calib" : (isCalibrateRateVol ? "RVol Calib" : "No Vol Calib")),
		(spreadLocalSqueeze && rateLocalSqueeze) ? "S+R Sqz" : (spreadLocalSqueeze ? "S Sqz" : (rateLocalSqueeze ? "R Sqz" : "No Sqz")));
fclose(f1);
****/

				}
				else
				{	
					/// Past reset w.r.t call date => set all adjustments to 0
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(KADJ_RATE_FLT,evalTimes[j],resetTimes[i],0.0);

					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_SPREAD,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_RATE,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_SPREAD_FLT,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(ADJ_RATE_FLT,evalTimes[j],resetTimes[i],0.0);

					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_SPREAD_FLT,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(VOL_RATE_FLT,evalTimes[j],resetTimes[i],0.0);

					GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(CORREL_SPREAD_RATE,evalTimes[j],resetTimes[i],0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(CORREL_SPREAD_RATE_FLT,evalTimes[j],resetTimes[i],0.0);
						
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DBLECOR_ADJ_PAY_FLT,evalTimes[j],payIndexResetTimes[periodIndex],0.0);
				} // if past reset flows w.r.t. call date

			} // for each call date

		} // if non null flow

		if( i == nbFlows-1 || ( i<nbFlows-1 && periodIndexes[i+1] != periodIndexes[i] ) )
		{
			if(spreadSqueeze)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : variance squeeze for spread condition in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				spreadSqueeze=false;
			}
			if(rateSqueeze)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : variance squeeze for rate condition in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				rateSqueeze=false;
			}
			if(correlSqueeze)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : covariance squeeze for spread/rate conditions in the period " << periodIndex;
				os << ", Max Squeeze = " << 100*maxCorrelSqueeze << "%" << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				correlSqueeze=false;
			}

			if(nullSpreadVolInput)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : strike too far from forward spread => null local vol set in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				nullSpreadVolInput=false;
			}
			if(nullRateVolInput)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : strike too far from forward rate => null local vol set in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				nullRateVolInput=false;
			}
			if(nullSpreadRateCorrelInput)
			{
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Double Corridor Option Local Calibration : null local vols or spread/rate variance squeeze => null local correl set in the period " << periodIndex << endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				nullSpreadRateCorrelInput=false;
			}
		}
	} // for each double condition reset date

	/// Save residual prices at the right place
	double price;
	size_t nbValues;
	bool isNewEvent;
	for(j=0;j<nbEvalTimes;++j)
	{
		nbValues = static_cast<const ARM_SurfaceListModelParam&>(GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation)).rows(DBLECOR_MKT_RESIDUAL_RA2_PRICE);
		isNewEvent = (nbValues<=j);
		for(i=0;i<nbFlows;++i)
		{
			if(resetTimes[i]>=evalTimes[j])
			{
				if(isNewEvent)
					price = 0.0;
				else
					price = GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).GetValue(DBLECOR_MKT_RESIDUAL_RA2_PRICE,evalTimes[j],resetTimes[i]);

				if( barrierType == DBLECOR_SDOWN_RDOWN || barrierType == DBLECOR_SUP_RUP ||
					((barrierType == DBLECOR_SDOWN_RUP || barrierType == DBLECOR_SUP_RDOWN) && spreadCallPut != rateCallPut) )
					price += residualPrices(j,i);
				else
					price -= residualPrices(j,i);
			}
			else
				price = 0.0;

			GetModelParams()->GetModelParam(ARM_ModelParamType::Correlation).SetValue(DBLECOR_MKT_RESIDUAL_RA2_PRICE,evalTimes[j],resetTimes[i],price);
		}
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : VolUnSqueezer
///	Returns : bool
///	Action  : Try to fit model vol covariance by shifting strike
////////////////////////////////////////////////////
bool ARM_Local_Normal_Model::VolUnSqueezer(
	double fwd, 
	double targetPrice,
	double strike,
	double mat, 
	int callPut,
	double targetVol,
	double &adj)
{
	double inf,sup;
	if (callPut==K_CAP)
	{
		inf=0.;
		sup=0.05;
	}
	else
	{
		inf=CC_Max<double>(-0.05,-strike+0.0001);
		sup=0.;
	}
	
	int iter;
	double tol=1e-6;
	double a = inf,b=sup,c=sup,d,e,min1,min2;

	double fa = VanillaOption_N(fwd,targetVol,strike+a,mat,callPut)-targetPrice;
	double fb = VanillaOption_N(fwd,targetVol,strike+b,mat,callPut)-targetPrice;

	double fc,p,q,r,s,tol1,xm;

	if((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) return false;
	fc=fb;
	for (iter=1;iter<=100;iter++) 
	{
		if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
		{
			c=a; 
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) 
		{
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*3e-8*fabs(b)+0.5*tol; 
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0)
		{
			adj=b;
			return true;
		}
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) 
		{
			s=fb/fa; 
			if (a == c) 
			{
				p=2.0*xm*s;
				q=1.0-s;
			} 
			else 
			{
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0) q = -q; 
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) 
			{
				e=d;
				d=p/q;
			} 
			else 
			{
				d=xm;
				e=d;
			}
		} 
		else 
		{ 
			d=xm;
			e=d;
		}
		a=b; 
		fa=fb;
		if (fabs(d) > tol1) b +=d;
		else
		b += (xm>=0 ? fabs(tol1):-fabs(tol1));
		fb=VanillaOption_N(fwd,targetVol,strike+b,mat,callPut)-targetPrice;
	}
	return false;
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CorrelUnSqueeze & associated err function
///	Returns : bool
///	Action  : Try to fit covariance by modifying input
///			  adjusted strikes
////////////////////////////////////////////////////
void NAG_CALL CorrelUnSqueezerFunc(Integer n,
	double* x, 
    double* fx,
	double* g,
    Nag_Comm *comm )
{
	ARM::ARM_Local_Normal_Model::CorrelUnSqueezerData* data = static_cast<ARM::ARM_Local_Normal_Model::CorrelUnSqueezerData*>(comm->p);
	*fx = data->ComputeErr(x);
}


bool ARM_Local_Normal_Model::CorrelUnSqueezer(double Tcall, double Tfix,
	double spreadTarget,double rateTarget,double spreadRateCorrelTarget,
	double spreadVarModel,double rateVarModel,double spreadRateCovModel,double sTfixTcall,
	double &newSpreadKAdj,double &newRateKAdj,double &newSpreadStdDev,double &newRateStdDev,
	double &newLocalSpreadStdDev,double &newLocalRateStdDev,double &newLocalSpreadRateCorrel,
	size_t& iter)
{
	double xS = (spreadTarget-newSpreadKAdj)/newSpreadStdDev;
	if(-K_NEW_DOUBLE_TOL<xS && xS<K_NEW_DOUBLE_TOL)
		return false; /// new fit is not possible !

	double xR = (rateTarget-newRateKAdj)/newRateStdDev;
	if(-K_NEW_DOUBLE_TOL<xR && xR<K_NEW_DOUBLE_TOL)
		return false; /// new fit is not possible !

	double localSpreadVarMin = LOCAL_SPREAD_VOL_MIN*sTfixTcall;
	localSpreadVarMin *= localSpreadVarMin;

	double spreadVarTarget0 = newSpreadStdDev*newSpreadStdDev;
	double localSpreadVar = spreadVarTarget0 - spreadVarModel;
	if(localSpreadVar < localSpreadVarMin)
		localSpreadVar = localSpreadVarMin;
	double spreadVarTarget		= spreadVarModel + localSpreadVar;
	double spreadVarTargetMin	= spreadVarModel + localSpreadVarMin;


	double localRateVarMin = LOCAL_RATE_VOL_MIN*sTfixTcall;
	localRateVarMin *= localRateVarMin;

	double rateVarTarget0 = newRateStdDev*newRateStdDev;
	double localRateVar = rateVarTarget0 - rateVarModel;
	if(localRateVar < localRateVarMin)
		localRateVar = localRateVarMin;
	double rateVarTarget	= rateVarModel + localRateVar;
	double rateVarTargetMin	= rateVarModel + localRateVarMin;

	newSpreadStdDev			= sqrt(spreadVarTarget);
	newRateStdDev			= sqrt(rateVarTarget);
	newLocalSpreadStdDev	= sqrt(localSpreadVar);
	newLocalRateStdDev		= sqrt(localRateVar);
	newLocalSpreadRateCorrel =  (newSpreadStdDev*newRateStdDev*spreadRateCorrelTarget - spreadRateCovModel)
								/ (newLocalSpreadStdDev*newLocalRateStdDev);

	if(fabs(newLocalSpreadRateCorrel)<=1.0)
	{
		/// No more squeeze : compute new strikes so that both binaries stay unchanged
		newSpreadKAdj			= spreadTarget - xS * newSpreadStdDev;
		newRateKAdj				= rateTarget - xR * newRateStdDev;
		return true;
	}


	/// Use NAG solver in 2D with constraints. The optimisation try to find
	/// a 
	double err;
    Nag_BoundType bound = Nag_Bounds;
    double x[2];
	double bl[2];
	double bu[2];
	double g[2];

	bl[0] = sqrt(spreadVarTargetMin);
	bl[1] = sqrt(rateVarTargetMin);
	
	/// Initialisation must be not too close to lower bounds !
	x[0] = newSpreadStdDev > bl[0]*1.1 ? newSpreadStdDev : bl[0]*1.1;
	x[1] = newRateStdDev > bl[1]*1.1 ? newRateStdDev : bl[1]*1.1;

	double localSpreadRateCorrelNagTarget=0.99;
	double localSpreadRateCorrelTarget=1.0;
	if(newLocalSpreadRateCorrel < -1)
	{
		localSpreadRateCorrelNagTarget=-0.99;
		localSpreadRateCorrelTarget=-1.0;
	}
	CorrelUnSqueezerData data(spreadVarModel,rateVarModel,spreadRateCovModel,spreadRateCorrelTarget,
		localSpreadRateCorrelNagTarget,x[0],x[1]);

	bu[0] = 1.0e+10; /// no upper bound
	bu[1] = 1.0e+10; /// no upper bound

    Nag_E04_Opt options;
	nag_opt_init(&options);
    Nag_Comm comm;
    comm.p = &data;

    static NagError fail,fail2;
	fail.print = FALSE;

/***
if(Tcall*365 == 3575 && Tfix*365 == 3589)
{
	strcpy(options.outfile, "c:\\temp\\dumpNAG.txt");
	options.list				= TRUE;
	options.print_level			= Nag_Soln_Iter;
}
***/

	options.list				= FALSE;
	options.print_level			= Nag_NoPrint;
	options.output_level		= Nag_NoOutput;
	options.minor_print_level	= Nag_NoPrint;
	options.print_deriv			= Nag_D_NoPrint;

	/// NAG optimisation with boundary using no derivatives (e04jbc)
	options.max_iter = 30;
	options.optim_tol = 0.0001*CC_Min(bl[0],bl[1]);
    nag_opt_bounds_no_deriv(2,&CorrelUnSqueezerFunc,bound,bl,bu,x, &err, g,&options, &comm, &fail);

	/// Output number of iteration
	iter = options.iter;

	/// free the option (e04xzc)
    nag_opt_free(&options,"all",&fail2);

	newSpreadStdDev				= x[0];
	newRateStdDev				= x[1];
	newLocalSpreadStdDev		= sqrt(x[0]*x[0] - spreadVarModel);
	newLocalRateStdDev			= sqrt(x[1]*x[1] - rateVarModel);
	newLocalSpreadRateCorrel	= (newSpreadStdDev*newRateStdDev*spreadRateCorrelTarget - spreadRateCovModel)
								  / (newLocalSpreadStdDev*newLocalRateStdDev);
	newSpreadKAdj				= spreadTarget - xS*newSpreadStdDev;
	newRateKAdj					= rateTarget - xR*newRateStdDev;

	/// Sucessful correction if local correlation under prices the theoretical target
	/// or over prices but less than 2%
	return localSpreadRateCorrelTarget*(newLocalSpreadRateCorrel-localSpreadRateCorrelTarget) <= 0.02;
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateCapFloor
///	Returns : void
///	Action  : Local model calibration to a 
///			  capl / floor
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateCapFloor(	ARM_CapFloor* capFloor,
													double targetPrice, 
													const ARM_GP_Vector& evalTimes)
{
	/// target price is ignored because we need the pv of all caplets / floorlets
	/// --> use GetCashFlowValues()
	size_t capSize = capFloor->GetResetDates()->GetSize();
	ARM_Vector* cashFlows = capFloor->GetCashFlowValues();

	for (size_t i(0); i<capSize; i++)
		CalibrateCapFloorLet(capFloor, i, cashFlows->Elt(i), evalTimes);
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Normal_Model
///	Routine : CalibrateCapFloorLet
///	Returns : void
///	Action  : Local model calibration to a 
///			  caplet / floorlet
////////////////////////////////////////////////////
void ARM_Local_Normal_Model::CalibrateCapFloorLet(	ARM_CapFloor* capFloor,
													size_t periodIdx,
													double targetPrice, 
													const ARM_GP_Vector& evalTimes)
{	
	// will generate an exception if numerical model is not supported for calibration
	Local_Normal_Model_Calibration_Helper helper (GetNumericalModel());

	size_t i, sizeEvalTimes = evalTimes.size();
	
	// get computed libor rate (potentially adjusted)
	double liborRateTarget = capFloor->GetSwapLeg()->GetFwdRates()->Elt(periodIdx);
	liborRateTarget *= 0.01;
	
	double resetTime, fwdStartTime, fwdEndTime, fwdPeriod, payTime, payPeriod, payNotional, strike;
	int callPut;
	
	GetCapFloorLetData(capFloor, periodIdx, resetTime, payTime, payPeriod, payNotional, fwdStartTime, fwdEndTime, fwdPeriod, strike, callPut);
	
	double dfPay = GetNumericalModel()->GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);
	double Tfix	 = resetTime / K_YEAR_LEN;
	
	targetPrice /= (payNotional * dfPay * payPeriod) ;
	double targetVol;

	
	bool success (true);
	if (!LiborAsShiftedLognormal())
		targetVol = VanillaImpliedVol_N (liborRateTarget, targetPrice, strike, Tfix, callPut, NULL, &success);
	else
		targetVol = VanillaImpliedVol_BS (	liborRateTarget + 1./fwdPeriod, 
											strike + 1./fwdPeriod, 
											Tfix, 
											targetPrice, 
											callPut,
											NULL,
											&success);
	
	///
	///  compute forward spread vol & forward adjustments
	///
	for (i=0; i<sizeEvalTimes; i++)
	{
		// calibration only if option maturity is after evalTime
		if (resetTime >= evalTimes[i] )
		{			
			/// adjusment on libor fwd
			double liborRateModel = helper.LiborRateFromNumericalModel(evalTimes[i], 
																		payTime, 
																		fwdStartTime, 
																		fwdEndTime, 
																		fwdPeriod, 
																		strike);
				
			// this adjust will be = 0 if standard libor
			double adjustment  = liborRateTarget  - liborRateModel;
						

			// set adjustments to model params
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(CAPLET_ADJ, evalTimes[i], resetTime, adjustment);
						
			/// compute libor forward volatility
			/// if the implied vol computation has failed, this means that the strike is too far from ATM
			/// in that case, a forward vol = 0 is inputed
			///
			double forwardVol;
			
			if (success)
			{
				double modelVol; /// model spread vol up to evalTimes[i] 
			
				if (!LiborAsShiftedLognormal())
					modelVol = helper.LiborNormalVolatilityFromNumericalModel (	evalTimes[i], 
																				fwdStartTime, 
																				fwdEndTime, 
																				fwdPeriod, 
																				strike);
				else
					modelVol = helper.FwdZcVolatilityFromNumericalModel  (  evalTimes[i], 
																			fwdStartTime, 
																			fwdEndTime, 
																			fwdPeriod, 
																			strike/* not used for HW, used for QGM1F*/);

				
				
				double Tcall = evalTimes[i] / K_YEAR_LEN;
				
				if (Tcall==Tfix)
					Tfix += 1./K_YEAR_LEN;
				
				double var = (	targetVol  * targetVol * Tfix 
							  -	modelVol   * modelVol  * Tcall) ;

				/// everything is fine
				if (var>=0)
					forwardVol = sqrt( var / (Tfix - Tcall) );
				
				/// oups: variance squeeze
				else 
				{
					forwardVol = 0.0;
					SetVarianceSqueezeStatus(true);
					CC_Ostringstream os;
					os << ARM_USERNAME << " : Local_Normal_Model / Caplet/Floorlet Calibration : variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(Tfix*10.0)/10.0 << " Y]"<< endl;
					ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				}
			}
			else
			{
				forwardVol = 0.0;
				CC_Ostringstream os;
				os << ARM_USERNAME << " : Local_Normal_Model / Caplet/Floorlet Calibration : strike seems to be to far from ATM. Null forward vol inputed..."<< endl;
				ARM_TheEventViewer.Instance()->AddToMessage(os.str());
			}

			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(CAPLET_VOL, evalTimes[i], resetTime, forwardVol);
		}
		else
		{	
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(CAPLET_ADJ, evalTimes[i], resetTime, 0.0);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(CAPLET_VOL, evalTimes[i], resetTime, 0.0);
		}
	}
}
			
///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ
////////////////////////////////////////////////////////
/// Local_Normal_Model_Helper implementation ///////////
////////////////////////////////////////////////////////
///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : constructor
///	Returns : void
///	Action  : contextual constructor
////////////////////////////////////////////////////
Local_Normal_Model_Calibration_Helper::Local_Normal_Model_Calibration_Helper (ARM_PricingModel* numericalModel)
:	itsNumericalModel (numericalModel),
	itsModelType (ARM_PricingModelType::HWM1F)
{
	if ( dynamic_cast<ARM_HullWhite1F*> (numericalModel) )
	{
		if ( dynamic_cast<ARM_ModelParamsHW1FStd*> (numericalModel->GetModelParams()) )
			itsModelType = ARM_PricingModelType::HWM1F;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model Calibration : numerical model is not supported" );
	}
	else if ( dynamic_cast<ARM_HullWhite2F*> (numericalModel) )
	{
		itsModelType = ARM_PricingModelType::HWM2F;
	}
	else if ( dynamic_cast<ARM_QGM1F*> (numericalModel) )
	{
		itsModelType = ARM_PricingModelType::QGM1F;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model Calibration : numerical model is not supported" );
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : constructor
///	Returns : void
///	Action  : copy constructor
////////////////////////////////////////////////////
Local_Normal_Model_Calibration_Helper::Local_Normal_Model_Calibration_Helper (const Local_Normal_Model_Calibration_Helper& rhs)
:	itsNumericalModel	(rhs.itsNumericalModel),
	itsModelType		(rhs.itsModelType)
{
}



////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : LiborRateFromNumericalModel
///	Returns : void
///	Action  : LIBOR rate under num. model, 
///			 fixed @ expiryTime & payed @ payTime 
///			 ( = E^Q(Tp)[L(Texp)] )
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::LiborRateFromNumericalModel(
																		double expiryTime,
																		double payTime, 
																		double fwdStartTime,
																		double fwdEndTime,
																		double fwdPeriod,
																		double strike ) const
{		
	if (itsModelType == ARM_PricingModelType::HWM1F || itsModelType == ARM_PricingModelType::HWM2F)
	{
		double dfStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdStartTime / K_YEAR_LEN);
		double dfEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdEndTime / K_YEAR_LEN);
		double cov;

		ARM_ModelParamsHW* params = (ARM_ModelParamsHW*)itsNumericalModel->GetModelParams();
		ARM_GP_Vector notUsed (2);
		cov   = params->FwdZcLocalCovariance(0, expiryTime, fwdStartTime, fwdEndTime, payTime, fwdEndTime, notUsed);
		
		return ( (dfStart/dfEnd) * exp(cov) + - 1.0 ) / fwdPeriod;
	}
	else if (itsModelType == ARM_PricingModelType::QGM1F)
	{		
		ARM_VectorPtr libor = ((ARM_QGM1F*)itsNumericalModel)->Libor("", 0/*as of*/, fwdStartTime, fwdEndTime, fwdPeriod, expiryTime, payTime, ARM_PricingStatesPtr( new ARM_PricingStates(1,1,0) ) );
		return (*libor)[0];
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::LiborRateFromNumericalModel : numerical model is not supported" );
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : LiborNormalVolatilityFromNumericalModel
///	Returns : void
///	Action  : normal vol of a libor on [0, expiryTime]
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::LiborNormalVolatilityFromNumericalModel(
																		double expiryTime,
																//		double payTime, // not used in HW case
																		double fwdStartTime,
																		double fwdEndTime,
																		double fwdPeriod,
																		double strike ) const
{		
	//
	// caution : gaussian approx on libor (which is lognorm. shifted under HW)
	//
	double dfStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdStartTime / K_YEAR_LEN);
	double dfEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdEndTime / K_YEAR_LEN);
	
	if (itsModelType == ARM_PricingModelType::HWM1F || itsModelType == ARM_PricingModelType::HWM2F)
	{	
		ARM_ModelParamsHW* params = (ARM_ModelParamsHW*)itsNumericalModel->GetModelParams();

		double vol = params->FwdZcLocalVariance(0, expiryTime, fwdStartTime, fwdEndTime);
		vol = sqrt( vol / (expiryTime/K_YEAR_LEN) );
		return (1./fwdPeriod) * vol * dfStart / dfEnd;
	}
	///
	/// To be modified ...
	/// this works only for QGM1F with Skew = 0 (i.e. HW1F)
	else if (itsModelType == ARM_PricingModelType::QGM1F)
	{
		// caplet price under QG
		double fwd		= (dfStart/dfEnd - 1.0)/fwdPeriod;
		// choose OTM option...
		int callPut		= (fwd>strike) ? -1 : 1 ;
		ARM_VectorPtr caplet = ((ARM_QGM1F*)itsNumericalModel)->VanillaCaplet("", 0 /*as of date*/, fwdEndTime, 1.0, 1.0, expiryTime, fwdStartTime, fwdEndTime, fwdPeriod, ARM_GP_Vector(1, strike), callPut, ARM_PricingStatesPtr( new ARM_PricingStates(1,1,0) ) );
		double price	= (*caplet)[0] / dfEnd;
		// corresponding implied vol for Libor under QGM1F
		double qgmVol	= VanillaImpliedVol_N  (fwd, 
												price,
												strike, 
												expiryTime / K_YEAR_LEN, 
												callPut);
		return qgmVol;

	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::LiborNormalVolatilityFromNumericalModel : numerical model is not supported" );
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : FwdZcVolatilityFromNumericalModel
///	Returns : void
///	Action  : lognormal vol of B(.,Tstart)/B(.,Tend) on [0, expiryTime]
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::FwdZcVolatilityFromNumericalModel(
																		double expiryTime,
																//		double payTime, // not used in HW case
																		double fwdStartTime,
																		double fwdEndTime,
																		double fwdPeriod,
																		double strike ) const
{
	if (itsModelType == ARM_PricingModelType::HWM1F || itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_ModelParamsHW* params = (ARM_ModelParamsHW*)itsNumericalModel->GetModelParams();

		double var = params->FwdZcLocalVariance(0, expiryTime, fwdStartTime, fwdEndTime);
		return  sqrt( var / (expiryTime/K_YEAR_LEN) );
	}
	///
	/// To be modified ...
	/// this works only for QGM1F with Skew = 0 (i.e. HW1F)
	else if (itsModelType == ARM_PricingModelType::QGM1F)
	{
		// caplet price under QG
		double dfStart	= itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdStartTime / K_YEAR_LEN);
		double dfEnd	= itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdEndTime / K_YEAR_LEN);
		double fwd		= (dfStart/dfEnd - 1.0)/fwdPeriod;
		// choose OTM option...
		int callPut		= (fwd>strike) ? -1 : 1 ;
		ARM_VectorPtr caplet = ((ARM_QGM1F*)itsNumericalModel)->VanillaCaplet("", 0 /*as of date*/, fwdEndTime, 1.0, 1.0, expiryTime, fwdStartTime, fwdEndTime, fwdPeriod, ARM_GP_Vector(1, strike), callPut, ARM_PricingStatesPtr( new ARM_PricingStates(1,1,0) ) );
		double price	= (*caplet)[0] / dfEnd;
		// corresponding implied vol for fwd discount factor
		double qgmVol	= VanillaImpliedVol_BS (fwd + 1./fwdPeriod, 
												strike + 1./fwdPeriod, 
												expiryTime / K_YEAR_LEN, 
												price, 
												callPut);
		return qgmVol;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::Local_Normal_Model_Calibration_Helper::FwdZcVolatilityFromNumericalModel : numerical model is not supported" );
}
	
////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : CmsRateFromNumericalModel
///	Returns : convexified CMS rate
///	Action  : CMS rate under num. model, 
///			 fixed @ expiryTime & payed @ payTime 
///			 ( = E^Q(Tp)[S(Texp)] )
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::CmsRateFromNumericalModel(
																	double expiryTime,
																	double payTime, 
																	double swapFloatStartTime,	// adjusted ...
																	double swapFloatEndTime,	// adjusted ...
																	const ARM_GP_Vector& swapFixPayTimes,
																	const ARM_GP_Vector& swapFixPayPeriods,
																	double* factor1, double* factor2) const
{
	double cms=0.0;

	size_t i, size = swapFixPayTimes.size();
		
	double dfStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFloatStartTime / K_YEAR_LEN);
	double dfEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFloatEndTime / K_YEAR_LEN);

	ARM_GP_Vector swapFixPayDfs (size);

	for (i=0; i<size; i++)
		swapFixPayDfs[i] = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFixPayTimes[i]/K_YEAR_LEN );
	
	double swapRate;	
	///
	/// HW 1F case (constant mean reversion)
	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_HullWhite1F* mod = (ARM_HullWhite1F*)itsNumericalModel;
		ARM_ModelParamsHW1FStd* params = (ARM_ModelParamsHW1FStd*)mod->GetModelParams();

		double swapVolFactor (0);
		double numeraireVolFactor (0);
		double lambda =	((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

		swapVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs, 
													&swapRate,
													payTime, 
													&numeraireVolFactor);
		if(factor1)
			/// Only output H&W vol factor
			*factor1 = swapVolFactor;
		else
		{
			double phi = params->StateLocalVariance(0.0, expiryTime, expiryTime);

			cms = swapRate + numeraireVolFactor * swapVolFactor * phi;
		}
	}
	///
	/// HW 2F case (constant mean reversion)
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_HullWhite2F* mod = (ARM_HullWhite2F*)itsNumericalModel;
		ARM_ModelParamsHW2FStd* params = (ARM_ModelParamsHW2FStd*)mod->GetModelParams();

		double swapVolFactor1, swapVolFactor2;
		double numeraireVolFactor1, numeraireVolFactor2;
		double lambda1 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
		double lambda2 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
		lambda2 += lambda1;

		swapVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs, 
													&swapRate,
													payTime, 
													&numeraireVolFactor1);
		
		swapVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs, 
													&swapRate,
													payTime, 
													&numeraireVolFactor2);
		if(factor1 && factor2)
		{
			/// Only output H&W vol factors
			*factor1 = swapVolFactor1;
			*factor2 = swapVolFactor2;
		}
		else
		{
			ARM_GP_TriangularMatrix* phi = params->StateLocalVariance(0.0, expiryTime);
			double phi11 = (*phi)(0, 0);
			double phi22 = (*phi)(1, 1);
			double phi12 = (*phi)(1, 0);
			delete phi;
			
			cms = swapRate;
			cms += numeraireVolFactor1 * swapVolFactor1 * phi11; 
			cms += numeraireVolFactor2 * swapVolFactor2 * phi22; 
			cms += numeraireVolFactor1 * swapVolFactor2 * phi12; 
			cms += numeraireVolFactor2 * swapVolFactor1 * phi12; 
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented yet" );
		
	return cms;
}


////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : SpreadNormalVolatilityFromNumericalModel
///	Returns : normal CMS spread volatility
///	Action  : normal vol of a CMS spread on [0, exprityTime]
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::SpreadNormalVolatilityFromNumericalModel(
																		double expiryTime,
																		double payTime, // not used in HW case
																		double coeffLong,
																		double coeffShort,
																		double strike, // not used in HW case
																		double swapLongFloatStartTime,	// adjusted ...
																		double swapLongFloatEndTime,	// adjusted ...
																		const ARM_GP_Vector& swapLongFixPayTimes,
																		const ARM_GP_Vector& swapLongFixPayPeriods,
																		double swapShortFloatStartTime,
																		double swapShortFloatEndTime,
																		const ARM_GP_Vector& swapShortFixPayTimes,
																		const ARM_GP_Vector& swapShortFixPayPeriods,
																		double* factor1, double* factor2) const
{
	double spreadVol=0.0;

	size_t i;

	// Computations for Long CMS
	size_t sizeLong = swapLongFixPayTimes.size();
		
	double dfStartLong = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapLongFloatStartTime / K_YEAR_LEN);
	double dfEndLong   = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapLongFloatEndTime / K_YEAR_LEN);

	ARM_GP_Vector swapLongFixPayDfs (sizeLong);

	for (i=0; i<sizeLong; i++)
		swapLongFixPayDfs[i] = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapLongFixPayTimes[i]/K_YEAR_LEN );
	
	// Computations for Short CMS
	size_t sizeShort = swapShortFixPayTimes.size();
		
	double dfStartShort = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapShortFloatStartTime / K_YEAR_LEN);
	double dfEndShort   = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapShortFloatEndTime / K_YEAR_LEN);

	ARM_GP_Vector swapShortFixPayDfs (sizeShort);

	for (i=0; i<sizeShort; i++)
		swapShortFixPayDfs[i] = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapShortFixPayTimes[i]/K_YEAR_LEN );
	
	///
	/// HW 1F case (constant mean reversion)
	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_HullWhite1F* mod = (ARM_HullWhite1F*)itsNumericalModel;
		ARM_ModelParamsHW1FStd* params = (ARM_ModelParamsHW1FStd*)mod->GetModelParams();
		double lambda =	((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

		double swapLongVolFactor  (0);
		double swapShortVolFactor (0);
		
		swapLongVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
													expiryTime, 
													swapLongFloatStartTime, 
													swapLongFloatEndTime, 
													dfStartLong, 
													dfEndLong, 
													swapLongFixPayTimes, 
													swapLongFixPayPeriods, 
													swapLongFixPayDfs);
		
		swapShortVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
													expiryTime, 
													swapShortFloatStartTime, 
													swapShortFloatEndTime, 
													dfStartShort, 
													dfEndShort, 
													swapShortFixPayTimes, 
													swapShortFixPayPeriods, 
													swapShortFixPayDfs);
													
		
		
		double spreadFactor = coeffLong * swapLongVolFactor - coeffShort * swapShortVolFactor;
		if(factor1)
			/// Only output H&W vol factor
			*factor1 = spreadFactor;
		else
		{

			double phi = params->StateLocalVariance(0.0, expiryTime, expiryTime);
			
			spreadVol = spreadFactor * spreadFactor * phi;
			if ( fabs(expiryTime)<=K_NEW_DOUBLE_TOL )
				spreadVol = 0.0;
			else
				spreadVol = sqrt (spreadVol / (expiryTime / K_YEAR_LEN));
		}
	}
	///
	/// HW 2F case (constant mean reversion)
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_HullWhite2F* mod = (ARM_HullWhite2F*)itsNumericalModel;
		ARM_ModelParamsHW2FStd* params = (ARM_ModelParamsHW2FStd*)mod->GetModelParams();
		
		double lambda1 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
		double lambda2 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
		lambda2 += lambda1;

		double swapLongVolFactor1, swapLongVolFactor2;
		double swapShortVolFactor1, swapShortVolFactor2;

		swapLongVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
													expiryTime, 
													swapLongFloatStartTime, 
													swapLongFloatEndTime, 
													dfStartLong, 
													dfEndLong, 
													swapLongFixPayTimes, 
													swapLongFixPayPeriods, 
													swapLongFixPayDfs);
		swapLongVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
													expiryTime, 
													swapLongFloatStartTime, 
													swapLongFloatEndTime, 
													dfStartLong, 
													dfEndLong, 
													swapLongFixPayTimes, 
													swapLongFixPayPeriods, 
													swapLongFixPayDfs);
		swapShortVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
													expiryTime, 
													swapShortFloatStartTime, 
													swapShortFloatEndTime, 
													dfStartShort, 
													dfEndShort, 
													swapShortFixPayTimes, 
													swapShortFixPayPeriods, 
													swapShortFixPayDfs);
		swapShortVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
													expiryTime, 
													swapShortFloatStartTime, 
													swapShortFloatEndTime, 
													dfStartShort, 
													dfEndShort, 
													swapShortFixPayTimes, 
													swapShortFixPayPeriods, 
													swapShortFixPayDfs);
		

		double spreadFactor1 = coeffLong * swapLongVolFactor1 - coeffShort * swapShortVolFactor1;
		double spreadFactor2 = coeffLong * swapLongVolFactor2 - coeffShort * swapShortVolFactor2;
		if(factor1 && factor2)
		{
			/// Only output H&W vol factors
			*factor1 = spreadFactor1;
			*factor2 = spreadFactor2;
		}
		else
		{
			ARM_GP_TriangularMatrix* phi = params->StateLocalVariance(0.0, expiryTime);
			double phi11 = (*phi)(0, 0);
			double phi22 = (*phi)(1, 1);
			double phi12 = (*phi)(1, 0);
			delete phi;
			
			
			spreadVol =				spreadFactor1 * spreadFactor1 * phi11
						+			spreadFactor2 * spreadFactor2 * phi22
						+	2.0  *	spreadFactor1 * spreadFactor2 * phi12 ;

			if ( fabs(expiryTime)<=K_NEW_DOUBLE_TOL )
				spreadVol = 0.0;
			else
				spreadVol = sqrt (spreadVol / (expiryTime / K_YEAR_LEN));
		}
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented" );
		
	return spreadVol;
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : FLTconvexifiedRateFromNumericalModel
///	Returns : convexified CMS rate w.r.t. the payment index numeraire
///	Action  : input : CMS rate under num. model, fixed @ expiryTime & payed @ payTime
///			  output : E^Q(Num)[CMS(Texp)] i.e. expectation under the numeraire Num(t)= B(t,Tp).Spay(t)
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::FLTconvexifiedRateFromNumericalModel(
														double cms, 
														double expiryTime,
														double swapFloatStartTime,	// adjusted ...
														double swapFloatEndTime,	// adjusted ...
														const ARM_GP_Vector& swapFixPayTimes,
														const ARM_GP_Vector& swapFixPayPeriods,
														/// -- payed rate
														double payIndexFloatStartTime,	// adjusted ...
														double payIndexFloatEndTime,	// adjusted ...
														const ARM_GP_Vector& payIndexFixPayTimes,
														const ARM_GP_Vector& payIndexFixPayPeriods) const
{	
	size_t i, size = swapFixPayTimes.size();
	size_t paySize = payIndexFixPayTimes.size();
		
	double dfStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFloatStartTime / K_YEAR_LEN);
	double dfEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFloatEndTime / K_YEAR_LEN);

	double dfPayStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(payIndexFloatStartTime / K_YEAR_LEN);
	double dfPayEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(payIndexFloatEndTime / K_YEAR_LEN);

	ARM_GP_Vector swapFixPayDfs (size);
	ARM_GP_Vector payIndexFixPayDfs (paySize);

	for (i=0; i<size; i++)
		swapFixPayDfs[i] = itsNumericalModel->GetZeroCurve()->DiscountPrice(swapFixPayTimes[i]/K_YEAR_LEN );
	
	for (i=0; i<paySize; i++)
		payIndexFixPayDfs[i] = itsNumericalModel->GetZeroCurve()->DiscountPrice(payIndexFixPayTimes[i]/K_YEAR_LEN );
	
	
	double payIndexFwd;	
	///
	/// HW 1F case (constant mean reversion)
	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_HullWhite1F* mod = (ARM_HullWhite1F*)itsNumericalModel;
		ARM_ModelParamsHW1FStd* params = (ARM_ModelParamsHW1FStd*)mod->GetModelParams();

		double swapVolFactor (0), payIndexVolFactor(0);
		double lambda =	((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];

		swapVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs);

		payIndexVolFactor = ARM_HullWhite::SwapRateVolFactor (lambda, 
													expiryTime, 
													payIndexFloatStartTime, 
													payIndexFloatEndTime, 
													dfPayStart, 
													dfPayEnd, 
													payIndexFixPayTimes, 
													payIndexFixPayPeriods, 
													payIndexFixPayDfs,
													&payIndexFwd);
													
							
		double phi = params->StateLocalVariance(0.0, expiryTime, expiryTime);
		
		/// (payIndexVolFactor / payIndexFwd) is an approx of lognormal vol of payIndex...
		cms += (payIndexVolFactor / payIndexFwd) * swapVolFactor * phi; // heu..., au signe près

		return cms;
	}
	///
	/// HW 2F case (constant mean reversion)
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{

		ARM_HullWhite2F* mod = (ARM_HullWhite2F*)itsNumericalModel;
		ARM_ModelParamsHW2FStd* params = (ARM_ModelParamsHW2FStd*)mod->GetModelParams();

		double swapVolFactor1, swapVolFactor2;
		double payIndexVolFactor1, payIndexVolFactor2;
		double lambda1 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversion)).GetCurve()->GetOrdinates()[0];
		double lambda2 =((ARM_CurveModelParam&)params->GetModelParam(ARM_ModelParamType::MeanReversionSpread)).GetCurve()->GetOrdinates()[0];
		lambda2 += lambda1;

		swapVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs);
		
		swapVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
													expiryTime, 
													swapFloatStartTime, 
													swapFloatEndTime, 
													dfStart, 
													dfEnd, 
													swapFixPayTimes, 
													swapFixPayPeriods, 
													swapFixPayDfs);

		payIndexVolFactor1 = ARM_HullWhite::SwapRateVolFactor (lambda1, 
													expiryTime, 
													payIndexFloatStartTime, 
													payIndexFloatEndTime, 
													dfPayStart, 
													dfPayEnd, 
													payIndexFixPayTimes, 
													payIndexFixPayPeriods, 
													payIndexFixPayDfs);

		payIndexVolFactor2 = ARM_HullWhite::SwapRateVolFactor (lambda2, 
													expiryTime, 
													payIndexFloatStartTime, 
													payIndexFloatEndTime, 
													dfPayStart, 
													dfPayEnd, 
													payIndexFixPayTimes, 
													payIndexFixPayPeriods, 
													payIndexFixPayDfs,
													&payIndexFwd);
													

		

		ARM_GP_TriangularMatrix* phi = params->StateLocalVariance(0.0, expiryTime);
		double phi11 = (*phi)(0, 0);
		double phi22 = (*phi)(1, 1);
		double phi12 = (*phi)(1, 0);
		delete phi;
		
		/// (payIndexVolFactor / payIndexFwd) is an approx of lognormal vol of payIndex...
		payIndexVolFactor1 /= payIndexFwd;
		payIndexVolFactor2 /= payIndexFwd;
		cms += payIndexVolFactor1 * swapVolFactor1 * phi11; 
		cms += payIndexVolFactor2 * swapVolFactor2 * phi22; 
		cms += payIndexVolFactor1 * swapVolFactor2 * phi12; 
		cms += payIndexVolFactor2 * swapVolFactor1 * phi12; 
		// heu..., au signe près

		return cms;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented yet" );
		
	return 0.0;
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : SpreadRateNormalCovarianceFromNumericalModel
///	Returns : void
///	Action  : input : CMS Spread = long*CMSLong - short*CMSShort & rate
///			  output : integrated normal covariance from 0 to expiryTime
////////////////////////////////////////////////////
double Local_Normal_Model_Calibration_Helper::SpreadRateNormalCovarianceFromNumericalModel(
														double expiryTime,
														double payTime,
														double coeffLong,
														double coeffShort,
														double strike,
														double swapLongFloatStartTime,
														double swapLongFloatEndTime,
														const ARM_GP_Vector& swapLongFixPayTimes,
														const ARM_GP_Vector& swapLongFixPayPeriods,
														double swapShortFloatStartTime,
														double swapShortFloatEndTime,
														const ARM_GP_Vector& swapShortFixPayTimes,
														const ARM_GP_Vector& swapShortFixPayPeriods,
														double rateStrike,
														double rateFloatStartTime,
														double rateFloatEndTime,
														const ARM_GP_Vector& rateFixPayTimes,
														const ARM_GP_Vector& rateFixPayPeriods,
														double rateCoeffLong,
														double rateCoeffShort,
														double rateShortFloatStartTime,
														double rateShortFloatEndTime,
														const ARM_GP_Vector& rateShortFixPayTimes,
														const ARM_GP_Vector& rateShortFixPayPeriods) const
{
	double spreadRateCovariance = 0.0;

	/// Compute CMS spread volatility factors
	double spreadFactor1,spreadFactor2;
	double unsedVol = SpreadNormalVolatilityFromNumericalModel(	expiryTime,
																payTime,
																coeffLong,
																coeffShort,
																strike,
																swapLongFloatStartTime,
																swapLongFloatEndTime,
																swapLongFixPayTimes,
																swapLongFixPayPeriods,
																swapShortFloatStartTime,
																swapShortFloatEndTime,
																swapShortFixPayTimes,
																swapShortFixPayPeriods,
																&spreadFactor1,&spreadFactor2);

	double rateFactor1,rateFactor2;
	if(rateShortFixPayTimes.size()==0)
	{
		/// Compute rate volatility factors
		double unusedRate = CmsRateFromNumericalModel(	expiryTime,
														payTime, 
														rateFloatStartTime,
														rateFloatEndTime,
														rateFixPayTimes,
														rateFixPayPeriods,
														&rateFactor1,&rateFactor2);
	}
	else
	{
		// Compute CMS spread volatility factors
		double unsedVol = SpreadNormalVolatilityFromNumericalModel(	expiryTime,
																	payTime,
																	rateCoeffLong,
																	rateCoeffShort,
																	rateStrike,
																	rateFloatStartTime,
																	rateFloatEndTime,
																	rateFixPayTimes,
																	rateFixPayPeriods,
																	rateShortFloatStartTime,
																	rateShortFloatEndTime,
																	rateShortFixPayTimes,
																	rateShortFixPayPeriods,
																	&rateFactor1,&rateFactor2);
	}


	/// HW 1F case (constant mean reversion)
	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_HullWhite1F* mod = (ARM_HullWhite1F*)itsNumericalModel;
		ARM_ModelParamsHW1FStd* params = (ARM_ModelParamsHW1FStd*)mod->GetModelParams();

		double phi = params->StateLocalVariance(0.0, expiryTime, expiryTime);

		spreadRateCovariance = spreadFactor1 * rateFactor1 * phi;
	}


	/// HW 2F case (constant mean reversion)
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{
		ARM_HullWhite2F* mod = (ARM_HullWhite2F*)itsNumericalModel;
		ARM_ModelParamsHW2FStd* params = (ARM_ModelParamsHW2FStd*)mod->GetModelParams();

		ARM_GP_TriangularMatrix* phi = params->StateLocalVariance(0.0, expiryTime);
		double phi11 = (*phi)(0, 0);
		double phi22 = (*phi)(1, 1);
		double phi12 = (*phi)(1, 0);
		delete phi;

		spreadRateCovariance =		spreadFactor1 * rateFactor1 * phi11
								+	spreadFactor2 * rateFactor2 * phi22
								+	(spreadFactor1 * rateFactor2 + spreadFactor2 * rateFactor1) * phi12;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "not implemented yet" );

	return spreadRateCovariance;
}

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

