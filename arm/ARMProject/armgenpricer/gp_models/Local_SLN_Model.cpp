/*!
 *
 * Copyright (c) IXIS CIB July 2005 Paris
 *
 *	\file Local_SLN_Model.cpp
 *
 *  \brief class for Shifted LogNormal local model
 *	\author  J-M Prié
 *	\version 1.0
 *	\date July 2005
 */


/// remove identified warning
#include "gpbase/removeidentifiedwarning.h"

/// gpmodels
#include "gpmodels/local_sln_model.h"
#include "gpmodels/Local_sln_modelparams.h"
#include "gpmodels/SFRM.h"
#include "gpmodels/ModelParamsSFRM.h"
#include "gpmodels/QModelAnalytics.h"
#include "gpmodels/q1f_fx.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/surfacelistmodelparam.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/pricingcontext.h"

/// gpbase
#include "gpbase/surface.h"
#include "gpbase/curve.h"
#include "gpbase/eventviewerfwd.h"
#include "gpbase/stringconvert.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"
#include "gpclosedforms/vanilla_bs.h"

/// gpnumlib
#include "gpnumlib/newtonraphson.h"

// kernel
#include "capfloor.h"
#include "option.h"
#include "corridorleg.h"
#include <util/fromto.h>

//gpcalculators
#include "gpcalculators/fxracalculator.h"


CC_BEGIN_NAMESPACE( ARM )

const size_t SURFACE_SIZE	=	10; //size of the Surface list 

const size_t FIRST_SURFACE	        = 0; // the first surface in the list
const size_t TARGET_PRICE_SURFACE   = 1; // to allow numerical calibration at pricing time
const size_t NUM_CALIB_SURFACE      = 2; // to save numerical calibration result at pricing time

// surfaces for the corridorlet at barrier DOWN or UP priced as call spread at strikes (K-epsilon) and (K+epsilon)
const size_t DOWN_MINUS_SURFACE	=	0; 
const size_t DOWN_PLUS_SURFACE	=	1;
const size_t UP_MINUS_SURFACE	=	2;
const size_t UP_PLUS_SURFACE	=	3;

const size_t DOWN_MINUS_SURFACE_PAYFLOAT	=	4; 
const size_t DOWN_PLUS_SURFACE_PAYFLOAT		=	5;
const size_t UP_MINUS_SURFACE_PAYFLOAT		=	6;
const size_t UP_PLUS_SURFACE_PAYFLOAT		=	7;

const double CALL_SPREAD_SHIFT  =	0.0001; //Call_Spread shift for corridor pricing =1bp

const size_t NB_STDDEV_NO_CALIB         = 8; 
const size_t NB_STDDEV_NO_CALIB_EQFX    = 6;

//the adjustment surfaces for corridor paying float Pricing 
//the first FIRST_SURFACE index is given to the ref_index  PayAdj   
const size_t DOWN_ADJ_SURFACE	=	1; //the covariance between the ref_idx et the pay_idx at with the variance smiled at the level Down
const size_t UP_ADJ_SURFACE		=	2; //the covariance between the ref_idx et the pay_idx at with the variance smiled at the level Up
const size_t PAY_ADJ_SURFACE	=	3;  // the Pay_index  PayAdj

//the shift surfaces for corridor paying float Pricing 
//the first FIRST_SURFACE index is given to the ref_index  shift   
const size_t PAY_IDX_SHIFT_SURFACE	=	1;

const size_t COL_NB_STATIC_RESULTS	=	10;


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Local_SLN_Model
////////////////////////////////////////////////////
ARM_Local_SLN_Model::ARM_Local_SLN_Model(const ARM_Local_SLN_ModelParams& params)
:	ARM_Local_Model( ARM_ZeroCurvePtr(NULL), params ),
itsResizeStatisticResult(true),
itsStatisticCalibResult( new ARM_GP_Matrix( 1,COL_NB_STATIC_RESULTS) )
{}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_Local_SLN_Model::ARM_Local_SLN_Model(const ARM_Local_SLN_Model& rhs)
:	ARM_Local_Model (rhs),
itsResizeStatisticResult(rhs.itsResizeStatisticResult),
itsStatisticCalibResult( ARM_GP_MatrixPtr(static_cast<ARM_GP_Matrix*>(rhs.itsStatisticCalibResult->Clone())))
{}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_Local_SLN_Model& ARM_Local_SLN_Model::operator = (const ARM_Local_SLN_Model& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_Local_SLN_Model();
		new (this) ARM_Local_SLN_Model (rhs);
	}
	return *this;
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: Destructor
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_Local_SLN_Model::~ARM_Local_SLN_Model()
{}
////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: Libor
///	Returns: ARM_VectorPtr
///	Action : computes the Libor rate in local normal model
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_SLN_Model::Libor( 
		const string& curveName,
        double evalTime,
		double fwdStartTime,
        double fwdEndTime,
		double period,
        double resetTime,
        double payTime,
        const ARM_PricingStatesPtr& states) const 
{
	// downcast to get functor
	ARM_PricingFunctionIR* irFunctor = dynamic_cast < ARM_PricingFunctionIR* > (GetNumericalModel());

	if( !irFunctor )
		ARM_THROW(  ERR_INVALID_ARGUMENT, "IR functions are not supported by numerical model" );

	
	// compute standard libor rate from model (payTime = fwdEndTime)
	ARM_VectorPtr liborRate	= irFunctor->Libor(curveName, evalTime,fwdStartTime, 
														 fwdEndTime, period, resetTime, payTime, states );

	// adjust forward (additive, same adjustment in all states)
	// libor adjustment supposed to be stored in first surface of the list
	double adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(FIRST_SURFACE, evalTime, resetTime);
	(*liborRate) *= adjustment;

	return liborRate;
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_Normal_Model
///	Routine: VanillaCaplet
///	Returns: ARM_VectorPtr
///	Action : computes the price of a caplet
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_SLN_Model::VanillaCaplet(
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
	// compute libor rate with previous method (additive adjustement on standard FRA)
	ARM_VectorPtr liborRate = ARM_Local_SLN_Model::Libor(curveName, evalTime,fwdStartTime, 
															 fwdEndTime, fwdPeriod, fwdResetTime, fwdEndTime, states );
	
	// compute discount factor with num model
	ARM_VectorPtr df	= GetNumericalModel()->DiscountFactor( curveName, evalTime, payTime, states );

	// compute option maturity
	double maturity	;
	if (fwdResetTime > evalTime)
		maturity = (fwdResetTime - evalTime) / K_YEAR_LEN;
	else if (fwdResetTime == evalTime)
		maturity = 1./K_YEAR_LEN;
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaCaplet: reset date is prior eval date !" );
	
	
	// get volatility (same vol in all states)
	// vols are supposed to be stored in the first surface of the list
	double volatility	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(FIRST_SURFACE, evalTime, fwdResetTime);
	double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(FIRST_SURFACE, evalTime, fwdResetTime);

	// compute caplet/floorlet for all states
    int i, nbStates = df->size();
    ARM_VectorPtr values (new std::vector<double>(nbStates));

	for (i=0; i<nbStates; i++)
		(*values)[i]  =		payNotional 
						 *	period 
						 *  (*df)[i] 
						 *  BS(	(*liborRate)[i]     + shift, 
								 strikesPerState[i] + shift,
								 maturity,
								 volatility,
								 capFloor );
	
		
				
    return values;
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: VanillaCorridorlet
///	Returns: a vector of Corridor let(t,L(R,S),K,S-E)
///	Action : Closed form formula (call spread) for standard
///          Corridor caplet/floorlet 
////////////////////////////////////////////////////

ARM_VectorPtr ARM_Local_SLN_Model::VanillaCorridorlet(
		const   string& curveName, 
		double  evalTime,
        double  payTime,
        double  resetTime,
        double  startTime,
        double  endTime,
        int     fwdPaymentType, 
        double  fwdPaymentPeriod,
        const std::vector<double>& RefIdxResettimes,
        const std::vector<double>& RefIdxStarttimes,
        const std::vector<double>& RefIdxEndtimes,
        const std::vector<double>& RefFwdPeriods,
        const std::vector<double>& RefIndexWeight,
        double  couponMargin,
        const vector<const std::vector<double>*> downBarrierPerState,
        const vector<const std::vector<double>*> upBarrierPerState,
        double  payNotional,
        int     capFloor,
        const   ARM_PricingStatesPtr& states) const 
{
	//init return values
	ARM_VectorPtr dfPay	= GetNumericalModel()->DiscountFactor( curveName, evalTime, payTime, states );
	int nbStates = dfPay->size();
	ARM_VectorPtr corridor (new std::vector<double>(nbStates));

	//pay index
	ARM_VectorPtr liborPay;
	double shift = 0.0;
	if(fwdPaymentType != IDXFIXED)
    {
		liborPay	= (dynamic_cast < ARM_PricingFunctionIR* > (GetNumericalModel()))->Libor(curveName, evalTime,startTime, 
														 endTime, fwdPaymentPeriod, resetTime, endTime, states );

		double adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(PAY_ADJ_SURFACE, evalTime, resetTime);
		(*liborPay) *= adjustment;

		shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(PAY_IDX_SHIFT_SURFACE, evalTime, resetTime);

		(*liborPay) += couponMargin;
	}
    else
    {
        liborPay = ARM_VectorPtr( new std::vector<double>(nbStates, couponMargin ) );
    }

	//compute proexercise
	std::vector<double> ExerciseProbability(nbStates,0.0); 
	
	for( int k = 0; k < RefIdxResettimes.size(); ++k)
    {
        double ref_resetTime = RefIdxResettimes[k];
        double ref_startTime = RefIdxStarttimes[k];
        double ref_endTime   = RefIdxEndtimes[k];
		double ref_FwdPeriod = RefFwdPeriods[k];

		ARM_VectorPtr ref_Fwd = ARM_Local_SLN_Model::Libor(curveName, evalTime,ref_startTime, 
														 ref_endTime, ref_FwdPeriod, ref_resetTime, ref_endTime, states );
		ARM_VectorPtr ref_Fwd_AdjPay;
		if(fwdPaymentType != IDXFIXED)
		{
			ref_Fwd_AdjPay = ARM_Local_SLN_Model::Libor(curveName, evalTime,ref_startTime, 
															 ref_endTime, ref_FwdPeriod, ref_resetTime, payTime, states );
		}

		double ref_shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(FIRST_SURFACE, evalTime, ref_resetTime);

        double ref_VolDownMinus	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(DOWN_MINUS_SURFACE, evalTime, ref_resetTime);
		double ref_VolDownPlus	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(DOWN_PLUS_SURFACE, evalTime, ref_resetTime);
		double ref_VolUpMinus	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(UP_MINUS_SURFACE, evalTime, ref_resetTime);
		double ref_VolUpPlus	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(UP_PLUS_SURFACE, evalTime, ref_resetTime);

		ARM_VectorPtr DownProbability(new std::vector<double>(nbStates,1.0));
        ARM_VectorPtr UpProbability(new std::vector<double>(nbStates,1.0));

		double ref_VolDownMinus_PayFloat, ref_VolDownPlus_PayFloat, ref_VolUpMinus_PayFloat,ref_VolUpPlus_PayFloat;
		if(fwdPaymentType != IDXFIXED)
		{
			ref_VolDownMinus_PayFloat	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(DOWN_MINUS_SURFACE_PAYFLOAT, evalTime, ref_resetTime);
			ref_VolDownPlus_PayFloat	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(DOWN_PLUS_SURFACE_PAYFLOAT, evalTime, ref_resetTime);
			ref_VolUpMinus_PayFloat	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(UP_MINUS_SURFACE_PAYFLOAT, evalTime, ref_resetTime);
			ref_VolUpPlus_PayFloat	= GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(UP_PLUS_SURFACE_PAYFLOAT, evalTime, ref_resetTime);
		}
		//maturity of the corridorlet
		double maturity	;
		if (ref_resetTime > evalTime)
			maturity = (ref_resetTime - evalTime) / K_YEAR_LEN;
		else if (ref_resetTime == evalTime)
			maturity = 1./K_YEAR_LEN;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "VanillaCorridorlet: reset date is prior eval date !" );
	
		
	    for(int i=0 ; i < nbStates; ++i)
        {
            /// Warning: To be coherent with kernel code
            if((*(downBarrierPerState[i]))[k] > 0.0001)
            {
               	double  strike_DownMinus=(*(downBarrierPerState[i]))[k]-CALL_SPREAD_SHIFT;
				double price_DownMinus = BS(	(*ref_Fwd)[i]+ ref_shift, 
												strike_DownMinus + ref_shift,
												maturity,
												ref_VolDownMinus,
												capFloor );

				double  strike_DownPlus=(*(downBarrierPerState[i]))[k]+CALL_SPREAD_SHIFT;
				double price_DownPlus = BS(	(*ref_Fwd)[i]+ ref_shift, 
												strike_DownPlus + ref_shift,
												maturity,
												ref_VolDownPlus,
												capFloor );

				(*DownProbability)[i] = (price_DownMinus-price_DownPlus)/(2*CALL_SPREAD_SHIFT);
            }

			if(((*(upBarrierPerState[i]))[k] > 0.0001) && ((*(upBarrierPerState[i]))[k] >=(*(downBarrierPerState[i]))[k]) )
			{
				double  strike_UpMinus=(*(upBarrierPerState[i]))[k]-CALL_SPREAD_SHIFT;
				double price_UpMinus = BS(	(*ref_Fwd)[i]+ ref_shift, 
												strike_UpMinus + ref_shift,
												maturity,
												ref_VolUpMinus,
												capFloor );

				double  strike_UpPlus=(*(upBarrierPerState[i]))[k]+CALL_SPREAD_SHIFT;
				double price_UpPlus = BS(	(*ref_Fwd)[i]+ ref_shift, 
												strike_UpPlus + ref_shift,
												maturity,
												ref_VolUpPlus,
												capFloor );

				(*UpProbability)[i] = (price_UpMinus-price_UpPlus)/(2*CALL_SPREAD_SHIFT);
			}
        }
        
        if(fwdPaymentType != IDXFIXED)
        {
			double ref_Fwd_AdjDown; 
			double ref_Fwd_AdjUp; 
			double adjustmentDown = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(DOWN_ADJ_SURFACE, evalTime, ref_resetTime);
			double adjustmentUp = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(UP_ADJ_SURFACE, evalTime, ref_resetTime);
			double adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(FIRST_SURFACE, evalTime, ref_resetTime);
			
			double rTime = resetTime > ref_resetTime ? ref_resetTime : resetTime;
			double toTime = rTime > ref_startTime ? ref_resetTime : rTime;
		   
			for(int i=0 ; i < nbStates; ++i)
			{
				double ref_FwdModel_Adj = (*ref_Fwd)[i]/adjustment;//cov*((*ref_Fwd_AdjPay)[i]/adjustment + ref_shift) - ref_shift;

				/// Warning: To be coherent with kernel code
				if((*(downBarrierPerState[i]))[k] > 0.0001)
				{
					ref_Fwd_AdjDown = ref_FwdModel_Adj*adjustmentDown;//(ref_FwdModel_Adj+ref_shift)*adjustmentDown-ref_shift;
               		double  strike_DownMinus=(*(downBarrierPerState[i]))[k]-CALL_SPREAD_SHIFT;
					double price_DownMinus = BS(	ref_Fwd_AdjDown+ ref_shift, 
													strike_DownMinus + ref_shift,
													maturity,
													ref_VolDownMinus_PayFloat,
													capFloor );

					double  strike_DownPlus=(*(downBarrierPerState[i]))[k]+CALL_SPREAD_SHIFT;
					double price_DownPlus = BS(	ref_Fwd_AdjDown+ ref_shift, 
													strike_DownPlus + ref_shift,
													maturity,
													ref_VolDownPlus_PayFloat,
													capFloor );

					//(*DownProbability)[i] = ((*DownProbability)[i]*(couponMargin-shift)+(price_DownMinus-price_DownPlus)/(2*CALL_SPREAD_SHIFT)*((*liborPay)[i]-couponMargin+shift))/((*liborPay)[i]);
					(*DownProbability)[i] = ((*DownProbability)[i]*(couponMargin)+(price_DownMinus-price_DownPlus)/(2*CALL_SPREAD_SHIFT)*((*liborPay)[i]-couponMargin))/((*liborPay)[i]);
				}

				if(((*(upBarrierPerState[i]))[k] > 0.0001) && ((*(upBarrierPerState[i]))[k] >=(*(downBarrierPerState[i]))[k]) )
				{
					ref_Fwd_AdjUp = ref_FwdModel_Adj*adjustmentUp;//(ref_FwdModel_Adj+ref_shift)*adjustmentUp-ref_shift;
					double  strike_UpMinus=(*(upBarrierPerState[i]))[k]-CALL_SPREAD_SHIFT;
					double price_UpMinus = BS(	ref_Fwd_AdjUp + ref_shift, 
													strike_UpMinus + ref_shift,
													maturity,
													ref_VolUpMinus_PayFloat,
													capFloor );

					double  strike_UpPlus=(*(upBarrierPerState[i]))[k]+CALL_SPREAD_SHIFT;
					double price_UpPlus = BS(	ref_Fwd_AdjUp+ ref_shift, 
													strike_UpPlus + ref_shift,
													maturity,
													ref_VolUpPlus_PayFloat,
													capFloor );

					
					//(*UpProbability)[i] =( (*UpProbability)[i]*(couponMargin-shift)+(price_UpMinus-price_UpPlus)/(2*CALL_SPREAD_SHIFT)*((*liborPay)[i]-couponMargin+shift) )/((*liborPay)[i]) ;
					(*UpProbability)[i] =( (*UpProbability)[i]*(couponMargin)+(price_UpMinus-price_UpPlus)/(2*CALL_SPREAD_SHIFT)*((*liborPay)[i]-couponMargin) )/((*liborPay)[i]) ;
				}
			}
		
		}

		for(i=0 ; i < nbStates; ++i)
		{
			double ref_IndexWeight = RefIndexWeight[k];
            ExerciseProbability[i] += ((*DownProbability)[i] - (*UpProbability)[i])* RefIndexWeight[k];
		}
	}
	for(int i=0; i<nbStates; i++)
        (*corridor)[i]	= payNotional*capFloor*(*dfPay)[i]*((*liborPay)[i])*ExerciseProbability[i];

    return corridor;
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Local_SLN_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Local Shifted LogNormal Model \n";
    os << indent << "-----------------------------\n";
    os << ARM_PricingModel::toString(indent);

	os << indent << "Statistic On Variance squeeze \n";
	int rowNb = itsStatisticCalibResult->GetRowsNb();
	int colNb = itsStatisticCalibResult->GetColsNb();
	
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "EvalDate";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbResetDate";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzDown-";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzDown+";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzUp-";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzUp+";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzPFDown-";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzPFDown+";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzPFUp-";
	os << CC_NS(std,setw)(16)<< CC_NS(std,left)<< "NbVarSqzPFUp+";
	os << "\n";

	for( size_t i=0; i<rowNb; ++i )
	{
        for( size_t j=0; j<colNb; ++j )
        {
            os << CC_NS(std,setw)(16)<< CC_NS(std,fixed) << CC_NS(std,setprecision)(0) << (*itsStatisticCalibResult)(i,j);
        }
		os << "\n";
	}
	return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Local_SLN_Model::ValidateModelParams(const ARM_ModelParams& params) const
{	
	if( !dynamic_cast<const ARM_Local_SLN_ModelParams*>(&params) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_Local_SLN_ModelParams" );
	return true;
}



////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CreateDefaultModelParams
///	Returns : ARM_ModelParams*
///	Action  : static method to build a default ModelParams 
///			  for ARM_Local_SLN_Model 
////////////////////////////////////////////////////
ARM_ModelParams* ARM_Local_SLN_Model::CreateDefaultModelParams()
{
    /// Create a absolute shift model by default
	ARM_ModelParam* fwdAdj  = CreateDefaultModelParam(ARM_ModelParamType::ForwardAdjustment);
	ARM_ModelParam* vol     = CreateDefaultModelParam(ARM_ModelParamType::Volatility);
	ARM_ModelParam* shift   = CreateDefaultModelParam(ARM_ModelParamType::Shift);

	ARM_ModelParamVector params(3);
	params[0] = fwdAdj;
	params[1] = vol;
	params[2] = shift;

    ARM_Local_SLN_ModelParams* slnModelParams = new ARM_Local_SLN_ModelParams(params);

    delete fwdAdj;
    delete vol;
    delete shift;

	return slnModelParams;
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CreateDefaultModelParam
///	Returns : ARM_ModelParam*
///	Action  : static method to build a default Volatility or
///			  ForwardAdjustment ModelParam 
////////////////////////////////////////////////////
ARM_ModelParam* ARM_Local_SLN_Model::CreateDefaultModelParam(ARM_ModelParamType::ParamNb paramType)
{
	ARM_InterpolType type = ARM_InterpolationType::linear_row_extrapoleCst;

    /// Becareful : interpolator is not cloned by the constructor
 	ARM_SurfaceWithInterpol surface (std::vector<double>(0),
									 std::vector<double>(0), 
									 ARM_GP_Matrix(0,0),
									 type,0.0);
	
	/// Build paramType model param : only one surface is used
    ARM_SurfacePtr surfacePtr;
	ARM_SurfacePtrVector surfaceList;
	std::vector<double> index;
	size_t i;
	for(i=0; i<SURFACE_SIZE; i++)
	{		
		surfacePtr =  ARM_SurfacePtr (static_cast <ARM_SurfaceWithInterpol*> (surface.Clone()));
		surfaceList.push_back(surfacePtr);
		index.push_back(i);
	}
		
	return new ARM_SurfaceListModelParam(paramType,index,surfaceList);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : ResetModelParams
///	Returns : void
///	Action  : resets model params (to be used before
///			  calibration)
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::ResetModelParams()
{
	/// Create default params
	ARM_ModelParam* newAdjParam = CreateDefaultModelParam(ARM_ModelParamType::ForwardAdjustment);
	GetModelParams()->SetModelParam(newAdjParam); // SetModelParam() already deletes if necessary
    delete newAdjParam;

    /// Free model params if any
	GetModelParams()->DeleteModelParam(ARM_ModelParamType::Volatility);
	GetModelParams()->DeleteModelParam(ARM_ModelParamType::Shift);
	GetModelParams()->DeleteModelParam(ARM_ModelParamType::QVol);
	GetModelParams()->DeleteModelParam(ARM_ModelParamType::QParameter);

    /// Create ans set new ones
    ARM_ModelParam *newVolParam,*newShiftParam;
    ARM_PricingModel* numericalModel = GetNumericalModel();
    if( numericalModel->GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Volatility) &&
        numericalModel->GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Shift))
    {
        newVolParam     = CreateDefaultModelParam(ARM_ModelParamType::Volatility);
        newShiftParam   = CreateDefaultModelParam(ARM_ModelParamType::Shift);
    }
    else if( numericalModel->GetModelParams()->DoesModelParamExist(ARM_ModelParamType::QVol) &&
             numericalModel->GetModelParams()->DoesModelParamExist(ARM_ModelParamType::QParameter))
    {
        newVolParam     = CreateDefaultModelParam(ARM_ModelParamType::QVol);
        newShiftParam   = CreateDefaultModelParam(ARM_ModelParamType::QParameter);
    }
    else
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_Local_SLN_Model only runs absolute or relative shifts (Q) models" );

	GetModelParams()->SetModelParam(newVolParam);
    delete newVolParam;
	GetModelParams()->SetModelParam(newShiftParam);
    delete newShiftParam;
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: Forward
///	Returns: ARM_VectorPtr
///	Action : Computes the forward
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_SLN_Model::Forward(
	const string& modelName, 
    double evalTime,
	double expiryTime,
	double settlementTime,
	double payTime,
    const ARM_PricingStatesPtr& states) const
{
	/// Get the numerical model and check its type
	ARM_PricingFunctionEquity* eqFxModel = dynamic_cast < ARM_PricingFunctionEquity* > (GetNumericalModel());

	/// Delegate the forward computation to its numerical model
	ARM_VectorPtr forward = eqFxModel->Forward(modelName,evalTime,expiryTime,settlementTime,payTime,states);

	/// Adjust forward in a deterministic multiplicative way
	/// Adjustment is stored in the first surface of the list
	double adjustment = GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).GetValue(FIRST_SURFACE,evalTime,expiryTime);
	(*forward) *= adjustment;

	return forward;
}


////////////////////////////////////////////////////
///	Class  : ARM_Local_SLN_Model
///	Routine: CallVectorial
///	Returns: ARM_VectorPtr
///	Action : Computes the price of a call/put option
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_SLN_Model::CallVectorial(
	const string& modelName,
    double evalTime,
	double expiryTime,
	double settlementTime,
	const std::vector<double>& strikesPerState,
    int callPut,
	double payTime,
	const ARM_PricingStatesPtr& states,
    ARM_PricingContext* context) const
{	
	/// Compute the forward (a local adjustment may be used)
	ARM_VectorPtr forward = Forward(modelName,evalTime,expiryTime,settlementTime,payTime,states);
	
	/// Compute the payment discount factor with its numerical model
    ARM_PricingModel* numericalModel = GetNumericalModel();
	ARM_VectorPtr zcPay = numericalModel->DiscountFactor(modelName,evalTime,payTime,states);

	/// Compute option maturity
	double maturity	;
	if (expiryTime > evalTime)
		maturity = (expiryTime - evalTime) / K_YEAR_LEN;
	else if (expiryTime == evalTime)
        maturity = 1./K_YEAR_LEN; /// for calibration purpose !
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " Fx reset date is prior eval date !" );
	
	/// Get the local determinsitic SLN volatility and shift
	/// Vols and shifts are supposed to be stored in the first surface of the list
    bool isAbsShift = GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Shift);

	// Compute call/put prices
    size_t i, nbStates = zcPay->size();
    ARM_VectorPtr values (new std::vector<double>(nbStates));

    if(isAbsShift)
    {
        /// Shifted lognormal pricing formula (<=> SLN & absolute shift)
	    double volatility = GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).GetValue(FIRST_SURFACE, evalTime, expiryTime);
	    double shift = GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).GetValue(FIRST_SURFACE, evalTime, expiryTime);
	    for (i=0; i<nbStates; i++)
		    (*values)[i]  =	(*zcPay)[i] *
                            BS((*forward)[i] + shift, 
							    strikesPerState[i] + shift,
							    maturity,
							    volatility,
							    callPut);
    }
    else
    {
        /// Q pricing formula (<=> SLN & relative shift)
	    double qvol = GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).GetValue(FIRST_SURFACE, evalTime, expiryTime);
	    double qparam = GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter).GetValue(FIRST_SURFACE, evalTime, expiryTime);
        const ARM_SurfaceListModelParam& qVolParam = GetModelParams()->GetModelParam(ARM_ModelParamType::QVol).ToSurfaceListModelParam();
        double fwd;
        if(qparam != 1.0 && qVolParam.cols(TARGET_PRICE_SURFACE)>0)
        {
            ARM_PricingStatesPtr dumStates(NULL);
	        double forward0 = (*(Forward(modelName,0.0,expiryTime,settlementTime,payTime,dumStates)))[0];
            if(fabs(qparam)<=K_NEW_DOUBLE_TOL)
                qparam = (qparam > 0 ? K_NEW_DOUBLE_TOL : -K_NEW_DOUBLE_TOL);

            /// Try a numerical calibration of the local Q vol 
		    int timeIdx		= numericalModel->GetNumMethod()->GetLastTimeIdx();
		    ARM_GP_VectorPtr arrowDebreuPrices = numericalModel->GetNumMethod()->GetArrowDebreuPrices( timeIdx, *numericalModel);

            if(arrowDebreuPrices == ARM_GP_VectorPtr(NULL))
		        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " Arrow-Debreu prices are missing for local Q model numerical calibration" );

            double targetPrice = qVolParam.GetValue(TARGET_PRICE_SURFACE, 0.0, expiryTime);

            /// Use a Newtow-Raphson to get closer
            double price,priceUp,newPrice,modelPrice;
            double volShift=0.002; // 0.2% de Q vol
            double volMove,maxVolMove;
            size_t iter=0,maxi=10;
            newPrice = 0.0;
	        for(i=0; i<nbStates; i++)
            {
                fwd = (*forward)[i];
		        modelPrice = QModelAnalytics::BSQFunction( fwd,
                                strikesPerState[i], qvol, maturity, qparam,
                                (*zcPay)[i], callPut, /*forward0*/ fwd );
                newPrice += (*arrowDebreuPrices)[i] * modelPrice;
                (*values)[i] = modelPrice;
            }
/****
FILE* f=fopen("c:\\temp\\dumpSLN.txt","a");
fprintf(f,"t=%7.2lf\tT=%7.2lf\ttarget=%15.10lf\tprice=%15.10lf\tqVol=%15.10lf\tqParam=%15.10lf\n",evalTime,expiryTime,targetPrice,newPrice,qvol,qparam);
fclose(f);
****/
			double PVAdj=1.0;
            while(iter<maxi && fabs(newPrice-targetPrice)>0.01)
            {
                priceUp = 0.0;
	            for(i=0; i<nbStates; i++)
                {
                    fwd = (*forward)[i];
		            modelPrice = QModelAnalytics::BSQFunction( fwd,
                                    strikesPerState[i], qvol+volShift, maturity, qparam,
                                    (*zcPay)[i], callPut, /*forward0*/ fwd );
                    priceUp += (*arrowDebreuPrices)[i] * modelPrice;
                }
            
                price   = newPrice;
                if(fabs(priceUp-price)>K_NEW_DOUBLE_TOL)
                {
                    volMove = -(price - targetPrice)/(priceUp-price)*volShift;
                    maxVolMove = 0.5*qvol;
                    if(fabs(volMove)>maxVolMove)
                        qvol += (volMove>0 ? maxVolMove : -maxVolMove);
                    else
                        qvol += volMove;
                }
                else
                    // No vega then quit
                    break;

                newPrice = 0.0;
	            for(i=0; i<nbStates; i++)
                {
                    fwd = (*forward)[i];
		            modelPrice = QModelAnalytics::BSQFunction( fwd,
                                    strikesPerState[i], qvol, maturity, qparam,
                                    (*zcPay)[i], callPut, /*forward0*/ fwd );
                    newPrice += (*arrowDebreuPrices)[i] * modelPrice;
                    (*values)[i] = modelPrice;
                }
				PVAdj = targetPrice/newPrice;

                ++iter;
/****
FILE* f=fopen("c:\\temp\\dumpSLN.txt","a");
fprintf(f,"     #%1d\tprice=%15.10lf\tpriceShift=%15.10lf\tnewprice=%15.10lf\n",iter,price,priceUp-price,newPrice);
fclose(f);
****/
            }
/**** Version with a final PV adjuster
			if(PVAdj!=1.0)
			{
				newPrice = 0.0;
				for(i=0; i<nbStates; i++)
				{
					(*values)[i] *= PVAdj;
                    newPrice += (*arrowDebreuPrices)[i] * (*values)[i];
				}
FILE* f=fopen("c:\\temp\\dumpSLN.txt","a");
fprintf(f,"--------------> PV Adjsuted price=%15.10lf\n",newPrice);
fclose(f);
			}
****/

            const_cast< ARM_SurfaceListModelParam& >(qVolParam).SetValue(NUM_CALIB_SURFACE, evalTime, expiryTime,qvol);
        }
        else
        {
            /// Use calibrated Q vol (analytical or numerical as computed above)
	        for (i=0; i<nbStates; i++)
            {
                fwd = (*forward)[i];
		        (*values)[i]  =	QModelAnalytics::BSQFunction( fwd,
                                    strikesPerState[i], qvol, maturity, qparam,
                                    (*zcPay)[i], callPut, fwd );
            }
        }
    }

			
    return values;
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : GetSettlementCalendar
///	Returns : string
///	Action  : Delegate to equity/fx type numerical model
////////////////////////////////////////////////////
string ARM_Local_SLN_Model::GetSettlementCalendar(const string& modelName) const
{
	ARM_PricingFunctionEquity* eqFxModel = dynamic_cast< ARM_PricingFunctionEquity* >(GetNumericalModel());
    if(!eqFxModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for GetSettlementCalendar()");

    return eqFxModel->GetSettlementCalendar(modelName);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : GetSettlementGap
///	Returns : double
///	Action  : Delegate to equity/fx type numerical model
////////////////////////////////////////////////////
double ARM_Local_SLN_Model::GetSettlementGap(const string& modelName) const
{
	ARM_PricingFunctionEquity* eqFxModel = dynamic_cast< ARM_PricingFunctionEquity* >(GetNumericalModel());
    if(!eqFxModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Equity or Fx compatible model required for GetSettlementGap()");

    return eqFxModel->GetSettlementGap(modelName);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration 
///			  (for 1 given ARM_Security)
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateLocalModel(const ARM_Security& security, double targetPrice, const std::vector<double>& evalTimes, size_t secIdx)
{	//
	// Switch on security type
	//
    ARM_Option* option=NULL;
    ARM_CapFloor* capFloor=NULL;
	ARM_CorridorLeg* corridor=NULL;
	if( (option = dynamic_cast< ARM_Option* >((ARM_Security*)&security)) )
	{			
		CalibrateEqFxOption(option,targetPrice,evalTimes);
	}
	else if( (capFloor = dynamic_cast< ARM_CapFloor* >((ARM_Security*)&security)) )
	{			
		CalibrateCapFloor(capFloor,targetPrice,evalTimes);
	}
	else if( (corridor = dynamic_cast< ARM_CorridorLeg* >((ARM_Security*)&security)) )
	{			
		CalibrateCorridorLet(corridor,targetPrice,evalTimes);
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local SLN Model calibration: instrument not supported" );
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration (given a calculator -> computation of target prices)
///				used for FX and Libor dual condition
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateLocalModel(const ARM_GenCalculator& calculator, const std::vector<double>& evalTimes) 
{
	if(ARM_FXRACalculator* fxraCalculator = dynamic_cast<ARM_FXRACalculator*>((ARM_GenCalculator*)&calculator))
	{
		fxraCalculator->Price();
		const ARM_Vector* allFixingPrices = fxraCalculator->GetvAllFixingPrices();

		const size_t EVALTIME=1;
		const size_t FXPRICEBDOWN=2;
		const size_t FXPROBABDOWN=3;
		const size_t FXPRICEBUP=4;
		const size_t FXPROBABUP=5;
		const size_t IRPRICEBDOWN=6;
		const size_t IRPROBABDOWN=7;
		const size_t IRPRICEBUP=8;
		const size_t IRPROBABUP=9;
		const size_t TOTAL=9;

		size_t LOCALBUPFX	=0;
		size_t LOCALBDOWNFX	=1;
		size_t LOCALBUPIR	=2;
		size_t LOCALBDOWNIR	=3;


		ARM_PricingModel* numericalModel = GetNumericalModel();
		ARM_QModel1F_Fx* eqFxModel = dynamic_cast< ARM_QModel1F_Fx* >( numericalModel );
		
		//Get the barriers from the calculator
		ARM_Curve fxDownBarrierCv	= fxraCalculator->GetFXDownBarrier();
		ARM_Curve fxUpBarrierCv		= fxraCalculator->GetFXUpBarrier();
		ARM_Curve irDownBarrierCv	= fxraCalculator->GetIRDownBarrier();
		ARM_Curve irUpBarrierCv		= fxraCalculator->GetIRUpBarrier();

		char* ccyName				= fxraCalculator->GetPayCcyName();

		string irIdxIT = fxraCalculator->GetIRIdxIT();
		double irIndexTerm		= StringMaturityToYearTerm(irIdxIT);
		string irIndexTermStr	= ConvertYearTermToStringMatu(irIndexTerm);
		ARM_INDEX_TYPE irIndexType;
		irIndexType = (ARM_INDEX_TYPE)FromIndexAndTermToIndexType(irIndexTermStr, ccyName);
		ARM_IRIndex irIndex(irIndexType);

		//Get the dates from the calculator using the date strip
		ARM_DateStripPtr dateStrip		= fxraCalculator->GetDateStrip() ;
		std::vector<double>& startDates		= dateStrip->GetFwdRateStartDates();
		std::vector<double>& endDates			= dateStrip->GetFwdRateEndDates();
		std::vector<double>& payDates			= dateStrip->GetPaymentDates();
		std::vector<double>& interestTerms	= dateStrip->GetInterestTerms() ;

		size_t nbEval = evalTimes.size();

		bool isFXDown=true, isIRDown=true; 
		bool isFXUpSqueezeVol=false, isFXDownSqueezeVol=false, isIRUpSqueezeVol=false, isIRDownSqueezeVol=false;

		//Target values : prices and barriers computed in the calculator
		double pricefxBDown, pricefxBUp, priceirBDown, priceirBUp;
		double fxBDownTgt, fxBUpTgt, irBDownTgt, irBUpTgt;

		double volfxBDownTgt, volfxBUpTgt, volirBDownTgt, volirBUpTgt;
		double varfxBDownTgt, varfxBUpTgt, varirBDownTgt, varirBUpTgt;
		double varfxModel, varirModel;
		double fxDownBarrier, fxUpBarrier, irDownBarrier, irUpBarrier;

		double irfwdStartTime, irfwdEndTime;
		
		double localBUpvarFX, localBDownvarFX, localBUpvarIR, localBDownvarIR;
		double localBUpvolFX, localBDownvolFX, localBUpvolIR, localBDownvolIR;
		double localBUpFX, localBDownFX, localBUpIR, localBDownIR;
		
		size_t size = allFixingPrices->size();

		double fixingTime;
		size_t nbFixIdx=0;
		size_t nbFixings=0; 
		
		for (size_t i=0; i<nbEval; i++)
		{	
			double evalTime = evalTimes[i];
				
			ARM_Date startDate = (*startDates)[i];
				
			ARM_Date endDate = (*endDates)[i];
			ARM_Date payDate = (*payDates)[i];
			fxraCalculator->FixingsStructure(startDate, endDate, ccyName, irIndex, irIndexTermStr, numericalModel);

			irDownBarrier	= irDownBarrierCv.Interpolate(evalTime);
			irUpBarrier		= irUpBarrierCv.Interpolate(evalTime);

			fxDownBarrier	= fxDownBarrierCv.Interpolate(evalTime);
			if (fxDownBarrier <= 0) isFXDown = false;
			fxUpBarrier		= fxUpBarrierCv.Interpolate(evalTime);

			//les dates sont calculées pour l'exercice courant par la méthode FixingsStructure()
			std::vector<double> irIndexStartTimes		= fxraCalculator->GetIRIndexStartTimes() ;
			std::vector<double> irIndexEndTimes		= fxraCalculator->GetIRIndexEndTimes() ;
			std::vector<double> irIndexTerms			= fxraCalculator->GetIRIndexTerms() ;
				
			if (i != 0 )
				nbFixIdx += nbFixings*TOTAL + 1;

			nbFixings=(*allFixingPrices)[nbFixIdx];
			for (size_t p=0; p<nbFixings; p++)
			{
				
				fixingTime = (*allFixingPrices)[nbFixIdx + p*TOTAL + EVALTIME];

				// IR PART
				double fwdstartTime = irIndexStartTimes[p];
				double fwdendTime = irIndexEndTimes[p];
				double fwdPeriod = irIndexTerms[p];
				double dfStart	= eqFxModel->GetZeroCurve()->DiscountPrice(fwdstartTime / K_YEAR_LEN);
				double dfEnd	= eqFxModel->GetZeroCurve()->DiscountPrice(fwdendTime / K_YEAR_LEN);
				double fwdLibor = (dfStart / dfEnd - 1.0) / fwdPeriod;

				//Model variances until evalTime
				//	-> volatility of fwd ZC
				ARM_VectorPtr fwdZCMP_Eval_Start = eqFxModel->ComputeFwdZCModelParam(0, evalTime, fwdstartTime, fwdendTime);
				double volirModel = (*fwdZCMP_Eval_Start)[0];
				
				//	-> transform to volatility of Libor (H&W framework) : approx normal vol of Libor
				volirModel =  (1./fwdPeriod) * volirModel * pow((dfStart / dfEnd),2);
				varirModel = volirModel * volirModel * evalTime / K_YEAR_LEN;

				ARM_VectorPtr fwdZCMP_Start_Start = eqFxModel->ComputeFwdZCModelParam(0, fwdstartTime, fwdstartTime, fwdendTime);
				double totalvolirModel = (*fwdZCMP_Start_Start)[0];

				totalvolirModel =  (1./fwdPeriod) * totalvolirModel * pow((dfStart / dfEnd), 2);

				localBDownvolIR = 0.0;
				localBDownIR = 0.0;
				if (isIRDown)
				{
					//target prices 
					priceirBDown		= (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPRICEBDOWN];

					try
					{
						//Compute target normal volatilities associated with theses prices
						volirBDownTgt	= VanillaImpliedVol_N(fwdLibor, priceirBDown, irDownBarrier, fixingTime/K_YEAR_LEN, 1);

						//Compute target normal variances until fixingTime
						varirBDownTgt	= volirBDownTgt * volirBDownTgt * fixingTime/K_YEAR_LEN; 

						localBDownvarIR = varirBDownTgt - varirModel;
					
						if(localBDownvarIR>0.0)
							localBDownvolIR = sqrt( localBDownvarIR / ((fixingTime - evalTime) / K_YEAR_LEN));
						else // Variance squeeze
							isIRDownSqueezeVol=true;
					}
					catch(Exception& x)
					{
						localBDownvolIR = 0.0;
						isIRDownSqueezeVol=true;
					}

					//Local Barrier Down calibrated on probas
					irBDownTgt = (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPROBABDOWN];
					double totalvarBDown=0.0;
					if (isIRDownSqueezeVol)
						totalvarBDown = totalvolirModel * totalvolirModel * fixingTime / K_YEAR_LEN;
					else 
						totalvarBDown = varirModel + localBDownvarIR;

					localBDownIR = fwdLibor + irBDownTgt * sqrt(totalvarBDown);	
					
				}

				try
				{
					priceirBUp			= (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPRICEBUP];
					volirBUpTgt		= VanillaImpliedVol_N(fwdLibor, priceirBUp, irUpBarrier, fixingTime/K_YEAR_LEN, 1, &totalvolirModel);
					//double test = VanillaImpliedVol_BS	(fwdLibor, irUpBarrier, fixingTime/K_YEAR_LEN, 
				//							             priceirBUp,1);
					varirBUpTgt		= volirBUpTgt * volirBUpTgt * fixingTime/K_YEAR_LEN;

					localBUpvarIR = varirBUpTgt - varirModel;
					localBUpvolIR = 0.0;
					if(localBUpvarIR>0.0)
						localBUpvolIR = sqrt( localBUpvarIR / ((fixingTime - evalTime) / K_YEAR_LEN));
					else // Variance squeeze
						isIRUpSqueezeVol=true;
				}
				catch(Exception& x)
				{
					localBUpvolIR = 0.0;
					isIRUpSqueezeVol=true;
				}

				//Local Barrier Down calibrated on probas
				irBUpTgt = (*allFixingPrices)[nbFixIdx + p*TOTAL + IRPROBABUP];
				double totalvarBUp=0.0;
				if (isIRUpSqueezeVol)
					totalvarBUp = totalvolirModel * totalvolirModel * fixingTime / K_YEAR_LEN;
				else 
					totalvarBUp = varirModel + localBUpvarIR;
					
				localBUpIR = fwdLibor + irBUpTgt * sqrt(totalvarBUp);
				

				//FX PART
				ARM_VectorPtr fx = eqFxModel->Forward(eqFxModel->GetModelName(), 0, fixingTime, fixingTime, fixingTime, ARM_PricingStatesPtr( new ARM_PricingStates(1,1,0) ));
				double fwdFX = (*fx)[0];

				ARM_VectorPtr fwdFXMP_Eval_Fix= eqFxModel->ComputeFwdFXModelParam(0, evalTime, fixingTime);
				double volfxModel = (*fwdFXMP_Eval_Fix)[1];

				ARM_VectorPtr fwdFXMP_Fix_Fix = eqFxModel->ComputeFwdFXModelParam(0, fixingTime, fixingTime);
				double totalvolfxModel = (*fwdFXMP_Fix_Fix)[1];

				varfxModel = volfxModel * volfxModel * evalTime/K_YEAR_LEN;

				//test
				ARM_VectorPtr t= eqFxModel->ComputeFwdFXModelParam(evalTime, fixingTime, fixingTime);
				double tt = (*t)[1];
				double stdtt = tt * sqrt((fixingTime - evalTime)/K_YEAR_LEN);
				double std1 = sqrt(varfxModel);
				double std2 = totalvolfxModel * sqrt(fixingTime/K_YEAR_LEN);

				//fin test
								
				localBDownvolFX = 0.0;
				localBDownFX = 0.0;
				if(isFXDown)
				{
					//Local Volatility calibrated on prices
					
					//target prices 
					pricefxBDown		= (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPRICEBDOWN];
						
					
					try
					{
						//Compute target normal volatilities associated with theses prices
						volfxBDownTgt	= VanillaImpliedVol_BS	(fwdFX, fxDownBarrier, fixingTime/K_YEAR_LEN, 
											             pricefxBDown,1);
					
						//Compute target lognormal variances until fixingTime
						varfxBDownTgt = volfxBDownTgt * volfxBDownTgt * fixingTime/K_YEAR_LEN;
					
						localBDownvarFX = varfxBDownTgt - varfxModel;
					
						if(localBDownvarFX>0.0)
							localBDownvolFX = sqrt(localBDownvarFX / ((fixingTime - evalTime) / K_YEAR_LEN));
						else // Variance squeeze
						{
							isFXDownSqueezeVol=true;
						}	
					}
					catch(Exception& x)
					{
						localBDownvolFX = 0.0;
						isFXDownSqueezeVol=true;

					}


					//Local Barrier Down calibrated on probas
					fxBDownTgt = (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPROBABDOWN];
					double totalvarBDown=0.0;
					if (isFXDownSqueezeVol)
						totalvarBDown = totalvolfxModel * totalvolfxModel * fixingTime / K_YEAR_LEN;
					else 
						totalvarBDown = varfxModel + localBDownvarFX;

					localBDownFX = fwdFX * exp(fxBDownTgt * sqrt(totalvarBDown) - 0.5 * totalvarBDown);	

				}

				pricefxBUp			= (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPRICEBUP];
				
				try
				{
					volfxBUpTgt	= VanillaImpliedVol_BS(fwdFX, fxUpBarrier, fixingTime/K_YEAR_LEN, pricefxBUp, 1);
					varfxBUpTgt	= volfxBUpTgt * volfxBUpTgt * fixingTime/K_YEAR_LEN;
					localBUpvarFX = varfxBUpTgt - varfxModel;

					localBUpvolFX = 0.0;
					if(localBUpvarFX>0.0)
						localBUpvolFX = sqrt(localBUpvarFX / ((fixingTime - evalTime) / K_YEAR_LEN));
					else // Variance squeeze
					{
						isFXUpSqueezeVol=true;
					}
					

				}
				
				catch(Exception& x)
				{
					//Si exception, pas de calibration de la vol locale
					localBUpvolFX = 0.0;
					isFXUpSqueezeVol=true;				
				}
				
				double test = localBUpvolFX*sqrt((fixingTime - evalTime) / K_YEAR_LEN);
				
				//Local Barrier Up calibrated on probas
				fxBUpTgt = (*allFixingPrices)[nbFixIdx + p*TOTAL + FXPROBABUP];
				totalvarBUp=0.0;
				if (isFXUpSqueezeVol)
					totalvarBUp = totalvolfxModel * totalvolfxModel * fixingTime / K_YEAR_LEN;
				else 
					totalvarBUp = varfxModel + localBUpvarFX;

				localBUpFX = fwdFX * exp(fxBUpTgt * sqrt(totalvarBUp) - 0.5 * totalvarBUp);

				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(LOCALBUPFX,evalTime,fixingTime,localBUpvolFX);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(LOCALBDOWNFX,evalTime,fixingTime,localBDownvolFX);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(LOCALBUPIR,evalTime,fixingTime,localBUpvolIR);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(LOCALBDOWNIR,evalTime,fixingTime,localBDownvolIR);

				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(LOCALBUPFX,evalTime,fixingTime,localBUpFX);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(LOCALBDOWNFX,evalTime,fixingTime,localBDownFX);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(LOCALBUPIR,evalTime,fixingTime,localBUpIR);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(LOCALBDOWNIR,evalTime,fixingTime,localBDownIR);

			}
				
		}	

		
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model calibration: calculator not supported" );
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateEqFxOption
///	Returns : void
///	Action  : Local model calibration to a Equity or
///			  Fx vanilla option
///           Numerical model is assumed to be
///           SLN compatible
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateEqFxOption(ARM_Option* option, double targetPrice, const std::vector<double>& evalTimes)
{	
	size_t i, sizeEvalTimes = evalTimes.size();

    /// Collect datas of ARM_Option
	double resetTime,settlementTime,payTime,strike;
	int callPut;
	GetEqFxOptionData(option,resetTime,settlementTime,payTime,strike,callPut);

    double invYearLen = 1.0/K_YEAR_LEN;
	double tFix	 = resetTime * invYearLen ;

	/// Calibration only if option maturity is after evalTime
    ARM_PricingModel* numericalModel = GetNumericalModel();
	ARM_PricingFunctionEquity* eqFxModel = dynamic_cast< ARM_PricingFunctionEquity* >( numericalModel );

    /// Absolute or relative shift model
    ARM_ModelParamType::ParamNb volType    = ARM_ModelParamType::Volatility;
    ARM_ModelParamType::ParamNb shiftType  = ARM_ModelParamType::Shift;
    bool isAbsShift = GetModelParams()->DoesModelParamExist(ARM_ModelParamType::Shift);
    if(!isAbsShift)
    {
        volType     = ARM_ModelParamType::QVol;
        shiftType   = ARM_ModelParamType::QParameter;

        /// Set targetPrice for further numerical calibration
        GetModelParams()->GetModelParam(volType).SetValue(TARGET_PRICE_SURFACE, 0.0, resetTime, targetPrice);
    }

	/// Get market forward rate (potentially adjusted if payment is not standard)
    /// and market (i.e. B&S vol) forward Fx vol previously computed
	double targetFwd = option->GetCalcFwd();
	double targetVol = 0.01*option->GetCalcVol();

    /// Get unit sensi B&S target price
    targetPrice /= option->GetCalcO1();
	
    std::vector<double> strikesPerState(1,strike);
    ARM_PricingCallContext callContext;
    double fwdAdj,volAdj,shift,modelVol,tCall,var,stdDev;
    ARM_GP_VectorPtr forward;
    bool success;
	for(i=0; i<sizeEvalTimes; i++)
	{
        fwdAdj=1.0;
		volAdj=0.0;
        shift = (isAbsShift ? 0.0 : 1.0);
        if(evalTimes[i] <= resetTime)
        {
		    tCall = evalTimes[i] * invYearLen;
            stdDev = sqrt(tFix-tCall)*targetVol*targetFwd; /// proxy of normal vol is LNVol*Fwd
            if(fabs(targetFwd-strike) < NB_STDDEV_NO_CALIB_EQFX*stdDev)
            {
            
		        /// Forward	& volatility adjusment + shift :
                /// price the call then restore pricing context to get datas
				// NDC
				/*
	            ARM_VectorPtr fxOption = eqFxModel->CallVectorial(
                    "",0.0,
                    evalTimes[i],
                    settlementTime,
                    strikesPerState,
                    callPut,
                    payTime,
                    ARM_PricingStatesPtr(NULL),
                    &callContext);
				*/
				size_t statesSize = numericalModel->ModelStatesSize();
				ARM_VectorPtr fxOption = eqFxModel->CallVectorial(
                    "",0.0,
                    std::vector<double>(statesSize,evalTimes[i]),
                    std::vector<double>(statesSize,settlementTime),
                    strikesPerState,
                    callPut,
                    std::vector<double>(statesSize,payTime),
                    ARM_PricingStatesPtr(NULL),
                    &callContext);

                forward = callContext.GetForward();
                fwdAdj  = targetFwd / (*forward)[0];


		        modelVol    = callContext.GetVol();
			    shift       = callContext.GetShift(); /// absolute or relative (Q)

	            if( ( isAbsShift && shift != 0.0) ||
                    (!isAbsShift && shift != 1.0) )
                {
                    /// The B&S targetVol must be converted to be SLN compatible
                    success = true;
                    if(isAbsShift)
                    {
                        /// Convert targetVol in an absolute shift world
	                    targetVol = VanillaImpliedVol_BS(targetFwd+shift,strike+shift,tFix,
											             targetPrice,callPut,NULL,&success);
                        if(!success)
                            /// Get a proxy for SLN vol
                            targetVol *= targetFwd/(targetFwd+shift);
                    }
                    else
                    {
                        /// Convert targetVol in a relative shift (Q) world
                        /// The standard call to VanillaImpliedVol_BS() is not used
                        //  to allow negative Q (handled by QPriceFunct object)
                        QPriceFunct qPriceFct(targetFwd,strike,tFix,shift,targetFwd,callPut);
                        UnaryFuncWithNumDerivative<double,double> fctToSolve( qPriceFct );

                        T_SmoothNewtonRaphsonSolver< UnaryFuncWithNumDerivative<double,double> > solver(fctToSolve,targetPrice,1.0e-8);
	                    solver.setInitialGuess(targetVol);

                        targetVol=solver.Solve();
                    }
                }

                if(evalTimes[i] == resetTime)
                    /// Shift by an artificial 1d to allow vol correction
                    tFix += invYearLen;

		        var = targetVol*targetVol*tFix - modelVol*modelVol*tCall;
		        if( var>=0 )
			        volAdj = sqrt( var / (tFix - tCall) );

		        else 
		        {
					SetVarianceSqueezeStatus(true);
					/// Variance squeeze : put a warning message
			        CC_Ostringstream os;
			        os << ARM_USERNAME << " : Local_SLN_Model & Equity/Fx option Calibration : variance squeeze [call time = " << evalTimes[i] << ", reset time = " << resetTime << " ]"<< endl;
			        ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					volAdj = 0.01; // 1% of qVol
		        }
            }
            else
            {
                /// No calibration => keep LN ATM vol
                volAdj = targetVol;

                /// Use shift interpolated at evalTime
                shift = numericalModel->GetModelParams()->GetModelParam(shiftType).ToCurveModelParam().GetCurve()->Interpolate(evalTimes[i]);

            }
        }

		// Set forward adjustment, volatility & shift parameter to model params
		GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(FIRST_SURFACE,evalTimes[i],resetTime,fwdAdj);
		GetModelParams()->GetModelParam(volType).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, volAdj);
		GetModelParams()->GetModelParam(shiftType).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, shift);

		GetModelParams()->GetModelParam(volType).SetValue(NUM_CALIB_SURFACE, evalTimes[i], resetTime, 0.0);
	}
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateCapFloor
///	Returns : void
///	Action  : Local model calibration to a 
///			  cap / floor
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateCapFloor(ARM_CapFloor* capFloor, double targetPrice, const std::vector<double>& evalTimes)
{	
	/// target price is ignored because we need the pv of all caplets / floorlets
	/// --> use GetCashFlowValues()
	size_t capSize = capFloor->GetResetDates()->GetSize();
	ARM_Vector* cashFlows = capFloor->GetCashFlowValues();

	for (size_t i(0); i<capSize; i++)
		CalibrateCapFloorLet(capFloor, i, cashFlows->Elt(i), evalTimes);
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateCapFloorLet
///	Returns : void
///	Action  : Local model calibration to a 
///			  caplet / floorlet
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateCapFloorLet(ARM_CapFloor* capFloor, size_t periodIdx, double targetPrice, const std::vector<double>& evalTimes)
{	
    // will generate an exception if numerical model is not supported for calibration
	Local_SLN_Model_Calibration_Helper helper (GetNumericalModel());

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
	double shift = helper.LiborSLNShiftFromNumericalModel(0.0,
															resetTime,
															fwdStartTime, 
															fwdEndTime, 
															fwdPeriod);
	targetVol = VanillaImpliedVol_BS (	liborRateTarget+shift, 
											strike+shift, 
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
																		resetTime,
																		fwdStartTime, 
																		fwdEndTime, 
																		fwdPeriod, 
																		strike);
				
			// this adjust will be = 0 if standard libor
			double adjustment  = liborRateTarget/liborRateModel;
						

			// set adjustments to model params
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, adjustment);
							
			GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, shift);

			/// compute libor forward volatility
			/// if the implied vol computation has failed, this means that the strike is too far from ATM
			/// in that case, a forward vol = 0 is inputed
			///
			double forwardVol;
			
			if (success)
			{
				double modelVol; /// model spread vol up to evalTimes[i] 
			
				modelVol = helper.LiborSLNVolatilityFromNumericalModel(	evalTimes[i],
																		resetTime,
																		fwdStartTime, 
																		fwdEndTime, 
																		fwdPeriod, 
																		strike,
																		shift);
							
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

			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, forwardVol);
		}
		else
		{	GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, 1.0);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, 0.0);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(FIRST_SURFACE, evalTimes[i], resetTime, 0.0);
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_SLN_Model
///	Routine : CalibrateCorridorLet
///	Returns : void
///	Action  : Local model calibration to a 
///			  corridorlet / floorlet
////////////////////////////////////////////////////
void ARM_Local_SLN_Model::CalibrateCorridorLet(ARM_CorridorLeg* corridor, double targetPrice, const std::vector<double>& evalTimes)
{	
    // will generate an exception if numerical model is not supported for calibration
	Local_SLN_Model_Calibration_Helper helper (GetNumericalModel());
	double AsOfDate = GetAsOfDate().GetJulian();
	
	size_t i, sizeEvalTimes = evalTimes.size();
	// some validation
	if (corridor->GetResetDates()->size() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Normal Model Calibration : ARM_corridor is required to have only one optionlet" );

	double resetTime, payTime, fwdPaymentPeriod, couponMargin, payNotional;
	std::vector<double> refIdxResettimes, refIdxStarttimes,	refIdxEndtimes,	refFwdPeriods,refIndexWeight, DownBarrier, UpBarrier;
	int callPut,indexPaymentType;

	if(itsResizeStatisticResult)
	{
		itsStatisticCalibResult->resize(sizeEvalTimes,COL_NB_STATIC_RESULTS);
		itsResizeStatisticResult = false;
	}

	GetCorridorLegData	(corridor,
						callPut,
						resetTime,
						payTime,
						indexPaymentType,
						fwdPaymentPeriod,
						refIdxResettimes ,
						refIdxStarttimes,
						refIdxEndtimes,
						refFwdPeriods,
						refIndexWeight,
						DownBarrier,
						UpBarrier,
						couponMargin,
						payNotional);

	double dfPay = GetNumericalModel()->GetZeroCurve()->DiscountPrice(payTime/K_YEAR_LEN);


	double ref_resetTime, ref_startTime, ref_endTime, ref_FwdPeriod;
	double ref_FwdTarget, ref_FwdTarget_AdjDown, ref_FwdTarget_AdjUp;
	int ref_Size = refIdxResettimes.size();
	std::vector<double> ref_TargetVol_DownMinus(ref_Size,0.0);
	std::vector<double> ref_TargetVol_DownPlus(ref_Size,0.0);
	std::vector<double> ref_TargetVol_UpMinus(ref_Size,0.0);
	std::vector<double> ref_TargetVol_UpPlus(ref_Size,0.0);

	std::vector<double> ref_TargetVol_DownMinus_PayFloat(ref_Size,0.0);
	std::vector<double> ref_TargetVol_DownPlus_PayFloat(ref_Size,0.0);
	std::vector<double> ref_TargetVol_UpMinus_PayFloat(ref_Size,0.0);
	std::vector<double> ref_TargetVol_UpPlus_PayFloat(ref_Size,0.0);

	std::vector<double> ref_Shift_Vector(ref_Size);
	double ref_Shift, ref_Tfix;

// FIXMEFRED: mig.vc8 (30/05/2007 16:14:42): no vector<bool> anymore
	std::deque<bool>	CalibrateDown(ref_Size,true);
	std::deque<bool>	CalibrateUp(ref_Size,true);

	std::deque<bool>	CalibrateDown_PayFloat(ref_Size,true);
	std::deque<bool>	CalibrateUp_PayFloat(ref_Size,true);
	
	double  barrier_Down, barrier_DownMinus, barrier_DownPlus, barrier_Up, barrier_UpMinus, barrier_UpPlus;
	double ref_TargetPrice_DownMinus, ref_TargetPrice_DownPlus, ref_TargetPrice_UpMinus, ref_TargetPrice_UpPlus;

	double ref_TargetPrice_DownMinus_PayFloat, ref_TargetPrice_DownPlus_PayFloat, ref_TargetPrice_UpMinus_PayFloat, ref_TargetPrice_UpPlus_PayFloat;

	double atm_volDown, atm_volUp, atm_stdDevDown, atm_stdDevUp;
	double atm_stdDevDown_PayFloat, atm_stdDevUp_PayFloat;

	//calib barrier UP ?
	double barrier_UP = UpBarrier[0];

	double liborPayTarget, liborPayModel, startTime, endTime, shift;

	//Target Pay if Float
	if(indexPaymentType != IDXFIXED)
	{		
		liborPayTarget = corridor->GetPaymentFwdRates()->Elt(0)-couponMargin*100.0;
		liborPayTarget *= 0.01;

		startTime = corridor->GetPaymentStartDates()->Elt(0)-AsOfDate;
		endTime = corridor->GetPaymentEndDates()->Elt(0)-AsOfDate;

		shift = helper.LiborSLNShiftFromNumericalModel(0.0,
														resetTime,
														startTime, 
														endTime, 
														fwdPaymentPeriod);
	}


	//compute target vol
	for( int k = 0; k < ref_Size; ++k)
	{	
		ref_FwdTarget = corridor->GetRefIndexForwards()->Elt(k);
		ref_FwdTarget*=0.01;

		ref_resetTime = refIdxResettimes[k];
		ref_startTime = refIdxStarttimes[k];
		ref_endTime   = refIdxEndtimes[k];
		ref_FwdPeriod = refFwdPeriods[k];
	
		ref_Tfix	 = ref_resetTime / K_YEAR_LEN;

		ref_Shift = helper.LiborSLNShiftFromNumericalModel(0.0,
																ref_resetTime,
																ref_startTime, 
																ref_endTime, 
																ref_FwdPeriod);
		ref_Shift_Vector[k] = ref_Shift;
		
		//calib barrier DOWN ?
		barrier_Down= DownBarrier[k];
		atm_volDown = corridor->GetRefATMVols()->Elt(k); //VolATM_FromCorridor
		atm_stdDevDown = atm_volDown*ref_FwdTarget*sqrt(ref_Tfix);			
		if (fabs(barrier_Down - ref_FwdTarget)> NB_STDDEV_NO_CALIB * atm_stdDevDown)
			CalibrateDown[k] = false;
	
		bool boolCalibrateDown (CalibrateDown[k]);
		if(boolCalibrateDown)
		{
			ref_TargetPrice_DownMinus = corridor->GetPricesDownMinus()->Elt(k)*0.01;
			ref_TargetPrice_DownPlus = corridor->GetPricesDownPlus()->Elt(k)*0.01;			
			
			barrier_DownMinus	=barrier_Down-CALL_SPREAD_SHIFT;
			barrier_DownPlus	=barrier_Down+CALL_SPREAD_SHIFT;

			ref_TargetVol_DownMinus[k] = VanillaImpliedVol_BS(ref_FwdTarget+ref_Shift, 
															barrier_DownMinus+ref_Shift, 
															ref_Tfix, 
															ref_TargetPrice_DownMinus, 
															callPut,
															&atm_volDown);

			ref_TargetVol_DownPlus[k] = VanillaImpliedVol_BS(ref_FwdTarget+ref_Shift, 
															barrier_DownPlus+ref_Shift, 
															ref_Tfix, 
															ref_TargetPrice_DownPlus, 
															callPut,
															&atm_volDown);
		}

		//calib barrier UP ?
		barrier_Up= UpBarrier[k];
		atm_volUp = corridor->GetRefATMVols()->Elt(k);//VolATM_FromCorridor
		atm_stdDevUp = atm_volUp*ref_FwdTarget*sqrt(ref_Tfix);			
		if (fabs(barrier_Up - ref_FwdTarget)> NB_STDDEV_NO_CALIB * atm_stdDevUp)
			CalibrateUp[k] = false;

		bool boolCalibrateUp (CalibrateUp[k]);
		if(boolCalibrateUp)
		{
			ref_TargetPrice_UpMinus = corridor->GetPricesUpMinus()->Elt(k)*0.01;
			ref_TargetPrice_UpPlus = corridor->GetPricesUpPlus()->Elt(k)*0.01;

			barrier_UpMinus		=barrier_Up-CALL_SPREAD_SHIFT;
			barrier_UpPlus		=barrier_Up+CALL_SPREAD_SHIFT;

			ref_TargetVol_UpMinus[k] = VanillaImpliedVol_BS(ref_FwdTarget+ref_Shift, 
															barrier_UpMinus+ref_Shift, 
															ref_Tfix, 
															ref_TargetPrice_UpMinus, 
															callPut,
															&atm_volUp);
			
			ref_TargetVol_UpPlus[k] = VanillaImpliedVol_BS(ref_FwdTarget+ref_Shift, 
															barrier_UpPlus+ref_Shift, 
															ref_Tfix, 
															ref_TargetPrice_UpPlus, 
															callPut,
															&atm_volUp);
		}

		///Float Pay Case
		if(indexPaymentType != IDXFIXED)
		{
			//Barrier down
			ref_FwdTarget_AdjDown = corridor->GetRefIndexForwardsAdjDown()->Elt(k);
			ref_FwdTarget_AdjDown *= 0.01;
			atm_stdDevDown_PayFloat = atm_volDown*ref_FwdTarget_AdjDown*sqrt(ref_Tfix);			
			if (fabs(barrier_Down - ref_FwdTarget_AdjDown)> NB_STDDEV_NO_CALIB * atm_stdDevDown_PayFloat)
				CalibrateDown_PayFloat[k] = false;
					
			if(CalibrateDown_PayFloat[k])
			{
				ref_TargetPrice_DownMinus_PayFloat = corridor->GetPricesDownMinus_PayFloat()->Elt(k)*0.01;
				ref_TargetPrice_DownPlus_PayFloat = corridor->GetPricesDownPlus_PayFloat()->Elt(k)*0.01;			
				
				barrier_DownMinus	=barrier_Down-CALL_SPREAD_SHIFT;
				barrier_DownPlus	=barrier_Down+CALL_SPREAD_SHIFT;

				ref_TargetVol_DownMinus_PayFloat[k] = VanillaImpliedVol_BS(ref_FwdTarget_AdjDown+ref_Shift, 
																barrier_DownMinus+ref_Shift, 
																ref_Tfix, 
																ref_TargetPrice_DownMinus_PayFloat, 
																callPut,
																&atm_volDown);

				ref_TargetVol_DownPlus_PayFloat[k] = VanillaImpliedVol_BS(ref_FwdTarget_AdjDown+ref_Shift, 
																barrier_DownPlus+ref_Shift, 
																ref_Tfix, 
																ref_TargetPrice_DownPlus_PayFloat, 
																callPut,
																&atm_volDown);
			}

			//Barrier Up
			ref_FwdTarget_AdjUp = corridor->GetRefIndexForwardsAdjUp()->Elt(k);
			ref_FwdTarget_AdjUp *= 0.01;
			atm_stdDevUp_PayFloat = atm_volUp*ref_FwdTarget_AdjUp*sqrt(ref_Tfix);			
			if (fabs(barrier_Up - ref_FwdTarget_AdjUp)> NB_STDDEV_NO_CALIB * atm_stdDevUp_PayFloat)
				CalibrateUp_PayFloat[k] = false;

			if(CalibrateUp_PayFloat[k])
			{
				ref_TargetPrice_UpMinus_PayFloat = corridor->GetPricesUpMinus_PayFloat()->Elt(k)*0.01;
				ref_TargetPrice_UpPlus_PayFloat = corridor->GetPricesUpPlus_PayFloat()->Elt(k)*0.01;

				barrier_UpMinus		=barrier_Up-CALL_SPREAD_SHIFT;
				barrier_UpPlus		=barrier_Up+CALL_SPREAD_SHIFT;

				ref_TargetVol_UpMinus_PayFloat[k] = VanillaImpliedVol_BS(ref_FwdTarget_AdjUp+ref_Shift, 
																barrier_UpMinus+ref_Shift, 
																ref_Tfix, 
																ref_TargetPrice_UpMinus_PayFloat, 
																callPut,
																&atm_volUp);
				
				ref_TargetVol_UpPlus_PayFloat[k] = VanillaImpliedVol_BS(ref_FwdTarget_AdjUp+ref_Shift, 
																barrier_UpPlus+ref_Shift, 
																ref_Tfix, 
																ref_TargetPrice_UpPlus_PayFloat, 
																callPut,
																&atm_volUp);
			}			
		}
	}
	for (i=0; i<sizeEvalTimes; i++)
	{
		// Adj_payIdx if pay Float
		if(indexPaymentType != IDXFIXED)
		{
			liborPayModel = helper.LiborRateFromNumericalModel(evalTimes[i], 
																			payTime, 
																			resetTime,
																			startTime, 
																			endTime, 
																			fwdPaymentPeriod, 
																			couponMargin);

			double adjustment  = liborPayTarget/liborPayModel;				
			GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(PAY_ADJ_SURFACE, evalTimes[i], resetTime, adjustment);
			GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(PAY_IDX_SHIFT_SURFACE, evalTimes[i], resetTime, shift);
		}

		(*itsStatisticCalibResult)(i,0) = evalTimes[i];
		int nbOfReset = 0;
		int nbOfResetMinusDown = 0;
		int nbOfResetPlusDown = 0;
		int nbOfResetMinusUp = 0;
		int nbOfResetPlusUp = 0;
		int nbOfResetPayFloatMinusDown = 0;
		int nbOfResetPayFloatPlusDown = 0;
		int nbOfResetPayFloatMinusUp = 0;
		int nbOfResetPayFloatPlusUp = 0;

		for( int k = 0; k < ref_Size; ++k)
		{
			// get computed libor rate (potentially adjusted)
			double ref_FwdTarget = corridor->GetRefIndexForwards()->Elt(k);
			ref_FwdTarget *= 0.01;

			ref_resetTime = refIdxResettimes[k];
			ref_startTime = refIdxStarttimes[k];
			ref_endTime   = refIdxEndtimes[k];
			ref_FwdPeriod = refFwdPeriods[k];
			ref_Shift = ref_Shift_Vector[k];

			ref_Tfix	 = ref_resetTime / K_YEAR_LEN;

			// calibration only if option maturity is after evalTime
			if (ref_resetTime >= evalTimes[i] )
			{	
				nbOfReset++;
				double ref_FwdModel = helper.LiborRateFromNumericalModel(evalTimes[i], 
																				payTime, 
																		ref_resetTime,
																			ref_startTime, 
																			ref_endTime, 
																			ref_FwdPeriod, 
																			couponMargin);

				
				double adjustment  = ref_FwdTarget/ref_FwdModel;				
				GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(FIRST_SURFACE, evalTimes[i], ref_resetTime, adjustment);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(FIRST_SURFACE, evalTimes[i], ref_resetTime, ref_Shift);
				
				// Adj_FwdRef_Down Adj_FwdRef_Up if payIdx Float
				if(indexPaymentType != IDXFIXED)
				{
					double ref_FwdTarget_AdjDown = corridor->GetRefIndexForwardsAdjDown()->Elt(k);
					ref_FwdTarget_AdjDown *= 0.01;

					double ref_FwdTarget_AdjUp = corridor->GetRefIndexForwardsAdjUp()->Elt(k);
					ref_FwdTarget_AdjUp *= 0.01;

					double rTime = resetTime > ref_resetTime ? ref_resetTime : resetTime;
					double toTime = rTime > ref_startTime ? ref_resetTime : rTime;
				   
					double evalTime= evalTimes[i]; 
					double cov = helper.CovarFromNumericalModel(0.0,
														evalTime,
														rTime,
														ref_startTime);			
				    
					cov *=((liborPayModel+shift)/liborPayModel);
					double adj = exp(cov);
					double ref_FwdModel_Adj = adj*(ref_FwdModel+ref_Shift)-ref_Shift;
					double adjustmentdown  = (ref_FwdTarget_AdjDown)/ref_FwdModel_Adj;				
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DOWN_ADJ_SURFACE, evalTimes[i], ref_resetTime, adjustmentdown);
					double adjustmentUp  = (ref_FwdTarget_AdjUp)/ref_FwdModel_Adj;				
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(UP_ADJ_SURFACE, evalTimes[i], ref_resetTime, adjustmentUp);
				}

				double Tcall = evalTimes[i] / K_YEAR_LEN;				
				if (Tcall==ref_Tfix)
					ref_Tfix += 1./K_YEAR_LEN;

				double ref_forwardVol_DownMinus;
				double ref_forwardVol_DownPlus;
				double ref_forwardVol_UpMinus;
				double ref_forwardVol_UpPlus;

				double ref_forwardVol_DownMinus_PayFloat;
				double ref_forwardVol_DownPlus_PayFloat;
				double ref_forwardVol_UpMinus_PayFloat;
				double ref_forwardVol_UpPlus_PayFloat;


				double var,ref_ModelVol;
				
				//DOWN
				bool boolCalibrateDown (CalibrateDown[k]);
				if (boolCalibrateDown) 
				{			
					//Down_Minus
					barrier_Down = DownBarrier[k];
					barrier_DownMinus	=barrier_Down-CALL_SPREAD_SHIFT;
					

					ref_ModelVol = helper.LiborSLNVolatilityFromNumericalModel(evalTimes[i],
																				ref_resetTime,
																				ref_startTime, 
																				ref_endTime, 
																				ref_FwdPeriod, 
																				barrier_DownMinus,
																				ref_Shift);			
					
					var = (	ref_TargetVol_DownMinus[k]  * ref_TargetVol_DownMinus[k] * ref_Tfix 
								  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

					/// everything is fine
					if (var>=0)
						ref_forwardVol_DownMinus = sqrt( var / (ref_Tfix - Tcall) );
					
					/// oups: variance squeeze
					else 
					{
						ref_forwardVol_DownMinus = 0.0;
						SetVarianceSqueezeStatus(true);
					
						nbOfResetMinusDown++;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Down minus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}

					//Down_Plus
					barrier_DownPlus	=barrier_Down+CALL_SPREAD_SHIFT;
					ref_ModelVol = helper.LiborSLNVolatilityFromNumericalModel(evalTimes[i],
																				ref_resetTime,
																				ref_startTime, 
																				ref_endTime, 
																				ref_FwdPeriod, 
																				barrier_DownPlus,
																				ref_Shift);			
					
					var = (	ref_TargetVol_DownPlus[k] * ref_TargetVol_DownPlus[k] * ref_Tfix 
								  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

					/// everything is fine
					if (var>=0)
						ref_forwardVol_DownPlus = sqrt( var / (ref_Tfix - Tcall) );
					
					/// oups: variance squeeze
					else 
					{
						ref_forwardVol_DownPlus = 0.0;
						SetVarianceSqueezeStatus(true);
					
						nbOfResetPlusDown++;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Down plus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}
				}
				else
				{
					ref_forwardVol_DownMinus = 0.0;
					ref_forwardVol_DownPlus  = 0.0;
					CC_Ostringstream os;
					os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Down seems to be to far from ATM. Null forward vol inputed..."<< int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
					ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				}

				//UP
				bool boolCalibrateUp (CalibrateUp[k]);
				if(boolCalibrateUp)
				{
					//Up_Minus
					barrier_Up = UpBarrier[k];
					barrier_UpMinus	=barrier_Up-CALL_SPREAD_SHIFT;
					

					ref_ModelVol = helper.LiborSLNVolatilityFromNumericalModel(evalTimes[i],
																				ref_resetTime,
																				ref_startTime, 
																				ref_endTime, 
																				ref_FwdPeriod, 
																				barrier_UpMinus,
																				ref_Shift);			
					
					var = (	ref_TargetVol_UpMinus[k]  * ref_TargetVol_UpMinus[k] * ref_Tfix 
								  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

					/// everything is fine
					if (var>=0)
						ref_forwardVol_UpMinus = sqrt( var / (ref_Tfix - Tcall) );
					
					/// oups: variance squeeze
					else 
					{
						ref_forwardVol_UpMinus = 0.0;
						SetVarianceSqueezeStatus(true);

						nbOfResetMinusUp++;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Up Minus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}

					//Up_Plus
					barrier_UpPlus	=barrier_Up+CALL_SPREAD_SHIFT;
					ref_ModelVol = helper.LiborSLNVolatilityFromNumericalModel(evalTimes[i],
																				ref_resetTime,
																				ref_startTime, 
																				ref_endTime, 
																				ref_FwdPeriod, 
																				barrier_UpPlus,
																				ref_Shift);			
					
					var = (	ref_TargetVol_UpPlus[k] * ref_TargetVol_UpPlus[k] * ref_Tfix 
								  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

					/// everything is fine
					if (var>=0)
						ref_forwardVol_UpPlus = sqrt( var / (ref_Tfix - Tcall) );
					
					/// oups: variance squeeze
					else 
					{
						ref_forwardVol_UpPlus = 0.0;
						SetVarianceSqueezeStatus(true);
					
						nbOfResetPlusUp++;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Up plus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}
				}
				else
				{
					ref_forwardVol_UpMinus = 0.0;
					ref_forwardVol_UpPlus  = 0.0;
					CC_Ostringstream os;
					os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Barrier Up seems to be to far from ATM. Null forward vol inputed..."<< int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
					ARM_TheEventViewer.Instance()->AddToMessage(os.str());
				}

				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_MINUS_SURFACE, evalTimes[i], ref_resetTime, ref_forwardVol_DownMinus);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_PLUS_SURFACE, evalTimes[i], ref_resetTime, ref_forwardVol_DownPlus);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_MINUS_SURFACE, evalTimes[i], ref_resetTime, ref_forwardVol_UpMinus);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_PLUS_SURFACE, evalTimes[i], ref_resetTime, ref_forwardVol_UpPlus);

				//if pay Float 
				if(indexPaymentType != IDXFIXED)
				{
					////////////////////////////////////////////////////
					//DOWN
					bool boolCalibrateDown_PayFloat(CalibrateDown_PayFloat[k]);
					if (boolCalibrateDown_PayFloat) 
					{				
						var = (	ref_TargetVol_DownMinus_PayFloat[k]  * ref_TargetVol_DownMinus_PayFloat[k] * ref_Tfix 
									  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

						/// everything is fine
						if (var>=0)
							ref_forwardVol_DownMinus_PayFloat = sqrt( var / (ref_Tfix - Tcall) );
						
						/// oups: variance squeeze
						else 
						{
							ref_forwardVol_DownMinus_PayFloat = 0.0;
							SetVarianceSqueezeStatus(true);
					
							nbOfResetPayFloatMinusDown ++;		
							CC_Ostringstream os;
							os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Down Minus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
							ARM_TheEventViewer.Instance()->AddToMessage(os.str());
						}

						//Down_Plus
						var = (	ref_TargetVol_DownPlus_PayFloat[k] * ref_TargetVol_DownPlus_PayFloat[k] * ref_Tfix 
									  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

						/// everything is fine
						if (var>=0)
							ref_forwardVol_DownPlus_PayFloat = sqrt( var / (ref_Tfix - Tcall) );
						
						/// oups: variance squeeze
						else 
						{
							ref_forwardVol_DownPlus_PayFloat = 0.0;
							SetVarianceSqueezeStatus(true);
					
							nbOfResetPayFloatPlusDown++;		
							CC_Ostringstream os;
							os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Down Plus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
							ARM_TheEventViewer.Instance()->AddToMessage(os.str());
						}
					}
					else
					{
						ref_forwardVol_DownMinus_PayFloat = 0.0;
						ref_forwardVol_DownPlus_PayFloat  = 0.0;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Barrier Down seems to be to far from ATM. Null forward vol inputed..."<< int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}

					//UP
					bool boolCalibrateUp_PayFloat (CalibrateUp_PayFloat[k]);
					if(boolCalibrateUp_PayFloat)
					{					
						var = (	ref_TargetVol_UpMinus_PayFloat[k]  * ref_TargetVol_UpMinus_PayFloat[k] * ref_Tfix 
									  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

						/// everything is fine
						if (var>=0)
							ref_forwardVol_UpMinus_PayFloat = sqrt( var / (ref_Tfix - Tcall) );
						
						/// oups: variance squeeze
						else 
						{
							ref_forwardVol_UpMinus_PayFloat = 0.0;
							nbOfResetPayFloatMinusUp++;
							SetVarianceSqueezeStatus(true);
					
		
							CC_Ostringstream os;
							os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Up Minus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
							ARM_TheEventViewer.Instance()->AddToMessage(os.str());
						}

						//Up_Plus
						var = (	ref_TargetVol_UpPlus_PayFloat[k] * ref_TargetVol_UpPlus_PayFloat[k] * ref_Tfix 
									  -	ref_ModelVol   * ref_ModelVol  * Tcall) ;

						/// everything is fine
						if (var>=0)
							ref_forwardVol_UpPlus_PayFloat = sqrt( var / (ref_Tfix - Tcall) );
						
						/// oups: variance squeeze
						else 
						{
							ref_forwardVol_UpPlus_PayFloat = 0.0;
							SetVarianceSqueezeStatus(true);
					
							nbOfResetPayFloatPlusUp++;
							CC_Ostringstream os;
							os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Up Plus variance squeeze [call date = " << int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
							ARM_TheEventViewer.Instance()->AddToMessage(os.str());
						}
					}
					else
					{
						ref_forwardVol_UpMinus_PayFloat = 0.0;
						ref_forwardVol_UpPlus_PayFloat  = 0.0;
						CC_Ostringstream os;
						os << ARM_USERNAME << " : Local_SLN_Model / corridorlet Calibration : Pay Float Barrier Up seems to be to far from ATM. Null forward vol inputed..."<< int(Tcall*10.0)/10.0 << " Y, reset date = " << int(ref_Tfix*10.0)/10.0 << " Y]"<< endl;
						ARM_TheEventViewer.Instance()->AddToMessage(os.str());
					}

					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_MINUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, ref_forwardVol_DownMinus_PayFloat);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_PLUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, ref_forwardVol_DownPlus_PayFloat);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_MINUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, ref_forwardVol_UpMinus_PayFloat);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_PLUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, ref_forwardVol_UpPlus_PayFloat);
				}//PayFloat
			}
			else
			{
				GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(FIRST_SURFACE, evalTimes[i], ref_resetTime, 1.0);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Shift).SetValue(FIRST_SURFACE, evalTimes[i], ref_resetTime, 0.0);

				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_MINUS_SURFACE, evalTimes[i], ref_resetTime, 0.0);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_PLUS_SURFACE, evalTimes[i], ref_resetTime, 0.0);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_MINUS_SURFACE, evalTimes[i], ref_resetTime, 0.0);
				GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_PLUS_SURFACE, evalTimes[i], ref_resetTime, 0.0);

				if(indexPaymentType != IDXFIXED)
				{
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(DOWN_ADJ_SURFACE, evalTimes[i], ref_resetTime, 1.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::ForwardAdjustment).SetValue(UP_ADJ_SURFACE, evalTimes[i], ref_resetTime, 1.0);

					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_MINUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(DOWN_PLUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_MINUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, 0.0);
					GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility).SetValue(UP_PLUS_SURFACE_PAYFLOAT, evalTimes[i], ref_resetTime, 0.0);

				}
			}
		}
		(*itsStatisticCalibResult)(i,1) +=nbOfReset;
		(*itsStatisticCalibResult)(i,2) +=nbOfResetMinusDown;
		(*itsStatisticCalibResult)(i,3) +=nbOfResetPlusDown;
		(*itsStatisticCalibResult)(i,4) +=nbOfResetMinusUp;
		(*itsStatisticCalibResult)(i,5) +=nbOfResetPlusUp;
		(*itsStatisticCalibResult)(i,6) +=nbOfResetPayFloatMinusDown;
		(*itsStatisticCalibResult)(i,7) +=nbOfResetPayFloatPlusDown;
		(*itsStatisticCalibResult)(i,8) +=nbOfResetPayFloatMinusUp;
		(*itsStatisticCalibResult)(i,9) +=nbOfResetPayFloatPlusUp;			
	}
}

///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ
////////////////////////////////////////////////////////
/// Local_Normal_Model_Helper implementation ///////////
////////////////////////////////////////////////////////
///µµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµµ


// set static attributes
int Local_SLN_Model_Calibration_Helper::SFRM  = 1;

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : constructor
///	Returns : void
///	Action  : contextual constructor
////////////////////////////////////////////////////
Local_SLN_Model_Calibration_Helper::Local_SLN_Model_Calibration_Helper (ARM_PricingModel* numericalModel)
:	itsNumericalModel (numericalModel),
	itsModelType (0)
{
	if ( dynamic_cast<ARM_SFRM*> (numericalModel) )
		itsModelType = SFRM;
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local SLN Model Calibration : numerical model is not supported" );
}

////////////////////////////////////////////////////
///	Class   : Local_Normal_Model_Calibration_Helper
///	Routine : constructor
///	Returns : void
///	Action  : copy constructor
////////////////////////////////////////////////////
Local_SLN_Model_Calibration_Helper::Local_SLN_Model_Calibration_Helper (const Local_SLN_Model_Calibration_Helper& rhs)
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
double Local_SLN_Model_Calibration_Helper::LiborRateFromNumericalModel(	double evalTime,
																		double payTime, 
																		double fwdResetTime,
																		double fwdStartTime,
																		double fwdEndTime,
																		double fwdPeriod,
																		double strike ) const
{		
	if (itsModelType == SFRM )
	{
		ARM_ModelParamsSFRM* params = (ARM_ModelParamsSFRM*)itsNumericalModel->GetModelParams();
		
		double dfStart = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdStartTime / K_YEAR_LEN);
		double dfEnd   = itsNumericalModel->GetZeroCurve()->DiscountPrice(fwdEndTime / K_YEAR_LEN);


		if (fabs(fwdEndTime - payTime) < 5.0 )
		{
			double libor = (dfStart/dfEnd-1.0)/fwdPeriod;
			return libor;
		}
		else
		{
			double libor =0.0;
			double var  = params->IntegratedVariance(0.0,evalTime,fwdResetTime);
			double M    = params->ShiftValue(fwdResetTime);
			double F,FM;
			if(payTime <= fwdStartTime + 5)
			{
				/// Use the in-arrears correction (exact if pay = start !)
				double payConvex = exp(var) - 1.0;
				F = (dfStart/dfEnd-1.0)/fwdPeriod;
				FM = F + M;
				libor = F + fwdPeriod*FM*FM*payConvex/(1+fwdPeriod*F);
				
			}
			else
			{
				/// Case start < pay <= end or end < pay
				/// Use the payment lag correction with the additive approximation
				/// It converges to the exact in-arrears correction if pay=start
				/// and obviously to no conversion if pay=end

				double ZcPay = itsNumericalModel->GetZeroCurve()->DiscountPrice(payTime / K_YEAR_LEN);
				
				bool isPayBeforeEnd = (payTime < fwdEndTime);

				// Compute the ccovariance between L(t,Start,End) & L(t,Start,Pay)
				double coVarSP=0.0;

				ARM_IRIndex* pIndex	= params->GetIRIndex();
				double indexPeriod = 1.0/pIndex->GetResetFrequency();

				int nbFwds = (payTime-fwdStartTime)/K_YEAR_LEN/indexPeriod;

				double startFwd, endFwd, fwdDelta, shift, cov, fwdVal;

				double ZcFwdStart;
				double ZcFwdEnd;

				for (int i = 0; i <= nbFwds; ++i)
				{
					startFwd = fwdStartTime+i*indexPeriod*K_YEAR_LEN;
					if (i < nbFwds)
						endFwd = startFwd+indexPeriod*K_YEAR_LEN;
					else
						endFwd = payTime;

					if (fabs(endFwd-startFwd) > K_DOUBLE_TOL)
					{
						shift = params->ShiftValue(startFwd);

						cov = params->IntegratedCovariance(0.0,evalTime,fwdStartTime,startFwd);;

						ZcFwdStart	=	itsNumericalModel->GetZeroCurve()->DiscountPrice(startFwd / K_YEAR_LEN);
						ZcFwdEnd	=	itsNumericalModel->GetZeroCurve()->DiscountPrice(endFwd / K_YEAR_LEN);

						fwdDelta = (endFwd-startFwd)/360;
						fwdVal= (ZcFwdStart/ZcFwdEnd-1.0)/fwdDelta;
						coVarSP += fwdDelta*(fwdVal+shift)/(1+fwdDelta*fwdVal)*cov;
					}
				}

				double MSP = M; // assuming L(t,Start,Pay) same shift as L(t,Start,End)
				double deltaSP = fwdPeriod*(payTime - fwdStartTime)/(fwdEndTime-fwdStartTime);

				double MEP  = M; // assuming L(t,Pay,End) or L(t,End,Pay) same shift as L(t,Start,End)
				double deltaEP,varEP=0.0;
				if(isPayBeforeEnd)
					deltaEP = fwdPeriod*(fwdEndTime-payTime)/(fwdEndTime-fwdStartTime);
				else
				{
					deltaEP = fwdPeriod*(payTime-fwdEndTime)/(fwdEndTime-fwdStartTime);
					varEP = params->IntegratedCovariance(0.0,evalTime,fwdEndTime,payTime,payTime);
				}

				double FSP,FMSP; // L(t,Start,Pay) & shifted libor 
				double FEP,FMEP; // L(t,Pay,End) or L(t,End,Pay) & shifted libor
				double coVar,volCoefEP,adjustEP=0.0;
				
				F = (dfStart/dfEnd-1.0)/fwdPeriod;
				FM = F + M;

				FSP = (dfStart/ZcPay-1.0)/deltaSP;
				FMSP = FSP + MSP;

				if(isPayBeforeEnd)
				{
					FEP = (ZcPay/dfEnd-1.0) / deltaEP;
					volCoefEP = deltaEP*(FEP+MEP)/(1 + deltaEP*FEP);
				}
				else
				{
					FEP = (dfEnd/ZcPay-1.0) / deltaEP;
					FMEP = FEP + MEP;
					volCoefEP = - deltaEP * FMEP / (1 + deltaEP*FEP);
					adjustEP = - volCoefEP * FMEP * (exp(varEP)-1.0);
					FEP += adjustEP;
					adjustEP *= F/FMEP;
				}
				FMEP = FEP + MEP;

				coVar = ( fwdPeriod*FM/(1+fwdPeriod*F)*var - coVarSP ) / volCoefEP;

				libor = F + volCoefEP * ( FM * (exp(coVar) - 1.0) + adjustEP );				
			}
			return libor;
		}		
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
double Local_SLN_Model_Calibration_Helper::LiborSLNVolatilityFromNumericalModel(
																		double expiryTime,
																		double fwdResetTime, 
																		double fwdStartTime,
																		double fwdEndTime,
																		double fwdPeriod,
																		double strike,
																		double shift) const
{		
	//
	// caution : gaussian approx on libor (which is lognorm. shifted under HW)
	//
	if (itsModelType == SFRM )
	{	
		ARM_ModelParamsSFRM* params = (ARM_ModelParamsSFRM*)itsNumericalModel->GetModelParams();

		//const double shift	= params->ShiftValue(fwdResetTime);
		const double stdDev	= params->IntegratedVol(0.0,expiryTime,fwdResetTime)/sqrt(expiryTime/K_YEAR_LEN);
			
		return stdDev;
	}
	///
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::LiborNormalVolatilityFromNumericalModel : numerical model is not supported" );
}


double Local_SLN_Model_Calibration_Helper::LiborSLNShiftFromNumericalModel(
									double expiryTime,
									double fwdResetTime , 
									double fwdStartTime,
									double fwdEndTime,
									double fwdPeriod) const
{
	if (itsModelType == SFRM )
	{	
		ARM_ModelParamsSFRM* params = (ARM_ModelParamsSFRM*)itsNumericalModel->GetModelParams();

		const double shift	= params->ShiftValue(fwdResetTime);
		return shift;
	}
	///
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::LiborNormalVolatilityFromNumericalModel : numerical model is not supported" );

}

double Local_SLN_Model_Calibration_Helper::CovarFromNumericalModel(
									double evalTime,
									double resetTime,
									double fwdResetTime1 , 
									double fwdResetTime2) const
{
	if (itsModelType == SFRM )
	{	
		ARM_ModelParamsSFRM* params = (ARM_ModelParamsSFRM*)itsNumericalModel->GetModelParams();	

		double Vol1		= params->IntegratedVol(evalTime,resetTime,fwdResetTime1);
		double Vol2		= params->IntegratedVol(evalTime,resetTime,fwdResetTime2);
            
					/// WARNING ASSUMPTION: correlation is set to 100% 
		double cov=(Vol1*Vol2);
		return cov;
	}
	///
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local_Normal_Model_Calibration_Helper::LiborNormalVolatilityFromNumericalModel : numerical model is not supported" );

}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

