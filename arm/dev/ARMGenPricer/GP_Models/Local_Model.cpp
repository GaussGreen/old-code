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

#include "gpmodels/local_model.h"
#include "gpmodels/multiassets.h"

/// gpinfra
#include "gpinfra/pricingstates.h"
#include "gpinfra/modelnamemap.h"


/// gpbase
#include "gpbase/surface.h"
#include "gpbase/datestrip.h"

/// closedforms
#include "gpclosedforms/normal.h"


// kernel
#include <inst/portfolio.h>
#include <inst/spreadoption.h>
#include <inst/capfloor.h>
#include <inst/option.h>
#include <inst/corridorleg.h>
#include <inst/corridordblcondition.h>


/// gpcalib
#include "gpcalib/kerneltogp.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/vanillacap.h"
#include "gpcalib/vanillaeqfxoption.h"
#include "gpcalib/vanillacorridor.h"
#include "gpcalib/densityfunctors.h"



CC_BEGIN_NAMESPACE( ARM )


////////////////////////////////////////////////////
///	Class  : ARM_Local_Model
///	Routine: Constructor
///	Returns: 
///	Action : builds the ARM_Local_Model
////////////////////////////////////////////////////
ARM_Local_Model::ARM_Local_Model(const ARM_ZeroCurvePtr& zc, const ARM_ModelParams& params)
:	ARM_PricingModel( zc, &params ),
	itsNumericalModel (NULL),
	itsVarianceSqueezInfo(),
	itsVarianceSqueezeStatus(false),
	itsResetCalib(true),
	itsCorrelUnSqueezer(false),
	itsVolUnSqueezer(false),
	itsFunctionals(0),
	itsGrid(NULL),
	itsFwds(0),
	itsVols(0),
	itsShifts(0),
	itsResetTimes(0),
	itsFunctionalFlag(false)
{}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routines: Copy constructor
///	Returns :
///	Action  : 
////////////////////////////////////////////////////
ARM_Local_Model::ARM_Local_Model(const ARM_Local_Model& rhs)
:	ARM_PricingModel (rhs),
	itsNumericalModel (rhs.itsNumericalModel),
	itsVarianceSqueezInfo(rhs.itsVarianceSqueezInfo),
	itsVarianceSqueezeStatus(rhs.itsVarianceSqueezeStatus),
	itsResetCalib(rhs.itsResetCalib),
	itsCorrelUnSqueezer(rhs.itsCorrelUnSqueezer),
	itsVolUnSqueezer(rhs.itsVolUnSqueezer),
	itsGrid(rhs.itsGrid),
	itsFwds(rhs.itsFwds),
	itsVols(rhs.itsVols),
	itsShifts(rhs.itsShifts),
	itsResetTimes(rhs.itsResetTimes),
	itsFunctionalFlag(rhs.itsFunctionalFlag)
{
	DuplicateCloneablePtrVectorInPlace<ARM_GP_Vector> (rhs.itsFunctionals, itsFunctionals);
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Model
///	Routine: GetNumericalModel
///	Returns: ARM_PricingModel*
///	Action : returns numerical model (if already setted)
////////////////////////////////////////////////////
ARM_PricingModel* ARM_Local_Model::GetNumericalModel() const
{
	if (itsNumericalModel) 
		return itsNumericalModel;
	else 
		ARM_THROW( ERR_INVALID_ARGUMENT, "Numerical Model has not been setted" );
}

////////////////////////////////////////////////////
///	Class  : ARM_Local_Model
///	Routine: SetNumericalModel
///	Returns: void
///	Action : set the numerical model
////////////////////////////////////////////////////
void ARM_Local_Model::SetNumericalModel (ARM_PricingModel* numericalModel)
{
	itsNumericalModel = numericalModel;
	SetZeroCurve(itsNumericalModel->GetZeroCurve());
}



////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routines: Update the links on the numerical model
///	Returns : void
///	Action  : set or update numerical model
////////////////////////////////////////////////////
void ARM_Local_Model::UpdateLinks( const ARM_MultiAssetsModel& multiAssetsModel )
{	
	const ARM_ModelNameMap* modelMap = multiAssetsModel.GetModelMap();

	ARM_PricingModel* newNumericalModel=NULL;

	/// Restore linked models
	ARM_IntVector itsOtherModels = (*modelMap)[GetModelName()]->OtherModelRefNb();
	
	if(itsOtherModels.size() > 0)
		/// By default the numerical model is the first in the list
		newNumericalModel = &*((*modelMap)[itsOtherModels[0]]->Model());
	else if(itsNumericalModel)
		newNumericalModel = itsNumericalModel;
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, " : no numerical model for LocalModel" );
	}

	SetNumericalModel(newNumericalModel); /// to update zc curve if changed
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : CalibrateLocalModel
///	Returns : void
///	Action  : Local model calibration (caplet / spread option)
///			  uses CalibrateLocalModel( ARM_Security, evalTimes) defined
///			  in specialized classes
////////////////////////////////////////////////////
void ARM_Local_Model::CalibrateLocalModel(const ARM_Portfolio& portfolio, const ARM_GP_Vector& evalTimes)
{	
	// reset model params
	ResetModelParams();
	
	// we should check that instruments expiry dates are all different
	int size = portfolio.GetSize();
	
	for (size_t i(0); i<size; i++)
	{
		ARM_Security* security = portfolio.GetAsset(i);
		if(security)
		{
			double targetPrice = portfolio.GetMktPrices()->Elt(i);
			CalibrateLocalModel(*security, targetPrice, evalTimes,i);
		}
	}
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : CalibrateLocalModelFunctional
///	Returns : void
///	Action  : Local model calibration on functionals
////////////////////////////////////////////////////
void ARM_Local_Model::CalibrateLocalModelFunctional(
			vector<ARM_Security*> securities,
			vector<ARM_DensityFunctor*> densities,
			int sizeGrid,
			double nbStdDev,
			bool rescaling)
{
	SetFunctionalFlag(true);
	int size = securities.size();

	ARM_NoVolDensityFunctor defaultDensity;

	if (densities.size()==size)
	{
		itsGrid = ARM_GP_VectorPtr(new ARM_GP_Vector(sizeGrid));
		for (int k=0;k<sizeGrid;k++)
		{
			(*itsGrid)[k]=(k-sizeGrid/2)*2.*nbStdDev/(sizeGrid-1);
		}
		
		for (int i=0;i<size;i++)
		{
			CalibrateLocalModel(*securities[i],densities[i] ? *densities[i] : defaultDensity,rescaling);
		}
	}
	else
	{
		for (int i=0;i<size;i++)
		{
			CalibrateLocalModel(*securities[i],defaultDensity,rescaling);
		}
	}
};

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : GetCapFloorData
///	Returns : void
///	Action  : collects all needed info from ARM_CapFloor
////////////////////////////////////////////////////
void ARM_Local_Model::GetCapFloorLetData (	ARM_CapFloor* capFloor,
											size_t periodIdx,
											// outputs
											double& resetTime,
											double& payTime,
											double& payPeriod,
											double& payNotional,
											double& liborStartTime,
											double& liborEndTime,
											double& liborPeriod,
											double& strike,
											int&	callPut) const
{
	// Validations ...
	size_t capSize = capFloor->GetFlowEndDates()->GetSize();

	if (periodIdx<0 || periodIdx>=capSize)
		ARM_THROW( ERR_INVALID_ARGUMENT, "GetCapFloorLetData : invalid periodIdx" );

	if (capFloor->GetSwapLeg()->GetSpread() != 0.0)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Model Calibration : ARM_CapFloor : spread not supported" );

	callPut	= capFloor->GetCapFloorType();

	// to be checked with JP Riaudel
	resetTime		= (*capFloor->GetResetDates())[periodIdx]; 
	payTime			= (*capFloor->GetPaymentDates())[periodIdx];
	payPeriod		= (*capFloor->GetSwapLeg()->GetInterestTerms())[periodIdx];
	payNotional		= capFloor->GetAmount()->CptReferenceValue(payTime);
	liborStartTime	= (*capFloor->GetSwapLeg()->GetFwdRateStartDates())[periodIdx];
	liborEndTime	= (*capFloor->GetSwapLeg()->GetFwdRateEndDates())[periodIdx];
			
	strike			= ( ( capFloor->GetStrikes()
						? capFloor->GetStrikes()->CptReferenceValue( resetTime)
						: capFloor->GetStrike() ) - capFloor->GetSwapLeg()->ComputeSpread(periodIdx) 
					  ) / CC_NS(ARM_Constants,rateBase);
			
	liborPeriod = CountYears ( capFloor->GetSwapLeg()->GetIRIndex()->GetDayCount(),
							   liborStartTime, 
                               liborEndTime); 

	
	resetTime		= GetTimeFromDate(resetTime);
	payTime			= GetTimeFromDate(payTime);
	liborStartTime	= GetTimeFromDate(liborStartTime);
	liborEndTime	= GetTimeFromDate(liborEndTime);
	
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : GetSpreadOptionLetData
///	Returns : void
///	Action  : collects all needed info from ARM_SpreadOption
////////////////////////////////////////////////////
void ARM_Local_Model::GetSpreadOptionLetData (	ARM_SpreadOption*	spreadOption,
												size_t				periodIdx,

												// outputs
												int&			callPut,
												double&			resetTime,
												double&			payTime,
												double&			payPeriod,
												double&			payNotional,
												double&			coeffLong,
												double&			coeffShort,
												double&			strike,
												double&			swapLongFloatStartTime,
												double&			swapLongFloatEndTime,
												ARM_GP_Vector&	swapLongFixPayTimes,
												ARM_GP_Vector&	swapLongFixPayPeriods,
												double&			swapShortFloatStartTime,
												double&			swapShortFloatEndTime,
												ARM_GP_Vector&	swapShortFixPayTimes,
												ARM_GP_Vector&	swapShortFixPayPeriods)  const 
{
	// Validations ...
	size_t soSize = spreadOption->GetResetDates()->GetSize();

	if (periodIdx<0 || periodIdx>=soSize)
		ARM_THROW( ERR_INVALID_ARGUMENT, "GetSpreadOptionLetData : invalid periodIdx" );

	
	// convert into vanilla arg
	ARM_VanillaArg* arg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(spreadOption, GetAsOfDate().GetJulian()) ;

	// at this step we are sure that arg is of type ARM_VanillaSpreadOptionArg
	// anyway, a dynamic_cast would be more secure...
	ARM_VanillaSpreadOptionArg* soArg = (ARM_VanillaSpreadOptionArg*)arg;


	/// compute schedules
	ARM_Currency* ccy = GetNumericalModel()->GetCurrency(GetNumericalModel()->GetModelName());
	
	soArg->ComputeIndexSchedulesAndAdjustDates(ccy, GetNumericalModel()->GetAsOfDate().GetJulian());
	
	payTime = (*soArg->GetPayTimes())[periodIdx];
		
	// compute target spread vol
	resetTime		= (*soArg->GetResetTimes())[periodIdx] ;
	coeffLong		= (*soArg->GetCoeffLong())[periodIdx];
	coeffShort		= (*soArg->GetCoeffShort())[periodIdx];
	strike			= (*soArg->GetStrikes())[periodIdx];
	payPeriod		= (*soArg->GetPayPeriods())[periodIdx];
	payNotional		= (*soArg->GetNotional())[periodIdx];
	callPut			= soArg->GetCallPut();
	

	swapLongFixPayTimes   = *soArg->GetSwapLongFixPayTimes()[periodIdx];
	swapLongFixPayPeriods = *soArg->GetSwapLongFixPayPeriods()[periodIdx];
	swapLongFloatStartTime = soArg->GetSwapLongFloatStartTime()->Elt(periodIdx);
	swapLongFloatEndTime   = soArg->GetSwapLongFloatEndTime()->Elt(periodIdx);


	swapShortFixPayTimes   = *soArg->GetSwapShortFixPayTimes()[periodIdx];
	swapShortFixPayPeriods = *soArg->GetSwapShortFixPayPeriods()[periodIdx];
	swapShortFloatStartTime = soArg->GetSwapShortFloatStartTime()->Elt(periodIdx);
	swapShortFloatEndTime   = soArg->GetSwapShortFloatEndTime()->Elt(periodIdx);

	/// don't forget to free memory
	delete soArg;

}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : GetCorridorSpreadOptionData
///	Returns : void
///	Action  : collects all needed info from ARM_SpreadOption which
/// is a corridor spread option in fact
////////////////////////////////////////////////////
void ARM_Local_Model::GetCorridorSpreadOptionData (	ARM_SpreadOption* spreadOption,
										// outputs
										int&			callPut,
										ARM_GP_Vector&	resetTimes,
										ARM_GP_Vector&	payTimes,
										ARM_GP_Vector&	payPeriods,
										ARM_GP_Vector&	payNotionals,
										ARM_GP_Vector&	coeffsLong,
										ARM_GP_Vector&	coeffsShort,
										ARM_GP_Vector&	strikes,
										ARM_GP_Vector&	swapsLongFloatStartTimes,
										ARM_GP_Vector&	swapsLongFloatEndTimes,
										vector<ARM_GP_Vector>&	swapsLongFixPayTimes,
										vector<ARM_GP_Vector>&	swapsLongFixPayPeriods,
										ARM_GP_Vector&	swapsShortFloatStartTimes,
										ARM_GP_Vector&	swapsShortFloatEndTimes,
										vector<ARM_GP_Vector>&	swapsShortFixPayTimes,
										vector<ARM_GP_Vector>&	swapsShortFixPayPeriods,
										ARM_GP_Vector& fixValues,
										double& spread1,
										double& spread2,
										int& payIndexType,
										ARM_GP_Vector& payIndexLeverages,
										ARM_GP_Vector& swapsPayFloatStartTimes,
										ARM_GP_Vector& swapsPayFloatEndTimes,
										vector<ARM_GP_Vector>& swapsPayFixPayTimes,
										vector<ARM_GP_Vector>& swapsPayFixPayPeriods,
										ARM_GP_Vector& payIndexResetTimes,
										ARM_IntVector& periodIndex,
										int& rateCallPut,
										ARM_GP_Vector& rateStrikes,
										ARM_GP_Vector&	rateFloatStartTimes,
										ARM_GP_Vector&	rateFloatEndTimes,
										vector<ARM_GP_Vector>&	rateFixPayTimes,
										vector<ARM_GP_Vector>&	rateFixPayPeriods,
										double& rateSpread1,double& rateSpread2)  const 
{
	// convert into vanilla arg
	ARM_VanillaArg* arg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(spreadOption, GetAsOfDate().GetJulian()) ;

	// at this step we are sure that arg is of type ARM_VanillaSpreadOptionArg
	// anyway, a dynamic_cast would be more secure...
	ARM_VanillaSpreadOptionArg* soArg = (ARM_VanillaSpreadOptionArg*)arg;

	size_t i,nbFlows = spreadOption->GetNumFlows();
	ARM_Currency* ccy = GetNumericalModel()->GetCurrency(GetNumericalModel()->GetModelName());

	/// Double condition case
	ARM_VanillaSpreadOptionArg* rateArg=NULL;
	bool isDbleCondition = spreadOption->IsCorridorDblCondition();
	int spreadCallPut=K_CALL;
	if(isDbleCondition)
	{
		ARM_CorridorDblCondition* dblCorridor=static_cast<ARM_CorridorDblCondition*>(spreadOption);
		rateArg = static_cast<ARM_VanillaSpreadOptionArg*>( ARM_ConverterFromKernel::ConvertSecuritytoArgObject(dblCorridor->GetSpreadDigital(), GetAsOfDate().GetJulian()) );
		rateArg->ComputeIndexSchedulesAndAdjustDates(ccy, GetNumericalModel()->GetAsOfDate().GetJulian());
		rateFloatStartTimes.resize(nbFlows);
		rateFloatEndTimes.resize(nbFlows);
		rateFixPayTimes.resize(nbFlows);
		rateFixPayPeriods.resize(nbFlows);
		rateStrikes.resize(nbFlows);
		for (i = 0; i < nbFlows; ++i)
		{
			rateFixPayTimes[i]		= *rateArg->GetSwapLongFixPayTimes()[i];
			rateFixPayPeriods[i]	= *rateArg->GetSwapLongFixPayPeriods()[i];
			rateFloatStartTimes[i]	= (*rateArg->GetSwapLongFloatStartTime())[i];
			rateFloatEndTimes[i]	= (*rateArg->GetSwapLongFloatEndTime())[i];
			rateStrikes[i]			= (*rateArg->GetStrikes())[i];
		}

		rateCallPut		= dblCorridor->GetDigitalCapFloor();
		spreadCallPut	= dblCorridor->GetSpreadCapFloor();

		rateSpread1 = rateArg->GetSpread1();
		rateSpread2 = rateArg->GetSpread2();

		delete rateArg;
	}

	/// compute schedules
	
	soArg->ComputeIndexSchedulesAndAdjustDates(ccy, GetNumericalModel()->GetAsOfDate().GetJulian());


	resetTimes.resize(nbFlows);
	payTimes.resize(nbFlows);
	payPeriods.resize(nbFlows);
	payNotionals.resize(nbFlows);
	coeffsLong.resize(nbFlows);
	coeffsShort.resize(nbFlows);
	strikes.resize(nbFlows);
	swapsLongFloatStartTimes.resize(nbFlows);
	swapsLongFloatEndTimes.resize(nbFlows);
	swapsLongFixPayTimes.resize(nbFlows);
	swapsLongFixPayPeriods.resize(nbFlows);
	swapsShortFloatStartTimes.resize(nbFlows);
	swapsShortFloatEndTimes.resize(nbFlows);
	swapsShortFixPayTimes.resize(nbFlows);
	swapsShortFixPayPeriods.resize(nbFlows);
	
	payIndexType = soArg->GetPayIndexType();

	size_t nbPeriods = soArg->GetNumPeriods();
	
	
	swapsPayFloatStartTimes.resize(nbPeriods);
	swapsPayFloatEndTimes.resize(nbPeriods);
	swapsPayFixPayTimes.resize(nbPeriods);
	swapsPayFixPayPeriods.resize(nbPeriods);
	payIndexLeverages.resize(nbPeriods);
	payIndexResetTimes.resize(nbPeriods);

	periodIndex.resize(nbFlows);
	

	fixValues.resize(nbFlows);

	/// Becareful for double condition rateArg & soArg are always set to K_CALL
	/// then C/F type is set in the DblCorridor object
	callPut	= (isDbleCondition ? spreadCallPut : soArg->GetCallPut());

	spread1 = soArg->GetSpread1();
	spread2 = soArg->GetSpread2();



	for (i = 0; i < nbFlows; ++i)
	{
		payTimes[i] = (*soArg->GetPayTimes())[i];
		
		// compute target spread vol
		resetTimes[i]		= (*soArg->GetResetTimes())[i] ;
		coeffsLong[i]		= (*soArg->GetCoeffLong())[i];
		coeffsShort[i]		= (*soArg->GetCoeffShort())[i];
		strikes[i]			= (*soArg->GetStrikes())[i];
		payPeriods[i]		= (*soArg->GetPayPeriods())[i];
		payNotionals[i]		= (*soArg->GetNotional())[i];

		swapsLongFixPayTimes[i]		= *soArg->GetSwapLongFixPayTimes()[i];
		swapsLongFixPayPeriods[i]	= *soArg->GetSwapLongFixPayPeriods()[i];
		swapsLongFloatStartTimes[i]	= soArg->GetSwapLongFloatStartTime()->Elt(i);
		swapsLongFloatEndTimes[i]	= soArg->GetSwapLongFloatEndTime()->Elt(i);


		swapsShortFixPayTimes[i]	= *soArg->GetSwapShortFixPayTimes()[i];
		swapsShortFixPayPeriods[i]	= *soArg->GetSwapShortFixPayPeriods()[i];
		swapsShortFloatStartTimes[i]= soArg->GetSwapShortFloatStartTime()->Elt(i);
		swapsShortFloatEndTimes[i]	= soArg->GetSwapShortFloatEndTime()->Elt(i);

		fixValues[i]				= soArg->GetFixValues()->Elt(i);

		periodIndex[i]				= soArg->GetPeriodIndex()->Elt(i);
	
	}

	for (i = 0; i < nbPeriods; ++i)
	{
		swapsPayFixPayTimes[i]		= *soArg->GetSwapPayFixPayTimes()[i];
		swapsPayFixPayPeriods[i]	= *soArg->GetSwapPayFixPayPeriods()[i];
		swapsPayFloatStartTimes[i]	= soArg->GetSwapPayFloatStartTime()->Elt(i);
		swapsPayFloatEndTimes[i]	= soArg->GetSwapPayFloatEndTime()->Elt(i);
		payIndexLeverages[i]		= soArg->GetPayIndexLeverages()->Elt(i);
		payIndexResetTimes[i]		= soArg->GetPayIndexResetTimes()->Elt(i);
	}

	/// don't forget to free memory
	delete soArg;
}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routine : GetCorridorlegData
///	Returns : void
///	Action  : collects all needed info from ARM_CorridorLeg
////////////////////////////////////////////////////
void ARM_Local_Model::GetCorridorLegData (	ARM_CorridorLeg* corridor,
											// outputs
											int&			callPut,
											double&			resetTime,
											double&			payTime,
											int&			indexPaymentType,
											double&			fwdPaymentPeriod,
											ARM_GP_Vector&	refIdxResettimes ,
											ARM_GP_Vector&	refIdxStarttimes,
											ARM_GP_Vector&	refIdxEndtimes,
											ARM_GP_Vector&	refFwdPeriods,
											ARM_GP_Vector&	refIndexWeight,
											ARM_GP_Vector&	DownBarrier,
											ARM_GP_Vector&	UpBarrier,
											double&			couponMargin,
											double&			payNotional)  const 
{
	// Validations ...	
	// check size == 1
	if (corridor->GetResetDates()->GetSize() != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Model Calibration : ARM_CorridorLeg is required to have only one optionlet" );

	
	// convert into vanilla arg
	ARM_VanillaArg* arg = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(corridor, GetAsOfDate().GetJulian()) ;

	// at this step we are sure that arg is of type ARM_VanillaSpreadOptionArg
	// anyway, a dynamic_cast would be more secure...
	ARM_VanillaCorridorLegArg* corridorArg = (ARM_VanillaCorridorLegArg*)arg;


	/// compute schedules
	ARM_Currency* ccy = GetNumericalModel()->GetCurrency(GetNumericalModel()->GetModelName());
	
	//soArg->ComputeIndexSchedulesAndAdjustDates(ccy, GetNumericalModel()->GetAsOfDate().GetJulian());


	// some validation...
	if (corridorArg->GetResetTimes()->size()  != 1)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Local Model Calibration : ARM_CorridorLeg is required to have only one optionlet" );
	
	payTime = (*corridorArg->GetPayTimes())[0];
		
	resetTime			= (*corridorArg->GetResetTimes())[0] ;
	couponMargin		= (*corridorArg->GetCouponMargin())[0];
	fwdPaymentPeriod	= (*corridorArg->GetFwdPaymentPeriod())[0];
	payNotional			= (*corridorArg->GetNominals())[0];
	callPut				= corridorArg->GetCallPut();
	indexPaymentType	= corridorArg->GetIndexPaymentType();
	
	refIdxResettimes	=*corridorArg->GetRefIdxResettimes()[0];
	refIdxStarttimes	=*corridorArg->GetRefIdxStarttimes()[0];
	refIdxEndtimes		=*corridorArg->GetRefIdxEndtimes()[0];
	refFwdPeriods		=*corridorArg->GetRefFwdPeriods()[0];
	refIndexWeight		=*corridorArg->GetRefIndexWeight()[0];
	DownBarrier			=*corridorArg->GetDownBarrier()[0];
	UpBarrier			=*corridorArg->GetUpBarrier()[0];

	/// don't forget to free memory
	delete corridorArg;

}

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routines: GetEqFxOptionData
///	Returns : void
///	Action  : collects all needed info from ARM_SpreadOption
////////////////////////////////////////////////////
void ARM_Local_Model::GetEqFxOptionData(    ARM_Option* option,
											// outputs
									        double& resetTime,
									        double& settlementTime,
									        double& payTime,
									        double& strike,
									        int&	callPut) const
{
	// convert into vanilla arg
	ARM_VanillaEqFxOption* eqFxOpt = dynamic_cast<ARM_VanillaEqFxOption*>( ARM_ConverterFromKernel::ConvertSecuritytoArgObject(option, GetAsOfDate().GetJulian()) );
    if(!eqFxOpt)
		ARM_THROW( ERR_INVALID_ARGUMENT, "Kernel to GP converter failed on an ARM_Option" );

    resetTime       = eqFxOpt->GetExpiry();
    settlementTime  = eqFxOpt->GetFwdTime();
    payTime         = eqFxOpt->GetPayTime();
    strike          = eqFxOpt->GetStrike();
    callPut         = eqFxOpt->GetCallPut();

	delete eqFxOpt;
}
 

////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routines: DiscountFactor
///	Returns : 
///	Action  : returns numerical model DF
////////////////////////////////////////////////////
ARM_VectorPtr ARM_Local_Model::DiscountFactor(	const string& curveName,
												double evalTime, 
												double maturityTime,
												const ARM_PricingStatesPtr& states) const
{	
	return itsNumericalModel->DiscountFactor(curveName, evalTime, maturityTime, states);
}


////////////////////////////////////////////////////
///	Class   : ARM_Local_Model
///	Routines: toString
///	Returns : 
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_Local_Model::toString(const string& indent,const string& nextIndent) const
{
    CC_Ostringstream os;
    os << "\n\n";
    os << indent << "Local Model \n";
    os << indent << "---------------------\n";
    os << ARM_PricingModel::toString(indent);
    return os.str();
}


////////////////////////////////////////////////////
///	Class   : ARM_BS_Model
///	Routine : ValidateModelParams
///	Returns : bool
///	Action  : validate model params
////////////////////////////////////////////////////
bool ARM_Local_Model::ValidateModelParams(const ARM_ModelParams& params) const
{	
	/*
	if( !dynamic_cast<const ARM_AnalyticModelParams*>(&params) )
		ARM_THROW( ERR_INVALID_ARGUMENT, "modelparams is not of type ARM_AnalyticModelParams" );
	*/
	return true;
}


CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

