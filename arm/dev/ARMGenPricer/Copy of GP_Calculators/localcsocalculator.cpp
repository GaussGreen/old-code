/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file localcsocalculator.cpp
 *
 *  \brief file for the calculator for Callable SpreadOption
 *
 *	\author  JP Riaudel, Y. Khlif, A. Chaix
 *	\version 1.0
 *	\date July 2005
 */


#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/localcsocalculator.h"


/// gpbase
#include "gpbase/env.h"
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/gpmatrixlinalg.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/modelnamemap.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/mktdatamanagerrep.h"

/// gpmodels
#include "gpmodels/hw1f.h"
#include "gpmodels/hw2f.h"
#include "gpmodels/hwfactory.h"
#include "gpmodels/QGM1F.h"
#include "gpmodels/ModelParamsQGM1F.h"
#include "gpmodels/modelparamshw1f.h"
#include "gpmodels/modelparamshw2f.h"
#include "gpmodels/local_normal_model.h"
#include "gpmodels/local_normal_modelparams.h"
#include "gpmodels/multiassets.h"
#include "gpmodels/marketirmodel.h"
#include "gpmodels/ForwardMarginBasis.h"
#include "gpmodels/ForwardMarginIR.h"
#include "gpmodels/forwardforex.h"


/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaspreadoption.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/stripper.h"
#include "gpcalib/modelfitterdes.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_normal.h"


/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

#include <inst/swaption.h>
#include "inst/spreadoption.h"
#include <util/fromto.h>
#include <inst/fixleg.h>
#include <inst/swap.h>
#include <inst/swapleg.h>
#include <mod/bssmiled.h>

using ARM::ARM_MarketData_ManagerRep;

CC_BEGIN_NAMESPACE( ARM )

/// Reference schedules for CSO date structure
const unsigned int EXERCISE_SCHED	= 0;
const unsigned int NB_CSO_SCHED		= 1;

const double DEFAULT_PRICE			= 1.0e+100;
const double DEFAULT_PRECISION		= 1.0e-6;
const double DEFAULT_WEIGHT			= 1.0;

const double SENSI_CMS				= 1.0e-5;

const double NO_FLOOR				= -0.99;
const double NO_CAP					=  0.99;

const double SIGMA_LOWER_BOUND      = 0.0005; // 5bps
const double SIGMA_UPPER_BOUND      = 0.5; // 500bps
const double SIGMA_DEFAULT_VALUE	= 0.01;

const double VOLRATIO_DEFAULT_VALUE			= 1.;
const double CORREL_DEFAULT_VALUE			= 0.;

/// Tree using by default 30 time steps per year
const int DEFAULT_TREE_NBSTEPS_PER_YEAR		= 15;

/// QGM Skew range [-15,40] with a 0 (<=>H&W) default value
const double SKEW_LOWER_BOUND      = 1.0e-10;
const double SKEW_UPPER_BOUND      = 150.0;
const double SKEW_DEFAULT_VALUE    = 25.0;

/// QGM Default Parameters for the calibration
const double SKEW_CALIB_PRECISION   = 1.0e-6;
const double VOL_CALIB_PRECISION    = 1.0e-6;
const size_t MAXITER				= 30;
const size_t NB_ITER_MAX  = 5;

const string LOCAL_FLOOR_MODEL_NAME  = "LOCFLOOR";
const string LOCAL_CAP_MODEL_NAME	 = "LOCCAP";

/// for Probas
const unsigned int NB_CALCULATED_PROBA      = 21;

// Vega to select swaption in the calibration portfolio
const double VEGA_MIN_TO_SELECT = 1e-12;

/// Default MDM key names


const string ARM_LocalCSOCalculator::LocalCSOColNamesTable [] =
{
    "EventDate",
	"StartDate",
	"EndDate",
	"FixFundLeg",
    "FloorFundLeg",
    "CapFundLeg",
    "FundingLeg",
    "FixLeg",
    "FloorLeg",
    "CapLeg",
    "SOLeg",
    "Fees",
	"Option",
	"CSO",
	"Funding",
    "Fix",
	"Floor",
	"Cap",
	"ExerSwapRate",
	"Frontier"
};

const int ARM_LocalCSOCalculator::LocalCSOProductToPriceColums [] =
{
    CSO,
    Funding,
    Fix,
    Floor,
	Cap
};

const string ARM_LocalCSOCalculator::LocalCSOProbaColNamesTable [] =
{
	"ResetPlus1",
	"FinalDate",
	"ExerciseIndicator",
	"Prob1",
	"Prob10",
	"Prob2",
	"Prob20",
	"Prob3",
	"Prob30",
	"Prob4",
	"Prob40",
	"Prob5",
	"Prob50",
	"Prob6",
	"Prob60",
	"Prob7",
	"Prob70",
	"Prob8",
	"Prob80",
	"Prob9",
	"Prob90",
	"Proba10",
	"Proba100",
	"Prob11",
	"Prob110",
	"Prob12",
	"Prob120",
	"Prob13",
	"Prob130",
	"Prob14",
	"Prob140",
	"Prob15",
	"Prob150",
	"Prob16",
	"Prob160",
	"Prob17",
	"Prob170",
	"Prob18",
	"Prob180",
	"Prob19",
	"Prob190",
	"Proba20",
	"Proba200",
	"Proba21",
	"Proba210",
};

const string ARM_LocalCSOCalculator::LocalCSOProbaPricedColNamesTable [] =
{
	"Prob10",
	"Prob20",
	"Prob30",
	"Prob40",
	"Prob50",
	"Prob60",
	"Prob70",
	"Prob80",
	"Prob90",
	"Proba100",
	"Prob110",
	"Prob120",
	"Prob130",
	"Prob140",
	"Prob150",
	"Prob160",
	"Prob170",
	"Prob180",
	"Prob190",
	"Proba200",
	"Proba210",
};

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: constructor (contextual)
///	Returns: void
///	Action : builds a ARM_LocalCSOCalculator
/////////////////////////////////////////////////////////////////
ARM_LocalCSOCalculator::ARM_LocalCSOCalculator(const ARM_Date&	startDate,
											   const ARM_Date&	endDate,
											   int CMSLong,
											   int CMSShort,
											   int cpnDayCount,
											   int cpnFreq,
											   int cpnResetTiming,
											   const ARM_Curve& cpnnominal,
											   const ARM_Curve& fixCoupon,
											   const ARM_Curve& leverageLong,
											   const ARM_Curve& leverageShort,
											   const ARM_Curve& cpnMin,
											   const ARM_Curve& cpnMax,
											   const ARM_Curve& strike,
											   int fundFreq,
											   int fundDaycount,
											   const ARM_Curve& fundSpread,
											   const ARM_Curve& fundLeverage,
											   int exerciseFreq,
											   int noticeGap,
											   int payRec,
											   const ARM_Curve& fees,
											   bool switchFlag,
											   int fundingType,
											   std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/	productsToPrice,
											   vector<double> ModelDatas,
											   const ARM_MarketData_ManagerRep& mktDataManager,
											   const ARM_StringVector& mdmKeys,
											   ARM_ModelType modelType,
											   CalibrationType	calibrationType,
											   vector<CalibStrikeType>	calibStrikeType,
											   ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod)											   
:	ARM_GenCSOCalculator  (startDate,
						   endDate,
						   CMSLong,
						   CMSShort,
						   cpnDayCount,
						   cpnFreq,
						   cpnResetTiming,
						   cpnnominal,
						   fixCoupon,
						   leverageLong,
						   leverageShort,
						   cpnMin,
						   cpnMax,
						   strike,
						   fundFreq,
						   fundDaycount,
						   fundSpread,
						   fundLeverage,
						   exerciseFreq,
						   noticeGap,
						   payRec,
						   fees,
						   productsToPrice,
						   ModelDatas,
						   mktDataManager,
						   mdmKeys ) ,

	itsModelType		(modelType),
	itsCalibrationType	(calibrationType),
	itsCalibStrikeType	(calibStrikeType),
	itsCpnLongSensi		(0),
	itsCpnShortSensi	(0),
	itsCpnLongSensiForSwitch	(0),
	itsCpnShortSensiForSwitch	(0),
	itsCpnValue			(0),
	itsCpnValueForSwitch		(0),
	itsSpreadOptionFloorPF (ARM_StdPortfolioPtr(NULL)),
	itsSpreadOptionCapPF   (ARM_StdPortfolioPtr(NULL)),
	itsSpreadOptionFloorPFforSwitch (ARM_StdPortfolioPtr(NULL)),
	itsSpreadOptionCapPFforSwitch   (ARM_StdPortfolioPtr(NULL)),
	itsSoPricesComputed (false),
	itsVnsPricingMethod (vnsPricingMethod),
	itsMoyenessLevel(1.0),
	itsComputeProba(false),
	itsFirstComputedProba(0),
	itsNbComputedProba(0),
	itsProba(ARM_VectorPtr(NULL))
{
	// switchable features
	itsSwitchFlag = switchFlag;
	if (IsSwitchable())
	{
		if (CMSLong==fundingType || CMSShort==fundingType)
			itsCMSFunding = fundingType;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : funding should be one of the 2 CMS paid");
	}

    /// Check input datas
    CheckDataAndTimeIt();
	
	DatesStructure();

	/// nb fix flows
	ARM_GP_Vector* resetDates = GetCalibDateStrip()->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int resetSize=resetDates->size();
	int i=0;
	while( (i < resetSize) && (const_cast< ARM_Curve& >(GetFixCpn()).Interpolate((*resetDates)[i]-asOfDate) > 0.0) )
        ++i;
	SetNbFixFlows(i);

	// probas
	int stdSize = itsProductsToPrice.size()-1;
	
	ARM_GP_Vector* callDates = GetExerciseDateStripUnadj()->GetResetDates();
	itsComputeProba		= productsToPrice[stdSize];
	itsNbComputedProba	= (NB_CALCULATED_PROBA<callDates->size() )?NB_CALCULATED_PROBA:callDates->size();
	itsProba			= ARM_VectorPtr(new ARM_GP_Vector(itsNbComputedProba,0.0));
	bool otherPayoffs	= itsComputeProba;

	// compute vectors from curve
	ComputeProductVectorsFromCurves();
	
	/// get all the product to price and price in multi-column
	ARM_StringVector pricedColumns;
	int prodSize	= MIN(NbProductsToPrice,stdSize);
	int size		= prodSize;

	if (itsComputeProba)
		size += itsNbComputedProba;

	bool isFrontierCalib = (itsCalibrationType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER);
	if(isFrontierCalib)
		++size;

	pricedColumns.reserve(size);

    for (i = 0; i < prodSize; ++i)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back( LocalCSOColNamesTable[ LocalCSOProductToPriceColums[i] ] );
	}

	if (itsComputeProba)
	{
		for (i = prodSize; i < prodSize+itsNbComputedProba; ++i)
		{
			pricedColumns.push_back( LocalCSOProbaColNamesTable[4+2*(i-prodSize) ] );
		}
	}

	if(isFrontierCalib)
		pricedColumns.push_back( LocalCSOColNamesTable[Frontier] );

	/// Create Cst manager
	ARM_CstManagerPtr cstMgr = CreateCstManager();

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();
	
	/// Create the Generic Security paid in coupon currency
	CreateAndSetDealDescriptionAndTimeIt(GetKeys()[GetCpnModelKey()], pricedColumns, cstMgr,false,otherPayoffs);
	
	/// Create the HW pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping
    CreateAndSetCalibrationAndTimeIt();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: constructor (summit constructor)
///	Returns: void
///	Action : builds a ARM_LocalCSOCalculator
///          without MARKET DATA and model/calib parameters
/////////////////////////////////////////////////////////////////
ARM_LocalCSOCalculator::ARM_LocalCSOCalculator(const ARM_Date&	asOfDate,
											   const ARM_Date&	startDate,
                                               const ARM_Date& fixEndDate,
											   const ARM_Date&	endDate,
											   const ARM_Currency& CpnCcy,
											   const ARM_Currency& FundCcy,
											   int CMSLong,
											   int CMSShort,
											   int cpnDayCount,
											   int cpnFreq,
											   int cpnResetTiming,
											   const ARM_Curve& cpnnominal,
											   const ARM_Curve& fixCoupon,
											   const ARM_Curve& leverageLong,
											   const ARM_Curve& leverageShort,
											   const ARM_Curve& cpnMin,
											   const ARM_Curve& cpnMax,
											   const ARM_Curve& strike,
											   int fundFreq,
											   int fundDaycount,
											   const ARM_Curve& fundNominal,
											   const ARM_Curve& fundSpread,
											   const ARM_Curve& fundLeverage,
											   int exerciseFreq,
											   int noticeGap,
											   int payRec,
											   const ARM_Curve& fees,
											   bool switchFlag,
											   int fundingType,
											   int nbNoCall)
:	ARM_GenCSOCalculator(asOfDate,
						 startDate,
                         fixEndDate, // fixEnDate in the case of Initial Fix!
						 endDate,
						 CMSLong,
						 CMSShort,
						 cpnDayCount,
						 cpnFreq,
						 cpnResetTiming,
						 cpnnominal,
						 fixCoupon,
						 leverageLong,
						 leverageShort,
						 cpnMin,
						 cpnMax,
						 strike,
						 fundFreq,
						 fundDaycount,
						 fundSpread,
						 fundLeverage,
						 exerciseFreq,
						 noticeGap,
						 payRec,
						 fees),
	                     itsCpnLongSensi	(0),
	                     itsCpnShortSensi	(0),
						 itsCpnLongSensiForSwitch	(0),
	                     itsCpnShortSensiForSwitch	(0),
	                     itsCpnValue		(0),
						 itsCpnValueForSwitch	(0),
	                     itsSpreadOptionFloorPF (ARM_StdPortfolioPtr(NULL)),
	                     itsSpreadOptionCapPF   (ARM_StdPortfolioPtr(NULL)),
						 itsSpreadOptionFloorPFforSwitch (ARM_StdPortfolioPtr(NULL)),
	                     itsSpreadOptionCapPFforSwitch   (ARM_StdPortfolioPtr(NULL)),
	                     itsSoPricesComputed (false),
	                     itsMoyenessLevel(1.0),
	                     itsComputeProba(false),
	                     itsFirstComputedProba(0),
	                     itsNbComputedProba(0),
	                     itsProba(ARM_VectorPtr(NULL))
{
	// switchable features
	itsSwitchFlag = switchFlag;
	if (IsSwitchable())
	{
		if (CMSLong==fundingType || CMSShort==fundingType)
			itsCMSFunding = fundingType;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : funding should be one of the 2 CMS paid");
	}

	///Set currencies
	SetCurrencyUnit(const_cast< ARM_Currency* >(&CpnCcy));
	SetFundingCcy(const_cast< ARM_Currency& >(FundCcy));
	
	/// nb fix flows
	for(int i(0); i <nbNoCall; ++i)
		const_cast<ARM_Curve&>(GetFees()).GetOrdinates()[i] = NON_CALL_FEE;

	SetFundNominal(fundNominal);

	string cpnCcyName(CpnCcy.GetCcyName());
	string fundingCcyName(FundCcy.GetCcyName());
	ARM_StringVector keys(NbKeys,UNKNOWN_KEY_NAME);
 
    keys[YcKey]			= YC_KEY_NAME + cpnCcyName;
    keys[MrsKey]		= MRS_KEY_NAME + cpnCcyName;
	keys[CorrelKey]		= CORREL_KEY_NAME + cpnCcyName;
    keys[CfModelKey]	= CFMODEL_KEY_NAME + cpnCcyName;
    keys[OswModelKey]	= OSWMODEL_KEY_NAME + cpnCcyName;
	keys[SoModelKey]	= SOMODEL_KEY_NAME + cpnCcyName;
	keys[VolRatioKey]	= VOLRATIO_KEY_NAME + cpnCcyName;
	keys[MrsSpreadKey]	= MRSSPREAD_KEY_NAME + cpnCcyName;
	keys[FundingKey]	= keys[YcKey];
	keys[BasisKey]		= keys[YcKey];
	keys[fundBasisKey]  = keys[YcKey];

	if (IsBasis())
	{
		keys[FundingKey]= YC_KEY_NAME + fundingCcyName;
		keys[BasisKey]= YC_BASIS_KEY_NAME + cpnCcyName;
		keys[fundBasisKey]= YC_BASIS_KEY_NAME + fundingCcyName;
		keys[ForexKey]= FOREX_KEY_NAME + cpnCcyName +"/"+ fundingCcyName;
	}

	SetKeys(keys);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: InitCSOFromSummit
///	Returns: void
///	Action : initialize a ARM_LocalCSOCalculator
///          with market data and model/calib parameters
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::InitCSOFromSummit(const ARM_MarketData_ManagerRep* mktDataManager,
											   std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
											   vector<double> modelDatas,
											   CalibrationType	calibrationType,
											   vector<CalibStrikeType>	calibStrikeType,
											   ARM_ModelType modelType,
											   ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod)
{

	InitializeMktDataManagerOnly(*mktDataManager);

	itsModelDatas       = modelDatas;
	itsCalibrationType  = calibrationType;
	itsCalibStrikeType  = calibStrikeType;
	itsModelType        = modelType;
	itsVnsPricingMethod = vnsPricingMethod;

	itsProductsToPrice  = productsToPrice;

    // Check input datas
    CheckDataAndTimeIt();
	
    // Schedules initialization
	DatesStructure();

    // fix period managing by fixEndDate
    int NbFixFlows = 0;

	while ( NbFixFlows < itsStructDateStrip->GetFlowEndDates()->size() 
            && 
			(*itsStructDateStrip->GetFlowEndDates())[NbFixFlows] < GetFixEndDate().GetJulian()+7.0 /* ARM_GlobalConstant::ARM_SEVENDAYS_LAG */)
		  NbFixFlows++;

    SetNbFixFlows(NbFixFlows);

	// probas
	int stdSize = itsProductsToPrice.size()-1;
	
	ARM_GP_Vector* callDates = GetExerciseDateStripUnadj()->GetResetDates();
	itsComputeProba		= productsToPrice[stdSize];
	itsNbComputedProba	= (NB_CALCULATED_PROBA<callDates->size() )?NB_CALCULATED_PROBA:callDates->size();
	itsProba			= ARM_VectorPtr(new ARM_GP_Vector(itsNbComputedProba,0.0));
	bool otherPayoffs	= itsComputeProba;

	// compute vectors from curve
	ComputeProductVectorsFromCurves();
	
	/// get all the product to price and price in multi-column
	ARM_StringVector pricedColumns;
	int prodSize	= MIN(NbProductsToPrice,stdSize);
	int size		= prodSize;

	if (itsComputeProba)
		size += itsNbComputedProba;

	bool isFrontierCalib = (itsCalibrationType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER);
	if(isFrontierCalib)
		++size;

	pricedColumns.reserve(size);

    int i;

    for (i = 0; i < prodSize; ++i)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back( LocalCSOColNamesTable[ LocalCSOProductToPriceColums[i] ] );
	}

	if (itsComputeProba)
	{
		for (i = prodSize; i < prodSize+itsNbComputedProba; ++i)
		{
			pricedColumns.push_back( LocalCSOProbaColNamesTable[4+2*(i-stdSize) ] );
		}
	}

	if(isFrontierCalib)
		pricedColumns.push_back( LocalCSOColNamesTable[Frontier] );

	/// Create Cst manager
	ARM_CstManagerPtr cstMgr = CreateCstManager();

	if (IsBasis())
	{
		/// Set Foreign, basis and domestic Currencies
		ARM_Forex* forex = static_cast< ARM_Forex* >(mktDataManager->GetData(GetKeys()[ForexKey]));
		SetForeignCcy(*forex->GetMainCurrency());
		SetDomesticCcy(*forex->GetMoneyCurrency());

		ARM_Currency BasisCcy = *static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(GetKeys()[BasisKey]))->GetCurrencyUnit();
		SetBasisCcy(BasisCcy);
	}

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();
	
	/// Create the Generic Security paid in coupon currency
	CreateAndSetDealDescriptionAndTimeIt(GetKeys()[GetCpnModelKey()], pricedColumns, cstMgr);
	
	/// Create the HW pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping
    CreateAndSetCalibrationAndTimeIt();
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: InitCSOFromSummit
///	Returns: void
///	Action : initialize a ARM_LocalCSOCalculator
///          with market data and model/calib parameters
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::InitCSOFromSummit(ARM_ZeroCurve* zc,
											   ARM_VolCurve* capVol,
											   ARM_VolCurve* swoptVol,
											   int SABRSigmaOrAlpha,
											   ARM_VolCurve* rhoCap,
											   ARM_VolCurve* nuCap,
											   ARM_VolCurve* betaCap,
											   ARM_VolCurve* rhoSwopt,
											   ARM_VolCurve* nuSwopt,
											   ARM_VolCurve* betaSwopt,
											   ARM_VolCurve* flatVol,
											   ARM_VolCurve* convAdjustVol,
											   ARM_ConvAdjustManager* convAdjustManager,
											   ARM_VolCube* correlDiagCap,
											   ARM_VolCube* correlDiagSwopt,
											   ARM_CorrelManager* correlCorr,
											   ARM_CurveModelParam* mrs,
											   ARM_CurveModelParam* correl,
											   ARM_CurveModelParam* volRatio,
											   ARM_CurveModelParam* mrsSpread,
											   std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
											   vector<double> modelDatas,
											   CalibrationType	calibrationType,
											   vector<CalibStrikeType>	calibStrikeType,
											   ARM_ModelType modelType,
											   ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
											   ARM_Forex* forex,
											   ARM_ZeroCurve* fundZc,
											   ARM_ZeroCurve* domBasisZc,
											   ARM_ZeroCurve* fundBasisZc)
{
    ARM_StringVector keys = GetKeys();

	string ccyName(zc->GetCurrencyUnit()->GetCcyName());

	//Market Data Manager ----------------------------------------------------------->
	ARM_MarketData_ManagerRep* marketDataManager = &*GetMktDataManager();
	
	ARM_Date asof = marketDataManager->GetAsOfDate();

	// SABR Model 
	ARM_BSSmiledModel* capSabrModel = NULL;
	ARM_VolCurve* ATMCapVol = (capVol->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)capVol)->GetATMVol() : capVol;

	if (rhoCap && nuCap)
	{
		int methodType;
	
        if (betaCap)
		{
           methodType = betaCap->IsEqualToOne()? K_SABR_ARITH:K_SABR_IMPLNVOL;

           capSabrModel = new ARM_BSSmiledModel(asof, 
											    0.0, 
											    zc, 
											    zc, 
											    ATMCapVol, 
											    K_YIELD, 
											    rhoCap, 
											    nuCap, 
											    methodType, 
											    betaCap,
                                                0.5,  // SABR Weight: Irrelevant
                                                SABRSigmaOrAlpha);
		}
        else
        {
           capSabrModel = new ARM_BSSmiledModel(asof, 
											    0.0, 
											    zc, 
											    zc, 
											    ATMCapVol, 
											    K_YIELD, 
											    rhoCap, 
											    nuCap, 
											    K_SABR_ARITH, 
											    NULL,
                                                0.5, // SABR Weight: Irrelevant
                                                SABRSigmaOrAlpha);
        }
	}

	// Cap BS Gen Model
	ARM_BSModel* capModel = new ARM_BSModel(asof, 
											zc, 
											flatVol, //spreadlock
											convAdjustVol,
											capVol,
											correlCorr,
											convAdjustManager,
											zc, //discount
											correlDiagCap,
											NULL, //cash vol
											NULL, //spread vol
											K_2LOG,
											K_COMPUTED,
											capSabrModel); 

	// Swopt Model
	ARM_BSSmiledModel* swoptSabrModel = NULL;
	ARM_VolCurve* ATMSwoptVol = (swoptVol->GetName() == ARM_VOL_CUBE) ? ((ARM_VolCube*)swoptVol)->GetATMVol() : swoptVol;

	if (rhoSwopt && nuSwopt)
	{
		int methodType;
	
        if (betaSwopt)
		{
            methodType = betaSwopt->IsEqualToOne()? K_SABR_ARITH : K_SABR_IMPLNVOL; 
    
            swoptSabrModel = new ARM_BSSmiledModel(asof, 
												   0.0, 
												   zc, 
												   zc, 
												   ATMSwoptVol, 
												   K_YIELD, 
												   rhoSwopt, 
												   nuSwopt, 
												   methodType, 
												   betaSwopt,
                                                   0.5, // SABR Weight: Irrelevant
                                                   SABRSigmaOrAlpha);
		
		}
        else
        {
		   swoptSabrModel = new ARM_BSSmiledModel(asof, 
												  0.0, 
												  zc, 
												  zc, 
												  ATMSwoptVol, 
												  K_YIELD, 
												  rhoSwopt, 
												  nuSwopt, 
												  K_SABR_ARITH,
                                                  NULL,
                                                  0.5, // SABR Weight: Irrelevant
                                                  SABRSigmaOrAlpha);
        }
	}

	ARM_BSModel* swoptModel = new ARM_BSModel(asof, 
											  zc, 
											  NULL, //spreadlock
											  convAdjustVol, //convAdjVol
											  swoptVol,
											  correlCorr,
											  convAdjustManager, //conv adjust manager
											  zc, //discount
											  correlDiagSwopt,
											  NULL, //cash vol
											  NULL, //spread vol
											  K_2LOG,
											  K_COMPUTED,
											  swoptSabrModel);

	// SO Model
	ARM_BSModel* soModel = (ARM_BSModel*) capModel->Clone();

	// Fill mktData Manager with models and curves
	marketDataManager->RegisterData(keys[YcKey], zc);
	marketDataManager->RegisterData(keys[CfModelKey], capModel);
	marketDataManager->RegisterData(keys[OswModelKey], swoptModel);
	marketDataManager->RegisterData(keys[SoModelKey], soModel);
	marketDataManager->RegisterData(keys[MrsKey], mrs);
	if (correl)
		marketDataManager->RegisterData(keys[CorrelKey], correl);
	if (volRatio)
		marketDataManager->RegisterData(keys[VolRatioKey], volRatio);
	if (mrsSpread)
		marketDataManager->RegisterData(keys[MrsSpreadKey], mrsSpread);

	if (IsBasis())
    {
		if (fundZc)
			marketDataManager->RegisterData(keys[FundingKey], fundZc);

		if (domBasisZc)
		{
			marketDataManager->RegisterData(keys[BasisKey], domBasisZc);
			/// set basis currency
			ARM_Currency BasisCcy = *(domBasisZc->GetCurrencyUnit());
			SetBasisCcy(BasisCcy);
		}

		if (fundBasisZc)
			marketDataManager->RegisterData(keys[fundBasisKey], fundBasisZc);

		if (forex)
		{
			marketDataManager->RegisterData(keys[ForexKey], forex);
			/// Set Foreign and domestic Currencies
			SetForeignCcy(*forex->GetMainCurrency());
			SetDomesticCcy(*forex->GetMoneyCurrency());
		}
	}
	
	InitCSOFromSummit(marketDataManager,
					  productsToPrice,
					  modelDatas,
					  calibrationType,
					  calibStrikeType,
					  modelType,
					  vnsPricingMethod);

	delete capModel;
	capModel = NULL;

	delete swoptModel;
	swoptModel = NULL;

	delete soModel;
	soModel = NULL;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: constructor (contextual)
///	Returns: void
///	Action : builds a ARM_LocalCSOCalculator for basis 
/////////////////////////////////////////////////////////////////
ARM_LocalCSOCalculator::ARM_LocalCSOCalculator(const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& CpnCcy,
		const ARM_Currency& FundCcy,
		int CMSLong,
		int CMSShort,
		int cpnDayCount,
		int cpnFreq,
		int cpnResetTiming,
		int stubRule,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& leverageLong,
		const ARM_Curve& leverageShort,
		const ARM_Curve& cpnMin,
		const ARM_Curve& cpnMax,
		const ARM_Curve& strikes,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		bool switchFlag,
		int fundingType,
		size_t nbNCall,
		const ARM_Curve& fees,
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
		vector<double> ModelDatas,
		const ARM_MarketData_ManagerRep& mktDataManager,
		ARM_ModelType		modelType,
		CalibrationType	calibrationType,
		vector<CalibStrikeType>	calibStrikeType,
		ARM_MarketIRModel::VnsPricingMethod vnsPricingMethod,
		double moyenessLevel)
 : 
	  ARM_GenCSOCalculator (startDate,fixEndDate,endDate,CpnCcy, FundCcy,CMSLong,CMSShort,
						cpnDayCount,cpnFreq,cpnResetTiming,stubRule,cpnResetCal,cpnPayCal,
						cpnnominal,leverageLong,leverageShort, cpnMin, cpnMax, strikes,fundFreq,
						fundDayCount,fundnominal,fundSpread, ARM_FlatCurve(1.0), exerciseFreq,
						noticeGap,payRec,nbNCall,fees,productsToPrice,ModelDatas,mktDataManager),
		itsModelType		(modelType),
		itsCalibrationType	(calibrationType),
		itsCalibStrikeType	(calibStrikeType),
		itsCpnLongSensi		(0),
		itsCpnShortSensi	(0),
		itsCpnLongSensiForSwitch	(0),
		itsCpnShortSensiForSwitch	(0),
		itsCpnValue			(0),
		itsCpnValueForSwitch		(0),
		itsSpreadOptionFloorPF (ARM_StdPortfolioPtr(NULL)),
		itsSpreadOptionCapPF   (ARM_StdPortfolioPtr(NULL)),
		itsSpreadOptionFloorPFforSwitch (ARM_StdPortfolioPtr(NULL)),
		itsSpreadOptionCapPFforSwitch   (ARM_StdPortfolioPtr(NULL)),
		itsSoPricesComputed (false),
		itsVnsPricingMethod (vnsPricingMethod),
		itsMoyenessLevel(moyenessLevel),
		itsComputeProba(false),
		itsFirstComputedProba(0),
		itsNbComputedProba(0),
		itsProba(ARM_VectorPtr(NULL))
{
	/// Register input objects but no more than the internal number of access keys
    size_t nbToReg = GetKeys().size();
    for(size_t i(0); i<nbToReg; ++i)
	{
		if(!mktDataManager.TestIfKeyMissing(GetKeys()[i]))
			GetMktDataManager()->RegisterData(GetKeys()[i],mktDataManager.GetData(GetKeys()[i]));
	}

	GetMktDataManager()->SetDetailMode(mktDataManager.GetDetailMode());

	// switchable features
	itsSwitchFlag = switchFlag;
	if (IsSwitchable())
	{
		if (CMSLong==fundingType || CMSShort==fundingType)
			itsCMSFunding = fundingType;
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : funding should be one of the 2 CMS paid");
	}

	/// Check input datas
    CheckDataAndTimeIt();

	DatesStructure();

	// probas
	int stdSize = itsProductsToPrice.size()-1;
	
	ARM_GP_Vector* callDates = GetExerciseDateStripUnadj()->GetResetDates();
	itsComputeProba		= productsToPrice[stdSize];
	itsNbComputedProba	= (NB_CALCULATED_PROBA<callDates->size() )?NB_CALCULATED_PROBA:callDates->size();
	itsProba			= ARM_VectorPtr(new ARM_GP_Vector(itsNbComputedProba,0.0));
	bool otherPayoffs	= itsComputeProba;


	// compute vectors from curve
	ComputeProductVectorsFromCurves();
	
	/// get all the product to price and price in multi-column
/***/
	ARM_StringVector pricedColumns;
	int prodSize	= MIN(NbProductsToPrice,stdSize);
	int size		= prodSize;

	if (itsComputeProba)
		size += itsNbComputedProba;

	bool isFrontierCalib = (itsCalibrationType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER);
	if(isFrontierCalib)
		++size;

	pricedColumns.reserve(size);

    for (i = 0; i < prodSize; ++i)
	{
		if (itsProductsToPrice[i])
			pricedColumns.push_back( LocalCSOColNamesTable[ LocalCSOProductToPriceColums[i] ] );
	}

	if (itsComputeProba)
	{
		for (i = prodSize; i < prodSize+itsNbComputedProba; ++i)
		{
			pricedColumns.push_back( LocalCSOProbaColNamesTable[4+2*(i-stdSize) ] );
		}
	}

	if(isFrontierCalib)
		pricedColumns.push_back( LocalCSOColNamesTable[Frontier] );

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();
	
	/// Create Cst manager
	ARM_CstManagerPtr cstMgr = CreateCstManager();
	/// Create the Generic Security paid in coupon currency
	CreateAndSetDealDescriptionAndTimeIt(GetKeys()[GetCpnModelKey()], pricedColumns, cstMgr);
	
	/// Create the HW pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping
    CreateAndSetCalibrationAndTimeIt();

	//ComputeAll();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: copy constructor 
///	Returns: void
///	Action : copies a ARM_LocalCSOCalculator
/////////////////////////////////////////////////////////////////
ARM_LocalCSOCalculator::ARM_LocalCSOCalculator (const ARM_LocalCSOCalculator& rhs)
:	ARM_GenCSOCalculator  (rhs),
	itsCalibrationType	(rhs.itsCalibrationType),
	itsCalibStrikeType	(rhs.itsCalibStrikeType),
	itsVnsPricingMethod (rhs.itsVnsPricingMethod),
	itsMoyenessLevel	(rhs.itsMoyenessLevel),
	itsModelType		(rhs.itsModelType),
	itsvFundSpread		(rhs.itsvFundSpread),
	itsvInitialFundSpread (rhs.itsvInitialFundSpread),
	itsvFundNominal		(rhs.itsvFundNominal),
	itsFundSize			(rhs.itsFundSize),
	itsvCpnNominal		(rhs.itsvCpnNominal),
	itsvCpnFloor		(rhs.itsvCpnFloor),
	itsvCpnCap			(rhs.itsvCpnCap),
	itsvCpnLeverageShort(rhs.itsvCpnLeverageShort),
    itsvCpnLeverageLong	(rhs.itsvCpnLeverageLong),
	itsvCpnStrike		(rhs.itsvCpnStrike),
	itsvCpnIsCapped		(rhs.itsvCpnIsCapped),
	itsvCpnIsFloored	(rhs.itsvCpnIsFloored),	
	itsCpnSize			(rhs.itsCpnSize),
	itsvExerFees		(rhs.itsvExerFees),
	itsvIsExerDate		(rhs.itsvIsExerDate),
	itsvFundIndex		(rhs.itsvFundIndex),
	itsvCpnIndex		(rhs.itsvCpnIndex),
	itsExerSize			(rhs.itsExerSize),
	itsFundingPrice		(rhs.itsFundingPrice),
	itsFixPrice			(rhs.itsFixPrice),
	itsFloorPrice		(rhs.itsFloorPrice),
	itsCapPrice			(rhs.itsCapPrice),
	itsSpreadOptionFloorPF	(CreateClonedPtr(&*rhs.itsSpreadOptionFloorPF)),
	itsSpreadOptionCapPF	(CreateClonedPtr(&*rhs.itsSpreadOptionCapPF)),
	itsSpreadOptionFloorPFforSwitch	(CreateClonedPtr(&*rhs.itsSpreadOptionFloorPFforSwitch)),
	itsSpreadOptionCapPFforSwitch	(CreateClonedPtr(&*rhs.itsSpreadOptionCapPFforSwitch)),
	itsCpnLongSensi     (rhs.itsCpnLongSensi),
	itsCpnShortSensi    (rhs.itsCpnShortSensi),
	itsCpnLongSensiForSwitch     (rhs.itsCpnLongSensiForSwitch),
	itsCpnShortSensiForSwitch    (rhs.itsCpnShortSensiForSwitch),
	itsCpnValue		    (rhs.itsCpnValue),
	itsCpnValueForSwitch	    (rhs.itsCpnValueForSwitch),
	itsSoPricesComputed (rhs.itsSoPricesComputed),
	itsComputeProba		(rhs.itsComputeProba),
	itsFirstComputedProba	(rhs.itsFirstComputedProba),
	itsNbComputedProba	(rhs.itsNbComputedProba),
	itsProba			(rhs.itsProba),
	itsSwitchFlag		(rhs.itsSwitchFlag),
	itsCMSFunding		(rhs.itsCMSFunding)
{
}

////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: operator=
///	Returns: 
///	Action : Assignment operator
////////////////////////////////////////////////////
ARM_LocalCSOCalculator& ARM_LocalCSOCalculator::operator = (const ARM_LocalCSOCalculator& rhs)
{	
	if (&rhs != this)
	{ 
		this->~ARM_LocalCSOCalculator();
		new (this) ARM_LocalCSOCalculator (rhs);
	}
	return *this;
}

////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: GetSwaptionPF
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::GetSwaptionPF() const
{ 

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Calib method is null");
#endif
    
	if (itsCalibrationType==DIAG_SPREAD_LONG || itsCalibrationType==DIAG_SPREAD_SHORT || itsCalibrationType==DIAG_LONG_CALIBRATION)
	{
		return GetCalibMethod()->GetPortfolio();
	}
	else if (itsCalibrationType==SHORT_LONG_SPREAD || itsCalibrationType==LONG_CALIBRATION)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no OSW portfolio available");
	}
	else
		return GetCalibMethod()->GetPortfolio();
}

////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: GetSOPortfolio()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::GetSOPortfolio() const
{ 

	if (itsCalibrationType==DIAG_SPREAD_LONG || itsCalibrationType==DIAG_SPREAD_SHORT)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibrationType==SHORT_LONG_SPREAD)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no SO portfolio available");

    
	return GetCalibMethod()->GetPortfolio();
}

////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: GetCMSLONGPortolio()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::GetCMSLONGPortfolio() const
{ 

	if (itsCalibrationType==DIAG_SPREAD_LONG)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibrationType==SHORT_LONG_SPREAD || itsCalibrationType==DIAG_LONG_CALIBRATION)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibrationType==LONG_CALIBRATION)
	{
		return GetCalibMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no CMSLONG portfolio available");

    
	return GetCalibMethod()->GetPortfolio();
}

////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: GetCMSLONGPortolio()
///	Returns: 
///	Action : 
////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::GetCMSSHORTPortfolio() const
{ 

	if (itsCalibrationType==DIAG_SPREAD_SHORT)
	{
		return GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio();
	}
	else if (itsCalibrationType==SHORT_LONG_SPREAD)
	{
		return GetCalibMethod()->GetPortfolio();
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator:: no CMSLONG portfolio available");

    
	return GetCalibMethod()->GetPortfolio();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: SetOSWPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of swaptions
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::SetOSWPortfolio(const ARM_StdPortfolio& port)
{
    if(!(port.size()))
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Portfolio to set is empty");

     if(port.GetAsset(0)->GetName() != ARM_SWAPTION || !(port.IsSameAssetsName()))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Only a portfolio of swaptions can be valid, please advise ");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Calib method is null");
#endif
	
	GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeProductVectorsFromCurves()
{
	size_t i;
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector* fundResetDates = GetFundDateStrip()->GetResetDates() ;
	itsFundSize = fundResetDates->size();

	itsvFundSpread.resize(itsFundSize);
	itsvInitialFundSpread.resize(itsFundSize);
	itsvFundLeverage.resize(itsFundSize);
	itsvInitialFundNominal.resize(itsFundSize);
	itsvFundNominal.resize(itsFundSize);

    for (i=0; i<itsFundSize; i++)
	{
		double lag = (*fundResetDates)[i]-asOfDate;
		itsvFundSpread[i]			= GetFundSpread().Interpolate(lag);
		itsvInitialFundSpread[i]	= GetFundSpread().Interpolate(lag);
		itsvFundLeverage[i]			= GetFundLeverage().Interpolate(lag);
		itsvInitialFundNominal[i]	= GetFundNominal().Interpolate(lag);
		itsvFundNominal[i]			= GetCpnNominal().Interpolate(lag);
	}

	ARM_GP_Vector* cpnResetDates = GetStructDateStrip()->GetResetDates() ;
	itsCpnSize = cpnResetDates->size();
	
	itsvCpnNominal.resize(itsCpnSize);
	itsvCpnIsFloored.resize(itsCpnSize);
	itsvCpnIsCapped.resize(itsCpnSize);
	itsvCpnLeverageShort.resize(itsCpnSize);
	itsvCpnLeverageLong.resize(itsCpnSize);
	itsvCpnStrike.resize(itsCpnSize);
	itsvCpnFloor.resize(itsCpnSize);
	itsvCpnCap.resize(itsCpnSize);

	for (i=0; i<itsCpnSize; i++)
	{
		double lag = (*cpnResetDates)[i]-asOfDate;
		itsvCpnNominal[i] = GetCpnNominal().Interpolate(lag);

		if (i<GetNbFixFlows())
		{
			itsvCpnIsFloored[i]		= false;
			itsvCpnIsCapped[i]		= false;
			itsvCpnLeverageShort[i] = 0.0;
			itsvCpnLeverageLong[i]	= 0.0;
			double strike = GetFixCpn().Interpolate(lag);
			itsvCpnStrike[i]		= -strike;
			itsvCpnFloor[i]			= strike;
			itsvCpnCap[i]			= strike;
		}
		else
		{	double floor = GetCpnFloor().Interpolate(lag) ;
			double cap   = GetCpnCap  ().Interpolate(lag) ;
			
			itsvCpnFloor[i] = floor;
			itsvCpnIsFloored[i] = (floor<NO_FLOOR) ? false : true;

			itsvCpnCap[i] = cap;
			itsvCpnIsCapped[i] = (cap>NO_CAP) ? false : true;

			itsvCpnLeverageShort[i]	= GetLeverageShort().Interpolate(lag);
			itsvCpnLeverageLong[i]	= GetLeverageLong().Interpolate(lag);
			itsvCpnStrike[i]		= GetStrike().Interpolate(lag);
		}
	}
	if(IsBasis())
		itsvFundSpread = ComputeDomesticBasis();

	ARM_GP_Vector* exerciseDates = GetExerciseDateStrip()->GetResetDates() ;
	itsExerSize = exerciseDates->size();

	/// past  no call
	size_t nbPastNoCall=0;
	while(nbPastNoCall < itsExerSize && (*exerciseDates)[nbPastNoCall] < asOfDate)
        ++nbPastNoCall;

	itsvExerFees.resize(itsExerSize);
	itsvIsExerDate.resize(itsExerSize);
	itsvFundIndex.resize(itsExerSize+1);
	itsvCpnIndex.resize(itsExerSize+1);

	for (i=0; i<itsExerSize; i++)
	{
		double lag = (*exerciseDates)[i] - asOfDate;
		double fee =  (i < nbPastNoCall) ? NON_CALL_FEE : GetFees().Interpolate(lag);
		itsvExerFees[i] = fee;

		itsvIsExerDate[i] = ( fabs((fee-NON_CALL_FEE))<K_NEW_DOUBLE_TOL )? false : true;
	}

	/// last: compute indexes to relate exer - funding - cpn
		
	ARM_GP_Vector* cpnStartDates  = GetStructDateStrip()->GetFlowStartDates() ;
	ARM_GP_Vector* fundStartDates = GetFundDateStrip()->GetFlowStartDates() ;
	
	size_t k;
	for (k=0, i = 0; k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*fundStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : no call date allowed between fixing dates and start dates");
		itsvFundIndex[k] = i;
	}
	for (k=0, i= 0; k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*cpnStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*cpnResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : no call date allowed between fixing dates and start dates");
		itsvCpnIndex[k] = i;
	}
	itsvFundIndex[itsExerSize] = fundStartDates->size();
	itsvCpnIndex[itsExerSize] = cpnStartDates->size();
		
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<itsExerSize; k++)
	{
		if (itsvFundIndex[k]>=itsvFundIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Funding) is forbidden");

		if (itsvCpnIndex[k]>=itsvCpnIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Coupon) is forbidden");
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	for (k=0; k<itsExerSize; k++) 
	{		
		if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*cpnStartDates)[itsvCpnIndex[k]])>0)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: IsFixCoupon
///	Returns: a boolean
///	Action : tells if the period isFix
/////////////////////////////////////////////////////////////////
bool ARM_LocalCSOCalculator::IsFixCoupon (size_t cpnIdx, double& fixRate) const
{
	if (cpnIdx<0 || cpnIdx>= itsCpnSize)
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::IsFixCoupon : invalid index");

	if (cpnIdx<GetNbFixFlows())
	{		
		fixRate = -itsvCpnStrike[cpnIdx]; // as computed in method ComputeProductVectorsFromCurves
		return true;
	}

	bool leverageShortIsNull =	(fabs(itsvCpnLeverageShort[cpnIdx])<K_NEW_DOUBLE_TOL);
	bool leverageLongIsNull  =	(fabs(itsvCpnLeverageLong[cpnIdx]) <K_NEW_DOUBLE_TOL);

	if (leverageShortIsNull && leverageLongIsNull)
	{
		fixRate = -itsvCpnStrike[cpnIdx] ;
		if (fixRate<itsvCpnFloor[cpnIdx])
			fixRate = itsvCpnFloor[cpnIdx];
		if (fixRate>itsvCpnCap[cpnIdx])
			fixRate = itsvCpnCap[cpnIdx];

		return true;
	}

	return false;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateCstManager
///	Returns: ARM_CstManagerPtr
///	Action : creates the cst manager
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_LocalCSOCalculator::CreateCstManager()
{
	vector<string> names(0);
	vector<ARM_GramFctorArg> values(0);

	ARM_GP_Vector cpnStrikePlusFloor = itsvCpnFloor + itsvCpnStrike;
	ARM_GP_Vector cpnStrikePlusCap   = itsvCpnCap + itsvCpnStrike;
	ARM_GP_Vector cpnFloorTimesNotio = itsvCpnFloor * itsvCpnNominal;
    ARM_GP_Vector spreadTimesNotio	 = itsvFundSpread * itsvFundNominal;
	ARM_GP_Vector leverageTimesNotio = itsvFundLeverage * itsvFundNominal;

	for (size_t i(0); i<itsExerSize; i++)
	{		
		/// define index
		CC_Ostringstream indexStr;
		/// sort by size
		i < 10 ? indexStr << 0 << i :  indexStr << i;
	
		/// names
		names.push_back("FundNominal"		+ indexStr.str()); /// first name
		names.push_back("LeverageTimesNotio"+ indexStr.str()); /// ....
        names.push_back("CpnNominal"		+ indexStr.str()); /// ....
		names.push_back("CpnLeverageLong"	+ indexStr.str());
		names.push_back("CpnLeverageShort"	+ indexStr.str());
		names.push_back("CpnFloorTimesNotio"+ indexStr.str());
		names.push_back("CpnStrikePlusFloor"+ indexStr.str());
		names.push_back("CpnStrikePlusCap"	+ indexStr.str());

		/// values
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(itsvFundNominal.begin()		+ itsvFundIndex[i], itsvFundNominal.end()))) );/// first value
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(leverageTimesNotio.begin()	+ itsvFundIndex[i],	leverageTimesNotio.end()))) );   //// ...
        values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(itsvCpnNominal.begin()		+ itsvCpnIndex[i],	itsvCpnNominal.end()))) );   //// ...
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(itsvCpnLeverageLong.begin()	+ itsvCpnIndex[i],	itsvCpnLeverageLong.end()))) );
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(itsvCpnLeverageShort.begin()	+ itsvCpnIndex[i],	itsvCpnLeverageShort.end()))) );
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnFloorTimesNotio.begin()	+ itsvCpnIndex[i],	cpnFloorTimesNotio.end()))) );
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnStrikePlusFloor.begin()	+ itsvCpnIndex[i],	cpnStrikePlusFloor.end()))) );
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnStrikePlusCap.begin()		+ itsvCpnIndex[i],	cpnStrikePlusCap.end()))) );
	}
	
	names.push_back("SpreadTimesNotio"); /// ....
	values.push_back( ARM_GramFctorArg(ARM_VectorPtr((ARM_GP_Vector*)spreadTimesNotio.Clone() ) ));   //// ...

	if (IsSwitchable())
	{
		int size = itsvFundSpread.size();
		ARM_GP_Vector auxVector(size,1000.);
		ARM_GP_Vector auxNull(size,0.);
		ARM_GP_Vector cpnStrikePlusFloorFund = auxNull-itsvFundSpread-auxVector;
		ARM_GP_Vector cpnStrikePlusCapFund   = auxNull-itsvFundSpread+auxVector;
		ARM_GP_Vector cpnFloorFundTimesNotio = itsvCpnNominal*cpnStrikePlusFloorFund;
		cpnStrikePlusFloorFund -= itsvFundSpread;
		cpnStrikePlusCapFund -= itsvFundSpread;

		for (size_t i(0); i<itsExerSize; i++)
		{		
			/// define index
			CC_Ostringstream indexStr;
			/// sort by size
			i < 10 ? indexStr << 0 << i :  indexStr << i;
		
			/// names
			names.push_back("CpnFloorFundTimesNotio"+ indexStr.str());
			names.push_back("CpnStrikePlusFloorFund"+ indexStr.str());
			names.push_back("CpnStrikePlusCapFund"	+ indexStr.str());

			/// values
			values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnFloorFundTimesNotio.begin()	+ itsvCpnIndex[i],	cpnFloorFundTimesNotio.end()))) );
			values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnStrikePlusFloorFund.begin()	+ itsvCpnIndex[i],	cpnStrikePlusFloorFund.end()))) );
			values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnStrikePlusCapFund.begin()		+ itsvCpnIndex[i],	cpnStrikePlusCapFund.end()))) );
		}
	}
	
	return  ARM_CstManagerPtr(new ARM_CstManager(names, values));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateSwaptionPortfolio
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::CreateSwaptionPortfolio(CalibrationType type, CalibStrikeType calibStrikeType,bool isUpdateStrike)
{	
	ARM_GP_Vector* resetDates = GetStructDateStrip()->GetResetDates();
	ARM_GP_Vector* startDates = GetStructDateStrip()->GetFlowStartDates();

	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	int i;

	list< ARM_Security* > swaptionList;
	list< ARM_Security* > fwdStartSwaptionList;

	bool USE_NORMAL_SENSIS = true;

	/// this needs to be computed...
	if (type == BASKET_CALIBRATION)
	{
		ComputeSOPortfolioPrices();

		if (USE_NORMAL_SENSIS)
		{
			ComputeSONormalSensitivities();
		}
		else
		{
			ComputeSOSensitivities();
		}
	}
	int k;
	for (i = 0,k=0; i < itsExerSize; ++i)
	{		
		if (itsvIsExerDate[i])
		{
			if (type == DIAG_CALIBRATION)
			{
				ARM_Date startDate ((*startDates)[itsvCpnIndex[i]]);
				ARM_Date expiryDate((*resetDates)[itsvCpnIndex[i]]);

				ARM_Swap stdSwap(	startDate,
									GetEndDate(),
									GetIndexType(),0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

				stdSwap.SetModel(CFBSModel);
        
				/// ATM
				/// if (itsCalibStrikeType == EQUIVALENT)
					/// ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator: EQUIVALENT strike not implemented for DIAG calibration");
					
				double equivStrike = stdSwap.PriceToRate((ARM_Date) asOfDate, 0.0);

				int RecOrPay = K_RCV;
				ARM_Swaption swaption((&stdSwap),RecOrPay,K_EUROPEAN,equivStrike,expiryDate);
				
				swaptionList.push_back(static_cast<ARM_Swaption*>(swaption.Clone()));
			}
			else if (type == BASKET_CALIBRATION)
			{
				ARM_Swaption* fwdStartSwaption = NULL;
				ARM_Swaption* swaption = NULL;

				swaption = CreateVarNotionalSwaptionAtExer(i, fwdStartSwaption, calibStrikeType,isUpdateStrike);

				swaptionList.push_back(swaption);

				if (fwdStartSwaption)
					fwdStartSwaptionList.push_back(fwdStartSwaption);

			}
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator: invalid calibration type");
		}
		else
			k++;
	}

	ARM_StdPortfolioPtr swaptionPF(new ARM_StdPortfolio(swaptionList));
	itsFwdStartSwaptionPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(fwdStartSwaptionList));
    	
	for(i=0;i<swaptionPF->size();++i)
	{
		swaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        swaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		swaptionPF->SetPrice(DEFAULT_PRICE,i);
	}

	for(i=0;i<itsFwdStartSwaptionPF->size();++i)
	{
		itsFwdStartSwaptionPF->SetPrecision(DEFAULT_PRECISION,i);
        itsFwdStartSwaptionPF->SetWeight(DEFAULT_WEIGHT,i);
		itsFwdStartSwaptionPF->SetPrice(DEFAULT_PRICE,i);
	}

	return swaptionPF;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateSOPortfolio
///	Returns: void
///	Action : create ARM_SpreadOption PF's 
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::CreateSOPortfolio()
{	
	///
	ARM_GP_Vector* startDates = GetStructDateStrip()->GetFlowStartDates();
	ARM_GP_Vector* resetDates = GetStructDateStrip()->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_ZeroCurve* zc = dynamic_cast< ARM_ZeroCurve* >( GetMktDataManager()->GetData(GetKeys()[YcKey]) );
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );

	itsVanillaArgSOVect.resize(0);
	list< ARM_Security* > spreadOptionFloorList;
	list< ARM_Security* > spreadOptionCapList;
		
	ARM_Vector* fixing1 = NULL;
	ARM_Vector* fixing2 = NULL;

	double notUsed;
	
	for (int i = 0; i < itsCpnSize; i++)
	{			
		ARM_Date startDate((*startDates)[i]);
		ARM_Date endDate((*startDates)[i]);
		endDate.AddMonths(12/GetCpnFreq());
	
		ARM_ReferenceValue vWeight1(itsvCpnLeverageShort[i]);
		ARM_ReferenceValue vWeight2(itsvCpnLeverageLong[i]);
		ARM_ReferenceValue aStrikeFloor( (itsvCpnFloor[i]+itsvCpnStrike[i])*100. );
		ARM_ReferenceValue aStrikeCap( (itsvCpnCap[i]+itsvCpnStrike[i])*100. );
				
		if (!IsFixCoupon(i, notUsed) && ((*resetDates)[i] - asOfDate > ARM_NumericConstants::ARM_DOUBLE_TOLERENCE))
		{	
			///-----------------------
			/// Floorlet
			///-----------------------
			ARM_SpreadOption soFloor ( startDate,
									   endDate,
									   K_CAP,
									   &aStrikeFloor,
									   (ARM_INDEX_TYPE)GetCMSShort(),
									   (ARM_INDEX_TYPE)GetCMSLong(),
									   &vWeight1,
									   &vWeight2, 
									   GetCpnDaycount(),
									   GetCpnFreq(),
									   GetCpnFreq(),
									   GetCpnResetTiming(),
									   K_ARREARS,
									   GetCurrencyUnit(),
									   1.0,
									   fixing1,
									   fixing2);

			spreadOptionFloorList.push_back(static_cast<ARM_SpreadOption*>(soFloor.Clone()));

		
			ARM_VanillaSpreadOptionArg* vanillaSO = CreateVanillaArgSO(&soFloor);
			itsVanillaArgSOVect.push_back(ARM_VanillaArgPtr(vanillaSO));
					

			///-----------------------
			/// Caplet
			///-----------------------
			ARM_SpreadOption soCap(startDate,
								   endDate,
								   K_CAP,
								   &aStrikeCap,
								   (ARM_INDEX_TYPE)GetCMSShort(),
								   (ARM_INDEX_TYPE)GetCMSLong(),
								   &vWeight1,
								   &vWeight2, 
								   GetCpnDaycount(),
								   GetCpnFreq(),
								   GetCpnFreq(),
								   GetCpnResetTiming(),
								   K_ARREARS,
								   GetCurrencyUnit(),
								   1.0,
								   fixing1,
								   fixing2);

			spreadOptionCapList.push_back(static_cast<ARM_SpreadOption*>(soCap.Clone()));
		}
		else
		{
			// add vanilla arg anyway (easier pour la suite)
			itsVanillaArgSOVect.push_back(ARM_VanillaArgPtr(NULL));
		}
	}

	itsSpreadOptionFloorPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionFloorList));
	itsSpreadOptionCapPF   = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionCapList));

	if (IsSwitchable())
	{
		itsVanillaArgSOVectForSwitch.resize(0);
		list< ARM_Security* > spreadOptionFloorListForSwitch;
		list< ARM_Security* > spreadOptionCapListForSwitch;

		for (int i = 0; i < itsCpnSize; i++)
		{			
			ARM_Date startDate((*startDates)[i]);
			ARM_Date endDate((*startDates)[i]);
			endDate.AddMonths(12/GetCpnFreq());
		
			ARM_ReferenceValue vWeight1(0.);
			ARM_ReferenceValue vWeight2(1.);
			ARM_ReferenceValue aStrikeFloor( (-1000-2.*itsvFundSpread[i])*100 );
			ARM_ReferenceValue aStrikeCap( (1000-2.*itsvFundSpread[i])*100 );
					
			if ( ((*resetDates)[i] - asOfDate > ARM_NumericConstants::ARM_DOUBLE_TOLERENCE))
			{	
				ARM_SpreadOption soFloor ( startDate,
										   endDate,
										   K_CAP,
										   &aStrikeFloor,
										   (ARM_INDEX_TYPE)GetCMSShort(),
										   (ARM_INDEX_TYPE)GetCMSFunding(),
										   &vWeight1,
										   &vWeight2, 
										   GetFundDaycount(),
										   GetCpnFreq(),
										   GetCpnFreq(),
										   GetCpnResetTiming(),
										   K_ARREARS,
										   GetCurrencyUnit(),
										   1.0,
										   fixing1,
										   fixing2);

				spreadOptionFloorListForSwitch.push_back(static_cast<ARM_SpreadOption*>(soFloor.Clone()));

				ARM_VanillaSpreadOptionArg* vanillaSO = CreateVanillaArgSO(&soFloor);
				itsVanillaArgSOVectForSwitch.push_back(ARM_VanillaArgPtr(vanillaSO));

				ARM_SpreadOption soCap(startDate,
									   endDate,
									   K_CAP,
									   &aStrikeCap,
									   (ARM_INDEX_TYPE)GetCMSShort(),
									   (ARM_INDEX_TYPE)GetCMSFunding(),
									   &vWeight1,
									   &vWeight2, 
									   GetFundDaycount(),
									   GetCpnFreq(),
									   GetCpnFreq(),
									   GetCpnResetTiming(),
									   K_ARREARS,
									   GetCurrencyUnit(),
									   1.0,
									   fixing1,
									   fixing2);

				spreadOptionCapListForSwitch.push_back(static_cast<ARM_SpreadOption*>(soCap.Clone()));
			}
			else
			{
				itsVanillaArgSOVectForSwitch.push_back(ARM_VanillaArgPtr(NULL));
			}
		}

		itsSpreadOptionFloorPFforSwitch = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionFloorListForSwitch));
		itsSpreadOptionCapPFforSwitch   = ARM_StdPortfolioPtr(new ARM_StdPortfolio(spreadOptionCapListForSwitch));
	}
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::CreateEmptyCalibration()
{		
	if ( itsCalibrationType == BASKET_CALIBRATION || itsCalibrationType == DIAG_CALIBRATION || itsCalibrationType == DIAG_BASKET_CALIBRATION)
	{
// OLD CALIBRATION
		ARM_StdPortfolioPtr swaptionPF; 
		
		/// create diagonal swaptions or var notio swaptions
		if ( itsCalibrationType == BASKET_CALIBRATION || itsCalibrationType == DIAG_CALIBRATION )
		{
			swaptionPF = CreateSwaptionPortfolio(itsCalibrationType,itsCalibStrikeType[0]);
			if(itsCalibStrikeType[0] == FRONTIER)
			{
				SetDiagStrikeToFrontier(swaptionPF);
			}
		}
		/// create diagonal swaptions
		else if (itsCalibrationType == DIAG_BASKET_CALIBRATION)
		{
			swaptionPF = CreateSwaptionPortfolio(DIAG_CALIBRATION,itsCalibStrikeType[0]);
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator::CreateEmptyCalibration: invalid calibration type");
			
		/// create vol param & vol calib method
		ARM_GP_Vector calibTimes(1, 0.00);
		ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
		ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
		ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
		ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
		ARM_ModelParamVector paramVec (0);
		paramVec.push_back(volParam);
		ARM_CalibMethodPtr volCalibMethod (new ARM_CalibMethod(swaptionPF, paramVec, ARM_CalibMethodType::Bootstrap1D) );
		delete volParam; 

		if ( itsCalibrationType == BASKET_CALIBRATION || itsCalibrationType == DIAG_CALIBRATION )
		{
			SetCalibMethod(volCalibMethod);
		}
		else if (itsCalibrationType == DIAG_BASKET_CALIBRATION)
		{		
			ARM_StdPortfolioPtr varNotioPF = CreateSwaptionPortfolio(BASKET_CALIBRATION,itsCalibStrikeType[0]);

			ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));

			if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
				ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid MeanReversion parameter");

			ARM_ModelParamVector mrsParamVect(1);
			mrsParamVect[0]= mrsParam ;

			/// fake prices to avoid exception in ARM_CalibMethod::Validate()
			for(int i=0;i<varNotioPF->GetSize();++i)
			{			
				varNotioPF->SetPrice (i, i);
			}
		

			/// MRS optimisation with an embedded volatility bootstrapping
			ARM_CalibMethod* mrsCalib = new ARM_CalibMethod (varNotioPF,
															mrsParamVect,
															ARM_CalibMethodType::Optimize,
															50, /// max nb iter
															ARM_CalibrationTarget::PriceTarget,
															&*volCalibMethod);
			SetCalibMethod(ARM_CalibMethodPtr(mrsCalib));
		}
	}
	else
	{
		CreateEmptyCalibrationFor2F();
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: SetDiagStrikeToFrontier
///	Returns: void
///	Action : compute exercise strikes and set them in the diagonal
///			 swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::SetDiagStrikeToFrontier(ARM_StdPortfolioPtr& swaptionPF)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

	/// Set ATM strike
	size_t i,nbSwaptions = swaptionPF->size();
	double atmStrike,price,vega;
	ARM_Swaption *swaption;
	ARM_VectorPtr lastFrontier(new ARM_GP_Vector(nbSwaptions));
	ARM_VectorPtr atmStrikes(new ARM_GP_Vector(nbSwaptions));
	for(i=0;i<nbSwaptions;++i)
	{
		swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
		swaption->SetModel(oswBSModel);
		atmStrike = swaption->PriceToRate((ARM_Date) asOfDate, 0.0);
		swaption->UpdateStrike(atmStrike);
		price = swaption->ComputePrice();
		swaptionPF->SetPrice(price,i);
		(*atmStrikes)[i]=atmStrike;
		(*lastFrontier)[i]=atmStrike;
	}

	/// Calibrate ATM diagonal swaption
	ARM_CalibMethodPtr saveCalibMethod = GetCalibMethod();

	ARM_GP_Vector calibTimes(1, 0.00);
	ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
	ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
	ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
	ARM_CurveModelParam volParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
	ARM_ModelParamVector volParamVect(1,&volParam);
	ARM_ModelFitterDes modelfitter(ARM_ModelFitterSolverType::Brent);
	ARM_CalibMethodPtr AtmVolCalib(new ARM_CalibMethod(swaptionPF,volParamVect,ARM_CalibMethodType::Bootstrap1D,&modelfitter));
	SetCalibMethod(AtmVolCalib);
	Calibrate(); /// calibrate model vol and local models !

	/// Compute exercise frontier projection w.r.t. swap rate
    ARM_GenSecurityPtr genSec = GetGenSecurity();
	bool savePayoffsFlag = genSec->GetOtherPayoffsFlag();
	genSec->SetOtherPayoffsFlag(true);
	ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );
	ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer);
	genPricer->Price();
	ARM_VectorPtr frontier = genPricer->GetPricerInfo()->GetContents(LocalCSOColNamesTable[Frontier]).GetData("Intermediateprices").GetVector();
	vector<ARM_GP_Vector> frontiers;

	/// Set diagonal swaption strike to frontier swap rates
	double frontierStrike,strike,strikeStep;
	size_t j,nbStrikeSteps=20;
	size_t k,bestStrikeIdx=0,nbFrontierLoop=3;
	double strikeChange,bestStrikeChange=1000.,strikeMinChange=0.05; // 5bp
	for(k=0;k<nbFrontierLoop;++k)
	{
		strikeChange=0;
		for(i=0;i<nbSwaptions;++i)
		{
			(*frontier)[i] *= 100.0;
			strikeChange=CC_Max(fabs((*lastFrontier)[i]-(*frontier)[i]),strikeChange);

			swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			swaption->UpdateStrike((*frontier)[i]);
			price = swaption->ComputePrice();
			swaptionPF->SetPrice(price,i);
		}

		frontiers.push_back(*frontier);

		if(strikeChange<bestStrikeChange)
		{
			bestStrikeChange = strikeChange;
			bestStrikeIdx = k;
		}

		if(strikeChange<strikeMinChange || k+1==nbFrontierLoop)
		{
			for(i=0;i<nbSwaptions;++i)
			{
				swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
				frontierStrike = frontiers[bestStrikeIdx][i];
				vega = swaption->ComputeSensitivity(K_VEGA)/100.;
				if(vega<VEGA_MIN_TO_SELECT)
				{
					/// Find a more convenient strike !
					atmStrike = (*atmStrikes)[i];
					strikeStep = (frontierStrike-atmStrike)/nbStrikeSteps;
					for(j=0,strike=frontierStrike;j<nbStrikeSteps;++j,strike-=strikeStep)
					{
						swaption->UpdateStrike(strike);
						price = swaption->ComputePrice();
						vega = swaption->ComputeSensitivity(K_VEGA)/100.;
						if(vega >= VEGA_MIN_TO_SELECT)
							break;
					}
					if(j>=nbStrikeSteps)
					{
						/// No chance back to ATM strike !
						swaption->UpdateStrike(atmStrike);
						price = swaption->ComputePrice();
					}
					swaptionPF->SetPrice(price,i);
				}
			}
			break;
		}

		/// New Frontier computation
		GetCalibMethod()->SetPortfolio(swaptionPF);
		Calibrate();
		genPricer->Price();
		lastFrontier = frontier;
		frontier = genPricer->GetPricerInfo()->GetContents(LocalCSOColNamesTable[Frontier]).GetData("Intermediateprices").GetVector();
	}

	/// Restore previous values
	genSec->SetOtherPayoffsFlag(savePayoffsFlag);
	SetCalibMethod(saveCalibMethod);

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateEmptyCalibrationFor2F
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::CreateEmptyCalibrationFor2F()
{
	ARM_GP_Vector calibTimes(1, 0.00);
	ARM_GP_Vector volatility(1, SIGMA_DEFAULT_VALUE);
	ARM_GP_Vector lowerVolatility(1, SIGMA_LOWER_BOUND);
	ARM_GP_Vector upperVolatility(1, SIGMA_UPPER_BOUND);
	ARM_CurveModelParam volParam(ARM_ModelParamType::Volatility, &volatility, &calibTimes ,"","STEPUPRIGHT",&lowerVolatility,&upperVolatility);
	ARM_ModelParamVector volParamVec(1,&volParam);

	ARM_CalibMethodPtr calibMethod(NULL);

	size_t N;

	ARM_StdPortfolioPtr pf1; 
	ARM_StdPortfolioPtr pf2; 
	ARM_StdPortfolioPtr pf3; 
	if (itsCalibrationType==DIAG_SPREAD_LONG)
	{
		pf1 = CreateSwaptionPortfolioFor2F(itsCalibStrikeType[0]);
		pf2 = CreateSOPortfolioFor2F(true);
		pf3 = CreateSOPortfolioFor2F(false,true);
		N   = 3;
	}
	else if (itsCalibrationType==DIAG_SPREAD_SHORT)
	{
		pf1 = CreateSwaptionPortfolioFor2F(itsCalibStrikeType[0]);
		pf2 = CreateSOPortfolioFor2F(true);
		pf3 = CreateSOPortfolioFor2F(false,false);
		N   = 3;
	}
	else if (itsCalibrationType==SHORT_LONG_SPREAD)
	{
		pf1 = CreateSOPortfolioFor2F(false,false);
		pf2 = CreateSOPortfolioFor2F(false,true);
		pf3 = CreateSOPortfolioFor2F(true);
		N   = 3;
	}
	else if (itsCalibrationType==LONG_CALIBRATION)
	{
		pf1 = CreateSOPortfolioFor2F(false,true);
		N   = 1;
	}
	else if (itsCalibrationType==DIAG_LONG_CALIBRATION)
	{
		pf1 = CreateSwaptionPortfolioFor2F(itsCalibStrikeType[0]);
		pf2 = CreateSOPortfolioFor2F(false,true);
		N   = 2;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator::CreateEmptyCalibration: invalid calibration type");

	if (N==3)
	{
		ARM_CalibMethod calibVolMethod(pf3,
							volParamVec,ARM_CalibMethodType::HW2FOnly,N,ARM_CalibrationTarget::UnknownTarget,
							NULL,NULL,false,0,1,false);

		calibTimes.push_back(1.0); // to switch H&W2F params to time dependent version
		ARM_GP_Vector volatilityRatio(2, VOLRATIO_DEFAULT_VALUE);
		ARM_CurveModelParam volRatioParam(ARM_ModelParamType::VolatilityRatio, &volatilityRatio, &calibTimes ,"","STEPUPRIGHT");
		ARM_ModelParamVector volRatioParamVec(1,&volRatioParam);
		ARM_CalibMethod calibVolRatioMethod(pf2,
							volRatioParamVec,ARM_CalibMethodType::HW2FOnly,N,ARM_CalibrationTarget::UnknownTarget,
							&calibVolMethod,NULL,false,0,1,false);

		ARM_GP_Vector correl(2, CORREL_DEFAULT_VALUE);
		ARM_CurveModelParam correlParam(ARM_ModelParamType::Correlation, &correl, &calibTimes ,"","STEPUPRIGHT");
		ARM_ModelParamVector correlParamVec(1,&correlParam);
		calibMethod = ARM_CalibMethodPtr ( new ARM_CalibMethod (pf1,
							correlParamVec,ARM_CalibMethodType::HW2FOnly,N,ARM_CalibrationTarget::UnknownTarget,
							&calibVolRatioMethod,NULL,false,0,1,false) );

		SetCalibMethod(calibMethod);
	}
	else if (N==2)
	{
		ARM_CalibMethod calibVolMethod(pf2,
							volParamVec,ARM_CalibMethodType::HW2FOnly,N,ARM_CalibrationTarget::UnknownTarget,
							NULL,NULL,false,0,1,false);

		calibTimes.push_back(1.0); // to switch H&W2F params to time dependent version
		ARM_GP_Vector volatilityRatio(2, VOLRATIO_DEFAULT_VALUE);
		ARM_CurveModelParam volRatioParam(ARM_ModelParamType::VolatilityRatio, &volatilityRatio, &calibTimes ,"","STEPUPRIGHT");
		ARM_ModelParamVector volRatioParamVec(1,&volRatioParam);
		calibMethod = ARM_CalibMethodPtr(new ARM_CalibMethod (pf1,
							volRatioParamVec,ARM_CalibMethodType::HW2FOnly,N,ARM_CalibrationTarget::UnknownTarget,
							&calibVolMethod,NULL,false,0,1,false));
		
		SetCalibMethod(calibMethod);
	}
	else if (N==1)
	{
		calibMethod = ARM_CalibMethodPtr(new ARM_CalibMethod (pf1,volParamVec,ARM_CalibMethodType::Bootstrap1D));
		SetCalibMethod(calibMethod);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateSOPortfolioFor2F
///	Returns: a portfolio
///	Action : create the list of CMS spread options. Only one
///			 option is created for each call date
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::CreateSOPortfolioFor2F(bool isSOCond,bool isLong)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_GP_Vector* startDates = GetExerciseDateStripUnadj()->GetFlowStartDates();
	ARM_GP_Vector* endDates = GetExerciseDateStripUnadj()->GetFlowEndDates();
	ARM_GP_Vector* cpnResetDates = GetStructDateStrip()->GetResetDates() ;

	size_t i,nbFlows = startDates->size();

	list< ARM_Security* > soList;
	ARM_SpreadOption* soSec;

	ARM_ReferenceValue strikes(0.0); // strike will be computed later
	
	ARM_INDEX_TYPE index1 = (ARM_INDEX_TYPE)GetCMSShort();
	ARM_INDEX_TYPE index2 = (ARM_INDEX_TYPE)GetCMSLong();
	ARM_ReferenceValue* unUsedFixing=NULL;
	int defaultResetGap = 10000;
	if(!isSOCond)
	{
		if (isLong)
		{
			index1 = (ARM_INDEX_TYPE)GetCMSLong();
			index2 = (ARM_INDEX_TYPE)GetCMSShort();
		}
		else
		{
			index1 = (ARM_INDEX_TYPE)GetCMSShort();
			index2 = (ARM_INDEX_TYPE)GetCMSLong();
		}
	}
	
	for(i=0;i<nbFlows;++i)
	{
		if (itsvIsExerDate[i])
		{
			double lag = (*cpnResetDates)[i]-asOfDate;

			/// Compute unadjusted flow end date
			ARM_Date startDate((*startDates)[i]);
			ARM_Date endDate((*endDates)[i]);
			
			ARM_ReferenceValue weight1(isLong?GetLeverageLong().Interpolate(lag):GetLeverageShort().Interpolate(lag));
			ARM_ReferenceValue weight2(isLong?GetLeverageShort().Interpolate(lag):GetLeverageLong().Interpolate(lag));
			if (!isSOCond)
				weight2=ARM_ReferenceValue(0.0);

			/// Create the CMS spread option with a single flow at payment frequency
			soSec = new ARM_SpreadOption(startDate, endDate, K_CAP, &strikes,
				index1,index2,&weight1,&weight2,
				GetCpnDaycount(), GetExerciseFreq(), GetExerciseFreq(), K_ADVANCE, K_ARREARS, 
				GetCurrencyUnit(),defaultResetGap,unUsedFixing, unUsedFixing);

			soList.push_back(soSec);
		}
	}

	ARM_StdPortfolioPtr soPF(new ARM_StdPortfolio(soList));
    	
	for(i=0;i<soPF->size();++i)
	{
		soPF->SetPrecision(DEFAULT_PRECISION,i);
        soPF->SetWeight(DEFAULT_WEIGHT,i);
		soPF->SetPrice(DEFAULT_PRICE,i);
	}
	return soPF;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateSwaptionPortfolioFor2F
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_LocalCSOCalculator::CreateSwaptionPortfolioFor2F(CalibStrikeType cst)
{
	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_GP_Vector* startDates = GetExerciseDateStrip()->GetFlowStartDates();
	size_t i,nbFlows = startDates->size();

	list< ARM_Security* > swoptList;

	ARM_ReferenceValue strikes(0.0); // strike will be computed later
	int defaultResetGap = 10000;

	ARM_Date endDate = GetEndDate();

	if (cst == EQUIVALENT)
	{
		ComputeSOPortfolioPrices();
		ComputeSONormalSensitivities();
	}
	for(i=0;i<nbFlows;++i)
	{
		if (itsvIsExerDate[i])
		{
			ARM_Date startDate ((*startDates)[i]);

			ARM_Swap stdSwap(	startDate,
								endDate,
								GetIndexType(),
								0.0,
								1.0,
								GetPayRec(),
								K_DEF_FREQ,
								K_DEF_FREQ,
								GetCurrencyUnit());

			ARM_Date expiryDate((*(stdSwap.GetFloatLeg()->GetResetDates()))[0]);
			stdSwap.SetModel(CFBSModel);

			double residualPv=0.0;
			if(cst == EQUIVALENT)
			{
				double floatLegPv,soLegPv;
				residualPv = ComputeResidualUnderlying(i,floatLegPv,soLegPv);

				ARM_SwapLeg* stdFloatLeg = stdSwap.Get1stLeg();
				if(stdFloatLeg->GetLegType() == K_FIXED_LEG)
					stdFloatLeg = stdSwap.Get2ndLeg();
				stdFloatLeg->SetModel(CFBSModel);

				double stdFloatNotio = floatLegPv / fabs(stdFloatLeg->ComputePrice());
				residualPv = GetPayRec() * (soLegPv-floatLegPv)/stdFloatNotio;
			}
					
			double equivStrike = stdSwap.PriceToRate((ARM_Date) asOfDate, residualPv);
						
			int RecOrPay = K_RCV;
			ARM_Swaption swoptSec((&stdSwap),RecOrPay,K_EUROPEAN,equivStrike,expiryDate);
			swoptList.push_back(static_cast<ARM_Swaption*>(swoptSec.Clone()));
		}
	}

	ARM_StdPortfolioPtr swoptPF(new ARM_StdPortfolio(swoptList));
    	
	for(i=0;i<swoptPF->size();++i)
	{
		swoptPF->SetPrecision(DEFAULT_PRECISION,i);
        swoptPF->SetWeight(DEFAULT_WEIGHT,i);
		swoptPF->SetPrice(DEFAULT_PRICE,i);
	}
	return swoptPF;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeSwaptionPortfolioPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSwaptionPortfolioPrices()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	   	
	ARM_Swaption* swaption;	
	double price;

	///
	/// standard swaptions --> use ARM kernel swaption pricing
	///
	if (itsCalibrationType == DIAG_CALIBRATION)
	{
		ARM_StdPortfolioPtr pf = GetCalibMethod()->GetPortfolio();

		for(int i=0;i<pf->GetSize();++i)
		{
			swaption = static_cast< ARM_Swaption* >(pf->GetAsset(i));
			swaption->SetModel(oswBSModel);
			price = swaption->ComputePrice();
			pf->SetPrice (price, i);
		}
	}
	///
	/// variable notio swaptions --> use GP pricing
	///
	else if (itsCalibrationType == BASKET_CALIBRATION)
	{
		ARM_StdPortfolioPtr pf = GetCalibMethod()->GetPortfolio();

		ARM_StringVector keys(2);
		keys[0] = GetKeys()[YcKey];
		keys[1] = GetKeys()[OswModelKey];
		ARM_MarketData_ManagerRep mktDataManager(asOfDate);
		mktDataManager.RegisterData(GetKeys()[YcKey],cpnCurve);
		mktDataManager.RegisterData(GetKeys()[OswModelKey],oswBSModel);
		ARM_MarketIRModel analyticModel(mktDataManager, keys, itsVnsPricingMethod,itsMoyenessLevel);
		
		for(int i=0;i<pf->GetSize();++i)
		{
			ARM_VanillaArg* swaption = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(pf->GetAsset(i), asOfDate, GetKeys()[YcKey]);
			price = swaption->Price(&analyticModel);
			pf->SetPrice (price, i);
			delete swaption;
		}
	}
	else if (itsCalibrationType == DIAG_BASKET_CALIBRATION)
	{	
		ARM_StdPortfolioPtr varNotioPf	   = GetCalibMethod()->GetPortfolio();
		ARM_StdPortfolioPtr diagSwaptionPf = GetCalibMethod()->GetlinkedMethod()->GetPortfolio();

		
		/// compute price of diag swaptions
		for(int i=0;i<diagSwaptionPf->GetSize();++i)
		{
			swaption = static_cast< ARM_Swaption* >(diagSwaptionPf->GetAsset(i));
			swaption->SetModel(oswBSModel);
			price = swaption->ComputePrice();
			diagSwaptionPf->SetPrice (price, i);
		}
		
		/// compute price of var notion swaption
		ARM_StringVector keys(2);
		keys[0] = GetKeys()[YcKey];
		keys[1] = GetKeys()[OswModelKey];
		ARM_MarketData_ManagerRep mktDataManager(asOfDate);
		mktDataManager.RegisterData(GetKeys()[YcKey],cpnCurve);
		mktDataManager.RegisterData(GetKeys()[OswModelKey],oswBSModel);
		ARM_MarketIRModel analyticModel(mktDataManager, keys, itsVnsPricingMethod);
		analyticModel.SetZeroCurveKey(GetKeys()[YcKey]);
		analyticModel.SetBsModelKey(GetKeys()[OswModelKey]);
		
		for(i=0;i<varNotioPf->GetSize();++i)
		{
			ARM_VanillaArg* swaption = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(varNotioPf->GetAsset(i), asOfDate, GetKeys()[YcKey]);
			price = swaption->Price(&analyticModel);
			varNotioPf->SetPrice (price, i);
			delete swaption;
		}
	
	}
	else if (itsCalibrationType == DIAG_SPREAD_LONG || itsCalibrationType == DIAG_SPREAD_SHORT || itsCalibrationType == SHORT_LONG_SPREAD || itsCalibrationType == DIAG_LONG_CALIBRATION || itsCalibrationType == LONG_CALIBRATION)
	{
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " ARM_LocalCSOCalculator::ComputeSwaptionPortfolioPrices: invalid calibration type");
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeSwaptionPortfolioPrices
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSOPortfolioPrices()
{
	if (itsSoPricesComputed)
		return;

	if (	itsSpreadOptionFloorPF == ARM_StdPortfolioPtr(NULL)
		||  itsSpreadOptionFloorPF == ARM_StdPortfolioPtr(NULL))
		ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOPortfolioPrices : SO portfolios have not been created");


    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	   	
	ARM_SpreadOption* so;
	double price;
	
	/// floors
	for(int i=0;i<itsSpreadOptionFloorPF->GetSize();++i)
	{
		so = static_cast< ARM_SpreadOption* >(itsSpreadOptionFloorPF->GetAsset(i));
		so->SetModel(SOBSModel);
		price = so->ComputePrice();
		itsSpreadOptionFloorPF->SetPrice (price, i);
	}

	/// capns
	for(i=0;i<itsSpreadOptionCapPF->GetSize();++i)
	{
		so = static_cast< ARM_SpreadOption* >(itsSpreadOptionCapPF->GetAsset(i));
		so->SetModel(SOBSModel);
		price = so->ComputePrice();
		itsSpreadOptionCapPF->SetPrice (price, i);
	}

	if (IsSwitchable())
	{
		if (	itsSpreadOptionFloorPFforSwitch == ARM_StdPortfolioPtr(NULL)
			||  itsSpreadOptionFloorPFforSwitch == ARM_StdPortfolioPtr(NULL))
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOPortfolioPrices : SO portfolios have not been created");

		for(int i=0;i<itsSpreadOptionFloorPFforSwitch->GetSize();++i)
		{
			so = static_cast< ARM_SpreadOption* >(itsSpreadOptionFloorPFforSwitch->GetAsset(i));
			so->SetModel(SOBSModel);
			price = so->ComputePrice();
			itsSpreadOptionFloorPFforSwitch->SetPrice (price, i);
		}

		for(i=0;i<itsSpreadOptionCapPFforSwitch->GetSize();++i)
		{
			so = static_cast< ARM_SpreadOption* >(itsSpreadOptionCapPFforSwitch->GetAsset(i));
			so->SetModel(SOBSModel);
			price = so->ComputePrice();
			itsSpreadOptionCapPFforSwitch->SetPrice (price, i);
		}
	}

	itsSoPricesComputed = true;
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeCalibPortfolioPricesFor2F
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeCalibPortfolioPricesFor2F()
{
	if (itsCalibrationType == DIAG_SPREAD_LONG || itsCalibrationType == DIAG_SPREAD_SHORT)
	{
		ComputeSwPricesFor2F(GetCalibMethod()->GetPortfolio());
		ComputeSoPricesFor2F(GetCalibMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[1]);
		ComputeSoPricesFor2F(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[2]);
	}
	else if (itsCalibrationType == SHORT_LONG_SPREAD)
	{
		ComputeSoPricesFor2F(GetCalibMethod()->GetPortfolio(),itsCalibStrikeType[0]);
		ComputeSoPricesFor2F(GetCalibMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[2]);
		ComputeSoPricesFor2F(GetCalibMethod()->GetlinkedMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[1]);
	}
	else if (itsCalibrationType == DIAG_LONG_CALIBRATION)
	{
		ComputeSwPricesFor2F(GetCalibMethod()->GetPortfolio());
		ComputeSoPricesFor2F(GetCalibMethod()->GetlinkedMethod()->GetPortfolio(),itsCalibStrikeType[1]);
	}
	else if (itsCalibrationType == LONG_CALIBRATION)
	{
		ComputeSoPricesFor2F(GetCalibMethod()->GetPortfolio(),itsCalibStrikeType[0]);
	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeCalibPortfolioPricesFor2F
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSoPricesFor2F(ARM_StdPortfolioPtr& pf,CalibStrikeType cst)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	ARM_BSModel* SOBSModel;
	if (itsVnsPricingMethod == ARM_MarketIRModel::MONEYNESS)
		SOBSModel= dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	else
		SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	ARM_SpreadOption* spreadOption;
	double price,strike;
	
	double soFwdRate1,soFwdRate2,weight1,weight2,soFwd;
	ARM_ReferenceValue strikes;

	for (size_t i = 0; i < pf->GetSize(); ++i)
	{
		spreadOption = static_cast< ARM_SpreadOption* >(pf->GetAsset(i));
		if(spreadOption->GetNumFlows()!=1)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : SO in calibration portfolio must have a single flow" );
		spreadOption->SetModel(SOBSModel);

		/// Compute ATM strike
		price=spreadOption->ComputePrice();
		soFwdRate1 = (*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()))[0];
		soFwdRate2 = (*(spreadOption->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()))[0];
		weight1 = spreadOption->GetWeight1();
		weight2 = spreadOption->GetWeight2();
		soFwd = weight2*soFwdRate2 - weight1*soFwdRate1;

		if (cst==ZERO)
			strike=0.;
		else if (cst==ATM)
			strike=soFwd;
		else if (cst==CAP)
			strike=100*GetCpnCap().Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
		else if (cst==FLOOR)
			strike=100*GetCpnFloor().Interpolate((*(spreadOption->GetSpreadLeg()->GetFirstLeg()->GetResetDates()))[0]-asOfDate);
		else
			strike=0.;

		spreadOption->SetStrike(strike);
		strikes = ARM_ReferenceValue(strike);
		spreadOption->SetStrikes(&strikes);

		/// Compute & save price
		price=spreadOption->ComputePrice();
		pf->SetPrice(price, i);
	}
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeCalibPortfolioPricesFor2F
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSwPricesFor2F(ARM_StdPortfolioPtr& pf)
{
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

	ARM_Swaption* swaption;	
	double price;

	for(int i=0;i<pf->GetSize();++i)
	{
		swaption = static_cast< ARM_Swaption* >(pf->GetAsset(i));
		swaption->SetModel(oswBSModel);
		price = swaption->ComputePrice();
		pf->SetPrice (price, i);
	}
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::CreateAndSetCalibration()
{
	itsSoPricesComputed = false;
	CreateSOPortfolio();
	ComputeSOPortfolioPrices();				// prices for local vol

	CreateEmptyCalibration();				// portfolios
	ComputeSwaptionPortfolioPrices();		// prices for standard calib method

	if(itsModelType == ARM_PricingModelType::HWM2F)
		ComputeCalibPortfolioPricesFor2F();		// prices for 2F
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::CreateAndSetModel()
{	
	/// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	///---------------------------------------------------------
	/// Build HW model
	///---------------------------------------------------------
	ARM_PricingModelPtr refModel;

	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility,SIGMA_DEFAULT_VALUE ,"SIGMA" );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));

	if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid MeanReversion parameter");

	if (itsModelType == ARM_PricingModelType::HWM1F)
	{
		ARM_ModelParamVector paramVector(0);
		paramVector.push_back( &volParam );
		paramVector.push_back( mrsParam );
		
		ARM_ModelParamsHW1FStd params(paramVector);
		/// params are cloned in ARM_PricingModel constructor...
		refModel = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_HullWhite1F( CreateClonedPtr( curve ), &params ) ));
	}
	else if (itsModelType == ARM_PricingModelType::HWM2F)
	{		
		ARM_CurveModelParam* volRatioParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[VolRatioKey]));
		
		if(!volRatioParam || volRatioParam->GetType() != ARM_ModelParamType::VolatilityRatio)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid VolatilityRatio parameter");

		ARM_CurveModelParam* mrsSpreadParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsSpreadKey]));
		
		if(!mrsSpreadParam || mrsSpreadParam->GetType() != ARM_ModelParamType::MeanReversionSpread)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid MeanReversionSpread parameter");

		ARM_CurveModelParam* correlParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
		
		if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : invalid Correlation parameter");
		
		ARM_ModelParamVector paramVector(0);
		paramVector.push_back( &volParam );
		paramVector.push_back( mrsParam );
		paramVector.push_back( volRatioParam );
		paramVector.push_back( mrsSpreadParam );
		paramVector.push_back( correlParam );

		 if(volRatioParam->GetCurve()->GetAbscisses().size() <= 1 && correlParam->GetCurve()->GetAbscisses().size() <= 1 )
		{
			/// Constant correlation and 2nd factor volatility
			ARM_ModelParamsHW2FStd params(paramVector);
			refModel = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_HullWhite2F( CreateClonedPtr( curve ), &params )));
		}
        else
		{
			/// Time dependent correlation and/or 2nd factor volatility
			ARM_ModelParamsHW2FExt params(paramVector);
			refModel = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_HullWhite2F( CreateClonedPtr( curve ), &params )));
		}
	}
	else if (itsModelType == ARM_PricingModelType::QGM1F)
	{
		ARM_ModelParamVector paramVector(0);
		paramVector.push_back( &volParam );
		paramVector.push_back( mrsParam );
		ARM_CurveModelParam skewParam( ARM_ModelParamType::Skew,SKEW_DEFAULT_VALUE ,"SKEW" );
        paramVector.push_back(&skewParam);
		refModel = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_QGM1F(CreateClonedPtr( curve ),paramVector)));
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "Only HWM1F, HWM2F or QGM1F are avaliable to price CSO" );

	/// We have to use the H&W intial zero curve  to build the local models
	///---------------------------------------------------------
	/// Build the 2 local models (cap + floor)
	///---------------------------------------------------------
	ARM_Local_Normal_ModelParams* defaultParams = ARM_Local_Normal_Model::CreateDefaultModelParams();
	ARM_Local_Model* LocalFloorModel = new ARM_Local_Normal_Model (ARM_ZeroCurvePtr((ARM_ZeroCurve*)curve->Clone()), *defaultParams);
	ARM_Local_Model* LocalCapModel   = new ARM_Local_Normal_Model (ARM_ZeroCurvePtr((ARM_ZeroCurve*)curve->Clone()), *defaultParams);
	delete defaultParams;

	///---------------------------------------------------------
	/// Build MultiAssets
	///---------------------------------------------------------
	ARM_StringVector names (3);
	vector<ARM_PricingModelPtr> models (3);
	ARM_StringVectorVector depends(3);

	names[myRefModel]			=  GetKeys()[YcKey];
	names[myLocalFloorModel]	=  LOCAL_FLOOR_MODEL_NAME;
	names[myLocalCapModel]		=  LOCAL_CAP_MODEL_NAME;

	models[myRefModel]			= ARM_PricingModelPtr( refModel );
	models[myLocalFloorModel]	= ARM_PricingModelPtr( LocalFloorModel );
	models[myLocalCapModel]		= ARM_PricingModelPtr( LocalCapModel );
		
	// depends[0] = ARM_StringVector (0);
	depends[myLocalFloorModel]	= ARM_StringVector (1,names[myRefModel]);
	depends[myLocalCapModel]	= ARM_StringVector (1,names[myRefModel]);
	
	if(IsBasis())
	{
		names.resize(NbModels);
		models.resize(NbModels);
        depends.resize(NbModels);
        names[myBasisMarginModel] = GetKeys()[BasisKey];

		// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
		models[myBasisMarginModel] = ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve))));
		depends[myBasisMarginModel] = ARM_StringVector (1,names[myRefModel]);
	}


	ARM_ModelNameMap modelMap (names, models, depends);
		
	///	modelMap & correls are cloned in multi assets model
	ARM_PricingModelPtr pricingModel (new ARM_MultiAssetsModel ( &modelMap ) );
	
				
	///---------------------------------------------------------
	/// Set numerical method
	///---------------------------------------------------------
	/// compute nb steps for tree
	ARM_GP_Vector* exerciseDates = GetExerciseDateStrip()->GetResetDates() ;
	double lastEventTime = (*exerciseDates)[exerciseDates->size()-1] - asOfDate;
	int nbSteps ;
	if (itsModelDatas.size()>0)
	{
		if (itsModelDatas[0] > 250)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : tree nb steps per year is not relevant...");

		if (itsModelDatas[0] < DEFAULT_TREE_NBSTEPS_PER_YEAR)
			nbSteps = static_cast<int>(floor(DEFAULT_TREE_NBSTEPS_PER_YEAR*lastEventTime/K_YEAR_LEN));
		else
			nbSteps = static_cast<int>(floor(itsModelDatas[0]*lastEventTime/K_YEAR_LEN));
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, "CSO calculator : ModelsDatas is Empty...");
	

	int schedulerType=ARM_SchedulerBase::ConstantVarianceMeanReverting;
	ARM_GP_Vector schedulerDatas(3);
	schedulerDatas[0] = nbSteps;
	schedulerDatas[1] = 1;
	schedulerDatas[2] = 1.0e-3;
	int samplerType=ARM_SamplerBase::MeanReverting;
	ARM_GP_Vector samplerDatas(1,1.0e-3);
	int truncatorType=ARM_TruncatorBase::StandardDeviation;
	ARM_GP_Vector truncatorDatas(1,5.0);
	int reconnectorType=ARM_ReconnectorBase::Mean;
	int smootherType=ARM_SmootherBase::DoNothing;
	size_t nbDim = (itsModelType == ARM_PricingModelType::HWM1F) ? 1 : 2;
	ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(nbDim,schedulerType,schedulerDatas,
		samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
	pricingModel->SetNumMethod( ARM_NumMethodPtr( tree ) );

	///---------------------------------------------------------
	/// Create a Numeraire and set it
	///---------------------------------------------------------
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    pricingModel->SetNumeraire(numeraire);

	///---------------------------------------------------------
	/// Set the model
	///---------------------------------------------------------
	SetPricingModel(pricingModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::UpdateCalibration(bool isUpdateStrike)
{	
	/* 
	///  would work but in that case all ARM_SpreadOption are rebuilded -> slow
	itsSoPricesComputed = false;
	CreateEmptyCalibration();
	ComputeSwaptionPortfolioPrices();
	ComputeSOPortfolioPrices();
	*/
	/*if(IsBasis()){
		itsvFundSpread = ComputeDomesticBasis();
		ARM_GP_Vector spreadTimesNotio	 = itsvFundSpread * itsvFundNominal;
		ARM_GP_VectorPtr basisMargin = ARM_GP_VectorPtr((ARM_GP_Vector*)(basisConveter.ComputeDomMargin()).Clone());
			for (size_t i(0); i<itsExerSize; i++)
	{		
		/// define index
		CC_Ostringstream indexStr;
		/// sort by size
		i < 10 ? indexStr << 0 << i :  indexStr << i;
	
		/// names
		names.push_back("FundNominal"		+ indexStr.str()); /// first name
		names.push_back("SpreadTimesNotio"	+ indexStr.str()); /// ....
			}
	GetGenSecurity()->GetCstManager()->insert("FundMargin",ARM_GramFctorArg(basisMargin));

	}*/
	
	/// The portfolio of SO will not change, 
	/// we just have to recompute prices
	itsSoPricesComputed = false;
	ComputeSOPortfolioPrices();

	/// If DIAG ATM, the swaption set will not be modified, 
	/// But in other cases, the calib PF will change (equivalent strike, variable notional...)
	/// --> re-create the swaption portfolio (even in DIAG ATM case)
	
	if (itsCalibrationType == DIAG_CALIBRATION || itsCalibrationType == BASKET_CALIBRATION)
	{
		ARM_StdPortfolioPtr	swaptionPF = CreateSwaptionPortfolio(itsCalibrationType,itsCalibStrikeType[0], isUpdateStrike);
		if(itsCalibStrikeType[0] == FRONTIER)
		{
			SetDiagStrikeToFrontier(swaptionPF);
		}
		GetCalibMethod()->SetPortfolio(swaptionPF);
		ComputeSwaptionPortfolioPrices();
	}
	else if (itsCalibrationType == DIAG_BASKET_CALIBRATION)
	{
		ARM_StdPortfolioPtr diagPF     = CreateSwaptionPortfolio(DIAG_CALIBRATION,itsCalibStrikeType[0]);
		ARM_StdPortfolioPtr varNotioPF = CreateSwaptionPortfolio(BASKET_CALIBRATION,itsCalibStrikeType[0]);
		GetCalibMethod()->SetPortfolio(varNotioPF);
		GetCalibMethod()->GetlinkedMethod()->SetPortfolio(diagPF);
		ComputeSwaptionPortfolioPrices();
	}
	else if(itsModelType == ARM_PricingModelType::HWM2F)
	{
		ComputeCalibPortfolioPricesFor2F();
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, "Unknown calibration type or model");
	}
	
	
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::UpdateModel()
{
	/// get zc curve
	ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// set it to pricing model
	GetPricingModel()->SetZeroCurve(CreateClonedPtr( zcCurve ));
	
	/// this takes into account the new mean reversion if its has been changed
	/// un peu bourrin mais bon on est pas à ca près
	CreateAndSetModel();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::Calibrate()
{	
	///-------------------------------------
	/// HW calibration
	///-------------------------------------
	ARM_MultiAssetsModel*  refModel = static_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
	ARM_PricingModelPtr StochasticModel = (*refModel->GetModelMap())[GetKeys()[YcKey]]->Model();
	GetCalibMethod()->Calibrate(&*StochasticModel);

	///-------------------------------------
	/// Local models calibration
	///-------------------------------------
	ARM_PricingModelPtr model;
	model = (*refModel->GetModelMap())[LOCAL_FLOOR_MODEL_NAME]->Model();
	ARM_Local_Model* LocalFloorModel = static_cast<ARM_Local_Model*> (&*model);
	model = (*refModel->GetModelMap())[LOCAL_CAP_MODEL_NAME]->Model();
	ARM_Local_Model* LocalCapModel = static_cast<ARM_Local_Model*> (&*model);

	/// eval times
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_GP_Vector* exerciseDates = GetExerciseDateStrip()->GetResetDates() ;
	ARM_GP_Vector evalTimes(0);
	
	for (size_t i(0); i<itsExerSize; i++)
	{
		if (itsvIsExerDate[i])
			evalTimes.push_back((*exerciseDates)[i] - asOfDate);
	}
		
	LocalFloorModel->CalibrateLocalModel(*itsSpreadOptionFloorPF, evalTimes);
	LocalCapModel->CalibrateLocalModel(*itsSpreadOptionCapPF, evalTimes);

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : set zero in priced columns (only the first line will
///			 be filled)
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[CSO] = zeroValue;
    rowTypeVec[CSO] = ARM_DOUBLE;

	rowDescVec[Funding] = zeroValue;
    rowTypeVec[Funding] = ARM_DOUBLE;

	rowDescVec[Fix] = zeroValue;
    rowTypeVec[Fix] = ARM_DOUBLE;

	rowDescVec[Floor] = zeroValue;
    rowTypeVec[Floor] = ARM_DOUBLE;

	rowDescVec[Cap] = zeroValue;
    rowTypeVec[Cap] = ARM_DOUBLE;

	rowDescVec[ExerSwapRate] = zeroValue;
	rowTypeVec[ExerSwapRate] = ARM_DOUBLE;
	
	rowDescVec[Frontier] = zeroValue;
	rowTypeVec[Frontier] = ARM_DOUBLE;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_LocalCSOCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    ARM_DateStripPtr exerDateStrip = datesStructure.GetDateStrip(EXERCISE_SCHED);
	ARM_DateStripPtr cpnDateStrip  = itsStructDateStrip;
	ARM_DateStripPtr fundDateStrip = itsFundDateStrip;

	size_t eventSize = exerDateStrip->GetResetDates()->size();
    size_t descSize = sizeof(LocalCSOColNamesTable)/sizeof(LocalCSOColNamesTable[0]);
	if (itsComputeProba)
		descSize += 3 + 2*itsNbComputedProba;
	
	for(size_t firstEventIdx(0); firstEventIdx<itsExerSize; firstEventIdx++)
	{
		if(itsvIsExerDate[firstEventIdx]) 
			break;
	}
		
	// bool isFirstEvent = (eventIdx==itsFirstEventIdx);
	bool isFirstEvent = (eventIdx==firstEventIdx);
	bool isLastEvent  = (eventIdx==(eventSize-1));

    string cpnmodelName = GetKeys()[GetCpnModelKey()];
	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();   

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

	/// --- build cst manager names
	CC_Ostringstream indexStr;
	string SpreadTimesNotioStr	("SpreadTimesNotio");
	eventIdx < 10 ? indexStr << 0 << eventIdx :  indexStr << eventIdx;
	string FundNominalStr		("FundNominal"			+ indexStr.str());
	string LeverageTimesNotioStr("LeverageTimesNotio"	+ indexStr.str());
	string CpnNominalStr		("CpnNominal"			+ indexStr.str());
	string CpnLeverageLongStr	("CpnLeverageLong"		+ indexStr.str());
	string CpnLeverageShortStr	("CpnLeverageShort"		+ indexStr.str());
	string CpnFloorTimesNotioStr("CpnFloorTimesNotio"	+ indexStr.str());
	string CpnStrikePlusFloorStr("CpnStrikePlusFloor"	+ indexStr.str());
	string CpnStrikePlusCapStr  ("CpnStrikePlusCap"		+ indexStr.str());
	string CpnFloorFundTimesNotioStr("CpnFloorFundTimesNotio"	+ indexStr.str());
	string CpnStrikePlusFloorFundStr("CpnStrikePlusFloorFund"	+ indexStr.str());
	string CpnStrikePlusCapFundStr  ("CpnStrikePlusCapFund"		+ indexStr.str());
	
	/// Index tenors & type
	string longIndexType, longIndexTenor; 
	string shortIndexType, shortIndexTenor;

	string fundIndexType, fundIndexTenor; 
	/// LIBOR case
	if (IsLiborIndex((ARM_INDEX_TYPE)GetCMSFunding()))
	{
		fundIndexType = "LIBOR";
		ARM_IRIndex index ((ARM_INDEX_TYPE)GetCMSFunding());
		double tenor = index.GetYearTerm();
		int nmonth = int(tenor * 12);
		CC_Ostringstream fundLiborTenorDesc;
		fundLiborTenorDesc  << nmonth  << "m" ;
		fundIndexTenor = fundLiborTenorDesc.str();
	}
	/// CMS case
	else
	{
		fundIndexType = "CMS";
		CC_Ostringstream fundCmsTenorDesc;
		fundCmsTenorDesc  << GetCMSFunding()  - K_CMS1 + 1  << "y" ;
		fundIndexTenor = fundCmsTenorDesc.str();
	}
	
	/// LIBOR case
	if (IsLiborIndex((ARM_INDEX_TYPE)GetCMSLong()))
	{
		longIndexType = "LIBOR";
		ARM_IRIndex index ((ARM_INDEX_TYPE)GetCMSLong());
		double tenor = index.GetYearTerm();
		int nmonth = int(tenor * 12);
		CC_Ostringstream longLiborTenorDesc;
		longLiborTenorDesc  << nmonth  << "m" ;
		longIndexTenor = longLiborTenorDesc.str();
	}
	/// CMS case
	else
	{
		longIndexType = "CMS";
		CC_Ostringstream longCmsTenorDesc;
		longCmsTenorDesc  << GetCMSLong()  - K_CMS1 + 1  << "y" ;
		longIndexTenor = longCmsTenorDesc.str();
	}

	/// LIBOR case
	if (IsLiborIndex((ARM_INDEX_TYPE)GetCMSShort()))
	{
		shortIndexType = "LIBOR";
		ARM_IRIndex index ((ARM_INDEX_TYPE)GetCMSShort());
		double tenor = index.GetYearTerm();
		int nmonth = int(tenor * 12);
		CC_Ostringstream shortLiborTenorDesc;
		shortLiborTenorDesc  << nmonth  << "m" ;
		shortIndexTenor = shortLiborTenorDesc.str();
	}
	/// CMS case
	else
	{
		shortIndexType = "CMS";
		CC_Ostringstream shortCmsTenorDesc;
		shortCmsTenorDesc << GetCMSShort() - K_CMS1 + 1  << "y" ;
		shortIndexTenor = shortCmsTenorDesc.str();
	}

	// frequencies & daycounts
	string fundingFreq		= ARM_ArgConvReverse_MatFrequency.GetString(GetFundFreq());
	string fundingDayCount	= ARM_ArgConvReverse_DayCount.GetString(GetFundDaycount());
	string cpnFreq			= ARM_ArgConvReverse_MatFrequency.GetString(GetCpnFreq());
	string cpnDaycount		= ARM_ArgConvReverse_DayCount.GetString(GetCpnDaycount());
	string cpnResetTiming	= ARM_ArgConvReverse_Timing.GetString(GetCpnResetTiming());
		
	
	// EventDate (i.e. call date)
    double eventDate=(*(exerDateStrip->GetResetDates()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();

    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	// Start Date
	// funding start date is assumed to be the same...
    double startDate=(*(itsStructDateStrip->GetFlowStartDates()))[itsvCpnIndex[eventIdx]];
    CC_Ostringstream startDateDesc;
    startDateDesc << CC_NS(std,fixed) << startDate;
    rowDescVec[StartDate] = startDateDesc.str();
    rowTypeVec[StartDate] = ARM_DATE_TYPE;

	// End date (unadjusted!)
    CC_Ostringstream endDateDesc;
    endDateDesc << CC_NS(std,fixed) << GetEndDate().GetJulian();
    rowDescVec[EndDate] = endDateDesc.str();
    rowTypeVec[EndDate] = ARM_DATE_TYPE;


	// Non standard funding
	CC_Ostringstream FixFundLegDesc;
	CC_Ostringstream FloorFundLegDesc;
	CC_Ostringstream CapFundLegDesc;
	if (IsSwitchable())
	{
		string LocalcapModelStr = LOCAL_CAP_MODEL_NAME;
		string LocalfloorModelStr = LOCAL_FLOOR_MODEL_NAME;
		FixFundLegDesc << "ANNUITY(" << cpnmodelName << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
		FixFundLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  cpnFreq << ", " << fundingDayCount << ", "<<  CpnFloorFundTimesNotioStr <<")";	
		FloorFundLegDesc << "SPREADOPTION(" << LocalfloorModelStr << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
		FloorFundLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  fundIndexTenor << ", " << shortIndexTenor << ", ";
		FloorFundLegDesc << 1<< ", " << 0 << ", ";
		FloorFundLegDesc << CpnStrikePlusFloorFundStr << ", CAP, " << cpnFreq << ", " << fundingDayCount << ",,," << CpnNominalStr << "," << cpnResetTiming << "," << fundIndexType << ", " << shortIndexType << ")";
		CapFundLegDesc << "SPREADOPTION(" << LocalcapModelStr << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
		CapFundLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  fundIndexTenor << ", " << shortIndexTenor << ", ";
		CapFundLegDesc << 1 << ", " << 0 << ", ";
		CapFundLegDesc << CpnStrikePlusCapFundStr << ", CAP, " << cpnFreq << ", " << fundingDayCount << ",,," << CpnNominalStr << "," << cpnResetTiming << "," << fundIndexType << ", " << shortIndexType << ")";
	}
	else
	{
		FixFundLegDesc << "0";
		FloorFundLegDesc << "0";
		CapFundLegDesc << "0";
	}
	rowDescVec[FixFundLeg] = FixFundLegDesc.str();
	rowTypeVec[FixFundLeg] = ARM_STRING;
	rowDescVec[FloorFundLeg] = FloorFundLegDesc.str();
	rowTypeVec[FloorFundLeg] = ARM_STRING;
	rowDescVec[CapFundLeg] = CapFundLegDesc.str();
	rowTypeVec[CapFundLeg] = ARM_STRING;

    /// Funding Flow
	CC_Ostringstream FundingFlowDesc;
	if (IsSwitchable())
	{
		FundingFlowDesc << LocalCSOColNamesTable[FixFundLeg] << "[i]+";
		FundingFlowDesc << LocalCSOColNamesTable[FloorFundLeg] << "[i]-";
		FundingFlowDesc << LocalCSOColNamesTable[CapFundLeg] << "[i]";
	}
	else
	{
		FundingFlowDesc << "SWAP(" << cpnmodelName << "," << LocalCSOColNamesTable[StartDate] << "[i], "; 
		FundingFlowDesc << LocalCSOColNamesTable[EndDate] << "[i], 0, PAY,,," ;
		FundingFlowDesc << fundingFreq << ", " << fundingDayCount << ",, " << LeverageTimesNotioStr << ")";
		FundingFlowDesc << " + ANNUITY(" << cpnmodelName << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
		FundingFlowDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  fundingFreq << ", " << fundingDayCount << ", "<<  SpreadTimesNotioStr <<",,";
		FundingFlowDesc << itsvFundIndex[eventIdx] << ")";
	}
	rowDescVec[FundingLeg] = FundingFlowDesc.str();
	rowTypeVec[FundingLeg] = ARM_STRING;

	// FixLeg
	CC_Ostringstream FixLegDesc;
	FixLegDesc << "ANNUITY(" << cpnmodelName << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
	FixLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  cpnFreq << ", " << cpnDaycount << ", "<<  CpnFloorTimesNotioStr <<")";	
	rowDescVec[FixLeg] = FixLegDesc.str();
	rowTypeVec[FixLeg] = ARM_STRING;

	// FloorLeg
	string LocalfloorModelStr = LOCAL_FLOOR_MODEL_NAME;
	CC_Ostringstream FloorLegDesc;
	FloorLegDesc << "SPREADOPTION(" << LocalfloorModelStr << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
	FloorLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  longIndexTenor << ", " << shortIndexTenor << ", ";
	FloorLegDesc << CpnLeverageLongStr<< ", " << CpnLeverageShortStr << ", ";
	FloorLegDesc << CpnStrikePlusFloorStr << ", CAP, " << cpnFreq << ", " << cpnDaycount << ",,," << CpnNominalStr << "," << cpnResetTiming << "," << longIndexType << ", " << shortIndexType << ")";
	rowDescVec[FloorLeg] = FloorLegDesc.str();
	rowTypeVec[FloorLeg] = ARM_STRING;

	// CapLeg
	string LocalcapModelStr = LOCAL_CAP_MODEL_NAME;
	CC_Ostringstream CapLegDesc;
	CapLegDesc << "SPREADOPTION(" << LocalcapModelStr << "," << LocalCSOColNamesTable[StartDate] << "[i],"; 
	CapLegDesc << LocalCSOColNamesTable[EndDate] << "[i], " <<  longIndexTenor << ", " << shortIndexTenor << ", ";
	CapLegDesc << CpnLeverageLongStr << ", " << CpnLeverageShortStr << ", ";
	CapLegDesc << CpnStrikePlusCapStr << ", CAP, " << cpnFreq << ", " << cpnDaycount << ",,," << CpnNominalStr << "," << cpnResetTiming << "," << longIndexType << ", " << shortIndexType << ")";
	rowDescVec[CapLeg] = CapLegDesc.str();
	rowTypeVec[CapLeg] = ARM_STRING;
	
	// SOLeg
	CC_Ostringstream SOLegDesc;
	SOLegDesc << LocalCSOColNamesTable[FixLeg]		<< "[i] + ";  
	SOLegDesc << LocalCSOColNamesTable[FloorLeg]	<< "[i] - ";  
	SOLegDesc << LocalCSOColNamesTable[CapLeg]		<< "[i]";  
	rowDescVec[SOLeg] = SOLegDesc.str();
	rowTypeVec[SOLeg] = ARM_STRING;

	// Fees
	CC_Ostringstream feesDesc;
	feesDesc << CC_NS(std,fixed) << itsvExerFees[eventIdx];
	rowDescVec[Fees] = feesDesc.str();
	rowTypeVec[Fees] = ARM_DOUBLE;

	// Option
	CC_Ostringstream OptionDesc;
	/// RECEIVER CASE
	if (GetPayRec() == K_RCV)
	{
		OptionDesc << "MAX(" ;
		OptionDesc << LocalCSOColNamesTable[SOLeg]		<< "[i] - ";  
		OptionDesc << LocalCSOColNamesTable[FundingLeg]	<< "[i] - ";  
		OptionDesc << LocalCSOColNamesTable[Fees]		<< "[i], ";  
        !isLastEvent ? OptionDesc << "PV(" << LocalCSOColNamesTable[Option] << "[i+1]))": OptionDesc << "0)";  
	}
	/// PAYER CASE
	else
	{
		OptionDesc << "MAX(" ;
		OptionDesc << LocalCSOColNamesTable[FundingLeg]		<< "[i] - ";  
		OptionDesc << LocalCSOColNamesTable[SOLeg]			<< "[i] - ";  
		OptionDesc << LocalCSOColNamesTable[Fees]			<< "[i], ";  
        !isLastEvent ? OptionDesc << "PV(" << LocalCSOColNamesTable[Option] << "[i+1]))": OptionDesc << "0)";  
	}
	rowDescVec[Option] = OptionDesc.str();
	rowTypeVec[Option] = ARM_STRING;


	// Priced columns: CSO / Funding / Fix / Floor / Cap 
	if (isFirstEvent)
	{
		CC_Ostringstream CSODesc, FundingDesc, FixDesc, FloorDesc, CapDesc;

		CSODesc		<< LocalCSOColNamesTable[Option]		<< "[i]";
		rowDescVec[CSO] = CSODesc.str();
		rowTypeVec[CSO] = ARM_STRING;

		FundingDesc << LocalCSOColNamesTable[FundingLeg]	<< "[i]";
		rowDescVec[Funding] = FundingDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		FixDesc		<< LocalCSOColNamesTable[FixLeg]		<< "[i]";
		rowDescVec[Fix] = FixDesc.str();
		rowTypeVec[Fix] = ARM_STRING;

		FloorDesc	<< LocalCSOColNamesTable[FloorLeg]		<< "[i]";
		rowDescVec[Floor] = FloorDesc.str();
		rowTypeVec[Floor] = ARM_STRING;
		
		CapDesc		<< LocalCSOColNamesTable[CapLeg]		<< "[i]";
		rowDescVec[Cap] = CapDesc.str();
		rowTypeVec[Cap] = ARM_STRING;
	}


////////////////////////////////////////////////////////////////
/// FRONTIER COMPUTATION
///////////////////////////////////////////////////////////////
	bool isComputeFrontier = (itsCalibrationType == DIAG_CALIBRATION && itsCalibStrikeType[0] == FRONTIER);
	if(isComputeFrontier)
	{
		// SWAPRATE for exercise frontier projection
		CC_Ostringstream swaprateDesc;
		swaprateDesc << "SWAPRATE(" << cpnmodelName << ",";
		swaprateDesc << LocalCSOColNamesTable[StartDate] << "[i],";
		swaprateDesc << LocalCSOColNamesTable[EndDate] << "[i])";
		
		rowDescVec[ExerSwapRate] = swaprateDesc.str();
		rowTypeVec[ExerSwapRate] = ARM_STRING;

		// FRONTIER
		CC_Ostringstream frontierDesc;
		if(isLastEvent)
		{
			frontierDesc << "FRONTIER(" << LocalCSOColNamesTable[SOLeg] << "[i]";
			frontierDesc << "-" << LocalCSOColNamesTable[FundingLeg] << "[i]";
			frontierDesc << "-" << LocalCSOColNamesTable[Fees] << "[i],";
			frontierDesc << "0," << LocalCSOColNamesTable[ExerSwapRate] << "[i])";
		}
		else
		{
			frontierDesc << "FRONTIER(" << LocalCSOColNamesTable[SOLeg] << "[i]";
			frontierDesc << "-" << LocalCSOColNamesTable[FundingLeg] << "[i]";
			frontierDesc << "-" << LocalCSOColNamesTable[Fees] << "[i],";
			frontierDesc << LocalCSOColNamesTable[Option] << "[i+1],";
			frontierDesc << LocalCSOColNamesTable[ExerSwapRate] << "[i])";
		}
		rowDescVec[Frontier] = frontierDesc.str();
		rowTypeVec[Frontier] = ARM_STRING;
	}


////////////////////////////////////////////////////////////////
/// EXERCISE PROBABILITIES
///////////////////////////////////////////////////////////////
	if (itsComputeProba)
	{
		size_t probaOffset = Frontier;
		size_t nbCalls = exerDateStrip->GetResetDates()->size();

	// Reset [i+1]

		CC_Ostringstream reset_Plus1_Desc;
		double noticeDate_Plus1;
		if( eventIdx < nbCalls-1)
			noticeDate_Plus1 = exerDateStrip->GetResetDates()->Elt(eventIdx+1) ;
		else
			noticeDate_Plus1 = GetEndDate().GetJulian();

		reset_Plus1_Desc << CC_NS(std,fixed) << noticeDate_Plus1;
		rowDescVec[probaOffset + 1] = reset_Plus1_Desc.str();
		rowTypeVec[probaOffset + 1] = ARM_DATE_TYPE;

	// Final Date

		CC_Ostringstream FinalDate_Desc;
		FinalDate_Desc<< CC_NS(std,fixed) << GetEndDate().GetJulian();
		rowDescVec[probaOffset  + 2] = FinalDate_Desc.str();
		rowTypeVec[probaOffset  + 2] = ARM_DATE_TYPE;

	// Exercise Indicator 

		CC_Ostringstream Indicator_Desc;
		if(eventIdx == nbCalls-1)
		{
			if (GetPayRec() == K_RCV)
			{
				Indicator_Desc<< "if(" << LocalCSOColNamesTable[SOLeg] << "[i]-";
				Indicator_Desc <<  LocalCSOColNamesTable[Fees] <<"[i]>0,1,0)";
			}
			else
			{
				Indicator_Desc<< "if(0-" << LocalCSOColNamesTable[SOLeg] << "[i]-";
				Indicator_Desc <<  LocalCSOColNamesTable[Fees] <<"[i]>0,1,0)";
			}
		}
		else
		{
			if (GetPayRec() == K_RCV)
			{
				Indicator_Desc << "if(" << LocalCSOColNamesTable[SOLeg] << "[i]-"; 
				Indicator_Desc << LocalCSOColNamesTable[FundingLeg]	<< "[i] - ";  
				Indicator_Desc << LocalCSOColNamesTable[Fees] << "[i]-";
				Indicator_Desc << "PV(" << LocalCSOColNamesTable[Option] << "[i+1])>0,1,0)";
			}
			else
			{
				Indicator_Desc << "if(" << LocalCSOColNamesTable[FundingLeg] << "[i]-"; 
				Indicator_Desc << LocalCSOColNamesTable[SOLeg]	<< "[i] - ";  
				Indicator_Desc << LocalCSOColNamesTable[Fees] << "[i]-";
				Indicator_Desc << "PV(" << LocalCSOColNamesTable[Option] << "[i+1])>0,1,0)";
			}
		}
		rowDescVec[probaOffset  + 3] = Indicator_Desc.str();
		rowTypeVec[probaOffset  + 3] = ARM_STRING;

	// Probabilities

		size_t currProb = 4;
		int nbCols = 3 + 2*itsNbComputedProba;
		for(int i = 4; i<=nbCols;++i)
		{
			CC_Ostringstream Prob0_Desc;
			size_t maxIndex = (i == nbCols-1)?(nbCalls-1) : currProb-4;
			/////// Probi0
			if(eventIdx<maxIndex)
			{
				Prob0_Desc<<"MAX(ExerciseIndicator[i], PV("<< (LocalCSOProbaColNamesTable[i-1].c_str())<<"[i+1])/DF("<<cpnmodelName<<",ResetPlus1[i]))";
				rowTypeVec[probaOffset  + i] = ARM_STRING;
			}
			else if (eventIdx==maxIndex)
			{
				Prob0_Desc<<"ExerciseIndicator[i]";
				rowTypeVec[probaOffset  + i] = ARM_STRING;
			}
			else
			{
				Prob0_Desc<<"0";
				rowTypeVec[probaOffset  + i] = ARM_DOUBLE;
			}
			rowDescVec[probaOffset  + i] = Prob0_Desc.str();

			CC_Ostringstream Prob1_Desc;
			if(eventIdx<=maxIndex)
			{
				Prob1_Desc<<(LocalCSOProbaColNamesTable[i-1].c_str())<<"[i]*DF("<<cpnmodelName<<",FinalDate[i])";
				i++;
				rowTypeVec[probaOffset  + i] = ARM_STRING;
			}
			else
			{
				Prob1_Desc<<"0";
				i++;
				rowTypeVec[probaOffset  + i] = ARM_DOUBLE;
			}
			rowDescVec[probaOffset  + i] = Prob1_Desc.str();
			currProb++;
		}
	}
	return ARM_RowInfo(rowDescVec,rowTypeVec);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CSOCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the CSO
/////////////////////////////////////////////////////////////////
double ARM_LocalCSOCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

	bool vasyFrankie (false);
	int size = MIN(NbProductsToPrice,itsProductsToPrice.size());
	for (size_t i(0); i < size; ++i)
		if (itsProductsToPrice[i]) {vasyFrankie = true; break;}

	if (vasyFrankie)
	{
		ARM_GenPricer* genPricer = new ARM_GenPricer( &*GetGenSecurity(),&*GetPricingModel());

		genPricer->Price();
	
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

		double price;

		for (size_t i(0); i < size; ++i)
		{
			if (itsProductsToPrice[i])
			{ 
				price	= genPricer->GetPricerInfo()->GetContents( LocalCSOColNamesTable[ LocalCSOProductToPriceColums[i] ] ).GetData("Price").GetDouble();
								
				if (i == CSOPrice)
				{
					itsCSOPrice = price;
				}
				else if (i == FundingPrice)
				{
					itsFundingPrice = price;
				}
				else if (i == FixPrice)
				{
					itsFixPrice = price;
				}
				else if (i == FloorPrice)
				{
					itsFloorPrice = price;
				}
				else if (i == CapPrice)
				{
					itsCapPrice = price;
				}
			}
		}

		if (itsComputeProba)
		{
			ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
			double DF = pCurve->DiscountPrice(((GetEndDate().GetJulian())-(GetMktDataManager()->GetAsOfDate().GetJulian()))/K_YEAR_LEN);
			for(size_t i = 0; i < itsProba->size(); i++)
			{
				itsProba->Elt(i) = genPricer->GetPricerInfo()->GetContents( (LocalCSOProbaPricedColNamesTable[i].c_str())).GetData("IntermediatePrices").GetVector()->Elt(0)/ DF;
			}
			for(i=itsProba->size()-2;i>0;i--)
				itsProba->Elt(i) = itsProba->Elt(i)-itsProba->Elt(i-1);
			GetPricingData()["Proba"] = itsProba;
		}
	}

	///
	/// to add the pricing of analytic products, see callable snowball calculator
	///

	itsHasBeenPriced = true;

    return itsCSOPrice;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_LocalCSOCalculator::ColumnNames() const
{
	size_t colNamesSize_std = sizeof(LocalCSOColNamesTable)/sizeof(LocalCSOColNamesTable[0]);
	size_t colNamesSize_proba = 3 + 2*itsNbComputedProba;
	size_t colNamesSize = colNamesSize_std;
	if (itsComputeProba)
		colNamesSize += colNamesSize_proba; 

	vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize_std; ++i)
        colNamesVec[i] = LocalCSOColNamesTable[i];
	for(i=colNamesSize_std;i<colNamesSize;++i)
		colNamesVec[i] = LocalCSOProbaColNamesTable[i-colNamesSize_std];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_LocalCSOCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_LocalCSOCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[CSOPrice])
		GetPricingData()[ "CSO" ] = itsCSOPrice;
	
	if (itsProductsToPrice[FundingPrice])
		GetPricingData()[ "Funding" ] = itsFundingPrice;		

	if (itsProductsToPrice[FixPrice])
		GetPricingData()[ "Fix" ] = itsFixPrice;		

	if (itsProductsToPrice[FloorPrice])
		GetPricingData()[ "Floor" ] = itsFloorPrice;		
	
	if (itsProductsToPrice[CapPrice])
		GetPricingData()[ "Cap" ] = itsCapPrice;		
}



///-------------------------------------------------------------------------
///------ BASKET CALIBRATION STUFF -----------------------------------------
///-------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateVanillaArgSO
///	Returns: ARM_VanillaSpreadOptionArg
///	Action : create  a VanillaArg for an ARM_SpreadOption and generate the fix leg and float leg ofthe CMSs 
/////////////////////////////////////////////////////////////////
ARM_VanillaSpreadOptionArg* ARM_LocalCSOCalculator::CreateVanillaArgSO(ARM_SpreadOption* spreadOption) 
{
	double asOfDate	  = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_Currency* ccy = GetCurrencyUnit();
	ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)ARM_ConverterFromKernel::ConvertSecuritytoArgObject(spreadOption,asOfDate);
	soVanillaArg->ComputeIndexSchedulesAndAdjustDates(ccy, asOfDate);
	return soVanillaArg;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeSOSensitivities
///	Returns: void
///	Action : compute so cap/floor cpn sensis w.r.t. long & short CMS
///			 NB: So portfolio prices are supposed to be already computed
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSOSensitivities()
{
	/// ---- NOTE ------
	/// itsCpnLongSensi[i] and itsCpnShortSensi[i] are the sensitivities
	/// of coupon[i] with respect to long & short CMS rates. 
	/// >>> IT DOES NOT TAKE INTO ACCOUNT notional, interest period and payment discount factor !!!
	/// Same thing for itsCpnValue
	/// 
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	ARM_ZeroCurve* curve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		
	ARM_SpreadOption* soFloor;
	ARM_SpreadOption* soCap;

	size_t sizeSO = itsSpreadOptionFloorPF->size();

	itsCpnLongSensi.resize(itsCpnSize);
	itsCpnShortSensi.resize(itsCpnSize);
	itsCpnValue.resize(itsCpnSize);

	size_t soIdx (0);
	double fixRate;

	const double LOCAL_SENSI_CMS = 100.0 * SENSI_CMS;
		
	for (size_t i (0); i<itsCpnSize; i++)
	{
		if (IsFixCoupon(i, fixRate))
		{
			itsCpnShortSensi[i] = 0.0;
			itsCpnLongSensi[i]  = 0.0;
			itsCpnValue[i]	    = fixRate;
		}
		else
		{
			if (soIdx>=sizeSO || ( (i==itsCpnSize-1) && (soIdx !=sizeSO-1) ) )
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOSensitivities : unresolved problem" );

			soFloor = static_cast< ARM_SpreadOption* >(itsSpreadOptionFloorPF->GetAsset(soIdx));
			soCap   = static_cast< ARM_SpreadOption* >(itsSpreadOptionCapPF->GetAsset(soIdx));
			soIdx ++;

			double shortFwd	 = soFloor->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(0);
			double longFwd   = soFloor->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(0);
			double correl	 = 0.01 * soFloor->GetCorrelVector()->Elt(0);
			double coeffShort= soFloor->GetWeight1();
			double coeffLong = soFloor->GetWeight2();
			// double optMat	 = (soFloor->GetMaturity().GetJulian()-asOfDate)/K_YEAR_LEN;
			double optMat = (GetStructDateStrip()->GetResetDates()->Elt(i) - asOfDate)/K_YEAR_LEN;

			double strikeFloor   = soFloor->GetStrikes()->GetDiscreteValues()->Elt(0);
			double shortVolFloor = 0.01 * soFloor->GetVol1Vector()->Elt(0);
			double longVolFloor  = 0.01 * soFloor->GetVol2Vector()->Elt(0);
			
			double strikeCap   = soCap->GetStrikes()->GetDiscreteValues()->Elt(0);
			double shortVolCap = 0.01 * soCap->GetVol1Vector()->Elt(0);
			double longVolCap  = 0.01 * soCap->GetVol2Vector()->Elt(0);

			double sensi;
			double floor0, cap0;

			///---------------------------------
			/// initial prices
			///---------------------------------
			/// it is safer to recompute them instead of taking the price
			/// already computed in the portfolio
			floor0 = SOBSModel->SpreadOptionPrice ( coeffShort*shortFwd,
													coeffLong*longFwd,
													shortVolFloor, 
													longVolFloor,
													correl, 
													strikeFloor, 
													optMat, 
													K_CAP);
			
			cap0  = SOBSModel->SpreadOptionPrice (  coeffShort*shortFwd,
													coeffLong*longFwd,
													shortVolCap, 
													longVolCap,
													correl, 
													strikeCap, 
													optMat, 
													K_CAP);

			///---------------------------------
			/// sensi w.r.t short CMS rate
			///---------------------------------
			sensi  = SOBSModel->SpreadOptionPrice ( coeffShort*(shortFwd+LOCAL_SENSI_CMS),
													coeffLong*longFwd,
													shortVolFloor, 
													longVolFloor,
													correl, 
													strikeFloor, 
													optMat, 
													K_CAP);

			sensi -= floor0;

			sensi += cap0;

			sensi -= SOBSModel->SpreadOptionPrice ( coeffShort*(shortFwd+LOCAL_SENSI_CMS),
													coeffLong*longFwd,
													shortVolCap, 
													longVolCap,
													correl, 
													strikeCap, 
													optMat, 
													K_CAP);

			sensi /= SENSI_CMS;

			itsCpnShortSensi[i] = sensi;


			///---------------------------------
			/// sensi w.r.t long CMS rate
			///---------------------------------
			sensi  = SOBSModel->SpreadOptionPrice ( coeffShort*shortFwd,
													coeffLong*(longFwd+LOCAL_SENSI_CMS),
													shortVolFloor, 
													longVolFloor,
													correl, 
													strikeFloor, 
													optMat, 
													K_CAP);

			sensi -= floor0;

			sensi += cap0;

			sensi -= SOBSModel->SpreadOptionPrice ( coeffShort*shortFwd,
													coeffLong*(longFwd+LOCAL_SENSI_CMS),
													shortVolCap, 
													longVolCap,
													correl, 
													strikeCap, 
													optMat, 
													K_CAP);

			sensi /= SENSI_CMS;

			itsCpnLongSensi[i] = sensi;

						
			// sacré kernel !
			itsCpnLongSensi[i]	*= 0.01;
			itsCpnShortSensi[i] *= 0.01;
			
			/// cpn value
			floor0 *= 0.01;
			cap0   *= 0.01;

			itsCpnValue[i] = itsvCpnFloor[i] + floor0 - cap0;

			
			
		}
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeSONormalSensitivities
///	Returns: void
///	Action : compute so cap/floor cpn sensis w.r.t. long & short CMS
///			 under the NORMAL ASSUMPTION on spread dynamics
///			 NB: So portfolio prices are supposed to be already computed
/////////////////////////////////////////////////////////////////
void ARM_LocalCSOCalculator::ComputeSONormalSensitivities()
{
	/// ---- NOTE ------
	/// itsCpnLongSensi[i] and itsCpnShortSensi[i] are the sensitivities
	/// of coupon[i] with respect to long & short CMS rates. 
	/// >>> IT DOES NOT TAKE INTO ACCOUNT notional, interest period and payment discount factor !!!
	/// Same thing for itsCpnValue
	/// 
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	ARM_ZeroCurve* curve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		
	ARM_SpreadOption* soFloor;
	ARM_SpreadOption* soCap;

	size_t sizeSO = itsSpreadOptionFloorPF->size();

	itsCpnLongSensi.resize(itsCpnSize);
	itsCpnShortSensi.resize(itsCpnSize);
	itsCpnValue.resize(itsCpnSize);

	size_t soIdx (0);
	double fixRate;

	const double NSTDEV_NO_IMPLIED_VOL = 6.0;
		
	for (size_t i (0); i<itsCpnSize; i++)
	{
		double optMat  = (GetStructDateStrip()->GetResetDates()->Elt(i) - asOfDate)/K_YEAR_LEN;
		if(optMat <- ARM_NumericConstants::ARM_DOUBLE_TOLERENCE) continue;
		if (IsFixCoupon(i, fixRate))
		{
			itsCpnShortSensi[i] = 0.0;
			itsCpnLongSensi[i]  = 0.0;
			itsCpnValue[i]	    = fixRate;
		}
		else
		{
			if (soIdx>=sizeSO || ( (i==itsCpnSize-1) && (soIdx !=sizeSO-1) ) )
				ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOSensitivities : unresolved problem" );

			soFloor = static_cast< ARM_SpreadOption* >(itsSpreadOptionFloorPF->GetAsset(soIdx));
			soCap   = static_cast< ARM_SpreadOption* >(itsSpreadOptionCapPF->GetAsset(soIdx));
			
			/// get corresponding vanilla arg
			ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVect[i]);
			
			double shortFwd	 = 0.01 * soFloor->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(0);
			double longFwd   = 0.01 * soFloor->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(0);
			double correl	 = 0.01 * soFloor->GetCorrelVector()->Elt(0);
			double coeffShort= soFloor->GetWeight1();
			double coeffLong = soFloor->GetWeight2();			
			
			double strikeFloor = 0.01 * soFloor->GetStrikes()->GetDiscreteValues()->Elt(0);
			double strikeCap   = 0.01 * soCap->GetStrikes()->GetDiscreteValues()->Elt(0);

			double longTenor  = (soVanillaArg->GetSwapLongFloatEndTime()->Elt(0) - soVanillaArg->GetSwapLongFloatStartTime()->Elt(0)) / K_YEAR_LEN;
			double shortTenor = (soVanillaArg->GetSwapShortFloatEndTime()->Elt(0) - soVanillaArg->GetSwapShortFloatStartTime()->Elt(0)) / K_YEAR_LEN;

			double atmShortVol  =  0.01 * soFloor->GetVol1ATM(optMat, shortTenor);
			double atmLongVol   =  0.01 * soFloor->GetVol2ATM(optMat, longTenor);
			atmShortVol	*= coeffShort * shortFwd ;
			atmLongVol	*= coeffLong  * longFwd;
			double atmSpreadVol =  sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
			double atmSpreadStdev = atmSpreadVol * sqrt(optMat);
			
			double forwardSpread = coeffLong * longFwd - coeffShort * shortFwd;

			double payDf		= curve->DiscountPrice(soVanillaArg->GetPayTimes()->Elt(0)/K_YEAR_LEN);
			double payNotional	= soVanillaArg->GetNotional()->Elt(0);
			double payPeriod	= soVanillaArg->GetPayPeriods()->Elt(0);

			/// take vols ATM
			double floorVol = atmSpreadVol;
			double capVol	= atmSpreadVol;

						
			/// overwrite only if not too far from ATM
			if (fabs(strikeFloor - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = soFloor->IsCap()-soFloor->IsFloor();
				double targetPrice_N = itsSpreadOptionFloorPF->GetMktPrices()->Elt(soIdx);
				targetPrice_N /= (payDf * payNotional * payPeriod);
				floorVol = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strikeFloor, optMat, callPut, &atmSpreadVol);
			}

			/// overwrite only if not too far from ATM
			if (fabs(strikeCap - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
			{
				int callPut = soCap->IsCap()-soCap->IsFloor();
				double targetPrice_N = itsSpreadOptionCapPF->GetMktPrices()->Elt(soIdx);
				targetPrice_N /= (payDf * payNotional * payPeriod);
				capVol = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strikeCap, optMat, callPut, &atmSpreadVol);
			}


			///---------------------------------
			/// initial prices
			///---------------------------------
			/// it is safer to recompute them instead of taking the price
			/// already computed in the portfolio
			double sensi;
			double floor0, cap0;
			
			floor0 = VanillaOption_N (forwardSpread, floorVol, strikeFloor, optMat, K_CALL);
			cap0   = VanillaOption_N (forwardSpread, capVol,   strikeCap,   optMat, K_CALL);
			
			
			///---------------------------------
			/// sensi w.r.t short CMS rate
			///---------------------------------
			sensi  = VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, floorVol, strikeFloor, optMat, K_CALL);
			sensi -= floor0;

			sensi += cap0;
			sensi -= VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, capVol,  strikeCap,   optMat, K_CALL);

			sensi /= SENSI_CMS;

			itsCpnShortSensi[i] = sensi;


			///---------------------------------
			/// sensi w.r.t long CMS rate
			///---------------------------------
			sensi  = VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, floorVol, strikeFloor, optMat, K_CALL);
			sensi -= floor0;

			sensi += cap0;
			sensi -= VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, capVol,  strikeCap,   optMat, K_CALL);

			sensi /= SENSI_CMS;

			itsCpnLongSensi[i] = sensi;


			///---------------------------------
			/// Cpn value
			///---------------------------------
			itsCpnValue[i] = itsvCpnFloor[i] + floor0 - cap0;


			///---------------------------------
			/// don't forget to soIdx ++
			///---------------------------------
			soIdx ++;
		}
	}
	if (IsSwitchable())
		ComputeExtraSONormalSensitivitiesForSwitch();
}

void ARM_LocalCSOCalculator::ComputeExtraSONormalSensitivitiesForSwitch()
{
	ARM_BSModel* SOBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[SoModelKey]) );
	ARM_ZeroCurve* curve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		
	ARM_SpreadOption* soFloor;
	ARM_SpreadOption* soCap;

	size_t sizeSO = itsSpreadOptionFloorPFforSwitch->size();

	itsCpnLongSensiForSwitch.resize(itsCpnSize);
	itsCpnShortSensiForSwitch.resize(itsCpnSize);
	itsCpnValueForSwitch.resize(itsCpnSize);

	size_t soIdx (0);

	const double NSTDEV_NO_IMPLIED_VOL = 6.0;
		
	for (size_t i (0); i<itsCpnSize; i++)
	{
		double optMat  = (GetStructDateStrip()->GetResetDates()->Elt(i) - asOfDate)/K_YEAR_LEN;
		if(optMat <- ARM_NumericConstants::ARM_DOUBLE_TOLERENCE) continue;
		
		if (soIdx>=sizeSO || ( (i==itsCpnSize-1) && (soIdx !=sizeSO-1) ) )
			ARM_THROW( ERR_INVALID_ARGUMENT, "ARM_LocalCSOCalculator::ComputeSOSensitivities : unresolved problem" );

		soFloor = static_cast< ARM_SpreadOption* >(itsSpreadOptionFloorPFforSwitch->GetAsset(soIdx));
		soCap   = static_cast< ARM_SpreadOption* >(itsSpreadOptionCapPFforSwitch->GetAsset(soIdx));
			
		/// get corresponding vanilla arg
		ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVectForSwitch[i]);
			
		double shortFwd	 = 0.01 * soFloor->GetSpreadLeg()->GetFirstLeg()->GetFwdRates()->Elt(0);
		double longFwd   = 0.01 * soFloor->GetSpreadLeg()->GetSecondLeg()->GetFwdRates()->Elt(0);
		double correl	 = 0.01 * soFloor->GetCorrelVector()->Elt(0);
		double coeffShort= soFloor->GetWeight1();
		double coeffLong = soFloor->GetWeight2();			
			
		double strikeFloor = 0.01 * soFloor->GetStrikes()->GetDiscreteValues()->Elt(0);
		double strikeCap   = 0.01 * soCap->GetStrikes()->GetDiscreteValues()->Elt(0);

		double longTenor  = (soVanillaArg->GetSwapLongFloatEndTime()->Elt(0) - soVanillaArg->GetSwapLongFloatStartTime()->Elt(0)) / K_YEAR_LEN;
		double shortTenor = (soVanillaArg->GetSwapShortFloatEndTime()->Elt(0) - soVanillaArg->GetSwapShortFloatStartTime()->Elt(0)) / K_YEAR_LEN;

		double atmShortVol  =  0.01 * soFloor->GetVol1ATM(optMat, shortTenor);
		double atmLongVol   =  0.01 * soFloor->GetVol2ATM(optMat, longTenor);
		atmShortVol	*= coeffShort * shortFwd ;
		atmLongVol	*= coeffLong  * longFwd;
		double atmSpreadVol =  sqrt( atmShortVol * atmShortVol + atmLongVol * atmLongVol - 2.* correl * atmShortVol * atmLongVol );
		double atmSpreadStdev = atmSpreadVol * sqrt(optMat);
			
		double forwardSpread = coeffLong * longFwd - coeffShort * shortFwd;

		double payDf		= curve->DiscountPrice(soVanillaArg->GetPayTimes()->Elt(0)/K_YEAR_LEN);
		double payNotional	= soVanillaArg->GetNotional()->Elt(0);
		double payPeriod	= soVanillaArg->GetPayPeriods()->Elt(0);

		/// take vols ATM
		double floorVol = atmSpreadVol;
		double capVol	= atmSpreadVol;

						
		/// overwrite only if not too far from ATM
		if (fabs(strikeFloor - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
		{
			int callPut = soFloor->IsCap()-soFloor->IsFloor();
			double targetPrice_N = itsSpreadOptionFloorPFforSwitch->GetMktPrices()->Elt(soIdx);
			targetPrice_N /= (payDf * payNotional * payPeriod);
			floorVol = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strikeFloor, optMat, callPut, &atmSpreadVol);
		}

		/// overwrite only if not too far from ATM
		if (fabs(strikeCap - forwardSpread)<NSTDEV_NO_IMPLIED_VOL * atmSpreadStdev)
		{
			int callPut = soCap->IsCap()-soCap->IsFloor();
			double targetPrice_N = itsSpreadOptionCapPFforSwitch->GetMktPrices()->Elt(soIdx);
			targetPrice_N /= (payDf * payNotional * payPeriod);
			capVol = VanillaImpliedVol_N (forwardSpread, targetPrice_N, strikeCap, optMat, callPut, &atmSpreadVol);
		}


		///---------------------------------
		/// initial prices
		///---------------------------------
		/// it is safer to recompute them instead of taking the price
		/// already computed in the portfolio
		double sensi;
		double floor0, cap0;
			
		floor0 = VanillaOption_N (forwardSpread, floorVol, strikeFloor, optMat, K_CALL);
		cap0   = VanillaOption_N (forwardSpread, capVol,   strikeCap,   optMat, K_CALL);
			
			
		///---------------------------------
		/// sensi w.r.t short CMS rate
		///---------------------------------
		sensi  = VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, floorVol, strikeFloor, optMat, K_CALL);
		sensi -= floor0;

		sensi += cap0;
		sensi -= VanillaOption_N (forwardSpread - coeffShort * SENSI_CMS, capVol,  strikeCap,   optMat, K_CALL);

		sensi /= SENSI_CMS;

		itsCpnShortSensiForSwitch[i] = sensi;


		///---------------------------------
		/// sensi w.r.t long CMS rate
		///---------------------------------
		sensi  = VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, floorVol, strikeFloor, optMat, K_CALL);
		sensi -= floor0;

		sensi += cap0;
		sensi -= VanillaOption_N (forwardSpread + coeffLong * SENSI_CMS, capVol,  strikeCap,   optMat, K_CALL);

		sensi /= SENSI_CMS;

		itsCpnLongSensiForSwitch[i] = sensi;


		///---------------------------------
		/// Cpn value
		///---------------------------------
		itsCpnValueForSwitch[i] = (-1000 - itsvFundSpread[i]) + floor0 - cap0;


		///---------------------------------
		/// don't forget to soIdx ++
		///---------------------------------
		soIdx ++;
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: CreateVarNotionalSwaptionAtExer
///	Returns: void
///	Action : create the variable notional swaption
///			 
/////////////////////////////////////////////////////////////////
ARM_Swaption* ARM_LocalCSOCalculator::CreateVarNotionalSwaptionAtExer (int exerIdx, ARM_Swaption*& fwdStartSwaption,CalibStrikeType calibStrikeType,bool isUpdateStrike)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_GP_Vector* exerDates = itsExerciseDateStrip->GetResetDates();	
	double ExerciseDate = (*exerDates)[exerIdx];	

	/// get defaults swaption params for currency
	ARM_Currency* ccy		 = GetCurrencyUnit();
	int spotDays			 = ccy->GetSpotDays(); 
	ARM_INDEX_TYPE indexType = ccy->GetVanillaIndexType();
	char* resetCalendar		 = ccy->GetResetCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdresetCalendar(resetCalendar);
	char* payCalendar		 = ccy->GetPayCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdpayCalendar(payCalendar);
	int stdFixFreq			 = ccy->GetFixedPayFreq();
	int stdFixDayCount		 = ccy->GetFixedDayCount();
	int resetGap			 = - spotDays;
	int fwdRule				 = K_MOD_FOLLOWING;	// for forward dates
    int intRule				 = K_ADJUSTED;		// for interest dates
    int stubRule			 = K_SHORTSTART;
	int resetTiming			 = K_ADVANCE;
	int payTiming			 = K_ARREARS;
		
	// compute start date
	ARM_Date startDate(ExerciseDate);
	startDate.GapBusinessDay(spotDays,resetCalendar);
	
		
	/// find (unadjusted) end date for variable notional swaption
	ARM_Date realEndDate ;
	ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVect[itsCpnSize - 1]);
	if (GetCMSLong()>GetCMSShort())
		realEndDate = soVanillaArg->GetSwapLongFloatEndTime()->Elt(0) + asOfDate;
	else
		realEndDate = soVanillaArg->GetSwapShortFloatEndTime()->Elt(0) + asOfDate;
	
	ARM_Date endDate(startDate), prevEndDate;

	while (endDate<realEndDate)
	{
		prevEndDate = endDate;
		endDate.AddMonths( 12/stdFixFreq ) ;
	}

	if (fabs(prevEndDate.GetJulian()-realEndDate.GetJulian()) < fabs(endDate.GetJulian()-realEndDate.GetJulian()) )
		endDate = prevEndDate;




	/// create date strip
	ARM_DateStrip dateStrip  (	startDate,endDate,stdFixFreq,stdFixDayCount,
								resetCalendar,
								fwdRule,intRule,stubRule,
								resetGap,stdFixFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);	

	/// instanciate stripper
	ARM_STRIPPER stripper(asOfDate, ARM_DateStripPtr(new ARM_DateStrip(dateStrip)));

	ARM_GP_Vector* resetDates	 = dateStrip.GetResetDates();
	ARM_GP_Vector* startDates	 = dateStrip.GetFlowStartDates();
	ARM_GP_Vector* endDates		 = dateStrip.GetFlowEndDates();
	ARM_GP_Vector* interestTerms = dateStrip.GetInterestTerms();


	/// compute start df + fwd swap rates to feed stripper
	double startDf = curve->DiscountPrice(((*startDates)[0]-asOfDate)/K_YEAR_LEN);
	
	int i, j, k, sizeSched = startDates->size(); 

	double level (0.0);
	double df;
	ARM_GP_Vector forwards(sizeSched);

	stripper.setStartDf(startDf);

	for (j=0; j<sizeSched; j++)
	{
		df = curve->DiscountPrice(((*endDates)[j]-asOfDate)/K_YEAR_LEN);
		level += df * (*interestTerms)[j];
		forwards[j] = (startDf - df) / level;
		stripper.setSwapRate (j, forwards[j]);
	}
	
	// strip, given start df and forward swap rates
	stripper.strip();
 
	///
	/// Fixed Leg notios
	///
	/// create notional
	ARM_Vector fixAbs(sizeSched);
	ARM_Vector fixNotios(sizeSched);

	// compute fixIndex so that we are as close as possible from product end date
	size_t fixIndex (0);

	if ((*endDates)[sizeSched-1]<=GetEndDate().GetJulian())
	{
		fixIndex = sizeSched - 1;
	}
	else
	{
		while ((*endDates)[fixIndex]<GetEndDate().GetJulian())
			fixIndex ++;


		if (	(fixIndex>0)
			&&  ( (*endDates)[fixIndex] - GetEndDate().GetJulian() >= GetEndDate().GetJulian()  - (*endDates)[fixIndex-1] )  )
			fixIndex = fixIndex - 1;
	}
	
	for (j=0; j<sizeSched; j++)
	{
		fixAbs[j] = (*endDates)[j];
		
		if (j>fixIndex)
			fixNotios[j] = 0.0;
		else
			/// rather approximative ...
			fixNotios[j] = GetCpnNominal().Interpolate((*resetDates)[j] - asOfDate);
	}

	// precomputation
	ARM_GP_Vector A(itsCpnSize), Bshort(itsCpnSize), Blong (itsCpnSize);
	ARM_GP_Vector Afund(itsCpnSize), Bfund(itsCpnSize);

	ARM_GP_Vector* cpnCoverages = GetStructDateStrip()->GetInterestTerms();
	ARM_GP_Vector* cpnPaydates  = GetStructDateStrip()->GetPaymentDates();
	double soLegPv (0.0),soLegPvForSwitch (0.0);
	for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
	{
		df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
		A[i]		= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsCpnValue[i];
		Blong[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsCpnLongSensi[i];
		Bshort[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsCpnShortSensi[i];
		soLegPv		+= A[i] * df;

		if (IsSwitchable())
		{
			Afund[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsCpnValueForSwitch[i];
			Bfund[i]	= (*cpnCoverages)[i] * itsvCpnNominal[i] * df * itsCpnLongSensiForSwitch[i];
			soLegPvForSwitch		+= Afund[i] * df;
		}
	}

	/// take exercise fees into account
	if ( GetPayRec()==K_RCV )
		soLegPv -= itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);
	else
		soLegPv += itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);

	ARM_GP_Vector df0(itsCpnSize), Slong0(itsCpnSize), Sshort0(itsCpnSize), Sfund0(itsCpnSize);
	double dfFee0;

	double Numeraire0 (0.0);
	ARM_GP_Vector NumerizedLiborFlow0 (sizeSched, 0.0);

	/// Notional used for normalization
	double fixNotional (1.0);
	for (j=0; j<sizeSched; j++)
	{
		if (fixNotios[j])
		{
			fixNotional = fixNotios[j];
			break;
		}
	}
	

	for (j=0; j<sizeSched; j++)
		Numeraire0 += (*interestTerms)[j] * (fixNotios[j]/fixNotional) * stripper.df((*endDates)[j]-asOfDate) ;
	
	for (j=0; j<sizeSched; j++)
		NumerizedLiborFlow0[j] =  (stripper.df((*startDates)[j]-asOfDate) -  stripper.df((*endDates)[j]-asOfDate)) / Numeraire0;

	/// initial values for dfs and swaprates
	double notUsed;
	for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
	{

		df0[i] = stripper.df((*cpnPaydates)[i]-asOfDate);

		if (!IsFixCoupon(i, notUsed))
		{
			ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVect[i]);
					
			double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(0);
			double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(0);
			ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[0];
			ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[0];
			Slong0[i] = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);

			double shortStartTime = soVanillaArg->GetSwapShortFloatStartTime()->Elt(0);
			double shortEndTime   = soVanillaArg->GetSwapShortFloatEndTime()->Elt(0);
			ARM_GP_Vector* shortPayTimes    = soVanillaArg->GetSwapShortFixPayTimes()[0];
			ARM_GP_Vector* shortPayPeriods  = soVanillaArg->GetSwapShortFixPayPeriods()[0];
			Sshort0[i] = stripper.swapRate(shortStartTime, shortEndTime, *shortPayTimes, *shortPayPeriods);
		}
		else
		{
			Slong0[i]  = 0.0;
			Sshort0[i] = 0.0;
		}

		if (IsSwitchable())
		{
			ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVectForSwitch[i]);
					
			double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(0);
			double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(0);
			ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[0];
			ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[0];
			Sfund0[i] = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);
		}
	
	}

	dfFee0 =  stripper.df(ExerciseDate-asOfDate);


	/// initial value for funding leg
	//size_t fundsize = itsFundSize-itsvFundIndex[exerIdx];
	ARM_GP_Vector fundingStarts (0) ;
	ARM_GP_Vector fundingEnds	(0) ;
	ARM_GP_Vector fundingPeriods(0) ;
	ARM_GP_Vector fundingNotios (0) ;
	ARM_GP_Vector fundingMargin (0) ;
	ARM_GP_Vector fundingLeverage (0) ;

	for (i=itsvFundIndex[exerIdx],j=0; i<itsFundSize; i++,++j)
	{
		fundingStarts.push_back(GetFundDateStrip()->GetFlowStartDates()->Elt(i) - asOfDate) ;
		fundingEnds.push_back(GetFundDateStrip()->GetFlowEndDates()->Elt(i) - asOfDate );
		fundingPeriods.push_back(GetFundDateStrip()->GetInterestTerms()->Elt(i));
		fundingMargin.push_back(itsvFundSpread[i]);
		fundingLeverage.push_back(itsvFundLeverage[i]);
		fundingNotios.push_back(itsvFundNominal[i]);
	}
	
	double floatLeg0 = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargin, &fundingLeverage);

	/// donc forget to store underlying PV for later strike determination
	double underlyingPv;
	if (!IsSwitchable())
		underlyingPv = floatLeg0 - soLegPv;
	else
		underlyingPv = soLegPvForSwitch - soLegPv;
	
	/// compute deal sensitivities w.r.t fwd swap rates
	ARM_GP_Vector sensitivities (sizeSched, 0.0);
	double Sshort, Slong, floatLeg, Numeraire,Sfund;
	ARM_GP_Vector invNumeraireSensis (sizeSched, 0.0);
	ARM_GP_Matrix toInvert(sizeSched, sizeSched);
		
	/// main loop on bumped forwards
	for (j=0; j<sizeSched; j++)
	{
		/// shift fwd swap rate #j
		stripper.setSwapRate(j, forwards[j] + SENSI_CMS);
		stripper.strip();

		/// (1) for spread option leg
		for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
		{
			df = stripper.df((*cpnPaydates)[i]-asOfDate);			
			sensitivities[j] -= A[i] * (df - df0[i]) / SENSI_CMS;

			if (!IsFixCoupon(i, notUsed))
			{
				ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVect[i]);
				double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(0);
				double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(0);
				ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[0];
				ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[0];
				Slong = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);

				double shortStartTime = soVanillaArg->GetSwapShortFloatStartTime()->Elt(0);
				double shortEndTime   = soVanillaArg->GetSwapShortFloatEndTime()->Elt(0);
				ARM_GP_Vector* shortPayTimes    = soVanillaArg->GetSwapShortFixPayTimes()[0];
				ARM_GP_Vector* shortPayPeriods  = soVanillaArg->GetSwapShortFixPayPeriods()[0];
				Sshort = stripper.swapRate(shortStartTime, shortEndTime, *shortPayTimes, *shortPayPeriods);

				sensitivities[j] -= Blong[i]  * (Slong  - Slong0[i])  / SENSI_CMS;
				sensitivities[j] -= Bshort[i] * (Sshort - Sshort0[i]) / SENSI_CMS;
			}

			if(IsSwitchable())
			{
				sensitivities[j] += Afund[i] * (df - df0[i]) / SENSI_CMS;
				ARM_VanillaSpreadOptionArg* soVanillaArg = (ARM_VanillaSpreadOptionArg*)(&*itsVanillaArgSOVectForSwitch[i]);
				double longStartTime = soVanillaArg->GetSwapLongFloatStartTime()->Elt(0);
				double longEndTime   = soVanillaArg->GetSwapLongFloatEndTime()->Elt(0);
				ARM_GP_Vector* longPayTimes    = soVanillaArg->GetSwapLongFixPayTimes()[0];
				ARM_GP_Vector* longPayPeriods  = soVanillaArg->GetSwapLongFixPayPeriods()[0];
				Sfund = stripper.swapRate(longStartTime, longEndTime, *longPayTimes, *longPayPeriods);

				sensitivities[j] += Bfund[i]  * (Sfund  - Sfund0[i])  / SENSI_CMS;
			}
		}

		/// take fee into account
		df = stripper.df(ExerciseDate-asOfDate);
		if ( GetPayRec()==K_RCV )
			sensitivities[j] += itsvExerFees[exerIdx] * (df - dfFee0);
		else
			sensitivities[j] -= itsvExerFees[exerIdx] * (df - dfFee0);
		

		if (!IsSwitchable())
		{
			/// (2) for funding leg + margin
			floatLeg = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargin,  &fundingLeverage);

			sensitivities[j] += (floatLeg - floatLeg0) / SENSI_CMS;
		}


		/// (3) for 1/Numeraire
		Numeraire = 0.0;
		for (k=0; k<sizeSched; k++)
		{	
			df = stripper.df((*endDates)[k]-asOfDate);
			Numeraire += (*interestTerms)[k] * (fixNotios[k]/fixNotional) * df;
		}

		for (k=0; k<sizeSched; k++)
		{
			toInvert(j,k)  =  (stripper.df((*startDates)[k]-asOfDate) -  stripper.df((*endDates)[k]-asOfDate)) / Numeraire;
			toInvert(j,k) -= NumerizedLiborFlow0[k];
			toInvert(j,k) /=  SENSI_CMS;
		}
		
		invNumeraireSensis[j] = (1./Numeraire - 1./Numeraire0) / SENSI_CMS;

		/// RAZ
		stripper.setSwapRate(j, forwards[j]);
	}

	/// Reset stripper
	stripper.strip();

	///
	/// build matrix to invert
	/// for the moment its contains the sensis of the Underlying w.r.t swap rates
	/// now, we want it to contain the sensi of Underlying/Numeraire w.r.t. swap rates
	for (j=0; j<sizeSched; j++)
	{
		sensitivities[j] *= 1./Numeraire0;
		sensitivities[j] += underlyingPv * invNumeraireSensis[j];
	}
	
	/// invert matrix
	LinSolve(&toInvert, &sensitivities);
	ARM_Vector floatNotios (sizeSched);
	for (j=0; j<sizeSched; j++)
		floatNotios[j] = sensitivities[j];
		
	/// compute strike
	double strike;
	if (isUpdateStrike)
	{
		double varFloatLeg = 0.0;

		for (j=0; j<sizeSched; j++)
		{			
			varFloatLeg += floatNotios[j] * (stripper.df((*startDates)[j]-asOfDate) - stripper.df((*endDates)[j]-asOfDate));
		}

		if (calibStrikeType == ATM)
			underlyingPv = 0.0;

		strike = (varFloatLeg - underlyingPv) / (Numeraire0 * fixNotional);
	}

	
	///  vector need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue floatNotionalRefVal ((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)floatNotios.Clone());
	floatNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);

	///  vectors need to be cloned (deleted in ARM_ReferenceValue destructor)
	ARM_ReferenceValue fixNotionalRefVal((ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)fixNotios.Clone());
	fixNotionalRefVal.SetCalcMethod(K_STEPUP_RIGHT);

		
	// index type: same frequency as fixed leg
	
	ARM_IRIndex irIndex (indexType, stdFixFreq, stdFixFreq, ccy);
	irIndex.SetTerm(stdFixFreq);
	irIndex.SetYearTerm(1.0/stdFixFreq);

	ARM_SwapLeg armFloatLeg( startDate, 
							 endDate, 
							 &irIndex, 
							 K_PAY,
							 0.0, 
							 K_SHORTSTART, 
							 K_COMP_PROP,
							 ccy,
							 ccy->GetLiborIndexDayCount());

	/// floatNotionalRefVal is cloned in SetAmount
	armFloatLeg.SetAmount(&floatNotionalRefVal);


	// create fixed leg with good strike
	ARM_FixLeg armFixLeg ( startDate,
						   endDate, 
						   strike * 100.0,
						   K_RCV, 
						   stdFixFreq,
						   stdFixDayCount,
						   K_COMP_PROP,
						   K_ARREARS, 
						   K_ADJUSTED,
						   K_SHORTSTART,
						   ccy);
	
	
	/// fixNotionalRefVal is cloned in SetAmount
	armFixLeg.SetAmount(&fixNotionalRefVal);
	
	
	/// create swap
	ARM_Swap swap(&armFixLeg, &armFloatLeg);
	ARM_Date expiryDate (ExerciseDate);
			
	/// create swaptionj
	ARM_Swaption* swaption = new ARM_Swaption((&swap), K_RCV,K_EUROPEAN, strike * 100.0, expiryDate);


	///------------------------------------------------------------
	/// compute associated forward start swaption if necessary
	///------------------------------------------------------------
	if (exerIdx > 0)
	{
		for (int prevExerIdx=exerIdx-1; prevExerIdx>=0; prevExerIdx--)
		{
			if (itsvIsExerDate[prevExerIdx])
				break;

			if (prevExerIdx==0)
			{
				fwdStartSwaption = NULL;
				return swaption;
			}
		}


		/// find start date
		ARM_Date prevExerDate = exerDates->Elt(prevExerIdx);
		prevExerDate -= 10; /// -10 days for security
		ARM_Date newStartDate(startDate);

		while (prevExerDate<newStartDate) 
			newStartDate.AddMonths(-3);

		newStartDate.AddMonths(3);

		/*
		if (newStartDate<prevExerDate)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : unresolved date problem" );
		if (newStartDate.GetJulian()-prevExerDate.GetJulian()>10)
			ARM_THROW( ERR_INVALID_ARGUMENT, "CreateVarNotionalSwaptionAtExer : unresolved date problem" );
		*/

		/// build new ref val
		ARM_Vector newFixAbs (sizeSched+1), newFixNotios(sizeSched+1), newFloatNotios(sizeSched+1);
		newFixAbs[0]		= asOfDate;
		newFixNotios[0]		= 0.0;
		newFloatNotios[0]	= 0.0;
		for (j=0; j<sizeSched; j++)
		{
			newFixAbs[j+1]		= fixAbs[j];
			newFixNotios[j+1]	= fixNotios[j];
			newFloatNotios[j+1] = floatNotios[j];
		}

		///  vectors need to be cloned (deleted in ARM_ReferenceValue destructor)
		ARM_ReferenceValue newFixNotionalRefVal ( (ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)newFixNotios.Clone());
		newFixNotionalRefVal.SetCalcMethod(K_STEPUP_LEFT);
		ARM_ReferenceValue newFloatNotionalRefVal ( (ARM_Vector*)fixAbs.Clone(), (ARM_Vector*)newFloatNotios.Clone());
		newFloatNotionalRefVal.SetCalcMethod(K_STEPUP_LEFT);

		/// build new legs
		ARM_SwapLeg newArmFloatLeg( newStartDate, 
							 endDate, 
							 &irIndex, 
							 K_PAY,
							 0.0, 
							 K_SHORTSTART, 
							 K_COMP_PROP,
							 ccy,
							 ccy->GetLiborIndexDayCount());
	
		/// ref value is cloned in SetAmount
		newArmFloatLeg.SetAmount(&newFloatNotionalRefVal);


		// create fixed leg with good strike
		ARM_FixLeg newArmFixLeg ( newStartDate,
							   endDate, 
							   strike * 100.0,
							   K_RCV, 
							   stdFixFreq,
							   stdFixDayCount,
							   K_COMP_PROP,
							   K_ARREARS, 
							   K_ADJUSTED,
							   K_SHORTSTART,
							   ccy);

		/// ref value is cloned in SetAmount
		newArmFixLeg.SetAmount(&newFixNotionalRefVal);


		/// create swap
		ARM_Swap newSwap(&newArmFixLeg, &newArmFloatLeg);

	
		/// associated expiry date
		ARM_Date newExpiryDate(newStartDate);
		newExpiryDate.GapBusinessDay(-spotDays,resetCalendar);
		
		/// create swaption
		fwdStartSwaption = new ARM_Swaption((&newSwap), K_RCV,K_EUROPEAN, strike * 100.0, newExpiryDate);

		return swaption;
	}
	else
	{
		fwdStartSwaption = NULL;
		return swaption;
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_LocalCSOCalculator
///	Routine: ComputeResidualUnderlying
///	Returns: void
///	Action : compute PV of the resiudal underlying
///			 at current the exercise date
/////////////////////////////////////////////////////////////////
double ARM_LocalCSOCalculator::ComputeResidualUnderlying(int exerIdx, double& floatLegPv, double& soLegPv)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
		
	ARM_GP_Vector* exerDates = itsExerciseDateStrip->GetResetDates();	
	double ExerciseDate = (*exerDates)[exerIdx];	
	
	/// Restore coupon datas
	ARM_GP_Vector* cpnCoverages = GetStructDateStrip()->GetInterestTerms();
	ARM_GP_Vector* cpnPaydates  = GetStructDateStrip()->GetPaymentDates();

	/// Compute exotic coupon leg
	size_t i;
	double df,soLegPvForSwitch=0.0;
	soLegPv=0.0;
	for (i = itsvCpnIndex[exerIdx]; i<itsCpnSize; i++)
	{
		df			= curve->DiscountPrice(((*cpnPaydates)[i]-asOfDate)/K_YEAR_LEN);
		soLegPv		+= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsCpnValue[i] * df;
		if (IsSwitchable())
			soLegPvForSwitch	+= (*cpnCoverages)[i] * itsvCpnNominal[i] * itsCpnValueForSwitch[i] * df;
	}

	/// take exercise fees into account
	if ( GetPayRec()==K_RCV )
		soLegPv -= itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);
	else
		soLegPv += itsvExerFees[exerIdx] * curve->DiscountPrice((ExerciseDate-asOfDate)/K_YEAR_LEN);

	/// Compute floating and fixed part of the funding leg
	double t,it,dfStart,dfEnd;
	floatLegPv = 0.0;
	for (i=itsvFundIndex[exerIdx]; i<itsFundSize; i++)
	{
		t = (*(GetFundDateStrip()->GetFlowStartDates()))[i];
		dfStart = curve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);
		t = (*(GetFundDateStrip()->GetFlowEndDates()))[i];
		dfEnd = curve->DiscountPrice((t-asOfDate)/K_YEAR_LEN);

		it = (*(GetFundDateStrip()->GetInterestTerms()))[i];

		if (!IsSwitchable()) 
			soLegPv -= it * itsvFundSpread[i]*dfEnd * itsvFundNominal[i];
		floatLegPv += (dfStart-dfEnd) * itsvFundNominal[i];
	}
	
	if (IsSwitchable())
		soLegPv = soLegPv - soLegPvForSwitch + floatLegPv;
		
	return GetPayRec() * (soLegPv - floatLegPv);
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
