/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file tarncalculator.cpp
 *
 *  \brief file for the calculator for TARN
 *
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date March 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/tarncalculator.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gptrigomatrix.h"

/// gpinfra
#include "gpinfra/argconvdefault.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numeraire.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/gramfunctorargdict.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/modelnamemap.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"
#include "gpcalib/modelparamsfactory.h"

/// gpmodels
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmdiag.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/ForwardMarginIR.h"
#include "gpmodels/ForwardMarginBasis.h"
#include "gpmodels/forwardforex.h"
#include "gpmodels/hybridbasisfwdir.h"
#include "gpmodels/smiledfrm.h"
#include "gpmodels/markovfunctional.h"
#include "gpmodels/modelparamsMF.h"

/// gpnumlib
#include "gpnumlib/argconvdefault.h"
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"

/// gpnummethods
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"

#include "gpclosedforms/sabrbdiff1.h"

/// kernel
//#include <inst/armdigital.h>
//#include <inst/swaption.h>
//#include <inst/portfolio.h>
//#include <inst/forex.h>
#include <glob/paramview.h>
//#include <mod/bssmiled.h>

/// STL
#include <iomanip> /// for setprecision()
#include <memory>

CC_BEGIN_NAMESPACE( ARM )

// SFRM default factors number
const int SFRM_NB_FACTORS			= 2;
const int SFRM_VOL_TYPE				= K_DIAG;

// MC default pricing values
const double MC_NB_INTER_STEPS		= 1;

/// SFRM sigma range [1bp,10000bp] with a 500bp default value
const double SIGMA_LOWER_BOUND      = 0.0001;
const double SIGMA_UPPER_BOUND      = 1.0;
const double SIGMA_DEFAULT_VALUE    = 0.20;

/// SFRM beta range [0%,130%] with a 100% default value
const double BETA_DEFAULT_VALUE    = 1;
const double BETA_LOWER_BOUND	   = 0.01;
const double BETA_UPPER_BOUND	   = 1.3;

/// SFRM MRS range [-15%,50%] with a 2% default value
const double MRS_LOWER_BOUND        = -0.15;
const double MRS_UPPER_BOUND        = 0.15;
const double MRS_DEFAULT_VALUE      = 0.02;

/// Defaults for SBGM
const double HUMP_DEFAULT_VALUE       = 1.;

const double BETACORREL_DEFAULT_VALUE = 0.;
const double BETACORREL_LOWER_BOUND = 0.;
const double BETACORREL_UPPER_BOUND = 1.;

// Theta (Correl)
const double THETA_DEFAULT_VALUE = 0.0;

// Minimum Strike Value for Cap Floor calibration
const double MIN_CF_STRIKE			= 0.0025;

/// 0.001bp of vega to be selected in portfolio for volatility bootstrapping 
const double VEGA_MIN_TO_SELECT=0.0000001;
const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

const double CF_DEFAULT_PRICE=1.0e+100;
const double CF_DEFAULT_WEIGHT=1.0;

const double DIGITAL_DEFAULT_WEIGHT=0.000001;

/// To indicate no link to any additional standard swaps for
/// equivalent strike/vol computation (diagonal swaption generation)
const int NO_OTHER_STDSWAPS=-1;

/// Equivalent datas computed for diagonal swaption
const unsigned int OSW_SWAP_RATE    = 0;
const unsigned int OSW_TARGET_PRICE = 1;
const unsigned int OSW_TARGET_VEGA  = 2;
const unsigned int OSW_NB_EQUIVDATA = 3;

/// Reference schedules for TARN date structure
const unsigned int RF_CPN_SCHED =0;
const unsigned int NB_TARN_SCHED =1;

// Digital spreads
const double DigitalSpread1 = -0.1;
const double DigitalSpread2 = +0.1;

// Digital smoothing spread
const double DigitalSmoothingSpread = 0.0001;

// Nb iterations max for strike calibration 
const int MaxIterCalibStrike = 10000;

// Nb iterations for control variable coeff calculation
const int ControlVariableNbIter = 500;

// Tolerance to replace swaption reset dates with model dates
const int SwaptionDateTolerance = 10;


// default PDE params for HK/SBGM
const int		PDE_NB_STEPS		= 1000;
const double	PDE_NB_STDEV		= 6;
const int		PDE_GRIDSIZE		= 901;
const int		SABR_GRIDSIZE		= 901;

// default params for SBGM
const int SBGM_NB_FACTORS				= 4;
const double SBGM_MC_NB_INTER_STEPS		= 0;


const string ARM_TARNCalculator::TARNColNamesTable [] =
{
    "EventDate",
	"NextEventDate",
    "StartDate",
    "EndDate",
    "PayDate",
    "IT",
    "IndexStartDate",
	"IndexEndDate",
	"ITStdSwap",
    "Strike",
    "Leverage",
	"CouponMin",
	"CouponMax",
	"LevPrev",
	"LifeTimeCapTarget",
	"LifeTimeFloorTarget",
    "Nominal",
	"Fees",
    "FundingSpread",
    "FundingStartDate",
    "FundingEndDate",
    "FundingAnnuity",
    "FundingVarFlow",
    "FundingSpreadFlow",
	"Funding",
	"FundingPaid",
	"DFPay",
    "CpnIndex",
    "CouponWithoutIT",
    "Coupon",
	"PaidNb",
	"CouponCF",
	"PaidCoupon",
	"SumCoupons",
	"DigitalUp",
	"DigitalDown",
	"IsAlived",
	"HasTriggered",
	"EstimExerciseStrike",
	"StdFixFlow",
	"StdFixRFFlow",
	"Swap",
	"TargetCap",
	"TargetFloor",
	"ControleVariable",
	"ControleVariable2",
	"ExerciseStrike",
	"ExerciseProba",
	"ExerciseTime",
	"LifeTimeCap",
	"LastLifeTimeCap",
	"LifeTimeFloor",
	"DigitalFunding",
	"Digital",
    "TARN"
};


const int ARM_TARNCalculator::ProductToPriceColumns[] =
{
	ARM_TARNCalculator::TARN,
	ARM_TARNCalculator::Swap,
	ARM_TARNCalculator::LifeTimeCap,
	ARM_TARNCalculator::LifeTimeFloor,
	ARM_TARNCalculator::DigitalFunding,
	ARM_TARNCalculator::Digital,
	ARM_TARNCalculator::FundingPaid,
	ARM_TARNCalculator::ExerciseStrike,
	ARM_TARNCalculator::ExerciseProba,
	ARM_TARNCalculator::ExerciseTime,
	ARM_TARNCalculator::LastLifeTimeCap
};




/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Constructor
///	Returns: void
///	Action : contructor for callable TARN SB
/////////////////////////////////////////////////////////////////

ARM_TARNCalculator::ARM_TARNCalculator(const ARM_Date& startDate,
    const ARM_Date& endDate,
    const ARM_Curve& strike,
    int payRec,
    int cpnDayCount,
    int cpnFreq,
    int cpnTiming,
    const string& cpnIndexTerm,
    int cpnIndexDayCount,
    const string& cpnResetCal,
    const string& cpnPayCal,
    int cpnResetGap,
	int intRule,
	int stubRule,
    const ARM_Curve& leverage,
	const ARM_Curve& cpnMin,
	const ARM_Curve& cpnMax,
	const ARM_Curve& levPrev,
	double lifeTimeCapTarget,
	bool globalCapFlag,
    double lifeTimeFloorTarget,
    const ARM_Curve& fundSpread,
    int fundFreq,
    int fundDayCount,
    const ARM_Curve& nominal,
	const ARM_Curve& fees,
	const ARM_Curve& fundNominal,
	const std::vector<double>& nbIterations,
	ARM_ModelType modelType,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool oswCalibFlag,
	bool controlVariableFlag,
	bool digitalSmoothingFlag,
	bool smiledFRMRescalling,
	const string& genType1,
	const string& genType2,
	int firstNbTimes,
	int firstNbDims,
	const string& pathScheme,
	const string& pathOrder,
	bool antithetic,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
    const ARM_MarketData_ManagerRep& mktDataManager)
 : ARM_GenCalculator(mktDataManager),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsStrike(strike),
	itsPayRec(payRec),
	itsCpnDayCount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnTiming(cpnTiming),
	itsCpnIndexTerm(cpnIndexTerm),
	itsCpnIndexDayCount(cpnIndexDayCount),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsCpnResetGap(cpnResetGap),
	itsIntRule(intRule),
	itsStubRule(stubRule),
	itsLeverage(leverage),
	itsCpnMin(cpnMin),
	itsCpnMax(cpnMax),
	itsLevPrev(levPrev),
	itsLifeTimeCapTarget(lifeTimeCapTarget),
	itsGlobalCapFlag(globalCapFlag),
	itsLifeTimeFloorTarget(lifeTimeFloorTarget),
	itsFundSpread(fundSpread),
	itsFundFreq(fundFreq),
	itsFundDayCount(fundDayCount),
	itsNominal(nominal),
	itsFees(fees),
	itsFundNominal(fundNominal),
	itsNbIterations(nbIterations),
	itsModelType(modelType),
	itsCapFloorCalibMode(capFloorCalibMode),
	itsDigitalCalibMode(digitalCalibMode),
	itsOSWCalibFlag(oswCalibFlag),
	itsControlVariableFlag(controlVariableFlag),
	itsDigitalSmoothingFlag(digitalSmoothingFlag),
	itsSmiledFRMRescallingFlag(smiledFRMRescalling),
	itsGenType1(genType1),
	itsGenType2(genType2),
	itsFirstNbTimes(firstNbTimes),
	itsFirstNbDims(firstNbDims),
	itsPathScheme(pathScheme),
	itsPathOrder(pathOrder),
	itsAntithetic(antithetic),
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsFirstEventIdx(0),
	itsFirstRFIdx(0),
	itsTARNPrice(0.0),
	itsTARNStdDev(0.0),
	itsTARNRA(NULL),
	itsSwapPrice(0.0),
	itsSwapStdDev(0.0),
	itsSwapRA(NULL),
	itsLifeTimeCapPrice(0.0),
	itsLifeTimeCapStdDev(0.0),
	itsLifeTimeCapRA(NULL),
	itsLifeTimeFloorPrice(0.0),
	itsLifeTimeFloorStdDev(0.0),
	itsLifeTimeFloorRA(NULL),
	itsDigitalFundingPrice(0.0),
	itsDigitalFundingStdDev(0.0),
	itsDigitalFundingRA(NULL),
	itsDigitalPrice(0.0),
	itsDigitalStdDev(0.0),
	itsDigitalRA(NULL),
	itsFundingPrice(0.0),
	itsFundingStdDev(0.0),
	itsFundingRA(NULL),
	itsExerciseStrikes(0),
	itsExerciseProbas(0),
	itsExerciseTimeAverage(0.0),
	itsCVRefPrices(NULL),
	itsFundDateStrip(0),
	itsStructDateStrip(0),
	itsvFundNominal(0),
	itvCpnsNominal(0),
	itsvInitialFundSpread(0)
{
	SetName(ARM_TARN);

	// By default there is no fees in the TARN SB.
	std::vector<double> dates(1,0.0);
	std::vector<double> values(1,0.0);
	itsFees = ARM_Curve(dates, values, new ARM::ARM_LinInterpCstExtrapolDble);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Constructor
///	Returns: void
///	Action : TARN ETK
// 
// N.B : This version is dedicated to SFRM pricing
//
/////////////////////////////////////////////////////////////////

ARM_TARNCalculator::ARM_TARNCalculator(const ARM_Date& asOfDate,
	const ARM_Date& startDate,
    const ARM_Date& endDate,
    const ARM_Curve& strike,
    int payRec,
    int cpnDayCount,
    int cpnFreq,
    int cpnTiming,
    const string& cpnIndexTerm,
    int cpnIndexDayCount,
    const string& cpnResetCal,
    const string& cpnPayCal,
    int cpnResetGap,
	int intRule,
    const ARM_Curve& leverage,
	const ARM_Curve& cpnMin,
	const ARM_Curve& cpnMax,
	double lifeTimeCapTarget,
	bool globalCapFlag,
    double lifeTimeFloorTarget,
    const ARM_Curve& fundSpread,
    int fundFreq,
    int fundDayCount,
    const ARM_Curve& nominal,
	const ARM_Curve& fees,
	const std::vector<double>& nbIterations,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool oswCalibFlag,
	bool controlVariableFlag,
	bool digitalSmoothingFlag,
	bool smiledFRMRescallingFlag,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
    const ARM_Currency& cpnCcy,
    const ARM_Currency& fundingCcy,
    const ARM_Currency& basisCcy,
    const ARM_Currency& domesticCcy,
    const ARM_Currency& foreignCcy,
	const ARM_Curve& fundnominal,
    const ARM_StringVector *mdmKeys,
    const ARM_MarketData_ManagerRep *mktDataManager,
	bool isCustomResetFlag,
	std::vector<double>& customResetDates
	) : ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsStrike(strike),
	itsPayRec(payRec),
	itsCpnDayCount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnTiming(cpnTiming),
	itsCpnIndexTerm(cpnIndexTerm),
	itsCpnIndexDayCount(cpnIndexDayCount),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsCpnResetGap(cpnResetGap),
	itsIntRule(intRule),
	itsStubRule(K_SHORTSTART),
	itsLeverage(leverage),
	itsCpnMin(cpnMin),
	itsCpnMax(cpnMax),
	itsLifeTimeCapTarget(lifeTimeCapTarget),
	itsGlobalCapFlag(globalCapFlag),
	itsLifeTimeFloorTarget(lifeTimeFloorTarget),
	itsFundSpread(fundSpread),
	itsFundFreq(fundFreq),
	itsFundDayCount(fundDayCount),
	itsNominal(nominal),
	itsFees(fees),
	itsNbIterations(nbIterations),
	itsModelType(ARM_PricingModelType::SFRM2F),
	itsCapFloorCalibMode(capFloorCalibMode),
	itsDigitalCalibMode(digitalCalibMode),
	itsOSWCalibFlag(oswCalibFlag),
	itsControlVariableFlag(controlVariableFlag),
	itsDigitalSmoothingFlag(digitalSmoothingFlag),
	itsSmiledFRMRescallingFlag(smiledFRMRescallingFlag),
	itsGenType1("MRGK5"),
	itsGenType2("Sobol"),
	itsFirstNbTimes(0),
	itsFirstNbDims(0),
	itsPathScheme("BrownianBridge"),
	itsPathOrder("PathOrder"),
	itsAntithetic(true),
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsFirstEventIdx(0),
	itsFirstRFIdx(0),
	itsTARNPrice(0.0),
	itsTARNStdDev(0.0),
	itsTARNRA(0),
	itsSwapPrice(0.0),
	itsSwapStdDev(0.0),
	itsSwapRA(NULL),
	itsLifeTimeCapPrice(0.0),
	itsLifeTimeCapStdDev(0.0),
	itsLifeTimeCapRA(NULL),
	itsLifeTimeFloorPrice(0.0),
	itsLifeTimeFloorStdDev(0.0),
	itsLifeTimeFloorRA(NULL),
	itsDigitalFundingPrice(0.0),
	itsDigitalFundingStdDev(0.0),
	itsDigitalFundingRA(NULL),
	itsDigitalPrice(0.0),
	itsDigitalStdDev(0.0),
	itsDigitalRA(NULL),
	itsFundingPrice(0.0),
	itsFundingStdDev(0.0),
	itsFundingRA(NULL),
	itsExerciseStrikes(0),
	itsExerciseProbas(0),
	itsExerciseTimeAverage(0.0),
	itsCVRefPrices(NULL),
	itsFundNominal(fundnominal),
	itsCustomResetFlag(isCustomResetFlag),
	itsCustomResetDates(customResetDates),
	itsFundDateStrip(0),
	itsStructDateStrip(0),
	itsvFundNominal(0),
	itvCpnsNominal(0),
	itsvInitialFundSpread(0)
{
	SetName(ARM_TARN);

    SetCurrencyUnit(const_cast< ARM_Currency* >(&cpnCcy));
    SetFundingCcy(const_cast< ARM_Currency& >(fundingCcy));
    SetBasisCcy(const_cast< ARM_Currency& >( basisCcy));
    SetDomesticCcy(const_cast< ARM_Currency& > (domesticCcy));
    SetForeignCcy(const_cast< ARM_Currency& > (foreignCcy));

	// By default there is no lev prev in the TARN.
	std::vector<double> dates(1,0.0);
	std::vector<double> values(1,0.0);
	itsLevPrev = ARM_Curve(dates, values, new ARM::ARM_LinInterpCstExtrapolDble);

    /// Set keys for MDM datas access
	if ( mdmKeys == NULL || mdmKeys->size() == 0 )
	{
		string cpnCcyName(cpnCcy.GetCcyName());
	    string fundingCcyName(fundingCcy.GetCcyName());

		// create internal keys
		ARM_StringVector keys(NbKeys);
		keys[YcKey]				= YC_KEY_NAME			+ cpnCcyName;
		keys[OswModelKey]		= OSWMODEL_KEY_NAME		+ cpnCcyName;
		keys[CfModelKey]		= CFMODEL_KEY_NAME		+ cpnCcyName;
		keys[MrsKey]			= MRS_KEY_NAME			+ cpnCcyName;
		keys[BetaKey]			= BETA_KEY_NAME			+ cpnCcyName;
		keys[CorrelKey]			= CORREL_KEY_NAME		+ cpnCcyName;
		keys[HumpKey]			= HUMP_KEY_NAME			+ cpnCcyName;
		keys[BetaCorrelKey]		= BETACORREL_KEY_NAME   + cpnCcyName;
 
		if( fundingCcyName == cpnCcyName )
		{
			keys.resize(NbKeys-1);
			// no basis, no forex but keep compatibility !
			keys[FundingKey]	= keys[YcKey];
			keys[BasisKey]		= keys[YcKey];
			keys[FundBasisKey]   = keys[YcKey];

		}
		else
		{
			keys[FundingKey] = YC_KEY_NAME       + fundingCcyName;
			keys[BasisKey]   = YC_BASIS_KEY_NAME + string(basisCcy.GetCcyName());
			keys[FundBasisKey]   = YC_BASIS_KEY_NAME + string(fundingCcy.GetCcyName());
			keys[ForexKey]   = FOREX_KEY_NAME + string(foreignCcy.GetCcyName()) + "_" + string(domesticCcy.GetCcyName());
		}

        SetKeys(keys);
	}
	else
	{
		if(mdmKeys->size() < (NbKeys-5))
		{
			/// To be compatible with SBGM version, should have hump, betacorrel, fund and bs keys
			ARM_StringVector newMdMKeys(*mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[HumpKey]			= newMdMKeys[YcKey];
			newMdMKeys[BetaCorrelKey]   = newMdMKeys[YcKey];
			newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]	= newMdMKeys[YcKey];
			SetKeys(newMdMKeys);
		}
		else if(mdmKeys->size() < (NbKeys-1))
		{
			/// To be compatible with non basis version, should at least have fund and bs keys
			ARM_StringVector newMdMKeys(*mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]	= newMdMKeys[YcKey];
			SetKeys(newMdMKeys);
		}
		else
			SetKeys(*mdmKeys);
	}

	// 12M is 
	if (itsCpnIndexTerm == "12M")
	{
		itsCpnIndexTerm = "1Y";
	}
	/// Check input datas
    CheckDataAndTimeIt();

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

	//Summit customised reset dates
	if (isCustomResetFlag == true)
	{
		SetCustomResetFlag(isCustomResetFlag);
		SetCustomResetDates(customResetDates);
	}

	bool needOtherPayoff = 
		itsProductsToPrice[ExerciseStrikesPrice] 
		|| itsProductsToPrice[ExerciseProbasPrice]
		|| controlVariableFlag;


    /// Create the Generic Security
	CreateAndSetDealDescription((GetKeys()[itsCpnModelKey]), CreateColumnNames(), ARM_CstManagerPtr(NULL),false,needOtherPayoff);
	if (mktDataManager)
	{
		Update(mktDataManager);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Constructor
///	Action : builds the TARN
// 
// N.B : This constructor is dedicated to SBGM or SFRM pricing
//
///	Action : TARN Addin
/////////////////////////////////////////////////////////////////

ARM_TARNCalculator::ARM_TARNCalculator(const ARM_Date& asOfDate,
	const ARM_Date& startDate,
    const ARM_Date& endDate,
    const ARM_Curve& strike,
    int payRec,
    int cpnDayCount,
    int cpnFreq,
    int cpnTiming,
    const string& cpnIndexTerm,
    int cpnIndexDayCount,
    const string& cpnResetCal,
    const string& cpnPayCal,
    int cpnResetGap,
	int intRule,
    const ARM_Curve& leverage,
	const ARM_Curve& cpnMin,
	const ARM_Curve& cpnMax,
	double lifeTimeCapTarget,
	bool globalCapFlag,
    double lifeTimeFloorTarget,
    const ARM_Curve& fundSpread,
    int fundFreq,
    int fundDayCount,
    const ARM_Curve& nominal,
	const ARM_Curve& fees,
	const std::vector<double>& nbIterations,
	ARM_ModelType modelType,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool oswCalibFlag,
	bool controlVariableFlag,
	bool digitalSmoothingFlag,
	bool smiledFRMRescallingFlag,
	const string& genType1,
	const string& genType2,
	int firstNbTimes,
	int firstNbDims,
	const string& pathScheme,
	const string& pathOrder,
	bool antithetic,
	std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
    const ARM_Currency& cpnCcy,
    const ARM_Currency& fundingCcy,
    const ARM_Currency& basisCcy,
    const ARM_Currency& domesticCcy,
    const ARM_Currency& foreignCcy,
	const ARM_Curve& fundnominal,
    const ARM_StringVector *mdmKeys,
    const ARM_MarketData_ManagerRep *mktDataManager,
	bool isCustomResetFlag,
	std::vector<double>& customResetDates
	) : ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsStrike(strike),
	itsPayRec(payRec),
	itsCpnDayCount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnTiming(cpnTiming),
	itsCpnIndexTerm(cpnIndexTerm),
	itsCpnIndexDayCount(cpnIndexDayCount),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsCpnResetGap(cpnResetGap),
	itsIntRule(intRule),
	itsStubRule(K_SHORTSTART),
	itsLeverage(leverage),
	itsCpnMin(cpnMin),
	itsCpnMax(cpnMax),
	itsLifeTimeCapTarget(lifeTimeCapTarget),
	itsGlobalCapFlag(globalCapFlag),
	itsLifeTimeFloorTarget(lifeTimeFloorTarget),
	itsFundSpread(fundSpread),
	itsFundFreq(fundFreq),
	itsFundDayCount(fundDayCount),
	itsNominal(nominal),
	itsFees(fees),
	itsNbIterations(nbIterations),
	itsModelType(modelType),
	itsCapFloorCalibMode(capFloorCalibMode),
	itsDigitalCalibMode(digitalCalibMode),
	itsOSWCalibFlag(oswCalibFlag),
	itsControlVariableFlag(controlVariableFlag),
	itsDigitalSmoothingFlag(digitalSmoothingFlag),
	itsSmiledFRMRescallingFlag(smiledFRMRescallingFlag),
	itsGenType1(genType1),
	itsGenType2(genType2),
	itsFirstNbTimes(firstNbTimes),
	itsFirstNbDims(firstNbDims),
	itsPathScheme(pathScheme),
	itsPathOrder(pathOrder),
	itsAntithetic(antithetic),
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsFirstEventIdx(0),
	itsFirstRFIdx(0),
	itsTARNPrice(0.0),
	itsTARNStdDev(0.0),
	itsTARNRA(0),
	itsSwapPrice(0.0),
	itsSwapStdDev(0.0),
	itsSwapRA(NULL),
	itsLifeTimeCapPrice(0.0),
	itsLifeTimeCapStdDev(0.0),
	itsLifeTimeCapRA(NULL),
	itsLifeTimeFloorPrice(0.0),
	itsLifeTimeFloorStdDev(0.0),
	itsLifeTimeFloorRA(NULL),
	itsDigitalFundingPrice(0.0),
	itsDigitalFundingStdDev(0.0),
	itsDigitalFundingRA(NULL),
	itsDigitalPrice(0.0),
	itsDigitalStdDev(0.0),
	itsDigitalRA(NULL),
	itsFundingPrice(0.0),
	itsFundingStdDev(0.0),
	itsFundingRA(NULL),
	itsExerciseStrikes(0),
	itsExerciseProbas(0),
	itsExerciseTimeAverage(0.0),
	itsCVRefPrices(NULL),
	itsFundNominal(fundnominal),
	itsCustomResetFlag(isCustomResetFlag),
	itsCustomResetDates(customResetDates),
	itsFundDateStrip(0),
	itsStructDateStrip(0),
	itsvFundNominal(0),
	itvCpnsNominal(0),
	itsvInitialFundSpread(0)
{
	SetName(ARM_TARN);

    SetCurrencyUnit(const_cast< ARM_Currency* >(&cpnCcy));
    SetFundingCcy(const_cast< ARM_Currency& >(fundingCcy));
    SetBasisCcy(const_cast< ARM_Currency& >( basisCcy));
    SetDomesticCcy(const_cast< ARM_Currency& > (domesticCcy));
    SetForeignCcy(const_cast< ARM_Currency& > (foreignCcy));

	// By default there is no lev prev in the TARN.
	std::vector<double> dates(1,0.0);
	std::vector<double> values(1,0.0);
	itsLevPrev = ARM_Curve(dates, values, new ARM::ARM_LinInterpCstExtrapolDble);


    /// Set keys for MDM datas access
	if ( mdmKeys == NULL || mdmKeys->size() == 0 )
	{
		string cpnCcyName(cpnCcy.GetCcyName());
	    string fundingCcyName(fundingCcy.GetCcyName());

		// create internal keys
		ARM_StringVector keys(NbKeys);
		keys[YcKey]				= YC_KEY_NAME			+ cpnCcyName;
		keys[OswModelKey]		= OSWMODEL_KEY_NAME		+ cpnCcyName;
		keys[CfModelKey]		= CFMODEL_KEY_NAME		+ cpnCcyName;
		keys[MrsKey]			= MRS_KEY_NAME			+ cpnCcyName;
		keys[BetaKey]			= BETA_KEY_NAME			+ cpnCcyName;
		keys[CorrelKey]			= CORREL_KEY_NAME		+ cpnCcyName;
		keys[HumpKey]			= HUMP_KEY_NAME			+ cpnCcyName;
		keys[BetaCorrelKey]		= BETACORREL_KEY_NAME   + cpnCcyName;
 
		if( fundingCcyName == cpnCcyName )
		{
			keys.resize(NbKeys-1);
			// no basis, no forex but keep compatibility !
			keys[FundingKey] = keys[YcKey];
			keys[BasisKey]   = keys[YcKey];
			keys[FundBasisKey]   = keys[YcKey];
		}
		else
		{
			keys[FundingKey]	= YC_KEY_NAME       + fundingCcyName;
			keys[BasisKey]		= YC_BASIS_KEY_NAME + string(basisCcy.GetCcyName());
			keys[FundBasisKey]  = YC_BASIS_KEY_NAME + string(fundingCcy.GetCcyName());
			keys[ForexKey]		= FOREX_KEY_NAME + string(foreignCcy.GetCcyName()) + "_" + string(domesticCcy.GetCcyName());
		}

        SetKeys(keys);
	}
	else
	{
		if(mdmKeys->size() < (NbKeys-5))
		{
			/// To be compatible with SBGM version, should have hump, betacorrel, fund and bs keys
			ARM_StringVector newMdMKeys(*mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[HumpKey]			= newMdMKeys[YcKey];
			newMdMKeys[BetaCorrelKey]   = newMdMKeys[YcKey];
			newMdMKeys[ReCorrelKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]	= newMdMKeys[YcKey];
			SetKeys(newMdMKeys);
		}
		else if(mdmKeys->size() < (NbKeys-1))
		{
			/// To be compatible with non basis version, should at least have fund and bs keys
			ARM_StringVector newMdMKeys(*mdmKeys);
			newMdMKeys.resize(NbKeys-1);
			newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
			newMdMKeys[BasisKey]		= newMdMKeys[YcKey];
			newMdMKeys[FundBasisKey]	= newMdMKeys[YcKey];

			SetKeys(newMdMKeys);
		}
		else
			SetKeys(*mdmKeys);
	}

	// 12M is 
	if (itsCpnIndexTerm == "12M")
	{
		itsCpnIndexTerm = "1Y";
	}
	/// Check input datas
    CheckDataAndTimeIt();

    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

	//Summit customised reset dates
	if (isCustomResetFlag == true)
	{
		SetCustomResetFlag(isCustomResetFlag);
		SetCustomResetDates(customResetDates);
	}

	bool needOtherPayoff = 
		itsProductsToPrice[ExerciseStrikesPrice] 
		|| itsProductsToPrice[ExerciseProbasPrice]
		|| controlVariableFlag;


    /// Create the Generic Security
	CreateAndSetDealDescription("", CreateColumnNames(), ARM_CstManagerPtr(NULL),false,needOtherPayoff);
	if (mktDataManager)
	{
		Update(mktDataManager);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_TARNCalculator::ARM_TARNCalculator(const ARM_TARNCalculator& rhs)
: ARM_GenCalculator(rhs),
	itsStartDate(rhs.itsStartDate),
	itsEndDate(rhs.itsEndDate),
	itsStrike(rhs.itsStrike),
	itsPayRec(rhs.itsPayRec),
	itsCpnDayCount(rhs.itsCpnDayCount),
	itsCpnFreq(rhs.itsCpnFreq),
	itsCpnTiming(rhs.itsCpnTiming),
	itsCpnIndexTerm(rhs.itsCpnIndexTerm),
	itsCpnIndexDayCount(rhs.itsCpnIndexDayCount),
	itsCpnResetCal(rhs.itsCpnResetCal),
	itsCpnPayCal(rhs.itsCpnPayCal),
	itsCpnResetGap(rhs.itsCpnResetGap),
	itsIntRule(rhs.itsIntRule),
	itsStubRule(rhs.itsStubRule),
	itsLeverage(rhs.itsLeverage),
	itsCpnMin(rhs.itsCpnMin),
	itsCpnMax(rhs.itsCpnMax),
	itsLifeTimeCapTarget(rhs.itsLifeTimeCapTarget),
	itsGlobalCapFlag(rhs.itsGlobalCapFlag),
	itsLifeTimeFloorTarget(rhs.itsLifeTimeFloorTarget),
	itsFundSpread(rhs.itsFundSpread),
	itsFundFreq(rhs.itsFundFreq),
	itsFundDayCount(rhs.itsFundDayCount),
	itsNominal(rhs.itsNominal),
	itsFees(rhs.itsFees),
	itsNbIterations(rhs.itsNbIterations),
	itsProductsToPrice(rhs.itsProductsToPrice),
	itsHasBeenPriced(rhs.itsHasBeenPriced),
	itsModelType(rhs.itsModelType),
	itsCapFloorCalibMode(rhs.itsCapFloorCalibMode),
	itsDigitalCalibMode(rhs.itsDigitalCalibMode),
	itsOSWCalibFlag(rhs.itsOSWCalibFlag),
	itsControlVariableFlag(rhs.itsControlVariableFlag),
	itsDigitalSmoothingFlag(rhs.itsDigitalSmoothingFlag),
	itsSmiledFRMRescallingFlag(rhs.itsSmiledFRMRescallingFlag),
	itsGenType1(rhs.itsGenType1),
	itsGenType2(rhs.itsGenType2),
	itsFirstNbTimes(rhs.itsFirstNbTimes),
	itsFirstNbDims(rhs.itsFirstNbDims),
	itsPathScheme(rhs.itsPathScheme),
	itsPathOrder(rhs.itsPathOrder),
	itsAntithetic(rhs.itsAntithetic),
	itsFirstEventIdx(rhs.itsFirstEventIdx),
	itsFirstRFIdx(rhs.itsFirstRFIdx),
	itsTARNPrice(rhs.itsTARNPrice),
	itsTARNStdDev(rhs.itsTARNStdDev),
	itsTARNRA(rhs.itsTARNRA),
	itsSwapPrice(rhs.itsSwapPrice),
	itsSwapStdDev(rhs.itsSwapStdDev),
	itsSwapRA(rhs.itsSwapRA),
	itsLifeTimeCapPrice(rhs.itsLifeTimeCapPrice),
	itsLifeTimeCapStdDev(rhs.itsLifeTimeCapStdDev),
	itsLifeTimeCapRA(rhs.itsLifeTimeCapRA),
	itsLastLifeTimeCapPrice(rhs.itsLastLifeTimeCapPrice),
	itsLastLifeTimeCapStdDev(rhs.itsLastLifeTimeCapStdDev),
	itsLastLifeTimeCapRA(rhs.itsLastLifeTimeCapRA),
	itsLifeTimeFloorPrice(rhs.itsLifeTimeFloorPrice),
	itsLifeTimeFloorStdDev(rhs.itsLifeTimeFloorStdDev),
	itsLifeTimeFloorRA(rhs.itsLifeTimeFloorRA),
	itsDigitalFundingPrice(rhs.itsDigitalFundingPrice),
	itsDigitalFundingStdDev(rhs.itsDigitalFundingStdDev),
	itsDigitalFundingRA(rhs.itsDigitalFundingRA),
	itsDigitalPrice(rhs.itsDigitalPrice),
	itsDigitalStdDev(rhs.itsDigitalStdDev),
	itsDigitalRA(rhs.itsDigitalRA),
	itsFundingPrice(rhs.itsFundingPrice),
	itsFundingStdDev(rhs.itsFundingStdDev),
	itsFundingRA(rhs.itsFundingRA),
	itsExerciseStrikes(rhs.itsExerciseStrikes),
	itsExerciseProbas(rhs.itsExerciseProbas),
	itsExerciseTimeAverage(rhs.itsExerciseTimeAverage),
	itsCpnModelKey(rhs.itsCpnModelKey),
	itsFundingModelKey(rhs.itsFundingModelKey),
	itsBasisRefModelKey(rhs.itsBasisRefModelKey),
	itsCVRefPrices(rhs.itsCVRefPrices),
	itsFundNominal(rhs.itsFundNominal),
	itsCustomResetFlag(rhs.itsCustomResetFlag),
	itsCustomResetDates(rhs.itsCustomResetDates),
	itsFundDateStrip(rhs.itsFundDateStrip),
	itsStructDateStrip(rhs.itsStructDateStrip),
	itsvFundNominal(rhs.itsvFundNominal),
	itvCpnsNominal(rhs.itvCpnsNominal),
	itsvInitialFundSpread(rhs.itsvInitialFundSpread)
{
	// Deep copy of caplets and swaps
	size_t i;

	// Deep copy of caplets and swaps
	itsCaplets.resize(rhs.itsCaplets.size());
	for(i=0;i<rhs.itsCaplets.size();++i)
    {
        itsCaplets[i].first=static_cast< ARM_CapFloor* >(const_cast< ARM_TARNCalculator& >(rhs).itsCaplets[i].first->Clone());
        itsCaplets[i].second=rhs.itsCaplets[i].second;
    }

	itsStdSwaps.resize(rhs.itsStdSwaps.size());
	for(i=0;i<rhs.itsStdSwaps.size();++i)
    {
        itsStdSwaps[i].first=static_cast< ARM_Swap* >(const_cast< ARM_TARNCalculator& >(rhs).itsStdSwaps[i].first->Clone());
        itsStdSwaps[i].second=rhs.itsStdSwaps[i].second;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_TARNCalculator::~ARM_TARNCalculator()
{
	size_t i;

	for(i=0;i<itsCaplets.size();++i)
    {
        delete itsCaplets[i].first;
    }
	for(i=0;i<itsStdSwaps.size();++i)
    {
        delete itsStdSwaps[i].first;
    }

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: InitTARNForSummit
///	Returns: void
///	Action : Init Tarn for Summit construction
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::InitTARNForSummit(ARM_ZeroCurve    *zcCpn, 
										   ARM_VolCurve     *swoptVC, 
										   ARM_VolCurve     *capVC, 
										   ARM_VolLInterpol *capRo,
										   ARM_VolLInterpol *capNu,
										   ARM_VolLInterpol *capBeta,
										   ARM_VolLInterpol *swoptRo,
										   ARM_VolLInterpol *swoptNu,
										   ARM_VolLInterpol *swoptBeta,
										   ARM_ZeroCurve    *zcFund,
										   ARM_ZeroCurve    *zcCpnBasis,
										   ARM_ZeroCurve    *zcFundBasis,
										   double fxSpot,
										   int hedgeUpdate,
                                           ARM_ModelType modelType,
                                           double theBetaCorrel,
                                           double theHump,
                                           int SABRSigmaOrAlpha)
{
    itsModelType = modelType;

    ARM_MarketData_ManagerRepPtr mktDataManager = GetMktDataManager();

	ARM_Date asof;
    
    // Do it systematicly    if (mktDataManager.IsNull())
    {
       asof = zcCpn->GetAsOfDate();

       // TMP: May be this have to be changed in the future!

       ARM_Currency* cpnCcy   = zcCpn->GetCurrencyUnit();

       
 	   ARM_Currency* fundCcy  = cpnCcy;
       
       if (zcFund)
          fundCcy = zcFund->GetCurrencyUnit();

       ARM_Currency* basisCcy = fundCcy;
       
       if (zcFundBasis)
 	      basisCcy = zcFundBasis->GetCurrencyUnit();

 	   ARM_Currency* forCcy   = fundCcy;       
 	
       ARM_Currency* domCcy   = zcCpn->GetCurrencyUnit();
    
 
       SetCurrencyUnit(const_cast< ARM_Currency* >(cpnCcy));
       
       SetFundingCcy(const_cast< ARM_Currency& >(*fundCcy));
       SetBasisCcy(const_cast< ARM_Currency& >(*basisCcy));

       SetDomesticCcy(const_cast< ARM_Currency& > (*domCcy));
       SetForeignCcy(const_cast< ARM_Currency& > (*forCcy));
    }
/*
    else
    {
       asof = GetMktDataManager()->GetAsOfDate();
    }
*/
	///////////////////////
	// SWOPT Vol Model ////
	///////////////////////

	ARM_BSModel* swoptBSmod = NULL;	

    // SABR with Beta = 1: 
	if (swoptRo && swoptNu && !swoptBeta)
	{
		swoptBSmod = new ARM_BSSmiledModel(asof, 0,
										   zcCpn,
										   zcCpn,
										   swoptVC,
										   K_YIELD,
										   swoptRo,
										   swoptNu,
										   K_SABR_ARITH,
                                           NULL, // swoptBeta: beta == 1
                                           0.5,  // SabrWeight = 0.5
                                           SABRSigmaOrAlpha);
	}
	//Complete SABR :
	else if (swoptRo && swoptNu && swoptBeta)
	{
		int methodType;

		if (swoptBeta->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;

		swoptBSmod = new ARM_BSSmiledModel(asof, 0,
										   zcCpn,
										   zcCpn,
										   swoptVC,
										   K_YIELD,
										   swoptRo,
										   swoptNu,
										   methodType,
										   swoptBeta,
                                           0.5,  // SabrWeight = 0.5
                                           SABRSigmaOrAlpha);
	}
	//No SABR --> Vol Cube
	else if (!swoptRo && !swoptNu && !swoptBeta)
	{
		swoptBSmod = new ARM_BSModel(asof,
									 0.0,
									 zcCpn,
									 zcCpn,
									 swoptVC,
									 K_YIELD);

	}

	///////////////////////
	// IRG Vol Model //////
	///////////////////////

	ARM_BSModel* capBSmod = NULL;
	if (capRo && capNu && !capBeta)
	{
		capBSmod = new ARM_BSSmiledModel(asof, 0,
										 zcCpn,
										 zcCpn,
										 capVC,
										 K_YIELD,
										 capRo,
										 capNu,
										 K_SABR_ARITH,
                                         NULL, // Beta == NULL means beta == 1
                                         0.5,  // SabrWeight = 0.5
                                         SABRSigmaOrAlpha);
	}
	//Complete SABR :
	else if (capRo && capNu && capBeta)
	{
		int methodType;

		if (capBeta->IsEqualToOne())
			methodType = K_SABR_ARITH;
		else
			methodType = K_SABR_IMPLNVOL;

		capBSmod = new ARM_BSSmiledModel(asof, 0,
										 zcCpn,
										 zcCpn,
										 capVC,
										 K_YIELD,
										 capRo,
										 capNu,
										 methodType,
										 capBeta,
                                         0.5, // SabrWeight = 0.5
                                         SABRSigmaOrAlpha);
	}
	//No SABR --> Vol Cube
	else if (!capRo && !capNu && !capBeta)
	{

		capBSmod = new ARM_BSModel(asof,
		    					   0.0,
								   zcCpn,
								   zcCpn,
								   capVC,
								   K_YIELD);

	}


	// Mean Reversion
    /// Give a default value to MRS if not set.
    ARM_CurveModelParam mr;

    if (GetMktDataManager()->TestIfKeyMissing(GetKeys()[MrsKey]))
    {
	    std::vector<double>       values(1, 0.0);
	    std::vector<double>       times (1, 0.0);
		mr = ARM_CurveModelParam(	ARM_ModelParamType::MeanReversion, 
									&values, 
									&times, 
									"MRS");
    }
    else
	{
	   mr = * (ARM_CurveModelParam *)(&GetMRS());
	}

	// Beta
    /// Give a default value to beta if not set.
    ARM_CurveModelParam beta;

    if (GetMktDataManager()->TestIfKeyMissing(GetKeys()[BetaKey]))
    {
	    std::vector<double>       values(1, BETA_DEFAULT_VALUE);
	    std::vector<double>       times (1, 0.0);
	
        beta = ARM_CurveModelParam(	ARM_ModelParamType::Beta, 
									&values, 
									&times, 
									"BETA");
    }
    else
	{
		beta = * (ARM_CurveModelParam *)(&GetBeta());
	}

	// Correl
	/// Give a default value to Correl if not set.
	bool deleteCorrel = false;
	ARM_CorrelMatParam* correl = NULL;

	if ( GetMktDataManager()->TestIfKeyMissing(GetKeys()[CorrelKey]) )
    {
		deleteCorrel = true;

		// Buid a DateStrip to have business days in times
	    ARM_DateStrip resetLegSched(itsStartDate, 
									itsEndDate,
									itsCpnFreq,
									itsCpnDayCount,
									GetCurrencyUnit()->GetCcyName(),
									K_MOD_FOLLOWING,
									itsIntRule,
									K_SHORTSTART,
									itsCpnResetGap,
									itsCpnFreq,
									GETDEFAULTVALUE,
									"DEFAULT",
									itsCpnTiming);

		// Conversion from ARM_Vector To ARM_RealVector
		std::vector<double>& times = NULL;
		times = resetLegSched.GetResetDates();
		*times -= asof.GetJulian(); // times will be deleted in DateStrip destructor

		// Conversion from ARM_GP_Matrix To ARM_RealVector
		int n = times->size();
		ARM_GP_Matrix* trigoMatrix = TrigoMatrix(n,THETA_DEFAULT_VALUE);
        std::vector<double>& values = new std::vector<double>(n*n);

		size_t i = 0,j = 0;
		for (i = 0; i < n; ++i)
			for (j = 0; j < n; ++j)
				(*values)[i*n+j] = (*trigoMatrix)(i,j);
        

        correl = static_cast<ARM_CorrelMatParam*>(ARM_ModelParamFactory.Instance()->CreateModelParam(
			                                      ARM_ModelParamType::Correlation,
                                                  values,
                                                  times));

		if (trigoMatrix)
		{
			delete trigoMatrix;
			trigoMatrix = NULL;
		}
	
        if (values)
		{
			delete values;
			values = NULL;
		}
    }
    else
    {
	    correl = (ARM_CorrelMatParam *)(&GetCorrel());
    }

    //*------> create mktData Manager & fill it with curves
	vector <ARM_Object*> marketDatas;

	marketDatas.push_back(zcCpn);
	marketDatas.push_back(swoptBSmod);
	marketDatas.push_back(capBSmod);
	marketDatas.push_back(&mr);
	marketDatas.push_back(&beta);
	marketDatas.push_back(correl);

    ARM_CurveModelParam hump;
    ARM_CurveModelParam betaCorrel;

    if ( modelType == ARM_PricingModelType::SBGM )
    {
        // Retrieve the Hump

        if ( theHump != -10000.0 )
        {
           std::vector<double>       values(1, theHump);
	       std::vector<double>       times (1, 0.0);

		   hump = ARM_CurveModelParam(ARM_ModelParamType::Hump, 
								      &values, 
								      &times, 
								      "HUMP");
        }
        else
        {
           if (GetMktDataManager()->TestIfKeyMissing(GetKeys()[HumpKey]))
           {
	          std::vector<double>       values(1, HUMP_DEFAULT_VALUE);
	          std::vector<double>       times (1, 0.0);

		      hump = ARM_CurveModelParam(ARM_ModelParamType::Hump, 
								         &values, 
								         &times, 
								         "HUMP");
           }
           else
           {
              hump = *(static_cast< ARM_CurveModelParam* >(GetMktDataManager()->GetData(GetKeys()[HumpKey])));
           }
        }
        
        // Beta Correl
       
        if ( theBetaCorrel != -10000.0 )
        {
	       std::vector<double>       values(1, theBetaCorrel);
	       std::vector<double>       times (1, 0.0);

		   betaCorrel = ARM_CurveModelParam(ARM_ModelParamType::BetaCorrelation, 
									        &values, 
									        &times, 
									        "BETACORREL");
        }
        else
        {
            if (GetMktDataManager()->TestIfKeyMissing(GetKeys()[BetaCorrelKey]))
            {
	            std::vector<double>       values(1, BETACORREL_DEFAULT_VALUE);
	            std::vector<double>       times (1, 0.0);

		        betaCorrel = ARM_CurveModelParam(ARM_ModelParamType::BetaCorrelation, 
									             &values, 
									             &times, 
									             "BETACORREL");
            }
            else
            {
               betaCorrel = *(static_cast< ARM_CurveModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaCorrelKey])));
            }
        }

        // Add Mkt datas relatively to SBGM

        marketDatas.push_back(&hump);
        marketDatas.push_back(&betaCorrel);
    }

	// Add funding, basis and forex for Basis Tarn
    
    ARM_ZeroCurve* zcForDomBS = NULL;
	
    if (IsBasis())
	{
		ARM_Forex forex(&GetForeignCcy(), &GetDomesticCcy(), fxSpot);

		
		ARM_CRV_TERMS psMatu;

		zcForDomBS = GenerateTwoCcyBSAdjusted(zcFund, zcFundBasis, 
                                              zcCpn, zcCpnBasis,
                                              0, psMatu);

		marketDatas.push_back(zcFund);
		marketDatas.push_back(zcForDomBS);
		marketDatas.push_back(&forex);
    }

	// Update Mkt data
    if (hedgeUpdate) // Just update correctly the calculator (in Hedge context)
    {
       Update(marketDatas);
    }
    else
    {
       Init(marketDatas); // Creation of Market data manager
    }

    if (zcForDomBS)
	   delete zcForDomBS;
	zcForDomBS = NULL;

	if (swoptBSmod)
	{
	   delete swoptBSmod;
	   swoptBSmod = NULL;
	}

	if (capBSmod)
	{
	   delete capBSmod;
	   capBSmod = NULL;
	}

	if (deleteCorrel)
	{
	   delete correl;
	   correl = NULL;
	}
}	


					   
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Set???ToPrice() and Is???ToPrice()
///	Returns: void/boolean
///	Action : Set the flag to know which product to price.
///          Test the product currently able to be priced
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::SetTARNToPrice(bool toPrice)					{itsProductsToPrice[TARNPrice]				=   toPrice;}
bool ARM_TARNCalculator::IsTARNToPrice() const							{return itsProductsToPrice[TARNPrice];}

void ARM_TARNCalculator::SetSwapToPrice(bool toPrice)					{itsProductsToPrice[SwapPrice]				=   toPrice;}
bool ARM_TARNCalculator::IsSwapToPrice() const							{return itsProductsToPrice[SwapPrice];}

void ARM_TARNCalculator::SetLifeTimeCapToPrice(bool toPrice)			{itsProductsToPrice[LifeTimeCapPrice]	=   toPrice;}
bool ARM_TARNCalculator::IsLifeTimeCapToPrice() const					{return itsProductsToPrice[LifeTimeCapPrice];}

void ARM_TARNCalculator::SetLifeTimeFloorToPrice(bool toPrice)			{itsProductsToPrice[LifeTimeFloorPrice]		=   toPrice;}
bool ARM_TARNCalculator::IsLifeTimeFloorToPrice() const					{return itsProductsToPrice[LifeTimeFloorPrice];}

void ARM_TARNCalculator::SetDigitalFundingToPrice(bool toPrice)			{itsProductsToPrice[DigitalFundingPrice]		=	toPrice;}
bool ARM_TARNCalculator::IsDigitalFundingToPrice() const				{return itsProductsToPrice[DigitalFundingPrice];}

void ARM_TARNCalculator::SetDigitalToPrice(bool toPrice)				{itsProductsToPrice[DigitalPrice]		=	toPrice;}
bool ARM_TARNCalculator::IsDigitalToPrice() const						{return itsProductsToPrice[DigitalPrice];}

void ARM_TARNCalculator::SetFundingToPrice(bool toPrice)				{itsProductsToPrice[FundingPrice]			=	toPrice;}
bool ARM_TARNCalculator::IsFundingToPrice() const						{return itsProductsToPrice[FundingPrice];}

void ARM_TARNCalculator::SetExerciseStrikesToPrice(bool toPrice) 		{itsProductsToPrice[ExerciseStrikesPrice]			=	toPrice;}
bool ARM_TARNCalculator::IsExerciseStrikesToPrice() const				{return itsProductsToPrice[ExerciseStrikesPrice];}

void ARM_TARNCalculator::SetExercisesProbasToPrice(bool toPrice) 		{itsProductsToPrice[ExerciseProbasPrice]			=	toPrice;}
bool ARM_TARNCalculator::IsExercisesProbasToPrice() const				{return itsProductsToPrice[ExerciseProbasPrice];}
	
void ARM_TARNCalculator::SetExercisesTimeAverageToPrice(bool toPrice)	{itsProductsToPrice[ExerciseTimeAverage]			=	toPrice;}
bool ARM_TARNCalculator::IsExercisesTimeAverageToPrice() const			{return itsProductsToPrice[ExerciseTimeAverage];}

bool ARM_TARNCalculator::IsResetPricing() const							{return itsHasBeenPriced;}
void ARM_TARNCalculator::SetHasToPrice(bool toReset)					{itsHasBeenPriced = toReset;}

//Set the product to price. Flags is a Y/N string vector.
//Used for a Summit tarn object where only the Tarn price is set.
void ARM_TARNCalculator::SetProductToPrice(vector<string>& flags)
{
    //TARN Price
	if(flags[0] == "Y")
	    SetTARNToPrice(true);
	else if (flags[0] == "N")
	    SetTARNToPrice(false);

	//SWAP Price
	if(flags[1] == "Y")
        SetSwapToPrice(true);
	else if (flags[1] == "N")
        SetSwapToPrice(false);

	//Life Time Cap Price
	if(flags[2] == "Y")
        SetLifeTimeCapToPrice(true);
	else if (flags[2] == "N")
        SetLifeTimeCapToPrice(false);

	//Life Time Floor Price
	if(flags[3] == "Y")
        SetLifeTimeFloorToPrice(true);
	else if (flags[3] == "N")
        SetLifeTimeFloorToPrice(false);
	
	//Digital Funding Price
	if(flags[4] == "Y")
        SetDigitalFundingToPrice(true);
	else if (flags[4] == "N")
        SetDigitalFundingToPrice(false);

	//Digital Price
    if(flags[5] == "Y")
        SetDigitalToPrice(true);
	else if (flags[5] == "N")
        SetDigitalToPrice(false);
	
	//Funding Price
    if(flags[6] == "Y")
        SetFundingToPrice(true);
	else if (flags[6] == "N")
		SetFundingToPrice(false);

	//Exercise Strike Price
    if(flags[7] == "Y")
        SetExerciseStrikesToPrice(true);
	else if (flags[7] == "N")
		SetExerciseStrikesToPrice(false);

	//Exercise Probas Price
    if(flags[8] == "Y")
        SetExercisesProbasToPrice(true);
	else if (flags[8] == "N")
        SetExercisesProbasToPrice(false);

	//Exercise Time Average
    if(flags[9] == "Y")
        SetExercisesTimeAverageToPrice(true);
	else if (flags[9] == "N")   
		SetExercisesTimeAverageToPrice(false);
    
	return;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: SetModelKeys
///	Returns: void
///	Action : set the coupon, the funding & the basis reference
///          model names alias
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::SetModelKeys()
{
	// Get the coupon & funding curve ccy
    string cpnCcyName(GetCurrencyUnit()->GetCcyName());
    string fundingCcyName(GetFundingCcy().GetCcyName());
    string basisCcyName(GetBasisCcy().GetCcyName());

    if(cpnCcyName == fundingCcyName)
    {
        // no basis effect
        itsBasisRefModelKey = YcKey;
        itsCpnModelKey      = YcKey;
        itsFundingModelKey  = YcKey;
    }
    else
    {
        if(cpnCcyName == basisCcyName)
        {
            itsBasisRefModelKey     = YcKey;        // Cpn forecast = YcKey (diffused)
            itsCpnModelKey          = BasisKey;     // Cpn discount = BasisKey
            itsFundingModelKey      = FundingKey;   // Funding forecast & discount = FundingKey
        }
        else
        {
            itsCpnModelKey          = YcKey;        // Cpn forecast & discount = YcKey (diffused)
            itsBasisRefModelKey     = FundingKey;   // Funding forecast = FundingKey
            itsFundingModelKey      = BasisKey;     // Funding discount = BasisKey
        }
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_TARNCalculator::GetMRS() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[MrsKey])));
}

void ARM_TARNCalculator::SetMRS(ARM_ModelParam* mrsParam)
{
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an MRS Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[MrsKey],static_cast< ARM_Object* >(mrsParam));
}

// Hump setting

void ARM_TARNCalculator::SetHump(ARM_ModelParam* humpParam)
{
    if (!humpParam || humpParam->GetType() != ARM_ModelParamType::Hump)
	   ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : a Hump Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[HumpKey],static_cast< ARM_Object* >(humpParam));
}


// BetaCorrel Setting
void ARM_TARNCalculator::SetBetaCorrel(ARM_ModelParam* betaCorrelParam)
{
    if (!betaCorrelParam || betaCorrelParam->GetType() != ARM_ModelParamType::BetaCorrelation)
	   ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : a Beta correl Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaCorrelKey],static_cast< ARM_Object* >(betaCorrelParam));
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetBeta & SetBeta
///	Returns: 
///	Action : get & set the Beta param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_TARNCalculator::GetBeta() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaKey])));
}

void ARM_TARNCalculator::SetBeta(ARM_ModelParam* betaParam)
{
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an Beta Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaKey],static_cast< ARM_Object* >(betaParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetCorrel & SetCorrel
///	Returns: 
///	Action : get & set the Correlation param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_TARNCalculator::GetCorrel() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[CorrelKey])));
}

void ARM_TARNCalculator::SetCorrel(ARM_ModelParam* correlParam)
{
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : a Correl Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[CorrelKey],static_cast< ARM_Object* >(correlParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetCFCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for cap floor
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_TARNCalculator::GetCFCalibMethod(bool calibSwaption) const
{
    if(itsCapFloorCalibMode == NOCalib)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (calibSwaption && GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

    if(calibSwaption)
        return GetCalibMethod()->GetlinkedMethod();
    else
        return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_TARNCalculator::GetOSWCalibMethod() const
{
    if(!itsOSWCalibFlag)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetCFPortfolio
///	Returns: ARM_Portfolio
///	Action : get the cap floor calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_TARNCalculator::GetCFPortfolio() const
{
    if(itsCapFloorCalibMode == NOCalib)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (itsOSWCalibFlag && (GetCalibMethod()->GetlinkedMethod() == NULL ||
            GetCalibMethod()->GetlinkedMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL))) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method or portfolio not found");
#endif
    if(itsOSWCalibFlag)
        return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
    else
        return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: SetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of cap floors
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::SetCFPortfolio(const ARM_StdPortfolio& port)
{
    if(itsCapFloorCalibMode == NOCalib)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio settable because calibration is off");

    int pfSize=port.GetSize();
    if(pfSize<1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Not any asset in the given portfolio");


    /// Test for portfolio consistency
    for(int i=0;i<pfSize;++i)
    {
        if(port.GetAsset(i)->GetName() != ARM_CAPFLOOR)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Diagonal swaption portfolio is not only made of swaptions");
    }

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        (itsOSWCalibFlag && GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

    /// Update portfolio
    if(itsOSWCalibFlag)
        GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
    else
        GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_TARNCalculator::GetOSWPortfolio() const
{

    if(!itsOSWCalibFlag)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");
#endif

    return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: SetOSWPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of digonal 
///			 swaptions.
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::SetOSWPortfolio(const ARM_StdPortfolio& port)
{
    if(!itsOSWCalibFlag)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio settable because calibration is off");

    int pfSize=port.GetSize();
    if(pfSize<1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Not any asset in the given portfolio");


#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    /// Update portfolio
    GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CheckData & CheckMktData
///	Returns: void
///	Action : check if TARN datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CheckData()
{
	// TARN parameters checking

	if ( itsEndDate < itsStartDate )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Start Date of the deal has to be before the end date.");
	}

	if (!(itsStrike > K_NEW_DOUBLE_TOL))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Strike has to be strictly positive.");
	}

	if (!(itsLeverage > K_NEW_DOUBLE_TOL))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Leverage has to be strictly positive.");
	}

	if (itsLifeTimeCapTarget < K_NEW_DOUBLE_TOL )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Life Time Cap should be positive.");
	}

	if (itsLifeTimeFloorTarget < -K_NEW_DOUBLE_TOL )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Life Time Floor should be positive.");
	}

	if ((itsLifeTimeCapTarget-itsLifeTimeFloorTarget) < -K_NEW_DOUBLE_TOL )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Life Time Cap should be greater than Life Time Floor.");
	}

}

void ARM_TARNCalculator::CheckMktData()
{
    /// MdM datas checking
	ARM_ZeroCurve* cpnCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
    if(!cpnCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcKey] + " is expected in the Market Data Manager");
    string cpnCcy(cpnCurve->GetCurrencyUnit()->GetCcyName());


	ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
    if(!fundingCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[FundingKey] + " is expected in the Market Data Manager");
    string fundingCcy(fundingCurve->GetCurrencyUnit()->GetCcyName());


	ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
    if(!basisCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[BasisKey] + " is expected in the Market Data Manager");
    string basisCcy(basisCurve->GetCurrencyUnit()->GetCcyName());

    if(cpnCcy != fundingCcy && basisCcy != cpnCcy && basisCcy != fundingCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Basis Curve currency should be " + cpnCcy + " or " + fundingCcy);


	ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
    if(!cfBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : cap/floor B&S model for key=" + GetKeys()[CfModelKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
    if(!oswBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : swaption B&S model for key=" + GetKeys()[OswModelKey] + " is expected in the Market Data Manager");


	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Beta Param for key=" + GetKeys()[BetaKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Correl Param for key=" + GetKeys()[CorrelKey] + " is expected in the Market Data Manager");

	if(cpnCcy != fundingCcy)
    {
        if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[ForexKey]))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : a Forex is required for different coupon & funding currencies");

	    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        if(!forex)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Forex for key=" + GetKeys()[ForexKey] + " is expected in the Market Data Manager");

        string domesticCcy(forex->GetMainCurrency()->GetCcyName());
        if(domesticCcy != cpnCcy && domesticCcy != fundingCcy)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic currency should be " + cpnCcy + " or " + fundingCcy);

        string foreignCcy(forex->GetMoneyCurrency()->GetCcyName());
        if(foreignCcy != cpnCcy && foreignCcy != fundingCcy)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign currency should be " + cpnCcy + " or " + fundingCcy);
    }

	// TARN parameters checking

	ARM_Date firstResetDate(itsStartDate);

	if (itsCpnTiming == K_ARREARS)
	{
		firstResetDate.AddPeriod(itsCpnFreq);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(TARNColNamesTable)/sizeof(TARNColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = TARNColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    /// Here only exercise events are generated. At each date the coupon
    /// value is known with its analytical formula because
    /// keywords LIBOR & CAPLET were extended to handle the in arrears case

    size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
    size_t descSize = sizeof(TARNColNamesTable)/sizeof(TARNColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// EventDate description
    double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
	double lastEventDate=(*(datesStructure.GetMergeData()))[eventSize-1];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

    /// Convention conversion : K_QUATERLY => 3M for instance
    string fundMatFreq=ARM_ArgConvReverse_MatFrequency.GetString(itsFundFreq);
    string fundDayCount=ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount);
    string cpnIndexDayCount=ARM_ArgConvReverse_DayCount.GetString(itsCpnIndexDayCount);
    string cpnDayCount=ARM_ArgConvReverse_DayCount.GetString(itsCpnDayCount);

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

    /// Get the model names for coupon & funding legs descritpion
    string cpnModelName = GetKeys()[YcKey];
    string basisModelName = GetKeys()[BasisKey];

    bool isArrearsCpn = (itsCpnTiming==K_ARREARS);
	bool isFirstEvent = (eventIdx==itsFirstEventIdx);
	bool isFirstExer = (eventIdx==itsFirstRFIdx);
	
	int lastCallIdx = eventSize+1;
    bool isLastExer = (eventIdx+1>=eventSize);

    /// Get the number of passed exercises
    size_t i,nbExerBefore=0;
    for(i=0;i<eventIdx;++i)
    {
        if((*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetResetDates()))[i] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
            ++nbExerBefore;
    }

    /// Define the right flow to describe the row
    string capletResetTiming("ADV");
    string liborPayTiming("ARR");
    if(isArrearsCpn)
    {
        capletResetTiming="ARR";
        liborPayTiming="ADV";
    }


	// In the IN ARREARS case there is no funding for the last event date
	if (!isArrearsCpn || !isLastExer)
	{
		/// EventDate description
		double nextEventDate;
		if (isArrearsCpn)
		{
			nextEventDate=(*(datesStructure.GetMergeData()))[eventIdx+1];
			CC_Ostringstream nextEventDateDesc;
			nextEventDateDesc << CC_NS(std,fixed) << nextEventDate;
			rowDescVec[NextEventDate] = nextEventDateDesc.str();
			rowTypeVec[NextEventDate] = ARM_DATE_TYPE;
		}

		int fundingIdx = eventIdx;
		if (isArrearsCpn)
		{
			fundingIdx = eventIdx+1;
		}

		/// Funding Flow description
		double  flowStartDate=(*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetFlowStartDates()))[fundingIdx];
		CC_Ostringstream fundStartDesc;
		fundStartDesc << CC_NS(std,fixed) << flowStartDate;
		rowDescVec[FundingStartDate] = fundStartDesc.str();
		rowTypeVec[FundingStartDate] = ARM_DATE_TYPE;

		double  flowEndDate=(*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetFlowEndDates()))[fundingIdx];
		CC_Ostringstream fundEndDesc;
		fundEndDesc << CC_NS(std,fixed) << flowEndDate;
		rowDescVec[FundingEndDate] = fundEndDesc.str();
		rowTypeVec[FundingEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream fundSpreadDesc;
		fundSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundSpread).Interpolate(flowStartDate-asOfDate);
		rowDescVec[FundingSpread] = fundSpreadDesc.str();
		rowTypeVec[FundingSpread] = ARM_DOUBLE;

		/// Cpn Nominal description
		CC_Ostringstream nominalDesc;
		double payDate=(*(datesStructure.GetDateStrip(RF_CPN_SCHED)->GetPaymentDates()))[fundingIdx];
		nominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNominal).Interpolate(payDate-asOfDate);
		rowDescVec[Nominal] = nominalDesc.str();
		rowTypeVec[Nominal] = ARM_DOUBLE;

		CC_Ostringstream fundCommonDesc;
		fundCommonDesc << "(" << basisModelName << "," << TARNColNamesTable[FundingStartDate] << "[i],";
		fundCommonDesc << TARNColNamesTable[FundingEndDate] << "[i],";
		fundCommonDesc << fundMatFreq << ",";
		fundCommonDesc << fundDayCount << ")";
		CC_Ostringstream fundAnnuityDesc;
		fundAnnuityDesc << "ANNUITY" << fundCommonDesc.str() << "*" << TARNColNamesTable[Nominal] << "[i]";
		rowDescVec[FundingAnnuity] = fundAnnuityDesc.str();
		rowTypeVec[FundingAnnuity] = ARM_STRING;


		CC_Ostringstream fundVarFlowDesc;
		fundVarFlowDesc << "SWAPRATE" << fundCommonDesc.str() << "*" << TARNColNamesTable[FundingAnnuity] << "[i]";
		rowDescVec[FundingVarFlow] = fundVarFlowDesc.str();
		rowTypeVec[FundingVarFlow] = ARM_STRING;

		CC_Ostringstream fundSpreadFlowDesc;
		fundSpreadFlowDesc << TARNColNamesTable[FundingSpread] << "[i]*" << TARNColNamesTable[FundingAnnuity] << "[i]";
		rowDescVec[FundingSpreadFlow] = fundSpreadFlowDesc.str();
		rowTypeVec[FundingSpreadFlow] = ARM_STRING;

		CC_Ostringstream fundingDesc;
		fundingDesc << "(" << TARNColNamesTable[FundingVarFlow] << "[i]+" << TARNColNamesTable[FundingSpreadFlow]  << "[i])";

		// In the IN ARREAR case we carry the funding flow at the next event date
		if (isArrearsCpn)
		{
			fundingDesc << "/DF(" << basisModelName << "," << TARNColNamesTable[NextEventDate] << "[i])";
		}
		rowDescVec[Funding] = fundingDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		CC_Ostringstream fundingPaidDesc;
		fundingPaidDesc << TARNColNamesTable[FundingVarFlow] << "[i]+" << TARNColNamesTable[FundingSpreadFlow]  << "[i]";
		rowDescVec[FundingPaid] = fundingPaidDesc.str();
		rowTypeVec[FundingPaid] = ARM_STRING;

	}

	// In the IN ARREARS case there is no RF coupon for the first event date
	if (!isArrearsCpn || eventIdx != 0)
	{
		string fundingRef= isArrearsCpn? "[i-1]" : "[i]";
 		
		/// Flow Date descriptions	
		int curSched =  RF_CPN_SCHED;
		double flowStartDate = (*(datesStructure.GetDateStrip(curSched)->GetFlowStartDates()))[eventIdx];
		CC_Ostringstream startDateDesc;
		startDateDesc << CC_NS(std,fixed) << flowStartDate;
		rowDescVec[StartDate] = startDateDesc.str();
		rowTypeVec[StartDate] = ARM_DATE_TYPE;

		CC_Ostringstream endDateDesc;
		endDateDesc << CC_NS(std,fixed) << (*(datesStructure.GetDateStrip(curSched)->GetFlowEndDates()))[eventIdx];
		rowDescVec[EndDate] = endDateDesc.str();
		rowTypeVec[EndDate] = ARM_DATE_TYPE;

		CC_Ostringstream payDateDesc;
		double payDate=(*(datesStructure.GetDateStrip(curSched)->GetPaymentDates()))[eventIdx];
		payDateDesc << CC_NS(std,fixed) << payDate;
		rowDescVec[PayDate] = payDateDesc.str();
		rowTypeVec[PayDate] = ARM_DATE_TYPE;

		/// Cpn Nominal description
		CC_Ostringstream nominalDesc;
		nominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNominal).Interpolate(payDate-asOfDate);
		rowDescVec[Nominal] = nominalDesc.str();
		rowTypeVec[Nominal] = ARM_DOUBLE;

		/// Interest Period description 
		CC_Ostringstream itDesc;
		itDesc << CC_NS(std,fixed) << CC_NS(std,setprecision)(8) << (*(datesStructure.GetDateStrip(curSched)->GetInterestTerms()))[eventIdx];
		rowDescVec[IT] = itDesc.str();
		rowTypeVec[IT] = ARM_DOUBLE;		

		/// Fees description
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFees).Interpolate(flowStartDate-asOfDate);
		rowDescVec[Fees] = feesDesc.str();
		rowTypeVec[Fees] = ARM_DOUBLE;

		/// Strike description
		CC_Ostringstream strikeDesc;
		double strike = const_cast< ARM_Curve& >(itsStrike).Interpolate(flowStartDate-asOfDate);
		strikeDesc << CC_NS(std,fixed) << strike;
		rowDescVec[Strike] = strikeDesc.str();
		rowTypeVec[Strike] = ARM_DOUBLE;

		/// Payment DF description
		CC_Ostringstream dfPayDesc;
		dfPayDesc << "DF(" << basisModelName << "," << TARNColNamesTable[PayDate] << "[i])";
		rowDescVec[DFPay] = dfPayDesc.str();
		rowTypeVec[DFPay] = ARM_STRING;

		/// Descriptions related to the RF coupon event

		/// Cpn Index Date and interest term descriptions
		CC_Ostringstream indexStartDateDesc;
		
		ARM_Date indexStartDate(eventDate);
		
		int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
		indexStartDate.GapBusinessDay(defaultResetGap, const_cast<char *>(itsCpnResetCal.c_str()));
		
		indexStartDateDesc << CC_NS(std,fixed) << indexStartDate.GetJulian();
		rowDescVec[IndexStartDate] = indexStartDateDesc.str();
		rowTypeVec[IndexStartDate] = ARM_DATE_TYPE;


		CC_Ostringstream indexEndDateDesc;
		
		ARM_Date indexEndDate(indexStartDate);
		
		indexEndDate.AddPeriod(itsCpnIndexTerm);
		
		indexEndDateDesc << CC_NS(std,fixed) << indexEndDate.GetJulian();
		rowDescVec[IndexEndDate] = indexEndDateDesc.str();
		rowTypeVec[IndexEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream itStdSwapDesc;
		
		itStdSwapDesc << "DCF(" << TARNColNamesTable[IndexStartDate] << "[i],";
		itStdSwapDesc << TARNColNamesTable[IndexEndDate] << "[i],";
		itStdSwapDesc << fundDayCount << ")";
		
		rowDescVec[ITStdSwap] = itStdSwapDesc.str();
		rowTypeVec[ITStdSwap] = ARM_STRING;


		/// Cpn Index description
		CC_Ostringstream cpnIndexDesc;
		cpnIndexDesc << "LIBOR(" << cpnModelName << "," << TARNColNamesTable[IndexStartDate] << "[i],";
		cpnIndexDesc << itsCpnIndexTerm << "," << cpnIndexDayCount << ",";
		cpnIndexDesc << ",0," << liborPayTiming << ")";
		rowDescVec[CpnIndex] = cpnIndexDesc.str();
		rowTypeVec[CpnIndex] = ARM_STRING;

		/// Coupon Min, Max & Leverage descriptions, LifeTimeCapTarget, LifeTimeFloorTarget
		CC_Ostringstream leverageDesc;
		double leverage = const_cast< ARM_Curve& >(itsLeverage).Interpolate(flowStartDate-asOfDate);
		leverageDesc << CC_NS(std,fixed) << leverage;
		rowDescVec[Leverage] = leverageDesc.str();
		rowTypeVec[Leverage] = ARM_DOUBLE;

		CC_Ostringstream couponMinDesc;
		double couponMin = const_cast< ARM_Curve& >(itsCpnMin).Interpolate(flowStartDate-asOfDate);
		couponMinDesc << CC_NS(std,fixed) << couponMin;
		rowDescVec[CouponMin] = couponMinDesc.str();
		rowTypeVec[CouponMin] = ARM_DOUBLE;

		CC_Ostringstream couponMaxDesc;
		double couponMax = const_cast< ARM_Curve& >(itsCpnMax).Interpolate(flowStartDate-asOfDate);
		couponMaxDesc << CC_NS(std,fixed) << couponMax;
		rowDescVec[CouponMax] = couponMaxDesc.str();
		rowTypeVec[CouponMax] = ARM_DOUBLE;

		CC_Ostringstream levPrevDesc;
		double levPrev = const_cast< ARM_Curve& >(itsLevPrev).Interpolate(flowStartDate-asOfDate);
		levPrevDesc << CC_NS(std,fixed) << levPrev;
		rowDescVec[LevPrev] = levPrevDesc.str();
		rowTypeVec[LevPrev] = ARM_DOUBLE;

		CC_Ostringstream lifeTimeCapTargetDesc;
		lifeTimeCapTargetDesc << CC_NS(std,fixed) << itsLifeTimeCapTarget;
		rowDescVec[LifeTimeCapTarget] = lifeTimeCapTargetDesc.str();
		rowTypeVec[LifeTimeCapTarget] = ARM_DOUBLE;

		CC_Ostringstream lifeTimeFloorTargetDesc;
		lifeTimeFloorTargetDesc << CC_NS(std,fixed) << itsLifeTimeFloorTarget;
		rowDescVec[LifeTimeFloorTarget] = lifeTimeFloorTargetDesc.str();
		rowTypeVec[LifeTimeFloorTarget] = ARM_DOUBLE;

		/// Coupon description
		string sCouponDesc = CreateCpnDescription(isFirstExer);
		rowDescVec[CouponWithoutIT] = sCouponDesc;
		rowTypeVec[CouponWithoutIT] = ARM_STRING;

		CC_Ostringstream couponDesc;
		couponDesc << TARNColNamesTable[CouponWithoutIT] << "[i]*"<<TARNColNamesTable[IT] << "[i]";

		rowDescVec[Coupon] = couponDesc.str();
		rowTypeVec[Coupon] = ARM_STRING;

		/// Coupon Paid description
		CC_Ostringstream couponPaidDesc;
		couponPaidDesc << TARNColNamesTable[Coupon] << "[i]*" << TARNColNamesTable[DFPay] << "[i]*";
		couponPaidDesc << TARNColNamesTable[Nominal] << "[i]";

		rowDescVec[PaidCoupon] = couponPaidDesc.str();
		rowTypeVec[PaidCoupon] = ARM_STRING;

		rowDescVec[PaidNb] = "0";
		rowTypeVec[PaidNb] = ARM_STRING;
	

		/// coupon CF for the control variate
		string couponCFDesc = CreateCFCpnDdescription(isFirstExer,cpnModelName);
		rowDescVec[CouponCF] = couponCFDesc;
		rowTypeVec[CouponCF] = ARM_STRING;

		/// Standard fixed Flow (estimated strikes)
		CC_Ostringstream stdFixFlowDesc;
		
		stdFixFlowDesc << TARNColNamesTable[EstimExerciseStrike] << "[i]*";
		stdFixFlowDesc << "DF(" << basisModelName << "," << TARNColNamesTable[IndexEndDate] << "[i])*";
		stdFixFlowDesc << TARNColNamesTable[ITStdSwap] << "[i]";

		rowDescVec[StdFixFlow] = stdFixFlowDesc.str();
		rowTypeVec[StdFixFlow] = ARM_STRING;

		/// Standard fixed Flow (reverse floater strike)
		CC_Ostringstream stdFixRFFlowDesc;
		
		stdFixRFFlowDesc << strike/leverage << "*";
		stdFixRFFlowDesc << "DF(" << basisModelName << "," << TARNColNamesTable[IndexEndDate] << "[i])*";
		stdFixRFFlowDesc << TARNColNamesTable[ITStdSwap] << "[i]";

		rowDescVec[StdFixRFFlow] = stdFixRFFlowDesc.str();
		rowTypeVec[StdFixRFFlow] = ARM_STRING;

		// Sum Coupons description
		CC_Ostringstream sumCouponsDesc;
		if (isFirstExer)
		{
			sumCouponsDesc << TARNColNamesTable[Coupon] << "[i]";
		}
		else
		{
			sumCouponsDesc << TARNColNamesTable[SumCoupons] << "[i-1]+" <<  TARNColNamesTable[Coupon] << "[i]";
		}

		rowDescVec[SumCoupons] = sumCouponsDesc.str();
		rowTypeVec[SumCoupons] = ARM_STRING;

		// Digital Funding

		CC_Ostringstream digitalUpDesc;
		CC_Ostringstream digitalDownDesc;

		digitalUpDesc << "MAX(" << TARNColNamesTable[Coupon] << "[i]-" << TARNColNamesTable[TargetCap] << "[i]" << ",0)";
		digitalDownDesc << "MAX(" << TARNColNamesTable[Coupon] << "[i]-" << TARNColNamesTable[TargetCap] << "[i]+" << DigitalSmoothingSpread << ",0)";

		rowDescVec[DigitalUp] = digitalUpDesc.str();
		rowTypeVec[DigitalUp] = ARM_STRING;
		rowDescVec[DigitalDown] = digitalDownDesc.str();
		rowTypeVec[DigitalDown] = ARM_STRING;

		// Is Alive description
		CC_Ostringstream isAliveDesc;
		if(itsDigitalSmoothingFlag)
		{
			if (!isFirstExer)
			{
				isAliveDesc << "if(" << TARNColNamesTable[IsAlived] << "[i-1]==1,";
				isAliveDesc << "1-(" << TARNColNamesTable[DigitalDown] << "[i]";
				isAliveDesc << "-" << TARNColNamesTable[DigitalUp]  << "[i])/" << DigitalSmoothingSpread;
				isAliveDesc << "," << TARNColNamesTable[IsAlived] << "[i-1])";
			}
			else
			{
				isAliveDesc << "1-(" << TARNColNamesTable[DigitalDown] << "[i]";
				isAliveDesc << "-" << TARNColNamesTable[DigitalUp]  << "[i])/" << DigitalSmoothingSpread;
			}
		}
		else
		{
			isAliveDesc << "if(" << TARNColNamesTable[SumCoupons];
 			isAliveDesc << "[i] < "  << TARNColNamesTable[LifeTimeCapTarget]  << "[i],1,0)";
		}

		rowDescVec[IsAlived] = isAliveDesc.str();
		rowTypeVec[IsAlived] = ARM_STRING;

		// Has Triggered description
		CC_Ostringstream hasTriggeredDesc;
		if (!isFirstExer)
		{
			hasTriggeredDesc << TARNColNamesTable[IsAlived] << "[i-1]";
			hasTriggeredDesc << "-" << TARNColNamesTable[IsAlived] << "[i]";
		}
		else
		{
			hasTriggeredDesc << "1-" << TARNColNamesTable[IsAlived] << "[i]";
		}

		rowDescVec[HasTriggered] = hasTriggeredDesc.str();
		rowTypeVec[HasTriggered] = ARM_STRING;

		// Reverse floater swap
		CC_Ostringstream rfSwapDesc;

		if(itsPayRec==K_PAY)
		{
			if (!isFirstEvent || !isArrearsCpn)
				rfSwapDesc <<  TARNColNamesTable[Funding] << fundingRef;

			rfSwapDesc << "-" <<  TARNColNamesTable[PaidCoupon] << "[i]" ;
		}
		else
		{
			rfSwapDesc <<  TARNColNamesTable[PaidCoupon] << "[i]";
			if (!isFirstEvent || !isArrearsCpn)
				rfSwapDesc << "-" << TARNColNamesTable[Funding] << fundingRef;
		}

		rowDescVec[Swap] = rfSwapDesc.str();
		rowTypeVec[Swap] = ARM_STRING;

		// Target Cap
		CC_Ostringstream targetCapDesc;

		if (isFirstExer)
		{
			targetCapDesc << TARNColNamesTable[LifeTimeCapTarget] << "[i]";
		}
		else
		{
			targetCapDesc << "MAX(" << TARNColNamesTable[TargetCap] << "[i-1]-" << TARNColNamesTable[Coupon]  << "[i-1]" << ",0)";
		}

		rowDescVec[TargetCap] = targetCapDesc.str();
		rowTypeVec[TargetCap] = ARM_STRING;

		// Target Floor
		CC_Ostringstream targetFloorDesc;

		if (isFirstExer)
		{
			targetFloorDesc << TARNColNamesTable[LifeTimeFloorTarget] << "[i]";
		}
		else
		{
			targetFloorDesc << "MAX(" << TARNColNamesTable[TargetFloor] << "[i-1]-" << TARNColNamesTable[Coupon]  << "[i-1]" << ",0)";
		}

		rowDescVec[TargetFloor] = targetFloorDesc.str();
		rowTypeVec[TargetFloor] = ARM_STRING;

		// Control Variable

		CC_Ostringstream controlVariableDesc;

		controlVariableDesc << "DIGITAL(" << basisModelName << "," ;
		controlVariableDesc	<< TARNColNamesTable[IndexStartDate] << "[i],";
		controlVariableDesc	<< TARNColNamesTable[IndexEndDate] << "[i],";
		controlVariableDesc <<  TARNColNamesTable[EstimExerciseStrike] << "[i],F)";
		
		rowDescVec[ControlVariable] = controlVariableDesc.str();
		rowTypeVec[ControlVariable] = ARM_STRING;

		/// Control variable2 (the std swap)
		CC_Ostringstream controlVariable2Desc;
		controlVariable2Desc << TARNColNamesTable[FundingPaid] << "[i]-" << TARNColNamesTable[CouponCF] << "[i]";
		rowDescVec[ControlVariable2] = controlVariable2Desc.str();
		rowTypeVec[ControlVariable2] = ARM_STRING;

		// Exercise Strike
		string sExStrikeDesc = CreateExStrikeDescription(isFirstExer, cpnModelName);
		rowDescVec[ExerciseStrike] = sExStrikeDesc;
		rowTypeVec[ExerciseStrike] = ARM_STRING;

		// Estimate Exercise Strike
		CC_Ostringstream estimateExerciseStrikeDesc;
		estimateExerciseStrikeDesc << 0;
		rowDescVec[EstimExerciseStrike] = estimateExerciseStrikeDesc.str();
		rowTypeVec[EstimExerciseStrike] = ARM_STRING;

		// Exercise Probas
		CC_Ostringstream exerciseProbaDesc;
		exerciseProbaDesc << "UNPAY(if(" << TARNColNamesTable[HasTriggered] << "[i],"; 
		exerciseProbaDesc << "1,0))";
		rowDescVec[ExerciseProba] = exerciseProbaDesc.str();
		rowTypeVec[ExerciseProba] = ARM_STRING;

		// Exercise Times
		CC_Ostringstream exerciseTimeDesc;
		double time = (payDate - itsStartDate.GetJulian()) / K_YEAR_LEN;

		if (isLastExer)
		{
			exerciseTimeDesc << "UNPAY(if(" << TARNColNamesTable[IsAlived] << "[i],"; 
			exerciseTimeDesc << time << ",0))";
		}
		else
		{
			exerciseTimeDesc << "UNPAY(if(" << TARNColNamesTable[HasTriggered] << "[i],"; 
			exerciseTimeDesc << time << ",0))";
		}

		rowDescVec[ExerciseTime] = exerciseTimeDesc.str();
		rowTypeVec[ExerciseTime] = ARM_STRING;

		// Life Time Cap
		CC_Ostringstream lifeTimeCapDesc;

		lifeTimeCapDesc << "MAX(" << TARNColNamesTable[Coupon] << "[i]-" << TARNColNamesTable[TargetCap] << "[i],0)*";
		lifeTimeCapDesc << TARNColNamesTable[DFPay] << "[i]*";
		lifeTimeCapDesc << TARNColNamesTable[Nominal] << "[i]";

		rowDescVec[LifeTimeCap] = lifeTimeCapDesc.str();
		rowTypeVec[LifeTimeCap] = ARM_STRING;

		// Last Life Time Cap
		CC_Ostringstream lastLifeTimeCapDesc;

		lastLifeTimeCapDesc << "IF(" << TARNColNamesTable[TargetCap] << "[i]>0,";
		lastLifeTimeCapDesc << TARNColNamesTable[LifeTimeCap] << "[i],0)";

		rowDescVec[LastLifeTimeCap] = lastLifeTimeCapDesc.str();
		rowTypeVec[LastLifeTimeCap] = ARM_STRING;

		// Life Time Floor

		if (isLastExer)
		{
			CC_Ostringstream lifeTimeFloorDesc;

			lifeTimeFloorDesc << "MAX(" << TARNColNamesTable[TargetFloor] << "[i]-" << TARNColNamesTable[Coupon] << "[i],0)*";
			lifeTimeFloorDesc << TARNColNamesTable[DFPay] << "[i]*";
			lifeTimeFloorDesc << TARNColNamesTable[Nominal] << "[i]";

			rowDescVec[LifeTimeFloor] = lifeTimeFloorDesc.str();
			rowTypeVec[LifeTimeFloor] = ARM_STRING;
		}

		CC_Ostringstream digitalFundingDesc;

		if (!isFirstEvent || !isArrearsCpn)
		{
			digitalFundingDesc << "if(!" << TARNColNamesTable[HasTriggered] << "[i],";
			digitalFundingDesc << "(1-" <<TARNColNamesTable[IsAlived] << "[i])*";
			digitalFundingDesc << TARNColNamesTable[Funding] << fundingRef;
			digitalFundingDesc << ",0)";
		}
		else
		{
			digitalFundingDesc << 0;
		}
		
		rowDescVec[DigitalFunding] = digitalFundingDesc.str();
		rowTypeVec[DigitalFunding] = ARM_STRING;

		// Digital Flow description
		CC_Ostringstream digitalDesc;

		if (!isLastExer)
		{
			digitalDesc << "if(" << TARNColNamesTable[HasTriggered] << "[i],";
			digitalDesc << TARNColNamesTable[Fees] << "[i]*" << TARNColNamesTable[DFPay] << "[i]";
			digitalDesc << ",0)";
		}
		else
		{
			digitalDesc << "if(" << TARNColNamesTable[IsAlived] << "[i]||"<< TARNColNamesTable[HasTriggered] << "[i],";
			digitalDesc << TARNColNamesTable[Fees] << "[i]*" << TARNColNamesTable[DFPay] << "[i]";
			digitalDesc << ",0)";
		}
		
		rowDescVec[Digital] = digitalDesc.str();
		rowTypeVec[Digital] = ARM_STRING;

		/// TARN Flow description
		CC_Ostringstream tarnDesc;
		
		if(itsPayRec==K_RCV)
		{
			tarnDesc << TARNColNamesTable[LifeTimeFloor] << "[i]-" << TARNColNamesTable[LifeTimeCap] << "[i]+" << TARNColNamesTable[DigitalFunding] << "[i]";
			if (!itsGlobalCapFlag)
				tarnDesc << "+" << TARNColNamesTable[LastLifeTimeCap] << "[i]";
		}
		else
		{
			tarnDesc << TARNColNamesTable[LifeTimeCap] << "[i]-" << TARNColNamesTable[LifeTimeFloor] << "[i]-" << TARNColNamesTable[DigitalFunding] << "[i]";
			if (!itsGlobalCapFlag)
				tarnDesc << "-" << TARNColNamesTable[LastLifeTimeCap]<< "[i]";
		}

		rowDescVec[TARN] = tarnDesc.str();
		rowTypeVec[TARN] = ARM_STRING;
	}

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}


string ARM_TARNCalculator::CreateCpnDescription(bool isFirstExer) const
{
	string couponDesc;
	couponDesc = "MIN(MAX(";
	couponDesc += TARNColNamesTable[Strike] + "[i]-";
	couponDesc += TARNColNamesTable[Leverage] + "[i]*" + TARNColNamesTable[CpnIndex] + "[i],"+ TARNColNamesTable[CouponMin] + "[i])," + TARNColNamesTable[CouponMax] + "[i])";

	return couponDesc;
}


string ARM_TARNCalculator::CreateCFCpnDdescription( bool isFirstExer, const string& cpnModelName) const
{
	CC_Ostringstream CouponCFDesc;
	string basisModelName = GetKeys()[BasisKey];
	CouponCFDesc << TARNColNamesTable[Leverage] << "[i]*(";
	CouponCFDesc << "Caplet(" << basisModelName << "," << TARNColNamesTable[IndexStartDate] << "[i]," + itsCpnIndexTerm  << ",(" << TARNColNamesTable[Strike] << "[i]-" << TARNColNamesTable[CouponMin] << "[i])/" << TARNColNamesTable[Leverage] << "[i],C)-";
	CouponCFDesc << "Caplet(" << basisModelName << "," << TARNColNamesTable[IndexStartDate] << "[i]," + itsCpnIndexTerm  << ",(" << TARNColNamesTable[Strike] << "[i]-" << TARNColNamesTable[CouponMax] << "[i])/" << TARNColNamesTable[Leverage] << "[i],C)";
	CouponCFDesc << ")*" << TARNColNamesTable[Nominal] << "[i]";
	return CouponCFDesc.str();
}



string ARM_TARNCalculator::CreateExStrikeDescription(bool isFirstExer, const string& cpnModelName) const
{
	string exerciseStrikeDesc;
	string basisModelName = GetKeys()[BasisKey];

	char minCF[15];
	sprintf(minCF,"%lf",MIN_CF_STRIKE);

	exerciseStrikeDesc = "MAX((" + TARNColNamesTable[Strike] + "[i]-" + TARNColNamesTable[TargetCap] + "[i]/" + TARNColNamesTable[IT] + "[i])/" + TARNColNamesTable[Leverage] + "[i]," + minCF + ")*";
	exerciseStrikeDesc += "DF(" + basisModelName + "," + TARNColNamesTable[IndexEndDate] + "[i])";

	return exerciseStrikeDesc;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[Swap] = zeroValue;
    rowTypeVec[Swap] = ARM_DOUBLE;

    rowDescVec[FundingPaid] = zeroValue;
    rowTypeVec[FundingPaid] = ARM_DOUBLE;

	rowDescVec[ExerciseStrike] = zeroValue;
    rowTypeVec[ExerciseStrike] = ARM_DOUBLE;

	rowDescVec[ExerciseProba] = zeroValue;
    rowTypeVec[ExerciseProba] = ARM_DOUBLE;

	rowDescVec[ExerciseTime] = zeroValue;
    rowTypeVec[ExerciseTime] = ARM_DOUBLE;

	rowDescVec[ControlVariable] = zeroValue;
    rowTypeVec[ControlVariable] = ARM_DOUBLE;

	rowDescVec[ControlVariable2] = zeroValue;
    rowTypeVec[ControlVariable2] = ARM_DOUBLE;

	rowDescVec[LifeTimeCap] = zeroValue;
    rowTypeVec[LifeTimeCap] = ARM_DOUBLE;

	rowDescVec[LastLifeTimeCap] = zeroValue;
    rowTypeVec[LastLifeTimeCap] = ARM_DOUBLE;

	rowDescVec[LifeTimeFloor] = zeroValue;
    rowTypeVec[LifeTimeFloor] = ARM_DOUBLE;

	rowDescVec[DigitalFunding] = zeroValue;
    rowTypeVec[DigitalFunding] = ARM_DOUBLE;

	rowDescVec[Digital] = zeroValue;
    rowTypeVec[Digital] = ARM_DOUBLE;

	rowDescVec[TARN] = zeroValue;
    rowTypeVec[TARN] = ARM_DOUBLE;

	rowDescVec[PaidNb] = "1";
    rowTypeVec[TARN] = ARM_DOUBLE;
	
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the TARN. The
///          DateStripCombiner merges event dates of each
///          legs
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_TARNCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the TARN. The
///          DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_TARNCalculator::DatesStructure() const
{
    /// Get reset & payment calendars
    const char* resetCalendar   = itsCpnResetCal.c_str();
    const char* payCalendar     = itsCpnPayCal.c_str();

    int fwdRule		=   K_MOD_FOLLOWING;	// for forward dates
    int resetTiming =   K_ADVANCE;
	int payGap      =   GETDEFAULTVALUE;
    int payTiming   =   K_ARREARS;
	
    ARM_DateStrip RFLegSched(itsStartDate,itsEndDate,itsCpnFreq,itsCpnDayCount,resetCalendar,fwdRule,itsIntRule,itsStubRule,
							itsCpnResetGap,	itsCpnFreq,payGap,payCalendar, itsCpnTiming,payTiming);


	char* ccyName = GetCurrencyUnit()->GetCcyName();

	// For the in arrears case we create a fake start date to get the first 
	// reset date of the funding leg
	if (itsCpnTiming == K_ARREARS)
	{
		ARM_Date  firstFundingResetDate(itsStartDate);
// FIXMEFRED: mig.vc8 (25/05/2007 14:32:36):fabs(int) doesnt exist
		firstFundingResetDate.PreviousBusinessDay(abs(itsCpnResetGap),ccyName);
		RFLegSched.InsertDate(0,0,0,0,0,firstFundingResetDate.GetJulian(),0,0,0);
	}	

	//Summit customised reset dates
	//We replace the reset, fwdstart and fwdend dates.
	if (GetCustomResetFlag() == true)
	{
		int spotDays = GetCurrencyUnit()->GetSpotDays();
	
		if (itsCpnTiming == K_ARREARS)
		{
			std::vector<double> custSmallResetDates = GetCustomResetDates();
			std::vector<double> currResetDates = *(RFLegSched.GetResetDates());
			
			std::vector<double> customResetDates = std::vector<double>(custSmallResetDates.size()+1);
			//The first reset date is not cutstomised.
			customResetDates.Elt(0) = currResetDates.Elt(0);
			//The other reset dates are customised.
			for (int i=1; i< currResetDates.size();i++)
			{
				customResetDates.Elt(i) = custSmallResetDates.Elt(i-1);
			}
			//Then, we build the forward schedule.
			std::vector<double> customFwdStartDates = customResetDates;
			std::vector<double> customFwdEndDates = std::vector<double>(customResetDates.size());

			for(i=0; i<customResetDates.size(); ++i )
			{
				ARM_Date tmpStartDate(customFwdStartDates.Elt(i));
				tmpStartDate.GapBusinessDay(spotDays,(char *) resetCalendar);
				customFwdStartDates.Elt(i) = tmpStartDate.GetJulian();
			
				ARM_Date tmpEndDate(customFwdStartDates.Elt(i));
				tmpEndDate.AddPeriod(itsCpnFreq, (char*)resetCalendar);
				tmpEndDate.GoodBusinessDay(fwdRule*itsIntRule, (char*)resetCalendar);
				customFwdEndDates.Elt(i) = tmpEndDate.GetJulian();
			}

			RFLegSched.SetResetDates(&customResetDates);
			RFLegSched.SetFwdRateStartDates(&customFwdStartDates);
			RFLegSched.SetFwdRateEndDates(&customFwdEndDates);
		}
		else
		{
			//For the standard case, all the reset dates are customised.
			std::vector<double> customResetDates = GetCustomResetDates();
		
			//Sanity check
			if (((RFLegSched.GetResetDates())->size()) != (customResetDates.size()))
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Your custom reset dates are invalid, please try again.");

			std::vector<double> customFwdStartDates = customResetDates;
			std::vector<double> customFwdEndDates = std::vector<double>(customResetDates.size());

			for(int i=0; i<customResetDates.size(); ++i )
			{
				ARM_Date tmpStartDate(customFwdStartDates.Elt(i));
				tmpStartDate.GapBusinessDay(spotDays,(char *) resetCalendar);
				customFwdStartDates.Elt(i) = tmpStartDate.GetJulian();
			
				ARM_Date tmpEndDate(customFwdStartDates.Elt(i));
				tmpEndDate.AddPeriod(itsCpnFreq, (char*)resetCalendar);
				tmpEndDate.GoodBusinessDay(fwdRule*itsIntRule, (char*)resetCalendar);
				customFwdEndDates.Elt(i) = tmpEndDate.GetJulian();
			}
			RFLegSched.SetResetDates(&customResetDates);
			RFLegSched.SetFwdRateStartDates(&customFwdStartDates);
			RFLegSched.SetFwdRateEndDates(&customFwdEndDates);
		}
	}

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(NB_TARN_SCHED,NULL);
    SchedVect[RF_CPN_SCHED]     = &RFLegSched;

    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();
	std::vector<double>& rfResetDates = EventSchedule.GetDateStrip(RF_CPN_SCHED)->GetResetDates();

	CC_MUTABLE( ARM_TARNCalculator, itsFirstEventIdx ) = 0;
	while(itsFirstEventIdx < eventDates->size() && (*eventDates)[itsFirstEventIdx] < asOfDate)
        ++CC_MUTABLE( ARM_TARNCalculator, itsFirstEventIdx );

	// Check if already start
	CC_MUTABLE( ARM_TARNCalculator, itsFirstRFIdx ) = 0;
	if (itsFirstEventIdx != 0)
		CC_MUTABLE( ARM_TARNCalculator, itsFirstRFIdx ) = itsFirstEventIdx;
	else
	{
		if (itsCpnTiming==K_ARREARS)
			CC_MUTABLE( ARM_TARNCalculator, itsFirstRFIdx ) = itsFirstEventIdx+1;
	}

	if(IsBasis()){
		
		ARM_DateStrip structSched(itsStartDate,itsEndDate,itsCpnFreq,itsCpnDayCount,resetCalendar,fwdRule,itsIntRule,itsStubRule,
			itsCpnResetGap,	itsCpnFreq,payGap,payCalendar, resetTiming,payTiming);
		
		ARM_DateStrip FundingSched(itsStartDate,itsEndDate,itsFundFreq,itsFundDayCount,resetCalendar,fwdRule,itsIntRule,itsStubRule,
			GETDEFAULTVALUE,itsFundFreq, payGap, payCalendar, resetTiming, payTiming);

	const_cast< ARM_TARNCalculator* >(this)->itsStructDateStrip = ARM_DateStripPtr(new ARM_DateStrip(structSched));	
	const_cast< ARM_TARNCalculator* >(this)->itsFundDateStrip   = ARM_DateStripPtr(new ARM_DateStrip(FundingSched));
		std::vector<double>& fundResetDates = FundingSched.GetResetDates() ;
		double lag, margin;
		for (size_t i(0); i<fundResetDates->size(); i++)
		{
			lag = (*fundResetDates)[i]-asOfDate;
			margin = itsFundSpread.Interpolate(lag);
			const_cast< ARM_TARNCalculator* >(this)->itsvInitialFundSpread.push_back(margin );
			const_cast< ARM_TARNCalculator* >(this)->itsvFundNominal.push_back(itsFundNominal.Interpolate(lag));
		}

		fundResetDates = itsStructDateStrip->GetResetDates();
		for (i= 0; i<fundResetDates->size(); i++){
			lag = (*fundResetDates)[i]-asOfDate;
			const_cast< ARM_TARNCalculator* >(this)->itvCpnsNominal.push_back(itsNominal.Interpolate(lag));
		}
		
		const_cast< ARM_TARNCalculator* >(this)->itsFundSpread.SetOrdinates(ComputeDomesticBasis());	
	}
    return EventSchedule;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetModel_SFRM2F
///	Returns: ARM_PricingModel*
///	Action : creates the model
/////////////////////////////////////////////////////////////////
ARM_PricingModel*  ARM_TARNCalculator::GetModel_SFRM2F(ARM_ZeroCurve* pCurve)
{
	ARM_PricingModel* pModel = 0;
	/// Creates default values for volatility & mean reversion
	/// (calibration will be called at pricing time) and set them
	/// in the reference model
	std::vector<double> defaultTimes;
	std::vector<double> defaultSigmas;
	std::vector<double> defaultBetas;

	/// Get asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam newBetaParam( ARM_ModelParamType::Beta, &defaultBetas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* betaParam = NULL;	
	betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
	if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Beta Param for key=" << GetKeys()[BetaKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
	if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Correl Param for key=" << GetKeys()[CorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	ARM_ModelParamVector paramVector(4);
	paramVector[0] = &volParam;
	paramVector[1] = mrsParam;
	if (itsDigitalCalibMode != NOCalib)
	{
		paramVector[2] = &newBetaParam;
	}
	else
	{
		paramVector[2] = betaParam;
	}

	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,SFRM_VOL_TYPE);

	/// Build the default stochastic model of the calculator : SFRM 2F
	pModel = new ARM_SFRM( CreateClonedPtr( pCurve ), *pSFRMModelParams);
	// Delte the model params because it is cloned in the model
	delete pSFRMModelParams;
	return pModel;

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetModel_SBGM
///	Returns: ARM_PricingModel*
///	Action : creates the model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_TARNCalculator::GetModel_SBGM(ARM_ZeroCurve* pCurve)
{
	int nbSteps   = (itsNbIterations.size()>4)?itsNbIterations[4]:PDE_NB_STEPS;
	int gridSize  = (itsNbIterations.size()>5)?itsNbIterations[5]:PDE_GRIDSIZE;
	double nbStdDevs = (itsNbIterations.size()>6)?itsNbIterations[6]:PDE_NB_STDEV;

	ARM_PricingModel* pModel = new ARM_SmiledFRM(CreateClonedPtr( pCurve ),NULL,nbSteps,gridSize,nbStdDevs,true, false, ARM_ModelParamsSmiled::LocalVolatilityWithRescaling, itsSmiledFRMRescallingFlag);

	ARM_ModelParam* humpParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[HumpKey]));
	if(!humpParam || humpParam->GetType() != ARM_ModelParamType::Hump)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Hump Param for key=" << GetKeys()[HumpKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* betaCorrel = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaCorrelKey]));
	if(!betaCorrel || betaCorrel->GetType() != ARM_ModelParamType::BetaCorrelation)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : BetaCorrelation Param for key=" << GetKeys()[BetaCorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}

	ARM_ModelParam* reCorrel  = NULL;
	
	if ((ReCorrelKey < GetKeys().size()) && !GetMktDataManager()->TestIfKeyMissing(GetKeys()[ReCorrelKey]))
		reCorrel = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[ReCorrelKey]));
		
	ARM_ModelParamVector modelParams(2);
	modelParams[0] = humpParam;
	modelParams[1] = betaCorrel;

	if (reCorrel)
		modelParams.push_back(reCorrel);
	    
	ARM_ModelParamsSmiled* modelParamsSmiledFRM;
	if (itsNbIterations.size()>7)
		modelParamsSmiledFRM = new ARM_ModelParamsSmiled(modelParams,itsNbIterations[7]);
	else
		modelParamsSmiledFRM = new ARM_ModelParamsSmiled(modelParams,SBGM_NB_FACTORS);

	pModel->SetModelParams(*modelParamsSmiledFRM);
	delete modelParamsSmiledFRM;
	return pModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetModel_HK
///	Returns: ARM_PricingModel*
///	Action : creates the model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_TARNCalculator::GetModel_HK(ARM_ZeroCurve* pCurve)
{
	std::vector<double> defaultTimes(1,0.);
	std::vector<double> defaultSigmas(1,0.01);
		
	ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}	

	ARM_ModelParamVector modelParams(2);
	modelParams[0] = volParam;
	modelParams[1] = mrsParam;
	   
	/// Create Markov Functional
	ARM_PricingModel*pModel = new ARM_MarkovFunctional( CreateClonedPtr( pCurve ) );
	/// Sets model params
    pModel->SetModelParams(ARM_ModelParamsMF(modelParams));
	return pModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetModel()
{
	/// Get yield curve
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_PricingModel* pModel = 0;
	ARM_MCMethod* pMCMethod = 0;

	/// Base Generator 1
	ARM_RandomGeneratorPtr  pBaseRandomGen1( ARM_RandGenFactory.Instance()->CreateRandGen( 
		(ARM_RandGenFactoryImp::BaseGenType)	 ARM_ArgConv_BaseGenAlgoType.GetNumber(itsGenType2),
		ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	/// Base Generator 2
	ARM_RandomGeneratorPtr  pBaseRandomGen2(ARM_RandGenFactory.Instance()->CreateRandGen( 
		(ARM_RandGenFactoryImp::BaseGenType)	 ARM_ArgConv_BaseGenAlgoType.GetNumber(itsGenType1),
		ARM_RandGenFactoryImp::UnknownTransformAlgo ) );

	ARM_RandomGeneratorPtr normRandGen1(NULL); 
	ARM_RandomGeneratorPtr normRandGen2(NULL); 

	// Normal Rand Gen
	// Inv Norm Cum for the quasi random
	normRandGen1 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::InvNormCum,
		pBaseRandomGen1 ) );

	// Box Muller for the pseudo random
	normRandGen2 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::BoxMuller,
		pBaseRandomGen2 ) );

	/// tranposed random gen
	ARM_RandomGeneratorPtr transpRandGen2( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::Transposer,
		ARM_RandomGeneratorPtr(normRandGen2),
		ARM_RandomGeneratorPtr(NULL),
		-1,
		1,
		1,
		ARM_GP_T_Vector<size_t>(1,10),
		4.0,
		0,
		0,
		(ARM_RandGenFactoryImp::RandGenOrder) ARM_ArgConv_RandGenOrder.GetNumber(itsPathOrder)));

	ARM_RandomGeneratorPtr antiRandGen1, antiRandGen2;
	if (itsAntithetic)
	{
		antiRandGen1 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen( 
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::AntitheticOne,
				normRandGen1 ) );

		antiRandGen2 = ARM_RandomGeneratorPtr( ARM_RandGenFactory.Instance()->CreateRandGen(
				ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
				ARM_RandGenFactoryImp::AntitheticOne,
				transpRandGen2 ) );
	}
	else
	{
		antiRandGen1 = normRandGen1;
		antiRandGen2 = transpRandGen2;
	}

	ARM_RandomGeneratorPtr mixteRandGen( ARM_RandGenFactory.Instance()->CreateRandGen( 
		ARM_RandGenFactoryImp::UnknownBaseGenAlgorithm,
		ARM_RandGenFactoryImp::MixteGen,
		ARM_RandomGeneratorPtr(antiRandGen1),
		ARM_RandomGeneratorPtr(antiRandGen2),
		-1,
		1,
		1,
		ARM_GP_T_Vector<size_t>(1,10),
		4.0,
		itsFirstNbTimes,
		itsFirstNbDims));
		

	int nbInterSteps   = (itsNbIterations.size()>3)?itsNbIterations[3]:MC_NB_INTER_STEPS; 
	ARM_TimeStepPerYearScheduler scheduler(nbInterSteps);

	// Path Scheme
	int pathSchemeType = ARM_ArgConv_PathSchemeType.GetNumber(itsPathScheme);
	ARM_PathSchemePtr pathScheme(ARM_PathSchemeFactory.Instance()->CreatePathScheme(pathSchemeType));

	if (itsModelType==ARM_PricingModelType::SFRM2F)
	{
		pModel = GetModel_SFRM2F(pCurve);
		ARM_NormalCentredSamplerND sampler(&scheduler);
		pMCMethod = new ARM_MCMethod(itsNbIterations[0],mixteRandGen,&sampler,itsNbIterations[2],ARM_ImpSamplerPtr(NULL),pathScheme);


	}
	else if (itsModelType==ARM_PricingModelType::SBGM)
	{
		int nbInterStepsSBGM   = (itsNbIterations.size()>3)?itsNbIterations[3]:SBGM_MC_NB_INTER_STEPS; 
		ARM_TimeStepPerYearScheduler schedulerSBGM(nbInterStepsSBGM);
		pModel = GetModel_SBGM(pCurve);
		ARM_NormalCentredSamplerND sampler(&schedulerSBGM);
		pMCMethod = new ARM_MCMethod(itsNbIterations[0],mixteRandGen,&sampler,itsNbIterations[2],ARM_ImpSamplerPtr(NULL),pathScheme);

	}
	else if (itsModelType==ARM_PricingModelType::HK)
	{
		pModel = GetModel_HK(pCurve);
		ARM_MeanRevertingSamplerND sampler(&scheduler);
		pMCMethod = new ARM_MCMethod(itsNbIterations[0],mixteRandGen,&sampler,itsNbIterations[2],ARM_ImpSamplerPtr(NULL),pathScheme);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" ARM_TARNCalculator::CreateAndSetModel : only SFRM2F, SBGM and HK are supported !" );
	}
	
	ARM_PricingModelPtr refModel( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(pModel)));	
	ARM_PricingModelPtr hybridmodel;	 
	if(IsBasis())
	{
		ARM_StringVector names (1);
		vector<ARM_PricingModelPtr> models (1);
		ARM_StringVectorVector depends(1);
		
		// Key names
		names[0] = GetKeys()[YcKey];

		/// models
		models[0] = refModel;
		
		/// basis => build a Hybrid BFIR
		/// Push back model pricing names
        names.push_back(GetKeys()[BasisKey]);

		// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
		models.push_back(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve)))));
		depends.push_back(ARM_StringVector (1,names[0]));

		ARM_ModelNameMap modelMap (names, models, depends);

		///	modelMap & correls are cloned in multi assets model
		 hybridmodel =ARM_PricingModelPtr (new ARM_MultiAssetsModel ( &modelMap) );

	}
	else{
		hybridmodel = refModel;
	}
	
	ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    hybridmodel->SetNumeraire(numeraire);
	hybridmodel->SetNumMethod(ARM_NumMethodPtr( pMCMethod ) );

	/// Set the model
	SetPricingModel(hybridmodel);	
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CleanAfterPricing
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CleanAfterPricing(void)
{
    ARM_PricingModelPtr pModel = GetPricingModel();

    if (!(pModel.IsNull()))
       pModel->SetNumMethod(ARM_NumMethodPtr(NULL));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateEmptyCalibration(
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool swaptionCalibFlag)
{
	/// Create cap floor portfolio (sigma calibration)
    ARM_StdPortfolioPtrVector capFloorAndDigitalPFs(CreateTARNCapletFloorletAndDigital());

    /// Build a volatility bootstrap calibration on the latter portfolio
    /// The CalibMethod object will be cloned because by default it is not shared

	ARM_CalibMethod* volCalibMethod = NULL;
	ARM_CalibMethod* lastCalibMethod = NULL;

	int nbFlows = capFloorAndDigitalPFs.size()/2;

	for(int i = 0; i < nbFlows; ++i)
	{
		ARM_CalibMethod* newCFCalibMethod = new ARM_CalibMethod(
			capFloorAndDigitalPFs[2*i],
			ARM_ModelParamVector(),
			ARM_CalibMethodType::Bootstrap1D,
			ARM_MAX_ITER,
			ARM_CalibrationTarget::PriceTarget,
			NULL,
			NULL);

		if ((digitalCalibMode!=NOCalib) && ((i != nbFlows-1) || capFloorCalibMode == ExerStrike))
		{
			ARM_CalibMethod* newDigitCalibMethod = NULL;
			newDigitCalibMethod = new ARM_CalibMethod(
				capFloorAndDigitalPFs[2*i+1],
				ARM_ModelParamVector(),
				ARM_CalibMethodType::Optimize,
				ARM_MAX_ITER,
				ARM_CalibrationTarget::PriceTarget);
			newDigitCalibMethod->SetlinkedMethod(newCFCalibMethod);
			newDigitCalibMethod->SetPreviousMethod(lastCalibMethod);
			lastCalibMethod = newDigitCalibMethod;
			
		}
		else
		{
			newCFCalibMethod->SetPreviousMethod(lastCalibMethod);
			lastCalibMethod = newCFCalibMethod;
		}	
	}

	volCalibMethod = lastCalibMethod;

    if(swaptionCalibFlag)
    {
        /// Create standard diagonal swaption portfolio (MRS calibration)
        ARM_StdPortfolioPtr diagonalSwaptionPF(CreateDiagonalSwaption());

        /// Build a MRS optimisation and embed the volatility bootstrapping
        /// The CalibMethod object will be cloned because by default it is not shared
        ARM_CalibMethod* mrsCalib = new ARM_CalibMethod(diagonalSwaptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
                                 ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);
		mrsCalib->SetlinkedMethod(volCalibMethod);
        SetCalibMethod( ARM_CalibMethodPtr( mrsCalib ) );
    }
    else
        SetCalibMethod( ARM_CalibMethodPtr( volCalibMethod ) );
}

///////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: FindValueInVector
///	Returns: size_t
///	Action : 
////////////////////////////////////////////////////
double ARM_TARNCalculator::FindValueInVector( const std::vector<double>& timeVector, double time, double daysNb ) const
{
	bool found = false;
	int idx,N = timeVector.size();
	for (int i=0;(i<N && (!found));i++)
	{
		if ( abs( time - timeVector[i] ) <= daysNb )
		{
			found = true;
			idx=i;
		}
	}
	if (found)
		return timeVector[idx];
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + "ARM_TARNCalculator::FindValueInVector: date not found in datestrip!");
}

///////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: FindIdxInVector
///	Returns: size_t
///	Action : 
////////////////////////////////////////////////////
bool ARM_TARNCalculator::FindIdxInVector( const std::vector<double>& timeVector, double time, double daysNb, size_t& idx) const
{
	bool found = false;
	int N = timeVector.size();
	for (int i=0;(i<N && (!found));i++)
		if ( abs( time - timeVector[i] ) <= daysNb )
		{
			found = true;
			idx=i;
		}
	return found;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetCalibration()
{
	switch (itsModelType)
	{
	case ARM_PricingModelType::SFRM2F:
		CreateAndSetCalibration_SFRM2F();
		break;
	case ARM_PricingModelType::SBGM:
		CreateAndSetCalibration_SBGM();
		break;
	case ARM_PricingModelType::HK:
		CreateAndSetCalibration_HK();
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::CreateAndSetCalibration : wrong model type !" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetCalibration_SFRM2F
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetCalibration_SFRM2F()
{
	bool isFreezeWeights = false;					// reinit non null weights in the portfolio
	bool isInitSigma = true;						// init sigma param for calibration (bounds & init guess)
	bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

	if ((itsCapFloorCalibMode == ExerStrike) || (itsDigitalCalibMode == ExerStrike))
	{
		CreateEmptyCalibration(
			RFStrike,
			NOCalib,
			false);

		ComputeTARNCapletFloorletPrices(
				isFreezeWeights,
				isInitSigma,
				isArrearsCpn,
				RFStrike,
				NOCalib,
				false);

		ComputeMCExerciseStrikes();
	}

	CreateEmptyCalibration(
		itsCapFloorCalibMode,
		itsDigitalCalibMode,
		itsOSWCalibFlag);

	if(itsCapFloorCalibMode != NOCalib)
	{
		ComputeTARNCapletFloorletPrices(
			isFreezeWeights,
			isInitSigma,
			isArrearsCpn,
			itsCapFloorCalibMode,
			itsDigitalCalibMode,
			itsOSWCalibFlag);
	}

	/// Compute diagonal swaption target prices
	if(itsOSWCalibFlag)
	{
		ComputeDiagonalSwaptionPrice(isFreezeWeights,isArrearsCpn);
		PrepareMeanRevForSwaptionCalib();
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetCalibration_SBGM
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetCalibration_SBGM()
{
	CreateAndSetDensityCalibration(ARM_StdPortfolioPtr(NULL));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetCalibration_HK
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetCalibration_HK()
{
 	bool isFreezeWeights = false;					// reinit non null weights in the portfolio
	bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

	if(itsOSWCalibFlag)
	{
		ARM_StdPortfolioPtr pfPtr = CreateDiagonalSwaption();
		CreateAndSetDensityCalibration(pfPtr);
		ComputeDiagonalSwaptionPrice(isFreezeWeights,isArrearsCpn);
	}
	else
		CreateAndSetDensityCalibration(ARM_StdPortfolioPtr(NULL));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateAndSetNumCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::CreateAndSetDensityCalibration(ARM_StdPortfolioPtr pfPtr)
{
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_BSSmiledModel* cfBSSmiledModel = dynamic_cast< ARM_BSSmiledModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	if (!cfBSSmiledModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" ARM_TARNCalculator::CreateAndSetNumCalibration: no BSSMiledModel available !" );
		
	ARM_VolCurve* pSigma	= cfBSSmiledModel->GetVolatility();

	ARM_VolCurve* pAlpha	= cfBSSmiledModel->FromSigmaToAlpha(pSigma);
	ARM_VolCurve* pBeta		= cfBSSmiledModel->GetBeta();
	ARM_VolCurve* pRho		= cfBSSmiledModel->GetRho();
	ARM_VolCurve* pNu		= cfBSSmiledModel->GetNu();

	std::vector<double>& pStrike	= NULL;
	ARM_DateStrip* pDS		= GetCalibSchedule(pStrike);	////// A CHANGER!!!!!!!!!!!!!!!!!

	std::vector<double>& pRD		= pDS->GetResetDates();	//not cloned, do not delete
	
	ARM_VanillaSecDensityPtrVector vanillaSecDensities;
	std::vector<double>::iterator resetD	= pRD->begin();
	std::vector<double>::iterator startD	= pDS->GetFwdRateStartDates()->begin();
	std::vector<double>::iterator endD	= pDS->GetFwdRateEndDates()->begin();
	std::vector<double>::iterator strike	= pStrike->begin();

	for(; (resetD != pRD->end()) && (strike != pStrike->end()); ++resetD,++startD,++endD,++strike)
	{
		ARM_VanillaSecurityDensity* pVSD	= GetMarketDensity(
												*resetD, *startD, *endD,
												pSigma, pBeta, pRho, pNu, cfBSSmiledModel->GetSABRFlag(),pCurve,*strike);

		vanillaSecDensities.push_back(ARM_VanillaSecDensityPtr(pVSD));
	}
		
	string MethodTypeStr = "Numerical";
	ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);

	string MktTargetStr = "UNKNOWN_TAR";
	ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);

	ARM_CalibMethod* pMethod = new ARM_CalibMethod ( pfPtr,
										ARM_ModelParamVector(0),
										methodType,
										100,
										mktTargetType,
										NULL, 
										NULL, 
										false, 
										0,
										1,
										true, 
										ARM_DateStripPtr(pDS),
										vanillaSecDensities);

	SetCalibMethod(ARM_CalibMethodPtr(pMethod));	////// A CHANGER!!!!!!!!!!!!!!!!!
	delete pStrike;
}


///////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetCalibSchedule
///	Returns: ARM_DateStrip*
///	Action : 
////////////////////////////////////////////////////
ARM_DateStrip* ARM_TARNCalculator::GetCalibSchedule(std::vector<double>&& pStrike) const
{
	bool isArrearsCpn					= (itsCpnTiming == K_ARREARS);
	int cpnIndexFreq					= ARM_ArgConv_MatFrequency.GetNumber(itsCpnIndexTerm);
	const ARM_DealDescription dealDesc	= GetGenSecurity()->GetDealDescription();
	size_t nbEvents						= dealDesc.GetRowsNb();
		
	size_t firstIndex					= (isArrearsCpn?2:1);
	size_t nbIndex						= (nbEvents-firstIndex);

	int maxFreq							= cpnIndexFreq<itsCpnFreq?itsCpnFreq:cpnIndexFreq;
	
	ARM_Date Date1( atof( dealDesc.GetElem(firstIndex,StartDate).c_str() ));
	ARM_Date Date2( itsEndDate);
	size_t nbmonth ;
	if		(cpnIndexFreq == K_ANNUAL)		nbmonth	= 12;
	else if (cpnIndexFreq == K_SEMIANNUAL)	nbmonth	= 6;
	else if (cpnIndexFreq == K_QUARTERLY)	nbmonth	= 3;
	else if (cpnIndexFreq == K_MONTHLY)		nbmonth	= 1;
	else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::GetCalibSchedule : index type not supported" );

	ARM_DateStrip* pDS;
	if (isArrearsCpn)
		Date2.AddMonths(nbmonth);
	else
		if (cpnIndexFreq < itsCpnFreq)
		{
			size_t aux ;
			if		(itsCpnFreq == K_ANNUAL)		aux	= 12;
			else if (itsCpnFreq == K_SEMIANNUAL)	aux	= 6;
			else if (itsCpnFreq == K_QUARTERLY)		aux	= 3;
			else if (itsCpnFreq == K_MONTHLY)		aux	= 1;
			else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::GetCalibSchedule : cpn freq not supported" );
			Date2.AddMonths(-aux);
			Date2.AddMonths(nbmonth);
		}

	int gap;
	if (GetCustomResetFlag() == true)
	{
		ARM_Date reset( atof(dealDesc.GetElem(1,EventDate).c_str()) );
		ARM_Date start( atof(dealDesc.GetElem(firstIndex,StartDate).c_str()) );
		gap = reset.GetJulian()-start.GetJulian();
		gap = gap - 2 * gap / 7;
	}
	else
		gap = itsCpnResetGap;


	pDS = new ARM_DateStrip( Date1, Date2, maxFreq,
							itsCpnDayCount, itsCpnResetCal.c_str(),	K_MOD_FOLLOWING,
							itsIntRule, K_SHORTSTART, gap,
							maxFreq, GETDEFAULTVALUE, itsCpnPayCal.c_str(),
							K_ADVANCE, K_ARREARS);
		
	size_t size = pDS->GetResetDates()->size();
	std::vector<double> julianResetDates(*pDS->GetResetDates());
	std::vector<double> julianStartDates(*pDS->GetFlowStartDates());
	std::vector<double> julianEndDates(*pDS->GetFlowEndDates());
	std::vector<double> julianFwdStartDates(*pDS->GetFwdRateStartDates());
	std::vector<double> julianFwdEndDates(*pDS->GetFwdRateEndDates());

	vector<int> vFreq(size,-1);
	pStrike = new std::vector<double>(size,-1.);


	size_t nper = maxFreq / cpnIndexFreq;
	
	/// 1. - adjust schedule for TARN dates
	size_t LOOK_FOR_DATE = 10;

	size_t idx;
	for (size_t eventIdx = 1 ; eventIdx < nbEvents ; ++eventIdx)
	{
		ARM_Date indexResetDate( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );
		idx = (eventIdx - 1) * maxFreq / itsCpnFreq;
		julianResetDates.Elt(idx)  = indexResetDate.GetJulian();
		if (!(eventIdx==1 && isArrearsCpn))
		{
			ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );
			ARM_Date indexEndDate(   atof(dealDesc.GetElem(eventIdx,IndexEndDate).c_str()) );
			julianFwdStartDates.Elt(idx) = indexStartDate.GetJulian();
			julianFwdEndDates.Elt(idx)   = indexEndDate.GetJulian();
			vFreq[idx]						= cpnIndexFreq;
			(*pStrike)[idx]					= atof(dealDesc.GetElem(eventIdx,Strike).c_str())
											/ atof(dealDesc.GetElem(eventIdx,Leverage).c_str());
		}
	}

	/// 2. - add extra dates if needed
	idx = 0;
	for (std::vector<double>::iterator iter = julianResetDates.begin() ; iter != julianResetDates.end() ; ++iter,++idx)
	{
		ARM_Date indexStartDate(julianResetDates.Elt(idx));
		int defaultResetGap = GetCurrencyUnit()->GetSpotDays();
		indexStartDate.GapBusinessDay(defaultResetGap, const_cast<char *>(itsCpnResetCal.c_str()));
	
		ARM_Date indexEndDate(indexStartDate);
		indexEndDate.AddPeriod(ARM_ArgConvReverse_MatFrequency.GetString(maxFreq));

		julianStartDates.Elt(idx)	= indexStartDate.GetJulian();
		julianEndDates.Elt(idx)		= indexEndDate.GetJulian();
			
		if ( vFreq[idx]==-1)
		{
			if ( idx % (maxFreq / itsCpnFreq) ==1 )
				julianStartDates.Elt(idx) = julianEndDates.Elt(idx-1);
			julianFwdStartDates.Elt(idx)	= julianStartDates.Elt(idx);
			julianFwdEndDates.Elt(idx)		= julianEndDates.Elt(idx);
			vFreq[idx]						= maxFreq;
			(*pStrike)[idx]					= -1;
		}
	}

	/// 3. - recombine schedule
	for ( idx = 0 ; idx < size-1 ; ++idx)
		julianEndDates.Elt(idx) = julianStartDates.Elt(idx+1);
		
	/// 4. - adjust dates 
	for ( idx = 0 ; idx < size ; ++idx)
	{
		julianFwdStartDates.Elt(idx) = julianStartDates.Elt(idx);
		julianFwdEndDates.Elt(idx) = julianEndDates.Elt(CC_Min(idx + nper - 1, size - 1));
	}
        
	pDS->SetResetDates(&julianResetDates);
	pDS->SetFlowStartDates(&julianStartDates);
	pDS->SetFlowEndDates(&julianEndDates);
	pDS->SetFwdRateStartDates(&julianFwdStartDates);
	pDS->SetFwdRateEndDates(&julianFwdEndDates);
	
	return pDS;
}

///////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetMarketDensity
///	Returns: ARM_VanillaSecurityDensity*
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensity* ARM_TARNCalculator::GetMarketDensity( 
								double julianRD, double julianSD, double julianED,
								ARM_VolCurve* pSigma, ARM_VolCurve* pBeta, ARM_VolCurve* pRho, ARM_VolCurve* pNu,
								int sabrFlag, ARM_ZeroCurve* pCurve, double strike) const
{
	int sabrGridSize  = (itsNbIterations.size()>8)?itsNbIterations[8]:SABR_GRIDSIZE;

	const ARM_Date asOfDate				= GetMktDataManager()->GetAsOfDate();
	double julianAsof					= asOfDate.GetJulian();

	ARM_VanillaSecDensityPtrVector vanillaSecDensities;
	
	double mat		= (julianRD-julianAsof)/K_YEAR_LEN;
	ARM_Currency* pCcy	= pCurve->GetCurrencyUnit();
	int floatDayCount = pCcy->GetLiborIndexDayCount();

	double delta	= CountYears( floatDayCount, julianSD, julianED);
	double theta	= CC_Round( (julianED-julianSD) / 30 ) / 12.;
	int freq		= ARM_ArgConv_MatFrequency.GetNumber(itsCpnIndexTerm);
	double fwd      = 1. / delta * (pCurve->DiscountPrice((julianSD-julianAsof)/K_YEAR_LEN)/pCurve->DiscountPrice((julianED-julianAsof)/K_YEAR_LEN)-1.);
			
	ARM_DensityFunctor* pDF;

	if ( mat < K_DOUBLE_TOL )
	{
		pDF	= new ARM_ShiftedLNDensityFunctor(0.2,0.0);
	}
	else
	{
		double sigma	= pSigma->ComputeVolatility(mat,theta)/100.0;
		double beta		= pBeta->ComputeVolatility(mat,theta);
		double rho		= pRho->ComputeVolatility(mat,theta);
		double nu		= pNu->ComputeVolatility(mat,theta);
		double alpha	= ComputeAlpha(fwd, fwd, mat, sigma,
												rho, nu, beta, 0, sabrFlag);

		ARM_SABRDensityFunctor* pSDF		= new ARM_SABRDensityFunctor(
												alpha, 
												beta, 
												rho, 
												nu, 
												sabrFlag,sabrGridSize);

		if (itsCapFloorCalibMode == RFStrike)
		{
    		if (itsDigitalCalibMode == BlackShift || itsDigitalCalibMode == NOCalib)
			{
				pDF = pSDF->toShiftedLN_1strike1shift(fwd,mat,(strike==-1?fwd:strike),0.);
				delete pSDF;
			}
			else if (itsDigitalCalibMode == GaussShift)
			{
				pDF = pSDF->toShiftedLN_1strike1shift(fwd,mat,(strike==-1?fwd:strike),1./delta);
				delete pSDF;
			}
			else if (itsDigitalCalibMode == RFStrike)
			{
				pDF = pSDF->toShiftedLN_2strikes(fwd,mat,(strike==-1?fwd:strike)-1e-4,(strike==-1?fwd:strike)+1e-4,fwd*(1./1.3-1),4./delta);
				delete pSDF;
			}
			else
				pDF = pSDF;
		}
		else if (itsCapFloorCalibMode == ATMStrike)
		{
			if (itsDigitalCalibMode == BlackShift || itsDigitalCalibMode == NOCalib)
			{
				pDF = pSDF->toShiftedLN_1strike1shift(fwd,mat,fwd,0.);
				delete pSDF;
			}
			else if (itsDigitalCalibMode == GaussShift)
			{
				pDF = pSDF->toShiftedLN_1strike1shift(fwd,mat,fwd,1./delta);
				delete pSDF;
			}
			else
				pDF = pSDF;

		}
		else
			pDF=pSDF;
	}

	ARM_VanillaSecurityDensity* pVSD	= new ARM_VanillaSecurityDensity(
											julianRD,
											julianSD,
											julianED,
											ARM_DensityFunctorPtr(pDF), 
											freq,
											itsCpnIndexDayCount,
											GETDEFAULTVALUE,
											1.,0.,1.,
											CreateClonedPtr(pCurve));

	return pVSD;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_TARNCalculator::GetIndexType()
{
    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string cpnIndexTerm(itsCpnIndexTerm);
    if(cpnIndexTerm=="12M")
        cpnIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += cpnIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Routine: ReplaceSwaptionDates
///	Returns: 
///	Action : Function to replace swaption reset dates with model dates
/////////////////////////////////////////////////////////////////
void ReplaceSwaptionDates(ARM_Vector* swaptionResetDates, const std::vector<double>& modelDates)
{
	size_t i;

	size_t modelIdx = 0;

	for (i = 0; i < swaptionResetDates->size(); ++i)
	{
		while ((modelIdx < modelDates->size()-1 ) && ((*modelDates)[modelIdx] < (*swaptionResetDates)[i]))
			modelIdx++;

		if ((modelIdx > 0) && (fabs((*modelDates)[modelIdx-1]-(*swaptionResetDates)[i]) < SwaptionDateTolerance))
			(*swaptionResetDates)[i] = (*modelDates)[modelIdx-1];

		if ((modelIdx > 0) && (fabs((*modelDates)[modelIdx]-(*swaptionResetDates)[i]) < SwaptionDateTolerance))
			(*swaptionResetDates)[i] = (*modelDates)[modelIdx];
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_TARNCalculator::CreateDiagonalSwaption()
{
	// Get the model dates
	const ARM_DateStripCombiner& dateStructure = DatesStructure();
	std::vector<double>& modelDates = dateStructure.GetDateStrip(0)->GetResetDates();

    /// Clean standard swap sets
    size_t i;
    for(i=0;i<itsStdSwaps.size();++i)
        delete itsStdSwaps[i].first;
    itsStdSwaps.resize(0);

    ARM_Swap* stdSwap;
    ARM_Swaption* swaption;

    list < ARM_Security* > swaptionList;
    double defaultStrike=0.0; // equivalent strike will be calculated later


    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names


    /// Find end date of diagonal swaps
    ARM_Date swapEndDate( atof(dealDesc.GetElem(nbEvents-1,IndexStartDate).c_str()) );
    swapEndDate.AddPeriod(itsCpnIndexTerm);

	ARM_INDEX_TYPE indexType = (ARM_INDEX_TYPE) FromIndexAndTermToIndexType(itsCpnIndexTerm, GetCurrencyUnit()->GetCcyName());
	int resetFreq = 1.0/StringMaturityToYearTerm(itsCpnIndexTerm);

    /// Analyse event dates
    for(size_t exerIdx=1;exerIdx<nbEvents-1;++exerIdx)
    {
		if (dealDesc.GetElemFormat(exerIdx,IndexStartDate) != ARM_MISSING_TYPE)
		{
			ARM_Date swapStartDate( atof(dealDesc.GetElem(exerIdx,IndexStartDate).c_str()) );

            /// The equivalent standard swap ends at the index standard swap end
            /// (except the case where pure fixed coupons appear at the end of the leg)

			// To Prevent stubs
			ARM_Date mathcSwapEndDate(swapStartDate);

			int nbMonths = floor((swapEndDate.GetJulian()-swapStartDate.GetJulian())/30.5+0.5);

			mathcSwapEndDate.AddMonths(nbMonths);

            stdSwap  = new ARM_Swap(
				swapStartDate,
                mathcSwapEndDate,
                indexType,0.0,1.0,K_RCV,resetFreq,resetFreq,GetCurrencyUnit());

			// Replace the swap reset dates with the model dates
			ReplaceSwaptionDates(stdSwap->GetResetDates(), modelDates);

            itsStdSwaps.push_back(pair<ARM_Swap*,int>(stdSwap,exerIdx));

            /// Build the market european swaption and set it
            /// (equivalent strike/vol will be computed latter)
            ARM_Date expiryDate((*(stdSwap->GetFloatLeg()->GetResetDates()))[0]);
            swaption = new ARM_Swaption(stdSwap,K_RCV,K_EUROPEAN,defaultStrike,expiryDate);
            swaptionList.push_back(static_cast< ARM_Security* >(swaption));
		}

    } /// for event date

    /// Built OSW portfolio and set default weight
    /// A default price is set to pass the validation test in CalibMethod constructor
    ARM_StdPortfolio* port = new ARM_StdPortfolio(swaptionList);
    for(i=0;i<port->size();++i)
    {
        port->SetWeight(OSW_DEFAULT_WEIGHT,i);
        port->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(port);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateTARNCapletFloorlet
///	Returns: nothing
///	Action : create the list of the TARN caplet & floorlet
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtrVector ARM_TARNCalculator::CreateTARNCapletFloorletAndDigital()
{
    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

	ARM_StdPortfolioPtrVector capFloorAndDigitalPortfolios;
	ARM_Digital* digital = NULL;
	
	bool sameLifeTimeCapFloor = (itsLifeTimeCapTarget == itsLifeTimeFloorTarget);

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Get the model name for coupon leg to get a caplet
    /// deal description with the right model name
    string cpnModelName = GetKeys()[itsCpnModelKey];


    size_t nbCF=0;
    double leverage,strike,capStrike;

	size_t i;

	for (i = 0; i < itsCaplets.size(); ++i)
		delete itsCaplets[i].first;
	itsCaplets.resize(0);

	std::vector<double> initBreakpointTimes;
	std::vector<double> initVolatility;
	std::vector<double> initBeta;

	bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

	int cpnIndexFreq=ARM_ArgConv_MatFrequency.GetNumber(itsCpnIndexTerm);

	int step = 1;

	if (cpnIndexFreq < itsCpnFreq)
	{
		step = itsCpnFreq/cpnIndexFreq;
	}

    /// Build a caplet/floorlet security using the coupon currency
    ARM_INDEX_TYPE liborType = GetIndexType();
    for(size_t eventIdx=(isArrearsCpn?2:1);eventIdx<nbEvents;eventIdx+=step)
    {
		leverage    =   atof( dealDesc.GetElem(eventIdx,Leverage).c_str() );
        strike      =   atof( dealDesc.GetElem(eventIdx,Strike).c_str() );
/*        cpnMin      =   atof( dealDesc.GetElem(eventIdx,CpnMin).c_str() );
        cpnMax      =   atof( dealDesc.GetElem(eventIdx,CpnMax).c_str() );
*/
        ARM_Date indexStartDate( atof(dealDesc.GetElem(eventIdx,IndexStartDate).c_str()) );

         /// Compute not adjusted index end date to avoid problem in caplet building
        ARM_Date indexEndDate(indexStartDate);
        indexEndDate.AddPeriod(itsCpnIndexTerm,itsCpnPayCal.c_str());

        /// Build the cap (ARM_Security & VanillaArg versions)
        capStrike = (strike>0 ? 100.0*strike : 0.0)/leverage;
        if(capStrike > K_NEW_DOUBLE_TOL)
        {	
            ARM_CapFloor capFloor(indexStartDate,indexEndDate,K_CAP,capStrike,liborType,
                0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

			double expiry = (capFloor.GetExpiryDate().GetJulian()-asOfDate.GetJulian());

			initBreakpointTimes.push_back(expiry);

			ARM_SwapLeg swapLeg(
				indexStartDate,
				indexEndDate,
				liborType,
				(int)K_RCV,
				(double)0.0,
				(int)K_DEF_FREQ,
				(int)K_DEF_FREQ,
				(int)K_ADVANCE,
				(int)K_ARREARS,
				GetCurrencyUnit()
				);

			ARM_ReferenceValue refValStrike(capStrike);
			ARM_ReferenceValue refValPayOff(1.0);

			if (eventIdx != nbEvents-step)
			{
				digital = new ARM_Digital(
					&swapLeg,
					K_CAP,
					&refValStrike,
					DigitalSpread1,
					DigitalSpread2,
					&refValPayOff);
			}

			itsCaplets.push_back(pair<ARM_CapFloor*,int>(static_cast<ARM_CapFloor*>(capFloor.Clone()),eventIdx));

            if(capFloor.GetResetDates()->size() != 1)
            {
	            CC_Ostringstream os;
	            os << ARM_USERNAME << " : Implied Caplet/Floorlet #" << eventIdx << " has more than one flow";
	            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
            }

			list< ARM_Security* > oneCapFloor;
			oneCapFloor.push_back(static_cast<ARM_CapFloor*>(capFloor.Clone()));
			capFloorAndDigitalPortfolios.push_back(ARM_StdPortfolioPtr(new ARM_StdPortfolio(oneCapFloor)));
			capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetWeight(CF_DEFAULT_WEIGHT,0);
			capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetPrice(((capFloorAndDigitalPortfolios.size()-1)/2+1)*CF_DEFAULT_PRICE,0);

			if (eventIdx != nbEvents-step)
			{
				list< ARM_Security* > oneDigital;
				oneDigital.push_back(static_cast<ARM_Security*>(digital));
				capFloorAndDigitalPortfolios.push_back(ARM_StdPortfolioPtr(new ARM_StdPortfolio(oneDigital)));
				capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetWeight(CF_DEFAULT_WEIGHT,0);
				capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetPrice(((capFloorAndDigitalPortfolios.size()-2)/2+1)*CF_DEFAULT_PRICE,0);
			}
			else
			{
				list< ARM_Security* > oneSecondCapFloor;
				ARM_CapFloor* secondCapFloor = static_cast<ARM_CapFloor*>(capFloor.Clone());
				// To prevent twice the same cap floor in the calib method
				secondCapFloor->SetStrike(0.0);
				oneSecondCapFloor.push_back(secondCapFloor);
				capFloorAndDigitalPortfolios.push_back(ARM_StdPortfolioPtr(new ARM_StdPortfolio(oneSecondCapFloor)));
				capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetWeight(CF_DEFAULT_WEIGHT,0);
				capFloorAndDigitalPortfolios[capFloorAndDigitalPortfolios.size()-1]->SetPrice(((capFloorAndDigitalPortfolios.size()-2)/2+2)*CF_DEFAULT_PRICE,0);
			}
			++nbCF;
        }

    } /// for event date

	ARM_CurveModelParam modelparamVol(ARM_ModelParamType::Volatility, SIGMA_DEFAULT_VALUE,"SIGMA",SIGMA_LOWER_BOUND,SIGMA_UPPER_BOUND );
	GetPricingModel()->GetModelParams()->SetModelParam(&modelparamVol);
	if (itsDigitalCalibMode != NOCalib)
	{
		ARM_CurveModelParam modelparamBeta(ARM_ModelParamType::Beta, BETA_DEFAULT_VALUE,"BETA",BETA_LOWER_BOUND,BETA_UPPER_BOUND );
		GetPricingModel()->GetModelParams()->SetModelParam(&modelparamBeta);
	}

	ARM_Portfolio* port=NULL;

	ARM_SFRM* refModel = GetRefModel();

	refModel->ConvertToShiftorBetaParam(*port);

	return capFloorAndDigitalPortfolios;

    
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetSubPrice
///	Returns: double
///	Action : compute the analytical price of a sub selection of
///          the deal description
/////////////////////////////////////////////////////////////////
double ARM_TARNCalculator::GetSubPrice(int startRowIdx,int endRowIdx,TARNColAlias columnName,
                                      const string& evalDateStr,
                                      const ARM_DealDescription& dealDesc) const
{
    /// Loop for each row between start and end row to set asOfDate because the
    /// GP can't price a GS with all event dates set to a unique date (asOf here)
    double price=0.0;
    ARM_DealDescriptionPtr genDesc;
    ARM_GenSecurity* genLeg;
    ARM_GenPricer* genPricer;
    for(size_t i=startRowIdx;i<=endRowIdx;++i)
    {
	    genDesc=dealDesc.GetSubDescription(i,i,columnName+1);
        genDesc->SetElem(1,0,evalDateStr,ARM_DATE_TYPE);
        genLeg = new ARM_GenSecurity(genDesc);
	    genPricer = new ARM_GenPricer( genLeg, &*GetPricingModel() );

        ARM_AutoCleaner< ARM_GenSecurity > HoldGS(genLeg);
        ARM_AutoCleaner< ARM_GenPricer > HoldGP(genPricer);

        /// Compute the analytical leg price (eventDate was changed to asOfDate)
        price += genPricer->Price();
    }
    return price;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputeEquivalentDatas
///	Returns: double
///	Action : compute the equivalent datas of a diagonal swaption
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::ComputeEquivalentDatas(const ARM_DealDescription& dealDesc,
    const pair<ARM_Swap*,int>& stdSwap,
    ARM_Swaption* swaption, const string& evalDateStr, std::vector<double>& equivDatas)
{
    ARM_Swap* swap=stdSwap.first;
    int startIdx = stdSwap.second;

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    const ARM_DealDescription dealDescr = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDescr.GetRowsNb(); /// + 1 because of column names

    // Compute the annuity price
    double annuityPrice = swap->GetFixedLeg()->ComputePrice();

	double swapRate  = swap->CptMarketSwapRate();
	double fixedRate = swapRate;

	if ( itsModelType==ARM_PricingModelType::SFRM2F )
	{
		/// Compute the analytical (eventDate is changed to asOfDate)
		/// fixed leg price (strike & funding spread)

		TARNColAlias columnNb;

		if (itsCapFloorCalibMode == ExerStrike)
		{
			columnNb = StdFixFlow;
		}
		else
		{
			columnNb = StdFixRFFlow;
		}

		double fixLegPrice = GetSubPrice(startIdx,nbEvents-1,columnNb,evalDateStr,dealDesc);

		fixedRate = fixLegPrice / annuityPrice * 100;
	}
	
	/// Take care : use UpdateStrike because SetStrike() only set strike and not
    /// the decomp strike used pricing !!
    swaption->UpdateStrike(fixedRate);

    /// Compute target price & vega
    double price=swaption->ComputePrice();
    double vega=swaption->ComputeSensitivity(K_VEGA);

    double optMat = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
    double swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;

    equivDatas[OSW_SWAP_RATE]=swapRate/100.0;
    equivDatas[OSW_TARGET_PRICE]=price;
    equivDatas[OSW_TARGET_VEGA]=vega;

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::Calibrate()
{
	switch (itsModelType)
	{
	case ARM_PricingModelType::SFRM2F:
		Calibrate_SFRM2F();
		break;
	case ARM_PricingModelType::SBGM:
		Calibrate_SBGM();
		break;
	case ARM_PricingModelType::HK:
		Calibrate_HK();
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::Calibrate : wrong model type !" );
	}

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::Calibrate_SFRM2F()
{
	/// Volatility bootstrapping with possible MRS optimisation
    if(itsCapFloorCalibMode != NOCalib)
    {
        //1.- functional calibration
		ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());
		if(model)
			GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
		else
			GetCalibMethod()->Calibrate(&*GetPricingModel());

    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Calibrate_HK
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::Calibrate_HK()
{
	ARM_NumMethodPtr mcMethod = GetPricingModel()->GetNumMethod();
	
	int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber("CN1F");
	ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType);

	int nbSteps   = (itsNbIterations.size()>4)?itsNbIterations[4]:PDE_NB_STEPS;
	int gridSize  = (itsNbIterations.size()>5)?itsNbIterations[5]:PDE_GRIDSIZE;
	double nbStdDevs = (itsNbIterations.size()>6)?itsNbIterations[6]:PDE_NB_STDEV;

	ARM_PDEMethod* pPDEMethod = new ARM_PDEMethod(numScheme,nbSteps, gridSize, 2*nbStdDevs);
	GetPricingModel()->SetNumMethod( ARM_NumMethodPtr( pPDEMethod ) );
	
	ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());
	if(model)
		GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
	else
		GetCalibMethod()->Calibrate(&*GetPricingModel());
    
	GetPricingModel()->SetNumMethod(mcMethod);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Calibrate_SBGM
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::Calibrate_SBGM()
{
	//1.- functional calibration
	ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());
	if(model)
		GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
	else
		GetCalibMethod()->Calibrate(&*GetPricingModel());

	//2.- correl calibration
	if(itsOSWCalibFlag)
	{
		bool isFreezeWeights = false;
		bool isArrearsCpn = (itsCpnTiming==K_ARREARS);

		ARM_StdPortfolioPtr pfPtr(CreateDiagonalSwaption());

		ARM_CalibMethod* correlCalib = new ARM_CalibMethod(pfPtr,ARM_ModelParamVector(),ARM_CalibMethodType::Bootstrap1D,
                                 ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

		SetCalibMethod( ARM_CalibMethodPtr( correlCalib ) );
		ComputeDiagonalSwaptionPrice(isFreezeWeights,isArrearsCpn);

		/// Replace the beta param with new initialisation
		std::vector<double> initTimes(1,0.);
		std::vector<double> initBetas(1,0.1);
		std::vector<double> betaLowerBound(1,0.);
		std::vector<double> betaUpperBound(1,4.);

		ARM_CurveModelParam* beta = NULL;
			
		beta = new ARM_CurveModelParam(
				ARM_ModelParamType::BetaCorrelation,
				&initBetas,
				&initTimes,
				"BETA",
				"STEPUPRIGHT",
				&betaLowerBound,
				&betaUpperBound);
			
		GetCalibMethod()->GetCalibParams().push_back(beta);
		
		ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());
		if(model)
			GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
		else
			GetCalibMethod()->Calibrate(&*GetPricingModel());
	}
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the TARN deal.
/////////////////////////////////////////////////////////////////
double ARM_TARNCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());
	ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);

	ARM_GenPricer* genPricer = NULL;
	if( itsControlVariableFlag )
	{
		ARM_StringVector cvColumnNames( 2);
		cvColumnNames[0] = TARNColNamesTable[ ControlVariable ];
		cvColumnNames[1] = TARNColNamesTable[ ControlVariable2 ];
		std::vector<double> columnPrices( 2 );
		ComputeControlVariable( columnPrices );
		itsCVRefPrices = ARM_GP_VectorPtr( new std::vector<double>(columnPrices) );
		string refPriceColumn = TARNColNamesTable[ TARN ];
		genPricer = new ARM_GenPricer( &* GetGenSecurity(),&*GetPricingModel(), cvColumnNames, columnPrices, refPriceColumn );
	}
	else
	{
		genPricer = new ARM_GenPricer( &*GetGenSecurity(),&*GetPricingModel() );
	}

	ARM_StringVector columnNames = CreateColumnNames();

	if (columnNames.size() > 0)
	{
		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
		genPricer->Price();
		double price,stdDev;
		ARM_GP_VectorPtr ra;

		for (int i = 0; i < NbProductsToPrice; ++i)
		{
			if (itsProductsToPrice[i])
			{
				price	= genPricer->GetPricerInfo()->GetContents( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
				stdDev	= genPricer->GetPricerInfo()->GetContents( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] ).GetData("StdDev").GetDouble();
				ra = genPricer->GetPricerInfo()->GetContents( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] ).GetData("RunningAvg").GetVector();

 				if(i == TARNPrice)
				{
					if( itsControlVariableFlag )
					{
						price			= genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("Price").GetDouble();
						stdDev			= genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("StdDev").GetDouble();
						ra	= genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("RunningAvg" ).GetVector();
					}
					else
					{
						ra	= genPricer->GetPricerInfo()->GetContents( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] ).GetData("RunningAvg" ).GetVector();
					}
					itsTARNPrice	=  price;
					itsTARNStdDev	= stdDev;
					itsTARNRA		= ra;
				}
				else if(i == SwapPrice)
				{
					itsSwapPrice = price;
					itsSwapStdDev = stdDev;
					itsSwapRA = ra;
				}
				else if(i == LifeTimeCapPrice)
				{
					itsLifeTimeCapPrice = price;
					itsLifeTimeCapStdDev = stdDev;
					itsLifeTimeCapRA = ra;
				}
				else if(i == LastLifeTimeCapPrice)
				{
					itsLastLifeTimeCapPrice = price;
					itsLastLifeTimeCapStdDev = stdDev;
					itsLastLifeTimeCapRA = ra;
				}
				else if(i == LifeTimeFloorPrice)
				{
					itsLifeTimeFloorPrice = price;
					itsLifeTimeFloorStdDev = stdDev;
					itsLifeTimeFloorRA = ra;
				}
				else if(i == DigitalFundingPrice)
				{
					itsDigitalFundingPrice = price;
					itsDigitalFundingStdDev = stdDev;
					itsDigitalFundingRA = ra;
				}
				else if(i == DigitalPrice)
				{
					itsDigitalPrice = price;
					itsDigitalStdDev = stdDev;
					itsDigitalRA = ra;
				}
				else if(i == FundingPrice)
				{
					itsFundingPrice = price;
					itsFundingStdDev = stdDev;
					itsFundingRA = ra;
				}
				else if(i == ExerciseTimeAverage)
				{
					itsExerciseTimeAverage = price;
				}
				else if(i == ExerciseStrikesPrice || i == ExerciseProbasPrice)
				{
					ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] ).GetData("IntermediatePrices");
					const ARM_DealDescription& newDealDesc = GetGenSecurity()->GetDealDescription();

					bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

					ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

					size_t newItermPricesSize = itermPrices->size();
					if( newItermPricesSize && isArrearsCpn ) 
						--newItermPricesSize;
					ARM_VectorPtr  newItermPrices( new std::vector<double>(newItermPricesSize ) );

					for (size_t j = (isArrearsCpn?1:0); j < itermPrices->size(); ++j)
					{
						double df = 1.0;

						if (i == ExerciseStrikesPrice)
						{
							double indexEndDate = atof( newDealDesc.GetElem(j+1,IndexEndDate).c_str() );
							df = zeroCurve->DiscountPrice((indexEndDate-asOfDate.GetJulian())/K_YEAR_LEN);
						}

						(*newItermPrices)[j-(isArrearsCpn?1:0)] = (*itermPrices)[j]/df;
					}
					
					itermPrices = newItermPrices;

					if (i == ExerciseStrikesPrice)
						itsExerciseStrikes = itermPrices;
					else if	(i == ExerciseProbasPrice)
						itsExerciseProbas = itermPrices;
				}
			}
		}
		itsHasBeenPriced = true;
	}


    return(itsTARNPrice);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_TARNCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_TARNCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[TARNPrice])
	{
		GetPricingData()[ "TARNPrice"		] = itsTARNPrice;
		GetPricingData()[ "TARNStdDev"		] = itsTARNStdDev;
		GetPricingData()[ "TARNRunningAvg"	] = itsTARNRA;
	}
	if (itsProductsToPrice[SwapPrice])
	{
		GetPricingData()[ "SwapPrice"		] = itsSwapPrice;
		GetPricingData()[ "SwapStdDev"		] = itsSwapStdDev;
		GetPricingData()[ "SwapRunningAvg"	] = itsSwapRA;
	}
	if (itsProductsToPrice[LifeTimeCapPrice])
	{
		GetPricingData()[ "LifeTimeCapPrice"		] = itsLifeTimeCapPrice;
		GetPricingData()[ "LifeTimeCapStdDev"	] = itsLifeTimeCapStdDev;
		GetPricingData()[ "LifeTimeCapRunningAvg"	] = itsLifeTimeCapRA;
	}
	if (itsProductsToPrice[LastLifeTimeCapPrice])
	{
		GetPricingData()[ "LastLifeTimeCapPrice"		] = itsLastLifeTimeCapPrice;
		GetPricingData()[ "LastLifeTimeCapStdDev"	] = itsLastLifeTimeCapStdDev;
		GetPricingData()[ "LastLifeTimeCapRunningAvg"	] = itsLastLifeTimeCapRA;
	}
	if (itsProductsToPrice[LifeTimeFloorPrice])
	{
		GetPricingData()[ "LifeTimeFloorPrice"	] = itsLifeTimeFloorPrice;
		GetPricingData()[ "LifeTimeFloorStdDev"	] = itsLifeTimeFloorStdDev;
		GetPricingData()[ "LifeTimeFloorRunningAvg"	] = itsLifeTimeFloorRA;
	}
	if (itsProductsToPrice[DigitalFundingPrice])
	{
		GetPricingData()[ "DigitalFundingPrice"	] = itsDigitalFundingPrice;
		GetPricingData()[ "DigitalFundingStdDev"	] = itsDigitalFundingStdDev;
		GetPricingData()[ "DigitalFundingRunningAvg"] = itsDigitalFundingRA;
	}
	if (itsProductsToPrice[DigitalPrice])
	{
		GetPricingData()[ "DigitalPrice"	] = itsDigitalPrice;
		GetPricingData()[ "DigitalStdDev"	] = itsDigitalStdDev;
		GetPricingData()[ "DigitalRunningAvg"	] = itsDigitalRA;
	}
	if (itsProductsToPrice[FundingPrice])
	{
		GetPricingData()[ "FundingPrice"			] = itsFundingPrice;
		GetPricingData()[ "FundingStdDev"		] = itsFundingStdDev;
		GetPricingData()[ "FundingRunningAvg"		] = itsFundingRA;
	}
	if (itsProductsToPrice[ExerciseStrikesPrice])
		GetPricingData()[ "ExerciseStrikes"		] = itsExerciseStrikes;
	if (itsProductsToPrice[ExerciseProbasPrice])
		GetPricingData()[ "ExerciseProbas"		] = itsExerciseProbas;
	if (itsProductsToPrice[ExerciseTimeAverage])
		GetPricingData()[ "ExerciseTimeAverage"	] = itsExerciseTimeAverage;

	if( itsCVRefPrices != ARM_GP_VectorPtr(NULL) )
		GetPricingData()[ "CVRefPrices"	] = itsCVRefPrices;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputeTARNCapletFloorletPrices
///	Returns: nothing
///	Action : compute market target prices of the C/F portfolio
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::ComputeTARNCapletFloorletPrices(
	bool isFreezeWeights, 
	bool isInitParam,
	bool isArrearCpn,
	CalibrationMode capFloorCalibMode,
	CalibrationMode digitalCalibMode,
	bool computeSwaptionPrice)
{
	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* cfBSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    double capFloorPrice,vega, digitalPrice, secondCapFloorPrice;
	double expiry, tenor, fwdRate, vol, lifeTimeCap, strike, cfStrike, leverage, theta;
    ARM_CapFloor* capFloor;
	ARM_CapFloor* secondCapFloor;
	ARM_Digital* digital;
	ARM_CalibMethod* cfCalibMethod = GetCFCalibMethod(computeSwaptionPrice);
	ARM_CalibMethod* digitCalibMethod = NULL;
	ARM_StdPortfolioPtr capFloorPF, digitPF;
	
	ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();

	int eventIdx;
	int i = itsCaplets.size()-1;

	// Set the exercise strikes
    while (cfCalibMethod != NULL)
    {
		eventIdx = itsCaplets[i].second;

		if (digitalCalibMode != NOCalib)
		{
			digitCalibMethod = cfCalibMethod;;
		
			if (cfCalibMethod->GetlinkedMethod())
			{
				cfCalibMethod = cfCalibMethod->GetlinkedMethod();
				digitPF = digitCalibMethod->GetPortfolio();
				if ((i != itsCaplets.size()-1))
				{
					digital=static_cast< ARM_Digital* >(digitPF->GetAsset(0));
					digital->SetModel(cfBSModel);
				}
				else if (capFloorCalibMode == ExerStrike)
				{
					secondCapFloor=static_cast< ARM_CapFloor* >(digitPF->GetAsset(0));	
					secondCapFloor->SetModel(cfBSModel);
				}
			}
		}

		capFloorPF = cfCalibMethod->GetPortfolio();
		capFloor=static_cast< ARM_CapFloor* >(capFloorPF->GetAsset(0));
		capFloor->SetModel(cfBSModel);
		

		expiry = (capFloor->GetExpiryDate().GetJulian()-asOfDate.GetJulian());
		lifeTimeCap      =   atof( dealDesc.GetElem(eventIdx,LifeTimeCapTarget).c_str() );
		strike      =   atof( dealDesc.GetElem(eventIdx,Strike).c_str() );
		leverage      =   atof( dealDesc.GetElem(eventIdx,Leverage).c_str() );
		theta = atof( dealDesc.GetElem(eventIdx,IT).c_str() );

		cfStrike		=   atof( dealDesc.GetElem(eventIdx,EstimExerciseStrike).c_str() );

		// The cap floor strike should not be less than MIN_CF_STRIKE to prevent calibration
		// instability
		double rfStrike;

		// For Tarn SnowBall
		if (i == 0)
			rfStrike = (strike - GetCoefReverse() * GetC0())/leverage;
		else
			rfStrike = strike/leverage;

		if (capFloorCalibMode == RFStrike)
		{
			capFloor->SetStrike((rfStrike<MIN_CF_STRIKE ? MIN_CF_STRIKE : rfStrike)*100.0);
		}
		else if (capFloorCalibMode == ExerStrike)
		{
			capFloor->SetStrike((cfStrike<MIN_CF_STRIKE ? MIN_CF_STRIKE : cfStrike)*100.0);
		}
		else if (capFloorCalibMode == ATMStrike)
		{
			capFloor->SetStrike(capFloor->GetSwapLeg()->GetFwdRates()->Elt(0)<MIN_CF_STRIKE/100. ? MIN_CF_STRIKE*100. : capFloor->GetSwapLeg()->GetFwdRates()->Elt(0));
		}
		if (digitalCalibMode != NOCalib)
		{
			if ((i != itsCaplets.size()-1))
			{
				if (digitalCalibMode == RFStrike)
				{
					digital->SetConstantStrike((rfStrike<MIN_CF_STRIKE ? MIN_CF_STRIKE : rfStrike)*100.0);
				}
				else if (digitalCalibMode == ExerStrike)
				{
					digital->SetConstantStrike((cfStrike<MIN_CF_STRIKE ? MIN_CF_STRIKE : cfStrike)*100.0);
				}
				else if (digitalCalibMode == ATMStrike)
				{
					digital->SetConstantStrike(digital->GetSwapLeg()->GetFwdRates()->Elt(0)<MIN_CF_STRIKE*100. ? MIN_CF_STRIKE : capFloor->GetSwapLeg()->GetFwdRates()->Elt(0));
				}
				digitalPrice=digital->ComputePrice();
			}
			else if (capFloorCalibMode == ExerStrike)
			{
				secondCapFloor->SetStrike((rfStrike<MIN_CF_STRIKE ? MIN_CF_STRIKE : rfStrike)*100.0);
				secondCapFloorPrice = secondCapFloor->ComputePrice();
			}
		}

        vega=capFloor->ComputeSensitivity(K_VEGA);
        vega=fabs(vega);
		capFloorPrice=capFloor->ComputePrice();


		if (isInitParam)
		{
			/// Test if volatility must be initialise. In this case,
			/// its schedule contents the expiry date of each swaption product
			size_t sigmaIdx,betaIdx,paramSize=cfCalibMethod->GetCalibParams().size();
			for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
				if(cfCalibMethod->GetCalibParam(sigmaIdx) &&
				(cfCalibMethod->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
					break;
			if (digitalCalibMode != NOCalib)
			{
				for(betaIdx=0;betaIdx<paramSize;++betaIdx)
					if(digitCalibMethod->GetCalibParam(betaIdx) &&
					(digitCalibMethod->GetCalibParams())[betaIdx]->GetType() == ARM_ModelParamType::Beta)
						break;
			}

			tenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();
			if (capFloor->GetSwapLeg()->GetFwdRates() 
				&& capFloor->GetSwapLeg()->GetFwdRates()->GetSize() >= 1)
			{
				fwdRate = capFloor->GetSwapLeg()->GetFwdRates()->Elt(0);
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : cap floor fwd rates are not available !" );
			}
			strike = capFloor->GetStrike();
			vol = cfBSModel->ComputeVol(expiry/K_YEAR_LEN,tenor,fwdRate,strike);

			/// Replace the sigma param with new initialisations
			std::vector<double> initTimes(1,expiry);
			std::vector<double> initSigmas(1,vol/100.0);
			std::vector<double> sigmaLowerBound(1,SIGMA_LOWER_BOUND);
			std::vector<double> sigmaUpperBound(1,SIGMA_UPPER_BOUND);

			ARM_CurveModelParam* sigma = new ARM_CurveModelParam(
				ARM_ModelParamType::Volatility,
				&initSigmas,
				&initTimes,
				"SIGMA",
				"STEPUPRIGHT",
				&sigmaLowerBound,
				&sigmaUpperBound);

			/// Replace the sigma param with new initialisations
			std::vector<double> initBetas(1,BETA_DEFAULT_VALUE);
			std::vector<double> betaLowerBound(1,BETA_LOWER_BOUND);
			std::vector<double> betaUpperBound(1,BETA_UPPER_BOUND);

			ARM_CurveModelParam* beta = NULL;
			
			if ((digitalCalibMode != NOCalib) && ((i != itsCaplets.size()-1) || (capFloorCalibMode == ExerStrike)))
			{
				beta = new ARM_CurveModelParam(
					ARM_ModelParamType::Beta,
					&initBetas,
					&initTimes,
					"BETA",
					"STEPUPRIGHT",
					&betaLowerBound,
					&betaUpperBound);
			}

			if(sigmaIdx >= paramSize || paramSize == 0)
				cfCalibMethod->GetCalibParams().push_back(sigma);
			else
			{
				delete cfCalibMethod->GetCalibParam(sigmaIdx);
				(cfCalibMethod->GetCalibParams())[sigmaIdx] = sigma;
			}

			if ((digitalCalibMode != NOCalib) && ((i != itsCaplets.size()-1) || (capFloorCalibMode == ExerStrike)))
			{
				if(betaIdx >= paramSize || paramSize == 0)
					digitCalibMethod->GetCalibParams().push_back(beta);
				else
				{
					delete digitCalibMethod->GetCalibParam(betaIdx);
					(digitCalibMethod->GetCalibParams())[betaIdx] = beta;
				}
			}
		}

        capFloorPF->SetPrecision(vega,0);
        capFloorPF->SetPrice(capFloorPrice,0);

		if (digitalCalibMode != NOCalib)
		{
			if (i != itsCaplets.size()-1)
			{
				digitPF->SetPrecision(DIGITAL_DEFAULT_WEIGHT,0);
				digitPF->SetPrice(digitalPrice,0);
			}
			else if (capFloorCalibMode == ExerStrike)
			{
				digitPF->SetPrecision(DIGITAL_DEFAULT_WEIGHT,0);
				digitPF->SetPrice(secondCapFloorPrice,0);
			}
		}
	
		if (digitalCalibMode != NOCalib)
		{
			cfCalibMethod = digitCalibMethod->GetPreviousMethod();
		}
		else
		{
			cfCalibMethod = cfCalibMethod->GetPreviousMethod();
		}

		i--;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: GetRefModel
///	Returns: ARM_PricingModelPtr
///	Action : return the reference model of the TARN
/////////////////////////////////////////////////////////////////
ARM_SFRM* ARM_TARNCalculator::GetRefModel()
{
	ARM_SFRM* refModel;
    ARM_HybridBasisFwdIR* bfirModel=NULL;
    if( string(GetCurrencyUnit()->GetCcyName()) == string(GetFundingCcy().GetCcyName()) )
    {
        /// Allow update only for default model (H&W1F)
        refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());
        if( !refModel )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only H&W1F model is allowed for updating");
    }
    else
    {
        bfirModel = dynamic_cast< ARM_HybridBasisFwdIR* >(&*GetPricingModel());
        refModel = dynamic_cast< ARM_SFRM* >(&*(bfirModel->GetRefModel()));
        if( !bfirModel || !refModel)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only Hybrid BFIR model with a H&W1F model as reference is allowed for updating");
    }

	return refModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputeMCExerciseStrikes
///	Returns: nothing
///	Action : compute exercise strikes with Monte Carlo
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::ComputeMCExerciseStrikes()
{
	/// Get the zero curve of the calculator
    ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	/// Make a copy of the initial description to update strike spreads
    ARM_DealDescriptionPtr dealDesc( (ARM_DealDescription*) const_cast< ARM_DealDescription& >(GetGenSecurity()->GetDealDescription()).Clone() );

    /// Select the right column that describes the product
    TARNColAlias prodIdx;
    prodIdx = ExerciseTime;

    /// Build the sub-generic security
    size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
    size_t endRowIdx = GetGenSecurity()->GetDealDescription().GetRowsNb() - 1;
	ARM_StringVector refColumns(2);
	refColumns[0] = TARNColNamesTable[ExerciseStrike];
	refColumns[1] = TARNColNamesTable[ExerciseTime];
    ARM_DealDescriptionPtr subDealDesc = GetGenSecurity()->GetDealDescription().GetSubDescription(startRowIdx,endRowIdx,prodIdx+1,refColumns);
    ARM_GenSecurityPtr genSec(new ARM_GenSecurity(subDealDesc));

	ARM_GenPricer* genPricer = new ARM_GenPricer( &*genSec,&*GetPricingModel() );

    ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );
	CalibrateAndTimeIt();

	ARM_SFRM* refModel = GetRefModel();

	// For Exercise strike we use 10000 iterations maximum
	int nbIter = itsNbIterations[1] > MaxIterCalibStrike? MaxIterCalibStrike: itsNbIterations[1];
	ARM_NumMethodPtr oldNumMethod = refModel->GetNumMethod();
	ARM_MCMethod* oldMcMethod = dynamic_cast<ARM_MCMethod*>( &*oldNumMethod );

	if( !oldMcMethod )
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : TARN should be priced with MC.");

	ARM_TimeStepPerYearScheduler scheduler(MC_NB_INTER_STEPS);
	ARM_NormalCentredSamplerND sampler(&scheduler);

	ARM_MCMethod* newMcMethod = new ARM_MCMethod(nbIter,oldMcMethod->GetRandGen(),&sampler);
	refModel->SetNumMethod(ARM_NumMethodPtr( newMcMethod ) );
    double price = genPricer->Price();
	
	/// restore old numerical method
	refModel->SetNumMethod( oldNumMethod );

	ARM_PricerInfo* pricerInfo = genPricer->GetPricerInfo();
	ARM_VectorPtr intermediatePrices = pricerInfo->GetContents(TARNColNamesTable[ExerciseStrike]).GetData("IntermediatePrices").GetVector();

	// In the in arrears case we should not care about the first line
	size_t i;
	for (i = 0; i < intermediatePrices->size(); ++i)
	{
		double df = 1.0;

		/// Caplet/Floorlet selection
        if( dealDesc->GetElemFormat(i+1,IndexEndDate) != ARM_MISSING_TYPE )
		{
			double indexEndDate = atof( dealDesc->GetElem(i+1,IndexEndDate).c_str() );
			df = zeroCurve->DiscountPrice((indexEndDate-asOfDate.GetJulian())/K_YEAR_LEN);
		}
    
		CC_Ostringstream exerStrikeStr;
		exerStrikeStr << (*intermediatePrices)[i]/df;
		dealDesc->SetElem(i+1,EstimExerciseStrike,exerStrikeStr.str(),ARM_DOUBLE);
	}

	/// last thing, remove the swap after the 
	if( dealDesc->GetRowsNb()>=3 )
	{
		double exerciseTime = pricerInfo->GetContents(TARNColNamesTable[ExerciseTime]).GetData("Price").GetDouble();
		double firstEvent = atof( dealDesc->GetElem(1,0).c_str() );
		double period = (atof( dealDesc->GetElem(2,0).c_str() ) - firstEvent)/K_YEAR_LEN;

		for( i=1; i<dealDesc->GetRowsNb()-1; ++i )
		{
			if( ( atof( dealDesc->GetElem(i+1,0).c_str() ) - firstEvent )/K_YEAR_LEN > exerciseTime + period/2 )
				break;
		}
		size_t knockOutIndex= i;
		SetCVSwapMiddleText( dealDesc, knockOutIndex );
	}
	
	/// Rebuild the GenSec to update nodes of the syntaxic tree
	SetGenSecurity( ARM_GenSecurityPtr( new ARM_GenSecurity( dealDesc, GetGenSecurity()->GetPayModelName() ) ) );
}


void ARM_TARNCalculator::SetCVSwapMiddleText( ARM_DealDescriptionPtr& dealDesc, size_t knockOutIndex ) const
{
	for( size_t i=knockOutIndex; i<dealDesc->GetRowsNb()-1; ++i )
		dealDesc->SetElem(i+1,ControlVariable2,"0",ARM_DOUBLE);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputeControlVariable
///	Returns: double
///	Action : computes the price of the control variable (before compute the estimated strikes)
/////////////////////////////////////////////////////////////////

void ARM_TARNCalculator::ComputeControlVariable(std::vector<double>& cvprices)
{
	/// get the exercise strikes
	ComputeMCExerciseStrikes();

	/// Get the zero curve of the calculator
	bool isArrearsCpn = (itsCpnTiming==K_ARREARS);
	const ARM_DealDescription& dealDes = GetGenSecurity()->GetDealDescription();
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());
	size_t endRowIdx = dealDes.GetRowsNb() - 1;
	double controlVariablePrice  = GetSubPrice((isArrearsCpn?2:1),endRowIdx,ControlVariable,evalDateStr,dealDes );
	double controlVariable2Price = GetSubPrice((isArrearsCpn?2:1),endRowIdx,ControlVariable2,evalDateStr,dealDes );
	cvprices[0] = controlVariablePrice;
	cvprices[1] = controlVariable2Price;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: ComputeDiagonalSwaptionPrice
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::ComputeDiagonalSwaptionPrice(bool isFreezeWeights, bool isArrearCpn)
{
    /// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

    const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();

	/// Get asOfDate
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());


    ARM_VolCurve* oswVolCurve = oswBSModel->GetVolatility();
    ARM_YCModel* ycModel = oswBSModel->GetYCModel();

    ARM_StdPortfolioPtr swaptionPF = GetOSWPortfolio();
    size_t i,nbOSW=swaptionPF->GetSize();

	double price,vega,nominal,swapRate,weight;
    ARM_Swaption* swaption;
    std::vector<double> equivDatas(OSW_NB_EQUIVDATA);
    bool isNotSelected;
    for(i=0;i<nbOSW;++i)
    {
        swaption=static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
	    swaption->SetModel(oswBSModel);
        itsStdSwaps[i].first->SetModel(ycModel);
      
		ComputeEquivalentDatas(dealDesc,itsStdSwaps[i],swaption,evalDateStr,equivDatas);

        swapRate    = equivDatas[OSW_SWAP_RATE];
        price       = equivDatas[OSW_TARGET_PRICE];
        vega        = equivDatas[OSW_TARGET_VEGA];

        nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight=(*(swaptionPF->GetWeights()))[i];

        isNotSelected = vega < VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        swaptionPF->SetWeight(isNotSelected ? 0.0 : weight,i);
        swaptionPF->SetPrecision(0.0001,i);
        swaptionPF->SetPrice(price,i);
		swaption->SetModel(NULL);
    }
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: PrepareMeanRevForSwaptionCalib
///	Returns: nothing
///	Action : 
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::PrepareMeanRevForSwaptionCalib()
{	
	ARM_CalibMethod* oswCalibMethod = GetOSWCalibMethod();
    std::vector<double> initTimes(1,0.0);
	std::vector<double> initMRS(1,MRS_DEFAULT_VALUE);
	std::vector<double> mrsLowerBound(1,MRS_LOWER_BOUND);
	std::vector<double> mrsUpperBound(1,MRS_UPPER_BOUND);

	ARM_CurveModelParam* mrs = new ARM_CurveModelParam(
		ARM_ModelParamType::MeanReversion,
		&initMRS,
		&initTimes,
		"MEANREVERSION",
		"STEPUPRIGHT",
		&mrsLowerBound,
		&mrsUpperBound);

	size_t mrsIdx,paramSize=oswCalibMethod->GetCalibParams().size();
    for(mrsIdx=0;mrsIdx<paramSize;++mrsIdx)
    {
        if( oswCalibMethod->GetCalibParam(mrsIdx) &&
            (oswCalibMethod->GetCalibParams())[mrsIdx]->GetType() == ARM_ModelParamType::MeanReversion )
            break;
    }

	if(mrsIdx >= paramSize || paramSize == 0)
		oswCalibMethod->GetCalibParams().push_back(mrs);
    else
    {
        delete oswCalibMethod->GetCalibParam(mrsIdx);
        (oswCalibMethod->GetCalibParams())[mrsIdx] = mrs;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateModel_SFRM2F
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateModel()
{
	itsHasBeenPriced = false;
	switch (itsModelType)
	{
	case ARM_PricingModelType::SFRM2F:
		UpdateModel_SFRM2F();
		break;
	case ARM_PricingModelType::SBGM:
		UpdateModel_SBGM();
		break;
	case ARM_PricingModelType::HK:
		UpdateModel_HK();
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::UpdateModel : wrong model type !" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateModel_SBGM
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateModel_SBGM()
{
	CreateAndSetModel();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateModel_HK
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateModel_HK()
{
	CreateAndSetModel();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateModel_SFRM2F
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateModel_SFRM2F()
{

	ARM_ModelParam* mrsParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }
	ARM_ModelParam* betaParam = NULL;	
	betaParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaKey]));
	if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Beta Param for key=" << GetKeys()[BetaKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}
	ARM_ModelParam* correlParam = dynamic_cast<ARM_ModelParam*>(GetMktDataManager()->GetData(GetKeys()[CorrelKey]));
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : Correl Param for key=" << GetKeys()[CorrelKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

    /// Check for model consistency
    ARM_SFRM* refModel;
    ARM_HybridBasisFwdIR* bfirModel=NULL;
    if( string(GetCurrencyUnit()->GetCcyName()) == string(GetFundingCcy().GetCcyName()) )
    {
        /// Allow update only for default model (H&W1F)
        refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());
        if( !refModel )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only SFRM model is allowed for updating");
    }
    else
    {
        bfirModel = dynamic_cast< ARM_HybridBasisFwdIR* >(&*GetPricingModel());
        refModel = dynamic_cast< ARM_SFRM* >(&*(bfirModel->GetRefModel()));
        if( !bfirModel || !refModel)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only Hybrid BFIR model with a H&W1F model as reference is allowed for updating");
    }

	ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
	ARM_ZeroCurve* basisCurve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));

    /// Update yield curves
    if(bfirModel)
    {
	    ARM_Forex* forex   = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        bfirModel->UpdateCurves( CreateClonedPtr( cpnCurve) ,CreateClonedPtr( fundingCurve ),CreateClonedPtr( basisCurve ),*forex);
    }
    else
	    refModel->SetZeroCurve(CreateClonedPtr( cpnCurve ) );

	ARM_ModelParamVector paramVector(4);

	paramVector[0] = static_cast<ARM_ModelParam*>(&refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility));

    /// Update the MRS
	paramVector[1] = mrsParam;
	/// Update the Beta
	if (itsDigitalCalibMode != NOCalib)
	{
		paramVector[2] = static_cast<ARM_ModelParam*>(&refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Beta));
	}
	else
	{
		paramVector[2] = betaParam;
	}
	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetIndexType();

	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,SFRM_VOL_TYPE);

	refModel->SetModelParams(*pSFRMModelParams);
	delete pSFRMModelParams;

	ARM_Portfolio* port=NULL;
	refModel->ConvertToShiftorBetaParam(*port);

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateCalibration(bool isUpdateStrike)
{
	itsHasBeenPriced = false;
	if(IsBasis())
		itsFundSpread.SetOrdinates(ComputeDomesticBasis());

    switch (itsModelType)
	{
	case ARM_PricingModelType::SFRM2F:
		UpdateCalibration_SFRM2F();
		break;
	case ARM_PricingModelType::SBGM:
		UpdateCalibration_SBGM();
		break;
	case ARM_PricingModelType::HK:
		UpdateCalibration_HK();
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_TARNCalculator::UpdateCalibration : wrong model type !" );
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateCalibration_SFRM2F
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateCalibration_SFRM2F(bool isUpdateStrike)
{
    /// Update cap floor prices
    bool isFreezeWeights = true;					// keep the portfolio size to avoid hedge jumps
    bool isInitSigma = true;						// keep current sigma param to init calibration
    bool isArrearsCpn = (itsCpnTiming==K_ARREARS);	// In Arrears Case

    if(itsCapFloorCalibMode != NOCalib)
	{
		ComputeTARNCapletFloorletPrices(
			isFreezeWeights,
			isInitSigma,
			isArrearsCpn,
			itsCapFloorCalibMode,
			itsDigitalCalibMode,
			itsOSWCalibFlag);
	}

	/// Compute diagonal swaption target prices
    if(itsOSWCalibFlag)
    {
        ComputeDiagonalSwaptionPrice(isFreezeWeights,isArrearsCpn);
		PrepareMeanRevForSwaptionCalib();
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateCalibration_SBGM
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateCalibration_SBGM(bool isUpdateStrike)
{
 	CreateAndSetCalibration_SBGM();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: UpdateCalibration_HK
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_TARNCalculator::UpdateCalibration_HK(bool isUpdateStrike)
{
  	CreateAndSetCalibration_HK();
}

////////////////////////////////////////////////////
///	Class   : ARM_TARNCalculator
///	Routines: Clone,View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////
ARM_Object* ARM_TARNCalculator::Clone() const
{
	return new ARM_TARNCalculator(*this);
}

void ARM_TARNCalculator::View(char* id, FILE* ficOut) const
{
    FILE* fOut;
    char fOutName[200];
	
    if ( ficOut == NULL )
    {
		ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);
		(void) unlink(fOutName);
		fOut = fopen(fOutName, "w");
    }
    else
		fOut = ficOut;

    /// TARN Calculator specific datas viewing
    fprintf(fOut,"\n\n =======> TARN REVERSE FLOATER CALCULATOR <====== \n");

    CC_Ostringstream tarnData;
    tarnData << "\nStartDate = " <<  itsStartDate.toString() << "\n";
    tarnData << "EndDate = " << itsEndDate.toString() << "\n";
    tarnData << "Pay/Rec = " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << "\n";


    tarnData << "\nReverse Leg Datas :\n";
    tarnData << "StartDate = " <<  itsFixEndDate.toString() << "\n";
    tarnData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsCpnDayCount ) << "\n";
    tarnData << "Frequency (Reset & Pay) = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    tarnData << "Pay Calendar = " << itsCpnPayCal << "\n";
    tarnData << "Reset Timing = " << ARM_ParamView::GetMappingName(S_TIMING_MOD, itsCpnTiming ) << "\n";
	
	//Customised reset dates.
	if (GetCustomResetFlag() == true)
	{
		std::vector<double> customResetDates = GetCustomResetDates();
		ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
		for (int i=0; i< customResetDates.size();i++)
		{
			if (asOfDate.GetJulian() < customResetDates.Elt(i))
			{
				tarnData << "ResetDate[" << i+1 << "]  :" << ARM_Date(customResetDates.Elt(i)).toString() << "\n";
			}
		}
	}
	else
	{
		tarnData << "Reset Gap = " << itsCpnResetGap << "\n";
	}
	tarnData << "Reset Calendar = " << itsCpnResetCal << "\n";
    tarnData << "Cpn Index Term = " << itsCpnIndexTerm << "\n";
	tarnData << "Cpn Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    tarnData << "Cpn Index Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT, itsCpnIndexDayCount )  << "\n";

    tarnData << "\nFunding Leg Datas :\n";
    tarnData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFundDayCount ) << "\n";
    tarnData << "Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsFundFreq ) << "\n";
    

    if( GetGenSecurity() != ARM_GenSecurityPtr(NULL))
    {
        const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
        size_t i,nbRows = dealDesc.GetRowsNb();

        tarnData << "\nReverse Coupon Profile :\n";
        tarnData << "StartDate \tStrike    \tLeverage  \tLifeTimeCapTarget   \tLifeTimeFloorTarget \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,StartDate) != ARM_MISSING_TYPE)
            {
                tarnData << ARM_Date(atof(dealDesc.GetElem(i,StartDate).c_str())).toString() << "\t";
                tarnData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Strike).c_str());
                if(dealDesc.GetElemFormat(i,Leverage) != ARM_MISSING_TYPE)
                {
                    tarnData << "\t" << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Leverage).c_str()) << "\t";
                    tarnData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,LifeTimeCapTarget).c_str()) << "\t";
                    tarnData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,LifeTimeFloorTarget).c_str()) << "\n";
                }
                else
                    tarnData << "\n";
            }
        }

        tarnData << "\nFunding Spread Profile :\n";
        tarnData << "StartDate \tSpread   \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,FundingStartDate) != ARM_MISSING_TYPE)
            {
                tarnData << ARM_Date(atof(dealDesc.GetElem(i,FundingStartDate).c_str())).toString() << "\t";
                tarnData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,FundingSpread).c_str()) << "\n";
            }
        }

        tarnData << "\nNominal Profile :\n";
        tarnData << "PayDate   \tNominal   \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,PayDate) != ARM_MISSING_TYPE)
            {
                tarnData << ARM_Date(atof(dealDesc.GetElem(i,PayDate).c_str())).toString() << "\t";
                tarnData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Nominal).c_str()) << "\n";
            }
        }

        
    }

    	/// part common to gencalculator
	tarnData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s",tarnData.str().c_str());

	string cfTxt("\n\nAuto-calibrated caplet & floorlet portfolio (Sigma Calibration)\n");
    cfTxt += "Cap Floor Calibration=";
    cfTxt += ARM_ArgConvReverse_TARNCalibMode.GetString(itsCapFloorCalibMode);

	string dgtTxt("\n\nAuto-calibrated digital portfolio (Beta Calibration)\n");
	dgtTxt += "Digital Calibration=";
	dgtTxt += ARM_ArgConvReverse_TARNCalibMode.GetString(itsDigitalCalibMode);

	string capFloorPFStr;
	string digitalPFStr;

    if( itsCapFloorCalibMode != NOCalib && itsModelType==ARM_PricingModelType::SFRM2F )
    {
		ARM_CalibMethod* cfCalibMethod = GetCFCalibMethod(itsOSWCalibFlag);
		ARM_CalibMethod* digitCalibMethod = NULL;
		ARM_StdPortfolioPtr capFloorPF;
		ARM_StdPortfolioPtr digitPF;

		int size = itsCaplets.size();
		size_t i = size;

		while (cfCalibMethod != NULL)
		{
			bool withDigitCalibMethod = false;

			if (itsDigitalCalibMode!=NOCalib)
			{
				digitCalibMethod = NULL;
				if (cfCalibMethod->GetlinkedMethod())
				{
					withDigitCalibMethod = true;
					digitCalibMethod = cfCalibMethod;
					cfCalibMethod = cfCalibMethod->GetlinkedMethod();
					digitPF = digitCalibMethod->GetPortfolio();
				}
			
			}

			capFloorPF = cfCalibMethod->GetPortfolio();

			capFloorPFStr = capFloorPF->toStringByAsset(0, i) + capFloorPFStr;

			if ((itsDigitalCalibMode!=NOCalib) && digitCalibMethod)
			{
					digitalPFStr = digitPF->toStringByAsset(0, i) + digitalPFStr;
			}

			if ((itsDigitalCalibMode!=NOCalib) && withDigitCalibMethod)
				cfCalibMethod = digitCalibMethod->GetPreviousMethod();
			else
				cfCalibMethod = cfCalibMethod->GetPreviousMethod();
			--i;
		}

		cfCalibMethod = GetCFCalibMethod(itsOSWCalibFlag);
		capFloorPF = cfCalibMethod->GetPortfolio();
		capFloorPFStr = capFloorPF->toStringHeader(size) + capFloorPFStr;
		cfTxt += capFloorPFStr;
		if (itsDigitalCalibMode!=NOCalib)
		{
			digitalPFStr = digitPF->toStringHeader(size) + digitalPFStr;
			dgtTxt += digitalPFStr;
		}
    }

	fprintf(fOut,"%s\n",cfTxt.c_str());
	fprintf(fOut,"%s\n",dgtTxt.c_str());
	
	string oswTxt("\n\nAuto-calibrated diagonal swaption portfolio (MRS calibration)\n");
    oswTxt += "Calibration=";
    oswTxt += itsOSWCalibFlag ? "On" : "Off";
    fprintf(fOut,"%s\n",oswTxt.c_str());
    if(itsOSWCalibFlag)
    {
        ARM_StdPortfolioPtr diagonalSwaptionPF=GetOSWPortfolio();
        if(diagonalSwaptionPF!=ARM_StdPortfolioPtr(NULL))
            diagonalSwaptionPF->View(id,fOut);
    }

    string prodTxt("Product to be priced = ");
    if(IsTARNToPrice())						prodTxt += "Trigger Accrual Range Note";
	else if(IsLifeTimeCapToPrice())			prodTxt += "Life Time Cap";
	else if(IsLifeTimeFloorToPrice())		prodTxt += "Life Time Floor";
	else if(IsDigitalFundingToPrice())		prodTxt += "Digital Funding";
	else if(IsFundingToPrice())				prodTxt += "DigitalFunding";
    else prodTxt += "Unknown !";
    fprintf(fOut,"\n\n %s\n",prodTxt.c_str());

    /// Common viewing
    ARM_GenCalculator::View(id,fOut);

    if ( ficOut == NULL )
       fclose(fOut);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNCalculator
///	Routine: CreateColumnNames
///	Returns: 
///	Action :  Create columnNames
/////////////////////////////////////////////////////////////////
ARM_StringVector ARM_TARNCalculator::CreateColumnNames()
{
	ARM_StringVector columnNames;
	columnNames.reserve(NbProductsToPrice);
	size_t i;

    for (i = 0; i < NbProductsToPrice; ++i)
	{
		if (itsProductsToPrice[i])
			columnNames.push_back( TARNColNamesTable[ ARM_TARNCalculator::ProductToPriceColumns[i] ] );
	}

	if( itsControlVariableFlag )
	{
		columnNames.push_back( TARNColNamesTable[ ControlVariable ] );
		columnNames.push_back( TARNColNamesTable[ ControlVariable2 ] );
	}

	return columnNames;
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

