/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file callablesnowballcalculator.h
 *
 *  \base callable snowball
 *	\author  JP Riaudel, Richard Guillemot
 *	\version 1.0
 *	\date March 2005
 */


/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/callablesnowballcalculator.h"
#include "gpcalculators/basisconverter.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/singleton.h"
#include "gpbase/autocleaner.h"
#include "gpbase/gpmatrixlinalg.h"
#include "gpbase/stringconvert.h"
#include "gpbase/gplinalgconvert.h"

/// gpinfra
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/gensecmanipulator.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/correlmatparam.h"
#include "gpinfra/modelnamemap.h"


/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillapricer.h"
#include "gpcalib/modelfitterdes.h"

/// gpmodels
#include "gpmodels/multiassets.h"
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/smiledfrm.h"
#include "gpmodels/modelparamsmf.h"
#include "gpmodels/markovfunctional.h"
#include "gpmodels/ForwardMarginBasis.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/argconvdefault.h"

/// gpnummethods
#include "gpnummethods/amcmethod.h"
#include "gpnummethods/amc_andersen.h"
#include "gpnummethods/amc_ls.h"
#include "gpnummethods/cfmethod.h"
#include "gpnummethods/normalcentredsampler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/argconvdefault.h"
#include "gpnummethods/pdemethod.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"

/// kernel
#include <glob/paramview.h>
//#include <inst/swaption.h>
//#include <inst/sumopt.h>
//#include <mod/bssmiled.h>
//#include <inst/forex.h>


CC_BEGIN_NAMESPACE( ARM )


/// Reference schedules for CRF date structure
const unsigned int CPN_SCHED=0;
const unsigned int EXER_SCHED   =1;
const unsigned int NB_CALLABLESB_SCHED =2;


const double CF_DEFAULT_PRICE=1.0e+100;
const double CF_DEFAULT_WEIGHT=1.0;

const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

/// SFRM sigma range [1bp,50000bp] with a 2000bp default value
const double SIGMA_DEFAULT_VALUE    = 0.20;
const double SIGMA_LOWER_BOUND      = 0.0001;
const double SIGMA_UPPER_BOUND      = 0.5;


/// SFRM MRS range [-15%,50%] with a 2% default value
const double BETA_LOWER_BOUND        = 0.1;
const double BETA_UPPER_BOUND        = 1.5;
const double BETA_DEFAULT_VALUE      = 1;

const double MIN_STRIKE_VALUE		 = 0.0025;

const double FILTER_STRIKE_TOL		= 1e-8;
const double ATM_VOL_FLAG			= -1;

const double NON_CALL_FEES			= 1e12;

/// SBGM Beta Correl
const double BETA_CORREL_LOWER_BOUND        = 0;
const double BETA_CORREL_UPPER_BOUND        = 4;
const double BETA_CORREL_DEFAULT_VALUE      = 0.1;

// DEFAULT for Summit/Calypco
const double MRS_DEFAULT_VALUE		= 0.0;
const double CORREL_DEFAULT_VALUE	= 0.0;
const double HUMP_DEFAULT_VALUE		= 1.0;
const double RECORREL_DEFAULT_VALUE = 0.0;
	
// default PDE params for HK/SBGM
const int PDE_NB_STEPS		= 800;
const int PDE_NB_STDEV		= 5;
const int PDE_GRIDSIZE		= 901;

// default params for SBGM
const int SBGM_NB_FACTORS	= 5;


/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string BETA_KEY_NAME          = "BETA_";
const string CORREL_KEY_NAME        = "CORREL_";
const string HUMP_KEY_NAME          = "HUMP_";
const string BETA_CORREL_KEY_NAME   = "BETACORREL_";
const string YC_BASIS_KEY_NAME      = "YC_BASIS_";
const string FOREX_KEY_NAME         = "FOREX_";
const string UNKNOWN_KEY_NAME       = "UNKNOWN";

const string ARM_CallableSnowBallCalculator::CallableSnowBallColNamesTable [] =
{
        "EventDate",
		"Index",
		"SumOptStartDate",
		"EndDate",
        "CpnFwdStartDate",
        "CpnFwdEndDate",
        "CpnPayDate",
        "CpnIT",
        "FundingStartDate",
        "FundingEndDate",
        "NextCpnFwdStartDate",
        "NextCpnFwdEndDate",
        "NextCpnPayDate",
        "NextCpnIT",
        "PrevCpnFwdStartDate",
        "PrevCpnFwdEndDate",
		"ProductEndDate",
		"DFPay",
		"DFNum",
        "DFNextPay",
		"Nominal",
        "NextNominal",
        "Fees",
		"Const",
		"LevPrevCpn",
		"LevNewOpt",
		"StrikeOpt",
		"MinCpn",
		"MaxCpn",
		"NextConst",
		"NextLevPrevCpn",
		"NextLevNewOpt",
		"NextStrikeOpt",
		"NextMinCpn",
		"NextMaxCpn",
		"PrevConst",
		"PrevLevPrevCpn",
		"PrevLevNewOpt",
		"PrevStrikeOpt",
		"PrevMinCpn",
		"PrevMaxCpn",
		"CpnIndex",
		"Coupon",
		"CouponSimple",
		"CouponOption",
		"StrikeSum",
		"PaidCoupon",
		"LBCoupon",
		"LBCouponPaid",
		"UBCoupon",
		"UBCouponPaid",
		"CouponOpt",
		"CouponOptPaid",
		"LBCouponAvg",
		"LBCouponAvgPaid",
		"AnaLBCouponAvg",
		"AnaLBCoupon",
		"CouponOptAvg",
		"CouponOptAvgPaid",
		"PrevCpnIndex",
		"PrevCoupon",
		"PrevCouponSimple",
		"PrevCouponOpt",
		"NextCpnIndex",
		"NextSimpleCoupon",
		"NextCpnIndexWithoutAdj",
		"NextSimpleCouponWithoutAdj",
		"StrikeCap",
		"NextCouponCap",
		"StrikeFloor",
		"NextCouponFloor",
		"NextCouponOpt",
		"NextCoupon",
		"FundingCoeff",
		"FundingMargin",
		"Funding",
		"FundingPaid",
		"CF",
		"Option",
		"Bermuda",
		"CF1",
		"CF2",
		"AMCIndex",
		"Option2",
		"Bermuda2"
};


const int ARM_CallableSnowBallCalculator::Product1ToPriceColumns[] =
{
	CF,
	FundingPaid,
	Bermuda,
	LBCouponPaid,
	UBCouponPaid,
	CouponOptPaid,
	LBCouponAvgPaid,
	CouponOptAvgPaid
};

const int ARM_CallableSnowBallCalculator::Product2ToPriceColumns[] =
{
	CF1,
	FundingPaid,
	Bermuda2,
	LBCouponPaid,
	UBCouponPaid,
	CouponOptPaid,
	LBCouponAvgPaid,
	CouponOptAvgPaid
};

const int ARM_CallableSnowBallCalculator::AnaProductToPriceColumns[] =
{
	AnaLBCouponAvg,
	AnaLBCoupon
};

// SFRM default factors number
const int SFRM_NB_FACTORS			= 2;
const int SFRM_VOL_TYPE				= K_DIAG;


/////////////////////////////////////////////////////////////////
///	Class  : 
///	Routine: FindIdxInVector
///	Returns: bool
///	Action : finds idx in vector
/////////////////////////////////////////////////////////////////
bool FindIdxInVector( const std::vector<double>& timeVector, double time, double daysNb, size_t& idx)
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
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Constructor
///	Returns: void
///	Action : ...
/////////////////////////////////////////////////////////////////
ARM_CallableSnowBallCalculator::ARM_CallableSnowBallCalculator(const ARM_Date& startDate,
															   const ARM_Date& endDate,
															   int cpnFreq,
															   int cpnDaycount,
															   int cpnIndexTiming,
															   const string& cpnIdxTerm,
															   int cpnIndexDaycount,
															   const string& cpnResetCal,
															   const string& cpnPayCal,
															   int cpnIntRule,
															   int cpnIndexResetLag,
															   int payRec,
															   int CF,
															   const ARM_Curve& notionalProfile,
															   const ARM_Curve& fundnotionalProfile,
															   const ARM_Curve& constProfile,
															   const ARM_Curve& lPrevCpnProfile,
															   const ARM_Curve& lNewOptProfile,
															   const ARM_Curve& strikeOptProfile,
															   const ARM_Curve& minCpnProfile,
															   const ARM_Curve& maxCpnProfile,
															   int fundFreq,
															   int fundDaycount,
															   const ARM_Curve& fundCoeffProfile,
															   const ARM_Curve& fundMarginProfile,
															   int NotifDays,
															   int NonCall,
															   const string& exerciseCal,
															   bool callSBOrStacky,
															   const ARM_Curve& feesProfile,
															   int NbPathBounding,
															   int NbPathPricing,
															   int NbMaxBucket,
															   int fixSteps,
															   const string& USMethod,
															   CalibMode calibMode,
															   bool betaCalib,
															   bool fixBoundary,
															   bool fixBeta,
															   ControlVariableMode controlVariableFlag,
															   bool fixControlVariable,
															   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
															   const ARM_MarketData_ManagerRep& mktDataManager,
															   const ARM_StringVector& mdmKeys,
															   const string& genType1,
															   const string& genType2,
															   int firstNbTimes,
															   int firstNbDims,
															   const string& pathScheme,
															   const string& pathOrder,
															   bool	antithetic,
															   ARM_ModelType modelType,
															   TriggerMode triggerMode,
															   int calibSwoptFreq,
															   const string& regressors,
															   const vector< double >& hkData)
:	ARM_GenCalculator(mktDataManager),
		itsStartDate(startDate),
		itsEndDate(endDate),
		itsCpnFreq(cpnFreq),
		itsCpnDayCount(cpnDaycount),
		itsCpnIndexTiming(cpnIndexTiming),
		itsCpnIdxTerm(cpnIdxTerm),
		itsCpnIndexDayCount(cpnIndexDaycount),
		itsCpnResetCal(cpnResetCal),
		itsCpnPayCal(cpnPayCal),
		itsCpnIntRule(cpnIntRule),
		itsCpnResetLag(-abs(cpnIndexResetLag)),
		itsPayRec(payRec),
		itsCapOrFloor(CF),
		itsNotionalProfile(notionalProfile),
		itsFundNotionalProfile(fundnotionalProfile),
		itsConstProfile(constProfile),
		itsLPrevCpnProfile(lPrevCpnProfile),
		itsLNewOptProfile(lNewOptProfile),
		itsStrikeOptProfile(strikeOptProfile),
		itsMinCpnProfile(minCpnProfile),
		itsMaxCpnProfile(maxCpnProfile),
		itsFundFreq(fundFreq),
		itsFundDayCount(fundDaycount),
		itsFundCoeffProfile(fundCoeffProfile),
		itsFundMarginProfile(fundMarginProfile),
		itsInitFundMarginProfile(fundMarginProfile),
		itsNoticeDays(-abs(NotifDays)),
		itsNbNonCall(NonCall),
		itsExerCal(exerciseCal),
		itsCallSBOrStacky(callSBOrStacky),
		itsFeesProfile(feesProfile),
		itsNbPathBounding(NbPathBounding),
		itsNbPathPricing(NbPathPricing),
		itsNbMaxBucket(NbMaxBucket),
		itsFixSteps(fixSteps),
		itsUSMethod(USMethod),
		itsCalibMode(calibMode),
		itsBetaCalib(betaCalib),
		itsFixBoundary(fixBoundary),
		itsFixBeta(fixBeta),
		itsControlVariableFlag(controlVariableFlag),
		itsFixControlVariable(fixControlVariable),
		itsProductsToPrice(productsToPrice),
		itsHasBeenPriced(false),
		itsFirstCallIdx(0),
		itsSnowBallSJPrice(0.0),
		itsSnowBallSJStdDev(0.0),
		itsSnowBallSJRA(0),
		itsFundingPrice(0.0),
		itsFundingStdDev(0.0),
		itsFundingRA(0),
		itsSnowBallCallablePrice(0.0),
		itsSnowBallCallableStdDev(0.0),
		itsSnowBallCallableRA(0),
		itsNoOptSBPrice(0.0),
		itsNoOptSBStdDev(0.0),
		itsNoOptSBRA(0),
		itsStackySBPrice(0.0),
		itsStackySBStdDev(0.0),
		itsStackySBRA(0),
		itsCouponOptPrice(0.0),
		itsCouponOptStdDev(0.0),
		itsCouponOptRA(0),
		itsLBCouponAvgPrice(0.0),
		itsLBCouponAvgStdDev(0),
		itsLBCouponAvgRA(0),
	    itsCouponOptAvgPrice(0.0),
		itsCouponOptAvgStdDev(0.0),
		itsCouponOptAvgRA(0),
		itsIsCpnFlow(0),
		itsResetDatesVec(0),
		itsStartDatesVec(0),
		itsEndDatesVec(0),
		itsPeriodsVec(0),
		itsPaymentDatesVec(0),
		itsConstVec(0),
		itsLevPrevVec(0),
		itsLevNewVec(0),
		itsStrikeOptVec(0),
		itsCpnMinVec(0),
		itsCpnMaxVec(0),
		itsStrikeSum(0),
		itsFeesVec(0),
		itsvCpnNominal(0),
		itsvInitialFundNominal(0),
		itsvInitialFundSpread(0),
		itsGenType1(genType1),
		itsGenType2(genType2),
		itsFirstNbTimes(firstNbTimes),
		itsFirstNbDims(firstNbDims),
		itsPathScheme(pathScheme),
		itsPathOrder(pathOrder),
		itsAntithetic(antithetic),
		itsModelType (modelType),
		itsTriggerMode(triggerMode),
		itsHkData (hkData),
		itsCalibSwoptFreq(calibSwoptFreq),
		itsRegressors(regressors),
		itsFixBetaTimes(0),
		itsFixBetaValues(0),
		itsFundDateStrip(0),
		itsStructDateStrip(0)

{
	SetName(ARM_CALLABLE_SNOWBALL);

	 /// Set keys for MDM datas access
    if(mdmKeys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(mdmKeys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey]		  = newMdMKeys[YcKey];
        newMdMKeys[YcBasisDomKey]     = newMdMKeys[YcKey];
		newMdMKeys[YcBasisFundKey]    = newMdMKeys[YcKey];
        newMdMKeys[ForexKey]          = UNKNOWN_KEY_NAME;
        SetKeys(newMdMKeys);
    }
    else
        SetKeys(mdmKeys);
    
    /// Check input datas
    CheckDataAndTimeIt();

    /// Set the coupon/payment currency (inherited from ARM_Security)
    ARM_Currency* cpnCcy = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]))->GetCurrencyUnit();
    SetCurrencyUnit(cpnCcy);

	/// Set funding & basis currencies
    ARM_Currency fundCcy   = *static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[FundingKey]))->GetCurrencyUnit();
    SetFundingCcy(fundCcy);

	// Precalcul sur le datestrip
	itsDateStructure = InitDateStrip();

	// We need this price for the control variate
	if (itsControlVariableFlag)
	{
		itsProductsToPrice[NbProductsToPrice+AnaLBCouponPrice] = true;
		itsProductsToPrice[NoOptSBPrice] = true;

		itsCVBetas = ARM_VectorPtr(new std::vector<double>(0));
	}


    /// Create the Generic Security paid in coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcBasisDomKey],CreateColumnNames(),CreateCstManager(),itsFixBoundary,(itsControlVariableFlag!=NO));


    /// Create the SFRM pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping and
    /// strike spread adjustments
    CreateAndSetCalibrationAndTimeIt();


}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Constructor
///	Returns: void
///	Action : ...
/////////////////////////////////////////////////////////////////
ARM_CallableSnowBallCalculator::ARM_CallableSnowBallCalculator(ARM_Currency& ccy,
															   const ARM_Date& startDate,
															   const ARM_Date& endDate,
															   int cpnFreq,
															   int cpnDaycount,
															   int cpnIndexTiming,
															   const string& cpnIdxTerm,
															   int cpnIndexDaycount,
															   const string& cpnResetCal,
															   const string& cpnPayCal,
															   int cpnIntRule,
															   int cpnIndexResetLag,
															   int payRec,
															   int CF,
															   const ARM_Curve& notionalProfile,
															   const ARM_Curve& constProfile,
															   const ARM_Curve& lPrevCpnProfile,
															   const ARM_Curve& lNewOptProfile,
															   const ARM_Curve& strikeOptProfile,
															   const ARM_Curve& minCpnProfile,
															   const ARM_Curve& maxCpnProfile,
															   int fundFreq,
															   int fundDaycount,
															   const ARM_Curve& fundCoeffProfile,
															   const ARM_Curve& fundMarginProfile,
															   int NotifDays,
															   int NonCall,
															   const string& exerciseCal,
															   bool callSBOrStacky,
															   const ARM_Curve& feesProfile,
															   int NbPathBounding,
															   int NbPathPricing,
															   int NbMaxBucket,
															   int fixSteps,
															   const string& USMethod,
															   CalibMode calibMode,
															   bool betaCalib,
															   bool fixBoundary,
															   bool fixBeta,
															   ControlVariableMode controlVariableFlag,
															   bool fixControlVariable,
															   const std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/& productsToPrice,
															   const string& genType1,
															   const string& genType2,
															   int firstNbTimes,
															   int firstNbDims,
															   const string& pathScheme,
															   const string& pathOrder,
															   bool	antithetic,
															   ARM_ModelType modelType,
															   TriggerMode triggerMode,
															   int calibSwoptFreq,
															   const string& regressors,
															   const vector< double >& hkData)
:	ARM_GenCalculator(),
		itsStartDate(startDate),
		itsEndDate(endDate),
		itsCpnFreq(cpnFreq),
		itsCpnDayCount(cpnDaycount),
		itsCpnIndexTiming(cpnIndexTiming),
		itsCpnIdxTerm(cpnIdxTerm),
		itsCpnIndexDayCount(cpnIndexDaycount),
		itsCpnResetCal(cpnResetCal),
		itsCpnPayCal(cpnPayCal),
		itsCpnIntRule(cpnIntRule),
		itsCpnResetLag(-abs(cpnIndexResetLag)),
		itsPayRec(payRec),
		itsCapOrFloor(CF),
		itsNotionalProfile(notionalProfile),
		itsFundNotionalProfile(),
		itsConstProfile(constProfile),
		itsLPrevCpnProfile(lPrevCpnProfile),
		itsLNewOptProfile(lNewOptProfile),
		itsStrikeOptProfile(strikeOptProfile),
		itsMinCpnProfile(minCpnProfile),
		itsMaxCpnProfile(maxCpnProfile),
		itsFundFreq(fundFreq),
		itsFundDayCount(fundDaycount),
		itsFundCoeffProfile(fundCoeffProfile),
		itsFundMarginProfile(fundMarginProfile),
		itsInitFundMarginProfile(fundMarginProfile),
		itsNoticeDays(-abs(NotifDays)),
		itsNbNonCall(NonCall),
		itsExerCal(exerciseCal),
		itsCallSBOrStacky(callSBOrStacky),
		itsFeesProfile(feesProfile),
		itsNbPathBounding(NbPathBounding),
		itsNbPathPricing(NbPathPricing),
		itsNbMaxBucket(NbMaxBucket),
		itsFixSteps(fixSteps),
		itsUSMethod(USMethod),
		itsCalibMode(calibMode),
		itsBetaCalib(betaCalib),
		itsFixBoundary(fixBoundary),
		itsFixBeta(fixBeta),
		itsControlVariableFlag(controlVariableFlag),
		itsFixControlVariable(fixControlVariable),
		itsProductsToPrice(productsToPrice),
		itsHasBeenPriced(false),
		itsFirstCallIdx(0),
		itsSnowBallSJPrice(0.0),
		itsSnowBallSJStdDev(0.0),
		itsSnowBallSJRA(0),
		itsFundingPrice(0.0),
		itsFundingStdDev(0.0),
		itsFundingRA(0),
		itsSnowBallCallablePrice(0.0),
		itsSnowBallCallableStdDev(0.0),
		itsSnowBallCallableRA(0),
		itsNoOptSBPrice(0.0),
		itsNoOptSBStdDev(0.0),
		itsNoOptSBRA(0),
		itsStackySBPrice(0.0),
		itsStackySBStdDev(0.0),
		itsStackySBRA(0),
		itsCouponOptPrice(0.0),
		itsCouponOptStdDev(0.0),
		itsCouponOptRA(0),
		itsLBCouponAvgPrice(0.0),
		itsLBCouponAvgStdDev(0),
		itsLBCouponAvgRA(0),
	    itsCouponOptAvgPrice(0.0),
		itsCouponOptAvgStdDev(0.0),
		itsCouponOptAvgRA(0),
		itsIsCpnFlow(0),
		itsResetDatesVec(0),
		itsStartDatesVec(0),
		itsEndDatesVec(0),
		itsPeriodsVec(0),
		itsPaymentDatesVec(0),
		itsConstVec(0),
		itsLevPrevVec(0),
		itsLevNewVec(0),
		itsStrikeOptVec(0),
		itsCpnMinVec(0),
		itsCpnMaxVec(0),
		itsStrikeSum(0),
		itsFeesVec(0),
		itsvCpnNominal(0),
		itsvInitialFundNominal(0),
		itsvInitialFundSpread(0),
		itsGenType1(genType1),
		itsGenType2(genType2),
		itsFirstNbTimes(firstNbTimes),
		itsFirstNbDims(firstNbDims),
		itsPathScheme(pathScheme),
		itsPathOrder(pathOrder),
		itsAntithetic(antithetic),
		itsModelType (modelType),
		itsTriggerMode(triggerMode),
		itsHkData (hkData),
		itsCalibSwoptFreq(calibSwoptFreq),
		itsRegressors(regressors),
		itsFixBetaTimes(0),
		itsFixBetaValues(0),
		itsFundDateStrip(0),
		itsStructDateStrip(0)

{
	SetName(ARM_CALLABLE_SNOWBALL);

    /// Set keys for MDM datas access
    /// Check input datas
    CheckDataAndTimeIt();

    /// Set the coupon/payment currency (inherited from ARM_Security)
	SetCurrencyUnit(&ccy);

	ARM_StringVector newMdMKeys(NbKeys);
	string ccyName				= (string) (ccy.GetCcyName());
	newMdMKeys[YcKey]			= YC_KEY_NAME + ccyName;
	newMdMKeys[OswModelKey]		= OSWMODEL_KEY_NAME + ccyName;
	newMdMKeys[CfModelKey]		= CFMODEL_KEY_NAME + ccyName;
	newMdMKeys[MrsKey]			= MRS_KEY_NAME + ccyName;
	newMdMKeys[BetaKey]			= BETA_KEY_NAME + ccyName;
	newMdMKeys[CorrelKey]		= CORREL_KEY_NAME + ccyName;
	newMdMKeys[HumpKey]			= HUMP_KEY_NAME + ccyName;
	newMdMKeys[BetaCorrelKey]	= BETA_CORREL_KEY_NAME + ccyName;
    newMdMKeys[FundingKey]		= newMdMKeys[YcKey];
    newMdMKeys[YcBasisDomKey]   = newMdMKeys[YcKey];
	newMdMKeys[YcBasisFundKey]  = newMdMKeys[YcKey];
    newMdMKeys[ForexKey]		= UNKNOWN_KEY_NAME;
	SetKeys(newMdMKeys);


	// Precalcul sur le datestrip
	itsDateStructure = InitDateStrip();


	// We need this price for the control variate
	if (itsControlVariableFlag)
	{
		itsProductsToPrice[NbProductsToPrice+AnaLBCouponPrice] = true;
		itsProductsToPrice[NoOptSBPrice] = true;

		itsCVBetas = ARM_VectorPtr(new std::vector<double>(0));
	}


    /// Create the Generic Security paid in coupon currency
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcKey],CreateColumnNames(),CreateCstManager(),itsFixBoundary,(itsControlVariableFlag!=NO));


    /// Create the SFRM pricing model with its default parameters
//    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping and
    /// strike spread adjustments
//    CreateAndSetCalibrationAndTimeIt();


}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CallableSnowBallCalculator::~ARM_CallableSnowBallCalculator()
{
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CallableSnowBallCalculator::ARM_CallableSnowBallCalculator(const ARM_CallableSnowBallCalculator& rhs)
:	ARM_GenCalculator(rhs),
	itsStartDate(rhs.itsStartDate),
	itsEndDate(rhs.itsEndDate),
    itsPayRec(rhs.itsPayRec),
	itsCpnDayCount(rhs.itsCpnDayCount),
	itsCpnFreq(rhs.itsCpnFreq),
	itsCpnIndexTiming(rhs.itsCpnIndexTiming),
	itsCpnIdxTerm(rhs.itsCpnIdxTerm),
	itsCpnIndexDayCount(rhs.itsCpnIndexDayCount),
	itsCpnResetCal(rhs.itsCpnResetCal),
	itsCpnPayCal(rhs.itsCpnPayCal),
	itsCpnResetLag(rhs.itsCpnResetLag),
	itsCpnIntRule(rhs.itsCpnIntRule),
	itsNotionalProfile(rhs.itsNotionalProfile),
	itsFundNotionalProfile(rhs.itsFundNotionalProfile),
	itsConstProfile(rhs.itsConstProfile),
	itsLPrevCpnProfile(rhs.itsLPrevCpnProfile),
	itsLNewOptProfile(rhs.itsLNewOptProfile),
	itsStrikeOptProfile(rhs.itsStrikeOptProfile),
	itsMinCpnProfile(rhs.itsMinCpnProfile),
	itsMaxCpnProfile(rhs.itsMaxCpnProfile),
	itsCapOrFloor(rhs.itsCapOrFloor),
	itsFundMarginProfile(rhs.itsFundMarginProfile),
	itsInitFundMarginProfile(rhs.itsInitFundMarginProfile),
	itsFundFreq(rhs.itsFundFreq),
	itsFundDayCount(rhs.itsFundDayCount),
	itsNoticeDays(rhs.itsNoticeDays),
	itsExerCal(rhs.itsExerCal),
	itsFeesProfile(rhs.itsFeesProfile),
	itsNbNonCall(rhs.itsNbNonCall),
	itsSBRank2Pos(rhs.itsSBRank2Pos),
	itsSBPos2Rank(rhs.itsSBPos2Rank),
	itsExerRank2Pos(rhs.itsExerRank2Pos),
	itsExerPos2Rank(rhs.itsExerPos2Rank),
	itsDateStructure(rhs.itsDateStructure),
	itsNbCoupons(rhs.itsNbCoupons),
	itsNbPathPricing(rhs.itsNbPathPricing),
	itsNbPathBounding(rhs.itsNbPathBounding),
	itsNbMaxBucket(rhs.itsNbMaxBucket),
	itsFixSteps(rhs.itsFixSteps),
	itsUSMethod(rhs.itsUSMethod),
	itsCalibMode(rhs.itsCalibMode),
	itsBetaCalib(rhs.itsBetaCalib),
	itsFixBoundary(rhs.itsFixBoundary),
	itsFixBeta(rhs.itsFixBeta),
	itsControlVariableFlag(rhs.itsControlVariableFlag),
	itsFixControlVariable(rhs.itsFixControlVariable),
	itsProductsToPrice(rhs.itsProductsToPrice),
	itsHasBeenPriced(rhs.itsHasBeenPriced),
	itsFirstCallIdx(rhs.itsFirstCallIdx),
	itsSnowBallSJPrice(rhs.itsSnowBallSJPrice),
	itsSnowBallSJStdDev(rhs.itsSnowBallSJStdDev),
	itsSnowBallSJRA(rhs.itsSnowBallSJRA),
	itsFundingPrice(rhs.itsFundingPrice),
	itsFundingStdDev(rhs.itsFundingStdDev),
	itsFundingRA(rhs.itsFundingRA),
	itsSnowBallCallablePrice(rhs.itsSnowBallCallablePrice),
	itsSnowBallCallableStdDev(rhs.itsSnowBallCallableStdDev),
	itsSnowBallCallableRA(rhs.itsSnowBallCallableRA),
	itsNoOptSBPrice(rhs.itsNoOptSBPrice),
	itsNoOptSBStdDev(rhs.itsNoOptSBStdDev),
	itsNoOptSBRA(rhs.itsNoOptSBRA),
	itsStackySBPrice(rhs.itsStackySBPrice),
	itsStackySBStdDev(rhs.itsStackySBStdDev),
	itsStackySBRA(rhs.itsStackySBRA),
	itsCouponOptPrice(rhs.itsCouponOptPrice),
	itsCouponOptStdDev(rhs.itsCouponOptStdDev),
	itsCouponOptRA(rhs.itsCouponOptRA),
	itsLBCouponAvgPrice(rhs.itsLBCouponAvgPrice),
	itsLBCouponAvgStdDev(rhs.itsLBCouponAvgStdDev),
	itsLBCouponAvgRA(rhs.itsLBCouponAvgRA),
	itsCouponOptAvgPrice(rhs.itsCouponOptAvgPrice),
	itsCouponOptAvgStdDev(rhs.itsCouponOptAvgStdDev),
	itsCouponOptAvgRA(rhs.itsCouponOptAvgRA),
	itsCVBetas(rhs.itsCVBetas),
	itsIsCpnFlow(rhs.itsIsCpnFlow),
	itsResetDatesVec(rhs.itsResetDatesVec),
	itsStartDatesVec(rhs.itsStartDatesVec),
	itsEndDatesVec(rhs.itsEndDatesVec),
	itsPeriodsVec(rhs.itsPeriodsVec),
	itsPaymentDatesVec(rhs.itsPaymentDatesVec),
	itsConstVec(rhs.itsConstVec),
	itsLevPrevVec(rhs.itsLevNewVec),
	itsLevNewVec(rhs.itsLevNewVec),
	itsStrikeOptVec(rhs.itsStrikeOptVec),
	itsCpnMinVec(rhs.itsCpnMinVec),
	itsCpnMaxVec(rhs.itsCpnMaxVec),
	itsStrikeSum(rhs.itsStrikeSum),
	itsFeesVec(rhs.itsFeesVec),
	itsvCpnNominal(rhs.itsvCpnNominal),
	itsvInitialFundNominal(rhs.itsvInitialFundNominal),
	itsvInitialFundSpread(rhs.itsvInitialFundSpread),
	itsGenType1(rhs.itsGenType1),
	itsGenType2(rhs.itsGenType2),
	itsFirstNbTimes(rhs.itsFirstNbTimes),
	itsFirstNbDims(rhs.itsFirstNbDims),
	itsPathScheme(rhs.itsPathScheme),
	itsPathOrder(rhs.itsPathOrder),
	itsAntithetic(rhs.itsAntithetic),
	itsModelType(rhs.itsModelType),
	itsTriggerMode(rhs.itsTriggerMode),
	itsHkData(rhs.itsHkData),
	itsCalibSwoptFreq(rhs.itsCalibSwoptFreq),
	itsRegressors(rhs.itsRegressors),
	itsFixBetaTimes(rhs.itsFixBetaTimes),
	itsFixBetaValues(rhs.itsFixBetaValues),
	itsFundDateStrip(rhs.itsFundDateStrip),
	itsStructDateStrip(rhs.itsStructDateStrip)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CRFCalculator
///	Routine: InitCRFForSummit
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::InitCSBFromSummit(ARM_ZeroCurve* zcCpn,
													   ARM_VolCurve* swoptVC,
													   ARM_VolCurve* capVC,
													   ARM_VolLInterpol* capRho,
													   ARM_VolLInterpol* capNu,
													   ARM_VolLInterpol* capBeta,
													   ARM_VolLInterpol* SwoptRho,
													   ARM_VolLInterpol* SwoptNu,
													   ARM_VolLInterpol* SwoptBeta,
													   double theHump,
													   double theBetaCorrel,
													   double theReCorrel,
													   int SABRSigmaOrAlpha)// SABR Sigma or Alpha
																			// By default=1(Sigma))
{
	vector <ARM_Object*> marketDatas;
	marketDatas.push_back(zcCpn);

	ARM_BSModel* capBSmod = NULL;
	ARM_BSModel* SwoptBSmod = NULL;
	int SABR_Flag;

	if (capRho && capNu)
	{
		if(capBeta)
		{		
			SABR_Flag = capBeta->IsEqualToOne()? K_SABR_ARITH:K_SABR_IMPLNVOL;
			capBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 0, zcCpn, zcCpn, capVC, K_YIELD, capRho,
											 capNu, SABR_Flag, capBeta, 0.5, SABRSigmaOrAlpha);
		}
		else
		{
			capBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 0, zcCpn, zcCpn, capVC, K_YIELD, capRho,
											 capNu, K_SABR_ARITH, NULL, 0.5, SABRSigmaOrAlpha);
		}
	}
	else
	{
		capBSmod = new ARM_BSModel(zcCpn->GetAsOfDate(), 0.0, zcCpn, zcCpn, capVC, K_YIELD);
	}

	if( SwoptRho && SwoptNu )
	{
		if(SwoptBeta)
		{
			SABR_Flag = SwoptBeta->IsEqualToOne()? K_SABR_ARITH : K_SABR_IMPLNVOL;
			SwoptBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 0.0,zcCpn,zcCpn,swoptVC,K_YIELD,
										       SwoptRho,SwoptNu,SABR_Flag, SwoptBeta, 0.5, SABRSigmaOrAlpha);
		}
		else
		{
			SwoptBSmod = new ARM_BSSmiledModel(zcCpn->GetAsOfDate(), 0.0,zcCpn,zcCpn,swoptVC,K_YIELD,
										       SwoptRho,SwoptNu,K_SABR_ARITH, NULL, 0.5, SABRSigmaOrAlpha);
		}
	}
	else
	{
		SwoptBSmod = new ARM_BSModel(zcCpn->GetAsOfDate(), 0.0, zcCpn, zcCpn, swoptVC, K_YIELD);
	}


	marketDatas.push_back(SwoptBSmod);
	marketDatas.push_back(capBSmod);

	std::vector<double> times(1,0.0);

	ARM_CurveModelParam tMr;
	std::vector<double> values(1,MRS_DEFAULT_VALUE);
	tMr = ARM_CurveModelParam(ARM_ModelParamType::MeanReversion,&values,&times,"MRS");
	marketDatas.push_back(&tMr);

	ARM_CurveModelParam tBeta;
	values = std::vector<double>(1,BETA_DEFAULT_VALUE);
    tBeta = ARM_CurveModelParam(ARM_ModelParamType::Beta,&values,&times,"BETA");
	marketDatas.push_back(&tBeta);

	ARM_CurveModelParam tCorrel;
	values = std::vector<double>(1,CORREL_DEFAULT_VALUE);
    tCorrel = ARM_CurveModelParam(ARM_ModelParamType::Correlation,&values,&times,"CORRELATION");
	marketDatas.push_back(&tCorrel);

	if ( theHump != -10000.0 )
	{
		values= std::vector<double>(1, theHump);
    }
    else
    {
		values= std::vector<double>(1, HUMP_DEFAULT_VALUE);
	}
	ARM_CurveModelParam hump = ARM_CurveModelParam(ARM_ModelParamType::Hump,
												  &values,
												  &times,
												  "HUMP");
	marketDatas.push_back(&hump);

	if ( theBetaCorrel != -10000.0 )
	{
		values= std::vector<double>(1, theBetaCorrel);
    }
    else
    {
		values= std::vector<double>(1, BETA_CORREL_DEFAULT_VALUE);
	}
	ARM_CurveModelParam betacorrel = ARM_CurveModelParam(ARM_ModelParamType::BetaCorrelation,
												  &values,
												  &times,
												  "BETACORREL");
	marketDatas.push_back(&betacorrel);

	if ( theReCorrel != -10000.0 )
	{
		values= std::vector<double>(1, theReCorrel);
    }
    else
    {
		values= std::vector<double>(1, RECORREL_DEFAULT_VALUE);
	}
	ARM_CurveModelParam recorrel = ARM_CurveModelParam(ARM_ModelParamType::ReCorrelation,
												  &values,
												  &times,
												  "RECORRELATION");
	marketDatas.push_back(&recorrel);

	Init(marketDatas); // Creation of Market data manager

	if (capBSmod)
	   delete capBSmod;
	capBSmod = NULL;

	if( SwoptBSmod )
		delete SwoptBSmod;
	SwoptBSmod = NULL;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: InitDateStrip
///	Returns: itsDateStructure
///	Action : precompute on the DateStrip
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CallableSnowBallCalculator::InitDateStrip()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	// Build the exercise schedule
	ARM_DateStrip ExerDateStrip(itsStartDate,itsEndDate, itsCpnFreq, itsCpnDayCount, itsCpnResetCal.c_str(), K_MOD_FOLLOWING, itsCpnIntRule,
					K_SHORTSTART, -abs(itsCpnResetLag),itsCpnFreq,GETDEFAULTVALUE,itsCpnPayCal.c_str(),itsCpnIndexTiming, K_ARREARS);

	/// past  no call
	size_t nbPastNoCall=0;
	std::vector<double>& dates = ExerDateStrip.GetResetDates() ;
	size_t size = ExerDateStrip.size();
	while(nbPastNoCall < size && (*dates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
		++nbPastNoCall;	
	ExerDateStrip.ResizeAndBuilt(nbPastNoCall,size);	
	double firstCancelDate = (*ExerDateStrip.GetResetDates())[0];

	// builds the Cpn schedule
	ARM_DateStrip CpnDateStrip(itsStartDate,itsEndDate, itsCpnFreq, itsCpnDayCount, itsCpnResetCal.c_str(), K_MOD_FOLLOWING, 
				itsCpnIntRule,K_SHORTSTART, -abs(itsCpnResetLag),itsCpnFreq,GETDEFAULTVALUE,itsCpnPayCal.c_str(),itsCpnIndexTiming, K_ARREARS);

	/// past  no call
	nbPastNoCall=0;
	dates = CpnDateStrip.GetResetDates() ;
	size = CpnDateStrip.size();
	while(nbPastNoCall < size && (*dates)[nbPastNoCall] < firstCancelDate - K_NEW_DOUBLE_TOL)
		++nbPastNoCall;	
	CpnDateStrip.ResizeAndBuilt(nbPastNoCall,size);

	// builds the Cpn schedule
	ARM_DateStrip FundingSched(itsStartDate,itsEndDate, itsFundFreq, itsFundDayCount, itsCpnResetCal.c_str(), K_MOD_FOLLOWING, 
					itsCpnIntRule, K_SHORTSTART, GETDEFAULTVALUE,itsFundFreq, GETDEFAULTVALUE, itsCpnPayCal.c_str(), itsCpnIndexTiming, K_ARREARS);

	/// past  no call
	nbPastNoCall=0;
	dates = FundingSched.GetResetDates() ;
	size = FundingSched.size();
	while(nbPastNoCall < size && (*dates)[nbPastNoCall] < firstCancelDate - K_NEW_DOUBLE_TOL)
		++nbPastNoCall;
	
	FundingSched.ResizeAndBuilt(nbPastNoCall,size);

	itsStructDateStrip = ARM_DateStripPtr(new ARM_DateStrip(CpnDateStrip));	
	itsFundDateStrip   = ARM_DateStripPtr(new ARM_DateStrip(FundingSched));	

	itsNbCoupons = CpnDateStrip.GetFlowStartDates()->size();
	
    ARM_DateStripVector SchedVect(NB_CALLABLESB_SCHED,NULL);
    SchedVect[CPN_SCHED] = &CpnDateStrip;
    SchedVect[EXER_SCHED] = &ExerDateStrip;

	ARM_DateStripCombiner EventSchedule (SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	CC_MUTABLE( ARM_CallableSnowBallCalculator, itsFirstCallIdx ) = 0;
	while(itsFirstCallIdx < eventDates->size() && (*eventDates)[itsFirstCallIdx] <= asOfDate)
        ++CC_MUTABLE( ARM_CallableSnowBallCalculator, itsFirstCallIdx );

	ARM_VectorPtr resetDates = EventSchedule.GetMergeData();
	std::vector<double>& cpnResetDates = EventSchedule.GetDateStrip(CPN_SCHED)->GetResetDates();
	std::vector<double>& exerDates = EventSchedule.GetDateStrip(EXER_SCHED)->GetResetDates();

	itsSBRank2Pos.resize(cpnResetDates->size(),ARM_DateStripCombiner::DateStripCombiner_BlankData);
	itsSBPos2Rank.resize(resetDates->size(),ARM_DateStripCombiner::DateStripCombiner_BlankData);
	itsExerRank2Pos.resize(cpnResetDates->size(),ARM_DateStripCombiner::DateStripCombiner_BlankData);
	itsExerPos2Rank.resize(resetDates->size(),ARM_DateStripCombiner::DateStripCombiner_BlankData);

	int curSBRank = 0;
	int curExerRank = 0;

	//FIXMEFRED: no vector<bool> anymore
	itsIsCpnFlow = std::deque<bool>(cpnResetDates->size(),false);

	int i;
	for (i = 0; i < resetDates->size(); i++)
	{
		if ((*cpnResetDates)[i] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
		{
			itsSBRank2Pos[curSBRank] = i;
			itsSBPos2Rank[i] = curSBRank;
			curSBRank++;

			itsIsCpnFlow[i] = 1;
		}
		if ((*exerDates)[i] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
		{
			itsExerRank2Pos[curExerRank] = i;
			itsExerPos2Rank[i] = curExerRank;
			curExerRank++;
		}
	}

	itsResetDatesVec = *(CpnDateStrip.GetResetDates());
	itsStartDatesVec = *(CpnDateStrip.GetFwdRateStartDates());
	itsEndDatesVec = *(CpnDateStrip.GetFwdRateEndDates());
	itsPaymentDatesVec = *(CpnDateStrip.GetPaymentDates());

	itsPeriodsVec.resize(itsResetDatesVec.size());
	itsConstVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsLevPrevVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsLevNewVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsStrikeOptVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsCpnMinVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsCpnMaxVec = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsvCpnNominal= std::vector<double>(itsResetDatesVec.size(),0.0);
	itsStrikeSum = std::vector<double>(itsResetDatesVec.size(),0.0);
	itsFeesVec = std::vector<double>(exerDates->size(),0.0);

	std::vector<double>& flowStartDates = CpnDateStrip.GetFlowStartDates();

	for (i = 0; i < itsStartDatesVec.size(); ++i)
	{
		itsPeriodsVec[i] = CountYears( itsCpnIndexDayCount,itsStartDatesVec[i],itsEndDatesVec[i]);
		itsConstVec[i] = itsConstProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsLevPrevVec[i] = itsLPrevCpnProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsLevNewVec[i] = itsLNewOptProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsStrikeOptVec[i] = itsStrikeOptProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsCpnMaxVec[i] = itsMaxCpnProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsCpnMinVec[i] = itsMinCpnProfile.Interpolate((*flowStartDates)[i]-asOfDate);
		itsvCpnNominal[i] = itsNotionalProfile.Interpolate((*flowStartDates)[i]-asOfDate);
	}

	std::vector<double>& exerciseDates = ExerDateStrip.GetResetDates();
	for (i = 0; i < exerciseDates->size(); ++i)
	{
		itsFeesVec[i] = itsFeesProfile.Interpolate((*exerciseDates)[i]-asOfDate);
	}

	std::vector<double>& fundResetDates = FundingSched.GetResetDates() ;
	double lag, margin;

	if(IsBasis()){
		for (i=0; i<fundResetDates->size(); i++)
		{
			lag = (*fundResetDates)[i]-asOfDate;
			margin = itsInitFundMarginProfile.Interpolate(lag);
			itsvInitialFundSpread.push_back(margin );
			itsvInitialFundNominal.push_back(itsFundNotionalProfile.Interpolate(lag));
		}
		itsFundMarginProfile.SetOrdinates(ComputeDomesticBasis());
	}

	// Compute the sum option coeffs
	ComputeSumOptCoeffs();

    return EventSchedule;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if Caption datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CheckData()
{
	if (((itsModelType == ARM_PricingModelType::HK) || (itsModelType == ARM_PricingModelType::SBGM)) 
		&& ((itsCpnIndexTiming == K_ARREARS) || (itsNoticeDays != itsCpnResetLag))
		&& (itsTriggerMode == TCoupon)
		&& (itsProductsToPrice[SnowBallCallablePrice]))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT,"You cannot price a CSB ARREARS with the Trigger Type 'Coupon' and the models HK & SBGM");
	}

	if (((itsModelType == ARM_PricingModelType::HK) || (itsModelType == ARM_PricingModelType::SBGM)) 
		&& (itsNoticeDays != itsCpnResetLag) 
		&& (itsProductsToPrice[SnowBallCallablePrice]))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT,"You cannot price a CSB with reset gap <> notice gap and the models HK & SBGM");
	}

	if (itsCalibSwoptFreq > itsCpnFreq)
	{
		ARM_THROW( ERR_INVALID_ARGUMENT,"The frequency of the calbration swaption cannot be greater than the cpn index frequency.");
	}

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the Caption deal.
/////////////////////////////////////////////////////////////////
double ARM_CallableSnowBallCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

	ARM_StringVector columnNames = CreateColumnNames();

	size_t i;

	ARM_StringVector anaColumnNames;
	anaColumnNames.reserve(NbAnaProductsToPrice);

	for (i = 0; i < NbAnaProductsToPrice; ++i)
	{
		if (itsProductsToPrice[NbProductsToPrice+i])
			anaColumnNames.push_back( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::AnaProductToPriceColumns[i] ] );
	}

	if (anaColumnNames.size() > 0)
	{
		// Use the CF Num Method
		ARM_PricingModelPtr anaModel = GetPricingModel();
		ARM_NumMethodPtr mcNumMethod = anaModel->GetNumMethod();
		ARM_NumMethodPtr cfNumMethod = ARM_NumMethodPtr(new ARM_CFMethod);
		anaModel->SetNumMethod(cfNumMethod);

		ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
		ARM_DealDescriptionPtr anaDealDesc = ARM_DealDescriptionPtr( new ARM_DealDescription( dealDesc.GetText(), dealDesc.GetFormat(), dealDesc.GetRowsNb(), dealDesc.GetColsNb(), anaColumnNames ) );
		ARM_GenSecurityPtr anaGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(anaDealDesc,"",CreateCstManager(),true,false));
		ARM_GenPricer* anaGenPricer = new ARM_GenPricer( &*anaGenSec, &*anaModel);
		anaGenPricer->Price();
		ARM_AutoCleaner<ARM_GenPricer> HoldAnaGP( anaGenPricer );

		double price;

		for ( i = 0; i < NbAnaProductsToPrice; ++i)
		{
			if (itsProductsToPrice[NbProductsToPrice+i])
			{
				price	= anaGenPricer->GetPricerInfo()->GetContents( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::AnaProductToPriceColumns[i] ] ).GetData("Price").GetDouble();
				
				if (i == AnaLBCouponAvgPrice)
				{
					itsAnaLBCouponAvgPrice = price;
				}
				else if (i == AnaLBCouponPrice)
				{
					itsAnaLBCouponPrice = price;
				}
			}
		}

		anaModel->SetNumMethod(mcNumMethod);
	}

	if (columnNames.size())
	{
		ARM_GenPricer* genPricer;

		if( itsControlVariableFlag != NO)
		{
			ARM_StringVector cvColumnNames(1);
			cvColumnNames[0] = CallableSnowBallColNamesTable[ LBCouponPaid ];
			std::vector<double> columnPrices(1);
			columnPrices[0] = itsAnaLBCouponPrice;

			string refPriceColumn;

			if (itsControlVariableFlag == SB)
				refPriceColumn = CallableSnowBallColNamesTable[ CF ];
			else if (itsControlVariableFlag == CSB)
				refPriceColumn = CallableSnowBallColNamesTable[ Bermuda ];

			if (!itsFixControlVariable)
				itsCVBetas->resize(0);

			genPricer = new ARM_GenPricer( &* GetGenSecurity(),&*GetPricingModel(), cvColumnNames, columnPrices, refPriceColumn, *itsCVBetas);
		}
		else
		{
			genPricer = new ARM_GenPricer( &*GetGenSecurity(),&*GetPricingModel());
		}

		genPricer->Price();

		if (itsFixBoundary)
		{
			ARM_GenSecManipulator genSecManipulator;

			ARM_DealDescriptionPtr dealDes= GetGenSecurity()->GetDealDescriptionPtr();

			genSecManipulator.ChangeAmericanIntoTrigger( *GetGenSecurity(), dealDes);
		
			itsFixBoundary = false;
		}

		if (itsControlVariableFlag != NO)
			itsCVBetas = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("Beta").GetVector();

		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

		double price, stdDev;
		ARM_GP_VectorPtr ra;

		for ( i = 0; i < NbProductsToPrice; ++i)
		{
			if (itsProductsToPrice[i])
			{
				price	= genPricer->GetPricerInfo()->GetContents( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::ProductToPriceColumns()[i] ] ).GetData("Price").GetDouble();
				stdDev	= genPricer->GetPricerInfo()->GetContents( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::ProductToPriceColumns()[i] ] ).GetData("StdDev").GetDouble();
				ra		= genPricer->GetPricerInfo()->GetContents( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::ProductToPriceColumns()[i] ] ).GetData("RunningAvg").GetVector();
				
				if (i == SnowBallSJPrice)
				{
					if( itsControlVariableFlag == SB)
					{
						itsSnowBallSJPrice = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("Price").GetDouble();
						itsSnowBallSJStdDev = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("StdDev").GetDouble();
						itsSnowBallSJRA = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("RunningAvg").GetVector();
					}
					else
					{
						itsSnowBallSJPrice = price;
						itsSnowBallSJStdDev = stdDev;
						itsSnowBallSJRA = ra;
					}
				}
				else if (i == FundingPrice)
				{
					itsFundingPrice = price;
					itsFundingStdDev = stdDev;
					itsFundingRA = ra;
				}
				else if (i == SnowBallCallablePrice)
				{
					if( itsControlVariableFlag == CSB)
					{
						itsSnowBallCallablePrice = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("Price").GetDouble();
						itsSnowBallCallableStdDev = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("StdDev").GetDouble();
						itsSnowBallCallableRA = genPricer->GetPricerInfo()->GetContents( ARM_GenPricer::CVTag ).GetData("RunningAvg").GetVector();
					}
					else
					{
						itsSnowBallCallablePrice = price;
						itsSnowBallCallableStdDev = stdDev;
						itsSnowBallCallableRA = ra;
					}
				}
				else if (i == NoOptSBPrice)
				{
					itsNoOptSBPrice = price;
					itsNoOptSBStdDev = stdDev;
					itsNoOptSBRA = ra;
				}
				else if (i == StackySBPrice)
				{
					itsStackySBPrice = price;
					itsStackySBStdDev = stdDev;
					itsStackySBRA = ra;
				}
				else if (i == CouponOptPrice)
				{
					itsCouponOptPrice = price;
					itsCouponOptStdDev = stdDev;
					itsCouponOptRA = ra;
				}
				else if (i == LBCouponAvgPrice)
				{
					itsLBCouponAvgPrice = price;
					itsLBCouponAvgStdDev = stdDev;
					itsLBCouponAvgRA = ra;
				}
				else if ( i == CouponOptAvgPrice)
				{
					itsCouponOptAvgPrice = price;
					itsCouponOptAvgStdDev = stdDev;
					itsCouponOptAvgRA = ra;
				}
			}
		}
	}

	itsHasBeenPriced = true;
    
    return itsSnowBallCallablePrice;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_CallableSnowBallCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CallableSnowBallCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[SnowBallCallablePrice])
	{
		GetPricingData()[ "CSBPrice"		] = itsSnowBallCallablePrice;
		GetPricingData()[ "CSBStdDev"		] = itsSnowBallCallableStdDev;
		GetPricingData()[ "CSBRunningAvg"	] = itsSnowBallCallableRA;
	}
	if (itsProductsToPrice[SnowBallSJPrice])
	{
		GetPricingData()[ "SBPrice"		] = itsSnowBallSJPrice;
		GetPricingData()[ "SBStdDev"	] = itsSnowBallSJStdDev;
		GetPricingData()[ "SBRunningAvg"] = itsSnowBallSJRA;
	}
	if (itsProductsToPrice[FundingPrice])
	{
		GetPricingData()[ "FundingPrice"		] = itsFundingPrice;
		GetPricingData()[ "FundingStdDev"		] = itsFundingStdDev;
		GetPricingData()[ "FundingRunningAvg"	] = itsFundingRA;
	}
	if (itsProductsToPrice[NoOptSBPrice])
	{
		GetPricingData()[ "NoOptSBPrice"		] = itsNoOptSBPrice;
		GetPricingData()[ "NoOptSBStdDev"		] = itsNoOptSBStdDev;
		GetPricingData()[ "NoOptSBRunningAvg"	] = itsNoOptSBRA;
	}
	if (itsProductsToPrice[StackySBPrice])
	{
		GetPricingData()[ "StackySBPrice"		] = itsStackySBPrice;
		GetPricingData()[ "StackySBStdDev"		] = itsStackySBStdDev;
		GetPricingData()[ "StackySBRunningAvg"	] = itsStackySBRA;
	}
	if(itsProductsToPrice[CouponOptPrice])
	{
		GetPricingData()[ "CouponOptPrice"		] = itsCouponOptPrice;
		GetPricingData()[ "CouponOptStdDev"		] = itsCouponOptStdDev;
		GetPricingData()[ "CouponOptRunningAvg"	] = itsCouponOptRA;
	}
	if(itsProductsToPrice[LBCouponAvgPrice])
	{
		GetPricingData()[ "LBCouponAvgPrice"		] = itsLBCouponAvgPrice;
		GetPricingData()[ "LBCouponAvgStdDev"		] = itsLBCouponAvgStdDev;
		GetPricingData()[ "LBCouponAvgRunningAvg"	] = itsLBCouponAvgRA;
	}
	if(itsProductsToPrice[CouponOptAvgPrice])
	{
		GetPricingData()[ "CouponOptAvgPrice"		] = itsCouponOptAvgPrice;
		GetPricingData()[ "CouponOptAvgStdDev"		] = itsCouponOptAvgStdDev;
		GetPricingData()[ "CouponOptAvgRunningAvg"	] = itsCouponOptAvgRA;
	}
	if(itsProductsToPrice[NbProductsToPrice+AnaLBCouponAvgPrice])
	{
		GetPricingData()[ "AnaLBCouponAvgPrice"		] = itsAnaLBCouponAvgPrice;
	}
	if(itsProductsToPrice[NbProductsToPrice+AnaLBCouponPrice])
	{
		GetPricingData()[ "AnaLBCouponPrice"		] = itsAnaLBCouponPrice;
	}
	if (itsControlVariableFlag != NO)
	{
		GetPricingData()[ "CVBetas"					] = itsCVBetas;
	}

}




////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::Calibrate_SFRM2F()
{
	//1.- functional calibration
	ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());

	GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::Calibrate_SBGM()
{
	//1.- functional calibration
	ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());

	GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());

	//2.- correl calibration
	if(itsCalibMode == CAPSWAPTION)
	{
		ARM_StdPortfolioPtr pfPtr(CreateSwoptPortfolio());

		ARM_CalibMethod* correlCalib = new ARM_CalibMethod(pfPtr,ARM_ModelParamVector(),ARM_CalibMethodType::Bootstrap1D,
															ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);

		SetCalibMethod( ARM_CalibMethodPtr( correlCalib ) );

		ComputePortPrices(false, SWAPTION, &*GetCalibMethod());

		ARM_ModelParam* initBetaCorrel = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[BetaCorrelKey]));
		
		ARM_CurveModelParam* beta = NULL;
			
		if (!itsFixBeta || !itsFixBetaTimes.size() || !itsFixBetaValues.size())
		{
			if (!initBetaCorrel)
			{
				/// Replace the beta param with new initialisation
				std::vector<double> initTimes(1,0.0);
				std::vector<double> initBetas(1,BETA_CORREL_DEFAULT_VALUE);
				std::vector<double> betaLowerBound(1,BETA_CORREL_LOWER_BOUND);
				std::vector<double> betaUpperBound(1,BETA_CORREL_UPPER_BOUND);
				beta = new ARM_CurveModelParam(
						ARM_ModelParamType::BetaCorrelation,
						&initBetas,
						&initTimes,
						"BETA",
						"STEPUPRIGHT",
						&betaLowerBound,
						&betaUpperBound);
			}
			else
			{
				beta = static_cast<ARM_CurveModelParam*>(initBetaCorrel->Clone());
			}
		}
		else
		{
			beta = new ARM_CurveModelParam(
					ARM_ModelParamType::BetaCorrelation,
					&itsFixBetaValues,
					&itsFixBetaTimes,
					"BETA",
					"STEPUPRIGHT",
					&itsFixBetaValues,
					&itsFixBetaValues);
		}
			
		GetCalibMethod()->GetCalibParams().push_back(beta);
		
		GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());

		if (itsFixBeta)
			FixBetaModelParam();
	
	}
}



////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::Calibrate_HK()
{
	ARM_NumMethodPtr mcMethod = GetPricingModel()->GetNumMethod();
	
	int numSchemeType = ARM_ArgConv_PDENumSchemeType.GetNumber("CN1F");
	ARM_PDENumericalScheme* numScheme = ARM_PDENumericalScheme::getNumericalSchemeInstanceById(numSchemeType);

	int nbSteps   = (itsHkData.size()>0)?itsHkData[0]:PDE_NB_STEPS;
	int gridSize  = (itsHkData.size()>1)?itsHkData[1]:PDE_GRIDSIZE;
	int nbStdDevs = (itsHkData.size()>2)?itsHkData[2]:PDE_NB_STDEV;

	ARM_PDEMethod* pPDEMethod = new ARM_PDEMethod(numScheme,nbSteps, gridSize, 2*nbStdDevs);
	GetPricingModel()->SetNumMethod( ARM_NumMethodPtr( pPDEMethod ) );

	//1.- functional calibration
	ARM_MultiAssetsModel* model= dynamic_cast< ARM_MultiAssetsModel* >(&*GetPricingModel());

	GetCalibMethod()->Calibrate(&*(*model->GetModelMap())[GetKeys()[YcKey]]->Model());
	
	//GetCalibMethod()->Calibrate(&*GetPricingModel());
	    
	GetPricingModel()->SetNumMethod(mcMethod);
}


////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::Calibrate()
{
	itsHasBeenPriced = false;
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
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CSBCalculator::Calibrate : wrong model type !" );
	}
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateCalibration_SFRM2F(bool isUpdateStrike)
{
	ARM_CalibMethod* volCalibMethod	= GetVolCalibMethod();
	ComputePortPrices(true,itsCalibMode,volCalibMethod);
	 
	if((itsCalibMode == SUMOPT) && itsBetaCalib)
	{
		ARM_CalibMethod* betaCalibMethod	= GetBetaCalibMethod();
		ARM_StdPortfolioPtr betaPF		= betaCalibMethod->GetPortfolio();
		UpdateBetaCapNominal(betaCalibMethod);
		ComputePortPrices(false,CAP,betaCalibMethod);
	}
	itsHasBeenPriced=false;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateCalibration_SBGM(bool isUpdateStrike)
{
	CreateAndSetCalibration();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateCalibration_HK(bool isUpdateStrike)
{
	/// todo
	ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::UpdateCalibration_HK : not implemented !" );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateCalibration(bool isUpdateStrike)
{
	itsHasBeenPriced = false;
		/// re-update the funding spread if necessairy
	if(IsBasis()){	
		itsFundMarginProfile.SetOrdinates(ComputeDomesticBasis());
	/// Fake FIX FIX FIX
	    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[YcBasisDomKey],CreateColumnNames(),CreateCstManager(),itsFixBoundary,(itsControlVariableFlag!=NO));
	}

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
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::UpdateCalibration : wrong model type !" );
	}
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateModel_SFRM2F()
{
	ARM_SFRM* refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());

	ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	refModel->SetZeroCurve(CreateClonedPtr( zcCurve ));

	/// create model params
	ARM_ModelParamsSFRM* pSFRMModelParams = CreateSFRMModelParams( itsCalibMode );

	refModel->SetModelParams(*pSFRMModelParams);

	delete pSFRMModelParams;

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateModel_SBGM()
{
	/// un peu bourrin
	CreateAndSetModel();
	CreateAndSetCalibration();
	itsHasBeenPriced=false;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateModel_HK()
{
	/// un peu bourrin
	CreateAndSetModel();
	CreateAndSetCalibration();
	itsHasBeenPriced=false;

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateModel()
{
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
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CSBCalculator::UpdateModel : wrong model type !" );
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CallableSnowBallCalculator::CreateEmptyCalibration( CalibMode calibMode, bool betaCalib)
{
	/// Create cap floor portfolio (sigma calibration)
    ARM_StdPortfolioPtr PF;

	if (calibMode == CAP )
		PF = CreateCapPortfolio(true);
	else if (calibMode == SWAPTION )
		PF = CreateSwoptPortfolio();
	else if (calibMode == SUMOPT)
		PF = CreateSumOptionPortfolio();

	ARM_CalibMethod* calibMethod = new ARM_CalibMethod(
				PF,
				ARM_ModelParamVector(),
				ARM_CalibMethodType::Bootstrap1D,
				ARM_MAX_ITER,
				ARM_CalibrationTarget::PriceTarget,
				NULL,
				NULL);

	if (betaCalib)
	{
		if (calibMode == SUMOPT )
		{
			ARM_StdPortfolioPtr capPF = CreateBetaCapPortfolio();

			double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();

			std::vector<double> initTimes(1,itsResetDatesVec[0]-asOf);
			std::vector<double> initBeta(1,BETA_DEFAULT_VALUE);
			std::vector<double> BetaLB(1,BETA_LOWER_BOUND);
			std::vector<double> BetaUB(1,BETA_UPPER_BOUND);

			ARM_CurveModelParam* beta = new ARM_CurveModelParam(
				ARM_ModelParamType::Beta,
				&initBeta,
				&initTimes,
				"Beta",
				"STEPUPRIGHT",
				&BetaLB,
				&BetaUB);

			ARM_ModelFitterDes betaModelFitter(ARM_ModelFitterSolverType::NewtonRaphsonWithDichotomy);

			ARM_ModelParamVector betaVec(1,beta);

			ARM_CalibMethod* betaCalibMethod =
				new ARM_CalibMethod(
					capPF,
					betaVec,
					ARM_CalibMethodType::Bootstrap1D,
					&betaModelFitter,
					ARM_CalibrationTarget::PriceTarget,
					NULL,
					NULL,
					false,
					0,
					1
					);

			betaCalibMethod->SetlinkedMethod(calibMethod);

			calibMethod = betaCalibMethod;
		}
		else
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Beta Calibration is just possible with volatility calibrated on CAP.");
		}
	}
	
	return calibMethod;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetCalibration_SFRM2F		
///	Returns: void
///	Action : create the calibration for SFRM
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetCalibration_SFRM2F()
{
	CheckCalibMode();

	ARM_CalibMethod* calibMethod = CreateEmptyCalibration(itsCalibMode,itsBetaCalib);
	SetCalibMethod(ARM_CalibMethodPtr(calibMethod));

	ARM_CalibMethod* volCalibMethod	= GetVolCalibMethod();
	ComputePortPrices(true,itsCalibMode,volCalibMethod);
	 
	if((itsCalibMode == SUMOPT) && itsBetaCalib)
	{
		/// 1) calibrate the model ATM
		ARM_CalibMethod* betaCalibMethod	= GetBetaCalibMethod();
		ARM_StdPortfolioPtr betaPF		= betaCalibMethod->GetPortfolio();
		ComputePortPrices(false,CAP,betaCalibMethod);
	}

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetCalibration_SBGM	
///	Returns: void
///	Action : create the calibration for SBGM
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetCalibration_SBGM()
{
	if (itsCalibMode == CAP || itsCalibMode == CAPSWAPTION)
		CreateAndSetDensityCalibration(ARM_StdPortfolioPtr(NULL));
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" ARM_CallableSnowBallCalculator::CreateAndSetCalibration: SBGM supports only CAP and CAPSWOPT calib mode!" );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetCalibration_HK
///	Returns: void
///	Action : create the calibration for HK
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetCalibration_HK()
{	
	if (itsCalibMode == CAPSWAPTION )
	{
		ARM_StdPortfolioPtr pf = CreateSwoptPortfolio();
		
		CreateAndSetDensityCalibration( pf );

		/// à faire là, ou bien plus tard ?  
		ComputePortPrices(false, SWAPTION, &*GetCalibMethod());
	}
	else if (itsCalibMode == CAP)
		CreateAndSetDensityCalibration( ARM_StdPortfolioPtr(NULL) );
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +" ARM_CallableSnowBallCalculator::CreateAndSetCalibration: HK supports only CAP and CAPSWOPT calib mode!" );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetDensityCalibration
///	Returns: void
///	Action : create the calibration for HK / SBGM
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetDensityCalibration(ARM_StdPortfolioPtr pfPtr)
{
	ARM_ZeroCurve* pCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	ARM_BSSmiledModel* cfBSSmiledModel = dynamic_cast< ARM_BSSmiledModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	if (!cfBSSmiledModel)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" ARM_CallableSnowBallCalculator::CreateAndSetNumCalibration: no BSSMiledModel available !" );
		
	ARM_VolCurve* pSigma	= cfBSSmiledModel->GetVolatility();

	ARM_VolCurve* pAlpha	= cfBSSmiledModel->FromSigmaToAlpha(pSigma);
	ARM_VolCurve* pBeta		= cfBSSmiledModel->GetBeta();
	ARM_VolCurve* pRho		= cfBSSmiledModel->GetRho();
	ARM_VolCurve* pNu		= cfBSSmiledModel->GetNu();

	ARM_DateStrip* pDS		= GetCalibSchedule();

	std::vector<double>& pRD		= pDS->GetResetDates();	//not cloned, do not delete
	
	ARM_VanillaSecDensityPtrVector vanillaSecDensities;
	std::vector<double>::iterator resetD	= pRD->begin();
	std::vector<double>::iterator startD	= pDS->GetFlowStartDates()->begin();
	std::vector<double>::iterator endD	= pDS->GetFlowEndDates()->begin();

	for(; resetD != pRD->end(); ++resetD,++startD,++endD)
	{
		ARM_VanillaSecurityDensity* pVSD = GetMarketDensity( *resetD, *startD, *endD,
															 pAlpha, pBeta, pRho, pNu, 
															 cfBSSmiledModel->GetSABRFlag(), pCurve);

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

	SetCalibMethod(ARM_CalibMethodPtr(pMethod));	
}

///////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetMarketDensity
///	Returns: ARM_VanillaSecurityDensity*
///	Action : 
////////////////////////////////////////////////////
ARM_VanillaSecurityDensity* 
ARM_CallableSnowBallCalculator::GetMarketDensity(	double julianRD, double julianSD, double julianED,
													ARM_VolCurve* pAlpha, ARM_VolCurve* pBeta, ARM_VolCurve* pRho, ARM_VolCurve* pNu,
													int sabrFlag, ARM_ZeroCurve* pCurve) const
{
	const ARM_Date asOfDate				= GetMktDataManager()->GetAsOfDate();
	double julianAsof					= asOfDate.GetJulian();

	ARM_VanillaSecDensityPtrVector vanillaSecDensities;
	
	double mat		= (julianRD - julianAsof) / K_YEAR_LEN;
	double delta	= (julianED - julianSD  ) / K_YEAR_LEN;
	int freq		=  CC_Round( 1./delta ) ;

	ARM_DensityFunctor* pSDF = NULL;

	if (mat < K_DOUBLE_TOL)
	{
		pSDF = new ARM_ShiftedLNDensityFunctor(0.2,0.0);
	}
	else
	{
		pSDF = new ARM_SABRDensityFunctor(	pAlpha->ComputeVolatility(mat,delta)/100.0, 
																	pBeta->ComputeVolatility(mat,delta), 
																	pRho->ComputeVolatility(mat,delta), 
																	pNu->ComputeVolatility(mat,delta), 
																	sabrFlag );
	}

	ARM_VanillaSecurityDensity* pVSD = new ARM_VanillaSecurityDensity(	julianRD,
																		julianSD,
																		julianED,
																		ARM_DensityFunctorPtr(pSDF), 
																		freq,
																		itsCpnIndexDayCount,
																		GETDEFAULTVALUE,
																		1.,0.,1.,
																		CreateClonedPtr(pCurve));
	return pVSD;
}

///////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetCalibSchedule
///	Returns: ARM_DateStrip*
///	Action : 
////////////////////////////////////////////////////
ARM_DateStrip* ARM_CallableSnowBallCalculator::GetCalibSchedule() const
{	
	if (itsCpnResetLag != itsNoticeDays)
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::GetCalibSchedule : Nb of notice days is required to match reset lag." );

	bool isArrearsCpn	= (itsCpnIndexTiming == K_ARREARS);
	int cpnIndexFreq	= ARM_ArgConv_MatFrequency.GetNumber(itsCpnIdxTerm);
	int maxFreq			= (cpnIndexFreq<itsCpnFreq) ? itsCpnFreq : cpnIndexFreq;
	
	ARM_Date Date1( itsStartDate );
	ARM_Date Date2( itsEndDate);
	
	size_t nbmonth ;
	if		(cpnIndexFreq == K_ANNUAL)		nbmonth	= 12;
	else if (cpnIndexFreq == K_SEMIANNUAL)	nbmonth	= 6;
	else if (cpnIndexFreq == K_QUARTERLY)	nbmonth	= 3;
	else if (cpnIndexFreq == K_MONTHLY)		nbmonth	= 1;
	else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::GetCalibSchedule : index type not supported" );

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
			else ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::GetCalibSchedule : cpn freq not supported" );
			
			Date2.AddMonths(-aux);
			Date2.AddMonths(nbmonth);
		}

	pDS = new ARM_DateStrip(Date1, Date2, maxFreq,
							itsCpnDayCount, itsCpnResetCal.c_str(),	K_MOD_FOLLOWING,
							itsCpnIntRule, K_SHORTSTART, -abs(itsCpnResetLag),
							maxFreq, GETDEFAULTVALUE, itsCpnPayCal.c_str(),
							K_ADVANCE, K_ARREARS);

	/// 
	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	std::vector<double>& resetDates = pDS->GetResetDates();

	
	/// if start date is in the past, regenerate schedule
	int idx = 0;
	
	while ((*resetDates)[idx]<asOf)
		idx++;

	if (idx)
	{
		Date1 = pDS->GetFlowStartDates()->Elt(idx);
		
		delete pDS;

		pDS = new ARM_DateStrip(Date1, Date2, maxFreq,
							itsCpnDayCount, itsCpnResetCal.c_str(),	K_MOD_FOLLOWING,
							itsCpnIntRule, K_SHORTSTART, itsCpnResetLag,
							maxFreq, GETDEFAULTVALUE, itsCpnPayCal.c_str(),
							K_ADVANCE, K_ARREARS);
	}
	
	size_t nper = maxFreq / cpnIndexFreq;

	std::vector<double>& fwdEndDates    = pDS->GetFwdRateEndDates();
	std::vector<double>& fwdStartDates  = pDS->GetFwdRateStartDates();
	std::vector<double>& flowStartDates = pDS->GetFlowStartDates();
	std::vector<double>& flowEndDates   = pDS->GetFlowEndDates();

	for (size_t i (0); i<fwdEndDates->size(); i++)
	{
		fwdEndDates->Elt(i)   = flowEndDates->Elt(CC_Min(i + nper - 1, flowEndDates->size() - 1));
		fwdStartDates->Elt(i) = flowStartDates->Elt(i);
	}

	// to be sure the reset date of numerical calibration 
	// and the last eval date are equal in case in arrears
	if(isArrearsCpn){
		std::vector<double>& resetDates    = pDS->GetResetDates();
		std::vector<double>& eventDates    = itsDateStructure.GetDateStrip(EXER_SCHED)->GetResetDates();
		(*resetDates)[resetDates->size()-1] = (*eventDates)[eventDates->size()-1];
	}
	return pDS;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetCalibration()
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
		ARM_THROW( ERR_INVALID_ARGUMENT," ARM_CallableSnowBallCalculator::CreateAndSetCalibration : wrong model type !" );
	}
	
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CheckCalibMode
///	Returns: void
///	Action : Check the calib mode flag based on the strike sum
/// values
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CheckCalibMode()
{
	for (size_t i = 0; i < itsStrikeSum.size(); ++i)
	{
		if (itsStrikeSum[i] > MIN_STRIKE_VALUE)
			return;
	}

	itsCalibMode = CAP;
	itsBetaCalib = false;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetVolCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for volatility
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CallableSnowBallCalculator::GetVolCalibMethod() const
{
	ARM_CalibMethod* volCalibMethod;

	if( itsCalibMode == SUMOPT && itsBetaCalib)
	{
		volCalibMethod = GetCalibMethod()->GetlinkedMethod();
	}
	else
	{
		volCalibMethod = &(*GetCalibMethod());
	}

#ifdef __GP_STRICT_VALIDATION
    if( volCalibMethod == NULL )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Volatility Bootstrap calib method not found");
#endif

	return volCalibMethod;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetCorrelCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for correlation
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CallableSnowBallCalculator::GetBetaCalibMethod() const
{
	ARM_CalibMethod* betaCalibMethod = NULL;

	if( (itsCalibMode == SUMOPT) && itsBetaCalib)
	{
		betaCalibMethod = &(*GetCalibMethod());
#ifdef __GP_STRICT_VALIDATION
    if( !betaCalibMethod  )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : correl Calib Method not found");
#endif
	}
	return betaCalibMethod;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: ComputePortPrices
///	Returns: nothing
///	Action : compute market target prices of a portfolio
/////////////////////////////////////////////////////////////////

void ARM_CallableSnowBallCalculator::ComputePortPrices(bool isInitParam, CalibMode calibMode, ARM_CalibMethod* calibMethod)
{
	ARM_StdPortfolioPtr PF = calibMethod->GetPortfolio();
	double asOf = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_Date asOfDate(asOf);

	/// Get the B&S model and vol curves from the market data manager
    ARM_BSModel* BSModel = NULL;
	
	if ((calibMode == CAP) || (calibMode == SUMOPT))
	{
		BSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	}
	else if (calibMode == SWAPTION)
	{
		BSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT, " error on calibMode type" );

	int nbProducts = PF->GetSize();
	double expiry, tenor, fwdRate, vol, vega, price, strike;

	std::vector<double> initTimes(nbProducts);
	std::vector<double> initSigmas(nbProducts);
	std::vector<double> sigmaLB(nbProducts);
	std::vector<double> sigmaUB(nbProducts);

	if (calibMode == CAP)
	{
		// Set the exercise strikes
		for (size_t i = 0; i < nbProducts; ++i)
		{
			ARM_CapFloor* capFloor=static_cast< ARM_CapFloor* >(PF->GetAsset(i));
			capFloor->SetModel(BSModel);

			expiry = (capFloor->GetExpiryDate().GetJulian()-asOfDate.GetJulian());
			if (capFloor->GetSwapLeg()->GetFwdRates() 
					&& capFloor->GetSwapLeg()->GetFwdRates()->GetSize() >= 1)
			{
				fwdRate = capFloor->GetSwapLeg()->GetFwdRates()->Elt(0);
			}
			else
			{
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : cap floor fwd rates are not available !" );
			}

			// To maintain the kernel compatibility
			if (capFloor->GetStrike() != K_MARKET_RATE)
				strike = capFloor->GetStrike();
			else
				strike = fwdRate;

			tenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();
			vol = BSModel->ComputeVol(
			expiry/K_YEAR_LEN,
			tenor,
			fwdRate,
			strike);

			initTimes[i]	= expiry;
			initSigmas[i]	= vol/100;
			sigmaLB[i]		= SIGMA_LOWER_BOUND;
			sigmaUB[i]		= SIGMA_UPPER_BOUND;

			vega=capFloor->ComputeSensitivity(K_VEGA);
			vega=fabs(vega);
			price=capFloor->ComputePrice();

			PF->SetPrice(price, i);
			PF->SetPrecision(vega*0.0001, i);
		}
	}
	else if (calibMode == SWAPTION)
	{
		// Set the exercise strikes
		for (size_t i = 0; i < nbProducts; ++i)
		{
			ARM_Swaption* swaption=static_cast< ARM_Swaption* >(PF->GetAsset(i));
			swaption->SetModel(BSModel);

			double atmStrike = swaption->PriceToRate(asOfDate, 0.0);

			expiry = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian());

			if (expiry > 0)
			{
				if (expiry > 0)
				{
					vol = BSModel->ComputeVol(
						expiry/K_YEAR_LEN,
						(itsEndDate.GetJulian()-swaption->GetStartDate().GetJulian())/K_YEAR_LEN,
						atmStrike,
						swaption->GetStrike());
				}
				else
				{
					vol = 0;
				}

				initTimes[i] = expiry;
				initSigmas[i] = vol/100;
				sigmaLB[i] = SIGMA_LOWER_BOUND;
				sigmaUB[i] = SIGMA_UPPER_BOUND;

			/// Compute target price & vega
				double price=swaption->ComputePrice();
				double vega=swaption->ComputeSensitivity(K_VEGA);
				vega=fabs(vega);
				price=swaption->ComputePrice();

				PF->SetPrice(price,i);
				PF->SetPrecision(vega, i);
			}
		}
	}
	else if (calibMode == SUMOPT)
	{
		CreateSumOptionModel();

		// Set the exercise strikes
		for (size_t i = 0; i < nbProducts; ++i)
		{
			ARM_SumOpt* sumOpt=static_cast< ARM_SumOpt* >(PF->GetAsset(i));

			expiry = (sumOpt->GetResetDates()->Elt(0)-asOfDate.GetJulian());

			if (expiry > 0)
			{
				vol = BSModel->ComputeVol(
						expiry/K_YEAR_LEN,
						StringMaturityToYearTerm(itsCpnIdxTerm),
						0.05,
						0.05);
			}
			else
			{
				vol = 0;
			}

			initTimes[i] = expiry;
			initSigmas[i] = vol/100;
			sigmaLB[i] = SIGMA_LOWER_BOUND;
			sigmaUB[i] = SIGMA_UPPER_BOUND;

			UpdateSumOptionStrike(sumOpt->GetAvgStrike());
			double price=ComputeSumOptionPrice(sumOpt);

			PF->SetPrice(price,i);
			PF->SetPrecision(1.0, i);
		}
	}

	if (isInitParam)
	{
		/// Test if volatility must be initialise. In this case,
		/// its schedule contents the expiry date of each swaption product
		size_t sigmaIdx,paramSize=calibMethod->GetCalibParams().size();
		for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
			if(calibMethod->GetCalibParam(sigmaIdx) &&
			(calibMethod->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
				break;
		ARM_CurveModelParam* sigma = new ARM_CurveModelParam(
			ARM_ModelParamType::Volatility,
			&initSigmas,
			&initTimes,
			"SIGMA",
			"STEPUPRIGHT",
			&sigmaLB,
			&sigmaUB);

		if(sigmaIdx >= paramSize || paramSize == 0)
			calibMethod->GetCalibParams().push_back(sigma);
		else
		{
			delete calibMethod->GetCalibParam(sigmaIdx);
			(calibMethod->GetCalibParams())[sigmaIdx] = sigma;
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateSFRMModelParams
///	Returns: ARM_ModelParamsSFRM*
///	Action : creates the model params
/////////////////////////////////////////////////////////////////
ARM_ModelParamsSFRM* ARM_CallableSnowBallCalculator::CreateSFRMModelParams( CalibMode calibMode )
{
	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	std::vector<double> defaultTimes;
	std::vector<double> defaultSigmas;
	std::vector<double> defaultBetas;

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
	paramVector[2] = betaParam;
	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams;
	if ( (calibMode == CAP) || (calibMode == SUMOPT))
		pSFRMModelParams =
			ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,K_DIAG);
	else
		pSFRMModelParams =
			ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,SFRM_NB_FACTORS,K_ROW);

	return pSFRMModelParams;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: BuildModel_SFRM2F
///	Returns: void
///	Action : creates the SFRM2F model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_CallableSnowBallCalculator::BuildModel_SFRM2F(ARM_ZeroCurve* curve)
{
	/// create model params
	ARM_ModelParamsSFRM* pSFRMModelParams = CreateSFRMModelParams( itsCalibMode );

	/// Build the default stochastic model of the calculator : SFRM 2F
	ARM_SFRM* sfrm = new ARM_SFRM( CreateClonedPtr( curve ), *pSFRMModelParams);

	// Delte the model params because it is cloned in the model
	delete pSFRMModelParams;

	// Fix the SFRM Model on the start of the first forward
	ARM_Date startDate(itsStartDatesVec[0]);
	if (itsCpnIndexTiming == K_ARREARS)
	{	
		startDate.AddPeriod(-itsCpnFreq,const_cast<char*>(itsCpnPayCal.c_str()));
	}
	sfrm->SetFixStartDate(startDate);
	
	return sfrm;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: BuildModel_SBGM
///	Returns: void
///	Action : creates the SBGM model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_CallableSnowBallCalculator::BuildModel_SBGM(ARM_ZeroCurve* pCurve)
{
	int nbSteps      = int( (itsHkData.size()>0)?itsHkData[0]:PDE_NB_STEPS );
	int gridSize	 = int( (itsHkData.size()>1)?itsHkData[1]:PDE_GRIDSIZE );
	double nbStdDevs =      (itsHkData.size()>2)?itsHkData[2]:PDE_NB_STDEV ;

	ARM_PricingModel* pModel = new ARM_SmiledFRM(CreateClonedPtr( pCurve ),NULL,nbSteps,gridSize,nbStdDevs,true);

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
	if (itsHkData.size()>3)
		modelParamsSmiledFRM = new ARM_ModelParamsSmiled(modelParams, int(itsHkData[3]));
	else
		modelParamsSmiledFRM = new ARM_ModelParamsSmiled(modelParams,SBGM_NB_FACTORS);

	pModel->SetModelParams(*modelParamsSmiledFRM);
	delete modelParamsSmiledFRM;
	return pModel;

}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: BuildModel_HK
///	Returns: void
///	Action : creates the HK model
/////////////////////////////////////////////////////////////////
ARM_PricingModel* ARM_CallableSnowBallCalculator::BuildModel_HK(ARM_ZeroCurve* pCurve)
{
	std::vector<double> defaultTimes(1,0.);
	std::vector<double> defaultSigmas(1,0.01);
		
	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam* mrsParam = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsKey]));
	if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
	{
		CC_Ostringstream os;
		os << ARM_USERNAME << " : MRS Param for key=" << GetKeys()[MrsKey] << " is expected in the Market Data Manager";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
	}	

	ARM_ModelParamVector modelParams(2);
	modelParams[0] = &volParam;
	modelParams[1] = mrsParam;
	   
	/// Create Markov Functional
	ARM_PricingModel*pModel = new ARM_MarkovFunctional( CreateClonedPtr( pCurve ) );
	
	/// Sets model params
	pModel->SetModelParams(ARM_ModelParamsMF(modelParams));
	
	return pModel;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateAndSetModel()
{
	 /// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	/// model to be created
	ARM_PricingModel* model = NULL;
	
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

	/// antithetic variates!
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

	ARM_ExerciseBoundaryCalc* exerBoundCal;
	if (itsUSMethod == "LS")
	{
		exerBoundCal = new ARM_AMCLongstaffSchwartz(itsNbPathBounding);
	}
	else
	{
		exerBoundCal = new ARM_AMCAndersen(itsNbPathBounding, true);
	}

	ARM_TimeStepPerYearScheduler scheduler(itsFixSteps);

	ARM_AMCMethod* amcMethod = NULL;

	// Path Scheme
	int pathSchemeType = ARM_ArgConv_PathSchemeType.GetNumber(itsPathScheme);
	ARM_PathSchemePtr pathScheme(ARM_PathSchemeFactory.Instance()->CreatePathScheme(pathSchemeType));
		
	/// build model
	if (itsModelType == ARM_PricingModelType::SFRM2F)
	{
		model = BuildModel_SFRM2F(curve);
		ARM_NormalCentredSamplerND sampler(&scheduler);
		amcMethod = new ARM_AMCMethod (itsNbPathPricing,mixteRandGen,&sampler,exerBoundCal,itsNbMaxBucket,ARM_ImpSamplerPtr(NULL),pathScheme);
	}
	else if (itsModelType == ARM_PricingModelType::SBGM)
	{
		model = BuildModel_SBGM(curve);
		ARM_NormalCentredSamplerND sampler(&scheduler);
		amcMethod = new ARM_AMCMethod (itsNbPathPricing,mixteRandGen,&sampler,exerBoundCal,itsNbMaxBucket,ARM_ImpSamplerPtr(NULL),pathScheme);
	}
	else if (itsModelType == ARM_PricingModelType::HK)
	{
		model = BuildModel_HK(curve);
		ARM_MeanRevertingSamplerND sampler(&scheduler);
		amcMethod = new ARM_AMCMethod (itsNbPathPricing,mixteRandGen,&sampler,exerBoundCal,itsNbMaxBucket,ARM_ImpSamplerPtr(NULL),pathScheme);
	}
	else
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" ARM_CallableSnowBallCalculator::CreateAndSetModel : only SFRM2F, SBGM and HK are supported !" );
	}

	/// don't forgeteu the deleteu
	delete exerBoundCal;

	/// ptr-ize model
	ARM_PricingModelPtr refModel( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(model)) );
	
	ARM_StringVector names (1);
	vector<ARM_PricingModelPtr> models (1);
	ARM_StringVectorVector depends(1);

	 // Key names
    names[0] = GetKeys()[YcKey];

	/// models
    models[0] = refModel;

	if(IsBasis())
	{
		/// Push back model pricing names
        names.push_back(GetKeys()[YcBasisDomKey]);

		// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
		models.push_back(ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr(basisCurve)))));
		depends.push_back(ARM_StringVector (1,names[0]));
	}
	ARM_ModelNameMap modelMap (names, models, depends);

	///	modelMap & correls are cloned in multi assets model
	ARM_PricingModelPtr hybridmodel (new ARM_MultiAssetsModel ( &modelMap) );


    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    hybridmodel->SetNumeraire(numeraire);
	hybridmodel->SetNumMethod(ARM_NumMethodPtr( amcMethod ) );

	/// Set the model
	SetPricingModel(hybridmodel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CallableSnowBallCalculator::GetIndexType()
{
    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");

	string cpnIndexTerm=itsCpnIdxTerm;

	if(cpnIndexTerm=="12M")
        cpnIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M

	int cpnIndexFreq=ARM_ArgConv_MatFrequency.GetNumber(cpnIndexTerm);
	int freq = CC_Max(cpnIndexFreq, itsCpnFreq);

    string term = ARM_ArgConvReverse_MatFrequency.GetString(freq);
    liborTypeName += term;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CallableSnowBallCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CallableSnowBallCalculator::DatesStructure() const
{
	return itsDateStructure;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateCapPortfolio
///	Returns: nothing
///	Action : Create caplets portfolio
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CallableSnowBallCalculator::CreateCapPortfolio(bool atmOrStacky)
{	
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	string cpnIndexTerm(itsCpnIdxTerm);
	if(cpnIndexTerm=="12M")
        cpnIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
	int cpnIndexFreq=ARM_ArgConv_MatFrequency.GetNumber(cpnIndexTerm);
	
	// builds the Cpn schedule
	ARM_DateStrip dateStrip(itsStartDate,
							   itsEndDate,
							   cpnIndexFreq,
							   itsCpnDayCount,
							   itsCpnResetCal.c_str(),
							   K_MOD_FOLLOWING,
							   itsCpnIntRule, 
							   K_SHORTSTART,
							   itsCpnResetLag,
							   -1,
							   GETDEFAULTVALUE,
							   itsCpnPayCal.c_str(),
							   itsCpnIndexTiming,
							   K_ARREARS);


	ARM_INDEX_TYPE indexType = (ARM_INDEX_TYPE) FromIndexAndTermToIndexType(itsCpnIdxTerm, GetCurrencyUnit()->GetCcyName());
	int resetFreq = 1.0/StringMaturityToYearTerm(itsCpnIdxTerm);

	std::vector<double>& startDates = dateStrip.GetFwdRateStartDates();
	int nbProducts = startDates->size();

	vector<ARM_Security*> securities;
	vector<double> mktPrices;
	vector<double> weights;
	vector<double> precisions;

	double strike = 0;
	for(size_t i=0;i<nbProducts;i++)
	{
		strike = 0;
		if (!atmOrStacky)
		{
			strike += itsStrikeOptVec[i];
			strike *= 100;
		}
		else
		{
			strike = -1.0;
		}

		ARM_Date indexStartDate((*startDates)[i]);

		 /// Compute not adjusted index end date to avoid problem in caplet building
		ARM_Date indexEndDate(indexStartDate);
		indexEndDate.AddPeriod(itsCpnIdxTerm,itsCpnPayCal.c_str());

		/// Build the cap (ARM_Security & VanillaArg versions)
		ARM_Security* sec = new  ARM_CapFloor(
			indexStartDate,
			indexEndDate,
			K_CAP,
			strike,
			indexType,
			0.0,
			K_DEF_FREQ,
			K_DEF_FREQ,
			GetCurrencyUnit());

		if(sec->GetResetDates()->size() != 1)
		{
			CC_Ostringstream os;
			os << ARM_USERNAME << " : Implied Caplet/Floorlet #" << i << " has more than one flow";
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
		}

		// Removed products already fixed
		if (sec->GetResetDates()->Elt(0) > asOfDate.GetJulian())
		{
			securities.push_back(sec);
			mktPrices.push_back(i+1);
			weights.push_back(CF_DEFAULT_WEIGHT);
			precisions.push_back(CF_DEFAULT_WEIGHT);
		}
		else
		{
			delete sec;
		}
	}

	ARM_StdPortfolioPtr port(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

	for(i=0;i<securities.size();i++)
		delete securities[i];

	return port;
}

ARM_StdPortfolioPtr ARM_CallableSnowBallCalculator::CreateBetaCapPortfolio()
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	vector<ARM_Security*> securities;
	vector<double> mktPrices;
	vector<double> weights;
	vector<double> precisions;

	int nbFlows = itsResetDatesVec.size();

	ARM_Vector* strikeValues = new ARM_Vector(nbFlows);
	ARM_Vector* nomValues = new ARM_Vector(nbFlows);

	double nominal = 0.0;
	double df;

	for (int i = nbFlows-1; i >= 0; --i)
	{
		if (fabs(itsLevPrevVec[i]) < K_DOUBLE_TOL)
			nominal = 0.0;

		df = GetPricingModel()->GetZeroCurve()->DiscountPrice((itsPaymentDatesVec[i]-asOfDate.GetJulian())/K_YEAR_LEN);
		if (i == nbFlows-1)
			nominal = df*itsLevNewVec[i];
		else
			nominal = df*itsLevNewVec[i] + itsLevPrevVec[i+1]*nominal;

		(*strikeValues)[i] = itsStrikeOptVec[i]*100;
		(*nomValues)[i] = nominal*100/df;
	}

	ARM_Date capStartDate(itsStartDate);
	
	if (itsCpnIndexTiming == K_ARREARS)
	{
		capStartDate.AddPeriod(itsCpnFreq,const_cast<char*>(itsCpnPayCal.c_str()));
	}

	ARM_Date capEndDate(capStartDate);
	capEndDate.AddPeriodMult(itsCpnFreq,nbFlows,const_cast<char*>(itsCpnPayCal.c_str()));

	ARM_INDEX_TYPE indexType = (ARM_INDEX_TYPE) FromIndexAndTermToIndexType(itsCpnIdxTerm, GetCurrencyUnit()->GetCcyName());

	/// Build the cap (ARM_Security & VanillaArg versions)
	ARM_CapFloor* sec = new  ARM_CapFloor(
		capStartDate,
		capEndDate,
		itsCapOrFloor,
		0.0,
		indexType,
		0.0,
		itsCpnFreq,
		itsCpnFreq,
		GetCurrencyUnit());

	ARM_Vector* resetDates = static_cast<ARM_Vector*>(sec->GetSwapLeg()->GetResetDates()->Clone());
	ARM_Vector* payDates = static_cast<ARM_Vector*>(sec->GetSwapLeg()->GetPaymentDates()->Clone());

	ARM_ReferenceValue strikes(resetDates, strikeValues);
	ARM_ReferenceValue nominals(payDates, nomValues);

	sec->SetStrikes(&strikes);
	sec->SetAmount(&nominals);

	securities.push_back(sec);
	mktPrices.push_back(1);
	weights.push_back(CF_DEFAULT_WEIGHT);
	precisions.push_back(CF_DEFAULT_WEIGHT);

	ARM_StdPortfolioPtr port(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

	delete sec;

	return port;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: UpdateBetaCapNominal
///	Returns: nothing
///	Action : Create sum option portfolio
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::UpdateBetaCapNominal(ARM_CalibMethod* calibMethod)
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	ARM_StdPortfolioPtr PF = calibMethod->GetPortfolio();

	ARM_CapFloor* capFloor=static_cast< ARM_CapFloor* >(PF->GetAsset(0));

	int nbFlows = itsResetDatesVec.size();

	ARM_Vector* nomValues = new ARM_Vector(nbFlows);

	double nominal = 0.0;
	double df;

	ARM_ZeroCurve* zcCurve = static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	for (int i = nbFlows-1; i >= 0; --i)
	{
		if (fabs(itsLevPrevVec[i]) < K_DOUBLE_TOL)
			nominal = 0.0;

		df = zcCurve->DiscountPrice((itsPaymentDatesVec[i]-asOfDate.GetJulian())/K_YEAR_LEN);
		if (i == nbFlows-1)
			nominal = df*itsLevNewVec[i];
		else
			nominal = df*itsLevNewVec[i] + itsLevPrevVec[i+1]*nominal;
		(*nomValues)[i] = nominal*100/df;
	}

	ARM_Vector* payDates = static_cast<ARM_Vector*>(capFloor->GetSwapLeg()->GetPaymentDates()->Clone());
	ARM_ReferenceValue nominals(payDates, nomValues);

	capFloor->SetAmount(&nominals);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateSumOptionPortfolio
///	Returns: nothing
///	Action : Create sum option portfolio
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_CallableSnowBallCalculator::CreateSumOptionPortfolio()
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	size_t nbProducts = itsStartDatesVec.size();
	
	size_t idx;

	std::vector<double> coeffs;

	vector<ARM_Security*> securities;
	vector<double> mktPrices;
	vector<double> weights;
	vector<double> precisions;

	size_t firstIdx = 0;

	double strike = 0;

	for (idx = 0; idx < nbProducts; ++idx)
	{
		
		if (fabs(itsLevPrevVec[idx]) < K_DOUBLE_TOL)
		{
			firstIdx = idx;
		}

		if (fabs(itsLevNewVec[idx]) > K_DOUBLE_TOL)
		{
			strike += itsLevNewVec[idx]*itsStrikeOptVec[idx];

			ARM_SumOpt* sec = new ARM_SumOpt(
				itsResetDatesVec,
				itsStartDatesVec,
				itsEndDatesVec,
				itsPaymentDatesVec,
				itsPeriodsVec,
				firstIdx,
				idx+1,
				-1,
				itsStrikeSum[idx],
				itsCapOrFloor,
				*(itsSumOptionCoeffs.GetRow(idx)));

			// Removed products already fixed
			if ((sec->GetResetDates()->Elt(0) > asOfDate.GetJulian()) && (itsStrikeSum[idx] > MIN_STRIKE_VALUE))
			{
				securities.push_back(sec);
				weights.push_back(CF_DEFAULT_WEIGHT);
				precisions.push_back(CF_DEFAULT_WEIGHT);
				mktPrices.push_back(CF_DEFAULT_PRICE);

			}
			else
			{
				delete sec;
			}
		}
	}

	ARM_StdPortfolioPtr port(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

	for (idx = 0; idx < securities.size(); ++idx)
		delete securities[idx];

	return port;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateSwoptPortfolio
///	Returns: nothing
///	Action : Create swaption portfolio
/////////////////////////////////////////////////////////////////

ARM_StdPortfolioPtr ARM_CallableSnowBallCalculator::CreateSwoptPortfolio()
{
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	int stdIndexFreq = GetCurrencyUnit()->GetFixedPayFreq();

	
	// builds the Cpn schedule
	ARM_DateStrip swoptDateStrip(itsStartDate,
							   itsEndDate,
							   (itsCalibSwoptFreq!=-1?itsCalibSwoptFreq:stdIndexFreq),
							   itsCpnDayCount,
							   itsCpnResetCal.c_str(),
							   K_MOD_FOLLOWING,
							   itsCpnIntRule, 
							   K_SHORTSTART,
							   itsCpnResetLag,
							   -1,
							   GETDEFAULTVALUE,
							   itsCpnPayCal.c_str(),
							   K_ADVANCE,
							   K_ARREARS);


	// We take the payment frequency linked to the currency
	ARM_INDEX_TYPE indexType = (ARM_INDEX_TYPE) FromIndexAndTermToIndexType(itsCpnIdxTerm, GetCurrencyUnit()->GetCcyName());
	int resetFreq = 1.0/StringMaturityToYearTerm(itsCpnIdxTerm);

	std::vector<double>& startDates = swoptDateStrip.GetFwdRateStartDates();
	int nbProducts = startDates->size();

	vector<ARM_Security*> securities;
	vector<double> mktPrices;
	vector<double> weights;
	vector<double> precisions;

	int flowIdx = 0;

	for(size_t i=0;i<nbProducts;i++)
	{
		ARM_Date indexStartDate((*startDates)[i]);
		ARM_Date indexEndDate(itsEndDate);

		/// Build the swaption (ARM_Security & VanillaArg versions)
		ARM_Swap swap(indexStartDate, indexEndDate,
			indexType,0.0,1.0,K_RCV,resetFreq,resetFreq,GetCurrencyUnit());

		ARM_Date expiryDate((*(swap.GetFloatLeg()->GetResetDates()))[0]);

		ARM_Security* sec = new ARM_Swaption(&swap,K_RCV,K_EUROPEAN,-1.0,expiryDate);

		// The std swaption reset freq could be different to deal reset freq
		while ((flowIdx<itsResetDatesVec.size()) && (itsResetDatesVec[flowIdx] < expiryDate.GetJulian() + K_NEW_DOUBLE_TOL))
			flowIdx++;

		// Removed products already fixed
		if ((flowIdx >= itsNbNonCall) && (flowIdx<nbProducts)
			&& (itsFeesVec[flowIdx] < NON_CALL_FEES)
			&& (expiryDate.GetJulian() > asOfDate.GetJulian()))
		{
			securities.push_back(sec);
			mktPrices.push_back(i+1);
			weights.push_back(CF_DEFAULT_WEIGHT);
			precisions.push_back(CF_DEFAULT_WEIGHT);
		}
		else
		{
			delete sec;
		}
	}

	ARM_StdPortfolioPtr port(new ARM_StdPortfolio(securities, weights, mktPrices, precisions));

	for (i = 0; i < securities.size(); ++i)
		delete securities[i];

	return port;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: GetRefModel
///	Returns: ARM_PricingModelPtr
///	Action : return the reference model of the Callable SnowBall
/////////////////////////////////////////////////////////////////
ARM_SFRM* ARM_CallableSnowBallCalculator::GetRefModel()
{
	ARM_SFRM* refModel;

	/// Allow update only for default model (H&W1F)
    refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());
    if( !refModel )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : only H&W1F model is allowed for updating");

	return refModel;
}




/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[FundingPaid] = zeroValue;
    rowTypeVec[FundingPaid] = ARM_DOUBLE;

    rowDescVec[CF] = zeroValue;
    rowTypeVec[CF] = ARM_DOUBLE;

	rowDescVec[Bermuda] = zeroValue;
    rowTypeVec[Bermuda] = ARM_DOUBLE;

	rowDescVec[LBCouponPaid] = zeroValue;
    rowTypeVec[LBCouponPaid] = ARM_DOUBLE;

	rowDescVec[UBCouponPaid] = zeroValue;
    rowTypeVec[UBCouponPaid] = ARM_DOUBLE;

	rowDescVec[PaidCoupon] = zeroValue;
    rowTypeVec[PaidCoupon] = ARM_DOUBLE;

	rowDescVec[NextCoupon] = zeroValue;
    rowTypeVec[NextCoupon] = ARM_DOUBLE;

	rowDescVec[CouponSimple] = zeroValue;
    rowTypeVec[CouponSimple] = ARM_DOUBLE;

	rowDescVec[CouponOption] = zeroValue;
    rowTypeVec[CouponOption] = ARM_DOUBLE;

	rowDescVec[CouponOptPaid] = zeroValue;
    rowTypeVec[CouponOptPaid] = ARM_DOUBLE;

	rowDescVec[LBCouponAvgPaid] = zeroValue;
    rowTypeVec[LBCouponAvgPaid] = ARM_DOUBLE;

	rowDescVec[AnaLBCouponAvg] = zeroValue;
    rowTypeVec[AnaLBCouponAvg] = ARM_DOUBLE;

	rowDescVec[AnaLBCoupon] = zeroValue;
    rowTypeVec[AnaLBCoupon] = ARM_DOUBLE;

	rowDescVec[CouponOptAvgPaid] = zeroValue;
    rowTypeVec[CouponOptAvgPaid] = ARM_DOUBLE;

	rowDescVec[NextSimpleCoupon] = zeroValue;
    rowTypeVec[NextSimpleCoupon] = ARM_DOUBLE;

	rowDescVec[NextCouponCap] = zeroValue;
    rowTypeVec[NextCouponCap] = ARM_DOUBLE;

	rowDescVec[NextCouponFloor] = zeroValue;
    rowTypeVec[NextCouponFloor] = ARM_DOUBLE;

	rowDescVec[NextCouponOpt] = zeroValue;
    rowTypeVec[NextCouponOpt] = ARM_DOUBLE;

	rowDescVec[Coupon] = zeroValue;
    rowTypeVec[Coupon] = ARM_DOUBLE;

	rowDescVec[CF1] = zeroValue;
    rowTypeVec[CF1] = ARM_DOUBLE;

	rowDescVec[CF2] = zeroValue;
    rowTypeVec[CF2] = ARM_DOUBLE;

	rowDescVec[AMCIndex] = zeroValue;
    rowTypeVec[AMCIndex] = ARM_DOUBLE;

	rowDescVec[Option2] = zeroValue;
    rowTypeVec[Option2] = ARM_DOUBLE;

	rowDescVec[Bermuda2] = zeroValue;
    rowTypeVec[Bermuda2] = ARM_DOUBLE;

}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CallableSnowBallCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
    size_t descSize = sizeof(CallableSnowBallColNamesTable)/sizeof(CallableSnowBallColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// Get the model names for coupon & funding legs descritpion
    string cpnModelName = GetKeys()[YcKey];
	string payModelName = GetKeys()[YcBasisDomKey];

    /// Convention conversion : K_QUATERLY => 3M for instance
	string fundFreq=ARM_ArgConvReverse_MatFrequency.GetString(itsFundFreq);
    string fundDayCount=ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount);

	string cpnFreq=ARM_ArgConvReverse_MatFrequency.GetString(itsCpnFreq);
    string cpnIndexDayCount=ARM_ArgConvReverse_DayCount.GetString(itsCpnIndexDayCount);
	string cpnIndexTiming =  itsCpnIndexTiming==K_ADVANCE ? "ADV":"ARR";
	string capFloor = (itsCapOrFloor==1 ?"C":"F");

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

	ARM_VectorPtr resetDates = datesStructure.GetMergeData();
	std::vector<double>& cpnResetDates = datesStructure.GetDateStrip(CPN_SCHED)->GetResetDates();

	std::vector<double>& cpnStartDates = datesStructure.GetDateStrip(CPN_SCHED)->GetFlowStartDates();
	std::vector<double>& cpnEndDates = datesStructure.GetDateStrip(CPN_SCHED)->GetFlowEndDates();
	std::vector<double>& cpnFwdStartDates = datesStructure.GetDateStrip(CPN_SCHED)->GetFwdRateStartDates();
	std::vector<double>& cpnFwdEndDates = datesStructure.GetDateStrip(CPN_SCHED)->GetFwdRateEndDates();
	std::vector<double>& cpnPayDates = datesStructure.GetDateStrip(CPN_SCHED)->GetPaymentDates();
	std::vector<double>& cpnIT = datesStructure.GetDateStrip(CPN_SCHED)->GetInterestTerms();

	std::vector<double>& fundingStartDates = datesStructure.GetDateStrip(EXER_SCHED)->GetFlowStartDates();
	std::vector<double>& fundingEndDates = datesStructure.GetDateStrip(EXER_SCHED)->GetFlowEndDates();

	std::vector<double>& exerDates = datesStructure.GetDateStrip(EXER_SCHED)->GetResetDates();

    /// EventDate description
    double eventDate=(*resetDates)[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;

	CC_Ostringstream indexDesc;
    indexDesc << CC_NS(std,fixed) << eventIdx;
    rowDescVec[Index] = indexDesc.str();
    rowTypeVec[Index] = ARM_DOUBLE;

    double cpnFwdStartDate=(*cpnFwdStartDates)[eventIdx];
	double cpnStartDate=(*cpnStartDates)[eventIdx];
	double cpnEndDate=(*cpnEndDates)[eventIdx];
	double fundingStartDate=(*fundingStartDates)[eventIdx];
	

	if (itsSBPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
    {
		// DealStartDate, EndDate
		CC_Ostringstream sumOptStartDateDesc;
		sumOptStartDateDesc << CC_NS(std,fixed) << itsStartDate.GetJulian();
		rowDescVec[SumOptStartDate] = sumOptStartDateDesc.str();
		rowTypeVec[SumOptStartDate] = ARM_DATE_TYPE;

		CC_Ostringstream endDateDesc;
		endDateDesc << CC_NS(std,fixed) << cpnEndDate;
		rowDescVec[EndDate] = endDateDesc.str();
		rowTypeVec[EndDate] = ARM_DATE_TYPE;

		// Cpn Fwd StartDate desc
		CC_Ostringstream cpnFwdStartDateDesc;
		cpnFwdStartDateDesc << CC_NS(std,fixed) << cpnFwdStartDate;
		rowDescVec[CpnFwdStartDate] = cpnFwdStartDateDesc.str();
		rowTypeVec[CpnFwdStartDate] = ARM_DATE_TYPE;

		ARM_Date dFwdEndDate((*cpnFwdStartDates)[eventIdx]);
		dFwdEndDate.AddPeriod(itsCpnIdxTerm);
		CC_Ostringstream cpnFwdEndDateDesc;
		cpnFwdEndDateDesc << CC_NS(std,fixed) << dFwdEndDate.GetJulian();
		rowDescVec[CpnFwdEndDate] = cpnFwdEndDateDesc.str();
		rowTypeVec[CpnFwdEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream cpnPayDateDesc;
		cpnPayDateDesc << CC_NS(std,fixed) << (*cpnPayDates)[eventIdx];
		rowDescVec[CpnPayDate] = cpnPayDateDesc.str();
		rowTypeVec[CpnPayDate] = ARM_DATE_TYPE;
	
		CC_Ostringstream cpnITDesc;
		cpnITDesc << CC_NS(std,fixed) << (*cpnIT)[eventIdx];
		rowDescVec[CpnIT] = cpnITDesc.str();
		rowTypeVec[CpnIT] = ARM_DOUBLE;

		CC_Ostringstream ConstDesc;
		ConstDesc << CC_NS(std,fixed) << itsConstVec[itsSBPos2Rank[eventIdx]];

		rowDescVec[Const] = ConstDesc.str();
		rowTypeVec[Const] = ARM_STRING;

		CC_Ostringstream LPrevCpnDesc;
		LPrevCpnDesc << CC_NS(std,fixed) << itsLevPrevVec[itsSBPos2Rank[eventIdx]];
		rowDescVec[LevPrevCpn] = LPrevCpnDesc.str();
		rowTypeVec[LevPrevCpn] = ARM_DOUBLE;

		CC_Ostringstream LNewOptDesc;
		LNewOptDesc << CC_NS(std,fixed) << itsLevNewVec[itsSBPos2Rank[eventIdx]];
		rowDescVec[LevNewOpt] = LNewOptDesc.str();
		rowTypeVec[LevNewOpt] = ARM_DOUBLE;

		CC_Ostringstream LStrikeOptDesc;
		LStrikeOptDesc << CC_NS(std,fixed) << itsStrikeOptVec[itsSBPos2Rank[eventIdx]];
		rowDescVec[StrikeOpt] = LStrikeOptDesc.str();
		rowTypeVec[StrikeOpt] = ARM_DOUBLE;	

		CC_Ostringstream minCpnDesc;
		minCpnDesc << CC_NS(std,fixed) << itsCpnMinVec[itsSBPos2Rank[eventIdx]];
		rowDescVec[MinCpn] = minCpnDesc.str();
		rowTypeVec[MinCpn] = ARM_DOUBLE;	

		CC_Ostringstream maxCpnDesc;
		maxCpnDesc << CC_NS(std,fixed) << itsCpnMaxVec[itsSBPos2Rank[eventIdx]];
		rowDescVec[MaxCpn] = maxCpnDesc.str();
		rowTypeVec[MaxCpn] = ARM_DOUBLE;	
	}

	if (itsExerPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
	{
		CC_Ostringstream fundingStartDateDesc;
		fundingStartDateDesc << CC_NS(std,fixed) << fundingStartDate;
		rowDescVec[FundingStartDate] = fundingStartDateDesc.str();
		rowTypeVec[FundingStartDate] = ARM_DATE_TYPE;

		CC_Ostringstream fundingEndDateDesc;
		fundingEndDateDesc << CC_NS(std,fixed) << (*fundingEndDates)[eventIdx];
		rowDescVec[FundingEndDate] = fundingEndDateDesc.str();
		rowTypeVec[FundingEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream productEndDateDesc;
		productEndDateDesc << CC_NS(std,fixed) << itsEndDate.GetJulian();
		rowDescVec[ProductEndDate] = productEndDateDesc.str();
		rowTypeVec[ProductEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream dfNumDesc;
		dfNumDesc << "DF(" << payModelName<<  ","<< CallableSnowBallColNamesTable[ProductEndDate] << "[i])";
		rowDescVec[DFNum] = dfNumDesc.str();
		rowTypeVec[DFNum] = ARM_STRING;

		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << itsFeesVec[itsExerPos2Rank[eventIdx]];
		rowDescVec[Fees] = feesDesc.str();
		rowTypeVec[Fees] = ARM_DOUBLE;
	
		CC_Ostringstream nextNominalDesc;
		nextNominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNotionalProfile).Interpolate((*cpnPayDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextNominal] = nextNominalDesc.str();
		rowTypeVec[NextNominal] = ARM_DOUBLE;

		CC_Ostringstream dfNextPayDesc;
		dfNextPayDesc << "DF(" << payModelName << "," << CallableSnowBallColNamesTable[NextCpnPayDate] << "[i])";
		rowDescVec[DFNextPay] = dfNextPayDesc.str();
		rowTypeVec[DFNextPay] = ARM_STRING;

		CC_Ostringstream nextCpnFwdStartDateDesc;
		nextCpnFwdStartDateDesc << CC_NS(std,fixed) << (*cpnFwdStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]];
		rowDescVec[NextCpnFwdStartDate] = nextCpnFwdStartDateDesc.str();
		rowTypeVec[NextCpnFwdStartDate] = ARM_DATE_TYPE;
	
		ARM_Date dFwdEndDate((*cpnFwdStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]);
		dFwdEndDate.AddPeriod(itsCpnIdxTerm);
		CC_Ostringstream nextCpnFwdEndDateDesc;
		nextCpnFwdEndDateDesc << CC_NS(std,fixed) << dFwdEndDate.GetJulian();
		rowDescVec[NextCpnFwdEndDate] = nextCpnFwdEndDateDesc.str();
		rowTypeVec[NextCpnFwdEndDate] = ARM_DATE_TYPE;
		
		CC_Ostringstream nextCpnPayDateDesc;
		nextCpnPayDateDesc << CC_NS(std,fixed) << (*cpnPayDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]];
		rowDescVec[NextCpnPayDate] = nextCpnPayDateDesc.str();
		rowTypeVec[NextCpnPayDate] = ARM_DATE_TYPE;
	
		CC_Ostringstream nextCpnITDesc;
		nextCpnITDesc << CC_NS(std,fixed) << (*cpnIT)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]];
		rowDescVec[NextCpnIT] = nextCpnITDesc.str();
		rowTypeVec[NextCpnIT] = ARM_DOUBLE;

		CC_Ostringstream nextConstDesc;
		nextConstDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsConstProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);

		rowDescVec[NextConst] = nextConstDesc.str();
		rowTypeVec[NextConst] = ARM_STRING;

		CC_Ostringstream nextLevPrevCpnDesc;
		nextLevPrevCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsLPrevCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextLevPrevCpn] = nextLevPrevCpnDesc.str();
		rowTypeVec[NextLevPrevCpn] = ARM_DOUBLE;

		CC_Ostringstream nextNewOptDesc;
		nextNewOptDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsLNewOptProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextLevNewOpt] = nextNewOptDesc.str();
		rowTypeVec[NextLevNewOpt] = ARM_DOUBLE;

		CC_Ostringstream nextStrikeOptDesc;
		nextStrikeOptDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsStrikeOptProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextStrikeOpt] = nextStrikeOptDesc.str();
		rowTypeVec[NextStrikeOpt] = ARM_DOUBLE;

		CC_Ostringstream nextMinCpnDesc;
		nextMinCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsMinCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextMinCpn] = nextMinCpnDesc.str();
		rowTypeVec[NextMinCpn] = ARM_DOUBLE;

		CC_Ostringstream nextMaxCpnDesc;
		nextMaxCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsMaxCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate);
		rowDescVec[NextMaxCpn] = nextMaxCpnDesc.str();
		rowTypeVec[NextMaxCpn] = ARM_DOUBLE;

		CC_Ostringstream nextCpnIndexDesc;
		nextCpnIndexDesc << "LIBOR(" << cpnModelName << "," << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i]," << cpnIndexDayCount << ",," << CallableSnowBallColNamesTable[NextCpnPayDate] << "[i])";
		rowDescVec[NextCpnIndex]	= nextCpnIndexDesc.str();
		rowTypeVec[NextCpnIndex]= ARM_STRING;

		CC_Ostringstream nextCpnIndexWithoutAdjDesc;
		nextCpnIndexWithoutAdjDesc << "LIBOR(" << cpnModelName << "," << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i]," << cpnIndexDayCount << ")";
		rowDescVec[NextCpnIndexWithoutAdj]	= nextCpnIndexWithoutAdjDesc.str();
		rowTypeVec[NextCpnIndexWithoutAdj]= ARM_STRING;

		if (itsCapOrFloor == K_CAP)
		{

			CC_Ostringstream nextCouponOptDesc;
			nextCouponOptDesc << CallableSnowBallColNamesTable[NextCouponFloor] << "[i]-" << CallableSnowBallColNamesTable[NextCouponCap] << "[i]";

			rowDescVec[NextCouponOpt] = nextCouponOptDesc.str();
			rowTypeVec[NextCouponOpt] = ARM_STRING;

			if (itsExerPos2Rank[eventIdx] > 0)
			{
				CC_Ostringstream nextSimpleCpnDesc;

				nextSimpleCpnDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]+";
				
				nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextCpnIndex]  << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCoupon] = nextSimpleCpnDesc.str();
				rowTypeVec[NextSimpleCoupon] = ARM_STRING;

				CC_Ostringstream nextSimpleCpnWithoutAdjDesc;
				nextSimpleCpnWithoutAdjDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]+";
				
				nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj]  << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCouponWithoutAdj] = nextSimpleCpnWithoutAdjDesc.str();
				rowTypeVec[NextSimpleCouponWithoutAdj] = ARM_STRING;

				CC_Ostringstream strikeCapDesc;
				strikeCapDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeCapDesc << "-(" << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeCapDesc << "+" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeCap] = strikeCapDesc.str();
				rowTypeVec[StrikeCap] = ARM_STRING;

				CC_Ostringstream strikeFloorDesc;

				strikeFloorDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";
				
				if (itsCallSBOrStacky)
					strikeFloorDesc << "-(" << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeFloorDesc << "+" << CallableSnowBallColNamesTable[NextMinCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeFloor] = strikeFloorDesc.str();
				rowTypeVec[StrikeFloor] = ARM_STRING;
			}
			else
			{
				CC_Ostringstream nextSimpleCpnDesc;
				nextSimpleCpnDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*("  << CallableSnowBallColNamesTable[NextCpnIndex]	 << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*("  << CallableSnowBallColNamesTable[NextCpnIndex]	 << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCoupon] = nextSimpleCpnDesc.str();
				rowTypeVec[NextSimpleCoupon] = ARM_STRING;

				CC_Ostringstream nextSimpleWithoutAdjCpnDesc;
				nextSimpleWithoutAdjCpnDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleWithoutAdjCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*("  << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj]	 << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleWithoutAdjCpnDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*("  << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj]	 << "[i]-" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCouponWithoutAdj] = nextSimpleWithoutAdjCpnDesc.str();
				rowTypeVec[NextSimpleCouponWithoutAdj] = ARM_STRING;

				CC_Ostringstream strikeCapDesc;

				strikeCapDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeCapDesc << "-(" << CallableSnowBallColNamesTable[NextConst] << "[i]-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeCapDesc << "+" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeCap] = strikeCapDesc.str();
				rowTypeVec[StrikeCap] = ARM_STRING;

				CC_Ostringstream strikeFloorDesc;

				strikeFloorDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeFloorDesc << "-(" << CallableSnowBallColNamesTable[NextConst] << "[i]-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeFloorDesc << "+" << CallableSnowBallColNamesTable[NextMinCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeFloor] = strikeFloorDesc.str();
				rowTypeVec[StrikeFloor] = ARM_STRING;
			}
		}
		else
		{
			CC_Ostringstream nextCouponOptDesc;
			nextCouponOptDesc << CallableSnowBallColNamesTable[NextCouponCap] << "[i]-" << CallableSnowBallColNamesTable[NextCouponFloor] << "[i]";
			rowDescVec[NextCouponOpt] = nextCouponOptDesc.str();
			rowTypeVec[NextCouponOpt] = ARM_STRING;

			if (itsExerPos2Rank[eventIdx] > 0)
			{
				CC_Ostringstream nextSimpleCpnDesc;

				nextSimpleCpnDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndex] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndex] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCoupon] = nextSimpleCpnDesc.str();
				rowTypeVec[NextSimpleCoupon] = ARM_STRING;

				CC_Ostringstream nextSimpleCpnWithoutAdjDesc;

				nextSimpleCpnWithoutAdjDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCouponWithoutAdj] = nextSimpleCpnWithoutAdjDesc.str();
				rowTypeVec[NextSimpleCouponWithoutAdj] = ARM_STRING;

				CC_Ostringstream strikeCapDesc;

				strikeCapDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeCapDesc << "+(" << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeCapDesc << "-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeCap] = strikeCapDesc.str();
				rowTypeVec[StrikeCap] = ARM_STRING;

				CC_Ostringstream strikeFloorDesc;

				strikeFloorDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeFloorDesc << "+(" << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeFloorDesc << "-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeFloor] = strikeFloorDesc.str();
				rowTypeVec[StrikeFloor] = ARM_STRING;
			}
			else
			{
				CC_Ostringstream nextSimpleCpnDesc;

				nextSimpleCpnDesc << "(";
				if (itsCallSBOrStacky)
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndex] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleCpnDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndex] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCoupon] = nextSimpleCpnDesc.str();
				rowTypeVec[NextSimpleCoupon] = ARM_STRING;

				CC_Ostringstream nextSimpleCpnWithoutAdjDesc;
				nextSimpleCpnWithoutAdjDesc << "(";

				if (itsCallSBOrStacky)
					nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
				else
					nextSimpleCpnWithoutAdjDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[NextCpnIndexWithoutAdj] << "[i]))*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";

				rowDescVec[NextSimpleCouponWithoutAdj] = nextSimpleCpnWithoutAdjDesc.str();
				rowTypeVec[NextSimpleCouponWithoutAdj] = ARM_STRING;

				CC_Ostringstream strikeCapDesc;

				strikeCapDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeCapDesc << "+(" << CallableSnowBallColNamesTable[NextConst] << "[i]-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeCapDesc << "-" << CallableSnowBallColNamesTable[NextMinCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeCap] = strikeCapDesc.str();
				rowTypeVec[StrikeCap] = ARM_STRING;

				CC_Ostringstream strikeFloorDesc;

				strikeFloorDesc << CallableSnowBallColNamesTable[NextStrikeOpt] << "[i]";

				if (itsCallSBOrStacky)
					strikeFloorDesc << "+(" << CallableSnowBallColNamesTable[NextConst] << "[i]-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i])/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";
				else
					strikeFloorDesc << "-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]/" << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]";

				rowDescVec[StrikeFloor] = strikeFloorDesc.str();
				rowTypeVec[StrikeFloor] = ARM_STRING;
			}
		}

		CC_Ostringstream nextCpnDesc;

		if (!itsCallSBOrStacky)
				nextCpnDesc << CallableSnowBallColNamesTable[NextConst] << "[i]+" << CallableSnowBallColNamesTable[NextLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[PrevCoupon] << "[i]+";

		nextCpnDesc << "(" << CallableSnowBallColNamesTable[NextSimpleCoupon] << "[i]+" << CallableSnowBallColNamesTable[NextCouponOpt] << "[i])*" << CallableSnowBallColNamesTable[NextCpnIT] << "[i]";

		rowDescVec[NextCoupon] = nextCpnDesc.str();
		rowTypeVec[NextCoupon] = ARM_STRING;


		if (const_cast< ARM_Curve& >(itsLNewOptProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]]]-asOfDate) > 0.0)
		{
			CC_Ostringstream nextCpnCapDesc;

			nextCpnCapDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*Caplet(" << payModelName << "," << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i]," << CallableSnowBallColNamesTable[StrikeCap] << "[i],C,A360,," << CallableSnowBallColNamesTable[NextCpnPayDate] << "[i],,,," << cpnIndexDayCount << ")/DCF(" << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i],A360)";
			rowDescVec[NextCouponCap] = nextCpnCapDesc.str();
			rowTypeVec[NextCouponCap] = ARM_STRING;

			CC_Ostringstream nextCpnFloorDesc;
			nextCpnFloorDesc << CallableSnowBallColNamesTable[NextLevNewOpt] << "[i]*Caplet(" << payModelName << "," << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i]," << CallableSnowBallColNamesTable[StrikeFloor] << "[i],F,A360,," << CallableSnowBallColNamesTable[NextCpnPayDate] << "[i],,,," << cpnIndexDayCount << ")/DCF(" << CallableSnowBallColNamesTable[NextCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[NextCpnFwdEndDate] << "[i],A360)";
			rowDescVec[NextCouponFloor] = nextCpnFloorDesc.str();
			rowTypeVec[NextCouponFloor] = ARM_STRING;
		}
		else
		{
			CC_Ostringstream nextCpnCapDesc;
			nextCpnCapDesc << CC_NS(std,fixed) << 0.0;
			rowDescVec[NextCouponCap] = nextCpnCapDesc.str();
			rowTypeVec[NextCouponCap] = ARM_DOUBLE;

			CC_Ostringstream nextCpnFloorDesc;
			nextCpnFloorDesc << CC_NS(std,fixed) << 0.0;
			rowDescVec[NextCouponFloor] = nextCpnFloorDesc.str();
			rowTypeVec[NextCouponFloor] = ARM_DOUBLE;
		}

		if (itsExerPos2Rank[eventIdx] > 0)
		{
			if (itsSBRank2Pos[itsExerPos2Rank[eventIdx] - 1] >= eventIdx)
			{
				CC_Ostringstream prevCpnFwdStartDateDesc;
				prevCpnFwdStartDateDesc << CC_NS(std,fixed) << (*cpnFwdStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]];
				rowDescVec[PrevCpnFwdStartDate] = prevCpnFwdStartDateDesc.str();
				rowTypeVec[PrevCpnFwdStartDate] = ARM_DATE_TYPE;

				ARM_Date dFwdEndDate((*cpnFwdStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]);
				dFwdEndDate.AddPeriod(itsCpnIdxTerm);
				CC_Ostringstream prevCpnFwdEndDateDesc;
				prevCpnFwdEndDateDesc << CC_NS(std,fixed) << dFwdEndDate.GetJulian();
				rowDescVec[PrevCpnFwdEndDate] = prevCpnFwdEndDateDesc.str();
				rowTypeVec[PrevCpnFwdEndDate] = ARM_DATE_TYPE;

				CC_Ostringstream prevConstDesc;
				prevConstDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsConstProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);

				rowDescVec[PrevConst] = prevConstDesc.str();
				rowTypeVec[PrevConst] = ARM_STRING;

				CC_Ostringstream prevLevPrevCpnDesc;
				prevLevPrevCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsLPrevCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);
				rowDescVec[PrevLevPrevCpn] = prevLevPrevCpnDesc.str();
				rowTypeVec[PrevLevPrevCpn] = ARM_DOUBLE;

				CC_Ostringstream prevNewOptDesc;
				prevNewOptDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsLNewOptProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);
				rowDescVec[PrevLevNewOpt] = prevNewOptDesc.str();
				rowTypeVec[PrevLevNewOpt] = ARM_DOUBLE;

				CC_Ostringstream prevStrikeOptDesc;
				prevStrikeOptDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsStrikeOptProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);
				rowDescVec[PrevStrikeOpt] = prevStrikeOptDesc.str();
				rowTypeVec[PrevStrikeOpt] = ARM_DOUBLE;

				CC_Ostringstream prevMinCpnDesc;
				prevMinCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsMinCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);
				rowDescVec[PrevMinCpn] = prevMinCpnDesc.str();
				rowTypeVec[PrevMinCpn] = ARM_DOUBLE;

				CC_Ostringstream prevMaxCpnDesc;
				prevMaxCpnDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsMaxCpnProfile).Interpolate((*cpnStartDates)[itsSBRank2Pos[itsExerPos2Rank[eventIdx]-1]]-asOfDate);
				rowDescVec[PrevMaxCpn] = prevMaxCpnDesc.str();
				rowTypeVec[PrevMaxCpn] = ARM_DOUBLE;

				CC_Ostringstream prevCpnIndexDesc;
				prevCpnIndexDesc << "LIBOR(" << cpnModelName << "," << CallableSnowBallColNamesTable[PrevCpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[PrevCpnFwdEndDate] << "[i]," << cpnIndexDayCount << ")";

				rowDescVec[PrevCpnIndex] = prevCpnIndexDesc.str();
				rowTypeVec[PrevCpnIndex] = ARM_STRING;

				CC_Ostringstream prevCpnSimpleDesc;
				prevCpnSimpleDesc << CallableSnowBallColNamesTable[PrevConst] << "[i]+";

				// Switch betwween call the SB and call the stacky
				if ((itsExerPos2Rank[eventIdx] > 1) && itsCallSBOrStacky)
				{
					int decalage = eventIdx - itsSBRank2Pos[itsExerPos2Rank[eventIdx] - 2];

					prevCpnSimpleDesc << CallableSnowBallColNamesTable[PrevLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[Coupon] << "[i-" << CC_NS(std,fixed) << decalage << "]+";
				}

				if (itsCapOrFloor == K_CAP)
				{
					prevCpnSimpleDesc << CallableSnowBallColNamesTable[PrevLevNewOpt] << "[i]*("<< CallableSnowBallColNamesTable[PrevCpnIndex] <<"[i]-" << CallableSnowBallColNamesTable[PrevStrikeOpt] << "[i])";
				}
				else
				{
					prevCpnSimpleDesc << CallableSnowBallColNamesTable[PrevLevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[PrevStrikeOpt] << "[i]-"<< CallableSnowBallColNamesTable[PrevCpnIndex] <<"[i])";
				}

				rowDescVec[PrevCouponSimple] = prevCpnSimpleDesc.str();
				rowTypeVec[PrevCouponSimple] = ARM_STRING;

				CC_Ostringstream prevCpnOptDesc;

				prevCpnOptDesc << "Max(" <<   CallableSnowBallColNamesTable[PrevMinCpn] << "[i]";
				prevCpnOptDesc << "-" << CallableSnowBallColNamesTable[PrevCouponSimple] << "[i],0)";
				prevCpnOptDesc << "+Max(" <<   CallableSnowBallColNamesTable[PrevCouponSimple] << "[i]-" << CallableSnowBallColNamesTable[PrevMaxCpn] << "[i],0)";

				rowDescVec[PrevCouponOpt] = prevCpnOptDesc.str();
				rowTypeVec[PrevCouponOpt] = ARM_STRING;

				CC_Ostringstream prevCpnDesc;

				// Switch betwween call the SB and call the stacky
				if ((itsExerPos2Rank[eventIdx] > 1) && !itsCallSBOrStacky)
				{
					int decalage = eventIdx - itsSBRank2Pos[itsExerPos2Rank[eventIdx] - 2];

					prevCpnDesc << CallableSnowBallColNamesTable[PrevLevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[Coupon] << "[i-" << CC_NS(std,fixed) << decalage << "]+";
				}

				prevCpnDesc << CallableSnowBallColNamesTable[PrevCouponSimple] << "[i]+" << CallableSnowBallColNamesTable[PrevCouponOpt] << "[i]";

				rowDescVec[PrevCoupon] = prevCpnDesc.str();
				rowTypeVec[PrevCoupon] = ARM_STRING;
			}
			else
			{
				int Decalage = eventIdx - itsSBRank2Pos[itsExerPos2Rank[eventIdx] - 1];

				CC_Ostringstream prevCpnDesc;
				prevCpnDesc << CC_NS(std,fixed) << CallableSnowBallColNamesTable[Coupon] << "[i-" << CC_NS(std,fixed) << Decalage << "]";
				rowDescVec[PrevCoupon] = prevCpnDesc.str();
				rowTypeVec[PrevCoupon] = ARM_STRING;
			}
		}
		else
		{
			CC_Ostringstream prevCpnDesc;
			prevCpnDesc << CC_NS(std,fixed) << 0;
			rowDescVec[PrevCoupon] = prevCpnDesc.str();
			rowTypeVec[PrevCoupon] = ARM_STRING;
		}

		CC_Ostringstream fundingCoeffDesc;
		fundingCoeffDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundCoeffProfile).Interpolate(fundingStartDate-asOfDate);
		rowDescVec[FundingCoeff] = fundingCoeffDesc.str();
		rowTypeVec[FundingCoeff] = ARM_DOUBLE;

		CC_Ostringstream fundingMarginDesc;
		fundingMarginDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundMarginProfile).Interpolate(fundingStartDate-asOfDate);
		rowDescVec[FundingMargin] = fundingMarginDesc.str();
		rowTypeVec[FundingMargin] = ARM_DOUBLE;

		CC_Ostringstream fundingDesc;
		fundingDesc << CallableSnowBallColNamesTable[FundingCoeff] << "[i]*Swap(" << payModelName << "," << CallableSnowBallColNamesTable[FundingStartDate] << "[i]," << CallableSnowBallColNamesTable[FundingEndDate] << "[i],0,P,,," << fundFreq << "," << fundDayCount << ",0)+";
		fundingDesc << CallableSnowBallColNamesTable[FundingMargin] << "[i]*Annuity(" << payModelName << "," << CallableSnowBallColNamesTable[FundingStartDate] << "[i]," << CallableSnowBallColNamesTable[FundingEndDate] << "[i]," << fundFreq << "," << fundDayCount+ ")";
		rowDescVec[Funding] = fundingDesc.str();
		rowTypeVec[Funding] = ARM_STRING;

		CC_Ostringstream fundingPaidDesc;
		fundingPaidDesc << CallableSnowBallColNamesTable[Funding] << "[i]*" << CallableSnowBallColNamesTable[NextNominal] << "[i]";
		
		
		rowDescVec[FundingPaid] = fundingPaidDesc.str();
		rowTypeVec[FundingPaid] = ARM_STRING;

		if (itsPayRec == K_PAY)
		{
			CC_Ostringstream CFDesc;
			CFDesc << "(" << CallableSnowBallColNamesTable[Funding] << "[i]-" << CallableSnowBallColNamesTable[NextCoupon] << "[i])*" << CallableSnowBallColNamesTable[NextNominal] << "[i]";
			rowDescVec[CF] = CFDesc.str();
			rowTypeVec[CF] = ARM_STRING;
		}
		else
		{
			CC_Ostringstream CFDesc;
			CFDesc << "(" << CallableSnowBallColNamesTable[NextCoupon] << "[i]-" << CallableSnowBallColNamesTable[Funding] << "[i])*" << CallableSnowBallColNamesTable[NextNominal] << "[i]";
			rowDescVec[CF] = CFDesc.str();
			rowTypeVec[CF] = ARM_STRING;
		}

		if (itsExerPos2Rank[eventIdx] <= itsNbNonCall - 1)
		{
			CC_Ostringstream bermudaDesc;
			bermudaDesc << CallableSnowBallColNamesTable[CF] << "[i]";
			rowDescVec[Bermuda] = bermudaDesc.str();
			rowTypeVec[Bermuda] = ARM_STRING;
		}
		else if (itsExerPos2Rank[eventIdx] == itsNbNonCall)
		{
			CC_Ostringstream bermudaDesc;
			bermudaDesc << CallableSnowBallColNamesTable[Option] << "[i]";
			rowDescVec[Bermuda] = bermudaDesc.str();
			rowTypeVec[Bermuda] = ARM_STRING;
		}
		else
		{
			CC_Ostringstream bermudaDesc;
			bermudaDesc << CC_NS(std,fixed) << 0.0;
			rowDescVec[Bermuda] = bermudaDesc.str();
			rowTypeVec[Bermuda] = ARM_DOUBLE;
		}

		if (itsExerPos2Rank[eventIdx] > itsNbNonCall - 1)
		{
			if (itsExerPos2Rank[eventIdx] < itsNbCoupons - 1)
			{
				int decalage = itsExerRank2Pos[itsExerPos2Rank[eventIdx] + 1] -eventIdx;

				CC_Ostringstream optionDesc;
				optionDesc << "Exercise(" << CallableSnowBallColNamesTable[CF] << "[i],-" << CallableSnowBallColNamesTable[CF] << "[i]-" << CallableSnowBallColNamesTable[Fees] << "[i]," << CallableSnowBallColNamesTable[Option] << "[i+" << CC_NS(std,fixed) << decalage << "]";

				if (itsRegressors.size() != 0)
					optionDesc << itsRegressors;

				optionDesc <<")";
				rowDescVec[Option] = optionDesc.str();
				rowTypeVec[Option] = ARM_STRING;
			}
			else
			{
				CC_Ostringstream optionDesc;
				optionDesc << CallableSnowBallColNamesTable[CF] << "[i]+Max(-" << CallableSnowBallColNamesTable[CF] << "[i]-" << CallableSnowBallColNamesTable[Fees] << "[i],0)";
				rowDescVec[Option] = optionDesc.str();
				rowTypeVec[Option] = ARM_STRING;
			}
		}
	}

	if (itsSBPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
    {
		CC_Ostringstream dfPayDesc;
		dfPayDesc << "DF(" << payModelName << "," << CallableSnowBallColNamesTable[CpnPayDate] << "[i])";
		rowDescVec[DFPay] = dfPayDesc.str();
		rowTypeVec[DFPay] = ARM_STRING;

		CC_Ostringstream nominalDesc;
		nominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNotionalProfile).Interpolate((*cpnPayDates)[eventIdx]-asOfDate);
		rowDescVec[Nominal] = nominalDesc.str();
		rowTypeVec[Nominal] = ARM_DOUBLE;

		CC_Ostringstream couponIndexDesc;
		couponIndexDesc << "LIBOR(" << cpnModelName << "," << CallableSnowBallColNamesTable[CpnFwdStartDate] << "[i]," << CallableSnowBallColNamesTable[CpnFwdEndDate] << "[i]," << cpnIndexDayCount << ")";
		rowDescVec[CpnIndex] = couponIndexDesc.str();
		rowTypeVec[CpnIndex] = ARM_STRING;

		CC_Ostringstream couponIncDesc;

		if (itsCapOrFloor == K_CAP)
		{
			couponIncDesc << CallableSnowBallColNamesTable[LevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[CpnIndex] << "[i]-" << CallableSnowBallColNamesTable[StrikeOpt] << "[i])";
		}
		else if (itsCapOrFloor == K_FLOOR)
		{
			couponIncDesc << CallableSnowBallColNamesTable[LevNewOpt] << "[i]*(" << CallableSnowBallColNamesTable[StrikeOpt] << "[i]-" << CallableSnowBallColNamesTable[CpnIndex] << "[i])";
		}

		if (itsSBPos2Rank[eventIdx] > 0)
		{
			int decalage = eventIdx - itsSBRank2Pos[itsSBPos2Rank[eventIdx] - 1];
			
			CC_Ostringstream couponSimpleDesc;

			if (itsCallSBOrStacky)
			{
				couponSimpleDesc << CallableSnowBallColNamesTable[Const] << "[i]+" << CallableSnowBallColNamesTable[LevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[Coupon] 
					<< "[i-" << CC_NS(std,fixed) << decalage << "]+";
			}

			couponSimpleDesc << couponIncDesc.str();
			
			rowDescVec[CouponSimple] = couponSimpleDesc.str();
			rowTypeVec[CouponSimple] = ARM_STRING;

			CC_Ostringstream majoCouponDesc;
			majoCouponDesc << CallableSnowBallColNamesTable[MinCpn]<<"[i]+"<<CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]*("<<CallableSnowBallColNamesTable[UBCoupon]<<"[i-"<<CC_NS(std,fixed) << decalage <<"]-"<<CallableSnowBallColNamesTable[MinCpn]<<"[i-" << CC_NS(std,fixed) << decalage << "])+"
				<<"Min(Max("<< CallableSnowBallColNamesTable[Const] <<"[i]+" <<couponIncDesc.str() <<"+" << CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]*"<<CallableSnowBallColNamesTable[MinCpn] << "[i-"<< CC_NS(std,fixed) << decalage <<"]-"<<CallableSnowBallColNamesTable[MinCpn] << "[i], 0)," << CallableSnowBallColNamesTable[MaxCpn] << "[i])";
			rowDescVec[UBCoupon] = majoCouponDesc.str();
			rowTypeVec[UBCoupon] = ARM_STRING;

			CC_Ostringstream minoCouponDesc;
			minoCouponDesc << CallableSnowBallColNamesTable[LevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[LBCoupon] << "[i-" << CC_NS(std,fixed) << decalage << "]+" <<  CallableSnowBallColNamesTable[Const] << "[i]+" << couponIncDesc.str();
			rowDescVec[LBCoupon] = minoCouponDesc.str();
			rowTypeVec[LBCoupon] = ARM_STRING;

		}
		else
		{
			CC_Ostringstream couponSimpleDesc;
			couponSimpleDesc << CallableSnowBallColNamesTable[Const] << "[i]+" << couponIncDesc.str();

			rowDescVec[CouponSimple] = couponSimpleDesc.str();
			rowTypeVec[CouponSimple] = ARM_STRING;

			CC_Ostringstream majoCouponDesc;
			majoCouponDesc<< "Min(Max("<< CallableSnowBallColNamesTable[Const] <<"[i]+" <<couponIncDesc.str() <<"-"<<CallableSnowBallColNamesTable[MinCpn] << "[i], 0)," << CallableSnowBallColNamesTable[MaxCpn] << "[i])";
			rowDescVec[UBCoupon] = majoCouponDesc.str();
			rowTypeVec[UBCoupon] = ARM_STRING;

			CC_Ostringstream minoCouponDesc;
			minoCouponDesc <<  CallableSnowBallColNamesTable[Const] << "[i]+" << couponIncDesc.str();
			rowDescVec[LBCoupon] = minoCouponDesc.str();
			rowTypeVec[LBCoupon] = ARM_STRING;
		}

		CC_Ostringstream couponOptionDesc;

		if ((itsSBPos2Rank[eventIdx] > 0) && !itsCallSBOrStacky)
		{
			couponOptionDesc << CallableSnowBallColNamesTable[Const] << "[i]+";
			if (itsSBPos2Rank[eventIdx] > 0)
			{
				int decalage = eventIdx - itsSBRank2Pos[itsSBPos2Rank[eventIdx] - 1];
				couponOptionDesc << CallableSnowBallColNamesTable[LevPrevCpn] << "[i]*" << CallableSnowBallColNamesTable[Coupon];
				couponOptionDesc << "[i-" << CC_NS(std,fixed) << decalage << "]+";
			}
		}

		couponOptionDesc << "Max(" <<   CallableSnowBallColNamesTable[MinCpn] << "[i]";
		couponOptionDesc << "-" << CallableSnowBallColNamesTable[CouponSimple] << "[i],0)";
		couponOptionDesc << "-Max(" <<   CallableSnowBallColNamesTable[CouponSimple] << "[i]-" << CallableSnowBallColNamesTable[MaxCpn] << "[i],0)";

		rowDescVec[CouponOption] = couponOptionDesc.str();
		rowTypeVec[CouponOption] = ARM_STRING;

		CC_Ostringstream couponDesc;
		couponDesc << CallableSnowBallColNamesTable[CouponSimple] << "[i]+" << CallableSnowBallColNamesTable[CouponOption] << "[i]";

		rowDescVec[Coupon] = couponDesc.str();
		rowTypeVec[Coupon] = ARM_STRING;



		const_cast< std::vector<double>& >(itsStrikeSum)[itsSBPos2Rank[eventIdx]] = itsStrikeOptVec[itsSBPos2Rank[eventIdx]]*itsLevNewVec[itsSBPos2Rank[eventIdx]]-itsCapOrFloor*itsConstVec[itsSBPos2Rank[eventIdx]];
		if ((itsSBPos2Rank[eventIdx]>0))
		{
			const_cast< std::vector<double>& >(itsStrikeSum)[itsSBPos2Rank[eventIdx]]+=itsStrikeSum[itsSBPos2Rank[eventIdx]-1]*itsLevPrevVec[itsSBPos2Rank[eventIdx]];
		}
			

		CC_Ostringstream strikeSumDesc;
		strikeSumDesc << itsStrikeSum[itsSBPos2Rank[eventIdx]];

		rowDescVec[StrikeSum] = strikeSumDesc.str();
		rowTypeVec[StrikeSum] = ARM_DOUBLE;


		//LBCouponAvg
		CC_Ostringstream minoCouponAvgDesc;
		minoCouponAvgDesc<< "Max(" << CallableSnowBallColNamesTable[LBCoupon] << "[i],"<<CallableSnowBallColNamesTable[MinCpn]<<"[i])";
		rowDescVec[LBCouponAvg] = minoCouponAvgDesc.str();
		rowTypeVec[LBCouponAvg] = ARM_STRING;

		//Ana LB Coupon Avg
		CC_Ostringstream anaLBCouponAvgDesc;
		if (itsPayRec == K_PAY)
		{
			anaLBCouponAvgDesc << "-";
		}
		else
		{
			anaLBCouponAvgDesc << "+";
		}

		if ((itsSBPos2Rank[eventIdx] == 0) || (fabs(itsLevPrevVec[itsSBPos2Rank[eventIdx]]) < K_DOUBLE_TOL))
		{
			const_cast<int&>(itsLastReset) = itsSBPos2Rank[eventIdx];
		}

		anaLBCouponAvgDesc <<  "(" << CallableSnowBallColNamesTable[MinCpn] << "[i]*DFPay[i]+"; 
		anaLBCouponAvgDesc <<  "SumOption(" << payModelName << "," << CallableSnowBallColNamesTable[SumOptStartDate] << "[i],"; 
		anaLBCouponAvgDesc <<  CallableSnowBallColNamesTable[EndDate] << "[i]," << CallableSnowBallColNamesTable[CpnPayDate] << "[i],";
		anaLBCouponAvgDesc << CallableSnowBallColNamesTable[StrikeSum] << "[i]";
		if (itsCapOrFloor)
			anaLBCouponAvgDesc << "-";
		else
			anaLBCouponAvgDesc << "+";
		anaLBCouponAvgDesc << CallableSnowBallColNamesTable[MinCpn] << "[i]";
		anaLBCouponAvgDesc << "," << capFloor << "," << cpnFreq << "," << itsCpnResetLag << ","; 
		anaLBCouponAvgDesc << cpnIndexTiming << "," << itsCpnIdxTerm << "," << cpnIndexDayCount <<  ",Coeff" << itsSBPos2Rank[eventIdx] << "," << itsLastReset << ")";
		anaLBCouponAvgDesc << ")*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		rowDescVec[AnaLBCouponAvg] = anaLBCouponAvgDesc.str();
		rowTypeVec[AnaLBCouponAvg] = ARM_STRING;


		//Ana LB Coupon
		CC_Ostringstream anaLBCouponDesc;
		int sign = (itsPayRec == K_RCV?1:-1);
		sign *= itsCapOrFloor;

		if (sign < 0)
		{
			anaLBCouponDesc << "-";
		}

		anaLBCouponDesc <<  "(SumOption(" << payModelName << "," << CallableSnowBallColNamesTable[SumOptStartDate] << "[i],";
		anaLBCouponDesc <<  CallableSnowBallColNamesTable[EndDate] << "[i]," << CallableSnowBallColNamesTable[CpnPayDate] << "[i],";
		anaLBCouponDesc << "0,C," << cpnFreq << "," << itsCpnResetLag << ","; 
		anaLBCouponDesc << cpnIndexTiming << "," << itsCpnIdxTerm << "," << cpnIndexDayCount <<  ",Coeff" << itsSBPos2Rank[eventIdx] << "," << itsLastReset << ")";
		anaLBCouponDesc << "-" << CallableSnowBallColNamesTable[StrikeSum] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i])";
		anaLBCouponDesc << "*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		rowDescVec[AnaLBCoupon] = anaLBCouponDesc.str();
		rowTypeVec[AnaLBCoupon] = ARM_STRING;


		//CouponOpt and CouponOptAvg
		if (itsSBPos2Rank[eventIdx] > 0)
		{
			int decalage = eventIdx - itsSBRank2Pos[itsSBPos2Rank[eventIdx] - 1];

			CC_Ostringstream couponOptDesc;	
			couponOptDesc << "Max(-"<< couponIncDesc.str() <<"-"<<CallableSnowBallColNamesTable[Const]<<"[i]+"<<CallableSnowBallColNamesTable[MinCpn]<<"[i]-"<<CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]*"<<CallableSnowBallColNamesTable[Coupon]<<"[i-"<< CC_NS(std,fixed) << decalage <<"], 0)"
						  << "+" << CallableSnowBallColNamesTable[CouponOpt] << "[i-"<< CC_NS(std,fixed) << decalage <<"]*"<<CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]";
			rowDescVec[CouponOpt] = couponOptDesc.str();
			rowTypeVec[CouponOpt] = ARM_STRING;

			CC_Ostringstream couponOptAvgDesc;
			couponOptAvgDesc << "Max(-"<< couponIncDesc.str() <<"-"<<CallableSnowBallColNamesTable[Const]<<"[i]+"<<CallableSnowBallColNamesTable[MinCpn]<<"[i]-"<<CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]*"<<CallableSnowBallColNamesTable[LBCoupon]<<"[i-"<< CC_NS(std,fixed) << decalage <<"], 0)"
			  << "+" << CallableSnowBallColNamesTable[CouponOptAvg] << "[i-"<< CC_NS(std,fixed) << decalage <<"]*"<<CallableSnowBallColNamesTable[LevPrevCpn]<<"[i]";
			rowDescVec[CouponOptAvg] = couponOptAvgDesc.str();
			rowTypeVec[CouponOptAvg] = ARM_STRING;

		}
		else
		{		
			CC_Ostringstream couponOptDesc;	
			couponOptDesc << "Max(-"<< couponIncDesc.str() <<"-"<<CallableSnowBallColNamesTable[Const]<<"[i]+"<<CallableSnowBallColNamesTable[MinCpn]<<"[i], 0)";
			rowDescVec[CouponOpt] = couponOptDesc.str();
			rowTypeVec[CouponOpt] = ARM_STRING;

			CC_Ostringstream couponOptAvgDesc;
			couponOptAvgDesc<<  "Max(-"<< couponIncDesc.str() <<"-"<<CallableSnowBallColNamesTable[Const]<<"[i]+"<<CallableSnowBallColNamesTable[MinCpn]<<"[i],0)";
			rowDescVec[CouponOptAvg] = couponOptAvgDesc.str();
			rowTypeVec[CouponOptAvg] = ARM_STRING;
		}

		//CouponOptPaid
		CC_Ostringstream couponOptPaidDesc;

		if (itsPayRec == K_PAY)
		{
			couponOptPaidDesc << "-" << CallableSnowBallColNamesTable[CouponOpt] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		else
		{
			couponOptPaidDesc << CallableSnowBallColNamesTable[CouponOpt] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}

		rowDescVec[CouponOptPaid] = couponOptPaidDesc.str();
		rowTypeVec[CouponOptPaid] = ARM_STRING;


		CC_Ostringstream LBCouponPaidDesc;

		if (itsPayRec == K_PAY)
		{
			LBCouponPaidDesc << "-" << CallableSnowBallColNamesTable[LBCoupon] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		else
		{
			LBCouponPaidDesc << CallableSnowBallColNamesTable[LBCoupon] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}

		rowDescVec[LBCouponPaid] = LBCouponPaidDesc.str();
		rowTypeVec[LBCouponPaid] = ARM_STRING;

		CC_Ostringstream UBCouponPaidDesc;

		if (itsPayRec == K_PAY)
		{
			UBCouponPaidDesc << "-" << CallableSnowBallColNamesTable[UBCoupon] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		else
		{
			UBCouponPaidDesc << CallableSnowBallColNamesTable[UBCoupon] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		
		rowDescVec[UBCouponPaid] = UBCouponPaidDesc.str();
		rowTypeVec[UBCouponPaid] = ARM_STRING;

		CC_Ostringstream couponPaidDesc;
		couponPaidDesc << CallableSnowBallColNamesTable[Coupon] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]";
		
		rowDescVec[PaidCoupon] = couponPaidDesc.str();
		rowTypeVec[PaidCoupon] = ARM_STRING;

		//CouponOptAvgPaid
		CC_Ostringstream couponOptAvgPaidDesc;
		if (itsPayRec == K_PAY)
		{
			couponOptAvgPaidDesc << "-" << CallableSnowBallColNamesTable[CouponOptAvg] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		else
		{
			couponOptAvgPaidDesc << CallableSnowBallColNamesTable[CouponOptAvg] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
				
		rowDescVec[CouponOptAvgPaid] = couponOptAvgPaidDesc.str();
		rowTypeVec[CouponOptAvgPaid] = ARM_STRING;

		//LBCouponAvgPaid
		CC_Ostringstream LBCouponAvgPaidDesc;

		if (itsPayRec == K_PAY)
		{
			LBCouponAvgPaidDesc << "-" << CallableSnowBallColNamesTable[LBCouponAvg] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		else
		{
			LBCouponAvgPaidDesc << CallableSnowBallColNamesTable[LBCouponAvg] << "[i]*" << CallableSnowBallColNamesTable[DFPay] << "[i]*" << CallableSnowBallColNamesTable[CpnIT] << "[i]*" << CallableSnowBallColNamesTable[Nominal] << "[i]";
		}
		
		rowDescVec[LBCouponAvgPaid] = LBCouponAvgPaidDesc.str();
		rowTypeVec[LBCouponAvgPaid] = ARM_STRING;
	}

	MiddleRowsExt(eventIdx, datesStructure, rowDescVec, rowTypeVec);

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: MiddleRowsExt
///	Returns: void
///	Action : create a row of a deal description
/// This function seems to be use less but I have some compiler crashs 
/// and its the only way I found to fix it.
/// I think the MiddleRows function is really too long
/////////////////////////////////////////////////////////////////

void ARM_CallableSnowBallCalculator::MiddleRowsExt(int eventIdx, const ARM_DateStripCombiner& datesStructure,  vector< string >& rowDescVec, vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
	size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();

	CC_Ostringstream CF1Desc;
	CC_Ostringstream CF2Desc;
	CC_Ostringstream Option2Desc;
	CC_Ostringstream AMCIndexDesc;

	int rank = 0;

	if ( (itsExerPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData) &&
		 (itsSBPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData))
    {
		rank = itsExerPos2Rank[eventIdx];
		if (itsCpnIndexTiming == K_ADVANCE)
		{
			if (itsPayRec == K_PAY)
			{
				CF1Desc << "(" << CallableSnowBallColNamesTable[FundingPaid] << "[i]-";
				CF1Desc << CallableSnowBallColNamesTable[PaidCoupon] << "[i]" << ")";

				CF2Desc << "(" << CallableSnowBallColNamesTable[PaidCoupon] << "[i]-";
				CF2Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]" << ")";
			}
			else
			{
				CF1Desc << "(" << CallableSnowBallColNamesTable[PaidCoupon] << "[i]-";
				CF1Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]" << ")";

				CF2Desc << "(" << CallableSnowBallColNamesTable[FundingPaid] << "[i]-";
				CF2Desc << CallableSnowBallColNamesTable[PaidCoupon] << "[i]" << ")";
			}
		}
		else
		{
			if (itsPayRec == K_PAY)
			{
				CF1Desc << "(" << CallableSnowBallColNamesTable[FundingPaid] << "[i]-";
				CF1Desc << CallableSnowBallColNamesTable[PaidCoupon] << "[i]" << ")";

				CF2Desc << "-" << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
			}
			else
			{
				CF1Desc << "(" << CallableSnowBallColNamesTable[PaidCoupon] << "[i]-";
				CF1Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]" << ")";

				CF2Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
			}	
		}

		AMCIndexDesc << CallableSnowBallColNamesTable[Funding] << "[i]-";
		AMCIndexDesc << "(" << CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i]+";
		AMCIndexDesc << "(Max(" <<   CallableSnowBallColNamesTable[NextMinCpn] << "[i]*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
		AMCIndexDesc << "-" << CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i],0)";
		AMCIndexDesc << "+Max(" <<   CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i]";
		AMCIndexDesc << "-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]*" << CallableSnowBallColNamesTable[DFNextPay] << "[i],0)))*";
		AMCIndexDesc << CallableSnowBallColNamesTable[NextCpnIT] << "[i]";
		if ( eventIdx < eventSize-1 )
		{
			Option2Desc << "Exercise(" << CallableSnowBallColNamesTable[CF1] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[CF2] << "[i]";
			Option2Desc << "-" << CallableSnowBallColNamesTable[Fees] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[Option2] << "[i+1],";
			if (itsRegressors.size() == 0)
			{
				Option2Desc << CallableSnowBallColNamesTable[AMCIndex] << "[i],1,";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],2),";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],3),";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],4)";
			}
			else
			{
				Option2Desc << itsRegressors;
			}

			Option2Desc << ")";
		}
		else
		{
			Option2Desc << "Max(" << CallableSnowBallColNamesTable[CF1] << "[i],0)";
		}
	}
	else if (itsExerPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
	{
		rank = itsExerPos2Rank[eventIdx];
		if (itsPayRec == K_PAY)
		{
			CF1Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
			CF2Desc << "-" << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
		}
		else
		{
			CF1Desc << "-" << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
			CF2Desc << CallableSnowBallColNamesTable[FundingPaid] << "[i]";
		}

		AMCIndexDesc << "(" << CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i]+";
		AMCIndexDesc << "(Max(" <<   CallableSnowBallColNamesTable[NextMinCpn] << "[i]*" << CallableSnowBallColNamesTable[DFNextPay] << "[i]";
		AMCIndexDesc << "-" << CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i],0)";
		AMCIndexDesc << "+Max(" <<   CallableSnowBallColNamesTable[NextSimpleCouponWithoutAdj] << "[i]";
		AMCIndexDesc << "-" << CallableSnowBallColNamesTable[NextMaxCpn] << "[i]*" << CallableSnowBallColNamesTable[DFNextPay] << "[i],0)))*";
		AMCIndexDesc << CallableSnowBallColNamesTable[NextCpnIT] << "[i]";

		if ( eventIdx < eventSize-1 )
		{
			Option2Desc << "Exercise(" << CallableSnowBallColNamesTable[CF1] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[CF2] << "[i]";
			Option2Desc << "-" << CallableSnowBallColNamesTable[Fees] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[Option2] << "[i+1],";
			if (itsRegressors.size() == 0)
			{
				Option2Desc << CallableSnowBallColNamesTable[AMCIndex] << "[i],1,";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],2),";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],3),";
				Option2Desc << "POW(" << CallableSnowBallColNamesTable[AMCIndex] << "[i],4)";
			}
			else
			{
				Option2Desc << itsRegressors;
			}

			Option2Desc << ")";
		}
		else
		{
			Option2Desc << "Max(" << CallableSnowBallColNamesTable[CF1] << "[i],0)";
		}
	}
	else if (itsSBPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData)
	{
		rank = itsSBPos2Rank[eventIdx];
		if (itsPayRec == K_PAY)
			CF1Desc << "-" << CallableSnowBallColNamesTable[PaidCoupon] << "[i]";
		else
			CF1Desc << CallableSnowBallColNamesTable[PaidCoupon] << "[i]";
			
		CF2Desc << "-10000000";

		AMCIndexDesc << CallableSnowBallColNamesTable[CpnIndex] << "[i]";

		if ( eventIdx < eventSize-1 )
		{
			Option2Desc << "Exercise(" << CallableSnowBallColNamesTable[CF1] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[CF2] << "[i],";
			Option2Desc << CallableSnowBallColNamesTable[Option2] << "[i+1],";
			Option2Desc << CallableSnowBallColNamesTable[AMCIndex] << "[i])";
		}
		else
		{
			Option2Desc << CallableSnowBallColNamesTable[CF1] << "[i]";
		}
	}

	rowDescVec[CF1] = CF1Desc.str();
	rowTypeVec[CF1] = ARM_STRING;

	rowDescVec[CF2] = CF2Desc.str();
	rowTypeVec[CF2] = ARM_STRING;

	rowDescVec[AMCIndex] = AMCIndexDesc.str();
	rowTypeVec[AMCIndex] = ARM_STRING;

	rowDescVec[Option2] = Option2Desc.str();
	rowTypeVec[Option2] = ARM_STRING;

	CC_Ostringstream Bermuda2Desc;

	if (rank <= itsNbNonCall - 1)
		Bermuda2Desc << CallableSnowBallColNamesTable[CF1] << "[i]";
	else if ((rank == itsNbNonCall) && (itsExerPos2Rank[eventIdx] != ARM_DateStripCombiner::DateStripCombiner_BlankData))
		Bermuda2Desc << CallableSnowBallColNamesTable[Option2] << "[i]";
	else
		Bermuda2Desc << CC_NS(std,fixed) << 0.0;

	rowDescVec[Bermuda2] = Bermuda2Desc.str();
	rowTypeVec[Bermuda2] = ARM_STRING;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CallableSnowBallCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(CallableSnowBallColNamesTable)/sizeof(CallableSnowBallColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = CallableSnowBallColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: Clone, View
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_Object* ARM_CallableSnowBallCalculator::Clone() const
{
	return new ARM_CallableSnowBallCalculator(*this);
}


void ARM_CallableSnowBallCalculator::View(char* id, FILE* ficOut) const
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

    /// CSB Calculator specific datas viewing
    fprintf(fOut,"\n\n =======> CALLABLE SNOWBALL CALCULATOR <====== \n");

    CC_Ostringstream callableSbData;
    callableSbData << "\nStartDate = " <<  itsStartDate.toString() << "\n";
    callableSbData << "EndDate = " << itsEndDate.toString() << "\n";
    callableSbData << "Pay/Rec = " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << "\n";


    callableSbData << "\nCoupon Datas :\n";
    callableSbData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsCpnDayCount ) << "\n";
    callableSbData << "Frequency (Reset & Pay) = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    callableSbData << "Pay Calendar = " << itsCpnPayCal << "\n";
    callableSbData << "Reset Timing = " << ARM_ParamView::GetMappingName(S_TIMING_MOD, itsCpnIndexTiming ) << "\n";
    callableSbData << "Reset Gap = " << itsCpnResetLag << "\n";
    callableSbData << "Reset Calendar = " << itsCpnResetCal << "\n";
    callableSbData << "Cpn Index Term = " << itsCpnIdxTerm << "\n";
	callableSbData << "Cpn Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsCpnFreq ) << "\n";
    callableSbData << "Cpn Index Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT, itsCpnIndexDayCount )  << "\n";

    callableSbData << "\nFunding Leg Datas :\n";
    callableSbData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFundDayCount ) << "\n";
    callableSbData << "Frequency = " << ARM_ParamView::GetMappingName(S_FREQUENCY, itsFundFreq ) << "\n\n" ;

    /// part common to gencalculator
	//callableSbData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
	fprintf(fOut,"%s",callableSbData.str().c_str());

	CC_Ostringstream sumOptCoeffData;

	size_t i, j;

	sumOptCoeffData << "Sum Option Coeffs: " << endl;

	for (i = 0; i < itsSumOptionCoeffs.GetRowsNb(); ++i)
	{
		for (j = 0; j < itsSumOptionCoeffs.GetColsNb(); ++j)
		{
			sumOptCoeffData << setw(5);
			sumOptCoeffData << itsSumOptionCoeffs(i,j) << "\t";
		}

		sumOptCoeffData << endl;
	}

	fprintf(fOut,"%s",sumOptCoeffData.str().c_str());

	string capFloorPFStr;

	ARM_CalibMethod* volCalibMethod = GetVolCalibMethod();
	ARM_StdPortfolioPtr volPF = volCalibMethod->GetPortfolio();

	if (itsCalibMode == SUMOPT)
	{
		CC_Ostringstream os;

		os << "Volatility Portfolio Calibration : " << endl;

		os << std::setiosflags(std::ios::left);
		os << std::setiosflags(std::ios::showpoint);

		os << setw(5);
		os << "Nb\t";
		os << setw(11);
		os << "AssetName\t";
		os << setw(11);
		os << "StartDates\t";
		os << setw(11);
		os << "EndDates\t";
		os << setw(11);
		os << "PayDates\t";
		os << setw(11);
		os << "Strike(%)\t";
		os << setw(11);
		os << "SumFwd(%)\t";
		os << setw(11);
		os << "SumVol(%)\t";
		os << setw(11);
		os << "MktPrecision\t";
		os << setw(11);
		os << "MktPrice(%)\t";
		os << setw(11);
		os << "Weight" << endl << endl;

		int size = volPF->GetSize();

		ARM_Vector* precisions = volPF->GetPrecision();
		ARM_Vector*  weights = volPF->GetWeights();
		ARM_Vector* mktPrices =  volPF->GetMktPrices();

		for(size_t i = 0; i < size; ++i)
		{
			ARM_SumOpt* sumOpt = dynamic_cast<ARM_SumOpt*>(volPF->GetAsset(i));
			if (sumOpt)
			{
				ARM_Date startDate(sumOpt->GetStartDate());
				ARM_Date endDate(sumOpt->GetEndDate());
				ARM_Date payDate(sumOpt->GetPayDate());

				os << setw(5) << setfill(' ');
				os << i <<"\t";
				os << setw(11);
				os << "SUMOPTION\t";
				os << setw(11);
				os << startDate.toString() << "\t";
				os << setw(11);
				os << endDate.toString() << "\t";
				os << setw(11);
				os << payDate.toString() << "\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << sumOpt->GetStrike()*100 << "\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << sumOpt->GetSumFwd()*100 << "\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << sumOpt->GetSumVol()*100 <<"\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << (*precisions)[i] << "\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << (*mktPrices)[i]*100 << "\t";
				os << setfill('0') << setw(11) << setprecision(8);
				os << (*weights)[i] << endl;
			}
		}

		os << endl;

		fprintf(fOut, "%s", os.str().c_str());
	}
	else
	{
		if (!volPF.IsNull())
			volPF->View(id,fOut);
	}

	if ((itsCalibMode == SUMOPT) && itsBetaCalib)
	{
		ARM_CalibMethod* betaCalibMethod = GetBetaCalibMethod();
		ARM_StdPortfolioPtr betaPF = betaCalibMethod->GetPortfolio();
		ARM_CapFloor* capFloor = static_cast<ARM_CapFloor*>(betaPF->GetAsset(0));
		ARM_SwapLeg* swapLeg = capFloor->GetSwapLeg();

		ARM_BSModel* BSModel = dynamic_cast< ARM_BSModel* >(GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );

		capFloor->SetModel(BSModel);
		capFloor->ComputePrice();

		ARM_Vector* flows = capFloor->GetCashFlowValues();

		ARM_Vector* startDates = swapLeg->GetFwdRateStartDates();
		ARM_Vector* endDates = swapLeg->GetFwdRateEndDates();
		ARM_Vector* payDates = swapLeg->GetPaymentDates();
		ARM_Vector* resetDates = swapLeg->GetResetDates();

		ARM_ReferenceValue* strikes = capFloor->GetStrikes();
		ARM_ReferenceValue* nominals = capFloor->GetAmount();

		double strike, nominal;

		CC_Ostringstream os;

		os << "Beta Portfolio Calibration : " << endl;

		os << std::setiosflags(std::ios::left);
		os << std::setiosflags(std::ios::showpoint);

		os << setw(5);
		os << "Nb\t";
		os << "ResetDates\t";
		os << setw(11);
		os << "StartDates\t";
		os << setw(11);
		os << "EndDates\t";
		os << setw(11);
		os << "PayDates\t";
		os << setw(11);
		os << "Strike(%)\t";
		os << setw(11);
		os << "Nominal(%)\t";
		os << setw(11);
		os << "Price(%)" << endl << endl;

		int size = startDates->GetSize();

		for(size_t i = 0; i < size; ++i)
		{
			ARM_Date resetDate((*resetDates)[i]);
			ARM_Date startDate((*startDates)[i]);
			ARM_Date endDate((*endDates)[i]);
			ARM_Date payDate((*payDates)[i]);

			strike = strikes->Interpolate((*resetDates)[i]);
			nominal = nominals->Interpolate((*payDates)[i]);


			os << setw(5) << setfill(' ');
			os << i <<"\t";
			os << setw(11);
			os << resetDate.toString() << "\t";
			os << setw(11);
			os << startDate.toString() << "\t";
			os << setw(11);
			os << endDate.toString() << "\t";
			os << setw(11);
			os << payDate.toString() << "\t";
			os << setfill('0') << setw(11) << setprecision(8);
			os << strike << "\t";
			os << setfill('0') << setw(11) << setprecision(8);
			os << nominal << "\t";
			os << setfill('0') << setw(11) << setprecision(8);
			os << (*flows)[i] << endl;
		}

		os << endl;

		fprintf(fOut, "%s", os.str().c_str());
	}

    /// Common viewing
    ARM_GenCalculator::View(id,fOut);

    if ( ficOut == NULL )
       fclose(fOut);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateSumOptionModel
///	Returns: 
///	Action :  Create an analytic lognormal (Beta=100%,MRS=0) SFRM 
///			  model to valuate sum option
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::CreateSumOptionModel()
{
	// Create a simple 1F SFRM to valuate sum option
	std::vector<double> defaultTimes;
	std::vector<double> defaultSigmas;
	
	ARM_CurveModelParam volParam( ARM_ModelParamType::Volatility, &defaultSigmas, &defaultTimes );
	ARM_CurveModelParam mrsParam( ARM_ModelParamType::MeanReversion, 0 );
	ARM_CurveModelParam betaParam( ARM_ModelParamType::Beta, 1 );

	ARM_ModelParamVector paramVector(3);
	paramVector[0] = &volParam;
	paramVector[1] = &mrsParam;
	paramVector[2] = &betaParam;

	ARM_INDEX_TYPE liborType = GetIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams =
			ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,1,K_DIAG);

	itsSumOptModel = ARM_PricingModelPtr(new ARM_SFRM( GetPricingModel()->GetZeroCurve(), *pSFRMModelParams));

	delete pSFRMModelParams;

	// This analytic model is calibrated on the CAP volatility for strike of the sum option.
	itsSumOptCalibMethod = ARM_CalibMethodPtr(CreateEmptyCalibration(CAP,false));
}

void ARM_CallableSnowBallCalculator::UpdateSumOptionStrike(double Strike)
{
	ARM_StdPortfolioPtr PF = itsSumOptCalibMethod->GetPortfolio();

	size_t nbProducts = PF->GetSize();

	for (size_t i = 0; i < nbProducts; ++i)
	{
		ARM_CapFloor* capFloor=static_cast< ARM_CapFloor* >(PF->GetAsset(i));

		if (Strike > K_DOUBLE_TOL)
			capFloor->SetStrike(Strike*100);
		else
			capFloor->SetStrike(-1.0);
	}
}

double ARM_CallableSnowBallCalculator::ComputeSumOptionPrice(ARM_SumOpt* sumOpt)
{
	ComputePortPrices(true,CAP,&*itsSumOptCalibMethod);
	itsSumOptCalibMethod->Calibrate(&*itsSumOptModel);

	return ARM_VanillaPricer::Price(sumOpt,&*itsSumOptModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: ComputeSumOptCoeffs
///	Returns: 
///	Action : Compute the sum option coeffs, 
/// When the previous coeff is different to 0 or 1, we need to take
/// its effect in the sumoption coeff
/////////////////////////////////////////////////////////////////
void ARM_CallableSnowBallCalculator::ComputeSumOptCoeffs()
{
	size_t nbFlows = itsLevNewVec.size();

	itsSumOptionCoeffs.resize(nbFlows,nbFlows);

	size_t i, j, firstIdx = 0;

	for (i = 0; i < nbFlows; ++i)
	{
		if (fabs(itsLevPrevVec[i]) < K_DOUBLE_TOL)
		{
			firstIdx = i;
		}

		if ((i == 0) || (i == firstIdx))
		{
			itsSumOptionCoeffs(i,i) = itsLevNewVec[i];
		}
		else
		{
			for (j = 0; j < nbFlows; ++j)
			{
				if ((j >= firstIdx) && (j < i))
					itsSumOptionCoeffs(i,j) = itsLevPrevVec[i]*itsSumOptionCoeffs(i-1,j);
				else if (j == i)
					itsSumOptionCoeffs(i,j) = itsLevNewVec[i];
				else
					itsSumOptionCoeffs(i,j) = 0.0;

			}
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateCstManager
///	Returns: 
///	Action :  Create the constant manager linked to the genere
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_CallableSnowBallCalculator::CreateCstManager()
{
	size_t nbFlows = itsLevNewVec.size();

	size_t i;

	vector<string> names;
	vector<ARM_GramFctorArg> values;


	double coeffPrev = 1.0;

	for (i = 0; i < nbFlows; ++i)
	{
		ostringstream os;

		os << "Coeff" << i;
		names.push_back(os.str());
		values.push_back(ARM_GramFctorArg(ARM_VectorPtr(static_cast<ARM_GP_Vector*>(itsSumOptionCoeffs.GetRow(i)->Clone()))));
	}

	return ARM_CstManagerPtr(new ARM_CstManager(names, values));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: CreateColumnNames
///	Returns: 
///	Action :  Create columnNames
/////////////////////////////////////////////////////////////////
ARM_StringVector ARM_CallableSnowBallCalculator::CreateColumnNames()
{
	ARM_StringVector columnNames;
	columnNames.reserve(NbProductsToPrice);
	size_t i;

    for (i = 0; i < NbProductsToPrice; ++i)
	{
		if (itsProductsToPrice[i])
		{
			columnNames.push_back( CallableSnowBallColNamesTable[ ARM_CallableSnowBallCalculator::ProductToPriceColumns()[i] ] );
		}
	}

	return columnNames;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CallableSnowBallCalculator
///	Routine: ProductToPriceColumns
///	Returns: 
///	Action :  Return the columns to be priced
/// Based on the trigger mode
/////////////////////////////////////////////////////////////////
const int* ARM_CallableSnowBallCalculator::ProductToPriceColumns() const
{
	if (itsTriggerMode == TCoupon)
	{
		return Product1ToPriceColumns;
	}
	else if (itsTriggerMode == TApproxCoupon)
	{
		return Product2ToPriceColumns;
	}
	else
		ARM_THROW( ERR_INVALID_ARGUMENT," CSBCalculator : wrong trigger mode !" );
}

void ARM_CallableSnowBallCalculator::FixBetaModelParam()
{
	const ARM_ModelParam& betaModelParam = GetPricingModel()->GetModelParams()->GetModelParam(ARM_ModelParamType::BetaCorrelation);

	itsFixBetaTimes.resize(0);
	itsFixBetaTimes.resize(0);

	for (size_t i = 0; i < betaModelParam.size(); ++i)
	{
		itsFixBetaTimes.push_back(betaModelParam.GetTimeAtPoint(i));
		itsFixBetaValues.push_back(betaModelParam.GetValueAtPoint(i));
	}

}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
