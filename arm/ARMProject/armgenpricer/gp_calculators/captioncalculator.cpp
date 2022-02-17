/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file captioncalculator.h
 *
 *  \base class for Caption calculators
 *	\author  y khlif
 *	\version 1.0
 *	\date November 2004
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/captioncalculator.h"
#include "gpcalculators/argconvdefault.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/autocleaner.h"
#include "gpbase/singleton.h"
#include "gpbase/datestrip.h"
#include "gpbase/stringconvert.h"

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

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaarg.h"
#include "gpcalib/kerneltogp.h"

/// gpmodels
#include "gpmodels/modelparamssfrm.h"
#include "gpmodels/modelparamssfrmdiag.h"
#include "gpmodels/modelparamssfrmfactory.h"
#include "gpmodels/sfrm.h"
#include "gpmodels/ForwardMarginBasis.h"
#include "gpmodels/ForwardMarginIR.h"
#include "gpmodels/forwardforex.h"
#include "gpmodels/hybridbasisfwdir.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/compositegen.h"
#include "gpnumlib/antitheticgen.h"

/// gpnummethods
#include "gpnummethods/treemethod.h"
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

/// gpbase (kernel)
//#include <inst/swaption.h>
//#include <inst/portfolio.h>
//#include <inst/forex.h>
#include <glob/paramview.h>
//#include <inst/fixleg.h>
//#include <mod/y2cmodel.h>

/// STL
#include <iomanip> /// for setprecision()
#include <memory>

CC_BEGIN_NAMESPACE( ARM )

const string ARM_CaptionCalculator::CaptionColNamesTable [] =
{
    "EventDate",
	"PayDate",    
    "FWDStartDate",
	"FWDEndDate",
	"IT",
	"Strike",
    "CpnSpread",
    "Notional",
    "FundStartDate",
	"FundEndDate",
	"FundSpread",
    "FundNotional",
	"NotionalExchange",
    "CapIndex",
    "CpnNPVIntermediate",
	"CpnNPV",
    "FundFragmentNPV",
	"FundNPV",
	"Caplet",
	"CapIntermediate",
    "Cap",
	"FeesPayDate",
	"Fees",
    "TotalNPV",
	"ExerciseFees",
	"BermudaProfile",
	"Caption",
	"CaptionStrike",
	"NumeraireDate",
	"ExerciseCondition",
	"ExerciseConditionOfIndex",
	"ProbaOfExercise"	
};

// SFRM default factors number
const int SFRM_NB_FACTORS			= 1;
const int SFRM_VOL_TYPE				= K_DIAG;

const double NON_CALL_FEE=1000000000000000.0;

// Tree default pricing values
const int TREE_NBSTEPS_PER_YEAR	=30;
const double STD_DEV_RATIO		=5.0;
const double MIN_STD_DEV		=0.001;

/// SFRM sigma range [1bp,10000bp] with a 500bp default value
const double SIGMA_LOWER_BOUND      = 0.0001;
const double SIGMA_UPPER_BOUND      = 1.0;
const double SIGMA_DEFAULT_VALUE    = 0.20;

/// SFRM beta range [0%,130%] with a 100% default value
const double BETA_DEFAULT_VALUE    = 1;
const double BETA_LOWER_BOUND	   = 0.01;
const double BETA_UPPER_BOUND	   = 1.3;

/// SFRM MRS range [-15%,50%] with a 2% default value
const double MRS_LOWER_BOUND        = -0.20;
const double MRS_UPPER_BOUND        = 0.20;
const double MRS_DEFAULT_VALUE      = 0.02;

// Minimum Strike Value for Cap Floor calibration
const double MIN_CF_STRIKE			= 0.0025;

/// 0.001bp of vega to be selected in portfolio for volatility bootstrapping 
const double VEGA_MIN_TO_SELECT=0.0000001;
const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

const double CF_DEFAULT_PRICE=1.0e+100;
const double CF_DEFAULT_PRECISION=0.0000001;
const double CF_DEFAULT_WEIGHT=1.0;

const double CF_PORTFOLIO_IDX=0;
const double DIGIT_PORTFOLIO_IDX=1;

const double DIGITAL_DEFAULT_WEIGHT=0.000001;

/// Default MDM key names
const string YC_KEY_NAME            = "YC_";
const string YC_BASIS_KEY_NAME      = "YC_BASIS_";
const string OSWMODEL_KEY_NAME      = "OSWMOD_";
const string CFMODEL_KEY_NAME       = "CFMOD_";
const string MRS_KEY_NAME           = "MRS_";
const string BETA_KEY_NAME          = "BETA_";
const string CORREL_KEY_NAME        = "CORREL_";
const string FOREX_KEY_NAME         = "FOREX_";

/// To indicate no link to any additional standard swaps for
/// equivalent strike/vol computation (diagonal swaption generation)
const int NO_OTHER_STDSWAPS=-1;

/// Equivalent datas computed for diagonal swaption
const unsigned int OSW_SWAP_RATE    = 0;
const unsigned int OSW_TARGET_VOL	= 1;
const unsigned int OSW_TARGET_PRICE = 2;
const unsigned int OSW_TARGET_VEGA  = 3;
const unsigned int OSW_NB_EQUIVDATA = 4;

/// Reference schedules for Caption date structure
const unsigned int Principal_SCHED =0;
const unsigned int NB_Caption_SCHED =1;




//  1)        create Calculator object
//  2)        Create generic security
//  3)        Model Create
//  4)        Calibration create
//  5)        Pricing
//  6)		  Get data & Update


///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///                                            1)        create Calculator object
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////

ARM_CaptionCalculator::ARM_CaptionCalculator(const ARM_Date& startDate,
			const ARM_Date& endDate,
            const string& cpnIdxTerm,
			int payRec,
			int CF,
			const ARM_Curve& couponProfile,
			const string& FundIdxTerm,
			int NotifDays,
			int NonCall,
			const ARM_Curve& exerciseProfile,		
			const ARM_Curve& notionalProfile,
			int cpnDayCount,
			int cpnResetTiming,
			const string& cpnResetCal,
			const string& cpnPayCal,
			const ARM_Curve& cpnSpreadProfile,
			int fundDayCount,
			const string& fundResetCal,
			const string& fundPayCal,
			const ARM_Curve& fundSpreadProfile,
			int factorNb,
			int SFRMVolType,
			CalibrationMode SWOPTCalibMode,
			CalibrationMode BETACalibFlag,
			std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice,
			const ARM_MarketData_ManagerRep& mktDataManager,
			const ARM_StringVector& keys) : ARM_GenCalculator(mktDataManager),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsCpnIdxTerm(cpnIdxTerm),
	itsPayRec(payRec),
	itsCapOrFloor(CF),
	itsCoupon(couponProfile),  
	itsFundIdxTerm(FundIdxTerm),
	itsNotifDays(NotifDays),
	itsNonCall(NonCall),	
	itsNotional(notionalProfile),
	itsCpnDayCount(cpnDayCount),
	itsCpnResetTiming(cpnResetTiming),
	itsCpnResetCal(cpnResetCal),         
    itsCpnPayCal(cpnPayCal),           
	itsCpnSpread(cpnSpreadProfile),
	itsFundDayCount(fundDayCount),
	itsFundResetCal(fundResetCal),         
    itsFundPayCal(fundPayCal),
	itsFundSpread(fundSpreadProfile),
	itsPayRecFund(-payRec),
	itsFactorNb(factorNb),
	itsSFRMVolType(SFRMVolType),
	itsSwoptCalibMode(SWOPTCalibMode),		
	itsBetaCalibFlag(BETACalibFlag),	
	itsProductsToPrice(productsToPrice),
	itsHasBeenPriced(false),
	itsFirstCallIdx(0),
	itsCaptionPrice(0.0),
	itsCapPrice(0),
	itsCouponLegPrice(0),
	itsFundingLegPrice(0),
	itsCaptionStrikes(0),
	itsExerciseProbas(0.0),
	itsNbSteps(100),
	itsCpnDateStrip(NULL),
	itsFundDateStrip(NULL),
	itsCalibPorfolioDateStrip(NULL),
	itsPrincipalDateStrip(NULL),
	itsStrike(couponProfile),
	itsCapFloorPF(NULL),
	itsSwaptionPF(NULL),
	itsNbdDeedResets(0)
{
	ARM_Timer timer;

	SetName(ARM_CAPTION);

    /// Set keys for MDM datas access
    if(keys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(keys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey] = newMdMKeys[YcKey];
        newMdMKeys[BasisKey]   = newMdMKeys[YcKey];
        SetKeys(newMdMKeys);
    }
    else
        SetKeys(keys);

	// 12M is 
	if (itsCpnIdxTerm == "12M")
	{
		itsCpnIdxTerm = "1Y";
	}
	itsCpnFreq  = ARM_ArgConv_MatFrequency.GetNumber(itsCpnIdxTerm);
	itsFundFreq = ARM_ArgConv_MatFrequency.GetNumber(itsFundIdxTerm);
    /// Check input datas
    CheckDataAndTimeIt();

    /// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit( static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]))->GetCurrencyUnit() );

    /// Set funding & basis currencies
    itsFundingCcy = *(static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[FundingKey]))->GetCurrencyUnit());
    itsBasisCcy = *(static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[BasisKey]))->GetCurrencyUnit());

	itsCpnIndexType = GetCpnIndexType();
	ARM_INDEX_TYPE fundIndexType = GetFundIndexType();

	itsNotifDays = (NotifDays == GETDEFAULTVALUE ? GetCurrencyUnit()->GetSpotDays(): NotifDays );
    itsCpnDayCount = (cpnDayCount == GETDEFAULTVALUE ? GetCurrencyUnit()->GetMMDayCount(): cpnDayCount );
	itsCpnResetCal= (cpnResetCal == GETDEFAULTVALUESTR ? GetCurrencyUnit()->GetResetCalName(itsCpnIndexType) : cpnResetCal );         
    itsCpnPayCal = (cpnPayCal == GETDEFAULTVALUESTR ? GetCurrencyUnit()->GetPayCalName(itsCpnIndexType) : cpnPayCal );  
	itsCpnResetTiming = (cpnResetTiming == GETDEFAULTVALUE ? K_ADVANCE : cpnResetTiming );

	itsFundDayCount = (fundDayCount == GETDEFAULTVALUE ? itsFundingCcy.GetMMDayCount(): fundDayCount );
	itsFundResetCal= (fundResetCal == GETDEFAULTVALUESTR ? itsFundingCcy.GetResetCalName(fundIndexType) : fundResetCal );         
    itsFundPayCal = (fundPayCal == GETDEFAULTVALUESTR ? itsFundingCcy.GetPayCalName(fundIndexType) : fundPayCal );  
	
	GenerateCapStrikes();
	GenerateEquivalentFundSpreads();
	
	GenerateExerciseSyle(exerciseProfile);
	
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

    /// Create the Generic Security
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[itsCpnModelKey]);	 

    /// Create the SFRM 1F pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping and
    /// mean reversion optimization
	std::vector<double>& resetDates = itsPrincipalDateStrip->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int resetSize=resetDates->size();
	int i=0;
	while(i < resetSize && (*resetDates)[i] <= asOfDate)
        ++i;
	itsNbdDeedResets = i;

    CreateAndSetCalibrationAndTimeIt();
}


ARM_CaptionCalculator::ARM_CaptionCalculator( const ARM_Date& asOfDate,
			const ARM_Date& startDate,
			const ARM_Date& endDate,
			const string& cpnIdxTerm,
			int payRec,
			int CF,
			const ARM_Curve& couponProfile,
			const string& FundIdxTerm,
			int NotifDays,
			int NonCall,
			const ARM_Curve& exerciseProfile,		
			const ARM_Curve& notionalProfile,
			int cpnDayCount,
			int cpnResetTiming,
			const string& cpnResetCal,
			const string& cpnPayCal,
			const ARM_Curve& cpnSpreadProfile,
			int fundDayCount,
			const string& fundResetCal,
			const string& fundPayCal,
			const ARM_Curve& fundSpreadProfile):
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsCpnIdxTerm(cpnIdxTerm),
	itsPayRec(payRec),
	itsCapOrFloor(CF),
	itsCoupon(couponProfile),  
	itsFundIdxTerm(FundIdxTerm),
	itsNotifDays(NotifDays),
	itsNonCall(NonCall),	
	itsNotional(notionalProfile),
	itsCpnDayCount(cpnDayCount),
	itsCpnResetTiming(cpnResetTiming),
	itsCpnResetCal(cpnResetCal),         
    itsCpnPayCal(cpnPayCal),           
	itsCpnSpread(cpnSpreadProfile),
	itsFundDayCount(fundDayCount),
	itsFundResetCal(fundResetCal),         
    itsFundPayCal(fundPayCal),
	itsFundSpread(fundSpreadProfile),
	itsPayRecFund(-payRec),
	itsFactorNb(0),
	itsSFRMVolType(0),
	itsSwoptCalibMode(NOCalib),		
	itsBetaCalibFlag(NOCalib),	
	itsProductsToPrice(NULL),
	itsHasBeenPriced(false),
	itsFirstCallIdx(0),
	itsCaptionPrice(0.0),
	itsCapPrice(0),
	itsCouponLegPrice(0),
	itsFundingLegPrice(0),
	itsCaptionStrikes(0),
	itsExerciseProbas(0.0),
	itsNbSteps(100),
	itsCpnDateStrip(NULL),
	itsFundDateStrip(NULL),
	itsCalibPorfolioDateStrip(NULL),
	itsPrincipalDateStrip(NULL),
	itsStrike(couponProfile),
	itsCapFloorPF(NULL),
	itsSwaptionPF(NULL),
	itsNbdDeedResets(0)
{
	SetName(ARM_CAPTION);

	// 12M is 
	if (itsCpnIdxTerm == "12M")
	{
		itsCpnIdxTerm = "1Y";
	}
	itsCpnFreq  = ARM_ArgConv_MatFrequency.GetNumber(itsCpnIdxTerm);
	itsFundFreq = ARM_ArgConv_MatFrequency.GetNumber(itsFundIdxTerm);

	itsCpnIndexType = GetCpnIndexType();
	ARM_INDEX_TYPE fundIndexType = GetFundIndexType();

	GenerateExerciseSyle(exerciseProfile);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_BermudaSwaptionCalculator
///	Routine: InitBermudaSwaptionForSummit
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::InitCaptionForSummit(ARM_StringVector& mdmKeys,
											ARM_MarketData_ManagerRep* mktDataManager,
											int factorNb,
											int SFRMVolType,
											CalibrationMode SWOPTCalibMode,
											CalibrationMode BETACalibFlag,
											std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice)
{
	ARM_Timer timer;

    /// Set keys for MDM datas access
    if(mdmKeys.size() < NbKeys)
    {
        /// To be compatible with non basis version
        ARM_StringVector newMdMKeys(mdmKeys);
        newMdMKeys.resize(NbKeys);
        newMdMKeys[FundingKey] = newMdMKeys[YcKey];
        newMdMKeys[BasisKey]   = newMdMKeys[YcKey];
        SetKeys(newMdMKeys);
    }
    else
        SetKeys(mdmKeys);

	InitializeMktDataManagerOnly(*mktDataManager);

    /// Check input datas
    CheckDataAndTimeIt();

	itsFactorNb = factorNb;
	itsSFRMVolType = SFRMVolType;
	itsSwoptCalibMode = SWOPTCalibMode;
	itsBetaCalibFlag = BETACalibFlag;
	itsProductsToPrice = productsToPrice;

    /// Set the coupon/payment currency (inherited from ARM_Security)
    SetCurrencyUnit( static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[YcKey]))->GetCurrencyUnit() );

    /// Set funding & basis currencies
    itsFundingCcy = *(static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[FundingKey]))->GetCurrencyUnit());
    itsBasisCcy = *(static_cast< ARM_ZeroCurve* >(GetMktDataManager()->GetData(GetKeys()[BasisKey]))->GetCurrencyUnit());

	if (itsNotifDays == GETDEFAULTVALUE)
		itsNotifDays = GetCurrencyUnit()->GetSpotDays();

    if (itsCpnDayCount == GETDEFAULTVALUE)
		itsCpnDayCount = GetCurrencyUnit()->GetMMDayCount();

	if (itsCpnResetCal == GETDEFAULTVALUESTR)
		itsCpnResetCal = GetCurrencyUnit()->GetResetCalName(itsCpnIndexType);

    if (itsCpnPayCal == GETDEFAULTVALUESTR)
		itsCpnPayCal = GetCurrencyUnit()->GetPayCalName(itsCpnIndexType);

	if (itsCpnResetTiming == GETDEFAULTVALUE)
		itsCpnResetTiming = K_ADVANCE;

	if (itsFundDayCount == GETDEFAULTVALUE)
		itsFundDayCount = itsFundingCcy.GetMMDayCount();

	itsCpnIndexType = GetCpnIndexType();
	ARM_INDEX_TYPE fundIndexType = GetFundIndexType();

	if (itsFundResetCal == GETDEFAULTVALUESTR)
		itsFundResetCal = itsFundingCcy.GetResetCalName(fundIndexType);
	
    if (itsFundPayCal == GETDEFAULTVALUESTR)
		itsFundPayCal = itsFundingCcy.GetPayCalName(fundIndexType);  
	
	GenerateCapStrikes();
	GenerateEquivalentFundSpreads();
	
    /// Set the coupon & funding model name alias i.e. keys to access to
    /// the right model for each GP keyword (1st argument)
    SetModelKeys();

    /// Create the Generic Security
    CreateAndSetDealDescriptionAndTimeIt(GetKeys()[itsCpnModelKey]);	 

    /// Create the SFRM 1F pricing model with its default parameters
    CreateAndSetModelAndTimeIt();

    /// Create the calibration set for volatility bootstapping and
    /// mean reversion optimization
	std::vector<double>& resetDates = itsPrincipalDateStrip->GetResetDates();
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int resetSize = resetDates->size();
	int i = 0;
	while (i < resetSize && (*resetDates)[i] <= asOfDate)
        ++i;
	itsNbdDeedResets = i;

    CreateAndSetCalibrationAndTimeIt();
}
																  

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_CaptionCalculator::ARM_CaptionCalculator(const ARM_CaptionCalculator& rhs)
: ARM_GenCalculator(rhs),
	itsStartDate(rhs.itsStartDate),
	itsEndDate(rhs.itsEndDate),
	itsCpnIdxTerm(rhs.itsCpnIdxTerm),
	itsPayRec(rhs.itsPayRec),
	itsCapOrFloor(rhs.itsCapOrFloor),
	itsCoupon(rhs.itsCoupon),  
	itsFundIdxTerm(rhs.itsFundIdxTerm),
	itsNotifDays(rhs.itsNotifDays),
	itsNonCall(rhs.itsNonCall),
	itsExerStyle(rhs.itsExerStyle),
	itsNotional(rhs.itsNotional),
	itsCpnDayCount(rhs.itsCpnDayCount),
	itsCpnResetTiming(rhs.itsCpnResetTiming),
	itsCpnResetCal(rhs.itsCpnResetCal),         
    itsCpnPayCal(rhs.itsCpnPayCal),           
	itsCpnSpread(rhs.itsCpnSpread),
	itsFundDayCount(rhs.itsFundDayCount),
	itsFundResetCal(rhs.itsFundResetCal),         
    itsFundPayCal(rhs.itsFundPayCal),
	itsFundSpread(rhs.itsFundSpread),
	itsPayRecFund(rhs.itsPayRecFund),
	itsFactorNb(rhs.itsFactorNb),
	itsSFRMVolType(rhs.itsSFRMVolType),
	itsSwoptCalibMode(rhs.itsSwoptCalibMode),		
	itsBetaCalibFlag(rhs.itsBetaCalibFlag),	
	itsProductsToPrice(rhs.itsProductsToPrice),
	itsHasBeenPriced(rhs.itsHasBeenPriced),
	itsFirstCallIdx(rhs.itsFirstCallIdx),
	itsCaptionPrice(rhs.itsCaptionPrice),
	/*itsCapPrice(rhs.itsCapPrice),
	itsCouponLegPrice(rhs.itsCouponLegPrice),
	itsFundingLegPrice(rhs.itsFundingLegPrice),
	itsCaptionStrikes(rhs.itsCaptionStrikes),*/
	itsCapPrice(ARM_VectorPtr((rhs.itsCapPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCapPrice->Clone() : NULL)),
	itsCouponLegPrice(ARM_VectorPtr((rhs.itsCouponLegPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCouponLegPrice->Clone() : NULL)),
	itsFundingLegPrice(ARM_VectorPtr((rhs.itsFundingLegPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsFundingLegPrice->Clone() : NULL)),
	itsCaptionStrikes(ARM_VectorPtr((rhs.itsCaptionStrikes != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCaptionStrikes->Clone() : NULL)),
	itsExerciseProbas(rhs.itsExerciseProbas),
	itsNbSteps(rhs.itsNbSteps),
	itsCpnFreq(rhs.itsCpnFreq),
	itsFundFreq(rhs.itsFundFreq),
	itsFundingCcy(rhs.itsFundingCcy),
	itsBasisCcy(rhs.itsBasisCcy),
	itsCpnModelKey(rhs.itsCpnModelKey),
	itsFundingModelKey(rhs.itsFundingModelKey),
	itsBasisRefModelKey(rhs.itsBasisRefModelKey),
	itsCpnDateStrip(rhs.itsCpnDateStrip),
	itsFundDateStrip(rhs.itsFundDateStrip),
	itsCalibPorfolioDateStrip(rhs.itsCalibPorfolioDateStrip),
	itsPrincipalDateStrip(rhs.itsPrincipalDateStrip),
	itsStrike(rhs.itsStrike),
	/*itsCapFloorPF(rhs.itsCapFloorPF),
	itsSwaptionPF(rhs.itsSwaptionPF),*/
	itsCapFloorPF(ARM_StdPortfolioPtr((rhs.itsCapFloorPF != ARM_StdPortfolioPtr(NULL)) ? (ARM_StdPortfolio*) rhs.itsCapFloorPF->Clone() : NULL)),
	itsSwaptionPF(ARM_StdPortfolioPtr((rhs.itsSwaptionPF != ARM_StdPortfolioPtr(NULL)) ? (ARM_StdPortfolio*) rhs.itsSwaptionPF->Clone() : NULL)),
	itsNbdDeedResets(rhs.itsNbdDeedResets),
	itsCpnIndexType(rhs.itsCpnIndexType)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: clone
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////

ARM_Object* ARM_CaptionCalculator::Clone() const
{
	return new ARM_CaptionCalculator(*this);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: assignment operator
///	Returns: 
///	Action : 
/////////////////////////////////////////////////////////////////
ARM_CaptionCalculator& ARM_CaptionCalculator::operator=(const ARM_CaptionCalculator& rhs)
{
	if (this != &rhs)
	{
		ARM_GenCalculator::operator=(rhs);

		itsStartDate				= rhs.itsStartDate;
		itsEndDate					= rhs.itsEndDate;
		itsCpnIdxTerm				= rhs.itsCpnIdxTerm;
		itsPayRec					= rhs.itsPayRec;
		itsCapOrFloor				= rhs.itsCapOrFloor;
		itsCoupon					= rhs.itsCoupon;
		itsFundIdxTerm				= rhs.itsFundIdxTerm;
		itsNotifDays				= rhs.itsNotifDays;
		itsNonCall					= rhs.itsNonCall;
		itsExerStyle				= rhs.itsExerStyle;
		itsNotional					= rhs.itsNotional;
		itsCpnDayCount				= rhs.itsCpnDayCount;
		itsCpnResetTiming			= rhs.itsCpnResetTiming;
		itsCpnResetCal				= rhs.itsCpnResetCal;     
		itsCpnPayCal				= rhs.itsCpnPayCal;         
		itsCpnSpread				= rhs.itsCpnSpread;
		itsFundDayCount				= rhs.itsFundDayCount;
		itsFundResetCal				= rhs.itsFundResetCal;         
		itsFundPayCal				= rhs.itsFundPayCal;
		itsFundSpread				= rhs.itsFundSpread;
		itsPayRecFund				= rhs.itsPayRecFund;
		itsFactorNb					= rhs.itsFactorNb;
		itsSFRMVolType				= rhs.itsSFRMVolType;
		itsSwoptCalibMode			= rhs.itsSwoptCalibMode;		
		itsBetaCalibFlag			= rhs.itsBetaCalibFlag;	
		itsProductsToPrice			= rhs.itsProductsToPrice;
		itsHasBeenPriced			= rhs.itsHasBeenPriced;
		itsFirstCallIdx				= rhs.itsFirstCallIdx;
		itsCaptionPrice				= rhs.itsCaptionPrice;
		/*itsCapPrice					= rhs.itsCapPrice;
		itsCouponLegPrice			= rhs.itsCouponLegPrice;
		itsFundingLegPrice			= rhs.itsFundingLegPrice;
		itsCaptionStrikes			= rhs.itsCaptionStrikes;*/
		itsExerciseProbas			= rhs.itsExerciseProbas;
		itsNbSteps					= rhs.itsNbSteps;
		itsCpnModelKey				= rhs.itsCpnModelKey;
		itsFundingCcy				= rhs.itsFundingCcy;
		itsBasisCcy					= rhs.itsBasisCcy;
		itsFundingModelKey			= rhs.itsFundingModelKey;
		itsBasisRefModelKey			= rhs.itsBasisRefModelKey;
		itsFundFreq					= rhs.itsFundFreq;
		itsCpnFreq					= rhs.itsCpnFreq;
		itsCpnDateStrip				= rhs.itsCpnDateStrip;
		itsFundDateStrip			= rhs.itsFundDateStrip;
		itsCalibPorfolioDateStrip	= rhs.itsCalibPorfolioDateStrip;
		itsPrincipalDateStrip		= rhs.itsPrincipalDateStrip,
		itsStrike					= rhs.itsStrike;
		itsCapFloorPF				= rhs.itsCapFloorPF;
		itsSwaptionPF				= rhs.itsSwaptionPF;
		itsNbdDeedResets			= rhs.itsNbdDeedResets;
		itsCpnIndexType				= rhs.itsCpnIndexType;
		itsCapPrice					= ARM_VectorPtr((rhs.itsCapPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCapPrice->Clone() : NULL);
		itsCouponLegPrice			= ARM_VectorPtr((rhs.itsCouponLegPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCouponLegPrice->Clone() : NULL);
		itsFundingLegPrice			= ARM_VectorPtr((rhs.itsFundingLegPrice != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsFundingLegPrice->Clone() : NULL);
		itsCaptionStrikes			= ARM_VectorPtr((rhs.itsCaptionStrikes != ARM_VectorPtr(NULL)) ? (std::vector<double>&) rhs.itsCaptionStrikes->Clone() : NULL);


	}

	return *this;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_CaptionCalculator::~ARM_CaptionCalculator()
{	
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if Caption datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CheckData()
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

	// Caption parameters checking

	if ( itsEndDate < itsStartDate )
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Start Date of the deal has to be before the end date.");
	}

	if (!(itsCoupon > K_NEW_DOUBLE_TOL))
	{
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The Strike has to be strictly positive.");
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the TARN. The
///          DateStripCombiner merges event dates of each
///          legs
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CaptionCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event dates of the Caption. The
///          DateStripCombiner merges event dates of each
///          legs
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_CaptionCalculator::DatesStructure() const
{

	/// The Caption is mono-currency at the moment : the currency
    /// is then the payement one & it is saved at the ARM_Security level

    /// Get reset & payment calendars
	const char* cpnResetCalendar   = itsCpnResetCal.c_str();
	const char* cpnPayCalendar     = itsCpnPayCal.c_str();

	int payFreq     =   itsCpnFreq;
	int payGap      =   GETDEFAULTVALUE;
	int payTiming   =   K_ARREARS;
	
	// principal datestrip
	int principalResetTiming =   K_ADVANCE;
	ARM_Date principalEndDate = itsEndDate;
	if (itsCpnResetTiming == K_ARREARS)
	{
		char* ccyName=GetCurrencyUnit()->GetCcyName();
		ARM_Date tmpDate(principalEndDate);
		tmpDate.AddPeriod(itsCpnFreq,ccyName); 
		principalEndDate = tmpDate;
	}

	ARM_DateStrip PrincipalDateStrip(itsStartDate,
									principalEndDate,
									itsCpnFreq,						
									itsCpnDayCount,					
									cpnResetCalendar,			
									K_MOD_FOLLOWING,
									K_ADJUSTED, 
									K_SHORTSTART,						
									-itsNotifDays,						
									payFreq,						
									payGap,						
									cpnPayCalendar,			
									principalResetTiming,					
									payTiming);


	if( const_cast< ARM_CaptionCalculator* >(this)->itsPrincipalDateStrip == ARM_DateStripPtr(NULL))
		const_cast< ARM_CaptionCalculator* >(this)->itsPrincipalDateStrip = ARM_DateStripPtr(new ARM_DateStrip(PrincipalDateStrip));	

   
    
	//Coupon datestrip 
	int cpnResetGap = GetCurrencyUnit()->GetSpotDays();
	ARM_DateStrip CouponDateStrip(itsStartDate,
									itsEndDate,
									itsCpnFreq,						
									itsCpnDayCount,					
									cpnResetCalendar,			
									K_MOD_FOLLOWING,
									K_ADJUSTED, 
									K_SHORTSTART,						
									-cpnResetGap,						
									payFreq,						
									payGap,						
									cpnPayCalendar,			
									itsCpnResetTiming,					
									payTiming);

	if( const_cast< ARM_CaptionCalculator* >(this)->itsCpnDateStrip == ARM_DateStripPtr(NULL))
		const_cast< ARM_CaptionCalculator* >(this)->itsCpnDateStrip = ARM_DateStripPtr(new ARM_DateStrip(CouponDateStrip));	
	
	//Funding datestrip 
	int fundResetGap = itsFundingCcy.GetSpotDays();
	const char* fundResetCalendar   = itsFundResetCal.c_str();
	const char* fundPayCalendar     = itsFundPayCal.c_str();

	ARM_DateStrip FundDateStrip(itsStartDate,
								  itsEndDate,
								  itsCpnFreq,						
								  itsFundDayCount,					
								  fundResetCalendar,			
								  K_MOD_FOLLOWING,
								  K_ADJUSTED, 
								  K_SHORTSTART,						
								  -fundResetGap,						
								  payFreq,						
								  payGap,						
								  fundPayCalendar,			
								  principalResetTiming,					
								  payTiming);

	if( const_cast< ARM_CaptionCalculator* >(this)->itsFundDateStrip == ARM_DateStripPtr(NULL))
		const_cast< ARM_CaptionCalculator* >(this)->itsFundDateStrip = ARM_DateStripPtr( new ARM_DateStrip(FundDateStrip));

	char* stdResetCalendar   = GetCurrencyUnit()->GetResetCalName(itsCpnIndexType);

	char* stdPayCalendar   = GetCurrencyUnit()->GetPayCalName(itsCpnIndexType);


	//calib portfolio datestrip 
	ARM_DateStrip CalibPorfolioDateStrip(itsStartDate,
									principalEndDate,
									itsCpnFreq,						
									itsCpnDayCount,					
									stdResetCalendar,			
									K_MOD_FOLLOWING,
									K_ADJUSTED, 
									K_SHORTSTART,						
									-cpnResetGap,						
									payFreq,						
									payGap,						
									stdPayCalendar,			
									principalResetTiming,					
									payTiming);

	if( const_cast< ARM_CaptionCalculator* >(this)->itsCalibPorfolioDateStrip == ARM_DateStripPtr(NULL))
		const_cast< ARM_CaptionCalculator* >(this)->itsCalibPorfolioDateStrip = ARM_DateStripPtr(new ARM_DateStrip(CalibPorfolioDateStrip));	

	
	/// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(NB_Caption_SCHED,NULL);
    SchedVect[Principal_SCHED] = &PrincipalDateStrip;
	
	ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

    /// Initialise the memorisation of the first call line
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	CC_MUTABLE( ARM_CaptionCalculator, itsFirstCallIdx ) = 0;
	while(itsFirstCallIdx < eventDates->size() && (*eventDates)[itsFirstCallIdx] <= asOfDate)
        ++CC_MUTABLE( ARM_CaptionCalculator, itsFirstCallIdx );

	delete stdResetCalendar;
	delete stdPayCalendar;

    return EventSchedule;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetCpnIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CaptionCalculator::GetCpnIndexType()
{
    string liborTypeName(string(GetCurrencyUnit()->GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string cpnIndexTerm(itsCpnIdxTerm);
    if(cpnIndexTerm=="12M")
        cpnIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += cpnIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetFundIndexType
///	Returns: the index type
///	Action : get the ARM index type of the coupon index
/////////////////////////////////////////////////////////////////
ARM_INDEX_TYPE ARM_CaptionCalculator::GetFundIndexType()
{
    string liborTypeName(string(itsFundingCcy.GetCcyName()) == string("EUR") ? "EURIBOR" : "LIBOR");
    string fundIndexTerm(itsFundIdxTerm);
    if(fundIndexTerm=="12M")
        fundIndexTerm="1Y"; // because we cant convert EURIBOR1Y and not EURIBOR12M
    liborTypeName += fundIndexTerm;
    return static_cast< ARM_INDEX_TYPE > (ARM_ArgConv_IndexType.GetNumber(liborTypeName));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GenerateExerciseSyle
///	Returns: void
///	Action : generate Exercise Schedule
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::GenerateExerciseSyle(const ARM_Curve& exerciseProfile)
{
	std::vector<double> Abscisses = const_cast< ARM_Curve& >(exerciseProfile).GetAbscisses();
	std::vector<double> Ordinates = const_cast< ARM_Curve& >(exerciseProfile).GetOrdinates();
	int size = Abscisses.size();
	int newSize = size+2;
	std::vector<double> newAbscisses(newSize); 
	std::vector<double> newOrdinates(newSize);

	newAbscisses[0] = Abscisses[0]-1.0;	
	newOrdinates[0] = NON_CALL_FEE;	
    for(int j = 1; j<newSize-1;j++)
	{
		newAbscisses[j] = Abscisses[j-1];	
		newOrdinates[j] = Ordinates[j-1];
	}
	
	newAbscisses[newSize-1] = Abscisses[size-1]+1.0;	
	newOrdinates[newSize-1] = NON_CALL_FEE;	

	const_cast< ARM_Curve& >(itsExerStyle).SetAbscisses(newAbscisses);
	const_cast< ARM_Curve& >(itsExerStyle).SetOrdinates(newOrdinates);
	//ARM_Interpolator<double, double>* interpolator = new ARM_StepUpRightOpenCstExtrapolDble ;
	const_cast< ARM_Curve& >(itsExerStyle).SetInterpolator(new ARM_StepUpRightOpenCstExtrapolDble);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GenerateCapStrikes
///	Returns: void
///	Action : generate Cap strikes (coupons-cpnSpreads)  Schedule
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::GenerateCapStrikes(void)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		
	std::vector<double> CouponAbscisses = const_cast< ARM_Curve& >(itsStrike).GetAbscisses();
	int size = CouponAbscisses.size();

	std::vector<double> CouponOrdinates(size);
	
	for(int i = 0; i<size;i++)
	{
		double abscisseDate= CouponAbscisses[i];
		double coupon = const_cast< ARM_Curve& >(itsStrike).Interpolate(abscisseDate);
		double spread = const_cast< ARM_Curve& >(itsCpnSpread).Interpolate(abscisseDate);
		
		CouponOrdinates[i] = coupon-spread;		
	}

	const_cast< ARM_Curve& >(itsStrike).SetOrdinates(CouponOrdinates);

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GenerateEquivalentFundSpreads
///	Returns: void
///	Action : generate equivalent fundSpreads structure with freq freqCpn 
///          in the case where freqFund>freqCpn)
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::GenerateEquivalentFundSpreads(void)
{
	if(itsFundFreq>itsCpnFreq)
	{
		string fundingCcyName(itsFundingCcy.GetCcyName());
		string basisCcyName(itsBasisCcy.GetCcyName());

		double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
		ARM_ZeroCurve* fundZeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[FundingKey]);
		ARM_ZeroCurve* basisZeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[BasisKey]);	

		ARM_ZeroCurve* basisFundCurve = NULL; 

		if(basisCcyName==fundingCcyName)
		{
			basisFundCurve = basisZeroCurve;
		}
		else
		{
			basisFundCurve = fundZeroCurve;
		}

		ARM_Y2CModel* fundYCModel = new ARM_Y2CModel(fundZeroCurve,basisFundCurve); 
		ARM_INDEX_TYPE fundLiborType = GetFundIndexType();

		
		//Funding datestrip 
		int fundResetGap = itsFundingCcy.GetSpotDays();
		const char* fundResetCalendar   = itsFundResetCal.c_str();
		const char* fundPayCalendar     = itsFundPayCal.c_str();

		int payFreq     =   itsCpnFreq;
		int payGap      =   GETDEFAULTVALUE;
		int payTiming   =   K_ARREARS;
		int ResetTiming =   K_ADVANCE;
		

		ARM_DateStrip FundDateStrip(itsStartDate,
									  itsEndDate,
									  itsCpnFreq,						
									  itsFundDayCount,					
									  fundResetCalendar,			
									  K_MOD_FOLLOWING,
									  K_ADJUSTED, 
									  K_SHORTSTART,						
									  -fundResetGap,						
									  payFreq,						
									  payGap,						
									  fundPayCalendar,			
									  ResetTiming,					
									  payTiming);

		std::vector<double>& fundStartDates = FundDateStrip.GetFlowStartDates();
		std::vector<double>& fundEndDates = FundDateStrip.GetFlowEndDates();
		std::vector<double>& resetDates = FundDateStrip.GetResetDates();
		int newSize = resetDates->size();
		
		std::vector<double> SpreadAbscisses(newSize);
		std::vector<double> spreadOrdinates(newSize);

		//transfor fundSpread to refvalue
		std::vector<double> FundSpreadAbscisses = const_cast< ARM_Curve& >(itsFundSpread).GetAbscisses();
		std::vector<double> FundSpreadOrdinates = const_cast< ARM_Curve& >(itsFundSpread).GetOrdinates();
		int fundSize = FundSpreadAbscisses.size();
		ARM_Vector SpreadRefvalueAbs(fundSize);
		ARM_Vector spreadRefvalueOrd(fundSize);
		for(int j = 0; j<fundSize;j++)
		{
			SpreadRefvalueAbs[j] = FundSpreadAbscisses[j]+asOfDate;	
			spreadRefvalueOrd[j] = FundSpreadOrdinates[j]*100.0;
		}		

		string interpolatorNameSpread = const_cast< ARM_Curve& >(itsFundSpread).GetInterpolator()->toString();
		int interpolMethodSpread  = ARM_ArgConv_interpolCurveType.GetNumber(interpolatorNameSpread);		
		ARM_ReferenceValue* refValueFundSpread = new ARM_ReferenceValue((&SpreadRefvalueAbs), (&spreadRefvalueOrd));
		refValueFundSpread->SetCalcMethod(interpolMethodSpread);

		//transfor notional to refvalue
		std::vector<double> notionalAbscisses = const_cast< ARM_Curve& >(itsNotional).GetAbscisses();		
		std::vector<double> notionalOrdinates = const_cast< ARM_Curve& >(itsNotional).GetOrdinates();
		int notionalSize = notionalAbscisses.size();
		ARM_Vector notionalRefvalueAbs(notionalSize);
		ARM_Vector notionalRefvalueOrd(notionalSize);
		for(j = 0; j<notionalSize;j++)
		{
			notionalRefvalueAbs[j] = notionalAbscisses[j]+asOfDate;	
			notionalRefvalueOrd[j] = notionalOrdinates[j];
		}		
		
		string interpolatorNameNotional = const_cast< ARM_Curve& >(itsNotional).GetInterpolator()->toString();
		int interpolMethodNotional  = K_STEPUP_RIGHT;
		ARM_ReferenceValue* refValueNotional = new ARM_ReferenceValue((&notionalRefvalueAbs), (&notionalRefvalueOrd));
		refValueNotional->SetCalcMethod(interpolMethodNotional);
		
		for(int i =0 ; i<newSize; i++)
		{
			double dStartDate = (*fundStartDates)[i];
			ARM_Date fundStartDate(dStartDate);

			double dEndDate = (*fundEndDates)[i];
			ARM_Date fundEndDate(dEndDate);

			ARM_FixLeg fundFixLeg(fundStartDate,
						   fundEndDate, 
						   refValueFundSpread,
						   K_RCV, 
						   itsFundFreq,
						   itsFundDayCount,
						   K_COMP_PROP,
						   K_ARREARS, 
						   K_ADJUSTED,
						   K_SHORTSTART,
						   &itsFundingCcy);

			fundFixLeg.SetAmount(refValueNotional);
			fundFixLeg.SetModel(fundYCModel);
			double fundspreadNPV = fundFixLeg.ComputePrice();

			double spread_1 =1.0;
			ARM_FixLeg fundFixLegAtSpread_1(fundStartDate,
						   fundEndDate, 
						   spread_1,
						   K_RCV, 
						   itsFundFreq,
						   itsFundDayCount,
						   K_COMP_PROP,
						   K_ARREARS, 
						   K_ADJUSTED,
						   K_SHORTSTART,
						   &itsFundingCcy);

			fundFixLegAtSpread_1.SetAmount(refValueNotional);

			
			fundFixLegAtSpread_1.SetModelVariable(NULL);
			fundFixLegAtSpread_1.SetModel(fundYCModel);
			double fundspread_1NPV = fundFixLegAtSpread_1.ComputePrice();

			double eqSpread = fundspreadNPV/fundspread_1NPV;

			SpreadAbscisses[i] = (*resetDates)[i]-asOfDate;	
			spreadOrdinates[i] = eqSpread/100.0;
		}

		const_cast< ARM_Curve& >(itsFundSpread).SetAbscisses(SpreadAbscisses);
		const_cast< ARM_Curve& >(itsFundSpread).SetOrdinates(spreadOrdinates);	
	}
		
}

			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///                                                    2)      Deal Description
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: InitPriceableColumns
///	Returns: nothing
///	Action : Feed the MDM with input market objects
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

	rowDescVec[Caption] = zeroValue;
    rowTypeVec[Caption] = ARM_DOUBLE;

	rowDescVec[Cap] = zeroValue;
    rowTypeVec[Cap] = ARM_DOUBLE;

    rowDescVec[CpnNPV] = zeroValue;
    rowTypeVec[CpnNPV] = ARM_DOUBLE;

	rowDescVec[FundNPV] = zeroValue;
    rowTypeVec[FundNPV] = ARM_DOUBLE;

	rowDescVec[CaptionStrike] = zeroValue;
    rowTypeVec[CaptionStrike] = ARM_DOUBLE;

	rowDescVec[ProbaOfExercise] = zeroValue;
    rowTypeVec[ProbaOfExercise] = ARM_DOUBLE;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CaptionCalculator::ColumnNames() const
{
	size_t colNamesSize = sizeof(CaptionColNamesTable)/sizeof(CaptionColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = CaptionColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_CaptionCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
    ARM_DateStripPtr principalDateStrip = datesStructure.GetDateStrip(Principal_SCHED);
	ARM_DateStripPtr cpnDateStrip = itsCpnDateStrip;
	ARM_DateStripPtr fundDateStrip = itsFundDateStrip;
	
    size_t eventSize = principalDateStrip->GetResetDates()->size();
    size_t descSize = sizeof(CaptionColNamesTable)/sizeof(CaptionColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();   

    /// Convention conversion : K_QUATERLY => 3M for instance
    string fundMatFreq=ARM_ArgConvReverse_MatFrequency.GetString(itsFundFreq);
    string fundDayCount=ARM_ArgConvReverse_DayCount.GetString(itsFundDayCount);	
	string liborPayTiming("ARR");
	if(itsCpnResetTiming != K_ADVANCE)
	{
		liborPayTiming = "ADV";
	}
    
    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);

	string cpnModelName = GetKeys()[itsCpnModelKey];
    string fundingModelName = GetKeys()[itsFundingModelKey];


	double spotAsOfValue = 1.0;
    string spotFutureValue("");
    if(fundingModelName != cpnModelName)
    {
		ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        string foreignCcy(forex->GetMoneyCurrency()->GetCcyName());
		//string domesticCcy(forex->GetMoneyCurrency()->GetCcyName());

		//string(GetCurrencyUnit()->GetCcyName())
        spotAsOfValue = static_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]))->GetMarketPrice();
        if( foreignCcy == string(itsFundingCcy.GetCcyName()))
            spotFutureValue = "/SPOT(" + GetKeys()[ForexKey] + ")";
        else
        {
            spotAsOfValue = 1.0/spotAsOfValue;
            spotFutureValue = "*SPOT(" + GetKeys()[ForexKey] + ")";
        }
    }

	bool isArrearsCpn = (itsCpnResetTiming==K_ARREARS);
	bool isFirstEvent = (eventIdx==itsFirstCallIdx);
	
	// EventDate 
    double eventDate=(*(principalDateStrip->GetResetDates()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[EventDate] = eventDateDesc.str();
    rowTypeVec[EventDate] = ARM_DATE_TYPE;
   

	double	fundflowPayDate;	
	double  fundStartDate;
	double  fundEndDate;
	double  fundResetDate;

	double cpnResetDate;
	double cpnflowPayDate;
	double fwdStartDate;
	double fwdEndDate;
	double InterestTerm;

	double zeroFlow = 0.0;
			

	if (isArrearsCpn)
	{
		bool isFundLastEvent = (eventIdx==(eventSize-2));
		bool isCpnLastEvent = (eventIdx==(eventSize-1));

		if(isFirstEvent)
		{
			fundflowPayDate=(*(cpnDateStrip->GetPaymentDates()))[eventIdx];
			fundStartDate=(*(fundDateStrip->GetFlowStartDates()))[eventIdx];
			fundEndDate=(*(fundDateStrip->GetFlowEndDates()))[eventIdx];
			fundResetDate=(*(fundDateStrip->GetResetDates()))[eventIdx];

			//PayDate
			//FSD
			//FED
			//It
			//Strike
			//Spread
			
			//Notional		
					
			//FundStartDate		
			CC_Ostringstream fundStartDesc;
			fundStartDesc << CC_NS(std,fixed) << fundStartDate;
			rowDescVec[FundStartDate] = fundStartDesc.str();
			rowTypeVec[FundStartDate] = ARM_DATE_TYPE;

			//FundEndDate		
			CC_Ostringstream fundEndDesc;
			fundEndDesc << CC_NS(std,fixed) << fundEndDate;
			rowDescVec[FundEndDate] = fundEndDesc.str();
			rowTypeVec[FundEndDate] = ARM_DATE_TYPE;

			//FundSpread		
			CC_Ostringstream fundSpreadDesc;
			fundSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundSpread).Interpolate(fundResetDate-asOfDate);
			rowDescVec[FundSpread] = fundSpreadDesc.str();
			rowTypeVec[FundSpread] = ARM_DOUBLE;

			//FundNotional
			CC_Ostringstream fundNominalDesc;
			double nominalValue = const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate);
			nominalValue *= spotAsOfValue;
			fundNominalDesc << CC_NS(std,fixed) << nominalValue << spotFutureValue;
			rowDescVec[FundNotional] = fundNominalDesc.str();
			rowTypeVec[FundNotional] = ARM_STRING;
			
			//NotionalExchange
			CC_Ostringstream fundNominalExchangeDesc;
			double nominalValueExchange =  const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate);
			fundNominalExchangeDesc << "(-DF(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
			fundNominalExchangeDesc << "+DF(" << fundingModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
			fundNominalExchangeDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
			fundNominalExchangeDesc << "+(DF(" << cpnModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
			fundNominalExchangeDesc << "-DF(" << cpnModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
			fundNominalExchangeDesc << "*" << CC_NS(std,fixed) << nominalValueExchange;
			rowDescVec[NotionalExchange] = fundNominalExchangeDesc.str();
			rowTypeVec[NotionalExchange] = ARM_STRING;

			//CapIndex
			//CpnNPVIntermediate

			//CpnNPV
			CC_Ostringstream cpnNPVDesc;
			cpnNPVDesc << "PV(" << CaptionColNamesTable[CpnNPVIntermediate] << "[i+1])";
			rowDescVec[CpnNPV] = cpnNPVDesc.str();
			rowTypeVec[CpnNPV] = ARM_STRING;

			//FundFragmentNPV
			double strikeFund =0.0;
			string fundSwap_POrR = "P";//ARM_ArgConvReverse_PayRec.GetString(-itsPayRecFund);
			CC_Ostringstream fundSwapDesc;
			fundSwapDesc << "Swap(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i],";
			fundSwapDesc << CaptionColNamesTable[FundEndDate] << "[i],";
			fundSwapDesc << CC_NS(std,fixed) << strikeFund << ",";
			fundSwapDesc << fundSwap_POrR << ",";
			fundSwapDesc << "," << ",";
			fundSwapDesc << fundMatFreq << ",";
			fundSwapDesc << fundDayCount << ",";
			fundSwapDesc <<CaptionColNamesTable[FundSpread] << "[i])";
			fundSwapDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
			fundSwapDesc << "+" << CaptionColNamesTable[NotionalExchange] << "[i]";
			rowDescVec[FundFragmentNPV] = fundSwapDesc.str();
			rowTypeVec[FundFragmentNPV] = ARM_STRING;


			// FundNPV
			CC_Ostringstream fundNPVDesc;
			fundNPVDesc << CaptionColNamesTable[FundFragmentNPV] << "[i]";
			if(!isFundLastEvent) 
				fundNPVDesc << "+PV(" << CaptionColNamesTable[FundNPV] << "[i+1]" << ")";
			
			rowDescVec[FundNPV] = fundNPVDesc.str();
			rowTypeVec[FundNPV] = ARM_STRING;

			//Caplet

			//Cap
			CC_Ostringstream capDesc;
			capDesc << "PV(" << CaptionColNamesTable[CapIntermediate] << "[i+1])";
			rowDescVec[Cap] = capDesc.str();
			rowTypeVec[Cap] = ARM_STRING;	

			//TotalNPV  LS_Times_PayRec*(  (C/F)*LeCapOuFloor - CpnNPV + FundNPV - PayRecFees_Times_PayRec*NPVFees)
			CC_Ostringstream totalNPVDesc;
			int  LS_Times_PayRec = -itsPayRec;
			totalNPVDesc << "(" << CC_NS(std,fixed) << LS_Times_PayRec << ")*(";
			totalNPVDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*";
			totalNPVDesc << CaptionColNamesTable[Cap] << "[i]";
			totalNPVDesc << "-" << CaptionColNamesTable[CpnNPV] << "[i]";
			totalNPVDesc << "+" << CaptionColNamesTable[FundNPV] << "[i])";
			rowDescVec[TotalNPV] = totalNPVDesc.str();
			rowTypeVec[TotalNPV] = ARM_STRING;			

			//FeesPayDate
			CC_Ostringstream feesPayDesc;
			feesPayDesc << CC_NS(std,fixed) << fundflowPayDate;
			rowDescVec[FeesPayDate] = feesPayDesc.str();
			rowTypeVec[FeesPayDate] = ARM_DATE_TYPE;

			//Fees		
			CC_Ostringstream feesDesc;
			feesDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsExerStyle).Interpolate(eventDate-asOfDate);
			rowDescVec[Fees] = feesDesc.str();
			rowTypeVec[Fees] = ARM_DOUBLE;

			//ExerciseFees
			CC_Ostringstream exerciseFeesDesc;
			exerciseFeesDesc << CaptionColNamesTable[Fees] << "[i]*";
			exerciseFeesDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[FeesPayDate] << "[i])";
			rowDescVec[ExerciseFees] = exerciseFeesDesc.str();
			rowTypeVec[ExerciseFees] = ARM_STRING;

			//BermudaProfile
			CC_Ostringstream bermudaDesc;		
			bermudaDesc << "MAX(";
			bermudaDesc << CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i],";
			if(!isFundLastEvent)
				bermudaDesc << "PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1]))";
			else
				bermudaDesc << "0)";
			
			rowDescVec[BermudaProfile] = bermudaDesc.str();
			rowTypeVec[BermudaProfile] = ARM_STRING;

			//Caption
			CC_Ostringstream captionDesc;		
			captionDesc << "MAX(";
			captionDesc << CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i],";
			if(!isFundLastEvent)
				captionDesc << "PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1]))";
			else
				captionDesc << "0)";
			
			rowDescVec[Caption] = captionDesc.str();
			rowTypeVec[Caption] = ARM_STRING;

			//CaptionStrike	(C/F)*(   CpnNPV - FundNPV + coeffBeforeFees*ExerciseFees  )
			CC_Ostringstream captionStrikeDesc;	
			int coeffBeforeFees = -itsPayRec;
			captionStrikeDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*(";
			captionStrikeDesc <<  CaptionColNamesTable[CpnNPV] << "[i]";
			captionStrikeDesc << "-" << CaptionColNamesTable[FundNPV] << "[i]";
			captionStrikeDesc << "+(" << CC_NS(std,fixed) << coeffBeforeFees << ")*"; 
			captionStrikeDesc << CaptionColNamesTable[ExerciseFees] << "[i])";		
			rowDescVec[CaptionStrike] = captionStrikeDesc.str();
			rowTypeVec[CaptionStrike] = ARM_STRING;

			//NumeraireDate
			double numeraireDate = (*(cpnDateStrip->GetFwdRateEndDates()))[(*(cpnDateStrip->GetFwdRateEndDates())).size()-1];
			CC_Ostringstream numeraireDateDesc;
			numeraireDateDesc << CC_NS(std,fixed) << numeraireDate;
			rowDescVec[NumeraireDate] = numeraireDateDesc.str();
			rowTypeVec[NumeraireDate] = ARM_DATE_TYPE;

			//ExerciseCondition
			CC_Ostringstream ExerciseConditionDesc;		
			ExerciseConditionDesc << "IF(";
			ExerciseConditionDesc << "("<< CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i]";
			if(!isFundLastEvent)
				ExerciseConditionDesc << "-PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1])";
			ExerciseConditionDesc << ")>=0";
			ExerciseConditionDesc << ",1,0)";
			rowDescVec[ExerciseCondition] = ExerciseConditionDesc.str();
			rowTypeVec[ExerciseCondition] = ARM_STRING;

			//ExerciseConditionOfIndex
			CC_Ostringstream ExerciseConditionOfIndexDesc;
			if(eventIdx<=itsNonCall)
			{					
				ExerciseConditionOfIndexDesc << "MAX(";
				ExerciseConditionOfIndexDesc << CaptionColNamesTable[ExerciseCondition] << "[i]*";
				ExerciseConditionOfIndexDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[NumeraireDate] << "[i])";
				if((!isFundLastEvent)&&(eventIdx != itsNonCall))
					ExerciseConditionOfIndexDesc << ",PV(" << CaptionColNamesTable[ExerciseConditionOfIndex] << "[i+1]))";
				else
					ExerciseConditionOfIndexDesc << ",0)";				
			}
			else
				ExerciseConditionOfIndexDesc <<CC_NS(std,fixed) << zeroFlow;

			rowDescVec[ExerciseConditionOfIndex] = ExerciseConditionOfIndexDesc.str();
			rowTypeVec[ExerciseConditionOfIndex] = ARM_STRING;

			//ProbaOfExercise
			CC_Ostringstream ProbaOfExerciseDesc;		
			ProbaOfExerciseDesc << "MAX(";
			ProbaOfExerciseDesc << CaptionColNamesTable[ExerciseCondition] << "[i]*";
			ProbaOfExerciseDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[NumeraireDate] << "[i])";
			if(!isFundLastEvent)
				ProbaOfExerciseDesc << ",PV(" << CaptionColNamesTable[ExerciseConditionOfIndex] << "[i+1]))";
			else
				ProbaOfExerciseDesc << ",0)";
			
			rowDescVec[ProbaOfExercise] = ProbaOfExerciseDesc.str();
			rowTypeVec[ProbaOfExercise] = ARM_STRING;

		}
		else //if isn't the first event
		{
			cpnflowPayDate=(*(cpnDateStrip->GetPaymentDates()))[eventIdx-1];
			cpnResetDate=(*(cpnDateStrip->GetResetDates()))[eventIdx-1];
			fwdStartDate = (*(cpnDateStrip->GetFwdRateStartDates()))[eventIdx-1];
			fwdEndDate = (*(cpnDateStrip->GetFwdRateEndDates()))[eventIdx-1];
			InterestTerm = (*(cpnDateStrip->GetInterestTerms()))[eventIdx-1];
			
			
			//PayDate
			CC_Ostringstream payDateDesc;
			payDateDesc << CC_NS(std,fixed) << cpnflowPayDate;
			rowDescVec[PayDate] = payDateDesc.str();
			rowTypeVec[PayDate] = ARM_DATE_TYPE;
			
			//FWDStartDate
			CC_Ostringstream FSDDesc;
			FSDDesc << CC_NS(std,fixed) << fwdStartDate;
			rowDescVec[FWDStartDate] = FSDDesc.str();
			rowTypeVec[FWDStartDate] = ARM_DATE_TYPE;

			//FWDEndDate
			CC_Ostringstream FEDDesc;
			FEDDesc << CC_NS(std,fixed) << fwdEndDate;
			rowDescVec[FWDEndDate] = FEDDesc.str();
			rowTypeVec[FWDEndDate] = ARM_DATE_TYPE;

			//IT
			CC_Ostringstream ITDesc;
			ITDesc << CC_NS(std,fixed) << InterestTerm;
			rowDescVec[IT] = ITDesc.str();
			rowTypeVec[IT] = ARM_DOUBLE;


			//Strike
			CC_Ostringstream strikeDesc;
			strikeDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsStrike).Interpolate(cpnResetDate-asOfDate);
			rowDescVec[Strike] = strikeDesc.str();
			rowTypeVec[Strike] = ARM_DOUBLE;
			
			//CpnSpread
			CC_Ostringstream cpnSpreadDesc;
			cpnSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsCpnSpread).Interpolate(cpnResetDate-asOfDate);
			rowDescVec[CpnSpread] = cpnSpreadDesc.str();
			rowTypeVec[CpnSpread] = ARM_DOUBLE;

			//Notional		
			CC_Ostringstream cpnNominalDesc;
			cpnNominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNotional).Interpolate(cpnflowPayDate-asOfDate);
			rowDescVec[Notional] = cpnNominalDesc.str();
			rowTypeVec[Notional] = ARM_DOUBLE;
			
			if(!isCpnLastEvent) //nothing to do for funding leg if is the last event
			{
				fundflowPayDate=(*(cpnDateStrip->GetPaymentDates()))[eventIdx];
				fundStartDate=(*(fundDateStrip->GetFlowStartDates()))[eventIdx];
				fundEndDate=(*(fundDateStrip->GetFlowEndDates()))[eventIdx];
				fundResetDate=(*(fundDateStrip->GetResetDates()))[eventIdx];
			
				//FundStartDate		
				CC_Ostringstream fundStartDesc;
				fundStartDesc << CC_NS(std,fixed) << fundStartDate;
				rowDescVec[FundStartDate] = fundStartDesc.str();
				rowTypeVec[FundStartDate] = ARM_DATE_TYPE;

				//FundEndDate		
				CC_Ostringstream fundEndDesc;
				fundEndDesc << CC_NS(std,fixed) << fundEndDate;
				rowDescVec[FundEndDate] = fundEndDesc.str();
				rowTypeVec[FundEndDate] = ARM_DATE_TYPE;

				//FundSpread		
				CC_Ostringstream fundSpreadDesc;
				fundSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundSpread).Interpolate(fundResetDate-asOfDate);
				rowDescVec[FundSpread] = fundSpreadDesc.str();
				rowTypeVec[FundSpread] = ARM_DOUBLE;

				//FundNotional
				CC_Ostringstream fundNominalDesc;
				double nominalValue = const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate);
				nominalValue *= spotAsOfValue;
				fundNominalDesc << CC_NS(std,fixed) << nominalValue << spotFutureValue;
				rowDescVec[FundNotional] = fundNominalDesc.str();
				rowTypeVec[FundNotional] = ARM_STRING;
					
				//NotionalExchange
				CC_Ostringstream fundNominalExchangeDesc;
				double nominalValueExchange =  const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate);
				fundNominalExchangeDesc << "(-DF(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
				fundNominalExchangeDesc << "+DF(" << fundingModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
				fundNominalExchangeDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
				fundNominalExchangeDesc << "+(DF(" << cpnModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
				fundNominalExchangeDesc << "-DF(" << cpnModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
				fundNominalExchangeDesc << "*" << CC_NS(std,fixed) << nominalValueExchange;

				rowDescVec[NotionalExchange] = fundNominalExchangeDesc.str();
				rowTypeVec[NotionalExchange] = ARM_STRING;
			}

			//CapIndex
			CC_Ostringstream cpnIndexDesc;
			cpnIndexDesc << "LIBOR(" << cpnModelName << "," << CaptionColNamesTable[FWDStartDate] << "[i],";
			cpnIndexDesc << CaptionColNamesTable[FWDEndDate] << "[i],,,,";
			cpnIndexDesc << liborPayTiming << ")";
			rowDescVec[CapIndex] = cpnIndexDesc.str();
			rowTypeVec[CapIndex] = ARM_STRING;
			

			//CpnNPVIntermediate
			CC_Ostringstream CpnNPVIntermediateDesc;
			CpnNPVIntermediateDesc<< "(" << CaptionColNamesTable[CapIndex] << "[i]";
			CpnNPVIntermediateDesc << "+" << CaptionColNamesTable[CpnSpread] << "[i])";
			CpnNPVIntermediateDesc << "*" << CaptionColNamesTable[IT] << "[i]";
			CpnNPVIntermediateDesc << "*" << CaptionColNamesTable[Notional] << "[i]"; 
			CpnNPVIntermediateDesc << "*DF(" << cpnModelName << "," << CaptionColNamesTable[PayDate] << "[i])";
			if(!isCpnLastEvent)
				CpnNPVIntermediateDesc << "+" <<"PV(" << CaptionColNamesTable[CpnNPVIntermediate] << "[i+1])";
					
			rowDescVec[CpnNPVIntermediate] = CpnNPVIntermediateDesc.str();
			rowTypeVec[CpnNPVIntermediate] = ARM_STRING;
					
			if(!isCpnLastEvent) //nothing to do for funding leg if is the last event
			{
				//CpnNPV
				CC_Ostringstream cpnNPVDesc;
				cpnNPVDesc << "PV(" << CaptionColNamesTable[CpnNPVIntermediate] << "[i+1])";
				rowDescVec[CpnNPV] = cpnNPVDesc.str();
				rowTypeVec[CpnNPV] = ARM_STRING;
						
				//FundFragmentNPV
				double strikeFund =0.0;
				string fundSwap_POrR = "P";
				CC_Ostringstream fundSwapDesc;
				fundSwapDesc << "Swap(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i],";
				fundSwapDesc << CaptionColNamesTable[FundEndDate] << "[i],";
				fundSwapDesc << CC_NS(std,fixed) << strikeFund << ",";
				fundSwapDesc << fundSwap_POrR << ",";
				fundSwapDesc << "," << ",";
				fundSwapDesc << fundMatFreq << ",";
				fundSwapDesc << fundDayCount << ",";
				fundSwapDesc <<CaptionColNamesTable[FundSpread] << "[i])";
				fundSwapDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
				fundSwapDesc << "+" << CaptionColNamesTable[NotionalExchange] << "[i]";
				rowDescVec[FundFragmentNPV] = fundSwapDesc.str();
				rowTypeVec[FundFragmentNPV] = ARM_STRING;

				// FundNPV
				CC_Ostringstream fundNPVDesc;
				fundNPVDesc << CaptionColNamesTable[FundFragmentNPV] << "[i]";
				if(!isFundLastEvent) 
					fundNPVDesc << "+PV(" << CaptionColNamesTable[FundNPV] << "[i+1]" << ")";
				
				rowDescVec[FundNPV] = fundNPVDesc.str();
				rowTypeVec[FundNPV] = ARM_STRING;
			}			

			//Caplet
			CC_Ostringstream CapletDesc;
			CapletDesc << "MAX(";
			CapletDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*";
			CapletDesc << "("  << CaptionColNamesTable[CapIndex] << "[i]";
			CapletDesc << "-" << CaptionColNamesTable[Strike] << "[i])";
			CapletDesc << ",0)";
			CapletDesc << "*" << CaptionColNamesTable[IT] << "[i]";
			CapletDesc << "*" << CaptionColNamesTable[Notional] << "[i]"; 
			CapletDesc << "*DF(" << cpnModelName << "," << CaptionColNamesTable[PayDate] << "[i])";

			rowDescVec[Caplet] = CapletDesc.str();
			rowTypeVec[Caplet] = ARM_STRING;
			
			//CapIntermediate
			CC_Ostringstream capIntermediateDesc;
			capIntermediateDesc << CaptionColNamesTable[Caplet] << "[i]";
			if(!isCpnLastEvent)
				capIntermediateDesc << "+PV(" << CaptionColNamesTable[CapIntermediate] << "[i+1])";
			
			rowDescVec[CapIntermediate] = capIntermediateDesc.str();
			rowTypeVec[CapIntermediate] = ARM_STRING;

			if(!isCpnLastEvent) //no capflow, no fees,  at the last the last event
			{
				//Cap
				CC_Ostringstream capDesc;
				capDesc << "PV(" << CaptionColNamesTable[CapIntermediate] << "[i+1])";
				rowDescVec[Cap] = capDesc.str();
				rowTypeVec[Cap] = ARM_STRING;
			
				//TotalNPV  LS_Times_PayRec*(  (C/F)*LeCapOuFloor - CpnNPV + FundNPV - PayRecFees_Times_PayRec*NPVFees)
				CC_Ostringstream totalNPVDesc;
				int  LS_Times_PayRec = -itsPayRec;
				totalNPVDesc << "(" << CC_NS(std,fixed) << LS_Times_PayRec << ")*(";
				totalNPVDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*";
				totalNPVDesc << CaptionColNamesTable[Cap] << "[i]";
				totalNPVDesc << "-" << CaptionColNamesTable[CpnNPV] << "[i]";
				totalNPVDesc << "+" << CaptionColNamesTable[FundNPV] << "[i])";
				rowDescVec[TotalNPV] = totalNPVDesc.str();
				rowTypeVec[TotalNPV] = ARM_STRING;			

				//FeesPayDate
				CC_Ostringstream feesPayDesc;
				feesPayDesc << CC_NS(std,fixed) << fundflowPayDate;
				rowDescVec[FeesPayDate] = feesPayDesc.str();
				rowTypeVec[FeesPayDate] = ARM_DATE_TYPE;

				//Fees		
				CC_Ostringstream feesDesc;
				feesDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsExerStyle).Interpolate(eventDate-asOfDate);
				rowDescVec[Fees] = feesDesc.str();
				rowTypeVec[Fees] = ARM_DOUBLE;

				//ExerciseFees
				CC_Ostringstream exerciseFeesDesc;
				exerciseFeesDesc << CaptionColNamesTable[Fees] << "[i]*";
				exerciseFeesDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[FeesPayDate] << "[i])";
				rowDescVec[ExerciseFees] = exerciseFeesDesc.str();
				rowTypeVec[ExerciseFees] = ARM_STRING;

				//CaptionStrike	(C/F)*(   CpnNPV - FundNPV + coeffBeforeFees*ExerciseFees  )
				CC_Ostringstream captionStrikeDesc;	
				int coeffBeforeFees = -itsPayRec;
				captionStrikeDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*(";
				captionStrikeDesc <<  CaptionColNamesTable[CpnNPV] << "[i]";
				captionStrikeDesc << "-" << CaptionColNamesTable[FundNPV] << "[i]";
				captionStrikeDesc << "+(" << CC_NS(std,fixed) << coeffBeforeFees << ")*"; 
				captionStrikeDesc << CaptionColNamesTable[ExerciseFees] << "[i])";		
				rowDescVec[CaptionStrike] = captionStrikeDesc.str();
				rowTypeVec[CaptionStrike] = ARM_STRING;				

				//NumeraireDate
				double numeraireDate = (*(cpnDateStrip->GetFwdRateEndDates()))[(*(cpnDateStrip->GetFwdRateEndDates())).size()-1];
				CC_Ostringstream numeraireDateDesc;
				numeraireDateDesc << CC_NS(std,fixed) << numeraireDate;
				rowDescVec[NumeraireDate] = numeraireDateDesc.str();
				rowTypeVec[NumeraireDate] = ARM_DATE_TYPE;

				//ExerciseCondition
				CC_Ostringstream ExerciseConditionDesc;		
				ExerciseConditionDesc << "IF(";
				ExerciseConditionDesc << "("<< CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i]";
				if(!isFundLastEvent)
					ExerciseConditionDesc << "-PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1])";
				ExerciseConditionDesc << ")>=0";
				ExerciseConditionDesc << ",1,0)";
				rowDescVec[ExerciseCondition] = ExerciseConditionDesc.str();
				rowTypeVec[ExerciseCondition] = ARM_STRING;

				//ExerciseConditionOfIndex
				CC_Ostringstream ExerciseConditionOfIndexDesc;
				if(eventIdx<=itsNonCall)
				{					
					ExerciseConditionOfIndexDesc << "MAX(";
					ExerciseConditionOfIndexDesc << CaptionColNamesTable[ExerciseCondition] << "[i]*";
					ExerciseConditionOfIndexDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[NumeraireDate] << "[i])";
					if((!isFundLastEvent)&&(eventIdx != itsNonCall))
						ExerciseConditionOfIndexDesc << ",PV(" << CaptionColNamesTable[ExerciseConditionOfIndex] << "[i+1]))";
					else
						ExerciseConditionOfIndexDesc << ",0)";
				}
				else
					ExerciseConditionOfIndexDesc <<CC_NS(std,fixed) << zeroFlow;
				
					
				rowDescVec[ExerciseConditionOfIndex] = ExerciseConditionOfIndexDesc.str();
				rowTypeVec[ExerciseConditionOfIndex] = ARM_STRING;

			}			

			//BermudaProfile			
			CC_Ostringstream bermudaDesc;
			if(!isCpnLastEvent)
			{
				bermudaDesc << "MAX(";
				bermudaDesc << CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i],";
				if(!isFundLastEvent)
					bermudaDesc << "PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1]))";
				else
					bermudaDesc << "0)";
			}
			else
				bermudaDesc <<CC_NS(std,fixed) << zeroFlow;
			
			
			rowDescVec[BermudaProfile] = bermudaDesc.str();
			rowTypeVec[BermudaProfile] = ARM_STRING;	
		}
	}
	else //ADV case
	{
		bool isLastEvent = (eventIdx==(eventSize-1));

		
		fundflowPayDate=(*(cpnDateStrip->GetPaymentDates()))[eventIdx];
		fundStartDate=(*(fundDateStrip->GetFlowStartDates()))[eventIdx];
		fundEndDate=(*(fundDateStrip->GetFlowEndDates()))[eventIdx];
		fundResetDate=(*(fundDateStrip->GetResetDates()))[eventIdx];

		cpnflowPayDate=(*(cpnDateStrip->GetPaymentDates()))[eventIdx];
		cpnResetDate=(*(cpnDateStrip->GetResetDates()))[eventIdx];
		fwdStartDate = (*(cpnDateStrip->GetFwdRateStartDates()))[eventIdx];
		fwdEndDate = (*(cpnDateStrip->GetFwdRateEndDates()))[eventIdx];
		InterestTerm = (*(cpnDateStrip->GetInterestTerms()))[eventIdx];
		
		
		//PayDate
		CC_Ostringstream payDateDesc;
		payDateDesc << CC_NS(std,fixed) << cpnflowPayDate;
		rowDescVec[PayDate] = payDateDesc.str();
		rowTypeVec[PayDate] = ARM_DATE_TYPE;
		
		//FWDStartDate
		CC_Ostringstream FSDDesc;
		FSDDesc << CC_NS(std,fixed) << fwdStartDate;
		rowDescVec[FWDStartDate] = FSDDesc.str();
		rowTypeVec[FWDStartDate] = ARM_DATE_TYPE;

		//FWDEndDate
		CC_Ostringstream FEDDesc;
		FEDDesc << CC_NS(std,fixed) << fwdEndDate;
		rowDescVec[FWDEndDate] = FEDDesc.str();
		rowTypeVec[FWDEndDate] = ARM_DATE_TYPE;

		//IT
		CC_Ostringstream ITDesc;
		ITDesc << CC_NS(std,fixed) << InterestTerm;
		rowDescVec[IT] = ITDesc.str();
		rowTypeVec[IT] = ARM_DOUBLE;

		//Strike
		CC_Ostringstream strikeDesc;
		strikeDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsStrike).Interpolate(cpnResetDate-asOfDate);
		rowDescVec[Strike] = strikeDesc.str();
		rowTypeVec[Strike] = ARM_DOUBLE;
		
		//CpnSpread
		CC_Ostringstream cpnSpreadDesc;
		cpnSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsCpnSpread).Interpolate(cpnResetDate-asOfDate);
		rowDescVec[CpnSpread] = cpnSpreadDesc.str();
		rowTypeVec[CpnSpread] = ARM_DOUBLE;

		//Notional		
		CC_Ostringstream cpnNominalDesc;
		cpnNominalDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsNotional).Interpolate(cpnflowPayDate-asOfDate);
		rowDescVec[Notional] = cpnNominalDesc.str();
		rowTypeVec[Notional] = ARM_DOUBLE;	
		
		//FundStartDate		
		CC_Ostringstream fundStartDesc;
		fundStartDesc << CC_NS(std,fixed) << fundStartDate;
		rowDescVec[FundStartDate] = fundStartDesc.str();
		rowTypeVec[FundStartDate] = ARM_DATE_TYPE;

		//FundEndDate		
		CC_Ostringstream fundEndDesc;
		fundEndDesc << CC_NS(std,fixed) << fundEndDate;
		rowDescVec[FundEndDate] = fundEndDesc.str();
		rowTypeVec[FundEndDate] = ARM_DATE_TYPE;

		//FundSpread		
		CC_Ostringstream fundSpreadDesc;
		fundSpreadDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsFundSpread).Interpolate(fundResetDate-asOfDate);
		rowDescVec[FundSpread] = fundSpreadDesc.str();
		rowTypeVec[FundSpread] = ARM_DOUBLE;

		//FundNotional
		CC_Ostringstream fundNominalDesc;
		double nominalValue = const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate);
		nominalValue *= spotAsOfValue;
		fundNominalDesc << CC_NS(std,fixed) << nominalValue << spotFutureValue;
		rowDescVec[FundNotional] = fundNominalDesc.str();
		rowTypeVec[FundNotional] = ARM_STRING;
		
		//NotionalExchange
		double nominalValueExchange =  const_cast< ARM_Curve& >(itsNotional).Interpolate(fundflowPayDate-asOfDate); 
		CC_Ostringstream fundNominalExchangeDesc;
		fundNominalExchangeDesc << "(-DF(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
		fundNominalExchangeDesc << "+DF(" << fundingModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
		fundNominalExchangeDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
		fundNominalExchangeDesc << "+(DF(" << cpnModelName << "," << CaptionColNamesTable[FundStartDate] << "[i])";
		fundNominalExchangeDesc << "-DF(" << cpnModelName << "," << CaptionColNamesTable[FundEndDate] << "[i]))";
		fundNominalExchangeDesc << "*" << CC_NS(std,fixed) << nominalValueExchange;
		rowDescVec[NotionalExchange] = fundNominalExchangeDesc.str();
		rowTypeVec[NotionalExchange] = ARM_STRING;

		//CapIndex
		CC_Ostringstream cpnIndexDesc;
		cpnIndexDesc << "LIBOR(" << cpnModelName << "," << CaptionColNamesTable[FWDStartDate] << "[i],";
		cpnIndexDesc << CaptionColNamesTable[FWDEndDate] << "[i],,,,";
		cpnIndexDesc << liborPayTiming << ")";
		rowDescVec[CapIndex] = cpnIndexDesc.str();
		rowTypeVec[CapIndex] = ARM_STRING;
		

		//CpnNPVIntermediate
		CC_Ostringstream CpnNPVIntermediateDesc;
		CpnNPVIntermediateDesc<< "(" << CaptionColNamesTable[CapIndex] << "[i]";
		CpnNPVIntermediateDesc << "+" << CaptionColNamesTable[CpnSpread] << "[i])";
		CpnNPVIntermediateDesc << "*" << CaptionColNamesTable[IT] << "[i]";
		CpnNPVIntermediateDesc << "*" << CaptionColNamesTable[Notional] << "[i]"; 
		CpnNPVIntermediateDesc << "*DF(" << cpnModelName << "," << CaptionColNamesTable[PayDate] << "[i])";
		if(!isLastEvent)
			CpnNPVIntermediateDesc << "+" <<"PV(" << CaptionColNamesTable[CpnNPVIntermediate] << "[i+1])";
				
		rowDescVec[CpnNPVIntermediate] = CpnNPVIntermediateDesc.str();
		rowTypeVec[CpnNPVIntermediate] = ARM_STRING;
		
		
	
		//CpnNPV
		CC_Ostringstream cpnNPVDesc;
		cpnNPVDesc << CaptionColNamesTable[CpnNPVIntermediate] << "[i]";
		rowDescVec[CpnNPV] = cpnNPVDesc.str();
		rowTypeVec[CpnNPV] = ARM_STRING;
				
		//FundFragmentNPV
		double strikeFund =0.0;
		string fundSwap_POrR = "P";
		CC_Ostringstream fundSwapDesc;
		fundSwapDesc << "Swap(" << fundingModelName << "," << CaptionColNamesTable[FundStartDate] << "[i],";
		fundSwapDesc << CaptionColNamesTable[FundEndDate] << "[i],";
		fundSwapDesc << CC_NS(std,fixed) << strikeFund << ",";
		fundSwapDesc << fundSwap_POrR << ",";
		fundSwapDesc << "," << ",";
		fundSwapDesc << fundMatFreq << ",";
		fundSwapDesc << fundDayCount << ",";
		fundSwapDesc <<CaptionColNamesTable[FundSpread] << "[i])";
		fundSwapDesc << "*" << CaptionColNamesTable[FundNotional] << "[i]";
		fundSwapDesc << "+" << CaptionColNamesTable[NotionalExchange] << "[i]";
		rowDescVec[FundFragmentNPV] = fundSwapDesc.str();
		rowTypeVec[FundFragmentNPV] = ARM_STRING;

		// FundNPV
		CC_Ostringstream fundNPVDesc;
		fundNPVDesc << CaptionColNamesTable[FundFragmentNPV] << "[i]";
		if(!isLastEvent) 
			fundNPVDesc << "+PV(" << CaptionColNamesTable[FundNPV] << "[i+1]" << ")";
		
		rowDescVec[FundNPV] = fundNPVDesc.str();
		rowTypeVec[FundNPV] = ARM_STRING;
					

		//Caplet
		CC_Ostringstream CapletDesc;
		CapletDesc << "MAX(";
		CapletDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*";
		CapletDesc << "("  << CaptionColNamesTable[CapIndex] << "[i]";
		CapletDesc << "-" << CaptionColNamesTable[Strike] << "[i])";
		CapletDesc << ",0)";
		CapletDesc << "*" << CaptionColNamesTable[IT] << "[i]";
		CapletDesc << "*" << CaptionColNamesTable[Notional] << "[i]"; 
		CapletDesc << "*DF(" << cpnModelName << "," << CaptionColNamesTable[PayDate] << "[i])";

		rowDescVec[Caplet] = CapletDesc.str();
		rowTypeVec[Caplet] = ARM_STRING;
		
		//CapIntermediate
		CC_Ostringstream capIntermediateDesc;
		capIntermediateDesc << CaptionColNamesTable[Caplet] << "[i]";
		if(!isLastEvent)
			capIntermediateDesc << "+PV(" << CaptionColNamesTable[CapIntermediate] << "[i+1])";
		
		rowDescVec[CapIntermediate] = capIntermediateDesc.str();
		rowTypeVec[CapIntermediate] = ARM_STRING;

		//Cap
		CC_Ostringstream capDesc;
		capDesc << CaptionColNamesTable[CapIntermediate] << "[i]";
		rowDescVec[Cap] = capDesc.str();
		rowTypeVec[Cap] = ARM_STRING;
	
		//TotalNPV  LS_Times_PayRec*(  (C/F)*LeCapOuFloor - CpnNPV + FundNPV - PayRecFees_Times_PayRec*NPVFees)
		CC_Ostringstream totalNPVDesc;
		int  LS_Times_PayRec = -itsPayRec;
		totalNPVDesc << "(" << CC_NS(std,fixed) << LS_Times_PayRec << ")*(";
		totalNPVDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*";
		totalNPVDesc << CaptionColNamesTable[Cap] << "[i]";
		totalNPVDesc << "-" << CaptionColNamesTable[CpnNPV] << "[i]";
		totalNPVDesc << "+" << CaptionColNamesTable[FundNPV] << "[i])";
		rowDescVec[TotalNPV] = totalNPVDesc.str();
		rowTypeVec[TotalNPV] = ARM_STRING;			

		//FeesPayDate
		CC_Ostringstream feesPayDesc;
		feesPayDesc << CC_NS(std,fixed) << fundflowPayDate;
		rowDescVec[FeesPayDate] = feesPayDesc.str();
		rowTypeVec[FeesPayDate] = ARM_DATE_TYPE;

		//Fees		
		CC_Ostringstream feesDesc;
		feesDesc << CC_NS(std,fixed) << const_cast< ARM_Curve& >(itsExerStyle).Interpolate(eventDate-asOfDate);
		rowDescVec[Fees] = feesDesc.str();
		rowTypeVec[Fees] = ARM_DOUBLE;

		//ExerciseFees
		CC_Ostringstream exerciseFeesDesc;
		exerciseFeesDesc << CaptionColNamesTable[Fees] << "[i]*";
		exerciseFeesDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[FeesPayDate] << "[i])";
		rowDescVec[ExerciseFees] = exerciseFeesDesc.str();
		rowTypeVec[ExerciseFees] = ARM_STRING;

		//CaptionStrike	(C/F)*(   CpnNPV - FundNPV + coeffBeforeFees*ExerciseFees  )
		CC_Ostringstream captionStrikeDesc;	
		int coeffBeforeFees = -itsPayRec;
		captionStrikeDesc << "(" << CC_NS(std,fixed) << itsCapOrFloor << ")*(";
		captionStrikeDesc <<  CaptionColNamesTable[CpnNPV] << "[i]";
		captionStrikeDesc << "-" << CaptionColNamesTable[FundNPV] << "[i]";
		captionStrikeDesc << "+(" << CC_NS(std,fixed) << coeffBeforeFees << ")*"; 
		captionStrikeDesc << CaptionColNamesTable[ExerciseFees] << "[i])";		
		rowDescVec[CaptionStrike] = captionStrikeDesc.str();
		rowTypeVec[CaptionStrike] = ARM_STRING;	

		//BermudaProfile
		double zeroFlow = 0.0;
		CC_Ostringstream bermudaDesc;
		bermudaDesc << "MAX(";
		bermudaDesc << CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i],";
		if(!isLastEvent)
			bermudaDesc << "PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1]))";
		else
			bermudaDesc << "0)";
		
		rowDescVec[BermudaProfile] = bermudaDesc.str();
		rowTypeVec[BermudaProfile] = ARM_STRING;

		//Caption		
		CC_Ostringstream captionDesc;
		if(isFirstEvent)
		{
			captionDesc << "MAX(";
			captionDesc << CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i],";
			if(!isLastEvent)
				captionDesc << "PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1]))";
			else
				captionDesc << "0)";
			
			rowDescVec[Caption] = captionDesc.str();
			rowTypeVec[Caption] = ARM_STRING;
		}

		//NumeraireDate
		double numeraireDate = (*(cpnDateStrip->GetFwdRateEndDates()))[(*(cpnDateStrip->GetFwdRateEndDates())).size()-1];
		CC_Ostringstream numeraireDateDesc;
		numeraireDateDesc << CC_NS(std,fixed) << numeraireDate;
		rowDescVec[NumeraireDate] = numeraireDateDesc.str();
		rowTypeVec[NumeraireDate] = ARM_DATE_TYPE;

		//ExerciseCondition
		CC_Ostringstream ExerciseConditionDesc;		
		ExerciseConditionDesc << "IF(";
		ExerciseConditionDesc << "("<< CaptionColNamesTable[TotalNPV] << "[i]-" << CaptionColNamesTable[ExerciseFees] << "[i]";
		if(!isLastEvent)
			ExerciseConditionDesc << "-PV(" << CaptionColNamesTable[BermudaProfile] << "[i+1])";
		ExerciseConditionDesc << ")>=0";
		ExerciseConditionDesc << ",1,0)";
		rowDescVec[ExerciseCondition] = ExerciseConditionDesc.str();
		rowTypeVec[ExerciseCondition] = ARM_STRING;

		//ExerciseConditionOfIndex
		CC_Ostringstream ExerciseConditionOfIndexDesc;
		if(eventIdx<=itsNonCall)
		{					
			ExerciseConditionOfIndexDesc << "MAX(";
			ExerciseConditionOfIndexDesc << CaptionColNamesTable[ExerciseCondition] << "[i]*";
			ExerciseConditionOfIndexDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[NumeraireDate] << "[i])";
			if((!isLastEvent)&&(eventIdx != itsNonCall))
				ExerciseConditionOfIndexDesc << ",PV(" << CaptionColNamesTable[ExerciseConditionOfIndex] << "[i+1]))";
			else
				ExerciseConditionOfIndexDesc << ",0)";
		}
		else
			ExerciseConditionOfIndexDesc <<CC_NS(std,fixed) << zeroFlow;		
			
		rowDescVec[ExerciseConditionOfIndex] = ExerciseConditionOfIndexDesc.str();
		rowTypeVec[ExerciseConditionOfIndex] = ARM_STRING;

		//ProbaOfExercise
		CC_Ostringstream ProbaOfExerciseDesc;
		if(isFirstEvent)
		{
			ProbaOfExerciseDesc << "MAX(";
			ProbaOfExerciseDesc << CaptionColNamesTable[ExerciseCondition] << "[i]*";
			ProbaOfExerciseDesc << "DF(" << cpnModelName << "," << CaptionColNamesTable[NumeraireDate] << "[i])";
			if(!isLastEvent)
				ProbaOfExerciseDesc << ",PV(" << CaptionColNamesTable[ExerciseConditionOfIndex] << "[i+1]))";
			else
				ProbaOfExerciseDesc << ",0)";
			
			rowDescVec[ProbaOfExercise] = ProbaOfExerciseDesc.str();
			rowTypeVec[ProbaOfExercise] = ARM_STRING;
		}
		
	}

    return ARM_RowInfo(rowDescVec,rowTypeVec);
}

			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///                                                3)    Model Create
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: SetModelKeys
///	Returns: void
///	Action : set the coupon, the funding & the basis reference
///          model names alias
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::SetModelKeys()
{
	// Get the coupon & funding curve ccy
    string cpnCcyName(GetCurrencyUnit()->GetCcyName());
    string fundingCcyName(itsFundingCcy.GetCcyName());
    string basisCcyName(itsBasisCcy.GetCcyName());

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
///	Class  : ARM_CaptionCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : creates the model
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CreateAndSetModel()
{
    /// Get yield curve
	ARM_ZeroCurve* curve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));

	const ARM_DateStripCombiner& dateStructure = DatesStructure();

	/// Creates default values for volatility & mean reversion
    /// (calibration will be called at pricing time) and set them
    /// in the reference model
	std::vector<double> defaultTimes(1,0.0);
	std::vector<double> defaultSigmas(1,SIGMA_DEFAULT_VALUE);
	std::vector<double> defaultBetas(1,BETA_DEFAULT_VALUE);

	/// Get asOfDate
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	bool isArrearsCpn = (itsCpnResetTiming==K_ARREARS);	// In Arrears Case

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
	if (itsBetaCalibFlag != NOCalib)
	{
		paramVector[2] = &newBetaParam;
	}
	else
	{
		paramVector[2] = betaParam;
	}

	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetCpnIndexType();
	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,itsFactorNb,itsSFRMVolType);

	/// Build the default stochastic model of the calculator : SFRM 2F
    ARM_PricingModelPtr refModel( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(new ARM_SFRM( CreateClonedPtr( curve ), *pSFRMModelParams))) );

	// Delte the model params because it is cloned in the model
	delete pSFRMModelParams;

	////////////////////// Create a TreeMethod with a default step number per year and set it
    const ARM_DealDescription dealDescr = GetGenSecurity()->GetDealDescription();
    double lastEventTime = atof(dealDescr.GetElem(dealDescr.GetRowsNb()-1,EventDate).c_str()) - asOfDate;
    int nbSteps=static_cast<int>(floor(TREE_NBSTEPS_PER_YEAR*lastEventTime/K_YEAR_LEN));
	int initialStep =static_cast<int> (nbSteps/10.0);

    if(isTree2G)
    {
        int schedulerType=ARM_SchedulerBase::ConstantVariance;
        std::vector<double> schedulerDatas(3);
        schedulerDatas[0] = nbSteps;
        schedulerDatas[1] = initialStep;
        int samplerType=ARM_SamplerBase::NormalCentred;
        std::vector<double> samplerDatas;
        int truncatorType=ARM_TruncatorBase::StandardDeviation;
        std::vector<double> truncatorDatas(1,STD_DEV_RATIO);
        int reconnectorType=ARM_ReconnectorBase::Mean;
        int smootherType=ARM_SmootherBase::DoNothing;
        ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(1,schedulerType,schedulerDatas,
            samplerType,samplerDatas,truncatorType,truncatorDatas,false,reconnectorType,smootherType);
        refModel->SetNumMethod(ARM_NumMethodPtr( tree ) );
    }
    else
    {
        ARM_TreeMethod* tree = new ARM_TreeMethod(STD_DEV_RATIO, MIN_STD_DEV, initialStep);
        tree->SetNbSteps(nbSteps);
	    refModel->SetNumMethod(ARM_NumMethodPtr( tree ) );
    }
	///////////////////////////////////////////////

    /// Create a Numeraire and set it
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::TerminalZc ) );
    refModel->SetNumeraire(numeraire);
	
    /// Build the model with or without basis effect
    ARM_PricingModelPtr model;
    if(string(GetCurrencyUnit()->GetCcyName()) == string(itsFundingCcy.GetCcyName()))
        // no basis
        model = refModel;
    else
    {
        /// basis => build a Hybrid BFIR
        vector< ARM_PricingModelPtr > models(ARM_HybridBasisFwdIR::NbModels);
        ARM_StringVector names(ARM_HybridBasisFwdIR::NbModels);

        /// Reference model
        models[ARM_HybridBasisFwdIR::RefModel]=refModel;
        names[ARM_HybridBasisFwdIR::RefModel]=GetKeys()[YcKey];

        /// Build the IR forward margin model
	    ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
        models[ARM_HybridBasisFwdIR::IrMarginModel] = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(new ARM_ForwardMarginIR( CreateClonedPtr( fundingCurve ) ) ) );
        names[ARM_HybridBasisFwdIR::IrMarginModel]=GetKeys()[FundingKey];

        /// Build the Basis forward margin model
	    ARM_ZeroCurve* basisCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));
        ARM_HybridBasisFwdIR::modelsAlias refBasisModelIdx(ARM_HybridBasisFwdIR::RefModel);
        if(itsBasisRefModelKey != YcKey)
            refBasisModelIdx = ARM_HybridBasisFwdIR::IrMarginModel;

        models[ARM_HybridBasisFwdIR::BasisMarginModel] = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(new ARM_ForwardMarginBasis( CreateClonedPtr( basisCurve ) ) ) ); 
        names[ARM_HybridBasisFwdIR::BasisMarginModel]=GetKeys()[BasisKey];


        /// Get domestic & foreign models
        string fundingCcy(itsFundingCcy.GetCcyName());
        string basisCcy(itsBasisCcy.GetCcyName());

	    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        string domesticCcy(forex->GetMoneyCurrency()->GetCcyName());
        string foreignCcy(forex->GetMainCurrency()->GetCcyName());

        ARM_HybridBasisFwdIR::modelsAlias domesticModelIdx,foreignModelIdx,modelIdx;
        if(domesticCcy == basisCcy)
        {
            domesticModelIdx = ARM_HybridBasisFwdIR::BasisMarginModel;
            foreignModelIdx  = (foreignCcy == fundingCcy) ? ARM_HybridBasisFwdIR::IrMarginModel
                                                          : ARM_HybridBasisFwdIR::RefModel;
            modelIdx         = foreignModelIdx;
        }
        else
        {
            domesticModelIdx = (domesticCcy == fundingCcy) ? ARM_HybridBasisFwdIR::IrMarginModel
                                                           : ARM_HybridBasisFwdIR::RefModel;
            foreignModelIdx  = ARM_HybridBasisFwdIR::BasisMarginModel;
            modelIdx         = domesticModelIdx;
        }


        /// Build the Forward Forex model
        models[ARM_HybridBasisFwdIR::ForexModel] = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(
            new ARM_ForwardForex(*forex,
            modelIdx == ARM_HybridBasisFwdIR::RefModel ? refModel->GetZeroCurve() : CreateClonedPtr( fundingCurve ),
                CreateClonedPtr( basisCurve) ) ) ); 
        names[ARM_HybridBasisFwdIR::ForexModel]=GetKeys()[ForexKey];


        /// Build the hybrid BFIR
        model = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(new ARM_HybridBasisFwdIR(names,models,
            refBasisModelIdx,domesticModelIdx,foreignModelIdx)) );
    }

	/// Set the model
	SetPricingModel(model);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_CaptionCalculator::GetMRS() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[MrsKey])));
}

void ARM_CaptionCalculator::SetMRS(ARM_ModelParam* mrsParam)
{
    if(!mrsParam || mrsParam->GetType() != ARM_ModelParamType::MeanReversion)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an MRS Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[MrsKey],static_cast< ARM_Object* >(mrsParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetMRS & SetMRS
///	Returns: 
///	Action : get & set the MRS param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_CaptionCalculator::GetBeta() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[BetaKey])));
}

void ARM_CaptionCalculator::SetBeta(ARM_ModelParam* betaParam)
{
    if(!betaParam || betaParam->GetType() != ARM_ModelParamType::Beta)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : an Beta Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[BetaKey],static_cast< ARM_Object* >(betaParam));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetCorrel & SetCorrel
///	Returns: 
///	Action : get & set the Correlation param for SUMMIT interface
/////////////////////////////////////////////////////////////////
const ARM_ModelParam& ARM_CaptionCalculator::GetCorrel() const
{
    return *(static_cast< ARM_ModelParam* >(GetMktDataManager()->GetData(GetKeys()[CorrelKey])));
}

void ARM_CaptionCalculator::SetCorrel(ARM_ModelParam* correlParam)
{
    if(!correlParam || correlParam->GetType() != ARM_ModelParamType::Correlation)
		    ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : a Correl Param is expected !" );

    GetMktDataManager()->RegisterData(GetKeys()[CorrelKey],static_cast< ARM_Object* >(correlParam));
}


			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///                                                4)    Calibration create
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CreateCaptionCapletFloorlet
///	Returns: nothing
///	Action : create the list of the Caption caplet & floorlet
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CreateCaptionCapletFloorletAndDigital()
{
    std::vector<double>& resetDates = itsCalibPorfolioDateStrip->GetResetDates();
	std::vector<double>& startDates = itsCalibPorfolioDateStrip->GetFlowStartDates();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();


	int resetSize=resetDates->size();

	int portfolioSize = resetSize-itsNbdDeedResets;
	
	std::vector<double> initBreakpointTimes(portfolioSize);
	std::vector<double> initVolatility(portfolioSize,SIGMA_DEFAULT_VALUE);	

	list< ARM_Security* > capFloorList;

	for (int i = 0; i < portfolioSize; ++i)
	{
		double resetDate = (*resetDates)[i+itsNbdDeedResets];
		initBreakpointTimes[i]=resetDate-asOfDate;
		ARM_Date startDate((*startDates)[i+itsNbdDeedResets]); 
		ARM_Date endDate(startDate);
		endDate.AddPeriod(itsCpnIdxTerm,itsCpnPayCal.c_str());

		double Strike =  const_cast< ARM_Curve& >(itsStrike).Interpolate(resetDate-asOfDate)*100.0;

		ARM_CapFloor capFloor(startDate,endDate,K_CAP,Strike,itsCpnIndexType,
                    0.0,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

		capFloorList.push_back(static_cast<ARM_CapFloor*>(capFloor.Clone()));
	}

	itsCapFloorPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(capFloorList));
	ARM_SFRM* refModel;
	ARM_HybridBasisFwdIR* bfirModel=NULL;
    for(i=0;i<itsCapFloorPF->size();++i)
	{
        itsCapFloorPF->SetWeight(CF_DEFAULT_WEIGHT,i);
		itsCapFloorPF->SetPrice(CF_DEFAULT_PRICE,i);
	}

	if( string(GetCurrencyUnit()->GetCcyName()) == string(itsFundingCcy.GetCcyName()) )
    {
        refModel = dynamic_cast< ARM_SFRM* >(&*GetPricingModel());
	}
	else
    {
        bfirModel = dynamic_cast< ARM_HybridBasisFwdIR* >(&*GetPricingModel());
        refModel = dynamic_cast< ARM_SFRM* >(&*(bfirModel->GetRefModel()));
	}

	//((ARM_CurveModelParam&) refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility)).SetValuesAndTimes(&initBreakpointTimes,&initVolatility); 
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CreateDiagonalSwaption()
{
	std::vector<double>& resetDates = itsCalibPorfolioDateStrip->GetResetDates();
	std::vector<double>& startDates = itsCalibPorfolioDateStrip->GetFlowStartDates();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	bool isArrearsCpn = (itsCpnResetTiming==K_ARREARS);
	int resetSize = resetDates->size();
	int portfolioSize=((!isArrearsCpn)? resetSize : (resetSize-1));

	std::vector<double> equivStrikes(portfolioSize);
	ComputeEquivalentSwoptStrikes(equivStrikes);
	int i;

	list< ARM_Security* > swaptionList;
	
    ARM_Swaption swaption;

	if(itsNbdDeedResets>=portfolioSize)
    {
		CC_Ostringstream os;
		os << ARM_USERNAME << " : No Swaption To calibrate For this deal ";
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, os.str() );
    }

	portfolioSize = portfolioSize-itsNbdDeedResets;

	for (i = 0; i < portfolioSize; ++i)
	{
		double resetDate = (*resetDates)[i+itsNbdDeedResets];
		ARM_Date startDate((*startDates)[i+itsNbdDeedResets]);
		ARM_Date expiryDate(resetDate);

		ARM_Swap stdSwap(startDate,
                itsEndDate,
                itsCpnIndexType,0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

		double equivStrike = equivStrikes[i+itsNbdDeedResets];

		int RecOrPay = itsPayRec;
		ARM_Swaption swaption((&stdSwap),RecOrPay,K_EUROPEAN,equivStrike,expiryDate);
		
		swaptionList.push_back(static_cast<ARM_Swaption*>(swaption.Clone()));
	}

	itsSwaptionPF = ARM_StdPortfolioPtr(new ARM_StdPortfolio(swaptionList));
    
	
	for(i=0;i<itsSwaptionPF->size();++i)
	{
        itsSwaptionPF->SetWeight(CF_DEFAULT_WEIGHT,i);
		itsSwaptionPF->SetPrice(CF_DEFAULT_PRICE,i);
	}
}
			/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CreateEmptyCalibration
///	Returns: void
///	Action : create an empty calibration structure
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CreateEmptyCalibration()
	
{
	/// Create cap floor portfolio (sigma calibration)
        /// Build a volatility bootstrap calibration on the latter portfolio
    /// The CalibMethod object will be cloned because by default it is not shared

	ARM_CalibMethod* volCalibMethod = NULL;
	
	CreateCaptionCapletFloorletAndDigital();
	volCalibMethod = new ARM_CalibMethod(itsCapFloorPF,ARM_ModelParamVector(),ARM_CalibMethodType::Bootstrap1D);

	
	if(itsSwoptCalibMode!=NOCalib)
    {
        CreateDiagonalSwaption();

		/// Create standard diagonal swaption portfolio (MRS calibration)
        ARM_StdPortfolioPtr swaptionPF(itsSwaptionPF);

        /// Build a MRS optimisation and embed the volatility bootstrapping
        /// The CalibMethod object will be cloned because by default it is not shared

		ARM_CalibMethod* mrsCalibMethod = new ARM_CalibMethod(swaptionPF,ARM_ModelParamVector(),ARM_CalibMethodType::Optimize,
                                 ARM_MAX_ITER,ARM_CalibrationTarget::PriceTarget);
		mrsCalibMethod->SetlinkedMethod(volCalibMethod);
        SetCalibMethod( ARM_CalibMethodPtr( mrsCalibMethod ) );
    }
    else
        SetCalibMethod( ARM_CalibMethodPtr( volCalibMethod ) );
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::CreateAndSetCalibration()
{
	bool isFreezeWeights = false;   // reinit non null weights in the portfolio
    bool isInitSigma = true;        // init sigma param for calibration (bounds & init guess)
    bool isUpdateStrike = true;     // force equivalent strike computation

	bool isArrearsCpn = (itsCpnResetTiming==K_ARREARS);	// In Arrears Case

	CreateEmptyCalibration();

	ComputeCaptionCapletFloorletPrices(
				isFreezeWeights,
				isInitSigma,
				isArrearsCpn,
				SWOPT_At_Exer,
				NOCalib,
				false);
	
	if(itsSwoptCalibMode!=NOCalib)
    {
		ComputeDiagonalSwaptionPrice(isFreezeWeights, isArrearsCpn);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: ComputeCaptionCapletFloorletPrices
///	Returns: nothing
///	Action : compute market target prices of the C/F portfolio
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::ComputeCaptionCapletFloorletPrices(
	bool isFreezeWeights, 
	bool isInitParam,
	bool isArrearCpn,
	CalibrationMode capFloorCalibMode,
	CalibrationMode BetaCalibFlag,
	bool computeSwaptionPrice)
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	ARM_VolCurve* CFVolCurve = CFBSModel->GetVolatility();
    ARM_YCModel* ycModel = CFBSModel->GetYCModel();
    ARM_CalibMethod* volCalibMethod = GetCFCalibMethod();
	std::vector<double>& resetDates = itsCalibPorfolioDateStrip->GetResetDates();
	int resetSize=resetDates->size();

	ARM_StdPortfolioPtr CapFloorPF = itsCapFloorPF;
   
	size_t nbCapFloor =CapFloorPF->GetSize();	
	ARM_CapFloor* capFloor;
	double price,fwdRate,strike,tenor,vol,expiry;

	//size_t NbdDeedResets = resetSize-nbCapFloor;
	std::vector<double> initTimes(nbCapFloor);
	std::vector<double> initSigmas(nbCapFloor);
	
   
	/// Then fills rows one by one
	for(int i =itsNbdDeedResets; i<resetSize; ++i )
    {
		double resetDate = (*resetDates)[i];	
	
		capFloor = static_cast< ARM_CapFloor* >(CapFloorPF->GetAsset(i-itsNbdDeedResets));
	    capFloor->SetModel(CFBSModel);
		price=capFloor->ComputePrice();
		CapFloorPF->SetPrice(price,i-itsNbdDeedResets);
		CapFloorPF->SetPrecision(CF_DEFAULT_PRECISION,i-itsNbdDeedResets);


		expiry = (capFloor->GetExpiryDate().GetJulian()-asOfDate);
		tenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();

		if(fabs(expiry-(resetDate-asOfDate))>K_NEW_DOUBLE_TOL)
		{
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : caplet reset date is not the same than the reset date of the coupon schedule  !" );
		}

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
		vol = CFBSModel->ComputeVol(expiry/K_YEAR_LEN,tenor,fwdRate,strike);		
				
			
		initTimes[i-itsNbdDeedResets]    = expiry; //capFloor->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
		initSigmas[i-itsNbdDeedResets]   = vol/100.0;// (*swapRate; for H&W)
	}

	size_t sigmaIdx,paramSize=volCalibMethod->GetCalibParams().size();
    for(sigmaIdx=0;sigmaIdx<paramSize;++sigmaIdx)
        if(volCalibMethod->GetCalibParam(sigmaIdx) &&
           (volCalibMethod->GetCalibParams())[sigmaIdx]->GetType() == ARM_ModelParamType::Volatility)
            break;	

	std::vector<double> sigmaLowerBound(nbCapFloor,SIGMA_LOWER_BOUND);
    std::vector<double> sigmaUpperBound(nbCapFloor,SIGMA_UPPER_BOUND);
    ARM_CurveModelParam* sigma = new ARM_CurveModelParam(ARM_ModelParamType::Volatility,&initSigmas,&initTimes,
        "SIGMA","STEPUPRIGHT",&sigmaLowerBound,&sigmaUpperBound);
    
	if(sigmaIdx >= paramSize || paramSize == 0)
        volCalibMethod->GetCalibParams().push_back(sigma);
    else
    {
        delete volCalibMethod->GetCalibParam(sigmaIdx);
        (volCalibMethod->GetCalibParams())[sigmaIdx] = sigma;
    }

	SetCFPortfolio((*CapFloorPF));

	//volCalibMethod->GetCalibParams().push_back(sigma);
}

////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: ComputeDiagonalSwaptionPrice
///	Returns: nothing
///	Action : compute market target prices of the swaption portfolio
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::ComputeDiagonalSwaptionPrice(bool isFreezeWeights, bool isArrearCpn)
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );
	ARM_VolCurve* OswVolCurve = oswBSModel->GetVolatility();
    ARM_YCModel* ycModel = oswBSModel->GetYCModel();	

	ARM_CalibMethod* oswCalibMethod = GetSWOPTCalibMethod();
    ARM_StdPortfolioPtr swaptionPF = itsSwaptionPF;
    int i,nbOSW=itsSwaptionPF->GetSize();

	ARM_Swaption* swaption;
	double price;
	
	if(itsSwoptCalibMode==SWOPT_At_Exer) 
	{
		const char* cpnResetCalendar   = itsCpnResetCal.c_str();
		const char* cpnPayCalendar     = itsCpnPayCal.c_str();

		int payFreq     =   itsCpnFreq;
		int payGap      =   GETDEFAULTVALUE;
		int payTiming   =   K_ARREARS;
		
		// principal datestrip
		/*int principalResetTiming =   K_ADVANCE;
		ARM_Date principalEndDate = itsEndDate;
		if (itsCpnResetTiming == K_ARREARS)
		{
			char* ccyName=GetCurrencyUnit()->GetCcyName();
			ARM_Date tmpDate(principalEndDate);
			tmpDate.AddPeriod(itsCpnFreq,ccyName); 
			principalEndDate = tmpDate;
		}
		ARM_DateStrip PrincipalDateStrip(itsStartDate,
										principalEndDate,
										itsCpnFreq,						
										itsCpnDayCount,					
										cpnResetCalendar,			
										K_MOD_FOLLOWING,
										K_ADJUSTED, 
										K_SHORTSTART,						
										-itsNotifDays,						
										payFreq,						
										payGap,						
										cpnPayCalendar,			
										principalResetTiming,					
										payTiming);*/
		
		std::vector<double>& NotifDates = itsPrincipalDateStrip->GetResetDates();
		
		for(i=0;i<nbOSW;++i)
		{
			double notifDate = (*NotifDates)[i+itsNbdDeedResets];	
			double exerciseFees = const_cast< ARM_Curve& >(itsExerStyle).Interpolate(notifDate-asOfDate);
			swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			
			if(fabs((exerciseFees-NON_CALL_FEE))<K_NEW_DOUBLE_TOL)
			{
				swaptionPF->SetPrice(0.0,i);
				swaptionPF->SetPrecision(CF_DEFAULT_PRECISION,i);
				swaptionPF->SetWeight(0.0,i);
			}
			else
			{
				swaption->SetModel(oswBSModel);
				price=swaption->ComputePrice();
				swaptionPF->SetPrice(price,i);
				swaptionPF->SetPrecision(CF_DEFAULT_PRECISION,i);
			}
		}

	}

	else
    {
		for(i=0;i<nbOSW;++i)
		{
			swaption = static_cast< ARM_Swaption* >(swaptionPF->GetAsset(i));
			swaption->SetModel(oswBSModel);
			price=swaption->ComputePrice();
			swaptionPF->SetPrice(price,i);
			swaptionPF->SetPrecision(CF_DEFAULT_PRECISION,i);
		}
	}

	size_t mrsIdx,paramSize=oswCalibMethod->GetCalibParams().size();
    for(mrsIdx=0;mrsIdx<paramSize;++mrsIdx)
    {
        if( oswCalibMethod->GetCalibParam(mrsIdx) &&
            (oswCalibMethod->GetCalibParams())[mrsIdx]->GetType() == ARM_ModelParamType::MeanReversion )
            break;
    }        

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

	if(mrsIdx >= paramSize || paramSize == 0)
		oswCalibMethod->GetCalibParams().push_back(mrs);
    else
    {
        delete oswCalibMethod->GetCalibParam(mrsIdx);
        (oswCalibMethod->GetCalibParams())[mrsIdx] = mrs;
    }

	SetSWOPTPortfolio((*swaptionPF));
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GenerateEquivalentSwoptStrikes
///	Returns: std::vector<double>&
///	Action : generate equivalent swaptionStrikes on fund 
/////////////////////////////////////////////////////////////////

void ARM_CaptionCalculator::ComputeEquivalentSwoptStrikes(std::vector<double>& equivStrikes)
{
	string fundingCcyName(itsFundingCcy.GetCcyName());
    string basisCcyName(itsBasisCcy.GetCcyName());

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* CpnZeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);
	ARM_ZeroCurve* fundZeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[FundingKey]);
	ARM_ZeroCurve* basisZeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[BasisKey]);	

	ARM_ZeroCurve* basisFundCurve = NULL; 
	ARM_ZeroCurve* basisCpnCurve = NULL;
	if(basisCcyName==fundingCcyName)
	{
		basisFundCurve = basisZeroCurve;
		basisCpnCurve = CpnZeroCurve;	
	}
	else
	{
		basisFundCurve = fundZeroCurve;
		basisCpnCurve = basisZeroCurve;	
	}

	ARM_BSModel* CFBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[CfModelKey]) );
	ARM_Y2CModel* fundYCModel = new ARM_Y2CModel(fundZeroCurve,basisFundCurve); 


	std::vector<double>& resetDates = itsCalibPorfolioDateStrip->GetResetDates();
	std::vector<double>& startDates = itsCalibPorfolioDateStrip->GetFlowStartDates();
	std::vector<double>& payDates = itsCalibPorfolioDateStrip->GetPaymentDates();

	std::vector<double>& fundStartDates = itsFundDateStrip->GetFlowStartDates();
	std::vector<double>& fundEndDates = itsFundDateStrip->GetFlowEndDates();
	ARM_INDEX_TYPE fundLiborType = GetFundIndexType();


	//transform CouponSpread in refvalue
	std::vector<double> CouponAbscisses = const_cast< ARM_Curve& >(itsCoupon).GetAbscisses();
	std::vector<double> CouponOrdinates = const_cast< ARM_Curve& >(itsCoupon).GetOrdinates();
	int CouponSize = CouponAbscisses.size();
	ARM_Vector couponRefvalueAbs(CouponSize);
	ARM_Vector couponRefvalueOrd(CouponSize);
	for(int j = 0; j<CouponSize;j++)
	{
		couponRefvalueAbs[j] = CouponAbscisses[j]+asOfDate;	
		couponRefvalueOrd[j] = CouponOrdinates[j]*100.0;
	}
	string interpolatorNameCoupon = const_cast< ARM_Curve& >(itsCoupon).GetInterpolator()->toString();
	int interpolMethodCoupon  = ARM_ArgConv_interpolCurveType.GetNumber(interpolatorNameCoupon);
	ARM_ReferenceValue* refValueCoupon= new ARM_ReferenceValue((&couponRefvalueAbs), (&couponRefvalueOrd));
	refValueCoupon->SetCalcMethod(interpolMethodCoupon);


	//transform FundSpread in refvalue
	std::vector<double> FundSpreadAbscisses = const_cast< ARM_Curve& >(itsFundSpread).GetAbscisses();	
	std::vector<double> FundSpreadOrdinates = const_cast< ARM_Curve& >(itsFundSpread).GetOrdinates();
	int fundSize = FundSpreadAbscisses.size();
	ARM_Vector fundSpreadRefvalueAbs(fundSize);
	ARM_Vector fundSpreadRefvalueOrd(fundSize);
	for(j = 0; j<fundSize;j++)
	{
		fundSpreadRefvalueAbs[j] = FundSpreadAbscisses[j]+asOfDate;	
		fundSpreadRefvalueOrd[j] = FundSpreadOrdinates[j]*100.0;
	}
	string interpolatorNameSpread = const_cast< ARM_Curve& >(itsFundSpread).GetInterpolator()->toString();
	int interpolMethodSpread  = ARM_ArgConv_interpolCurveType.GetNumber(interpolatorNameSpread);
	ARM_ReferenceValue* refValueFundSpread = new ARM_ReferenceValue((&fundSpreadRefvalueAbs), (&fundSpreadRefvalueOrd));
	refValueFundSpread->SetCalcMethod(interpolMethodSpread);
	
	//transform notional in refvalue
	std::vector<double> notionalAbscisses = const_cast< ARM_Curve& >(itsNotional).GetAbscisses();
	std::vector<double> notionalOrdinates = const_cast< ARM_Curve& >(itsNotional).GetOrdinates();
	int notionalSize = notionalAbscisses.size();
	ARM_Vector notionalRefvalueAbs(notionalSize);
	ARM_Vector notionalRefvalueOrd(notionalSize);
	for(j = 0; j<notionalSize;j++)
	{
		notionalRefvalueAbs[j] = notionalAbscisses[j]+asOfDate;	
		notionalRefvalueOrd[j] = notionalOrdinates[j];
	}
	string interpolatorNameNotional = const_cast< ARM_Curve& >(itsNotional).GetInterpolator()->toString();
	int interpolMethodNotional  = K_STEPUP_RIGHT;
	ARM_ReferenceValue* refValueNotional = new ARM_ReferenceValue((&notionalRefvalueAbs), (&notionalRefvalueOrd));
	refValueNotional->SetCalcMethod(interpolMethodNotional);


	ARM_Currency* cpnCcy = GetCurrencyUnit();
	int size = equivStrikes.size();
	
	ARM_FixLeg couponFixLeg(itsStartDate,
                       itsEndDate, 
                       refValueCoupon,
                       K_RCV, 
					   itsCpnFreq,
                       itsCpnDayCount,
                       K_COMP_PROP,
                       K_ARREARS, 
					   K_ADJUSTED,
                       K_SHORTSTART,
                       cpnCcy);

	couponFixLeg.SetAmount(refValueNotional);

	ARM_FixLeg fundFixLeg(itsStartDate,
                       itsEndDate, 
                       refValueFundSpread,
                       K_RCV, 
					   itsFundFreq,
                       itsFundDayCount,
                       K_COMP_PROP,
                       K_ARREARS, 
					   K_ADJUSTED,
                       K_SHORTSTART,
                       &itsFundingCcy);

	fundFixLeg.SetAmount(refValueNotional);

	ARM_SwapLeg fundSwapLeg(itsStartDate,
                       itsEndDate, 
                       fundLiborType,
                       K_RCV, 
					   0,
					   itsFundFreq,
					   itsFundFreq,
					   K_ADVANCE,
					   K_ARREARS,
					   &itsFundingCcy);
	
	fundSwapLeg.SetDayCount(itsFundDayCount);
	fundSwapLeg.SetAmount(refValueNotional);

	ARM_Swap stdCpnSwap(itsStartDate,
                itsEndDate,
                itsCpnIndexType,0.0,1.0,K_RCV,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	//datestrip for notif dates
	const char* cpnResetCalendar   = itsCpnResetCal.c_str();
		const char* cpnPayCalendar     = itsCpnPayCal.c_str();

		int payFreq     =   itsCpnFreq;
		int payGap      =   GETDEFAULTVALUE;
		int payTiming   =   K_ARREARS;
		
		// principal datestrip
		int principalResetTiming =   K_ADVANCE;
		ARM_Date principalEndDate = itsEndDate;
		if (itsCpnResetTiming == K_ARREARS)
		{
			char* ccyName=GetCurrencyUnit()->GetCcyName();
			ARM_Date tmpDate(principalEndDate);
			tmpDate.AddPeriod(itsCpnFreq,ccyName); 
			principalEndDate = tmpDate;
		}

		ARM_DateStrip principalDateStrip(itsStartDate,
										principalEndDate,
										itsCpnFreq,						
										itsCpnDayCount,					
										cpnResetCalendar,			
										K_MOD_FOLLOWING,
										K_ADJUSTED, 
										K_SHORTSTART,						
										-itsNotifDays,						
										payFreq,						
										payGap,						
										cpnPayCalendar,			
										principalResetTiming,					
										payTiming);


	double notionalExchange = 0.0;

	for(int i = size-1; i>=0; i--)
	{
		//cpn fixleg At Coupon
		double startDate = (*startDates)[i];
		ARM_Date cpnStartdate(startDate);

		couponFixLeg.SetStartDate(cpnStartdate);
        couponFixLeg.CptCashFlowDates();
		couponFixLeg.SetModelVariable(NULL);
        couponFixLeg.SetModel(CFBSModel);
		double couponlegNPV = couponFixLeg.ComputePrice();

		//fund fixleg At spread
		double dStartDate = (*fundStartDates)[i];
		ARM_Date fundStartDate(dStartDate);

		fundFixLeg.SetStartDate(fundStartDate);
        fundFixLeg.CptCashFlowDates();
		fundFixLeg.SetModelVariable(NULL);
        fundFixLeg.SetModel(fundYCModel);
		double fundspreadLegNPV = fundFixLeg.ComputePrice();

		//Fees
		//double resetDate = (*resetDates)[i];		
		double eventDate=(*(principalDateStrip.GetResetDates()))[i];
		double fees = const_cast< ARM_Curve& >(itsExerStyle).Interpolate(eventDate-asOfDate);
		double paydate = (*payDates)[i];
		double FeesNPV =0.0;
		if(fabs((fees-NON_CALL_FEE))>K_NEW_DOUBLE_TOL)
		{
			FeesNPV = fees*basisCpnCurve->DiscountPrice((paydate-asOfDate)/K_YEAR_LEN);
		}

		//notionalExchange
		double dEndDate = (*fundEndDates)[i];
		double notional_i = const_cast< ARM_Curve& >(itsNotional).Interpolate(paydate-asOfDate);
		double notionalexchangeFragment = -(basisFundCurve->DiscountPrice((dStartDate-asOfDate)/K_YEAR_LEN)-basisFundCurve->DiscountPrice((dEndDate-asOfDate)/K_YEAR_LEN)); 		
		notionalexchangeFragment += (basisCpnCurve->DiscountPrice((dStartDate-asOfDate)/K_YEAR_LEN)-basisCpnCurve->DiscountPrice((dEndDate-asOfDate)/K_YEAR_LEN));
		notionalexchangeFragment*=notional_i;
		notionalExchange+= notionalexchangeFragment;

		//total NPV
		int positionOnFees = -itsPayRec;
		double totalNPV = couponlegNPV-fundspreadLegNPV-notionalExchange+positionOnFees*FeesNPV; 

		//fundFloat
		fundSwapLeg.SetStartDate(fundStartDate);
        fundSwapLeg.CptCashFlowDates();
		fundSwapLeg.SetModelVariable(NULL);
        fundSwapLeg.SetModel(fundYCModel);
		double fundLegNPV = fundSwapLeg.ComputePrice();

		//std swap rate 
		(stdCpnSwap.Get1stLeg())->SetStartDate(cpnStartdate);
		(stdCpnSwap.Get2ndLeg())->SetStartDate(cpnStartdate);
		(stdCpnSwap.Get1stLeg())->CptCashFlowDates();
		(stdCpnSwap.Get2ndLeg())->CptCashFlowDates();
		stdCpnSwap.CptCashFlowDates();
		stdCpnSwap.SetModelVariable(NULL);
        stdCpnSwap.SetModel(CFBSModel);
		double stdSwapNPV = stdCpnSwap.ComputePrice();
		double stdSwapRate = stdCpnSwap.CptMarketSwapRate();

		equivStrikes[i]=totalNPV*stdSwapRate/fundLegNPV; 		
	}
	delete fundYCModel;
	/*delete refValueCoupon;
	delete refValueNotional;
	delete refValueFundSpread;*/

		/*delete refValueFundSpread;
		delete refValueNotional;
		delete refValueCoupon;*/
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetCFCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for cap floor
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CaptionCalculator::GetCFCalibMethod() const
{
    
#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        ((itsSwoptCalibMode!=NOCalib) && GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

    if(itsSwoptCalibMode!=NOCalib)
        return GetCalibMethod()->GetlinkedMethod();
    else
        return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetOSWCalibMethod
///	Returns: ARM_CalibMethod
///	Action : get calibration method for diagonal swaptions
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_CaptionCalculator::GetSWOPTCalibMethod() const
{
	if(itsSwoptCalibMode == NOCalib)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No method available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL))
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method not found");
#endif

    return &(*GetCalibMethod());
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetCFPortfolio
///	Returns: ARM_Portfolio
///	Action : get the cap floor calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CaptionCalculator::GetCFPortfolio() const
{
#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        ((itsSwoptCalibMode != NOCalib) && (GetCalibMethod()->GetlinkedMethod() == NULL ||
            GetCalibMethod()->GetlinkedMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL))) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method or portfolio not found");
#endif
    if(itsSwoptCalibMode != NOCalib)
        return GetCalibMethod()->GetlinkedMethod()->GetPortfolio();
    else
        return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: SetCFPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of cap floors
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::SetCFPortfolio(const ARM_StdPortfolio& port)
{
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
        ((itsSwoptCalibMode != NOCalib) && GetCalibMethod()->GetlinkedMethod() == NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Bootstrap calib method not found");
#endif

    /// Update portfolio
    if(itsSwoptCalibMode != NOCalib)
        GetCalibMethod()->GetlinkedMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
    else
        GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(port).Clone())) );
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_Portfolio
///	Action : get the diagonal swaptions calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_CaptionCalculator::GetSWOPTPortfolio() const
{

    if(itsSwoptCalibMode == NOCalib)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : No portfolio available because calibration is off");

#ifdef __GP_STRICT_VALIDATION
    if( GetCalibMethod() == ARM_CalibMethodPtr(NULL) ||
        GetCalibMethod()->GetPortfolio() == ARM_StdPortfolioPtr(NULL) )
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME+" : Optimisation calib method or portfolio not found");
#endif

    return GetCalibMethod()->GetPortfolio();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: SetOSWPortfolio
///	Returns: void
///	Action : affectation of a calibration portfolio made of digonal 
///			 swaptions.
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::SetSWOPTPortfolio(const ARM_StdPortfolio& port)
{
    if(itsSwoptCalibMode == NOCalib)
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
///	Class  : ARM_CaptionCalculator
///	Routine: GetSubPrice
///	Returns: double
///	Action : compute the analytical price of a sub selection of
///          the deal description
/////////////////////////////////////////////////////////////////
double ARM_CaptionCalculator::GetSubPrice(int startRowIdx,int endRowIdx,CaptionColAlias columnName,
                                      const string& evalDateStr,
                                      const ARM_DealDescription& dealDesc)
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


////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Calibrate
///	Returns: void
///	Action : Calibrate the model of the calculator.
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::Calibrate()
{
	/// Volatility bootstrapping with possible MRS optimisation
    /// FIX FIX : ugly test waiting for the generic calibration extended for hybrid models
    /// with an equivalent modelMap on modelParam
    ARM_HybridBasisFwdIR* bfirModel = dynamic_cast< ARM_HybridBasisFwdIR* >(&*GetPricingModel());
    if(bfirModel)
        /// Delegate to the reference model
		GetCalibMethod()->Calibrate(&*(bfirModel->GetRefModel()));
    else
		GetCalibMethod()->Calibrate(&*GetPricingModel());   
}


///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///                                                5)    Pricing
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------
///------------------------------------------------------------------------------------------------------------------------------------

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Set???ToPrice() and Is???ToPrice()
///	Returns: void/boolean
///	Action : Set the flag to know which product to price.
///          Test the product currently able to be priced
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::SetCaptionToPrice(bool toPrice)			{itsProductsToPrice[CaptionPrice]				=   toPrice;}
bool ARM_CaptionCalculator::IsCaptionToPrice() const					{return itsProductsToPrice[CaptionPrice];}

void ARM_CaptionCalculator::SetCapToPrice(bool toPrice)				{itsProductsToPrice[CapPrice]				=   toPrice;}
bool ARM_CaptionCalculator::IsCapToPrice() const					{return itsProductsToPrice[CapPrice];}

void ARM_CaptionCalculator::SetCouponLegToPrice(bool toPrice)		{itsProductsToPrice[CouponLegPrice]	=   toPrice;}
bool ARM_CaptionCalculator::IsCouponLegToPrice() const				{return itsProductsToPrice[CouponLegPrice];}

void ARM_CaptionCalculator::SetFundingLegToPrice(bool toPrice)		{itsProductsToPrice[FundingLegPrice]		=   toPrice;}
bool ARM_CaptionCalculator::IsFundingLegToPrice() const				{return itsProductsToPrice[FundingLegPrice];}

void ARM_CaptionCalculator::SetCaptionStrikesToPrice(bool toPrice)	{itsProductsToPrice[CaptionStrikes]		=	toPrice;}
bool ARM_CaptionCalculator::IsCaptionStrikesToPrice() const		{return itsProductsToPrice[CaptionStrikes];}

void ARM_CaptionCalculator::SetExerciseProbasToPrice(bool toPrice)		{itsProductsToPrice[ExerciseProbas]			=	toPrice;}
bool ARM_CaptionCalculator::IsExerciseProbasToPrice() const				{return itsProductsToPrice[ExerciseProbas];}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the Caption deal.
/////////////////////////////////////////////////////////////////
double ARM_CaptionCalculator::Price()
{
	CalibrateAndTimeIt();

    /// Price the implicit product according to internal flag
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	CC_Ostringstream evalDateDesc;
    evalDateDesc << CC_NS(std,fixed) << asOfDate.GetJulian();
    string evalDateStr(evalDateDesc.str());

	ARM_ZeroCurve* zeroCurve = (ARM_ZeroCurve*) GetMktDataManager()->GetData(GetKeys()[YcKey]);

    for (size_t i = 0; i < NbProductsToPrice; ++i)
	{
		if (itsProductsToPrice[i])
		{
			/// Select the right column that describes the product
			CaptionColAlias prodIdx;

			if(i == CaptionPrice)
			{
				prodIdx = Caption;				
			}
			
			else if(i == CapPrice)
				prodIdx = Cap;

			else if(i == CouponLegPrice)
				prodIdx = CpnNPV;

			else if(i == CaptionStrikes)
				prodIdx = CaptionStrike;

			else if(i == FundingLegPrice)
				prodIdx = FundNPV;

			else if(i == ExerciseProbas)
				prodIdx = ProbaOfExercise;
			
			ARM_GenSecurityPtr genSec = GetGenSecurity();
			ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();

			/// Build the sub-generic security
			size_t startRowIdx = 1; /// start with the first actual row (=> skip column names)
			size_t endRowIdx = dealDesc.GetRowsNb() - 1;
			ARM_DealDescriptionPtr subDealDesc = dealDesc.GetSubDescription(startRowIdx,endRowIdx,prodIdx+1);
			ARM_GenSecurityPtr subGenSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,genSec->GetPayModelName()));

			//genSec = ARM_GenSecurityPtr(new ARM_GenSecurity(subDealDesc,genSec->GetPayModelName()));

			ARM_GenPricer* genPricer = new ARM_GenPricer( &*subGenSec,&*GetPricingModel());

			ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

			double price = genPricer->Price();
			
			if(i == CaptionPrice)
			{
				itsCaptionPrice = price;				
			}
			else if(i == CapPrice)
			{
				ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");

				ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

				ARM_VectorPtr  newItermPrices(new std::vector<double>(itermPrices->size()));

				for (size_t j =0; j < itermPrices->size(); ++j)
				{
					double df = 1.0;
					
					(*newItermPrices)[j] = (*itermPrices)[j];
				}
				
				itermPrices = newItermPrices;

				itsCapPrice = itermPrices;
			}
			else if(i == CouponLegPrice)
			{
				ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");

				ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

				ARM_VectorPtr  newItermPrices(new std::vector<double>(itermPrices->size()));

				for (size_t j =0; j < itermPrices->size(); ++j)
				{
					double df = 1.0;
					
					(*newItermPrices)[j] = (*itermPrices)[j];
				}
				
				itermPrices = newItermPrices;

				itsCouponLegPrice = itermPrices;
			}
			else if(i == CaptionStrikes)
			{
				ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");

				ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

				ARM_VectorPtr  newItermPrices(new std::vector<double>(itermPrices->size()));

				for (size_t j =0; j < itermPrices->size(); ++j)
				{
					double df = 1.0;
					
					(*newItermPrices)[j] = (*itermPrices)[j];
				}
				
				itermPrices = newItermPrices;

				itsCaptionStrikes = itermPrices;
			}
			else if(i == FundingLegPrice)
			{
				ARM_GramFctorArg gramFunctorArg = genPricer->GetPricerInfo()->GetContents().GetData("IntermediatePrices");

				ARM_VectorPtr itermPrices = gramFunctorArg.GetVector();

				ARM_VectorPtr  newItermPrices(new std::vector<double>(itermPrices->size()));

				for (size_t j =0; j < itermPrices->size(); ++j)
				{
					double df = 1.0;
					
					(*newItermPrices)[j] = (*itermPrices)[j];
				}
				
				itermPrices = newItermPrices;

				itsFundingLegPrice = itermPrices;
			}
			else if(i == ExerciseProbas)
			{
				double numeraireDate = (*(itsCpnDateStrip->GetFwdRateEndDates()))[(*(itsCpnDateStrip->GetFwdRateEndDates())).size()-1];

				double numeraireDF = zeroCurve->DiscountPrice((numeraireDate-asOfDate.GetJulian())/K_YEAR_LEN);
				itsExerciseProbas = price/numeraireDF;
			}
		}
	}

	itsHasBeenPriced = true;
    
    return itsCaptionPrice*GetPorS();
}

			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///                                                6)    Get data & Update
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------
			///------------------------------------------------------------------------------------------------------------------------------------


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: ComputePricingData
///	Returns: a ARM_MultiTypeDict
///	Action : get pricing data of the calculaltor
/////////////////////////////////////////////////////////////////

void ARM_CaptionCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_CaptionCalculator*>(this)->PriceAndTimeIt();

	if (itsProductsToPrice[CaptionPrice])
		GetPricingData()[ "CaptionPrice" ] = itsCaptionPrice;
	
	if (itsProductsToPrice[CapPrice])
		GetPricingData()[ "CapPrice" ] = itsCapPrice;		

	if (itsProductsToPrice[CouponLegPrice])	
		GetPricingData()[ "CouponLegPrice" ] = itsCouponLegPrice;		
	
	if (itsProductsToPrice[FundingLegPrice])
		GetPricingData()[ "FundingLegPrice"	] = itsFundingLegPrice;	
	
	if (itsProductsToPrice[CaptionStrikes])
		GetPricingData()[ "CaptionStrikes" ] = itsCaptionStrikes;

	if (itsProductsToPrice[ExerciseProbas])
		GetPricingData()[ "ExerciseProbas" ] = itsExerciseProbas;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: GetRefModel
///	Returns: ARM_PricingModelPtr
///	Action : return the reference model of the Caption
/////////////////////////////////////////////////////////////////
ARM_SFRM* ARM_CaptionCalculator::GetRefModel()
{
	ARM_SFRM* refModel;
    ARM_HybridBasisFwdIR* bfirModel=NULL;
    if( string(GetCurrencyUnit()->GetCcyName()) == string(itsFundingCcy.GetCcyName()) )
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
///	Class  : ARM_CaptionCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::UpdateModel()
{
	const ARM_DateStripCombiner& dateStructure = DatesStructure();

	std::vector<double>& resetDates = dateStructure.GetDateStrip(0)->GetResetDates();

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
    if( string(GetCurrencyUnit()->GetCcyName()) == string(itsFundingCcy.GetCcyName()) )
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

	ARM_ZeroCurve* cpnCurve     = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcKey]));
	ARM_ZeroCurve* fundingCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[FundingKey]));
	ARM_ZeroCurve* basisCurve   = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[BasisKey]));

    /// Update yield curves
    if(bfirModel)
    {
	    ARM_Forex* forex   = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
        bfirModel->UpdateCurves(CreateClonedPtr( cpnCurve ),CreateClonedPtr( fundingCurve ),CreateClonedPtr( basisCurve ),*forex);
    }
    else
	    refModel->SetZeroCurve( CreateClonedPtr( cpnCurve ) );

	ARM_ModelParamVector paramVector(4);

	paramVector[0] = (&refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility));

    /// Update the MRS
	paramVector[1] = mrsParam;
	/// Update the Beta
	if (itsBetaCalibFlag != NOCalib)
	{
		paramVector[2] = (&refModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Beta));
	}
	else
	{
		paramVector[2] = betaParam;
	}
	paramVector[3] = correlParam;

	ARM_INDEX_TYPE liborType = GetCpnIndexType();

	ARM_IRIndex IRIndex(liborType,K_DEF_FREQ,K_DEF_FREQ,GetCurrencyUnit());

	/// use the factory class to get the model params!
	ARM_ModelParamsSFRM* pSFRMModelParams = 
		ARM_ModelParamsSFRMFactory.Instance()->CreateModelParamsSFRM(paramVector,&IRIndex,itsFactorNb,SFRM_VOL_TYPE);

	refModel->SetModelParams(*pSFRMModelParams);
	delete pSFRMModelParams;

	ARM_Portfolio* port=NULL;
	refModel->ConvertToShiftorBetaParam(*port);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_CaptionCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_CaptionCalculator::UpdateCalibration(bool isUpdateStrike)
{/// Get the current market model for swaption and its associated YC model
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswModelKey]) );

    /// Update cap floor prices
    bool isFreezeWeights = true; // keep the portfolio size to avoid hedge jumps
    bool isInitSigma = true;       // keep current sigma param to init calibration
    bool isArrearsCpn = (itsCpnResetTiming==K_ARREARS);	// In Arrears Case

    ComputeCaptionCapletFloorletPrices(
			isFreezeWeights,
			isInitSigma,
			isArrearsCpn,
			itsSwoptCalibMode,
			itsBetaCalibFlag,
			isInitSigma);

	/// Compute diagonal swaption target prices
    if(itsSwoptCalibMode != NOCalib)
    {
        ComputeDiagonalSwaptionPrice(isFreezeWeights,isArrearsCpn);
    }

}


////////////////////////////////////////////////////
///	Class   : ARM_CaptionCalculator
///	Routines: View
///	Returns :
///	Action  : Standard ARM object support
////////////////////////////////////////////////////

void ARM_CaptionCalculator::View(char* id, FILE* ficOut) const
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

    /// Caption Calculator specific datas viewing
    fprintf(fOut,"\n\n =======> Caption CALCULATOR <====== \n");

    CC_Ostringstream CaptionData;
    CaptionData << "\nStartDate = " <<  itsStartDate.toString() << "\n";
    CaptionData << "EndDate = " << itsEndDate.toString() << "\n";
	CaptionData << "Notification days = " << itsNotifDays << "\n";
	CaptionData << "Pay/Rec option = " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, ((ARM_GenCalculator*)this)->GetPorS()) << "\n";
	    
    CaptionData << "\n Cpn Leg Datas :\n";
	CaptionData << "Cpn Index Term = " << itsCpnIdxTerm << "\n";
	CaptionData << "Pay/Rec = " << ARM_ParamView::GetMappingName(S_RECEIVE_PAY, itsPayRec) << "\n";
	CaptionData << "Cap/Floor = " << ARM_ParamView::GetMappingName(S_CAPLIKE_TYPE, itsCapOrFloor) << "\n";
    CaptionData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsCpnDayCount ) << "\n";
	CaptionData << "Reset Timing = " << ARM_ParamView::GetMappingName(S_TIMING_MOD, itsCpnResetTiming ) << "\n";
    CaptionData << "Reset Calendar = " << itsCpnResetCal << "\n";
    CaptionData << "Pay Calendar = " << itsCpnPayCal << "\n";     
	
    CaptionData << "\nFunding Leg Datas :\n";
	CaptionData << "Fund Index Term = " << itsFundIdxTerm << "\n";
    CaptionData << "Day Count = " << ARM_ParamView::GetMappingName( S_DAYCOUNT,  itsFundDayCount ) << "\n";
	CaptionData << "Reset Calendar = " << itsFundResetCal << "\n";
    CaptionData << "Pay Calendar = " << itsFundPayCal << "\n";    
    
    if( GetGenSecurity() != ARM_GenSecurityPtr(NULL))
    {
        const ARM_DealDescription dealDesc = GetGenSecurity()->GetDealDescription();
        size_t i,nbRows = dealDesc.GetRowsNb();

        CaptionData << "\nCoupon Leg Profile :\n";
        CaptionData << "StartDate   \tStrike    \tSpread  \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,PayDate) != ARM_MISSING_TYPE)
            {
                CaptionData << ARM_Date(atof(dealDesc.GetElem(i,EventDate).c_str())).toString() << "\t";
                //CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Strike).c_str());
                if(dealDesc.GetElemFormat(i,Strike) != ARM_MISSING_TYPE)
                {
                    CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Strike).c_str()) << "\t";
                    CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,CpnSpread).c_str()) << "\n";
                    
                }
                else
                    CaptionData << "\n";
            }
        }

        CaptionData << "\nFunding Spread Profile :\n";
        CaptionData << "StartDate \tSpread   \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,FundStartDate) != ARM_MISSING_TYPE)
            {
                CaptionData << ARM_Date(atof(dealDesc.GetElem(i,EventDate).c_str())).toString() << "\t";
                CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,FundSpread).c_str()) << "\n";
            }
        }

        CaptionData << "\nNominal Profile :\n";
        CaptionData << "PayDate   \tNominal   \n";
        for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,PayDate) != ARM_MISSING_TYPE)
            {
                CaptionData << ARM_Date(atof(dealDesc.GetElem(i,PayDate).c_str())).toString() << "\t";
                CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Notional).c_str()) << "\n";
            }
        }

		CaptionData << "\nExerciseFees Profile :\n";
        CaptionData << "EventDate   \tPayDate   \tFees   \n";
		for(i=1;i<nbRows;++i)
        {
            if(dealDesc.GetElemFormat(i,ExerciseFees) != ARM_MISSING_TYPE)
            {
                CaptionData << ARM_Date(atof(dealDesc.GetElem(i,EventDate).c_str())).toString() << "\t";
				CaptionData << ARM_Date(atof(dealDesc.GetElem(i,FeesPayDate).c_str())).toString() << "\t";
                CaptionData << CC_NS(std,fixed) << atof(dealDesc.GetElem(i,Fees).c_str()) << "\n";
            }
        }      
    }

	/// part common to gencalculator
	CaptionData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString() << "\n\n";
    fprintf(fOut,"%s",CaptionData.str().c_str());

	string cfTxt("\n\n Cap Floor Calibration For sigma Curve\n");
    fprintf(fOut,"%s\n",cfTxt.c_str());
      
    ARM_StdPortfolioPtr capFloorPF = GetCFPortfolio();
        if(capFloorPF!=ARM_StdPortfolioPtr(NULL))
            capFloorPF->View(id,fOut);

	
	
	string oswTxt("\n\nDiagonal swaption portfolio (MRS calibration)\n");
    oswTxt += "Calibration=";
	oswTxt += ARM_ArgConvReverse_CaptionCalibMode.GetString(itsSwoptCalibMode);

    fprintf(fOut,"%s\n",oswTxt.c_str());
    if(itsSwoptCalibMode != NOCalib)
    {
        ARM_StdPortfolioPtr diagonalSwaptionPF=GetSWOPTPortfolio();
        if(diagonalSwaptionPF!=ARM_StdPortfolioPtr(NULL))
            diagonalSwaptionPF->View(id,fOut);
    }

    if ( ficOut == NULL )
       fclose(fOut);
}

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

