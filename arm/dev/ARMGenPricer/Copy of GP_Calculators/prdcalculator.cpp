/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prdccalculator.cpp
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  JP.Prie & E.Ezzine
 *	\version 1.0
 *	\date February 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/prdcalculator.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/utilityport.h"  
#include "gpbase/curveconvert.h"  
#include "gpbase/gplinalgconvert.h"
#include "gpbase/datestripconvert.h"
#include "gpbase/globalconstant.h"
#include "gpbase/surface.h"


#include "gpcalib/modelfitterdes.h"

/// gpinfra
#include "gpinfra/mktdatamanagerrep.h"
#include "gpinfra/argconvdefault.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/gensecurity.h"
#include "gpinfra/dealdescription.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpinfra/genpricer.h"
#include "gpinfra/numerairefactory.h"
#include "gpinfra/curvemodelparam.h"
#include "gpinfra/pricerinfo.h"
#include "gpinfra/cstmanager.h"
#include "gpinfra/discretisationscheme.h"
#include "gpinfra/surfacemodelparam.h"

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/VanillaSwaption.h"
#include "gpcalib/VanillaFxOption.h"
#include "gpcalib/VanillaIrFxSwaption.h"
#include "gpcalib/KernelToGP.h"
#include "gpcalib/stripper.h"

/// gpmodels
#include "gpmodels/2irfxModel.h"
#include "gpmodels/q1f.h"
#include "gpmodels/q1f_fx.h"
#include "gpmodels/modelparamsq1f.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"
#include "gpmodels/local_sln_model.h"
#include "gpmodels/local_sln_modelparams.h"
#include "gpmodels/MarketIRModel.h"
#include "gpmodels/MarketHybridModel.h"
#include "gpmodels/BS_ModelParams.h"
#include "gpmodels/BS_Model.h"

/// gpnummethods
#include "gpnummethods/treebase.h"
#include "gpnummethods/treefactory.h"
#include "gpnummethods/scheduler.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/truncator.h"
#include "gpnummethods/reconnector.h"
#include "gpnummethods/smoother.h"

/// gpclosedforms
#include "gpclosedforms/vanilla_bs.h"

/// kernel
#include <crv/volflat.h>
#include <inst/powrev.h>
#include <mod/xbsfx.h>
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/fixleg.h>
#include <inst/swapleg.h>
#include <inst/swap.h>
#include <util/fromto.h>
#include <inst/optionportfolio.h>


/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

/// BoosterData column names
const unsigned int NOTICE_DATE              = 0;
const unsigned int FX_RESET_DATE            = 1;
const unsigned int FX_PAY_DATE              = 2;
const unsigned int FUNDING_START_DATE       = 3;
const unsigned int FUNDING_END_DATE         = 4;
const unsigned int FUNDING_PAY_DATE         = 5;
const unsigned int FUNDING_SPREAD           = 6;
const unsigned int FX_NOTIO_MULTIPLIER      = 7;
const unsigned int FX_STRIKE                = 8;
const unsigned int FX_CAP                   = 9;
const unsigned int FX_ACCRUAL_BASIS         = 10;
const unsigned int FUNDING_ACCRUAL_BASIS    = 11;
const unsigned int COUPON_DOMESTIC          = 12;

const double FOREX_CPN = -1.0;
const double FIXED_CPN = +1.0;

const double NO_FLOOR				= -0.99;
const double NO_CAP					=  99;
const double NON_CALL_FEE			= 1.0e15;

/// Strike maximum value to allow call pricing
const double MAX_FX_STRIKE = 5000.0;


/// RedemptionData names
const unsigned int REDEMPTION_PAY_DATE      = 2;


/// H&W vol range [10bp,500bp]
const double HWVOL_LOWER_BOUND      = 0.001;
const double HWVOL_UPPER_BOUND      = 0.05;

/// Q vol range [2.5%,100%]
const double QVOL_LOWER_BOUND       = 0.025;
const double QVOL_UPPER_BOUND       = 1.0;

/// H&W MRS range [-15%,50%] with a 0% default value
const double MRS_LOWER_BOUND        = -0.15;
const double MRS_UPPER_BOUND        = 0.5;
const double MRS_DEFAULT_VALUE      = 0.0;

/// 10-3 bp of vega to be selected in portfolio for IR volatility bootstrapping 
const double IR_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double OSW_DEFAULT_WEIGHT     = 1.0;
const double OSW_DEFAULT_PRICE      = 1.0e+10;

/// 10-3 bp of vega to be selected in portfolio for FX volatility bootstrapping 
const double FX_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double FX_DEFAULT_WEIGHT      = 1.0;
const double FX_DEFAULT_PRICE       = 1.0e+30;

/// For ATS minimum volatility search
const double FX_MIN_MONEYNESS       = 0.5;
const double FX_MAX_MONEYNESS       = 1.75;
const double FX_STEP_MONEYNESS      = 0.05;

/// Basket IR/FX constants
const size_t VNS_STRIKE_MAX_ITER	= 100;
const double VNS_STRIKE_TOL			= 0.00001; // 0.1bp
const double VNS_SWAP_BUMP			= 0.00001; // 0.1bp

/// Default MDM key names
const string YC_KEY_NAME                    = "YC_";
const string YC_BASIS_KEY_NAME              = "YC_BASIS_";
const string FOREX_KEY_NAME                 = "FOREX_";
const string OSWMODEL_KEY_NAME              = "OSWMOD_";
const string FXMODEL_KEY_NAME               = "FXMOD_";
const string CORREL_KEY_NAME                = "CORREL_";
const string MRS_KEY_NAME                   = "MRS_";
const string Q_KEY_NAME                     = "Q_";
const string FLOORED_FOREX_KEY_NAME         = "FLOORED_FOREX_";
const string CAPPED_FOREX_KEY_NAME          = "CAPPED_FOREX_";
const string REDEMPTION_FOREX_KEY_NAME      = "REDEMPTION_FOREX_";
const string UNKNOWN_KEY_NAME               = "UNKNOWN";
const string LOCAL_FXMODEL_KEY_NAME         = "LOCAL_FXMOD_";
const string MARKET_IRMODEL_KEY_NAME		= "MARKET_IRMOD_";
/// ugly ugly, we have to fix fix fix 
const string BSFXVol_KEY_NAME				= "BSFXVol_";

const string ARM_PRDCalculator::PRDColNamesTable [] =
{
    "EventDate",
    "FxResetDate", 
	"StartDate",
	"FundingEndDate",
	"FundingPayDate",
	"EndDate",
    "CpnPayDate", 
	"CpnIT",
	"CpnLeverage",
	"FxLowStrike",
    "FxHighStrike",
	"FundingSpread",    
    "FundingLeg",    
	"FxLowCall",
	"FxHighCall",
    "FixedLeg",
    "FxLowCallStrip",
    "FxHighCallStrip",
    "RedemptionResetDate",
    "RedemptionPayDate",
    "RedemptionStrike",
	"Nominal",
    "Redemption",
    "CpnLeg",
    "PRDSwap",
	"FundingFlow",
	"FixFlow",
    "FxStrip",
	"PRDCFirstSwap",
    "PRDCFirstEuropean",

};

const string ARM_PRDCalculator::PRDProfileNamesTable [] =
{
	"FxDateStripProfile",
	"CpnDateStripProfile",
	"FundingDateStripProfile",
	"NominalProfile",
	"FxLowStrikeProfile",
    "FxHighStrikeProfile",
	"FxLeverageProfile",
    "MarginProfile",
	"MarginTimesNotional",
	"FundLvgeTimesNotional",
	"CpnMinTimesNotional",
};


/// Maximum calibration errors in % (IR) and real terms (FX)
const double MAX_CALIB_ERR[] = {    0.01,
                                    0.01,
                                    0.0001};

const string CALIB_ERR_MSGE[] = {   ": domestic IR calibration failed",
                                    ": foreign IR calibration failed",
                                    ": forex calibration failed",
                                    ": hybrid IR/FX calibration failed"};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_PRDCalculator::ARM_PRDCalculator( const ARM_PRDCalculator& rhs )
:	ARM_GenCalculator( rhs ),
	
	itsStartDate(rhs.itsStartDate), 
	itsFixedEndDate(rhs.itsFixedEndDate),
	itsEndDate(rhs.itsEndDate), 
	itsRefDate(rhs.itsRefDate),

	itsDomesticCcy(rhs.itsDomesticCcy),
	itsForeignCcy(rhs.itsForeignCcy),
	itsFundingCcy(rhs.itsFundingCcy),

	itsCpnDaycount(rhs.itsCpnDaycount),      
	itsCpnFreq(rhs.itsCpnFreq),          
	itsFxResetGap(rhs.itsFxResetGap),
	itsCpnResetCal(rhs.itsCpnResetCal),      
	itsCpnPayCal(rhs.itsCpnPayCal),        
	itsStubRule(rhs.itsStubRule),  
	itsResetTiming(rhs.itsResetTiming),

	itsFundFreq(rhs.itsFundFreq),
	itsFundDaycount(rhs.itsFundDaycount),

	itsExerciseFreq(rhs.itsExerciseFreq),
	itsNoticeGap(rhs.itsNoticeGap),
	itsPayRec(rhs.itsPayRec),
	itsFees(rhs.itsFees),
	itsNbNoCall(rhs.itsNbNoCall),

	itsCpnnominalCv(rhs.itsCpnnominalCv),
	itsDomesticCpnCv(rhs.itsDomesticCpnCv),
	itsForeignCpnCv(rhs.itsForeignCpnCv),
	itsMinCpnCv(rhs.itsMinCpnCv),
	itsMaxCpnCv(rhs.itsMaxCpnCv),
	itsInitialFxCv(rhs.itsInitialFxCv),
	itsFundnominalCv(rhs.itsFundnominalCv),
	itsFundSpreadCv(rhs.itsFundSpreadCv),
	itsFundlevrageCv(rhs.itsFundlevrageCv),
	
	itsNbFixFlows(rhs.itsNbFixFlows),
	itsRedemptionType(rhs.itsRedemptionType),
	itsRedemptionStrike(rhs.itsRedemptionStrike),
	itsRedemptionGap(rhs.itsRedemptionGap),
	itsResetRedemptionDate(rhs.itsResetRedemptionDate),

	itsDomesticCpn(rhs.itsDomesticCpn),
	itsForeignCpn(rhs.itsForeignCpn),
	itsMinCpn(rhs.itsMinCpn),
	itsMaxCpn(rhs.itsMaxCpn),
	itsInitialFx(rhs.itsInitialFx),

	itsvFundSpread(rhs.itsvFundSpread),
	itsfundMargin(rhs.itsfundMargin),
	itsBasisType(rhs.itsBasisType),
	itsvInitialFundSpread(rhs.itsvInitialFundSpread),
	itsvFundNominal(rhs.itsvFundNominal),
	itsFundSize(rhs.itsFundSize),
	itsvFundIndex(rhs.itsvFundIndex),
	itsvInitialFundNominal(rhs.itsvInitialFundNominal),
	itsvFundLeverage(rhs.itsvFundLeverage),


	itsvCpnNominal(rhs.itsvCpnNominal),
	itsvLeverage(rhs.itsvLeverage),
	itsvFixCpn(rhs.itsvFixCpn),	
	itsvLowStrikeCpn(rhs.itsvLowStrikeCpn),	
    itsvHighStrikeCpn(rhs.itsvHighStrikeCpn),	
	itsvCpnIsCapped(rhs.itsvCpnIsCapped),
	itsvCpnIsFloored(rhs.itsvCpnIsFloored),
	itsvCpnIsFixed(rhs.itsvCpnIsFixed),
	itsCpnSize(rhs.itsCpnSize),
	itsvCpnIndex(rhs.itsvCpnIndex),

	itsvIsExerDate(rhs.itsvIsExerDate),
	itsExerSize(rhs.itsExerSize),

	itsFundDateStripFrom(rhs.itsFundDateStripFrom),
	itsFundDateStrip(rhs.itsFundDateStrip),
	itsStructDateStrip(rhs.itsStructDateStrip),
	itsForexDateStrip(rhs.itsForexDateStrip),
	itsExerciseDateStrip(rhs.itsExerciseDateStrip),

    itsColumnsToPrice( rhs.itsColumnsToPrice),
    itsSchedulerDatas( rhs.itsSchedulerDatas),
    itsTruncatorDatas( rhs.itsTruncatorDatas),

    itsDomModelHWFlag( rhs.itsDomModelHWFlag),
    itsForModelHWFlag( rhs.itsForModelHWFlag),

	itsFxResetTimeBefore(rhs.itsFxResetTimeBefore),
    itsMarkovianDriftSamplerFlag( rhs.itsMarkovianDriftSamplerFlag),

	itsFlooredFxOptionPF( CreateClone(rhs.itsFlooredFxOptionPF)),
	itsCappedFxOptionPF( CreateClone(rhs.itsCappedFxOptionPF)),
	itsRedemptionFxOptionPF( CreateClone(rhs.itsRedemptionFxOptionPF)),

    itsFxLocalModelFlag( rhs.itsFxLocalModelFlag),
    itsCalibType( rhs.itsCalibType),
    itsCalibDatas( rhs.itsCalibDatas),

	itsATMFwdFxVol(CreateClone(rhs.itsATMFwdFxVol)),
    itsBasisIRCalibFlag ( rhs.itsBasisIRCalibFlag),
	itsATMDoubleCalibMethod ( CreateClone(rhs.itsATMDoubleCalibMethod)),

	itsAnalyticalModel(CreateClone(rhs.itsAnalyticalModel)),
	itsPowerReverseSwap(CreateClone(rhs.itsPowerReverseSwap)),
	itsMktHybridModel(CreateClone(rhs.itsMktHybridModel)),
	itsHybridBasketCalibMethod(CreateClone(rhs.itsHybridBasketCalibMethod)),

	itsHasBeenComputed(rhs.itsHasBeenComputed),
	itsTradeFromDataBase(rhs.itsTradeFromDataBase)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_PRDCalculator::~ARM_PRDCalculator()
{
	delete itsPowerReverseSwap;
	itsPowerReverseSwap = NULL;

	delete itsAnalyticalModel;
	itsAnalyticalModel = NULL;

    delete itsFlooredFxOptionPF;
    itsFlooredFxOptionPF=NULL;

    delete itsCappedFxOptionPF;
    itsCappedFxOptionPF=NULL;

    delete itsRedemptionFxOptionPF;
    itsRedemptionFxOptionPF=NULL;

    delete itsATMFwdFxVol;
    itsATMFwdFxVol=NULL;

    delete itsATMDoubleCalibMethod;
    itsATMDoubleCalibMethod=NULL;

	delete itsMktHybridModel;
	itsMktHybridModel=NULL;

	delete itsHybridBasketCalibMethod;
	itsHybridBasketCalibMethod=NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDCalculator::ARM_PRDCalculator(const ARM_Date& asofDate,
		const ARM_Date& startDate,
		const ARM_Date& fixEndDate,
		const ARM_Date& endDate,
		const ARM_Currency& domCcy,
		const ARM_Currency& forCcy,
		const ARM_Currency& fundCcy,
		int cpnDayCount,
		int cpnFreq,
		int FxResetGap,
		int stubRule,
		int resetTiming,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const ARM_Curve& cpnnominal,
		const ARM_Curve& domesticCpn,
		const ARM_Curve& foreignCpn,
		const ARM_Curve& minCpn,
		const ARM_Curve& maxCpn,
		const ARM_Curve& initialFx,
		int fundFreq,
		int fundDayCount,
		const ARM_Curve& fundnominal,
		const ARM_Curve& fundSpread,
		int exerciseFreq,
		int noticeGap,
		int payRec,
		size_t nbNCall,
		const ARM_Curve& fees,
		int redemptionGap ,
		double redemptionStrike,
		const ARM_RedemptionType& redemptionType,
		const ARM_StringVector& columnsToPrice,
		bool fxLocalModelFlag,
		bool basisIRCalibFlag,
		const ARM_Curve& fundlevrage)
:	
	ARM_GenCalculator(asofDate),

	itsStartDate(startDate),        
	itsEndDate(endDate),
    itsRefDate(startDate),
	itsFixedEndDate(fixEndDate),

	itsDomesticCcy(domCcy),
	itsForeignCcy(forCcy),
	itsFundingCcy(fundCcy),

	itsCpnDaycount(cpnDayCount),      
	itsCpnFreq(cpnFreq),          
	itsFxResetGap(FxResetGap),
	itsCpnResetCal(cpnResetCal),      
	itsCpnPayCal(cpnPayCal),        
	itsStubRule(stubRule), 
	itsResetTiming(resetTiming),

	itsFundFreq(fundFreq),
	itsFundDaycount(fundDayCount),

	itsExerciseFreq(exerciseFreq),
	itsNoticeGap(noticeGap),
	itsPayRec(payRec),
	itsFees(fees),
	itsNbNoCall(nbNCall),


	itsCpnnominalCv(cpnnominal),
	itsDomesticCpnCv(domesticCpn),
	itsForeignCpnCv(foreignCpn),
	itsMinCpnCv(minCpn),
	itsMaxCpnCv(maxCpn),
	itsInitialFxCv(initialFx),
	itsFundnominalCv(fundnominal),
	itsFundSpreadCv(fundSpread),
	itsFundlevrageCv(fundlevrage),

	itsRedemptionType(redemptionType),
	itsRedemptionStrike(redemptionStrike),
	itsRedemptionGap(redemptionGap),
	itsResetRedemptionDate(),

	itsDomesticCpn(0),
	itsForeignCpn(0),
	itsMinCpn(0),
	itsMaxCpn(0),
	itsInitialFx(0),

	itsvFundSpread(0),
	itsBasisType(ARM_PRCSBasisType::flowByflow),
	itsfundMargin(),
	itsvInitialFundSpread(0),
	itsvCpnNominal(0),
	itsvFundNominal(0),
	itsvInitialFundNominal(0),
	itsFundSize(0),
	itsvFundIndex(0),
	itsvLeverage(0),
	itsvFundLeverage(0),
	itsvFixCpn(0),
	itsvCpnIsCapped(0),
	itsvCpnIsFloored(0),
	itsvCpnIsFixed(0), 
	itsCpnSize(0),
	itsvCpnIndex(0),
	itsFundDateStrip(0),
	itsFundDateStripFrom(0),
	itsStructDateStrip(0),
	itsForexDateStrip(0),
	itsExerciseDateStrip(0),
	itsNbFixFlows(0),
    itsColumnsToPrice(columnsToPrice),
    itsDomModelHWFlag(false),
    itsForModelHWFlag(false),
	itsFxResetTimeBefore(0),
    itsFlooredFxOptionPF(NULL),
    itsCappedFxOptionPF(NULL),
    itsRedemptionFxOptionPF(NULL),
	itsFxLocalModelFlag(fxLocalModelFlag),   
    itsATMFwdFxVol(NULL),
    itsBasisIRCalibFlag(basisIRCalibFlag),
    itsATMDoubleCalibMethod(NULL),
	itsPowerReverseSwap(NULL),
	itsAnalyticalModel(NULL),
	itsMktHybridModel(NULL),
	itsHybridBasketCalibMethod(NULL),
	itsHasBeenComputed(false),
	itsTradeFromDataBase(false)

{
	DatesStructure();

	/// fix period managing by fixEndDate
	while(itsNbFixFlows < itsStructDateStrip->GetFlowEndDates()->size() && 
		(*itsStructDateStrip->GetFlowEndDates())[itsNbFixFlows] < fixEndDate.GetJulian()+ ARM_GlobalConstant::ARM_SEVENDAYS_LAG )
		itsNbFixFlows++;

	/// Check input datas
    CheckDataAndTimeIt();

    /// Set the domestic=coupon=payment currency
    SetCurrencyUnit(&const_cast<ARM_Currency&>(domCcy));
	SetFundingCcy(const_cast<ARM_Currency&> (fundCcy));
    SetDomesticCcy(const_cast<ARM_Currency&> (domCcy));
	SetForeignCcy(const_cast<ARM_Currency&> (forCcy));

    string domCcyName( domCcy.GetCcyName() );
    string forCcyName( forCcy.GetCcyName() );
	string fundCcyName( fundCcy.GetCcyName() );
    string fxName(forCcyName + "/" + domCcyName);
	string fxfundName(fundCcyName + "/" + domCcyName);
    string domForFxName("(" + domCcyName + "," + forCcyName + "," + fxName + ")");

    ARM_StringVector mdmKeys(NbKeys);
    mdmKeys[YcDomKey]               = YC_KEY_NAME               + domCcyName;
    mdmKeys[YcForKey]               = YC_KEY_NAME               + forCcyName;
	mdmKeys[YcFundKey]              = YC_KEY_NAME               + fundCcyName;
    mdmKeys[ForexKey]               = FOREX_KEY_NAME            + fxName;
	mdmKeys[FunForexKey]            = FOREX_KEY_NAME			+ fxfundName;
    mdmKeys[YcBasisDomKey]          = YC_BASIS_KEY_NAME         + domCcyName;
    mdmKeys[YcBasisForKey]          = YC_BASIS_KEY_NAME         + forCcyName;
	mdmKeys[YcBasisFundKey]         = YC_BASIS_KEY_NAME         + fundCcyName;
    mdmKeys[OswDomModelKey]         = OSWMODEL_KEY_NAME         + domCcyName;
    mdmKeys[OswForModelKey]         = OSWMODEL_KEY_NAME         + forCcyName;
    mdmKeys[FxModelKey]             = FXMODEL_KEY_NAME          + fxName;
    mdmKeys[CorrelMatrixKey]        = CORREL_KEY_NAME           + domForFxName;
    mdmKeys[MrsDomKey]              = MRS_KEY_NAME              + domCcyName;
    mdmKeys[MrsForKey]              = MRS_KEY_NAME              + forCcyName;
    mdmKeys[QFxKey]                 = Q_KEY_NAME                + fxName;
    mdmKeys[QDomKey]                = Q_KEY_NAME                + domCcyName;
    mdmKeys[QForKey]                = Q_KEY_NAME                + forCcyName;
    mdmKeys[FlooredFxModelKey]      = FLOORED_FOREX_KEY_NAME    + fxName;
    mdmKeys[CappedFxModelKey]       = CAPPED_FOREX_KEY_NAME     + fxName;
    mdmKeys[RedemptionFxModelKey]   = REDEMPTION_FOREX_KEY_NAME + fxName;
	mdmKeys[BSFxVol]				= BSFXVol_KEY_NAME			+ fxName;
    mdmKeys[LocalFxModelKey]		= LOCAL_FXMODEL_KEY_NAME    + fxName;

    SetKeys(mdmKeys);
	
	ComputeProductVectorsFromCurves();   	
	
}
	/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Init
///	Returns: void
///	Action : Initialize the object with MKT
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::Init(const ARM_MarketData_ManagerRep& mktDataManager,
		const ARM_GP_Vector& schedulerDatas,
		const ARM_GP_Vector& truncatorDatas,
		bool markovianDriftSamplerFlag,
		ARM_PRDCCalibType calibType,
		const ARM_GP_Vector& calibDatas)
{

	itsSchedulerDatas = schedulerDatas;
	itsTruncatorDatas = truncatorDatas;
	itsMarkovianDriftSamplerFlag = markovianDriftSamplerFlag;
	itsCalibType = calibType;
	itsCalibDatas = calibDatas;

	/// Register input objects but no more than the internal number of access keys
    size_t nbToReg = GetKeys().size();
    for(size_t i(0); i<nbToReg; ++i)
	{
		if(!mktDataManager.TestIfKeyMissing(GetKeys()[i]))
			GetMktDataManager()->RegisterData(GetKeys()[i],mktDataManager.GetData(GetKeys()[i]));
	}
	GetMktDataManager()->SetDetailMode(mktDataManager.GetDetailMode());

	//Set the convexityAdjustment model

	/// Check market datas
    CheckMktDataAndTimeIt();

	
	/// Convert le basis from funding Ccy to domestic Ccy
	if(IsBasis()){
		ComputeDomesticBasis();
		ARM_GramFctorArg cst(ARM_VectorPtr((ARM_GP_Vector*)itsvFundSpread.Clone()));
		GetGenSecurity()->GetCstManager()->insert(PRDProfileNamesTable[MarginProfile],cst );
	}

	/// Create a 2IR+FX model
    CreateAndSetModelAndTimeIt();
	
/// Create calibration sets
    CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRDCalculator::ARM_PRDCalculator(
                ARM_PowerReverse* powRev,
                ARM_DFBSModel* model,
                const ARM_ObjectVector& otherMktDatas,
                const ARM_GP_Vector& schedulerDatas,
                const ARM_GP_Vector& truncatorDatas,
                const ARM_StringVector& columnsToPrice,
                bool markovianDriftSamplerFlag,
                bool fxLocalModelFlag,
                ARM_PRDCCalibType calibType,
                const ARM_GP_Vector& calibDatas,
                bool basisIRCalibFlag,
				ARM_BasisType basisType)
:	
	ARM_GenCalculator(model ? model->GetStartDate() : ARM_Date()),
	itsPowerReverseSwap( CreateClone(powRev) ),
	itsAnalyticalModel( CreateClone(model) ),
    itsSchedulerDatas(schedulerDatas),
    itsTruncatorDatas(truncatorDatas),
    itsColumnsToPrice(columnsToPrice),
    itsDomModelHWFlag(false),
    itsForModelHWFlag(false),
    itsMarkovianDriftSamplerFlag(markovianDriftSamplerFlag),
    itsFlooredFxOptionPF(NULL),
    itsCappedFxOptionPF(NULL),
    itsRedemptionFxOptionPF(NULL),
    itsFxLocalModelFlag(fxLocalModelFlag),
    itsCalibType(calibType),
    itsCalibDatas(calibDatas),
    itsATMFwdFxVol(NULL),
    itsBasisIRCalibFlag(basisIRCalibFlag),
    itsATMDoubleCalibMethod(NULL),
	itsMktHybridModel(NULL),
	itsHybridBasketCalibMethod(NULL),
	itsStartDate(powRev->GetItsRealFundLeg()->GetStartDateNA()),        
	itsEndDate(powRev->GetItsRealFundLeg()->GetEndDateNA()),  
	itsRefDate(powRev->GetInitFixedLeg() && powRev->GetInitFixedLeg()->GetRefDate() ? (ARM_Date)powRev->GetInitFixedLeg()->GetRefDate(): itsStartDate),
	itsFixedEndDate(powRev->GetInitFixedLeg() ? powRev->GetInitFixedLeg()->GetEndDateNA(): itsStartDate),
	itsDomesticCcy(*powRev->GetItsFxNumLeg()->GetCurrencyUnit()),
	itsForeignCcy(*powRev->GetItsFxUnderLeg()->GetCurrencyUnit()),
	itsFundingCcy(*powRev->GetItsRealFundLeg()->GetCurrencyUnit()),
	itsCpnDaycount(powRev->GetItsFxUnderLeg()->GetDayCount()),      
	itsCpnFreq(powRev->GetItsFxUnderLeg()->GetIRIndex()->GetResetFrequency()),          
	itsFxResetGap(),
	itsCpnResetCal(powRev->GetItsRealFundLeg()->GetResetCalName()),      
	itsCpnPayCal(powRev->GetItsRealFundLeg()->GetPayCalName()),        
	itsStubRule(powRev->GetItsFxUnderLeg()->GetStubMeth()), 
	itsResetTiming(powRev->GetItsFxUnderLeg()->GetIRIndex()->GetResetTiming()),
	itsFundFreq(powRev->GetItsFxUnderLeg()->GetIRIndex()->GetResetFrequency()),
	itsFundDaycount(powRev->GetItsRealFundLeg()->GetDayCount()),
	itsExerciseFreq(),
	itsNoticeGap(),
	itsPayRec(),
	itsFees(),
	itsCpnnominalCv(),
	itsDomesticCpnCv(),
	itsForeignCpnCv(),
	itsMinCpnCv(),
	itsMaxCpnCv(),
	itsInitialFxCv(),
	itsFundnominalCv(),
	itsFundSpreadCv(),
	itsFundlevrageCv(ARM_FlatCurve(1.0)),
	itsNbNoCall(0),
	itsFxResetTimeBefore(0),
	itsRedemptionType((ARM_RedemptionType)powRev->GetDualOptionFlag()),
	itsRedemptionStrike(powRev->GetDualOptionStrike()),
	itsResetRedemptionDate(powRev->GetRedempNoticeDate()),
	itsRedemptionGap(),
	itsDomesticCpn(0),
	itsForeignCpn(0),
	itsMinCpn(0),
	itsMaxCpn(0),
	itsInitialFx(0),
	itsvFundSpread(0),
	itsBasisType(basisType),
	itsfundMargin(),
	itsvInitialFundSpread(0),
	itsvCpnNominal(0),
	itsvFundNominal(0),
	itsvInitialFundNominal(0),
	itsFundSize(0),
	itsvFundIndex(0),
	itsvLeverage(0),
	itsvFundLeverage(0),
	itsvFixCpn(0),
	itsvCpnIsCapped(0),
	itsvCpnIsFloored(0),
	itsvCpnIsFixed(0), 
	itsCpnSize(0),
	itsvCpnIndex(0),
	itsFundDateStrip(0),
	itsFundDateStripFrom(0),
	itsStructDateStrip(0),
	itsForexDateStrip(0),
	itsExerciseDateStrip(0),
	itsNbFixFlows(0),
	itsHasBeenComputed(false),
	itsTradeFromDataBase(true)
{

    powRev->SetAsOfEqNoticeFlag(true);

    /// Set default values for scheduler & truncator datas if necessary
    if(itsSchedulerDatas.size()<6)
    {
        itsSchedulerDatas.resize(6);
        itsSchedulerDatas[0]=4;     // Nb min steps before 1st notice
        itsSchedulerDatas[1]=20;    // Nb max steps before 1st notice
        itsSchedulerDatas[2]=4;     // Nb steps per year before 1st notice
        itsSchedulerDatas[3]=7.0;   // Optimal date in year fraction
        itsSchedulerDatas[4]=2;     // Nb steps after 1st notice & before optimal
        itsSchedulerDatas[5]=1;     // Nb steps after optimal
    }

    if(itsTruncatorDatas.size()<4)
    {
        itsTruncatorDatas.resize(4);
        itsTruncatorDatas[0]=4;         // Nb max std dev
        itsTruncatorDatas[1]=25;        // Nb max space steps (for each half space)
        itsTruncatorDatas[2]=1.0e-4;    // Arrow-Debreu limit price for truncation
        itsTruncatorDatas[3]=4.0;       // Booster date in year fraction that allow AD strategy truncation
    }

    /// Check input datas
    CheckDataAndTimeIt();

    /// Set the domestic=coupon=payment currency
    SetCurrencyUnit(&itsDomesticCcy);
	SetFundingCcy(itsFundingCcy);
    SetDomesticCcy(itsDomesticCcy);
	SetForeignCcy(itsForeignCcy);

    string domCcyName( itsDomesticCcy.GetCcyName() );
    string forCcyName( itsForeignCcy.GetCcyName() );
	string fundCcyName( itsFundingCcy.GetCcyName() );
    string fxName(forCcyName + "/" + domCcyName);
	string fxfundName(fundCcyName + "/" + domCcyName);
    string domForFxName("(" + domCcyName + "," + forCcyName + "," + fxName + ")");

    ARM_StringVector mdmKeys(NbKeys);
    mdmKeys[YcDomKey]               = YC_KEY_NAME               + domCcyName;
    mdmKeys[YcForKey]               = YC_KEY_NAME               + forCcyName;
	mdmKeys[YcFundKey]              = YC_KEY_NAME               + fundCcyName;
    mdmKeys[ForexKey]               = FOREX_KEY_NAME            + fxName;
	mdmKeys[FunForexKey]            = FOREX_KEY_NAME			+ fxfundName;
    mdmKeys[YcBasisDomKey]          = YC_BASIS_KEY_NAME         + domCcyName;
    mdmKeys[YcBasisForKey]          = YC_BASIS_KEY_NAME         + forCcyName;
	mdmKeys[YcBasisFundKey]         = YC_BASIS_KEY_NAME         + fundCcyName;
    mdmKeys[OswDomModelKey]         = OSWMODEL_KEY_NAME         + domCcyName;
    mdmKeys[OswForModelKey]         = OSWMODEL_KEY_NAME         + forCcyName;
    mdmKeys[FxModelKey]             = FXMODEL_KEY_NAME          + fxName;
    mdmKeys[CorrelMatrixKey]        = CORREL_KEY_NAME           + domForFxName;
    mdmKeys[MrsDomKey]              = MRS_KEY_NAME              + domCcyName;
    mdmKeys[MrsForKey]              = MRS_KEY_NAME              + forCcyName;
    mdmKeys[QFxKey]                 = Q_KEY_NAME                + fxName;
    mdmKeys[QDomKey]                = Q_KEY_NAME                + domCcyName;
    mdmKeys[QForKey]                = Q_KEY_NAME                + forCcyName;
    mdmKeys[FlooredFxModelKey]      = FLOORED_FOREX_KEY_NAME    + fxName;
    mdmKeys[CappedFxModelKey]       = CAPPED_FOREX_KEY_NAME     + fxName;
    mdmKeys[RedemptionFxModelKey]   = REDEMPTION_FOREX_KEY_NAME + fxName;

    /// Market FX model to compute target prices for local FX model calibration
    mdmKeys[LocalFxModelKey]        = LOCAL_FXMODEL_KEY_NAME    + fxName;

	/// Market IR model to compute market price of variable notional domestic swaption
	/// part of the new hybrid basket option calibration stuff
    mdmKeys[MarketIrModelKey]		= MARKET_IRMODEL_KEY_NAME	+ domCcyName;

    SetKeys(mdmKeys);	
	
	/// Collect market datas from the input model and intialize the MdM
    FillMarketDataManager(otherMktDatas);

	/// Compute redemption datas & booster table then reaffect
    /// notice dates to relevant flows
    CreateDataFromUnderlying();

    /// Check market datas
    CheckMktDataAndTimeIt();

    /// Create a 2IR+FX model
    CreateAndSetModelAndTimeIt();
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: FillMarketDataManager
///	Returns: ARM_ObjectVector
///	Action : Collect market datas from the input model
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::FillMarketDataManager(const ARM_ObjectVector& mktDatas)
{
    /// Collect market datas from the input model
    /// -----------------------------------------

    /// Domestic vanilla swaptions are computed using a standard non basis BS model
    ARM_BSModel* oswDomModel = itsAnalyticalModel->GetDBSModel();
    GetMktDataManager()->RegisterData(GetKeys()[OswDomModelKey],static_cast< ARM_Object* >(oswDomModel));

    ARM_ZeroCurve* ycDomCurve = itsAnalyticalModel->GetDBSModel()->GetZeroCurve();
    GetMktDataManager()->RegisterData(GetKeys()[YcDomKey],static_cast< ARM_Object* >(ycDomCurve));

    /// Foreign vanilla swaptions are computed using a standard non basis BS model
    ARM_BSModel* oswForModel = itsAnalyticalModel->GetFBSModel();
    GetMktDataManager()->RegisterData(GetKeys()[OswForModelKey],static_cast< ARM_Object* >(oswForModel));

    ARM_ZeroCurve* ycForCurve = itsAnalyticalModel->GetFBSModel()->GetZeroCurve();
    GetMktDataManager()->RegisterData(GetKeys()[YcForKey],static_cast< ARM_Object* >(ycForCurve));

    //ARM_BSModel* fxModel = itsAnalyticalModel->CreateFxBSModel();

	/// To support new FX GP models !
    ARM_DFBSModel* fxModel = itsAnalyticalModel;
    GetMktDataManager()->RegisterData(GetKeys()[FxModelKey],static_cast< ARM_Object* >(fxModel));

    double fxSpot = itsAnalyticalModel->GetFxSpot();
    ARM_Forex forex(ycForCurve->GetCurrencyUnit(),ycDomCurve->GetCurrencyUnit(), fxSpot);
    GetMktDataManager()->RegisterData(GetKeys()[ForexKey],static_cast< ARM_Object* >(&forex));

	///Funding Zero Coupon Yield Curve if we have to do
	if(IsBasis())
	{
		ARM_ZeroCurve* ycFunCurve = itsAnalyticalModel->GetFundCrv();
		GetMktDataManager()->RegisterData(GetKeys()[YcFundKey],static_cast< ARM_Object* >(ycFunCurve));
		///Funding Zero Coupon Yield Curve
		ARM_ZeroCurve* ycBasisFunCurve = itsAnalyticalModel->GetFundBSCrv();
		GetMktDataManager()->RegisterData(GetKeys()[YcBasisFundKey],static_cast< ARM_Object* >(ycBasisFunCurve));
		double fundfxSpot = itsAnalyticalModel->GetFundFxSpot();
		ARM_Forex fundforex(&GetFundingCcy(),GetCurrencyUnit(), fundfxSpot);
		GetMktDataManager()->RegisterData(GetKeys()[FunForexKey],static_cast< ARM_Object* >(&fundforex));
	}

    ARM_ZeroCurve* ycBasisDomCurve = itsAnalyticalModel->GetDBsCrv();
    GetMktDataManager()->RegisterData(GetKeys()[YcBasisDomKey],static_cast< ARM_Object* >(ycBasisDomCurve));

    ARM_ZeroCurve* ycBasisForCurve = itsAnalyticalModel->GetFBsCrv();
    GetMktDataManager()->RegisterData(GetKeys()[YcBasisForKey],static_cast< ARM_Object* >(ycBasisForCurve));

    ARM_VolFlat* domFxCorrel = dynamic_cast<ARM_VolFlat*>(itsAnalyticalModel->GetDomFxCorrel());
    ARM_VolFlat* forFxCorrel = dynamic_cast<ARM_VolFlat*>(itsAnalyticalModel->GetForeignFxCorrel());
    if(!domFxCorrel || !forFxCorrel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't restore correlations for input model");

    ARM_GP_Matrix correlMatrix(3,3,1.0);
    correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel) = domFxCorrel->GetVolatility();
    correlMatrix(ARM_2IRFXModel::FxModel,ARM_2IRFXModel::DomModel) = correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
    correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel) = forFxCorrel->GetVolatility();
    correlMatrix(ARM_2IRFXModel::FxModel,ARM_2IRFXModel::ForModel) = correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel);

    correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel) = itsAnalyticalModel->GetRatesCorr()/100.0;
    correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::DomModel) = correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel);
    GetMktDataManager()->RegisterData(GetKeys()[CorrelMatrixKey],static_cast< ARM_Object* >(&correlMatrix));

    /// Get other market datas using the input object vector
    /// ----------------------------------------------------

    /// Mean Reversions (on FX not input because it doesn't matter !)
    ARM_CurveModelParam* mrsDom = dynamic_cast<ARM_CurveModelParam*>(mktDatas[MrsDomKey]);
    if(!mrsDom)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic MRS is missing");
    GetMktDataManager()->RegisterData(GetKeys()[MrsDomKey],static_cast< ARM_Object* >(mrsDom));

    ARM_CurveModelParam* mrsFor = dynamic_cast<ARM_CurveModelParam*>(mktDatas[MrsForKey]);
    if(!mrsFor)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign MRS is missing");
    GetMktDataManager()->RegisterData(GetKeys()[MrsForKey],static_cast< ARM_Object* >(mrsFor));

    /// Q Parameters
    ARM_CurveModelParam* qFx = dynamic_cast<ARM_CurveModelParam*>(mktDatas[QFxKey]);
    /// If no Q parameter is input, model will be pure log-normal (Q=1)
    qFx ? GetMktDataManager()->RegisterData(GetKeys()[QFxKey],static_cast< ARM_Object* >(qFx)) : "  pure lognormal model";

    ARM_CurveModelParam* qDom = dynamic_cast<ARM_CurveModelParam*>(mktDatas[QDomKey]);
	/// If no Q parameter is input, model will be degenerated in H&W
	qDom ? GetMktDataManager()->RegisterData(GetKeys()[QDomKey],static_cast< ARM_Object* >(qDom)) : " Hull & White model";

    ARM_CurveModelParam* qFor = dynamic_cast<ARM_CurveModelParam*>(mktDatas[QForKey]);
	/// If no Q parameter is input, model will be degenerated in H&W
	qFor ? GetMktDataManager()->RegisterData(GetKeys()[QForKey],static_cast< ARM_Object* >(qFor)) : " Hull & White model";

    if(itsFxLocalModelFlag)
    {
        ARM_DFBSModel* localDFBSModel = dynamic_cast<ARM_DFBSModel*>(mktDatas[LocalFxModelKey]);
        if(localDFBSModel)
        {
			/// To support new FX GP models !
            GetMktDataManager()->RegisterData(GetKeys()[LocalFxModelKey],static_cast< ARM_Object* >(localDFBSModel));
        }
        else
        {
            /// Use the default one computed by the dual currencies BS model
            GetMktDataManager()->RegisterData(GetKeys()[LocalFxModelKey],static_cast< ARM_Object* >(fxModel));
        }
    }

	ARM_MarketIRModel* marketIrModel = dynamic_cast<ARM_MarketIRModel*>(mktDatas[MarketIrModelKey]);
	if(itsCalibType == ARM_PRCSCalibTypes::HybridBasketCalib && marketIrModel)
	{
		/// If market IR model is missing, VNS evaluation will be degenerated
		GetMktDataManager()->RegisterData(GetKeys()[MarketIrModelKey],static_cast< ARM_Object* >(marketIrModel));
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeDomesticBasis
///	Returns: void
///	Action : Convert le basis from funding Ccy to domestic Ccy
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeDomesticBasis()
{
	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* ycfundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundKey]));
	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisFundKey]));
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[FunForexKey]));

	ARM_BasisConverter basisConveter(itsDomesticCcy,
					itsFundingCcy,
					itsStructDateStrip,
					itsFundDateStripFrom,
					itsFundDateStrip,
					ycDomCurve,	
					ycfundCurve,
					ycBasisDomCurve,
					basisFundCurve,
					*forex,
					itsvCpnNominal,
					itsvInitialFundNominal,
					itsvInitialFundSpread);

	ARM_GP_Vector vdomMargin = basisConveter.ComputeDomMargin();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();	
	ARM_GP_Vector* fundResetDates = itsFundDateStrip->GetResetDates() ;
	ARM_GP_Vector* cpnResetDates = itsStructDateStrip->GetResetDates();
	ARM_GP_Vector timeLags(itsCpnSize);
	for (int i(0); i<itsCpnSize; i++)
		timeLags[i] = (*cpnResetDates)[i]-asOfDate;

	/// we create a local curve to resize.....
	ARM_Curve cvMargin(timeLags,vdomMargin, new ARM_StepUpLeftOpenCstExtrapolDble);

	for (i=0; i<itsFundSize; i++)
	{
		double lag = (*fundResetDates)[i]-asOfDate;
		itsvFundSpread[i] = (itsBasisType == ARM_PRCSBasisType::average) ? itsfundMargin.Interpolate(lag) :  cvMargin.Interpolate(lag);
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeProductVectorsFromCurves()
{
	size_t i;
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();	
	
	/// Convert ARM_Curves into relevant vectors
	/// this will be much easier for computations
	ARM_GP_Vector* fundResetDates = itsFundDateStrip->GetResetDates() ;
	
	double lag,margin, fundLvge;

	ARM_GP_Vector* originalfundResetDates = itsFundDateStripFrom->GetResetDates() ;
	for (i=0; i<originalfundResetDates->size(); i++)
	{
		lag = (*originalfundResetDates)[i]-asOfDate;
		margin = itsFundSpreadCv.Interpolate(lag);
		fundLvge = itsFundlevrageCv.Interpolate(lag);
		if(IsBasis()){
			itsvInitialFundSpread.push_back(margin/fundLvge);
			itsvInitialFundNominal.push_back(itsFundnominalCv.Interpolate(lag));
		}
	}
	for (i=0; i<itsFundSize; i++)
	{
		lag = (*fundResetDates)[i]-asOfDate;
		fundLvge =  IsBasis() ? itsFundSpreadCv.Interpolate(lag) : 1.0 ;
		margin = itsFundSpreadCv.Interpolate(lag)/fundLvge;
		itsvFundSpread.push_back(margin );
		itsvFundNominal.push_back(itsCpnnominalCv.Interpolate(lag));
		itsvFundLeverage.push_back(itsFundlevrageCv.Interpolate(lag));
	}

	ARM_GP_Vector* cpnResetDates = itsStructDateStrip->GetResetDates() ;
	double fixRate,minCpn,maxCpn,domesticCpn,foreignCpn,
		initialFx,leverage,lowStrike,highStrike;
	
	for (i=0; i<itsCpnSize; i++)
	{
		lag = (*cpnResetDates)[i]-asOfDate;
		itsvCpnNominal.push_back(itsCpnnominalCv.Interpolate(lag) );

		if (i < itsNbFixFlows)
		{
			itsvCpnIsFloored.push_back(false);
			itsvCpnIsCapped.push_back(false);
			itsvCpnIsFixed.push_back(true),
			itsvLeverage.push_back(0.0);
			fixRate = itsDomesticCpnCv.Interpolate(lag);
			itsvFixCpn.push_back(fixRate);

			itsDomesticCpn.push_back(0.0);
			itsForeignCpn.push_back(0.0);
			itsInitialFx.push_back(1.0);
			itsMinCpn.push_back(fixRate);
			itsMaxCpn.push_back(fixRate);

			itsvLowStrikeCpn.push_back(fixRate);
			itsvHighStrikeCpn.push_back(fixRate);
		}
		else
		{	
			minCpn = itsMinCpnCv.Interpolate(lag);
			(minCpn < NO_FLOOR) ? itsvCpnIsFloored.push_back(false) : itsvCpnIsFloored.push_back(true);
			itsvFixCpn.push_back(minCpn);
			itsMinCpn.push_back(minCpn);
				
			maxCpn = itsMaxCpnCv.Interpolate(lag);
			(maxCpn > NO_CAP) ? itsvCpnIsCapped.push_back(false) : itsvCpnIsCapped.push_back(true);
			itsMaxCpn.push_back(maxCpn);

			domesticCpn = itsDomesticCpnCv.Interpolate(lag);
			itsDomesticCpn.push_back(domesticCpn);
			foreignCpn  = itsForeignCpnCv.Interpolate(lag);
			itsForeignCpn.push_back(foreignCpn);
			initialFx   = itsInitialFxCv.Interpolate(lag);
			itsInitialFx.push_back(initialFx);

			if(initialFx <K_NEW_DOUBLE_TOL)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ": Invalid intial Fx value,make sur you don't put a no positive value");
			leverage = foreignCpn/initialFx;
			if(leverage <K_NEW_DOUBLE_TOL){
				itsvCpnIsFloored[i] = false;
				itsvCpnIsCapped[i]  = false;
				itsMaxCpn[i]  = minCpn;
				itsDomesticCpn[i] = 0.0;
				itsForeignCpn[i]= 0.0;
				itsInitialFx[i]= 1.0;
				itsvCpnIsFixed.push_back(true),
				itsvLeverage.push_back(0.0);				
				itsvLowStrikeCpn.push_back(minCpn);
				itsvHighStrikeCpn.push_back(minCpn);
			}
			else{
				itsvCpnIsFixed.push_back(false),
				itsvLeverage.push_back(leverage);
				lowStrike = (domesticCpn + minCpn)/leverage;
				highStrike  = (domesticCpn + maxCpn)/leverage; 
				itsvLowStrikeCpn.push_back(lowStrike);
				itsvHighStrikeCpn.push_back(highStrike);
			}
		}
	}
		
	
	ARM_GP_Vector* exerciseDates = itsExerciseDateStrip->GetResetDates() ;
	if(!itsFees.empty())
	{
		double fee; 
		for (i=0; i<itsExerSize; i++)
		{
			lag = (*exerciseDates)[i] - asOfDate;
			fee =  itsFees.Interpolate(lag);
			((fee > NON_CALL_FEE ) || (i < itsNbNoCall)) ? itsvIsExerDate.push_back(false) : itsvIsExerDate.push_back(true);
		}
		//itsNbNoCall = CC_NS(std,distance)(itsvIsExerDate.begin(),itsvIsExerDate.find(true));
	}

	/// last: compute indexes to relate exer - funding - cpn		
	ARM_GP_Vector* cpnStartDates  = itsStructDateStrip->GetFlowStartDates() ;
	ARM_GP_Vector* fundStartDates = itsFundDateStrip->GetFlowStartDates() ;
	
	i = 0;
	for (size_t k(0); k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*fundStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "PRDC calculator : no call date allowed between fixing dates and start dates");
		itsvFundIndex.push_back(i);
	}
	for (k=0, i=0; k<itsExerSize; k++)
	{
		while ((*exerciseDates)[k]>(*cpnStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*cpnResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "PRDC calculator : no call date allowed between fixing dates and start dates");
		itsvCpnIndex.push_back(i);
	}
	itsvFundIndex.push_back(fundStartDates->size());
	itsvCpnIndex.push_back(cpnStartDates->size());
		
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<itsExerSize; k++)
	{
		if (itsvFundIndex[k]>=itsvFundIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "PRDC calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Funding) is forbidden");

		if (itsvCpnIndex[k]>=itsvCpnIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "PRDC calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Coupon) is forbidden");
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	for (k=0; k<itsExerSize; k++) 
	{		
		if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*cpnStartDates)[itsvCpnIndex[k]])> ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			ARM_THROW( ERR_INVALID_ARGUMENT, "PRDC calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeProductVectorsFromCurves
///	Returns: void
///	Action : compute vector attributes from ARM_Curve 
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CreateDataFromUnderlying()
{
	/// Compute redemption datas & booster table then reaffect
    /// notice dates to relevant flows
	ComputeBoosterDatas();

	double firstEventDate = (*itsExerciseDateStrip->GetResetDates())[0];
	size_t exersize = itsExerciseDateStrip->size();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	/// Get DateStrip From SwapLeg.
	/// 2- original funding schedule (sometimes, it's not synchnonised with the structured leg)
	ARM_DateStrip* OriginalFundingSched = (itsBasisType == ARM_PRCSBasisType::average) ? SwapLegToDateStrip(*itsPowerReverseSwap->GetItsFictivFundLeg()):
																              SwapLegToDateStrip(*itsPowerReverseSwap->GetItsRealFundLeg());
	CC_NS(std,auto_ptr)<ARM_DateStrip> HoldFundingSched(OriginalFundingSched);

	size_t nbPastNoCall=0;
	ARM_GP_Vector* exerciseDates = OriginalFundingSched->GetResetDates() ;
	size_t size = OriginalFundingSched->size();
	while(nbPastNoCall < size && (*exerciseDates)[nbPastNoCall] < firstEventDate + K_NEW_DOUBLE_TOL)
		++nbPastNoCall;
	OriginalFundingSched->ResizeAndBuilt(nbPastNoCall,size);

	itsFundDateStripFrom = ARM_DateStripPtr(new ARM_DateStrip(*OriginalFundingSched));
	if(itsBasisType == ARM_PRCSBasisType::average) 
		itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(*OriginalFundingSched));

	/// 3- structured schedule
	ARM_DateStrip* FxNumSched = SwapLegToDateStrip(*itsPowerReverseSwap->GetItsFxNumLeg());
	CC_NS(std,auto_ptr)<ARM_DateStrip> HoldFxNumSched(FxNumSched);

	ARM_SwapLeg* fixedLeg = itsPowerReverseSwap->GetInitFixedLeg();
	ARM_DateStrip* StructureSched = fixedLeg ? SwapLegToDateStrip(*fixedLeg) : NULL;
	CC_NS(std,auto_ptr)<ARM_DateStrip> HoldStructureSched(StructureSched);

	if(StructureSched)
		StructureSched->fill(*FxNumSched);
	else
		StructureSched = FxNumSched;

	nbPastNoCall=0;
	exerciseDates = StructureSched->GetResetDates() ;
	size = StructureSched->size();
	while(nbPastNoCall < size && (*exerciseDates)[nbPastNoCall] < firstEventDate + K_NEW_DOUBLE_TOL)
		++nbPastNoCall;
	StructureSched->ResizeAndBuilt(nbPastNoCall,size);

	itsStructDateStrip = ARM_DateStripPtr(new ARM_DateStrip(*StructureSched));
	size = itsStructDateStrip->size();
	itsStructDateStrip->ResizeAndBuilt(size-exersize,size);

		/// to set the real dates ( start, end,it,...) to exerciseStripDate
	itsExerciseDateStrip->SetFlowStartDates(itsStructDateStrip->GetFlowStartDates());
	itsExerciseDateStrip->SetFlowEndDates(itsStructDateStrip->GetFlowEndDates());
	itsExerciseDateStrip->SetFwdRateStartDates(itsStructDateStrip->GetFwdRateStartDates());
	itsExerciseDateStrip->SetFwdRateEndDates(itsStructDateStrip->GetFwdRateEndDates());
	itsExerciseDateStrip->SetPaymentDates(itsStructDateStrip->GetPaymentDates());
	itsExerciseDateStrip->SetInterestDays(itsStructDateStrip->GetInterestDays());
	itsExerciseDateStrip->SetInterestTerms(itsStructDateStrip->GetInterestTerms());

	/// 4- forex schedule
	ARM_DateStrip* FxUnderSched = SwapLegToDateStrip(*itsPowerReverseSwap->GetItsFxUnderLeg());
	CC_NS(std,auto_ptr)<ARM_DateStrip> HoldFxUnderSched(FxUnderSched);

	ARM_DateStrip* ForexSched = fixedLeg ? SwapLegToDateStrip(*fixedLeg) : NULL;
	CC_NS(std,auto_ptr)<ARM_DateStrip> HoldForexSched(ForexSched);

	if(ForexSched)
		ForexSched->fill(*FxUnderSched);
	else
		ForexSched = FxUnderSched;

	ForexSched->ResizeAndBuilt(nbPastNoCall,ForexSched->size());
	itsForexDateStrip = ARM_DateStripPtr(new ARM_DateStrip(*ForexSched));

	DatesStructure();

	double refDate = asOfDate + ARM_GlobalConstant::ARM_SEVENDAYS_LAG;

	/// Lower and upper strikes
	itsCpnnominalCv	= RefValueToCurve(itsPowerReverseSwap->GetItsFxNumLeg()->GetAmount(),refDate,1.0,K_STEPUP_LEFT);
	itsDomesticCpnCv	= RefValueToCurve(itsPowerReverseSwap->GetItsFxNumLeg()->GetVarCoupons(),refDate,100,K_STEPUP_LEFT);
	ARM_GP_Vector abscisses = itsDomesticCpnCv.GetAbscisses();
	itsForeignCpnCv	= RefValueToCurve(itsPowerReverseSwap->GetItsFxUnderLeg()->GetVarCoupons(),asOfDate,100,K_STEPUP_LEFT);
	itsForeignCpnCv.SetAbscisses(abscisses);
	itsMinCpnCv		= RefValueToCurve(itsPowerReverseSwap->GetItsFloor(),asOfDate,100,K_STEPUP_LEFT);
	itsMinCpnCv.SetAbscisses(abscisses);
	itsMaxCpnCv		= RefValueToCurve(itsPowerReverseSwap->GetItsCap(),asOfDate,100,K_STEPUP_LEFT);
	itsMaxCpnCv.SetAbscisses(abscisses);
	itsInitialFxCv		= RefValueToCurve(itsPowerReverseSwap->GetItsFxVariable(),asOfDate,1.0,K_STEPUP_LEFT);
	itsInitialFxCv.SetAbscisses(abscisses);

	itsFundnominalCv	= RefValueToCurve(itsPowerReverseSwap->GetItsRealFundLeg()->GetAmount(),refDate,1.0,K_STEPUP_LEFT);
	itsFundSpreadCv	= RefValueToCurve(itsPowerReverseSwap->GetItsRealFundLeg()->GetSpreads(),refDate,100.0,K_STEPUP_LEFT);

	/// fix period managing by fixEndDate
	while(itsNbFixFlows < itsStructDateStrip->GetFlowEndDates()->size() && 
		(*itsStructDateStrip->GetFlowEndDates())[itsNbFixFlows] < itsFixedEndDate.GetJulian()+ ARM_GlobalConstant::ARM_SEVENDAYS_LAG )
		itsNbFixFlows++;

	double fixRate = fixedLeg ? fixedLeg->GetFixedRate()/100 : 0.0;

	itsDomesticCpnCv.insert(0.0, fixRate);

	//ARM_Curve fixRate = RefValueToCurve(fixedLeg ? fixedLeg->GetFixedRates() : NULL ,asOfDate,100,K_STEPUP_LEFT );

	ComputeProductVectorsFromCurves();
	
	if(IsBasis())
		ComputeDomesticBasis();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeCstManager
///	Returns: void
///	Action : Build a constant manager with Fx option date strip &
///          profiles
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_PRDCalculator::ComputeCstManager()
{
    ARM_GP_Vector fxExpiryDates;
    ARM_GP_Vector fxPaymentDates;
	 ARM_GP_VectorPtr vnominal(new ARM_GP_Vector());
    ARM_GP_VectorPtr fxLowStrikes(new ARM_GP_Vector());
    ARM_GP_VectorPtr fxHighStrikes(new ARM_GP_Vector());
    ARM_GP_VectorPtr fxLeverages(new ARM_GP_Vector());
    ARM_GP_Vector fundingPaymentDates;
    ARM_GP_VectorPtr fundingSpreads(new ARM_GP_Vector());

	for(size_t eventIdx = 0;eventIdx<itsExerSize;++eventIdx)
    {
		/// Coupon pay date
		double payDate=(*(itsStructDateStrip->GetPaymentDates()))[itsvCpnIndex[eventIdx]];
        fundingPaymentDates.push_back( payDate);
		fxPaymentDates.push_back( payDate);

		/// funding marging in domestic currecy
		double fundMargin = itsvFundSpread[itsvFundIndex[eventIdx]];
        fundingSpreads->push_back(fundMargin);

		/// forex reset date
		double fxResetDate = (*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[eventIdx]];
        fxExpiryDates.push_back(  fxResetDate);
       
		/// low strike of call spread
		double lowStrike = itsvLowStrikeCpn[itsvCpnIndex[eventIdx]];
        fxLowStrikes->push_back(lowStrike);

		/// high strike of call spread
		double highStrike = itsvHighStrikeCpn[itsvCpnIndex[eventIdx]];
        fxHighStrikes->push_back(highStrike);
		
		/// levrage times interestperiod for callStrip Option
		double term=(*(itsStructDateStrip->GetInterestTerms()))[itsvCpnIndex[eventIdx]];
		double leverage = itsvLeverage[itsvCpnIndex[eventIdx]];
        fxLeverages->push_back(leverage*term);

		/// domestic nominal
		double nominal = itsvCpnNominal[itsvCpnIndex[eventIdx]];
		vnominal->push_back(nominal);
    }

    ARM_StringVector names(0);
    vector <ARM_GramFctorArg> values(0);

    /// Funding  stripDates
    names.push_back(PRDProfileNamesTable[FundingDateStripProfile]);
    values.push_back(ARM_GramFctorArg(itsFundDateStrip));

	/// Cpn stripDates
    names.push_back(PRDProfileNamesTable[CpnDateStripProfile]);
    values.push_back(ARM_GramFctorArg(itsStructDateStrip));
	
	/// funding margin times notional
	ARM_GP_VectorPtr spreadTimesNotioptr = ARM_GP_VectorPtr(new ARM_GP_Vector(itsvFundSpread * itsvFundNominal) );
    names.push_back(PRDProfileNamesTable[MarginTimesNotional]);
    values.push_back(ARM_GramFctorArg(spreadTimesNotioptr));

	//ARM_GP_VectorPtr  vMarginptr ((ARM_GP_Vector*)itsvFundSpread.Clone());
	ARM_GP_VectorPtr  vMarginptr (new ARM_GP_Vector(itsvFundSpread));

    names.push_back(PRDProfileNamesTable[MarginProfile]);
    values.push_back(ARM_GramFctorArg(vMarginptr));


    /// Fx option datas
    ARM_DateStripPtr fxDateStrip(new ARM_DateStrip());
    fxDateStrip->SetResetDates(&fxExpiryDates);
    fxDateStrip->SetPaymentDates(&fxPaymentDates);

    names.push_back(PRDProfileNamesTable[FxDateStripProfile]);
    values.push_back(ARM_GramFctorArg(fxDateStrip));

    names.push_back(PRDProfileNamesTable[FxLowStrikeProfile]);
    values.push_back(ARM_GramFctorArg(fxLowStrikes));

    names.push_back(PRDProfileNamesTable[NominalProfile]);
    values.push_back(ARM_GramFctorArg(vnominal));

	names.push_back(PRDProfileNamesTable[FxHighStrikeProfile]);
    values.push_back(ARM_GramFctorArg(fxHighStrikes));

    names.push_back(PRDProfileNamesTable[FxLeverageProfile]);
    values.push_back(ARM_GramFctorArg(fxLeverages));


	ARM_GP_Vector cpnFloorTimesNotio = itsvFixCpn * itsvCpnNominal;

	ARM_GP_VectorPtr fundLvgeTimesNotioptr = ARM_GP_VectorPtr(new ARM_GP_Vector(itsvFundLeverage * itsvFundNominal) );
	names.push_back(PRDProfileNamesTable[FundLvgeTimesNotional]);
	values.push_back(ARM_GramFctorArg(fundLvgeTimesNotioptr));
	

	for (size_t i(0); i<itsExerSize; i++)
	{		
		/// define index
		CC_Ostringstream indexStr;
		/// sort by size
		i < 10 ? indexStr << 0 << i :  indexStr << i;
	
		/// names
		names.push_back("CpnMinTimesNotio"					+ indexStr.str());
		
		/// values
		values.push_back( ARM_GramFctorArg(ARM_VectorPtr(new ARM_GP_Vector(cpnFloorTimesNotio.begin()	  + itsvCpnIndex[i],	 cpnFloorTimesNotio.end()))) );
	}
	
	return  ARM_CstManagerPtr(new ARM_CstManager(names, values));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if PRDC datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CheckData()
{
    /// Check if input columns to price are present in the deal description
    size_t colNamesSize = sizeof(PRDColNamesTable)/sizeof(PRDColNamesTable[0]);
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        for(size_t j=0;j<colNamesSize;++j)
            if(itsColumnsToPrice[i] == PRDColNamesTable[j])
                break;
        if(j==itsColumnsToPrice.size())
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't price this unknown column name " + itsColumnsToPrice[i]);  
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if PRDC market datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CheckMktData()
{
    /// Market datas checking

	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
    if(!domCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcDomKey] + " is expected in the Market Data Manager");
    string domCcy(domCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
    if(!basisDomCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcBasisDomKey] + " is expected in the Market Data Manager");
    string basisDomCcy(basisDomCurve->GetCurrencyUnit()->GetCcyName());

    if(domCcy != basisDomCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic(=coupon) Basis Curve currency should consistent with reference curve");

	ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]));
    if(!forCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcForKey] + " is expected in the Market Data Manager");
    string forCcy(forCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    if(!basisForCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcBasisForKey] + " is expected in the Market Data Manager");
    string basisForCcy(basisForCurve->GetCurrencyUnit()->GetCcyName());

    if(forCcy != basisForCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Basis Curve currency should consistent with reference curve");

	ARM_ZeroCurve* fundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcFundKey]));
    if(!fundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcFundKey] + " is expected in the Market Data Manager");
    string fundCcy(fundCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisFundKey]));
    if(!basisFundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + GetKeys()[YcBasisFundKey] + " is expected in the Market Data Manager");
    string basisFundCcy(basisFundCurve->GetCurrencyUnit()->GetCcyName());

    if(fundCcy != basisFundCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Funding Basis Curve currency should consistent with reference curve");

	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
	if(!forex )
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + GetKeys()[ForexKey] + " is expected in the Market Data Manager");

	ARM_Forex* fundforex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[FunForexKey]));
	if(!fundforex && IsBasis())
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + GetKeys()[FunForexKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    if(!oswDomBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + GetKeys()[OswDomModelKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    if(!oswForBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + GetKeys()[OswForModelKey] + " is expected in the Market Data Manager");

    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    if(!fxBSModel)
	{
		/// Build it from market data manager and replace it
		bool useMktDataMgrCorrels=true;
		fxBSModel = CreateAnalyticalModel(FxModelKey,useMktDataMgrCorrels);
		GetMktDataManager()->RegisterData(GetKeys()[FxModelKey],static_cast< ARM_Object* >(fxBSModel));
		delete fxBSModel; /// cloned by RegisterData
	}

	ARM_ModelParam* mrsDomParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsDomKey]) );
    if(!mrsDomParam || mrsDomParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + GetKeys()[MrsDomKey] + " is expected in the Market Data Manager");

	ARM_ModelParam* mrsForParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(GetKeys()[MrsForKey]) );
    if(!mrsForParam || mrsForParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign MRS Param for key=" + GetKeys()[MrsForKey] + " is expected in the Market Data Manager");

    /// Q parameters are not tested because if missing,
    /// the 2IRFXModel will be degenerated

    if(itsFxLocalModelFlag)
    {
		/// To support new FX GP models !
		if(!GetMktDataManager()->TestIfKeyMissing(GetKeys()[LocalFxModelKey]))
		{
			ARM_DFBSModel* localFxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[LocalFxModelKey]) );
			if(!localFxBSModel)
			{
				/// Build it from market data manager and replace it
				bool useMktDataMgrCorrels=true;
				localFxBSModel = CreateAnalyticalModel(LocalFxModelKey,useMktDataMgrCorrels);
				GetMktDataManager()->RegisterData(GetKeys()[LocalFxModelKey],static_cast< ARM_Object* >(localFxBSModel));
				delete localFxBSModel; /// cloned by RegisterData
			}
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : local fx B&S model for key=" + GetKeys()[LocalFxModelKey] + " is expected in the Market Data Manager");
    }

	if(itsCalibType == ARM_PRCSCalibTypes::HybridBasketCalib)
	{
		/// If market IR is missing IR part in the hybrid basket
		/// will be priced by an ATM calibrated model
		ARM_MarketIRModel* marketIrModel=NULL;
		if(!GetMktDataManager()->TestIfKeyMissing(GetKeys()[MarketIrModelKey]))
		{
			marketIrModel = dynamic_cast<ARM_MarketIRModel*>( GetMktDataManager()->GetData(GetKeys()[MarketIrModelKey]) );
			if(!marketIrModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : market IR model is not of good type for hybrid basket calibration");
		}
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDCalculator::ColumnNames() const
{
    size_t colNamesSize = sizeof(PRDColNamesTable)/sizeof(PRDColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = PRDColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC.
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_PRDCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeBoosterDatas
///	Returns: void
///	Action : Build booster table & redemption vector then
///          associate notice dates to the right flow
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeBoosterDatas()
{	
    /// Compute the PRDC analytical price, to generate fictive funding spreads,
    /// low & high strikes, leverage etc... saved in THE booster matrix !
    itsPowerReverseSwap->SetModel(itsAnalyticalModel);
    itsPowerReverseSwap->SetFXSmileAdjFlag(true); // to allow market vol smile adj to be saved for further using
    double prdcSwapPrice = itsPowerReverseSwap->ComputePrice();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// Get the BoosterDatas table
    ARM_Matrix* dataMatrix = itsPowerReverseSwap->CalcReformattedPRCSDataMatrix();
	CC_NS(std,auto_ptr)<ARM_Matrix> HolddataMatrix(dataMatrix);

    if(!dataMatrix)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't get the PRDC data matrix");

	ARM_GP_Matrix BoosterDatas(To_ARM_GP_Matrix(*dataMatrix));
	size_t nbFlows = BoosterDatas.rows();

    /// Initialise exercise flag vector to say if a notice is actual or not
    ARM_BoolVector exerciseFlag(nbFlows,false);

    /// Affect notice dates to relevant flows and
    /// activate exercise flag for an actual exercise
    double noticeDate,defaultNoticeDate;
    ARM_GP_Vector noticeDates(nbFlows,0.0);

    int lastNoticeIdx=-1;
    int flowIdx=0;
    int firstNoticeIdx=-1;
	itsNbNoCall = 0;
    for(size_t noticeIdx=0;noticeIdx<nbFlows && (noticeDate=BoosterDatas(noticeIdx,NOTICE_DATE)) > 0;++noticeIdx)
    {
        while(flowIdx < nbFlows && (defaultNoticeDate=CC_Min(BoosterDatas(flowIdx,FX_RESET_DATE),BoosterDatas(flowIdx,FUNDING_START_DATE))) < noticeDate)
        {
            /// No call date => notice date is set to minimum value between
            /// fx reset date & funding start date
            noticeDates[flowIdx] = defaultNoticeDate - ARM_GlobalConstant::ARM_SEVENDAYS_LAG;
			++flowIdx;
        }
        if(flowIdx < nbFlows)
        {
            noticeDates[flowIdx] = noticeDate;
            exerciseFlag[flowIdx] = true;

            /// Keep record of the 1st notice position after the as of date
            if(firstNoticeIdx==-1 && floor(noticeDate) > asOfDate + K_NEW_DOUBLE_TOL)
                firstNoticeIdx=flowIdx;

            lastNoticeIdx = flowIdx;
            ++flowIdx;
        }
    }

    /// Update NoticeDate column in the booster table
    for(flowIdx=0;flowIdx <= lastNoticeIdx;++flowIdx)
        BoosterDatas(flowIdx,NOTICE_DATE) = noticeDates[flowIdx];

    /// After last notice date, just set notice to minimum value between
    /// fx reset date & funding start date
    for(;flowIdx < nbFlows;++flowIdx)
        BoosterDatas(flowIdx,NOTICE_DATE) = CC_Min(BoosterDatas(flowIdx,FX_RESET_DATE),BoosterDatas(flowIdx,FUNDING_START_DATE)) - ARM_GlobalConstant::ARM_SEVENDAYS_LAG;


	/// Save all futur Fx reset times associated to a notice date before the 1st one
    double fxResetTime;
    for(flowIdx=0;flowIdx < firstNoticeIdx;++flowIdx)
    {
        if((fxResetTime=BoosterDatas(flowIdx,FX_RESET_DATE)-asOfDate) >= 0)
            itsFxResetTimeBefore.push_back(fxResetTime);
    }

    /// Until a generic security can't have different lines with same event date
    /// (here asOfDate will be the best), it will not possible to price the deal
    /// before the 1st notice date without adding no callable fx reset dates !
    /// Here to get the same diffusion schedule as Tree3F, this is not posssible
    /// to add such artificial event dates
    /// Do the same for exercise flags that begin necessary to a "true"
	ARM_GP_Vector realNoticeDates;
    if(firstNoticeIdx > -1)
    {
        ARM_GP_Matrix initialBoosterDatas(BoosterDatas);
        BoosterDatas.resize(initialBoosterDatas.rows()-firstNoticeIdx,initialBoosterDatas.cols());

        for(flowIdx=0; flowIdx<BoosterDatas.rows(); ++flowIdx)
        {
            itsvIsExerDate.push_back(exerciseFlag[firstNoticeIdx+flowIdx]);
            for(int j=0;j<BoosterDatas.cols();++j)
                BoosterDatas(flowIdx,j)=initialBoosterDatas(firstNoticeIdx+flowIdx,j);
			realNoticeDates.push_back(floor(BoosterDatas(flowIdx,NOTICE_DATE)));
        }
    }
	
	ARM_DateStrip ExerSched(&realNoticeDates,&realNoticeDates,&realNoticeDates,&realNoticeDates,&realNoticeDates,
                                &realNoticeDates,&realNoticeDates,&realNoticeDates);

	ARM_GP_Vector* margin = BoosterDatas.GetColumn(FUNDING_SPREAD);
	CC_NS(std,auto_ptr)<ARM_GP_Vector> holdmargin(margin);

	ARM_StepUpLeftOpenCstExtrapolDble* interpolator = new ARM::ARM_StepUpLeftOpenCstExtrapolDble();
	itsfundMargin = ARM_Curve(realNoticeDates-asOfDate, *margin,interpolator);

	itsExerciseDateStrip = ARM_DateStripPtr(new ARM_DateStrip(ExerSched));	
	itsExerSize = ExerSched.size();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_PRDCalculator::DatesStructure() const
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	
	int fwdRule			= K_MOD_FOLLOWING;	// for forward dates
	int intRule			= K_ADJUSTED;			// for interest dates
	int resetTiming		= K_ADVANCE;
	int payTiming		= K_ARREARS;
	
	ARM_INDEX_TYPE indexType = ((ARM_Currency*)GetCurrencyUnit())->GetVanillaIndexType();
	char* DefaultresetCalendar = GetCurrencyUnit()->GetResetCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdresetCalendar(DefaultresetCalendar);
	char* DefaultpayCalendar  = GetCurrencyUnit()->GetPayCalName(indexType);
	CC_NS(std,auto_ptr)<char> holdpayCalendar(DefaultpayCalendar);

	const char* resetCalendar	= itsCpnResetCal == "" ? DefaultresetCalendar : itsCpnResetCal.c_str();
	const char* payCalendar		= itsCpnPayCal == "" ? DefaultpayCalendar: itsCpnPayCal.c_str();

	size_t nbPastNoCall;
	ARM_GP_Vector* exerciseDates = NULL;
	int indexFirstNotice =0;
	size_t size;
	/// 1- exercise schedule
	/// generate schedules only once ...
	if (itsExerciseDateStrip == ARM_DateStripPtr(NULL))
	{
		ARM_DateStrip ExerSched(itsStartDate,itsEndDate,itsExerciseFreq,itsCpnDaycount,resetCalendar,fwdRule,intRule,itsStubRule,
			-fabs(itsNoticeGap),itsExerciseFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);
	
		/// past  no call
		nbPastNoCall=0;
	    exerciseDates = ExerSched.GetResetDates() ;
		size = ExerSched.size();
		while(nbPastNoCall < size && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		
		indexFirstNotice = CC_Max(nbPastNoCall,itsNbNoCall);
		const_cast< ARM_PRDCalculator* >(this)->itsNbNoCall= itsNbNoCall > nbPastNoCall ? itsNbNoCall - nbPastNoCall : 0;

		/// Keep only the futur notification dates
		ExerSched.ResizeAndBuilt(nbPastNoCall,size);

		const_cast< ARM_PRDCalculator* >(this)->itsExerciseDateStrip = ARM_DateStripPtr(new ARM_DateStrip(ExerSched));	
		const_cast< ARM_PRDCalculator* >(this)->itsExerSize = ExerSched.size();		
	}
	double firstNoticeDate = (*itsExerciseDateStrip->GetResetDates())[0];
	
		/// 3- structured schedule
	/// generate schedules only once ...
	if (itsStructDateStrip == ARM_DateStripPtr(NULL))
	{
		/// 4- structured schedule
		ARM_DateStrip StructureSched(itsStartDate,itsEndDate,itsCpnFreq,itsCpnDaycount,resetCalendar,fwdRule,intRule,
									 itsStubRule,GETDEFAULTVALUE,itsCpnFreq, GETDEFAULTVALUE, payCalendar, resetTiming, payTiming);

		nbPastNoCall=0;
		exerciseDates = StructureSched.GetResetDates() ;
		size = StructureSched.size();
		while(nbPastNoCall < size && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		//while( nbPastNoCall < size  && (*exerciseDates)[nbPastNoCall] < firstNoticeDate - ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
		//	++nbPastNoCall;

		StructureSched.ResizeAndBuilt(nbPastNoCall,size);
		const_cast< ARM_PRDCalculator* >(this)->itsStructDateStrip = ARM_DateStripPtr( new ARM_DateStrip(StructureSched));

		if(fabs(itsFxResetGap) > fabs(itsNoticeGap) + K_NEW_DOUBLE_TOL)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Fx Reset gap must be less than notification gap, please advise");

		/// 4'- forex schedule
		ARM_DateStrip ForexSched(itsStartDate,itsEndDate,itsCpnFreq,itsCpnDaycount,resetCalendar,fwdRule,intRule,itsStubRule,
								-fabs(itsFxResetGap) ,itsCpnFreq, GETDEFAULTVALUE, payCalendar, itsResetTiming, payTiming);
		for(int i(0); i<indexFirstNotice; ++i)
		{
			double fxResetTime = (*ForexSched.GetResetDates())[i];
			if( asOfDate - K_NEW_DOUBLE_TOL < fxResetTime)
				const_cast< ARM_PRDCalculator* >(this)->itsFxResetTimeBefore.push_back(fxResetTime - asOfDate);
		}

		ForexSched.ResizeAndBuilt(nbPastNoCall,ForexSched.size());
		const_cast< ARM_PRDCalculator* >(this)->itsForexDateStrip = ARM_DateStripPtr( new ARM_DateStrip(ForexSched));

		const_cast< ARM_PRDCalculator* >(this)->itsResetRedemptionDate =  itsEndDate;
		
		/// 5- Reset date creating
		const_cast< ARM_PRDCalculator* >(this)->itsResetRedemptionDate.PreviousBusinessDay(fabs(itsRedemptionGap),const_cast<char*> (resetCalendar));
	}
	double firstCancelDate = (*itsStructDateStrip->GetResetDates())[0];
		
	/// 2- funding schedule
	/// generate schedules only once ...
	if (itsFundDateStrip == ARM_DateStripPtr(NULL))
	{
		char* refDateStr = itsRefDate.PrtStrDate();
		CC_NS(std,auto_ptr)<char> HoldrefDateStr(refDateStr);
		ARM_DateStrip FundingSched(	itsStartDate,itsEndDate,itsFundFreq,itsFundDaycount,resetCalendar,fwdRule,intRule,
								    itsStubRule,GETDEFAULTVALUE,itsFundFreq, GETDEFAULTVALUE, payCalendar, resetTiming,
									payTiming,true,refDateStr);

		nbPastNoCall=0;

		exerciseDates = FundingSched.GetResetDates() ;
		size = FundingSched.size();
		while( nbPastNoCall < size  && (*exerciseDates)[nbPastNoCall] < asOfDate + K_NEW_DOUBLE_TOL)
			++nbPastNoCall;
		while( nbPastNoCall < size  && (*exerciseDates)[nbPastNoCall] < CC_Min(firstNoticeDate,firstCancelDate)  - ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			++nbPastNoCall;
		FundingSched.ResizeAndBuilt(nbPastNoCall,size);

		const_cast< ARM_PRDCalculator* >(this)->itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(FundingSched));

		/// 2'- Original funding schedule
		if (itsFundDateStripFrom == ARM_DateStripPtr(NULL))
			const_cast< ARM_PRDCalculator* >(this)->itsFundDateStripFrom = ARM_DateStripPtr(new ARM_DateStrip(FundingSched));
	}
	
	const_cast< ARM_PRDCalculator* >(this)->itsFundSize = itsFundDateStrip->size();	
	const_cast< ARM_PRDCalculator* >(this)->itsCpnSize = itsStructDateStrip->size();

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(1);
    SchedVect[0] = &*itsExerciseDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect,"ResetDate");

	return EventSchedule;
}
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");

    rowDescVec[FundingLeg] = zeroValue;
    rowTypeVec[FundingLeg] = ARM_DOUBLE;

    rowDescVec[FixedLeg] = zeroValue;
    rowTypeVec[FixedLeg] = ARM_DOUBLE;

    rowDescVec[FxLowCallStrip] = zeroValue;
    rowTypeVec[FxLowCallStrip] = ARM_DOUBLE;

    rowDescVec[FxHighCallStrip] = zeroValue;
    rowTypeVec[FxHighCallStrip] = ARM_DOUBLE;

    rowDescVec[CpnLeg] = zeroValue;
    rowTypeVec[CpnLeg] = ARM_DOUBLE;

    rowDescVec[Redemption] = zeroValue;
    rowTypeVec[Redemption] = ARM_DOUBLE;


    rowDescVec[PRDCFirstEuropean] = zeroValue;
    rowTypeVec[PRDCFirstEuropean] = ARM_DOUBLE;

	rowDescVec[FundingFlow] = zeroValue;
    rowTypeVec[FundingFlow] = ARM_DOUBLE;

	rowDescVec[FixFlow] = zeroValue;
    rowTypeVec[FixFlow] = ARM_DOUBLE;

    rowDescVec[FxStrip] = zeroValue;
    rowTypeVec[FxStrip] = ARM_DOUBLE;

	rowDescVec[PRDCFirstSwap] = zeroValue;
    rowTypeVec[PRDCFirstSwap] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRDCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
	size_t descSize = sizeof(PRDColNamesTable)/sizeof(PRDColNamesTable[0]);

	vector< string > rowDescVec(descSize);
	vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 

	/// EventDate description
	double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
	CC_Ostringstream eventDateDesc;
	eventDateDesc << CC_NS(std,fixed) << eventDate;
	rowDescVec[EventDate] = eventDateDesc.str();
	rowTypeVec[EventDate] = ARM_DATE_TYPE;

	/// --- build cst manager names
	CC_Ostringstream indexStr;
	eventIdx < 10 ? indexStr << 0 << eventIdx :  indexStr << eventIdx;
	string CpnMinTimesNotioStr		("CpnMinTimesNotio"				+ indexStr.str());

	/// Set default 0 value for each column to be able to sum it
	ARM_PRDCalculator::InitPriceableColumns(rowDescVec,rowTypeVec);

	/// Get the model names for domestic and forex models
	string basisDomModelName = GetKeys()[YcBasisDomKey];
	string fxModelName = GetKeys()[ForexKey];

	string nextExerIdx("[i+1]");
	bool isLastEvent = (eventIdx+1>=eventSize);

	/// Funding description : in the booster table the spread is already converted
	/// in the domestic ccy then the leg is described as a simple (but fictive)
	/// domestic floating leg

	double startDate=(*(itsStructDateStrip->GetFlowStartDates()))[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fundStartDesc;
	fundStartDesc << CC_NS(std,fixed) << startDate;
	rowDescVec[StartDate] = fundStartDesc.str();
	rowTypeVec[StartDate] = ARM_DATE_TYPE;

	double endDate=(*(itsStructDateStrip->GetFlowEndDates()))[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fundEndDesc;
	fundEndDesc << CC_NS(std,fixed) << endDate;
	rowDescVec[FundingEndDate] = fundEndDesc.str();
	rowTypeVec[FundingEndDate] = ARM_DATE_TYPE;

	double payDate=(*(itsStructDateStrip->GetPaymentDates()))[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fundPayDesc;
	fundPayDesc << CC_NS(std,fixed) << payDate;
	rowDescVec[FundingPayDate] = fundPayDesc.str();
	rowTypeVec[FundingPayDate] = ARM_DATE_TYPE;

	CC_Ostringstream fundSpreadDesc;
	double fundMargin = itsvFundSpread[itsvFundIndex[eventIdx]];
	fundSpreadDesc << CC_NS(std,fixed) << fundMargin;
	rowDescVec[FundingSpread] = fundSpreadDesc.str();
	rowTypeVec[FundingSpread] = ARM_DOUBLE;

	CC_Ostringstream cpnITDesc;
	double it=(*(itsStructDateStrip->GetInterestTerms()))[itsvCpnIndex[eventIdx]];
	cpnITDesc << CC_NS(std,fixed) << it;
	rowDescVec[CpnIT] = cpnITDesc.str();
	rowTypeVec[CpnIT] = ARM_DOUBLE;

	/// to convert from long to string
	string fundingDayCount	= ARM_ArgConvReverse_DayCount.GetString(itsFundDaycount);
	string fundingFreq		= ARM_ArgConvReverse_MatFrequency.GetString(itsFundFreq);
	
	/// to convert from long to string
	string cpnDayCount	= ARM_ArgConvReverse_DayCount.GetString(itsCpnDaycount);
	string cpnFreq		= ARM_ArgConvReverse_MatFrequency.GetString(itsCpnFreq);


	CC_Ostringstream EndDateDesc;
	EndDateDesc << CC_NS(std,fixed) << itsEndDate.GetJulian();
	rowDescVec[EndDate] = EndDateDesc.str();
	rowTypeVec[EndDate] = ARM_DATE_TYPE;

	/// Spread profile in input through notional argument
	int doubleNotional = (itsBasisType == ARM_PRCSBasisType::average) ? 1 : 0;
	string freq = (itsBasisType == ARM_PRCSBasisType::average) ? cpnFreq:fundingFreq;
	CC_Ostringstream fundingDesc;
	fundingDesc << "SWAP(" << basisDomModelName << "," << PRDColNamesTable[StartDate] << "[i], ";
	fundingDesc << PRDColNamesTable[EndDate] << "[i], 0, PAY,,,"<< freq << ", " << fundingDayCount << ",";
	itsBasisType == ARM_PRCSBasisType::average ? fundingDesc << "0,":fundingDesc  << PRDProfileNamesTable[MarginProfile] << "," ;
	fundingDesc << PRDProfileNamesTable[FundLvgeTimesNotional] << ",," << doubleNotional << ",";
	fundingDesc << PRDProfileNamesTable[FundingDateStripProfile] << ",," << itsvFundIndex[eventIdx] << ")";
	if (itsBasisType == ARM_PRCSBasisType::average){
		fundingDesc << " + ANNUITY(" << basisDomModelName << "," << PRDColNamesTable[StartDate] << "[i],"; 
		fundingDesc << PRDColNamesTable[EndDate] << "[i], " <<  freq << ", ";
		fundingDesc << fundingDayCount << ", "<<  PRDProfileNamesTable[MarginTimesNotional] <<",";
		fundingDesc << PRDProfileNamesTable[FundingDateStripProfile] << "," << itsvFundIndex[eventIdx] << ")";
	}

	rowDescVec[FundingLeg] = fundingDesc.str();
	rowTypeVec[FundingLeg] = ARM_STRING;

	CC_Ostringstream cpnPayDesc;
	cpnPayDesc << CC_NS(std,fixed) << payDate;
	rowDescVec[CpnPayDate] = cpnPayDesc.str();
	rowTypeVec[CpnPayDate] = ARM_DATE_TYPE;

	double leverage = itsvLeverage[itsvCpnIndex[eventIdx]];
	CC_Ostringstream cpnLeverageDesc;
	cpnLeverageDesc << CC_NS(std,fixed) << leverage;
	rowDescVec[CpnLeverage] = cpnLeverageDesc.str();
	rowTypeVec[CpnLeverage] = ARM_DOUBLE;

	/// Coupon description : flow is always described but residual coupons
	/// are described only for local model case
	string FlooredFxModelName(fxModelName);
	string CappedFxModelName(fxModelName);
	string RedemptionFxModelName(fxModelName);
	if(itsFxLocalModelFlag)
	{
		FlooredFxModelName      = GetKeys()[FlooredFxModelKey];
		CappedFxModelName       = GetKeys()[CappedFxModelKey];
		RedemptionFxModelName   = GetKeys()[RedemptionFxModelKey];
	}
	CC_Ostringstream cpnFlowDesc;
	CC_Ostringstream cpnDesc;

	double fxResetDate = (*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fxResetDesc;
	fxResetDesc << CC_NS(std,fixed) << fxResetDate;
	rowDescVec[FxResetDate] = fxResetDesc.str();
	rowTypeVec[FxResetDate] = ARM_DATE_TYPE;

	double lowStrike = itsvLowStrikeCpn[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fxLowStrikeDesc;
	fxLowStrikeDesc << CC_NS(std,fixed) << lowStrike;
	rowDescVec[FxLowStrike] = fxLowStrikeDesc.str();
	rowTypeVec[FxLowStrike] = ARM_DOUBLE;

	double highStrike = itsvHighStrikeCpn[itsvCpnIndex[eventIdx]];
	CC_Ostringstream fxHighStrikeDesc;
	fxHighStrikeDesc << CC_NS(std,fixed) << highStrike;
	rowDescVec[FxHighStrike] = fxHighStrikeDesc.str();
	rowTypeVec[FxHighStrike] = ARM_DOUBLE;

	CC_Ostringstream FixLegDesc;
	FixLegDesc << "ANNUITY(" << basisDomModelName << "," << PRDColNamesTable[StartDate] << "[i],"; 
	FixLegDesc << PRDColNamesTable[EndDate] << "[i], " <<  cpnFreq << ", " ;
	FixLegDesc << cpnDayCount << ", "<<  CpnMinTimesNotioStr << ",";
	FixLegDesc << PRDProfileNamesTable[CpnDateStripProfile] << "," << itsvCpnIndex[eventIdx] << ")";

	rowDescVec[FixedLeg] = FixLegDesc.str();
	rowTypeVec[FixedLeg] = ARM_STRING;

	/// Current FX option description
    CC_Ostringstream fxLowCapDesc;
    fxLowCapDesc << CC_NS(std,fixed) << PRDColNamesTable[CpnLeverage] << "[i]*";
    fxLowCapDesc << PRDColNamesTable[CpnIT] << "[i]*";
    fxLowCapDesc << "CALL(" << FlooredFxModelName << "," << PRDColNamesTable[FxResetDate] << "[i],";
    fxLowCapDesc << PRDColNamesTable[FxLowStrike] << "[i],C,,";
    fxLowCapDesc << PRDColNamesTable[CpnPayDate] << "[i])";
    rowDescVec[FxLowCall] = fxLowCapDesc.str();
    rowTypeVec[FxLowCall] = ARM_STRING;

	if(itsvCpnIsCapped[itsvCpnIndex[eventIdx]] || itsvCpnIsFixed[itsvCpnIndex[eventIdx]] )
	{
		CC_Ostringstream fxHighCapDesc;
		fxHighCapDesc << CC_NS(std,fixed) << PRDColNamesTable[CpnLeverage] << "[i]*";
		fxHighCapDesc << PRDColNamesTable[CpnIT] << "[i]*";
		fxHighCapDesc << "CALL(" << CappedFxModelName << "," << PRDColNamesTable[FxResetDate] << "[i],";
		fxHighCapDesc << PRDColNamesTable[FxHighStrike] << "[i],C,,";
		fxHighCapDesc << PRDColNamesTable[CpnPayDate] << "[i])";
		rowDescVec[FxHighCall] = fxHighCapDesc.str();
		rowTypeVec[FxHighCall] = ARM_STRING;
	}

	/// Residual FX option sum description with CALLSTRIP keyword
	CC_Ostringstream fxLowCallStripDesc;
	fxLowCallStripDesc << CC_NS(std,fixed) << "CALLSTRIP(" << FlooredFxModelName << ",";
	fxLowCallStripDesc << PRDColNamesTable[FxResetDate] << "[i],";
	fxLowCallStripDesc << PRDColNamesTable[FxResetDate] << "[i],";
	fxLowCallStripDesc << PRDProfileNamesTable[FxLowStrikeProfile] << ",C,,,,"<< PRDProfileNamesTable[NominalProfile]<< ",,";
	fxLowCallStripDesc << PRDProfileNamesTable[FxLeverageProfile] << ",";
	fxLowCallStripDesc << PRDProfileNamesTable[FxDateStripProfile] << ",";
	fxLowCallStripDesc << eventIdx << ")";
	rowDescVec[FxLowCallStrip] = fxLowCallStripDesc.str();
	rowTypeVec[FxLowCallStrip] = ARM_STRING;

	/// Check if there is a significative residual high strike            
	for(size_t idx=eventIdx;idx<eventSize;++idx)
	{
		if(itsvHighStrikeCpn[itsvCpnIndex[idx]] < MAX_FX_STRIKE)
		{
			CC_Ostringstream fxHighCallStripDesc;
			fxHighCallStripDesc << CC_NS(std,fixed) << "CALLSTRIP(" << CappedFxModelName << ",";
			fxHighCallStripDesc << PRDColNamesTable[FxResetDate] << "[i],";
			fxHighCallStripDesc << PRDColNamesTable[FxResetDate] << "[i],";
			fxHighCallStripDesc << PRDProfileNamesTable[FxHighStrikeProfile] << ",C,,,,"<< PRDProfileNamesTable[NominalProfile]<< ",,";
			fxHighCallStripDesc << PRDProfileNamesTable[FxLeverageProfile] << ",";
			fxHighCallStripDesc << PRDProfileNamesTable[FxDateStripProfile] << ",";
			fxHighCallStripDesc << eventIdx<< ")";
			rowDescVec[FxHighCallStrip] = fxHighCallStripDesc.str();
			rowTypeVec[FxHighCallStrip] = ARM_STRING;
			break;
		}
	}
	cpnDesc << PRDColNamesTable[FxLowCallStrip] << "[i]-" << PRDColNamesTable[FxHighCallStrip] << "[i]";

	if((isLastEvent || itsFxLocalModelFlag) && itsRedemptionType != ARM_PRCSRedemptionType::standard)
	{
		/// Redemption flow description
		CC_Ostringstream redResetDesc;
		redResetDesc << CC_NS(std,fixed) << itsResetRedemptionDate.GetJulian();
		rowDescVec[RedemptionResetDate] = redResetDesc.str();
		rowTypeVec[RedemptionResetDate] = ARM_DATE_TYPE;

		CC_Ostringstream redPayDesc;
		redPayDesc << CC_NS(std,fixed) << itsEndDate.GetJulian();
		rowDescVec[RedemptionPayDate] = redPayDesc.str();
		rowTypeVec[RedemptionPayDate] = ARM_DATE_TYPE;

		CC_Ostringstream redStrikeDesc;
		redStrikeDesc << CC_NS(std,fixed) << itsRedemptionStrike;
		rowDescVec[RedemptionStrike] = redStrikeDesc.str();
		rowTypeVec[RedemptionStrike] = ARM_DOUBLE;

		CC_Ostringstream redNominalDesc;
		redNominalDesc << CC_NS(std,fixed) << itsvCpnNominal[itsvCpnIndex[eventIdx]];
		rowDescVec[Nominal] = redNominalDesc.str();
		rowTypeVec[Nominal] = ARM_DOUBLE;

		CC_Ostringstream redemptionDesc;
		if(itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)
		{
			redemptionDesc << PRDColNamesTable[Nominal] << "[i]*";
			redemptionDesc << "(FWD(" << RedemptionFxModelName << "," << PRDColNamesTable[RedemptionResetDate] << "[i],,";
			redemptionDesc << PRDColNamesTable[RedemptionPayDate] << "[i])/";
			redemptionDesc << PRDColNamesTable[RedemptionStrike] << "[i]" << "-1)*";
			redemptionDesc << "DF(" << basisDomModelName << "," << PRDColNamesTable[RedemptionPayDate] << "[i])";

		}
		else if(itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption)
		{
			redemptionDesc << "-" << PRDColNamesTable[Nominal] << "[i]*";
			redemptionDesc << "CALL(" << RedemptionFxModelName << "," << PRDColNamesTable[RedemptionResetDate] << "[i],";
			redemptionDesc << PRDColNamesTable[RedemptionStrike] << "[i]" << ",P,,";
			redemptionDesc << PRDColNamesTable[RedemptionPayDate] << "[i])/";
			redemptionDesc << PRDColNamesTable[RedemptionStrike] << "[i]";

		}
		else
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown redemption type");

		if(isLastEvent)
		{
			rowDescVec[Redemption] = redemptionDesc.str();
			rowTypeVec[Redemption] = ARM_STRING;
			cpnFlowDesc << "+" << PRDColNamesTable[Redemption] << "[i]";
		}
		if(itsFxLocalModelFlag)
		{
			/// Local model => always add redemption flow
			/// (because of local convexity adjustment and/or local volatility)
			isLastEvent ? cpnDesc << "+" << PRDColNamesTable[Redemption] << "[i]" : cpnDesc << "+" << redemptionDesc.str();
		}
	}

	/// Swap = residual PRD swap
	rowDescVec[CpnLeg] = cpnDesc.str();
	rowTypeVec[CpnLeg] = ARM_STRING;

	/// PRD swap description
	CC_Ostringstream prdcSwapDesc;
	prdcSwapDesc << PRDColNamesTable[FixedLeg] << "[i]+" <<PRDColNamesTable[CpnLeg] << "[i]-" << PRDColNamesTable[FundingLeg] << "[i]";
	rowDescVec[PRDSwap] = prdcSwapDesc.str();
	rowTypeVec[PRDSwap] = ARM_STRING;

	/// The single PRDC option a first line to get the price
	bool isFirstEvent = (eventIdx==itsNbNoCall);
	if(isFirstEvent)
	{
		CC_Ostringstream prdcFundingDesc;
		prdcFundingDesc << PRDColNamesTable[FundingLeg] << "[i]";
		rowDescVec[FundingFlow] = prdcFundingDesc.str();
		rowTypeVec[FundingFlow] = ARM_STRING;

		CC_Ostringstream prdcFixDesc;
		prdcFixDesc << PRDColNamesTable[FixedLeg] << "[i]";
		rowDescVec[FixFlow] = prdcFixDesc.str();
		rowTypeVec[FixFlow] = ARM_STRING;

		CC_Ostringstream FxStripDesc;
		FxStripDesc << PRDColNamesTable[CpnLeg] << "[i]";
		rowDescVec[FxStrip] = FxStripDesc.str();
		rowTypeVec[FxStrip] = ARM_STRING;

		CC_Ostringstream prdcFirstSwapDesc;
        prdcFirstSwapDesc << PRDColNamesTable[PRDSwap] << "[i]";
        rowDescVec[PRDCFirstSwap] = prdcFirstSwapDesc.str();
        rowTypeVec[PRDCFirstSwap] = ARM_STRING;

		CC_Ostringstream prdcFirstEuropeanDesc;
		prdcFirstEuropeanDesc << "EXERCISE(0," << PRDColNamesTable[PRDSwap] << "[i],0)";
		rowDescVec[PRDCFirstEuropean] = prdcFirstEuropeanDesc.str();
		rowTypeVec[PRDCFirstEuropean] = ARM_STRING;

	}

	return ARM_RowInfo(rowDescVec,rowTypeVec);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CreateAndSetModel()
{
    /// Create the 2IR+FX model

    /// Q vol without any curve because will be bootstrapped latter
    ARM_CurveModelParam volParam( ARM_ModelParamType::QVol,0.0,"QVOL");
    ARM_ModelParamVector modelParams(3);
    modelParams[0]  = &volParam;

    /// Multi-assets avec modèles locaux sur le Fx si nécessaire
    int nbModels = ARM_2IRFXModel::NbModels;
    if(!itsFxLocalModelFlag)
        nbModels -= ARM_2IRFXModel::NbLocalModels;
    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    names[ARM_2IRFXModel::DomModel]        = GetKeys()[YcDomKey];
    names[ARM_2IRFXModel::ForModel]        = GetKeys()[YcForKey];
    names[ARM_2IRFXModel::FxModel]         = GetKeys()[ForexKey];	
    names[ARM_2IRFXModel::DomBasisModel]   = GetKeys()[YcBasisDomKey];
    names[ARM_2IRFXModel::ForBasisModel]   = GetKeys()[YcBasisForKey];

	depends[ARM_2IRFXModel::FxModel]       = ARM_StringVector(2);
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::DomModel]	   = names[ARM_2IRFXModel::DomBasisModel];
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::ForModel]	   = names[ARM_2IRFXModel::ForBasisModel];
	depends[ARM_2IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::DomModel]);
	depends[ARM_2IRFXModel::ForBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::ForModel]);
	
    vector< ARM_PricingModelPtr > models(nbModels);
    /// Create the Q1F model for domestic IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsDomKey]));
	ARM_CurveModelParam QmodelParam;

    itsDomModelHWFlag=false;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[QDomKey]))
    {
        /// Degenerate to H&W
        itsDomModelHWFlag = true;
		QmodelParam = ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
        modelParams[2] = &QmodelParam;
    }
    else
	    modelParams[2] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[QDomKey]));

	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
    models[ARM_2IRFXModel::DomModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycDomCurve),ARM_ModelParamsQ1F(modelParams),itsDomModelHWFlag) );
	
    /// Create the Q1F model for foreign IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[MrsForKey]));
	
    itsForModelHWFlag=false;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[QForKey]))
    {
        /// Degenerate to H&W
		itsForModelHWFlag = true;
        QmodelParam = ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
		modelParams[2] =  &QmodelParam;
    }
    else
	    modelParams[2] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[QForKey]));

	ARM_ZeroCurve* ycForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]));
    models[ARM_2IRFXModel::ForModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycForCurve),ARM_ModelParamsQ1F(modelParams),itsForModelHWFlag) );

    /// Create the Q1F model for FX market (with a default MRS=0 because it doesn't matter !)
	ARM_CurveModelParam mrsFxParam( ARM_ModelParamType::MeanReversion,0.0,"FXMRS");
    modelParams[1]   = &mrsFxParam;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[QFxKey]))
    {
        /// Degenerate to pure LN
		QmodelParam = ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
        modelParams[2] =  &QmodelParam;
    }
    else
	    modelParams[2] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(GetKeys()[QFxKey]));

	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]));

    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
    models[ARM_2IRFXModel::FxModel] = ARM_PricingModelPtr( ARM_EqFx_ModelFactory.Instance()->CreateModel( CreateClonedPtr(ycBasisDomCurve), 
									modelParams, forex->GetMarketPrice(), 
									CreateClonedPtr(ycBasisForCurve),
									*correlMatrix,
									ARM_EqFx_ModelFactoryImp::Q1F_Model ) );

    /// Create both domestic & foreign forward margin models
    models[ARM_2IRFXModel::DomBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve)) );
    models[ARM_2IRFXModel::ForBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisForCurve)) );

    /// Create local Fx models if necessary
    if(itsFxLocalModelFlag)
    {
        ARM_Local_SLN_ModelParams* defaultLocalParams = static_cast< ARM_Local_SLN_ModelParams* >( ARM_Local_SLN_Model::CreateDefaultModelParams() );
		CC_NS(std,auto_ptr)<ARM_Local_SLN_ModelParams> holdLocalParams(defaultLocalParams);
        models[ARM_2IRFXModel::FlooredFxLocalModel]     = ARM_PricingModelPtr( new ARM_Local_SLN_Model(*defaultLocalParams) );
        models[ARM_2IRFXModel::CappedFxLocalModel]      = ARM_PricingModelPtr( new ARM_Local_SLN_Model(*defaultLocalParams) );
        models[ARM_2IRFXModel::RedemptionFxLocalModel]  = ARM_PricingModelPtr( new ARM_Local_SLN_Model(*defaultLocalParams) );

		names[ARM_2IRFXModel::FlooredFxLocalModel]      = GetKeys()[FlooredFxModelKey];
        names[ARM_2IRFXModel::CappedFxLocalModel]       = GetKeys()[CappedFxModelKey];
        names[ARM_2IRFXModel::RedemptionFxLocalModel]   = GetKeys()[RedemptionFxModelKey];

		depends[ARM_2IRFXModel::FlooredFxLocalModel]	= ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
		depends[ARM_2IRFXModel::CappedFxLocalModel]		= ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
		depends[ARM_2IRFXModel::RedemptionFxLocalModel] = ARM_StringVector(1,names[ARM_2IRFXModel::FxModel]);
    }
	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	/// Finally, We create and set the brand new 2IR+FX Ferrari model !
    ARM_PricingModelPtr hybridModel = ARM_PricingModelPtr( new ARM_2IRFXModel(modelMap,*correlMatrix) );

    /// Create a 3 factors version of the tree ND then set it
    /// The tree is paramerised to replicate the Kernel version
    int schedulerType=ARM_SchedulerBase::MultiRegime;
    int samplerType = itsMarkovianDriftSamplerFlag ? ARM_SamplerBase::MarkovianDrift: ARM_SamplerBase::DriftedMeanReverting;
    int truncatorType=ARM_TruncatorBase::ArrowDebreu;
    int reconnectorType=ARM_ReconnectorBase::Mean;
    int smootherType=ARM_SmootherBase::Linear;
    ARM_TreeBase* tree = ARM_TreeFactory.Instance()->CreateTreeND(3,schedulerType,itsSchedulerDatas,
						samplerType,ARM_GP_Vector(0),truncatorType,itsTruncatorDatas,false,reconnectorType,smootherType);
    hybridModel->SetNumMethod( ARM_NumMethodPtr( tree ) );

    /// Create and set the cash numeraire
    ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( ARM_Numeraire::Cash ) );
    hybridModel->SetNumeraire(numeraire);

	/// Set the model
	SetPricingModel(hybridModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CreateAndSetCalibration()
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;

    bool isFreezeWeights = false;   // reinit non null weights in the portfolio
    bool isInitVol = true;          // init vol param for calibration (bounds & init guess)


	/// Hybrid basket swaptions (may be also used for domestic calibration)
    /// -------------------------------------------------------------------

	if(itsCalibType == ARM_PRCSCalibTypes::HybridBasketCalib)
	{
		ARM_StdPortfolioPtr hybridBasketPF = CreateHybridBasketOption();

		/// Build the hybrid market model with input market IR & FX models
		ARM_ObjectVector mktDatas(ARM_MarketHybridModel::NbKeys,NULL);
		if(! GetMktDataManager()->TestIfKeyMissing(GetKeys()[MarketIrModelKey]))
			mktDatas[ARM_MarketHybridModel::MarketIrModelKey] = GetMktDataManager()->GetData(GetKeys()[MarketIrModelKey]);

		mktDatas[ARM_MarketHybridModel::MarketEqFxModelKey] = GetMktDataManager()->GetData(GetKeys()[FxModelKey]);

		/// Default correls are not internally computed but
		/// are the same as the 2IRFX ones (instantaneous correlations !)
		/// They could be modified by a first calibration which
		/// yield an integrated correlation structure w.r.t. notice dates
		ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]));
		ARM_GP_VectorPtr hybridDatas( new ARM_GP_Vector(ARM_MarketHybridModel::NbDatas,0.0) );
		(*hybridDatas)[ARM_MarketHybridModel::Time]				=0.0;
		(*hybridDatas)[ARM_MarketHybridModel::DomForCor]		=(*correlMatrix)(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel);
		(*hybridDatas)[ARM_MarketHybridModel::DomFxCor]			=(*correlMatrix)(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
		(*hybridDatas)[ARM_MarketHybridModel::ForFxCor]			=(*correlMatrix)(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel);
		(*hybridDatas)[ARM_MarketHybridModel::IrFxSpreadCor]	=(*correlMatrix)(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
		(*hybridDatas)[ARM_MarketHybridModel::DomO1LogNorVol]	=0.0;
		(*hybridDatas)[ARM_MarketHybridModel::ForO1LogNorVol]	=0.0;
		(*hybridDatas)[ARM_MarketHybridModel::IrNorVol]			=0.0;
		size_t nbBaskets=hybridBasketPF->GetSize();
		ARM_GP_T_Vector< ARM_VectorPtr > hybridDatasStruct(nbBaskets);
		for(size_t i=0;i<nbBaskets;++i)
			hybridDatasStruct[i] = hybridDatas; /// same hybrid datas are shared by default w.r.t. notice schedule
		mktDatas[ARM_MarketHybridModel::HybridDatasKey]	= static_cast< ARM_Object* >(&hybridDatasStruct);

		/// Hybrid basket prices will be computed latter because
		/// Zc models must be calibrated first
		itsMktHybridModel =  new ARM_MarketHybridModel(GetMktDataManager()->GetAsOfDate(),mktDatas);
		const ARM_ModelNameMap& modelMap = * (static_cast< ARM_2IRFXModel* > (&*GetPricingModel()))->GetModelMap();	
		itsMktHybridModel->SetZcModels(modelMap[ARM_2IRFXModel::DomBasisModel]->Model(),modelMap[ARM_2IRFXModel::ForBasisModel]->Model());

		/// Build & save the associated calib method
		itsHybridBasketCalibMethod = CreateFxCalibration(hybridBasketPF,ARM_2IRFXModel::FxModel);
	}

    /// Diagonal domestic swaptions
    /// ---------------------------

    /// Create standard diagonal swaption portfolios (domestic & foreign sigmas calibration)
    pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > diagonalSwaptionPF(CreateDiagonalSwaption());


    /// Build an empty calib method for latter domestic volatility bootstrapping
    ARM_Currency* domCcy = GetCurrencyUnit();
    ARM_CalibMethod* domCalib = CreateIRCalibration(diagonalSwaptionPF.first,domCcy,
                                    itsBasisIRCalibFlag ? ARM_2IRFXModel::DomBasisModel
                                                        : ARM_2IRFXModel::DomModel);

    /// Compute domestic target prices
    ComputeIROptionPrices(domCalib,OswDomModelKey,itsDomModelHWFlag,isFreezeWeights,isInitVol);


    /// Diagonal foreign swaptions
    /// --------------------------
    /// Build an empty calib method for latter foreign volatility bootstrapping
	ARM_Currency* forCcy = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]))->GetCurrencyUnit();
    ARM_CalibMethod* forCalib = CreateIRCalibration(diagonalSwaptionPF.second,forCcy,
                                    itsBasisIRCalibFlag ? ARM_2IRFXModel::ForBasisModel
                                                        : ARM_2IRFXModel::ForModel);

    /// Compute foreign target prices
    ComputeIROptionPrices(forCalib,OswForModelKey,itsForModelHWFlag,isFreezeWeights,isInitVol);


    /// ATM Forex options
    /// -----------------

    /// Create a Fx option portfolio and compute Fx option target prices
    ARM_StdPortfolioPtr fxOptionPF;
    ARM_CalibMethod* fxCalib;
    if(itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib)
    {
        /// Create a standard ATM Fx option portfolio for first calibration step
        itsCalibType    = ARM_PRCSCalibTypes::ATMCalib;
        fxOptionPF      = CreateFxOption();
        fxCalib         = CreateFxCalibration(fxOptionPF,ARM_2IRFXModel::FxModel);
        ComputeFxOptionPrices(fxCalib,isFreezeWeights,isInitVol);
        itsCalibType    = ARM_PRCSCalibTypes::ATMDoubleCalib;
    }
    else
    {
        fxOptionPF      = CreateFxOption();
        fxCalib         = CreateFxCalibration(fxOptionPF,ARM_2IRFXModel::FxModel);
        ComputeFxOptionPrices(fxCalib,isFreezeWeights,isInitVol);
    }


    /// Set next methods to handle calibration with a single method
    domCalib->SetNextMethod(forCalib);
    forCalib->SetNextMethod(fxCalib);

    SetCalibMethod(ARM_CalibMethodPtr(domCalib));

    if(itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib)
    {
        /// Create an ATS moneyness = 100% like Fx option portfolio
        /// for second calibration step
        fxOptionPF = CreateFxOption();
        itsATMDoubleCalibMethod = CreateFxCalibration(fxOptionPF,ARM_2IRFXModel::FxModel);

        /// Compute Fx option target prices (+ ATM strike +...)
        ComputeFxOptionPrices(itsATMDoubleCalibMethod,isFreezeWeights,isInitVol);
    }


    /// Forex options for floored, capped & terminal redemption coupon
    /// --------------------------------------------------------------
    if(itsFxLocalModelFlag)
    {
        /// Create Fx option portfolios used for local Fx model calibration
        CreateLocalFxOption();

        /// Compute prices of floored,capped & terminal redemption Fx option portfolio
        /// Only used for precision because calibration is done directly on volatilities
        /// and not on prices
        ComputeLocalFxOptionPrices();
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions for domestic market
/////////////////////////////////////////////////////////////////
pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > ARM_PRDCalculator::CreateDiagonalSwaption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// Get market swaption standard expiries for domestic & foreign volatilities
    size_t i,nbExp;
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    ARM_VolCurve* oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector domStdExp(nbExp+1);
    domStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        domStdExp[i+1] = asOfDate + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

    oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector forStdExp(nbExp+1);
    forStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        forStdExp[i+1] = asOfDate + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

    ARM_Swaption* swaption;

    list < ARM_Security* > domSwaptionList;
    list < ARM_Security* > forSwaptionList;


    size_t nbEvents = itsvIsExerDate.size();

    /// Standard index of the domestic & foreign currencies
    ARM_Currency* domCcy = GetCurrencyUnit();
    int domSpotDays = domCcy->GetSpotDays();
    ARM_INDEX_TYPE domLiborIndex = domCcy->GetVanillaIndexType();
    char* domResetCal = domCcy->GetResetCalName(domLiborIndex);

	ARM_Currency* forCcy = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]))->GetCurrencyUnit();
    int forSpotDays = forCcy->GetSpotDays();
    ARM_INDEX_TYPE forLiborIndex = forCcy->GetVanillaIndexType();
    char* forResetCal = forCcy->GetResetCalName(forLiborIndex);

    /// Diagonal swaptions are based on funding start/end dates but with the right currency
    ARM_Date swapStartDate,expiryDate;
    ARM_Date swapEndDate((*itsStructDateStrip->GetFlowEndDates())[itsvCpnIndex[nbEvents-1]]);	
		
    bool isDomSwap,isForSwap;
    ARM_Swap domSwap,forSwap;

    /// The 1st notice date is always added then get the very last expiry
    double lastNotice,lastInsertedNotice;
    isDomSwap=false;
    size_t firstNoticeIdx = 0;
    for(size_t eventIdx=0; eventIdx<nbEvents; ++eventIdx)
    {
        if(itsvIsExerDate[eventIdx])
        {
			expiryDate = ARM_Date((*itsExerciseDateStrip->GetResetDates())[eventIdx]);
            if(!isDomSwap)
            {
                swapStartDate   = expiryDate;
                swapStartDate.GapBusinessDay(domSpotDays,domResetCal);

                domSwap = ARM_Swap(swapStartDate,swapEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                 K_DEF_FREQ,K_DEF_FREQ,domCcy);
                swaption = new ARM_Swaption(&domSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                domSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                isDomSwap=true;

                swapStartDate   = expiryDate;
                swapStartDate.GapBusinessDay(forSpotDays,forResetCal);

                forSwap = ARM_Swap(swapStartDate,swapEndDate,forLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                 K_DEF_FREQ,K_DEF_FREQ,forCcy);
                swaption = new ARM_Swaption(&forSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                forSwaptionList.push_back(static_cast< ARM_Security* >(swaption));

                firstNoticeIdx = eventIdx;
                lastInsertedNotice = expiryDate.GetJulian();
            }
            lastNotice = expiryDate.GetJulian();
        }
    }
    if(!isDomSwap)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " creation of calibration diagonal swaption not possible");


    /// Between each standard market expiries, only one swaption is created corresponding to
    /// the earlier notice date
    /// To be consistent with Tree3F code (may be it is a bug !), the very last notice date
    /// is excluded
    size_t domExpIdx=0,forExpIdx=0;
    for(eventIdx = firstNoticeIdx+1;eventIdx < nbEvents; ++eventIdx)
    {
		expiryDate = ARM_Date((*itsExerciseDateStrip->GetResetDates())[eventIdx]);

        if(itsvIsExerDate[eventIdx] && expiryDate.GetJulian() < lastNotice)  /// to be consistent to Tree3F code !!
        {
            isDomSwap=false;
            while(!isDomSwap && domExpIdx < domStdExp.size() && domStdExp[domExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
            {
                if( ( domExpIdx+1 == domStdExp.size() ||
                      (domExpIdx+1 < domStdExp.size() && expiryDate.GetJulian() < domStdExp[domExpIdx+1] + K_NEW_DOUBLE_TOL) )
                    && lastInsertedNotice + K_NEW_DOUBLE_TOL < domStdExp[domExpIdx] )
                {
                    /// Create a domestic underlying swap
                    isDomSwap=true;
                    swapStartDate   = expiryDate;
                    swapStartDate.GapBusinessDay(domSpotDays,domResetCal);
                    domSwap = ARM_Swap(swapStartDate,swapEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                     K_DEF_FREQ,K_DEF_FREQ,domCcy);
                }

                ++domExpIdx;
            }

            isForSwap=false;
            while(!isForSwap && forExpIdx < forStdExp.size() && forStdExp[forExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
            {
                if( ( forExpIdx+1 == forStdExp.size() ||
                      (forExpIdx+1 < forStdExp.size() && expiryDate.GetJulian() < forStdExp[forExpIdx+1] + K_NEW_DOUBLE_TOL) )
                    && lastInsertedNotice + K_NEW_DOUBLE_TOL < forStdExp[forExpIdx] )
                {
                    /// Create a Forestic underlying swap
                    isForSwap=true;
                    swapStartDate   = expiryDate;
                    swapStartDate.GapBusinessDay(forSpotDays,forResetCal);
                    forSwap = ARM_Swap(swapStartDate,swapEndDate,forLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                     K_DEF_FREQ,K_DEF_FREQ,forCcy);
                }

                ++forExpIdx;
            }


            if(isDomSwap)
            {
                /// Domestic swaption
                swaption = new ARM_Swaption(&domSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                domSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                lastInsertedNotice = expiryDate.GetJulian();
            }

            /// Foreign swaption
            if(isForSwap)
            {
                swaption = new ARM_Swaption(&forSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                forSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
                lastInsertedNotice = expiryDate.GetJulian();
            }
        }
    }

    ARM_StdPortfolio* domPort = new ARM_StdPortfolio(domSwaptionList);
    for(i=0;i<domPort->size();++i)
    {
        domPort->SetWeight(OSW_DEFAULT_WEIGHT,i);
        domPort->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }
    ARM_StdPortfolio* forPort = new ARM_StdPortfolio(forSwaptionList);
    for(i=0;i<forPort->size();++i)
    {
        forPort->SetWeight(OSW_DEFAULT_WEIGHT,i);
        forPort->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    /// Don't forget to free char* calendars !
    delete domResetCal;
    delete forResetCal;


    return pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr >(ARM_StdPortfolioPtr(domPort),ARM_StdPortfolioPtr(forPort));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateIRCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRDCalculator::CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
            ARM_Currency* ccy, int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;


    /// Build an empty bootstrap calib method (filled latter for domestic volatility calibration)
    ARM_CalibMethod* volCalib = new ARM_CalibMethod(diagonalSwaptionPF,emptyCalibParam,
                                    ARM_CalibMethodType::Bootstrap1D,ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,noPreviousMethod,
                                    isCalibShared,modelIdx);


    return volCalib;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeIROptionPrice
///	Returns: nothing
///	Action : compute market target prices for calibration purpose
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeIROptionPrices(ARM_CalibMethod* calibMethod, 
	          mdmKeysAlias oswModelIdx, 
			  bool isModelHW, 
			  bool isFreezeWeights, 
			  bool isInitVolParam)
{
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[oswModelIdx]) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Restore calibration portfolios
    ARM_CalibMethod* volCalibMethod;
    ARM_StdPortfolioPtr oswPortfolio;

    double price,vega,nominal,weight;
    ARM_Swaption* swaption;
    size_t i;
    bool isNotSelected;

    volCalibMethod = calibMethod;


    /// Update datas for q volatility calibration
    /// -----------------------------------------
    oswPortfolio = volCalibMethod->GetPortfolio();
    size_t nbOSW = oswPortfolio->GetSize();
    size_t volParamSize=volCalibMethod->GetCalibParams().size();
    bool isInitVol = isInitVolParam || volParamSize == 0;

    ARM_GP_Vector initTimes(nbOSW);
    ARM_GP_Vector initVols(nbOSW);
    double optMat,swapMat,volATM,strike,swapRate;

    for(i=0;i<nbOSW;++i)
    {
        swaption=static_cast< ARM_Swaption* >(oswPortfolio->GetAsset(i));
	    swaption->SetModel(oswBSModel);
        price=swaption->ComputePrice();
        vega=swaption->ComputeSensitivity(K_VEGA);

        nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(oswPortfolio->GetWeights()))[i];

        isNotSelected = vega < IR_VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        oswPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        oswPortfolio->SetPrecision(0.001*vega,i);
        oswPortfolio->SetPrice(price,i);

        if(isInitVol)
        {
            /// Vol initialisation
            optMat  = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
            swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
            strike  = swaption->GetStrike();
            volATM  = oswBSModel->ComputeVol(optMat,swapMat,strike,strike)/100.0;

            initTimes[i]    = swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
            if(isModelHW)
            {
                swapRate = swaption->CptMarketSwapRate();
                initVols[i] = volATM * swapRate;
            }
            else
                initVols[i] = volATM;
        }
    }

    if(isInitVol)
    {
        /// Replace the sigma param with new initialisations
        ARM_GP_Vector volLowerBound(nbOSW,isModelHW ? HWVOL_LOWER_BOUND : QVOL_LOWER_BOUND);
        ARM_GP_Vector volUpperBound(nbOSW,isModelHW ? HWVOL_UPPER_BOUND : QVOL_UPPER_BOUND);
        ARM_CurveModelParam* vol = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
            "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
        if(volParamSize==0)
            volCalibMethod->GetCalibParams().push_back(vol);
        else if(volParamSize==1)
        {
            delete volCalibMethod->GetCalibParam(0);
            (volCalibMethod->GetCalibParams())[0] = vol;
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : diagonal swaption calibrate only volatilities");
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeFxEquivStrike
///	Returns: a doubke
///	Action : compute an equivalent strike for FX options
/////////////////////////////////////////////////////////////////
double ARM_PRDCalculator::ComputeFxEquivStrike(size_t eventIdx)
{
    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    double spotFx       = fxBSModel->GetFxSpot(); //fxBSModel->GetSpot();
    double spotFxShift  = 0.1;

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	string settlCal  = forex->ComputeSettlementCalendar();
	int settlGap = forex->ComputeSettlementGap();

    ARM_Date fundStartDate( atof(dealDesc.GetElem(eventIdx,StartDate).c_str()) );
    double fundingLeg = basisDomCurve->DiscountPrice(fundStartDate);

    ARM_Date payDate,fxResetDate,fxSetDate;

    ARM_Option* fxOption;

    double fxStrike,price,delta,period,zcPay,zcDomSet,zcForSet,zcFwd,fundingIT,fundingSpread;
    double firstZcFwd=1.0;
    double fxSensi=0.0,strikeLeg=0.0;
    for(size_t i=eventIdx;i<nbEvents;++i)
    {
        /// Compute low strike delta call (by shift because coded BS delta is
        /// not accurate, Fx settlement date not used in d1 of B&S...)
        fxResetDate = ARM_Date( atof(dealDesc.GetElem(i,FxResetDate).c_str()) );
        payDate     = ARM_Date( atof(dealDesc.GetElem(i,CpnPayDate).c_str()) );
        fxStrike    = atof(dealDesc.GetElem(i,FxLowStrike).c_str());
        fxOption    = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
        fxOption->SetPayDate(payDate);
	    fxOption->SetModel(fxBSModel);
        price       = fxOption->ComputePrice();
        fxBSModel->SetFxSpot(spotFx+spotFxShift); //fxBSModel->SetSpot(spotFx+spotFxShift);
	    fxOption->SetModel(fxBSModel);
        delta       = (fxOption->ComputePrice()-price)/spotFxShift;
        fxBSModel->SetFxSpot(spotFx); //fxBSModel->SetSpot(spotFx);

        /// Update fxSensi (taking into account only low strike calls)
        /// Convexity adjustment if settl != pay is assumed to be
        /// the same on PRD swap and strike equivalent FX flow
        period      = delta * atof(dealDesc.GetElem(i,CpnLeverage).c_str());
        zcPay       = basisDomCurve->DiscountPrice(payDate);
        fxSetDate   = fxResetDate;
        fxSetDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));
        zcDomSet    = basisDomCurve->DiscountPrice(fxSetDate);
        zcForSet    = basisForCurve->DiscountPrice(fxSetDate);
        zcFwd       = zcForSet/zcDomSet;
        fxSensi     += period*zcFwd;

        if(i==eventIdx)
            firstZcFwd=zcFwd;

        /// Update low strike leg price
        strikeLeg += fxStrike*period;

        /// Update funding leg price
        payDate         = ARM_Date( atof(dealDesc.GetElem(i,FundingPayDate).c_str()) );
		fundingIT       = (*(itsStructDateStrip->GetInterestTerms()))[itsvCpnIndex[i-1]];

        fundingSpread   = atof(dealDesc.GetElem(i,FundingSpread).c_str());
		fundingSpread  =  itsvFundSpread[itsvFundIndex[i-1]];
        zcPay           = basisDomCurve->DiscountPrice(payDate);
        fundingLeg      += fundingSpread*fundingIT*zcPay;

        /// Free memory
        delete fxOption;

    }

    /// Terminal condition flow
    if(itsRedemptionType != ARM_PRCSRedemptionType::standard)
    {
        payDate   = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,RedemptionPayDate).c_str()) );
        fxResetDate = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,RedemptionResetDate).c_str()) );
        fxStrike    = atof(dealDesc.GetElem(nbEvents-1,RedemptionStrike).c_str());
        fxSetDate   = fxResetDate;
        fxSetDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));
        zcDomSet    = basisDomCurve->DiscountPrice(fxSetDate);
        zcForSet    = basisForCurve->DiscountPrice(fxSetDate);
        zcPay       = basisDomCurve->DiscountPrice(payDate);
        zcFwd       = zcForSet/zcDomSet;
        if(itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)
        {
            /// Mandatory termination => Fwd/K of the nominal is paid back
            /// <=> deltaPut=-1 in dual option case
            fxSensi     += zcFwd*zcPay/fxStrike;
            strikeLeg   -= -zcPay;
        }
        else
        {
            /// Dual option termination : compute delta put
            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
            fxOption->SetPayDate(payDate);
	        fxOption->SetModel(fxBSModel);
            price       = fxOption->ComputePrice();
            fxBSModel->SetFxSpot(spotFx+spotFxShift); //fxBSModel->SetSpot(spotFx+spotFxShift);
	        fxOption->SetModel(fxBSModel);
            delta       = (fxOption->ComputePrice()-price)/spotFxShift;
            fxBSModel->SetFxSpot(spotFx); //fxBSModel->SetSpot(spotFx);

            fxSensi     -= delta*zcFwd/fxStrike;
            strikeLeg   -= delta;

            /// Free memory
            delete fxOption;
        }

    }

    fxSensi     /= firstZcFwd;
    fundingLeg  -= zcPay;

    return (strikeLeg+fundingLeg)/fxSensi;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_PRDCalculator::CreateFxOption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    
    ARM_Option* fxOption;
    list < ARM_Security* > fxOptionList;

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	
	ARM_GP_Vector* exerResetDates = itsExerciseDateStrip->GetResetDates();
    if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
    {
        /// Portfolio = market standard FX options (ATM)

		/// Compute standard vol expiries
		ARM_BSModel* bsModel = fxBSModel->CreateFxBSModel();
		ARM_VolCurve* fxBSVol = bsModel->GetVolatility();
		size_t nbExp = fxBSVol->GetExpiryTerms()->GetSize();
		ARM_GP_Vector fxStdExp(nbExp);
		for(int i=0;i<nbExp;++i)
			fxStdExp[i] = asOfDate + K_YEAR_LEN * (*(fxBSVol->GetExpiryTerms()))[i];
		delete bsModel;

        /// Check if the 1st notice date may replace the 1st standard fx option expiry
        ARM_Date firstNoticeDate;
        size_t fxExpIdx=0;
        for(size_t exerIdx=0; exerIdx<itsExerSize; ++exerIdx)
        {
            if(itsvIsExerDate[exerIdx] && !itsvCpnIsFixed [exerIdx])
            {
                firstNoticeDate = ARM_Date( (*itsExerciseDateStrip->GetResetDates())[exerIdx]);
                if(firstNoticeDate.GetJulian() >= fxStdExp[fxExpIdx]-K_NEW_DOUBLE_TOL)
                {
                    /// Replace the 1st standard expiry by the 1st notice
                    fxOption    = new ARM_Option(forex,firstNoticeDate,K_MARKET_RATE,K_CALL,K_EUROPEAN);
                    fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
                    ++fxExpIdx;
                    break;
                }
            }
        }

        /// Create a fx option for each following standard expiry greater than the 1st notice
        /// but lower or equal to 30y
        double maxFxOptionExpiry = asOfDate + 11111; /// 11111/365 is > 30y!!
        for(;fxExpIdx < nbExp && fxStdExp[fxExpIdx] < maxFxOptionExpiry;++fxExpIdx)
        {
            if(fxStdExp[fxExpIdx] > firstNoticeDate.GetJulian() + K_NEW_DOUBLE_TOL)
            {
                fxOption    = new ARM_Option(forex,ARM_Date(fxStdExp[fxExpIdx]),K_MARKET_RATE,K_CALL,K_EUROPEAN);
                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }
    else
    {
        /// Portfolio = FX options of the PRD leg, calibration strike depends on calibType
        ARM_Date fxResetDate,fxPayDate;
        double fxStrike;
        size_t strikeIdx=0,nbStrikes = itsCalibDatas.size();
        for(size_t exerIdx=0; exerIdx<itsExerSize; ++exerIdx)
        {
			if(itsvIsExerDate[exerIdx] && !itsvCpnIsFixed [exerIdx])
            //if(itsvCpnIsFloored[itsvCpnIndex[exerIdx]])
            {
                /// Here is a Fx option for a floored coupon
				fxResetDate = ARM_Date((*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[exerIdx]]);
                fxPayDate   =ARM_Date((*(itsStructDateStrip->GetPaymentDates()))[itsvCpnIndex[exerIdx]]);
				
                if(itsCalibType == ARM_PRCSCalibTypes::ATSFxEquivCalib)
                {
                    fxStrike = ComputeFxEquivStrike(exerIdx+1);
                }
                else if(itsCalibType == ARM_PRCSCalibTypes::ATSFxProfileCalib)
                {
                    fxStrike = itsCalibDatas[strikeIdx];
                    if(strikeIdx+1 < nbStrikes)
                        ++strikeIdx;
                }
                else
                {
					fxStrike    = itsvLowStrikeCpn[itsvCpnIndex[exerIdx]];
                    if(itsvCpnIsCapped[itsvCpnIndex[exerIdx]])
                        /// Here is a Fx option for a capped coupon
                        fxStrike = 0.5*( fxStrike + itsvHighStrikeCpn[itsvCpnIndex[exerIdx]]);
                }

                fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
                fxOption->SetPayDate(fxPayDate);

                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }


    ARM_StdPortfolio* fxPort = new ARM_StdPortfolio(fxOptionList);
    for(int i=0;i<fxPort->size();++i)
    {
        fxPort->SetWeight(FX_DEFAULT_WEIGHT,i);
        fxPort->SetPrice((i+1)*FX_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(fxPort);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateFxCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRDCalculator::CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;


    /// Build an empty bootstrap calib method (filled latter for fx volatility calibration)
    ARM_CalibMethod* fxCalib = new ARM_CalibMethod(fxOptionPF,emptyCalibParam,
                                    ARM_CalibMethodType::Bootstrap1D,ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,noPreviousMethod,
                                    isCalibShared,modelIdx);

    return fxCalib;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeFxOptionPrice
///	Returns: nothing
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, bool isFreezeWeights, bool isInitVolParam)
{
    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    ARM_StdPortfolioPtr fxPortfolio = calibMethod->GetPortfolio();
    size_t nbFx = fxPortfolio->GetSize();
    size_t volParamSize=calibMethod->GetCalibParams().size();
    bool isInitVol = isInitVolParam || volParamSize == 0;

    ARM_GP_Vector initTimes(nbFx);
    ARM_GP_Vector initVols(nbFx);

    ARM_Option* fxOption;
    double optTime,vol,strike,price,vega=0.0,nominal,weight,fwd, atmVol;
    bool isNotSelected;
    size_t strikeIdx=0;
    for(size_t i(0);i<nbFx;++i)
    {
        fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(i));
	    fxOption->SetModel(fxBSModel);
        price=fxOption->ComputePrice();

        switch(itsCalibType)
        {
        case ARM_PRCSCalibTypes::ATSFxMixedCalib : /// Choose ATM/ATS strike with minimum volatility
            vol = fxOption->GetCalcVol();
            strike = fxOption->GetStrike();
            fxOption->SetStrike(K_MARKET_RATE);
            price=fxOption->ComputePrice();
            atmVol = fxOption->GetCalcVol();
            if(vol < atmVol)
            {
                fxOption->SetStrike(strike);
                price=fxOption->ComputePrice();
            }
            break;

        case ARM_PRCSCalibTypes::ATMDoubleCalib :
            /// strike = fwdFx
            fxOption->SetStrike(fxOption->GetCalcFwd());
            price=fxOption->ComputePrice();
            break;

        case ARM_PRCSCalibTypes::ATSFxMoneynessCalib :
        case ARM_PRCSCalibTypes::ATSFxShiftedCalib :
		case ARM_PRCSCalibTypes::HybridBasketCalib :
            if(itsCalibType == ARM_PRCSCalibTypes::ATSFxMoneynessCalib)                 
				///strike = moneyness * fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] *fxOption->GetCalcFwd());
			else if(itsCalibType == ARM_PRCSCalibTypes::HybridBasketCalib)
				 /// First calibration is ATSFxMoneyness=100% like
                fxOption->SetStrike(fxOption->GetCalcFwd());
            else
                /// strike = shift + fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] + fxOption->GetCalcFwd());

            price=fxOption->ComputePrice();
            if(strikeIdx+1 < itsCalibDatas.size())
                ++strikeIdx;
            break;

        case ARM_PRCSCalibTypes::ATSFxMinVolCalib : /// Choose ATS with mimimum volatility
            if(itsCalibDatas.size()==0 || i <= itsCalibDatas[0])
            {
                double minMoneyness = FX_MIN_MONEYNESS;
                double maxMoneyness = FX_MAX_MONEYNESS;
                double stepMoneyness= FX_STEP_MONEYNESS;
                if(itsCalibDatas.size()>1)
                    minMoneyness = itsCalibDatas[1];
                if(itsCalibDatas.size()>2)
                    maxMoneyness = itsCalibDatas[2];
                if(itsCalibDatas.size()>3)
                    stepMoneyness = itsCalibDatas[3];

                double moneyness = minMoneyness;
                double minVolStrike = fxOption->GetStrike();
                double minVol = fxOption->GetCalcVol();
                while (moneyness <= maxMoneyness)
                {
                    fxOption->SetStrike(moneyness * fxOption->GetCalcFwd());
                    price=fxOption->ComputePrice();
                    if(fxOption->GetCalcVol() < minVol)
                    {
                        minVolStrike = fxOption->GetStrike();
                        minVol = fxOption->GetCalcVol();
                    }
                    moneyness += stepMoneyness;
                }
                fxOption->SetStrike(minVolStrike);
                price=fxOption->ComputePrice();
            }
            break;

		case ARM_PRCSCalibTypes::ATSFxBarrierMoneyness : /// Only for PRDKO, ATS at Barrier
			UpdateFxOption(i,fxOption);
			price=fxOption->ComputePrice();
        }

        if(!isFreezeWeights)
            vega=fxOption->ComputeSensitivity(K_VEGA);

        nominal = fxOption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(fxPortfolio->GetWeights()))[i];

        isNotSelected = vega < FX_VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        fxPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        fxPortfolio->SetPrecision(0.001*vega,i);
        fxPortfolio->SetPrice(price,i);

        if(isInitVol)
        {
            /// Vol initialisation (may be too high because forward FX vol is used)
            optTime = fxOption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();
            strike  = fxOption->GetStrike();
            fwd = fxOption->GetCalcFwd();
            vol = fxOption->GetCalcVol()/100.0;

            initTimes[i] = optTime;
            initVols[i]  = vol;
        }
    }

    if(isInitVol)
    {
        /// Replace the sigma param with new initialisations
        ARM_GP_Vector volLowerBound(nbFx,QVOL_LOWER_BOUND);
        ARM_GP_Vector volUpperBound(nbFx,QVOL_UPPER_BOUND);
        ARM_CurveModelParam* vol = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
            "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
        if(volParamSize==0)
            calibMethod->GetCalibParams().push_back(vol);
        else if(volParamSize==1)
        {
            delete calibMethod->GetCalibParam(0);
            (calibMethod->GetCalibParams())[0] = vol;
        }
        else
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx options calibrate only volatilities at the moment");
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateBasketOption
///	Returns: a portfolio
///	Action : Build a portfolio made of spread between a 
///			 variable notional Fx leg and a IR swap 
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_PRDCalculator::CreateHybridBasketOption()
{
/// Restore models
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
    //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

    ARM_Option*				fxOption;
	ARM_ReferenceValue*		fixRefVal;
	ARM_ReferenceValue*		varRefVal;
	ARM_SwapLeg*			fixLeg;
	ARM_SwapLeg*			varLeg;
    ARM_Swap*				irSwap;
    ARM_Swaption*			irSwaption;
	ARM_StdPortfolio*		irFxBasket;
	ARM_OptionPortfolio*	irFxBasketOption;
	ARM_ExerciseStyle		xStyle;
	ARM_ReferenceValue		xStrike;

    list < ARM_Security* > irFxBasketList;
    list < ARM_Security* > basketList;
	list < ARM_Security* >::iterator begin,pos;


    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t eventIdx,nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
	int flowIdx;


    /// Restore forex
	ARM_Forex* forex	= static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	string settlCal		= forex->ComputeSettlementCalendar();
	int settlGap		= forex->ComputeSettlementGap();

	int nbFunds = nbEvents-1;
	ARM_GP_Vector fundingStarts(nbFunds),fundingEnds(nbFunds),fundingPeriods(nbFunds);
	ARM_GP_Vector fundingMargins(nbFunds),fundingNotios(nbFunds,1.0);


	/// Compute sensitivities to forex (i.e. to FwdFX.ZcDom(Setll)=SptFx.ZcFor(Settl))
	/// and sensivities to domestic Zc at payment dates
    string columnText,fxCallText("CALL");
    ARM_Date eventDate,fxResetDate,fxSetDate,payDate;
    double fxStrike,price,fx,stdDev,d1;
	ARM_GP_Vector fxSensis(nbEvents),zcSensis(nbEvents),fxPrices(nbEvents),cpnPrices(nbEvents);
	ARM_GP_Vector fwdFxs(nbEvents),zcFxSets(nbEvents),zcCpnPays(nbEvents); /// for dump !
	vector< ARM_Date > fxResetDates(nbEvents),fxSetDates(nbEvents); /// for dump !
	vector< ARM_Date > payDates(nbEvents);
	double fxCoef,cpnCoef,zcPay,zcSet;
	double fxLeg=0.0;
	double strikeLeg=0.0;
	double fundingLeg=0.0;
    for(eventIdx=0,flowIdx=0; eventIdx<itsExerSize; ++eventIdx,++flowIdx)
    {
        eventDate = ARM_Date((*itsExerciseDateStrip->GetResetDates())[eventIdx]);

		/// Save funding datas
		fundingStarts[flowIdx]	= (*itsStructDateStrip->GetFlowStartDates())[itsvCpnIndex[eventIdx]]- asOfDate;
		fundingEnds[flowIdx]	= (*itsStructDateStrip->GetFlowEndDates())[itsvCpnIndex[eventIdx]] - asOfDate;
		double term=(*(itsStructDateStrip->GetInterestTerms()))[itsvCpnIndex[eventIdx]];
		fundingPeriods[flowIdx] = term;
		fundingMargins[flowIdx] = itsvFundSpread[itsvFundIndex[eventIdx]];

		/// Compute coupon sensitivities
        payDates[flowIdx] = ARM_Date( (*itsStructDateStrip->GetPaymentDates())[itsvCpnIndex[eventIdx]] );
		cpnCoef = term * itsvLeverage[itsvCpnIndex[eventIdx]];

        if(itsvCpnIsFloored[itsvCpnIndex[eventIdx]])
        {
            /// Here is a Fx option for a floored coupon
            fxResetDate = ARM_Date( (*itsForexDateStrip->GetFlowStartDates())[itsvCpnIndex[eventIdx]]);
            fxStrike    = itsvLowStrikeCpn[itsvCpnIndex[eventIdx]];

			/// Compute dom Zc at settlement date
			fxSetDate   = fxResetDate;
			fxSetDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));

			/// Compute FwdFX.ZcDom(Setll) and ZcDom(Pay) sensitivities
			/// (ZcDom(Pay)/ZcDom(Settl) sensitivity w.r.t. ZcDom(Pay) is neglected)
            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
            fxOption->SetPayDate(payDates[flowIdx]);
			fxOption->SetModel(fxBSModel);
			price	= fxOption->ComputePrice(); // only to compute fwd, vol...
			fx		= fxOption->GetCalcFwd();
			stdDev	= 0.01*fxOption->GetCalcVol()*sqrt(fxOption->GetCalcMat());
			d1		= log(fx/fxStrike)/stdDev + 0.5*stdDev;

			zcPay	= fxOption->GetCalcO1();
			zcSet	= basisDomCurve->DiscountPrice(fxSetDate);
			fxCoef	= zcPay / zcSet;

			fwdFxs[flowIdx]			= fx;
			zcFxSets[flowIdx]		= zcSet;
			fxResetDates[flowIdx]	= fxResetDate;
			fxSetDates[flowIdx]		= fxSetDate;
			zcCpnPays[flowIdx]		= zcPay;


			fxSensis[flowIdx]	= cdfNormal(d1);
			zcSensis[flowIdx]	= -fxStrike*cdfNormal(d1-stdDev);
			cpnPrices[flowIdx]	= price;

            if(itsvCpnIsFloored[itsvCpnIsCapped[eventIdx]])
			{
                /// Here is a Fx option for a capped coupon
				fxStrike    = itsvHighStrikeCpn[itsvCpnIndex[eventIdx]];
				fxOption->SetStrike(fxStrike);
				price	= fxOption->ComputePrice();
				fx		= fxOption->GetCalcFwd();
				stdDev	= 0.01*fxOption->GetCalcVol()*sqrt(fxOption->GetCalcMat());
				d1		= log(fx/fxStrike)/stdDev + 0.5*stdDev;

				fxSensis[flowIdx]	-= cdfNormal(d1);
				zcSensis[flowIdx]	+= fxStrike*cdfNormal(d1-stdDev);
				cpnPrices[flowIdx]	-= price;
			}

			fxSensis[flowIdx] *= fxCoef * cpnCoef;

			/// Save as an option on fwd FX to build the FX portfolio
			fxOption->SetMaturity(eventDate);
            fxOption->SetPayDate(fxSetDate);
            irFxBasketList.push_back(static_cast< ARM_Security* >(fxOption));
        }
		else
		{
			/// Only fixed flow
			zcSensis[flowIdx]	= itsvFixCpn[itsvCpnIndex[eventIdx]];
			zcPay				= basisDomCurve->DiscountPrice(payDate);
			cpnPrices[flowIdx]	= zcPay * zcSensis[flowIdx];
		}

		zcSensis[flowIdx]	*= cpnCoef;
		cpnPrices[flowIdx]	*= cpnCoef;
		fxPrices[flowIdx]	= cpnPrices[flowIdx] - zcPay * zcSensis[flowIdx];

		/// Legs price test
		fxLeg += fxSensis[flowIdx] * fx*basisDomCurve->DiscountPrice(fxSetDate);
		strikeLeg += -zcSensis[flowIdx] * zcPay;
		fundingLeg += basisDomCurve->DiscountPrice(fundingStarts[flowIdx]/K_YEAR_LEN)
					 +(fundingPeriods[flowIdx]*fundingMargins[flowIdx]-1)*
						basisDomCurve->DiscountPrice(fundingEnds[flowIdx]/K_YEAR_LEN);
    }


	/// Case of non standard termination
	flowIdx = itsExerSize-1;
	ARM_Date lastCpnPayDate(payDate);
    if(itsRedemptionType != ARM_PRCSRedemptionType::standard)
    {
		eventDate = ARM_Date((*itsExerciseDateStrip->GetResetDates())[flowIdx]);

        payDates[flowIdx] = ARM_Date( atof(dealDesc.GetElem(flowIdx,RedemptionPayDate).c_str()) );
		if(lastCpnPayDate < payDate) lastCpnPayDate = payDate;
        fxResetDate = itsResetRedemptionDate;
        fxStrike    = itsRedemptionStrike;

        fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
        fxOption->SetPayDate(payDates[flowIdx]);
	    fxOption->SetModel(fxBSModel);
        price	= fxOption->ComputePrice(); // only to compute fwd, vol...
		fx		= fxOption->GetCalcFwd();

        fxSetDate = fxResetDate;
        fxSetDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));

		zcPay	= fxOption->GetCalcO1();
		zcSet	= basisDomCurve->DiscountPrice(fxSetDate);
		fxCoef	= zcPay / zcSet;

		fwdFxs[flowIdx]			= fx;
		zcFxSets[flowIdx]		= zcSet;
		fxResetDates[flowIdx]	= fxResetDate;
		fxSetDates[flowIdx]		= fxSetDate;
		zcCpnPays[flowIdx]		= zcPay;

        if(itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)
        {
            /// Mandatory termination => Fwd/K of the nominal is paid back
            /// <=> deltaPut=-1 in dual option case
            fxSensis[flowIdx]	= 1.0;
            zcSensis[flowIdx]	= -1.0;
            cpnPrices[flowIdx]	= zcPay * (fx-fxStrike);
        }
        else
        {
            /// Dual option termination
            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
            fxOption->SetPayDate(payDate);
	        fxOption->SetModel(fxBSModel);
            price	= fxOption->ComputePrice();
			stdDev	= 0.01*fxOption->GetCalcVol()*sqrt(fxOption->GetCalcMat());
			d1		= log(fx/fxStrike)/stdDev + 0.5*stdDev;

			fxSensis[flowIdx]	= cdfNormal(-d1);
			zcSensis[flowIdx]	= -cdfNormal(-d1+stdDev);
			cpnPrices[flowIdx]	= -price;

        }
        fxSensis[flowIdx]	= fxCoef / fxStrike;
		cpnPrices[flowIdx]	/= fxStrike;
		fxPrices[flowIdx]	= cpnPrices[flowIdx] - zcPay * zcSensis[flowIdx];

		fxOption->SetMaturity(eventDate);
        fxOption->SetPayDate(fxSetDate);
        irFxBasketList.push_back(static_cast< ARM_Security* >(fxOption));

		/// Legs price test
		fxLeg += fxSensis[flowIdx] * fx*basisDomCurve->DiscountPrice(fxSetDate);
		strikeLeg += -zcSensis[flowIdx] * zcPay;

    }
	else
	{
		/// No special termination
        payDates[flowIdx]	= ARM_Date( atof(dealDesc.GetElem(flowIdx,CpnPayDate).c_str()) );
        fxSensis[flowIdx]	= 0.0;
		zcSensis[flowIdx]	= 0.0;
		cpnPrices[flowIdx]	= 0.0;
		fxPrices[flowIdx]	= 0.0;

		fxResetDates[flowIdx]	= fxResetDate;
		fxSetDates[flowIdx]		= fxSetDate;
	}

	/// Get standard domestic curves parameters
	ARM_Currency* ccy			= basisDomCurve->GetCurrencyUnit();
	int spotDays				= ccy->GetSpotDays();
	ARM_INDEX_TYPE indexType	= ccy->GetVanillaIndexType();
	char* resetCal				= ccy->GetResetCalName(indexType);
	char* payCal				= ccy->GetPayCalName(indexType);
	int stdFixFreq				= ccy->GetFixedPayFreq();
	int stdFixDayCount			= ccy->GetFixedDayCount();
	int resetGap				= - spotDays;
	int fwdRule					= K_MOD_FOLLOWING;
    int intRule					= K_ADJUSTED;
    int stubRule				= K_SHORTSTART;
	int resetTiming				= K_ADVANCE;
	int payTiming				= K_ARREARS;

	ARM_STRIPPER stripper;

	ARM_IRIndex irIndex (indexType, stdFixFreq, stdFixFreq, ccy);
	irIndex.SetTerm(stdFixFreq);
	irIndex.SetYearTerm(1.0/stdFixFreq);

	ARM_Date startDate,endDate;
	ARM_Date prdcEndDate = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,FundingEndDate).c_str()) );
	if(prdcEndDate < lastCpnPayDate)
		prdcEndDate = lastCpnPayDate;

	ARM_DateStripPtr stdCurveSched;
	ARM_GP_Vector *startDates,*endDates,*interestTerms;
	ARM_GP_Vector dfs,swaps,stdSwapSensis,stdZcSensis,df0s(nbEvents);
	ARM_GP_Matrix dZc_dSwp;
	size_t nbDates;
	int i,j;
	double startDf,df,O1,fxPrice,swapPrice,varLegPrice,basketPrice,funding0,funding;
	double basketStrike,x,fixRate,vnsStrike;
	ARM_Vector *legDates,varLegNotios,fixLegNotios;
    for(eventIdx=0;eventIdx<nbEvents;++eventIdx)
    {
		/// Compute a standard domestic swap curve schedule as of current event date
		/// that covers the very last payment date
        eventDate = ARM_Date((*itsExerciseDateStrip->GetResetDates())[eventIdx]);
		startDate = eventDate; /// becareful : startDate = ARM_Date(eventDate) gives strange result in Release !!

		startDate.GapBusinessDay(spotDays,resetCal);
		endDate = startDate;
		while(endDate<prdcEndDate)
			endDate.AddMonths( 12/stdFixFreq );


		stdCurveSched = ARM_DateStripPtr( new ARM_DateStrip(startDate,endDate,stdFixFreq,stdFixDayCount,
			resetCal,fwdRule,intRule,stubRule,resetGap,stdFixFreq,GETDEFAULTVALUE,payCal,
			resetTiming,payTiming) );

		startDates		= stdCurveSched->GetFlowStartDates();
		endDates		= stdCurveSched->GetFlowEndDates();
		interestTerms	= stdCurveSched->GetInterestTerms();

		/// Initialise a new stripper based on the previous schedule
		stripper = ARM_STRIPPER(asOfDate,stdCurveSched);

		startDf = basisDomCurve->DiscountPrice( ((*startDates)[0]-asOfDate)/K_YEAR_LEN );
		stripper.setStartDf(startDf);


		nbDates = endDates->size(); 
		dfs.resize(nbDates+1);
		dfs[0]=startDf;
	
		O1=0.0; /// last one will be also the numeraire of the basket
		swaps.resize(nbDates);
		legDates = new ARM_Vector(nbDates);
		dZc_dSwp.resize(nbDates,nbDates);  /// dZc_dSwp(i,j) = dZc(j)/dS(i) <> 0 if i<=j
		for(i=0; i<nbDates; ++i)
		{
			(*legDates)[i]=(*endDates)[i];
			df = basisDomCurve->DiscountPrice(((*endDates)[i]-asOfDate)/K_YEAR_LEN);
			O1 += df * (*interestTerms)[i];
			dfs[i+1] = df;
			swaps[i] = (startDf - df) / O1;
			stripper.setSwapRate (i, swaps[i]);
			for(j=0; j<=i;++j)
				dZc_dSwp(j,i) = df;
		}
		stripper.strip();


		/// Compute non FX reference price
		fxPrice = 0.0;
		basketPrice = 0.0;
		for(flowIdx=eventIdx-1;flowIdx<nbEvents;++flowIdx)
		{
			df0s[flowIdx] = stripper.df(payDates[flowIdx].GetJulian()-asOfDate);
			fxPrice += fxPrices[flowIdx];
			basketPrice += cpnPrices[flowIdx];
		}
		funding0 = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargins);
		basketPrice -= funding0;


		/// Compute non FX flows sensitivities and standard Zc
		/// sensitivities to standard swap rates
		stdSwapSensis.resize(nbDates);
		for (i=0; i<nbDates; ++i)
		{
			stripper.setSwapRate(i, swaps[i] + VNS_SWAP_BUMP);
			stripper.strip();

			stdSwapSensis[i]=0.0;
			for(flowIdx=eventIdx-1;flowIdx<nbEvents;++flowIdx)
			{
				df = stripper.df(payDates[flowIdx].GetJulian()-asOfDate);
				stdSwapSensis[i] += zcSensis[flowIdx] * (df - df0s[flowIdx]);
			}
			funding = stripper.stdFundingLeg(fundingStarts, fundingEnds, fundingPeriods, fundingNotios, fundingMargins);
			stdSwapSensis[i] = (stdSwapSensis[i] - (funding - funding0)) / VNS_SWAP_BUMP;

			for(j=i; j<nbDates; ++j)
			{
				df = stripper.df((*endDates)[j]-asOfDate);
				dZc_dSwp(i,j) = (df - dZc_dSwp(i,j)) / VNS_SWAP_BUMP;
			}

			stripper.setSwapRate(i, swaps[i]);
		}


		/// Convert non FX PV std swap sensis to std Zc sensi by solving
		/// the implied sensitivies equation i.e. find  [C(k) k=1,N] such that :
		/// [dZc(k)/dS(i) i,k=1,N] * [C(k) k=1,N] = [dPV/dS(i) i=1,N]
		/// where [dZc(k)/dS(i)] is an upper triangular matrix
		stdZcSensis.resize(nbDates);
		for(i=nbDates-1;i>=0; --i)
		{
			stdZcSensis[i]=stdSwapSensis[i];
			for(j=i+1;j<nbDates;++j)
				stdZcSensis[i] -= dZc_dSwp(i,j)*stdZcSensis[j];
			stdZcSensis[i] /= dZc_dSwp(i,i);
		}

		/// Convert non FX PV to a variable notional swap w.r.t.the given
		/// numeraire O1 = sum{k=1,N a(k)Zc(k)} and matching std swap sensis
		/// Sensi to first Zc is not fulfilled (could be done by adding another
		/// degre of freedom with a variable notional on the fixed leg)
		fixLegNotios = ARM_Vector(nbDates,1.0);
		fixRefVal = new ARM_ReferenceValue( static_cast<ARM_Vector*>(legDates->Clone()),
										    static_cast<ARM_Vector*>(fixLegNotios.Clone()) );
		fixRefVal->SetCalcMethod(K_STEPUP_RIGHT);

		/// Becareful IR is now opposite in basket price i.e. Basket = FX - IR
		/// with IR = VNSVarLeg - FixRate.O1 then Zc sensis are used with negatively
		swapPrice = fxPrice - basketPrice;
		for (i=0; i<nbDates; ++i)
			stdZcSensis[i] = - stdZcSensis[i];
		fixRate=0.0;
		varLegNotios = ARM_Vector(nbDates);
		for(i=0;i<VNS_STRIKE_MAX_ITER;++i)
		{
			/// Compute floating leg notional
			j = nbDates-1;
			varLegNotios[j] = - stdZcSensis[j] - fixRate * (*interestTerms)[j];
			varLegPrice = varLegNotios[j] * (dfs[j]-dfs[j+1]);
			for(j=nbDates-2;j>=0;--j)
			{
				varLegNotios[j] = varLegNotios[j+1] - stdZcSensis[j] - fixRate * (*interestTerms)[j];
				varLegPrice += varLegNotios[j] * (dfs[j]-dfs[j+1]);
			}

			/// Compute strike to fit PV
			x=fixRate;
			fixRate=(varLegPrice-swapPrice)/O1;
			if(fabs(x-fixRate)<VNS_STRIKE_TOL)
				break;
		}
		if(i>=VNS_STRIKE_MAX_ITER)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't find an equivalent strike for IR/FX basket");

		basketStrike= fixRate; /// such that Basket = FX - VNSVarLeg + basketStrike
		vnsStrike	= fxPrice/O1 + basketStrike;

		/// Build the option on this variable notional swap
		varRefVal = new ARM_ReferenceValue( static_cast<ARM_Vector*>(legDates->Clone()),
											static_cast<ARM_Vector*>(varLegNotios.Clone()) );
		varRefVal->SetCalcMethod(K_STEPUP_RIGHT);

		varLeg = new ARM_SwapLeg( startDate, 
								 endDate, 
								 &irIndex, 
								 K_RCV,
								 0.0, 
								 K_SHORTSTART, 
								 K_COMP_PROP,
								 ccy,
								 ccy->GetLiborIndexDayCount());

		varLeg->SetAmount(varRefVal); // cloned

		fixLeg = new ARM_FixLeg( startDate,
							   endDate, 
							   vnsStrike * 100.0,
							   K_PAY, 
							   stdFixFreq,
							   stdFixDayCount,
							   K_COMP_PROP,
							   K_ARREARS, 
							   K_ADJUSTED,
							   K_SHORTSTART,
							   ccy);
		
		fixLeg->SetAmount(fixRefVal); // cloned
				
		/// Create the variable notional receiver swaption
		/// Swaption expiry is the basket option one (=current notice date)
		/// Swaption strike is the equivalent IR strike = FXSwapRate + basketStrike
		irSwap = new ARM_Swap(fixLeg,varLeg);
		irSwaption = new ARM_Swaption(irSwap,K_RCV,K_EUROPEAN,vnsStrike * 100.0,eventDate);

        irFxBasketList.push_back(static_cast< ARM_Security* >(irSwaption));

		/// Build the portfolio made of a basket on N fwd FX + 1 variable notional swaption
		irFxBasket = new ARM_StdPortfolio(irFxBasketList);
		for(i=0,flowIdx=eventIdx-1;flowIdx<nbEvents;++i,++flowIdx)
			irFxBasket->SetWeight(fxSensis[flowIdx],i);
		irFxBasket->SetWeight(1.0,irFxBasket->size()-1);

		/// Build an option on this portfolio to describe the IR/FX hybrid basket option
		/// Becareful constructor clones portfolio & its assets
		/// Call=max(FX-IR,0) type is set (put=(max(IR-FX,0)))
		xStrike = ARM_ReferenceValue(basketStrike * 100.0);
        xStyle = ARM_ExerciseStyle(eventDate);
		irFxBasketOption = new ARM_OptionPortfolio(irFxBasket,&xStyle,&xStrike,K_CALL);

		basketList.push_back(static_cast< ARM_Security* >(irFxBasketOption));


		/// Prepare next set of FX options by dropping first one and the last
		/// variable notional swaption in the list
		/// The list is cloned because the IR/FX basket will be deleted with
		/// its linked assets
		irFxBasketList.pop_front(); // current FX option
		irFxBasketList.pop_back();	// IR swaption
		for(pos=irFxBasketList.begin();pos!=irFxBasketList.end();++pos)
			*pos = static_cast<ARM_Security*>((*pos)->Clone());


		/// One step forward by dropping first element of
		/// coupon and funding datas
		fundingStarts.erase(fundingStarts.begin());
		fundingEnds.erase(fundingEnds.begin());
		fundingPeriods.erase(fundingPeriods.begin());
		fundingMargins.erase(fundingMargins.begin());
		fundingNotios.erase(fundingNotios.begin());


		/// Free memory
		delete legDates;
		delete varRefVal;
		delete fixRefVal;
		delete varLeg;
		delete fixLeg;
		delete irSwap;
		delete irFxBasket;

	} // for eventIdx


	/// Build the calibration portfolio
	ARM_StdPortfolio* basketPort = new ARM_StdPortfolio(basketList);
	for(i=0;i<basketList.size();++i)
		basketPort->SetWeight(1.0,i);


	/// Free memory
	for(pos=irFxBasketList.begin();pos!=irFxBasketList.end();++pos)
		delete *pos;
	delete resetCal;
	delete payCal;


    return ARM_StdPortfolioPtr(basketPort);

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRSCCalculator
///	Routine: ComputeHybridBasketPrices
///	Returns: nothing
///	Action : Compute prices of hybrid baskets of the input portfolio
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeHybridBasketPrices()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    size_t volParamSize = itsHybridBasketCalibMethod->GetCalibParams().size();

	/// Restore portfolio of IR/FX hybrid swaption
	ARM_Option* option;
	ARM_Swaption* swaption;
	ARM_OptionPortfolio* optPort;

	ARM_StdPortfolioPtr hybridBasketPF = itsHybridBasketCalibMethod->GetPortfolio();
	size_t i,nbBaskets = hybridBasketPF->size();

	/// Convert to GP IR/FX swaptions
	vector< ARM_VanillaIrFxSwaption* > hybridBaskets(nbBaskets);
    for(i=0;i<nbBaskets;++i)
		hybridBaskets[i] = static_cast<ARM_VanillaIrFxSwaption*>( ARM_ConverterFromKernel::ConvertSecuritytoArgObject(hybridBasketPF->GetAsset(i),asOfDate) );

	/// Idx = -1 => spot residual FX option used
	/// Idx = 0  => 1st residual FX option used
	/// Idx = 1  => 2nd residual FX option used etc...
	itsMktHybridModel->SetRefEqFxOptionIdx(-1);

	/// The previously calibrated model may be used
	/// to compute implicit correlations and replace
	/// the hybrid market model ones
	ComputeHybridBasketCorrelations(0.0,hybridBaskets);
	
	/// Restore previously calibrated FX model
	const ARM_ModelNameMap& modelMap = * (static_cast< ARM_2IRFXModel* > (&*GetPricingModel()))->GetModelMap();	
    ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[ARM_2IRFXModel::FxModel]->Model() );
    const ARM_Curve& fxVol = * (static_cast< const ARM_CurveModelParam& >(fxModel.GetModelParams()->GetModelParam(ARM_ModelParamType::QVol))).GetCurve();


    ARM_GP_Vector initTimes(nbBaskets);
    ARM_GP_Vector initVols(nbBaskets);
	ARM_GP_T_Vector< ARM_GP_VectorPtr > correls(nbBaskets);

	/// Reset market hybrid model to lognormal rates, reset limitation strike status
	/// and set reference FX option to the 1st one in the current strip
	itsMktHybridModel->SetIsLogNorRates(true); /// Start with a LN for VNX & VNS rates
	itsMktHybridModel->ResetIrStrikeStatus();

	double price;
    for(i=0;i<nbBaskets;++i)
    {
		price = hybridBaskets[i]->Price(itsMktHybridModel);
		hybridBasketPF->SetPrice (price, i);

		/// Save computed strike for further toString
		optPort = static_cast<ARM_OptionPortfolio*>(hybridBasketPF->GetAsset(i));
		option = static_cast<ARM_Option*>(optPort->GetPtf()->GetAsset(0));
		option->SetStrike(itsMktHybridModel->GetEqFxStrike());
		swaption = static_cast<ARM_Swaption*>(optPort->GetPtf()->GetAsset(optPort->GetPtf()->size()-1));
		swaption->SetStrike(itsMktHybridModel->GetIrStrike());

        /// Vol is initialised by the first calibration of FX vol
        initTimes[i] = hybridBaskets[i]->GetExpiry();
        initVols[i]  = fxVol.Interpolate(initTimes[i]);
	}

	/// Check limitation strike status for VNS
	if( itsMktHybridModel->IsLogNorRates() &&
		itsMktHybridModel->GetIrStrikeStatus() > 0.2*nbBaskets)
	{
		/// Due to lognormal assumption 20% of VNS in hybrid basket option
		/// led to strike limitation => switch to normal assumption and price again
		itsMktHybridModel->SetIsLogNorRates(false);
		itsMktHybridModel->ResetIrStrikeStatus();
		for(i=0;i<nbBaskets;++i)
		{
			price = hybridBaskets[i]->Price(itsMktHybridModel);
			hybridBasketPF->SetPrice (price, i);

			/// Save computed strike for further toString
			optPort = static_cast<ARM_OptionPortfolio*>(hybridBasketPF->GetAsset(i));
			option = static_cast<ARM_Option*>(optPort->GetPtf()->GetAsset(0));
			option->SetStrike(itsMktHybridModel->GetEqFxStrike());
			swaption = static_cast<ARM_Swaption*>(optPort->GetPtf()->GetAsset(optPort->GetPtf()->size()-1));
			swaption->SetStrike(itsMktHybridModel->GetIrStrike());
		}
	}

	/// Free memory
    DeletePointorVector<ARM_VanillaIrFxSwaption>(hybridBaskets);

    /// Build a new calibParam for spot FX vols
    ARM_GP_Vector volLowerBound(nbBaskets,QVOL_LOWER_BOUND);
    ARM_GP_Vector volUpperBound(nbBaskets,QVOL_UPPER_BOUND);
    ARM_CurveModelParam* vol = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
        "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
    if(volParamSize==0)
        itsHybridBasketCalibMethod->GetCalibParams().push_back(vol);
    else if(volParamSize==1)
    {
        delete itsHybridBasketCalibMethod->GetCalibParam(0);
        (itsHybridBasketCalibMethod->GetCalibParams())[0] = vol;
    }
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME +
			" :  IR/FX hybrid swaptions only calibrate spot FX volatilities");

	/// Reset the Fx vols in the model (else vol schedule will be merged with initTimes[]
	fxModel.GetModelParams()->DeleteModelParam(ARM_ModelParamType::QVol);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRSCCalculator
///	Routine: ComputeHybridBasketCorrelations
///	Returns: nothing 
///	Action : Compute implied correlation between FX and IR underlyings
///			 If domestic variable notional swap are assumed to be
///			 market instrument at a given expiry the hybrid basket PV is :
///
///				Basket(T) = FixO1dom(T)*[EqFxStrip(T)/FixO1dom(T) - (VNSRate(T)-Strike)]
///						  = FixO1dom(T)*[VNXRate(T) - (VNSRate(T)-Strike)]
///						  = FixO1dom(T)*[X(T,T0).Bdom(T,T0)/Bfor(T,T0).O1for(T)/FixO1dom(T) - VNSRate(T)-(-Strike)]
///
///			 If VNS is not a market instrument, domestic IR part is
///			 a simple Zc strip and its volatility/correlations are estimated
///			 by an ATM calibrated model. Hybrid basket PV is :
///
///				Basket(T) = Bdom(T,T0)*[EqFxStrip(T)/Bdom(T,T0) - DomZcStrip(T)/Bdom(T,T0)]
///						  = Bdom(T,T0)*[X(T,T0).O1for(T)/Bfor(T,T0) - O1dom(T)/Bdom(T,T0)]
///
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeHybridBasketCorrelations(double evalTime,
														 const vector< ARM_VanillaIrFxSwaption* >& hybridBaskets)
{
    ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());
	const ARM_ModelNameMap& modelMap = * (hybridModel->GetModelMap());

	ARM_PricingModelPtr domZcModel = modelMap[ARM_2IRFXModel::DomModel]->Model();
    const ARM_ModelParamsQ1F* const domZcParams   = static_cast< const ARM_ModelParamsQ1F* const >(domZcModel->GetModelParams());

	ARM_PricingModelPtr forZcModel = modelMap[ARM_2IRFXModel::ForModel]->Model();
    const ARM_ModelParamsQ1F* const forZcParams   = static_cast< const ARM_ModelParamsQ1F* const >(forZcModel->GetModelParams());

	ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[ARM_2IRFXModel::FxModel]->Model() );
    const ARM_ModelParamsQ1F_Fx& fxParams = static_cast< const ARM_ModelParamsQ1F_Fx&>(*(fxModel.GetModelParams()));
    const ARM_CurveModelParam& fxVolParam = static_cast< const ARM_CurveModelParam& >(fxParams.GetModelParam(ARM_ModelParamType::QVol));

	ARM_ZeroCurve* domZcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* forZcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));

	/// Restore model correlations
    ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey]));
	double domForCor	= (*correlMatrix)(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel);
	double domFxCor		= (*correlMatrix)(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
	double forFxCor		= (*correlMatrix)(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel);

	size_t idx,nbBaskets = hybridBaskets.size();
	size_t i,j,nbForCoefs,nbFloatFlows,nbFixFlows;
	ARM_GP_Vector domFloatZcCoefs,domFloatZcTimes,domFixZcCoefs,forZcCoefs;
	double expiryTime,expiryYf,refFxOptionSetTime;
	double forO1,domFloatO1,domFixO1,zcForSet,zcDomSet;
	double forBetaSum,forO1Var,forO1StdDev,domFloatBetaSum;
	double fwdFxForO1Covar;
	double fwdFxVar,fwdFxStdDev;
	double domFixBetaSum,domFixO1Var,domFixO1StdDev,domForFixO1Covar,fwdFxDomFixO1Covar;
	double vnxRateStdDev,vnsRateStdDev,vnxVnsRateCovar;
	double sfsf,sdsd,sdsf,xsf,sfzf,sfzd,xsd,sdzf,sdzd,xx,xzf,xzd,zfzf,zfzd,zdzd;
	double fwdFx,zcPay,beta,invSqrtt;


	ARM_GP_T_Vector< ARM_VectorPtr > hybridDatasStruct(nbBaskets);
	ARM_GP_Vector* hybridDatas;

	ARM_VanillaIrFxSwaption* irFxSwaption;
	for(idx=0;idx<nbBaskets;++idx)
	{
		irFxSwaption	= hybridBaskets[idx];

		expiryTime	= irFxSwaption->GetExpiry();
		expiryYf	= expiryTime/K_YEAR_LEN;
		invSqrtt	= 1.0/sqrt(expiryYf);

		/// Reference Fx option is the first one of the current strip (previously is was
		/// the FX option which expiries at hybrid option expiry i.e. current notice)


		/// Compute foreign Zc coefs
		const ARM_GP_Vector& fxNotionals = irFxSwaption->GetFxNominals();
		const ARM_GP_Vector& fxSettlementTimes = irFxSwaption->GetFxSettlementTimes();
		const ARM_GP_Vector& fxPayTimes = irFxSwaption->GetFxPayTimes();
		nbForCoefs = fxNotionals.size();

		if(itsMktHybridModel->GetRefEqFxOptionIdx()<0)
		{
			ARM_Date expirySetDate(fxParams.ComputeSettlementDate(domZcCurve,forZcCurve,domZcCurve->GetAsOfDate()+expiryTime));
			refFxOptionSetTime = expirySetDate.GetJulian()-domZcCurve->GetAsOfDate().GetJulian();
		}
		else
		{
			int refFxIdx = itsMktHybridModel->GetRefEqFxOptionIdx();
			refFxOptionSetTime = fxSettlementTimes[refFxIdx < nbForCoefs ? refFxIdx : nbForCoefs-1];
		}


		forZcCoefs.resize(nbForCoefs);

		forO1 = 0;
		for(i=0;i<nbForCoefs;++i)
		{
			fwdFx = fxParams.Forward(fxSettlementTimes[i]);
			zcPay = domZcCurve->DiscountPrice(fxPayTimes[i]/K_YEAR_LEN);

			zcForSet = forZcCurve->DiscountPrice(fxSettlementTimes[i]/K_YEAR_LEN);
			if(fxSettlementTimes[i] != fxPayTimes[i])
			{
				zcDomSet		= domZcCurve->DiscountPrice(fxSettlementTimes[i]/K_YEAR_LEN);
				forZcCoefs[i]	= zcForSet*zcPay*fxNotionals[i]/zcDomSet;
			}
			else
				forZcCoefs[i]	= zcForSet*fxNotionals[i];

			forO1	+= forZcCoefs[i];
		}


		/// Foreign fwd O1 fwd variance
		sfsf=forZcParams->StateLocalVariance(evalTime,expiryTime,refFxOptionSetTime);
		forBetaSum=0.0;
		for(i=0;i<nbForCoefs;++i)
		{
			beta		= forZcParams->BetatT(refFxOptionSetTime,fxSettlementTimes[i]) * forZcCoefs[i];
			forBetaSum	+= beta;
		}
		forBetaSum = -forBetaSum/forO1;
		forO1Var = forBetaSum * forBetaSum * sfsf;
		forO1StdDev = sqrt(forO1Var);


		/// Compute domestic Zc coef IR part = 
		///		FloatLeg(t)-FixLeg(t) = O1Fix(t).(VNSRate(t)-K)
		///		 where VNSRate(t) = FloatLeg(t)/O1Fix(t)
		const ARM_GP_Vector& floatNotionals = *(irFxSwaption->GetIrSwaption().GetFloatNotional());
		const ARM_GP_Vector& floatStartTimes = *(irFxSwaption->GetIrSwaption().GetFloatStartTimes());
		const ARM_GP_Vector& floatEndTimes = *(irFxSwaption->GetIrSwaption().GetFloatEndTimes());
		const ARM_GP_Vector& fixPayPeriods = *(irFxSwaption->GetIrSwaption().GetFixPayPeriods());
		const ARM_GP_Vector& fixPayTimes = *(irFxSwaption->GetIrSwaption().GetFixPayTimes());
		const ARM_GP_Vector& fixNotionals = *(irFxSwaption->GetIrSwaption().GetFixNotional());
		nbFixFlows = fixNotionals.size();
		nbFloatFlows = nbFixFlows + 1;

		/// Fixed leg strip
		domFixZcCoefs.resize(nbFixFlows);
		domFixO1=0.0;
		for(i=0;i<nbFixFlows;++i)
		{
			domFixZcCoefs[i] = fixPayPeriods[i]*fixNotionals[i]*
							domZcCurve->DiscountPrice(fixPayTimes[i]/K_YEAR_LEN);

			domFixO1 += domFixZcCoefs[i];
		}

		/// Domestic Fix O1 variance & std dev
		sdsd=domZcParams->StateLocalVariance(evalTime,expiryTime,refFxOptionSetTime);
		domFixBetaSum=0.0;
		for(i=0;i<nbFixFlows;++i)
		{
			beta		= domZcParams->BetatT(refFxOptionSetTime,fixPayTimes[i]) * domFixZcCoefs[i];
			domFixBetaSum	+= beta;
		}
		domFixBetaSum = -domFixBetaSum/domFixO1;
		domFixO1Var = domFixBetaSum * domFixBetaSum * sdsd;
		domFixO1StdDev = sqrt(domFixO1Var);

		/// Float leg strip
		domFloatZcTimes.resize(nbFloatFlows);
		domFloatZcCoefs.resize(nbFloatFlows);
		domFloatZcTimes[0] = floatStartTimes[0];
		domFloatZcCoefs[0] = floatNotionals[0] * domZcCurve->DiscountPrice(domFloatZcTimes[0]/K_YEAR_LEN);
		domFloatO1 = domFloatZcCoefs[0];
		for(j=0,i=1;i<nbFixFlows;++i,++j)
		{
			domFloatZcTimes[i] = floatEndTimes[j]; // should be = fixPayTimes[j]
			domFloatZcCoefs[i] = (floatNotionals[i] - floatNotionals[j])*
							domZcCurve->DiscountPrice(domFloatZcTimes[i]/K_YEAR_LEN);

			domFloatO1	+= domFloatZcCoefs[i];
		}
		j=nbFixFlows-1;
		domFloatZcTimes[nbFixFlows] = floatEndTimes[j];
		domFloatZcCoefs[nbFixFlows] = -floatNotionals[j]*
								 domZcCurve->DiscountPrice(domFloatZcTimes[nbFixFlows]/K_YEAR_LEN);
		domFloatO1	+= domFloatZcCoefs[nbFixFlows];

		/// Domestic Float O1 variance beta coefficient
		domFloatBetaSum=0.0;
		for(i=0;i<nbFloatFlows;++i)
		{
			beta		= domZcParams->BetatT(refFxOptionSetTime,domFloatZcTimes[i]) * domFloatZcCoefs[i];
			domFloatBetaSum	+= beta;
		}
		domFloatBetaSum = -domFloatBetaSum/domFloatO1;


		/// Compute model lognormal correlations between :
		///		O1for(t) & FixO1dom(t)
		///		X(t,T0) & O1for(t)
		///		X(t,T0) & FixO1dom(t)
		/// and lognormal or normal correlations between
		///		VNXRate(t)  & VNSRate(t)
		///		where VNXRate(t)=X(t,T0)*O1for(t)/FixO1dom(t) and VNSRate(t)=VNS(t)/FixO1dom(t)
		///
		/// In addition model lognormal volatilities are saved for :
		///		FixO1Dom(t)
		///		O1For(t)
		/// and normal volatility for :
		///		VNSRate(t)

		hybridDatas = new ARM_GP_Vector(ARM_MarketHybridModel::NbDatas);
		(*hybridDatas)[ARM_MarketHybridModel::Time] = expiryTime;


		/// Domestic Fix O1 / foreign O1 covariance
		sdsf = ARM_ModelParamsHW1F::HW1FStateCovariance(domZcParams,forZcParams,evalTime,expiryTime,refFxOptionSetTime);
		domForFixO1Covar = domForCor * forBetaSum * domFixBetaSum * sdsf;

		/// Fwd Fx(expirySet) / foreign O1 Zc covariance
		xsf		= ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,forZcParams,evalTime,expiryTime,refFxOptionSetTime);
		sfzf	= ARM_ModelParamsHW1F::HW1FStateZcCovariance(forZcParams,forZcParams,evalTime,expiryTime,refFxOptionSetTime,refFxOptionSetTime);
		sfzd	= ARM_ModelParamsHW1F::HW1FStateZcCovariance(forZcParams,domZcParams,evalTime,expiryTime,refFxOptionSetTime,refFxOptionSetTime);
		fwdFxForO1Covar = (forFxCor*xsf + sfzf - domForCor*sfzd) * forBetaSum;

		/// Fwd Fx((expirySet) / domestic Fix O1 covariance
		xsd		= ARM_ModelParamsHW1F::HW1FEqFxStateCovariance(fxVolParam,domZcParams,evalTime,expiryTime,refFxOptionSetTime);
		sdzf	= ARM_ModelParamsHW1F::HW1FStateZcCovariance(domZcParams,forZcParams,evalTime,expiryTime,refFxOptionSetTime,refFxOptionSetTime);
		sdzd	= ARM_ModelParamsHW1F::HW1FStateZcCovariance(domZcParams,domZcParams,evalTime,expiryTime,refFxOptionSetTime,refFxOptionSetTime);
		fwdFxDomFixO1Covar = (domFxCor*xsd + domForCor*sdzf - sdzd) * domFixBetaSum;

		/// Fwd Fx(expirySet) variance
		xx	= fxParams.StateLocalVariance(evalTime,expiryTime,expiryTime);
		xzf	= ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,forZcParams,evalTime,expiryTime,refFxOptionSetTime);
		xzd	= ARM_ModelParamsHW1F::HW1FEqFxZcCovariance(fxVolParam,domZcParams,evalTime,expiryTime,refFxOptionSetTime);
		zfzf = ARM_ModelParamsHW1F::HW1FZcCovariance(forZcParams,forZcParams,evalTime,expiryTime,refFxOptionSetTime);
		zfzd = ARM_ModelParamsHW1F::HW1FZcCovariance(forZcParams,domZcParams,evalTime,expiryTime,refFxOptionSetTime);
		zdzd = ARM_ModelParamsHW1F::HW1FZcCovariance(domZcParams,domZcParams,evalTime,expiryTime,refFxOptionSetTime);

		fwdFxVar = xx + zfzf + zdzd + 2*(xzf*forFxCor - xzd*domFxCor - zfzd*domForCor);
		fwdFxStdDev = sqrt(fwdFxVar);
		

		/// VNX rate std dev
		vnxRateStdDev = xx + (zfzf+forO1Var+2*forBetaSum*sfzf) + (zdzd+domFixO1Var+2*domFixBetaSum*sdzd)
						+ 2*( forFxCor*(xzf + forBetaSum*xsf) - domFxCor*(xzd + domFixBetaSum*xsd)
							  - domForCor*(zfzd+sdzf+sfzd+sdsf) );
		vnxRateStdDev = sqrt(vnxRateStdDev);

		///  VNSRate std dev
		vnsRateStdDev = fabs(domFloatBetaSum - domFixBetaSum)*sqrt(sdsd);

		/// VNX rate / VNS rate covariance
		vnxVnsRateCovar = domFxCor*xsd + domForCor*(sdzf+forBetaSum*sdsf)
						  - (sdzd+domFixBetaSum*sdsd);
		vnxVnsRateCovar *= (domFloatBetaSum - domFixBetaSum);


		/// Lognormal or normal correlations are identical (due to domestic & foreign
		/// curve approximations used to compute variances and covariances)
		(*hybridDatas)[ARM_MarketHybridModel::DomForCor]	= domForFixO1Covar / (domFixO1StdDev*forO1StdDev);
		(*hybridDatas)[ARM_MarketHybridModel::DomFxCor]		= fwdFxDomFixO1Covar / (fwdFxStdDev*domFixO1StdDev);
		(*hybridDatas)[ARM_MarketHybridModel::ForFxCor]		= fwdFxForO1Covar / (fwdFxStdDev*forO1StdDev);
		(*hybridDatas)[ARM_MarketHybridModel::IrFxSpreadCor]= vnxVnsRateCovar / (vnxRateStdDev*vnsRateStdDev);

		/// Fix O1dom(t), O1for(t) & VNSRate(t) volatilities
		(*hybridDatas)[ARM_MarketHybridModel::ForO1LogNorVol]	= forO1StdDev * invSqrtt;
		(*hybridDatas)[ARM_MarketHybridModel::DomO1LogNorVol]	= domFixO1StdDev * invSqrtt;
		(*hybridDatas)[ARM_MarketHybridModel::IrNorVol]			= fabs(domFloatO1/domFixO1)*vnsRateStdDev * invSqrtt;


		/// Save hybrid datas
		hybridDatasStruct[idx] = ARM_VectorPtr(hybridDatas);
	}

	/// Update hybrid datas structure (=time dependent hybrid datas) in the hybrid market model
	itsMktHybridModel->SetHybridDatas(hybridDatasStruct);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateLocalFxOption
///	Returns: nothing
///	Action : Build portfolios used for local Fx model calibration
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CreateLocalFxOption()
{
    ARM_Option* fxOption;

    list < ARM_Security* > fxFlooredOptionList;
    list < ARM_Security* > fxCappedOptionList;
    list < ARM_Security* > fxRedemptionOptionList;

    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );

    /// Create Fx options for floored & capped coupon
    ARM_Date fxResetDate,fxPayDate;
    double fxStrike;
    for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
        if(itsvCpnIsFloored[itsvCpnIndex[eventIdx-1]])
        {
            /// Here is a Fx option for a floored coupon
            //fxResetDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,FxResetDate).c_str()) );
			fxResetDate = ARM_Date((*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[eventIdx-1]]);

            //fxPayDate   = ARM_Date( atof(dealDesc.GetElem(eventIdx,CpnPayDate).c_str()) );
			fxPayDate = ARM_Date((*(itsStructDateStrip->GetPaymentDates()))[itsvCpnIndex[eventIdx-1]]);

            //fxStrike    = atof(dealDesc.GetElem(eventIdx,FxLowStrike).c_str());
			fxStrike    = itsvLowStrikeCpn[itsvCpnIndex[eventIdx-1]];

            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
            fxOption->SetPayDate(fxPayDate);

            fxFlooredOptionList.push_back(static_cast< ARM_Security* >(fxOption));

            if(itsvCpnIsCapped[itsvCpnIndex[eventIdx-1]])
            {
                /// Here is a Fx option for a capped coupon
                //fxStrike    = atof(dealDesc.GetElem(eventIdx,FxHighStrike).c_str());

				fxStrike = itsvHighStrikeCpn[itsvCpnIndex[eventIdx-1]];

                fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
                fxOption->SetPayDate(fxPayDate);

                fxCappedOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }

    if(fxFlooredOptionList.size()>0)
    {
        delete itsFlooredFxOptionPF;
        itsFlooredFxOptionPF = new ARM_StdPortfolio(fxFlooredOptionList);
        if(fxCappedOptionList.size()>0)
        {
            delete itsCappedFxOptionPF;
            itsCappedFxOptionPF = new ARM_StdPortfolio(fxCappedOptionList);
        }
    }

    /// Redemption terminal coupon
    if(itsRedemptionType != ARM_PRCSRedemptionType::standard)
    {
		fxPayDate   = itsEndDate;
		fxResetDate = itsResetRedemptionDate;
		/// Mandatory case : only forward is used but an ATM option is built
           /// to compute convexity adjustment through local model
        if(itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)            
            fxStrike = K_MARKET_RATE;
        else
            /// Dual option case
			fxStrike = itsRedemptionStrike;

        fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
        fxOption->SetPayDate(fxPayDate);

        fxRedemptionOptionList.push_back(static_cast< ARM_Security* >(fxOption));
        delete itsRedemptionFxOptionPF;
        itsRedemptionFxOptionPF = new ARM_StdPortfolio(fxRedemptionOptionList);
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeLocalFxOptionPrices
///	Returns: nothing
///	Action : Compute prices of Fx options in portfolios used for
///          local Fx model calibration
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputeLocalFxOptionPrices()
{
	/// To support new FX GP models !
    ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[LocalFxModelKey]) );

    size_t nbFloored=0,nbCapped=0,nbRedemption=0;
    if(itsFlooredFxOptionPF) nbFloored          = itsFlooredFxOptionPF->size();
    if(itsCappedFxOptionPF) nbCapped            = itsCappedFxOptionPF->size();
    if(itsRedemptionFxOptionPF) nbRedemption    = itsRedemptionFxOptionPF->size();
    size_t nbFx = nbFloored + nbCapped + nbRedemption;

    /// Compute market prices of Fx options
    ARM_Option* fxOption;
    double price;
    ARM_StdPortfolio* fxPortfolio = itsFlooredFxOptionPF;
    size_t i,assetIdx=0;
    for(i=0;i<nbFx;++i)
    {
        fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(assetIdx));
	    fxOption->SetModel(fxBSModel);
        price=fxOption->ComputePrice();

        fxPortfolio->SetWeight(1.0,assetIdx); // Not used anyway
        fxPortfolio->SetPrecision(K_NEW_DOUBLE_TOL,assetIdx); // Not used anyway

        fxPortfolio->SetPrice(price,assetIdx);

        if(i+1 == nbFloored+nbCapped)
        {
            fxPortfolio = itsRedemptionFxOptionPF;
            assetIdx=0;
        }
        else if(i+1==nbFloored)
        {
            fxPortfolio = itsCappedFxOptionPF;
            assetIdx=0;
        }
        else
            ++assetIdx;
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CalibrateLocalFxModel
///	Returns: nothing
///	Action : Calibrate local Fx Models using associated Fx option portfolios
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CalibrateLocalFxModel()
{
    if(itsFxLocalModelFlag)
    {
        if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
        {
            /// Compute and save ATM forward FX vols from previously
            /// bootstrapped spot FX and Zc vols. Done only if market
            /// standard expiries are calibrated to avoid interpolation
            /// errors on targets for local FX calibration (but smile
            /// adjustment is interpolated from standard expiries)
            if(itsATMFwdFxVol)
                delete itsATMFwdFxVol;
            itsATMFwdFxVol = ComputeATMFxOptionVols();
        }

	    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

        /// Collect all event dates
        const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
        size_t i,nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
        size_t resetIdx,nbResets=nbEvents-1;
        ARM_GP_Vector evalTimes(nbResets);
        ARM_GP_Vector resetTimes(nbResets,0.0);
        for(i=1;i<nbEvents;++i)
        {
            evalTimes[i-1]      = ARM_Date(atof(dealDesc.GetElem(i,EventDate).c_str())).GetJulian() - asOfDate;
            if(dealDesc.GetElemFormat(i,FxResetDate) != ARM_MISSING_TYPE)
                resetTimes[i-1]  = ARM_Date(atof(dealDesc.GetElem(i,FxResetDate).c_str())).GetJulian() - asOfDate;
        }
        
        ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());
        ARM_Local_Model* localFxModel;
        ARM_Option* fxOption;
        double atmFwdFxvol,yfExpiry,price,resetTime;

        ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[LocalFxModelKey]) );
		double saveStrike,unsedPrice,smileAdjVol;
        if(itsFlooredFxOptionPF)
        {
            /// Calibration of the local Fx model for floored coupon
            localFxModel = static_cast< ARM_Local_Model* >( &* hybridModel->GetModel(ARM_2IRFXModel::FlooredFxLocalModel) );
            localFxModel->ResetModelParams();

            ARM_GP_Vector noticeTimes(0);
            resetIdx=0;

            /// Change target vol (interpolated from smiled fwd FX vol)
            /// to 2IRFX ATM fwd FX vol + market smile adj vol and update
            /// target prices consequently
            for(i=0;i<itsFlooredFxOptionPF->size();++i)
            {
                fxOption    = static_cast< ARM_Option* >(itsFlooredFxOptionPF->GetAsset(i));
                resetTime   = fxOption->GetExpiry().GetJulian() - asOfDate;
                if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
                {
                    /// ATM vol is interpolated from bootstrapped vols but
                    // smile adjustment is interpolated from standard expiries
                    yfExpiry    = resetTime/K_YEAR_LEN;
                    atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());

					/// Get market smile adjustment
                    //fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());

					/// Stuff to support new FX GP models
					saveStrike = fxOption->GetStrike();
					fxOption->SetStrike(fxOption->GetCalcFwd());
					unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
					smileAdjVol = - fxOption->GetCalcVol(); /// ATM vol

					fxOption->SetStrike(saveStrike);
					unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
					smileAdjVol += fxOption->GetCalcVol(); /// Smiled vol

                    fxOption->SetCalcVol(atmFwdFxvol+smileAdjVol);

                    price = fxOption->GetCalcO1()*BS(fxOption->GetCalcFwd(),fxOption->GetStrike(), 
								                     yfExpiry,0.01*fxOption->GetCalcVol(),
                                                     fxOption->IsCall() ? K_CALL : K_PUT );
                    itsFlooredFxOptionPF->SetPrice(price,i);
                }
                else
                    /// No interpolation pb because calibrated options = FX options of the underlying
                    price = (*(itsFlooredFxOptionPF->GetMktPrices()))[i];

                /// Get the relevant notice dates (i.e. that may cancelled the fx option)
                while(resetIdx<nbResets && resetTimes[resetIdx] <= resetTime)
                {
                    noticeTimes.push_back(evalTimes[resetIdx]);
                    ++resetIdx;
                }
                localFxModel->CalibrateLocalModel(*fxOption,price,noticeTimes);

            }
        }

        if(itsCappedFxOptionPF)
        {
            /// Calibration of the local Fx model for capped coupon
            localFxModel = static_cast< ARM_Local_Model* >( &* hybridModel->GetModel(ARM_2IRFXModel::CappedFxLocalModel) );
            localFxModel->ResetModelParams();

            ARM_GP_Vector noticeTimes(0);
            resetIdx=0;

            for(i=0;i<itsCappedFxOptionPF->size();++i)
            {
                fxOption    = static_cast< ARM_Option* >(itsCappedFxOptionPF->GetAsset(i));
                resetTime   = fxOption->GetExpiry().GetJulian() - asOfDate; 
                if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
                {
                    yfExpiry    = resetTime/K_YEAR_LEN;
                    atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());

					/// Get market smile adjustment
                    //fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());

					/// Stuff to support new FX GP models
					saveStrike = fxOption->GetStrike();
					fxOption->SetStrike(fxOption->GetCalcFwd());
					unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
					smileAdjVol = - fxOption->GetCalcVol(); /// ATM vol

					fxOption->SetStrike(saveStrike);
					unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
					smileAdjVol += fxOption->GetCalcVol(); /// Smiled vol

                    fxOption->SetCalcVol(atmFwdFxvol+smileAdjVol);

                    price = fxOption->GetCalcO1()*BS(fxOption->GetCalcFwd(),fxOption->GetStrike(), 
								                     yfExpiry,0.01*fxOption->GetCalcVol(),
                                                     fxOption->IsCall() ? K_CALL : K_PUT );
                    itsCappedFxOptionPF->SetPrice(price,i);
                }
              else
                  price = (*(itsCappedFxOptionPF->GetMktPrices()))[i];

                /// Get the relevant notice dates (i.e. that may cancelled the fx option)
                while(resetIdx<nbResets && resetTimes[resetIdx] <= resetTime)
                {
                    noticeTimes.push_back(evalTimes[resetIdx]);
                    ++resetIdx;
                }
                localFxModel->CalibrateLocalModel(*fxOption,price,noticeTimes);
            }
        }

        if(itsRedemptionFxOptionPF)
        {
            /// Calibration of the local Fx model for terminal redemption coupon
            localFxModel = static_cast< ARM_Local_Model* >( &* hybridModel->GetModel(ARM_2IRFXModel::RedemptionFxLocalModel) );
            localFxModel->ResetModelParams();

            fxOption    = static_cast< ARM_Option* >(itsRedemptionFxOptionPF->GetAsset(0));
            resetTime   = fxOption->GetExpiry().GetJulian() - asOfDate;
			if(itsCalibType == ARM_PRCSCalibTypes::ATMCalib)
            {
                yfExpiry    = resetTime/K_YEAR_LEN;
                atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());
                fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());

                //fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());

				/// Stuff to support new FX GP models
				saveStrike = fxOption->GetStrike();
				fxOption->SetStrike(fxOption->GetCalcFwd());
				unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
				smileAdjVol = - fxOption->GetCalcVol(); /// ATM vol

				fxOption->SetStrike(saveStrike);
				unsedPrice = PriceFxOptionWithFxVolModel(fxOption,fxBSModel);
				smileAdjVol += fxOption->GetCalcVol(); /// Smiled vol

                fxOption->SetCalcVol(atmFwdFxvol+smileAdjVol);

                price = fxOption->GetCalcO1()*BS(fxOption->GetCalcFwd(),fxOption->GetStrike(), 
								                 yfExpiry,0.01*fxOption->GetCalcVol(),
                                                 fxOption->IsCall() ? K_CALL : K_PUT );
                itsRedemptionFxOptionPF->SetPrice(price,0);
            }
          else
              price = (*(itsRedemptionFxOptionPF->GetMktPrices()))[0];

            /// All notice dates may cancelled the fx flow
            localFxModel->CalibrateLocalModel(*fxOption,price,evalTimes);
        }
    }
}

////////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCCalculator
///	Routine: ComputePricingData
///	Returns: 
///	Action : pricing function called from the addin GetPricingData()
////////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_PRDCalculator*>(this)->PriceAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Calibraton
///	Returns: a double
///	Action : to calibrate the 2IRFxModel. At beginning, we run an auto-calibration 
///          by bootstrapping of the volatility on diagonal swaptions
///          on both domestic & foreign market. Secondly, the spot FX volatility
////         is boostrapped on FX options portfolio. finnaly, we calibrat
///          the local models 
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::Calibrate()
{ 
	ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());

    /// Calibrate stochatic models and check errors
    GetCalibMethod()->Calibrate(hybridModel);

	/// to check Calibration to avoid any variance squeeze
    CheckCalibErrors();

    /// Second step if necessary : calibration of PRD Fx options with
    /// bootstrapped ATM Fx vols
    if(itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib)
    {
        //ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
        ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

        /// Save market Fx vols
		//ARM_VolCurve* marketFxVol = static_cast<ARM_VolCurve*>(fxBSModel->GetVolatility());
		ARM_VolCurve* marketFxVol = static_cast<ARM_VolCurve*>(fxBSModel->GetFxVol());

        /// Compute previously bootstrapped ATM vols to give new targets
	    itsATMFwdFxVol = ComputeATMFxOptionVols();
		//fxBSModel->SetVolatility(itsATMFwdFxVol);
		fxBSModel->SetFxVol(itsATMFwdFxVol);
        bool isFreezeWeights=true;
        bool isInitVolParam=true;
        ComputeFxOptionPrices(itsATMDoubleCalibMethod,isFreezeWeights,isInitVolParam);


        /// Replace previous ATM calib method by additional one
        ARM_CalibMethod* atmFxCalib = GetFxCalib();
        SetFxCalib(itsATMDoubleCalibMethod);

        // Calibrate PRD leg Fx options and with bootstrapped ATM vols 
        GetCalibMethod()->Calibrate(hybridModel);

        CheckCalibErrors();

        /// Restore market Fx vols in the market data manager (side effect !)
		//fxBSModel->SetVolatility(marketFxVol);
		fxBSModel->SetFxVol(marketFxVol);

        // Restore ATM calib method
        SetFxCalib(atmFxCalib);

    }

	CalibrateLocalFxModel();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the PRDC option. An auto-calibration is done before
///          by bootstrapping the volatility on diagonal swaptions
///          on both domestic & foreign market then the spot FX volatility
////         is boostrapped on FX options portfolio
/////////////////////////////////////////////////////////////////
double ARM_PRDCalculator::Price()
{
    /// Calibrate
    CalibrateAndTimeIt();
	
	ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());

    ARM_PricingModelPtr initialModel;

	double firstColumnPrice = 0.0;
	if (itsColumnsToPrice.size() > 0)
	{
		if(itsMarkovianDriftSamplerFlag)
			/// Save calibrated model because of interpolation over tree schedule
			initialModel = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(hybridModel->Clone()) );

		ARM_GenPricer* genPricer = new ARM_GenPricer( &*GetGenSecurity(),hybridModel);

		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

			/// URGLY URGLY URGLY mais bon....
		double nominal = itsTradeFromDataBase ? itsvCpnNominal[0] : 1.0;

		firstColumnPrice = genPricer->Price() / nominal;

		string columnName;
		for(size_t i=0;i<itsColumnsToPrice.size();++i)
		{
			columnName = itsColumnsToPrice[i];
			GetPricingData()[columnName] = genPricer->GetPricerInfo()->GetContents(columnName).GetData("Price").GetDouble() / nominal;
		}
		
		/// Restore calibrated model
		if(itsMarkovianDriftSamplerFlag)
			SetPricingModel(initialModel);
	}
	itsHasBeenComputed = true;

    return firstColumnPrice;       
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateUnderlying
///	Returns: 
///	Action : Create a Power Reverse Swap
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::CreateUnderlying()
{
	///
	double dayLag = ARM_GlobalConstant::ARM_SEVENDAYS_LAG;
	ARM_IRIndex irIndex(itsDomesticCcy.GetCcyName(),itsCpnDaycount);

	irIndex.SetResetFrequency(itsCpnFreq);
	irIndex.SetPayFrequency(itsCpnFreq);
	int resetTiming		= itsResetTiming;
	int payTiming		= K_ARREARS;
	irIndex.SetResetTiming(resetTiming);
	irIndex.SetPayTiming(payTiming);

	irIndex.SetCompMeth(K_COMP_PROP);
	irIndex.SetFwdRule(K_MOD_FOLLOWING);

	irIndex.SetYearTerm(1./(double)irIndex.GetResetFrequency());
	irIndex.SetTerm(irIndex.GetResetFrequency());

	int intRule	= K_UNADJUSTED;
	irIndex.SetIntRule(intRule);

	ARM_SwapLeg FxUnderTmpLeg(itsStartDate,
							  itsEndDate,
							  &irIndex,							  						 
							  K_PAY,
							  0.0,
							  itsStubRule,
							  K_COMP_PROP,
							  &itsForeignCcy,
							  itsCpnDaycount,
							  -fabs(itsFxResetGap),
							  const_cast<char*>(itsCpnResetCal.c_str()),
							  const_cast<char*>(itsCpnPayCal.c_str()),
							  1);
	

	ARM_FixLeg* fixUndLeg = new  ARM_FixLeg;
	ARM_AutoCleaner<ARM_FixLeg> HoldFIXUND(fixUndLeg );

	/// foreign coupon from Vector to refrence value
	ARM_GP_Vector* fxResetDates = itsForexDateStrip->GetResetDates();
	ARM_ReferenceValue  rForeignCpn(To_pARM_Vector(&(*fxResetDates- dayLag)), To_pARM_Vector(&(itsForeignCpn*100.0)));
	rForeignCpn.SetCalcMethod(K_STEPUP_LEFT);

	fixUndLeg->SetVarCoupons(&rForeignCpn);

	ARM_SwapLeg* FxUnderLeg = ((ARM_SwapLeg *) fixUndLeg);
	*FxUnderLeg = FxUnderTmpLeg;

	ARM_GP_Vector* payDates = itsStructDateStrip->GetPaymentDates();

	/// domestic nominal by initial Fx to refrence value
	ARM_GP_Vector cpnNominalbyInitialFx = itsvCpnNominal/itsInitialFx;
	ARM_ReferenceValue  rcpnNominalbyInitialFx(To_pARM_Vector(&(*payDates- dayLag) ), To_pARM_Vector(&cpnNominalbyInitialFx));
	rcpnNominalbyInitialFx.SetCalcMethod(K_STEPUP_LEFT);
	FxUnderLeg->SetAmount(&rcpnNominalbyInitialFx);

	resetTiming		= K_ADVANCE;
	irIndex.SetResetTiming(resetTiming);
	ARM_SwapLeg FxNumTmpLeg(itsStartDate,
							itsEndDate,
							&irIndex,						 
							K_RCV,
							0.0,
							itsStubRule,
							K_COMP_PROP,
							&itsDomesticCcy,
							itsCpnDaycount,
							10000,
							const_cast<char*>(itsCpnResetCal.c_str()),
							const_cast<char*>(itsCpnPayCal.c_str()),
							1);
	
	ARM_FixLeg* fixNumLeg = new ARM_FixLeg;
	ARM_AutoCleaner<ARM_FixLeg> HoldFIXNUM(fixNumLeg );

	ARM_GP_Vector* cpnStartDates = itsStructDateStrip->GetFlowStartDates();
	/// domestic coupon from Vector to refrence value
	ARM_ReferenceValue  rDomesticCpn(To_pARM_Vector(&(*cpnStartDates -dayLag)), To_pARM_Vector(&(itsDomesticCpn*100.0)));
	rDomesticCpn.SetCalcMethod(K_STEPUP_LEFT);
	fixNumLeg->SetVarCoupons(&rDomesticCpn);

	ARM_SwapLeg* FxNumLeg = ((ARM_SwapLeg *) fixNumLeg);
	*FxNumLeg = FxNumTmpLeg;

	/// domestic nominal from Vector to refrence value
	ARM_ReferenceValue  rNominal(To_pARM_Vector(&(*payDates- dayLag) ), To_pARM_Vector(&itsvCpnNominal));
	rNominal.SetCalcMethod(K_STEPUP_LEFT);
	FxNumLeg->SetAmount(&rNominal);
	
	/// initial forex from Vector to refrence value
	ARM_ReferenceValue rFx(To_pARM_Vector(&(*fxResetDates- dayLag)), To_pARM_Vector(&itsInitialFx));
	rFx.SetCalcMethod(K_STEPUP_LEFT);
	/// coupon min coupon from Vector to refrence value
	ARM_ReferenceValue rFloor(To_pARM_Vector(&(*fxResetDates- dayLag)), To_pARM_Vector(&(itsMinCpn*100.0)));
	rFloor.SetCalcMethod(K_STEPUP_LEFT);
	///  coupon max from Vector to refrence value
	ARM_ReferenceValue rCap(To_pARM_Vector(&(*fxResetDates- dayLag)), To_pARM_Vector(&(itsMaxCpn*100.0)));
	rCap.SetCalcMethod(K_STEPUP_LEFT);

	/// Funding Leg
	int indexType = !strcmp("EUR",itsDomesticCcy.GetCcyName()) ? FromFrequencyToEuriborType(itsFundFreq) :
																 FromFrequencyToLiborType(itsFundFreq);
	ARM_SwapLeg fundinLeg(itsStartDate,
							itsEndDate,
							(ARM_INDEX_TYPE)indexType,
							K_RCV,
							0.0,
							itsFundFreq,
							itsFundFreq,
							K_ADVANCE,
							K_ARREARS,
							&itsDomesticCcy,
							K_ADJUSTED,
							10000,
							const_cast<char*>(itsCpnResetCal.c_str()),
							const_cast<char*>(itsCpnPayCal.c_str()),
							1,
							K_NX_NONE,
							itsStubRule);

	ARM_GP_Vector* fundStartDates = itsFundDateStrip->GetFlowStartDates() ;
	/// funding spread  from Vector to refrence value
	ARM_ReferenceValue  rspread(To_pARM_Vector(&(*fundStartDates - dayLag)), To_pARM_Vector(&(itsvFundSpread*100.0)));
	rspread.SetCalcMethod(K_STEPUP_LEFT);
	fundinLeg.SetVariableSpread(&rspread);

	ARM_GP_Vector* fundPayDates = itsFundDateStrip->GetPaymentDates() ;
	/// funding nominal from Vector to refrence value
	ARM_ReferenceValue  rfundNotional(To_pARM_Vector(&(*fundPayDates - dayLag)), To_pARM_Vector(&(itsvFundNominal)));
	rfundNotional.SetCalcMethod(K_STEPUP_LEFT);
	fundinLeg.SetAmount(&rfundNotional);
	
	ARM_GP_Vector* exerResetDates = itsExerciseDateStrip->GetResetDates();	
	ARM_Vector* noticeDates = To_pARM_Vector(&ARM_GP_Vector(exerResetDates->begin() + itsNbNoCall,exerResetDates->end()));
	ARM_AutoCleaner<ARM_Vector> HoldNotice(noticeDates );
	ARM_Vector* cancelDates = To_pARM_Vector(&ARM_GP_Vector(cpnStartDates->begin() + itsNbNoCall,cpnStartDates->end()));
	ARM_AutoCleaner<ARM_Vector> HoldCANCEL(cancelDates );

	int dualType = itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption ? K_OPT : 
				   itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption ? K_OPT_FWD : K_OPT_NO;

	ARM_SwapLeg* initFixedLeg = NULL;
	itsPowerReverseSwap = new ARM_PowerReverse(initFixedLeg,
		&fundinLeg,
		FxUnderLeg,
		FxNumLeg,
		noticeDates,
		cancelDates,
		&rFx,
		&rCap,
		&rFloor,
		dualType,
		itsRedemptionStrike,
		itsResetRedemptionDate);

	itsPowerReverseSwap->SetInCentsFactor(false);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateAnalyticalModel
///	Returns: 
///	Action : Create a Black&Scholes Model multi currencies to price
///////////  a power reverse funding in Ccy1 and the structured leg 
///////////  is a forex option on ccy2 et Ccy3
/////////////////////////////////////////////////////////////////
ARM_BS_Model* ARM_PRDCalculator::CreateConvAdjstModel()
{
	ARM_BSModel* DBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );

	/// Convexity model adjustment
	ARM_VolCurve* bsVol = ((ARM_VolLInterpol*)DBSModel->GetVolatility())->ConvertToBSVol(DBSModel->GetZeroCurve());

	ARM_GP_Matrix x3Mat(To_ARM_GP_Matrix(*bsVol->GetVolatilities()));
	ARM_GP_Vector x1Vec(To_ARM_GP_Vector(*bsVol->GetExpiryTerms()));
	ARM_GP_Vector x2Vec(To_ARM_GP_Vector(*(bsVol->GetStrikes())));

	ARM_InterpolType type = ARM_InterpolationType::linear_column_row_extrapoleCst_column_row;

	x1Vec*=365;
	ARM_SurfaceWithInterpol surface( x1Vec, x2Vec, x3Mat,type);

	ARM_SurfaceModelParam modelParam(ARM_ModelParamType::Volatility, &surface);

	ARM_ModelParamVector modelparams(1);
	modelparams[0] = &modelParam;
	ARM_BS_ModelParams BSModelParams(modelparams);

	ARM_BS_Model* model = new ARM_BS_Model( CreateClonedPtr( DBSModel->GetZeroCurve() ), BSModelParams );

	return model;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateAnalyticalModel
///	Returns: 
///	Action : Create a Black&Scholes Model multi currencies to price
///////////  a power reverse funding in Ccy1 and the structured leg 
///////////  is a forex option on ccy2 et Ccy3
/////////////////////////////////////////////////////////////////
ARM_DFBSModel* ARM_PRDCalculator::CreateAnalyticalModel(mdmKeysAlias key,bool useMktDataMgrCorrels)
{
	/// As of Date
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
	/// Black & Scholes Models
	ARM_BSModel* DBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
	ARM_BSModel* FBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );

	/// Correlation
	ARM_GP_Matrix correlMatrix;
	if(useMktDataMgrCorrels)
		correlMatrix = (*dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(GetKeys()[CorrelMatrixKey])));
	else
	{
		ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());
		correlMatrix = hybridModel->GetCorrelMatrix()->Interpolate(0);
	}

	double irCorrelValue	= correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel);
	double domFxcorrelValue	= correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
	double forFxcorrelValue	= correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel);
	ARM_VolFlat dfxCorr(asOfDate,domFxcorrelValue);
	ARM_VolFlat ffxCorr(asOfDate,forFxcorrelValue);

	/// Forex Volatility (compatible with new FX GP models !)
	ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[key]) );
	ARM_VolCurve* fxVol=NULL;
	ARM_EqFxBase* fxVolModel=NULL;
	if(fxBSModel)
		/// No FX GP model possible
		fxVol = fxBSModel->GetVolatility();
	else
	{
		ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[key]) );
		if(fxBSModel)
		{
			fxVol		= fxBSModel->GetFxVol();
			fxVolModel	= dynamic_cast< ARM_EqFxBase* >(fxBSModel->GetFxVolModel()); // may be NULL or contains a FX GP model
		}
		else
		{
			/// A FX GP model is input
			fxVolModel = dynamic_cast< ARM_EqFxBase* >( GetMktDataManager()->GetData(GetKeys()[key]));
			if(!fxVolModel)
				ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + ":  FX model is missing in the market data manager and it shoudl be ARM_EqFxBase");

			ARM_BS_Model* model = CreateConvAdjstModel(); // Don't delete it, will be removed by fxModel
			fxVolModel->SetConvAdjustModel(model);
			fxVolModel->SetRho(domFxcorrelValue);

			fxVol = dynamic_cast< ARM_VolCurve* >( GetMktDataManager()->GetData(GetKeys()[BSFxVol]) );
		}
	}

	/// Zc Curves
	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
	/// discount curves
	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));

	/// Fx Spot
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(GetKeys()[ForexKey]));
	double spoValue = forex->GetMarketPrice();
	
	return  new ARM_DFBSModel(DBSModel,
							 FBSModel,
                             &dfxCorr, 
                             &ffxCorr, 
                             fxVol,
                             irCorrelValue, 
                             basisDomCurve,
                             basisForCurve,
                             spoValue,
							 domCurve,
							 basisDomCurve,
							 1.0,NULL,fxVolModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: ComputeAll
///	Returns: an ARM_Vector
///	Action : call price() function and store results in a vector
/////////////////////////////////////////////////////////////////
ARM_Vector* ARM_PRDCalculator::ComputeAll()
{
	ARM_Vector* results = new ARM_Vector(12, 0.0);
	ARM_VolCurve* originalFXVol = NULL;
	ARM_VolLInterpol* latticeATMFXVol = NULL;
	if(!itsAnalyticalModel)
		itsAnalyticalModel = CreateAnalyticalModel();
	if(!itsPowerReverseSwap)
		CreateUnderlying();	
	
	/// Option price if necessairy
	double optionPrice = 0.0;
	if(!itsHasBeenComputed && itsTradeFromDataBase)
		optionPrice = Price();

	double undAfterFirstNoticeInLattice = 0.0;
	double europeanCall = 0.0;
	if(itsTradeFromDataBase){
		/// get the option to cancel price
		optionPrice = GetPricingData().GetData("PRDCOPTION").GetDouble(); 
		// Get from calculator pricing results
        undAfterFirstNoticeInLattice = GetPricingData().GetData("PRDCFIRSTSWAP").GetDouble() * (-1.0); 
		/// Get first PRD call
		europeanCall = GetPricingData().GetData("PRDCFIRSTEUROPEAN").GetDouble();
	}

	itsPowerReverseSwap->SetModel(itsAnalyticalModel);
	itsPowerReverseSwap->SetFXSmileAdjFlag(true); // to allow market vol smile adj to be saved for further using
	double prdcSwapPrice = itsPowerReverseSwap->ComputePrice();

    /// "Before" prices fully use market vols (unsmiled and smiled)
	double undPriceBeforeFirstNotice           = itsPowerReverseSwap->GetUndPriceBeforeFirstNoticeDate();
	double fxStripTotalPriceBeforeFirstNotice  = itsPowerReverseSwap->GetFxStripPriceBeforeNotice();
	double fundingTotalPriceBeforeFirstNotice  = itsPowerReverseSwap->GetFundingPriceBeforeFirstNotice();
	double FXSmileAdjustmentBeforeNotice       = itsPowerReverseSwap->GetFXSmileAdjBeforeNotice();
	double ExNotionalAtBeginning               = itsPowerReverseSwap->GetExNotionalAtBeginning();
	double NotionalExAtCancelDate              = itsPowerReverseSwap->GetNotionalExAtCancelDate();

    double FXSmileAdjustment;
	if(!itsPowerReverseSwap->IsFXSmileAdjFlag())
        FXSmileAdjustment    = itsPowerReverseSwap->GetFXSmileAdjustment();

	// save FXVol before modify it
	originalFXVol = static_cast<ARM_VolCurve*>(itsAnalyticalModel->GetFxVol()->Clone());
	CC_NS(std,auto_ptr)<ARM_VolCurve> HOLDOriginalFXVol(originalFXVol);
	        
	/// Not already computed then do it !
    ARM_PricingModel* theFxVolModel = itsAnalyticalModel->GetFxVolModel();

    if ((itsCalibType == ARM_PRCSCalibTypes::ATMCalib || itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib)
        || ( theFxVolModel == NULL ))
    {
       if (!itsATMFwdFxVol)      itsATMFwdFxVol  = ComputeATMFxOptionVols();
    }
    else
    {
        delete itsATMFwdFxVol;
        itsATMFwdFxVol = itsPowerReverseSwap->ComputeATMFxVolFromFxVolModel();
    }

    latticeATMFXVol = itsATMFwdFxVol;

    // AS WE ARE GOING TO PRICE WITH ATM FX VOL ===> We don't need an FX Vol Model
    // like Mixture model so we need to set it to NULL
    
    

    if (theFxVolModel)
    {
       theFxVolModel = (ARM_PricingModel *) itsAnalyticalModel->GetFxVolModel()->Clone();
       itsAnalyticalModel->SetFxVolModel(NULL);
    }

	itsAnalyticalModel->SetFxVol(latticeATMFXVol);
	itsPowerReverseSwap->SetModelVariable(NULL);
	itsPowerReverseSwap->SetModel(itsAnalyticalModel);

    /// Compute unsmiled underlying price (using ATM fwd FX vols given by
    /// bootstrapped spot FX & Zc vols)
	double undPrice        = itsPowerReverseSwap->ComputePrice();
	double dualOptionPrice = itsPowerReverseSwap->GetDualValue();
	double fxStripTotalPrice = itsPowerReverseSwap->GetfxStripTotalPrice();
	double fundingTotalPrice = itsPowerReverseSwap->GetFundingTotalPrice();
	double fixPrice = itsPowerReverseSwap->GetFixPrice();
	double floorCallPrice= itsPowerReverseSwap->GetFloorCallPrice();
	double capCallPrice= itsPowerReverseSwap->GetCapCallPrice();

    // Re-Set the fxVol model if necessary

    if (theFxVolModel)
    {
       itsAnalyticalModel->SetFxVolModel(theFxVolModel);
       delete theFxVolModel;
       theFxVolModel = NULL;
    }

    if( itsPowerReverseSwap->IsFXSmileAdjFlag())
    {
        /// Get new smile adjustment using model ATM fwd FX vols
        FXSmileAdjustment = itsPowerReverseSwap->GetFXSmileAdjustment();
        /// "Before" prices still inchanged else inconsistency !
		double inCentsFactor = itsTradeFromDataBase ? 1.0: itsvCpnNominal[0]/1.0e+4;
		double tolerance = 1e-3*inCentsFactor;
		bool isPriceConsistent = fabs(undPriceBeforeFirstNotice-itsPowerReverseSwap->GetUndPriceBeforeFirstNoticeDate()) > tolerance;
		bool isAdjsConsistent = fabs(FXSmileAdjustmentBeforeNotice-itsPowerReverseSwap->GetFXSmileAdjBeforeNotice()) > tolerance;

		/// reset original FXVol
		itsAnalyticalModel->SetFxVol(originalFXVol);
		itsPowerReverseSwap->SetModel(itsAnalyticalModel);

    }

	// Unavailable data
	double totalFXDeltaAfterFirstNotice = -9999;
	double underlyingPrice = fundingTotalPrice-(fxStripTotalPrice - FXSmileAdjustment - dualOptionPrice);
	if(itsTradeFromDataBase){
		results->Elt(0) = undAfterFirstNoticeInLattice;
		results->Elt(1) = optionPrice;
		results->Elt(2) = europeanCall;
		results->Elt(3) = totalFXDeltaAfterFirstNotice;
		results->Elt(4) = undPrice;
		results->Elt(5) = undPriceBeforeFirstNotice;
		results->Elt(6) = FXSmileAdjustment;
		results->Elt(7) = FXSmileAdjustmentBeforeNotice;
		results->Elt(8) = dualOptionPrice;
		results->Elt(9) = fxStripTotalPrice;

		results->Elt(10) = ExNotionalAtBeginning;
		results->Elt(11) = NotionalExAtCancelDate;
	}
	else{
	results->Elt(0) = underlyingPrice;
	results->Elt(1) = fxStripTotalPrice;
	results->Elt(2) = fundingTotalPrice;
	results->Elt(3) = fixPrice;
	results->Elt(4) = floorCallPrice;
	results->Elt(5) = capCallPrice;
	results->Elt(6) = FXSmileAdjustment;
	results->Elt(7) = dualOptionPrice;	

	results->Elt(8) = fxStripTotalPriceBeforeFirstNotice;
	results->Elt(9)  = fundingTotalPriceBeforeFirstNotice;
	results->Elt(10) = FXSmileAdjustmentBeforeNotice;
	}

	return results;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::UpdateModel()
{
	delete itsAnalyticalModel;
	itsAnalyticalModel = NULL;
	/// get zc curve
	ARM_ZeroCurve* zcCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcDomKey]));
    /// set it to pricing model
	GetPricingModel()->SetZeroCurve(CreateClonedPtr( zcCurve ));

	/// this takes into account the new mean reversion if its has been changed
	/// un peu bourrin mais bon on est pas à ca près
	CreateAndSetModel();

	itsHasBeenComputed = false;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::UpdateCalibration(bool isUpdateStrike)
{
	delete itsPowerReverseSwap;
	itsPowerReverseSwap = NULL;
	delete itsATMFwdFxVol;
	itsATMFwdFxVol = NULL;
	
	/// re-update the funding spread if necessairy
	if(IsBasis())
	{
		ComputeDomesticBasis();
		ARM_GramFctorArg cst(ARM_VectorPtr((ARM_GP_Vector*)itsvFundSpread.Clone()));
		GetGenSecurity()->GetCstManager()->insert(PRDProfileNamesTable[MarginProfile],cst );
	}
	
	bool isInitVol = true;
	bool isFreezeWeights = false;
    /// Compute domestic target prices
    ComputeIROptionPrices(&*GetCalibMethod(),OswDomModelKey,itsDomModelHWFlag,isFreezeWeights,isInitVol);

	/// Compute foreign target prices
    ComputeIROptionPrices(GetCalibMethod()->GetNextMethod(),OswForModelKey,itsForModelHWFlag,isFreezeWeights,isInitVol);

	if(itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib)
    {
        /// Create a standard ATM Fx option portfolio for first calibration step
        itsCalibType    = ARM_PRCSCalibTypes::ATMCalib;
        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVol);
        itsCalibType    = ARM_PRCSCalibTypes::ATMDoubleCalib;

		/// Compute Fx option target prices (+ ATM strike +...)
        ComputeFxOptionPrices(itsATMDoubleCalibMethod,isFreezeWeights,isInitVol);
    }
    else
    {
        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVol);
    }
	if(itsFxLocalModelFlag)
    {
        /// Compute prices of floored,capped & terminal redemption Fx option portfolio
        /// Only used for precision because calibration is done directly on volatilities
        /// and not on prices
        ComputeLocalFxOptionPrices();
    }

	itsHasBeenComputed = false;
}


////////////////////////////////////////////////////
///	Class   : ARM_PRDCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PRDCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;

    /// Additional ATM calib method
    if(itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib && itsATMDoubleCalibMethod)
    {
        prdcData << indent << "\nAdditional Calib Method :\n";
        prdcData << indent << itsATMDoubleCalibMethod->toString(indent,nextIndent);
    }

    /// GenCalculator general datas dump
	prdcData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";

    return prdcData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_PRDCCalculator
///	Routines: HybridBasketDump
///	Returns :
///	Action  : Dump hybrid IR/FX swaption portfolio as
///			  basket on CATU & domestic Zc
////////////////////////////////////////////////////
string ARM_PRDCalculator::HybridBasketDump()
{
    CC_Ostringstream dump;
	ARM_StdPortfolioPtr hybridBaskets;
	ARM_VanillaIrFxSwaption* irFxSwaption;
	ARM_OptionPortfolio* hybridBasketOption;

	double  asOfDate = GetPricingModel()->GetAsOfDate().GetJulian();

	if( itsHybridBasketCalibMethod &&
		(hybridBaskets = itsHybridBasketCalibMethod->GetPortfolio()) != ARM_StdPortfolioPtr(NULL) )
	{
		size_t i,nbBaskets=hybridBaskets->GetSize();
		double zcCoef,strike;
		for(i=0;i<nbBaskets;++i)
		{
			hybridBasketOption = static_cast<ARM_OptionPortfolio*>(hybridBaskets->GetAsset(i));
			irFxSwaption = static_cast<ARM_VanillaIrFxSwaption*>( ARM_ConverterFromKernel::ConvertSecuritytoArgObject(hybridBasketOption,asOfDate) );
			dump << "Basket #" << i << endl;
			dump << "Expiry\tC/P\tStrike\n";
			dump << irFxSwaption->GetExpiry() << "\t";
			dump << irFxSwaption->GetCallPut() << "\t" << irFxSwaption->GetStrike() << endl;
			dump << "FxReset\tFxSettl\tFxPay\tFxNotio\n";
			for(i=0;i<irFxSwaption->GetFxResetTimes().size();++i)
			{
				dump << irFxSwaption->GetFxResetTimes()[i] << "\t";
				dump << irFxSwaption->GetFxSettlementTimes()[i] << "\t";
				dump << irFxSwaption->GetFxPayTimes()[i] << "\t";
				dump << irFxSwaption->GetFxNominals()[i] << endl;
			}
			dump << "IrPay\tIrZcCoef\n";
			strike = irFxSwaption->GetStrike();
			ARM_GP_Vector floatNotio(*(irFxSwaption->GetIrSwaption().GetFloatNotional()));
			ARM_GP_Vector fixNotio(*(irFxSwaption->GetIrSwaption().GetFixNotional()));
			ARM_GP_Vector floatEndTimes(*(irFxSwaption->GetIrSwaption().GetFloatEndTimes()));
			ARM_GP_Vector fixPeriods(*(irFxSwaption->GetIrSwaption().GetFixPayPeriods()));
			zcCoef = - floatNotio[0];
			dump << (*(irFxSwaption->GetIrSwaption().GetFloatStartTimes()))[0] << "\t" << zcCoef << endl;
			for(i=0;i+1<floatNotio.size();++i)
			{
				zcCoef = -floatNotio[i+1] + floatNotio[i] + strike*fixPeriods[i]*fixNotio[i];
				dump << floatEndTimes[i] << "\t" << zcCoef<< endl;
			}
			i = floatNotio.size()-1;
			zcCoef = floatNotio[i] + strike*fixPeriods[i]*fixNotio[i];
			dump << floatEndTimes[i] << "\t" << zcCoef<< endl << endl;
		}
			
	}

	return dump.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_PRDCalculator
///	Routines: ComputeATMFxOptionVols
///	Returns :
///	Action  : Compute the term fx volatility for each
///           Fx option of the PRDC after the 1st notice
////////////////////////////////////////////////////
ARM_VolLInterpol* ARM_PRDCalculator::ComputeATMFxOptionVols()
{
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    

    ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());
    ARM_Date asOfDate(hybridModel->GetAsOfDate());
    double asOfJulian=asOfDate.GetJulian();

    bool isATMCalib = itsCalibType == ARM_PRCSCalibTypes::ATMCalib || itsCalibType == ARM_PRCSCalibTypes::ATMDoubleCalib;


    if(!isATMCalib)
    {
        /// Clone to avoid side effects !
        hybridModel = static_cast< ARM_2IRFXModel* > (GetPricingModel()->Clone());

        /// Save Fx option strikes and set them ATM
        ARM_StdPortfolioPtr fxPortfolio = GetFxPortfolio();
        size_t i,nbFx = fxPortfolio->GetSize();
        ARM_GP_Vector prevStrike(nbFx);
        ARM_Option* fxOption;
        for(i=0;i<nbFx;++i)
        {
            fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(i));
            prevStrike[i] = fxOption->GetStrike();
            fxOption->SetStrike(K_MARKET_RATE);
        }

        /// Compute new target ATM prices
        bool isFreezeWeights=true;
        bool isInitVolParam=true;

        ARM_PRDCCalibType prevCalibType = itsCalibType;
        itsCalibType = ARM_PRCSCalibTypes::ATMCalib;

        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVolParam);

        /// Calibrate stochatic models ATM
        GetCalibMethod()->Calibrate(hybridModel);

        CheckCalibErrors();

        /// Restore original calibration strikes & prices
        for(i=0;i<nbFx;++i)
            (static_cast< ARM_Option* >(fxPortfolio->GetAsset(i)))->SetStrike( prevStrike[i] );
        itsCalibType = prevCalibType;
        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVolParam);
    }

    /// Get the Fx model
    const ARM_ModelNameMap& modelMap = * hybridModel->GetModelMap();
    ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[ARM_2IRFXModel::FxModel]->Model() );


    /// Compute Fx options then get saved forward Fx vols
    double evalTime=0.0;
    bool isDualOption = itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption;
	double lastfxResetDate = (*(itsForexDateStrip->GetResetDates()))[itsCpnSize-1];
    bool isDiffLastExpiry = isDualOption && (itsResetRedemptionDate.GetJulian() != lastfxResetDate);

    double fxResetTime,settlementTime,payTime;
    int callPut=K_CALL;
    ARM_GP_Vector strikePerStates(1);
    ARM_GP_VectorPtr ATMStrike;
	size_t nbResetsBefore = itsFxResetTimeBefore.size();
    size_t nbResets = nbResetsBefore -itsNbNoCall + itsExerSize + (isDualOption ? 1 : 0);
    ARM_GP_Vector fxResetYfs(nbResets);
    ARM_GP_Matrix fwdFxVols(nbResets,4,0.0);
    ARM_GP_VectorPtr price;

    /// Compute market ATM FX vols
    ARM_VolLInterpol* marketATMFxVols = NULL;

    if (itsPowerReverseSwap->IsFxVolModelExisting())
    {
       marketATMFxVols = itsPowerReverseSwap->ComputeATMFxVolFromFxVolModel();
    }
    else
    {
		//ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
		//ARM_FXVolCurve* marketFxVol = static_cast<ARM_FXVolCurve*>(fxBSModel->GetVolatility());

		/// To support new FX GP models !
		ARM_DFBSModel* fxBSModel = dynamic_cast< ARM_DFBSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
		ARM_FXVolCurve* marketFxVol = static_cast<ARM_FXVolCurve*>(fxBSModel->GetFxVol());
		if(!marketFxVol)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : market FX vol is missing");
       marketATMFxVols = marketFxVol->ComputeATMFxVol();
    }

    /// The first notice date is always at the first line of the itsExerciseDateStrip
    double firstNoticeTime = (*itsExerciseDateStrip->GetResetDates())[0]-asOfJulian;

    bool prevIsImpliedVolCalc=fxModel.GetIsImpliedVolCalc();
    fxModel.SetIsImpliedVolCalc(true); // to activate B&S implied vol calculation

    /// Add resets and vols saved before the 1st notice
	double fxResetDate;
    for(size_t resetIdx(0);resetIdx < nbResetsBefore;++resetIdx)
    {
        /// Use market ATM forward FX vols
		fxResetTime = itsFxResetTimeBefore[resetIdx];
        fxResetYfs[resetIdx] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol will be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol will be used

	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        fwdFxVols(resetIdx,0) = (*ATMStrike)[0];
        fwdFxVols(resetIdx,1) = 0.01*marketATMFxVols->ComputeVolatility(fxResetYfs[resetIdx],(*ATMStrike)[0]);
    }

    /// Add resets and vols after the 1st notice
    for(resetIdx = itsNbNoCall; resetIdx < itsExerSize;++resetIdx)
    {
		int index = resetIdx-itsNbNoCall;
		fxResetDate = (*(itsForexDateStrip->GetResetDates()))[itsvCpnIndex[resetIdx]];
        fxResetTime = fxResetDate - asOfJulian;
        fxResetYfs[nbResetsBefore+index] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol will be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol will be used

        /// Use boostrapped forward FX vols
	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        price = fxModel.CallVectorial(fxModel.GetModelName(),evalTime,fxResetTime,settlementTime,*ATMStrike,callPut,payTime,dumStates);
        fwdFxVols(nbResetsBefore+index,0) = (*ATMStrike)[0];
        fwdFxVols(nbResetsBefore+index,1) = (fxModel.GetFwdFxVol())[0];
    }

    /// Get Dual redemption vol if any
    callPut=K_PUT;
    if(isDualOption)
    {
        fxResetTime = itsResetRedemptionDate.GetJulian()-asOfJulian;
        fxResetYfs[nbResetsBefore+resetIdx] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol is be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol is be used

        /// Use boostrapped forward FX vols
	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        strikePerStates[0]= (*ATMStrike)[0];

        price = fxModel.CallVectorial(fxModel.GetModelName(),evalTime,fxResetTime,settlementTime,strikePerStates,callPut,payTime,dumStates);
        fwdFxVols(nbResetsBefore+resetIdx,0) = (*ATMStrike)[0];
        fwdFxVols(+nbResetsBefore+resetIdx,1) = (fxModel.GetFwdFxVol())[0];
    }

    /// Get FX option strikes if necessary
    ARM_GP_Vector strikes(1,0.0);
    size_t nbStrikes=strikes.size();

    nbResets -= (isDualOption && !isDiffLastExpiry ? 1 : 0);
    fxResetYfs.resize(nbResets);
    ARM_Matrix* vols = new ARM_Matrix(nbResets,nbStrikes,0.0);

    /// Take care Kernel is in 100 base
    for(resetIdx=0; resetIdx<nbResetsBefore+itsExerSize-itsNbNoCall; ++resetIdx)
        vols->Elt(resetIdx,0) = 100.0*fwdFxVols(resetIdx,1);

    if(isDiffLastExpiry)
        /// Last Fx option & dual option reset at different times
        vols->Elt(resetIdx,0) = 100.0*fwdFxVols(resetIdx,1);

    fxModel.SetIsImpliedVolCalc(prevIsImpliedVolCalc);

    /// Convert ARM_GP_Vector to ARM_Vector
    ARM_Vector* theResets = To_pARM_Vector(&fxResetYfs);
    ARM_Vector* theStrikes = To_pARM_Vector(&strikes);

    /// Free memory if necessary
    if(!isATMCalib){
        delete hybridModel;
		/// Ugly Ugly but it's necessairy to avoid the PRDC Calculator crashing
		 /// Calibrate stochatic models ATM
        GetCalibMethod()->Calibrate(&*GetPricingModel());
	}

	if (marketATMFxVols)
		delete marketATMFxVols;

    return new ARM_VolLInterpol(asOfDate,theResets,theStrikes,vols);
}

////////////////////////////////////////////////////
///	Class   : ARM_PRDCalculator
///	Routines: CheckCalibErrors
///	Returns :
///	Action  : Check calibration errors and throw an
///           exception if errors are too large
////////////////////////////////////////////////////
void ARM_PRDCalculator::CheckCalibErrors()
{
    /// Recall that Dom -> For -> Fx with "next" links

    size_t i,nbErrors,calibIdx=0;
    double err;
    ARM_CalibMethod* calib = &*(GetCalibMethod());
    ARM_ModelFitterPtr modelFitter;


    while( calib != NULL &&
           (modelFitter=calib->GetModelFitter()) != ARM_ModelFitterPtr(NULL) )
    {
        modelFitter->SetUpError();

        nbErrors = modelFitter->GetError()->rows();
        for(i=0; i<nbErrors; ++i)
        {
            err=(*(modelFitter->GetError()))(i,0);
            if(fabs(err) > MAX_CALIB_ERR[calibIdx])
            {
                ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + CALIB_ERR_MSGE[calibIdx]);
            }
        }

        calib=calib->GetNextMethod();
        ++calibIdx;
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a diagonal swaption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRDCalculator::GetOSWPortfolio(mdmKeysAlias oswModelKey) const
{
    /// Recall that Dom -> For -> Fx with "next" links
    switch(oswModelKey)
    {
    case OswDomModelKey:
        if(GetCalibMethod() != ARM_CalibMethodPtr(NULL))
            return GetCalibMethod()->GetPortfolio();

    case OswForModelKey:
        if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
            GetCalibMethod()->GetNextMethod() != NULL )
            return GetCalibMethod()->GetNextMethod()->GetPortfolio();
    }

    return ARM_StdPortfolioPtr(NULL);
}

void ARM_PRDCalculator::SetOSWPortfolio(const ARM_StdPortfolio& pf, mdmKeysAlias oswModelKey)
{
	/// Recall that Dom -> For -> Fx with "next" links
    switch(oswModelKey)
    {
    case OswDomModelKey:
        if(GetCalibMethod() != ARM_CalibMethodPtr(NULL))
			GetCalibMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(pf).Clone())) );
		break;
    case OswForModelKey:
        if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) && GetCalibMethod()->GetNextMethod() != NULL )
			GetCalibMethod()->GetNextMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(pf).Clone())) );
		break;
	default:
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : portfolio type is not valid");

    }

}
void ARM_PRDCalculator::SetFxPortfolio(const ARM_StdPortfolio& pf)
{
	/// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL && 
        GetCalibMethod()->GetNextMethod()->GetNextMethod() != NULL ) 
		GetCalibMethod()->GetNextMethod()->GetNextMethod()->SetPortfolio( ARM_StdPortfolioPtr((ARM_StdPortfolio*)(const_cast< ARM_StdPortfolio& >(pf).Clone())) );
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: GetFxPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a fx option calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRDCalculator::GetFxPortfolio() const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL && 
        GetCalibMethod()->GetNextMethod()->GetNextMethod() != NULL ) 
        return GetCalibMethod()->GetNextMethod()->GetNextMethod()->GetPortfolio();
    else
        return ARM_StdPortfolioPtr(NULL);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: GetExtraFxPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get the fx option calibration portfolio created for ATDoubleCalib
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRDCalculator::GetExtraFxPortfolio() const
{
    if( itsATMDoubleCalibMethod ) 
        return itsATMDoubleCalibMethod->GetPortfolio();
    else if(itsHybridBasketCalibMethod)
        return itsHybridBasketCalibMethod->GetPortfolio();
	else
        return ARM_StdPortfolioPtr(NULL);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: GetFxCalib
///	Returns: ARM_CalibMethod*
///	Action : Get the fx calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRDCalculator::GetFxCalib() const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        return GetCalibMethod()->GetNextMethod()->GetNextMethod();
    else
        return NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: SetFxCalib
///	Returns: void
///	Action : Set the fx calibration method
/////////////////////////////////////////////////////////////////
void ARM_PRDCalculator::SetFxCalib(ARM_CalibMethod* fxCalib) const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        GetCalibMethod()->GetNextMethod()->SetNextMethod(fxCalib);
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't set the Fx calib method");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: GetExtraFxCalib
///	Returns: ARM_CalibMethod*
///	Action : Get the extra calibration method for ATMDoubleCalib
///			 or HybridBasketCalib types
/////////////////////////////////////////////////////////////////
const ARM_CalibMethod* ARM_PRDCalculator::GetExtraFxCalib() const
{
    if( itsATMDoubleCalibMethod ) 
        return itsATMDoubleCalibMethod;
    else if(itsHybridBasketCalibMethod)
        return itsHybridBasketCalibMethod;
	else
        return NULL;
}



CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

