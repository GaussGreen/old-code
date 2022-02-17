/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file prdccalculator.cpp
 *
 *  \brief file for the Power Reverse Dual Currencies Calculator
 *	\author  J-M Prié
 *	\version 1.0
 *	\date February 2005
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/prcscalculator.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/autocleaner.h"
#include "gpbase/ostringstream.h"
#include "gpbase/datestripcombiner.h"
#include "gpbase/datestrip.h"
#include "gpbase/singleton.h"
#include "gpbase/utilityport.h"  /// for CC_Min

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
#include "gpinfra/modelnamemap.h"
#include "gpinfra/timeinfo.h"
#include "gpinfra/discretisationscheme.h"

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
#include <inst/powrev.h>
#include <mod/xbsfx.h>
#include <crv/volflat.h>
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/fixleg.h>
#include <inst/swapleg.h>
#include <inst/swap.h>
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

/// Structured coupon type
const double FOREX_CPN = -1.0;
const double FIXED_CPN = +1.0;

/// Strike maximum value to allow call pricing
const double MAX_FX_STRIKE = 5000.0;


/// RedemptionData names
const unsigned int REDEMPTION_TYPE          = 0;
const unsigned int REDEMPTION_STRIKE        = 1;
const unsigned int REDEMPTION_PAY_DATE      = 2;
const unsigned int REDEMPTION_RESET_DATE    = 3;

/// Redemption clause
const unsigned STANDARD_REDEMPTION          = K_OPT_NO;
const unsigned MANDATORY_REDEMPTION         = K_OPT_FWD;
const unsigned DUAL_OPTION_REDEMPTION       = K_OPT;


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
const string LOCAL_FXMODEL_KEY_NAME         = "LOCAL_FXMOD_";
const string MARKET_IRMODEL_KEY_NAME		= "MARKET_IRMOD_";

const string UNKNOWN_KEY_NAME               = "UNKNOWN";

/// Flag to replace standard domestic diagonal swaption calibration by a
/// a calibration on domestic variable notional swaptions
/// The case "true" doesn't work at the moment because
/// Q1F model doesn't support variable notional swaption pricing
const bool IS_DOM_BASKET_CALIB	= false;

const string ARM_PRCSCalculator::PRDCColNamesTable [] =
{
    "EventDate",
    "FxResetDate",
    "FxLowStrike",
    "FxHighStrike",
    "CpnIT",
    "CpnPayDate",
    "FundingStartDate",
    "FundingEndDate",
    "FundingLastEndDate",
    "FundingPayDate",
    "FundingSpread",
    "FundingIT",
    "FundingFlow",
    "Funding",
    "CpnLeverage",
    "FixedCpn",
    "FxLowCall",
    "FxHighCall",
    "FxLowCallStrip",
    "FxHighCallStrip",
    "RedemptionResetDate",
    "RedemptionPayDate",
    "RedemptionStrike",
    "Redemption",
    "CpnFlow",
    "Cpn",
    "PRDCFlow",
    "PRDCSwap",
    "PRDCFirstSwap",
    "PRDCFirstEuropean",
    "PRDCBermuda",
    "PRDCOption"
};

const string ARM_PRCSCalculator::ExtPRDCColNamesTable [] =
{
    "FundingSum",
    "CouponSum",
    "FundingSum1",
    "CouponSum1",
    "FundingSum2",
    "CouponSum2",
    "F"
};


const string ARM_PRCSCalculator::PRDCProfileNamesTable [] =
{
    "FxDateStripProfile",
    "FxLowStrikeProfile",
    "FxHighStrikeProfile",
    "FxLeverageProfile",
    "FundingDateStripProfile",
    "FundingSpreadProfile"
};


/// Maximum calibration errors in % (IR) and real terms (FX)
const double MAX_CALIB_ERR[] = {    0.01,
                                    0.01,
                                    0.0001,
									0.0001};

const string CALIB_ERR_MSGE[] = {   ": domestic IR calibration failed",
                                    ": foreign IR calibration failed",
                                    ": forex calibration failed",
                                    ": hybrid IR/FX calibration failed"};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_PRCSCalculator::ARM_PRCSCalculator( const ARM_PRCSCalculator& rhs )
:	ARM_GenCalculator( rhs ),
    itsInputPowRev( CreateClone(rhs.itsInputPowRev) ),
    itsInputModel( CreateClone(rhs.itsInputModel) ),
    itsFxResetTimeBefore( rhs.itsFxResetTimeBefore),
    itsBoosterDatas( rhs.itsBoosterDatas),
    itsRedemptionDatas( rhs.itsRedemptionDatas),
    itsExerciseFlag( rhs.itsExerciseFlag),
    itsColumnsToPrice( rhs.itsColumnsToPrice),
    itsSchedulerDatas( rhs.itsSchedulerDatas),
    itsTruncatorDatas( rhs.itsTruncatorDatas),
    itsDomModelHWFlag( rhs.itsDomModelHWFlag),
    itsForModelHWFlag( rhs.itsForModelHWFlag),
    itsFirstFutureEventIdx( rhs.itsFirstFutureEventIdx),
    itsMarkovianDriftSamplerFlag( rhs.itsMarkovianDriftSamplerFlag),
	itsFlooredFxOptionPF( CreateClone(rhs.itsFlooredFxOptionPF)),
	itsCappedFxOptionPF( CreateClone(rhs.itsCappedFxOptionPF)),
	itsRedemptionFxOptionPF( CreateClone(rhs.itsRedemptionFxOptionPF)),
    itsFxLocalModelFlag( rhs.itsFxLocalModelFlag),
    itsCalibType    ( rhs.itsCalibType),
    itsCalibDatas   ( rhs.itsCalibDatas),
	itsATMFwdFxVol (CreateClone(rhs.itsATMFwdFxVol)),
    itsBasisIRCalibFlag ( rhs.itsBasisIRCalibFlag),
	itsATMDoubleCalibMethod ( CreateClone(rhs.itsATMDoubleCalibMethod)),
	itsMktHybridModel(CreateClone(rhs.itsMktHybridModel)),
	itsHybridBasketCalibMethod(CreateClone(rhs.itsHybridBasketCalibMethod)),
	itsHasBeenPriced(false)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_PRCSCalculator::~ARM_PRCSCalculator()
{
    delete itsInputPowRev;
	delete itsInputModel;

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
///	Class  : ARM_PRCSCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_PRCSCalculator::ARM_PRCSCalculator(
                ARM_PowerReverse* powRev,
                ARM_Model* model,
                const ARM_ObjectVector& otherMktDatas,
                const ARM_GP_Vector& schedulerDatas,
                const ARM_GP_Vector& truncatorDatas,
                const ARM_StringVector& columnsToPrice,
                bool markovianDriftSamplerFlag,
                bool fxLocalModelFlag,
                PRDCCalibType calibType,
                const ARM_GP_Vector& calibDatas,
                bool basisIRCalibFlag)
:	
	ARM_GenCalculator(model ? model->GetZeroCurve()->GetAsOfDate() : ARM_Date()),
	itsInputPowRev( CreateClone(powRev) ),
	itsInputModel( CreateClone(dynamic_cast< ARM_DFBSModel* >(model)) ),
    itsFxResetTimeBefore(0),
	itsBoosterDatas(0,0),
    itsSchedulerDatas(schedulerDatas),
    itsTruncatorDatas(truncatorDatas),
    itsExerciseFlag(0),
    itsColumnsToPrice(columnsToPrice),
    itsDomModelHWFlag(false),
    itsForModelHWFlag(false),
    itsFirstFutureEventIdx(0),
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
	itsHasBeenPriced(false)
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
    ARM_Currency* domCcy = itsInputModel->GetDBSModel()->GetZeroCurve()->GetCurrencyUnit();
    SetCurrencyUnit(domCcy);
    string domCcyName( domCcy->GetCcyName() );

    /// Set foreign currency (of forex)
    string forCcyName( itsInputModel->GetFBSModel()->GetZeroCurve()->GetCurrencyUnit()->GetCcyName() );

    string fxName(forCcyName + "/" + domCcyName);
    string domForFxName("(" + domCcyName + "," + forCcyName + "," + fxName + ")");

    ARM_StringVector mdmKeys(NbKeys);
    mdmKeys[YcDomKey]               = YC_KEY_NAME               + domCcyName;
    mdmKeys[YcForKey]               = YC_KEY_NAME               + forCcyName;
    mdmKeys[ForexKey]               = FOREX_KEY_NAME            + fxName;
    mdmKeys[YcBasisDomKey]          = YC_BASIS_KEY_NAME         + domCcyName;
    mdmKeys[YcBasisForKey]          = YC_BASIS_KEY_NAME         + forCcyName;
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

    /// Compute redemption datas & booster table then reaffect
    /// notice dates to relevant flows
    ComputeBoosterDatas();

    ARM_CstManagerPtr fxDataManager(NULL);
    if(itsFxLocalModelFlag)
    {
        /// Build an internal cst manager to save
        /// Fx options date strip + strikes, nominals & leverage profiles
        fxDataManager=ComputeCstManager();
    }

    /// Create the Generic Security paid in domestic=coupon currency
    CreateAndSetDealDescriptionAndTimeIt(mdmKeys[YcBasisDomKey],itsColumnsToPrice,fxDataManager);

    /// To boost time, disable intermediatePayoffs & snapshots computation
    GetGenSecurity()->SetOtherPayoffsFlag(false);


    /// Collect market datas from the input model and intialize the MdM
    FillMarketDataManager(otherMktDatas);

    /// Check market datas
    CheckMktDataAndTimeIt();

    /// Create a 2IR+FX model
    CreateAndSetModelAndTimeIt();

    /// Create calibration sets
    CreateAndSetCalibrationAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: GetOSWPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a diagonal swaption calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRCSCalculator::GetOSWPortfolio(mdmKeysAlias oswModelKey) const
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


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: GetFxPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get a fx option calibration portfolio
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRCSCalculator::GetFxPortfolio() const
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
///	Class  : ARM_PRCSCalculator
///	Routine: GetExtraFxPortfolio
///	Returns: ARM_StdPortfolioPtr
///	Action : Get the extra calibration portfolio created for ATMDoubleCalib
///			 or HybridBasketCalib
/////////////////////////////////////////////////////////////////
const ARM_StdPortfolioPtr ARM_PRCSCalculator::GetExtraFxPortfolio() const
{
    if( itsATMDoubleCalibMethod ) 
        return itsATMDoubleCalibMethod->GetPortfolio();
    else if(itsHybridBasketCalibMethod)
        return itsHybridBasketCalibMethod->GetPortfolio();
	else
        return ARM_StdPortfolioPtr(NULL);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: GetFxCalib
///	Returns: ARM_CalibMethod*
///	Action : Get the fx calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRCSCalculator::GetFxCalib() const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        return GetCalibMethod()->GetNextMethod()->GetNextMethod();
    else
        return NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: SetFxCalib
///	Returns: void
///	Action : Set the fx calibration method
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::SetFxCalib(ARM_CalibMethod* fxCalib) const
{
    /// Recall that Dom -> For -> Fx with "next" links
    if( GetCalibMethod() != ARM_CalibMethodPtr(NULL) &&
        GetCalibMethod()->GetNextMethod() != NULL ) 
        GetCalibMethod()->GetNextMethod()->SetNextMethod(fxCalib);
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't set the Fx calib method");
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: GetExtraFxCalib
///	Returns: ARM_CalibMethod*
///	Action : Get the extra calibration method for ATMDoubleCalib
///			 or HybridBasketCalib types
/////////////////////////////////////////////////////////////////
const ARM_CalibMethod* ARM_PRCSCalculator::GetExtraFxCalib() const
{
    if( itsATMDoubleCalibMethod ) 
        return itsATMDoubleCalibMethod;
    else if(itsHybridBasketCalibMethod)
        return itsHybridBasketCalibMethod;
	else
        return NULL;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: FillMarketDataManager
///	Returns: ARM_ObjectVector
///	Action : Collect market datas from the input model
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::FillMarketDataManager(const ARM_ObjectVector& mktDatas)
{
    /// Collect market datas from the input model
    /// -----------------------------------------

    /// Domestic vanilla swaptions are computed using a standard non basis BS model
    ARM_BSModel* oswDomModel = itsInputModel->GetDBSModel();
    GetMktDataManager()->RegisterData(GetKeys()[OswDomModelKey],static_cast< ARM_Object* >(oswDomModel));

    ARM_ZeroCurve* ycDomCurve = itsInputModel->GetDBSModel()->GetZeroCurve();
    GetMktDataManager()->RegisterData(GetKeys()[YcDomKey],static_cast< ARM_Object* >(ycDomCurve));

    /// Foreign vanilla swaptions are computed using a standard non basis BS model
    ARM_BSModel* oswForModel = itsInputModel->GetFBSModel();
    GetMktDataManager()->RegisterData(GetKeys()[OswForModelKey],static_cast< ARM_Object* >(oswForModel));

    ARM_ZeroCurve* ycForCurve = itsInputModel->GetFBSModel()->GetZeroCurve();
    GetMktDataManager()->RegisterData(GetKeys()[YcForKey],static_cast< ARM_Object* >(ycForCurve));

    ARM_BSModel* fxModel = itsInputModel->CreateFxBSModel();
    GetMktDataManager()->RegisterData(GetKeys()[FxModelKey],static_cast< ARM_Object* >(fxModel));

    double fxSpot = itsInputModel->GetFxSpot();
    ARM_Forex forex(ycDomCurve->GetCurrencyUnit(), ycForCurve->GetCurrencyUnit(), fxSpot);
    GetMktDataManager()->RegisterData(GetKeys()[ForexKey],static_cast< ARM_Object* >(&forex));

    ARM_ZeroCurve* ycBasisDomCurve = itsInputModel->GetDBsCrv();
    GetMktDataManager()->RegisterData(GetKeys()[YcBasisDomKey],static_cast< ARM_Object* >(ycBasisDomCurve));

    ARM_ZeroCurve* ycBasisForCurve = itsInputModel->GetFBsCrv();
    GetMktDataManager()->RegisterData(GetKeys()[YcBasisForKey],static_cast< ARM_Object* >(ycBasisForCurve));

    ARM_VolFlat* domFxCorrel = dynamic_cast<ARM_VolFlat*>(itsInputModel->GetDomFxCorrel());
    ARM_VolFlat* forFxCorrel = dynamic_cast<ARM_VolFlat*>(itsInputModel->GetForeignFxCorrel());
    if(!domFxCorrel || !forFxCorrel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't restore correlations for input model");

    ARM_GP_Matrix correlMatrix(3,3,1.0);
    correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel) = domFxCorrel->GetVolatility();
    correlMatrix(ARM_2IRFXModel::FxModel,ARM_2IRFXModel::DomModel) = correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::FxModel);
    correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel) = forFxCorrel->GetVolatility();
    correlMatrix(ARM_2IRFXModel::FxModel,ARM_2IRFXModel::ForModel) = correlMatrix(ARM_2IRFXModel::ForModel,ARM_2IRFXModel::FxModel);

    correlMatrix(ARM_2IRFXModel::DomModel,ARM_2IRFXModel::ForModel) = itsInputModel->GetRatesCorr()/100.0;
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
    if(qFx)
    {
        /// If no Q parameter is input, model will be pure log-normal (Q=1)
        GetMktDataManager()->RegisterData(GetKeys()[QFxKey],static_cast< ARM_Object* >(qFx));
    }

    ARM_CurveModelParam* qDom = dynamic_cast<ARM_CurveModelParam*>(mktDatas[QDomKey]);
    if(qDom)
    {
        /// If no Q parameter is input, model will be degenerated in H&W
        GetMktDataManager()->RegisterData(GetKeys()[QDomKey],static_cast< ARM_Object* >(qDom));
    }

    ARM_CurveModelParam* qFor = dynamic_cast<ARM_CurveModelParam*>(mktDatas[QForKey]);
    if(qFor)
    {
        /// If no Q parameter is input, model will be degenerated in H&W
        GetMktDataManager()->RegisterData(GetKeys()[QForKey],static_cast< ARM_Object* >(qFor));
    }

    if(itsFxLocalModelFlag)
    {
        ARM_DFBSModel* localDFBSModel = dynamic_cast<ARM_DFBSModel*>(mktDatas[LocalFxModelKey]);
        if(localDFBSModel)
        {
            /// Use the specific market Fx model given in input
            ARM_BSModel* localFxModel = localDFBSModel->CreateFxBSModel();
            GetMktDataManager()->RegisterData(GetKeys()[LocalFxModelKey],static_cast< ARM_Object* >(localFxModel));

            delete localFxModel;  // because cloned by RegisterData()
        }
        else
        {
            /// Use the default one computed by the dual currencies BS model
            GetMktDataManager()->RegisterData(GetKeys()[LocalFxModelKey],static_cast< ARM_Object* >(fxModel));
        }
    }

	ARM_MarketIRModel* marketIrModel = dynamic_cast<ARM_MarketIRModel*>(mktDatas[MarketIrModelKey]);
	if(itsCalibType == HybridBasketCalib && marketIrModel)
	{
		/// If market IR model is missing, VNS evaluation will be degenerated
		GetMktDataManager()->RegisterData(GetKeys()[MarketIrModelKey],static_cast< ARM_Object* >(marketIrModel));
	}

    delete fxModel; // because cloned by RegisterData()


}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeBoosterDatas
///	Returns: void
///	Action : Build booster table & redemption vector then
///          associate notice dates to the right flow
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputeBoosterDatas()
{
    /// Compute the PRDC analytical price, to generate fictive funding spreads,
    /// low & high strikes, leverage etc... saved in THE booster matrix !
    itsInputPowRev->SetModel(itsInputModel);
    itsInputPowRev->SetFXSmileAdjFlag(true); // to allow market vol smile adj to be saved for further using
    double prdcSwapPrice = itsInputPowRev->ComputePrice();

    double asOfDate = itsInputModel->GetZeroCurve()->GetAsOfDate().GetJulian();


    /// Get the redemption datas and save them only a special redemption clause 
    ARM_Vector* redemptionDatas = itsInputPowRev->GetDualOptVect();
    size_t i,j,nbRedDatas = redemptionDatas ? redemptionDatas->GetSize() : 0;
    if(nbRedDatas>0 && (*redemptionDatas)[0] != STANDARD_REDEMPTION)
    {
        itsRedemptionDatas.resize(nbRedDatas);
        for(i=0;i<nbRedDatas;++i)
            itsRedemptionDatas[i] = (*redemptionDatas)[i];
    }
    /// Call to GetDualOptVect() returns a new ARM_Vector then don't forget to delete it !
    delete redemptionDatas;

    /// Get the BoosterDatas table
    ARM_Matrix* dataMatrix1 = itsInputPowRev->GetPRCSDataMatrix();
    if(!dataMatrix1)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't get the PRDC data matrix");
    ARM_Matrix* dataMatrix2 = itsInputPowRev->CalcReformattedPRCSDataMatrixFromDataMatrix(dataMatrix1);
    if(!dataMatrix2)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't get the PRDC data matrix");

    itsBoosterDatas.resize(dataMatrix2->GetNumLines(),dataMatrix2->GetNumCols());
    size_t nbFlows = itsBoosterDatas.rows();
    for(i=0;i<nbFlows;++i)
        for(j=0;j<itsBoosterDatas.cols();++j)
            itsBoosterDatas(i,j)= dataMatrix2->Elt(i,j);

    delete dataMatrix1;
    delete dataMatrix2;


    /// Initialise exercise flag vector to say if a notice is actual or not
    itsExerciseFlag.resize(nbFlows,false);

    /// Affect notice dates to relevant flows and
    /// activate exercise flag for an actual exercise
    double noticeDate,defaultNoticeDate;
    ARM_GP_Vector noticeDates(nbFlows,0.0);

    int lastNoticeIdx=-1;
    int flowIdx=0;
    int firstNoticeIdx=-1;
    for(size_t noticeIdx=0;noticeIdx<nbFlows && (noticeDate=itsBoosterDatas(noticeIdx,NOTICE_DATE)) > 0;++noticeIdx)
    {
        while(flowIdx < nbFlows && (defaultNoticeDate=CC_Min(itsBoosterDatas(flowIdx,FX_RESET_DATE),itsBoosterDatas(flowIdx,FUNDING_START_DATE))) < noticeDate)
        {
            /// No call date => notice date is set to minimum value between
            /// fx reset date & funding start date
            noticeDates[flowIdx] = defaultNoticeDate;
            ++flowIdx;
        }
        if(flowIdx < nbFlows)
        {
            noticeDates[flowIdx] = noticeDate;
            itsExerciseFlag[flowIdx] = true;

            /// Keep record of the 1st notice position after the as of date
            if(firstNoticeIdx==-1 && floor(noticeDate) > asOfDate + K_NEW_DOUBLE_TOL)
                firstNoticeIdx=flowIdx;

            lastNoticeIdx = flowIdx;
            ++flowIdx;
        }
    }


    /// Update NoticeDate column in the booster table
    for(flowIdx=0;flowIdx <= lastNoticeIdx;++flowIdx)
        itsBoosterDatas(flowIdx,NOTICE_DATE) = noticeDates[flowIdx];

    /// After last notice date, just set notice to minimum value between
    /// fx reset date & funding start date
    for(;flowIdx < nbFlows;++flowIdx)
        itsBoosterDatas(flowIdx,NOTICE_DATE) = CC_Min(itsBoosterDatas(flowIdx,FX_RESET_DATE),itsBoosterDatas(flowIdx,FUNDING_START_DATE));

/**** Dump check ****
FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
fprintf(f,"\n\nAll structure\n");
for(flowIdx=0;flowIdx < itsBoosterDatas.rows();++flowIdx)
{
    fprintf(f,"Exer ?=%1d\t",(int)itsExerciseFlag[flowIdx]);
    fprintf(f,"Ev=%11.3lf\t",itsBoosterDatas(flowIdx,NOTICE_DATE)-asOfDate);
    fprintf(f,"FxR=%7.0lf\t",itsBoosterDatas(flowIdx,FX_RESET_DATE)-asOfDate);
    fprintf(f,"FxP=%7.0lf\t",itsBoosterDatas(flowIdx,FX_PAY_DATE)-asOfDate);
    fprintf(f,"FuS=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_START_DATE)-asOfDate);
    fprintf(f,"FuE=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_END_DATE)-asOfDate);
    fprintf(f,"FuP=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_PAY_DATE)-asOfDate);
    fprintf(f,"FuSp=%10.7lf\t",itsBoosterDatas(flowIdx,FUNDING_SPREAD));
    fprintf(f,"Lev=%10.7lf\t",itsBoosterDatas(flowIdx,FX_NOTIO_MULTIPLIER));
    fprintf(f,"LoK=%10.5lf\t",itsBoosterDatas(flowIdx,FX_STRIKE));
    fprintf(f,"HiK=%10.5lf\t",itsBoosterDatas(flowIdx,FX_CAP));
    fprintf(f,"FxI=%10.7lf\t",itsBoosterDatas(flowIdx,FX_ACCRUAL_BASIS));
    fprintf(f,"FuI=%10.7lf\t",itsBoosterDatas(flowIdx,FUNDING_ACCRUAL_BASIS));
    fprintf(f,"Typ=%2d\n",(int)itsBoosterDatas(flowIdx,COUPON_DOMESTIC));
}
fclose(f);
**** Dump check ****/

    /// Save all futur Fx reset times associated to a notice date before the 1st one
    double fxResetTime;
    for(flowIdx=0;flowIdx < firstNoticeIdx;++flowIdx)
    {
        if((fxResetTime=itsBoosterDatas(flowIdx,FX_RESET_DATE)-asOfDate) >= 0)
            itsFxResetTimeBefore.push_back(fxResetTime);
    }


    /// Truncate the booster data to begin at the 1st notice date
    /// Until a generic security can't have different lines with same event date
    /// (here asOfDate will be the best), it will not possible to price the deal
    /// before the 1st notice date without adding no callable fx reset dates !
    /// Here to get the same diffusion schedule as Tree3F, this is not posssible
    /// to add such artificial event dates
    /// Do the same for exercise flags that begin necessary to a "true"
    if(firstNoticeIdx > -1)
    {
        ARM_BoolVector initialExerciseFlag(itsExerciseFlag);
        itsExerciseFlag.resize(initialExerciseFlag.size()-firstNoticeIdx);

        ARM_GP_Matrix initialBoosterDatas(itsBoosterDatas);
        itsBoosterDatas.resize(initialBoosterDatas.rows()-firstNoticeIdx,initialBoosterDatas.cols());

        for(flowIdx=0;flowIdx<itsBoosterDatas.rows();++flowIdx)
        {
            itsExerciseFlag[flowIdx]=initialExerciseFlag[firstNoticeIdx+flowIdx];
            for(j=0;j<itsBoosterDatas.cols();++j)
                itsBoosterDatas(flowIdx,j)=initialBoosterDatas(firstNoticeIdx+flowIdx,j);
        }

    }

/**** Dump check ****
FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
fprintf(f,"\n\nFuture Fx Resets before 1st notice\n");
for(flowIdx=0;flowIdx < itsFxResetTimeBefore.size();++flowIdx)
    fprintf(f,"FxR=%7.0lf\n",itsFxResetTimeBefore[flowIdx]);
fprintf(f,"\n\nStructure from 1st notice\n");
for(flowIdx=0;flowIdx < itsBoosterDatas.rows();++flowIdx)
{
    fprintf(f,"Exer ?=%1d\t",(int)itsExerciseFlag[flowIdx]);
    fprintf(f,"Ev=%11.3lf\t",itsBoosterDatas(flowIdx,NOTICE_DATE)-asOfDate);
    fprintf(f,"FxR=%7.0lf\t",itsBoosterDatas(flowIdx,FX_RESET_DATE)-asOfDate);
    fprintf(f,"FxP=%7.0lf\t",itsBoosterDatas(flowIdx,FX_PAY_DATE)-asOfDate);
    fprintf(f,"FuS=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_START_DATE)-asOfDate);
    fprintf(f,"FuE=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_END_DATE)-asOfDate);
    fprintf(f,"FuP=%7.0lf\t",itsBoosterDatas(flowIdx,FUNDING_PAY_DATE)-asOfDate);
    fprintf(f,"FuSp=%10.7lf\t",itsBoosterDatas(flowIdx,FUNDING_SPREAD));
    fprintf(f,"Lev=%10.7lf\t",itsBoosterDatas(flowIdx,FX_NOTIO_MULTIPLIER));
    fprintf(f,"LoK=%10.5lf\t",itsBoosterDatas(flowIdx,FX_STRIKE));
    fprintf(f,"HiK=%10.5lf\t",itsBoosterDatas(flowIdx,FX_CAP));
    fprintf(f,"FxI=%10.7lf\t",itsBoosterDatas(flowIdx,FX_ACCRUAL_BASIS));
    fprintf(f,"FuI=%10.7lf\t",itsBoosterDatas(flowIdx,FUNDING_ACCRUAL_BASIS));
    fprintf(f,"Typ=%2d\n",(int)itsBoosterDatas(flowIdx,COUPON_DOMESTIC));
}
fclose(f);
**** Dump check ****/

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeCstManager
///	Returns: void
///	Action : Build a constant manager with Fx option date strip &
///          profiles
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_PRCSCalculator::ComputeCstManager()
{
    size_t eventIdx=0,nbEvents = itsBoosterDatas.rows();

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    ARM_GP_Vector fxExpiryDates;
    ARM_GP_Vector fxPaymentDates;
    ARM_GP_VectorPtr fxLowStrikes(new ARM_GP_Vector());
    ARM_GP_VectorPtr fxHighStrikes(new ARM_GP_Vector());
    ARM_GP_VectorPtr fxLeverages(new ARM_GP_Vector());
    ARM_GP_Vector fundingPaymentDates;
    ARM_GP_Vector fundingInterestTerms;
    ARM_GP_VectorPtr fundingSpreads(new ARM_GP_Vector());

    while(floor(itsBoosterDatas(eventIdx,NOTICE_DATE)) <= asOfDate + K_NEW_DOUBLE_TOL)
        ++eventIdx;
    for(;eventIdx<nbEvents;++eventIdx)
    {
        fundingPaymentDates.push_back(  itsBoosterDatas(eventIdx,FUNDING_PAY_DATE));
        fundingInterestTerms.push_back( itsBoosterDatas(eventIdx,FUNDING_ACCRUAL_BASIS));
        fundingSpreads->push_back(      itsBoosterDatas(eventIdx,FUNDING_SPREAD));

        if(itsBoosterDatas(eventIdx,COUPON_DOMESTIC) == FOREX_CPN)
        {
            fxExpiryDates.push_back(  itsBoosterDatas(eventIdx,FX_RESET_DATE));
            fxPaymentDates.push_back( itsBoosterDatas(eventIdx,FX_PAY_DATE));
            fxLowStrikes->push_back(  itsBoosterDatas(eventIdx,FX_STRIKE));
            fxHighStrikes->push_back( itsBoosterDatas(eventIdx,FX_CAP));
            fxLeverages->push_back(   itsBoosterDatas(eventIdx,FX_NOTIO_MULTIPLIER) *
                                      itsBoosterDatas(eventIdx,FX_ACCRUAL_BASIS));
        }
    }


    ARM_StringVector profileNames;
    vector <ARM_GramFctorArg> profileValues;


    /// Funding datas
    ARM_DateStripPtr fundingDateStrip(new ARM_DateStrip());
    fundingDateStrip->SetPaymentDates(&fundingPaymentDates);
    fundingDateStrip->SetInterestTerms(&fundingInterestTerms);

    profileNames.push_back(PRDCProfileNamesTable[FundingDateStripProfile]);
    profileValues.push_back(ARM_GramFctorArg(fundingDateStrip));

    profileNames.push_back(PRDCProfileNamesTable[FundingSpreadProfile]);
    profileValues.push_back(ARM_GramFctorArg(fundingSpreads));


    /// Fx option datas
    ARM_DateStripPtr fxDateStrip(new ARM_DateStrip());
    fxDateStrip->SetResetDates(&fxExpiryDates);
    fxDateStrip->SetPaymentDates(&fxPaymentDates);

    profileNames.push_back(PRDCProfileNamesTable[FxDateStripProfile]);
    profileValues.push_back(ARM_GramFctorArg(fxDateStrip));

    profileNames.push_back(PRDCProfileNamesTable[FxLowStrikeProfile]);
    profileValues.push_back(ARM_GramFctorArg(fxLowStrikes));

    profileNames.push_back(PRDCProfileNamesTable[FxHighStrikeProfile]);
    profileValues.push_back(ARM_GramFctorArg(fxHighStrikes));

    profileNames.push_back(PRDCProfileNamesTable[FxLeverageProfile]);
    profileValues.push_back(ARM_GramFctorArg(fxLeverages));

    return ARM_CstManagerPtr(new ARM_CstManager(profileNames,profileValues));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if PRDC datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CheckData()
{
    if(!itsInputPowRev)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : a ARM_PowerReverse is needed to initialise PRDC calculator");

    if(!itsInputModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : a ARM_BFBSModel is needed to initialise PRDC calculator");

    /// Check if input columns to price are present in the deal description
    size_t colNamesSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        for(size_t j=0;j<colNamesSize;++j)
            if(itsColumnsToPrice[i] == PRDCColNamesTable[j])
                break;
        if(j==itsColumnsToPrice.size())
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : can't price this unknown column name " + itsColumnsToPrice[i]);  
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if PRDC market datas are consistent
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CheckMktData()
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



    ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    if(!oswDomBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + GetKeys()[OswDomModelKey] + " is expected in the Market Data Manager");

    ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    if(!oswForBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + GetKeys()[OswForModelKey] + " is expected in the Market Data Manager");

    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    if(!fxBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx B&S model for key=" + GetKeys()[FxModelKey] + " is expected in the Market Data Manager");



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
		ARM_BSModel* localFxBSModel=NULL;
		if(!GetMktDataManager()->TestIfKeyMissing(GetKeys()[LocalFxModelKey]))
			localFxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[LocalFxModelKey]) );
        if(!localFxBSModel)
            ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : local fx B&S model for key=" + GetKeys()[LocalFxModelKey] + " is expected in the Market Data Manager");
    }

	if(itsCalibType == HybridBasketCalib)
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
///	Class  : ARM_PRCSCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRCSCalculator::ColumnNames() const
{
    size_t colNamesSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);
    vector< string > colNamesVec(colNamesSize);
    vector< ARM_GP_VALUE_TYPE > colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0;i<colNamesSize; ++i)
        colNamesVec[i] = PRDCColNamesTable[i];

    ARM_RowInfo rowInfo(colNamesVec,colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_PRCSCalculator::DatesStructure() const
{
    /// Built a pseudo date strip combiner just to save notice dates (actual or not)
    size_t i,nbFlows = itsBoosterDatas.rows();

    /// Take care : notice dates in BoosterDatas are not integer
    /// but shifted by +0.625 then go back to the floor integer !
    ARM_GP_Vector noticeDates(nbFlows);
    ARM_GP_Vector interestTerms(nbFlows,1.0);
    for(i=0;i<nbFlows;++i)
        noticeDates[i] = floor(itsBoosterDatas(i,NOTICE_DATE));

    ARM_DateStrip pseudoSched(&noticeDates,&noticeDates,&noticeDates,&noticeDates,&noticeDates,
                                &noticeDates,&interestTerms,&interestTerms);

    /// Create a pseudo combiner to give a date structure to the generic security builder
    ARM_DateStripVector schedVect(1);
    schedVect[0] = &pseudoSched;
    ARM_DateStripCombiner EventSchedule(schedVect,"ResetDate");

    /// Save the first event line strictly posterior to as of date
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();
    itsFirstFutureEventIdx = 0;
	while(itsFirstFutureEventIdx < eventDates->size() && (*eventDates)[itsFirstFutureEventIdx] < asOfDate)
        ++itsFirstFutureEventIdx;

    return EventSchedule;
}



/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the PRDC
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_PRCSCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::InitPriceableColumns(vector< string >& rowDescVec,vector< ARM_GP_VALUE_TYPE >& rowTypeVec) const
{
    string zeroValue("0");


    rowDescVec[PRDCColAlias::FundingFlow] = zeroValue;
    rowTypeVec[PRDCColAlias::FundingFlow] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::Funding] = zeroValue;
    rowTypeVec[PRDCColAlias::Funding] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::FixedCpn] = zeroValue;
    rowTypeVec[PRDCColAlias::FixedCpn] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::FxLowCall] = zeroValue;
    rowTypeVec[PRDCColAlias::FxLowCall] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::FxHighCall] = zeroValue;
    rowTypeVec[PRDCColAlias::FxHighCall] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::FxLowCallStrip] = zeroValue;
    rowTypeVec[PRDCColAlias::FxLowCallStrip] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::FxHighCallStrip] = zeroValue;
    rowTypeVec[PRDCColAlias::FxHighCallStrip] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::CpnFlow] = zeroValue;
    rowTypeVec[PRDCColAlias::CpnFlow] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::Cpn] = zeroValue;
    rowTypeVec[PRDCColAlias::Cpn] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::Redemption] = zeroValue;
    rowTypeVec[PRDCColAlias::Redemption] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::PRDCFlow] = zeroValue;
    rowTypeVec[PRDCColAlias::PRDCFlow] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::PRDCFirstEuropean] = zeroValue;
    rowTypeVec[PRDCColAlias::PRDCFirstEuropean] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::PRDCFirstSwap] = zeroValue;
    rowTypeVec[PRDCColAlias::PRDCFirstSwap] = ARM_DOUBLE;

    rowDescVec[PRDCColAlias::PRDCOption] = zeroValue;
    rowTypeVec[PRDCColAlias::PRDCOption] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_PRCSCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{

    size_t eventSize = datesStructure.GetDateStrip(0)->GetResetDates()->size();
    size_t descSize = sizeof(PRDCColNamesTable)/sizeof(PRDCColNamesTable[0]);

    vector< string > rowDescVec(descSize);
    vector< ARM_GP_VALUE_TYPE > rowTypeVec(descSize, ARM_MISSING_TYPE); 


    /// EventDate description
    double eventDate=(*(datesStructure.GetMergeData()))[eventIdx];
    CC_Ostringstream eventDateDesc;
    eventDateDesc << CC_NS(std,fixed) << eventDate;
    rowDescVec[PRDCColAlias::EventDate] = eventDateDesc.str();
    rowTypeVec[PRDCColAlias::EventDate] = ARM_DATE_TYPE;


    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec,rowTypeVec);


    /// Get the model names for domestic and forex models
    string basisDomModelName = GetKeys()[YcBasisDomKey];
    string fxModelName = GetKeys()[ForexKey];


    string nextExerIdx("[i+1]");

    bool isLastEvent = (eventIdx+1>=eventSize);
    bool isFirstEvent = (eventIdx==itsFirstFutureEventIdx);

    /// Funding description : in the booster table the spread is already converted
    /// in the domestic ccy then the leg is described as a simple (but fictive)
    /// domestic floating leg

    CC_Ostringstream fundStartDesc;
    fundStartDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FUNDING_START_DATE);
    rowDescVec[PRDCColAlias::FundingStartDate] = fundStartDesc.str();
    rowTypeVec[PRDCColAlias::FundingStartDate] = ARM_DATE_TYPE;

    CC_Ostringstream fundEndDesc;
    fundEndDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FUNDING_END_DATE);
    rowDescVec[PRDCColAlias::FundingEndDate] = fundEndDesc.str();
    rowTypeVec[PRDCColAlias::FundingEndDate] = ARM_DATE_TYPE;

    CC_Ostringstream fundPayDesc;
    fundPayDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FUNDING_PAY_DATE);
    rowDescVec[PRDCColAlias::FundingPayDate] = fundPayDesc.str();
    rowTypeVec[PRDCColAlias::FundingPayDate] = ARM_DATE_TYPE;

    CC_Ostringstream fundSpreadDesc;
    fundSpreadDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FUNDING_SPREAD);
    rowDescVec[PRDCColAlias::FundingSpread] = fundSpreadDesc.str();
    rowTypeVec[PRDCColAlias::FundingSpread] = ARM_DOUBLE;

    CC_Ostringstream fundITDesc;
    fundITDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FUNDING_ACCRUAL_BASIS);
    rowDescVec[PRDCColAlias::FundingIT] = fundITDesc.str();
    rowTypeVec[PRDCColAlias::FundingIT] = ARM_DOUBLE;

    CC_Ostringstream fundingFlowDesc;
    fundingFlowDesc << "DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::FundingStartDate] << "[i])";
    fundingFlowDesc << "-DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::FundingEndDate] << "[i])+";
    fundingFlowDesc << PRDCColNamesTable[PRDCColAlias::FundingSpread] << "[i]*" << PRDCColNamesTable[PRDCColAlias::FundingIT] << "[i]";
    fundingFlowDesc << "*DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::FundingPayDate] << "[i])";

    rowDescVec[PRDCColAlias::FundingFlow] = fundingFlowDesc.str();
    rowTypeVec[PRDCColAlias::FundingFlow] = ARM_STRING;

    if(itsFxLocalModelFlag)
    {
        CC_Ostringstream fundEndDesc;
        fundEndDesc << CC_NS(std,fixed) << itsBoosterDatas(eventSize-1,FUNDING_END_DATE);
        rowDescVec[PRDCColAlias::FundingLastEndDate] = fundEndDesc.str();
        rowTypeVec[PRDCColAlias::FundingLastEndDate] = ARM_DATE_TYPE;

        /// Spread profile in input through notional argument
        CC_Ostringstream fundingDesc;
        fundingDesc << "DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::FundingStartDate] << "[i])";
        fundingDesc << "-DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::FundingLastEndDate] << "[i])";
        fundingDesc << "+ANNUITY(" << basisDomModelName << ",";
        fundingDesc << PRDCColNamesTable[PRDCColAlias::FundingStartDate] << "[i],";
        fundingDesc << PRDCColNamesTable[PRDCColAlias::FundingStartDate] << "[i],,,";
        fundingDesc << PRDCProfileNamesTable[PRDCProfileType::FundingSpreadProfile] << ",";
        fundingDesc << PRDCProfileNamesTable[PRDCProfileType::FundingDateStripProfile] << ",";
        fundingDesc << eventIdx << ")";

        rowDescVec[PRDCColAlias::Funding] = fundingDesc.str();
        rowTypeVec[PRDCColAlias::Funding] = ARM_STRING;
    }


    CC_Ostringstream cpnITDesc;
    cpnITDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_ACCRUAL_BASIS);
    rowDescVec[PRDCColAlias::CpnIT] = cpnITDesc.str();
    rowTypeVec[PRDCColAlias::CpnIT] = ARM_DOUBLE;

    CC_Ostringstream cpnPayDesc;
    cpnPayDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_PAY_DATE);
    rowDescVec[PRDCColAlias::CpnPayDate] = cpnPayDesc.str();
    rowTypeVec[PRDCColAlias::CpnPayDate] = ARM_DATE_TYPE;

    CC_Ostringstream cpnLeverageDesc;
    cpnLeverageDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_NOTIO_MULTIPLIER);
    rowDescVec[PRDCColAlias::CpnLeverage] = cpnLeverageDesc.str();
    rowTypeVec[PRDCColAlias::CpnLeverage] = ARM_DOUBLE;


    /// Coupon description : flow is always described but residual coupons
    /// are described only for local model case

    /// FIX FIX : if there is a fixed coupon part, local model case should be adapted
    /// to described a fixed leg (ANNUITY) then a fx option strip (CALLSTRIP)
    if(itsBoosterDatas(eventIdx,COUPON_DOMESTIC) == FIXED_CPN && itsFxLocalModelFlag)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mixing of fixed & fx cpn with local vol model not supported");

    CC_Ostringstream cpnFlowDesc;
    CC_Ostringstream cpnDesc;
    if(itsBoosterDatas(eventIdx,COUPON_DOMESTIC) == FIXED_CPN)
    {
        /// Fixed flow description.
        CC_Ostringstream fixedCpnDesc;
        fixedCpnDesc << CC_NS(std,fixed) << PRDCColNamesTable[PRDCColAlias::FixedCpn] << "[i]*";
        fixedCpnDesc << PRDCColNamesTable[PRDCColAlias::CpnLeverage] << "[i]*";
        fixedCpnDesc << PRDCColNamesTable[PRDCColAlias::CpnIT] << "[i]*";
        fixedCpnDesc << "DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::CpnPayDate] << ")";
        rowDescVec[PRDCColAlias::FixedCpn] = fixedCpnDesc.str();
        rowTypeVec[PRDCColAlias::FixedCpn] = ARM_STRING;

        cpnFlowDesc << PRDCColNamesTable[FixedCpn] << "[i]";
    }
    else if(itsBoosterDatas(eventIdx,COUPON_DOMESTIC) == FOREX_CPN)
    {
        CC_Ostringstream fxResetDesc;
        fxResetDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_RESET_DATE);
        rowDescVec[PRDCColAlias::FxResetDate] = fxResetDesc.str();
        rowTypeVec[PRDCColAlias::FxResetDate] = ARM_DATE_TYPE;

        CC_Ostringstream fxLowStrikeDesc;
        fxLowStrikeDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_STRIKE);
        rowDescVec[PRDCColAlias::FxLowStrike] = fxLowStrikeDesc.str();
        rowTypeVec[PRDCColAlias::FxLowStrike] = ARM_DOUBLE;

        CC_Ostringstream fxHighStrikeDesc;
        fxHighStrikeDesc << CC_NS(std,fixed) << itsBoosterDatas(eventIdx,FX_CAP);
        rowDescVec[PRDCColAlias::FxHighStrike] = fxHighStrikeDesc.str();
        rowTypeVec[PRDCColAlias::FxHighStrike] = ARM_DOUBLE;

        string FlooredFxModelName(fxModelName);
        string CappedFxModelName(fxModelName);
        string RedemptionFxModelName(fxModelName);
        if(itsFxLocalModelFlag)
        {
            FlooredFxModelName      = GetKeys()[FlooredFxModelKey];
            CappedFxModelName       = GetKeys()[CappedFxModelKey];
            RedemptionFxModelName   = GetKeys()[RedemptionFxModelKey];
        }

        /// Current FX option description
        CC_Ostringstream fxLowCapDesc;
        fxLowCapDesc << CC_NS(std,fixed) << PRDCColNamesTable[PRDCColAlias::CpnLeverage] << "[i]*";
        fxLowCapDesc << PRDCColNamesTable[PRDCColAlias::CpnIT] << "[i]*";
        fxLowCapDesc << "CALL(" << FlooredFxModelName << "," << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
        fxLowCapDesc << PRDCColNamesTable[PRDCColAlias::FxLowStrike] << "[i],C,,";
        fxLowCapDesc << PRDCColNamesTable[PRDCColAlias::CpnPayDate] << "[i])";
        rowDescVec[PRDCColAlias::FxLowCall] = fxLowCapDesc.str();
        rowTypeVec[PRDCColAlias::FxLowCall] = ARM_STRING;

        if(itsBoosterDatas(eventIdx,FX_CAP) < MAX_FX_STRIKE)
        {
            CC_Ostringstream fxHighCapDesc;
            fxHighCapDesc << CC_NS(std,fixed) << PRDCColNamesTable[PRDCColAlias::CpnLeverage] << "[i]*";
            fxHighCapDesc << PRDCColNamesTable[PRDCColAlias::CpnIT] << "[i]*";
            fxHighCapDesc << "CALL(" << CappedFxModelName << "," << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
            fxHighCapDesc << PRDCColNamesTable[PRDCColAlias::FxHighStrike] << "[i],C,,";
            fxHighCapDesc << PRDCColNamesTable[PRDCColAlias::CpnPayDate] << "[i])";
            rowDescVec[PRDCColAlias::FxHighCall] = fxHighCapDesc.str();
            rowTypeVec[PRDCColAlias::FxHighCall] = ARM_STRING;
        }

        cpnFlowDesc << PRDCColNamesTable[PRDCColAlias::FxLowCall] << "[i]-" << PRDCColNamesTable[PRDCColAlias::FxHighCall] << "[i]";

        if(itsFxLocalModelFlag)
        {
            /// Residual FX option sum description with CALLSTRIP keyword
            CC_Ostringstream fxLowCallStripDesc;
            fxLowCallStripDesc << CC_NS(std,fixed) << "CALLSTRIP(" << FlooredFxModelName << ",";
            fxLowCallStripDesc << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
            fxLowCallStripDesc << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
            fxLowCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxLowStrikeProfile] << ",C,,,,,,";
            fxLowCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxLeverageProfile] << ",";
            fxLowCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxDateStripProfile] << ",";
            fxLowCallStripDesc << eventIdx << ")";
            rowDescVec[PRDCColAlias::FxLowCallStrip] = fxLowCallStripDesc.str();
            rowTypeVec[PRDCColAlias::FxLowCallStrip] = ARM_STRING;

            /// Check if there is a significative residual high strike            
            for(size_t idx=eventIdx;idx<eventSize;++idx)
            {
                if(itsBoosterDatas(idx,FX_CAP) < MAX_FX_STRIKE)
                {
                    CC_Ostringstream fxHighCallStripDesc;
                    fxHighCallStripDesc << CC_NS(std,fixed) << "CALLSTRIP(" << CappedFxModelName << ",";
                    fxHighCallStripDesc << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
                    fxHighCallStripDesc << PRDCColNamesTable[PRDCColAlias::FxResetDate] << "[i],";
                    fxHighCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxHighStrikeProfile] << ",C,,,,,,";
                    fxHighCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxLeverageProfile] << ",";
                    fxHighCallStripDesc << PRDCProfileNamesTable[PRDCProfileType::FxDateStripProfile] << ",";
                    fxHighCallStripDesc << eventIdx << ")";
                    rowDescVec[PRDCColAlias::FxHighCallStrip] = fxHighCallStripDesc.str();
                    rowTypeVec[PRDCColAlias::FxHighCallStrip] = ARM_STRING;
                    break;
                }
            }
            cpnDesc << PRDCColNamesTable[PRDCColAlias::FxLowCallStrip] << "[i]-" << PRDCColNamesTable[PRDCColAlias::FxHighCallStrip] << "[i]";
        }


        if((isLastEvent || itsFxLocalModelFlag) && itsRedemptionDatas.size()>0)
        {
            /// Redemption flow description
            CC_Ostringstream redResetDesc;
            redResetDesc << CC_NS(std,fixed) << itsRedemptionDatas[REDEMPTION_RESET_DATE];
            rowDescVec[PRDCColAlias::RedemptionResetDate] = redResetDesc.str();
            rowTypeVec[PRDCColAlias::RedemptionResetDate] = ARM_DATE_TYPE;

            CC_Ostringstream redPayDesc;
            redPayDesc << CC_NS(std,fixed) << itsRedemptionDatas[REDEMPTION_PAY_DATE];
            rowDescVec[PRDCColAlias::RedemptionPayDate] = redPayDesc.str();
            rowTypeVec[PRDCColAlias::RedemptionPayDate] = ARM_DATE_TYPE;

            CC_Ostringstream redStrikeDesc;
            redStrikeDesc << CC_NS(std,fixed) << itsRedemptionDatas[REDEMPTION_STRIKE];
            rowDescVec[PRDCColAlias::RedemptionStrike] = redStrikeDesc.str();
            rowTypeVec[PRDCColAlias::RedemptionStrike] = ARM_DOUBLE;

            CC_Ostringstream redemptionDesc;
            if(itsRedemptionDatas[REDEMPTION_TYPE] == MANDATORY_REDEMPTION)
            {
                redemptionDesc << "(FWD(" << RedemptionFxModelName << "," << PRDCColNamesTable[PRDCColAlias::RedemptionResetDate] << "[i],,";
                redemptionDesc << PRDCColNamesTable[PRDCColAlias::RedemptionPayDate] << "[i])/";
                redemptionDesc << PRDCColNamesTable[PRDCColAlias::RedemptionStrike] << "[i]" << "-1)*";
                redemptionDesc << "DF(" << basisDomModelName << "," << PRDCColNamesTable[PRDCColAlias::RedemptionPayDate] << "[i])";

            }
            else if(itsRedemptionDatas[REDEMPTION_TYPE] == DUAL_OPTION_REDEMPTION)
            {
                redemptionDesc << "-CALL(" << RedemptionFxModelName << "," << PRDCColNamesTable[PRDCColAlias::RedemptionResetDate] << "[i],";
                redemptionDesc << PRDCColNamesTable[PRDCColAlias::RedemptionStrike] << "[i]" << ",P,,";
                redemptionDesc << PRDCColNamesTable[PRDCColAlias::RedemptionPayDate] << "[i])/";
                redemptionDesc << PRDCColNamesTable[PRDCColAlias::RedemptionStrike] << "[i]";

            }
            else
                ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : unknown redemption type");

            if(isLastEvent)
            {
                rowDescVec[PRDCColAlias::Redemption] = redemptionDesc.str();
                rowTypeVec[PRDCColAlias::Redemption] = ARM_STRING;
                cpnFlowDesc << "+" << PRDCColNamesTable[PRDCColAlias::Redemption] << "[i]";
            }
            if(itsFxLocalModelFlag)
            {
                /// Local model => always add redemption flow
                /// (because of local convexity adjustment and/or local volatility)
                if(isLastEvent)
                    cpnDesc << "+" << PRDCColNamesTable[PRDCColAlias::Redemption] << "[i]";
                else
                    cpnDesc << "+" << redemptionDesc.str();
            }
        }
    }
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : PRDC coupon must be either a fixed rate or a FX option");

    rowDescVec[PRDCColAlias::CpnFlow] = cpnFlowDesc.str();
    rowTypeVec[PRDCColAlias::CpnFlow] = ARM_STRING;


    /// PRD flow description
    CC_Ostringstream prdcFlowDesc;
    prdcFlowDesc << PRDCColNamesTable[PRDCColAlias::CpnFlow] << "[i]-" << PRDCColNamesTable[PRDCColAlias::FundingFlow] << "[i]";
    rowDescVec[PRDCColAlias::PRDCFlow] = prdcFlowDesc.str();
    rowTypeVec[PRDCColAlias::PRDCFlow] = ARM_STRING;


    /// PRD swap description
    CC_Ostringstream prdcSwapDesc;
    if(itsFxLocalModelFlag)
    {
        /// Swap = residual PRD swap
        rowDescVec[PRDCColAlias::Cpn] = cpnDesc.str();
        rowTypeVec[PRDCColAlias::Cpn] = ARM_STRING;

        prdcSwapDesc << PRDCColNamesTable[PRDCColAlias::Cpn] << "[i]-" << PRDCColNamesTable[PRDCColAlias::Funding] << "[i]";
    }
    else
    {
        /// Swap = PV of PRDC flows
        if(!isLastEvent) 
            prdcSwapDesc << "PV(" << PRDCColNamesTable[PRDCColAlias::PRDCSwap] << nextExerIdx << ")+";
        prdcSwapDesc << PRDCColNamesTable[PRDCColAlias::PRDCFlow] << "[i]";
    }
    rowDescVec[PRDCColAlias::PRDCSwap] = prdcSwapDesc.str();
    rowTypeVec[PRDCColAlias::PRDCSwap] = ARM_STRING;


    /// PRDC option description
    CC_Ostringstream prdcBermudaDesc;
    if(itsExerciseFlag[eventIdx])
    {
        /// Actual exercise at this event date
        prdcBermudaDesc << "EXERCISE(0," << PRDCColNamesTable[PRDCColAlias::PRDCSwap] << "[i],";
        if(isLastEvent)
            /// No residual option
            prdcBermudaDesc << "0)";
        else
            /// Only one description to get the bermuda price
            prdcBermudaDesc << PRDCColNamesTable[PRDCColAlias::PRDCBermuda] << nextExerIdx << ")";
    }
    else
    {
        /// No exercise allowed at this event date
        if(isLastEvent)
            /// No residual option
            prdcBermudaDesc << "0";
        else
            /// Just actualise residual option
            prdcBermudaDesc << "PV(" << PRDCColNamesTable[PRDCColAlias::PRDCBermuda] << nextExerIdx << ")";
    }

    rowDescVec[PRDCColAlias::PRDCBermuda] = prdcBermudaDesc.str();
    rowTypeVec[PRDCColAlias::PRDCBermuda] = ARM_STRING;

    /// The single PRDC option a first line to get the price
    if(isFirstEvent)
    {
        CC_Ostringstream prdcFirstSwapDesc;
        prdcFirstSwapDesc << PRDCColNamesTable[PRDCColAlias::PRDCSwap] << "[i]";
        rowDescVec[PRDCColAlias::PRDCFirstSwap] = prdcFirstSwapDesc.str();
        rowTypeVec[PRDCColAlias::PRDCFirstSwap] = ARM_STRING;

        CC_Ostringstream prdcFirstEuropeanDesc;
        prdcFirstEuropeanDesc << "EXERCISE(0," << PRDCColNamesTable[PRDCColAlias::PRDCSwap] << "[i],0)";
        rowDescVec[PRDCColAlias::PRDCFirstEuropean] = prdcFirstEuropeanDesc.str();
        rowTypeVec[PRDCColAlias::PRDCFirstEuropean] = ARM_STRING;

        CC_Ostringstream prdcOptionDesc;
        prdcOptionDesc << PRDCColNamesTable[PRDCColAlias::PRDCBermuda] << "[i]";
        rowDescVec[PRDCColAlias::PRDCOption] = prdcOptionDesc.str();
        rowTypeVec[PRDCColAlias::PRDCOption] = ARM_STRING;
    }


    return ARM_RowInfo(rowDescVec,rowTypeVec);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CreateAndSetModel()
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

    itsDomModelHWFlag=false;
    if(GetMktDataManager()->TestIfKeyMissing(GetKeys()[QDomKey]))
    {
        /// Degenerate to H&W
        itsDomModelHWFlag = true;
        modelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
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
        modelParams[2] =  new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
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
        modelParams[2] =  new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 1.0,"Q");
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
///	Class  : ARM_PRCSCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create the calibration
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CreateAndSetCalibration()
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=false;

    bool isFreezeWeights = false;   // reinit non null weights in the portfolio
    bool isInitVol = true;          // init vol param for calibration (bounds & init guess)


	/// Hybrid basket swaptions (may be also used for domestic calibration)
    /// -------------------------------------------------------------------

	if(itsCalibType == HybridBasketCalib)
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
    ComputeIROptionPrices(domCalib,OswDomModelKey,itsDomModelHWFlag,isFreezeWeights,isInitVol,
		!(IS_DOM_BASKET_CALIB && (itsCalibType == HybridBasketCalib)));


    /// Diagonal foreign swaptions
    /// --------------------------

    /// Build an empty calib method for latter foreign volatility bootstrapping
	ARM_Currency* forCcy = static_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcForKey]))->GetCurrencyUnit();
    ARM_CalibMethod* forCalib = CreateIRCalibration(diagonalSwaptionPF.second,forCcy,
                                    itsBasisIRCalibFlag ? ARM_2IRFXModel::ForBasisModel
                                                        : ARM_2IRFXModel::ForModel);

    /// Compute foreign target prices
    ComputeIROptionPrices(forCalib,OswForModelKey,itsForModelHWFlag,isFreezeWeights,isInitVol);


    /// Forex calibration
    /// -----------------


    /// Create a Fx option portfolio and compute Fx option target prices
    ARM_StdPortfolioPtr fxOptionPF;
    ARM_CalibMethod* fxCalib;
	if(itsCalibType == ATMDoubleCalib)
	{
		/// Create a standard ATM Fx option portfolio for first calibration step
		itsCalibType    = ATMCalib;
		fxOptionPF      = CreateFxOption();
		fxCalib         = CreateFxCalibration(fxOptionPF,ARM_2IRFXModel::FxModel);
		ComputeFxOptionPrices(fxCalib,isFreezeWeights,isInitVol);
		itsCalibType    = ATMDoubleCalib;
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

    if(itsCalibType == ATMDoubleCalib)
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
///	Class  : ARM_PRCSCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions for domestic market
/////////////////////////////////////////////////////////////////
pair< ARM_StdPortfolioPtr, ARM_StdPortfolioPtr > ARM_PRCSCalculator::CreateDiagonalSwaption()
{
	bool isDomDiagCalib = !(IS_DOM_BASKET_CALIB && (itsCalibType == HybridBasketCalib));

    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    /// Get market swaption standard expiries for domestic & foreign volatilities
    size_t i,nbExp;
    ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswDomModelKey]) );
    ARM_VolCurve* oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector domStdExp(nbExp+1);
    domStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        domStdExp[i+1] = asOfDate + 365.0 * (*(oswBSVol->GetExpiryTerms()))[i];

    oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[OswForModelKey]) );
    oswBSVol = oswBSModel->GetVolatility();
    nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector forStdExp(nbExp+1);
    forStdExp[0] = asOfDate;
    for(i=0;i<nbExp;++i)
        forStdExp[i+1] = asOfDate + 365.0 * (*(oswBSVol->GetExpiryTerms()))[i];


    ARM_Swaption* swaption;

    list < ARM_Security* > domSwaptionList;
    list < ARM_Security* > forSwaptionList;


    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
    size_t nbFutureExer=itsExerciseFlag.size() - itsFirstFutureEventIdx;
    if(nbFutureExer+1 != nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mismatch between exercise flags & notice dates in the deal description");

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
    ARM_Date swapEndDate = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,FundingEndDate).c_str()) );
    bool isDomSwap,isForSwap;
    ARM_Swap domSwap,forSwap;

    /// The 1st notice date is always added then get the very last expiry
    double lastNotice,lastInsertedNotice;
    double is1stNoticeSwap=false;
    size_t firstNoticeIdx=1;
    for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
        if(itsExerciseFlag[eventIdx-1+itsFirstFutureEventIdx])
        {
            expiryDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );
            if(!is1stNoticeSwap)
            {
				if(isDomDiagCalib)
				{
					/// Create the 1st domestic diagonal swaption
					swapStartDate   = expiryDate;
					swapStartDate.GapBusinessDay(domSpotDays,domResetCal);

					domSwap = ARM_Swap(swapStartDate,swapEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
									 K_DEF_FREQ,K_DEF_FREQ,domCcy);
					swaption = new ARM_Swaption(&domSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
					domSwaptionList.push_back(static_cast< ARM_Security* >(swaption));
				}

				/// Create the 1st foreign diagonal swaption
                swapStartDate   = expiryDate;
                swapStartDate.GapBusinessDay(forSpotDays,forResetCal);

                forSwap = ARM_Swap(swapStartDate,swapEndDate,forLiborIndex,0.0,K_MARKET_RATE,K_RCV,
                                 K_DEF_FREQ,K_DEF_FREQ,forCcy);
                swaption = new ARM_Swaption(&forSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);
                forSwaptionList.push_back(static_cast< ARM_Security* >(swaption));

                firstNoticeIdx = eventIdx;
                lastInsertedNotice = expiryDate.GetJulian();

                is1stNoticeSwap=true;
            }
            lastNotice = expiryDate.GetJulian();
        }
    }
    if(!is1stNoticeSwap)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " creation of calibration diagonal swaption not possible");


    /// Between each standard market expiries, only one swaption is created corresponding to
    /// the earlier notice date
    /// To be consistent with Tree3F code (may be it is a bug !), the very last notice date
    /// is excluded
    size_t domExpIdx=0,forExpIdx=0;
    for(eventIdx = firstNoticeIdx+1;eventIdx<nbEvents;++eventIdx)
    {
        expiryDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );

        if(itsExerciseFlag[eventIdx-1+itsFirstFutureEventIdx] && expiryDate.GetJulian() < lastNotice)  /// to be consistent to Tree3F code !!
        {
			isDomSwap=false;
			if(isDomDiagCalib)
			{
				while(!isDomSwap && domExpIdx < domStdExp.size() && domStdExp[domExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
				{
					if( ( domExpIdx+1 == domStdExp.size() ||
						  (domExpIdx+1 < domStdExp.size() && expiryDate.GetJulian() < domStdExp[domExpIdx+1] + K_NEW_DOUBLE_TOL) )
						&& lastInsertedNotice + K_NEW_DOUBLE_TOL < domStdExp[domExpIdx] )
					{
						/// Create a domestic diagonal swap
						isDomSwap=true;
						swapStartDate   = expiryDate;
						swapStartDate.GapBusinessDay(domSpotDays,domResetCal);
						domSwap = ARM_Swap(swapStartDate,swapEndDate,domLiborIndex,0.0,K_MARKET_RATE,K_RCV,
										 K_DEF_FREQ,K_DEF_FREQ,domCcy);
					}

					++domExpIdx;
				}
			}

            isForSwap=false;
            while(!isForSwap && forExpIdx < forStdExp.size() && forStdExp[forExpIdx] + K_NEW_DOUBLE_TOL < expiryDate.GetJulian())
            {
                if( ( forExpIdx+1 == forStdExp.size() ||
                      (forExpIdx+1 < forStdExp.size() && expiryDate.GetJulian() < forStdExp[forExpIdx+1] + K_NEW_DOUBLE_TOL) )
                    && lastInsertedNotice + K_NEW_DOUBLE_TOL < forStdExp[forExpIdx] )
                {
                    /// Create a foreign diagonal swap
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

	if(!isDomDiagCalib)
	{
		/// Restore and clone domestic variable notional swaptions of hybrid IR/FX swaptions
		ARM_Portfolio* hybridBasket;
		ARM_StdPortfolioPtr hybridPort = itsHybridBasketCalibMethod->GetPortfolio();
		for(i=0;i<hybridPort->size();++i)
		{
			/// VN swaption is the last security of the portfolio
			hybridBasket = (static_cast<ARM_OptionPortfolio*>(hybridPort->GetAsset(i)))->GetPtf();
			domSwaptionList.push_back(static_cast<ARM_Security*>(hybridBasket->GetAsset(hybridBasket->GetSize()-1)->Clone()));
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
///	Class  : ARM_PRCSCalculator
///	Routine: CreateIRCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRCSCalculator::CreateIRCalibration(const ARM_StdPortfolioPtr& diagonalSwaptionPF,
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
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeIROptionPrice
///	Returns: nothing
///	Action : compute market target prices for calibration purpose
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputeIROptionPrices(ARM_CalibMethod* calibMethod, 
	          mdmKeysAlias oswModelIdx, 
			  bool isModelHW, 
			  bool isFreezeWeights, 
			  bool isInitVolParam,
			  bool isDiagCalib)
{
    ARM_BSModel* oswBSModel	= dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[oswModelIdx]) );
	ARM_MarketIRModel* oswMktModel = NULL;
	if(!GetMktDataManager()->TestIfKeyMissing(GetKeys()[MarketIrModelKey]))
		oswMktModel = dynamic_cast< ARM_MarketIRModel* >( GetMktDataManager()->GetData(GetKeys()[MarketIrModelKey]) );
	bool isBasketCalib = (!isDiagCalib && oswMktModel);

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Restore calibration portfolios
    ARM_CalibMethod* volCalibMethod;
    ARM_StdPortfolioPtr oswPortfolio;

    double price,vega,nominal,weight;
    ARM_Swaption* swaption;
	ARM_VanillaArg* vnSwaption;
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

		if(isBasketCalib)
		{
			/// Only GP is able to compute the price
			vnSwaption = ARM_ConverterFromKernel::ConvertSecuritytoArgObject(swaption, asOfDate.GetJulian());
			price = vnSwaption->Price(oswMktModel);
		}

		nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(oswPortfolio->GetWeights()))[i];

        isNotSelected = vega < IR_VEGA_MIN_TO_SELECT * nominal && !isFreezeWeights;

        oswPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        oswPortfolio->SetPrecision(0.001*vega,i);
        oswPortfolio->SetPrice(price,i);

        if(isInitVol)
        {
            /// Vol initialisation
            initTimes[i]    = swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian();

            optMat  = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
            swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
            strike  = swaption->GetStrike();
			volATM	= isBasketCalib ? oswMktModel->GetVnsBasketVol() :
						oswBSModel->ComputeVol(optMat,swapMat,strike,strike)/100.0;

            if(isModelHW)
            {
				swapRate = isBasketCalib ? oswMktModel->GetVnsForward() : swaption->CptMarketSwapRate();
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
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeFxEquivStrike
///	Returns: a doubke
///	Action : compute an equivalent strike for FX options
/////////////////////////////////////////////////////////////////
double ARM_PRCSCalculator::ComputeFxEquivStrike(size_t eventIdx)
{
    string columnText,fxCallText("CALL");

    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names

	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
	ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisForKey]));
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    double spotFx       = fxBSModel->GetSpot();

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );
	string settlCal  = forex->ComputeSettlementCalendar();
	int settlGap = forex->ComputeSettlementGap();

    ARM_Date fundStartDate( atof(dealDesc.GetElem(eventIdx,FundingStartDate).c_str()) );
    double fundingLeg = basisDomCurve->DiscountPrice(fundStartDate);

    ARM_Date payDate,fxResetDate,fxSetDate;

    ARM_Option* fxOption;

	/// Compute spot df ratio
	ARM_Date spotDate(basisDomCurve->GetAsOfDate());
    spotDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));
    double zcDomSet	= basisDomCurve->DiscountPrice(spotDate);
    double zcForSet	= basisForCurve->DiscountPrice(spotDate);
	double spotDfRatio = zcDomSet/zcForSet;

    double fxStrike,price,delta,period,zcPay,zcFwd,fundingIT,fundingSpread;
    double firstZcFwd=1.0;
    double fxSensi=0.0,strikeLeg=0.0;
	double fx,stdDev,d1,strikePrice;
    for(size_t i=eventIdx;i<nbEvents;++i)
    {
        columnText = dealDesc.GetElem(eventIdx,FxLowCall);
        if(columnText.find(fxCallText) != string::npos)
		{
			period      = atof(dealDesc.GetElem(i,CpnIT).c_str()) * atof(dealDesc.GetElem(i,CpnLeverage).c_str());

			/// Compute low strike delta
			fxResetDate = ARM_Date( atof(dealDesc.GetElem(i,FxResetDate).c_str()) );
			payDate     = ARM_Date( atof(dealDesc.GetElem(i,CpnPayDate).c_str()) );
			fxStrike    = atof(dealDesc.GetElem(i,FxLowStrike).c_str());
			fxOption    = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
			fxOption->SetPayDate(payDate);
			fxOption->SetModel(fxBSModel);
			price       = fxOption->ComputePrice();

			fxSetDate   = fxResetDate;
			fxSetDate.NextBusinessDay(settlGap,(char *)(settlCal.c_str()));
			zcDomSet    = basisDomCurve->DiscountPrice(fxSetDate);
			zcForSet    = basisForCurve->DiscountPrice(fxSetDate);
			zcFwd		= zcForSet/zcDomSet;

			fx			= fxOption->GetCalcFwd();
			stdDev		= 0.01*fxOption->GetCalcVol()*sqrt(fxOption->GetCalcMat());
			d1			= log(fx/fxStrike)/stdDev + 0.5*stdDev;
			zcPay		= fxOption->GetCalcO1();
			delta		= cdfNormal(d1);
			strikePrice = cdfNormal(d1-stdDev)*fxStrike;

			columnText = dealDesc.GetElem(eventIdx,FxHighCall);
			if(columnText.find(fxCallText) != string::npos)
			{
				/// Compute high strike delta
				fxStrike	= atof(dealDesc.GetElem(i,FxHighStrike).c_str());
				d1			= log(fx/fxStrike)/stdDev + 0.5*stdDev;
				delta		-= cdfNormal(d1);
				strikePrice -= cdfNormal(d1-stdDev)*fxStrike;
			}
			delta		*= period*spotDfRatio*zcFwd*zcPay;
			strikePrice *= period*zcPay;

			/// Update Sensi to 1st fwd FX
			/// Convexity adjustment if settl != pay is assumed to be
			/// the same on PRD swap and strike equivalent FX flow
			fxSensi     += delta;
			strikeLeg	+= strikePrice;

			if(i==eventIdx)
				firstZcFwd = zcFwd;

			/// Update funding leg price
			payDate         = ARM_Date( atof(dealDesc.GetElem(i,FundingPayDate).c_str()) );
			fundingIT       = atof(dealDesc.GetElem(i,FundingIT).c_str());
			fundingSpread   = atof(dealDesc.GetElem(i,FundingSpread).c_str());
			zcPay           = basisDomCurve->DiscountPrice(payDate);
			fundingLeg      += fundingSpread*fundingIT*zcPay;

			/// Free memory
			delete fxOption;
		}
    }

    /// Terminal condition flow
    if(itsRedemptionDatas.size()>0)
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
        if(itsRedemptionDatas[REDEMPTION_TYPE] == MANDATORY_REDEMPTION)
        {
            /// Mandatory termination => Fwd/K of the nominal is paid back
            /// <=> deltaPut=-1 in dual option case
            fxSensi     += spotDfRatio*zcFwd*zcPay/fxStrike;
            strikeLeg   += zcPay;
        }
        else
        {
            /// Dual option termination : compute delta put
            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
            fxOption->SetPayDate(payDate);
	        fxOption->SetModel(fxBSModel);
            price       = fxOption->ComputePrice();

			fx		= fxOption->GetCalcFwd();
			stdDev	= 0.01*fxOption->GetCalcVol()*sqrt(fxOption->GetCalcMat());
			d1		= log(fx/fxStrike)/stdDev + 0.5*stdDev;
			zcPay	= fxOption->GetCalcO1();
			delta	= cdfNormal(-d1)*spotDfRatio*zcFwd*zcPay;

            fxSensi		+= delta/fxStrike;
            strikePrice	= cdfNormal(-d1+stdDev)*zcPay;
            strikeLeg   += strikePrice;

            /// Free memory
            delete fxOption;
        }

    }

	/// Leg tests
	double fxLeg = fxSensi*spotFx;

    fxSensi     /= (firstZcFwd*spotDfRatio);
    fundingLeg  -= zcPay;

    return (strikeLeg+fundingLeg)/fxSensi;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_PRCSCalculator::CreateFxOption()
{
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    ARM_VolCurve* fxBSVol = fxBSModel->GetVolatility();
    size_t i, nbExp = fxBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector fxStdExp(nbExp);
    for(i=0;i<nbExp;++i)
        fxStdExp[i] = asOfDate + K_YEAR_LEN * (*(fxBSVol->GetExpiryTerms()))[i];
    
    ARM_Option* fxOption;

    list < ARM_Security* > fxOptionList;


    const ARM_DealDescription& dealDesc = GetGenSecurity()->GetDealDescription();
    size_t nbEvents=dealDesc.GetRowsNb(); // + 1 because of column names
    size_t nbFutureExer=itsExerciseFlag.size() - itsFirstFutureEventIdx;
    if(nbFutureExer+1 != nbEvents)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : mismatch between exercise flags & notice dates in the deal description");

    /// Restore forex
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(GetKeys()[ForexKey]) );


    if(itsCalibType == ATMCalib)
    {
        /// Portfolio = market standard FX options (ATM)

        /// Check if the 1st notice date may replace the 1st standard fx option expiry
        ARM_Date firstNoticeDate;
        size_t fxExpIdx=0;
        for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
        {
            if(itsExerciseFlag[eventIdx-1+itsFirstFutureEventIdx])
            {
                firstNoticeDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );
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
        string columnText,fxCallText("CALL");
        ARM_Date fxResetDate,fxPayDate;
        double fxStrike;
        size_t strikeIdx=0,nbStrikes = itsCalibDatas.size();
        for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
        {
            columnText = dealDesc.GetElem(eventIdx,FxLowCall);
            if(columnText.find(fxCallText) != string::npos)
            {
                /// Here is a Fx option for a floored coupon
                fxResetDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,FxResetDate).c_str()) );
                fxPayDate   = ARM_Date( atof(dealDesc.GetElem(eventIdx,CpnPayDate).c_str()) );
                if(itsCalibType == ATSFxEquivCalib)
                {
                    fxStrike = ComputeFxEquivStrike(eventIdx);
                }
                else if(itsCalibType == ATSFxProfileCalib)
                {
                    fxStrike = itsCalibDatas[strikeIdx];
                    if(strikeIdx+1 < nbStrikes)
                        ++strikeIdx;
                }
                else
                {
                    fxStrike    = atof(dealDesc.GetElem(eventIdx,FxLowStrike).c_str());

                    columnText = dealDesc.GetElem(eventIdx,FxHighCall);
                    if(columnText.find(fxCallText) != string::npos)
                        /// Here is a Fx option for a capped coupon
                        fxStrike = 0.5*( fxStrike + atof(dealDesc.GetElem(eventIdx,FxHighStrike).c_str()) );
                }

                fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
                fxOption->SetPayDate(fxPayDate);

                fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
            }
        }
    }


    ARM_StdPortfolio* fxPort = new ARM_StdPortfolio(fxOptionList);
    for(i=0;i<fxPort->size();++i)
    {
        fxPort->SetWeight(FX_DEFAULT_WEIGHT,i);
        fxPort->SetPrice((i+1)*FX_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(fxPort);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: CreateFxCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_PRCSCalculator::CreateFxCalibration(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx)
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
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeFxOptionPrice
///	Returns: nothing
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputeFxOptionPrices(ARM_CalibMethod* calibMethod, bool isFreezeWeights, bool isInitVolParam)
{
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    ARM_StdPortfolioPtr fxPortfolio = calibMethod->GetPortfolio();
    size_t nbFx = fxPortfolio->GetSize();
    size_t volParamSize=calibMethod->GetCalibParams().size();
    bool isInitVol = isInitVolParam || volParamSize == 0;

    ARM_GP_Vector initTimes(nbFx);
    ARM_GP_Vector initVols(nbFx);

    size_t i;
    ARM_Option* fxOption;
    double optTime,vol,strike,price,vega=0.0,nominal,weight,fwd, atmVol;
    bool isNotSelected;
    size_t strikeIdx=0;
    for(i=0;i<nbFx;++i)
    {
        fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(i));
	    fxOption->SetModel(fxBSModel);
        price=fxOption->ComputePrice();

        switch(itsCalibType)
        {
        case ATSFxMixedCalib : /// Choose ATM/ATS strike with minimum volatility
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

        case ATMDoubleCalib :
            /// strike = fwdFx
            fxOption->SetStrike(fxOption->GetCalcFwd());
            price=fxOption->ComputePrice();
            break;

        case ATSFxMoneynessCalib :
        case ATSFxShiftedCalib :
		case HybridBasketCalib :
            if(itsCalibType == ATSFxMoneynessCalib)
                /// strike = moneyness * fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] * fxOption->GetCalcFwd());
			else if(itsCalibType == HybridBasketCalib)
				 /// First calibration is ATSFxMoneyness=100% like
                fxOption->SetStrike(fxOption->GetCalcFwd());
            else
                /// strike = shift + fwdFx
                fxOption->SetStrike(itsCalibDatas[strikeIdx] + fxOption->GetCalcFwd());

            price=fxOption->ComputePrice();
            if(strikeIdx+1 < itsCalibDatas.size())
                ++strikeIdx;
            break;

        case ATSFxMinVolCalib : /// Choose ATS with mimimum volatility
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
///	Class  : ARM_PRCSCalculator
///	Routine: CreateHybridBasketOption
///	Returns: a portfolio
///	Action : Build a portfolio made of spread between a 
///			 variable notional Fx leg and a IR swap 
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_PRCSCalculator::CreateHybridBasketOption()
{
	/// Restore models
    double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(GetKeys()[YcBasisDomKey]));
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );
    
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
    for(eventIdx=1,flowIdx=0;eventIdx<nbEvents;++eventIdx,++flowIdx)
    {
        eventDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );

		/// Save funding datas
		fundingStarts[flowIdx]	= ARM_Date( atof(dealDesc.GetElem(eventIdx,FundingStartDate).c_str()) ).GetJulian() - asOfDate;
		fundingEnds[flowIdx]	= ARM_Date( atof(dealDesc.GetElem(eventIdx,FundingEndDate).c_str()) ).GetJulian() - asOfDate;
		fundingPeriods[flowIdx] = atof(dealDesc.GetElem(eventIdx,FundingIT).c_str());
		fundingMargins[flowIdx] = atof(dealDesc.GetElem(eventIdx,FundingSpread).c_str());

		/// Compute coupon sensitivities
        payDates[flowIdx] = ARM_Date( atof(dealDesc.GetElem(eventIdx,CpnPayDate).c_str()) );
		cpnCoef = atof(dealDesc.GetElem(eventIdx,CpnIT).c_str()) * atof(dealDesc.GetElem(eventIdx,CpnLeverage).c_str());

        columnText = dealDesc.GetElem(eventIdx,FxLowCall);
        if(columnText.find(fxCallText) != string::npos)
        {
            /// Here is a Fx option for a floored coupon
            fxResetDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,FxResetDate).c_str()) );
            fxStrike    = atof(dealDesc.GetElem(eventIdx,FxLowStrike).c_str());

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

            columnText = dealDesc.GetElem(eventIdx,FxHighCall);
            if(columnText.find(fxCallText) != string::npos)
			{
                /// Here is a Fx option for a capped coupon
				fxStrike    = atof(dealDesc.GetElem(eventIdx,FxHighStrike).c_str());
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
			zcSensis[flowIdx]	= atof(dealDesc.GetElem(eventIdx,FixedCpn).c_str());
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
	flowIdx = nbEvents-1;
	bool isRedemption;
	ARM_Date lastCpnPayDate(payDate);
    if((isRedemption=(itsRedemptionDatas.size()>0)))
    {
        eventDate = ARM_Date( atof(dealDesc.GetElem(flowIdx,EventDate).c_str()) );

        payDates[flowIdx] = ARM_Date( atof(dealDesc.GetElem(flowIdx,RedemptionPayDate).c_str()) );
		if(lastCpnPayDate < payDate)
			lastCpnPayDate = payDate;
        fxResetDate = ARM_Date( atof(dealDesc.GetElem(flowIdx,RedemptionResetDate).c_str()) );
        fxStrike    = atof(dealDesc.GetElem(nbEvents-1,RedemptionStrike).c_str());

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

        if(itsRedemptionDatas[REDEMPTION_TYPE] == MANDATORY_REDEMPTION)
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

FILE* f=NULL;
/****
f=fopen("c:\\temp\\dumpIRFXCalib.txt","a");
for(size_t ii=0;ii<fxSensis.size();++ii)
{
	double fundPrice = 0.0;
	if(ii<fundingStarts.size())
		fundPrice = basisDomCurve->DiscountPrice(fundingStarts[ii]/K_YEAR_LEN)
					+(fundingPeriods[ii]*fundingMargins[ii]-1)*
					basisDomCurve->DiscountPrice(fundingEnds[ii]/K_YEAR_LEN);
	fprintf(f,"#%3d\tFxReset=\t%10.1lf\tFxSet=\t%10.1lf\tFxSensi=\t%13.8lf\tFwdFx=\t%13.8lf\tZcFxSet=\t%13.8lf\tFxPrice=\t%13.8lf\tZcPay=\t%10.1lf\tZcSensi=\t%13.8lf\tZcCpnPay=\t%13.8lf\tCpnPrice=\t%13.8lf\tFundPrice=\t%13.8lf\n",
		ii,fxResetDates[ii].GetJulian()-asOfDate,fxSetDates[ii].GetJulian()-asOfDate,
		fxSensis[ii],fwdFxs[ii],zcFxSets[ii],fxPrices[ii],
		payDates[ii].GetJulian()-asOfDate,zcSensis[ii],zcCpnPays[ii],
		cpnPrices[ii],fundPrice);
}
fprintf(f,"The End\n\n");
fclose(f);
****/

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
    for(eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
		/// Compute a standard domestic swap curve schedule as of current event date
		/// that covers the very last payment date
        eventDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,EventDate).c_str()) );
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

/****
f=fopen("c:\\temp\\dumpIRFXCalib.txt","a");
fprintf(f,"T0=\t%10.1lf\tBasketStrike=\t%13.8lf\n",startDate.GetJulian()-asOfDate,basketStrike);
for(size_t ij=0;ij<nbDates;++ij)
{
	fprintf(f,"#%3d\tT=\t%10.1lf\tVarNotio=\t%13.8lf\tFixIT=\t%13.8lf\n",
		ij,(*legDates)[ij]-asOfDate,varLegNotios[ij],(*interestTerms)[ij]);
}
fprintf(f,"\n\n");
fclose(f);
****/

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
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeHybridBasketPrices
///	Returns: nothing
///	Action : Compute prices of hybrid baskets of the input portfolio
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputeHybridBasketPrices()
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
/****
FILE* f=fopen("c:\\temp\\dumpIRFXCalib.txt","a");
	hybridBaskets[i]->SetCurveName(modelMap[ARM_2IRFXModel::FxModel]->ModelName());
	fprintf(f,"#%3d\tMarket=\t%13.8lf\t2IRFX=\t%13.8lf\n",
		i,price,hybridBaskets[i]->Price(&*(GetPricingModel())));
	fprintf(f,"\n\n");
fclose(f);
****/
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
///	Class  : ARM_PRCSCalculator
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
void ARM_PRCSCalculator::ComputeHybridBasketCorrelations(double evalTime,const vector< ARM_VanillaIrFxSwaption* >& hybridBaskets)
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
///	Class  : ARM_PRCSCalculator
///	Routine: CreateLocalFxOption
///	Returns: nothing
///	Action : Build portfolios used for local Fx model calibration
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CreateLocalFxOption()
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
    string columnText,fxCallText("CALL");
    ARM_Date fxResetDate,fxPayDate;
    double fxStrike;
    for(size_t eventIdx=1;eventIdx<nbEvents;++eventIdx)
    {
        columnText = dealDesc.GetElem(eventIdx,FxLowCall);
        if(columnText.find(fxCallText) != string::npos)
        {
            /// Here is a Fx option for a floored coupon
            fxResetDate = ARM_Date( atof(dealDesc.GetElem(eventIdx,FxResetDate).c_str()) );
            fxPayDate   = ARM_Date( atof(dealDesc.GetElem(eventIdx,CpnPayDate).c_str()) );
            fxStrike    = atof(dealDesc.GetElem(eventIdx,FxLowStrike).c_str());

            fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_CALL,K_EUROPEAN);
            fxOption->SetPayDate(fxPayDate);

            fxFlooredOptionList.push_back(static_cast< ARM_Security* >(fxOption));

            columnText = dealDesc.GetElem(eventIdx,FxHighCall);
            if(columnText.find(fxCallText) != string::npos)
            {
                /// Here is a Fx option for a capped coupon
                fxStrike    = atof(dealDesc.GetElem(eventIdx,FxHighStrike).c_str());

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
    if(itsRedemptionDatas.size()>0)
    {
        fxPayDate   = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,RedemptionPayDate).c_str()) );
        fxResetDate = ARM_Date( atof(dealDesc.GetElem(nbEvents-1,RedemptionResetDate).c_str()) );
        if(itsRedemptionDatas[REDEMPTION_TYPE] == MANDATORY_REDEMPTION)
            /// Mandatory case : only forward is used but an ATM option is built
            /// to compute convexity adjustment through local model
            fxStrike = K_MARKET_RATE;
        else
            /// Dual option case
            fxStrike = atof(dealDesc.GetElem(nbEvents-1,RedemptionStrike).c_str());

        fxOption = new ARM_Option(forex,fxResetDate,fxStrike,K_PUT,K_EUROPEAN);
        fxOption->SetPayDate(fxPayDate);

        fxRedemptionOptionList.push_back(static_cast< ARM_Security* >(fxOption));
        delete itsRedemptionFxOptionPF;
        itsRedemptionFxOptionPF = new ARM_StdPortfolio(fxRedemptionOptionList);
    }
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeLocalFxOptionPrices
///	Returns: nothing
///	Action : Compute prices of Fx options in portfolios used for
///          local Fx model calibration
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputeLocalFxOptionPrices()
{
    ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[LocalFxModelKey]) );

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
///	Class  : ARM_PRCSCalculator
///	Routine: CalibrateLocalFxModel
///	Returns: nothing
///	Action : Calibrate local Fx Models using associated Fx option portfolios
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::CalibrateLocalFxModel()
{
    if(itsFxLocalModelFlag)
    {
        if(itsCalibType == ATMCalib)
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

//FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
//fprintf(f,"AsOf=%10.1lf\n",asOfTime);
//fclose(f);

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
                if(itsCalibType == ATMCalib)
                {
                    /// ATM vol is interpolated from bootstrapped vols but
                    // smile adjustment is interpolated from standard expiries
                    yfExpiry    = resetTime/K_YEAR_LEN;
                    atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());
                    fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());
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

//FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
//fprintf(f,"%2d\t%10.1lf\t%10.7lf\t%10.7lf\t%10.7lf\n",i,fxOption->GetExpiry().GetJulian(),fxOption->GetCalcFwd(),atmFwdFxvol+fxOption->GetCalcSmileAdjVol(),atmFwdFxvol);
//fclose(f);

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
                if(itsCalibType == ATMCalib)
                {
                    yfExpiry    = resetTime/K_YEAR_LEN;
                    atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());
                    fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());
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
          if(itsCalibType == ATMCalib)
            {
                yfExpiry    = resetTime/K_YEAR_LEN;
                atmFwdFxvol = itsATMFwdFxVol->ComputeVolatility(yfExpiry,fxOption->GetCalcFwd());
                fxOption->SetCalcVol(atmFwdFxvol+fxOption->GetCalcSmileAdjVol());
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
///	Class  : ARM_PRCSCalculator
///	Routine: ComputePricingData
///	Returns: 
///	Action : pricing function called from the addin GetPricingData()
////////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::ComputePricingData() const
{
	if (!itsHasBeenPriced)
		const_cast<ARM_PRCSCalculator*>(this)->PriceAndTimeIt();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: Price
///	Returns: a double
///	Action : price the PRDC option. An auto-calibration is done before
///          by bootstrapping the volatility on diagonal swaptions
///          on both domestic & foreign market then the spot FX volatility
////         is boostrapped on FX options portfolio
/////////////////////////////////////////////////////////////////
double ARM_PRCSCalculator::Price()
{
    ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());

    /// Calibrate stochatic models and check errors
    GetCalibMethod()->Calibrate(hybridModel);

    CheckCalibErrors();

    if(itsCalibType == ATMDoubleCalib)
    {
		/// Second step : calibration of PRD Fx options with bootstrapped ATM Fx vols
        ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(GetKeys()[FxModelKey]) );

        /// Save market Fx vols
		ARM_VolCurve* marketFxVol = static_cast<ARM_VolCurve*>(itsInputModel->GetFxVol());

        /// Compute previously bootstrapped ATM vols to give new targets
	    itsATMFwdFxVol = ComputeATMFxOptionVols();
		fxBSModel->SetVolatility(itsATMFwdFxVol);
        bool isFreezeWeights=true;
        bool isInitVolParam=true;
        ComputeFxOptionPrices(itsATMDoubleCalibMethod,isFreezeWeights,isInitVolParam);

		/// Reset the Fx vols in the model (else vol schedule will be merged with initTimes[]
		const ARM_ModelNameMap& modelMap = * (hybridModel->GetModelMap());	
		ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[ARM_2IRFXModel::FxModel]->Model() );
		fxModel.GetModelParams()->DeleteModelParam(ARM_ModelParamType::QVol);

        /// Replace previous ATM calib method by additional one
        ARM_CalibMethod* atmFxCalib = GetFxCalib();
        SetFxCalib(itsATMDoubleCalibMethod);

        // Calibrate PRD leg Fx options and with bootstrapped ATM vols 
        GetCalibMethod()->Calibrate(hybridModel);

        CheckCalibErrors();

        /// Restore market Fx vols in the market data manager (side effect !)
		fxBSModel->SetVolatility(marketFxVol);

        // Restore ATM calib method
        SetFxCalib(atmFxCalib);

    }
	else if(itsCalibType == HybridBasketCalib)
	{
		/// Second step : FX calibration with IR/FX hybrid basket options

		/// Compute IR/FX hybrid basket market prices using
		/// the previouly calibrated model
		/// Correlation may be computed at this stage
		ComputeHybridBasketPrices();

        // Calibrate FX vols on IR/FX hybrid basket options
        itsHybridBasketCalibMethod->Calibrate(hybridModel);

        CheckCalibErrors(itsHybridBasketCalibMethod,3);
	}

    /// Calibrate local Fx models
    CalibrateAndTimeIt();


    ARM_PricingModelPtr initialModel;
    if(itsMarkovianDriftSamplerFlag)
        /// Save calibrated model because of interpolation over tree schedule
        initialModel = ARM_PricingModelPtr( static_cast< ARM_PricingModel* >(hybridModel->Clone()) );

	ARM_GenPricer* genPricer = new ARM_GenPricer( &*GetGenSecurity(),hybridModel);

	ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

	double firstColumnPrice = genPricer->Price();

    string columnName;
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
        columnName = itsColumnsToPrice[i];
        GetPricingData()[columnName] = genPricer->GetPricerInfo()->GetContents(columnName).GetData("Price").GetDouble();
    }

    if(itsMarkovianDriftSamplerFlag)
        /// Restore calibrated model
        SetPricingModel(initialModel);

    itsHasBeenPriced = true;

    return firstColumnPrice;       
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: ComputeAll
///	Returns: an ARM_Vector
///	Action : call price() function and store results in a vector
/////////////////////////////////////////////////////////////////
ARM_Vector* ARM_PRCSCalculator::ComputeAll()
{
	ARM_Vector* results = new ARM_Vector(9, 0.0);
	ARM_VolCurve* originalFXVol = NULL;
	ARM_VolLInterpol* latticeATMFXVol = NULL;

	try
	{
		// Pricing
		double optionPrice = Price();

		// Get from calculator pricing results
        double undAfterFirstNoticeInLattice;
        undAfterFirstNoticeInLattice = GetPricingData().GetData("PRDCFIRSTSWAP").GetDouble() * (-1.0); 

        /// "Before" prices fully use market vols (unsmiled and smiled)
		double undPriceBeforeFirstNotice        = itsInputPowRev->GetUndPriceBeforeFirstNoticeDate();
		double FXSmileAdjustmentBeforeNotice    = itsInputPowRev->GetFXSmileAdjBeforeNotice();

        double FXSmileAdjustment;
		if( ! itsInputPowRev->IsFXSmileAdjFlag())
            FXSmileAdjustment    = itsInputPowRev->GetFXSmileAdjustment();

		// save FXVol before modify it
		originalFXVol = static_cast<ARM_VolCurve*>(itsInputModel->GetFxVol()->Clone());

        if(!itsATMFwdFxVol)
        {
            /// Not already computed then do it !
		    itsATMFwdFxVol  = ComputeATMFxOptionVols();
        }
        latticeATMFXVol = itsATMFwdFxVol;

/****
FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
fprintf(f,"\n\nTerm FX vols (2IRFX) :\n");
ARM_Matrix* vols = latticeATMFXVol->GetVolatilities();
ARM_Vector* terms = latticeATMFXVol->GetExpiryTerms();
for(size_t ii=0;ii<terms->size();++ii)
    fprintf(f,"%10.7lf\t%7.3lf\n",(*terms)[ii],vols->Elt(ii,0));
fclose(f);
****/


		itsInputModel->SetFxVol(latticeATMFXVol);

		itsInputPowRev->SetModelVariable(NULL);
		itsInputPowRev->SetModel(itsInputModel);

        /// Compute unsmiled underlying price (using ATM fwd FX vols given by
        /// bootstrapped spot FX & Zc vols)
		double undPrice             = itsInputPowRev->ComputePrice();
		double dualOptionPrice      = itsInputPowRev->GetDualValue();

        if( itsInputPowRev->IsFXSmileAdjFlag())
        {
            /// Get new smile adjustment using model ATM fwd FX vols
            FXSmileAdjustment = itsInputPowRev->GetFXSmileAdjustment();

            /// "Before" prices still inchanged else inconsistency !
            if( fabs(undPriceBeforeFirstNotice-itsInputPowRev->GetUndPriceBeforeFirstNoticeDate()) > 1e-3 ||
                fabs(FXSmileAdjustmentBeforeNotice-itsInputPowRev->GetFXSmileAdjBeforeNotice()) > 1e-2 )
            {
		        if (originalFXVol)
		        {
			        itsInputModel->SetFxVol(originalFXVol);
			        itsInputPowRev->SetModel(itsInputModel);
			        delete originalFXVol;
		        }
                ARM_THROW( ERR_SWAP_CALC_PB, ARM_USERNAME + " : inconsistency for Before underlying prices");
            }
        }

		// reset original FXVol
		itsInputModel->SetFxVol(originalFXVol);
		itsInputPowRev->SetModel(itsInputModel);

        delete originalFXVol;

        /// Get first PRD call
		double europeanCall = GetPricingData().GetData("PRDCFIRSTEUROPEAN").GetDouble();

		// Unavailable data
		double totalFXDeltaAfterFirstNotice = -9999;
		
		results->Elt(0) = undAfterFirstNoticeInLattice;
		results->Elt(1) = optionPrice;
		results->Elt(2) = europeanCall;
		results->Elt(3) = totalFXDeltaAfterFirstNotice;
		results->Elt(4) = undPrice;
		results->Elt(5) = undPriceBeforeFirstNotice;
		results->Elt(6) = FXSmileAdjustment;
		results->Elt(7) = FXSmileAdjustmentBeforeNotice;
		results->Elt(8) = dualOptionPrice;		
	}
	catch (...)
	{
		if (originalFXVol)
		{
			itsInputModel->SetFxVol(originalFXVol);
			itsInputPowRev->SetModel(itsInputModel);
			delete originalFXVol;
		}

        ARM_THROW( ERR_SWAP_CALC_PB, ARM_USERNAME + " : error in ARM_PRCSCalculator::ComputeAll()");
	}

	return results;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: UpdateModel
///	Returns: void
///	Action : update the model datas w.r.t. MarketDataManager
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::UpdateModel()
{
    /// To do...
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRCSCalculator
///	Routine: UpdateCalibration
///	Returns: void
///	Action : update the calibration datas w.r.t. MarketDataManager
///          The context is an hedge ratios computation
/////////////////////////////////////////////////////////////////
void ARM_PRCSCalculator::UpdateCalibration(bool isUpdateStrike)
{
    /// To do...
}


////////////////////////////////////////////////////
///	Class   : ARM_PRCSCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_PRCSCalculator::toString(const string& indent, const string& nextIndent) const
{
    CC_Ostringstream prdcData;

    /// PRDC Calculator specific datas dump
    prdcData << indent <<"\n\n =======> POWER REVERSE DUAL CURRENCIES CALCULATOR <====== \n\n";

    prdcData << indent << "Columns to Price :\n";
    for(size_t i=0;i<itsColumnsToPrice.size();++i)
    {
            prdcData << indent << itsColumnsToPrice[i] << "\n";
    }

    /// Additional ATM calib method
    if(itsCalibType == ATMDoubleCalib && itsATMDoubleCalibMethod)
    {
        prdcData << indent << "\nAdditional Calib Method :\n";
        prdcData << indent << itsATMDoubleCalibMethod->toString(indent,nextIndent);
    }

    /// GenCalculator general datas dump
	prdcData << "\n\nCommon Calculator Part\n" << ARM_GenCalculator::toString(indent,nextIndent) << "\n\n";

    return prdcData.str();
}

////////////////////////////////////////////////////
///	Class   : ARM_PRCSCalculator
///	Routines: ComputeATMFxOptionVols
///	Returns :
///	Action  : Compute the term fx volatility for each
///           Fx option of the PRDC after the 1st notice
////////////////////////////////////////////////////
ARM_VolLInterpol* ARM_PRCSCalculator::ComputeATMFxOptionVols()
{
    ARM_PricingStatesPtr dumStates( new ARM_PricingStates(1,1,0) );    

    ARM_2IRFXModel* hybridModel = static_cast< ARM_2IRFXModel* > (&*GetPricingModel());
    ARM_Date asOfDate(hybridModel->GetAsOfDate());
    double asOfJulian=asOfDate.GetJulian();

    bool isATMCalib = itsCalibType == ATMCalib || itsCalibType == ATMDoubleCalib;

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

        PRDCCalibType prevCalibType = itsCalibType;
        itsCalibType = ATMCalib;

        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVolParam);

        /// Calibrate stochatic models ATM
        GetCalibMethod()->Calibrate(hybridModel);

        CheckCalibErrors();

        /// Restore original calibration strikes & prices
        for(i=0;i<nbFx;++i)
            (static_cast< ARM_Option* >(fxPortfolio->GetAsset(i)))->SetStrike( prevStrike[i] );
        itsCalibType = prevCalibType;
        ComputeFxOptionPrices(GetFxCalib(),isFreezeWeights,isInitVolParam);

/****
        /// Collect time infos from generic security and get the event schedule
	    ARM_PricingAdviserPtr pricingAdviser = GetGenSecurity()->GetPricingAdviser();
	    pricingAdviser->DecidePricing(hybridModel);
	    ARM_TimeInfoPtrVector timeInfos	= pricingAdviser->GetTimeInfos(hybridModel);
	    ARM_EventTime discretisationScheme;
	    ARM_GP_Vector* eventTimes =  discretisationScheme.ModelTimesFromTimesInfo(timeInfos,*hybridModel);
	    hybridModel->GetNumMethod()->SetTimeSteps(*eventTimes);
****/
    }

    /// Get the Fx model
    const ARM_ModelNameMap& modelMap = * hybridModel->GetModelMap();
    ARM_QModel1F_Fx& fxModel = static_cast<ARM_QModel1F_Fx&>( * modelMap[ARM_2IRFXModel::FxModel]->Model() );


    /// Compute Fx options then get saved forward Fx vols
    double evalTime=0.0;
    bool isDualOption = (itsRedemptionDatas.size()>0 && itsRedemptionDatas[REDEMPTION_TYPE] == DUAL_OPTION_REDEMPTION);
    bool isDiffLastExpiry = isDualOption && (itsRedemptionDatas[REDEMPTION_RESET_DATE] != itsBoosterDatas(itsBoosterDatas.rows()-1,FX_RESET_DATE));

    double fxResetTime,settlementTime,payTime;
    int callPut=K_CALL;
    ARM_GP_Vector strikePerStates(1);
    ARM_GP_VectorPtr ATMStrike;
    size_t nbResetsBefore = itsFxResetTimeBefore.size();
    size_t resetIdx,nbResetsAfter = itsBoosterDatas.rows() + (isDualOption ? 1 : 0);
    size_t nbResets = nbResetsBefore + nbResetsAfter;
    ARM_GP_Vector fxResetYfs(nbResets);
    ARM_GP_Matrix fwdFxVols(nbResets,4,0.0);
    ARM_GP_VectorPtr price;

    /// Compute market ATM FX vols
    ARM_FXVolCurve * marketFxVols = dynamic_cast< ARM_FXVolCurve *>(itsInputModel->GetFxVol());
    if(!marketFxVols)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : market FX vol is missing");

    ARM_VolLInterpol* marketATMFxVols = marketFxVols->ComputeATMFxVol();
    const ARM_Vector& atmYfExpiries = *(marketATMFxVols->GetExpiryTerms());


    /// The first notice date is always at the first line of the booster data
    /// (because originally truncated in ComputeBoosterDatas() to begin with this date !)
    double firstNoticeTime = itsBoosterDatas(0,NOTICE_DATE)-asOfJulian;

    bool prevIsImpliedVolCalc=fxModel.GetIsImpliedVolCalc();
    fxModel.SetIsImpliedVolCalc(true); // to activate B&S implied vol calculation

    /// Add resets and vols saved before the 1st notice
    for(resetIdx=0;resetIdx < nbResetsBefore;++resetIdx)
    {
        /// Use market ATM forward FX vols
        fxResetTime = itsFxResetTimeBefore[resetIdx];
        fxResetYfs[resetIdx] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol is be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol is be used

	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        fwdFxVols(resetIdx,0) = (*ATMStrike)[0];
        fwdFxVols(resetIdx,1) = 0.01*marketATMFxVols->ComputeVolatility(fxResetYfs[resetIdx],(*ATMStrike)[0]);
    }

    /// Add resets and vols after the 1st notice
    for(resetIdx=0;resetIdx < itsBoosterDatas.rows();++resetIdx)
    {
        fxResetTime = itsBoosterDatas(resetIdx,FX_RESET_DATE)-asOfJulian;
        fxResetYfs[nbResetsBefore+resetIdx] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol is be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol is be used

        /// Use boostrapped forward FX vols
	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        price = fxModel.CallVectorial(fxModel.GetModelName(),evalTime,fxResetTime,settlementTime,*ATMStrike,callPut,payTime,dumStates);
        fwdFxVols(nbResetsBefore+resetIdx,0) = (*ATMStrike)[0];
        fwdFxVols(nbResetsBefore+resetIdx,1) = (fxModel.GetFwdFxVol())[0];
    }

    /// Get Dual redemption vol if any
    callPut=K_PUT;
    if(isDualOption)
    {
        fxResetTime = itsRedemptionDatas[REDEMPTION_RESET_DATE]-asOfJulian;
        fxResetYfs[nbResetsBefore+resetIdx] = fxResetTime/K_YEAR_LEN;
        settlementTime=fxResetTime; // not relevant because only fwd fx vol is be used
        payTime = fxResetTime;      // not relevant because only fwd fx vol is be used

        /// Use boostrapped forward FX vols
	    ATMStrike = fxModel.Forward(fxModel.GetModelName(), evalTime, fxResetTime, settlementTime, payTime, dumStates );
        strikePerStates[0]= (*ATMStrike)[0];

        price = fxModel.CallVectorial(fxModel.GetModelName(),evalTime,fxResetTime,settlementTime,strikePerStates,callPut,payTime,dumStates);
        fwdFxVols(nbResetsBefore+resetIdx,0) = (*ATMStrike)[0];
        fwdFxVols(nbResetsBefore+resetIdx,1) = (fxModel.GetFwdFxVol())[0];
    }

    /// Get FX option strikes if necessary
    ARM_GP_Vector strikes(1,0.0);
    size_t strikeIdx,nbStrikes=strikes.size();

    nbResets -= (isDualOption && !isDiffLastExpiry ? 1 : 0);
    fxResetYfs.resize(nbResets);
    ARM_Matrix* vols = new ARM_Matrix(nbResets,nbStrikes,0.0);

    /// Take care Kernel is in 100 base
    for(resetIdx=0;resetIdx<nbResetsBefore+itsBoosterDatas.rows();++resetIdx)
        vols->Elt(resetIdx,0) = 100.0*fwdFxVols(resetIdx,1);

    if(isDiffLastExpiry)
        /// Last Fx option & dual option reset at different times
        vols->Elt(resetIdx,0) = 100.0*fwdFxVols(resetIdx,1);

    fxModel.SetIsImpliedVolCalc(prevIsImpliedVolCalc);

/****
FILE* f=fopen("c:\\temp\\dumpPRDC.txt","a");
fprintf(f,"Fx option resets :\n");
fprintf(f,"1stNotice t=%6d, yf=%10.6lf\n",(int)(firstNoticeTime),firstNoticeTime/K_YEAR_LEN);
for(resetIdx=0;resetIdx<nbResets;++resetIdx)
    fprintf(f,"t=%6d, yf=%10.6lf\n",(int)(fxResetYfs[resetIdx]*K_YEAR_LEN),fxResetYfs[resetIdx]);
fprintf(f,"Fx option strikes :\n");
for(strikeIdx=0;strikeIdx<nbStrikes;++strikeIdx)
    fprintf(f,"%10.4lf\n",strikes[strikeIdx]);
fprintf(f,"\nFx option implied vol (reset,strike) :\n\n");
for(resetIdx=0;resetIdx<nbResets;++resetIdx)
{
    for(strikeIdx=0;strikeIdx+1<nbStrikes;++strikeIdx)
        fprintf(f,"%7.3lf\t",vols->Elt(resetIdx,strikeIdx));
    fprintf(f,"%7.3lf\n",vols->Elt(resetIdx,strikeIdx));
}
fclose(f);
****/


    /// Convert ARM_GP_Vector to ARM_Vector
    ARM_Vector* theResets = new ARM_Vector(nbResets);
    for(resetIdx=0;resetIdx<nbResets;++resetIdx)
        (*theResets)[resetIdx] = fxResetYfs[resetIdx];
    ARM_Vector* theStrikes = new ARM_Vector(nbStrikes);
    for(strikeIdx=0;strikeIdx<nbStrikes;++strikeIdx)
        (*theStrikes)[strikeIdx] = strikes[strikeIdx];


    /// Free memory if necessary
    if(!isATMCalib)
        delete hybridModel;


    return new ARM_VolLInterpol(asOfDate,theResets,theStrikes,vols);

}

////////////////////////////////////////////////////
///	Class   : ARM_PRCSCalculator
///	Routines: CheckCalibErrors
///	Returns :
///	Action  : Check calibration errors and throw an
///           exception if errors are too large
////////////////////////////////////////////////////
void ARM_PRCSCalculator::CheckCalibErrors(ARM_CalibMethod* theCalib, size_t theCalibIdx)
{
    /// Recall that Dom -> For -> Fx with "next" links

    size_t i,nbErrors,calibIdx=(theCalib ? theCalibIdx : 0);
    double err;
    ARM_CalibMethod* calib;
	if(theCalib)
		calib = theCalib;
	else
		calib = &*(GetCalibMethod());

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

////////////////////////////////////////////////////
///	Class   : ARM_PRCSCalculator
///	Routines: HybridBasketDump
///	Returns :
///	Action  : Dump hybrid IR/FX swaption portfolio as
///			  basket on CATU & domestic Zc
////////////////////////////////////////////////////
string ARM_PRCSCalculator::HybridBasketDump()
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

CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

