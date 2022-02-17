/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file argconvdefault.cpp
 *
 *  \brief
 *
 *	\author  R. GUILLEMOT
 *	\version 1.0
 *	\date October 2004
 */

#include "gpbase/removeidentifiedwarning.h"

/// gpcalculator
#include "gpcalculators/argconvdefault.h"
#include "gpcalculators/tarncalculator.h"
#include "gpcalculators/captioncalculator.h"
#include "gpcalculators/callablesnowballcalculator.h"
#include "gpcalculators/csocalculator.h"
#include "gpcalculators/craspreadcalculator.h"
#include "gpcalculators/prdccalculator.h"
#include "gpcalculators/localcsocalculator.h"
#include "gpcalculators/tarnfxcalculator.h"
#include "gpcalculators/enumgencalculator.h"

/// gpinfra
#include "gpinfra/gensecurity.h"
#include "gpinfra/pricingmodel.h"
#include "gpinfra/pricingadviser.h"
#include "gpcalib/calibmethod.h"
#include "gpcalib/vanillaarg.h"


CC_BEGIN_NAMESPACE( ARM )

ARGConvTable CalculatorCcyTable[] =
{	
    /// Type Name       /// number
    "DOMESTIC",			ARM_GenCalculatorCcyType::DomCurrency,
    "FOREIGN",			ARM_GenCalculatorCcyType::ForCurrency,
    "FUNDING",			ARM_GenCalculatorCcyType::FundCurrency,
    "COUPON",			ARM_GenCalculatorCcyType::CpnCurrency,
    "UNKNOWN",			ARM_GenCalculatorCcyType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_GenCalculatorCcyType( CalculatorCcyTable, "GenCalculatorCcy" );
const ARM_ArgConvReverse ARM_ArgConvReverse_GenCalculatorCcyType( CalculatorCcyTable, "GenCalculatorCcy" );

ARGConvTable MRSCalibTable[] =
{	
    /// Type Name       /// number
    "DIAG_START_FWD",   ARM_MRSCalibrationType::diagstartfwd,
    "STM_FIRST_COLUMN", ARM_MRSCalibrationType::stmfirstcolumn,
    "UNKNOWN",			ARM_MRSCalibrationType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MRSCalibType( MRSCalibTable, "MRSCalibration" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MRSCalibType( MRSCalibTable, "MRSCalibration" );

ARGConvTable MRSStrikeCalibTable[] =
{	
    /// Type Name       /// number
    "EQUIVALENT",       ARM_MRSStrikeCalibrationType::strikeEquivalent,
    "ATM",				ARM_MRSStrikeCalibrationType::strikeATM,
    "MOYENESS",			ARM_MRSStrikeCalibrationType::strikeKeepStdMoyenessCst,
    "UNKNOWN",			ARM_MRSStrikeCalibrationType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MRSStrikeCalibType( MRSStrikeCalibTable, "MRSStrikeCalibration" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MRSStrikeCalibType( MRSStrikeCalibTable, "MRSStrikeCalibration" );

ARGConvTable SigmaCalibTable[] =
{	
    /// Type Name       /// number
    "EQUIVALENT",        ARM_SigmaCalibrationType::strikeEquivalent,
    "ATM",				 ARM_SigmaCalibrationType::strikeATM,
	"DELTAPPROXI",		 ARM_SigmaCalibrationType::StrikeDeltaApproxi,
    "UNKNOWN",			 ARM_SigmaCalibrationType::Unknown,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};
const ARM_ArgConv ARM_ArgConv_SigmaCalibType( SigmaCalibTable, "SigmaStrikeCalibration" );
const ARM_ArgConvReverse ARM_ArgConvReverse_SigmaCalibType( SigmaCalibTable, "SigmaStrikeCalibration" );


ARGConvTable TARNCalibModeTable[] =
{	
    /// Type Name       /// number
    "N",				ARM_TARNCalculator::NOCalib,
	"RF",				ARM_TARNCalculator::RFStrike,
    "EX",				ARM_TARNCalculator::ExerStrike,
    "ATM",				ARM_TARNCalculator::ATMStrike,
    "BLACK",			ARM_TARNCalculator::BlackShift,
    "GAUSS",			ARM_TARNCalculator::GaussShift,
    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TARNCalibMode( TARNCalibModeTable, "TARNCalibMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TARNCalibMode( TARNCalibModeTable, "TARNCalibMode" );

////////////////////////////////////////////////////////////////////////////////
////////////////////Caption Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ARGConvTable CaptionCalibModeTable[] =
{	
    /// Type Name       /// number
    "EXSWOPT",		    ARM_CaptionCalculator::SWOPT_At_Exer,
	"ALLSWOPT",			ARM_CaptionCalculator::AllSWOPTION,
    "N",				ARM_CaptionCalculator::NOCalib,
	"Y",				ARM_CaptionCalculator::Calib,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CaptionCalibMode( CaptionCalibModeTable, "CaptionCalibMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CaptionCalibMode( CaptionCalibModeTable, "CaptionCalibMode" );

////////////////////////////////////////////////////////////////////////////////
////////////////////CSB Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ARGConvTable CSBCalibModeTable[] =
{	
    /// Type Name       /// number
    "CAP",				ARM_CallableSnowBallCalculator::CAP,
	"SWOPT",			ARM_CallableSnowBallCalculator::SWAPTION,
	"SUMOPT",			ARM_CallableSnowBallCalculator::SUMOPT,
    "NOCALIB",			ARM_CallableSnowBallCalculator::NOCALIB,
	"CAPSWOPT",			ARM_CallableSnowBallCalculator::CAPSWAPTION,
    
	/// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSBCalibMode( CSBCalibModeTable, "CSBCalibMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSBCalibMode( CSBCalibModeTable, "CSBCalibMode" );


ARGConvTable CSBControlVariableModeTable[] =
{	
    /// Type Name       /// number
    "N",				ARM_CallableSnowBallCalculator::NO,
	"SB",				ARM_CallableSnowBallCalculator::SB,
	"CSB",				ARM_CallableSnowBallCalculator::CSB,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSBControlVariableMode( CSBControlVariableModeTable, "CSBControlVariableMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSBControlVariableMode( CSBControlVariableModeTable, "CSBControlVariableMode" );

ARGConvTable CSBTriggerModeTable[] =
{	
    /// Type Name       /// number
	
    "Coupon",				ARM_CallableSnowBallCalculator::TCoupon,
	"ApproxCoupon",			ARM_CallableSnowBallCalculator::TApproxCoupon,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSBTriggerMode( CSBTriggerModeTable, "CSBTriggerMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSBTriggerMode( CSBTriggerModeTable, "CSBTriggerMode" );


////////////////////////////////////////////////////////////////////////////////
////////////////////CSO Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ARGConvTable CSOCalibModeTable[] =
{	
    /// Type Name       /// number
    "BESTFIT",			ARM_CSOCalculator::BestFitCalib,
	"VEGABESTFIT",		ARM_CSOCalculator::VegaBestFitCalib,
	"BOOTSTRAP",		ARM_CSOCalculator::BootstrapCalib,
    "NOCALIB",			ARM_CSOCalculator::NOCalib,
    "1D",				ARM_CSOCalculator::ONEDCalib,
    "2D",				ARM_CSOCalculator::TWODCalib,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOCalibMode( CSOCalibModeTable, "CSOCalibMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibMode( CSOCalibModeTable, "CSOCalibMode" );

////////////////////////////////////////////////////////////////////////////////
////////////////////CRA SO Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

ARGConvTable CRASpreadCalibrationTypeTable[] =
{	
    /// Type Name       /// number
    "DIAG",			    ARM_CRASpreadCalculator::DIAG_CALIBRATION,
	"BASKET",			ARM_CRASpreadCalculator::BASKET_CALIBRATION,
	"BASKET_SIMPLE",	ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CRASpreadCalibrationType( CRASpreadCalibrationTypeTable, "CRASpreadCalibrationType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CRASpreadCalibrationType( CRASpreadCalibrationTypeTable, "CRASpreadCalibrationType" );

ARGConvTable CRASpreadCalibStrikeTypeTable[] =
{	
    /// Type Name       /// number
    "ATM",			    ARM_CRASpreadCalculator::ATM,
	"EQUIVALENT",		ARM_CRASpreadCalculator::EQUIVALENT,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CRASpreadCalibStrikeType( CRASpreadCalibStrikeTypeTable, "CRASpreadCalibStrikeType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CRASpreadCalibStrikeType( CRASpreadCalibStrikeTypeTable, "CRASpreadCalibStrikeType" );


////////////////////////////////////////////////////////////////////////////////
///////////////////////PRDC Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ARGConvTable PRDCCalibModeTable[] =
{	
    /// Type Name       /// number
    "ATM",              ARM_PRCSCalibTypes::ATMCalib,
    "ATSFX",            ARM_PRCSCalibTypes::ATSFxCalib,
    "ATSFXMIXED",       ARM_PRCSCalibTypes::ATSFxMixedCalib,
	"ATSFXPROFILE",     ARM_PRCSCalibTypes::ATSFxProfileCalib,
	"ATSFXMONEYNESS",   ARM_PRCSCalibTypes::ATSFxMoneynessCalib,
    "ATSFXSHIFTED",     ARM_PRCSCalibTypes::ATSFxShiftedCalib,
    "ATSFXMINVOL",      ARM_PRCSCalibTypes::ATSFxMinVolCalib,
    "ATMDOUBLE",        ARM_PRCSCalibTypes::ATMDoubleCalib,
    "ATSFXEQUIV",       ARM_PRCSCalibTypes::ATSFxEquivCalib,
    "HYBRIDBASKET",     ARM_PRCSCalibTypes::HybridBasketCalib,
	"ATSBARRIER",		ARM_PRCSCalibTypes::ATSFxBarrierMoneyness,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_PRDCCalibMode( PRDCCalibModeTable, "PRDCCalibMode" );
const ARM_ArgConvReverse ARM_ArgConvReverse_PRDCCalibMode( PRDCCalibModeTable, "PRDCCalibMode" );

ARGConvTable PRCSRedemptionTypeTable[] =
{	
    /// Type Name       /// number
    "MANDATORY",            ARM_PRCSRedemptionType::mandatoryRedemption,
    "DUAL_OPTION",          ARM_PRCSRedemptionType::dualOptionRedemption,
    "STANDARD",			    ARM_PRCSRedemptionType::standard,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_PRCSRedemptionType( PRCSRedemptionTypeTable, "PRSCRedemptionType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_PRCSRedemptionType( PRCSRedemptionTypeTable, "PRSCRedemptionType" );


ARGConvTable PRCSBasisTypeTable[] =
{	
    /// Type Name       /// number
    "FLOWBYFLOW",       ARM_PRCSBasisType::flowByflow,
    "AVERAGE",          ARM_PRCSBasisType::average,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_PRCSBasisType( PRCSBasisTypeTable, "PRSCBasisType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_PRCSBasisType( PRCSBasisTypeTable, "PRSCBasisType" );

////////////////////////////////////////////////////////////////////////////////
///////////////////////CSO Convertion//////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
ARGConvTable CSOCalibModeHWM2FTable[] =
{	
    /// Type Name			/// number
    "DIAG",					ARM_LocalCSOCalculator::DIAG_CALIBRATION,
    "BASKET",				ARM_LocalCSOCalculator::BASKET_CALIBRATION,
    "DIAG_BASKET",			ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION,
	"DIAG_LONG_SPREAD",     ARM_LocalCSOCalculator::DIAG_SPREAD_LONG,
	"DIAG_SHORT_SPREAD",	ARM_LocalCSOCalculator::DIAG_SPREAD_SHORT,
    "SHORT_LONG_SPREAD",    ARM_LocalCSOCalculator::SHORT_LONG_SPREAD,
    

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOCalibModeHWM2F( CSOCalibModeHWM2FTable, "CSOCalibModeHWM2F" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibModeHWM2F( CSOCalibModeHWM2FTable, "CSOCalibModeHWM2F" );

ARGConvTable CSOCalibModeHWM1FTable[] =
{	
    /// Type Name			/// number
    "DIAG",					ARM_LocalCSOCalculator::DIAG_CALIBRATION,
    "BASKET",				ARM_LocalCSOCalculator::BASKET_CALIBRATION,
    "DIAG_BASKET",			ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOCalibModeHWM1F( CSOCalibModeHWM1FTable, "CSOCalibModeHWM1F" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOCalibModeHWM1F( CSOCalibModeHWM1FTable, "CSOCalibModeHWM1F" );

ARGConvTable CSOStrike1TypeTable[] =
{	
    /// Type Name       /// number
    "ATM",			    ARM_LocalCSOCalculator::ATM,
	"EQUIVALENT",		ARM_LocalCSOCalculator::EQUIVALENT,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOStrike1Type( CSOStrike1TypeTable, "CSOStrike1Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike1Type( CSOStrike1TypeTable, "CSOStrike1Type" );

ARGConvTable CSOStrike2TypeTable[] =
{	
    /// Type Name       /// number
    "ATM",			    ARM_LocalCSOCalculator::ATM,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOStrike2Type( CSOStrike2TypeTable, "CSOStrike2Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike2Type( CSOStrike2TypeTable, "CSOStrike2Type" );

ARGConvTable CSOStrike3TypeTable[] =
{	
    /// Type Name       /// number
    "ATM",			    ARM_LocalCSOCalculator::ATM,
	"ZERO",		        ARM_LocalCSOCalculator::ZERO,
	"CAP",			    ARM_LocalCSOCalculator::CAP,
	"FLOOR",		    ARM_LocalCSOCalculator::FLOOR,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_CSOStrike3Type( CSOStrike3TypeTable, "CSOStrike3Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_CSOStrike3Type( CSOStrike3TypeTable, "CSOStrike3Type" );


ARGConvTable TARNFXModelTypeTable[] =
{	
    /// Type Name       /// number
    "2IRFX",			ARM_TARNFXCalculator::Model2IRFX,
	"1IRFX",		    ARM_TARNFXCalculator::Model1IRFX,
	"NP1IRNFX",			ARM_TARNFXCalculator::ModelNP1IRNFX,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TARNFXModelType( TARNFXModelTypeTable, "CSOStrike3Type" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TARNFXModelType( TARNFXModelTypeTable, "CSOStrike3Type" );


ARGConvTable FXVanillaTypeTable[] =
{	
    /// Type Name       /// number
    "VANILLA",			ARM_FXVanillaType::vanilla,
    "SPREAD",			ARM_FXVanillaType::spread,
    "BASKET",			ARM_FXVanillaType::basket,
	"DIGIT",			ARM_FXVanillaType::digit,
    "PERF",				ARM_FXVanillaType::perf,
    "QUOTIENT",			ARM_FXVanillaType::quotient,
    "DIGITSPREAD",		ARM_FXVanillaType::digitspread,
	"FXBALL",			ARM_FXVanillaType::FXBall,
	"FXBALLPERF",		ARM_FXVanillaType::FXBallPerf,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_FXVanillaType( FXVanillaTypeTable, "FXVanillaType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_FXVanillaType( FXVanillaTypeTable, "FXVanillaType" );


ARGConvTable FXBasketTypeTable[] =
{	
    /// Type Name       /// number
    "MAX",				ARM_FXBasketType::max,
    "MIN",				ARM_FXBasketType::min,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_FXBasketType( FXBasketTypeTable, "FXBasketType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_FXBasketType( FXBasketTypeTable, "FXBasketType" );

ARGConvTable MixPRDNoticeTypeTable[] =
{	
    /// Type Name       /// number
    "CALL",				ARM_MixPRDNoticeType::Call,
    "KO",				ARM_MixPRDNoticeType::KO,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_MixPRDNoticeType( MixPRDNoticeTypeTable, "MixPRDNoticeType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_MixPRDNoticeType( MixPRDNoticeTypeTable, "MixPRDNoticeType" );


ARGConvTable TARNFXPayoffTable[] =
{	
    /// Type Name       /// number
    "TARNFX",			ARM_TARNFXPayoffType::TARNFX,
    "CHOOSERFX",		ARM_TARNFXPayoffType::CHOOSERFX,
	"INDIANFX",			ARM_TARNFXPayoffType::INDIANFX,
	"PRDKO",			ARM_TARNFXPayoffType::PRDKO,
	"SWITCHER",			ARM_TARNFXPayoffType::SWITCHER,

    /// very important as it tells that this is the end
    ENDOFLINE_CHAR	
};

const ARM_ArgConv ARM_ArgConv_TARNFXPayoffType( TARNFXPayoffTable, "TARNFXPayoffType" );
const ARM_ArgConvReverse ARM_ArgConvReverse_TARNFXPayoffType( TARNFXPayoffTable, "TARNFXPayoffType" );

CC_END_NAMESPACE()


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

