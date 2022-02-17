
#include "firstToBeIncluded.h"

//#include "ARM_local_gp_calculators.h"
#include "ARM_local_gp_FXVanilla_calculators.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include "ARM_local_wrapper.h"
#include "ARM_local_class.h"
#include "ARM_local_gp_genericaddin.h"

/// gpbase
#include <gpbase\autocleaner.h>
#include <gpbase\datestripcombiner.h>
#include <gpbase\curve.h>
#include <gpbase\curvetypedef.h>
#include <gpbase\interpolator.h>
#include <gpbase\stringmanip.h>
#include <gpbase\datemanip.h>
#include <gpbase\typedef.h>
#include <gpbase\warning.h>
#include <gpbase\curveconvert.h>
#include <gpbase\argconvdefault.h>
#include <gpbase\cloneutilityfunc.h>
#include <gpbase\numericconstant.h>


/// gpinfra
#include <gpinfra\mktdatamanagerrep.h>
#include <gpinfra\gensecurity.h>
#include <gpinfra\pricingmodel.h>
#include <gpinfra\dealdescription.h>
#include <gpinfra\pricingadviser.h>
#include <gpinfra\modelnrefcall.h>
#include <gpinfra\gramnode.h>
#include <gpinfra\modelparam.h>
#include <gpinfra\cstmanager.h>
#include <gpinfra\curvemodelparam.h>
#include <gpinfra\argconvdefault.h>

/// gpmodels
#include <gpmodels\HW1F.h>
#include <gpmodels\HW2F.h>
#include <gpmodels\SFRM.h>
#include <gpmodels\typedef.h>
#include <gpmodels\argconvdefault.h>
#include <gpmodels\HybridBasisFwdIR.h>
#include <gpmodels\marketirmodel.h>

/// gpcalib
#include <gpcalib\kerneltogp.h>
#include <gpcalib\calibmethod.h>
#include <gpcalib\vanillamepi.h>
#include <gpcalib\vanillaswaption.h>
#include <gpcalib\vanillapricer.h>

/// gpcalculators
#include <GP_Calculators\gpcalculators\typedef.h>
#include <GP_Calculators\gpcalculators\argconvdefault.h>
#include <GP_Calculators\gpcalculators\gencalculator.h>
#include <GP_Calculators\gpcalculators\crfcalculator.h>
#include <GP_Calculators\gpcalculators\tarncalculator.h>
#include <GP_Calculators\gpcalculators\tarncalculatorsnowball.h>
#include <GP_Calculators\gpcalculators\maturitycapcalculator.h>
#include <GP_Calculators\gpcalculators\captioncalculator.h>
#include <GP_Calculators\gpcalculators\prdccalculator.h>
#include <GP_Calculators\gpcalculators\prdkocalculator.h>
#include <GP_Calculators\gpcalculators\callablesnowballcalculator.h>
#include <GP_Calculators\gpcalculators\gencsocalculator.h>
#include <GP_Calculators\gpcalculators\csocalculator.h>
#include <GP_Calculators\gpcalculators\localcsocalculator.h>
#include <GP_Calculators\gpcalculators\bermudaswaptioncalculator.h>
#include <GP_Calculators\gpcalculators\cracalculator.h>
#include <GP_Calculators\gpcalculators\craspreadcalculator.h>
#include <GP_Calculators\gpcalculators\cralocalcalculator.h>
#include <GP_Calculators\gpcalculators\globalcapcalculator.h>
#include <GP_Calculators\gpcalculators\snowrangecalculator.h>
#include <GP_Calculators\gpcalculators\hybridirfxcalculator.h>
#include <GP_Calculators\gpcalculators\tarnfxcalculator.h>
#include <GP_Calculators\gpcalculators\tarncalculatorindian.h>
#include <GP_Calculators\gpcalculators\ccscalculator.h>
#include <GP_Calculators\gpcalculators\fxvanillacalculator.h>
#include <GP_Calculators\gpcalculators\fxracalculator.h>
#include <GP_Calculators\gpcalculators\basisconverter.h>
#include <GP_Calculators\gpcalculators\craquantocalculator.h>

/// gpnumlib
#include <gpnumlib\argconvdefault.h>

/// gphelp
#include <GP_Help\gphelp\crmcookies.h>

/// ARM Kernel
#include "portfolio.h"
#include "powrev.h"
#include "model.h"
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/globalcap.h>
#include <mod/markovtree.h>
#include <mod/calibratorfrm.h>




//#include <ARM\libarm_frometk\arm_local_parsexml_util.h>

//// using the namespace directive to access ARM object!
using ARM::ARM_NumericConstants;
using ARM::ARM_ArgConv_LgNameDayCount;
using ARM::ARM_ArgConv_LgNameFrequency;
using ARM::ARM_ArgConv_CompoundingType;
using ARM::ARM_ArgConv_CompoundingFrequency;
using ARM::ARM_ArgConv_StubRules;
using ARM::ARM_ArgConv_LgTimingMod;
using ARM::ARM_ArgConv_RcvOrPay;
using ARM::ARM_DateStrip;
using ARM::ARM_DateStripPtr;
using ARM::ARM_DateStripVector;
using ARM::ARM_DateStripCombiner;
using ARM::ARM_ArgConv_PricingModelType;
using ARM::ARM_PricingModelType;
using ARM::ARM_ModelType;
using ARM::ARM_ArgConv_GenCalculatorCcyType;
using ARM::ARM_GenCalculatorCcyType;
using ARM::ARM_CcyType;
using ARM::ARM_ArgConv_MRSCalibType;
using ARM::ARM_ArgConvReverse_MRSCalibType;
using ARM::ARM_MRSCalibrationType;
using ARM::ARM_MRSCalibType;
using ARM::ARM_ArgConv_SigmaCalibType;
using ARM::ARM_ArgConvReverse_SigmaCalibType;
using ARM::ARM_SigmaCalibrationType;
using ARM::ARM_SigmaCalibType;

using ARM::ARM_RedemptionType;
using ARM::ARM_ArgConv_PRCSRedemptionType;

using ARM::ARM_BasisType;
using ARM::ARM_ArgConv_PRCSBasisType;

using ARM::ARM_VanillaType;
using ARM::ARM_ArgConv_FXVanillaType;

using ARM::ARM_BasketType;
using ARM::ARM_ArgConv_FXBasketType;

using ARM::ARM_DigitType;
using ARM::ARM_ArgConv_DigitType;

using ARM::ARM_ArgConv_CallPut;

using ARM::ARM_PRDCCalibType;

using ARM::ARM_ArgConv_MRSStrikeCalibType;
using ARM::ARM_ArgConvReverse_MRSStrikeCalibType;
using ARM::ARM_MRSStrikeCalibrationType;
using ARM::ARM_MRSStrikeCalibType;
using ARM::ARM_ArgConv_CRASpreadCalibrationType;
using ARM::ARM_ArgConv_CRASpreadCalibStrikeType;
using ARM::ARM_ArgConv_PRDCCalibMode;

using ARM::ARM_ArgConv_CSOCalibModeHWM2F;
using ARM::ARM_ArgConvReverse_CSOCalibModeHWM2F;
using ARM::ARM_ArgConv_CSOCalibModeHWM1F;
using ARM::ARM_ArgConvReverse_CSOCalibModeHWM1F;

using ARM::ARM_ArgConv_CSOStrike1Type;
using ARM::ARM_ArgConvReverse_CSOStrike1Type;
using ARM::ARM_ArgConv_CSOStrike2Type;
using ARM::ARM_ArgConvReverse_CSOStrike2Type;
using ARM::ARM_ArgConv_CSOStrike3Type;
using ARM::ARM_ArgConvReverse_CSOStrike3Type;

using ARM::ARM_GenCalculator;
using ARM::ARM_CRFCalculator;
using ARM::ARM_TARNCalculator;
using ARM::ARM_TARNCalculatorSnowBall;
using ARM::ARM_MaturityCapCalculator;
using ARM::ARM_CaptionCalculator;
using ARM::ARM_PRDCalculator;
using ARM::ARM_PRDCCalculator;
using ARM::ARM_PRDKOCalculator;
using ARM::ARM_CallableSnowBallCalculator;
using ARM::ARM_GenCSOCalculator;
using ARM::ARM_CSOCalculator;
using ARM::ARM_LocalCSOCalculator;
using ARM::ARM_BermudaSwaptionCalculator;
using ARM::ARM_CRACalculator;
using ARM::ARM_CRALocalCalculator;
using ARM::ARM_CraQuantoCalculator;
using ARM::ARM_CRASpreadCalculator;
using ARM::ARM_GlobalCapCalculator;
using ARM::ARM_SnowRangeCalculator;
using ARM::ARM_HybridIRFXCalculator;
using ARM::ARM_TARNFXCalculator;
using ARM::ARM_TARNCalculatorIndian;
using ARM::ARM_CCSCalculator;
using ARM::ARM_FXVanillaCalculator;
using ARM::ARM_FXRACalculator;

using ARM::ARM_MarketData_ManagerRep;
using ARM::ARM_GenSecurity;
using ARM::ARM_PricingModel;
using ARM::ARM_CalibMethod;
using ARM::ARM_GenSecurityPtr;
using ARM::ARM_PricingModelPtr;
using ARM::ARM_CalibMethodPtr;
using ARM::ARM_MarketData_ManagerRepPtr;
using ARM::ARM_CRMCookies;
using ARM::ARM_Curve;
using ARM::ARM_StepUpLeftOpenCstExtrapolDble;
using ARM::ARM_FlatCurve;
using ARM::ARM_HullWhite1F;
using ARM::ARM_HullWhite2F;
using ARM::ARM_SFRM;
using ARM::ARM_HybridBasisFwdIR;
using ARM::ARM_GP_Vector;
using ARM::ARM_GP_StrVector;
using ARM::ARM_StringVector;
using ARM::ARM_GP_Matrix;
using ARM::stringGetUpper;
using ARM::ARM_ObjectVector;
using ARM::ARM_BoolVector;
using ARM::ARM_ModelParam;
using ARM::ARM_VanillaMepi;
using ARM::ConvertXLDateToJulian;
using ARM::ARM_Warning;
using ARM::ARM_StdPortfolioPtr;
using ARM::ARM_CstManager;
using ARM::ARM_AutoCleaner;
using ARM::ARM_ArgConv_CSBTriggerMode;
using ARM::ARM_CurveModelParam;
using ARM::CurveToRefValue;
using ARM::ARM_MarketIRModel;
using ARM::ARM_ArgConv_VnsPricingMethod;
using ARM::ARM_VanillaSwaptionArg;
using ARM::ARM_BasisConverter;
using ARM::ARM_ZeroCurvePtr;
using ARM::CreateClonedPtr;
using ARM::ARM_ArgConv_Timing;
using ARM::ARM_ArgConv_IntRule;

using ARM::ARM_ArgConv_BaseGenAlgoType;
using ARM::ARM_ArgConv_TransformAlgoType;
using ARM::ARM_ArgConv_YesNo;
using ARM::ARM_ArgConv_TARNFXModelType;
using ARM::ARM_ArgConv_Numeraire;
using ARM::ARM_ArgConv_MMCorrelType;
using ARM::ARM_ArgConvReverse_LgNameDayCount;
using ARM::ARM_ArgConv_TARNFXPayoffType;
using ARM::ARM_FXTARNPayoffType;

//***********************************************************************
//**********  FXVanilla Calculator Creating  ****************************
//***********************************************************************

long ARM_FXVanillaCalculator_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_FXVanillaCalculator* calculator = NULL;

	ARM_Curve*	nominalCv	= NULL;
	ARM_Curve*	strikeCv	= NULL;	

	ARM_Curve*	leverageCv	= NULL;
	ARM_Curve*	alphaCv		= NULL;
	ARM_Curve*	betaCv		= NULL;
	ARM_Curve*	strike2Cv	= NULL;
	ARM_Curve*	couponMinCv	= new ARM_FlatCurve(0.0);
	ARM_Curve*	couponMaxCv	= new ARM_FlatCurve(ARM_NumericConstants::ARM_INFINITY);

	ARM_DateStrip* dateStrip = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXVanilla Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];
		
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);

		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &dateStrip, genericParams->GetParamValue("DateStrip").GetObjectId() ," DateStrip", result ) ) return ARM_KO;

		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);
		int expiryGap  = genericParams->GetParamValue("ExpiryGap").GetDouble();
		int setlmentGap = genericParams->GetParamValue("SetlmentGap").GetDouble();
		int frequency  = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("Frequency").GetString());
		int dayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("DayCount").GetString());
		string resetCal = genericParams->GetParamValue("ResetCal").GetString();
		string payCal = genericParams->GetParamValue("PayCal").GetString();
		
		string fx1Name = genericParams->GetParamValue("FX1name").GetString();
		string fx2Name = genericParams->GetParamValue("FX2name").GetString();
		ARM_Currency	payCcy(genericParams->GetParamValue("PayCcy").GetString().c_str());
		int callPut    = ARM_ArgConv_CallPut.GetNumber(genericParams->GetParamValue("CallPut").GetString());

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &nominalCv, genericParams->GetParamValue("Nominal").GetObjectId() ," Nominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &strikeCv, genericParams->GetParamValue("Strike").GetObjectId() ,"Strike", result ) ) return ARM_KO;

		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &alphaCv, genericParams->GetParamValue("Alpha").GetObjectId() ,"Alpha", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &betaCv, genericParams->GetParamValue("Beta").GetObjectId() ," Beta", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &strike2Cv, genericParams->GetParamValue("Strike2").GetObjectId() ,"Strike2", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &leverageCv, genericParams->GetParamValue("Leverage").GetObjectId() ,"Leverage", result ) ) return ARM_KO;
		
		// Coupon Min/Max
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &couponMinCv, genericParams->GetParamValue("CouponMin").GetObjectId() ,"Coupon Min", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &couponMaxCv, genericParams->GetParamValue("CouponMax").GetObjectId() ,"Coupon Max", result ) ) return ARM_KO;
		
		int callPut2    = ARM_ArgConv_CallPut.GetNumber(genericParams->GetParamValue("CallPut2").GetString());
		ARM_BasketType minMax  = (ARM_BasketType)ARM_ArgConv_FXBasketType.GetNumber(genericParams->GetParamValue("MinMax").GetString());
		ARM_DigitType digitType  = (ARM_DigitType)ARM_ArgConv_DigitType.GetNumber(genericParams->GetParamValue("DigitType").GetString());
		double epsilon  = genericParams->GetParamValue("Epsilon").GetDouble();		
		ARM_VanillaType vanillaType  = (ARM_VanillaType)ARM_ArgConv_FXVanillaType.GetNumber(genericParams->GetParamValue("OptionType").GetString());

        calculator = new ARM_FXVanillaCalculator(	asOfDate,
													*dateStrip,
													startDate,
													endDate,
													payCcy,
													expiryGap,
													setlmentGap,
													frequency,
													dayCount,
													resetCal,
													payCal,
													fx1Name,
													fx2Name,
													leverageCv ? *leverageCv: ARM_FlatCurve(1.0),
													*nominalCv,
													*strikeCv,
													callPut,
													vanillaType,
													alphaCv ? *alphaCv: ARM_FlatCurve(1.0),
													betaCv ? *betaCv: ARM_FlatCurve(1.0),
													digitType,
													epsilon,
													strike2Cv ? *strike2Cv: *strikeCv,
													callPut2,
													minMax,
													couponMinCv ? *couponMinCv : ARM_FlatCurve(0.0),
													couponMaxCv ? *couponMaxCv : ARM_FlatCurve(ARM_NumericConstants::ARM_INFINITY));


		// assign object
		return  !assignObject( calculator, result, objId ) ? ARM_KO: ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	calculator;

		x.DebugPrint();
		ARM_RESULT();
	}
}

//***********************************************************************
//**********  FXVanilla Calculator initializing  ************************
//***********************************************************************

long ARM_FXVanillaCalculator_InitFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_FXVanillaCalculator*	calculator = NULL;
	ARM_FXVanillaCalculator*	newCalculator = NULL;

	ARM_MarketData_ManagerRep*	mktDataManager = NULL;

	string message("Done! ");

	ARM_Curve*	Kmin1Cv = NULL;
	ARM_Curve*	Kmax1Cv = NULL;
	ARM_Curve*	Kmin2Cv = NULL;
	ARM_Curve*	Kmax2Cv = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXVanilla Calculator Init" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calculator, genericParams->GetParamValue("FXCalculatorId").GetObjectId() ," FXVanilla Calculator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, genericParams->GetParamValue("MktDataManager").GetObjectId() ," Mkt Data Manager", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &Kmin1Cv, genericParams->GetParamValue("Kmin1").GetObjectId() ,"Kmin1", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &Kmax1Cv, genericParams->GetParamValue("Kmax1").GetObjectId() ,"Kmax1", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &Kmin2Cv, genericParams->GetParamValue("Kmin2").GetObjectId() ," Kmin2", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &Kmax2Cv, genericParams->GetParamValue("Kmax2").GetObjectId() ,"Kmax2", result ) ) return ARM_KO;

		int Nleft1  = genericParams->GetParamValue("Nleft1").GetDouble();
		int Ncenter1  = genericParams->GetParamValue("Ncenter1").GetDouble();
		int Nright1  = genericParams->GetParamValue("Nright1").GetDouble();
		int Nleft2  = genericParams->GetParamValue("Nleft2").GetDouble();
		int Ncenter2  = genericParams->GetParamValue("Ncenter2").GetDouble();
		int Nright2  = genericParams->GetParamValue("Nright2").GetDouble();
	
		newCalculator = (ARM_FXVanillaCalculator*)calculator->Clone();

		newCalculator->Init(*mktDataManager,
							*Kmin1Cv,
							*Kmax1Cv,
							*Kmin2Cv,
							*Kmax2Cv,
							Nleft1,
							Ncenter1,
							Nright1,
							Nleft2,
							Ncenter2,
							Nright2);

		// assign object
		return  ( !assignObject( newCalculator, result, objId ) )  ?  ARM_KO : ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	newCalculator;

		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}

//***********************************************************************
//**********  FXRA Calculator Creating  ****************************
//***********************************************************************

long ARM_FXRACalculator_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_FXRACalculator* calculator = NULL;

	ARM_Curve*	nominalCv = NULL;
	ARM_Curve*	strikeCv = NULL;
	ARM_Curve*	alphaCv = NULL;
	ARM_Curve*	betaCv = NULL;
	ARM_Curve*	strike2Cv = new ARM_Curve();
	ARM_Curve*	leverageCv = NULL;
	
	ARM_Curve*	fxDownBarrierCv = NULL;
	ARM_Curve*	fxUpBarrierCv = NULL;
	ARM_Curve*	irDownBarrierCv = NULL;
	ARM_Curve*	irUpBarrierCv = NULL;

	ARM_DateStrip* dateStrip;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXVanilla Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];
		
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);

		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &dateStrip, genericParams->GetParamValue("DateStrip").GetObjectId() ," DateStrip", result ) ) return ARM_KO;

		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);
		int expiryGap  = genericParams->GetParamValue("ExpiryGap").GetDouble();
		int setlmentGap = genericParams->GetParamValue("SetlmentGap").GetDouble();
		int paymentGap = genericParams->GetParamValue("PaymentGap").GetDouble();
		int frequency  = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("Frequency").GetString());
		int dayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("DayCount").GetString());
		string resetCal = genericParams->GetParamValue("ResetCal").GetString();
		string payCal = genericParams->GetParamValue("PayCal").GetString();

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &nominalCv, genericParams->GetParamValue("Nominal").GetObjectId() ," Nominal", result ) ) return ARM_KO;
		string fx1Name = genericParams->GetParamValue("FX1name").GetString();
		string fx2Name = genericParams->GetParamValue("FX2name").GetString();
		ARM_Currency	payCcy(genericParams->GetParamValue("PayCcy").GetString().c_str());

		int callPut    = ARM_ArgConv_CallPut.GetNumber(genericParams->GetParamValue("CallPut").GetString());
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &strikeCv, genericParams->GetParamValue("Strike").GetObjectId() ,"Strike", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &alphaCv, genericParams->GetParamValue("Alpha").GetObjectId() ,"Alpha", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &betaCv, genericParams->GetParamValue("Beta").GetObjectId() ," Beta", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &strike2Cv, genericParams->GetParamValue("Strike2").GetObjectId() ,"Strike2", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &leverageCv, genericParams->GetParamValue("Leverage").GetObjectId() ,"Leverage", result ) ) return ARM_KO;

		int callPut2    = ARM_ArgConv_CallPut.GetNumber(genericParams->GetParamValue("CallPut2").GetString());
		ARM_BasketType minMax  = (ARM_BasketType)ARM_ArgConv_FXBasketType.GetNumber(genericParams->GetParamValue("MinMax").GetString());
		ARM_DigitType digitType  = (ARM_DigitType)ARM_ArgConv_DigitType.GetNumber(genericParams->GetParamValue("DigitType").GetString());
		double epsilon  = genericParams->GetParamValue("Epsilon").GetDouble();
		
		ARM_VanillaType vanillaType  = (ARM_VanillaType)ARM_ArgConv_FXVanillaType.GetNumber(genericParams->GetParamValue("OptionType").GetString());

		int stubType	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());
		int resetTiming	= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("ResetTiming").GetString());


		int fixingfrequency  = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("FixingFrequency").GetString());
		string payIdx = genericParams->GetParamValue("PAYidx").GetString();
		double payIdxSpread = genericParams->GetParamValue("PAYidxSpread").GetDouble();
		string payIdxIT = genericParams->GetParamValue("PAYidxIT").GetString();
		string irIdx = genericParams->GetParamValue("IRidx").GetString();
		string irIdxIT = genericParams->GetParamValue("IRidxIT").GetString();

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fxDownBarrierCv, genericParams->GetParamValue("FxDownBarrier").GetObjectId() ,"FxDownBarrier", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fxUpBarrierCv, genericParams->GetParamValue("FxUpBarrier").GetObjectId() ,"FxUpBarrier", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &irDownBarrierCv, genericParams->GetParamValue("IrDownBarrier").GetObjectId() ,"IrDownBarrier", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &irUpBarrierCv, genericParams->GetParamValue("IrUpBarrier").GetObjectId() ,"IrUpBarrier", result ) ) return ARM_KO;
		
		calculator = new ARM_FXRACalculator(	asOfDate,
												*dateStrip,
												startDate,
												endDate,
												expiryGap,
												setlmentGap,
												paymentGap,
												frequency,
												dayCount,
												resetCal,
												payCal,
												fx1Name,
												fx2Name,
												payCcy,
												*nominalCv,
												callPut,
												*strikeCv,
												*alphaCv,
												*betaCv,
												*strike2Cv,
												*leverageCv,
												callPut2,
												minMax,
												digitType,
												epsilon,
												vanillaType,
												intRule,
												stubType,
												resetTiming,
												fixingfrequency,
												payIdx,
												payIdxSpread,
												payIdxIT,
												irIdx,
												irIdxIT,
												*fxDownBarrierCv,
												*fxUpBarrierCv,
												*irDownBarrierCv,
												*irUpBarrierCv);
		// assign object
		return  !assignObject( calculator, result, objId ) ? ARM_KO: ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	calculator;

		x.DebugPrint();
		ARM_RESULT();
	}
}

//***********************************************************************
//**********  FXRA Calculator initializing  ************************
//***********************************************************************

long ARM_FXRACalculator_InitFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_FXRACalculator*	calculator = NULL;
	ARM_FXRACalculator*	newCalculator = NULL;

	ARM_MarketData_ManagerRep*	mktDataManager = NULL;

	string message("Done! ");


	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXVanilla Calculator Init" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calculator, genericParams->GetParamValue("FXCalculatorId").GetObjectId() ," FXVanilla Calculator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, genericParams->GetParamValue("MktDataManager").GetObjectId() ," Mkt Data Manager", result ) ) return ARM_KO;

		int nbPoints1  = genericParams->GetParamValue("NbPoints1").GetDouble();
		int nbPoints2  = genericParams->GetParamValue("NbPoints2").GetDouble();
	
		newCalculator = (ARM_FXRACalculator*)calculator->Clone();

		newCalculator->Init(*mktDataManager,
							nbPoints1,
							nbPoints2);

		// assign object
		return  ( !assignObject( newCalculator, result, objId ) )  ?  ARM_KO : ARM_OK; 
	}

	catch(Exception& x)
	{
		delete	newCalculator;

		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_OK;
}


/****************************************************************
	QUANTO
*****************************************************************/

long ARMLOCAL_CRAQuantoCalculator_Create(
		const ARM_Currency& ccyDom,
		const ARM_Currency& ccyFor,
		const double& startDate,
		const double& endDate,
		const int& payReceive,
		const int& callFreq,
		const int& callNotice,
		const string& callCal,
		const int& fundFreq,
		const int& fundDayCount,
		const int& cpnDayCount,
		const int& cpnPayFreq,
		const string& cpnResetCal,
		const string& cpnPayCal,
		const int& cpnResetFreq,
		const int& cpnResetTiming,
		const int& refIndex,
		const int& payIndex,
		const int& payIndexResetTiming,
		const double& notional,
		const long& notionalId,
		const long& callFeesId,
		const double& fundSpread,
		const long& fundSpreadId,
		const double& boostedFix,
		const long& boostedFixId,
		const double& bDown,
		const long& bDownId,
		const double& bUp,
		const long& bUpId,
		const vector<string>& pricingFlags,
		ARM_result&	result,
        long objId)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CraQuantoCalculator* craQuantoCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* callFeesProfile = NULL;
    bool isCallFees=false;
	
	ARM_ReferenceValue* fundSpreadProfile = NULL;
    bool isFundSpread=false;

	ARM_ReferenceValue* boostedFixProfile = NULL;
    bool isBoostedFix=false;
	
	ARM_ReferenceValue* bDownProfile = NULL;
    bool isbDown=false;
	
	ARM_ReferenceValue* bUpProfile = NULL;
    bool isbUp=false;

	try
	{
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable Range Accrual Quanto Calculator" );

		char myStartDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
	
		//Notional Curve
        if (notionalId != ARM_NULL_OBJECT)
        {
		    notionalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId));

		    if (!notionalProfile)
		    {
			    result.setMsg ("ARM_ERR: notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            notionalProfile = new ARM_ReferenceValue(notional);
            isNotional=true;
        }

		//Call Fees Curve 
		callFeesProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(callFeesId));
		if (!callFeesProfile)
		{
		    result.setMsg ("ARM_ERR: call fees are not of a good type");
		    return ARM_KO;
		}

		//Fund Spread Curve
        if (fundSpreadId != ARM_NULL_OBJECT)
        {
		    fundSpreadProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: fund spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundSpreadProfile = new ARM_ReferenceValue(fundSpread);
            isFundSpread = true;
        }

		//BoostedFix Curve
        if (boostedFixId != ARM_NULL_OBJECT)
        {
		    boostedFixProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(boostedFixId));

		    if (!boostedFixProfile)
		    {
			    result.setMsg ("ARM_ERR: boosted fix is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            boostedFixProfile = new ARM_ReferenceValue(boostedFix);
            isBoostedFix = true;
        }

		//Barrier Down Curve
        if (bDownId != ARM_NULL_OBJECT)
        {
		    bDownProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bDownId));

		    if (!bDownProfile)
		    {
				result.setMsg ("ARM_ERR: barrier down is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bDownProfile = new ARM_ReferenceValue(bDown);
            isbDown = true;
        }
 
		// Barrier Up Curve
        if (bUpId != ARM_NULL_OBJECT)
        {
		    bUpProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bUpId));

		    if (!bUpProfile)
		    {
				result.setMsg ("ARM_ERR: barrier up is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bUpProfile = new ARM_ReferenceValue(bUp);
            isbUp = true;
        }

        if (pricingFlags.size()>ARM::ARM_CraQuantoCalculator::NbProductsToPrice)
		{
			result.setMsg ("ARM_ERR: too much products to price");
			return ARM_KO;
		}
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(ARM::ARM_CraQuantoCalculator::NbProductsToPrice, false);

		for (int i = 0; i < pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

		// Create the calculator
        craQuantoCalculator = new ARM_CraQuantoCalculator(
												ccyDom,
												ccyFor,
												ARM_Date(myStartDate),
												ARM_Date(myEndDate),
												payReceive,
												*notionalProfile,
												callFreq,
												callNotice,
												callCal,
												*callFeesProfile,
												fundFreq,
												fundDayCount,
												*fundSpreadProfile,
												cpnDayCount,
												cpnPayFreq,
												cpnResetCal,
												cpnPayCal,	
												payIndex,
												payIndexResetTiming,
												*boostedFixProfile,
												*bDownProfile,
												*bUpProfile,
												cpnResetFreq,
												cpnResetTiming,
												refIndex,
												productsToPrice);
		 		
        // Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		if (isCallFees)
			delete callFeesProfile;
		callFeesProfile = NULL;
	
		if (isFundSpread)
			delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		if (isBoostedFix)
			delete boostedFixProfile;
		boostedFixProfile = NULL;

		if (isbDown)
			delete bDownProfile;
		bDownProfile = NULL;
	
		if (isbUp)
			delete bUpProfile;
		bUpProfile = NULL;
	
		// assign object
		if ( !assignObject( craQuantoCalculator, result, objId ) )
		{
			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}
	}
	
	catch(Exception& x)
	{
		if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		if (isCallFees)
			delete callFeesProfile;
		callFeesProfile = NULL;
	
		if (isFundSpread)
			delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		if (isBoostedFix)
			delete boostedFixProfile;
		boostedFixProfile = NULL;

		if (isbDown)
			delete bDownProfile;
		bDownProfile = NULL;
	
		if (isbUp)
			delete bUpProfile;
		bUpProfile = NULL;
	
		x.DebugPrint();
		ARM_RESULT();
	}
}

//===========================================================//
// FXVanillaCalculator										 //
//===========================================================//
long ARMLOCAL_FXVanillaCalculator_CreateFromSecurity(
		long FxSpreadId,
        string basketTypeStr,
		string digitTypeStr,
		string vanillaTypeStr,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
	{
		return ARM_KO;
	}

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


	/// Calculator parameters
	ARM_BasketType	minMax;
	ARM_DigitType	digitType;
	ARM_VanillaType	vanillaType;

	ARM_FxSpreadStripOption* FxStrip	 = NULL;
	ARM_FXVanillaCalculator* vanillaCalc = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "FXVanilla Calculator create from security" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &FxStrip, FxSpreadId ,"Fx spread strip option ", result ) ) return ARM_KO;


		minMax		= static_cast<ARM_BasketType> (ARM_ArgConv_FXBasketType.GetNumber(stringGetUpper(basketTypeStr)));
		digitType	= static_cast<ARM_DigitType>  (ARM_ArgConv_DigitType.GetNumber(stringGetUpper(digitTypeStr)));	
		vanillaType	= static_cast<ARM_VanillaType>(ARM_ArgConv_FXVanillaType.GetNumber(stringGetUpper(vanillaTypeStr)));


        /// Create the CRF calculator
        vanillaCalc = new ARM_FXVanillaCalculator((*FxStrip).GetAsOfDate(),
								*(FxStrip->GetSchedule()),
								(*FxStrip).GetStartDate(),
								(*FxStrip).GetExpiryDate(),
								(*FxStrip).GetPaymentCcy(),
								(*FxStrip).GetExpiryGap(),
								(*FxStrip).GetSetlmentGap(),
								(*FxStrip).GetFrequency(),
								(*FxStrip).GetDayCount(), 
								string((*FxStrip).GetSchedule()->GetResetCalendar()),
								string((*FxStrip).GetSchedule()->GetPayCalendar()),
								(*FxStrip).GetFx1Name(),
								(*FxStrip).GetFx2Name(),
								(*FxStrip).GetLeverage(),
								(*FxStrip).GetNotional(),
								(*FxStrip).GetStrikes(),
								(*FxStrip).GetOptionType(),
								vanillaType,
								(*FxStrip).GetAlpha(),
								(*FxStrip).GetBeta(),
								digitType,
								(*FxStrip).GetEpsilon());

		/// assign object
		return (assignObject( vanillaCalc, result, objId ) )? ARM_OK: ARM_KO;
	}
	
	catch(Exception& x)
	{
		delete vanillaCalc;
		x.DebugPrint();
		ARM_RESULT();
	}
}
