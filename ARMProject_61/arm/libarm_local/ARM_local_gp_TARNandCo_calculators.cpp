/*
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * $Log: ARM_local_gp_calculators.cpp,v $
 * Revision 1.1  2003/13/07 15:08:43  ebenhamou
 * Initial version
 *
 */


#include "firstToBeIncluded.h"

#include "ARM_local_gp_TARNandCo_calculators.h"
#include "ARM_local_gp_calculators_tools.h"
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

long ARM_TARNFXCalculator_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_TARNFXCalculator*	tarnFXCalculator = NULL;

	ARM_Curve*	cpnNominal = NULL;
	ARM_Curve*	domCoupon = NULL;
	ARM_Curve*	forCoupon = NULL;
	ARM_Curve*	initialFX = NULL;
	ARM_Curve*	minCoupon = NULL;
	ARM_Curve*	maxCoupon = NULL;
	ARM_Curve*	fundingNominal = NULL;
	ARM_Curve*	fundingSpread = NULL;
	ARM_Curve*	fees = NULL;
	ARM_Curve*  targetCv = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];

		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);

		ARM_Currency	domCurrency(genericParams->GetParamValue("CpnCcy").GetString().c_str());
		
		vector<ARM_Currency> forCurrencies;
		if (genericParams->GetParamValue("FgnCcy").GetType() == GA_OBJECT)
		{
			ARM_GP_StrVector* vector = dynamic_cast<ARM_GP_StrVector*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FgnCcy").GetObjectId()));
			int i=0;
			for (i = 0; i < vector->size(); i++)
			{
				string currency = vector->Elt(i);
				forCurrencies.push_back(ARM_Currency(currency.c_str()));
			}
		}
		else
			forCurrencies.push_back(ARM_Currency(genericParams->GetParamValue("FgnCcy").GetString().c_str()));

		ARM_Currency	fundCurrency(genericParams->GetParamValue("FundCcy").GetString().c_str());
		
		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(genericParams->GetParamValue("PayRec").GetString());
		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("CpnFreq").GetString());
		int cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("CpnDayCount").GetString());
		int stubRule	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString().c_str());
		int timing		= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("Timing").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());
		

		cpnNominal = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("CpnNominal").GetObjectId()));
		if (!cpnNominal)
		{
			result.setMsg ("ARM_ERR: coupon nominal should be a curve");
			return ARM_KO;
		}

		vector<ARM_Curve> domCoupons;
		ARM_Object* object = LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("DomCpn").GetObjectId());
		if (object->GetName() == ARM_GP_VECTOR)
		{
			ARM_GP_StrVector* vector = dynamic_cast<ARM_GP_StrVector*>(object);
			int i=0;
			for (i = 0; i < vector->size(); i++)
			{
				long objectId = LocalGetNumObjectId(vector->Elt(i).c_str());
				ARM_Curve* domCoupon = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(objectId));
				if (!domCoupon)
				{
					result.setMsg ("ARM_ERR: domestic coupon should be a curve");
					return ARM_KO;
				}
				domCoupons.push_back(*domCoupon);
			}
		}
		else
		{
			domCoupon = dynamic_cast<ARM_Curve*>(object);
			if (!domCoupon)
			{
				result.setMsg ("ARM_ERR: domestic coupon should be a curve");
				return ARM_KO;
			}
			domCoupons.push_back(*domCoupon);
		}

		vector<ARM_Curve> forCoupons;
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FgnCpn").GetObjectId());
		if (object->GetName() == ARM_GP_VECTOR)
		{
			ARM_GP_StrVector* vector = dynamic_cast<ARM_GP_StrVector*>(object);
			int i=0;
			for (i = 0; i < vector->size(); i++)
			{
				long objectId = LocalGetNumObjectId(vector->Elt(i).c_str());
				ARM_Curve* forCoupon = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(objectId));
				if (!forCoupon)
				{
					result.setMsg ("ARM_ERR: foreign coupon should be a curve");
					return ARM_KO;
				}
				forCoupons.push_back(*forCoupon);
			}
		}
		else
		{
			forCoupon = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FgnCpn").GetObjectId()));
			if (!forCoupon)
			{
				result.setMsg ("ARM_ERR: foreign coupon should be a curve");
				return ARM_KO;
			}
			forCoupons.push_back(*forCoupon);
		}

		minCoupon = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MinCpn").GetObjectId()));
		if (!minCoupon)
		{
			result.setMsg ("ARM_ERR: min coupon should be a curve");
			return ARM_KO;
		}
		maxCoupon = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MaxCpn").GetObjectId()));
		if (!maxCoupon)
		{
			result.setMsg ("ARM_ERR: max coupon should be a curve");
			return ARM_KO;
		}

		vector<ARM_Curve> initialFXs;
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("InitFX").GetObjectId());
		if (object->GetName() == ARM_GP_VECTOR)
		{
			ARM_GP_StrVector* vector = dynamic_cast<ARM_GP_StrVector*>(object);
			int i=0;
			for (i = 0; i < vector->size(); i++)
			{
				long objectId = LocalGetNumObjectId(vector->Elt(i).c_str());
				ARM_Curve* initialFX = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(objectId));
				if (!initialFX)
				{
					result.setMsg ("ARM_ERR: foreign coupon should be a curve");
					return ARM_KO;
				}
				initialFXs.push_back(*initialFX);
			}
		}
		else
		{
			initialFX = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("InitFX").GetObjectId()));
			if (!initialFX)
			{
				result.setMsg ("ARM_ERR: initial FX should be a curve");
				return ARM_KO;
			}
			initialFXs.push_back(*initialFX);
		}

		int fundFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("FundFreq").GetString());
		int fundDayCount= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("FundDayCount").GetString());


		fundingNominal = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FundNominal").GetObjectId()));
		if (!fundingNominal)
		{
			result.setMsg ("ARM_ERR: funding nominal should be a curve");
			return ARM_KO;
		}

		fundingSpread = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FundSpread").GetObjectId()));
		if (!fundingSpread)
		{
			result.setMsg ("ARM_ERR: funding spread should be a curve");
			return ARM_KO;
		}

		ARM_RedemptionType	redemptionType = (ARM_RedemptionType)ARM_ArgConv_PRCSRedemptionType.GetNumber(genericParams->GetParamValue("ReDeemType").GetString());

		ARM_GP_Vector redemptionStrikes;
		if (genericParams->GetParamValue("RedeemStrike").GetType() == GA_OBJECT)
		{
			ARM_GP_Vector* vector = dynamic_cast<ARM_GP_Vector*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("RedeemStrike").GetObjectId()));
			int i=0;
			for (i = 0; i < vector->size(); i++)
				redemptionStrikes.push_back(vector->Elt(i));
		}
		else
			redemptionStrikes.push_back(genericParams->GetParamValue("RedeemStrike").GetDouble());

		fees = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("Fees").GetObjectId()));
		if (!fees)
		{
			result.setMsg ("ARM_ERR: fees should be a curve");
			return ARM_KO;
		}
			
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &targetCv, genericParams->GetParamValue("Target").GetObjectId() ," target", result ) ) return ARM_KO;


		int intermPrices = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("IntermPrices").GetString());

		ARM_StringVector columnsToPrice;
		ARM_GP_StrVector* productsToPrice = dynamic_cast<ARM_GP_StrVector*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("ProductsToPrice").GetObjectId()));;
		if (!productsToPrice)
			columnsToPrice.push_back("TARN");
		else
			columnsToPrice = ARM_StringVector(productsToPrice->GetValues());

		ARM_FXTARNPayoffType payoffType = (ARM_FXTARNPayoffType) ARM_ArgConv_TARNFXPayoffType.GetNumber(genericParams->GetParamValue("PayoffType").GetString());

		ARM_FixingSched* pastFixings = NULL;

		if (genericParams->GetParamValue("FixingSched").GetType() == GA_OBJECT)
			pastFixings = dynamic_cast<ARM_FixingSched*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("FixingSched").GetObjectId()));

        tarnFXCalculator = new ARM_TARNFXCalculator(asOfDate,
													startDate,
													endDate,
													domCurrency,
													forCurrencies,
													fundCurrency,
													payRec,
													cpnDayCount,
													cpnFreq,
													genericParams->GetParamValue("CpnResetGap").GetDouble(),
													genericParams->GetParamValue("ResetCal").GetString(),
													genericParams->GetParamValue("PayCal").GetString(),
													stubRule,
													timing,
													intRule,
													*cpnNominal,
													domCoupons,
													forCoupons,
													*minCoupon,
													*maxCoupon,
													initialFXs,
													fundFreq,
													fundDayCount,
													*fundingNominal,
													*fundingSpread,
													*targetCv,
													genericParams->GetParamValue("FXChoice").GetString(),
													redemptionType,
													genericParams->GetParamValue("RedeemGap").GetDouble(),
													redemptionStrikes,
													*fees,
													intermPrices,
													columnsToPrice,
													payoffType,
													pastFixings);

		// assign object
		if ( !assignObject( tarnFXCalculator, result, objId ) )
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
		delete	tarnFXCalculator;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARM_TARNFXCalculator_InitFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_TARNFXCalculator*	newTarnFXCalculator = NULL;

	ARM_Curve*	cpnNominal = NULL;
	ARM_Curve*	domCoupon = NULL;
	ARM_Curve*	forCoupon = NULL;
	ARM_Curve*	initialFX = NULL;
	ARM_Curve*	minCoupon = NULL;
	ARM_Curve*	maxCoupon = NULL;
	ARM_Curve*	fundingNominal = NULL;
	ARM_Curve*	fundingSpread = NULL;
	ARM_Curve*	fees = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator Init" );

		
		ARM_HybridIRFXCalculator*	tarnFXCalculator = dynamic_cast<ARM_HybridIRFXCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("TARNFXId").GetObjectId()));
		if (!tarnFXCalculator)
		{
			result.setMsg ("ARM_ERR: TARN FX is not of a good type");
			return ARM_KO;
		}

		ARM_MarketData_ManagerRep*	mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("MDM").GetObjectId()));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		int randGenType1 = ARM_ArgConv_BaseGenAlgoType.GetNumber(genericParams->GetParamValue("RandGenType1").GetString());
		int randGenAlgo1 = ARM_ArgConv_TransformAlgoType.GetNumber(genericParams->GetParamValue("RandGenAlgo1").GetString());
		int randGenType2 = ARM_ArgConv_BaseGenAlgoType.GetNumber(genericParams->GetParamValue("RandGenType2").GetString());
		int randGenAlgo2 = ARM_ArgConv_TransformAlgoType.GetNumber(genericParams->GetParamValue("RandGenAlgo2").GetString());

		int skipPDE = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("SkipPDE").GetString());

		int modelType = ARM_ArgConv_TARNFXModelType.GetNumber(genericParams->GetParamValue("ModelType").GetString());

		int oneFactorFlag = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("OneFactor").GetString());

		int smileFlag = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("Smile").GetString());

		int mixCalib = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("MixCalib").GetString());
		int rescalling = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("Rescalling").GetString());
		int correltype = ARM_ArgConv_MMCorrelType.GetNumber(genericParams->GetParamValue("CorrelType").GetString());
		

        newTarnFXCalculator = static_cast<ARM_TARNFXCalculator*>(tarnFXCalculator->Clone());
			
		newTarnFXCalculator->Init(
								(int)genericParams->GetParamValue("NbSimul").GetDouble(),
								(int)genericParams->GetParamValue("BucketSize").GetDouble(),
								randGenType1,
								randGenAlgo1,
								randGenType2,
								randGenAlgo2,
								genericParams->GetParamValue("FirstNbDims").GetDouble(),
								genericParams->GetParamValue("FirstNbTimes").GetDouble(),
								genericParams->GetParamValue("FactorNb").GetDouble(),
								genericParams->GetParamValue("TimeStepNb").GetDouble(),
								genericParams->GetParamValue("SpaceStepNb").GetDouble(),
								genericParams->GetParamValue("StdDevNb").GetDouble(),
								skipPDE,
								rescalling,
								(ARM_TARNFXCalculator::ModelType)modelType,
								smileFlag,
								mixCalib,
								oneFactorFlag,
								correltype,
								*mktDataManager
								);

		// assign object
		if ( !assignObject( newTarnFXCalculator, result, objId ) )
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
		delete	newTarnFXCalculator;

		x.DebugPrint();
		ARM_RESULT();
	}
}

//***********************************************************************
//**********  Indian TARN Calculator Creating  **************************
//***********************************************************************

long ARM_TARNCalculatorIndian_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_TARNCalculatorIndian* TARNCalculatorIndian = NULL;

	ARM_Curve*	cpnNominal = NULL;
	ARM_Curve*	barrierUp = NULL;
	ARM_Curve*	barrierDown = NULL;
	ARM_Curve*	strike = NULL;
	ARM_Curve*	fees = NULL;
	ARM_Curve*	targetCv = NULL;
	
	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];

		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);

		ARM_Currency cpnCurrency(genericParams->GetParamValue("CpnCcy").GetString().c_str());
		
		vector<ARM_Currency>	forCurrencies;
		forCurrencies.push_back((genericParams->GetParamValue("FgnCcy").GetString().c_str()));

		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(genericParams->GetParamValue("PayRec").GetString());
		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("CpnFreq").GetString());
		int cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("CpnDayCount").GetString());
		int stubRule	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString().c_str());
		int timing		= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("Timing").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());
		
		int optionType = K_CALL;
		if ((genericParams->GetParamValue("OptionType").GetString() == "Put")||
			(genericParams->GetParamValue("OptionType").GetString() == "P"))
			optionType = K_PUT;

		int indianType;
		if (genericParams->GetParamValue("IndianType").GetString() == "DownUp")
			indianType = ARM_TARNCalculatorIndian::DownUp;
		else if (genericParams->GetParamValue("IndianType").GetString() == "Trigger")
			indianType = ARM_TARNCalculatorIndian::Trigger;
		else if (genericParams->GetParamValue("IndianType").GetString() == "DigitalTrigger")
			indianType = ARM_TARNCalculatorIndian::DigitalTrigger;
		else
		{
			result.setMsg ("ARM_ERR: Indian type should be DownUp, Trigger or DigitalTrigger");
			return ARM_KO;
		}

		barrierDown = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("BarrierDown").GetObjectId()));
		if (!barrierDown && (indianType == ARM_TARNCalculatorIndian::DownUp))
		{
			result.setMsg ("ARM_ERR: barrier down should be a curve");
			return ARM_KO;
		}

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &cpnNominal, genericParams->GetParamValue("CpnNominal").GetObjectId() ," CpnNominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &barrierUp, genericParams->GetParamValue("BarrierUp").GetObjectId() ," BarrierUp", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &strike, genericParams->GetParamValue("Strike").GetObjectId() ," Strike", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fees, genericParams->GetParamValue("Fees").GetObjectId() ," Fees", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &targetCv, genericParams->GetParamValue("Target").GetObjectId() ," target", result ) ) return ARM_KO;

		int intermPrices = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("IntermPrices").GetString());

		ARM_StringVector columnsToPrice;
		ARM_GP_StrVector* productsToPrice = dynamic_cast<ARM_GP_StrVector*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("ProductsToPrice").GetObjectId()));;
		if (!productsToPrice)
			columnsToPrice.push_back("TARN");
		else
			columnsToPrice = ARM_StringVector(productsToPrice->GetValues());

        TARNCalculatorIndian = new ARM_TARNCalculatorIndian(asOfDate,
															startDate,
															endDate,
															cpnCurrency,
															forCurrencies,
															payRec,
															cpnDayCount,
															cpnFreq,
															genericParams->GetParamValue("CpnResetGap").GetDouble(),
															genericParams->GetParamValue("ResetCal").GetString(),
															genericParams->GetParamValue("PayCal").GetString(),
															stubRule,
															timing,
															intRule,
															*cpnNominal,
															*barrierUp,
															*barrierDown,
															*strike,
															*targetCv,
															genericParams->GetParamValue("Epsilon").GetDouble(),
															*fees,
															intermPrices,
															columnsToPrice,
															optionType,
															indianType);

		// assign object
		if ( !assignObject( TARNCalculatorIndian, result, objId ) )
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
		delete	TARNCalculatorIndian;

		x.DebugPrint();
		ARM_RESULT();
	}
	/*
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_TARNCalculatorIndian* TARNCalculatorIndian = NULL;

	ARM_Curve*	cpnNominal = NULL;
	ARM_Curve*	barrierUp = NULL;
	ARM_Curve*	barrierDown = NULL;
	ARM_Curve*	strike = NULL;
	ARM_Curve*	fees = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN FX Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];

		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);

		ARM_Currency cpnCurrency(genericParams->GetParamValue("CpnCcy").GetString().c_str());
		
		vector<ARM_Currency>	forCurrencies;
		forCurrencies.push_back((genericParams->GetParamValue("FgnCcy").GetString().c_str()));

		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(genericParams->GetParamValue("PayRec").GetString());
		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("CpnFreq").GetString());
		int cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("CpnDayCount").GetString());
		int stubRule	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString().c_str());
		int timing		= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("Timing").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());
		
		int optionType = K_CALL;
		if ((genericParams->GetParamValue("OptionType").GetString() == "Put")||
			(genericParams->GetParamValue("OptionType").GetString() == "P"))
			optionType = K_PUT;

		int indianType;
		if (genericParams->GetParamValue("IndianType").GetString() == "DownUp")
			indianType = ARM_TARNCalculatorIndian::DownUp;
		else if (genericParams->GetParamValue("IndianType").GetString() == "Trigger")
			indianType = ARM_TARNCalculatorIndian::Trigger;
		else if (genericParams->GetParamValue("IndianType").GetString() == "DigitalTrigger")
			indianType = ARM_TARNCalculatorIndian::DigitalTrigger;
		else
		{
			result.setMsg ("ARM_ERR: Indian type should be DownUp, Trigger or DigitalTrigger");
			return ARM_KO;
		}

		barrierDown = dynamic_cast<ARM_Curve*>(LOCAL_PERSISTENT_OBJECTS->GetObject(genericParams->GetParamValue("BarrierDown").GetObjectId()));
		if (!barrierDown && (indianType == ARM_TARNCalculatorIndian::DownUp))
		{
			result.setMsg ("ARM_ERR: barrier down should be a curve");
			return ARM_KO;
		}

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &cpnNominal, genericParams->GetParamValue("CpnNominal").GetObjectId() ," CpnNominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &barrierUp, genericParams->GetParamValue("BarrierUp").GetObjectId() ," BarrierUp", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &strike, genericParams->GetParamValue("Strike").GetObjectId() ," Strike", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fees, genericParams->GetParamValue("Fees").GetObjectId() ," Fees", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &targetCv, genericParams->GetParamValue("Target").GetObjectId() ," target", result ) ) return ARM_KO;


		int intermPrices = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("IntermPrices").GetString());
		ARM_IntVector productsToPrice(ARM_TARNCalculatorIndian::NbProductsToPrice);

		productsToPrice[ARM_TARNCalculatorIndian::SpotFwdPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("SpotFwdFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::ShortPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("ShortFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::LongPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("LongFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::TwiceCallPutPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("TwiceCallPutFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::CallMinusCallPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("CallMinusCallFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::PutMinusPutPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("PutMinusPutFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::CallDigitalPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("CallDigitalFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::PutDigitalPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("PutDigitalFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::IsAliveValue] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("IsAliveFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::CouponPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("CouponFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::PaidCouponPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("PaidCouponFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::RealCouponPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("RealCouponFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::ProbaValue] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("ProbaFlag").GetString());
		productsToPrice[ARM_TARNCalculatorIndian::TARNPrice] = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("TARNFlag").GetString());

        TARNCalculatorIndian = new ARM_TARNCalculatorIndian(asOfDate,
													startDate,
													endDate,
													cpnCurrency,
													forCurrencies,
													payRec,
													cpnDayCount,
													cpnFreq,
													genericParams->GetParamValue("CpnResetGap").GetDouble(),
													genericParams->GetParamValue("ResetCal").GetString(),
													genericParams->GetParamValue("PayCal").GetString(),
													stubRule,
													timing,
													intRule,
													*cpnNominal,
													*barrierUp,
													*barrierDown,
													*strike,
													*targetCv,
													genericParams->GetParamValue("Epsilon").GetDouble(),
													*fees,
													intermPrices,
													productsToPrice,
													optionType,
													indianType);

		// assign object
		if ( !assignObject( TARNCalculatorIndian, result, objId ) )
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
		delete	TARNCalculatorIndian;

		x.DebugPrint();
		ARM_RESULT();
	}*/
}

