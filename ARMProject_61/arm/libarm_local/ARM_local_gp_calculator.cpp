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

#include "ARM_local_gp_calculators.h"
#include "ARM_local_glob.h"
#include <ARM\local_xlarm\ARM_local_interglob.h>
#include <ARM\local_xlarm\ARM_gp_local_interglob.h>
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
#include <GP_Calculators\gpcalculators\volbondcalculator.h>

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


///////////////////////////////////////////////
//// Generic Function to set data to a sprecific Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_Calculator_Set(
        const long& calculatorId,
        const long& dataId,
        const string& setType,
		const vector< string >& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_LocalCSOCalculator* localCSO=NULL;
	ARM_CRFCalculator* localCRF=NULL;
	ARM_PRDCalculator* localPRDC=NULL;
	
	try
	{
		if( GetObjectFromIdWithDynamicCastCheck( &localCSO, calculatorId  ) )
			return ARMLOCAL_CSO_Set(calculatorId,dataId,setType,keys,isUpdated,result,objId);

		else if (GetObjectFromIdWithDynamicCastCheck( &localCRF, calculatorId  ) )
			return ARMLOCAL_CRF_Set(calculatorId,dataId,setType,keys,isUpdated,result,objId);

		else if (GetObjectFromIdWithDynamicCastCheck( &localPRDC, calculatorId  ) )
			return ARMLOCAL_PRDC_Set(calculatorId,dataId,setType,keys,isUpdated,result,objId);
		else
		{
			return ARMLOCAL_GC_Set(calculatorId,dataId,keys,isUpdated,result,objId);
			//result.setMsg ("ARM_ERR: cpn notional is not of a good type");
			//return ARM_KO;
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}

}
