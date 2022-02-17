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

////////////////////////////////////////////
//// Function to create a date strip combiner
////////////////////////////////////////////
////////////////////////////////////////////
//// Function to create a CRF Calculator
////////////////////////////////////////////
bool CRFCalculatorCheckAutoCalFlags(const vector< string >& flags,size_t nbFlagsToCheck)
{
    size_t nb=(flags.size()>nbFlagsToCheck ? nbFlagsToCheck : flags.size());
    for(size_t i=0;i<nb;++i)
        if(flags[i] != "Y" && flags[i] != "YES" && flags[i] != "N" && flags[i] != "NO")
            return false;

    return true;
}

extern long ARMLOCAL_CRFCalculator_Create(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
        long payRec,
        double fixEndDate,
        long fixDayCount,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
		long stubRule,
        long cpnResetGap,
        double leverage,
        long leverageId,
        double cpnMin,
        long cpnMinId,
        double cpnMax,
        long cpnMaxId,
        double fundSpread,
        long fundSpreadId,
        long fundFreq,
        long fundDayCount,
        double nominal,
        long nominalId,
        long exerGap,
        long nbNonCall,
        double exerFee,
        long exerFeeId,
        const vector< string >& flags,
        long mktDataManagerId,
        const vector< string >& keys,
        double fundnominal,
        long fundnominalId,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* crfCalculator = NULL;

	ARM_ReferenceValue* strikeProfile=NULL;
    bool isStrike=false;

	ARM_ReferenceValue* leverageProfile=NULL;
    bool isLeverage=false;

	ARM_ReferenceValue* cpnMinProfile=NULL;
    bool isCpnMin=false;

	ARM_ReferenceValue* cpnMaxProfile=NULL;
    bool isCpnMax=false;

	ARM_ReferenceValue* fundSpreadProfile=NULL;
    bool isFundSpread=false;

	ARM_ReferenceValue* nominalProfile=NULL;
    bool isNominal=false;

	ARM_ReferenceValue* exerFeeProfile=NULL;
    bool isExerFee=false;

    ARM_ReferenceValue* fundnominalProfile=NULL;
    bool isFundNominal=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "CRF Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		char myFixEndDate[20];
		

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		Local_XLDATE2ARMDATE(fixEndDate, myFixEndDate);

        /// Convert refValue Id to object if possible

        if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            strikeProfile = new ARM_ReferenceValue(strike);
            isStrike=true;
        }

        if(leverageId != ARM_NULL_OBJECT)
        {
		    leverageProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId));

		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: leverage is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            leverageProfile = new ARM_ReferenceValue(leverage);
            isLeverage=true;
        }

        if(cpnMinId != ARM_NULL_OBJECT)
        {
		    cpnMinProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMinId));

		    if (!cpnMinProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon min is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnMinProfile = new ARM_ReferenceValue(cpnMin);
            isCpnMin=true;
        }

        if(cpnMaxId != ARM_NULL_OBJECT)
        {
		    cpnMaxProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMaxId));

		    if (!cpnMaxProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon max is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnMaxProfile = new ARM_ReferenceValue(cpnMax);
            isCpnMax=true;
        }

        if(fundSpreadId != ARM_NULL_OBJECT)
        {
		    fundSpreadProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: funding spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundSpreadProfile = new ARM_ReferenceValue(fundSpread);
            isFundSpread=true;
        }

        if(nominalId != ARM_NULL_OBJECT)
        {
		    nominalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(nominalId));

		    if (!nominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            nominalProfile = new ARM_ReferenceValue(nominal);
            isNominal=true;
        }

        if(fundnominalId != ARM_NULL_OBJECT)
        {
		    fundnominalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundnominalId));

		    if (!fundnominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else if (fundnominal!=ARM_MISSING_VALUE)
        {
            fundnominalProfile = new ARM_ReferenceValue(fundnominal);
            isFundNominal=true;
        }
        else
        {
            fundnominalProfile = new ARM_ReferenceValue();
            isFundNominal=true;
        }


        if(exerFeeId != ARM_NULL_OBJECT)
        {
		    exerFeeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(exerFeeId));

		    if (!exerFeeProfile)
		    {
			    result.setMsg ("ARM_ERR: exercise fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            exerFeeProfile = new ARM_ReferenceValue(exerFee); 
            isExerFee=true;
        }


        /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}


        /// Convert autocal flags
		string sigmaTypeStr = (flags[0]=="N" || flags[0]=="NO") ?  string("UNKNOWN"):
							    (flags[0]=="Y" || flags[0]=="YES")  ? string("EQUIVALENT"):flags[0];
		ARM_SigmaCalibType sigmaType = (ARM_SigmaCalibType) ARM_ArgConv_SigmaCalibType.GetNumber(sigmaTypeStr);

		string mrsTypeStr = (flags[1]=="Y" || flags[1]=="YES") ?  string("STM_FIRST_COLUMN"):
							(flags[1]=="N" || flags[1]=="NO")  ? string("UNKNOWN"):flags[1];
		ARM_MRSCalibType mrsType = (ARM_MRSCalibType) ARM_ArgConv_MRSCalibType.GetNumber(mrsTypeStr);

        bool capCalibFlag=(flags[2]=="Y" || flags[2]=="YES");
        bool floorCalibFlag=(flags[3]=="Y" || flags[3]=="YES");
		bool skewCalFlag = false;

		ARM_ModelType modelType = ARM_PricingModelType::HWM1F;
		if (flags.size()>5)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(flags[5]);
		if (flags.size()>6)
			skewCalFlag=(flags[6]=="Y" || flags[6]=="YES");

		ARM_MRSStrikeCalibType	mrsStrikeType  = ARM_MRSStrikeCalibrationType::strikeEquivalent;
		if (flags.size()>7)
			 mrsStrikeType = (ARM_MRSStrikeCalibType) ARM_ArgConv_MRSStrikeCalibType.GetNumber(flags[7]);

		long isFrontier = 0;
		if (flags.size()>8)
			isFrontier = atof(flags[8].c_str());

        /// Create the CRF calculator
        crfCalculator = new ARM_CRFCalculator(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
            *strikeProfile,
            payRec,
            ARM_Date(myFixEndDate),
            fixDayCount,
            cpnDayCount,
            cpnFreq,
            cpnTiming,
            cpnIndexTerm,
            cpnIndexDayCount,
            cpnResetCal,
            cpnPayCal,
			stubRule,
            cpnResetGap,
            *leverageProfile,
            *cpnMinProfile,
            *cpnMaxProfile,
            *fundSpreadProfile,
            fundFreq,
            fundDayCount,
            *nominalProfile,
            exerGap,
            nbNonCall,
            *exerFeeProfile,
            sigmaType,
            mrsType,
			mrsStrikeType,
            capCalibFlag,
            floorCalibFlag,
            *mktDataManager,
            keys,
            *fundnominalProfile,
			modelType,
			skewCalFlag,
			isFrontier);


        /// Set product to price flag
        if(!CRF_SetProductToPrice(crfCalculator,flags[4]))
        {
			result.setMsg ("ARM_ERR: Unknown product to be priced within CRF calculator");
			return ARM_KO;
        }


        /// Free memory
        if(isStrike)
            delete strikeProfile;
		strikeProfile = NULL;
        if(isLeverage)
            delete leverageProfile;
		leverageProfile = NULL;
        if(isCpnMin)
            delete cpnMinProfile;
		cpnMinProfile = NULL;
        if(isCpnMax)
            delete cpnMaxProfile;
		cpnMaxProfile = NULL;
        if(isFundSpread)
            delete fundSpreadProfile;
		fundSpreadProfile = NULL;
        if(isNominal)
            delete nominalProfile;
		nominalProfile = NULL;
        if(isFundNominal)
            delete fundnominalProfile;

        if(isExerFee)
            delete exerFeeProfile;
		exerFeeProfile = NULL;

		/// assign object
		if( !assignObject( crfCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if(isStrike)
            delete strikeProfile;
        if(isLeverage)
            delete leverageProfile;
        if(isCpnMin)
            delete cpnMinProfile;
        if(isCpnMax)
            delete cpnMaxProfile;
        if(isFundSpread)
            delete fundSpreadProfile;
        if(isNominal)
            delete nominalProfile;
        if(isFundNominal)
            delete fundnominalProfile;
        if(isExerFee)
            delete exerFeeProfile;

		delete crfCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}

// Create a CRF Calculator without flags and Market data

extern long ARMLOCAL_Crude_CRFCalculator_Create(
        const double & startDate,
        const double & endDate,
        const double & strike,
        const long & strikeId,
        const long & payRec,
        const double & fixEndDate,
        const long & fixDayCount,
        const long & cpnDayCount,
        const long & cpnFreq,
        const long & cpnTiming,
        const string & cpnIndexTerm,
        const long & cpnIndexDayCount,
        const string & cpnResetCal,
        const string & cpnPayCal,
		const long & stubRule,
        const long & cpnResetGap,
        const double & leverage,
        const long & leverageId,
        const double & cpnMin,
        const long & cpnMinId,
        const double & cpnMax,
        const long & cpnMaxId,
        const double & fundSpread,
        const long & fundSpreadId,
        const long & fundFreq,
        const long & fundDayCount,
        const double & nominal,
        const long & nominalId,
        const long & exerGap,
        const long & nbNonCall,
        const double & exerFee,
        const long & exerFeeId,
        const double & fundnominal,
        const long & fundnominalId,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* crfCalculator = NULL;

	ARM_ReferenceValue* strikeProfile=NULL;
    bool isStrike=false;

	ARM_ReferenceValue* leverageProfile=NULL;
    bool isLeverage=false;

	ARM_ReferenceValue* cpnMinProfile=NULL;
    bool isCpnMin=false;

	ARM_ReferenceValue* cpnMaxProfile=NULL;
    bool isCpnMax=false;

	ARM_ReferenceValue* fundSpreadProfile=NULL;
    bool isFundSpread=false;

	ARM_ReferenceValue* nominalProfile=NULL;
    bool isNominal=false;

	ARM_ReferenceValue* exerFeeProfile=NULL;
    bool isExerFee=false;

    ARM_ReferenceValue* fundnominalProfile=NULL;
    bool isFundNominal=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "CRF Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		char myFixEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		Local_XLDATE2ARMDATE(fixEndDate, myFixEndDate);

        /// Convert refValue Id to object if possible

        if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            strikeProfile = new ARM_ReferenceValue(strike);
            isStrike=true;
        }

        if(leverageId != ARM_NULL_OBJECT)
        {
		    leverageProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId));

		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: leverage is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            leverageProfile = new ARM_ReferenceValue(leverage);
            isLeverage=true;
        }

        if(cpnMinId != ARM_NULL_OBJECT)
        {
		    cpnMinProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMinId));

		    if (!cpnMinProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon min is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnMinProfile = new ARM_ReferenceValue(cpnMin);
            isCpnMin=true;
        }

        if(cpnMaxId != ARM_NULL_OBJECT)
        {
		    cpnMaxProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMaxId));

		    if (!cpnMaxProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon max is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnMaxProfile = new ARM_ReferenceValue(cpnMax);
            isCpnMax=true;
        }

        if(fundSpreadId != ARM_NULL_OBJECT)
        {
		    fundSpreadProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: funding spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundSpreadProfile = new ARM_ReferenceValue(fundSpread);
            isFundSpread=true;
        }

        if(nominalId != ARM_NULL_OBJECT)
        {
		    nominalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(nominalId));

		    if (!nominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            nominalProfile = new ARM_ReferenceValue(nominal);
            isNominal=true;
        }

        if(fundnominalId != ARM_NULL_OBJECT)
        {
		    fundnominalProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundnominalId));

		    if (!fundnominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else if (fundnominal!=ARM_MISSING_VALUE)
        {
            fundnominalProfile = new ARM_ReferenceValue(fundnominal);
            isFundNominal=true;
        }
        else
        {
            fundnominalProfile = new ARM_ReferenceValue();
            isFundNominal=true;
        }


        if(exerFeeId != ARM_NULL_OBJECT)
        {
		    exerFeeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(exerFeeId));

		    if (!exerFeeProfile)
		    {
			    result.setMsg ("ARM_ERR: exercise fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            exerFeeProfile = new ARM_ReferenceValue(exerFee); 
            isExerFee=true;
        }


        /// Create the CRF calculator
       /* crfCalculator = new ARM_CRFCalculator(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
            *strikeProfile,
            payRec,
            ARM_Date(myFixEndDate),
            fixDayCount,
            cpnDayCount,
            cpnFreq,
            cpnTiming,
            cpnIndexTerm,
            cpnIndexDayCount,
            cpnResetCal,
            cpnPayCal,
			stubRule,
            cpnResetGap,
            *leverageProfile,
            *cpnMinProfile,
            *cpnMaxProfile,
            *fundSpreadProfile,
            fundFreq,
            fundDayCount,
            *nominalProfile,
            exerGap,
            nbNonCall,
            *exerFeeProfile,
            *fundnominalProfile);*/


        /// Free memory
        if(isStrike)
            delete strikeProfile;
		strikeProfile = NULL;
        if(isLeverage)
            delete leverageProfile;
		leverageProfile = NULL;
        if(isCpnMin)
            delete cpnMinProfile;
		cpnMinProfile = NULL;
        if(isCpnMax)
            delete cpnMaxProfile;
		cpnMaxProfile = NULL;
        if(isFundSpread)
            delete fundSpreadProfile;
		fundSpreadProfile = NULL;
        if(isNominal)
            delete nominalProfile;
		nominalProfile = NULL;
        if(isFundNominal)
            delete fundnominalProfile;

        if(isExerFee)
            delete exerFeeProfile;
		exerFeeProfile = NULL;

		/// assign object
		if( !assignObject( crfCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if(isStrike)
            delete strikeProfile;
        if(isLeverage)
            delete leverageProfile;
        if(isCpnMin)
            delete cpnMinProfile;
        if(isCpnMax)
            delete cpnMaxProfile;
        if(isFundSpread)
            delete fundSpreadProfile;
        if(isNominal)
            delete nominalProfile;
        if(isFundNominal)
            delete fundnominalProfile;
        if(isExerFee)
            delete exerFeeProfile;

		delete crfCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}



///////////////////////////////////////////////
//// Function to get data from a CRF Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_CRF_Get(
        long crfId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* crfCalculator;
    ARM_Object* object=NULL;

	try
	{
		crfCalculator = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));

		if (!crfCalculator)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_OSW_PORT || 
			typeToGet == GC_ACCESS_SHORT_TERM_PORT ||
            typeToGet == GC_ACCESS_CF_PORT || 
			typeToGet == GC_ACCESS_OSW_SECOND_PORT || 
			typeToGet == GC_ACCESS_OSW_BSMODEL)
        {
            if(typeToGet == GC_ACCESS_OSW_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(crfCalculator->GetOSWPortfolio())).Clone();
            else if(typeToGet == GC_ACCESS_SHORT_TERM_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(crfCalculator->GetSTMPortfolio())).Clone();
			else if(typeToGet==GC_ACCESS_OSW_BSMODEL)
				object = static_cast<ARM_BSModel> (*(crfCalculator->GetOSWBSModel())).Clone();
            else if(typeToGet==GC_ACCESS_CF_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(crfCalculator->GetCFPortfolio())).Clone();
			else if(typeToGet==GC_ACCESS_OSW_SECOND_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(crfCalculator->GetSkewPortfolio())).Clone();

            /// Assign the object in ARM cache
		    if( !assignObject( object, result, objId ) )
            {
                delete object;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
        else
            return ARMLOCAL_GC_Get(crfId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_CRF_GetMRS(
        long crfId,
        ARM_result&	result,
		long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* crfCalculator;
    ARM_Object* object=NULL;
	ARM_Object* prevObject=NULL;
	long mrsId;

	try
	{
		crfCalculator = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		double asOfDate = crfCalculator->GetMktDataManager()->GetAsOfDate().GetJulian();

		if (!crfCalculator)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
		ARM_CurveModelParam* MRS = const_cast<ARM_CurveModelParam*>(&crfCalculator->GetMRS());
		object = CurveToRefValue(*MRS->GetCurve(), asOfDate);

		if(objId == -1)
		{
			CREATE_GLOBAL_OBJECT();

			mrsId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)object);

			if (mrsId == RET_KO)
			{
				if (object)
					delete object;
				object = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(mrsId);

			return ARM_OK;
		}
		else
		{
			prevObject = (ARM_ReferenceValue *) LOCAL_PERSISTENT_OBJECTS->GetObject(objId);

			if (LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(prevObject, ARM_REFERENCE_VALUE) == 1)
			{
				if (prevObject)
				{
					delete prevObject;
					prevObject = NULL;
				}

				LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)object, objId);

				return ARM_OK;
			}
			else
			{
				if (object)
					delete object;
				object = NULL;

				result.setMsg ("ARM_ERR: previous object is not of a good type");
				return ARM_KO;
			}
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a CRF Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_CRF_Set(
        long crfId,
        long dataId,
        const string& setType,
		const vector< string >& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* oldCRF=NULL;
    ARM_CRFCalculator* newCRF=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");


	try
	{
		oldCRF = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		if (!oldCRF)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

		ARM_ModelParam* param = NULL;
        ARM_StdPortfolio* port = NULL;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
            if(isUpdated)
                /// Just update the initial CRF
                newCRF=oldCRF;
            else
                /// Clone the initial CRF then set
                newCRF = dynamic_cast<ARM_CRFCalculator *>(oldCRF->Clone());

            string typeToSet(stringGetUpper(setType));
            bool isUpdateOswStrike=true;
			bool updateCalibration = true;
            if(typeToSet == GC_ACCESS_CF_PORT)
            {
                newCRF->SetCFPortfolio(*port);
                updateMsg += "C/F portfolio";
            }
            else if(typeToSet == GC_ACCESS_SHORT_TERM_PORT)
            {
                newCRF->SetSTMPortfolio(*port);
                updateMsg += "STM vanillas portfolio";
				isUpdateOswStrike = false;
            }
		    else
            {
                /// Default
                newCRF->SetOSWPortfolio(*port);
                updateMsg += "OSW portfolio";
                isUpdateOswStrike = false;
            }

            /// Keep calibration data consistency
            if(updateCalibration)
				newCRF->UpdateCalibrationAndTimeIt(isUpdateOswStrike);

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRF in ARM cache
		    if(!assignObject( newCRF, result, objId ) )
            {
                delete newCRF;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else if (param = dynamic_cast< ARM_ModelParam* >(object))
		{
			newCRF = dynamic_cast<ARM_CRFCalculator *>(oldCRF->Clone());
			newCRF->SetCalibParam(param);
			updateMsg += "Model Param";
			
			/// Assign the new CRF in ARM cache
		    if(!assignObject( newCRF, result, objId ) )
            {
                delete newCRF;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
		}
		else
            return ARMLOCAL_GC_Set(crfId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newCRF;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to manage CRF auto-calibrations
////////////////////////////////////////////
extern long ARMLOCAL_CRF_SetAutoCalFlags(
	long crfId,
    vector< string >& flags,
    ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


    ARM_CRFCalculator* crf=NULL;
	try
	{
		crf = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		if (!crf)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}

        /// Convert autocal flags
        if(!CRFCalculatorCheckAutoCalFlags(flags,4))
		{
			result.setMsg ("ARM_ERR: autocal flags are not of a good type");
			return ARM_KO;
		}
		string sigmaTypeStr = (flags[0]=="N" || flags[0]=="NO") ?  string("UNKNOWN"):
							    (flags[0]=="Y" || flags[0]=="YES")  ? string("EQUIVALENT"):flags[0];
		ARM_SigmaCalibType sigmaType = (ARM_SigmaCalibType) ARM_ArgConv_SigmaCalibType.GetNumber(sigmaTypeStr);
		crf->SetOSWCalibFlag(sigmaType);

		string mrsTypeStr = (flags[1]=="Y" || flags[1]=="YES") ?  string("STM_FIRST_COLUMN"):
							(flags[1]=="N" || flags[1]=="NO")  ? string("UNKNOWN"):flags[1];
		ARM_MRSCalibType mrsType = (ARM_MRSCalibType) ARM_ArgConv_MRSCalibType.GetNumber(mrsTypeStr);

        bool capCalibFlag=(flags[2]=="Y" || flags[2]=="YES");
        bool floorCalibFlag=(flags[3]=="Y" || flags[3]=="YES");
		
		
		crf->SetSTMAndUpdateCalibFlag(mrsType);
		crf->SetCapCalibFlag(capCalibFlag);
		crf->SetFloorCalibFlag(floorCalibFlag);
		
		string txt( "AutoCal Flags : Swaption=" );
		txt += ARM_ArgConvReverse_SigmaCalibType.GetString(crf->GetOSWCalibFlag());
        txt += ", STM Vanilla=";
		txt += ARM_ArgConvReverse_MRSCalibType.GetString(crf->GetSTMCalibFlag());
        txt += ", Cap=";
		txt += crf->GetCapCalibFlag() ? "On" : "Off";
        txt += ", Floor=";
		txt += crf->GetFloorCalibFlag() ? "On" : "Off";
		
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to manage CRF auto-calibrations
////////////////////////////////////////////
extern long ARMLOCAL_CRF_SetOneCalFlag(
	const long& crfId,
    const string& flagStr,
	const string& typeStr,
    ARM_result&	result, 
	long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


    ARM_CRFCalculator* crf=NULL;
	try
	{
		crf = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		if (!crf)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}


        bool falg=(flagStr=="Y" || flagStr=="YES");
		
		string txt( "Auto Calibration : Param =" );
		if(typeStr == string("SKEW"))
		{
		    crf->SetSkewReCalibFlag(falg);
			txt+="SKEW";
			 }
		else
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT, "invalid Paramter type" );
		
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}
////////////////////////////////////////////
//// Function to manage CRF auto-calibrations and returns a new CRF
////////////////////////////////////////////
extern long ARMLOCAL_CRF_SetAutoCalFlagsAndClone(
	const long& crfId,
    const vector< string >& flags,
    ARM_result&	result,
	long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_CRFCalculator* newCRF=NULL;
	ARM_CRFCalculator* crf=NULL;
	try
	{
		crf = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		if (!crf)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}
		/// Clone the initial CRF then set
		newCRF = dynamic_cast<ARM_CRFCalculator *>(crf->Clone());

        /// Convert autocal flags
        if(!CRFCalculatorCheckAutoCalFlags(flags,4))
		{
			result.setMsg ("ARM_ERR: autocal flags are not of a good type");
			return ARM_KO;
		}
	    bool oswCalibFlag=(flags[0]=="Y" || flags[0]=="YES");

		string sigmaTypeStr = (flags[0]=="N" || flags[0]=="NO") ?  string("UNKNOWN"):
							    (flags[0]=="Y" || flags[0]=="YES")  ? string("EQUIVALENT"):flags[0];
		ARM_SigmaCalibType sigmaType = (ARM_SigmaCalibType) ARM_ArgConv_SigmaCalibType.GetNumber(sigmaTypeStr);
		newCRF->SetOSWCalibFlag(sigmaType);

        string mrsTypeStr = (flags[1]=="Y" || flags[1]=="YES") ?  string("UNKNOWN"):
							(flags[1]=="N" || flags[1]=="NO")  ? string("STM_FIRST_COLUMN"):flags[1];
		ARM_MRSCalibType mrsType = (ARM_MRSCalibType) ARM_ArgConv_MRSCalibType.GetNumber(mrsTypeStr);
        bool capCalibFlag=(flags[2]=="Y" || flags[2]=="YES");
        bool floorCalibFlag=(flags[3]=="Y" || flags[3]=="YES");
				
		newCRF->SetSTMAndUpdateCalibFlag(mrsType);
		newCRF->SetCapCalibFlag(capCalibFlag);
		newCRF->SetFloorCalibFlag(floorCalibFlag);
		
		/// Assign the new CRF in ARM cache
		if(!assignObject( newCRF, result, objId ) )
        {
			delete newCRF;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}
////////////////////////////////////////////
//// Function to manage product to price
//// within the CRF deal description
////////////////////////////////////////////
bool CRF_SetProductToPrice(ARM_Object* crf, const string& prodToPrice)
{
    ARM_CRFCalculator* theCRF = static_cast< ARM_CRFCalculator* >(crf);
    if(prodToPrice == "CRF")
        theCRF->SetCRFToPrice();
    else if(prodToPrice == "CAP")
        theCRF->SetCapToPrice();
    else if(prodToPrice == "FLOOR")
        theCRF->SetFloorToPrice();
    else if(prodToPrice == "FUNDING")
        theCRF->SetFundingToPrice();
    else if(prodToPrice == "STDLEG")
        theCRF->SetStdLegToPrice();
    else if(prodToPrice == "RFLEG")
        theCRF->SetRFLegToPrice();
    else if(prodToPrice == "STDSWAP")
        theCRF->SetStdSwapToPrice();
    else if(prodToPrice == "RFSWAP")
        theCRF->SetRFSwapToPrice();
    else if(prodToPrice == "BERMUDA")
        theCRF->SetBermudaToPrice();
    else
        return false;

    return true;
}

extern long ARMLOCAL_CRF_SetProductToPrice(
	long crfId,
	string prodToPrice,
	ARM_result&	result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");


    ARM_CRFCalculator* crf=NULL;

	try
	{
		crf = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(crfId));
		if (!crf)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}

		string txt( "Product to be priced=" );
        if(!CRF_SetProductToPrice(crf,prodToPrice))
        {
			result.setMsg ("ARM_ERR: Unknown product to be priced within CRF calculator");
			return ARM_KO;
        }
		txt += prodToPrice;
		result.setString(txt.c_str());
		
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}


						////////////////////////////////////////////
						////          Caption Calculator
						////////////////////////////////////////////

extern long ARMLOCAL_CaptionCalculator_Create(
        double startDate,
        double endDate,
		const string& cpnIdxTerm,
		long payRec,
		long CF,
		double coupon,
        long couponId,
		const string& FundIdxTerm,
		long NotifDays,
		long NonCall,
		double exercise,
		long exerStyleId,
		double notional,
        long notionalId,
		long cpnDayCount,
		long cpnResetTiming,
		const string& cpnResetCal,
        const string& cpnPayCal,
		double cpnSpread,
        long cpnSpreadId,
		long fundDayCount,
		const string& fundResetCal,
        const string& fundPayCal,
		double fundSpread,
        long fundSpreadId,
		long factorNb,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CaptionCalculator* captionCalculator = NULL;
	
	ARM_Curve* couponProfile=NULL;
    bool isCoupon=false;

	ARM_Curve* exerciseProfile=NULL;
    bool isExercise=false;

	ARM_Curve* notionalProfile=NULL;
    bool isNotional=false;

	ARM_Curve* cpnSpreadProfile=NULL;
    bool isCpnSpread=false;

	ARM_Curve* fundSpreadProfile=NULL;
    bool isFundSpread=false;

	ARM_Curve* feesProfile=NULL;
    bool isFees=false;


	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Caption Calculator" );

		char myStartDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

        /// Convert curve Id to object if possible
		//coupon
        if(couponId != ARM_NULL_OBJECT)
        {
		    couponProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(couponId));

		    if (!couponProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,coupon);
            couponProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCoupon=true;
        }

		//exerStyleId
		if(exerStyleId != ARM_NULL_OBJECT)
        {
		    exerciseProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(exerStyleId));

		    if (!exerciseProfile)
		    {
			    result.setMsg ("ARM_ERR: exerciseStyle  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,exercise);
            exerciseProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isExercise=true;
        }

		//notional
		if(notionalId != ARM_NULL_OBJECT)
        {
		    notionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId));

		    if (!notionalProfile)
		    {
			    result.setMsg ("ARM_ERR: notional  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,notional);
            notionalProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isNotional=true;
        }

		//cpnSpread
		if(cpnSpreadId != ARM_NULL_OBJECT)
        {
		    cpnSpreadProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnSpreadId));

		    if (!cpnSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon Spread  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,cpnSpread);
            cpnSpreadProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCpnSpread=true;
        }


		//fundSpread
		if(fundSpreadId != ARM_NULL_OBJECT)
        {
		    fundSpreadProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: Funding Spread  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundSpread);
            fundSpreadProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFundSpread=true;
        }

	    /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		if(flags.size() > 5)
		{
			result.setMsg ("ARM_ERR: flags size has to be less than 14.");
			return ARM_KO;
		}

		long SFRMVolType;
		const char*  volType = CalibMod[0].c_str();
		if ((SFRMVolType = ARM_ConvShapeType (CCString(volType), result)) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: Vol type is not of a good type");
			return ARM_KO;
		}
		
		ARM_CaptionCalculator::CalibrationMode SWOPTCalibMode = (ARM_CaptionCalculator::CalibrationMode) ARM::ARM_ArgConv_CaptionCalibMode.GetNumber(CalibMod[1]); 
		ARM_CaptionCalculator::CalibrationMode BETACalibMode = (ARM_CaptionCalculator::CalibrationMode) ARM::ARM_ArgConv_CaptionCalibMode.GetNumber(CalibMod[2]);
	    
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(6,false);
		productsToPrice[0] = true;

		for (size_t i = 0; i < flags.size(); ++i)
		{
			productsToPrice[i+1] = ((flags[i] == "Y") || (flags[i] == "YES"));
		}

        /// Create the Caption calculator
		captionCalculator = new ARM_CaptionCalculator(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
            cpnIdxTerm,
			payRec,
			CF,
			*couponProfile,
			FundIdxTerm,
			NotifDays,
			NonCall,
			*exerciseProfile,		
			*notionalProfile,
			cpnDayCount,
			cpnResetTiming,
			cpnResetCal,
			cpnPayCal,
			*cpnSpreadProfile,
			fundDayCount,
			fundResetCal,
			fundPayCal,
			*fundSpreadProfile,
			factorNb,
			SFRMVolType,
			SWOPTCalibMode,
			BETACalibMode,
			productsToPrice,
			*mktDataManager,
			keys);

	    /// Free memory
        if(isCoupon)
            delete couponProfile;
		couponProfile = NULL;

        if(isExercise)
            delete exerciseProfile;
		exerciseProfile = NULL;

        if(isNotional)
            delete notionalProfile;
		notionalProfile = NULL;

        if(isCpnSpread)
            delete cpnSpreadProfile;
		cpnSpreadProfile = NULL;

        if(isFundSpread)
            delete fundSpreadProfile;
		fundSpreadProfile = NULL;
       
		/// assign object
		if( !assignObject( captionCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isCoupon)
            delete couponProfile;
         if(isExercise)
            delete exerciseProfile;
        if(isNotional)
            delete notionalProfile;
        if(isCpnSpread)
            delete cpnSpreadProfile;
        if(isFundSpread)
            delete fundSpreadProfile;
	
		delete captionCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////////////////
//// Function to get data from a Caption Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_Caption_Get(
        long captionId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CaptionCalculator* captionCalculator;
    ARM_Object* object=NULL;

	try
	{
		captionCalculator = dynamic_cast<ARM_CaptionCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(captionId));

		if (!captionCalculator)
		{
			result.setMsg ("ARM_ERR: Caption Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_OSW_PORT ||
            typeToGet == GC_ACCESS_CF_PORT)
        {
            if(typeToGet == GC_ACCESS_OSW_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(captionCalculator->GetSWOPTPortfolio())).Clone();
            else if (typeToGet == GC_ACCESS_CF_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(captionCalculator->GetCFPortfolio())).Clone();

            /// Assign the object in ARM cache
		    if( !assignObject( object, result, objId ) )
            {
                delete object;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
        else
            return ARMLOCAL_GC_Get(captionId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a Caption Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_Caption_Set(
        long captionId,
        long dataId,
        const string& setType,
		const vector< string >& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CaptionCalculator* oldCaption=NULL;
    ARM_CaptionCalculator* newCaption=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");


	try
	{
		oldCaption = dynamic_cast<ARM_CaptionCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(captionId));
		if (!oldCaption)
		{
			result.setMsg ("ARM_ERR: Caption Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
            if(isUpdated)
                /// Just update the initial Caption
                newCaption=oldCaption;
            else
                /// Clone the initial Caption then set
                newCaption = dynamic_cast<ARM_CaptionCalculator *>(oldCaption->Clone());

            string typeToSet(stringGetUpper(setType));
            if(typeToSet == GC_ACCESS_CF_PORT)
            {
                newCaption->SetCFPortfolio(*port);
                updateMsg += "C/F portfolio";
            }
            else if(typeToSet == GC_ACCESS_OSW_PORT)
            {
                /// Default
                newCaption->SetSWOPTPortfolio(*port);
                updateMsg += "OSW portfolio";
            }

            /// Keep calibration data consistency
            newCaption->UpdateCalibrationAndTimeIt();

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new Caption in ARM cache
		    if(!assignObject( newCaption, result, objId ) )
            {
                delete newCaption;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
        else
            return ARMLOCAL_GC_Set(captionId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newCaption;
		x.DebugPrint();
		ARM_RESULT();
	}
}


/////////////////////////////////////////////////////////////////////
/// CSO calculator creation : standard version
///	--> leverageLong = leverageShort
/// --> strike = 0.0
/// --> will work for HW1F + basket calib + local vol and for SFRM
/////////////////////////////////////////////////////////////////////
long ARMLOCAL_CSOCalculator_Create(double startDate,
								   double endDate,
								   long CMS1Type,
								   long CMS2Type,
								   long cpnFreq,
								   long cpnDaycount,
								   long fundFreq,
								   long fundDaycount,
								   long exerFreq,
								   long NotifDays,
								   double notional,
								   long notionalId,
								   double minCpn,
								   long minCpnId,
								   double maxCpn,
								   long maxCpnId,
								   double leverage,
								   long leverageId,
								   double fundMargin,
								   long fundMarginId,
								   long fixCpnId,
								   long feesId,
								   const vector< string >& CalibMod,
								   const vector< string >& ProductFlags,
								   const vector< double >& ModelDatas,
								   long mktDataManagerId,
								   const vector< string >& keys,
								   ARM_result&	result, 
								   long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_GenCSOCalculator* csoCalculator = NULL;
	
	ARM_Curve* notionalProfile=NULL;
    bool isNotional=false;

	ARM_Curve* minCpnProfile=NULL;
    bool isMinCpn=false;

	ARM_Curve* maxCpnProfile=NULL;
    bool isMaxCpn=false;

	ARM_Curve* leverageProfile=NULL;
    bool isLeverage=false;

	ARM_Curve* fundMarginProfile=NULL;
    bool isFundMargin=false;

	ARM_Curve* feesProfile=NULL;
	ARM_Curve* fixCpnProfile=NULL;

	ARM_Curve* strikeProfile=NULL;
	

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable SpreadOption Calculator" );

		char myStartDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

        /// Convert curve Id to object if possible
		//notional
        if(notionalId != ARM_NULL_OBJECT)
        {
		    notionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(notionalId));

		    if (!notionalProfile)
		    {
			    result.setMsg ("ARM_ERR: notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            notionalProfile = new ARM_FlatCurve(notional);
            isNotional=true;
        }

		//minCpn
		if(minCpnId != ARM_NULL_OBJECT)
        {
		    minCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(minCpnId));

		    if (!minCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: minCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,minCpn);
            minCpnProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isMinCpn=true;
        }

		//maxCpn
		if(maxCpnId != ARM_NULL_OBJECT)
        {
		    maxCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(maxCpnId));

		    if (!maxCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: maxCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,maxCpn);
            maxCpnProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isMaxCpn=true;
        }

		//fundMargin
		if(fundMarginId != ARM_NULL_OBJECT)
        {
		    fundMarginProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundMarginId));

		    if (!fundMarginProfile)
		    {
			    result.setMsg ("ARM_ERR: fundMargin  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundMargin);
            fundMarginProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFundMargin=true;
        }

		//fundMargin
		if(leverageId != ARM_NULL_OBJECT)
        {
		    leverageProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId));

		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: leverage  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,leverage);
            leverageProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isLeverage=true;
        }

		feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));

		if (!feesProfile)
		{
			result.setMsg ("ARM_ERR: fees is not of a good type");
			return ARM_KO;
		}

		fixCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fixCpnId));

		if (!fixCpnProfile)
		{
			result.setMsg ("ARM_ERR: fix Coupon is not of a good type");
			return ARM_KO;
		}


		/// strike profile = 0 for the moment
		double value(0.0);
        strikeProfile = new ARM_FlatCurve(value);
        
	    /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}
		
		bool switchToHWwithLocalVol (false);
		
		if (CalibMod.size() == 4)
			switchToHWwithLocalVol = false;	
		else if (CalibMod.size() == 2)
			switchToHWwithLocalVol = true;
		else
		{
			result.setMsg ("ARM_ERR: calib flags size = 2 or 4");
			return ARM_KO;
		}


		//--------------------------------------------
		// HW + local vol 
		//--------------------------------------------
		if (switchToHWwithLocalVol)
		{
			/// set calibration type and model type
			ARM_LocalCSOCalculator::CalibrationType calibType ;
			ARM_LocalCSOCalculator::CalibStrikeType calibStrikeType ;
			
			if (CalibMod[0] == "DIAG") 
				calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;

			else if (CalibMod[0] == "BASKET") 
				calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
			else
			{
				result.setMsg ("ARM_ERR: invalid calib type");
				return ARM_KO;
			}

			if (CalibMod[1] == "ATM") 
				calibStrikeType = ARM_LocalCSOCalculator::ATM;

			else if (CalibMod[1] == "EQUIVALENT") 
				calibStrikeType = ARM_LocalCSOCalculator::EQUIVALENT;

			else if (calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
					 CalibMod[1] == "FRONTIER") 
				calibStrikeType = ARM_LocalCSOCalculator::FRONTIER;
			else
			{
				result.setMsg ("ARM_ERR: invalid calib strike type");
				return ARM_KO;
			}

			ARM_ModelType modelType = ARM_PricingModelType::HWM1F;

			/// set products to price
			std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(5,false);
			productsToPrice[0] = true;
			for (size_t i = 0; i < ProductFlags.size(); ++i)
			{
				productsToPrice[i] = ((ProductFlags[i] == "Y") || (ProductFlags[i] == "YES"));
			}
			vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(1,calibStrikeType);

            csoCalculator = new ARM_LocalCSOCalculator (
									ARM_Date(myStartDate),
									ARM_Date(myEndDate),
									CMS1Type,
									CMS2Type,
									cpnDaycount,
									cpnFreq,
									K_ADVANCE,
									*notionalProfile,
									*fixCpnProfile,
									*leverageProfile,
									*leverageProfile,
									*minCpnProfile,
									*maxCpnProfile,
									*strikeProfile,
									fundFreq,
									fundDaycount,
									*fundMarginProfile,
									ARM_FlatCurve(1.0), /// funding leverage (not interfaced)
									exerFreq,
									NotifDays,
									K_RCV,
									*feesProfile,
									false,
									0,
									productsToPrice,
									ModelDatas,
									*mktDataManager,
									keys,
									modelType,
									calibType,
									calibStrikeTypeVec
									);

		}

		//--------------------------------------------
		// SFRM
		//--------------------------------------------
		else
		{
			ARM_CSOCalculator::CalibrationMode sigmaCalibMode = (ARM_CSOCalculator::CalibrationMode) ARM::ARM_ArgConv_CSOCalibMode.GetNumber(CalibMod[0]); 
			ARM_CSOCalculator::CalibrationMode NDCalib = (ARM_CSOCalculator::CalibrationMode) ARM::ARM_ArgConv_CSOCalibMode.GetNumber(CalibMod[1]);

			bool MRSCalib;
			if (CalibMod[2] == "Y")
				MRSCalib = true;
			else if (CalibMod[2] == "N")
				MRSCalib = false;
			else
			{
				result.setMsg ("ARM_ERR: mean rev calib flag should be Y or N.");
				return ARM_KO;
			}

			bool thetaCalib;
			if (CalibMod[3] == "Y")
				thetaCalib = true;
			else if (CalibMod[3] == "N")
				thetaCalib = false;
			else
			{
				result.setMsg ("ARM_ERR: correl calib flag should be Y or N.");
				return ARM_KO;
			}

			std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(8,false);
			productsToPrice[0] = true;

			for (size_t i = 0; i < ProductFlags.size(); ++i)
			{
				productsToPrice[i] = ((ProductFlags[i] == "Y") || (ProductFlags[i] == "YES"));
			}

			/// Create the Callable SnowBall calculator
			csoCalculator = new ARM_CSOCalculator(
				ARM_Date(myStartDate),
				ARM_Date(myEndDate),
				CMS1Type,
				CMS2Type,
				cpnDaycount,
				cpnFreq,
				K_ADVANCE,
				*notionalProfile,
				*fixCpnProfile,
				*leverageProfile,
				*leverageProfile,
				*minCpnProfile,
				*maxCpnProfile,
				*strikeProfile,
				fundFreq,
				fundDaycount,
				*fundMarginProfile,
				exerFreq,
				NotifDays,
				K_RCV,
				*feesProfile,
				sigmaCalibMode,
				NDCalib,
				MRSCalib,
				thetaCalib,
				productsToPrice,
				ModelDatas,
				*mktDataManager,
				keys
				);

			
		}
		/// Free memory
		if(isNotional)
			delete notionalProfile;
		notionalProfile = NULL;
   
		if(isLeverage)
			delete leverageProfile;
		leverageProfile = NULL;

		if(isMinCpn)
			delete minCpnProfile;
		minCpnProfile = NULL;
   
		if(isMaxCpn)
			delete maxCpnProfile;
		maxCpnProfile = NULL;
   
		if(isFundMargin)
			delete fundMarginProfile;
		fundMarginProfile = NULL;

		if (strikeProfile)
			delete strikeProfile;
		strikeProfile = NULL;

		/// assign object
		if( !assignObject( csoCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
		
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isNotional)
            delete notionalProfile;
		notionalProfile = NULL;
       
        if(isLeverage)
            delete leverageProfile;
		leverageProfile = NULL;

        if(isMinCpn)
            delete minCpnProfile;
		minCpnProfile = NULL;
       
        if(isMaxCpn)
            delete maxCpnProfile;
		maxCpnProfile = NULL;
       
        if(isFundMargin)
            delete fundMarginProfile;
		fundMarginProfile = NULL;

		if (csoCalculator)
			delete csoCalculator;
		csoCalculator = NULL;

		if (strikeProfile)
			delete strikeProfile;
		strikeProfile = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_KO;
}

/////////////////////////////////////////////////////////////////////
/// CSO calculator creation : extended version
///	--> leverageLong != leverageShort possible
/// --> strike != 0.0 possible
///	--> will work only for HW + local vol + baske calib. SFRM impossible
/////////////////////////////////////////////////////////////////////
long ARMLOCAL_ExtendedCSOCalculator_Create (
        const double&		startDate,
        const double&		endDate,
        const long &		CMS1Type,
        const long &		CMS2Type,
        const long &		cpnFreq,
        const long &		cpnDaycount,
        const long &		cpnResetTiming,
        const long &		fundFreq,
        const long &		fundDaycount,
        const long &		fundResetTiming,
        const long &		exerFreq,
        const long &		NotifDays,
        const long &		payRec,
        const double &		cpnnotional,
        const long &		cpnnotionalId,
        const double &		minCpn,
        const long &		minCpnId,
        const double &		maxCpn,
        const long &		maxCpnId,
        const double &		leverageLong,
        const long &		leverageLongId,
        const double &		leverageShort,
        const long &		leverageShortId,
        const double &		strike,
        const long &		strikeId,
        const double &		fundMargin,
        const long &		fundMarginId,
		const double &		fundLeverage,
        const long &		fundLeverageId,
        const long &		fixCpnId,
        const long &		feesId,
		const bool &		switchFlag,
		const long &		fundingType,
        const vector< string >&	CalibMod,
        const vector< string >&	ProductFlags,
        const vector< double >&	ModelDatas,
        const long &		mktDataManagerId,
        const vector< string >&	keys,
        ARM_result&			result, 
        long				objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_GenCSOCalculator* csoCalculator = NULL;

	ARM_Curve* cpnnotionalProfile=NULL;
	ARM_Curve* minCpnProfile=NULL;
	ARM_Curve* maxCpnProfile=NULL;
	ARM_Curve* leverageLongProfile=NULL;
	ARM_Curve* leverageShortProfile=NULL;
	ARM_Curve* strikeProfile=NULL;
	ARM_Curve* fundMarginProfile=NULL;
	ARM_Curve* fundLeverageProfile=NULL;
	ARM_Curve* feesProfile=NULL;
	ARM_Curve* fixCpnProfile=NULL;

	ARM_AutoCleaner< ARM_Curve > HoldCNP(cpnnotionalProfile);
	ARM_AutoCleaner< ARM_Curve > HoldMCP(minCpnProfile);
	ARM_AutoCleaner< ARM_Curve > HoldMAXCP(maxCpnProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFMC(fundMarginProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFLC(fundLeverageProfile);
	ARM_AutoCleaner< ARM_Curve > HoldLLP(leverageLongProfile);
	ARM_AutoCleaner< ARM_Curve > HoldLSP(leverageShortProfile);
	ARM_AutoCleaner< ARM_Curve > HoldSP(strikeProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFEESC(feesProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFCP(fixCpnProfile);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable SpreadOption Calculator" );

		char myStartDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

        /// Convert curve Id to object if possible
		//notional
        if(cpnnotionalId != ARM_NULL_OBJECT)
        {
		    cpnnotionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnnotionalId));
		    if (!cpnnotionalProfile)
		    {
			    result.setMsg ("ARM_ERR: cpn notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            cpnnotionalProfile = new ARM_FlatCurve(cpnnotional);
			HoldCNP.setPtr(cpnnotionalProfile);
		}

		//minCpn
		if(minCpnId != ARM_NULL_OBJECT)
        {
		    minCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(minCpnId));

		    if (!minCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: minCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            minCpnProfile = new ARM_FlatCurve(minCpn);
			HoldMCP.setPtr(minCpnProfile);
		}

		//maxCpn
		if(maxCpnId != ARM_NULL_OBJECT)
        {
		    maxCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(maxCpnId));
		    if (!maxCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: maxCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            maxCpnProfile = new ARM_FlatCurve(maxCpn);
			HoldMAXCP.setPtr(maxCpnProfile);
		}

		//fundMargin
		if(fundMarginId != ARM_NULL_OBJECT)
        {
		    fundMarginProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundMarginId));

		    if (!fundMarginProfile)
		    {
			    result.setMsg ("ARM_ERR: fundMargin  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            fundMarginProfile = new ARM_FlatCurve(fundMargin);
			HoldFMC.setPtr(fundMarginProfile);
		}

		//fundLeverage
		if(fundLeverageId != ARM_NULL_OBJECT)
        {
		    fundLeverageProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundLeverageId));

		    if (!fundLeverageProfile)
		    {
			    result.setMsg ("ARM_ERR: fundLeverage  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            fundLeverageProfile = new ARM_FlatCurve(fundLeverage);
			HoldFLC.setPtr(fundLeverageProfile);
		}

		/// Leverage long
		if(leverageLongId != ARM_NULL_OBJECT)
        {
		    leverageLongProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageLongId));
		    if (!leverageLongProfile)
		    {
			    result.setMsg ("ARM_ERR: leverageLong  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            leverageLongProfile = new ARM_FlatCurve(leverageLong);
			HoldLLP.setPtr(leverageLongProfile);
		}

		/// leverageShort short
		if(leverageShortId != ARM_NULL_OBJECT)
        {
		    leverageShortProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageShortId));
		    if (!leverageShortProfile)
		    {
			    result.setMsg ("ARM_ERR: leverageShort  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            leverageShortProfile = new ARM_FlatCurve(leverageShort);
			HoldLSP.setPtr(leverageShortProfile);
		}
	
		/// strike
		if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));
		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            strikeProfile = new ARM_FlatCurve(strike);
			HoldSP.setPtr(strikeProfile);
		}

		/// fees
		feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));
		if (!feesProfile)
		{
			result.setMsg ("ARM_ERR: fees is not of a good type");
			return ARM_KO;
		}

		/// fixed coupon 
		fixCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fixCpnId));
		if (!fixCpnProfile)
		{
			result.setMsg ("ARM_ERR: fix Coupon is not of a good type");
			return ARM_KO;
		}
		
	    /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		//--------------------------------------------
		// HW + local vol 
		//--------------------------------------------
				
		/// set calibration type and model type

		// defaults
		ARM_LocalCSOCalculator::CalibrationType calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
		ARM_LocalCSOCalculator::CalibStrikeType calibStrikeType = ARM_LocalCSOCalculator::ATM;
		ARM_ModelType modelType = ARM_PricingModelType::HWM1F;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;

		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		if (CalibMod.size()>4) //switch to new calib hw2f
		{
			if (CalibMod.size()==6)
			{
				//MODEL
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(CalibMod[0]);

				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (CalibMod[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_LONG;
					else if (CalibMod[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_LocalCSOCalculator::DIAG_SPREAD_SHORT;
					else if (CalibMod[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_LocalCSOCalculator::SHORT_LONG_SPREAD;
					else if (CalibMod[1] == "DIAG")
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET")
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "DIAG_BASKET")
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else if (CalibMod[1] == "LONG")
						calibType = ARM_LocalCSOCalculator::LONG_CALIBRATION;
					else if (CalibMod[1] == "DIAG_LONG")
						calibType = ARM_LocalCSOCalculator::DIAG_LONG_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (CalibMod[1] == "DIAG") 
						calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET") 
						calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "DIAG_BASKET") 
						calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}

				if (CalibMod[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (CalibMod[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);
				else if (modelType==ARM_PricingModelType::HWM1F &&
						 calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 CalibMod[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (CalibMod[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else if (CalibMod[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ZERO);
				else if (CalibMod[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::CAP);
				else if (CalibMod[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (CalibMod[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}
					
			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[5]);
		}
		else
		{
			if (CalibMod[0] == "DIAG") 
				calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;

			else if (CalibMod[0] == "BASKET") 
				calibType = ARM_LocalCSOCalculator::BASKET_CALIBRATION;

			else if (CalibMod[0] == "DIAG_BASKET") 
				calibType = ARM_LocalCSOCalculator::DIAG_BASKET_CALIBRATION;

			else
			{
				result.setMsg ("ARM_ERR: invalid calib type");
				return ARM_KO;
			}

			if (CalibMod.size()>1)
			{
				if (CalibMod[1] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::ATM);

				else if (CalibMod[1] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::EQUIVALENT);

				else if (calibType == ARM_LocalCSOCalculator::DIAG_CALIBRATION &&
						 CalibMod[1] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_LocalCSOCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (CalibMod.size()>2)
				modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(CalibMod[2]);

			if (CalibMod.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[3]);
		}

		/// set products to price
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(6,false);
		productsToPrice[0] = true;

		for (size_t i = 0; i < ProductFlags.size(); ++i)
			productsToPrice[i] = ((ProductFlags[i] == "Y") || (ProductFlags[i] == "YES"));

		
		csoCalculator = new ARM_LocalCSOCalculator (
								ARM_Date(myStartDate),
								ARM_Date(myEndDate),
								CMS1Type,
								CMS2Type,
								cpnDaycount,
								cpnFreq,
								cpnResetTiming,
								*cpnnotionalProfile,
								*fixCpnProfile,
								*leverageLongProfile,
								*leverageShortProfile,
								*minCpnProfile,
								*maxCpnProfile,
								*strikeProfile,
								fundFreq,
								fundDaycount,
								*fundMarginProfile,
								*fundLeverageProfile,
								exerFreq,
								NotifDays,
								payRec,
								*feesProfile,
								switchFlag,
								fundingType,
								productsToPrice,
								ModelDatas,
								*mktDataManager,
								keys,
								modelType,
								calibType,
								calibStrikeTypeVec,
								vnsMethod
								);

		/// assign object
		if( !assignObject( csoCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
		
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_KO;
}

/////////////////////////////////////////////////////////////////////
/// CSO calculator creation : extended version
///	--> leverageLong != leverageShort possible
/// --> strike != 0.0 possible
///	--> will work only for HW + local vol + baske calib. SFRM impossible
/////////////////////////////////////////////////////////////////////
long ARMLOCAL_BasicCSOCalculator_Create(const double&		asOfDate,
                                        const double&		startDate,
                                        const double&		inFixEndDate,
                                        const double&		endDate,
                                        const long &		CMS1Type,
                                        const long &		CMS2Type,
                                        const long &		cpnFreq,
                                        const long &		cpnDaycount,
                                        const long &		cpnResetTiming,
                                        const long &		fundFreq,
                                        const long &		fundDaycount,
                                        const long &		fundResetTiming,
                                        const long &		exerFreq,
                                        const long &		NotifDays,
                                        const long &		payRec,
                                        const double &		cpnnotional,
                                        const long &		cpnnotionalId,
                                        const double &		minCpn,
                                        const long &		minCpnId,
                                        const double &		maxCpn,
                                        const long &		maxCpnId,
                                        const double &		leverageLong,
                                        const long &		leverageLongId,
                                        const double &		leverageShort,
                                        const long &		leverageShortId,
                                        const double &		strike,
                                        const long &		strikeId,
                                        const double &		fundNotional,
                                        const long &		fundNotionalId,
                                        const double &		fundMargin,
                                        const long &		fundMarginId,
		                                const double &		fundLeverage,
                                        const long &		fundLeverageId,
                                        const long &		fixCpnId,
                                        const long &		feesId,
                                        const CCString&		CpnCcy,
                                        const CCString&		FundCcy,
                                        const long&			nbNoCall,
                                        ARM_result&			result, 
                                        long				objId)
{
	/// input checks
	if (!GlobalPersistanceOk(result))
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg("");

    ARM_GenCSOCalculator* csoCalculator = NULL;

	ARM_Curve* cpnnotionalProfile   = NULL;
	ARM_Curve* minCpnProfile        = NULL;
	ARM_Curve* maxCpnProfile        = NULL;
	ARM_Curve* leverageLongProfile  = NULL;
	ARM_Curve* leverageShortProfile = NULL;
	ARM_Curve* strikeProfile        = NULL;
	ARM_Curve* fundNotionalProfile  = NULL;
	ARM_Curve* fundMarginProfile    = NULL;
	ARM_Curve* fundLeverageProfile  = NULL;
	ARM_Curve* feesProfile          = NULL;
	ARM_Curve* fixCpnProfile        = NULL;

	ARM_AutoCleaner< ARM_Curve > HoldCNP(cpnnotionalProfile);
	ARM_AutoCleaner< ARM_Curve > HoldMCP(minCpnProfile);
	ARM_AutoCleaner< ARM_Curve > HoldMAXCP(maxCpnProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFNP(fundNotionalProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFMC(fundMarginProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFLC(fundLeverageProfile);
	ARM_AutoCleaner< ARM_Curve > HoldLLP(leverageLongProfile);
	ARM_AutoCleaner< ARM_Curve > HoldLSP(leverageShortProfile);
	ARM_AutoCleaner< ARM_Curve > HoldSP(strikeProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFEESC(feesProfile);
	ARM_AutoCleaner< ARM_Curve > HoldFCP(fixCpnProfile);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService(ARM::ARM_USERNAME, "Callable SpreadOption Calculator" );

		char myAsOfDate[20];
		char myStartDate[20];
        char myFixEndDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(asOfDate, myAsOfDate);
		Local_XLDATE2ARMDATE(startDate, myStartDate);

        if ( inFixEndDate <= 0.0 )
        {
           // No initial FIX, so use the start date as fix End Date
           Local_XLDATE2ARMDATE(startDate, myFixEndDate);
        }
        else
        {
           Local_XLDATE2ARMDATE(inFixEndDate, myFixEndDate);
        }

		Local_XLDATE2ARMDATE(endDate, myEndDate);

        /// Convert curve Id to object if possible
		// notional
        if(cpnnotionalId != ARM_NULL_OBJECT)
        {
		    cpnnotionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnnotionalId));
		    if (!cpnnotionalProfile)
		    {
			    result.setMsg ("ARM_ERR: cpn notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            cpnnotionalProfile = new ARM_FlatCurve(cpnnotional);
			HoldCNP.setPtr(cpnnotionalProfile);
		}

		//minCpn
		if(minCpnId != ARM_NULL_OBJECT)
        {
		    minCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(minCpnId));

		    if (!minCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: minCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            minCpnProfile = new ARM_FlatCurve(minCpn);
			HoldMCP.setPtr(minCpnProfile);
		}

		//maxCpn
		if(maxCpnId != ARM_NULL_OBJECT)
        {
		    maxCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(maxCpnId));
		    if (!maxCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: maxCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            maxCpnProfile = new ARM_FlatCurve(maxCpn);
			HoldMAXCP.setPtr(maxCpnProfile);
		}

		// funding notional
        if(fundNotionalId != ARM_NULL_OBJECT)
        {
		    fundNotionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundNotionalId));
		    if (!fundNotionalProfile)
		    {
			    result.setMsg ("ARM_ERR: funding notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            fundNotionalProfile = new ARM_FlatCurve(fundNotional);
			HoldFNP.setPtr(fundNotionalProfile);
		}

		//fundMargin
		if(fundMarginId != ARM_NULL_OBJECT)
        {
		    fundMarginProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundMarginId));

		    if (!fundMarginProfile)
		    {
			    result.setMsg ("ARM_ERR: fundMargin  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            fundMarginProfile = new ARM_FlatCurve(fundMargin);
			HoldFMC.setPtr(fundMarginProfile);
		}

		//fundLeverage
		if(fundLeverageId != ARM_NULL_OBJECT)
        {
		    fundLeverageProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundLeverageId));

		    if (!fundLeverageProfile)
		    {
			    result.setMsg ("ARM_ERR: fundLeverage  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            fundLeverageProfile = new ARM_FlatCurve(fundLeverage);
			HoldFLC.setPtr(fundLeverageProfile);
		}

		/// Leverage long
		if(leverageLongId != ARM_NULL_OBJECT)
        {
		    leverageLongProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageLongId));
		    if (!leverageLongProfile)
		    {
			    result.setMsg ("ARM_ERR: leverageLong  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            leverageLongProfile = new ARM_FlatCurve(leverageLong);
			HoldLLP.setPtr(leverageLongProfile);
		}

		/// leverageShort short
		if(leverageShortId != ARM_NULL_OBJECT)
        {
		    leverageShortProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageShortId));
		    if (!leverageShortProfile)
		    {
			    result.setMsg ("ARM_ERR: leverageShort  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            leverageShortProfile = new ARM_FlatCurve(leverageShort);
			HoldLSP.setPtr(leverageShortProfile);
		}
	
		/// strike
		if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));
		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
		{
            strikeProfile = new ARM_FlatCurve(strike);
			HoldSP.setPtr(strikeProfile);
		}

		/// fees
		feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));
		if (!feesProfile)
		{
			result.setMsg ("ARM_ERR: fees is not of a good type");
			return ARM_KO;
		}

		/// fixed coupon 
		fixCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fixCpnId));
		if (!fixCpnProfile)
		{
			result.setMsg ("ARM_ERR: fix Coupon is not of a good type");
			return ARM_KO;
		}
		const ARM_Currency currency;
        ARM_Currency theCpnCcy;
        ARM_Currency theFundCcy;

        if ( CpnCcy == "DEFAULT" )
        {
           ARM_Currency* ccy = ARM_DEFAULT_CURRENCY;

           theCpnCcy  = *ccy;
           theFundCcy = *ccy;
        }
        else
        {
           theCpnCcy = ARM_Currency((const char *) CpnCcy);

           if ( FundCcy == "DEFAULT" )
           {
              theFundCcy = theCpnCcy;
           }
           else
           {
              theFundCcy = ARM_Currency((const char *) FundCcy);
           }
        }

        ARM_Date fixEndDate = ARM_Date(myFixEndDate); // In the case of Initial Fix

		csoCalculator = new ARM_LocalCSOCalculator(ARM_Date(myAsOfDate),
								                   ARM_Date(myStartDate),
                                                   fixEndDate,
								                   ARM_Date(myEndDate),
								                   theCpnCcy,
								                   theFundCcy,
								                   CMS1Type,
								                   CMS2Type,
								                   cpnDaycount,
								                   cpnFreq,
								                   cpnResetTiming,
								                   *cpnnotionalProfile,
								                   *fixCpnProfile,
								                   *leverageLongProfile,
								                   *leverageShortProfile,
								                   *minCpnProfile,
								                   *maxCpnProfile,
								                   *strikeProfile,
								                   fundFreq,
								                   fundDaycount,
												   *fundNotionalProfile,
								                   *fundMarginProfile,
								                   *fundLeverageProfile,
								                   exerFreq,
								                   NotifDays,
								                   payRec,
								                   *feesProfile,
												   false,
												   0,
												   nbNoCall);

		/// assign object
		if ( !assignObject(csoCalculator, result, objId) )
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
		x.DebugPrint();
		ARM_RESULT();
	}
	
    return(ARM_KO);
}



/////////////////////////////////////////////////////////////////////
/// CSO calculator creation : extended version
///	--> leverageLong != leverageShort possible
/// --> strike != 0.0 possible
///	--> will work only for HW + local vol + baske calib. SFRM impossible
/////////////////////////////////////////////////////////////////////
long ARMLOCAL_BasisCSOCalculator_Create (
        const double&		startDate,
		const double&		fixEndDate,
        const double&		endDate,
		const string&       cpnCcy,
		const string&       fundCcy,
        const long &		CMS1Type,
        const long &		CMS2Type,
        const long &		cpnFreq,
        const long &		cpnDaycount,
        const long &		cpnResetTiming,
		const long &        stubRule,
        const string&       cpnResetCal,
        const string&       cpnPayCal,
		const long &		cpnnotionalId,
		const long &		strikeId,
		const long &		leverageLongId,
        const long &		leverageShortId,
        const long &		minCpnId,
        const long &		maxCpnId,
        const long &		fundFreq,
        const long &		fundDaycount,
        const long &		fundResetTiming,
		const long &		fundingnotionalId,
		const long &		fundMarginId,
        const long &		exerFreq,
        const long &		NotifDays,
        const long &		payRec,
		const long &        nbNCall,
        const long &		feesId,
        const vector< string >&	CalibDatas,
        const vector< string >&	ProductFlags,
        const vector< double >&	ModelDatas,
        const long &		mktDataManagerId,
        ARM_result&			result, 
        long				objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_GenCSOCalculator* csoCalculator = NULL;
	
	ARM_Curve* cpnnotionalProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldCNP(cpnnotionalProfile);

    ARM_Curve* fundingnominalprofile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFNP(fundingnominalprofile);

	ARM_Curve* minCpnProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMCP(minCpnProfile);

	ARM_Curve* maxCpnProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMAXCP(maxCpnProfile);

	ARM_Curve* leverageLongProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLLP(leverageLongProfile);

	ARM_Curve* leverageShortProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLSP(leverageShortProfile);
	
	ARM_Curve* strikeProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldSP(strikeProfile);

	ARM_Curve* fundMarginProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFMC(fundMarginProfile);

	ARM_Curve* feesProfile=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFEESC(feesProfile);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable SpreadOption Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		char myFixEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		Local_XLDATE2ARMDATE(fixEndDate, myFixEndDate);

        /// Convert curve Id to object if possible
		//notionals
		if( !GetObjectFromIdWithDynamicCastCheck( &cpnnotionalProfile, cpnnotionalId  ) )
		{
			result.setMsg ("ARM_ERR: cpn notional is not of a good type");
			return ARM_KO;
		}
		if( !GetObjectFromIdWithDynamicCastCheck( &fundingnominalprofile, fundingnotionalId  ) )
		{
			result.setMsg ("ARM_ERR: funding nominal is not of a good type");
			return ARM_KO;
		}

		//minCpn
		if( !GetObjectFromIdWithDynamicCastCheck( &minCpnProfile, minCpnId  ) )
		{
			result.setMsg ("ARM_ERR: minCpn  is not of a good type");
			return ARM_KO;
		}

		//maxCpn
		if( !GetObjectFromIdWithDynamicCastCheck( &maxCpnProfile, maxCpnId  ) )
		{
			result.setMsg ("ARM_ERR: maxCpn  is not of a good type");
			return ARM_KO;
		}

		//fundMargin
		if( !GetObjectFromIdWithDynamicCastCheck( &fundMarginProfile, fundMarginId  ) )
		{
			result.setMsg ("ARM_ERR: fundMargin  is not of a good type");
			return ARM_KO;
		}

		/// Leverage long
		if( !GetObjectFromIdWithDynamicCastCheck( &leverageLongProfile, leverageLongId  ) )
		{
			result.setMsg ("ARM_ERR: leverageLong  is not of a good type");
			return ARM_KO;
		}


		/// leverageShort short
		if( !GetObjectFromIdWithDynamicCastCheck( &leverageShortProfile, leverageShortId  ) )
		{
			result.setMsg ("ARM_ERR: leverageShort  is not of a good type");
			return ARM_KO;
		}

		/// strike
		if( !GetObjectFromIdWithDynamicCastCheck( &strikeProfile, strikeId  ) )
		{
			result.setMsg ("ARM_ERR: strike  is not of a good type");
			return ARM_KO;
		}

		/// fees
		if( !GetObjectFromIdWithDynamicCastCheck( &feesProfile, feesId  ) )
		{
			result.setMsg ("ARM_ERR: fees is not of a good type");
			return ARM_KO;
		}
		
	    /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager= NULL;
		if( !GetObjectFromIdWithDynamicCastCheck( &mktDataManager, mktDataManagerId  ) )
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}
		
	

		//--------------------------------------------
		// HW + local vol 
		//--------------------------------------------
				
		/// set calibration type and model type

		// defaults
		ARM_LocalCSOCalculator::CalibrationType calibType = ARM_LocalCSOCalculator::DIAG_CALIBRATION;
		ARM_LocalCSOCalculator::CalibStrikeType calibStrike1Type = ARM_LocalCSOCalculator::ATM;
		ARM_LocalCSOCalculator::CalibStrikeType calibStrike2Type = ARM_LocalCSOCalculator::ATM;
		ARM_LocalCSOCalculator::CalibStrikeType calibStrike3Type = ARM_LocalCSOCalculator::ATM;
		ARM_ModelType modelType = ARM_PricingModelType::HWM1F;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;
		double moyenessLevel = 1.0;

		/// Which model?
		if (CalibDatas.size()>2)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(CalibDatas[2]);

		if (CalibDatas.size()>1)
			calibStrike1Type = (ARM_LocalCSOCalculator::CalibStrikeType)ARM_ArgConv_CSOStrike1Type.GetNumber(CalibDatas[1]);
		
		if (CalibDatas.size()>3)
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibDatas[3]);
		
		vector<ARM_LocalCSOCalculator::CalibStrikeType> calibStrikeTypeVec(1,calibStrike1Type);

		if (modelType == ARM_PricingModelType::HWM1F){
			calibType = (ARM_LocalCSOCalculator::CalibrationType) ARM_ArgConv_CSOCalibModeHWM1F.GetNumber(CalibDatas[0]);		

		/*if (CalibDatas.size()>4)
			moyenessLevel = atof(CalibDatas[4].c_str()); */
		}
		else if (modelType == ARM_PricingModelType::HWM2F){
			calibType = (ARM_LocalCSOCalculator::CalibrationType) ARM_ArgConv_CSOCalibModeHWM2F.GetNumber(CalibDatas[0]);

			if (CalibDatas.size()>4)
				calibStrike2Type = (ARM_LocalCSOCalculator::CalibStrikeType)ARM_ArgConv_CSOStrike2Type.GetNumber(CalibDatas[4]);
			calibStrikeTypeVec.push_back(calibStrike2Type);

			if (CalibDatas.size()>5)
				calibStrike3Type = (ARM_LocalCSOCalculator::CalibStrikeType)ARM_ArgConv_CSOStrike3Type.GetNumber(CalibDatas[5]);
			calibStrikeTypeVec.push_back(calibStrike3Type);
		}


		/// set products to price
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(6,false);
		productsToPrice[0] = true;

		for (size_t i = 0; i < ProductFlags.size(); ++i)
			productsToPrice[i] = ((ProductFlags[i] == "Y") || (ProductFlags[i] == "YES"));

		ARM_Currency CpnCurrency(cpnCcy.c_str());
		ARM_Currency FundCurrency(fundCcy.c_str());

		
		csoCalculator = new ARM_LocalCSOCalculator (ARM_Date(myStartDate),
								ARM_Date(myFixEndDate),
								ARM_Date(myEndDate),
								CpnCurrency,
								FundCurrency,
								CMS1Type,
								CMS2Type,
								cpnDaycount,
								cpnFreq,
								cpnResetTiming,
								stubRule,
                                cpnResetCal,
                                cpnPayCal,
								*cpnnotionalProfile,
								*leverageLongProfile,
								*leverageShortProfile,
								*minCpnProfile,
								*maxCpnProfile,
								*strikeProfile,
								fundFreq,
								fundDaycount,
								*fundingnominalprofile,
								*fundMarginProfile,
								exerFreq,
								NotifDays,
								payRec,
								false,
								0,
								nbNCall,
								*feesProfile,
								productsToPrice,
								ModelDatas,
								*mktDataManager,
								modelType,
								calibType,
								calibStrikeTypeVec,
								vnsMethod,
								moyenessLevel
								);

		/// assign object

		return (assignObject( csoCalculator, result, objId ) )? ARM_OK: ARM_KO;
		
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_KO;
}




/////////////////////////////////////////////////////////////////////
/// CSO calculator get datas
/////////////////////////////////////////////////////////////////////
long ARMLOCAL_CSO_Get(
        const long& CSOId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CSOCalculator* CSOCalculator;
    ARM_LocalCSOCalculator* localCSOCalculator;
    ARM_Object* object=NULL;

	try
	{
		CSOCalculator = dynamic_cast<ARM_CSOCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(CSOId));

		if (!CSOCalculator)
		{
			localCSOCalculator = dynamic_cast<ARM_LocalCSOCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(CSOId));
			if (!localCSOCalculator)
			{
				result.setMsg ("ARM_ERR: Callable SpreadOption Calculator is not of a good type");
				return ARM_KO;
			}
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));

        if	(CSOCalculator &&
			 (
			  (typeToGet == GC_ACCESS_OSW_PORT) || (typeToGet == GC_ACCESS_SO_PORT) || 
			  (typeToGet == GC_ACCESS_CMSLONG_PORT) || (typeToGet == GC_ACCESS_CMSSHORT_PORT)
			 )
			)
		{
			if(typeToGet == GC_ACCESS_OSW_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(CSOCalculator->GetOSWPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_SO_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(CSOCalculator->GetSOPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CMSLONG_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(CSOCalculator->GetCMSLONGPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CMSSHORT_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(CSOCalculator->GetCMSSHORTPortfolio())).Clone();
			}
			
			/// Assign the object in ARM cache
			if(!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}

		if	(localCSOCalculator &&
			 (
			  (typeToGet == GC_ACCESS_OSW_PORT) || 
			  (typeToGet == GC_ACCESS_OSW_SECOND_PORT) || 
			  (typeToGet == GC_ACCESS_SO_PORT)  || 
			  (typeToGet == GC_ACCESS_SO_SECOND_PORT) ||
              (typeToGet == GC_ACCESS_SO_CSTMANAGER) ||
			  (typeToGet == GC_ACCESS_CSO_PORT) ||
			  (typeToGet == GC_ACCESS_CMSLONG_PORT) ||
			  (typeToGet == GC_ACCESS_CMSSHORT_PORT)
			 )
			)
		{
			if(typeToGet == GC_ACCESS_OSW_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetSwaptionPF())).Clone();
			}
			else if (typeToGet == GC_ACCESS_OSW_SECOND_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetFwdStartSwaptionPF())).Clone();
			}
			else if (typeToGet == GC_ACCESS_SO_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetSpreadOptionFloorPF())).Clone();
			}
			else if (typeToGet == GC_ACCESS_SO_SECOND_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetSpreadOptionCapPF())).Clone();
			}
            else if (typeToGet == GC_ACCESS_SO_CSTMANAGER)
			{
				object = const_cast< ARM_CstManager* >(&*(localCSOCalculator->CreateCstManager()))->Clone();
			}
			else if (typeToGet == GC_ACCESS_CSO_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetSOPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CMSLONG_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetCMSLONGPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CMSSHORT_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(localCSOCalculator->GetCMSSHORTPortfolio())).Clone();
			}
			
		
			/// Assign the object in ARM cache
			if(!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}
				
		
		// call to generic calculator
		return ARMLOCAL_GC_Get(CSOId,getType,result,objId);
		
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}
///////////////////////////////////////////////
//// Function to set data to a CSO Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_CSO_Set(
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
    ARM_LocalCSOCalculator* newCSO=NULL;
    ARM_Object* object=NULL;
    string updateMsg(" CSO updated by the ");

	try
	{
		GetObjectFromIdWithDynamicCastCheck( &localCSO, calculatorId  ); 
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

		ARM_ModelParam* param = NULL;
        ARM_StdPortfolio* port = NULL;
		bool isUpdateOswStrike = true;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
			/// Just update the initial CSO or Clone the initial CSO then set
			newCSO= isUpdated ? localCSO: dynamic_cast<ARM_LocalCSOCalculator *>(localCSO->Clone());

            string typeToSet(stringGetUpper(setType));
			if(typeToSet == GC_ACCESS_OSW_PORT)
            {
                /// Default
                newCSO->SetOSWPortfolio(*port);
                updateMsg += "OSW portfolio";
                isUpdateOswStrike = false;
            }
            else if(typeToSet == GC_ACCESS_SHORT_TERM_PORT)
            {
                updateMsg += "BZZZZZZZZZZZZZZZ";
            }
		    else
            {
                /// Default
                updateMsg += "BZZZZZZZZZZZZZZZ";
            }

            /// Keep calibration data consistency
			newCSO->UpdateCalibrationAndTimeIt(isUpdateOswStrike);

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CSO in ARM cache
		    if(!assignObject( newCSO, result, objId ) )
            {
                delete newCSO;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else if (param = dynamic_cast< ARM_ModelParam* >(object))
		{
			newCSO = dynamic_cast<ARM_LocalCSOCalculator *>(localCSO->Clone());
			updateMsg += "BZZZZZZZZZZZZZZZZZZZZZZ";
			
			/// Assign the new CSO in ARM cache
		    if(!assignObject( newCSO, result, objId ) )
            {
                delete newCSO;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
		}
		else
            return ARMLOCAL_GC_Set(calculatorId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newCSO;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////////////////
/// CRA FRM Markov Tree Calculator (first version)              
///////////////////////////////////////////////////////////
extern long ARMLOCAL_CRA_GetModel(
        long craId,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRAFRMMarkovTree* craCalculator;
    ARM_Object* object=NULL;

	try
	{
		craCalculator = dynamic_cast<ARM_CRAFRMMarkovTree *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));

		if (!craCalculator)
		{
			result.setMsg ("ARM_ERR: CRA Markov Tree Calculator is not of a good type");
			return ARM_KO;
		}

      	ARM_PricingModel* pmod =  static_cast<ARM_PricingModel *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetPricingModel()->Clone());
	    object = (ARM_Object*)pmod;
		
		/// Assign the object in ARM cache
		if(!assignObject(object, result, objId))
		{
			delete object;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_CRA_GetCalibData(
        long craId,
		string calibOrPortfolio,
		string dataType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRAFRMMarkovTree* craCalculator;
    ARM_Object* object=NULL;
	
	ARM_CalibMethod* calib		= NULL;
	ARM_StdPortfolio* portfolio	= NULL;

	try
	{
		// CRA CALCULATOR
		craCalculator = dynamic_cast<ARM_CRAFRMMarkovTree *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));

		if (!craCalculator)
		{
			result.setMsg ("ARM_ERR: CRA Markov Tree Calculator is not of a good type");
			return ARM_KO;
		}

		// CALIB METHOD
		if (calibOrPortfolio == "CALIB")
		{
			if (dataType == "SIGMA")
			{
				calib =  static_cast<ARM_CalibMethod *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetCalibMethodSigma()->Clone());
			}
			else if (dataType == "BETA") 
			{
				calib =  static_cast<ARM_CalibMethod *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetCalibMethodBeta()->Clone());
			}
			else if (dataType == "MRS")
			{
				calib =  static_cast<ARM_CalibMethod *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetCalibMethodMR()->Clone());
			}
			
			object = (ARM_Object*)calib;
	
		}
		// PORTFOLIO
		else if (calibOrPortfolio == "PORTFOLIO")
		{
			if (dataType == "SIGMA")
			{
				portfolio =  static_cast<ARM_StdPortfolio *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetSigmaPortfolio()->Clone());
			}
			else if (dataType == "BETA") 
			{
				portfolio =  static_cast<ARM_StdPortfolio *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetBetaPortfolio()->Clone());
			}
			else if (dataType == "MRS")
			{
				portfolio =  static_cast<ARM_StdPortfolio *> (((ARM_CalibratorSFRM *)(craCalculator->GetCalibratorSFRM()))->GetMRSPortfolio()->Clone());
			}

			object = (ARM_Object*)portfolio;
		}

		/// Assign the object in ARM cache
		if(!assignObject(object, result, objId))
		{
			delete object;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_CRA_GetUnderlying(
        long craId,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRAFRMMarkovTree* craCalculator;
    ARM_Object* object=NULL;
	
	ARM_Swap*			exoticSwap		= NULL;
	ARM_CorridorLeg*	corridorLeg		= NULL;
	ARM_SwapLeg*		fundingLeg      = NULL;
	
	try
	{
		// CRA CALCULATOR
		craCalculator = dynamic_cast<ARM_CRAFRMMarkovTree *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));

		if (!craCalculator)
		{
			result.setMsg ("ARM_ERR: CRA Markov Tree Calculator is not of a good type");
			return ARM_KO;
		}

		// FUNDING / CORRIDOR LEGS
		ARM_OptionPortfolio* optionPortfolio = (ARM_OptionPortfolio*)((craCalculator->GetCalibratorSFRM())->GetSecurity());
		
		corridorLeg = optionPortfolio->GetCorridorLeg();
		fundingLeg  = optionPortfolio->GetLiborLeg();
		
		// EXOTIC SWAP
		exoticSwap = new ARM_Swap((ARM_SwapLeg*)corridorLeg, fundingLeg, -K_HUGE_DOUBLE);
		object = (ARM_Object*)exoticSwap;

		delete corridorLeg;
		corridorLeg = NULL;
		delete fundingLeg;
		fundingLeg = NULL;

		/// Assign the object in ARM cache
		if(!assignObject(object, result, objId))
		{
			delete object;
			return ARM_KO;
		}
		else
		{
			return ARM_OK;
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////////////////
/// BERMUDA swaption
///////////////////////////////////////////////////////////
extern long ARMLOCAL_BERMUDASWAPTIONCalculator_Create(
		ARM_Currency ccy,
	    double startDate,
		double endDate,
		double notional,
		long notionalId,
		double strike,
		long strikeId,
		double fees,
		long feesId,
		double spread,
		long spreadId,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		double firstCallDate,
		double lastCallDate,
		int fixFreq,
		int fixBasis,
		string fixPayCal,
		int fixAdjRule,
		int fixRule,
		int fixPayGap,
		bool isZc,
		int varFreq,
		int varBasis,
		string varResetCal,
		string varPayCal,
		string varIndexTerm,
		int varAdjRule,
		int varRule,
		int varResetGap,
		int varPayGap,
		int stubRule,
		int genSecType,
		vector<int>* controlVariates,
		vector<double>* controlPrices,
		vector<string> mdmKeys,
		long mktDataManagerId,
		int modelType,
        vector<string> modelParams,
		vector<string> calibFlags,
		int numMethodType,
		int amcIter,
		int mcIter,
		int maxBucketSize,
		string genType1,
		string genType2,
		string pathOrder,
		string pathScheme,
		int firstNbDims,
		int treeSteps,
		vector<int> portfolioMode,
		bool boundaryFlag,
		bool approxMarginFlag,
		bool freezeBetasFlag,
		bool calculateProbaFlag,
		ARM_result&	result, 
        long        objId )
		  
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_BermudaSwaptionCalculator* bsCalculator = NULL;

	ARM_ReferenceValue* strikeProfile = NULL;
    bool isStrike=false;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* feesProfile = NULL;
    bool isFees=false;

	ARM_ReferenceValue* spreadProfile = NULL;
    bool isSpread=false;

	vector<ARM_ReferenceValue*>* modelParamsRef;

	int i;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Bermuda Swaption Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		char myFirstCallDate[20];
		char myLastCallDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		Local_XLDATE2ARMDATE(firstCallDate, myFirstCallDate);
		Local_XLDATE2ARMDATE(lastCallDate, myLastCallDate);
	
		//Strike Curve
        if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            strikeProfile = new ARM_ReferenceValue(strike);
            isStrike=true;
        }

		//Notional Curve
        if(notionalId != ARM_NULL_OBJECT)
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

		//Fees Curve
        if(feesId != ARM_NULL_OBJECT)
        {
		    feesProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));

		    if (!feesProfile)
		    {
			    result.setMsg ("ARM_ERR: fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            feesProfile = new ARM_ReferenceValue(fees);
            isFees=true;
        }
		
		//Spread Curve
        if(spreadId != ARM_NULL_OBJECT)
        {
		    spreadProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(spreadId));

		    if (!spreadProfile)
		    {
			    result.setMsg ("ARM_ERR: spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            spreadProfile = new ARM_ReferenceValue(spread);
            isSpread=true;
        }

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		ARM_Date adjLastCallDate = ARM_Date(myLastCallDate);

		//Calib Flags
		bool mrsCalibFlag = false;
		bool betaCalibFlag = false;
		if (calibFlags[0] == "Y")
		{
			mrsCalibFlag = true;
		}
		if (calibFlags[1] == "Y")
		{
			betaCalibFlag = true;
		}

		// model params : convert into Reference Values
		modelParamsRef = new vector<ARM_ReferenceValue*>(7);
		for (i=0; i<modelParams.size(); i++)
		{
			CCString objId(modelParams[i].c_str());
			if (LocalGetNumObjectId(objId) == ARM_KO)
			{
				double value = 0.0;
				if (modelParams[i] == "DEFAULT")
				{
					if (i == 5)
						value = -0.1;
					if (i == 6)
						value = 0.1;
				}
				else
					value = atof(modelParams[i].c_str());
				
				(*modelParamsRef)[i] = new ARM_ReferenceValue(value);
			}
			else
			{
				(*modelParamsRef)[i] = dynamic_cast<ARM_ReferenceValue*>(LOCAL_PERSISTENT_OBJECTS->GetObject(LocalGetNumObjectId(objId))->Clone());
			}
		}

        /// Create the Bermuda Swaption calculator
        bsCalculator = new 	ARM_BermudaSwaptionCalculator(ccy,
								  ARM_Date(myStartDate),
								  ARM_Date(myEndDate),
								  *notionalProfile,
								  *strikeProfile,
								  payReceive,
								  stubRule,
								  callFreq,
								  callNotice,
								  callCal,
								  ARM_Date(myFirstCallDate),
								  adjLastCallDate,
								  *feesProfile,
								  fixFreq,
								  fixBasis,
								  fixPayCal,
							      fixAdjRule,
								  fixRule,
								  fixPayGap,
								  isZc,
								  varFreq,
								  varBasis,
								  varResetCal,
								  varPayCal,
								  varIndexTerm,
								  varAdjRule,
								  varRule,
								  varResetGap,
								  varPayGap,
								  *spreadProfile,
								  genSecType,
								  controlVariates,
								  controlPrices,
								  mdmKeys,
								  *mktDataManager,
								  modelType,
								  modelParamsRef,
								  mrsCalibFlag,
								  betaCalibFlag,
								  numMethodType,
								  amcIter,
								  mcIter,
								  maxBucketSize,
								  genType1,
								  genType2,
								  pathScheme,
								  pathOrder,
								  firstNbDims,
								  treeSteps,
								  portfolioMode,
								  boundaryFlag,
								  approxMarginFlag,
								  freezeBetasFlag,
								  calculateProbaFlag);
		 		
        /// Free memory
        if(isStrike)
            delete strikeProfile;
		strikeProfile = NULL;
        if(isNotional)
            delete notionalProfile;
		notionalProfile = NULL;
        if(isFees)
            delete feesProfile;
		feesProfile = NULL;
		if(isSpread)
            delete spreadProfile;
		spreadProfile = NULL;

		for (i = 0; i < modelParamsRef->size(); i++)
			delete (*modelParamsRef)[i];

		delete modelParamsRef;

		// assign object
		if( !assignObject( bsCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		if(isStrike)
            delete strikeProfile;
		strikeProfile = NULL;
        if(isNotional)
            delete notionalProfile;
		notionalProfile = NULL;
        if(isFees)
            delete feesProfile;
		feesProfile = NULL;
		if(isSpread)
            delete spreadProfile;
		spreadProfile = NULL;

		for (int i = 0; i < modelParamsRef->size(); i++)
			delete (*modelParamsRef)[i];

		delete modelParamsRef;

		//delete bsCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}

}

extern long ARMLOCAL_BERMUDASWAPTION_Get(
        long bsId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_BermudaSwaptionCalculator* bsCalculator;
    ARM_Object* object=NULL;

	try
	{
		bsCalculator = dynamic_cast<ARM_BermudaSwaptionCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bsId));

		if (!bsCalculator)
		{
			result.setMsg ("ARM_ERR: Bermuda Swaption Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if ((typeToGet == GC_ACCESS_OSW_PORT) || (typeToGet == GC_ACCESS_SHORT_TERM_PORT) || (typeToGet == GC_ACCESS_STD_PORT) ||
			(typeToGet == GC_ACCESS_CALIB) || (typeToGet == GC_ACCESS_EXO_SWAPTION) || (typeToGet == GC_ACCESS_PRECALIB)||
			(typeToGet == GC_ACCESS_EXOSWAP))
		{
			if(typeToGet == GC_ACCESS_OSW_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(bsCalculator->GetOSWPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_SHORT_TERM_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(bsCalculator->GetSTMPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_STD_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(bsCalculator->GetStdPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CALIB)
			{
				object = const_cast< ARM_CalibMethod& >(*(bsCalculator->GetOSWCalibMethod())).Clone();
			}
			else if (typeToGet == GC_ACCESS_PRECALIB)
			{
				object = const_cast< ARM_CalibMethod& >(*(bsCalculator->GetPreCalibMethod())).Clone();
			}
			else if (typeToGet == GC_ACCESS_EXO_SWAPTION)
			{
				object = const_cast< ARM_Swaption& >(*(bsCalculator->GetExoSwaption())).Clone();
			}
			else if (typeToGet == GC_ACCESS_EXOSWAP)
			{
				object = const_cast< ARM_Swap& >(*(bsCalculator->GetUnderSwap())).Clone();
			}

			/// Assign the object in ARM cache
			if(!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}
		else
		{
			return ARMLOCAL_GC_Get(bsId,getType,result,objId);
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a Bermuda Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_BERMUDASWAPTION_Set(
        long bsId,
        long dataId,
        const string& setType,
		const vector< string >& keys,
        bool isUpdated,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_BermudaSwaptionCalculator* oldBS=NULL;
    ARM_BermudaSwaptionCalculator* newBS=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");

	try
	{
		oldBS = dynamic_cast<ARM_BermudaSwaptionCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bsId));
		if (!oldBS)
		{
			result.setMsg ("ARM_ERR: Bermuda Swaption Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
            if(isUpdated)
                /// Just update the initial Bermuda Swaption
                newBS = oldBS;
            else
                /// Clone the initial Bermuda Swaption then set
                newBS = dynamic_cast<ARM_BermudaSwaptionCalculator *>(oldBS->Clone());

            string typeToSet(stringGetUpper(setType));
            if(typeToSet == GC_ACCESS_OSW_PORT)
            {
                newBS->SetOSWPortfolio(*port);
                updateMsg += "OSW portfolio";
            }
            else if(typeToSet == GC_ACCESS_SHORT_TERM_PORT)
            {
                newBS->SetSTMPortfolio(*port);
                updateMsg += "STM vanillas portfolio";
            }
		
            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRF in ARM cache
		    if(!assignObject( newBS, result, objId ) )
            {
                delete newBS;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else
            return ARMLOCAL_GC_Set(bsId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newBS;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_ARM_BermudaRootMrs(
        long		bsId,
        double		targetPrice,
		double		fTolerance,
		int			maxIter,
        ARM_result&	result)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_BermudaSwaptionCalculator* clonedBs = NULL;

	try
	{
		clonedBs = dynamic_cast<ARM_BermudaSwaptionCalculator *>((LOCAL_PERSISTENT_OBJECTS->GetObject(bsId))->Clone());
	
		if (!clonedBs)
		{
			result.setMsg ("ARM_ERR: Bermuda Swaption Calculator is not of a good type");
			return ARM_KO;
		}
        
		double rootMrs = clonedBs->MeanReversionRoot(targetPrice, fTolerance, maxIter);
	
		/// save the result!
		result.setDouble(rootMrs);
	
		if (clonedBs)
		{
			delete clonedBs;
			clonedBs = NULL;
		}

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		if (clonedBs)
		{
			delete clonedBs;
			clonedBs = NULL;
		}
        x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_CRALocalCalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		int fundFreq,
		int fundDayCount,
		int cpnPayFreq,
		string cpnResetCal,
		string cpnPayCal,
		int boostedIndexType,
		string boostedVarTerm,
		int boostedResetGap,
		int boostedResetTiming,
		int boostedDayCount,
		int boostedAdjRule,
		int boostedRule,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnResetGap,
		int refIndexType,
		string  refTerm,
		int refDayCount,
		double refCoeff,
		double notional,
		long notionalId,
		long callFeesId,
		double fundSpread,
		long fundSpreadId,
		double cpnSpread,
		long cpnSpreadId,
		double boostedFix,
		long boostedFixId,
		double bDown,
		long bDownId,
		double bUp,
		long bUpId,
		bool localModel,
		int localModelType,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		int localResetFreq,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool isStdCalib,
		ARM_result&	result, 
        long objId)
		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CRALocalCalculator* cfcraCalculator = NULL;
	ARM_CRACalculator* craCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* callFeesProfile = NULL;
    bool isCallFees=false;
	
	ARM_ReferenceValue* fundSpreadProfile = NULL;
    bool isFundSpread=false;
	
	ARM_ReferenceValue* cpnSpreadProfile = NULL;
    bool isCpnSpread=false;

	ARM_ReferenceValue* boostedFixProfile = NULL;
    bool isBoostedFix=false;
	
	ARM_ReferenceValue* bDownProfile = NULL;
    bool isbDown=false;
	
	ARM_ReferenceValue* bUpProfile = NULL;
    bool isbUp=false;
	
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "CRA Local Calculator" );

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
		    fundSpreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId)->Clone()));
			(*fundSpreadProfile) /= 100;
		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: fund spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundSpreadProfile = new ARM_ReferenceValue(fundSpread);
			(*fundSpreadProfile) /= 100;
            isFundSpread = true;
        }

		//Cpn Spread Curve
        if (cpnSpreadId != ARM_NULL_OBJECT)
        {
		    cpnSpreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(cpnSpreadId))->Clone());
			(*cpnSpreadProfile) /= 100;
		    if (!cpnSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: Cpn Spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnSpreadProfile = new ARM_ReferenceValue(cpnSpread);
			(*cpnSpreadProfile) /= 100;
            isCpnSpread = true;
        }

		//BoostedFix Curve
        if (boostedFixId != ARM_NULL_OBJECT)
        {
		    boostedFixProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(boostedFixId))->Clone());
			(*boostedFixProfile) /= 100;
		    if (!boostedFixProfile)
		    {
			    result.setMsg ("ARM_ERR: boosted fix is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            boostedFixProfile = new ARM_ReferenceValue(boostedFix);
			(*boostedFixProfile) /= 100;
            isBoostedFix = true;
        }

		//Barrier Down Curve
        if (bDownId != ARM_NULL_OBJECT)
        {
		    bDownProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(bDownId))->Clone());
			(*bDownProfile) /= 100;
		    if (!bDownProfile)
		    {
			    result.setMsg ("ARM_ERR: barrier down is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bDownProfile = new ARM_ReferenceValue(bDown);
			(*bDownProfile) /= 100;
            isbDown = true;
        }

		//Barrier Up Curve
        if (bUpId != ARM_NULL_OBJECT)
        {
		    bUpProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(bUpId))->Clone());
			(*bUpProfile) /= 100;
		    if (!bUpProfile)
		    {
			    result.setMsg ("ARM_ERR: barrier up is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bUpProfile = new ARM_ReferenceValue(bUp);
			(*bUpProfile) /= 100;
            isbUp = true;
        }

		// MRS and Beta
        double meanRevMin = mrsBeta[0];
        double meanRevMax = mrsBeta[1];
        double betaMin = mrsBeta[2];
        double betaMax = mrsBeta[3];
		int reCalibMrs = -1;
		int reCalibBeta = -1;
		if (mrsBeta.size() == 5)
			reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);

		if (mrsBeta.size() == 6)
		{
			reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
			reCalibBeta = (abs(mrsBeta[5] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
		}

		//test de coherence
		if (reCalibMrs == 1 && meanRevMin == meanRevMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Mrs min = Mrs max");
			return ARM_KO;
		}
		if (reCalibBeta == 1 && betaMin == betaMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Beta min = Beta max");
			return ARM_KO;
		}

		ARM_Vector* calibSecParams = new ARM_Vector(calibSecPFParams);

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true; //Bermuda must be computed
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

        // Create the CRA calculator
		if (localModel == true)
		{
			cfcraCalculator = new ARM_CRALocalCalculator(	ccy,
															ARM_Date(myStartDate),
															ARM_Date(myEndDate),
															payReceive,
															*notionalProfile,
															callFreq,
															callNotice,
															callCal,
															*callFeesProfile,
															fundFreq,
															*fundSpreadProfile,
															fundDayCount,
															*cpnSpreadProfile,
															cpnPayFreq,
															cpnResetCal,
															cpnPayCal,
															boostedIndexType,
															*boostedFixProfile,
															boostedVarTerm,
															boostedResetGap,
															boostedResetTiming,
															boostedDayCount,
															boostedAdjRule,
															boostedRule,
															*bDownProfile,
															*bUpProfile,
															cpnResetFreq,
															cpnResetTiming,
															cpnResetGap,
															refIndexType,
															refTerm,
															refDayCount,
															refCoeff,
															localModelType,
															meanRevMin,
															meanRevMax,
															betaMin,
															betaMax,
															calibSecParams,
															nbSteps,
															flagToGenerateOSWATM,
															localResetFreq,
															mdmKeys,
															*mktDataManager,
															productsToPrice,
															reCalibMrs,
															reCalibBeta,
															isStdCalib);
		}
		else
		{
			craCalculator = new ARM_CRACalculator(	ccy,
													ARM_Date(myStartDate),
													ARM_Date(myEndDate),
													payReceive,
													*notionalProfile,
													callFreq,
													callNotice,
													callCal,
													*callFeesProfile,
													fundFreq,
													*fundSpreadProfile,
													fundDayCount,
													*cpnSpreadProfile,
													cpnPayFreq,
													cpnResetCal,
													cpnPayCal,
													boostedIndexType,
													*boostedFixProfile,
													boostedVarTerm,
													boostedResetGap,
													boostedResetTiming,
													boostedDayCount,
													boostedAdjRule,
													boostedRule,
													*bDownProfile,
													*bUpProfile,
													cpnResetFreq,
													cpnResetTiming,
													cpnResetGap,
													refIndexType,
													refTerm,
													refDayCount,
													refCoeff,
													meanRevMin,
													meanRevMax,
													betaMin,
													betaMax,
													calibSecParams,
													nbSteps,
													flagToGenerateOSWATM,
													mdmKeys,
													*mktDataManager,
													productsToPrice,
													reCalibMrs,
													reCalibBeta,
													isStdCalib);
		}
		// Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		if (isCallFees)
			delete callFeesProfile;
		callFeesProfile = NULL;
	
		delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		delete cpnSpreadProfile;
		cpnSpreadProfile = NULL;

		delete boostedFixProfile;
		boostedFixProfile = NULL;

		delete bDownProfile;
		bDownProfile = NULL;
	
		delete bUpProfile;
		bUpProfile = NULL;

		// assign object
		if (localModel == true)
		{
			if ( !assignObject( cfcraCalculator, result, objId ) )
			{
				return ARM_KO; 
			}
			else
			{
				return ARM_OK; 
			}
		}
		else
		{
			if ( !assignObject( craCalculator, result, objId ) )
			{
				return ARM_KO; 
			}
			else
			{
				return ARM_OK; 
			}
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
	
		delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		delete cpnSpreadProfile;
		cpnSpreadProfile = NULL;

		delete boostedFixProfile;
		boostedFixProfile = NULL;

		delete bDownProfile;
		bDownProfile = NULL;
	
		delete bUpProfile;
		bUpProfile = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}


extern long ARMLOCAL_LocalCRA_Set(long craId,
	  long dataId,
	  const string& setType,
	  const vector<string>& keys,
	  bool isUpdated,
	  ARM_result& result, 
	  long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRALocalCalculator* oldCRA = NULL;
    ARM_CRALocalCalculator* newCRA = NULL;
    ARM_Object* object = NULL;
    string updateMsg("Calculator updated by the ");

	try
	{
		oldCRA = dynamic_cast<ARM_CRALocalCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));
		if (!oldCRA)
		{
			result.setMsg ("ARM_ERR: Local CRA Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio*		port = NULL;
		ARM_OptionPortfolio* optPort = NULL;
        if ((port = dynamic_cast<ARM_StdPortfolio*>(object)))
        {
            if (isUpdated)
                /// Just update the initial CRA
                newCRA = oldCRA;
            else
                /// Clone the initial CRA then set
                newCRA = dynamic_cast<ARM_CRALocalCalculator*>(oldCRA->Clone());

            string typeToSet(stringGetUpper(setType));
			if (typeToSet == GC_ACCESS_LOCAL_PORT)
            {	
					newCRA->SetLocalPortfolio(ARM_StdPortfolioPtr(port));
					updateMsg += "Local portfolio";
			}
            if (isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRA in ARM cache
		    if (!assignObject(newCRA, result, objId))
            {
                delete newCRA;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else if ((optPort = dynamic_cast<ARM_OptionPortfolio*>(object)))
        {
            if (isUpdated)
                /// Just update the initial CRA
                newCRA = oldCRA;
            else
                /// Clone the initial CRA then set
                newCRA = dynamic_cast<ARM_CRALocalCalculator*>(oldCRA->Clone());

            string typeToSet(stringGetUpper(setType));
            if (typeToSet == GC_ACCESS_OPT_PORT)
            {
                newCRA->SetOptionPortfolio(optPort);
                updateMsg += "Option portfolio";
            }
            if (isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRA in ARM cache
		    if (!assignObject(newCRA, result, objId))
            {
                delete newCRA;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else
            return ARMLOCAL_GC_Set(craId, dataId, keys, isUpdated, result, objId);
	}
	
	catch (Exception& x)
	{
        if (!isUpdated)
            delete newCRA;
		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_LocalCRA_Get(long craId,
						   const string& getType,
						   ARM_result& result, 
						   long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRALocalCalculator* craCalculator;
    ARM_Object* object = NULL;

	try
	{
		craCalculator = dynamic_cast<ARM_CRALocalCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));

		if (!craCalculator)
		{
			result.setMsg ("ARM_ERR: Local CRA Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if ( (typeToGet == GC_ACCESS_OPT_PORT) ||
			 (typeToGet == GC_ACCESS_LOCAL_PORT))
		{
			if (typeToGet == GC_ACCESS_OPT_PORT)
			{
				object = const_cast< ARM_OptionPortfolio& >(*(craCalculator->GetOptionPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_LOCAL_PORT)
			{
//				bool isLocalModel = craCalculator->GetIsLocalModel();
//				if (isLocalModel)
//				{
					object = const_cast< ARM_StdPortfolio& >(*(craCalculator->GetLocalPortfolio())).Clone();
//				}
			}

			/// Assign the object in ARM cache
			if (!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}
		else
		{
			return ARMLOCAL_GC_Get(craId, getType, result, objId);
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}



long ARMLOCAL_CRACalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		int callFreq,
		int callNotice,
		string callCal,
		int fundFreq,
		int fundDayCount,
		int cpnPayFreq,
		string cpnResetCal,
		string cpnPayCal,
		int boostedIndexType,
		string boostedVarTerm,
		int boostedResetGap,
		int boostedResetTiming,
		int boostedDayCount,
		int boostedAdjRule,
		int boostedRule,
		int cpnResetFreq,
		int cpnResetTiming,
		int cpnResetGap,
		int refIndexType,
		string  refTerm,
		int refDayCount,
		double refCoeff,
		double notional,
		long notionalId,
		long callFeesId,
		double fundSpread,
		long fundSpreadId,
		double cpnSpread,
		long cpnSpreadId,
		double boostedFix,
		long boostedFixId,
		double bDown,
		long bDownId,
		double bUp,
		long bUpId,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool isStdCalib,
		ARM_result&	result, 
        long objId) 
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CRACalculator* craCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* callFeesProfile = NULL;
    bool isCallFees=false;
	
	ARM_ReferenceValue* fundSpreadProfile = NULL;
    bool isFundSpread=false;
	
	ARM_ReferenceValue* cpnSpreadProfile = NULL;
    bool isCpnSpread=false;

	ARM_ReferenceValue* boostedFixProfile = NULL;
    bool isBoostedFix=false;

	ARM_ReferenceValue* bDownProfile = NULL;
    bool isbDown=false;
	
	ARM_ReferenceValue* bUpProfile = NULL;
    bool isbUp=false;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable Range Accrual Calculator" );

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
		    fundSpreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId))->Clone());
			(*fundSpreadProfile) /= 100;

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: fund spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundSpreadProfile = new ARM_ReferenceValue(fundSpread);
			(*fundSpreadProfile) /= 100;
            isFundSpread = true;
        }

		//Cpn Spread Curve
        if (cpnSpreadId != ARM_NULL_OBJECT)
        {
		    cpnSpreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(cpnSpreadId))->Clone());
			(*cpnSpreadProfile) /= 100;

		    if (!cpnSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: Cpn Spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cpnSpreadProfile = new ARM_ReferenceValue(cpnSpread);
			(*cpnSpreadProfile) /= 100;
            isCpnSpread = true;
        }

		//BoostedFix Curve
        if (boostedFixId != ARM_NULL_OBJECT)
        {
		    boostedFixProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(boostedFixId))->Clone());
			(*boostedFixProfile) /= 100;

		    if (!boostedFixProfile)
		    {
			    result.setMsg ("ARM_ERR: boosted fix is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            boostedFixProfile = new ARM_ReferenceValue(boostedFix);
			(*boostedFixProfile) /=100;
            isBoostedFix = true;
        }

		//Barrier Down Curve
        if (bDownId != ARM_NULL_OBJECT)
        {
		    bDownProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(bDownId))->Clone());
			(*bDownProfile) /= 100;
		    if (!bDownProfile)
		    {
			    result.setMsg ("ARM_ERR: barrier down is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bDownProfile = new ARM_ReferenceValue(bDown);
			(*bDownProfile) /= 100;
            isbDown = true;
        }

		//Barrier Up Curve
        if (bUpId != ARM_NULL_OBJECT)
        {
		    bUpProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(bUpId))->Clone());
			(*bUpProfile) /= 100;
		    if (!bUpProfile)
		    {
			    result.setMsg ("ARM_ERR: barrier up is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bUpProfile = new ARM_ReferenceValue(bUp);
			(*bUpProfile) /= 100;
            isbUp = true;
        }

		// MRS and Beta
        double meanRevMin = mrsBeta[0];
        double meanRevMax = mrsBeta[1];
        double betaMin = mrsBeta[2];
        double betaMax = mrsBeta[3];
		// Recalibrate MRS and Beta
		int reCalibMrs = -1;
		int reCalibBeta = -1;
		if (mrsBeta.size() == 5)
			reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);

		if (mrsBeta.size() == 6)
		{
			reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
			reCalibBeta = (abs(mrsBeta[5] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
		}

		//test de coherence
		if (reCalibMrs == 1 && meanRevMin == meanRevMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Mrs min = Mrs max");
			return ARM_KO;
		}
		if (reCalibBeta == 1 && betaMin == betaMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Beta min = Beta max");
			return ARM_KO;
		}

		ARM_Vector* calibSecParams = new ARM_Vector(calibSecPFParams);

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true; //Bermuda must be computed
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

        // Create the CRA calculator
        craCalculator = new ARM_CRACalculator(ccy,
											  ARM_Date(myStartDate),
											  ARM_Date(myEndDate),
											  payReceive,
											  *notionalProfile,
											  callFreq,
											  callNotice,
											  callCal,
											  *callFeesProfile,
											  fundFreq,
											  *fundSpreadProfile,
											  fundDayCount,
											  *cpnSpreadProfile,
											  cpnPayFreq,
											  cpnResetCal,
											  cpnPayCal,
											  boostedIndexType,
											  *boostedFixProfile,
											  boostedVarTerm,
											  boostedResetGap,
											  boostedResetTiming,
											  boostedDayCount,
											  boostedAdjRule,
											  boostedRule,
											  *bDownProfile,
											  *bUpProfile,
											  cpnResetFreq,
											  cpnResetTiming,
											  cpnResetGap,
											  refIndexType,
											  refTerm,
											  refDayCount,
											  refCoeff,
											  meanRevMin,
											  meanRevMax,
											  betaMin,
											  betaMax,
											  calibSecParams,
											  nbSteps,
											  flagToGenerateOSWATM,
											  mdmKeys,
											  *mktDataManager,
											  productsToPrice,
											  reCalibMrs,
											  reCalibBeta,
											  isStdCalib);
		 		
        // Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		if (isCallFees)
			delete callFeesProfile;
		callFeesProfile = NULL;
	
		delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		delete cpnSpreadProfile;
		cpnSpreadProfile = NULL;

		delete boostedFixProfile;
		boostedFixProfile = NULL;

		delete bDownProfile;
		bDownProfile = NULL;
	
		delete bUpProfile;
		bUpProfile = NULL;
	
		// assign object
		if ( !assignObject( craCalculator, result, objId ) )
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
	
		delete fundSpreadProfile;
		fundSpreadProfile = NULL;

		delete cpnSpreadProfile;
		cpnSpreadProfile = NULL;

		delete boostedFixProfile;
		boostedFixProfile = NULL;

		delete bDownProfile;
		bDownProfile = NULL;
	
		delete bUpProfile;
		bUpProfile = NULL;
	
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_CRACalculator_Create(
		int OptionPfId,
		vector<double> mrsBeta,
		vector<double> calibSecPFParams,
		int nbSteps,
		int flagToGenerateOSWATM,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		bool localModel,
		int localModelType,
		int localResetFreq,
		bool isStdCalib,
		ARM_result&	result, 
        long objId)
{
	CCString msg ("");

	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_OptionPortfolio* pOptionPf = NULL;
	
	try
	{
		pOptionPf = (ARM_OptionPortfolio*) LOCAL_PERSISTENT_OBJECTS->GetObject(OptionPfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(pOptionPf, ARM_OPTIONPORTFOLIO) == 0) 
		{
			result.setMsg ("ARM_ERR: Option portfolio is not of a good type");
			return ARM_KO;
		}

        double meanRevMin, meanRevMax, betaMin, betaMax;
		int reCalibMrs = -1;
		int reCalibBeta = -1;
		ARM_Vector* craPricing = pOptionPf->GetCraPricing();

		if (mrsBeta.size() > 0)
		{
			meanRevMin = mrsBeta[0];
			meanRevMax = mrsBeta[1];
			betaMin = mrsBeta[2];
			betaMax = mrsBeta[3];
			if (mrsBeta.size() == 5)
				reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);

			if (mrsBeta.size() == 6)
			{
				reCalibMrs = (abs(mrsBeta[4] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
				reCalibBeta = (abs(mrsBeta[5] - 1.0) < K_DOUBLE_TOL ? 1 : 0);
			}
		}
		else
		{
			meanRevMin	= craPricing->Elt(2);
			meanRevMax	= craPricing->Elt(3);
			betaMin		= craPricing->Elt(4);
			betaMax		= craPricing->Elt(5);
		}

		//test de coherence
		if (reCalibMrs == 1 && meanRevMin == meanRevMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Mrs min = Mrs max");
			return ARM_KO;
		}
		if (reCalibBeta == 1 && betaMin == betaMax)
		{
			result.setMsg ("ARM_ERR: can not calibrate if Beta min = Beta max");
			return ARM_KO;
		}

		ARM_Vector* calibSecParams = NULL;
		if (calibSecPFParams.size() > 0)
		{
			calibSecParams = new ARM_Vector(calibSecPFParams);
		}
		else
		{
			calibSecParams = new ARM_Vector(8, 0.0);

			calibSecParams->Elt(0) = 12.0/pOptionPf->GetCorridorLeg()->GetRefIndex()->GetTerm();
			calibSecParams->Elt(1) = 12.0/pOptionPf->GetCorridorLeg()->GetCurrencyUnit()->GetFixedPayFreq();
			calibSecParams->Elt(2) = 0;
			calibSecParams->Elt(3) = 1;
			calibSecParams->Elt(4) = 1;
			calibSecParams->Elt(5) = craPricing->Elt(0);
			calibSecParams->Elt(6) = craPricing->Elt(1);
			calibSecParams->Elt(7) = 12;
		}

		if (flagToGenerateOSWATM == -1)
			flagToGenerateOSWATM = craPricing->Elt(7);

		if (abs(nbSteps + 1.0) < K_DOUBLE_TOL)
			nbSteps = craPricing->Elt(6);

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true; //Bermuda must be computed
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

        // Create the CRA calculator
		ARM_CRALocalCalculator* cfcraCalculator = NULL;
		ARM_CRACalculator* craCalculator = NULL;
        
		if (localModel == true)
		{
			cfcraCalculator = new ARM_CRALocalCalculator(	pOptionPf,
															meanRevMin,
															meanRevMax,
															betaMin,
															betaMax,
															calibSecParams,
															nbSteps,
															flagToGenerateOSWATM,
															mdmKeys,
															*mktDataManager,
															productsToPrice,
															reCalibMrs,
															reCalibBeta,
															localModelType,
															localResetFreq,
															isStdCalib);		
		}
		else
		{
			craCalculator = new ARM_CRACalculator(pOptionPf,
						                          meanRevMin,
							                      meanRevMax,
								                  betaMin,
									              betaMax,
										          calibSecParams,
											      nbSteps,
											      flagToGenerateOSWATM,
												  mdmKeys,
												  *mktDataManager,
												  productsToPrice,
												  reCalibMrs,
												  reCalibBeta,
												  isStdCalib);
		}

		if (calibSecParams)
		{
			delete calibSecParams;
			calibSecParams = NULL;
		}

		// assign object
		if (localModel == true)
		{
			if ( !assignObject( cfcraCalculator, result, objId ) )
			{
				return ARM_KO; 
			}
			else
			{
				return ARM_OK; 
			}
		}
		else
		{
			if ( !assignObject( craCalculator, result, objId ) )
			{
				return ARM_KO; 
			}
			else
			{
				return ARM_OK; 
			}
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_CRA_Get(long craId,
					 const string& getType,
					 ARM_result& result, 
					 long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRACalculator* craCalculator;
    ARM_Object* object = NULL;

	try
	{
		craCalculator = dynamic_cast<ARM_CRACalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));

		if (!craCalculator)
		{
			result.setMsg ("ARM_ERR: CRA Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if ((typeToGet == GC_ACCESS_OSW_PORT) || 
			(typeToGet == GC_ACCESS_SHORT_TERM_PORT) || 
			(typeToGet == GC_ACCESS_CALIB) ||
			(typeToGet == GC_ACCESS_OPT_PORT))
		{
			if (typeToGet == GC_ACCESS_OPT_PORT)
			{
				object = const_cast< ARM_OptionPortfolio& >(*(craCalculator->GetOptionPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_OSW_PORT)
			{
				object = const_cast< ARM_StdPortfolio& >(*(craCalculator->GetOSWPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_SHORT_TERM_PORT)
			{
//				object = const_cast< ARM_StdPortfolio& >(*(craCalculator->GetSTMPortfolio())).Clone();
			}
			else if (typeToGet == GC_ACCESS_CALIB)
			{
				object = const_cast< ARM_CalibMethod& >(*(craCalculator->GetOSWCalibMethod())).Clone();
			}

			/// Assign the object in ARM cache
			if (!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}
		else
		{
			return ARMLOCAL_GC_Get(craId, getType, result, objId);
		}
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a CRA Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_CRA_Set(long craId,
		 long dataId,
		 const string& setType,
		 const vector<string>& keys,
		 bool isUpdated,
		 ARM_result& result, 
		 long objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRACalculator* oldCRA = NULL;
    ARM_CRACalculator* newCRA = NULL;
    ARM_Object* object = NULL;
    string updateMsg("Calculator updated by the ");

	try
	{
		oldCRA = dynamic_cast<ARM_CRACalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));
		if (!oldCRA)
		{
			result.setMsg ("ARM_ERR: CRA Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if ((port = dynamic_cast<ARM_StdPortfolio*>(object)))
        {
            if (isUpdated)
                /// Just update the initial CRA
                newCRA = oldCRA;
            else
                /// Clone the initial CRA then set
                newCRA = dynamic_cast<ARM_CRACalculator*>(oldCRA->Clone());

            string typeToSet(stringGetUpper(setType));
            if (typeToSet == GC_ACCESS_OSW_PORT)
            {
//                newCRA->SetOSWPortfolio(*port);
                updateMsg += "OSW portfolio";
            }
            else if (typeToSet == GC_ACCESS_SHORT_TERM_PORT)
            {
//                newCRA->SetSTMPortfolio(*port);
                updateMsg += "STM vanillas portfolio";
            }
		
            if (isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRA in ARM cache
		    if (!assignObject(newCRA, result, objId))
            {
                delete newCRA;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else
            return ARMLOCAL_GC_Set(craId, dataId, keys, isUpdated, result, objId);
	}
	
	catch (Exception& x)
	{
        if (!isUpdated)
            delete newCRA;
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_CRA_SetRecalibFlags(long craId,
								  int mrsFlag,
								  int betaFlag,
								  ARM_result& result,
								  long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_CRACalculator* cra = dynamic_cast<ARM_CRACalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId));
		if (!cra)
		{
			result.setMsg ("ARM_ERR: CRA Calculator is not of a good type");
			return ARM_KO;
		}

		ARM_CRALocalCalculator* clonedCra = dynamic_cast<ARM_CRALocalCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId)->Clone());
		if (clonedCra)
		{
			clonedCra->SetReCalibFlags(mrsFlag, betaFlag);
			
			/// assign object
			if( !assignObject( clonedCra, result, objId ) ){
				return ARM_KO; }
			else{
				return ARM_OK; }
		}
		else
		{
			ARM_CRACalculator* clonedCra = dynamic_cast<ARM_CRACalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(craId)->Clone());
			clonedCra->SetReCalibFlags(mrsFlag, betaFlag);
			
			/// assign object
			if( !assignObject( clonedCra, result, objId ) ){
				return ARM_KO; }
			else{
				return ARM_OK; }
		}
	}
	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////// CRA Spread Calculator
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

long ARMLOCAL_CRASpreadCalculator_Create(
		const ARM_Currency& ccy,
		const ARM_Currency& fundCcy,
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
		const int& cpnResetGap,
		const int& refIndex1,
		const int& refIndex2,
		const int& payIndex,
		const int& payIndexResetTiming,
		const double& notional,
		const long& notionalId,
		const double& fundNotionaldble,
		const long& callFeesId,
		const double& fundSpread,
		const long& fundSpreadId,
		const double& boostedFix,
		const long& boostedFixId,
		const double& payIndexMult,
		const long& payIndexMultId,
		const double& bDown,
		const long& bDownId,
		const double& bUp,
		const long& bUpId,
		const double& coeff1,
		const double& coeff2,
		const int& refCond1Index,
		const double& bCond1Down,
		const long& bCond1DownId,
		const double& bCond1Up,
		const long& bCond1UpId,
		const long& boostedFix2Id,
		const long& bDown2Id,
		const long& bUp2Id,
		const long& boostedFix3Id,
		const long& bDown3Id,
		const long& bUp3Id,
		const vector<string>& tenorVec,
		const vector<string>& CalibMod,
		const vector<string>& mdmKeys,
		const long& mktDataManagerId,
		const vector<string>& pricingFlags,
		const vector<string>& localCalibFlags,
		const vector<double>& optimResetData,
		const vector<double>& exerProbas,
		ARM_result&	result, 
        long objId)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CRASpreadCalculator* craSpreadCalculator = NULL;
	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* callFeesProfile = NULL;
    bool isCallFees=false;
	
	ARM_ReferenceValue* fundSpreadProfile = NULL;
    bool isFundSpread=false;

	ARM_ReferenceValue* boostedFixProfile = NULL;
    bool isBoostedFix=false;
	
	ARM_ReferenceValue* payIndexMultProfile = NULL;
    bool isPayIndexMult=false;

	ARM_ReferenceValue* bDownProfile = NULL;
    bool isbDown=false;
	
	ARM_ReferenceValue* bUpProfile = NULL;
    bool isbUp=false;

	bool isDbleCondition = (refCond1Index != K_FIXED); // => double condition is activated

	ARM_ReferenceValue* bCond1DownProfile = NULL;
    bool isbCond1Down=false;
	
	ARM_ReferenceValue* bCond1UpProfile = NULL;
    bool isbCond1Up=false;

	bool isTripleRange = (boostedFix2Id != ARM_NULL_OBJECT); // triple range activated

	ARM_ReferenceValue* boostedFix2Profile = NULL;
    bool isBoostedFix2=false;
	
	ARM_ReferenceValue* bDown2Profile = NULL;
    bool isbDown2=false;
	
	ARM_ReferenceValue* bUp2Profile = NULL;
    bool isbUp2=false;

	ARM_ReferenceValue* boostedFix3Profile = NULL;
    bool isBoostedFix3=false;
	
	ARM_ReferenceValue* bDown3Profile = NULL;
    bool isbDown3=false;
	
	ARM_ReferenceValue* bUp3Profile = NULL;
    bool isbUp3=false;

	bool isVms = tenorVec.size()>0; // VMS activated

	ARM_ReferenceValue* refTenor1Profile = NULL;
    bool isRefTenor1=false;

	try
	{
		// crm tracing
		if(isDbleCondition)
			ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable Double Condition Range Accrual Calculator" );
		else if(isTripleRange)
			ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable on Triple Range Accrual Spread Calculator" );
		else if(isVms)
			ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable on VMS Range Accrual Calculator" );
		else
			ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable Range Accrual Spread Calculator" );

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

		ARM_ReferenceValue fundnotionalProfile(fundNotionaldble);

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

		//PayIndexMult Curve
        if (payIndexMultId != ARM_NULL_OBJECT)
        {
		    payIndexMultProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexMultId));

		    if (!payIndexMultProfile)
		    {
			    result.setMsg ("ARM_ERR: payIndexMult is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            payIndexMultProfile = new ARM_ReferenceValue(payIndexMult);
            isPayIndexMult = true;
        }

		//Barrier Down Curve
        if (bDownId != ARM_NULL_OBJECT)
        {
		    bDownProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bDownId));

		    if (!bDownProfile)
		    {
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier down for 2nd condition is not of a good type");
				else
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
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier up for 2nd condition is not of a good type");
				else
					result.setMsg ("ARM_ERR: barrier up is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            bUpProfile = new ARM_ReferenceValue(bUp);
            isbUp = true;
        }

		//Barrier Down Curve
		if (bCond1DownId != ARM_NULL_OBJECT)
		{
			bCond1DownProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bCond1DownId));

			if (!bCond1DownProfile)
			{
				result.setMsg ("ARM_ERR: barrier down for 1st condition is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			bCond1DownProfile = new ARM_ReferenceValue(bCond1Down);
			isbCond1Down = true;
		}



		// Barrier Up Curve
		if (bCond1UpId != ARM_NULL_OBJECT)
		{
			bCond1UpProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bCond1UpId));

			if (!bCond1UpProfile)
			{
				result.setMsg ("ARM_ERR: barrier up for 1st condition is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			bCond1UpProfile = new ARM_ReferenceValue(bCond1Up);
			isbCond1Up = true;
		}

// for triple range
		//BoostedFix Curve
        if (boostedFix2Id != ARM_NULL_OBJECT)
        {
		    boostedFix2Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(boostedFix2Id));

		    if (!boostedFix2Profile)
		    {
			    result.setMsg ("ARM_ERR: boosted Fix2 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			boostedFix2Profile = new ARM_ReferenceValue(0.);
            isBoostedFix2 = true;
        }

		if (boostedFix3Id != ARM_NULL_OBJECT)
        {
		    boostedFix3Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(boostedFix3Id));

		    if (!boostedFix3Profile)
		    {
			    result.setMsg ("ARM_ERR: boosted Fix3 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			boostedFix3Profile = new ARM_ReferenceValue(0.);
            isBoostedFix3 = true;
        }

		//Barrier Down Curve
        if (bDown2Id != ARM_NULL_OBJECT)
        {
		    bDown2Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bDown2Id));

		    if (!bDown2Profile)
		    {
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier Down2 for 2nd condition is not of a good type");
				else
					result.setMsg ("ARM_ERR: barrier Down2 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			bDown2Profile = new ARM_ReferenceValue(0.);
            isbDown2 = true;
        }

		if (bDown3Id != ARM_NULL_OBJECT)
        {
		    bDown3Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bDown3Id));

		    if (!bDown3Profile)
		    {
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier Down3 for 2nd condition is not of a good type");
				else
					result.setMsg ("ARM_ERR: barrier Down3 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			bDown3Profile = new ARM_ReferenceValue(0.);
            isbDown3 = true;
        }
 
		// Barrier Up Curve
        if (bUp2Id != ARM_NULL_OBJECT)
        {
		    bUp2Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bUp2Id));

		    if (!bUp2Profile)
		    {
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier Up2 for 2nd condition is not of a good type");
				else
					result.setMsg ("ARM_ERR: barrier Up2 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			bUp2Profile = new ARM_ReferenceValue(0.);
            isbUp2 = true;
        }
		
		if (bUp3Id != ARM_NULL_OBJECT)
        {
		    bUp3Profile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(bUp3Id));

		    if (!bUp3Profile)
		    {
				if(isDbleCondition)
					result.setMsg ("ARM_ERR: barrier Up3 for 2nd condition is not of a good type");
				else
					result.setMsg ("ARM_ERR: barrier Up3 is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			bUp3Profile = new ARM_ReferenceValue(0.);
            isbUp3 = true;
        }

// for triple range END

// for VMS

		if (isVms)
        {
		    refTenor1Profile = (ARM_ReferenceValue*)bUpProfile->Clone();

			if (refTenor1Profile->size()==tenorVec.size())
			{
				for (int i=0;i<refTenor1Profile->size();i++)
					refTenor1Profile->GetDiscreteValues()->Elt(i) = ARM_ConvIrType(tenorVec[i].c_str());
			}
			else
			{
				result.setMsg ("ARM_ERR: pls check notional and tenor size");
				return ARM_KO;
			}
        }
        else
        {
			refTenor1Profile = new ARM_ReferenceValue(0.);
            isRefTenor1 = true;
        }


// for VMS END

		ARM_CRASpreadCalculator::CalibrationType calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
		ARM_CRASpreadCalculator::CalibStrikeType calibStrikeType = ARM_CRASpreadCalculator::ATM;
		ARM_ModelType modelType = ARM_PricingModelType::HWM1F;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;

		vector<ARM_CRASpreadCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		
		if (CalibMod.size()>4) //switch to new calib hw2f
		{
			if (CalibMod.size()==6)
			{
				//MODEL
				if (CalibMod[0] == "HW1F")
					modelType = ARM_PricingModelType::HWM1F;
				else if (CalibMod[0] == "HW2F")
					modelType = ARM_PricingModelType::HWM2F;
				else
				{
					result.setMsg ("ARM_ERR: invalid model type");
					return ARM_KO;
				}
				
				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (CalibMod[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_LONG;
					else if (CalibMod[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_SHORT;
					else if (CalibMod[1] == "DIAG_INDEX_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_INDEX;
					else if (CalibMod[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::SHORT_LONG_SPREAD;
					else if (CalibMod[1] == "DIAG")
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "BASKET_SIMPLE")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (CalibMod[1] == "DIAG") 
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "BASKET_SIMPLE") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}

				if (CalibMod[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (CalibMod[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);
				else if (CalibMod[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (CalibMod[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (CalibMod[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ZERO);
				else if (CalibMod[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::CAP);
				else if (CalibMod[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (CalibMod[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if( calibType == ARM_CRASpreadCalculator::DIAG_CALIBRATION ||
						 calibType == ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED)
				{
					if (CalibMod[3] == "ZERO") 
						calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ZERO);
					else if (CalibMod[3] == "CAP") 
						calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::CAP);
					else if (CalibMod[3] == "FLOOR") 
						calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FLOOR);
					else
					{
						result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
						return ARM_KO;
					}
				}
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}

			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[5]);
		}
		else
		{
			if (CalibMod[0] == "HW1F")
				modelType = ARM_PricingModelType::HWM1F;
			else if (CalibMod[0] == "HW2F")
				modelType = ARM_PricingModelType::HWM2F;
			else
			{
				result.setMsg ("ARM_ERR: invalid model type");
				return ARM_KO;
			}

			if (CalibMod.size()>1)
			{
				if (CalibMod[1] == "DIAG") 
					calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;

				else if (CalibMod[1] == "BASKET") 
					calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;

				else if (CalibMod[1] == "BASKET_SIMPLE") 
					calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;

				else if (CalibMod[1] == "DIAG_LONG_SPREAD")
					calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_LONG;

				else if (CalibMod[1] == "DIAG_SHORT_SPREAD")
					calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_SHORT;

				else if (CalibMod[1] == "SHORT_LONG_SPREAD")
					calibType = ARM_CRASpreadCalculator::SHORT_LONG_SPREAD;

				else
				{
					result.setMsg ("ARM_ERR: invalid calib type");
					return ARM_KO;
				}
			}

			if (CalibMod.size()>2)
			{
				if (CalibMod[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);

				else if (CalibMod[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);

				else if (CalibMod[2] == "FRONTIER") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FRONTIER);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (CalibMod.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[3]);
		}


		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

        if (pricingFlags.size()>ARM::ARM_CRALocalCalculator::localCfNbProductsToPrice)
		{
			result.setMsg ("ARM_ERR: too much products to price");
			return ARM_KO;
		}
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(ARM::ARM_CRALocalCalculator::localCfNbProductsToPrice, false);

		for (int i = 0; i < pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ localCalibSelectors(localCalibFlags.size(), false);
		for (i = 0; i < localCalibFlags.size(); ++i)
		{
			localCalibSelectors[i] = ((localCalibFlags[i] == "Y") || (localCalibFlags[i] == "YES"));
		}

		ARM_GP_Vector optimResetDataGP(optimResetData.size());
		for(i = 0;i<optimResetData.size();++i) optimResetDataGP[i] = optimResetData[i];

	    // Create the CRA calculator
		ARM_CRALocalCalculator* cfcraCalculator = NULL;
		ARM_CRACalculator* craCalculator = NULL;

		        
        // Create the CRA calculator
        craSpreadCalculator = new ARM_CRASpreadCalculator(ccy,
											  ARM_Date(myStartDate),
											  ARM_Date(myEndDate),
											  payReceive,
											  *notionalProfile,
											  callFreq,
											  callNotice,
											  callCal,
											  *callFeesProfile,
											  fundFreq,
											  *fundSpreadProfile,
											  fundDayCount,
											  cpnDayCount,
											  cpnPayFreq,
											  cpnResetCal,
											  cpnPayCal,
											  payIndex,
											  payIndexResetTiming,
											  *boostedFixProfile,
											  *payIndexMultProfile,
											  *bDownProfile,
											  *bUpProfile,
											  cpnResetFreq,
											  cpnResetTiming,
											  cpnResetGap,
											  refIndex1,
											  coeff1,
											  refIndex2,
											  coeff2,
											  refCond1Index,
											  *bCond1DownProfile,
											  *bCond1UpProfile,
											  isTripleRange,
											  isVms,
											  *boostedFix2Profile,
											  *bDown2Profile,
											  *bUp2Profile,
											  *boostedFix3Profile,
											  *bDown3Profile,
											  *bUp3Profile,
											  *refTenor1Profile,
											  modelType,
											  calibType,
											  calibStrikeTypeVec,
											  vnsMethod,
											  mdmKeys,
											  *mktDataManager,
											  productsToPrice,
											  localCalibSelectors,
											  optimResetDataGP,
											  exerProbas[0]==1.0,
											  static_cast<size_t>(floor(exerProbas[1])),
											  const_cast<ARM_Currency*>(&fundCcy),
											  fundnotionalProfile);
		 		
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

		if (isPayIndexMult)
			delete payIndexMultProfile;
		payIndexMultProfile = NULL;

		if (isbDown)
			delete bDownProfile;
		bDownProfile = NULL;
	
		if (isbUp)
			delete bUpProfile;
		bUpProfile = NULL;
	
		if (isbCond1Down)
			delete bCond1DownProfile;
	
		if (isbCond1Up)
			delete bCond1UpProfile;

		if (isBoostedFix2)
			delete boostedFix2Profile;
		boostedFix2Profile = NULL;

		if (isbDown2)
			delete bDown2Profile;
		bDown2Profile = NULL;
	
		if (isbUp2)
			delete bUp2Profile;
		bUp2Profile = NULL;

		if (isBoostedFix3)
			delete boostedFix3Profile;
		boostedFix3Profile = NULL;

		if (isbDown3)
			delete bDown3Profile;
		bDown3Profile = NULL;
	
		if (isbUp3)
			delete bUp3Profile;
		bUp3Profile = NULL;

		if (isRefTenor1)
			delete refTenor1Profile;
		refTenor1Profile = NULL;
	
		// assign object
		if ( !assignObject( craSpreadCalculator, result, objId ) )
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

		if (isPayIndexMult)
			delete payIndexMultProfile;
		payIndexMultProfile = NULL;

		if (isbDown)
			delete bDownProfile;
		bDownProfile = NULL;
	
		if (isbUp)
			delete bUpProfile;
		bUpProfile = NULL;
	
		if (isbCond1Down)
			delete bCond1DownProfile;
	
		if (isbCond1Up)
			delete bCond1UpProfile;

		if (isBoostedFix2)
			delete boostedFix2Profile;
		boostedFix2Profile = NULL;

		if (isbDown2)
			delete bDown2Profile;
		bDown2Profile = NULL;
	
		if (isbUp2)
			delete bUp2Profile;
		bUp2Profile = NULL;

		if (isBoostedFix3)
			delete boostedFix3Profile;
		boostedFix3Profile = NULL;

		if (isbDown3)
			delete bDown3Profile;
		bDown3Profile = NULL;
	
		if (isbUp3)
			delete bUp3Profile;
		bUp3Profile = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_CRASpreadCalculator_Create(int OptionPfId,
										int refResetFreq,
										vector<string> mdmKeys,
										vector<string> CalibMod,
										double payIndexMultValue, 
										long payIndexMultId,
										long mktDataManagerId,
										vector<string> pricingFlags,
										ARM_result&	result, 
										long objId)
{
	CCString msg ("");

	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_OptionPortfolio* inOptionPf = NULL;
	ARM_OptionPortfolio* pOptionPf = NULL;
	ARM_ReferenceValue*	payIndexMult = NULL;
	
	try
	{
		inOptionPf = (ARM_OptionPortfolio*) LOCAL_PERSISTENT_OBJECTS->GetObject(OptionPfId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inOptionPf, ARM_OPTIONPORTFOLIO) == 0) 
		{
			result.setMsg ("ARM_ERR: Option portfolio is not of a good type");
			return ARM_KO;
		}

		pOptionPf = (ARM_OptionPortfolio*) inOptionPf->Clone();

		// set an inputed frequency
		// useful for tests
		if (refResetFreq != K_DEF_FREQ)
			pOptionPf->GetCorridorLeg()->GetRefIndex()->SetResetFrequency(refResetFreq);

        if (payIndexMultId != ARM_NULL_OBJECT)
		{
			payIndexMult = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexMultId));

		    if(payIndexMult == NULL)
		    {
			    result.setMsg ("ARM_ERR: payIndexMult is not of a good type");
			    return ARM_KO;
		    }
		}
		else
		{
			payIndexMult = new ARM_ReferenceValue(payIndexMultValue);
		}

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true; //Bermuda must be computed
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

        // Force Down Barrier to -100

        pOptionPf->GetCorridorLeg()->GetDownBarriers()->SetCalcMethod(K_STEPUP_RIGHT);

        ARM_Vector* BarDownValues = pOptionPf->GetCorridorLeg()->GetDownBarriers()->GetDiscreteValues();
        
        int bDownSz = BarDownValues->GetSize();

        for (i = 0; i < bDownSz; i++)
        {
            BarDownValues->Elt(i) = -10000.0;
        }

        // Force UP Barrier

		pOptionPf->GetCorridorLeg()->GetUpBarriers()->SetCalcMethod(K_STEPUP_RIGHT);
		*(pOptionPf->GetCorridorLeg()->GetUpBarriers()) *= 100.0;

		ARM_ReferenceValue* theSpreads = pOptionPf->GetLiborLeg()->GetSpreads();

		ARM_Vector* resetDates = pOptionPf->GetLiborLeg()->GetResetDates();
		ARM_Vector* stDates =  (ARM_Vector *) pOptionPf->GetLiborLeg()->GetFlowStartDates()->Clone();

		if (theSpreads)
		{
		   ARM_Vector* resSpread = new ARM_Vector(resetDates->GetSize());

		   for (int i = 0; i < resetDates->GetSize(); i++)
		   {
        	   resSpread->Elt(i) = theSpreads->CptReferenceValue(resetDates->Elt(i));
		   }
		   		    
           theSpreads->SetDiscreteDates(stDates);

		   theSpreads->SetDiscreteValues(resSpread);
		
		   theSpreads->SetCalcMethod(K_LINEAR);
		}
		else
		{
		    ARM_Vector* resSpread = new ARM_Vector(resetDates->GetSize(), pOptionPf->GetLiborLeg()->GetSpread());

		    ARM_ReferenceValue genSpreads(stDates, resSpread);
					
		    genSpreads.SetCalcMethod(K_LINEAR);

            pOptionPf->GetLiborLeg()->SetSpreads(&genSpreads);
            
            theSpreads = pOptionPf->GetLiborLeg()->GetSpreads();
		}
		*theSpreads *= 100.0;

        // Boosted Fix Rate

         pOptionPf->GetCorridorLeg()->GetSpreads()->SetCalcMethod(K_STEPUP_RIGHT);
        *(pOptionPf->GetCorridorLeg()->GetSpreads()) *= 100.0;

		ARM_CRASpreadCalculator::CalibrationType calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
		ARM_CRASpreadCalculator::CalibStrikeType calibStrikeType = ARM_CRASpreadCalculator::ATM;
		ARM_ModelType modelType = ARM_PricingModelType::HWM1F;
		ARM_MarketIRModel::VnsPricingMethod vnsMethod = ARM_MarketIRModel::MONEYNESS;

		vector<ARM_CRASpreadCalculator::CalibStrikeType> calibStrikeTypeVec(0);

		
		if (CalibMod.size()>4) //switch to new calib hw2f
		{
			if (CalibMod.size()==6)
			{
				//MODEL
				if (CalibMod[0] == "HW1F")
					modelType = ARM_PricingModelType::HWM1F;
				else if (CalibMod[0] == "HW2F")
					modelType = ARM_PricingModelType::HWM2F;
				else
				{
					result.setMsg ("ARM_ERR: invalid model type");
					return ARM_KO;
				}
				
				if (modelType==ARM_PricingModelType::HWM2F)
				{
					//CALIB PATTERN
					if (CalibMod[1] == "DIAG_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_LONG;
					else if (CalibMod[1] == "DIAG_SHORT_SPREAD")
						calibType = ARM_CRASpreadCalculator::DIAG_SPREAD_SHORT;
					else if (CalibMod[1] == "SHORT_LONG_SPREAD")
						calibType = ARM_CRASpreadCalculator::SHORT_LONG_SPREAD;
					else if (CalibMod[1] == "DIAG")
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "BASKET_SIMPLE")
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}
				else
				{
					if (CalibMod[1] == "DIAG") 
						calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;
					else if (CalibMod[1] == "BASKET") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;
					else if (CalibMod[1] == "BASKET_SIMPLE") 
						calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;
					else
					{
						result.setMsg ("ARM_ERR: invalid calib type");
						return ARM_KO;
					}
				}

				if (CalibMod[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (CalibMod[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);
				else
				{
					result.setMsg ("ARM_ERR: invalid 1st calib strike type");
					return ARM_KO;
				}

				if (CalibMod[4] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else if (CalibMod[4] == "ZERO") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ZERO);
				else if (CalibMod[4] == "CAP") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::CAP);
				else if (CalibMod[4] == "FLOOR") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::FLOOR);
				else
				{
					result.setMsg ("ARM_ERR: invalid 3rd calib strike type");
					return ARM_KO;
				}

				if (CalibMod[3] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);
				else
				{
					result.setMsg ("ARM_ERR: invalid 2nd calib strike type");
					return ARM_KO;
				}
			}
			else
			{
				result.setMsg ("ARM_ERR: invalid calibflags size, should be 6");
				return ARM_KO;
			}
			vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[5]);
		}
		else
		{
			if (CalibMod[0] == "HW1F")
				modelType = ARM_PricingModelType::HWM1F;
			else if (CalibMod[0] == "HW2F")
				modelType = ARM_PricingModelType::HWM2F;
			else
			{
				result.setMsg ("ARM_ERR: invalid model type");
				return ARM_KO;
			}

			if (CalibMod.size()>1)
			{
				if (CalibMod[1] == "DIAG") 
					calibType = ARM_CRASpreadCalculator::DIAG_CALIBRATION;

				else if (CalibMod[1] == "BASKET") 
					calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION;

				else if (CalibMod[1] == "BASKET_SIMPLE") 
					calibType = ARM_CRASpreadCalculator::BASKET_CALIBRATION_SIMPLIFIED;

				else
				{
					result.setMsg ("ARM_ERR: invalid calib type");
					return ARM_KO;
				}
			}

			if (CalibMod.size()>2)
			{
				if (CalibMod[2] == "ATM") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::ATM);

				else if (CalibMod[2] == "EQUIVALENT") 
					calibStrikeTypeVec.push_back(ARM_CRASpreadCalculator::EQUIVALENT);
				else
				{
					result.setMsg ("ARM_ERR: invalid calib strike type");
					return ARM_KO;
				}
			}

			if (CalibMod.size()>3)
				vnsMethod = (ARM_MarketIRModel::VnsPricingMethod)ARM_ArgConv_VnsPricingMethod.GetNumber(CalibMod[3]);
		}

        // Create the CRA Spread calculator
		ARM_CRASpreadCalculator* craSpreadCalculator = new ARM_CRASpreadCalculator(
														pOptionPf,
														mdmKeys,
													    modelType,
													    calibType,
													    calibStrikeTypeVec,
													    vnsMethod,
														*payIndexMult,
														*mktDataManager,
														productsToPrice);

		delete pOptionPf;
		if ( !assignObject( craSpreadCalculator, result, objId ) )
		{
			if(payIndexMultId == ARM_NULL_OBJECT)
				delete payIndexMult;

			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}

	}
	catch(Exception& x)
	{
		if(payIndexMultId == ARM_NULL_OBJECT)
			delete payIndexMult;

		delete pOptionPf;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_CRASpreadCalculator_Create(double asOfDate,
										int securityId,
										int refResetFreq,
										double payIndexMultValue,
										long payIndexMultId,
										ARM_result&	result, 
										long objId)
{
	CCString msg ("");

	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	ARM_OptionPortfolio* inOptionPf = NULL;
	ARM_OptionPortfolio* pOptionPf = NULL;

	ARM_Swaption* inSwaption = NULL;
	ARM_Swaption* pSwaption = NULL;

	ARM_ReferenceValue*	payIndexMult = NULL;
	
	ARM_CRASpreadCalculator* craSpreadCalculator = NULL;

	try
	{
		char myAsOfDate[20];
		Local_XLDATE2ARMDATE(asOfDate, myAsOfDate);

		inOptionPf = dynamic_cast<ARM_OptionPortfolio*>(LOCAL_PERSISTENT_OBJECTS->GetObject(securityId));
		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inOptionPf, ARM_OPTIONPORTFOLIO) == 1) 
		{
			pOptionPf = (ARM_OptionPortfolio*) inOptionPf->Clone();

			// set an inputed frequency
			// useful for tests
			if (refResetFreq != K_DEF_FREQ)
				pOptionPf->GetCorridorLeg()->GetRefIndex()->SetResetFrequency(refResetFreq);

			if (payIndexMultId != ARM_NULL_OBJECT)
			{
				payIndexMult = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexMultId));

				if(payIndexMult == NULL)
				{
					result.setMsg ("ARM_ERR: payIndexMult is not of a good type");
					return ARM_KO;
				}
			}
			else
			{
				payIndexMult = new ARM_ReferenceValue(payIndexMultValue);
			}
			
			craSpreadCalculator = new ARM_CRASpreadCalculator(ARM_Date(myAsOfDate),
															  pOptionPf,
															  *payIndexMult);

			delete pOptionPf;
		}
		else
		{
			inSwaption = dynamic_cast<ARM_Swaption*>(LOCAL_PERSISTENT_OBJECTS->GetObject(securityId));
			if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(inSwaption, ARM_SWAPTION) == 0) 
			{
				result.setMsg ("ARM_ERR: Security is not of a good type (option pf or swaption");
				return ARM_KO;
			}
				
			pSwaption = (ARM_Swaption*) inSwaption->Clone();

			if (payIndexMultId != ARM_NULL_OBJECT)
			{
				payIndexMult = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexMultId));

				if(payIndexMult == NULL)
				{
					result.setMsg ("ARM_ERR: payIndexMult is not of a good type");
					return ARM_KO;
				}
			}
			else
			{
				payIndexMult = new ARM_ReferenceValue(payIndexMultValue);
			}

			craSpreadCalculator = new ARM_CRASpreadCalculator(ARM_Date(myAsOfDate),
															  pSwaption,
															  *payIndexMult);

			delete pSwaption;
		}

		if ( !assignObject( craSpreadCalculator, result, objId ) )
		{
			if(payIndexMultId == ARM_NULL_OBJECT)
				delete payIndexMult;

			return ARM_KO; 
		}
		else
		{
			return ARM_OK; 
		}

	}
	catch(Exception& x)
	{
		if(payIndexMultId == ARM_NULL_OBJECT)
			delete payIndexMult;

		delete pOptionPf;

		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_CRASpreadCalculator_Create(const ARM_Currency& ccy,
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
										 const int& cpnResetGap,
										 const int& refIndex1,
										 const int& refIndex2,
										 const int& payIndex,
										 const int& payIndexResetTiming,
										 const double& notional,
										 const long& notionalId,
										 const long& callFeesId,
										 const double& fundSpread,
										 const long& fundSpreadId,
										 const double& boostedFix,
										 const long& boostedFixId,
										 const double& payIndexMult,
										 const long& payIndexMultId,
										 const double& bDown,
										 const long& bDownId,
										 const double& bUp,
										 const long& bUpId,
										 const double& coeff1,
										 const double& coeff2,
										 ARM_result&	result, 
										 long objId)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_CRASpreadCalculator* craSpreadCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* callFeesProfile = NULL;
    bool isCallFees=false;
	
	ARM_ReferenceValue* fundSpreadProfile = NULL;
    bool isFundSpread=false;

	ARM_ReferenceValue* boostedFixProfile = NULL;
    bool isBoostedFix=false;
	
	ARM_ReferenceValue* payIndexMultProfile = NULL;
    bool isPayIndexMult=false;

	ARM_ReferenceValue* bDownProfile = NULL;
    bool isbDown=false;
	
	ARM_ReferenceValue* bUpProfile = NULL;
    bool isbUp=false;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable Range Accrual Spread Calculator" );

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

		//PayIndexMult Curve
        if (payIndexMultId != ARM_NULL_OBJECT)
        {
		    payIndexMultProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(payIndexMultId));

		    if (!payIndexMultProfile)
		    {
			    result.setMsg ("ARM_ERR: payIndexMult is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            payIndexMultProfile = new ARM_ReferenceValue(payIndexMult);
            isPayIndexMult = true;
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

	    // Create the CRA calculator
		ARM_CRALocalCalculator* cfcraCalculator = NULL;
		ARM_CRACalculator* craCalculator = NULL;

		        
        // Create the CRA calculator
        craSpreadCalculator = new ARM_CRASpreadCalculator(ccy,
														  ARM_Date(myStartDate),
														  ARM_Date(myEndDate),
														  payReceive,
														  *notionalProfile,
														  callFreq,
														  callNotice,
														  callCal,
														  *callFeesProfile,
														  fundFreq,
														  *fundSpreadProfile,
														  fundDayCount,
														  cpnDayCount,
														  cpnPayFreq,
														  cpnResetCal,
														  cpnPayCal,
														  payIndex,
														  payIndexResetTiming,
														  *boostedFixProfile,
														  *payIndexMultProfile,
														  *bDownProfile,
														  *bUpProfile,
														  cpnResetFreq,
														  cpnResetTiming,
														  cpnResetGap,
														  refIndex1,
														  coeff1,
														  refIndex2,
														  coeff2);		 		
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

		if (isPayIndexMult)
			delete payIndexMultProfile;
		payIndexMultProfile = NULL;

		if (isbDown)
			delete bDownProfile;
		bDownProfile = NULL;
	
		if (isbUp)
			delete bUpProfile;
		bUpProfile = NULL;
		
		// assign object
		if ( !assignObject( craSpreadCalculator, result, objId ) )
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

		if (isPayIndexMult)
			delete payIndexMultProfile;
		payIndexMultProfile = NULL;

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


long ARMLOCAL_CRASpread_Get(
        const long& CRASpreadId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRASpreadCalculator* craSpreadCalculator;
    ARM_Object* object = NULL;

	try
	{
		craSpreadCalculator = dynamic_cast<ARM_CRASpreadCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(CRASpreadId));

		if (!craSpreadCalculator)
		{
			result.setMsg ("ARM_ERR: CRA Spread Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_OSW_PORT ||
            typeToGet == GC_ACCESS_CSO_PORT ||
			typeToGet == GC_ACCESS_SO_PORT ||
			typeToGet == GC_ACCESS_CMSLONG_PORT ||
			typeToGet == GC_ACCESS_CMSSHORT_PORT ||
            (craSpreadCalculator->IsDbleCorridor() && (typeToGet == GC_ACCESS_SO_PORT ||
			typeToGet == GC_ACCESS_SO_SECOND_PORT)) )
        {
            if(typeToGet == GC_ACCESS_OSW_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetSwoptPF())).Clone();
            else if (typeToGet == GC_ACCESS_CSO_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetCSOPF())).Clone();
            else if (typeToGet == GC_ACCESS_SO_PORT)
				/// SO porffolio related to CMS spread condition
                object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetSOPF(true))).Clone();
            else if (typeToGet == GC_ACCESS_SO_SECOND_PORT)
				/// SO portfolio related to CMS or Libor condition
                object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetSOPF(false))).Clone();
			else if (typeToGet == GC_ACCESS_CMSLONG_PORT)
				object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetCMSLONGPortfolio())).Clone();
			else if (typeToGet == GC_ACCESS_CMSSHORT_PORT)
				object = const_cast< ARM_StdPortfolio& >(*(craSpreadCalculator->GetCMSSHORTPortfolio())).Clone();

			/// Assign the object in ARM cache
			if (!assignObject(object, result, objId))
			{
				delete object;
				return ARM_KO;
			}
			else
			{
				return ARM_OK;
			}
		}
		else
		{
			return ARMLOCAL_GC_Get(CRASpreadId, getType, result, objId);
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}

}

long ARMLOCAL_CRASpread_Set(
        const long& CRASpreadId,
        const long& dataId,
        const string& setType,
		const vector<string>& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId )
{
	ARM_CRASpreadCalculator* craSpreadCalculator = dynamic_cast<ARM_CRASpreadCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(CRASpreadId));

	if (!craSpreadCalculator)
	{
		result.setMsg ("ARM_ERR: CRA Spread Calculator is not of a good type");
		return ARM_KO;
	}

    /// Check for specific accessor
    string typeToSet(stringGetUpper(setType));
	if(typeToSet == GC_ACCESS_CALIB_FLAG)
	{
		ARM_GP_Vector* flags = dynamic_cast<ARM_GP_Vector*>(LOCAL_PERSISTENT_OBJECTS->GetObject(dataId));
		if (!flags || flags->size()==0)
		{
			result.setMsg ("ARM_ERR: flag must be a vector and contain at least one value");
			return ARM_KO;
		}

		/// Reset price if any
		craSpreadCalculator->ResetHasBeenPriced();

		string txt;
		bool calibFlag = (*flags)[0] ? true : false;
		if(craSpreadCalculator->SetSwoptCalib(calibFlag))
		{
			if(calibFlag)
				txt += "Model calib enabled";
			else
			{
				craSpreadCalculator->SetCalibMethod(ARM_CalibMethodPtr(NULL));
				txt += "Model calib disabled";
			}
		}
		else
		{
			craSpreadCalculator->SetSwoptCalib(false);
			result.setMsg ("ARM_ERR: portfolio is missing, calibration disabled");
			return ARM_KO;
		}
		if(flags->size()>1)
		{
			calibFlag = (*flags)[1] ? true : false;
			craSpreadCalculator->SetVolUnSqueezer(calibFlag);
			if(calibFlag)
				txt += ", Vol UnSqueezer enabled";
			else
				txt += ", Vol UnSqueezer disabled";
		}
		if(flags->size()>2)
		{
			calibFlag = (*flags)[2] ? true : false;
			craSpreadCalculator->SetCorrelUnSqueezer(calibFlag);
			if(calibFlag)
				txt += ", Correl UnSqueezer enabled";
			else
				txt += ", Correl UnSqueezer disabled";
		}
		if(flags->size()>3)
		{
			calibFlag = (*flags)[3] ? true : false;
			craSpreadCalculator->SetPVAdjuster(calibFlag);
			if(calibFlag)
				txt += ", PV Adjuster enabled";
			else
				txt += ", PV Adjuster disabled";
		}
		result.setString(txt.c_str());
		return ARM_OK;
	}
	else
		return ARMLOCAL_GC_Set(CRASpreadId,dataId,keys,isUpdated,result,objId);
}

