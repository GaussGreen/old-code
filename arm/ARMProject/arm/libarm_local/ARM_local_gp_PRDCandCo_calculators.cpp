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

#include "ARM_local_gp_PRDCandCo_calculators.h"
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
using ARM::std::vector<double>;
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
//// Function to create a PRDC Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_PRDCCalculator_Create(
        const long& prdcId,
        const long& modelId,
        const vector< long >& otherMktDataIds,
        const vector< double >& schedulerDatas,
        const vector< double >& truncatorDatas,
        const vector< string >& columnsToPrice,
        const vector< bool >& prdcFlags,
        const string& calibType,
        const vector< double >& calibDatas,
		const string& basisTypeStr,
        ARM_result&	result, 
        long        objId )
{
    /// input checks
    if( !GlobalPersistanceOk( result ) )
	    return ARM_KO;

    /// used in the MACRO ARM_RESULT
    CCString msg ("");
    ARM_PowerReverse* powRev;
    ARM_DFBSModel* model;
    ARM_PRDCCalculator* calculator = NULL;

    try
    {
        /// crm tracing
        ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "PRDC Calculator" );

		/// Prdc  Opbget check
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &powRev, prdcId ,"PRDC", result ) ) return ARM_KO;

		///  bi-BsModel
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &model, modelId ,"Model", result ) ) return ARM_KO;

        if(otherMktDataIds.size()<2)
        {
	        result.setMsg ("ARM_ERR: MRS are required in other mkt datas");
	        return ARM_KO;
        }

        ARM_ObjectVector otherMktVector(ARM_PRDCCalculator::NbKeys);
        otherMktVector[ARM_PRDCCalculator::MrsDomKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[0]);
        otherMktVector[ARM_PRDCCalculator::MrsForKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[1]);
        if(otherMktDataIds.size()>=3)
            otherMktVector[ARM_PRDCCalculator::QFxKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[2]);
        if(otherMktDataIds.size()>=4)
            otherMktVector[ARM_PRDCCalculator::LocalFxModelKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[3]);
        if(otherMktDataIds.size()>=5)
            otherMktVector[ARM_PRDCCalculator::MarketIrModelKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[4]);
        if(otherMktDataIds.size()>=6)
            otherMktVector[ARM_PRDCCalculator::QDomKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[5]);
        if(otherMktDataIds.size()>=7)
            otherMktVector[ARM_PRDCCalculator::QForKey] = LOCAL_PERSISTENT_OBJECTS->GetObject(otherMktDataIds[6]);

        std::vector<double> schedulerVector(schedulerDatas);
        std::vector<double> truncatorVector(truncatorDatas);
		std::vector<double> prdcCalibDatas(calibDatas);

        bool markovianDriftSampler  = prdcFlags[0];
        bool fxLocalModelFlag       = prdcFlags[1];
        bool basisIRCalibFlag       = prdcFlags[2];

        ARM_PRDCCalibType prdcCalibType = static_cast< ARM_PRDCCalibType > (ARM_ArgConv_PRDCCalibMode.GetNumber(calibType));

		ARM_BasisType basisType  = (ARM_BasisType)ARM_ArgConv_PRCSBasisType.GetNumber(basisTypeStr);

        calculator = new ARM_PRDCCalculator(powRev,
			model,
			otherMktVector,
			schedulerVector,
			truncatorVector,
			columnsToPrice,
			markovianDriftSampler,
			fxLocalModelFlag,
			prdcCalibType,
			prdcCalibDatas,
			basisIRCalibFlag,
			basisType);

        /// assign object
        if( !assignObject( calculator, result, objId ) ){
	        return ARM_KO; }
        else{
	        return ARM_OK; }

    }

    catch(Exception& x)
    {
	    delete calculator;
	    x.DebugPrint();
	    ARM_RESULT();
    }
}




/////////////////////////////////////////////////////////////////////
/// PRCS calculator creation : extended version
/////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_PRCSCalculator_Create (
        const double&			startDate,															
		const double&			fixEndDate,
        const double&			endDate,
		const string&			domCcy,
		const string&			forCcy,
		const string&			fundCcy,
        const string&			cpnFreqStr,
        const string&			cpnDaycountStr,
		const long&             fxResetGap,  
		const string&			stubRuleStr,
		const string&           resetTimingStr,
        const string&			cpnResetCalStr,
        const string&			cpnPayCalStr,
		const long&				cpnnotionalId,
		const long&				domesticCpnId,
		const long&				foreignCpnId,
		const long&				initialFxId,
        const long&				minCpnId,
        const long&				maxCpnId,
        const string&			fundFreqStr,
        const string&			fundDaycountStr,
		const long&				fundingnotionalId,
		const long&				fundMarginId,
		const long&				fundLvgeId,
        const string&			exerFreqStr,
        const long&				NotifGap,
        const string&			payRecStr,
		const long&				nbNCall,
        const long&				feesId,
		const string&			redempTypeStr,
		const long&				redemptionGap,
		const double&			redemptionStrike, 
		const vector <string >& calibTypesStr,
        const vector< double >&	calibDatasDble,
        const vector< string >&	productsStr,
        const vector< double >&	schedulerDatasDble,
		const vector< double >&	truncatorDatasDble,
        const long &			mktDataManagerId,
        ARM_result&				result, 
        long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_PRDCCalculator* calculator = NULL;
	
	ARM_Curve* cpnnotionalCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldCNP(cpnnotionalCve);

    ARM_Curve* fundingnominalCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFNP(fundingnominalCve);

	ARM_Curve* minCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMCP(minCpnCve);

	ARM_Curve* maxCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMAXCP(maxCpnCve);

	ARM_Curve* domesticCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLLP(domesticCpnCve);

	ARM_Curve* foreignCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLSP(foreignCpnCve);
	
	ARM_Curve* initialFxCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldSP(initialFxCve);

	ARM_Curve* fundMarginCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFMC(fundMarginCve);

	ARM_Curve* fundLvgeCve=new ARM_FlatCurve(1.0);
    ARM_AutoCleaner< ARM_Curve > HoldFLV(fundLvgeCve);

	ARM_Curve* feesCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFEESC(feesCve);
		
	/// Restore Market Data Manager Rep
	ARM_MarketData_ManagerRep* mktDataManager= NULL;
	ARM_AutoCleaner< ARM_MarketData_ManagerRep > HoldMKTDT(mktDataManager);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "PRDC Calculator" );

        // Convert dates
		ARM_Date myStartDate	= ConvertToARMDATE(startDate);
		ARM_Date myFixEndDate	= ConvertToARMDATE(fixEndDate);
		ARM_Date myEndDate		= ConvertToARMDATE(endDate);

        /// Convert curve Id to object if possible
		//notionals
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &cpnnotionalCve, cpnnotionalId ,"cpn notional", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundingnominalCve, fundingnotionalId,"funding notional", result  ) ) return ARM_KO;

		//minCpn
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &minCpnCve, minCpnId ,"minCpn", result ) ) return ARM_KO;

		//maxCpn
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &maxCpnCve, maxCpnId ,"maxCpn", result ) ) return ARM_KO;

		//fundMargin
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundMarginCve, fundMarginId,"fundMargin",result ) ) return ARM_KO;

		//fundLevrage
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &fundLvgeCve, fundLvgeId,"fundLevrage",result ) ) return ARM_KO;

		/// domestic Coupon
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &domesticCpnCve, domesticCpnId, "domestic Coupon", result  ) ) return ARM_KO;

		/// foreign Coupon
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &foreignCpnCve, foreignCpnId , "foreign coupon" , result) ) return ARM_KO;

		/// initial forex 
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &initialFxCve, initialFxId , " initial forex", result ) )  return ARM_KO;

		/// fees
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &feesCve, feesId, "Fees" , result ) ) return ARM_KO;
		
	    /// mkt data manager
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, mktDataManagerId ,"Market Data Manager" , result) ) return ARM_KO;

		/// convert from string to int
		int cpnDaycount = ARM_ArgConv_LgNameDayCount.GetNumber(cpnDaycountStr);
		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(cpnFreqStr);
		int stubRule    = ARM_ArgConv_StubRules.GetNumber(stubRuleStr);
		int resetTiming = ARM_ArgConv_LgTimingMod.GetNumber(resetTimingStr); 


		int fundFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(fundFreqStr);
		int fundDaycount= ARM_ArgConv_LgNameDayCount.GetNumber(fundDaycountStr);
		int exerFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(exerFreqStr);
		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(payRecStr);

		ARM_RedemptionType redemptionType  = (ARM_RedemptionType)ARM_ArgConv_PRCSRedemptionType.GetNumber(redempTypeStr);

		std::vector<double> schedulerDatas(schedulerDatasDble);
		std::vector<double> truncatorDatas(truncatorDatasDble);
		ARM_StringVector productsToPrice(productsStr);
		std::vector<double> calibDatas(calibDatasDble);

		ARM_PRDCCalibType calibType = static_cast< ARM_PRDCCalibType > (ARM_ArgConv_PRDCCalibMode.GetNumber(calibTypesStr[0]));
		bool mkvDriftSplerFlag  = ((calibTypesStr[1] == "Y") || (calibTypesStr[1] == "YES"));
		bool fxLocalModelFlag   = ((calibTypesStr[2] == "Y") || (calibTypesStr[2] == "YES"));
		bool basisIRCalibFlag   = ((calibTypesStr[3] == "Y") || (calibTypesStr[3] == "YES"));

		// defaults
		ARM_Currency domCurrency(domCcy.c_str());
		ARM_Currency forCurrency(forCcy.c_str());
		ARM_Currency fundCurrency(fundCcy.c_str());

		int compFreq = K_COMP_PROP;
		int compType = K_COMP_NONE;

		calculator = new ARM_PRDCCalculator (mktDataManager->GetAsOfDate(),
			myStartDate,
			myFixEndDate,
			myEndDate,
			domCurrency,
			forCurrency,
			fundCurrency,
			cpnDaycount,
			cpnFreq,
			fxResetGap,
			stubRule,
			resetTiming,
            cpnResetCalStr,
            cpnPayCalStr,
			*cpnnotionalCve,
			*domesticCpnCve,
			*foreignCpnCve,
			*minCpnCve,
			*maxCpnCve,
			*initialFxCve,
			fundFreq,
			fundDaycount,
			compFreq,
			compType,
			*fundingnominalCve,
			*fundMarginCve,
			exerFreq,
			NotifGap,
			payRec,
			nbNCall,
			*feesCve,
			redemptionGap,
			redemptionStrike,
			redemptionType,
			productsToPrice,
			fxLocalModelFlag,
			basisIRCalibFlag);

		calculator->Init(*mktDataManager,
				schedulerDatas, 
				truncatorDatas, 
				mkvDriftSplerFlag, 
				calibType, 
				calibDatas);

		/// assign object

		return (assignObject( calculator, result, objId ) )? ARM_OK: ARM_KO;
		
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_KO;
}

/////////////////////////////////////////////////////////////////////
/// PRCS calculator creation : extended version
/////////////////////////////////////////////////////////////////////
extern long ARMLOCAL_PRDKOCalculator_Create (
        const double&			startDate,															
		const double&			fixEndDate,
		const double&			switchDate,
        const double&			endDate,
		const string&			domCcy,
		const string&			forCcy,
		const string&			fundCcy,
        const string&			cpnFreqStr,
        const string&			cpnDaycountStr,
		const long&             fxResetGap,  
		const string&			stubRuleStr,
		const string&           resetTimingStr, 
        const string&			cpnResetCalStr,
        const string&			cpnPayCalStr,
		const long&				cpnnotionalId,
		const long&				domesticCpnId,
		const long&				foreignCpnId,
		const long&				initialFxId,
        const long&				minCpnId,
        const long&				maxCpnId,
		const long&				barrierId,
        const string&			fundFreqStr,
        const string&			fundDaycountStr,
		const long&				fundingnotionalId,
		const long&				fundMarginId,
		const long&				fundLvgeId,
        const string&			exerFreqStr,
        const long&				NotifGap,
        const string&			payRecStr,
		const long&				nbNCall,
        const long&				feesId,
		const string&			redempTypeStr,
		const long&				redemptionGap,
		const double&			redemptionStrike, 
		const vector <string >& calibTypesStr,
        const vector< double >&	calibDatasDble,
        const vector< string >&	productsStr,
        const vector< double >&	schedulerDatasDble,
		const vector< double >&	truncatorDatasDble,
        const long &			mktDataManagerId,
        ARM_result&				result, 
        long					objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_PRDKOCalculator* calculator = NULL;
	
	ARM_Curve* cpnnotionalCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldCNP(cpnnotionalCve);

    ARM_Curve* fundingnominalCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFNP(fundingnominalCve);

	ARM_Curve* minCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMCP(minCpnCve);

	ARM_Curve* maxCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldMAXCP(maxCpnCve);

	ARM_Curve* domesticCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLLP(domesticCpnCve);

	ARM_Curve* foreignCpnCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldLSP(foreignCpnCve);
	
	ARM_Curve* initialFxCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldSP(initialFxCve);

	ARM_Curve* barrierCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldBR(barrierCve);

	ARM_Curve* fundMarginCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFMC(fundMarginCve);

	ARM_Curve* fundLvgeCve=new ARM_FlatCurve(1.0);
    ARM_AutoCleaner< ARM_Curve > HoldFLV(fundLvgeCve);

	ARM_Curve* feesCve=NULL;
    ARM_AutoCleaner< ARM_Curve > HoldFEESC(feesCve);
		
	/// Restore Market Data Manager Rep
	ARM_MarketData_ManagerRep* mktDataManager= NULL;
	ARM_AutoCleaner< ARM_MarketData_ManagerRep > HoldMKTDT(mktDataManager);

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "PRCS Calculator" );

        // Convert dates
		ARM_Date myStartDate	= ConvertToARMDATE(startDate);
		ARM_Date myFixEndDate	= ConvertToARMDATE(fixEndDate);
		ARM_Date mySwitchDate	= ConvertToARMDATE(switchDate);
		ARM_Date myEndDate		= ConvertToARMDATE(endDate);
		

        /// Convert curve Id to object if possible
		//notionals
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &cpnnotionalCve, cpnnotionalId ,"cpn notional", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundingnominalCve, fundingnotionalId,"funding notional", result  ) ) return ARM_KO;

		//minCpn
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &minCpnCve, minCpnId ,"minCpn", result ) ) return ARM_KO;

		//maxCpn
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &maxCpnCve, maxCpnId ,"maxCpn", result ) ) return ARM_KO;

		//fundMargin
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundMarginCve, fundMarginId,"fundMargin",result ) ) return ARM_KO;

		//fundLevrage
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &fundLvgeCve, fundLvgeId,"fundLevrage",result ) ) return ARM_KO;

		/// domestic Coupon
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &domesticCpnCve, domesticCpnId, "domestic Coupon", result  ) ) return ARM_KO;

		/// foreign Coupon
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &foreignCpnCve, foreignCpnId , "foreign coupon" , result) ) return ARM_KO;

		/// initial forex 
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &initialFxCve, initialFxId , " initial forex", result ) )  return ARM_KO;

		/// barrier 
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &barrierCve, barrierId , " barrier", result ) )  return ARM_KO;

		/// fees
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &feesCve, feesId, "Fees" , result ) ) return ARM_KO;
		
	    /// mkt data manager
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, mktDataManagerId ,"Market Data Manager" , result) ) return ARM_KO;

		/// convert from string to int
		int cpnDaycount = ARM_ArgConv_LgNameDayCount.GetNumber(cpnDaycountStr);
		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(cpnFreqStr);
		int stubRule    = ARM_ArgConv_StubRules.GetNumber(stubRuleStr);
		int resetTiming = ARM_ArgConv_LgTimingMod.GetNumber(resetTimingStr); 

		int fundFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(fundFreqStr);
		int fundDaycount= ARM_ArgConv_LgNameDayCount.GetNumber(fundDaycountStr);
		int exerFreq	= ARM_ArgConv_LgNameFrequency.GetNumber(exerFreqStr);
		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(payRecStr);

		ARM_RedemptionType redemptionType  = (ARM_RedemptionType)ARM_ArgConv_PRCSRedemptionType.GetNumber(redempTypeStr);

		std::vector<double> schedulerDatas(schedulerDatasDble);
		std::vector<double> truncatorDatas(truncatorDatasDble);
		ARM_StringVector productsToPrice(productsStr);
		std::vector<double> calibDatas(calibDatasDble);

		ARM_PRDCCalibType calibType = static_cast< ARM_PRDCCalibType > (ARM_ArgConv_PRDCCalibMode.GetNumber(calibTypesStr[0]));
		bool mkvDriftSplerFlag  = ((calibTypesStr[1] == "Y") || (calibTypesStr[1] == "YES"));
		bool fxLocalModelFlag   = ((calibTypesStr[2] == "Y") || (calibTypesStr[2] == "YES"));
		bool basisIRCalibFlag   = ((calibTypesStr[3] == "Y") || (calibTypesStr[3] == "YES"));

		// defaults
		ARM_Currency domCurrency(domCcy.c_str());
		ARM_Currency forCurrency(forCcy.c_str());
		ARM_Currency fundCurrency(fundCcy.c_str());

		int compFreq = K_COMP_PROP;
		int compType = K_COMP_NONE;
		calculator = new ARM_PRDKOCalculator (mktDataManager->GetAsOfDate(),
			myStartDate,
			myFixEndDate,
			myEndDate,
			domCurrency,
			forCurrency,
			fundCurrency,
			cpnDaycount,
			cpnFreq,
			fxResetGap,
			stubRule,
			resetTiming,
            cpnResetCalStr,
            cpnPayCalStr,
			*cpnnotionalCve,
			*domesticCpnCve,
			*foreignCpnCve,
			*minCpnCve,
			*maxCpnCve,
			*initialFxCve,
			*barrierCve,
			ARM_FlatCurve(1),
			fundFreq,
			fundDaycount,
			compFreq,
			compType,
			*fundingnominalCve,
			*fundMarginCve,
			exerFreq,
			NotifGap,
			payRec,
			nbNCall,
			*feesCve,
			redemptionGap,
			redemptionStrike,
			redemptionType,
			productsToPrice,
			fxLocalModelFlag,
			basisIRCalibFlag);

		calculator->Init(*mktDataManager, 
			schedulerDatas, 
			truncatorDatas, 
			mkvDriftSplerFlag, 
			calibType, 
			calibDatas);

		/// assign object

		return (assignObject( calculator, result, objId ) )? ARM_OK: ARM_KO;
		
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
	return ARM_KO;
}

//***********************************************************************
//**********  PRDC Calculator Creating  **********************************
//***********************************************************************

long ARM_PRDCCalculator_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_PRDKOCalculator* calculator = NULL;	
	ARM_Curve* cpnNotionalCv=NULL;
    ARM_Curve* fundNominalCv=NULL;
	ARM_Curve* minCpnCv=NULL;
	ARM_Curve* maxCpnCv=NULL;
	ARM_Curve* domesticCpnCv=NULL;
	ARM_Curve* foreignCpnCv=NULL;
	ARM_Curve* initialFxCv=NULL;
	ARM_Curve * barrierCv = new ARM_FlatCurve(1.0e+15);
	ARM_Curve* NoticeTypeCv = new ARM_FlatCurve(0.0);
	ARM_Curve* fundMarginCv=NULL;
	ARM_Curve* feesCv=NULL;
	ARM_GP_StrVector*  productsToPrice = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "PRDC Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char fixEndDate[20];
		char endDate[20];
		
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("FixEndDate").GetDouble(), fixEndDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);

		ARM_Currency	cpnCurrency(genericParams->GetParamValue("CpnCcy").GetString().c_str());
		ARM_Currency	forCurrency(genericParams->GetParamValue("FgnCcy").GetString().c_str());
		ARM_Currency	fundCurrency(genericParams->GetParamValue("FundCcy").GetString().c_str());

		int cpnFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("CpnFreq").GetString());
		int cpnDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("CpnDayCount").GetString());
		int cpnResetGap = genericParams->GetParamValue("CpnResetGap").GetDouble();
		string resetCal = genericParams->GetParamValue("ResetCal").GetString();
		string payCal	= genericParams->GetParamValue("PayCal").GetString();
		int stubRule	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString());
		int resetTiming	= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("ResetTiming").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &cpnNotionalCv, genericParams->GetParamValue("CpnNominal").GetObjectId() ," cpn Nominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &domesticCpnCv, genericParams->GetParamValue("DomCpn").GetObjectId() ," dom Cpn", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &foreignCpnCv, genericParams->GetParamValue("FgnCpn").GetObjectId() ," fgn Cpn", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &initialFxCv, genericParams->GetParamValue("InitFX").GetObjectId() ," initial FX", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &minCpnCv, genericParams->GetParamValue("MinCpn").GetObjectId() ,"  min Cpn", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &maxCpnCv, genericParams->GetParamValue("MaxCpn").GetObjectId() ," max Cpn", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &barrierCv, genericParams->GetParamValue("Barrier").GetObjectId() ," barrier", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &NoticeTypeCv, genericParams->GetParamValue("NoticeType").GetObjectId() ," Notice Type", result ) ) return ARM_KO;

		int fundFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("FundFreq").GetString());
		int fundDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("FundDayCount").GetString());

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundNominalCv, genericParams->GetParamValue("FundNominal").GetObjectId() ,"  fund Nominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &fundMarginCv, genericParams->GetParamValue("FundSpread").GetObjectId() ," fund Spread", result ) ) return ARM_KO;

		int noticeFreq = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("NoticeFreq").GetString());	
		int noticeGap = genericParams->GetParamValue("NoticeGap").GetDouble();
		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(genericParams->GetParamValue("PayRec").GetString());
		int nbNoCall = genericParams->GetParamValue("NbNoCall").GetDouble();

		ARM_RedemptionType	redemType = (ARM_RedemptionType)ARM_ArgConv_PRCSRedemptionType.GetNumber(genericParams->GetParamValue("RedemType").GetString());
		double  redempStrike = genericParams->GetParamValue("RedemStrike").GetDouble();
		int redemGap = genericParams->GetParamValue("RedemGap").GetDouble();

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &feesCv, genericParams->GetParamValue("Fees").GetObjectId() ," fees", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &productsToPrice, genericParams->GetParamValue("ProductsToPrice").GetObjectId() ," products to price", result ) ) return ARM_KO;

		ARM_StringVector columnsToPrice(productsToPrice->GetValues());

		int yesNo = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("FXLocalFlag").GetString());
		bool fxLocalModelFlag = yesNo ? true:false;
		yesNo = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("BasisIRFlag").GetString());
		bool basisIRCalibFlag = yesNo ? true:false;

		int compFreq = ARM_ArgConv_CompoundingFrequency.GetNumber(genericParams->GetParamValue("CompFreq").GetString());
		int compType = ARM_ArgConv_CompoundingType.GetNumber(genericParams->GetParamValue("CompType").GetString());

		calculator = new ARM_PRDKOCalculator (asOfDate,
			startDate,
			fixEndDate,
			endDate,
			cpnCurrency,
			forCurrency,
			fundCurrency,
			cpnDayCount,
			cpnFreq,
			cpnResetGap,
			stubRule,
			resetTiming,
            resetCal,
            payCal,
			*cpnNotionalCv,
			*domesticCpnCv,
			*foreignCpnCv,
			*minCpnCv,
			*maxCpnCv,
			*initialFxCv,
			*barrierCv,
			*NoticeTypeCv,
			fundFreq,
			fundDayCount,
			compFreq,
			compType,
			*fundNominalCv,
			*fundMarginCv,
			noticeFreq,
			noticeGap,
			payRec,
			nbNoCall,
			*feesCv,
			redemGap,
			redempStrike,
			redemType,
			columnsToPrice,
			fxLocalModelFlag,
			basisIRCalibFlag);
		
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
//**********  CCS Calculator initializing  **********************************
//***********************************************************************

long ARM_PRDCCalculator_InitFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_PRDKOCalculator*	calculator = NULL;
	ARM_PRDKOCalculator*	newCalculator = NULL;
	
	ARM_MarketData_ManagerRep*	mktDataManager = NULL;
	std::vector<double>&	scheduler = NULL;
	std::vector<double>&	truncator  = NULL;
	std::vector<double>&	calibDatas = NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "PRDC Calculator Init" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calculator, genericParams->GetParamValue("PRDCId").GetObjectId() ," PRDC Calculator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, genericParams->GetParamValue("MktDataMger").GetObjectId() ," Mkt Data Manager", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &scheduler, genericParams->GetParamValue("SchedulerTree").GetObjectId() ," Tree 3F Scheduler", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &truncator, genericParams->GetParamValue("TruncatorTree").GetObjectId() ," Tree 3F Truncator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calibDatas, genericParams->GetParamValue("CalibDatas").GetObjectId() ," Datas of calibration", result ) ) return ARM_KO;

		ARM_PRDCCalibType calibType = (ARM_PRDCCalibType)ARM_ArgConv_PRDCCalibMode.GetNumber(genericParams->GetParamValue("CalibType").GetString());

		int yesNo = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("MarkovianDrift").GetString());
		bool markovianDriftSamplerFlag = yesNo ? true:false;
	
		newCalculator = (ARM_PRDKOCalculator*)calculator->Clone();
		newCalculator->Init(*mktDataManager,
							*scheduler,
							*truncator,
							markovianDriftSamplerFlag,
							calibType,
							*calibDatas);

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

///////////////////////////////////////////////
//// Function to set data to a PRDC Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_PRDC_Set(
        long prdcId,
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
    ARM_PRDCalculator* oldPRDC=NULL;
    ARM_PRDCalculator* newPRDC=NULL;
    ARM_Object* object=NULL;
    string updateMsg("update is done!!! ");


	try
	{
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &oldPRDC, prdcId ," prdc ", result ) ) return ARM_KO;

		ARM_StdPortfolio* pf = NULL;
		if( GetObjectFromIdWithDynamicCastCheckwNullandMsge( &pf, dataId,"portfolio",result ) ) return ARM_KO;
        if(pf)
        {
            newPRDC = (isUpdated) ? oldPRDC : dynamic_cast<ARM_PRDCalculator *>(oldPRDC->Clone());

            string typeToSet(stringGetUpper(setType));
            bool isUpdateOswStrike=true;
			bool updateCalibration = true;
            if(typeToSet == GC_ACCESS_DOM_OSW_PORT)
            {
                newPRDC->SetOSWPortfolio(*pf,ARM_PRDCCalculator::OswDomModelKey);
            }
            else if(typeToSet == GC_ACCESS_FOR_OSW_PORT)
            {
                newPRDC->SetOSWPortfolio(*pf,ARM_PRDCCalculator::OswForModelKey);
            }
			else if(typeToSet == GC_ACCESS_FX_PORT)
            {
                newPRDC->SetFxPortfolio(*pf);
            }
		    else
                updateMsg += "Portfolio setting failed";

            /// Keep calibration data consistency
            if(updateCalibration)
				newPRDC->UpdateCalibrationAndTimeIt(isUpdateOswStrike);

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }
			
			return  !assignObject( newPRDC, result, objId ) ? ARM_KO: ARM_OK; 

        }
		
		return ARMLOCAL_GC_Set(prdcId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newPRDC;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to get data from a PRDC Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_PRDC_Get(
        const long& prdcId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_PRDCalculator* calculator;
    ARM_Object* object=NULL;

	try
	{
		calculator = dynamic_cast<ARM_PRDCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(prdcId));

		if (!calculator)
		{
			result.setMsg ("ARM_ERR: PRDC Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_DOM_OSW_PORT		||
			typeToGet == GC_ACCESS_FOR_OSW_PORT		||
            typeToGet == GC_ACCESS_FX_PORT			||
			typeToGet == GC_ACCESS_EXTRA_FX_PORT	||
            typeToGet == GC_ACCESS_FLOORED_FX_PORT	||
			typeToGet == GC_ACCESS_CAPPED_FX_PORT	|| 
			typeToGet == GC_ACCESS_REDEMPTION_FX_PORT ||
			typeToGet == GC_ACCESS_POWER_REVERSE      ||
			typeToGet == GC_ACCESS_DFBSMODEL          ||
			typeToGet == GC_ACCESS_EXTRA_FX_CALIB)
        {
            if(typeToGet == GC_ACCESS_DOM_OSW_PORT)
                object = calculator->GetOSWPortfolio(ARM_PRDCCalculator::OswDomModelKey)->Clone();

            else if(typeToGet == GC_ACCESS_FOR_OSW_PORT)
                object = calculator->GetOSWPortfolio(ARM_PRDCCalculator::OswForModelKey)->Clone();

            else if(typeToGet == GC_ACCESS_FX_PORT)
                object = calculator->GetFxPortfolio()->Clone();

            else if(typeToGet == GC_ACCESS_EXTRA_FX_PORT)
                object = calculator->GetExtraFxPortfolio()->Clone();

			else if(typeToGet == GC_ACCESS_POWER_REVERSE)
                object = calculator->GetPowerReverseSwap()->Clone();

			else if(typeToGet == GC_ACCESS_DFBSMODEL)
                object = calculator->GetAnalyticalModel()->Clone();


            else if(typeToGet == GC_ACCESS_FLOORED_FX_PORT)
            {
                if(calculator->GetFlooredFxPortfolio())
                    object = const_cast<ARM_StdPortfolio*>(calculator->GetFlooredFxPortfolio())->Clone();
                else
                {
			        result.setMsg ("ARM_ERR: No floored coupon");
			        return ARM_KO;
                }
            }

            else if(typeToGet == GC_ACCESS_CAPPED_FX_PORT)
            {
                if(calculator->GetCappedFxPortfolio())
                    object = const_cast<ARM_StdPortfolio*>(calculator->GetCappedFxPortfolio())->Clone();
                else
                {
			        result.setMsg ("ARM_ERR: No Fx capped coupon");
			        return ARM_KO;
                }
            }

            else if(typeToGet == GC_ACCESS_REDEMPTION_FX_PORT)
            {
                if(calculator->GetRedemptionFxPortfolio())
                    object = const_cast<ARM_StdPortfolio*>(calculator->GetRedemptionFxPortfolio())->Clone();
                else
                {
			        result.setMsg ("ARM_ERR: No Fx terminal redemption coupon");
			        return ARM_KO;
                }
            }

            else if(typeToGet == GC_ACCESS_EXTRA_FX_CALIB)
            {
                if(calculator->GetExtraFxCalib())
                    object = const_cast<ARM_CalibMethod*>(calculator->GetExtraFxCalib())->Clone();
                else
                {
			        result.setMsg ("ARM_ERR: No Extra Fx calib method");
			        return ARM_KO;
                }
            }
			else
            {
			    result.setMsg ("ARM_ERR: Unknown type to get in PRDC calculator");
			    return ARM_KO;
            }

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
            return ARMLOCAL_GC_Get(prdcId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//***********************************************************************
//**********  CCS Calculator Creating  **********************************
//***********************************************************************

long ARM_CCSCalculator_CreateFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_CCSCalculator*	calculator = NULL;

	ARM_Curve*	domNominalCv = NULL;
	ARM_Curve*	domSpreadCv  = NULL;
	ARM_Curve*	forNominalCv = NULL;
	ARM_Curve*	forSpreadCv  = NULL;
	ARM_Curve*	feesCv  = NULL;
	ARM_GP_StrVector*  productsToPrice = NULL;


	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "CCS Calculator Create" );

		char asOfDate[20];
		char startDate[20];
		char endDate[20];
		
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("AsOfDate").GetDouble(), asOfDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("StartDate").GetDouble(), startDate);
		Local_XLDATE2ARMDATE(genericParams->GetParamValue("EndDate").GetDouble(), endDate);

		ARM_Currency	domCurrency(genericParams->GetParamValue("DomCcy").GetString().c_str());
		int domFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("DomFreq").GetString());
		int domDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("DomDayCount").GetString());

		
		ARM_Currency	forCurrency(genericParams->GetParamValue("ForCcy").GetString().c_str());
		int forFreq     = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("ForFreq").GetString());
		int forDayCount	= ARM_ArgConv_LgNameDayCount.GetNumber(genericParams->GetParamValue("ForDayCount").GetString());

		int fxResetGap = genericParams->GetParamValue("FxResetGap").GetDouble();
		string resetCal = genericParams->GetParamValue("ResetCal").GetString();
		string payCal = genericParams->GetParamValue("PayCal").GetString();
		string forResetCal = genericParams->GetParamValue("ForResetCal").GetString();
		string forPayCal = genericParams->GetParamValue("ForPayCal").GetString();


		int payRec      = ARM_ArgConv_RcvOrPay.GetNumber(genericParams->GetParamValue("PayRec").GetString());
		int noticeFreq = ARM_ArgConv_LgNameFrequency.GetNumber(genericParams->GetParamValue("NoticeFreq").GetString());
		int noticeGap = genericParams->GetParamValue("NoticeGap").GetDouble();
		int nbNoCall = genericParams->GetParamValue("NbNoCall").GetDouble();

		int stubRule	= ARM_ArgConv_StubRules.GetNumber(genericParams->GetParamValue("StubType").GetString());
		int resetTiming	= ARM_ArgConv_Timing.GetNumber(genericParams->GetParamValue("ResetTiming").GetString());
		int intRule		= ARM_ArgConv_IntRule.GetNumber(genericParams->GetParamValue("IntRule").GetString());	

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &domNominalCv, genericParams->GetParamValue("DomNominal").GetObjectId() ," dom Nominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &domSpreadCv, genericParams->GetParamValue("DomSpread").GetObjectId() ," dom Spread", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &forNominalCv, genericParams->GetParamValue("ForNominal").GetObjectId() ," for Nominal", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &forSpreadCv, genericParams->GetParamValue("ForSpread").GetObjectId() ," for Spread", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &feesCv, genericParams->GetParamValue("Fees").GetObjectId() ," fees", result ) ) return ARM_KO;

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &productsToPrice, genericParams->GetParamValue("ProductsToPrice").GetObjectId() ," products to price", result ) ) return ARM_KO;

        calculator = new ARM_CCSCalculator(	asOfDate,
												startDate,
												startDate,
												endDate,
												domCurrency,
												forCurrency,
												domDayCount,
												domFreq,
												resetCal,
												payCal,
												forFreq,
												forDayCount,
												forResetCal,
												forPayCal,
												fxResetGap,
												stubRule,											
												*domSpreadCv,
												*domNominalCv,												
												*forSpreadCv,
												*forNominalCv,
												noticeFreq,
												noticeGap,
												payRec,
												nbNoCall,
												*feesCv,
												productsToPrice->GetValues());
		
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
//**********  CCS Calculator initializing  **********************************
//***********************************************************************

long ARM_CCSCalculator_InitFunctor::operator()( ARM_result& result, long objId )
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenericParams* genericParams = GetGenericParams();

	ARM_CCSCalculator*	calculator = NULL;
	ARM_CCSCalculator*	newCalculator = NULL;

	ARM_MarketData_ManagerRep*	mktDataManager = NULL;
	std::vector<double>&	scheduler = NULL;
	std::vector<double>&	truncator  = NULL;
	std::vector<double>&	calibDatas = NULL;
	string message("Done! ");

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "CCS Calculator Init" );

		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calculator, genericParams->GetParamValue("CCSCalculator").GetObjectId() ," CCS Calculator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &mktDataManager, genericParams->GetParamValue("MktDataMger").GetObjectId() ," Mkt Data Manager", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &scheduler, genericParams->GetParamValue("SchedulerTree").GetObjectId() ," Tree 3F Scheduler", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &truncator, genericParams->GetParamValue("TruncatorTree").GetObjectId() ," Tree 3F Truncator", result ) ) return ARM_KO;
		if( GetObjectFromIdWithDynamicCastCheckandMsge( &calibDatas, genericParams->GetParamValue("CalibDatas").GetObjectId() ," Datas of calibration", result ) ) return ARM_KO;

		ARM_PRDCCalibType calibType = (ARM_PRDCCalibType)ARM_ArgConv_PRDCCalibMode.GetNumber(genericParams->GetParamValue("CalibType").GetString());

		int yesNo = ARM_ArgConv_YesNo.GetNumber(genericParams->GetParamValue("MarkovianDrift").GetString());
		bool markovianDriftSamplerFlag = yesNo? true:false;
	
		newCalculator = (ARM_CCSCalculator*)calculator->Clone();

		newCalculator->Init(*mktDataManager,
							*scheduler,
							*truncator,
							markovianDriftSamplerFlag,
							calibType,
							*calibDatas);

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

