
#include "firstToBeIncluded.h"

//#include "ARM_local_gp_calculators.h"
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

////////////////////////////////////////////
//// Function to create a date strip combiner
////////////////////////////////////////////
extern long ARMLOCAL_DateStripCombiner_Create(
         const VECTOR<long>&    C_dateStripIds,
		 const string& funcToMerge,
         ARM_result&	result, 
         long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_DateStripCombiner* dateStripCombiner= NULL;

	try
	{
        ARM_DateStripVector dateStripVec;
		ARM_DateStrip* dateStrip;
        size_t dateStripSize = C_dateStripIds.size();

        for(size_t i=0; i<dateStripSize; ++i)
        {
	        if( !GetObjectFromId( &dateStrip, C_dateStripIds[i], ARM_DATESTRIP) )
	        {
   		        result.setMsg ("ARM_ERR: date strip is not of a good type");
		        return ARM_KO;
	        }
            dateStripVec.push_back(dateStrip);
        }

        dateStripCombiner = new ARM_DateStripCombiner(dateStripVec, funcToMerge );  

		/// assign object
		if( !assignObject( dateStripCombiner, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}

	catch(Exception& x)
	{
		delete dateStripCombiner;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_DateStripCombiner_GetData(
         long dateStripCombinerId,
	     long dataType,
		 size_t dateStripNb,
		 VECTOR<double>& Data,
		 ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_DateStripCombiner* dateStripCombiner= NULL;
		if( !GetObjectFromId( &dateStripCombiner, dateStripCombinerId, ARM_DATESTRIP_COMBINER ) )
		{
			result.setMsg ("ARM_ERR: dateStrip combiner is not of a good type");
			return ARM_KO;
		};

		if( dateStripNb >= dateStripCombiner->size() )
		{
			char msg[255];
			sprintf(msg,"ARM_ERR: trying to access a datestrip outside the range total size is %d while asked %d",dateStripCombiner->size(), dateStripNb );
			result.setMsg(CCString(msg));
			return ARM_KO;
		};

		ARM_GP_Vector* tmpData = dateStripCombiner->GetDateStrip(dateStripNb)->GetMemberData(dataType);
		size_t i;

		/// test whether we have dates in which case, we convert this!
		if(		dataType == K_START_DATES 
			||	dataType == K_END_DATES 
			||  dataType == K_RESET_DATES
			||	dataType == K_PAY_DATES 
			||	dataType == K_FWD_START_DATES
			||	dataType == K_FWD_END_DATES 
			)
		{
			for( i=0; i<tmpData->size(); ++i)
                if( (*tmpData)[i] != ARM_DateStripCombiner::DateStripCombiner_BlankData )
					Data.push_back( JulianToXLDate((*tmpData)[i]) );
				else
					Data.push_back( (*tmpData)[i] );
		}
		/// otherwise no treatment!
		else
			for( i=0; i<tmpData->size(); ++i)
				Data.push_back( (*tmpData)[i] );

		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}



extern long ARMLOCAL_DateStripCombiner_GetMergeData(
	long dateStripCombinerId,
	VECTOR<double>& Data,
    ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result )  )
		return ARM_KO;
	
	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	
	try
	{
		ARM_DateStripCombiner* dateStripCombiner= NULL;
		if( !GetObjectFromId( &dateStripCombiner, dateStripCombinerId, ARM_DATESTRIP_COMBINER ) )
		{
			result.setMsg ("ARM_ERR: dateStrip combiner is not of a good type");
			return ARM_KO;
		};

		ARM_GP_Vector tmpData = *dateStripCombiner->GetMergeData();
		Data.reserve(tmpData.size());

		for( size_t i=0; i<tmpData.size(); ++i)
			Data.push_back( JulianToXLDate(tmpData[i]) );
		return ARM_OK;
	}

	catch(Exception& x)
	{
		x.DebugPrint();
		
		ARM_RESULT();
	}
}




///////////////////////////////////////////////
//// Function to get Currencies Name from a Calculator
///////////////////////////////////////////////

long ARMLOCAL_ARM_GetCcyFromGenCalculator(long calcId,
        const CCString& C_ccyType,
        ARM_result& result, 
        long objId)
{
    CCString msg(""); // used in macro ARM_RESULT()
    ARM_GenCalculator* calculator = NULL;
    
	try
	{
		calculator = (ARM_GenCalculator*) LOCAL_PERSISTENT_OBJECTS->GetObject(calcId);
		if ( !(LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(calculator, ARM_CALLREVFLOATER) == 1 ))
        {
			result.setMsg("ARM_ERR: ARM_GetDomCcy: unknown generic calculator");
			return ARM_KO;
        }
         string ccyTypeStr = CCSTringToSTLString(C_ccyType);
         ARM_CcyType ccyType = (ARM_CcyType) ARM_ArgConv_GenCalculatorCcyType.GetNumber(ccyTypeStr);
         CCString CcyStr;
         switch(ccyType)
         {
         case ARM_GenCalculatorCcyType::DomCurrency:
             CcyStr = calculator->GetDomesticCcy().GetCcyName();
             break;
         case ARM_GenCalculatorCcyType::ForCurrency:
             CcyStr = calculator->GetForeignCcy().GetCcyName();
             break;
         case ARM_GenCalculatorCcyType::FundCurrency:
             CcyStr = calculator->GetFundingCcy().GetCcyName();
             break;
         case ARM_GenCalculatorCcyType::CpnCurrency:
             CcyStr = calculator->GetCurrencyUnit()->GetCcyName();
             break;
         default:
             {
                 result.setMsg("Ccy Type: Unknown type");
                 return ARM_KO;
             }
         }

		result.setString(CcyStr);
		return ARM_OK;
	}
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////////////////
//// To initialize a crf calculator with a mrkDataManager
///////////////////////////////////////////////
extern long ARMLOCAL_Calculator_Initialize(
        const long& gcId,
        const long& mktMangerId,
        const CCString& toCalSigma,
        const CCString& toCalMrs,
		const CCString& strikeTypeToCalMrs,
        const CCString& toAdjKcap,
        const CCString& toAdjKfloor,
        const CCString& modelTypeCCStr,
        const CCString& toCalKskew,
        const double&   kShift,
		const long&     frontier,  
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CRFCalculator* crfCalculator;
	ARM_CRFCalculator* NewcrfCalculator;
	ARM_MarketData_ManagerRep* mktManager;
    ARM_Object* object=NULL;

	try
	{
		crfCalculator = dynamic_cast<ARM_CRFCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(gcId));

		if (!crfCalculator)
		{
			result.setMsg ("ARM_ERR: CRF Calculator is not of a good type");
			return ARM_KO;
		}

		mktManager = dynamic_cast<ARM_MarketData_ManagerRep *>(LOCAL_PERSISTENT_OBJECTS->GetObject(mktMangerId));

		if (!mktManager)
		{
			result.setMsg ("ARM_ERR: Mkt Data Manager is not of a good type");
			return ARM_KO;
		}
		NewcrfCalculator = (ARM_CRFCalculator*)crfCalculator->Clone();

		CCString sigmaTypeCCStr = (toCalSigma=="N" || toCalSigma=="NO") ?  CCString("UNKNOWN"):
							    (toCalSigma=="Y" || toCalSigma=="YES")  ? CCString("EQUIVALENT"):toCalSigma;
		string sigmaTypeStr = CCSTringToSTLString(sigmaTypeCCStr);
		ARM_SigmaCalibType sigmaType = (ARM_SigmaCalibType) ARM_ArgConv_SigmaCalibType.GetNumber(sigmaTypeStr);


		CCString mrsTypeCCStr = (toCalMrs=="N" || toCalMrs=="NO") ?  CCString("UNKNOWN"):
							    (toCalMrs=="Y" || toCalMrs=="YES")  ? CCString("STM_FIRST_COLUMN"):toCalMrs;
		string mrsTypeStr = CCSTringToSTLString(mrsTypeCCStr);
		ARM_MRSCalibType mrsType = (ARM_MRSCalibType) ARM_ArgConv_MRSCalibType.GetNumber(mrsTypeStr);

		string mrsStrikeTypeStr = CCSTringToSTLString(strikeTypeToCalMrs);
		ARM_MRSStrikeCalibType mrsStrikeType = (ARM_MRSStrikeCalibType) ARM_ArgConv_MRSStrikeCalibType.GetNumber(mrsStrikeTypeStr);

		

        bool Bool_toAdjKcap     =(toAdjKcap=="Y" || toAdjKcap=="YES");
        bool Bool_toAdjKfloor   =(toAdjKfloor=="Y" || toAdjKfloor=="YES");
        bool Bool_toCalKskew    =(toCalKskew=="Y" || toCalKskew=="YES");
        string modelTypeStr = CCSTringToSTLString(modelTypeCCStr);
        ARM_ModelType modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(modelTypeStr);

        NewcrfCalculator->SetModelType(modelType);
        NewcrfCalculator->SetOSWCalibFlag(sigmaType);
		NewcrfCalculator->SetSTMCalibFlag(mrsType);
		NewcrfCalculator->SetMRSStrikeType(mrsStrikeType);
		NewcrfCalculator->SetCapCalibFlag(Bool_toAdjKcap);
		NewcrfCalculator->SetFloorCalibFlag(Bool_toAdjKfloor);
		NewcrfCalculator->SetSkewCalibFlag(Bool_toCalKskew);
		NewcrfCalculator->SetNbIterFrontier(frontier);

        
        /// it's necessairy to set all flags before in order to update
        /// model, calibrator and Gensecurity
		ARM_MarketData_ManagerRep* mktManagerCloned = static_cast <ARM_MarketData_ManagerRep *>(mktManager->Clone());
		CC_NS(std,auto_ptr)<ARM_MarketData_ManagerRep> mktManagerPtr(mktManagerCloned);
        NewcrfCalculator->Initialize(mktManagerCloned);

	
        /// Assign the object in the ARM cache
		if( !assignObject( NewcrfCalculator, result, objId ) )
        {
            delete object;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}

extern string GCGetTypeToClass(const string& typeToGet, long gcId)
{
	string typeToGetUpper = ARM::stringGetUpper(typeToGet);
    bool gensecType = (typeToGetUpper == GC_ACCESS_SECURITY || typeToGetUpper == GC_ACCESS_CAP    ||
        typeToGetUpper == GC_ACCESS_FLOOR    || typeToGetUpper ==    GC_ACCESS_FUNDING            ||
        typeToGetUpper == GC_ACCESS_STDSWAP  || typeToGetUpper == GC_ACCESS_STDLEG                ||
        typeToGetUpper == GC_ACCESS_RFSWAP   || typeToGetUpper == GC_ACCESS_RFLEG                 ||
        typeToGetUpper == GC_ACCESS_BERMUDA );

    if(gensecType)
        return LOCAL_GENSEC_CLASS;

    else if(typeToGetUpper == GC_ACCESS_MODEL)
    {
        /// Find the model name for interface
        ARM_GenCalculator* genCalculator = dynamic_cast<ARM_GenCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(gcId));
		if (!genCalculator)
			return LOCAL_ANY_CLASS;

		ARM_PricingModel* model = dynamic_cast<ARM_PricingModel *>(&(*genCalculator->GetPricingModel()));
		if (!model)
			return LOCAL_ANY_CLASS;

		if( dynamic_cast<ARM_HullWhite1F*>(model) )
			return LOCAL_HW1FMOD_CLASS;
		
		else if( dynamic_cast<ARM_HullWhite2F*>(model) )
        	return LOCAL_HW2FMOD_CLASS;
		
		else if( dynamic_cast<ARM_SFRM*>(model) )
            return LOCAL_SFRMMODEL_CLASS;
		
		else if( dynamic_cast<ARM_HybridBasisFwdIR*>(model) )
            return LOCAL_BASISFWDIRMOD_CLASS;

		else
            /// New style to get the exportable 
            return model->GetExportShortName();
    }
    else if(typeToGetUpper == GC_ACCESS_CALIB ||
			typeToGetUpper == GC_ACCESS_EXTRA_FX_CALIB)
        return LOCAL_CALIBMETHOD_CLASS;

    else if(typeToGetUpper == GC_ACCESS_OSW_PORT            ||
            typeToGetUpper == GC_ACCESS_SHORT_TERM_PORT     ||
            typeToGetUpper == GC_ACCESS_DOM_OSW_PORT        || 
            typeToGetUpper == GC_ACCESS_FOR_OSW_PORT        ||
            typeToGetUpper == GC_ACCESS_FX_PORT             ||
            typeToGetUpper == GC_ACCESS_FLOORED_FX_PORT     ||
            typeToGetUpper == GC_ACCESS_CAPPED_FX_PORT      || 
            typeToGetUpper == GC_ACCESS_REDEMPTION_FX_PORT  ||
            typeToGetUpper == GC_ACCESS_EXTRA_FX_PORT       ||
            typeToGetUpper == GC_ACCESS_CF_PORT             || 
            typeToGetUpper == GC_ACCESS_OSW_SECOND_PORT     ||
			typeToGetUpper == GC_ACCESS_SO_PORT             || 
            typeToGetUpper == GC_ACCESS_SO_SECOND_PORT      || 
            typeToGetUpper == GC_ACCESS_CMSLONG_PORT        || 
            typeToGetUpper == GC_ACCESS_CMSSHORT_PORT       ||
			typeToGetUpper == GC_ACCESS_SO_PORT             || 
            typeToGetUpper == GC_ACCESS_SO_SECOND_PORT      || 
            typeToGetUpper == GC_ACCESS_CMSLONG_PORT        || 
            typeToGetUpper == GC_ACCESS_CMSSHORT_PORT       ||
			typeToGetUpper == GC_ACCESS_CSO_PORT            || 
            typeToGetUpper == GC_ACCESS_STD_PORT)
        return LOCAL_PF_CLASS;
	
	else if(typeToGetUpper == GC_ACCESS_OPT_PORT)
		return LOCAL_OPTION_PORTFOLIO_CLASS;

	else if(typeToGetUpper == GC_ACCESS_POWER_REVERSE)
		return LOCAL_POWER_REVERSE_CLASS;

	else if(typeToGetUpper == GC_ACCESS_DFBSMODEL)
		return LOCAL_DFFXBS_CLASS;

	else if(typeToGetUpper == GC_ACCESS_LOCAL_PORT)
		return LOCAL_LOCAL_PORTFOLIO_CLASS;

    else  if(typeToGetUpper == GC_ACCESS_SO_CSTMANAGER)
		return LOCAL_CSTMANAGER_CLASS;

    else if(typeToGetUpper == GC_ACCESS_MDM)
        return LOCAL_MKTDATAMANAGER_CLASS;

	else if(typeToGetUpper == GC_ACCESS_OSW_BSMODEL)
		return LOCAL_BSMODEL_CLASS;

	else if(typeToGetUpper == GC_ACCESS_EXO_SWAPTION)
		return LOCAL_SWAPTION_CLASS;
	else if(typeToGetUpper == GC_ACCESS_EXOSWAP)
		return LOCAL_SWAP_CLASS;
    else
        return LOCAL_ANY_CLASS;
}


///////////////////////////////////////////////
//// Function to get data from a Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_GC_Get(
        long gcId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_GenCalculator* genCalculator;
    ARM_Object* object=NULL;

	try
	{
		genCalculator = dynamic_cast<ARM_GenCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(gcId));

		if (!genCalculator)
		{
			result.setMsg ("ARM_ERR: Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));

        bool gensecType = (typeToGet == GC_ACCESS_SECURITY || typeToGet == GC_ACCESS_CAP    ||
        typeToGet == GC_ACCESS_FLOOR    || typeToGet == GC_ACCESS_FUNDING                   ||
        typeToGet == GC_ACCESS_STDSWAP  || typeToGet == GC_ACCESS_STDLEG                    ||
        typeToGet == GC_ACCESS_RFSWAP   || typeToGet == GC_ACCESS_RFLEG                     ||
        typeToGet == GC_ACCESS_BERMUDA );

        if(gensecType)
        {
            ARM_GenSecurity* genSec = NULL;

            if( typeToGet == GC_ACCESS_SECURITY)
                    genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetGenSecurity()->Clone());

            else if( typeToGet == GC_ACCESS_CAP)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::Caplet));

            else if( typeToGet == GC_ACCESS_FLOOR)
                    genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::Floorlet));

            else if( typeToGet == GC_ACCESS_FUNDING)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::Funding));

            else if( typeToGet == GC_ACCESS_STDLEG)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::StdFlow));

             else if( typeToGet == GC_ACCESS_RFLEG)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::RFFlow));

            else if( typeToGet == GC_ACCESS_STDSWAP)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::StdSwaplet));

            else if( typeToGet == GC_ACCESS_RFSWAP)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::RFSwaplet));

            else if( typeToGet == GC_ACCESS_BERMUDA)
                genSec = static_cast< ARM_GenSecurity* >(genCalculator->GetSubGenSecurity(ARM_CRFCalculator::StdBermudaPrice));
            else
            {
                result.setMsg ("ARM_ERR: Data type is not of a good type");
			    return ARM_KO;
            }

            object = (ARM_Object*)genSec;
        }
        else if(typeToGet == GC_ACCESS_MODEL)
        {
            ARM_PricingModel* pmod = static_cast< ARM_PricingModel* >(genCalculator->GetPricingModel()->Clone());
            object = (ARM_Object*)pmod;
        }
        else if(typeToGet == GC_ACCESS_CALIB)
        {
            ARM_CalibMethod* calibMeth = static_cast< ARM_CalibMethod* >(genCalculator->GetCalibMethod()->Clone());
            object =  static_cast< ARM_Object* >(calibMeth);
        }
        else if(typeToGet == GC_ACCESS_MDM)
        {
            ARM_MarketData_ManagerRep* mdm = static_cast< ARM_MarketData_ManagerRep* >(genCalculator->GetMktDataManager()->Clone());
            object =  static_cast< ARM_Object* >(mdm);
        }
        else
		{
			result.setMsg ("ARM_ERR: Accessor type not recognized for a Calculator");
            delete object;
			return ARM_KO;
		}

        /// Assign the object in the ARM cache
		if( !assignObject( object, result, objId ) )
        {
            delete object;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data from a Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_GC_Set(
        long calculatorId,
        long dataId,
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
	ARM_GenCalculator* calculator = NULL;
    ARM_GenCalculator* newGC=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");

	try
	{
		if( !GetObjectFromIdWithDynamicCastCheck( &calculator, calculatorId  ) )
		{
			result.setMsg ("ARM_ERR: maxCpn  is not of a good type");
			return ARM_KO;
		}
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        if(isUpdated)
            /// Just update the initial GC
            newGC=calculator;
        else
            /// Clone the initial GC then set
            newGC = (ARM_GenCalculator*) calculator->Clone();

        ARM_GenSecurity* genSec;
        ARM_PricingModel* pmod;
        ARM_CalibMethod* calibMeth;
        ARM_MarketData_ManagerRep* mdm;

        if((genSec = dynamic_cast< ARM_GenSecurity* >(object)))
        {
            newGC->SetGenSecurity( ARM_GenSecurityPtr((ARM_GenSecurity*)genSec->Clone()) );

            /// Re-build auto-calibration portfolios
            newGC->CreateAndSetCalibrationAndTimeIt();

            updateMsg += "generic security";
        }
        else if((pmod = dynamic_cast< ARM_PricingModel* >(object)))
        {
            newGC->SetPricingModel( ARM_PricingModelPtr(static_cast< ARM_PricingModel* >(pmod->Clone())) );
            updateMsg += "pricing model";
        }
        else if((calibMeth = dynamic_cast< ARM_CalibMethod* >(object)))
        {
            newGC->SetCalibMethod( ARM_CalibMethodPtr(static_cast< ARM_CalibMethod* >(calibMeth->Clone())) );

            /// Update internal calibration datas
            newGC->UpdateCalibrationAndTimeIt();

            updateMsg += "calib method";
        }
        else if((mdm = dynamic_cast< ARM_MarketData_ManagerRep* >(object)))
        {
			if(isUpdated)
			{
				newGC->UpDateWithSubMktDataManager(*mdm);
				updateMsg += "market data manager by updating only the existing datas";
			}
			else
			{
				newGC->SetMktDataManager( ARM_MarketData_ManagerRepPtr(static_cast< ARM_MarketData_ManagerRep* >(mdm->Clone())) );
				if (keys.size() != 0)
					newGC->SetKeys(keys);
				 /// Update internal datas depending on MDM
				newGC->Update();
				updateMsg += "market data manager using clone of mkt";
			}
        }
        else
		{
			result.setMsg ("ARM_ERR: Objet not settable for a Calculator");
            delete newGC;
			return ARM_KO;
		}

        if(isUpdated)
        {
             result.setString(updateMsg.c_str());
			return ARM_OK;
        }

        /// Assign the new GC in the ARM cache
		else if(!assignObject( newGC, result, objId ) )
        {
            delete newGC;
			return ARM_KO;
        }
		else
			return ARM_OK;
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newGC;
		x.DebugPrint();
		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to get pricing data from a Gen Calculator
////////////////////////////////////////////

long ARMLOCAL_GC_GetPricingData(
	long gcId,
	const string& key,
	ARM_GramFctorArg& argResult,
	ARM_result& result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	try
	{
		ARM_GenCalculator* maturityCapCalc = NULL;

		ARM_Object* armObj =  LOCAL_PERSISTENT_OBJECTS->GetObject(gcId);
		ARM_GenCalculator* genCalc = NULL;
		if( !(genCalc = dynamic_cast< ARM_GenCalculator* >( armObj )) )
		{
			result.setMsg ("ARM_ERR: calculator is not of a good type");
			return ARM_KO;
		};
		
		genCalc->ComputePricingData();
		argResult = genCalc->GetPricingData().GetData(key);

		return ARM_OK;
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_CallOnMepiVanillaArgCreate(
		const string& CurveName,
		const string& EquityName,
		double startDate,
		double endDate,
		long resetFreq,
		double riskFactor,
		double strike,
		double maxBorrow,
		double protectionCurveStart,
		double protectionCurveEnd,
		double startingPortfolio,
		double startingCash,
		double minInvested,
		double leverageCost,
		double cashSpread,
		double fees,
		double alreadyAsianed,
		long asianDatesNb,
        ARM_result&  result,
        long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_GenSecurity * gensec = NULL;

	try
	{
		ARM_VanillaMepi vm(CurveName, EquityName, ConvertXLDateToJulian(startDate), ConvertXLDateToJulian(endDate), resetFreq, riskFactor, strike, maxBorrow, protectionCurveStart, protectionCurveEnd, startingPortfolio, startingCash, minInvested, leverageCost, cashSpread, fees, alreadyAsianed, asianDatesNb );
		vm.CreateAndSetGenSec(CurveName, ConvertXLDateToJulian(startDate));

		gensec = static_cast<ARM_GenSecurity *> ( vm.GetGenSecurity()->Clone() );

		/// assign object
		if( !assignObject( gensec, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }

	}
	catch(Exception& x)
	{
		delete gensec;
		x.DebugPrint();
		ARM_RESULT();
	}
}

int ARMLOCAL_IsWarning_OnObject( long C_ObjIdId,
        ARM_result&  result )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
   
	try
	{	
		ARM_Object* obj = LOCAL_PERSISTENT_OBJECTS->GetObject(C_ObjIdId);
		ARM_GenCalculator* genCalc		= dynamic_cast<ARM_GenCalculator*>( obj );

		if( genCalc )
			return genCalc->GetWarning() != NULL;
		else 
		{
			ARM_CalibMethod* calibMethod	= dynamic_cast<ARM_CalibMethod*>( obj );
			if( calibMethod )
				return calibMethod->GetWarning() != NULL;
			else return false;
		}
	}
	
	catch(Exception& x)
	{
		x.DebugPrint();
		ARM_RESULT();
	}
}


long ARMLOCAL_GetWarning_OnObject( const long& C_ObjIdId,
	ARM_result&  result,
	long objId
	)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_Warning* warning= NULL;
   
	try
	{	
		ARM_Object* obj = LOCAL_PERSISTENT_OBJECTS->GetObject(C_ObjIdId);
		ARM_GenCalculator* genCalc		= dynamic_cast<ARM_GenCalculator*>( obj );

		if( genCalc )
			warning = genCalc->GetWarning();
		else
		{
			ARM_CalibMethod* calibMethod	= dynamic_cast<ARM_CalibMethod*>( obj );
			if( calibMethod )
				warning = calibMethod->GetWarning();
			else
				ARM_THROW( ERR_INVALID_ARGUMENT, "object says there is warning but found no warning!" );
		}
		
		if( !warning )
			ARM_THROW( ERR_INVALID_ARGUMENT, "object says there is warning but found no warning!" );

		// assign object
		if( !assignObject( warning, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		delete warning;
		x.DebugPrint();
		ARM_RESULT();
	}
}

/****************************************************************
	General Swaption Converter to Variable Notional Swaption
*****************************************************************/
long ARMLOCAL_ConvertToVarNotionalSwaption(long yieldCurveId,long swaptionId,
										   ARM_result& result,long objId)		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_ZeroCurve* yieldCurve = NULL;
	ARM_Swaption *swaption = NULL;
	ARM_Swaption *vnSwaption = NULL;
	try
	{
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Swaption Converter" );

		/// Get yield curve
		yieldCurve = dynamic_cast<ARM_ZeroCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(yieldCurveId));
		if(!yieldCurve)
		{
			result.setMsg ("ARM_ERR: Yield Curve is not of a good type");
			return ARM_KO;
		}

		/// Get general ARM swaption and convert it to vanilla arg
		swaption = dynamic_cast<ARM_Swaption *>(LOCAL_PERSISTENT_OBJECTS->GetObject(swaptionId));
		if(!swaption)
		{
			result.setMsg ("ARM_ERR: Swaption is not of a good type");
			return ARM_KO;
		}

		double asOfDate = yieldCurve->GetAsOfDate().GetJulian();
		string modelName(ARM::ARM_VanillaPricer::GetDefaultModelName(swaption,NULL));
		const ARM_VanillaSwaptionArg& swaptionVanillaArg = static_cast<const ARM_VanillaSwaptionArg&>(* ARM::ARM_ConverterFromKernel::ConvertSecuritytoArgObject(swaption,asOfDate,modelName) );

		/// Convert to VNS
        vnSwaption = ARM_CRASpreadCalculator::ConvertToVarNotionalSwaption(yieldCurve,swaptionVanillaArg);
		 						
		/// Assign new bject
		if ( !assignObject( vnSwaption, result, objId ) )
			return ARM_KO; 
		else
			return ARM_OK; 
	}
	
	catch(Exception& x)
	{	
		x.DebugPrint();
		ARM_RESULT();
	}
}



/****************************************************************
	Basis converter
*****************************************************************/

long ARMLOCAL_BasisConverter(
	const double& asOfDate,
	const string& domCcyStr,
	const string& forCcyStr,
	const long& domDateStripId,
	const long& forDateStripId,
	const long& fundDateStripId,
	const string& domDayCount,
	const string& domFreq,
	const string& forDayCount,
	const string& forFreq,
	const long& domZcId,
	const long& forZcId,
	const long& domDiscZcId,
	const long& forDiscZcId,
	const long& forexId,
	const long& domNotionalId,
	const long& forNotionalId,
	const long& forSpreadId,
	VECTOR<double>& Margin,
	ARM_result& result)
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM_DateStrip* domDateStrip;
	ARM_DateStrip* forDateStrip;
	ARM_DateStrip* fundDateStrip;
	int domDayCountInt;
	int domFreqInt;
	int forDayCountInt;
	int forFreqInt;
	ARM_ZeroCurve* domZc;
	ARM_ZeroCurve* forZc;
	ARM_ZeroCurve* domDiscZc;
	ARM_ZeroCurve* forDiscZc;
	ARM_Forex* forex;
	ARM_Curve* domNotional;
	ARM_Curve* forNotional;
	ARM_Curve* forSpread;

	try
	{
		/// CRM tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Basis Converter" );

		ARM_Currency domCcy(domCcyStr.c_str());
		ARM_Currency forCcy;(forCcyStr.c_str());

		// Convert into Julian dates
		double asOfDateJulian = ConvertXLDateToJulian(asOfDate);

		/// Get Domestic Date Strip
		domDateStrip = dynamic_cast<ARM_DateStrip *>(LOCAL_PERSISTENT_OBJECTS->GetObject(domDateStripId));
		if(!domDateStrip)
		{
			result.setMsg ("ARM_ERR: Domestic date strip is not of a good type");
			return ARM_KO;
		}

		/// Get Foreign Date Strip
		forDateStrip = dynamic_cast<ARM_DateStrip *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forDateStripId));
		if(!forDateStrip)
		{
			result.setMsg ("ARM_ERR: Foreign date strip is not of a good type");
			return ARM_KO;
		}

		/// Get Funding Date Strip
		fundDateStrip = dynamic_cast<ARM_DateStrip *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundDateStripId));
		if(!fundDateStrip)
		{
			result.setMsg ("ARM_ERR: Funding date strip is not of a good type");
			return ARM_KO;
		}

		// Get the domestic day count and frequency
		domDayCountInt = ARM_ArgConv_LgNameDayCount.GetNumber(domDayCount);
		domFreqInt = ARM_ArgConv_LgNameFrequency.GetNumber(domFreq);

		// Get the foreign day count and frequency
		forDayCountInt = ARM_ArgConv_LgNameDayCount.GetNumber(forDayCount);
		forFreqInt = ARM_ArgConv_LgNameFrequency.GetNumber(forFreq);

		// Get the domestic curve
		domZc = dynamic_cast<ARM_ZeroCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(domZcId));
		if(!domZc)
		{
			result.setMsg ("ARM_ERR: Domestic Curve is not of a good type");
			return ARM_KO;
		}

		// Get the foreign curve
		forZc = dynamic_cast<ARM_ZeroCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forZcId));
		if(!forZc)
		{
			result.setMsg ("ARM_ERR: Foreign Curve is not of a good type");
			return ARM_KO;
		}

		// Get the domestic discount curve (basis)
		domDiscZc = dynamic_cast<ARM_ZeroCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(domDiscZcId));
		if(!domDiscZc)
		{
			result.setMsg ("ARM_ERR: Domestic Discount Curve is not of a good type");
			return ARM_KO;
		}

		// Get the foreign discount curve (basis)
		forDiscZc = dynamic_cast<ARM_ZeroCurve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forDiscZcId));
		if(!forDiscZc)
		{
			result.setMsg ("ARM_ERR: Foreign Discount Curve is not of a good type");
			return ARM_KO;
		}

		// Get the forex
		forex = dynamic_cast<ARM_Forex *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forexId));
		if(!forex)
		{
			result.setMsg ("ARM_ERR: Forex is not of a good type");
			return ARM_KO;
		}

		// Get the domestic notional
		domNotional = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(domNotionalId));
		if(!domNotional)
		{
			result.setMsg ("ARM_ERR: Domestic notional is not of a good type");
			return ARM_KO;
		}

		// Get the foreign notional
		forNotional = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forNotionalId));
		if(!forNotional)
		{
			result.setMsg ("ARM_ERR: Foreign notional is not of a good type");
			return ARM_KO;
		}

		// Get the foreign spread
		forSpread = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(forSpreadId));
		if(!forSpread)
		{
			result.setMsg ("ARM_ERR: Foreign spread is not of a good type");
			return ARM_KO;
		}

		size_t i;
		size_t domSize = domDateStrip->size();
		size_t forSize = forDateStrip->size();

		ARM_GP_Vector* domStartDates = domDateStrip->GetFlowStartDates();
		ARM_GP_Vector* forStartDates = forDateStrip->GetFlowStartDates();
		ARM_GP_Vector* fundStartDates = fundDateStrip->GetFlowStartDates();

		ARM_GP_Vector domNotionalVec(domSize);
		ARM_GP_Vector forNotionalVec(forSize), forSpreadVec(forSize);

		for (i = 0; i < domSize; ++i)
		{
			domNotionalVec[i] = domNotional->Interpolate((*domStartDates)[i]-asOfDateJulian);
		}

		for (i = 0; i < forSize; ++i)
		{
			forNotionalVec[i] = forNotional->Interpolate((*forStartDates)[i]-asOfDateJulian);
			forSpreadVec[i] = forSpread->Interpolate((*forStartDates)[i]-asOfDateJulian);
		}

		ARM_DateStripPtr domDateStripPtr(CreateClonedPtr(domDateStrip));
		ARM_DateStripPtr forDateStripPtr(CreateClonedPtr(forDateStrip));
		ARM_DateStripPtr fundDateStripPtr(CreateClonedPtr(fundDateStrip));

		ARM_BasisConverter basisConverter( 
		domCcy,
		forCcy,
		domDateStripPtr,
		forDateStripPtr,
		fundDateStripPtr,
		domZc,
		forZc,
		domDiscZc,
		forDiscZc,
		*forex, 
		domNotionalVec,
		forNotionalVec,
		forSpreadVec);

		ARM_GP_Vector margin = basisConverter.ComputeDomMargin();
		
		ARM_Curve marginCurve(*domStartDates, margin, new ARM_StepUpLeftOpenCstExtrapolDble);

		Margin.resize(fundStartDates->size());

		for (i = 0; i < Margin.size(); ++i)
		{
			Margin[i] = marginCurve.Interpolate((*fundStartDates)[i]);
		}
	}
	catch(Exception& x)
	{	
		x.DebugPrint();
		ARM_RESULT();
	}

	return ARM_OK;
}
