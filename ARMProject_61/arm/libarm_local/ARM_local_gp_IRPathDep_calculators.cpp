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

#include "ARM_local_gp_IRPathDep_calculators.h"
#include "ARM_local_gp_calculators_tools.h"
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

////////////////////////////////////////////
//// Function to create a date strip combiner
////////////////////////////////////////////
////////////////////////////////////////////
//// Function to create a TARN Calculator
////////////////////////////////////////////
bool TARNCalculatorCheckAutoCalFlags(const vector< string >& flags,size_t nbFlagsToCheck)
{
    size_t nb=(flags.size()>nbFlagsToCheck ? nbFlagsToCheck : flags.size());
    for(size_t i=0;i<nb-1;++i)
        if(flags[i] != "Y" && flags[i] != "YES" && flags[i] != "N" && flags[i] != "NO")
            return false;

	if(flags[i] != "CF" && flags[i] != "MC" && flags[i] != "N" && flags[i] != "NO")
            return false;

    return true;
}

////////////////////////////////////////////
//// Function to create a TARN Calculator
////////////////////////////////////////////

long ARMLOCAL_TARNCalculator_Create(double startDate,
                                    double endDate,
                                    double strike,
                                    long strikeId,
                                    long payRec,
                                    long cpnDayCount,
                                    long cpnFreq,
                                    long cpnTiming,
                                    const string& cpnIndexTerm,
                                    long cpnIndexDayCount,
                                    const string& cpnResetCal,
                                    const string& cpnPayCal,
                                    long cpnResetGap,
		                            long intRule,
                                    double leverage,
                                    long leverageId,
									double cpnMin,
									long cpnMinId,
									double cpnMax,
									long cpnMaxId,
		                            double lifeTimeCapTarget,
									bool globalCapFlag,
		                            double lifeTimeFloorTarget,
                                    double fundSpread,
                                    long fundSpreadId,
                                    long fundFreq,
                                    long fundDayCount,
                                    double nominal,
                                    long nominalId,
		                            double fees,
		                            long feesId,
                                    long fundNomId,
		                            const vector< double >& nbIterations,
					                const string& genType1,
		                            const string& genType2,
		                            int firstNbTimes,
		                            int firstNbDims,
		                            const string& pathScheme,
		                            const string& pathOrder,
		                            const string& antithetic,
                                    const vector< string >& calibFlags,
                                    const vector< string >& outputFlags,
                                    long mktDataManagerId,
                                    const vector< string >& keys,
                                    bool isCustomResetFlag,
		                            long resetDatesId,
                                    double AsOfDate,
		                            ARM_result&	result, 
                                    long objId)
{
	/// input checks
	if ( !GlobalPersistanceOk(result) )
	   return(ARM_KO);

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

    ARM_TARNCalculator* tarnCalculator = NULL;

	ARM_Curve* strikeProfile     = NULL;
    bool isStrike     = false;
	ARM_Curve* leverageProfile   = NULL;
    bool isLeverage   = false;
	ARM_Curve* cpnMinProfile   = NULL;
    bool isCpnMin   = false;
	ARM_Curve* cpnMaxProfile   = NULL;
    bool isCpnMax   = false;
	ARM_Curve* fundSpreadProfile = NULL;
    bool isFundSpread = false;
	ARM_Curve* nominalProfile    = NULL;
    bool isNominal    = false;
	ARM_Curve* feesProfile       = NULL;
    bool isFees		  = false;
	ARM_Curve* fundNomProfile    = NULL;
    bool isFundNom    = false;

	ARM_GP_Vector* customResetDates = NULL;
	ARM_GP_Vector* julianResetDates = new ARM_GP_Vector(0);
	

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN Calculator" );

		char myStartDate[20];
		char myEndDate[20];
        char myAsOfDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

        if ( AsOfDate > 0 )
        {
           Local_XLDATE2ARMDATE(AsOfDate, myAsOfDate);
        }

        /// Convert refValue Id to object if possible
		
        if ( strikeId != ARM_NULL_OBJECT )
        {
		    strikeProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,strike);
           
            strikeProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isStrike = true;
        }

        if ( leverageId != ARM_NULL_OBJECT )
        {
		    leverageProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId));

		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: leverage is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,leverage);
            
            leverageProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isLeverage = true;
        }

        if(cpnMinId != ARM_NULL_OBJECT)
        {
		    cpnMinProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMinId));

		    if (!cpnMinProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon min is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,cpnMin);
            cpnMinProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCpnMin=true;
        }

        if(cpnMaxId != ARM_NULL_OBJECT)
        {
		    cpnMaxProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMaxId));

		    if (!cpnMaxProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon max is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,cpnMax);
            cpnMaxProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCpnMax=true;
        }

        if ( fundSpreadId != ARM_NULL_OBJECT )
        {
		    fundSpreadProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			   result.setMsg("ARM_ERR: funding spread is not of a good type");
			   
               return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundSpread);
            
            fundSpreadProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFundSpread      = true;
        }

        if ( nominalId != ARM_NULL_OBJECT )
        {
		    nominalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(nominalId));

		    if (!nominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,nominal);
            
            nominalProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isNominal      = true;
        }

		if ( feesId != ARM_NULL_OBJECT )
        {
		    feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));

		    if (!feesProfile)
		    {
			    result.setMsg ("ARM_ERR: fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fees);

            feesProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFees      = true;
        }

        if ( fundNomId != ARM_NULL_OBJECT )
        {
		    fundNomProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundNomId));

		    if (!fundNomProfile)
		    {
			    result.setMsg ("ARM_ERR: funding nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundNomProfile = new ARM_Curve(*nominalProfile);
			isFundNom = true;
        }

        /// Restore Market Data Manager Rep

        ARM_MarketData_ManagerRep* mktDataManager = NULL;

        if ( mktDataManagerId != ARM_NULL_OBJECT )
        {
		   mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		
           if (!mktDataManager)
           {
			  result.setMsg("ARM_ERR: market data manager is not of a good type");
			  
              return(ARM_KO);
           }
        }

 		/// Get cpn, funding, basis, domestic and foreign CCYs

        ARM_Currency DefCCy("EUR");

 		ARM_Currency* cpnCcy = NULL;
        
        if ( mktDataManager != NULL )
        {
           cpnCcy = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::YcKey]))->GetCurrencyUnit();
        }
        else
        {
           cpnCcy = &DefCCy;
        }

 		ARM_Currency* fundCcy  = NULL;
 		ARM_Currency* basisCcy = NULL;
 		ARM_Currency* forCcy   = NULL;
 		ARM_Currency* domCcy   = NULL;
 		ARM_Forex*    forex    = NULL;

 		if (( mktDataManager != NULL )
            &&
            ( keys.size() != 0 )
            &&
            ( keys.size() == ARM_TARNCalculator::NbKeys )
           )
 		{
 			fundCcy  = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::FundingKey]))->GetCurrencyUnit();
 			basisCcy = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::BasisKey]))->GetCurrencyUnit();
 
 			forex    = static_cast< ARM_Forex* >(mktDataManager->GetData(keys[ARM_TARNCalculator::ForexKey]));
 			forCcy   = forex->GetMainCurrency();
 			domCcy   = forex->GetMoneyCurrency();
 		}
 		else
 		{
 			fundCcy  = cpnCcy;
 			basisCcy = cpnCcy;
 			domCcy   = cpnCcy;
 			forCcy   = cpnCcy;
 		}

		if ((calibFlags.size() < 5) && (calibFlags.size() > 8))
		{
			result.setMsg ("ARM_ERR: calibration flags size has to be greater 5 than and less than 8");
			return ARM_KO;
		}

		ARM_TARNCalculator::CalibrationMode capFloorCalibMode = (ARM_TARNCalculator::CalibrationMode) ARM::ARM_ArgConv_TARNCalibMode.GetNumber(calibFlags[0]); 
		ARM_TARNCalculator::CalibrationMode digitalCalibMode = (ARM_TARNCalculator::CalibrationMode) ARM::ARM_ArgConv_TARNCalibMode.GetNumber(calibFlags[1]);
	  
        bool oswCalibFlag=(calibFlags[2]=="Y" || calibFlags[2]=="YES");
		bool controlVariableFlag=(calibFlags[3]=="Y" || calibFlags[3]=="YES");
		bool digitalSmoothingFlag=(calibFlags[4]=="Y" || calibFlags[4]=="YES");

		ARM_ModelType modelType;
		if(calibFlags.size() >= 6)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibFlags[5]);
		else
			modelType = ARM_PricingModelType::SFRM2F;

		bool smiledFRMRescallingFlag;
		if(calibFlags.size() >= 7)
			smiledFRMRescallingFlag = (calibFlags[6]=="Y" || calibFlags[6]=="YES");
		else
			smiledFRMRescallingFlag = false;

		bool antitheticFlag;
		if (antithetic == "Y")
			antitheticFlag = true;
		else if (antithetic == "N")
			antitheticFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: antithetic flag should be Y or N.");
			return ARM_KO;
		}

		if ( outputFlags.size() > 11 )
		{
			result.setMsg ("ARM_ERR: output flags size has to be less than 10.");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(11,false);

		productsToPrice[0] = true;

		for (size_t i = 0; i < outputFlags.size(); ++i)
		{
			productsToPrice[i] = ((outputFlags[i] == "Y") || (outputFlags[i] == "YES"));
		}

		//Do we have custom reset dates
		if (isCustomResetFlag == true)
		{
			double xlToJulianDate = 2415019.0;
			customResetDates = dynamic_cast<ARM_GP_Vector *>(LOCAL_PERSISTENT_OBJECTS->GetObject(resetDatesId));
			julianResetDates->resize(customResetDates->size());
			for (int z=0; z < customResetDates->size(); z++)
				(*julianResetDates)[z] = (*customResetDates)[z] + xlToJulianDate;
		}
		
		/// Create the TARN calculator
        tarnCalculator = new ARM_TARNCalculator(mktDataManager ? mktDataManager->GetAsOfDate(): 
                                                                 ARM_Date(myAsOfDate),
			                                    ARM_Date(myStartDate),
			                                    ARM_Date(myEndDate),
			                                    *strikeProfile,
			                                    payRec,
			                                    cpnDayCount,
			                                    cpnFreq,
			                                    cpnTiming,
			                                    cpnIndexTerm,
			                                    cpnIndexDayCount,
			                                    cpnResetCal,
			                                    cpnPayCal,
			                                    cpnResetGap,
			                                    intRule,
			                                    *leverageProfile,
												*cpnMinProfile,
												*cpnMaxProfile,
			                                    lifeTimeCapTarget,
												globalCapFlag,
			                                    lifeTimeFloorTarget,
			                                    *fundSpreadProfile,
			                                    fundFreq,
			                                    fundDayCount,
			                                    *nominalProfile,
			                                    *feesProfile,
			                                    ARM_GP_Vector( nbIterations ),
			                                    modelType,
			                                    capFloorCalibMode,
			                                    digitalCalibMode,
			                                    oswCalibFlag,
			                                    controlVariableFlag,
			                                    digitalSmoothingFlag,
												smiledFRMRescallingFlag,
			                                    genType1,
			                                    genType2,
			                                    firstNbTimes,
			                                    firstNbDims,
			                                    pathScheme,
			                                    pathOrder,
			                                    antitheticFlag,
			                                    productsToPrice,
 			                                    *cpnCcy,
 			                                    *fundCcy,
 			                                    *basisCcy,
 			                                    *domCcy,
 			                                    *forCcy,
 			                                    *fundNomProfile,
 			                                    &keys,
 			                                    mktDataManager,
			                                    isCustomResetFlag,
			                                    *julianResetDates);
			
		

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
		if(isFees)
            delete feesProfile;
		feesProfile = NULL;
		if (isFundNom)
			delete fundNomProfile;
		fundNomProfile = NULL;
		
		delete julianResetDates;
		julianResetDates = NULL;

		/// assign object
		if( !assignObject( tarnCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isStrike)
            delete strikeProfile;
        if(isLeverage)
            delete leverageProfile;
		if(isCpnMin)
			delete cpnMinProfile;
		cpnMinProfile = NULL;
		if(isCpnMax)
			delete cpnMaxProfile;
		cpnMaxProfile = NULL;
        if(isFundSpread)
            delete fundSpreadProfile;
        if(isNominal)
            delete nominalProfile;
		if(isFees)
            delete feesProfile;
		feesProfile = NULL;
		if (isFundNom)
			delete fundNomProfile;

		if (isCustomResetFlag)
		{
			//delete customResetDates;
			//customResetDates = NULL;

			delete julianResetDates;
			julianResetDates = NULL;
		}

		delete tarnCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}


////////////////////////////////////////////
//// Function to create a TARN SnowBall Calculator
////////////////////////////////////////////

long ARMLOCAL_TARNSBCalculator_Create(
        double startDate,
        double endDate,
        double strike,
        long strikeId,
		double coupon0,
        long payRec,
        long cpnDayCount,
        long cpnFreq,
        long cpnTiming,
        const string& cpnIndexTerm,
        long cpnIndexDayCount,
        const string& cpnResetCal,
        const string& cpnPayCal,
        long cpnResetGap,
		long intRule,
		string intSubRule,
        double leverage,
        long leverageId,
		double cpnMin,
		long cpnMinId,
		double cpnMax,
		long cpnMaxId,
		double levPrev,
        long levPrevId,
		long reverse,
		double lifeTimeCapTarget,
		bool globalCapFlag,
		double lifeTimeFloorTarget,
        double fundSpread,
        long fundSpreadId,
        long fundFreq,
        long fundDayCount,
        double nominal,
        long nominalId,
		double fees,
		long feesId,
		double fundNominal,
		long fundNomId,
		const vector< double >& nbIterations,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
        const vector< string >& calibFlags,
        const vector< string >& outputFlags,
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
    ARM_TARNCalculatorSnowBall* tarnCalculator = NULL;

	ARM_Curve* strikeProfile=NULL;
    bool isStrike=false;
	ARM_Curve* leverageProfile=NULL;
    bool isLeverage=false;
	ARM_Curve* cpnMinProfile   = NULL;
    bool isCpnMin   = false;
	ARM_Curve* cpnMaxProfile   = NULL;
    bool isCpnMax   = false;
	ARM_Curve* levPrevProfile = NULL;
    bool isLevPrev=false;
	ARM_Curve* fundSpreadProfile=NULL;
    bool isFundSpread=false;
	ARM_Curve* nominalProfile=NULL;
    bool isNominal=false;
	ARM_Curve* feesProfile=NULL;
    bool isFees=false;
	ARM_Curve* fundNomProfile=NULL;
	bool isFundNom=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "TARN Snow Ball Calculator" );

		char myStartDate[20];
		char myEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);

        /// Convert refValue Id to object if possible
		
        if(strikeId != ARM_NULL_OBJECT)
        {
		    strikeProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId));

		    if (!strikeProfile)
		    {
			    result.setMsg ("ARM_ERR: strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,strike);
            strikeProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isStrike=true;
        }

        if(leverageId != ARM_NULL_OBJECT)
        {
		    leverageProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId));

		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: leverage is not of a good type");
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

		if(cpnMinId != ARM_NULL_OBJECT)
        {
		    cpnMinProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMinId));

		    if (!cpnMinProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon min is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,cpnMin);
            cpnMinProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCpnMin=true;
        }

        if(cpnMaxId != ARM_NULL_OBJECT)
        {
		    cpnMaxProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cpnMaxId));

		    if (!cpnMaxProfile)
		    {
			    result.setMsg ("ARM_ERR: coupon max is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,cpnMax);
            cpnMaxProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isCpnMax=true;
        }

		if(levPrevId != ARM_NULL_OBJECT)
        {
		    levPrevProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(levPrevId));

		    if (!levPrevProfile)
		    {
			    result.setMsg ("ARM_ERR: levprev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,levPrev);
            levPrevProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isLevPrev=true;
        }

        if(fundSpreadId != ARM_NULL_OBJECT)
        {
		    fundSpreadProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundSpreadId));

		    if (!fundSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: funding spread is not of a good type");
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

        if(nominalId != ARM_NULL_OBJECT)
        {
		    nominalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(nominalId));

		    if (!nominalProfile)
		    {
			    result.setMsg ("ARM_ERR: nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,nominal);
            nominalProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isNominal=true;
        }

		if(feesId != ARM_NULL_OBJECT)
        {
		    feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));

		    if (!feesProfile)
		    {
			    result.setMsg ("ARM_ERR: fees is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fees);
            feesProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFees=true;
        }

		if ( fundNomId != ARM_NULL_OBJECT )
        {
		    fundNomProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundNomId));

		    if (!fundNomProfile)
		    {
			    result.setMsg ("ARM_ERR: funding nominal is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundNominal);
            fundNomProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
			isFundNom = true;
        }

        /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		/// Get cpn, funding, basis, domestic and foreign CCYs

        ARM_Currency DefCCy("EUR");

 		ARM_Currency* cpnCcy = NULL;
        
        if ( mktDataManager != NULL )
        {
           cpnCcy = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::YcKey]))->GetCurrencyUnit();
        }
        else
        {
           cpnCcy = &DefCCy;
        }

 		ARM_Currency* fundCcy  = NULL;
 		ARM_Currency* basisCcy = NULL;
 		ARM_Currency* forCcy   = NULL;
 		ARM_Currency* domCcy   = NULL;
 		ARM_Forex*    forex    = NULL;

 		if (( mktDataManager != NULL )
            &&
            ( keys.size() != 0 )
            &&
            ( keys.size() == ARM_TARNCalculator::NbKeys )
           )
 		{
 			fundCcy  = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::FundingKey]))->GetCurrencyUnit();
 			basisCcy = static_cast< ARM_ZeroCurve* >(mktDataManager->GetData(keys[ARM_TARNCalculator::BasisKey]))->GetCurrencyUnit();
 
 			forex    = static_cast< ARM_Forex* >(mktDataManager->GetData(keys[ARM_TARNCalculator::ForexKey]));
 			forCcy   = forex->GetMainCurrency();
 			domCcy   = forex->GetMoneyCurrency();
 		}
 		else
 		{
 			fundCcy  = cpnCcy;
 			basisCcy = cpnCcy;
 			domCcy   = cpnCcy;
 			forCcy   = cpnCcy;
 		}

		if ((calibFlags.size() < 5) && (calibFlags.size() > 8))
		{
			result.setMsg ("ARM_ERR: calibration flags size has to be greater 5 than and less than 8");
			return ARM_KO;
		}

		ARM_TARNCalculator::CalibrationMode capFloorCalibMode = (ARM_TARNCalculator::CalibrationMode) ARM::ARM_ArgConv_TARNCalibMode.GetNumber(calibFlags[0]); 
		ARM_TARNCalculator::CalibrationMode digitalCalibMode = (ARM_TARNCalculator::CalibrationMode) ARM::ARM_ArgConv_TARNCalibMode.GetNumber(calibFlags[1]);
	    bool oswCalibFlag=(calibFlags[2]=="Y" || calibFlags[2]=="YES");
		bool controlVariableFlag=(calibFlags[3]=="Y" || calibFlags[3]=="YES");
		bool digitalSmoothingFlag=(calibFlags[4]=="Y" || calibFlags[4]=="YES");

		bool antitheticFlag;
		if (antithetic == "Y")
			antitheticFlag = true;
		else if (antithetic == "N")
			antitheticFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: antithetic flag should be Y or N.");
			return ARM_KO;
		}

		ARM_ModelType modelType;
		if(calibFlags.size() >= 6)
			modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(calibFlags[5]);
		else
			modelType = ARM_PricingModelType::SFRM2F;

		bool smiledFRMRescallingFlag;
		if(calibFlags.size() >= 7)
			smiledFRMRescallingFlag = (calibFlags[6]=="Y" || calibFlags[6]=="YES");
		else
			smiledFRMRescallingFlag = false;
		
		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(11,false);

		productsToPrice[0] = true;

		for (size_t i = 0; i < outputFlags.size(); ++i)
		{
			productsToPrice[i] = ((outputFlags[i] == "Y") || (outputFlags[i] == "YES"));
		}

		int stubRule  = ARM_ArgConv_StubRules.GetNumber(intSubRule);
		/// Create the TARN calculator
        tarnCalculator = new ARM_TARNCalculatorSnowBall(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
            *strikeProfile,
			coupon0,
            payRec,
            cpnDayCount,
            cpnFreq,
            cpnTiming,
            cpnIndexTerm,
            cpnIndexDayCount,
            cpnResetCal,
            cpnPayCal,
            cpnResetGap,
			intRule,
			stubRule,
			reverse,
            *leverageProfile,
			*cpnMinProfile,
			*cpnMaxProfile,
			*levPrevProfile,
			lifeTimeCapTarget,
			globalCapFlag,
			lifeTimeFloorTarget,
            *fundSpreadProfile,
            fundFreq,
            fundDayCount,
            *nominalProfile,
			*feesProfile,
			*fundNomProfile,
			ARM_GP_Vector( nbIterations ),
			modelType,
			capFloorCalibMode,
			digitalCalibMode,
			oswCalibFlag,
			controlVariableFlag,
			digitalSmoothingFlag,
			smiledFRMRescallingFlag,
			genType1,
			genType2,
			firstNbTimes,
			firstNbDims,
			pathScheme,
			pathOrder,
			antitheticFlag,
			productsToPrice,
			*cpnCcy,
 			*fundCcy,
 			*basisCcy,
 			*domCcy,
 			*forCcy,
            *mktDataManager,
            keys);

        /// Free memory
        if(isStrike)
            delete strikeProfile;
		isStrike = NULL;
        if(isLeverage)
            delete leverageProfile;
		if(isCpnMin)
			delete cpnMinProfile;
		cpnMinProfile = NULL;
		if(isCpnMax)
			delete cpnMaxProfile;
		cpnMaxProfile = NULL;
		leverageProfile = NULL;
		if(isLevPrev)
            delete levPrevProfile;
		levPrevProfile = NULL;
        if(isFundSpread)
            delete fundSpreadProfile;
		fundSpreadProfile = NULL;
        if(isNominal)
            delete nominalProfile;
		nominalProfile = NULL;
		if(isFees)
            delete feesProfile;
		feesProfile = NULL;
		if (isFundNom)
			delete fundNomProfile;
		fundNomProfile = NULL;
		
		/// assign object
		if( !assignObject( tarnCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isStrike)
            delete strikeProfile;
        if(isLeverage)
            delete leverageProfile;
		if(isCpnMin)
			delete cpnMinProfile;
		cpnMinProfile = NULL;
		if(isCpnMax)
			delete cpnMaxProfile;
		cpnMaxProfile = NULL;
		if(isLevPrev)
            delete levPrevProfile;
		levPrevProfile = NULL;
        if(isFundSpread)
            delete fundSpreadProfile;
        if(isNominal)
            delete nominalProfile;
		if(isFees)
            delete feesProfile;
		feesProfile = NULL;

		delete tarnCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to get data from a TARN Calculator
///////////////////////////////////////////////
long ARMLOCAL_TARN_Get(
        long tarnId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_TARNCalculator* tarnCalculator;
    ARM_Object* object=NULL;

	try
	{
		tarnCalculator = dynamic_cast<ARM_TARNCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(tarnId));

		if (!tarnCalculator)
		{
			result.setMsg ("ARM_ERR: TARN Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_OSW_PORT ||
            typeToGet == GC_ACCESS_CF_PORT)
        {
            if(typeToGet == GC_ACCESS_OSW_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(tarnCalculator->GetOSWPortfolio())).Clone();
            else if (typeToGet == GC_ACCESS_CF_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(tarnCalculator->GetCFPortfolio())).Clone();

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
            return ARMLOCAL_GC_Get(tarnId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a TARN Calculator
///////////////////////////////////////////////
long ARMLOCAL_TARN_Set(
        long tarnId,
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
    ARM_TARNCalculator* oldTARN=NULL;
    ARM_TARNCalculator* newTARN=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");


	try
	{
		oldTARN = dynamic_cast<ARM_TARNCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(tarnId));
		if (!oldTARN)
		{
			result.setMsg ("ARM_ERR: TARN Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
            if(isUpdated)
                /// Just update the initial TARN
                newTARN=oldTARN;
            else
                /// Clone the initial TARN then set
                newTARN = dynamic_cast<ARM_TARNCalculator *>(oldTARN->Clone());

            string typeToSet(stringGetUpper(setType));
            if(typeToSet == GC_ACCESS_CF_PORT)
            {
                newTARN->SetCFPortfolio(*port);
                updateMsg += "C/F portfolio";
            }
            else if(typeToSet == GC_ACCESS_OSW_PORT)
            {
                /// Default
                newTARN->SetOSWPortfolio(*port);
                updateMsg += "OSW portfolio";
            }

            /// Keep calibration data consistency
            newTARN->UpdateCalibrationAndTimeIt();

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new TARN in ARM cache
		    if(!assignObject( newTARN, result, objId ) )
            {
                delete newTARN;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
        else
            return ARMLOCAL_GC_Set(tarnId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newTARN;
		x.DebugPrint();
		ARM_RESULT();
	}
}

//Aim: Get a Tarn security and set the calculated outputs
long ARMLOCAL_TARN_SetProductToPrice(
	long tarnId,
	vector<string>& prodToPrice,
	ARM_result&	result,
	long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

    ARM_TARNCalculator* tarn=NULL;

	CCString msg(" ");
	CCString thisClass;

	try
	{
		//We get the object to modify
		tarn = dynamic_cast<ARM_TARNCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(tarnId));
		if (!tarn)
		{
			result.setMsg ("ARM_ERR: Check your TARN ticker");
			return ARM_KO;
		}
		
		ARM_TARNCalculator* newTarn = NULL; 
		//We create a new object.
		if (objId == ARM_NULL_OBJECT_ID)
		{	
			newTarn = (ARM_TARNCalculator*) tarn->Clone();
			newTarn->SetProductToPrice(prodToPrice);
			newTarn->SetHasToPrice(false);
			CREATE_GLOBAL_OBJECT();
			objId = LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newTarn);
			if (objId == RET_KO)
			{
				if (newTarn)
					delete newTarn;
				newTarn = NULL;

				result.setMsg ("ARM_ERR: Pb with inserting object");				
				return ARM_KO;
			}

			result.setLong(objId);

			return ARM_OK;
		}
		//We modify an existing object.
		else
		{
			newTarn = dynamic_cast<ARM_TARNCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(objId));
			delete newTarn;
			newTarn = (ARM_TARNCalculator*) tarn->Clone();
			newTarn->SetProductToPrice(prodToPrice);
			newTarn->SetHasToPrice(false);
			
			LOCAL_PERSISTENT_OBJECTS->SetPersistent((ARM_Object*)newTarn, objId);

			return ARM_OK;
		}
	}

	catch(Exception& x)
	{
		x.DebugPrint();

		ARM_RESULT();
	}
}

////////////////////////////////////////////
//// Function to create a Maturity Cap Calculator
////////////////////////////////////////////

long ARMLOCAL_MaturityCapCalculator_Create(
        double startDate,
        double endDate,
		double underlyingEndDate,
        long longShort,
        long capFloor,
        long resetFreq,
		long payFreq,
        const string& indexTerm,
        long dayCount,
		long intRule,
		double spread,
		double initNominal,
		double initTRI,
		double annuity,
		long maturityCapMode,
		double coeff,
		double amortizing,
		long amortizingId,
		long resetGap,
        const string& resetCal,
        const string& payCal,
		long calibrationMode,
		long nbIterations,
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
    ARM_MaturityCapCalculator* maturityCapCalculator = NULL;

	ARM_Curve* amortizingProfile=NULL;
    bool isAmortizing=false;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Maturity CAP Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		char myUnderlyingEndDate[20];

        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);
		Local_XLDATE2ARMDATE(underlyingEndDate, myUnderlyingEndDate);

        /// Convert refValue Id to object if possible
		
        if(amortizingId != ARM_NULL_OBJECT)
        {
		    amortizingProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(amortizingId));

		    if (!amortizingProfile)
		    {
			    result.setMsg ("ARM_ERR: amortizing is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,amortizing);
            amortizingProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isAmortizing=true;
        }

        /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(4,false);

		productsToPrice[0] = true;

		for (size_t i = 0; i < flags.size(); ++i)
		{
			productsToPrice[i] = ((flags[i] == "Y") || (flags[i] == "YES"));
		}

        /// Create the Maturity Cap calculator
        maturityCapCalculator = new ARM_MaturityCapCalculator(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
			ARM_Date(myUnderlyingEndDate),
            longShort,
            capFloor,
			resetFreq,
			payFreq,
            indexTerm,
            dayCount,
            intRule,
            spread,
            initNominal,
            initTRI,
			annuity,
			(ARM_MaturityCapCalculator::ProductMode)maturityCapMode,
			coeff,
			*amortizingProfile,
			resetGap,
			resetCal,
			payCal,
			(ARM_MaturityCapCalculator::CalibrationMode)calibrationMode,
			nbIterations,		
			productsToPrice,
            *mktDataManager,
            keys);

        /// Free memory
        if(isAmortizing)
            delete amortizingProfile;
		amortizing = NULL;
        
		/// assign object
		if( !assignObject( maturityCapCalculator, result, objId ) ){
			return ARM_KO; }
		else{
			return ARM_OK; }
	}
	
	catch(Exception& x)
	{
		/// Free memory
        if(isAmortizing)
            delete amortizingProfile;

		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////////////////
//// Function to get data from a Maturity Cap Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_MaturityCap_Get(
        long maturityCapId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_MaturityCapCalculator* maturityCapCalculator;
    ARM_Object* object=NULL;

	try
	{
		maturityCapCalculator = dynamic_cast<ARM_MaturityCapCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(maturityCapId));

		if (!maturityCapCalculator)
		{
			result.setMsg ("ARM_ERR: Maturity Cap Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
        if( typeToGet == GC_ACCESS_CF_PORT )
        {
            if (typeToGet == GC_ACCESS_CF_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(maturityCapCalculator->GetCFPortfolio())).Clone();

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
            return ARMLOCAL_GC_Get(maturityCapId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}


///////////////////////////////////////////////
//// Function to set data to a Maturity Cap Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_MaturityCap_Set(
        long maturityCapId,
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
    ARM_MaturityCapCalculator* oldMaturityCap=NULL;
    ARM_MaturityCapCalculator* newMaturityCap=NULL;
    ARM_Object* object=NULL;
    string updateMsg("Calculator updated by the ");


	try
	{
		oldMaturityCap = dynamic_cast<ARM_MaturityCapCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(maturityCapId));
		if (!oldMaturityCap)
		{
			result.setMsg ("ARM_ERR: Maturity Cap Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if((port = dynamic_cast< ARM_StdPortfolio* >(object)))
        {
            if(isUpdated)
                /// Just update the initial Maturity Cap
                newMaturityCap=oldMaturityCap;
            else
                /// Clone the initial Maturity Cap
                newMaturityCap = dynamic_cast<ARM_MaturityCapCalculator *>(oldMaturityCap->Clone());

            string typeToSet(stringGetUpper(setType));
            if(typeToSet == GC_ACCESS_CF_PORT)
            {
                newMaturityCap->SetCFPortfolio(*port);
                updateMsg += "C/F portfolio";
            }

            /// Keep calibration data consistency
            newMaturityCap->UpdateCalibrationAndTimeIt();

            if(isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new Maturity Cap in ARM cache
		    if(!assignObject( newMaturityCap, result, objId ) )
            {
                delete newMaturityCap;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
        else
            return ARMLOCAL_GC_Set(maturityCapId,dataId,keys,isUpdated,result,objId);
	}
	
	catch(Exception& x)
	{
        if(!isUpdated)
            delete newMaturityCap;
		x.DebugPrint();
		ARM_RESULT();
	}

}


long ARMLOCAL_CallableSBCalculator_Create(
        double startDate,
        double endDate,
		long cpnFreq,
		long cpnDaycount,
		long cpnIndexTiming,
		const string& cpnIdxTerm,
		long cpnIndexDaycount,
		const string& cpnResetCal,
		const string& cpnPayCal,
		long cpnIntRule,
		long cpnIndexResetLag,
		long payRec,
		long CF,
		double notional,
		long notionalId,
		double fundnotional,
		long fundnotionalId,
		double Const,
		long constId,
		double lPrevCpn,
		long lPrevCpnId,
		double lNewOpt,
		long lNewOptId,
		double strikeOpt,
		long strikeOptId,
		double minCpn,
		long minCpnId,
		double maxCpn,
		long maxCpnId,
		long fundFreq,
		long fundDaycount,
		double fundCoeff,
		long fundCoeffId,
		double fundMargin,
		long fundMarginId,
		long NotifDays,
		long NonCall,
		const string& exerciseCal,
		bool CallSBOrStacky,
		long feesId,
		long NbPathBounding,
		long NbPathPricing,
		long NbMaxBucket,
		long FixSteps,
		const string& USMethod,
		const vector< string >& CalibMod,
		const vector< string >& flags,
		long mktDataManagerId,
        const vector< string >& keys,
		const string& genType1,
		const string& genType2,
		int firstNbTimes,
		int firstNbDims,
		const string& pathScheme,
		const string& pathOrder,
		const string& antithetic,
		const string& ModelType,
		const string& TriggerMode,
		int calibSwoptFreq,
		const string& regressors,
		const vector< double >& hkDatas,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CallableSnowBallCalculator* callableSBCalculator = NULL;
	
	ARM_Curve* notionalProfile=NULL;
    bool isNotional=false;

	ARM_Curve* fundnotionalProfile=NULL;
    bool isfundNotional=false;

	ARM_Curve* constProfile=NULL;
    bool isConst=false;

	ARM_Curve* lPrevCpnProfile=NULL;
    bool isLPrevCpn=false;

	ARM_Curve* lNewOptProfile=NULL;
    bool isLNewOpt=false;

	ARM_Curve* strikeOptProfile=NULL;
    bool isStrikeOpt=false;

	ARM_Curve* minCpnProfile=NULL;
    bool isMinCpn=false;

	ARM_Curve* maxCpnProfile=NULL;
    bool isMaxCpn=false;

	ARM_Curve* fundMarginProfile=NULL;
    bool isFundMargin=false;

	ARM_Curve* fundCoeffProfile=NULL;
    bool isFundCoeff=false;

	ARM_Curve* feesProfile=NULL;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Callable SnowBall Calculator" );

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
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,notional);
            notionalProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isNotional=true;
        }

		if(fundnotionalId != ARM_NULL_OBJECT)
        {
		    fundnotionalProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundnotionalId));

		    if (!fundnotionalProfile)
		    {
			    result.setMsg ("ARM_ERR: notional is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundnotional);
            fundnotionalProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isfundNotional=true;
        }

		//const
		if(constId != ARM_NULL_OBJECT)
        {
		    constProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(constId));

		    if (!constProfile)
		    {
			    result.setMsg ("ARM_ERR: const  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,Const);
            constProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isConst=true;
        }

		//lPrevCpn
		if(lPrevCpnId != ARM_NULL_OBJECT)
        {
		    lPrevCpnProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(lPrevCpnId));

		    if (!lPrevCpnProfile)
		    {
			    result.setMsg ("ARM_ERR: lPrevCpn  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,lPrevCpn);
            lPrevCpnProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isLPrevCpn=true;
        }

		//lNewOpt
		if(lNewOptId != ARM_NULL_OBJECT)
        {
		    lNewOptProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(lNewOptId));

		    if (!lNewOptProfile)
		    {
			    result.setMsg ("ARM_ERR: lNewOpt  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,lNewOpt);
            lNewOptProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isLNewOpt=true;
        }


		//strikeOpt
		if(strikeOptId != ARM_NULL_OBJECT)
        {
		    strikeOptProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(strikeOptId));

		    if (!strikeOptProfile)
		    {
			    result.setMsg ("ARM_ERR: strikeOpt  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,strikeOpt);
            strikeOptProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isStrikeOpt=true;
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

		//fundCoeff
		if(fundCoeffId != ARM_NULL_OBJECT)
        {
		    fundCoeffProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundCoeffId));

		    if (!fundCoeffProfile)
		    {
			    result.setMsg ("ARM_ERR: fundMargin  is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
			ARM_GP_Vector times(1,0.0);
			ARM_GP_Vector values(1,fundCoeff);
            fundCoeffProfile = new ARM_Curve(times,values,new ARM::ARM_LinInterpCstExtrapolDble);
            isFundCoeff=true;
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

		feesProfile = dynamic_cast<ARM_Curve *>(LOCAL_PERSISTENT_OBJECTS->GetObject(feesId));

		if (!feesProfile)
		{
			result.setMsg ("ARM_ERR: fees is not of a good type");
			return ARM_KO;
		}

	    /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager=dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}
		

		if (CalibMod.size() != 6)
		{
			result.setMsg ("ARM_ERR: calib flags size = 5");
			return ARM_KO;
		}	

		ARM_CallableSnowBallCalculator::CalibMode calibMode = (ARM_CallableSnowBallCalculator::CalibMode) ARM::ARM_ArgConv_CSBCalibMode.GetNumber(CalibMod[0]); 

		bool betaCalib;
		if (CalibMod[1] == "Y")
			betaCalib = true;
		else if (CalibMod[1] == "N")
			betaCalib = false;
		else
		{
			result.setMsg ("ARM_ERR: beta calib flag should be Y or N.");
			return ARM_KO;
		}

		bool fixBoundary;
		if (CalibMod[2] == "Y")
			fixBoundary = true;
		else if (CalibMod[2] == "N")
			fixBoundary = false;
		else
		{
			result.setMsg ("ARM_ERR: fix boundary flag should be Y or N.");
			return ARM_KO;
		}

		ARM_CallableSnowBallCalculator::ControlVariableMode controlVariableFlag = (ARM_CallableSnowBallCalculator::ControlVariableMode) ARM::ARM_ArgConv_CSBControlVariableMode.GetNumber(CalibMod[3]);

		bool fixControlVariable;
		if (CalibMod[4] == "Y")
			fixControlVariable = true;
		else if (CalibMod[4] == "N")
			fixControlVariable = false;
		else
		{
			result.setMsg ("ARM_ERR: fix control variable flag should be Y or N.");
			return ARM_KO;
		}

		bool fixBeta;
		if (CalibMod[5] == "Y")
			fixBeta = true;
		else if (CalibMod[5] == "N")
			fixBeta = false;
		else
		{
			result.setMsg ("ARM_ERR: fix beta flag should be Y or N.");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(10,false);
		productsToPrice[0] = true;

		for (size_t i = 0; i < flags.size(); ++i)
		{
			productsToPrice[i] = ((flags[i] == "Y") || (flags[i] == "YES"));
		}


		bool antitheticFlag;
		if (antithetic == "Y")
			antitheticFlag = true;
		else if (antithetic == "N")
			antitheticFlag = false;
		else
		{
			result.setMsg ("ARM_ERR: antithetic flag should be Y or N.");
			return ARM_KO;
		}

		ARM_ModelType modelType = (ARM_ModelType) ARM_ArgConv_PricingModelType.GetNumber(ModelType);

		ARM_CallableSnowBallCalculator::TriggerMode triggerMode = (ARM_CallableSnowBallCalculator::TriggerMode) ARM_ArgConv_CSBTriggerMode.GetNumber(TriggerMode);

        /// Create the Callable SnowBall calculator
		callableSBCalculator = new ARM_CallableSnowBallCalculator(
            ARM_Date(myStartDate),
            ARM_Date(myEndDate),
			cpnFreq,
			cpnDaycount,
			cpnIndexTiming,
			cpnIdxTerm,
			cpnIndexDaycount,
			cpnResetCal,
			cpnPayCal,
			cpnIntRule,
			cpnIndexResetLag,
			payRec,
			CF,
			*notionalProfile,
			*fundnotionalProfile,
			*constProfile,
			*lPrevCpnProfile,
			*lNewOptProfile,
			*strikeOptProfile,
			*minCpnProfile,
			*maxCpnProfile,
			fundFreq,
			fundDaycount,
			*fundMarginProfile,
			*fundCoeffProfile,
			NotifDays,
			NonCall,
			exerciseCal,
			CallSBOrStacky,
			*feesProfile,
			NbPathBounding,
			NbPathPricing,
			NbMaxBucket,
			FixSteps,
			USMethod,
			calibMode,
			betaCalib,
			fixBoundary,
			fixBeta,
			controlVariableFlag,
			fixControlVariable,
			productsToPrice,
			*mktDataManager,
			keys,
			genType1,
			genType2,
			firstNbTimes,
			firstNbDims,
			pathScheme,
			pathOrder,
			antitheticFlag,
			modelType,
			triggerMode,
			calibSwoptFreq,
			regressors,
			hkDatas
			);

	    /// Free memory
        if(isNotional)
            delete notionalProfile;
		notionalProfile = NULL;

		if(isfundNotional)
            delete fundnotionalProfile;
		fundnotionalProfile = NULL;

        if(isConst)
            delete constProfile;
		constProfile = NULL;

        if(isLPrevCpn)
            delete lPrevCpnProfile;
		lPrevCpnProfile = NULL;

        if(isLNewOpt)
            delete lNewOptProfile;
		lNewOptProfile = NULL;

        if(isStrikeOpt)
            delete strikeOptProfile;
		strikeOptProfile = NULL;

		if (isFundCoeff)
			delete fundCoeffProfile;
		fundCoeffProfile = NULL;
       
        if(isMinCpn)
            delete minCpnProfile;
		minCpnProfile = NULL;
       
        if(isMaxCpn)
            delete maxCpnProfile;
		maxCpnProfile = NULL;
       
        if(isFundMargin)
            delete fundMarginProfile;
		fundMarginProfile = NULL;

		/// assign object
		if( !assignObject( callableSBCalculator, result, objId ) ){
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

        if(isConst)
            delete constProfile;
		constProfile = NULL;

        if(isLPrevCpn)
            delete lPrevCpnProfile;
		lPrevCpnProfile = NULL;

        if(isLNewOpt)
            delete lNewOptProfile;
		lNewOptProfile = NULL;

        if(isStrikeOpt)
            delete strikeOptProfile;
		strikeOptProfile = NULL;

		if (isFundCoeff)
			delete fundCoeffProfile;
		fundCoeffProfile = NULL;
       
        if(isMinCpn)
            delete minCpnProfile;
		minCpnProfile = NULL;
       
        if(isMaxCpn)
            delete maxCpnProfile;
		maxCpnProfile = NULL;
       
        if(isFundMargin)
            delete fundMarginProfile;
		fundMarginProfile = NULL;

		delete callableSBCalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_CALLABLESB_Get(
        long callableSBId,
        const string& getType,
        ARM_result&	result, 
        long        objId )
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");
    ARM_CallableSnowBallCalculator* callableSBCalculator;
    ARM_Object* object=NULL;

	try
	{
		callableSBCalculator = dynamic_cast<ARM_CallableSnowBallCalculator *>(LOCAL_PERSISTENT_OBJECTS->GetObject(callableSBId));

		if (!callableSBCalculator)
		{
			result.setMsg ("ARM_ERR: Callable SnowBall Calculator is not of a good type");
			return ARM_KO;
		}

        /// Extract the internal data
        string typeToGet(stringGetUpper(getType));
 /*       if( typeToGet == GC_ACCESS_OSW_PORT ||
            typeToGet == GC_ACCESS_CF_PORT)
        {
            if(typeToGet == GC_ACCESS_OSW_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(tarnCalculator->GetOSWPortfolio())).Clone();
            else if (typeToGet == GC_ACCESS_CF_PORT)
                object = const_cast< ARM_StdPortfolio& >(*(tarnCalculator->GetCFPortfolio())).Clone();

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
*/            return ARMLOCAL_GC_Get(callableSBId,getType,result,objId );
	}
	
	catch(Exception& x)
	{
		delete object;
		x.DebugPrint();
		ARM_RESULT();
	}
}

///////////////////////////////////////////////
//// Function to set data to a CSB Calculator
///////////////////////////////////////////////
extern long ARMLOCAL_CALLABLESB_Set(
        const long& csbId,
        const long& dataId,
        const string& setType,
		const vector< string >& keys,
        const bool& isUpdated,
        ARM_result&	result, 
        long        objId )
{
	return ARMLOCAL_GC_Set(csbId,dataId,keys,isUpdated,result,objId);
}

/******************************************************************************
							Global Cap
*******************************************************************************/
long ARMLOCAL_GlobalCapCalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		double notional,
		long notionalId,
		string fundIndexType,
		int fundDayCount,
		int fundFreq,
		int fundResetGap,
		int fundPayGap,
		int fundResetTiming,
		int fundPayTiming,
		string fundResetCal,
		string fundPayCal,
		int fundAdjRule,
		int fundIntRule,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		ARM_Vector* globalCapParams,
		long pastFixingsId,
		double capFixed,
		long capFixedId,
		double capStrike,
		long capStrikeId,
		double capSpread,
		long capSpreadId,
		int	nbSteps,
		vector<string> randGenerator,
		int samplerType,
		ARM_Vector* calibParams,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId)
		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GlobalCapCalculator* globalCapCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* fundLevProfile = NULL;
    bool isFundLev=false;
	
	ARM_ReferenceValue* capLevProfile = NULL;
    bool isCapLev=false;
	
	ARM_ReferenceValue* capFixedProfile = NULL;
    bool isCapFixed=false;

	ARM_ReferenceValue* capStrikeProfile = NULL;
    bool isCapStrike=false;

	ARM_ReferenceValue* capSpreadProfile = NULL;
    bool isCapSpread=false;

	ARM_ReferenceValue* pastFixings = NULL;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Global Cap Calculator" );

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

		//Fund Lev Curve
        if (fundLevId != ARM_NULL_OBJECT)
        {
		    fundLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundLevId));

		    if (!fundLevProfile)
		    {
			    result.setMsg ("ARM_ERR: fund lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundLevProfile = new ARM_ReferenceValue(fundLev);
            isFundLev = true;
        }

		//Cap Lev Curve
        if (capLevId != ARM_NULL_OBJECT)
        {
		    capLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(capLevId));

		    if (!capLevProfile)
		    {
			    result.setMsg ("ARM_ERR: cap lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capLevProfile = new ARM_ReferenceValue(capLev);
            isCapLev = true;
        }

		//Cap Fixed Curve
        if (capFixedId != ARM_NULL_OBJECT)
        {
		    capFixedProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(capFixedId))->Clone());

		    if (!capFixedProfile)
		    {
			    result.setMsg ("ARM_ERR: Cap Fixed is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capFixedProfile = new ARM_ReferenceValue(capFixed);

            isCapFixed = true;
        }

		//Cap Strike Curve
        if (capStrikeId != ARM_NULL_OBJECT)
        {
		    capStrikeProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(capStrikeId))->Clone());

		    if (!capStrikeProfile)
		    {
			    result.setMsg ("ARM_ERR: Cap Strike is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capStrikeProfile = new ARM_ReferenceValue(capStrike);

            isCapStrike = true;
        }

		//Cap Spread Curve
        if (capSpreadId != ARM_NULL_OBJECT)
        {
		    capSpreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(capSpreadId))->Clone());

		    if (!capSpreadProfile)
		    {
			    result.setMsg ("ARM_ERR: Cap Spread is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capSpreadProfile = new ARM_ReferenceValue(capSpread);

            isCapSpread = true;
        }

		// Global Cap Value
		globalCapParams->Elt(2)/=100;

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true;
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}
		
		if( pastFixingsId >= 0)
		{
			pastFixings = (ARM_ReferenceValue*) LOCAL_PERSISTENT_OBJECTS->GetObject(pastFixingsId);

			if(LocalPersistent::LOCAL_IS_OBJECT_ROOT_CLASS_OK(pastFixings, ARM_REFERENCE_VALUE) == 0)
			{
				result.setMsg ("ARM_ERR: Past fixings are not of a good type");
				return ARM_KO;
			}
		}

        // Create the Global Cap calculator
        globalCapCalculator = new ARM_GlobalCapCalculator(ccy,
											  ARM_Date(myStartDate),
											  ARM_Date(myEndDate),
											  payReceive,
											  *notionalProfile,
											  (ARM_INDEX_TYPE) ARM_ConvIrType(fundIndexType.c_str()),
											  fundDayCount,
											  fundFreq,
											  fundResetGap,
											  fundPayGap,
											  fundResetTiming,
											  fundPayTiming,
											  fundResetCal,
											  fundPayCal,
											  fundAdjRule,
											  fundIntRule,
											  *fundLevProfile,
											  *capLevProfile,
											  *capFixedProfile,
											  *capStrikeProfile,
											  *capSpreadProfile,
											  globalCapParams,
											  pastFixings,
											  nbSteps,
											  randGenerator,
											  samplerType,
											  calibParams,
											  mdmKeys,
											  *mktDataManager,
											  productsToPrice);
		 		
        // Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;

		delete capFixedProfile;
		capFixedProfile = NULL;

		delete capStrikeProfile;
		capStrikeProfile = NULL;

		delete capSpreadProfile;
		capSpreadProfile = NULL;

		delete pastFixings;
		pastFixings = NULL;
	
		// assign object
		if ( !assignObject( globalCapCalculator, result, objId ) )
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

		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;

		delete capFixedProfile;
		capFixedProfile = NULL;

		delete capStrikeProfile;
		capStrikeProfile = NULL;

		delete capSpreadProfile;
		capSpreadProfile = NULL;
	
		x.DebugPrint();
		ARM_RESULT();
	}
}

long ARMLOCAL_GlobalCapCalculator_Create(
		long globalCapId,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		int	nbSteps,
		vector<string> randGenerator,
		int samplerType,
		ARM_Vector* calibParams,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId)
		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GlobalCapCalculator* globalCapCalculator = NULL;

	ARM_ReferenceValue* fundLevProfile = NULL;
    bool isFundLev=false;
	
	ARM_ReferenceValue* capLevProfile = NULL;
    bool isCapLev=false;

	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Global Cap Calculator" );

		ARM_GlobalCap* globalCap = NULL;
		
		globalCap = (ARM_GlobalCap*) LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(globalCap, ARM_GLOBALCAP) == 0) 
		{
			result.setMsg ("ARM_ERR: global cap is not of a good type");
			return ARM_KO;
		}

		//Fund Lev Curve
        if (fundLevId != ARM_NULL_OBJECT)
        {
		    fundLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundLevId));

		    if (!fundLevProfile)
		    {
			    result.setMsg ("ARM_ERR: fund lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundLevProfile = new ARM_ReferenceValue(fundLev);
            isFundLev = true;
        }

		//Cap Lev Curve
        if (capLevId != ARM_NULL_OBJECT)
        {
		    capLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(capLevId));

		    if (!capLevProfile)
		    {
			    result.setMsg ("ARM_ERR: cap lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capLevProfile = new ARM_ReferenceValue(capLev);
            isCapLev = true;
        }

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true;
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}
		
        // Create the Global Cap calculator
        globalCapCalculator = new ARM_GlobalCapCalculator(globalCap,
														  *fundLevProfile,
														  *capLevProfile,
														  nbSteps,
														  randGenerator,
														  samplerType,
														  calibParams,
														  mdmKeys,
														  *mktDataManager,
														  productsToPrice);
		 
        // Free memory
		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;
	
		// assign object
		if ( !assignObject( globalCapCalculator, result, objId ) )
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
		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}

// Temp interface for Global Cap using securities from Summit
long ARMLOCAL_GlobalCapCalculator_Create_WithoutMktData(
		double asOfDate,
		long globalCapId,
		double fundLev,
		long fundLevId,
		double capLev,
		long capLevId,
		ARM_result&	result, 
        long objId)
		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_GlobalCapCalculator* globalCapCalculator = NULL;

	ARM_ReferenceValue* fundLevProfile = NULL;
    bool isFundLev=false;
	
	ARM_ReferenceValue* capLevProfile = NULL;
    bool isCapLev=false;

	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Global Cap Calculator" );
		
		char myAsOfDate[11];
		Local_XLDATE2ARMDATE(asOfDate,myAsOfDate);
		
		ARM_GlobalCap* globalCap = NULL;
		
		globalCap = (ARM_GlobalCap*) LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapId);

		if (LocalPersistent::LOCAL_IS_OBJECT_CLASS_OK(globalCap, ARM_GLOBALCAP) == 0) 
		{
			result.setMsg ("ARM_ERR: global cap is not of a good type");
			return ARM_KO;
		}

		//Fund Lev Curve
        if (fundLevId != ARM_NULL_OBJECT)
        {
		    fundLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(fundLevId));

		    if (!fundLevProfile)
		    {
			    result.setMsg ("ARM_ERR: fund lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fundLevProfile = new ARM_ReferenceValue(fundLev);
            isFundLev = true;
        }

		//Cap Lev Curve
        if (capLevId != ARM_NULL_OBJECT)
        {
		    capLevProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(capLevId));

		    if (!capLevProfile)
		    {
			    result.setMsg ("ARM_ERR: cap lev is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            capLevProfile = new ARM_ReferenceValue(capLev);
            isCapLev = true;
        }
		

        // Create the Global Cap calculator
        globalCapCalculator = new ARM_GlobalCapCalculator(
											  (ARM_Date) myAsOfDate,
											  globalCap,
											  *fundLevProfile,
											  *capLevProfile);
		 		
        // Free memory
		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;
	
		// assign object
		if ( !assignObject( globalCapCalculator, result, objId ) )
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
		if (isFundLev)
			delete fundLevProfile;
		fundLevProfile = NULL;

		if (isCapLev)
			delete capLevProfile;
		capLevProfile = NULL;

		x.DebugPrint();
		ARM_RESULT();
	}
}

extern long ARMLOCAL_GlobalCap_Set(long globalCapId,
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
    ARM_GlobalCapCalculator* oldGlobalCap = NULL;
    ARM_GlobalCapCalculator* newGlobalCap = NULL;
    ARM_Object* object = NULL;
    string updateMsg("Calculator updated by the ");

	try
	{
		oldGlobalCap = dynamic_cast<ARM_GlobalCapCalculator*>(LOCAL_PERSISTENT_OBJECTS->GetObject(globalCapId));
		if (!oldGlobalCap)
		{
			result.setMsg ("ARM_ERR: Global Cap Calculator is not of a good type");
			return ARM_KO;
		}
        
		object = LOCAL_PERSISTENT_OBJECTS->GetObject(dataId);

        ARM_StdPortfolio* port = NULL;
        if ((port = dynamic_cast<ARM_StdPortfolio*>(object)))
        {
            if (isUpdated)
                /// Just update the initial CRA
                newGlobalCap = oldGlobalCap;
            else
                /// Clone the initial CRA then set
                newGlobalCap = dynamic_cast<ARM_GlobalCapCalculator*>(oldGlobalCap->Clone());

            string typeToSet(stringGetUpper(setType));
            if (typeToSet == GC_ACCESS_OSW_PORT)
                updateMsg += "OSW portfolio";
            else if (typeToSet == GC_ACCESS_SHORT_TERM_PORT)
                updateMsg += "STM vanillas portfolio";
		
            if (isUpdated)
            {
                result.setString(updateMsg.c_str());
			    return ARM_OK;
            }

            /// Assign the new CRA in ARM cache
		    if (!assignObject(newGlobalCap, result, objId))
            {
                delete newGlobalCap;
			    return ARM_KO;
            }
		    else
			    return ARM_OK;
        }
		else
            return ARMLOCAL_GC_Set(globalCapId, dataId, keys, isUpdated, result, objId);
	}
	
	catch (Exception& x)
	{
        if (!isUpdated)
            delete newGlobalCap;
		x.DebugPrint();
		ARM_RESULT();
	}
}

/******************************************************************************
							Snow Range Calculator
*******************************************************************************/
long ARMLOCAL_SnowRangeCalculator_Create(
		ARM_Currency ccy,
		double startDate,
		double endDate,
		int payReceive,
		double notional,
		long notionalId,
		string fundingIndexTerm,
		int fundingDayCount,
		string couponIndexTerm,
		int couponDayCount,
		int resetFreq,
		int payFreq,
		int resetTiming,
		int payTiming,
		int resetGap,
		string resetCal,
		string payCal,
		int adjRule,
		int intRule,
		double spread,
		long spreadId,
		double strike,
		long strikeId,
		double ratchet,
		long ratchetId,
		double cashFlow,
		long cashFlowId,
		double fixedRate,
		long fixedRateId,
		double leverage,
		long leverageId,
		ARM_Vector* snowRangeParams,
		ARM_Vector* calibParams,
		string modelName,
		ARM_Vector* modelParams,
		int	nbSteps,
		string generatorType,
		string inversionMethod,
		bool antithetic,
		int samplerType,
		vector<string> mdmKeys,
		long mktDataManagerId,
		vector<string> pricingFlags,
		ARM_result&	result, 
        long objId)		  
{
	// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	// used in the MACRO ARM_RESULT
	CCString msg ("");
	ARM_SnowRangeCalculator* snowRangeCalculator = NULL;

	ARM_ReferenceValue* notionalProfile = NULL;
    bool isNotional=false;

	ARM_ReferenceValue* spreadProfile = NULL;
    bool isSpread=false;
	
	ARM_ReferenceValue* strikeProfile = NULL;
    bool isStrike=false;
	
	ARM_ReferenceValue* ratchetProfile = NULL;
    bool isRatchet=false;

	ARM_ReferenceValue* cashFlowProfile = NULL;
    bool isCashFlow=false;

	ARM_ReferenceValue* fixedRateProfile = NULL;
    bool isFixedRate=false;

	ARM_ReferenceValue* leverageProfile = NULL;
    bool isLeverage=false;
	
	try
	{
		// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "Global Cap Calculator" );

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

		//Spread Curve
		if (spreadId != ARM_NULL_OBJECT)
		{
			spreadProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(spreadId))->Clone());
			(*spreadProfile) /= 100;

			if (!spreadProfile)
			{
				result.setMsg ("ARM_ERR: spread is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			spreadProfile = new ARM_ReferenceValue(spread);
			(*spreadProfile) /= 100;
			isSpread = true;
		}

		//Strike Curve
		if (strikeId != ARM_NULL_OBJECT)
		{
			strikeProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(strikeId))->Clone());
			(*strikeProfile) /= 100;

			if (!strikeProfile)
			{
				result.setMsg ("ARM_ERR: strike is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			strikeProfile = new ARM_ReferenceValue(strike);
			(*strikeProfile) /= 100;
			isStrike = true;
		}

		//Ratchet Curve
		if (ratchetId != ARM_NULL_OBJECT)
		{
			ratchetProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(ratchetId));

			if (!ratchetProfile)
			{
				result.setMsg ("ARM_ERR: ratchet is not of a good type");
				return ARM_KO;
			}
		}
		else
		{
			ratchetProfile = new ARM_ReferenceValue(ratchet);
			isRatchet = true;
		}

		//Cash Flow Curve
        if (cashFlowId != ARM_NULL_OBJECT)
        {
		    cashFlowProfile = dynamic_cast<ARM_ReferenceValue *>(LOCAL_PERSISTENT_OBJECTS->GetObject(cashFlowId));

		    if (!cashFlowProfile)
		    {
			    result.setMsg ("ARM_ERR: Cash Flow is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            cashFlowProfile = new ARM_ReferenceValue(cashFlow);
            isCashFlow = true;
        }

		//Fixed Rate Curve
        if (fixedRateId != ARM_NULL_OBJECT)
        {
		    fixedRateProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(fixedRateId))->Clone());
			(*fixedRateProfile) /= 100;

		    if (!fixedRateProfile)
		    {
			    result.setMsg ("ARM_ERR: Fixed Rate is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            fixedRateProfile = new ARM_ReferenceValue(fixedRate);
			(*fixedRateProfile) /= 100;
            isFixedRate = true;
        }

		//Leverage Curve
        if (leverageId != ARM_NULL_OBJECT)
        {
		    leverageProfile = dynamic_cast<ARM_ReferenceValue *>((LOCAL_PERSISTENT_OBJECTS->GetObject(leverageId))->Clone());
			
		    if (!leverageProfile)
		    {
			    result.setMsg ("ARM_ERR: Leverage is not of a good type");
			    return ARM_KO;
		    }
        }
        else
        {
            leverageProfile = new ARM_ReferenceValue(leverage);
			isLeverage = true;
        }

        // Market Data Manager Object
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(mktDataManagerId));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		std::deque<bool> /* FIXMEFRED: mig.vc8 (23/05/2007 14:27:19): dont use std::vector<bool>. see item 18 in "Effective STL" from Scott Meyers*/ productsToPrice(pricingFlags.size(),false);
		productsToPrice[0] = true;
		for (size_t i=1; i<pricingFlags.size(); ++i)
		{
			productsToPrice[i] = ((pricingFlags[i] == "Y") || (pricingFlags[i] == "YES"));
		}

        snowRangeCalculator = new ARM_SnowRangeCalculator(ccy,
														  ARM_Date(myStartDate),
														  ARM_Date(myEndDate),
														  payReceive,
														  *notionalProfile,
														  fundingIndexTerm,
														  fundingDayCount,
														  couponIndexTerm,
														  couponDayCount,
														  resetFreq,
														  payFreq,
														  resetTiming,
														  payTiming,
														  resetGap,
														  resetCal,
														  payCal,
														  adjRule,
														  intRule,
														  *spreadProfile,
														  *strikeProfile,
														  *ratchetProfile,
														  *cashFlowProfile,
														  *fixedRateProfile,
														  *leverageProfile,
														  snowRangeParams,
														  calibParams,
														  modelName,
														  modelParams,
														  nbSteps,
														  generatorType,
														  inversionMethod,
														  antithetic,
														  samplerType,
														  mdmKeys,
														  *mktDataManager,
														  productsToPrice);
		 					
        // Free memory
       	if (isNotional)
			delete notionalProfile;
		notionalProfile = NULL;

		delete spreadProfile;
		spreadProfile = NULL;

		delete strikeProfile;
		strikeProfile = NULL;

		if (isRatchet)
			delete ratchetProfile;
		ratchetProfile = NULL;

		if (isCashFlow)
			delete cashFlowProfile;
		cashFlowProfile = NULL;

		delete fixedRateProfile;
		fixedRateProfile = NULL;
	
		delete leverageProfile;
		leverageProfile = NULL;
	
		// assign object
		if ( !assignObject( snowRangeCalculator, result, objId ) )
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

		delete spreadProfile;
		spreadProfile = NULL;

		delete strikeProfile;
		strikeProfile = NULL;

		if (isRatchet)
			delete ratchetProfile;
		ratchetProfile = NULL;

		if (isCashFlow)
			delete cashFlowProfile;
		cashFlowProfile = NULL;

		delete fixedRateProfile;
		fixedRateProfile = NULL;
	
		x.DebugPrint();
		ARM_RESULT();
	}
}



//===========================================================//
// VolBondCalculator										 //
//===========================================================//
long ARMLOCAL_VolBondCalculator_Create(
		const double nominal,
		const double startDate,
		const double endDate,
		const string payFreq_s,
		const string resetFreq_s,
		const string dayCount_s,
		const string tenor_s,
		const string intrule_s,
		const string stubrule_s,
		const double resetgap,
		const string paycalendar,
		const string resetcalendar,
		const string type_odesolver_s,
		const vector<double>& RK4parameters,
		const vector<double>& RK5parameters,
		const double mcnbsteps,
		const double mcbucket,
		const double mcnbstepsbyyear,
		const vector<string>& randomgenerator,
		const long marketdatamanager_id,
		const vector<string>& marketdatamanagerkeys,
		const string payofftype,
		const vector<string>& productstoprice,
		ARM_result& result,
        long objId)
{
	/// input checks
	if( !GlobalPersistanceOk( result ) )
		return ARM_KO;

	/// used in the MACRO ARM_RESULT
	CCString msg ("");

	ARM::ARM_VolBondCalculator* mycalculator = 0;

	try
	{
		/// crm tracing
		ARM_CRMCookies.Instance()->RegisterService( ARM::ARM_USERNAME, "VolBond Calculator" );

		char myStartDate[20];
		char myEndDate[20];
		
        // Convert dates
		Local_XLDATE2ARMDATE(startDate, myStartDate);
		Local_XLDATE2ARMDATE(endDate, myEndDate);


		long payFreq_l;
		if((payFreq_l = ARM_ConvFrequency (payFreq_s.c_str(), result)) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: payFreq is not of a good type");
			return ARM_KO;
		}

		long resetFreq_l;
		if((resetFreq_l = ARM_ConvFrequency (resetFreq_s.c_str(), result)) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: resetFreq is not of a good type");
			return ARM_KO;
		}

		long dayCount_l;
		if((dayCount_l = ARM_ConvDayCount (dayCount_s.c_str())) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: dayCount is not of a good type");
			return ARM_KO;
		}


		long intRule_l;
		if ((intRule_l = ARM_ConvIntRule (intrule_s.c_str())) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: interest rule (intrule) is not of a good type");
			return ARM_KO;
		}

		long stubRule_l;
		if ((stubRule_l = ARM_ConvStubRule (stubrule_s.c_str())) == ARM_DEFAULT_ERR)
		{
			result.setMsg ("ARM_ERR: stub rule is not of a good type");
			return ARM_KO;
		}

		long type_odesolver_l;		
		if( (type_odesolver_l = ARM_ConvGPODESolverType( type_odesolver_s.c_str(), result)) == ARM_DEFAULT_ERR )
		{
			result.setMsg ("ARM_ERR: type ode solver is not of a good type");
			return ARM_KO;
		}

		ARM_GP_Vector RKparameters;
		if (type_odesolver_s != "RK5Adaptative")		
			for (int i = 0; i < RK4parameters.size(); ++i)
				RKparameters.push_back(RK4parameters[i]);
		else
			for (int i = 0; i < RK5parameters.size(); ++i)
				RKparameters.push_back(RK5parameters[i]);


		deque<bool> productstoprice_b(productstoprice.size());
		for (int i = 0; i < productstoprice.size(); ++i)
		{
			if (productstoprice[i] == "Y" || productstoprice[i] == "Yes" || productstoprice[i] == "YES")
				productstoprice_b[i] = true;
			else
				productstoprice_b[i] = false;
		}
				
		int payofftype_i;

		if (payofftype == "TypeI")
			payofftype_i = ARM::ARM_VolBondCalculator::payoffType::TypeI;
		else if (payofftype == "TypeII")
			payofftype_i = ARM::ARM_VolBondCalculator::payoffType::TypeII;
		else if (payofftype == "LookBack")
			payofftype_i = ARM::ARM_VolBondCalculator::payoffType::LookBack;
		else
		{
			result.setMsg ("ARM_ERR: payoff is not of a good type");
			return ARM_KO;
		}


        /// Restore Market Data Manager Rep
		ARM_MarketData_ManagerRep* mktDataManager = dynamic_cast< ARM_MarketData_ManagerRep* >(LOCAL_PERSISTENT_OBJECTS->GetObject(marketdatamanager_id));
		if (!mktDataManager)
		{
			result.setMsg ("ARM_ERR: market data manager is not of a good type");
			return ARM_KO;
		}

		mycalculator = 	new ARM::ARM_VolBondCalculator(
							*mktDataManager, marketdatamanagerkeys,
							nominal,							
							myStartDate, myEndDate, 
							payFreq_l, resetFreq_l, dayCount_l, tenor_s, intRule_l, stubRule_l, resetgap, paycalendar, resetcalendar,
							type_odesolver_l, RKparameters, 
							mcnbsteps, mcnbstepsbyyear, mcbucket, randomgenerator,
							payofftype_i,
							productstoprice_b);

		/// assign object
		return (assignObject( mycalculator, result, objId ) )? ARM_OK: ARM_KO;
	}
	
	catch(Exception& x)
	{
		if (0 != mycalculator)
		delete mycalculator;
		x.DebugPrint();
		ARM_RESULT();
	}
}

