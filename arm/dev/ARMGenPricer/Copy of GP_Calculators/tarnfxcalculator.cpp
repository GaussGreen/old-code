/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file tarnfxcalculator.cpp
 *
 *  \brief file for the TARN FX Calculator
 *	\author  A. Lekrafi && P. Lam
 *	\version 1.0
 *	\date August 2006
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/tarnfxcalculator.h"
#include "gpcalculators/basisconverter.h"

/// gpbase
#include "gpbase/env.h"
#include "gpbase/argconvdefault.h"
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
#include "gpbase/globalconstant.h"

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

/// gpcalib
#include "gpcalib/calibmethod.h"
#include "gpcalib/modelfitter.h"
#include "gpcalib/vanillapricer.h"


/// gpmodels
#include "gpmodels/2irfxModel.h"
#include "gpmodels/1irfxModel.h"
#include "gpmodels/NP1IRNFX.h"
#include "gpmodels/q1f.h"
#include "gpmodels/q1f_fx.h"
#include "gpmodels/Smiled_fx.h"
#include "gpmodels/mixture_fx.h"
#include "gpmodels/modelparamsq1f.h"
#include "gpmodels/eqfx_modelfactory.h"
#include "gpmodels/ModelParams_EqFxBase.h"
#include "gpmodels/forwardmarginbasis.h"

/// gpnummethods
#include "gpnummethods/scheduler.h"
#include "gpnummethods/meanrevertingsampler.h"
#include "gpnummethods/pathscheme.h"
#include "gpnummethods/pathschemefactory.h"
#include "gpnummethods/mcmethod.h"
#include "gpnummethods/sampler.h"
#include "gpnummethods/schedulerfactory.h"
#include "gpnummethods/samplerfactory.h"

/// gpnumlib
#include "gpnumlib/randomgenfactory.h"
#include "gpnumlib/transposer.h"

/// kernel
#include <crv/volflat.h>
#include <inst/forex.h>
#include <inst/swaption.h>
#include <inst/option.h>
#include <inst/portfolio.h>
#include <inst/swap.h>
#include <util/fromto.h>


/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )

/// H&W vol range [10bp,500bp]
const double HWVOL_LOWER_BOUND      = 0.00001;
const double HWVOL_UPPER_BOUND      = 0.05;

/// Q vol range [2.5%,100%]
const double QVOL_LOWER_BOUND       = 0.025;
const double QVOL_UPPER_BOUND       = 1.0;

/// 0.001bp of vega to be selected in portfolio for volatility bootstrapping 
const double IR_VEGA_MIN_TO_SELECT=0.0000001;
const double OSW_DEFAULT_WEIGHT=1.0;
const double OSW_DEFAULT_PRICE=1.0e+100;

/// 10-3 bp of vega to be selected in portfolio for FX volatility bootstrapping 
const double FX_VEGA_MIN_TO_SELECT  = 1.0e-7;
const double FX_DEFAULT_WEIGHT      = 1.0;
const double FX_DEFAULT_PRICE       = 1.0e+30;

/// Default MDM key names
const string YC_KEY_NAME		= "YC_";
const string YC_BASIS_KEY_NAME	= "YC_BASIS_";
const string FOREX_KEY_NAME		= "FOREX_";
const string OSWMODEL_KEY_NAME	= "OSWMOD_";
const string FXMODEL_KEY_NAME	= "FXMOD_";
const string CORREL_KEY_NAME	= "CORREL_";
const string MRS_KEY_NAME		= "MRS_";
const string QFX_KEY_NAME		= "Q_";

/// Reference schedules for TARN date structure
const unsigned int RF_CPNFX_SCHED = 0;
const unsigned int NB_TARNFX_SCHED = 1;

const string ARM_TARNFXCalculator::TARNFXColNamesTable [] =
{
    "ResetDate",
	"NextResetDate",
	"StartDate",
	"EndDate",
	"FundingStartDate",
	"FundingEndDate",
	"IT",
	"RedemptionResetDate",
	"MinCpn",
	"MaxCpn",
	"Notional",
	"Target",
	"DF",
	"Offset",
	"Width",
	"Funding",
	"InitFX",
	"DomCpn",
	"FgnCpn",
	"FxSpot",
	"Coupon",
	"PaidCoupon",
	"SumCoupons",
	"Fees",
	"IsAlive",
	"RealCoupon",
	"DiscountFunding",
	"RealFunding",
	"Cap",
	"Floor",
	"RedemptionStrike",
	"Redemption",
	"RealRedemption",	
	"TARN",	
	"Duration",
	"Proba",
	"InitFX2",
	"DomCpn2",
	"FgnCpn2",
	"RedemptionStrike2",
	"InitFX3",
	"DomCpn3",
	"FgnCpn3",
	"RedemptionStrike3"
};


const int ARM_TARNFXCalculator::ProductToPriceColumns[] =
{
	ARM_TARNFXCalculator::Funding,
	ARM_TARNFXCalculator::Coupon,
	ARM_TARNFXCalculator::DF,
	ARM_TARNFXCalculator::PaidCoupon,
	ARM_TARNFXCalculator::RealCoupon,
	ARM_TARNFXCalculator::IsAlive,
	ARM_TARNFXCalculator::Cap,
	ARM_TARNFXCalculator::Floor,
	ARM_TARNFXCalculator::Redemption,
	ARM_TARNFXCalculator::RealRedemption,
	ARM_TARNFXCalculator::RealFunding,
	ARM_TARNFXCalculator::Duration,
	ARM_TARNFXCalculator::Proba,
	ARM_TARNFXCalculator::TARN,
};


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::ARM_TARNFXCalculator( const ARM_TARNFXCalculator& rhs )
	:
	ARM_HybridIRFXCalculator(rhs),
	itsvFundIndex(rhs.itsvFundIndex),
	itsPayOffName(rhs.itsPayOffName)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::~ARM_TARNFXCalculator()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
/*ARM_TARNFXCalculator::ARM_TARNFXCalculator(
							 const ARM_Date& asOfDate,
							 const ARM_Date& startDate,
							 const ARM_Date& endDate,
							 const ARM_Currency& DomCcy,
							 const ARM_Currency& ForCcy,
							 const ARM_Currency* FundCcy,
							 int payRec,
							 int cpnDayCount,
							 int cpnFreq,
							 int cpnResetGap,
							 const string& cpnResetCal,
							 const string& cpnPayCal,
							 int stubRule,
							 int cpnTiming,
							 int cpnIntRule,
							 const ARM_Curve* cpnNominal,
							 const ARM_Curve* domesticCpn,
							 const ARM_Curve* foreignCpn,
							 const ARM_Curve* MinCpn,
							 const ARM_Curve* MaxCpn,
							 const ARM_Curve* InitialFX,
							 int fundFreq,
							 int fundDayCount,
							 const ARM_Curve* fundNominal,
							 const ARM_Curve* fundSpread,
							 double targetRedemption,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 double redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_IntVector& productsToPrice)
:	
	ARM_HybridIRFXCalculator(asOfDate,
							 startDate,
							 endDate,
							 DomCcy,
							 ForCcy,
							 FundCcy,
							 payRec,
							 cpnDayCount,
							 cpnFreq,
							 cpnResetGap,
							 cpnResetCal,
							 cpnPayCal,
							 stubRule,
							 cpnTiming,
							 cpnIntRule,
							 cpnNominal,
							 domesticCpn,
							 foreignCpn,
							 MinCpn,
							 MaxCpn,
							 InitialFX,
							 fundFreq,
							 fundDayCount,
							 fundNominal,
							 fundSpread,
							 redemptionType,
							 redemptionGap,
							 redemptionStrike,
							 &fees),
	itsCpnIntRule(cpnIntRule),
	itsTargetRedemption(targetRedemption),
	itsFXChoice(FXChoice),
	itsIntermediatePrices(intermediatePrices),
	itsProductsToPrice(productsToPrice)
{
	itsNbFX = 1;
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

//	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	bool	needOtherPayoff = (itsIntermediatePrices?true:false);

	string payModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());

	itsColumnsToPrice = ProductsToPriceColumnNames();

    /// Create the Generic Security
	CreateAndSetDealDescription(payModelName, itsColumnsToPrice, CreateCstManager(), false, needOtherPayoff);
}*/
	
/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_TARNFXCalculator::ARM_TARNFXCalculator(
							 const ARM_Date& asOfDate,
							 const ARM_Date& startDate,
							 const ARM_Date& endDate,
							 const ARM_Currency& DomCcy,
							 const vector<ARM_Currency>& ForCcy,
							 const ARM_Currency& FundCcy,
							 int payRec,
							 int cpnDayCount,
							 int cpnFreq,
							 int cpnResetGap,
							 const string& cpnResetCal,
							 const string& cpnPayCal,
							 int stubRule,
							 int cpnTiming,
							 int cpnIntRule,
							 const ARM_Curve& cpnNominal,
							 const vector<ARM_Curve>& domesticCpn,
							 const vector<ARM_Curve>& foreignCpn,
							 const ARM_Curve& MinCpn,
							 const ARM_Curve& MaxCpn,
							 const vector<ARM_Curve>& InitialFX,
							 int fundFreq,
							 int fundDayCount,
							 const ARM_Curve& fundNominal,
							 const ARM_Curve& fundSpread,
							 double targetRedemption,
							 const string& FXChoice,
							 int redemptionType,
							 int redemptionGap,
							 ARM_GP_Vector& redemptionStrike,
							 const ARM_Curve& fees,
							 int intermediatePrices,
							 const ARM_IntVector& productsToPrice,
							 const ARM_FXTARNPayoffType& payOffName)
:	
	ARM_HybridIRFXCalculator(asOfDate,
							 startDate,
							 endDate,
							 DomCcy,
							 ForCcy,
							 FundCcy,
							 payRec,
							 cpnDayCount,
							 cpnFreq,
							 cpnResetGap,
							 cpnResetCal,
							 cpnPayCal,
							 stubRule,
							 cpnTiming,
							 cpnIntRule,
							 cpnNominal,
							 domesticCpn,
							 foreignCpn,
							 MinCpn,
							 MaxCpn,
							 InitialFX,
							 fundFreq,
							 fundDayCount,
							 fundNominal,
							 fundSpread,
							 redemptionType,
							 redemptionGap,
							 redemptionStrike,
							 fees),
	itsCpnIntRule(cpnIntRule),
	itsTargetRedemption(targetRedemption),
	itsFXChoice(FXChoice),
	itsvFundIndex(0),
	itsIntermediatePrices(intermediatePrices),
	itsProductsToPrice(productsToPrice),
	itsPayOffName(payOffName)
{
	itsNbFX = ForCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

//	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	bool	needOtherPayoff = (itsIntermediatePrices?true:false);

	string payModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());

	itsColumnsToPrice = ProductsToPriceColumnNames();

    /// Create the Generic Security
	CreateAndSetDealDescription(payModelName, itsColumnsToPrice, CreateCstManager(), false, needOtherPayoff);
}

//////////////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Init
///	Returns: void
///	Action : initialise calculator with market data and model parameters
//////////////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::Init(
			int nbSimul,
			int bucketSize,
			int randGenType1,
			int randGenAlgo1,
			int randGenType2,
			int randGenAlgo2,
			int firstNbDims,
			int firstNbTimes,
			int factorNb,
			int timeStepNb,
			int spaceStepNb,
			double stdDevNb,
			int skipPDE,
			int rescalling,
			int modelType,
			int smileFlag,
			int mixCalib,
			int oneFactorFlag,
			int correlType,
			const ARM_MarketData_ManagerRep& mktDataManager)
{
	SetMktDataManager(CreateClonedPtr(const_cast<ARM_MarketData_ManagerRep*>(&mktDataManager)));

	itsNbSimul = nbSimul;
	itsBucketSize = bucketSize;
	itsRandGenType1 = randGenType1;
	itsRandGenAlgo1 = randGenAlgo1;
	itsRandGenType2 = randGenType2;
	itsRandGenAlgo2 = randGenAlgo2;
	itsFirstNbDims = firstNbDims;
	itsFirstNbTimes = firstNbTimes;
	itsFactorsNb = factorNb;
	itsTimeStepNb = timeStepNb;
	itsSpaceStepNb = spaceStepNb;
	itsStdDevNb = stdDevNb;
	itsSkipPDE = skipPDE;
	itsRescalling = rescalling;
	
	if ((itsNbFX > 1) && (modelType != ModelNP1IRNFX))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : TARN FX Chooser : model type should be 'NP1IRNFX'");

	itsModelType = modelType;
	itsSmileFlag = smileFlag;
	itsMixCalib = mixCalib;
	itsOneFactorFlag = oneFactorFlag;
	itsCorrelType = correlType;

	/// Check market datas
    CheckMktDataAndTimeIt();

	Update(&mktDataManager);
}

//////////////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: Init
///	Returns: void
///	Action : initialise calculator with market data and model parameters
//////////////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::Init(
			int nbSimul,
			int bucketSize,
			int randGenType1,
			int randGenAlgo1,
			int randGenType2,
			int randGenAlgo2,
			int firstNbDims,
			int firstNbTimes,
			int factorNb,
			int timeStepNb,
			int spaceStepNb,
			double stdDevNb,
			int skipPDE,
			int rescalling,
			int modelType,
			int smileFlag,
			int mixCalib,
			int oneFactorFlag,
			int correlType,
			vector<ARM_ZeroCurve*> zeroCurves,
			vector<ARM_ZeroCurve*> basisCurves,
			vector<ARM_Forex*> forex,
			vector<ARM_VolCurve*> ATMVol, //for swopt BSGen
			vector<ARM_VolCurve*> fxVol, //for BS fx models
			vector<ARM_ParamsMixture_Fx*> mixtureParams, //for mixture fx models
			vector<ARM_CurveModelParam*> mrsParams,
			vector<ARM_CurveModelParam*> QParams,
			ARM_GP_Matrix* correlMatrix)
{
	itsNbSimul = nbSimul;
	itsBucketSize = bucketSize;
	itsRandGenType1 = randGenType1;
	itsRandGenAlgo1 = randGenAlgo1;
	itsRandGenType2 = randGenType2;
	itsRandGenAlgo2 = randGenAlgo2;
	itsFirstNbDims = firstNbDims;
	itsFirstNbTimes = firstNbTimes;
	itsFactorsNb = factorNb;
	itsTimeStepNb = timeStepNb;
	itsSpaceStepNb = spaceStepNb;
	itsStdDevNb = stdDevNb;
	itsSkipPDE = skipPDE;
	itsRescalling = rescalling;
	itsNbFX = zeroCurves.size()-1;

	if ((itsNbFX > 1) && (modelType != ModelNP1IRNFX))
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : TARN FX Chooser : model type should be 'NP1IRNFX'");

	itsModelType = modelType;
	itsSmileFlag = smileFlag;
	itsMixCalib = mixCalib;
	itsOneFactorFlag = oneFactorFlag;
	itsCorrelType = correlType;

	ARM_Date asof = zeroCurves[0]->GetAsOfDate();

	//Market Data Manager ----------------------------------------------------------->
	ARM_MarketData_ManagerRep* marketDataManager = &*GetMktDataManager();
	marketDataManager->SetAsOfDate(asof);

	ARM_StringVector keys; //variable size (depends on the number of currencies)

	int i;
	string domCcyName(zeroCurves[0]->GetCurrencyUnit()->GetCcyName());
	for (i=0; i<itsNbFX+1; i++) // 1 dom + N fgn
	{
		string ccyName(zeroCurves[i]->GetCurrencyUnit()->GetCcyName());

		keys.push_back(YC_KEY_NAME + ccyName);
		marketDataManager->RegisterData(YC_KEY_NAME + ccyName, zeroCurves[i]);

		keys.push_back(YC_BASIS_KEY_NAME + ccyName);
		marketDataManager->RegisterData(YC_BASIS_KEY_NAME + ccyName, basisCurves[i]);

		keys.push_back(OSWMODEL_KEY_NAME + ccyName);
		ARM_BSModel* swoptModel = new ARM_BSNorModel(asof, zeroCurves[i], NULL, ATMVol[i], ATMVol[i]);
		marketDataManager->RegisterData(OSWMODEL_KEY_NAME + ccyName, swoptModel);
		delete swoptModel;

		keys.push_back(MRS_KEY_NAME + ccyName);
		marketDataManager->RegisterData(MRS_KEY_NAME + ccyName, mrsParams[i]);

		if (i==0)
		{
			keys.push_back(FOREX_KEY_NAME + domCcyName + "/" + domCcyName);
			marketDataManager->RegisterData(FOREX_KEY_NAME + domCcyName + "/" + domCcyName, forex[i]);
		}
		else //fgn and fx part
		{
			keys.push_back(FOREX_KEY_NAME + ccyName + "/" + domCcyName);
			marketDataManager->RegisterData(FOREX_KEY_NAME + ccyName + "/" + domCcyName, forex[i]);

			keys.push_back(FXMODEL_KEY_NAME + ccyName + "/" + domCcyName);
			ARM_Object* Model_FX = NULL;
			if (mixtureParams.size() > 0)
			{
				// mixture
				Model_FX = new ARM_MixtureModel_Fx (CreateClonedPtr( basisCurves[0] ), 
													CreateClonedPtr( basisCurves[i] ), 
													forex[i]->GetMarketPrice(), //spot
													mixtureParams[i-1]);
			}
			else
			{
				// BS
				Model_FX = new ARM_BSModel (asof,
											forex[i]->GetMarketPrice(), //spot
											basisCurves[i],
											basisCurves[0],
											fxVol[i-1]);
			}
			marketDataManager->RegisterData(FXMODEL_KEY_NAME + ccyName + "/" + domCcyName, Model_FX);
			delete Model_FX;
			
			keys.push_back(QFX_KEY_NAME + ccyName + "/" + domCcyName);
			marketDataManager->RegisterData(QFX_KEY_NAME + ccyName + "/" + domCcyName, QParams[i-1]);
		}
	}
	keys.push_back(string("CORREL_DOMCCY_FORCCY"));
	marketDataManager->RegisterData(string("CORREL_DOMCCY_FORCCY"), correlMatrix);

	SetKeys(keys);

	Update(marketDataManager);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateCstManager
///	Returns: void
///	Action : create the const manager (static data).
/////////////////////////////////////////////////////////////////
ARM_CstManagerPtr ARM_TARNFXCalculator::CreateCstManager()
{
	ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

	vector<string>	cstNames(3,"");
	cstNames[0] = "FundDateStrip";
	cstNames[1] = "FundNominal";
	cstNames[2] = "FundMargin";

	vector<ARM_GramFctorArg> cstVector;
	
	cstVector.push_back(ARM_GramFctorArg(itsFundDateStrip));

	size_t fundSize = itsFundDateStrip->size();

	ARM_GP_VectorPtr fundNominal(new ARM_GP_Vector(fundSize));
	ARM_GP_VectorPtr fundMargin(new ARM_GP_Vector(fundSize));

	ARM_GP_Vector* resetDates = itsFundDateStrip->GetResetDates();

	size_t i;
	for (i = 0; i < fundSize; ++i)
	{
		(*fundNominal)[i] = itsCpnNominalCv.Interpolate((*resetDates)[i]-asOfDate.GetJulian());
		(*fundMargin)[i] = itsFundSpreadCv.Interpolate((*resetDates)[i]-asOfDate.GetJulian());
	}

	cstVector.push_back(ARM_GramFctorArg(fundNominal));
	cstVector.push_back(ARM_GramFctorArg(fundMargin));

	ARM_CstManagerPtr cstManagerPtr = ARM_CstManagerPtr(new ARM_CstManager(cstNames, cstVector));

	return cstManagerPtr;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CheckData
///	Returns: void
///	Action : check if TARN FX data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::CheckData()
{
	if (itsStartDate >= itsEndDate)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The start date frequency should be before the end date.");

	if (itsCpnFreq > itsFundFreq)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : The funding frequency should be greater or equal than the coupon frequency.");

	if (itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption)
		ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Dual Option is not supported.");
	
	itsvFundIndex.clear();
/// last: compute indexes to relate exer - funding - cpn		
	ARM_GP_Vector* exerciseDates = itsCpnDateStrip->GetResetDates() ;
	ARM_GP_Vector* fundStartDates = itsFundDateStrip->GetFlowStartDates();
	ARM_GP_Vector* fundResetDates = itsFundDateStrip->GetResetDates() ;
	
	size_t exSize = itsCpnDateStrip->size()-1;
	size_t i = 0;
	for (size_t k(0); k<exSize; k++)
	{
		while ((*exerciseDates)[k]>(*fundStartDates)[i]) i++;
		if((*exerciseDates)[k]>(*fundResetDates)[i])
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : no call date allowed between fixing dates and start dates");
		itsvFundIndex.push_back(i);
	}
	itsvFundIndex.push_back(fundStartDates->size());
	/* 2 call dates cannot belong to the same period */
	/* this may happen when freq (call) > freq (exotic) or freq (funding) */
	for (k=0; k<exSize; k++)
	{
		if (itsvFundIndex[k]>=itsvFundIndex[k+1])
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : 2 call dates found within the same funding leg period. Freq (call) > Freq (Funding) is forbidden");
	}

	/* after each call call date, we must have same start dates for exo and funding legs */
	ARM_GP_Vector* cpnStartDates  = itsCpnDateStrip->GetFlowStartDates() ;
	for (k=0; k<exSize; k++) 
	{		
		if (fabs((*fundStartDates)[itsvFundIndex[k]]-(*cpnStartDates)[k+1])> ARM_GlobalConstant::ARM_SEVENDAYS_LAG)
			ARM_THROW( ERR_INVALID_ARGUMENT, "TARN calculator : call date  defines a swap with mismatching exo and floating start dates (no broken period allowed)");
	}
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CheckMktData
///	Returns: void
///	Action : check if TARN FX market data are consistent
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::CheckMktData()
{
	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
    if(!domCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcDomKey + " is expected in the Market Data Manager");
    string domCcy(domCurve->GetCurrencyUnit()->GetCcyName());

	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
    if(!basisDomCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisDomKey + " is expected in the Market Data Manager");
    string basisDomCcy(basisDomCurve->GetCurrencyUnit()->GetCcyName());

    if(domCcy != basisDomCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Domestic(=coupon) Basis Curve currency should consistent with reference curve");

	string YcFundKey = YC_KEY_NAME + itsFundingCcy.GetCcyName();
	string YcBasisFundKey = YC_BASIS_KEY_NAME + itsFundingCcy.GetCcyName();

	ARM_ZeroCurve* fundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcFundKey));
    if(!fundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcFundKey + " is expected in the Market Data Manager");
    string fundCcy(fundCurve->GetCurrencyUnit()->GetCcyName());

	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisFundKey));
    if(!basisFundCurve)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisFundKey + " is expected in the Market Data Manager");
    string basisFundCcy(basisFundCurve->GetCurrencyUnit()->GetCcyName());

    if(fundCcy != basisFundCcy)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Funding Basis Curve currency should consistent with reference curve");

	string FundForexKey = FOREX_KEY_NAME + itsFundingCcy.GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_Forex* fundforex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(FundForexKey));
	if(!fundforex)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + FundForexKey + " is expected in the Market Data Manager");

	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	ARM_ModelParam* mrsDomParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsDomKey) );
    if(!mrsDomParam || mrsDomParam->GetType() != ARM_ModelParamType::MeanReversion)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : MRS Param for key=" + MrsDomKey + " is expected in the Market Data Manager");

	string OswDomModelKey = OSWMODEL_KEY_NAME + itsDomesticCcy.GetCcyName();
    ARM_BSModel* oswDomBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswDomModelKey) );
    if(!oswDomBSModel)
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : domestic swaption B&S model for key=" + OswDomModelKey + " is expected in the Market Data Manager");

	string YcForKey;
	string YcBasisForKey;
	string ForexKey;
	string MrsForKey;
	string QFxKey;
	string CorrelMatrixKey;
	string OswForModelKey;
	string FxModelKey;
	int i = 0;

	for (i=0; i<itsNbFX; i++)
	{
		YcForKey = YC_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey));
		if(!forCurve)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcForKey + " is expected in the Market Data Manager");
		string forCcy(forCurve->GetCurrencyUnit()->GetCcyName());

		YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
		if(!basisForCurve)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Yield Curve for key=" + YcBasisForKey + " is expected in the Market Data Manager");
		string basisForCcy(basisForCurve->GetCurrencyUnit()->GetCcyName());

		if(forCcy != basisForCcy)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Basis Curve currency should consistent with reference curve");

		ForexKey = FOREX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));
		if(!forex)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : forex=" + ForexKey + " is expected in the Market Data Manager");

		MrsForKey = MRS_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_ModelParam* mrsForParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(MrsForKey) );
		if(!mrsForParam || mrsForParam->GetType() != ARM_ModelParamType::MeanReversion)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign MRS Param for key=" + MrsForKey + " is expected in the Market Data Manager");

		QFxKey = QFX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_ModelParam* qFxParam = dynamic_cast< ARM_ModelParam* >( GetMktDataManager()->GetData(QFxKey) );
		if(!qFxParam || qFxParam->GetType() != ARM_ModelParamType::QParameter)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : Foreign Q Fx Param for key=" + QFxKey + " is expected in the Market Data Manager");

		CorrelMatrixKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "," + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")";
		
		OswForModelKey = OSWMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName();
		ARM_BSModel* oswForBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(OswForModelKey) );
		if(!oswForBSModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : foreign swaption B&S model for key=" + OswForModelKey + " is expected in the Market Data Manager");

		FxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(FxModelKey) );
		ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(FxModelKey) );

		if (!fxMixModel && (itsModelType==Model1IRFX))
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : 1IRFX model needs fx Mixture model for key.");

		if(!fxMixModel && !fxBSModel)
			ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx Mixture or BS model for key=" + FxModelKey + " is expected in the Market Data Manager");
	}
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNFXCalculator::ColumnNames() const
{
    size_t	colNamesSize = sizeof(TARNFXColNamesTable)/sizeof(TARNFXColNamesTable[0]);
    vector<string>				colNamesVec(colNamesSize);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0; i<colNamesSize; ++i)
        colNamesVec[i] = TARNFXColNamesTable[i];

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: ProductsToPriceColumnNames
///	Returns: ARM_StringVector
///	Action : Create Products to price column names
/////////////////////////////////////////////////////////////////
ARM_StringVector ARM_TARNFXCalculator::ProductsToPriceColumnNames()
{
	ARM_StringVector columnNames;
	columnNames.reserve(NbProductsToPrice);

    for(size_t i = 0; i < NbProductsToPrice; ++i)
		if ((itsProductsToPrice.size() > i) && itsProductsToPrice[i])
			columnNames.push_back( TARNFXColNamesTable[ ARM_TARNFXCalculator::ProductToPriceColumns[i] ] );

	return columnNames;
}

/*
//////////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the TARN FX.
//////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_TARNFXCalculator::DatesStructure() const
{
    const char*	resetCalendar = itsCpnResetCal.c_str();
    const char*	payCalendar   = itsCpnPayCal.c_str();

    int fwdRule		= K_MOD_FOLLOWING;	// for forward dates
    int intRule		= itsCpnIntRule;	// for interest dates
    int stubRule	= itsStubRule;
	int resetFreq   = itsCpnFreq;
    int resetTiming = itsCpnTiming;
	int payFreq     = itsCpnFreq;
	int payGap      = GETDEFAULTVALUE;
    int payTiming   = K_ARREARS;

	int spotDays = itsDomesticCcy.GetSpotDays();
	int fundResetTiming = K_ADVANCE;

	if (itsFundDateStrip.IsNull())
	{
		itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(itsStartDate, itsEndDate,itsFundFreq, itsFundDayCount,resetCalendar,fwdRule, intRule, 
								stubRule,-spotDays, itsFundFreq,payGap, payCalendar,fundResetTiming, payTiming));
	}
	
	if (itsCpnDateStrip.IsNull())
	{
		itsCpnDateStrip = ARM_DateStripPtr(new ARM_DateStrip(itsStartDate, itsEndDate, resetFreq, itsCpnDayCount,resetCalendar, fwdRule, intRule, 
								stubRule,-fabs(itsCpnResetGap), payFreq, payGap, payCalendar,	resetTiming, payTiming));

		// For the in arrears case we create a fake start date to get the first 
		// reset date of the funding leg
		if (itsCpnTiming == K_ARREARS)
		{
			double  firstFundingResetDate = (*itsFundDateStrip->GetResetDates())[0];
			itsCpnDateStrip->InsertDate(0, 0, 0, 0, 0, firstFundingResetDate, 0, 0, 0);
		}	
	}

	

    /// Merge schedules on "ResetDate"
    ARM_DateStripVector SchedVect(NB_TARNFX_SCHED, NULL);
    SchedVect[RF_CPNFX_SCHED] = &*itsCpnDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect, "ResetDate");

	// Memorize first event index
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	int	firstEventIdx = 0;
	while( (firstEventIdx < eventDates->size()) && ((*eventDates)[firstEventIdx] < asOfDate) )
		++firstEventIdx;

	const_cast<ARM_TARNFXCalculator*>(this)->itsFirstEventIdx = firstEventIdx;


    return	EventSchedule;
}
*/

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: InitPriceableColumns
///	Returns: void
///	Action : initialise to 0 column to be able to sum
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::InitPriceableColumns(vector<string>& rowDescVec, vector<ARM_GP_VALUE_TYPE>& rowTypeVec) const
{
    string	zeroValue("0");

	rowDescVec[ARM_TARNFXCalculator::Funding] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Funding] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Coupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Coupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::DF] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::DF] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::PaidCoupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::PaidCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealCoupon] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealCoupon] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::IsAlive] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::IsAlive] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Cap] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Cap] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Floor] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Floor] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Redemption] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Redemption] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealRedemption] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealRedemption] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::RealFunding] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::RealFunding] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Duration] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Duration] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::Proba] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::Proba] = ARM_DOUBLE;

	rowDescVec[ARM_TARNFXCalculator::TARN] = zeroValue;
    rowTypeVec[ARM_TARNFXCalculator::TARN] = ARM_DOUBLE;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_TARNFXCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_DateStripPtr	vDateStrip = datesStructure.GetDateStrip(0);
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	size_t	eventSize = vDateStrip->GetResetDates()->size();
	size_t	descSize = sizeof(TARNFXColNamesTable)/sizeof(TARNFXColNamesTable[0]);

	vector<string>				rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE>	rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec, rowTypeVec);
	int i = 0;

	bool isArrearsCpn = (itsCpnTiming == K_ARREARS);
	bool isRealFirstEvent = (eventIdx == 0);
	bool isFirstEvent = (eventIdx == itsFirstEventIdx);
	bool isFirstCpnEvent = ((itsCpnTiming == K_ARREARS)&&(itsFirstEventIdx==0)?(eventIdx == 1):(eventIdx == itsFirstEventIdx));
	bool isLastEvent = (eventIdx == (eventSize-1));

    CC_Ostringstream	resetDateDesc;
    resetDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetResetDates()))[eventIdx];
    rowDescVec[TARNFXColAlias::ResetDate] = resetDateDesc.str();
    rowTypeVec[TARNFXColAlias::ResetDate] = ARM_DATE_TYPE;

    CC_Ostringstream	nextResetDateDesc;
	if( !(isArrearsCpn && isLastEvent) )
	{
		if(!isArrearsCpn)
			nextResetDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetResetDates()))[eventIdx];
		else if(!isLastEvent)
			nextResetDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetResetDates()))[eventIdx+1];

		rowDescVec[TARNFXColAlias::NextResetDate] = nextResetDateDesc.str();
		rowTypeVec[TARNFXColAlias::NextResetDate] = ARM_DATE_TYPE;

		CC_Ostringstream	fundingStartDateDesc;
		if(!isArrearsCpn)
			fundingStartDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowStartDates()))[eventIdx];
		else
			fundingStartDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowStartDates()))[eventIdx+1];
		rowDescVec[TARNFXColAlias::FundingStartDate] = fundingStartDateDesc.str();
		rowTypeVec[TARNFXColAlias::FundingStartDate] = ARM_DATE_TYPE;

		CC_Ostringstream	fundingEndDateDesc;
		if(!isArrearsCpn)
			fundingEndDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowEndDates()))[eventIdx];
		else
			fundingEndDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowEndDates()))[eventIdx+1];
		rowDescVec[TARNFXColAlias::FundingEndDate] = fundingEndDateDesc.str();
		rowTypeVec[TARNFXColAlias::FundingEndDate] = ARM_DATE_TYPE;

		CC_Ostringstream	fundingDesc;
		fundingDesc << "SWAP(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< TARNFXColNamesTable[TARNFXColAlias::FundingStartDate] << "[i]" << ","
						<< TARNFXColNamesTable[TARNFXColAlias::FundingEndDate] << "[i]" << ","
						<< "0" << ","
						<< ARM_ArgConvReverse_RcvOrPay.GetString(itsPayRec) << ","
						<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsCpnFreq) << ","
						<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsCpnDayCount) << ","
						<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsFundFreq) << ","
						<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsFundDayCount) << ","
						<< "FundMargin,FundNominal,,,FundDateStrip,FundDateStrip," 
						<< TARNFXColNamesTable[TARNFXColAlias::Offset] << "[i]," 
						<< TARNFXColNamesTable[TARNFXColAlias::Offset] << "[i],"
						<< TARNFXColNamesTable[TARNFXColAlias::Width] << "[i],"
						<< TARNFXColNamesTable[TARNFXColAlias::Width] << "[i])";

		rowDescVec[TARNFXColAlias::Funding] = fundingDesc.str();
		rowTypeVec[TARNFXColAlias::Funding] = ARM_STRING;

		CC_Ostringstream	discountFundingDesc;
		discountFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::Funding] << "[i]";
		discountFundingDesc << "/" << "DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
							<< TARNFXColNamesTable[TARNFXColAlias::NextResetDate] << "[i])";

		rowDescVec[TARNFXColAlias::DiscountFunding] = discountFundingDesc.str();
		rowTypeVec[TARNFXColAlias::DiscountFunding] = ARM_STRING;

		CC_Ostringstream offsetDesc;
		offsetDesc << itsvFundIndex[eventIdx] << endl;

		rowDescVec[TARNFXColAlias::Offset] = offsetDesc.str();
		rowTypeVec[TARNFXColAlias::Offset] = ARM_DOUBLE;

		size_t width = itsvFundIndex[eventIdx+1] -itsvFundIndex[eventIdx];
		CC_Ostringstream widthDesc; 
		widthDesc << CC_NS(std, fixed) << width;

		rowDescVec[TARNFXColAlias::Width] = widthDesc.str();
		rowTypeVec[TARNFXColAlias::Width] = ARM_DOUBLE;
	}

	if( !(isArrearsCpn && isRealFirstEvent) )
	{
		CC_Ostringstream	startDateDesc;
		double	vFlowStartDate = (*(vDateStrip->GetFlowStartDates()))[eventIdx];
		startDateDesc << CC_NS(std, fixed) << vFlowStartDate;
		rowDescVec[TARNFXColAlias::StartDate] = startDateDesc.str();
		rowTypeVec[TARNFXColAlias::StartDate] = ARM_DATE_TYPE;

		CC_Ostringstream	endDateDesc;
		endDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowEndDates()))[eventIdx];
		rowDescVec[TARNFXColAlias::EndDate] = endDateDesc.str();
		rowTypeVec[TARNFXColAlias::EndDate] = ARM_DATE_TYPE;

		CC_Ostringstream	ITDesc;
		ITDesc << CC_NS(std, fixed) << (*(vDateStrip->GetInterestTerms()))[eventIdx];
		rowDescVec[TARNFXColAlias::IT] = ITDesc.str();
		rowTypeVec[TARNFXColAlias::IT] = ARM_DOUBLE;

		CC_Ostringstream	redemptionResetDateDesc;
		if(isLastEvent)
		{
			redemptionResetDateDesc << CC_NS(std, fixed) << itsRedemptionResetDate.GetJulian();
			rowDescVec[TARNFXColAlias::RedemptionResetDate] = redemptionResetDateDesc.str();
			rowTypeVec[TARNFXColAlias::RedemptionResetDate] = ARM_DATE_TYPE;
		}

		CC_Ostringstream	minCpnDesc;
		minCpnDesc << CC_NS(std, fixed) << itsMinCpnCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::MinCpn] = minCpnDesc.str();
		rowTypeVec[TARNFXColAlias::MinCpn] = ARM_DOUBLE;

		CC_Ostringstream	maxCpnDesc;
		maxCpnDesc << CC_NS(std, fixed) << itsMaxCpnCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::MaxCpn] = maxCpnDesc.str();
		rowTypeVec[TARNFXColAlias::MaxCpn] = ARM_DOUBLE;

		CC_Ostringstream	notionalDesc;
		notionalDesc << CC_NS(std, fixed) << itsCpnNominalCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::Notional] = notionalDesc.str();
		rowTypeVec[TARNFXColAlias::Notional] = ARM_DOUBLE;

		CC_Ostringstream	targetRedemptionDesc; 
		targetRedemptionDesc << CC_NS(std, fixed) << itsTargetRedemption;
		rowDescVec[TARNFXColAlias::Target] = targetRedemptionDesc.str();
		rowTypeVec[TARNFXColAlias::Target] = ARM_DOUBLE;

		CC_Ostringstream	feesDesc;
		feesDesc << CC_NS(std, fixed) << itsFees.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::Fees] = feesDesc.str();
		rowTypeVec[TARNFXColAlias::Fees] = ARM_DOUBLE;

		CC_Ostringstream	DFDesc;
		DFDesc	<< "DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
				<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")";

		rowDescVec[TARNFXColAlias::DF] = DFDesc.str();
		rowTypeVec[TARNFXColAlias::DF] = ARM_STRING;

		CC_Ostringstream	initialFXDesc;
		initialFXDesc << CC_NS(std, fixed) << itsInitialFXCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::InitFX] = initialFXDesc.str();
		rowTypeVec[TARNFXColAlias::InitFX] = ARM_DOUBLE;

		CC_Ostringstream	domesticCpnDesc;
		domesticCpnDesc << CC_NS(std, fixed) << itsDomesticCpnCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::DomCpn] = domesticCpnDesc.str();
		rowTypeVec[TARNFXColAlias::DomCpn] = ARM_DOUBLE;

		CC_Ostringstream	foreignCpnDesc;
		foreignCpnDesc << CC_NS(std, fixed) << itsForeignCpnCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[TARNFXColAlias::FgnCpn] = foreignCpnDesc.str();
		rowTypeVec[TARNFXColAlias::FgnCpn] = ARM_DOUBLE;


		CC_Ostringstream	fxSpotDesc;
		fxSpotDesc	<< "SPOT(FOREX_" << itsForeignCcy[0].GetCcyName() << "/" << itsDomesticCcy.GetCcyName() << ")";

		rowDescVec[TARNFXColAlias::FxSpot] = fxSpotDesc.str();
		rowTypeVec[TARNFXColAlias::FxSpot] = ARM_STRING;


		CC_Ostringstream	couponDesc;
		couponDesc	<< "Min(Max(" ;
		if ((itsNbFX == 1) || (itsFXChoice == "First") || itsPayOffName == ARM_TARNFXPayoffType::TARNFX)
		{
			couponDesc	<< TARNFXColNamesTable[TARNFXColAlias::FgnCpn] << "[i]" 
						<< "*SPOT(FOREX_" << itsForeignCcy[0].GetCcyName() << "/" 
						<< itsDomesticCcy.GetCcyName() << ")/"
						<< TARNFXColNamesTable[TARNFXColAlias::InitFX] << "[i]" << "-"
						<< TARNFXColNamesTable[TARNFXColAlias::DomCpn] << "[i]" << ",";
		}
		else if (itsFXChoice == "Second")
		{
			couponDesc	<< TARNFXColNamesTable[TARNFXColAlias::FgnCpn2] << "[i]" 
						<< "*SPOT(FOREX_" << itsForeignCcy[1].GetCcyName() << "/" 
						<< itsDomesticCcy.GetCcyName() << ")/" 
						<< TARNFXColNamesTable[TARNFXColAlias::InitFX2] << "[i]" << "-"
						<< TARNFXColNamesTable[TARNFXColAlias::DomCpn2] << "[i]" << ",";
		}
		else
		{
			couponDesc << (itsFXChoice == "Worst" ?	"Min(" : "Max(" );
			if (i==0)
			{
				couponDesc	<< TARNFXColNamesTable[TARNFXColAlias::FgnCpn] << "[i]" 
							<< "*SPOT(FOREX_" << itsForeignCcy[i].GetCcyName() << "/" 
							<< itsDomesticCcy.GetCcyName() << ")/"
							<< TARNFXColNamesTable[TARNFXColAlias::InitFX] << "[i]" << "-"
							<< TARNFXColNamesTable[TARNFXColAlias::DomCpn] << "[i],";
			}
			for (i=1; i<itsNbFX; i++)
			{
				couponDesc	<< TARNFXColNamesTable[TARNFXColAlias::FgnCpn2 + 4*(i-1)] << "[i]"
							<< "*SPOT(FOREX_" << itsForeignCcy[i].GetCcyName() << "/" 
							<< itsDomesticCcy.GetCcyName() << ")/"
							<< TARNFXColNamesTable[TARNFXColAlias::InitFX2 + 4*(i-1)] << "[i]" << "-"
							<< TARNFXColNamesTable[TARNFXColAlias::DomCpn2 + 4*(i-1)] << "[i]";

				couponDesc <<  (i == itsNbFX-1 ?  ")," : ",");
			}
		}
		couponDesc	<< TARNFXColNamesTable[TARNFXColAlias::MinCpn] << "[i]" << "),"
					<< TARNFXColNamesTable[TARNFXColAlias::MaxCpn] << "[i]" << ")*"
					<< TARNFXColNamesTable[TARNFXColAlias::IT] << "[i]";
	
		rowDescVec[TARNFXColAlias::Coupon] = couponDesc.str();
		rowTypeVec[TARNFXColAlias::Coupon] = ARM_STRING;

		CC_Ostringstream	paidCouponDesc;
		paidCouponDesc	<< TARNFXColNamesTable[TARNFXColAlias::Coupon] << "[i]"
						<< "*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")*"
						<< TARNFXColNamesTable[TARNFXColAlias::Notional] << "[i]";

		if(isLastEvent)
		{
			paidCouponDesc << "-" << TARNFXColNamesTable[TARNFXColAlias::RealRedemption] << "[i]";
		}

		rowDescVec[TARNFXColAlias::PaidCoupon] = paidCouponDesc.str();
		rowTypeVec[TARNFXColAlias::PaidCoupon] = ARM_STRING;

		CC_Ostringstream	sumCouponsDesc;
		sumCouponsDesc << TARNFXColNamesTable[TARNFXColAlias::Coupon] << "[i]";
		if (!isFirstCpnEvent)
		{
			sumCouponsDesc	<< "+" << TARNFXColNamesTable[TARNFXColAlias::SumCoupons] << "[i-1]";
		}
		rowDescVec[TARNFXColAlias::SumCoupons] = sumCouponsDesc.str();
		rowTypeVec[TARNFXColAlias::SumCoupons] = ARM_STRING;
		
		bool isPRDKO = itsPayOffName == ARM_TARNFXPayoffType::PRDKO;
		int index  = isPRDKO ? TARNFXColAlias::FxSpot : TARNFXColAlias::SumCoupons;

		CC_Ostringstream	isAliveDesc;
		isAliveDesc	<< "IF(" << TARNFXColNamesTable[index] << "[i]+" 
						<< TARNFXColNamesTable[TARNFXColAlias::Fees] << "[i]<"
						<< TARNFXColNamesTable[TARNFXColAlias::Target] << "[i],1,0)";
		if(!isFirstCpnEvent)
		{
			isAliveDesc	<< "*" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]";
		}

		rowDescVec[TARNFXColAlias::IsAlive] = isAliveDesc.str();
		rowTypeVec[TARNFXColAlias::IsAlive] = ARM_STRING;

		CC_Ostringstream	realCouponDesc;
		realCouponDesc	<< TARNFXColNamesTable[TARNFXColAlias::PaidCoupon] << "[i]";
		if(!isFirstCpnEvent)
		{
			realCouponDesc	<< "*" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]";
		}

		rowDescVec[TARNFXColAlias::RealCoupon] = realCouponDesc.str();
		rowTypeVec[TARNFXColAlias::RealCoupon] = ARM_STRING;

		CC_Ostringstream	realFundingDesc;

		if(!isFirstCpnEvent)
		{
			realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*";
		}

		if (isArrearsCpn)
		{
			if (itsFirstEventIdx < eventIdx)
				realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::DiscountFunding] << "[i-1]";
			else
				realFundingDesc	<< "0";
		}
		else
			realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::DiscountFunding] << "[i]";
/*
		if (!isArrearsCpn)
		{
			if(isFirstCpnEvent)
				realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::DiscountFunding] << "[i]";
			else
				realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*"
								<< TARNFXColNamesTable[TARNFXColAlias::DiscountFunding] << "[i]";
		}
		else
		{
			if(isFirstCpnEvent)
				realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::Funding] << "[i-1]";
			else
				realFundingDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*"
								<< TARNFXColNamesTable[TARNFXColAlias::DiscountFunding] << "[i-1]";
		}
*/
		rowDescVec[TARNFXColAlias::RealFunding] = realFundingDesc.str();
		rowTypeVec[TARNFXColAlias::RealFunding] = ARM_STRING;

		CC_Ostringstream	capDesc;

		if(!isFirstCpnEvent)
		{
			capDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]*";
		}

		capDesc	<< "Max("
				<< TARNFXColNamesTable[TARNFXColAlias::SumCoupons] << "[i]" << "-"
				<< TARNFXColNamesTable[TARNFXColAlias::Target] << "[i]"
				<< ",0)*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
				<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")*"
				<< TARNFXColNamesTable[TARNFXColAlias::Notional] << "[i]";

		rowDescVec[TARNFXColAlias::Cap] = capDesc.str();
		rowTypeVec[TARNFXColAlias::Cap] = ARM_STRING;

		CC_Ostringstream	floorDesc;		
		if(isLastEvent)
		{
			floorDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i]*"
						<< "Max("
						<< TARNFXColNamesTable[TARNFXColAlias::Target] << "[i]" << "-"
						<< TARNFXColNamesTable[TARNFXColAlias::SumCoupons] << "[i]"
						<< ",0)*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")*"
						<< TARNFXColNamesTable[TARNFXColAlias::Notional] << "[i]";
			rowDescVec[TARNFXColAlias::Floor] = floorDesc.str();
			rowTypeVec[TARNFXColAlias::Floor] = ARM_STRING;
		}
		
		CC_Ostringstream	redemptionStrikeDesc;
		double	vRedemptionStrike = 0.;
		if(isLastEvent)
			vRedemptionStrike = itsRedemptionStrike.Elt(0);
		redemptionStrikeDesc << CC_NS(std, fixed) << vRedemptionStrike;
		rowDescVec[TARNFXColAlias::RedemptionStrike] = redemptionStrikeDesc.str();
		rowTypeVec[TARNFXColAlias::RedemptionStrike] = ARM_DOUBLE;

		CC_Ostringstream	redemptionDesc;
		if(isLastEvent)
		{
			if (itsNbFX > 1) // NP1IRNFX
			{
				if (itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)
				{
					if (itsFXChoice == "First")
					{
						redemptionDesc	<< "(Spot(" << FOREX_KEY_NAME << itsForeignCcy[0].GetCcyName()
										<< "/" << itsDomesticCcy.GetCcyName()
										<< ")/" << TARNFXColNamesTable[TARNFXColAlias::RedemptionStrike] << "[i]" << "-1)";
					}
					else if (itsFXChoice == "Second")
					{
						redemptionDesc	<< "(Spot(" << FOREX_KEY_NAME << itsForeignCcy[1].GetCcyName()
										<< "/" << itsDomesticCcy.GetCcyName()
										<< ")/" << TARNFXColNamesTable[TARNFXColAlias::RedemptionStrike2] << "[i]" << "-1)";
					}
					else
					{
						redemptionDesc	<< "Min(Spot(" << FOREX_KEY_NAME << itsForeignCcy[0].GetCcyName()
										<< "/" << itsDomesticCcy.GetCcyName()
										<< ")/" << TARNFXColNamesTable[TARNFXColAlias::RedemptionStrike] << "[i]" << "-1,";
						
						for (i=1; i<itsNbFX; i++)
						{
							redemptionDesc	<< "(Spot(" << FOREX_KEY_NAME << itsForeignCcy[i].GetCcyName()
											<< "/" << itsDomesticCcy.GetCcyName()
											<< ")/" << TARNFXColNamesTable[TARNFXColAlias::RedemptionStrike2+4*(i-1)] << "[i]" << "-1)";

							if (i == itsNbFX-1)
								redemptionDesc	<< ")";
							else
								redemptionDesc	<< ",";
						}
					}
					redemptionDesc	<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
									<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i])*"
									<< TARNFXColNamesTable[TARNFXColAlias::Notional] << "[i]";
				}
				else if (itsRedemptionType == ARM_PRCSRedemptionType::dualOptionRedemption)
				{
				}
				else
				{
					redemptionDesc	<< "0";
				}
			}
			else
			{
				if (itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption)
				{
					redemptionDesc	<< "-(FWD(" << FOREX_KEY_NAME << itsForeignCcy[0].GetCcyName() //TMP !
									<< "/" << itsDomesticCcy.GetCcyName() << ","
									<< TARNFXColNamesTable[TARNFXColAlias::RedemptionResetDate] << "[i]" << ", ,"
									<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")"
									<< "/" << TARNFXColNamesTable[TARNFXColAlias::RedemptionStrike] << "[i]" << "-1)"
									<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
									<< TARNFXColNamesTable[TARNFXColAlias::EndDate] << "[i]" << ")"
									<< "*" << TARNFXColNamesTable[TARNFXColAlias::Notional] << "[i]";
				}
				else
				{
					redemptionDesc	<< "0";
				}
			}
			rowDescVec[TARNFXColAlias::Redemption] = redemptionDesc.str();
			rowTypeVec[TARNFXColAlias::Redemption] = ARM_STRING;
		}

		CC_Ostringstream	realRedemptionDesc;
		if(isLastEvent)
		{
			if ((itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption) || (itsNbFX > 1))
			{
				if (itsNbFX > 1)
					realRedemptionDesc	<< "-";
				
				realRedemptionDesc	<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i]*"
									<< TARNFXColNamesTable[TARNFXColAlias::Redemption] << "[i]";
			}
			else
			{
				realRedemptionDesc	<< "0";
			}

			rowDescVec[TARNFXColAlias::RealRedemption] = realRedemptionDesc.str();
			rowTypeVec[TARNFXColAlias::RealRedemption] = ARM_STRING;
		}
		

		CC_Ostringstream	durationDesc;
		if(isFirstCpnEvent)
		{
			durationDesc	<< "UNPAY(" << TARNFXColNamesTable[TARNFXColAlias::IT] << "[i]" << ")";
		}
		else
		{
			durationDesc	<< "UNPAY(" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*"
							<< TARNFXColNamesTable[TARNFXColAlias::IT] << "[i]" << ")";
		}
		rowDescVec[TARNFXColAlias::Duration] = durationDesc.str();
		rowTypeVec[TARNFXColAlias::Duration] = ARM_STRING;

		CC_Ostringstream	probaDesc;
		if(isFirstCpnEvent)
		{
			probaDesc	<< "UNPAY(1-" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i]" << ")";
		}
		else if(isLastEvent)
		{
			probaDesc	<< "UNPAY(" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*(1-"
						<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i]" << ")"
						<< "+" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i])";
		}
		else
		{
			probaDesc	<< "UNPAY(" << TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i-1]" << "*(1-"
							<< TARNFXColNamesTable[TARNFXColAlias::IsAlive] << "[i]" << "))";
		}
		rowDescVec[TARNFXColAlias::Proba] = probaDesc.str();
		rowTypeVec[TARNFXColAlias::Proba] = ARM_STRING;

		CC_Ostringstream	tarnDesc;
		if(itsPayRec == K_PAY)
		{
			tarnDesc	<< TARNFXColNamesTable[TARNFXColAlias::RealFunding] << "[i]-" << TARNFXColNamesTable[TARNFXColAlias::RealCoupon] << "[i]";
		}
		else
		{
			tarnDesc	<< TARNFXColNamesTable[TARNFXColAlias::RealCoupon] << "[i]-" << TARNFXColNamesTable[TARNFXColAlias::RealFunding] << "[i]";
		}
		rowDescVec[TARNFXColAlias::TARN] = tarnDesc.str();
		rowTypeVec[TARNFXColAlias::TARN] = ARM_STRING;

		int i;
		for (i=1; i<itsNbFX; i++)
		{
			CC_Ostringstream	initialFXDesc;
			initialFXDesc << CC_NS(std, fixed) << itsInitialFXCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[TARNFXColAlias::InitFX2 + 4*(i-1)] = initialFXDesc.str();
			rowTypeVec[TARNFXColAlias::InitFX2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	domesticCpnDesc;
			domesticCpnDesc << CC_NS(std, fixed) << itsDomesticCpnCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[TARNFXColAlias::DomCpn2 + 4*(i-1)] = domesticCpnDesc.str();
			rowTypeVec[TARNFXColAlias::DomCpn2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	foreignCpnDesc;
			foreignCpnDesc << CC_NS(std, fixed) << itsForeignCpnCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[TARNFXColAlias::FgnCpn2 + 4*(i-1)] = foreignCpnDesc.str();
			rowTypeVec[TARNFXColAlias::FgnCpn2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	redemptionStrikeDesc;
			double	vRedemptionStrike = 0.;
			if(isLastEvent)
				vRedemptionStrike = itsRedemptionStrike.Elt(i);
			redemptionStrikeDesc << CC_NS(std, fixed) << vRedemptionStrike;
			rowDescVec[TARNFXColAlias::RedemptionStrike2 + 4*(i-1)] = redemptionStrikeDesc.str();
			rowTypeVec[TARNFXColAlias::RedemptionStrike2 + 4*(i-1)] = ARM_DOUBLE;
		}
	}

	return ARM_RowInfo(rowDescVec, rowTypeVec);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateModelNP1IRNFX
///	Returns: ARM_PricingModelPtr
///	Action : create the 1IR+nFX model
/////////////////////////////////////////////////////////////////
ARM_NP1IRNFXModel* ARM_TARNFXCalculator::CreateModelNP1IRNFX()
{
	const double initIRVol = 0.01;
	const double initFXVol = 0.1;
	const double initIRQ = 0.0;
	const double initFXQ = 1.0;
	ARM_GP_Vector lowerBound(1, 0.0001);
	ARM_GP_Vector upperBound(1, 0.05);
	ARM_GP_Vector lowerQ(1, -10.0);
	ARM_GP_Vector upperQ(1, 10.0);

	/// models :
	/// Dom + DomBS + N x {For, ForBasis, FX}
	int nbModels = 2 + itsNbFX*3;

    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    vector< ARM_PricingModelPtr > models(nbModels);

	/// Domestic parameters
	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();

	/// foreign and FX parameters
	vector<string> YcForKey, YcBasisForKey, ForexKey, MrsForKey, QFxKey;
	int i = 0;
	for (i=0; i<itsNbFX; i++)
	{
		YcForKey.push_back(YC_KEY_NAME + itsForeignCcy[i].GetCcyName());
		YcBasisForKey.push_back(YC_BASIS_KEY_NAME + itsForeignCcy[i].GetCcyName());
		ForexKey.push_back(FOREX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName());
		MrsForKey.push_back(MRS_KEY_NAME + itsForeignCcy[i].GetCcyName());
		QFxKey.push_back(QFX_KEY_NAME + itsForeignCcy[i].GetCcyName() + "/" + itsDomesticCcy.GetCcyName());
	}

	/// Correl parameter
	string CorrelMatrixKey = CORREL_KEY_NAME + "DOMCCY_FORCCY";

	/// domestic model
    ARM_ModelParamVector domModelParams(3);
	ARM_GP_Vector breakPointTimes(1, 0.0);
	ARM_GP_Vector values(1, initIRVol );

    ARM_CurveModelParam HWVolParam( ARM_ModelParamType::QVol,
									&values,
									&breakPointTimes,
									"QVOL",
									"STEPUPRIGHT",
									&lowerBound,
									&upperBound);

    domModelParams[0] = &HWVolParam;
	domModelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsDomKey));
	values[0] = initIRQ;
    domModelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter,
												&values,
												&breakPointTimes,
												"Q",
												"STEPUPRIGHT",
												&lowerQ,
												&upperQ);
	bool degenerateInHW = true;

	ARM_ZeroCurve* domCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
	names[0]  = YcDomKey;
	models[0] = ARM_PricingModelPtr( new ARM_QModel1F(CreateClonedPtr(domCurve), ARM_ModelParamsQ1F(domModelParams), degenerateInHW) );

	ARM_ZeroCurve* basisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	names[2*itsNbFX+1]  = YcBasisDomKey;
	models[2*itsNbFX+1] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(basisDomCurve)) );

	/// foreign and FX models
	for (i=0; i<itsNbFX; i++)
	{
		ARM_ModelParamVector forModelParams(3);
		values[0] = initIRVol;
		ARM_CurveModelParam forHWVolParam ( ARM_ModelParamType::QVol,
											&values,
											&breakPointTimes,
											"QVOL",
											"STEPUPRIGHT",
											&lowerBound,
											&upperBound);
		forModelParams[0] = &forHWVolParam;
		forModelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsForKey[i]));
		values[0] = initIRQ;
		forModelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter,
													&values,
													&breakPointTimes,
													"Q",
													"STEPUPRIGHT",
													&lowerQ,
													&upperQ);
		
		ARM_ZeroCurve* forCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey[i]));
		names[i+1]  = YcForKey[i];
		models[i+1] = ARM_PricingModelPtr( new ARM_QModel1F(CreateClonedPtr(forCurve), ARM_ModelParamsQ1F(forModelParams), degenerateInHW) );

		ARM_ModelParamVector fxModelParams(3);
		values[0] = initFXVol;
		ARM_CurveModelParam fxHWVolParam(ARM_ModelParamType::QVol,
										&values,
										&breakPointTimes,
										"QVOL",
										"STEPUPRIGHT",
										&lowerQ,
										&upperQ);
		fxModelParams[0] = &fxHWVolParam;
		fxModelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsForKey[i]));
		values[0] = initFXQ;
		fxModelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter,
													&values,
													&breakPointTimes,
													"Q",
													"STEPUPRIGHT",
													&lowerQ,
													&upperQ);

		double spot = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey[i]))->GetMarketPrice();
		ARM_ZeroCurve* basisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey[i]));
		names[i+1+itsNbFX]  = ForexKey[i];
		models[i+1+itsNbFX] = ARM_PricingModelPtr(ARM_EqFx_ModelFactory.Instance()->CreateModel(CreateClonedPtr(basisDomCurve), 
																						fxModelParams, 
																						spot, 
																						CreateClonedPtr(basisForCurve)));

		names[i+1+2*itsNbFX+1]  = YcBasisForKey[i];
		models[i+1+2*itsNbFX+1] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(basisForCurve)) );
	}

	int nbCcy = itsNbFX+1;
	int BasisModelIdx = 2*nbCcy-1;
	depends[BasisModelIdx] = ARM_StringVector(1, names[0]);

	for (i = 0; i < nbCcy-1; ++i)
	{
		depends[nbCcy+i] = ARM_StringVector(2);
		if (nbModels > BasisModelIdx)
		{
			depends[nbCcy+i][0]	= names[BasisModelIdx];
			depends[nbCcy+i][1]	= names[BasisModelIdx+1+i];
			depends[BasisModelIdx+1+i] = ARM_StringVector(1,names[1+i]);
		}
		else
		{
			depends[nbCcy+i][0]	= names[0];
			depends[nbCcy+i][1]	= names[BasisModelIdx+1+i];
		}
	}

   	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(CorrelMatrixKey));

	/// We create and set the brand new 1 IR + n FX model !
    ARM_NP1IRNFXModel* hybridModel = new ARM_NP1IRNFXModel(modelMap, *correlMatrix);

	return hybridModel;
}

/*
/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions for domestic market
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_TARNFXCalculator::CreateDiagonalSwaption(ARM_Currency& ccy)
{
	// Get the standard osw expiry dates
	string oswModelKey = OSWMODEL_KEY_NAME + ccy.GetCcyName();
	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(oswModelKey) );

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();
    ARM_VolCurve* oswBSVol = oswBSModel->GetVolatility();
    size_t nbExp = oswBSVol->GetExpiryTerms()->GetSize();
    ARM_GP_Vector stdExp(nbExp+1);
	size_t i;
    stdExp[0] = asOfDate.GetJulian();
    for(i=0;i<nbExp;++i)
        stdExp[i+1] = asOfDate.GetJulian() + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

	// Get the model dates
	int nbFund = itsFundDateStrip->size();

	ARM_GP_Vector* resetDates = itsFundDateStrip->GetResetDates();
	ARM_GP_Vector* startDates = itsFundDateStrip->GetFlowStartDates();

	ARM_INDEX_TYPE indexType = ccy.GetVanillaIndexType();
	int resetFreq = ccy.GetFixedPayFreq();

	ARM_Swap* stdSwap;
	ARM_Swaption* swaption;

	list< ARM_Security* > swaptionList;

	int curStdExpIdx = 0;
	bool insertSwaption = true;

    /// Analyse event dates
    for(i=0;i<nbFund;++i)
    {
		while ( (curStdExpIdx < nbExp-1) && ( stdExp[curStdExpIdx] < (*resetDates)[i]+K_NEW_DOUBLE_TOL) )
		{
			curStdExpIdx++;
			insertSwaption = true;
		}

		if (insertSwaption)
		{
			ARM_Date swapStartDate( (*startDates)[i] );

			// To Prevent stubs
			ARM_Date mathcSwapEndDate((*startDates)[i]);
			int nbMonths = CC_Round((itsEndDate.GetJulian()-(*startDates)[i])/30.5);
			mathcSwapEndDate.AddMonths(nbMonths);

			stdSwap  = new ARM_Swap(
				swapStartDate,
				mathcSwapEndDate,
				indexType,
				0.0,
				K_MARKET_RATE,
				K_RCV,
				resetFreq,
				resetFreq,
				&ccy);

			ARM_Date expiryDate((*resetDates)[i]);
			swaption = new ARM_Swaption(stdSwap,K_RCV,K_EUROPEAN,K_MARKET_RATE,expiryDate);

			delete stdSwap;

			swaptionList.push_back(static_cast< ARM_Security* >(swaption));

			insertSwaption = false;
		}
    }

    /// Built OSW portfolio and set default weight
    /// A default price is set to pass the validation test in CalibMethod constructor
    ARM_StdPortfolio* port = new ARM_StdPortfolio(swaptionList);
    for(i=0;i<port->size();++i)
    {
        port->SetWeight(OSW_DEFAULT_WEIGHT,i);
        port->SetPrice((i+1)*OSW_DEFAULT_PRICE,i);
    }

    return ARM_StdPortfolioPtr(port);
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: ComputeSwaptionPrice
///	Returns: void
///	Action : Compute the swaption price
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::ComputeSwaptionPrice(ARM_Currency& ccy, ARM_CalibMethod* IRVolCalibMethod)
{
	string oswModelKey = OSWMODEL_KEY_NAME + ccy.GetCcyName();

	ARM_BSModel* oswBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(oswModelKey) );
	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    /// Restore calibration portfolios
    ARM_StdPortfolioPtr oswPortfolio = IRVolCalibMethod->GetPortfolio();

    double price,vega,weight;
    ARM_Swaption* swaption;
    size_t i;
    bool isNotSelected;

    size_t nbOSW = oswPortfolio->GetSize();


    ARM_GP_Vector initTimes(nbOSW);
    ARM_GP_Vector initVols(nbOSW);
    double optMat,swapMat,volATM,nominal,swapRate;

    for(i=0;i<nbOSW;++i)
    {
        swaption=static_cast< ARM_Swaption* >(oswPortfolio->GetAsset(i));
	    swaption->SetModel(oswBSModel);
		swapRate = swaption->CptMarketSwapRate();

        price=swaption->ComputePrice();
        vega=swaption->ComputeSensitivity(K_VEGA);

		nominal = swaption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(oswPortfolio->GetWeights()))[i];

        isNotSelected = vega < IR_VEGA_MIN_TO_SELECT*nominal;

        oswPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        oswPortfolio->SetPrecision(0.001*vega,i);
        oswPortfolio->SetPrice(price,i);

        /// Vol initialisation
        optMat  = (swaption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
        swapMat = (swaption->GetEndDate().GetJulian() - swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
        volATM  = oswBSModel->ComputeVol(optMat,swapMat,swapRate,swapRate)/100.0;

		initTimes[i] = optMat*K_YEAR_LEN;
		initVols[i] = volATM * swapRate;
    }

    /// Replace the sigma param with new initialisations
    ARM_GP_Vector volLowerBound(nbOSW,HWVOL_LOWER_BOUND);
    ARM_GP_Vector volUpperBound(nbOSW,HWVOL_UPPER_BOUND);
    ARM_CurveModelParam* vol = new ARM_CurveModelParam(
		ARM_ModelParamType::QVol,
		&initVols,
		&initTimes,
        "QVOL",
		"STEPUPRIGHT",
		&volLowerBound,
		&volUpperBound);

	size_t volParamSize=IRVolCalibMethod->GetCalibParams().size();
    if(volParamSize==0)
        IRVolCalibMethod->GetCalibParams().push_back(vol);
    else if(volParamSize==1)
    {
        delete IRVolCalibMethod->GetCalibParam(0);
        (IRVolCalibMethod->GetCalibParams())[0] = vol;
    }
    else
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : diagonal swaption calibrate only volatilities");

}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateIRCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_TARNFXCalculator::CreateIRCalibMethod(
			const ARM_StdPortfolioPtr& diagonalSwaptionPF,
            int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=true;

	ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());

    /// Build an empty bootstrap calib method (filled latter for domestic volatility calibration)
    ARM_CalibMethod* volCalib = new ARM_CalibMethod(
									diagonalSwaptionPF,
									emptyCalibParam, //calibParam,
                                    ARM_CalibMethodType::Bootstrap1D,
									ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,noPreviousMethod,
                                    isCalibShared,modelIdx);


    return volCalib;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_TARNFXCalculator::CreateFxOption(ARM_Currency& foreignCcy)
{   
    ARM_Option* fxOption;
    list < ARM_Security* > fxOptionList;

    /// Restore forex
	string forexKey = FOREX_KEY_NAME + foreignCcy.GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	ARM_Forex* forex = static_cast< ARM_Forex* >( GetMktDataManager()->GetData(forexKey) );
	
	ARM_GP_Vector* resetDates = itsCpnDateStrip->GetResetDates();
	size_t nbReset = resetDates->size();
    
	size_t i;
    for(i=(itsCpnTiming==K_ARREARS?1:0); i<nbReset; ++i)
    {
		ARM_Date resetDate((*resetDates)[i]);

        fxOption    = new ARM_Option(forex,resetDate,K_MARKET_RATE,K_CALL,K_EUROPEAN);
        fxOptionList.push_back(static_cast< ARM_Security* >(fxOption));
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
///	Class  : ARM_TARNFXCalculator
///	Routine: ComputeFxOptionPrice
///	Returns: nothing
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::ComputeFxOptionPrices(ARM_CalibMethod* FXVolCalibMethod, int Ccyi)
{
    string fxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[Ccyi].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[Ccyi].GetCcyName();
	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[Ccyi].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();

    ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(fxModelKey) );
	ARM_BSModel* fxBSModel = dynamic_cast< ARM_BSModel* >( GetMktDataManager()->GetData(fxModelKey) );
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));

	const ARM_Date asOfDate = GetMktDataManager()->GetAsOfDate();

    ARM_StdPortfolioPtr fxPortfolio = FXVolCalibMethod->GetPortfolio();
    size_t nbFx = fxPortfolio->GetSize();
    size_t volParamSize=FXVolCalibMethod->GetCalibParams().size();

    ARM_GP_Vector initTimes(nbFx);
    ARM_GP_Vector initVols(nbFx);

    ARM_Option* fxOption;
    double optTime,payTime,vol,price,vega=0.0,nominal,weight,fwd,df;
    bool isNotSelected;
    size_t strikeIdx=0;

	double fxSpot = forex->GetMarketPrice();

	int spotDays = itsDomesticCcy.GetSpotDays();

	ARM_Date spotDate(asOfDate);
	spotDate.NextBusinessDay(spotDays);

	double dfRatio = ycBasisDomCurve->DiscountPrice((spotDate.GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN)/ycBasisForCurve->DiscountPrice((spotDate.GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN);

    for(size_t i(0);i<nbFx;++i)
    {
		fxOption = static_cast< ARM_Option* >(fxPortfolio->GetAsset(i));
		optTime = (fxOption->GetExpiryDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;
		payTime = (fxOption->GetForwardDate().GetJulian()-asOfDate.GetJulian())/K_YEAR_LEN;

		fwd = dfRatio*ycBasisForCurve->DiscountPrice(payTime)/ycBasisDomCurve->DiscountPrice(payTime)*fxSpot;
		fxOption->SetStrike(fwd);

		if (fxMixModel)
		{
			price=ARM_VanillaPricer::Price(fxOption,fxMixModel);
		}
		else
		{
			fxOption->SetModel(fxBSModel);
			price=fxOption->ComputePrice();
		}

        nominal = fxOption->GetAmount()->CptReferenceValue(0.0); // Cst RefValue
        weight  = (*(fxPortfolio->GetWeights()))[i];
	
		df = ycBasisDomCurve->DiscountPrice(payTime);
		vol = VanillaImpliedVol_BS(fwd,fwd,optTime,price/df,K_CALL,NULL,NULL);
		BS(fwd,fwd,optTime,vol,K_CALL,NULL,NULL,&vega,NULL);
		vega /= 100;

        isNotSelected = vega < FX_VEGA_MIN_TO_SELECT * nominal;

        fxPortfolio->SetWeight(isNotSelected ? 0.0 : weight,i);
        fxPortfolio->SetPrecision(0.001*vega,i);
        fxPortfolio->SetPrice(price,i);

        initTimes[i] = optTime*K_YEAR_LEN;
        initVols[i]  = vol;
    }

    /// Replace the sigma param with new initialisations
    ARM_GP_Vector volLowerBound(nbFx,QVOL_LOWER_BOUND);
    ARM_GP_Vector volUpperBound(nbFx,QVOL_UPPER_BOUND);
    ARM_CurveModelParam* volParam = new ARM_CurveModelParam(ARM_ModelParamType::QVol,&initVols,&initTimes,
        "QVOL","STEPUPRIGHT",&volLowerBound,&volUpperBound);
    if(volParamSize==0)
        FXVolCalibMethod->GetCalibParams().push_back(volParam);
    else if(volParamSize==1)
    {
        delete FXVolCalibMethod->GetCalibParam(0);
        (FXVolCalibMethod->GetCalibParams())[0] = volParam;
    }
    else
	{
        ARM_THROW( ERR_INVALID_ARGUMENT, ARM_USERNAME + " : fx options calibrate only volatilities at the moment");
    }
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateFxCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_TARNFXCalculator::CreateFxCalibMethod(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx)
{
    ARM_ModelParamVector emptyCalibParam;
    ARM_CalibMethod* noLinkedMethod=NULL;
    ARM_CalibMethod* noPreviousMethod=NULL;
    bool isCalibShared=true;

    /// Build an empty bootstrap calib method (filled latter for fx volatility calibration)
    ARM_CalibMethod* fxCalib = new ARM_CalibMethod(
									fxOptionPF,
									emptyCalibParam, //calibParam,
                                    ARM_CalibMethodType::Bootstrap1D,
									ARM_MAX_ITER,
                                    ARM_CalibrationTarget::PriceTarget,
                                    noLinkedMethod,
									noPreviousMethod,
                                    isCalibShared,
									modelIdx);

    return fxCalib;
}
*/
ARM_CalibMethod* ARM_TARNFXCalculator::CreateNumFxCalibMethod()
{
	/// This function is used only when pricing with 1IRFX	
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	string fxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();

	ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(fxModelKey) );
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));

	ARM_ModelParam& volATM = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility);
	ARM_ModelParam& decVol = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Smile);
	ARM_ModelParam& shift = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Shift);
	ARM_ModelParam& lambda = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter);

	size_t nbReset = itsCpnDateStrip->size();

	ARM_GP_Vector* startDates = itsCpnDateStrip->GetFlowStartDates();
	ARM_GP_Vector* endDates = itsCpnDateStrip->GetFlowEndDates();
	ARM_GP_Vector* fwdStartDates = itsCpnDateStrip->GetFwdRateStartDates();
	ARM_GP_Vector* fwdEndDates = itsCpnDateStrip->GetFwdRateEndDates();
	ARM_GP_Vector* resetDates = itsCpnDateStrip->GetResetDates();
	ARM_GP_Vector* payDates = itsCpnDateStrip->GetPaymentDates();
	ARM_GP_Vector* interestDays = itsCpnDateStrip->GetInterestDays();
	ARM_GP_Vector* interestTerms = itsCpnDateStrip->GetInterestTerms();

	double volATMValue, decVolValue, shiftValue, lambdaValue;

	double fxSpot = forex->GetMarketPrice();
	double fwd, maturity;

	ARM_VanillaSecDensityPtrVector vanillaSecDensities;

	ARM_GP_Vector calibStartDates;
	ARM_GP_Vector calibEndDates;
	ARM_GP_Vector calibFwdStartDates;
	ARM_GP_Vector calibFwdEndDates;
	ARM_GP_Vector calibResetDates;
	ARM_GP_Vector calibPayDates;
	ARM_GP_Vector calibInterestDays;
	ARM_GP_Vector calibInterestTerms;

	size_t i;
	for (i = 0; i < nbReset; ++i)
	{
		maturity = ((*resetDates)[i]-asOfDate);

		if (maturity > 0)
		{
			calibStartDates.push_back((*startDates)[i]);
			calibEndDates.push_back((*endDates)[i]);
			calibFwdStartDates.push_back((*fwdStartDates)[i]);
			calibFwdEndDates.push_back((*fwdEndDates)[i]);
			calibResetDates.push_back((*resetDates)[i]);
			calibPayDates.push_back((*payDates)[i]);
			calibInterestDays.push_back((*interestDays)[i]);
			calibInterestTerms.push_back((*interestTerms)[i]);

			fwd = ycBasisForCurve->DiscountPrice(maturity/K_YEAR_LEN)/ycBasisDomCurve->DiscountPrice(maturity/K_YEAR_LEN)*fxSpot;
			
			volATMValue = volATM.GetValue(maturity);

			ARM_DensityFunctorPtr df;

			if (itsSmileFlag)
			{
				decVolValue = decVol.GetValue(maturity);
				shiftValue = shift.GetValue(maturity);
				lambdaValue = lambda.GetValue(maturity);

				df = ARM_DensityFunctorPtr(new ARM_MixtureDensityFunctor(
				fwd,
				maturity/K_YEAR_LEN,
				volATMValue,
				decVolValue,
				shiftValue,
				lambdaValue));
			}
			else
			{
				df = ARM_DensityFunctorPtr(new ARM_ShiftedLNDensityFunctor(
				volATMValue,
				0.0));
			}

			ARM_VanillaSecDensityPtr ds(new ARM_VanillaSecurityDensityFX(
				(*resetDates)[i],
				df, 
				CreateClonedPtr(ycBasisDomCurve),
				CreateClonedPtr(ycBasisForCurve),
				fxSpot));

			vanillaSecDensities.push_back(ds);
		}
	}

	ARM_DateStripPtr calibDateStrip(new ARM_DateStrip(
					&calibStartDates,
					&calibEndDates,
					&calibFwdStartDates,
					&calibFwdEndDates,
					&calibResetDates,
					&calibPayDates,
					&calibInterestDays,
					&calibInterestTerms));

	ARM_StdPortfolioPtr portfolioPtr(NULL);

	string MethodTypeStr = "Numerical";
	ARM_MethodType methodType = (ARM_MethodType) ARM_ArgConv_CalibMethod.GetNumber(MethodTypeStr);

	string MktTargetStr = "UNKNOWN_TAR";
	ARM_MktTargetType  mktTargetType = (ARM_MktTargetType)ARM_ArgConv_TargetFuncMethod.GetNumber(MktTargetStr);

	ARM_CalibMethod* pMethod = new ARM_CalibMethod( portfolioPtr,
													ARM_ModelParamVector(0),
													methodType,
													100,
													mktTargetType,
													NULL, 
													NULL, 
													false, 
													0,
													1,
													true, 
													calibDateStrip,
													vanillaSecDensities );

	return pMethod;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_TARNFXCalculator::CreateAndSetModel()
{
	ARM_PricingModelPtr numModel;

	if (itsModelType == Model2IRFX)
	{
		numModel = ARM_PricingModelPtr(CreateModel2IRFX());
	}
	else if (itsModelType == Model1IRFX)
	{
		numModel = ARM_PricingModelPtr(CreateModel1IRFX());
	}
	else if (itsModelType == ModelNP1IRNFX)
	{
		numModel = ARM_PricingModelPtr(CreateModelNP1IRNFX());
	}

	ARM_RandomGeneratorPtr  randGen(ARM_RandGenFactory.Instance()->CreateSimpleRandGen(
		itsRandGenType2,
		itsRandGenType1,
		itsRandGenAlgo2,
		itsRandGenAlgo1,
		ARM_Transposer::PathOrder,
		itsFirstNbTimes,
		itsFirstNbDims,
		true));

	ARM_TimeStepPerYearScheduler scheduler(0);
	ARM_MeanRevertingSamplerND sampler(&scheduler);
	ARM_PathSchemePtr pathScheme(ARM_PathSchemeFactory.Instance()->CreatePathScheme(ARM_PathScheme::BrownianBridge));
	ARM_MCMethod* mcMethod = new ARM_MCMethod(itsNbSimul,randGen,&sampler,itsBucketSize,ARM_ImpSamplerPtr(NULL),pathScheme);

	numModel->SetNumMethod(ARM_NumMethodPtr( mcMethod ) );
	ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( (itsModelType==ModelNP1IRNFX)||(itsModelType==Model2IRFX)||(itsCorrelType==ARM_ModelParamsSmiled::Fwd)?ARM_Numeraire::RollingCash:ARM_Numeraire::TerminalZc ) );
    numModel->SetNumeraire(numeraire);
	
	SetPricingModel(numModel);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateCalibration
///	Returns: ARM_CalibMethod*
///	Action : create the calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_TARNFXCalculator::CreateCalibration(ModelType modelType)
{
	ARM_CalibMethod* domCalib = NULL;
	vector<ARM_CalibMethod*> forCalib;
	vector<ARM_CalibMethod*> fxCalib;

	int modelIdx;

	/// Diagonal domestic swaptions
	/// --------------------------
	/// Create standard diagonal swaption portfolios (domestic & foreign sigmas calibration)
	if (!itsOneFactorFlag)
	{
		itsDomSwaptionPF = CreateDiagonalSwaption(itsDomesticCcy);

		if (modelType == Model2IRFX)
			modelIdx = ARM_2IRFXModel::DomBasisModel;
		else if (modelType == Model1IRFX)
			modelIdx = ARM_1IRFXModel::DomBasisModel;
		else
			modelIdx = 2*itsNbFX + 1;

		/// Build an empty calib method for latter domestic volatility bootstrapping
		domCalib = CreateIRCalibMethod(itsDomSwaptionPF, modelIdx);
		ComputeSwaptionPrice(itsDomesticCcy,domCalib);
	}

	/// Calibrate swaptions on every foreign currency
	/// ---------------------------------------------
	/// Diagonal foreign swaptions
	/// --------------------------
	int i;
	for (i = 0; i < itsNbFX; i++)
	{
		if ( ((modelType == Model2IRFX) || (modelType == ModelNP1IRNFX)) && !itsOneFactorFlag)
		{
			/// Create standard diagonal swaption portfolios (domestic & foreign sigmas calibration)
			itsForSwaptionPF.push_back(CreateDiagonalSwaption(itsForeignCcy[i]));

			if (modelType == Model2IRFX)
				modelIdx = ARM_2IRFXModel::ForBasisModel;
			else
				modelIdx = i+1+2*itsNbFX+1;

			/// Build an empty calib method for latter domestic volatility bootstrapping
			forCalib.push_back(CreateIRCalibMethod(itsForSwaptionPF[i], modelIdx));
			ComputeSwaptionPrice(itsForeignCcy[i],forCalib[i]);

			if (i == 0)
				domCalib->SetNextMethod(forCalib[i]);
			else
				forCalib[i-1]->SetNextMethod(forCalib[i]);
		}
	}

	/// ATM Forex options
	/// -----------------
	for (i = 0; i < itsNbFX; i++)
	{
		/// Create a Fx option portfolio and compute Fx option target prices
		if ((modelType == Model2IRFX) || (modelType == ModelNP1IRNFX))
		{
			itsFxOptionPF.push_back(CreateFxOption(itsForeignCcy[i]));

			if (modelType == Model2IRFX)
				modelIdx = ARM_2IRFXModel::FxModel;
			else
				modelIdx = i+1+itsNbFX;

			fxCalib.push_back(CreateFxCalibMethod(itsFxOptionPF[i], modelIdx));
			ComputeFxOptionPrices(fxCalib[i], i);
	
			/// Set next methods to handle calibration with a single method
			if (!itsOneFactorFlag)
			{
				if (i == 0)
					forCalib[itsNbFX-1]->SetNextMethod(fxCalib[i]);
				else
					fxCalib[i-1]->SetNextMethod(fxCalib[i]);
			}
		}
	}

	if (modelType == Model1IRFX)
		itsSmiledFxCalibMethod = ARM_CalibMethodPtr(CreateNumFxCalibMethod());

	if (!itsOneFactorFlag)
		return domCalib;
	else
		return fxCalib[0];
}

vector<ARM_DensityFunctor*> ARM_TARNFXCalculator::CreateDensityFunctor()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int nbReset = itsCpnDateStrip->size();

	int start_i = (itsCpnTiming == K_ADVANCE) ? 0 : 1;
	vector<ARM_DensityFunctor*> densitiesVector(itsNbFX*(nbReset-start_i));

	size_t k;
	for (k=0; k<itsNbFX; k++)
	{
		string fxModelKey = FXMODEL_KEY_NAME + itsForeignCcy[k].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
		string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
		string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[k].GetCcyName();
		string ForexKey = FOREX_KEY_NAME + itsForeignCcy[k].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();

		ARM_MixtureModel_Fx* fxMixModel = dynamic_cast< ARM_MixtureModel_Fx* >( GetMktDataManager()->GetData(fxModelKey) );
		ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
		ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
		ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));

		ARM_ModelParam& volATM = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Volatility);
		ARM_ModelParam& decVol = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Smile);
		ARM_ModelParam& shift = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::Shift);
		ARM_ModelParam& lambda = fxMixModel->GetModelParams()->GetModelParam(ARM_ModelParamType::QParameter);

		ARM_GP_Vector* startDates = itsCpnDateStrip->GetFlowStartDates();
		ARM_GP_Vector* endDates = itsCpnDateStrip->GetFlowEndDates();
		ARM_GP_Vector* fwdStartDates = itsCpnDateStrip->GetFwdRateStartDates();
		ARM_GP_Vector* fwdEndDates = itsCpnDateStrip->GetFwdRateEndDates();
		ARM_GP_Vector* resetDates = itsCpnDateStrip->GetResetDates();
		ARM_GP_Vector* payDates = itsCpnDateStrip->GetPaymentDates();
		ARM_GP_Vector* interestDays = itsCpnDateStrip->GetInterestDays();
		ARM_GP_Vector* interestTerms = itsCpnDateStrip->GetInterestTerms();

		double volATMValue, decVolValue, shiftValue, lambdaValue;

		double fxSpot = forex->GetMarketPrice();
		double fwd, maturity;

		ARM_GP_Vector calibStartDates;
		ARM_GP_Vector calibEndDates;
		ARM_GP_Vector calibFwdStartDates;
		ARM_GP_Vector calibFwdEndDates;
		ARM_GP_Vector calibResetDates;
		ARM_GP_Vector calibPayDates;
		ARM_GP_Vector calibInterestDays;
		ARM_GP_Vector calibInterestTerms;

		for (int i=start_i; i<nbReset; ++i)
		{
			maturity = ((*resetDates)[i]-asOfDate);

			if (maturity > 0)
			{
				calibStartDates.push_back((*startDates)[i]);
				calibEndDates.push_back((*endDates)[i]);
				calibFwdStartDates.push_back((*fwdStartDates)[i]);
				calibFwdEndDates.push_back((*fwdEndDates)[i]);
				calibResetDates.push_back((*resetDates)[i]);
				calibPayDates.push_back((*payDates)[i]);
				calibInterestDays.push_back((*interestDays)[i]);
				calibInterestTerms.push_back((*interestTerms)[i]);

				fwd = ycBasisForCurve->DiscountPrice(maturity/K_YEAR_LEN)/ycBasisDomCurve->DiscountPrice(maturity/K_YEAR_LEN)*fxSpot;

				volATMValue = volATM.GetValue(maturity);

				if (itsSmileFlag)
				{
					decVolValue = decVol.GetValue(maturity);
					shiftValue = shift.GetValue(maturity);
					lambdaValue = lambda.GetValue(maturity);

					densitiesVector[(i-start_i)*itsNbFX+k] = new ARM_MixtureDensityFunctor(
														fwd,
														maturity/K_YEAR_LEN,
														volATMValue,
														decVolValue,
														shiftValue,
														lambdaValue);
				}
				else
				{
					densitiesVector[(i-start_i)*itsNbFX+k] = new ARM_ShiftedLNDensityFunctor(volATMValue, 0.0);
				}
			}
		}
	}

	return densitiesVector;
}

void ARM_TARNFXCalculator::Calibrate()
{
	ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());

	ARM_PricingModelPtr CalibModel2IRFXPtr;
	ARM_2IRFXModelPtr CalibModel2IRFX;
	ARM_CalibMethod* CalibMethod2IRFX;

	if (itsModelType == Model1IRFX)
	{
		ARM_1IRFXModel* Model1IRFX = dynamic_cast<ARM_1IRFXModel*>(hybridModel);
		CalibModel2IRFX = Model1IRFX->Get2IRFXModel();
		CalibMethod2IRFX = CreateCalibration(Model2IRFX);
		CalibMethod2IRFX->Calibrate(&*CalibModel2IRFX);

		ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*(*hybridModel->GetModelMap())[ARM_1IRFXModel::FxModel]->Model() );
		fxModel->SetModel2IRFX(&*CalibModel2IRFX);
		itsSmiledFxCalibMethod->Calibrate(fxModel);
		
		Model1IRFX->CalibrateCorrelation(CalibModel2IRFX);
	}

	if (itsModelType == Model2IRFX)
	{
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(false);
	}

    /// Calibrate stochatic models and check errors
	if (!itsOneFactorFlag)
		GetCalibMethod()->Calibrate(hybridModel);

	if (itsModelType == ModelNP1IRNFX)
	{
		ARM_NP1IRNFXModel* ModelNP1IRNFX = dynamic_cast<ARM_NP1IRNFXModel*>(hybridModel);

		vector<ARM_DensityFunctor*> densities = CreateDensityFunctor();

		int i = (itsCpnTiming == K_ADVANCE) ? 0 : 1;
		ARM_GP_Vector newResetTimes((&*itsCpnDateStrip)->GetResetDates()->begin()+i, (&*itsCpnDateStrip)->GetResetDates()->end());
		ARM_GP_VectorPtr ResetTimes( (ARM_GP_Vector*)newResetTimes.Clone() );
		
		for (i=0; i<ResetTimes->size(); i++)
		{
			ARM_Date armdate ((*ResetTimes)[i]);
			(*ResetTimes)[i] = ModelNP1IRNFX->GetTimeFromDate(armdate);
		}

		ModelNP1IRNFX->CalibrateFunctional( *ResetTimes, 
											densities, 
											densities.size()/itsNbFX,
											itsNbFX,
											itsSpaceStepNb, 
											itsStdDevNb,
											false ); // rescaling
	}

	if (itsModelType == Model2IRFX)
	{
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(true);
	}
	
	if (itsModelType == Model1IRFX)
	{
		delete CalibMethod2IRFX;
	}

	ComputeDomesticBasis();
}

double ARM_TARNFXCalculator::Price()
{	/// Calibrate
    CalibrateAndTimeIt();
	
	ARM_PricingModel* hybridModel = &*GetPricingModel();

	double firstColumnPrice = 0.0;
	if (itsColumnsToPrice.size() > 0)
	{
		ARM_GenPricer* genPricer = new ARM_GenPricer( &*GetGenSecurity(),hybridModel);

		ARM_AutoCleaner<ARM_GenPricer> HoldGP(genPricer );

		firstColumnPrice = genPricer->Price();

		string columnName;
		for(size_t i=0;i<itsColumnsToPrice.size();++i)
		{
			columnName = itsColumnsToPrice[i];
			GetPricingData()[columnName] = genPricer->GetPricerInfo()->GetContents(columnName).GetData("Price").GetDouble();
			GetPricingData()[columnName+"StdDev"] = genPricer->GetPricerInfo()->GetContents(columnName).GetData("StdDev").GetDouble();
		}

		if ((itsProductsToPrice[ProbaValue]) && itsIntermediatePrices) 
		{
			GetPricingData()["Probas"] = genPricer->GetPricerInfo()->GetContents( TARNFXColNamesTable[ ProductToPriceColumns[ProbaValue] ]).GetData("IntermediatePrices").GetVector();
		}
	}
	itsHasBeenComputed = true;

    return firstColumnPrice;
}

ARM_Vector* ARM_TARNFXCalculator::ComputeAll()
{	return	NULL;}

void ARM_TARNFXCalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_TARNFXCalculator*>(this)->PriceAndTimeIt();
}

////////////////////////////////////////////////////
///	Class   : ARM_TARNFXCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_TARNFXCalculator::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream	os;
	
	if (!itsDomSwaptionPF.IsNull())
	{
		os << indent << "Domestic Swaption Portfolio" << endl;
		os << indent << itsDomSwaptionPF->toString() << endl;
	}

	size_t i;
	for (i=0; i<itsNbFX; i++)
	{
		char num[3];
		sprintf(num, "%d", i);

		if (!itsForSwaptionPF[i].IsNull())
		{
			os << indent  << "Foreign Swaption Portfolio " << string(num) << endl;
			os << indent  << itsForSwaptionPF[i]->toString() << endl;
		}

		if (!itsFxOptionPF[i].IsNull())
		{
			os << indent  << "Fx Option Portfolio" << string(num) << endl;
			os << indent  << itsFxOptionPF[i]->toString() << endl;
		}
	}

	// A hook to prevent crash in the view
	if (itsModelType == Model2IRFX)
	{
		ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(false);
	}

    os <<  ARM_GenCalculator::toString(indent,nextIndent);

	// A hook to prevent crash in the view
	if (itsModelType == Model2IRFX)
	{
		ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());
		ARM_QModel1F_Fx* fxModel = dynamic_cast<ARM_QModel1F_Fx*>(&*(*hybridModel->GetModelMap())[ARM_2IRFXModel::FxModel]->Model() );
		fxModel->SetIntegratedVersion(true);
	}

	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

