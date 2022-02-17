/*!
 *
 * Copyright (c) IXIS CIB August 2006 Paris
 *
 *	\file hybridirfxcalculator.cpp
 *
 *  \brief file for a generic IR FX Calculator (root for TARN FX, PRDC...)
 *	\author  P. Lam
 *	\version 1.0
 *	\date May 2007
 */

/// this header comes firts as it includes some preprocessor constants!
#include "gpbase/removeidentifiedwarning.h"
#include "gpcalculators/hybridirfxcalculator.h"
#include "gpcalculators/basisconverter.h"
#include "gpcalculators/tarnfxcalculator.h"

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
const unsigned int NB_HybridIRFX_SCHED = 1;

class ARM_TARNFXCalculator;

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Copy constructor
///	Returns: a built object
///	Action : builds the object
/////////////////////////////////////////////////////////////////
ARM_HybridIRFXCalculator::ARM_HybridIRFXCalculator( const ARM_HybridIRFXCalculator& rhs )
:	ARM_GenCalculator( rhs ),
	itsStartDate(rhs.itsStartDate),
	itsEndDate(rhs.itsEndDate),
	itsRedemptionResetDate(rhs.itsRedemptionResetDate),
	itsNbFX(rhs.itsNbFX),
	itsDomesticCcy(rhs.itsDomesticCcy),
	itsForeignCcy(rhs.itsForeignCcy),
	itsFundingCcy(rhs.itsFundingCcy),
	itsPayRec(rhs.itsPayRec),
	itsCpnDayCount(rhs.itsCpnDayCount),
	itsCpnFreq(rhs.itsCpnFreq),
	itsCpnResetGap(rhs.itsCpnResetGap),
	itsStubRule(rhs.itsStubRule),
	itsCpnResetCal(rhs.itsCpnResetCal),
	itsCpnPayCal(rhs.itsCpnPayCal),
	itsCpnTiming(rhs.itsCpnTiming),
	itsCpnNominalCv(rhs.itsCpnNominalCv),
	itsFundNominalCv(rhs.itsFundNominalCv),
	itsFundSpreadCv(rhs.itsFundSpreadCv),
	itsDomesticCpnCv(rhs.itsDomesticCpnCv),
	itsForeignCpnCv(rhs.itsForeignCpnCv),
	itsMinCpnCv(rhs.itsMinCpnCv),
	itsMaxCpnCv(rhs.itsMaxCpnCv),
	itsInitialFXCv(rhs.itsInitialFXCv),	
	itsFundFreq(rhs.itsFundFreq),
	itsFundDayCount(rhs.itsFundDayCount),
	itsRedemptionType(rhs.itsRedemptionType),
	itsRedemptionGap(rhs.itsRedemptionGap),
	itsRedemptionStrike(rhs.itsRedemptionStrike),
	itsFees(rhs.itsFees),
	itsNbSimul(rhs.itsNbSimul),
	itsBucketSize(rhs.itsBucketSize),
	itsRandGenType1(rhs.itsRandGenType1),
	itsRandGenAlgo1(rhs.itsRandGenAlgo1),
	itsRandGenType2(rhs.itsRandGenType2),
	itsRandGenAlgo2(rhs.itsRandGenAlgo2),
	itsFirstNbDims(rhs.itsFirstNbDims),
	itsFirstNbTimes(rhs.itsFirstNbTimes),
	itsFactorsNb(rhs.itsFactorsNb),
	itsTimeStepNb(rhs.itsTimeStepNb),
	itsSpaceStepNb(rhs.itsSpaceStepNb),
	itsStdDevNb(rhs.itsStdDevNb),
	itsSkipPDE(rhs.itsSkipPDE),
	itsRescalling(rhs.itsRescalling),
	itsModelType(rhs.itsModelType),
	itsSmileFlag(rhs.itsSmileFlag),
	itsMixCalib(rhs.itsMixCalib),
	itsOneFactorFlag(rhs.itsOneFactorFlag),
	itsCorrelType(rhs.itsCorrelType),
	itsFirstEventIdx(rhs.itsFirstEventIdx),
	itsCpnDateStrip(rhs.itsCpnDateStrip),
	itsFundDateStrip(rhs.itsFundDateStrip),
	itsColumnsToPrice(rhs.itsColumnsToPrice),
	itsHasBeenComputed(rhs.itsHasBeenComputed),
	itsSmiledFxCalibMethod(rhs.itsSmiledFxCalibMethod),
	itsDomSwaptionPF(rhs.itsDomSwaptionPF),
	itsForSwaptionPF(rhs.itsForSwaptionPF),
	itsFxOptionPF(rhs.itsFxOptionPF)
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Destructor
///	Returns: void
///	Action : destroys the object
/////////////////////////////////////////////////////////////////
ARM_HybridIRFXCalculator::~ARM_HybridIRFXCalculator()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
/////////////////////////////////////////////////////////////////
ARM_HybridIRFXCalculator::ARM_HybridIRFXCalculator(
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
							 int redemptionType,
							 int redemptionGap,
							 double redemptionStrike,
							 const ARM_Curve* fees)
:	
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsDomesticCcy(DomCcy),
	itsFundingCcy(),
	itsPayRec(payRec),
	itsCpnDayCount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsCpnResetGap(cpnResetGap),
	itsStubRule(stubRule),
	itsCpnTiming(cpnTiming),
	itsCpnIntRule(cpnIntRule),
	itsFundFreq(fundFreq),
	itsFundDayCount(fundDayCount),
	itsRedemptionType(redemptionType),
	itsRedemptionGap(redemptionGap),
	itsNbSimul(0),
	itsBucketSize(0),
	itsRandGenType1(0),
	itsRandGenAlgo1(0),
	itsRandGenType2(0),
	itsRandGenAlgo2(0),
	itsFirstNbDims(0),
	itsFirstNbTimes(0),
	itsFactorsNb(0),
	itsTimeStepNb(0),
	itsSpaceStepNb(0),
	itsStdDevNb(0),
	itsSkipPDE(0),
	itsRescalling(0),
	itsModelType(0),
	itsSmileFlag(0),
	itsMixCalib(0),
	itsOneFactorFlag(0),
	itsFirstEventIdx(0),
	itsCpnDateStrip(0),
	itsFundDateStrip(0),
	itsColumnsToPrice(0),
	itsHasBeenComputed(false),
	itsSmiledFxCalibMethod(0),
	itsDomSwaptionPF(0),
	itsForSwaptionPF(0),
	itsFxOptionPF(0)
{
	itsForeignCcy.push_back(ForCcy);

	if (FundCcy)
		itsFundingCcy = *FundCcy;

	if (cpnNominal)
		itsCpnNominalCv = *cpnNominal;

	if (domesticCpn)
		itsDomesticCpnCv.push_back(*domesticCpn);

	if (foreignCpn)
		itsForeignCpnCv.push_back(*foreignCpn);

	if (InitialFX)
		itsInitialFXCv.push_back(*InitialFX);

	if (fundNominal)
		itsFundNominalCv = *fundNominal;

	if (fundSpread)
		itsFundSpreadCv = *fundSpread;

	if (MinCpn)
		itsMinCpnCv = *MinCpn;

	if (MaxCpn)
		itsMaxCpnCv = *MaxCpn;

	if (fees)
		itsFees = *fees;

	itsRedemptionStrike.push_back(redemptionStrike);

	itsNbFX = itsForeignCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	string payModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Constructor
///	Returns: void
///	Action : builds the object (stand alone version)
///          case of several foreign currencies
/////////////////////////////////////////////////////////////////
ARM_HybridIRFXCalculator::ARM_HybridIRFXCalculator(
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
							 int redemptionType,
							 int redemptionGap,
							 ARM_GP_Vector& redemptionStrike,
							 const ARM_Curve& fees)
:	
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsDomesticCcy(DomCcy),
	itsForeignCcy(ForCcy),
	itsFundingCcy(FundCcy),
	itsPayRec(payRec),
	itsCpnDayCount(cpnDayCount),
	itsCpnFreq(cpnFreq),
	itsCpnResetCal(cpnResetCal),
	itsCpnPayCal(cpnPayCal),
	itsCpnResetGap(cpnResetGap),
	itsStubRule(stubRule),
	itsCpnTiming(cpnTiming),
	itsCpnIntRule(cpnIntRule),
	itsCpnNominalCv(cpnNominal),
	itsDomesticCpnCv(domesticCpn),
	itsForeignCpnCv(foreignCpn),
	itsMinCpnCv(MinCpn),
	itsMaxCpnCv(MaxCpn),
	itsInitialFXCv(InitialFX),
	itsFundFreq(fundFreq),
	itsFundDayCount(fundDayCount),
	itsFundNominalCv(fundNominal),
	itsFundSpreadCv(fundSpread),
	itsRedemptionType(redemptionType),
	itsRedemptionGap(redemptionGap),
	itsRedemptionStrike(redemptionStrike),
	itsFees(fees),
	itsNbSimul(0),
	itsBucketSize(0),
	itsRandGenType1(0),
	itsRandGenAlgo1(0),
	itsRandGenType2(0),
	itsRandGenAlgo2(0),
	itsFirstNbDims(0),
	itsFirstNbTimes(0),
	itsFactorsNb(0),
	itsTimeStepNb(0),
	itsSpaceStepNb(0),
	itsStdDevNb(0),
	itsSkipPDE(0),
	itsRescalling(0),
	itsModelType(0),
	itsSmileFlag(0),
	itsMixCalib(0),
	itsOneFactorFlag(0),
	itsFirstEventIdx(0),
	itsCpnDateStrip(0),
	itsFundDateStrip(0),
	itsColumnsToPrice(0),
	itsHasBeenComputed(false),
	itsSmiledFxCalibMethod(0),
	itsDomSwaptionPF(0),
	itsForSwaptionPF(0),
	itsFxOptionPF(0)
{
	itsNbFX = ForCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	string payModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());
}

//////////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Init
///	Returns: void
///	Action : initialise calculator with market data and model parameters
//////////////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::Init(
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

    /// Create a 2IR+FX model
//    CreateAndSetModelAndTimeIt();

	Update(&mktDataManager);
}

//////////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Init
///	Returns: void
///	Action : initialise calculator with market data and model parameters
//////////////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::Init(
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: ComputeDomesticBasis
///	Returns: void
///	Action : Convert le basis from funding Ccy to domestic Ccy
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::ComputeDomesticBasis()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcFundKey = YC_KEY_NAME + itsFundingCcy.GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisFundKey = YC_BASIS_KEY_NAME + itsFundingCcy.GetCcyName();
	string FundForexKey = FOREX_KEY_NAME + itsFundingCcy.GetCcyName() + "/" + itsDomesticCcy.GetCcyName();

	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycfundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcFundKey));
	ARM_ZeroCurve* basisFundCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisFundKey));
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(FundForexKey));

	size_t size = itsFundDateStrip->size();

	ARM_GP_Vector* resetDates = itsFundDateStrip->GetResetDates();

	ARM_GP_Vector nominal(size);
	ARM_GP_Vector fundNominal(size);
	ARM_GP_Vector fundSpread(size);

	size_t i;

	for (i = 0; i < size; ++i)
	{
		nominal[i] = itsCpnNominalCv.Interpolate((*resetDates)[i]-asOfDate);
		fundNominal[i] = itsFundNominalCv.Interpolate((*resetDates)[i]-asOfDate);
		fundSpread[i] = itsFundSpreadCv.Interpolate((*resetDates)[i]-asOfDate);
	}
	ARM_BasisConverter basisConveter(itsDomesticCcy,
					itsFundingCcy,
					itsFundDateStrip,
					itsFundDateStrip,
					itsFundDateStrip,
					ycDomCurve,	
					ycfundCurve,
					ycBasisDomCurve,
					basisFundCurve,
					*forex,
					nominal,
					fundNominal,
					fundSpread);

	ARM_GP_VectorPtr basisMargin = ARM_GP_VectorPtr((ARM_GP_Vector*)(basisConveter.ComputeDomMargin()).Clone());
	GetGenSecurity()->GetCstManager()->insert("FundMargin",ARM_GramFctorArg(basisMargin));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_HybridIRFXCalculator::ColumnNames() const
{
    vector<string>				colNamesVec(0);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(0); 

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: ProductsToPriceColumnNames
///	Returns: ARM_StringVector
///	Action : Create Products to price column names
/////////////////////////////////////////////////////////////////
ARM_StringVector ARM_HybridIRFXCalculator::ProductsToPriceColumnNames()
{
	return 0;
}

//////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the TARN FX.
//////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_HybridIRFXCalculator::DatesStructure() const
{
/*    const char*	resetCalendar = itsCpnResetCal.c_str();
    const char*	payCalendar   = itsCpnPayCal.c_str();

    int fwdRule		= K_MOD_FOLLOWING;	// for forward dates
    int intRule		= itsCpnIntRule;	// for interest dates
    int stubRule	= itsStubRule;
	int resetFreq   = itsCpnFreq;
    int resetTiming = itsCpnTiming;
	int payFreq     = itsCpnFreq;
	int payGap      = GETDEFAULTVALUE;
    int payTiming   = K_ARREARS;
	int cpnResetGap	= -itsCpnResetGap;

	int spotDays = itsDomesticCcy.GetSpotDays();
	
	if (itsCpnDateStrip.IsNull())
	{
		itsCpnDateStrip = ARM_DateStripPtr(new ARM_DateStrip(
								itsStartDate, 
								itsEndDate, 
								resetFreq, 
								itsCpnDayCount,
								resetCalendar, 
								fwdRule, 
								intRule, 
								stubRule,
								cpnResetGap, 
								payFreq, 
								payGap, 
								payCalendar,
								resetTiming, 
								payTiming));

		if (&itsFundingCcy == NULL)
		{
			char* fundCcyName = itsFundingCcy.GetCcyName();

			// For the in arrears case we create a fake start date to get the first 
			// reset date of the funding leg
			if (itsCpnTiming == K_ARREARS)
			{
				ARM_Date	firstFundingResetDate(itsStartDate);
				firstFundingResetDate.PreviousBusinessDay(fabs(spotDays), fundCcyName);

				itsCpnDateStrip->InsertDate(0, 0, 0, 0, 0, firstFundingResetDate.GetJulian(), 0, 0, 0);
			}
		}
	}

	int	firstEventIdx = 0;
	if (&itsFundingCcy != NULL)
	{
		int fundFreq = itsFundFreq;
		int fundDayCount = itsFundDayCount;
		int fundResetTiming = K_ADVANCE;

		if (itsFundDateStrip.IsNull())
		{
			itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(
									itsStartDate, 
									itsEndDate, 
									fundFreq, 
									fundDayCount,
									resetCalendar, 
									fwdRule, 
									intRule, 
									stubRule,
									-spotDays, 
									fundFreq, 
									payGap, 
									payCalendar,
									fundResetTiming, 
									payTiming));
		}
	}

	/// Merge schedules on "ResetDate"
	ARM_DateStripVector SchedVect(NB_HybridIRFX_SCHED, NULL);
	SchedVect[RF_CPNFX_SCHED] = &*itsCpnDateStrip;

	ARM_DateStripCombiner EventSchedule(SchedVect, "ResetDate");

	// Memorize first event index
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	while( (firstEventIdx < eventDates->size()) && ((*eventDates)[firstEventIdx] < asOfDate) )
		++firstEventIdx;

	const_cast<ARM_HybridIRFXCalculator*>(this)->itsFirstEventIdx = firstEventIdx;


    return	EventSchedule;*/

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
    ARM_DateStripVector SchedVect(NB_HybridIRFX_SCHED, NULL);
    SchedVect[RF_CPNFX_SCHED] = &*itsCpnDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect, "ResetDate");

	// Memorize first event index
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	int	firstEventIdx = 0;
	while( (firstEventIdx < eventDates->size()) && ((*eventDates)[firstEventIdx] < asOfDate) )
		++firstEventIdx;

	const_cast<ARM_HybridIRFXCalculator*>(this)->itsFirstEventIdx = firstEventIdx;

    return	EventSchedule;
}

/////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CustomDatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the product
///			customized Schedule Not implemented yet.
/////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_HybridIRFXCalculator::CustomDatesStructure(const ARM_DateStripVector& dateStrips) const
{
	return ARM_DateStripCombiner();
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_HybridIRFXCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	vector<string>				colNamesVec(0);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(0); 

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: Create2IRFXModel
///	Returns: ARM_PricingModelPtr
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
ARM_2IRFXModel* ARM_HybridIRFXCalculator::CreateModel2IRFX()
{
    /// Create the 2IR+FX model

    /// Q vol without any curve because will be bootstrapped latter
    ARM_CurveModelParam volParam( ARM_ModelParamType::QVol,HWVOL_LOWER_BOUND,"QVOL");
    ARM_ModelParamVector modelParams(3);
    modelParams[0]  = &volParam;

    /// Multi-assets avec modèles locaux sur le Fx si nécessaire
    int nbModels = ARM_2IRFXModel::ForBasisModel+1;

	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcForKey = YC_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string MrsForKey = MRS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string QFxKey = QFX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string CorrelMatrixKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")"; 

    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    names[ARM_2IRFXModel::DomModel]        = YcDomKey;
    names[ARM_2IRFXModel::ForModel]        = YcForKey;
    names[ARM_2IRFXModel::FxModel]         = ForexKey;	
    names[ARM_2IRFXModel::DomBasisModel]   = YcBasisDomKey;
    names[ARM_2IRFXModel::ForBasisModel]   = YcBasisForKey;

	depends[ARM_2IRFXModel::FxModel]       = ARM_StringVector(2);
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::DomModel]	   = names[ARM_2IRFXModel::DomBasisModel];
	depends[ARM_2IRFXModel::FxModel][ARM_2IRFXModel::ForModel]	   = names[ARM_2IRFXModel::ForBasisModel];
	depends[ARM_2IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::DomModel]);
	depends[ARM_2IRFXModel::ForBasisModel] = ARM_StringVector(1,names[ARM_2IRFXModel::ForModel]);
	
    vector< ARM_PricingModelPtr > models(nbModels);
    /// Create the Q1F model for domestic IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsDomKey));
    modelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
   
	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
    models[ARM_2IRFXModel::DomModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycDomCurve),ARM_ModelParamsQ1F(modelParams),true) );
	
    /// Create the Q1F model for foreign IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsForKey));
    modelParams[2] =  new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
  
	ARM_ZeroCurve* ycForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey));
    models[ARM_2IRFXModel::ForModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycForCurve),ARM_ModelParamsQ1F(modelParams),true) );

    /// Create the Q1F model for FX market (with a default MRS=0 because it doesn't matter !)
	ARM_CurveModelParam mrsFxParam( ARM_ModelParamType::MeanReversion,0.0,"FXMRS");
    modelParams[1]   = &mrsFxParam;
    if(GetMktDataManager()->TestIfKeyMissing(QFxKey))
    {
        /// Degenerate to pure LN
        modelParams[2] =  new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 1.0,"Q");
    }
    else
	    modelParams[2] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(QFxKey));

	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
    ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(CorrelMatrixKey));

    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));

	ARM_ModelParamsQ1F_Fx modelParamsFx( modelParams, CreateClonedPtr(ycBasisDomCurve), CreateClonedPtr(ycBasisForCurve), forex->GetMarketPrice() );

	ARM_QModel1F_Fx* fxModel = static_cast<ARM_QModel1F_Fx*>( 
									ARM_EqFx_ModelFactory.Instance()->CreateModel( CreateClonedPtr(ycBasisDomCurve), 
									modelParams, forex->GetMarketPrice(), 
									CreateClonedPtr(ycBasisForCurve),
									*correlMatrix,
									ARM_EqFx_ModelFactoryImp::Q1F_Model ));

	fxModel->SetIntegratedVersion(false);

    models[ARM_2IRFXModel::FxModel] = ARM_PricingModelPtr( fxModel );

    /// Create both domestic & foreign forward margin models
    models[ARM_2IRFXModel::DomBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve)) );
    models[ARM_2IRFXModel::ForBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisForCurve)) );

   	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	/// Finally, We create and set the brand new 2IR+FX Ferrari model !
    ARM_2IRFXModel* hybridModel = new ARM_2IRFXModel(modelMap,*correlMatrix);

	return hybridModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateModel1IRFX
///	Returns: ARM_PricingModelPtr
///	Action : create the 1IR+FX model
/////////////////////////////////////////////////////////////////
ARM_1IRFXModel* ARM_HybridIRFXCalculator::CreateModel1IRFX()
{
    /// Create the 1IR+FX model

    /// Q vol without any curve because will be bootstrapped latter
    ARM_CurveModelParam HWVolParam( ARM_ModelParamType::QVol,HWVOL_LOWER_BOUND,"QVOL");
    ARM_ModelParamVector modelParams(3);
    modelParams[0]  = &HWVolParam;

    /// Multi-assets avec modèles locaux sur le Fx si nécessaire
    int nbModels = ARM_1IRFXModel::DomBasisModel+1;

	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string MrsForKey = MRS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string QFxKey = QFX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string CorrelMatrixKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")"; 

    ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
    names[ARM_1IRFXModel::DomModel]        = YcDomKey;
    names[ARM_1IRFXModel::FxModel]         = ForexKey;	
    names[ARM_1IRFXModel::DomBasisModel]   = YcBasisDomKey;

	depends[ARM_1IRFXModel::FxModel]       = ARM_StringVector(1);
	depends[ARM_1IRFXModel::FxModel][ARM_1IRFXModel::DomModel]	   = names[ARM_1IRFXModel::DomBasisModel];
	depends[ARM_1IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_1IRFXModel::DomModel]);
	
    vector< ARM_PricingModelPtr > models(nbModels);
    /// Create the Q1F model for domestic IR market
	modelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsDomKey));
    modelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");
   
	ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
    models[ARM_1IRFXModel::DomModel] = ARM_PricingModelPtr( new ARM_QModel1F( CreateClonedPtr(ycDomCurve),ARM_ModelParamsQ1F(modelParams),true) );
	

	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
    ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));
	ARM_GP_Matrix* correlMatrix = dynamic_cast<ARM_GP_Matrix*>(GetMktDataManager()->GetData(CorrelMatrixKey));

	double FXSpot = forex->GetMarketPrice();

	ARM_CurveModelParam volParam( ARM_ModelParamType::QVol,0.0,"QVol");
	ARM_CurveModelParam humpParam( ARM_ModelParamType::Hump,1.0,"Hump");
	ARM_CurveModelParam correlParam( ARM_ModelParamType::BetaCorrelation,0.0,"Beta");
    modelParams[0]  = &volParam;
	modelParams[1]  = &humpParam;
	modelParams[2]  = &correlParam;

	ARM_ModelParamsSmiled_Fx modelParamsObj(
			modelParams,
			ARM_ZeroCurvePtr(CreateClonedPtr(ycBasisDomCurve)),
			ARM_ZeroCurvePtr(CreateClonedPtr(ycBasisForCurve)),
			FXSpot,
			(itsCorrelType==ARM_ModelParamsSmiled::Fwd?2:itsFactorsNb),
			(ARM_ModelParamsSmiled::CorrelType) itsCorrelType,
			0.0);

	ARM_PricingModelPtr mod( new ARM_SmiledModel_Fx(CreateClonedPtr( ycBasisDomCurve ),&modelParamsObj,itsTimeStepNb,itsSpaceStepNb,itsStdDevNb,(itsSkipPDE?true:false),(itsRescalling?true:false)));

    models[ARM_1IRFXModel::FxModel] = mod;

    /// Create both domestic & foreign forward margin models
    models[ARM_1IRFXModel::DomBasisModel] = ARM_PricingModelPtr( new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve)) );
   	/// Create a modelnamemap
	ARM_ModelNameMap modelMap( names, models, depends );

	ARM_GP_Matrix correl(itsFactorsNb+1,itsFactorsNb+1,0.0);
	for (size_t i = 0; i < itsFactorsNb+1; ++i)
		correl(i,i) = 1.0;

	ARM_2IRFXModelPtr Model2IRFX = ARM_2IRFXModelPtr(CreateModel2IRFX());

	/// Finally, We create and set the brand new 2IR+FX Ferrari model !
    ARM_1IRFXModel* hybridModel = new ARM_1IRFXModel(modelMap,correl,Model2IRFX);

	return hybridModel;
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the 2IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::CreateAndSetModel()
{
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_PRDCalculator
///	Routine: CreateDiagonalSwaption
///	Returns: a portfolio
///	Action : create the list of diagonal swaptions for domestic market
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_HybridIRFXCalculator::CreateDiagonalSwaption(ARM_Currency& ccy)
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: ComputeSwaptionPrice
///	Returns: void
///	Action : Compute the swaption price
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::ComputeSwaptionPrice(ARM_Currency& ccy, ARM_CalibMethod* IRVolCalibMethod)
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
ARM_CalibMethod* ARM_HybridIRFXCalculator::CreateIRCalibMethod(
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateFxOption
///	Returns: a portfolio
///	Action : create the list of ATM FX options
/////////////////////////////////////////////////////////////////
ARM_StdPortfolioPtr ARM_HybridIRFXCalculator::CreateFxOption(ARM_Currency& foreignCcy)
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: ComputeFxOptionPrice
///	Returns: nothing
///	Action : compute market target prices of the FX option
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::ComputeFxOptionPrices(ARM_CalibMethod* FXVolCalibMethod, int Ccyi)
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateFxCalibration
///	Returns: void
///	Action : create the calibration for IR parameters (vol & MRS)
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_HybridIRFXCalculator::CreateFxCalibMethod(const ARM_StdPortfolioPtr fxOptionPF, int modelIdx)
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


/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateAndSetCalibration
///	Returns: void
///	Action : create and set the calib Method
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::CreateAndSetCalibration()
{
	SetCalibMethod(ARM_CalibMethodPtr(CreateCalibration((ModelType)itsModelType)));
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_TARNFXCalculator
///	Routine: CreateCalibration
///	Returns: ARM_CalibMethod*
///	Action : create the calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_HybridIRFXCalculator::CreateCalibration(ModelType modelType)
{
	return NULL;
}

void ARM_HybridIRFXCalculator::UpdateModel()
{
	itsHasBeenComputed = false;
	CreateAndSetModel();
	itsHasBeenComputed = false;
}

void ARM_HybridIRFXCalculator::UpdateCalibration(bool isUpdateStrike)
{
	itsHasBeenComputed = false;
	ComputeDomesticBasis();
	CreateAndSetCalibration();
	itsHasBeenComputed = false;
}

void ARM_HybridIRFXCalculator::Calibrate()
{
	ARM_MultiAssetsModel* hybridModel = dynamic_cast<ARM_MultiAssetsModel*>(&*GetPricingModel());

	ARM_2IRFXModelPtr CalibModel2IRFX;
	ARM_CalibMethod* CalibMethod2IRFX;

	ARM_1IRFXModel* Model1IRFX = dynamic_cast<ARM_1IRFXModel*>(hybridModel);
	CalibModel2IRFX = Model1IRFX->Get2IRFXModel();
	CalibMethod2IRFX = CreateCalibration(Model2IRFX);
	CalibMethod2IRFX->Calibrate(&*CalibModel2IRFX);

	ARM_SmiledModel_Fx* fxModel = dynamic_cast<ARM_SmiledModel_Fx*>(&*(*hybridModel->GetModelMap())[ARM_1IRFXModel::FxModel]->Model() );
	fxModel->SetModel2IRFX(&*CalibModel2IRFX);
	itsSmiledFxCalibMethod->Calibrate(fxModel);
	
	Model1IRFX->CalibrateCorrelation(CalibModel2IRFX);

    /// Calibrate stochatic models and check errors
	if (!itsOneFactorFlag)
		GetCalibMethod()->Calibrate(hybridModel);

	delete CalibMethod2IRFX;

	ComputeDomesticBasis();
}

////////////////////////////////////////////////////
///	Class   : ARM_HybridIRFXCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
string ARM_HybridIRFXCalculator::toString(const string& indent, const string& nextIndent) const
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

    os <<  ARM_GenCalculator::toString(indent,nextIndent);

	return os.str();
}


CC_END_NAMESPACE()

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

