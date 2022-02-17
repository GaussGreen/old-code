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
////#include "gpbase/datestripconvert.h"
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
////#include <crv/volflat.h>
//#include <inst/forex.h>
//#include <inst/swaption.h>
//#include <inst/option.h>
//#include <inst/portfolio.h>
//#include <inst/swap.h>
#include <util/fromto.h>


/// STL
#include <iomanip> /// for setprecision()
#include <list>
CC_USING_NS(std,list)


CC_BEGIN_NAMESPACE( ARM )


/// Default MDM key names
const string YC_KEY_NAME		= "YC_";
const string YC_BASIS_KEY_NAME	= "YC_BASIS_";
const string FOREX_KEY_NAME		= "FOREX_";
const string OSWMODEL_KEY_NAME	= "OSWMOD_";
const string FXMODEL_KEY_NAME	= "FXMOD_";
const string CORREL_KEY_NAME	= "CORREL_";
const string MRS_KEY_NAME		= "MRS_";
const string QFX_KEY_NAME		= "Q_";


const string ARM_HybridIRFXCalculator::ColNamesTable [] =
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
	"FixedFunding",
	"Funding",
	"InitFX",
	"DomCpn",
	"FgnCpn",
	"FxSpot1",
	"FxSpot2",
	"FxSpot3",
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


const int ARM_HybridIRFXCalculator::ProductToPriceColumns[] =
{
	ARM_HybridIRFXCalculator::Funding,
	ARM_HybridIRFXCalculator::Coupon,
	ARM_HybridIRFXCalculator::DF,
	ARM_HybridIRFXCalculator::PaidCoupon,
	ARM_HybridIRFXCalculator::RealCoupon,
	ARM_HybridIRFXCalculator::IsAlive,
	ARM_HybridIRFXCalculator::Cap,
	ARM_HybridIRFXCalculator::Floor,
	ARM_HybridIRFXCalculator::Redemption,
	ARM_HybridIRFXCalculator::RealRedemption,
	ARM_HybridIRFXCalculator::RealFunding,
	ARM_HybridIRFXCalculator::Duration,
	ARM_HybridIRFXCalculator::Proba,
	ARM_HybridIRFXCalculator::TARN,
};


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
	itsReferenceDate(rhs.itsReferenceDate),
	itsNbFX(rhs.itsNbFX),
	itsDomesticCcy(rhs.itsDomesticCcy),
	itsForeignCcy(rhs.itsForeignCcy),
	itsCouponCcy(rhs.itsCouponCcy),
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
	itsNbPastFixings(rhs.itsNbPastFixings),
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
	itsIntermediatePrices(rhs.itsIntermediatePrices),
	itsColumnsToPrice(rhs.itsColumnsToPrice),
	itsHasBeenComputed(rhs.itsHasBeenComputed),
	itsSmiledFxCalibMethod(rhs.itsSmiledFxCalibMethod),
	itsDomSwaptionPF(rhs.itsDomSwaptionPF),
	itsForSwaptionPF(rhs.itsForSwaptionPF),
	itsFxOptionPF(rhs.itsFxOptionPF),
	itsvFundIndex(rhs.itsvFundIndex),
	itsPayOffName(rhs.itsPayOffName),
	itsIsNormal(rhs.itsIsNormal)
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
///          called by the children classes only
/////////////////////////////////////////////////////////////////
ARM_HybridIRFXCalculator::ARM_HybridIRFXCalculator(
							 const ARM_Date& asOfDate,
							 const ARM_Date& startDate,
							 const ARM_Date& endDate,
							 const ARM_Currency& DomCcy,
							 const vector<ARM_Currency>& ForCcy,
							 const ARM_Currency& CpnCcy,
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
							 const ARM_Curve* target,
							 const ARM_Curve* fees,
							 const string& FXChoice,
							 int intermediatePrices,
							 const ARM_StringVector& columnsToPrice,
							 const ARM_FXTARNPayoffType& payOffName,
							 ARM_FixingSched* pastFixings,
							 char* refDate,
							 ARM_Date& effDate)
:	
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsReferenceDate(refDate),
	itsEffectiveDate(effDate),
	itsDomesticCcy(DomCcy),
	itsForeignCcy(ForCcy),
	itsCouponCcy(CpnCcy),
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
	itsColumnsToPrice(columnsToPrice),
	itsHasBeenComputed(false),
	itsSmiledFxCalibMethod(0),
	itsDomSwaptionPF(0),
	itsForSwaptionPF(0),
	itsFxOptionPF(0),
	itsvFundIndex(0),
	itsFXChoice(FXChoice),
	itsPayOffName(payOffName),
	itsIntermediatePrices(intermediatePrices),
	itsNbPastFixings(0)
{
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

	if (target)
		itsTarget = *target;

	itsRedemptionStrike.push_back(redemptionStrike);

	itsNbFX = itsForeignCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

	itsIsNormal = (strcmp(itsCouponCcy.GetCcyName(), itsDomesticCcy.GetCcyName()) == 0);
	if (asOfDate > itsStartDate)
	{
		if (pastFixings)
		{
			SetFixings(pastFixings);
			itsNbPastFixings = pastFixings->GetNbFixings();
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "AsOf > Start : past fixings expected");
	}

	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	itsPayModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());
	if (((itsPayOffName == ARM_TARNFXPayoffType::INDIANFX)
		 ||
		 (itsPayOffName == ARM_TARNFXPayoffType::SWITCHER))
		&&
		(!itsIsNormal))
	{
		itsPayModelName = YC_BASIS_KEY_NAME + string(itsForeignCcy[0].GetCcyName());
	}
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
							 std::vector<double>& redemptionStrike,
							 const ARM_Curve& target,
							 const ARM_Curve& fees,
							 const string& FXChoice,
							 int intermediatePrices,
							 const ARM_StringVector& columnsToPrice,
							 const ARM_FXTARNPayoffType& payOffName,
							 ARM_FixingSched* pastFixings,
							 char* refDate,
							 ARM_Date& effDate)
:	
	ARM_GenCalculator(asOfDate),
	itsStartDate(startDate),
	itsEndDate(endDate),
	itsReferenceDate(refDate),
	itsEffectiveDate(effDate),
	itsDomesticCcy(DomCcy),
	itsForeignCcy(ForCcy),
	itsFundingCcy(FundCcy),
	itsCouponCcy(DomCcy),
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
	itsTarget(target),
	itsFees(fees),
	itsNbPastFixings(0),
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
	itsColumnsToPrice(columnsToPrice),
	itsHasBeenComputed(false),
	itsSmiledFxCalibMethod(0),
	itsDomSwaptionPF(0),
	itsForSwaptionPF(0),
	itsFxOptionPF(0),
	itsvFundIndex(0),
	itsFXChoice(FXChoice),
	itsPayOffName(payOffName),
	itsIntermediatePrices(intermediatePrices)
{
	if (asOfDate > itsStartDate)
	{
		if (pastFixings)
		{
			SetFixings(pastFixings);
			itsNbPastFixings = pastFixings->GetNbFixings();
		}
		else
			ARM_THROW( ERR_INVALID_ARGUMENT, "AsOf > Start : past fixings expected");
	}

	itsNbFX = ForCcy.size();
	itsRedemptionResetDate = itsEndDate;
	itsRedemptionResetDate.PreviousBusinessDay( fabs(itsRedemptionGap), const_cast<char*>(itsCpnResetCal.c_str()) );

	itsIsNormal = (strcmp(itsCouponCcy.GetCcyName(), itsDomesticCcy.GetCcyName()) == 0);

	DatesStructure();

	/// Check input datas
    CheckDataAndTimeIt();

	itsPayModelName = YC_BASIS_KEY_NAME + string(itsDomesticCcy.GetCcyName());
	if (((itsPayOffName == ARM_TARNFXPayoffType::INDIANFX)
		 ||
		 (itsPayOffName == ARM_TARNFXPayoffType::SWITCHER))
		&&
		(!itsIsNormal))
	{
		itsPayModelName = YC_BASIS_KEY_NAME + string(itsForeignCcy[0].GetCcyName());
	}
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
	
	// chooser : model must be NP1IRNFX
	if (itsNbFX > 1)
		modelType = ModelNP1IRNFX;

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

	// chooser : model must be NP1IRNFX
	if (itsNbFX > 1)
		modelType = ModelNP1IRNFX;

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
	/// if funding ccy is different from dom and every fgn,
	/// add YC, BASIS YC & FOREX for this funding ccy
	if (strcmp(itsFundingCcy.GetCcyName(), itsDomesticCcy.GetCcyName()) != 0)
	{
		bool found = false;
		for (i=0; i<itsNbFX; i++)
		{
			if (strcmp(itsFundingCcy.GetCcyName(), itsForeignCcy[i].GetCcyName()) == 0)
				found = true;
		}
		if (!found)
		{
			string ccyName(itsFundingCcy.GetCcyName());
		
			keys.push_back(YC_KEY_NAME + ccyName);
			marketDataManager->RegisterData(YC_KEY_NAME + ccyName, zeroCurves[itsNbFX+1]);

			keys.push_back(YC_BASIS_KEY_NAME + ccyName);
			marketDataManager->RegisterData(YC_BASIS_KEY_NAME + ccyName, basisCurves[itsNbFX+1]);

			keys.push_back(FOREX_KEY_NAME + ccyName + "/" + domCcyName);
			marketDataManager->RegisterData(FOREX_KEY_NAME + ccyName + "/" + domCcyName, forex[itsNbFX+1]);
		}
	}

	string correlKey;
	if (itsModelType == ModelNP1IRNFX)
		correlKey = string("CORREL_DOMCCY_FORCCY");
	else
		correlKey = CORREL_KEY_NAME +"(" + itsDomesticCcy.GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "," + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName() + ")"; 

	keys.push_back(correlKey);
	marketDataManager->RegisterData(correlKey, correlMatrix);

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

	std::vector<double>& resetDates = itsFundDateStrip->GetResetDates();

	std::vector<double> nominal(size);
	std::vector<double> fundNominal(size);
	std::vector<double> fundSpread(size);

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

	ARM_GP_VectorPtr basisMargin = ARM_GP_VectorPtr(static_cast<ARM_GP_Vector*>(basisConveter.ComputeDomMargin()).Clone());
	GetGenSecurity()->GetCstManager()->insert("FundMargin",ARM_GramFctorArg(basisMargin));
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: ColumnNames
///	Returns: ARM_RowInfo
///	Action : create the column names of the deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_HybridIRFXCalculator::ColumnNames() const
{
    size_t	colNamesSize = sizeof(ColNamesTable)/sizeof(ColNamesTable[0]);
    vector<string>				colNamesVec(colNamesSize);
    vector<ARM_GP_VALUE_TYPE>	colTypeVec(colNamesSize, ARM_STRING); 

    for(size_t i=0; i<colNamesSize; ++i)
        colNamesVec[i] = ColNamesTable[i];

    ARM_RowInfo	rowInfo(colNamesVec, colTypeVec);

    return rowInfo;
}

//////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: DatesStructure
///	Returns: ARM_DateStripCombiner
///	Action : create the list of all event/notice dates of the TARN FX.
//////////////////////////////////////////////////////////////////////
ARM_DateStripCombiner ARM_HybridIRFXCalculator::DatesStructure() const
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
	int nbPastNoCall, size;
	std::vector<double>& payDates = NULL;

	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	if (itsFundDateStrip.IsNull())
	{
		itsFundDateStrip = ARM_DateStripPtr(new ARM_DateStrip(itsStartDate, itsEndDate,itsFundFreq, itsFundDayCount,resetCalendar,fwdRule, intRule, 
								stubRule,-fabs(spotDays), itsFundFreq,payGap, payCalendar,fundResetTiming, payTiming, true, itsReferenceDate));
		// past  no call
		nbPastNoCall=0;
	    payDates = itsFundDateStrip->GetPaymentDates();
		size = itsFundDateStrip->size();
		while(nbPastNoCall < size && (*payDates)[nbPastNoCall] < asOfDate - K_NEW_DOUBLE_TOL)
			++nbPastNoCall;

		/// Keep only the futur notification dates
		itsFundDateStrip->ResizeAndBuilt(nbPastNoCall,size);
	}
	
	if (itsCpnDateStrip.IsNull())
	{
		itsCpnDateStrip = ARM_DateStripPtr(new ARM_DateStrip(itsStartDate, itsEndDate, resetFreq, itsCpnDayCount,resetCalendar, fwdRule, intRule, 
								stubRule,-fabs(itsCpnResetGap), payFreq, payGap, payCalendar, resetTiming, payTiming, true, itsReferenceDate));
		
		// past  no call
		nbPastNoCall=0;
	    payDates = itsCpnDateStrip->GetPaymentDates() ;
		size = itsCpnDateStrip->size();
		while(nbPastNoCall < size && (*payDates)[nbPastNoCall] < asOfDate - K_NEW_DOUBLE_TOL)
			++nbPastNoCall;

		/// Keep only the futur notification dates
		itsCpnDateStrip->ResizeAndBuilt(nbPastNoCall,size);

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
    SchedVect[CPNFX_SCHED] = &*itsCpnDateStrip;

    ARM_DateStripCombiner EventSchedule(SchedVect, "ResetDate");

	// Memorize first event index
	ARM_VectorPtr eventDates = EventSchedule.GetMergeData();

	int	firstEventIdx = 0;
	while( (firstEventIdx < eventDates->size()) && ((*eventDates)[firstEventIdx] < asOfDate) )
		++firstEventIdx;

	if ((firstEventIdx>0) && const_cast<ARM_HybridIRFXCalculator*>(this)->GetFixings() != NULL)
		firstEventIdx--;

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
///	Class  : ARM_TARNFXCalculator
///	Routine: FXValue
///	Returns: string
///	Action : fill the FX value(s) to pass to COUPON
/////////////////////////////////////////////////////////////////
string ARM_HybridIRFXCalculator::FXValue(ARM_DateStripPtr dateStrip, size_t eventIdx, int FXNb)
{
	string fxSpotDesc;
	double fixing = 0.0;
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();

	if (itsNbPastFixings)
	{
		if (itsCpnTiming == K_ARREARS)
			eventIdx--;

		double prevReset = (*(dateStrip->GetResetDates()))[eventIdx];
		fixing = GetFixings()->GetFxFixing(ARM_Date(prevReset), 
										   itsForeignCcy[FXNb].GetCcyName(), 
										   itsDomesticCcy.GetCcyName() );

		if ((prevReset <= asOfDate) && (fixing == -1))
			ARM_THROW( ERR_INVALID_ARGUMENT, "past FX fixing expected at " + ARM_Date(prevReset).toString());
	}

	char value[50];
	if (fixing > 0)
		sprintf(value, "%f", fixing);
	else
	{
		if (itsIsNormal)
			sprintf(value, "SPOT(FOREX_%s/%s)", itsForeignCcy[FXNb].GetCcyName(), itsDomesticCcy.GetCcyName());
		else
			sprintf(value, "SPOT(FOREX_%s/%s)", itsDomesticCcy.GetCcyName(), itsForeignCcy[FXNb].GetCcyName());
	}
	fxSpotDesc = string(value);

	return fxSpotDesc;
}


/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: MiddleRows
///	Returns: ARM_RowInfo
///	Action : create a row of a deal description
/////////////////////////////////////////////////////////////////
ARM_RowInfo ARM_HybridIRFXCalculator::MiddleRows( size_t eventIdx, const ARM_DateStripCombiner& datesStructure ) const
{
	ARM_DateStripPtr	vDateStrip = datesStructure.GetDateStrip(0);
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	size_t	eventSize = vDateStrip->GetResetDates()->size();
	size_t	descSize = sizeof(ColNamesTable)/sizeof(ColNamesTable[0]);

	vector<string>				rowDescVec(descSize);
	vector<ARM_GP_VALUE_TYPE>	rowTypeVec(descSize, ARM_MISSING_TYPE); 

    /// Set default 0 value for each column to be able to sum it
    InitPriceableColumns(rowDescVec, rowTypeVec);
	int i = 0;

	double resetDate = (*(vDateStrip->GetResetDates()))[eventIdx];

	/// WARNING : can have only one past fixing !
	if ((eventIdx < vDateStrip->GetResetDates()->size()-1)
		&&
		((*(vDateStrip->GetResetDates()))[eventIdx+1] < asOfDate) )
		return ARM_RowInfo();

	double payDate = (*(vDateStrip->GetPaymentDates()))[eventIdx];

	bool isArrearsCpn = (itsCpnTiming == K_ARREARS);
	bool isRealFirstEvent = (eventIdx == 0);
	bool isFirstEvent = (eventIdx == itsFirstEventIdx);
	bool isFirstCpnEvent = (isArrearsCpn ? (eventIdx == (itsFirstEventIdx+1)):(eventIdx == itsFirstEventIdx));
	bool isLastEvent = (eventIdx == (eventSize-1));

    CC_Ostringstream	resetDateDesc;
	double newResetDate = resetDate;
    resetDateDesc << CC_NS(std, fixed) << newResetDate;
    rowDescVec[ColAlias::ResetDate] = resetDateDesc.str();
    rowTypeVec[ColAlias::ResetDate] = ARM_DATE_TYPE;

    CC_Ostringstream	nextResetDateDesc;
	if( !(isArrearsCpn && isLastEvent) )
	{
		if(!isArrearsCpn)
			nextResetDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetResetDates()))[eventIdx];
		else if(!isLastEvent)
			nextResetDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetResetDates()))[eventIdx+1];

		rowDescVec[ColAlias::NextResetDate] = nextResetDateDesc.str();
		rowTypeVec[ColAlias::NextResetDate] = ARM_DATE_TYPE;

		CC_Ostringstream	fundingStartDateDesc;
		double fundingStart;
		if(!isArrearsCpn)
			fundingStart = (*(vDateStrip->GetFlowStartDates()))[eventIdx];
		else
			fundingStart = (*(vDateStrip->GetFlowStartDates()))[eventIdx+1];

		fundingStartDateDesc << CC_NS(std, fixed) << fundingStart;
		rowDescVec[ColAlias::FundingStartDate] = fundingStartDateDesc.str();
		rowTypeVec[ColAlias::FundingStartDate] = ARM_DATE_TYPE;

		CC_Ostringstream	fundingEndDateDesc;
		double fundingEnd;
		if(!isArrearsCpn)
			fundingEnd = (*(vDateStrip->GetFlowEndDates()))[eventIdx];
		else
			fundingEnd = (*(vDateStrip->GetFlowEndDates()))[eventIdx+1];

		fundingEndDateDesc << CC_NS(std, fixed) << fundingEnd;
		rowDescVec[ColAlias::FundingEndDate] = fundingEndDateDesc.str();
		rowTypeVec[ColAlias::FundingEndDate] = ARM_DATE_TYPE;

		double IT = 0.0;
		double pastValue = const_cast<ARM_HybridIRFXCalculator*>(this)->PastLiborValue(vDateStrip, IT, eventIdx, rowDescVec, rowTypeVec);

		CC_Ostringstream	fixedFundingDesc;
		if (pastValue > 0)
		{
			//don't fill IT if no past value
			CC_Ostringstream	ITDesc;
			ITDesc << CC_NS(std, fixed) << IT;
			rowDescVec[ColAlias::IT] = ITDesc.str();
			rowTypeVec[ColAlias::IT] = ARM_DOUBLE;

			fixedFundingDesc << pastValue;
		}
		else
			fixedFundingDesc << "0";

		rowDescVec[ColAlias::FixedFunding] = fixedFundingDesc.str();
		rowTypeVec[ColAlias::FixedFunding] = ARM_STRING;

		CC_Ostringstream	fundingDesc;
		if ((asOfDate == payDate)
			&&
			(const_cast<ARM_HybridIRFXCalculator*>(this)->GetDiscPricingMode() == ARM_DISC_ACCOUNTING_METH))
		{
			fundingDesc << "0";
		}
		else
		{
			if (isRealFirstEvent && (pastValue > 0))
				fundingDesc << "FixedFunding[i]";
			else
				fundingDesc << "FixedFunding[i]+SWAP(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
							<< ColNamesTable[ColAlias::FundingStartDate] << "[i]" << ","
							<< ColNamesTable[ColAlias::FundingEndDate] << "[i]" << ","
							<< "0" << ","
							<< ARM_ArgConvReverse_RcvOrPay.GetString(itsPayRec) << ","
							<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsCpnFreq) << ","
							<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsCpnDayCount) << ","
							<< ARM_ArgConvReverse_LgNameFrequency.GetString(itsFundFreq) << ","
							<< ARM_ArgConvReverse_LgNameDayCount.GetString(itsFundDayCount) << ","
							<< "FundMargin,FundNominal,,,FundDateStrip,FundDateStrip," 
							<< ColNamesTable[ColAlias::Offset] << "[i]," 
							<< ColNamesTable[ColAlias::Offset] << "[i],"
							<< ColNamesTable[ColAlias::Width] << "[i],"
							<< ColNamesTable[ColAlias::Width] << "[i])";
		}

		rowDescVec[ColAlias::Funding] = fundingDesc.str();
		rowTypeVec[ColAlias::Funding] = ARM_STRING;

		CC_Ostringstream	discountFundingDesc;
		discountFundingDesc	<< ColNamesTable[ColAlias::Funding] << "[i]";
		discountFundingDesc << "/" << "DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
							<< ColNamesTable[ColAlias::NextResetDate] << "[i])";

		rowDescVec[ColAlias::DiscountFunding] = discountFundingDesc.str();
		rowTypeVec[ColAlias::DiscountFunding] = ARM_STRING;

		CC_Ostringstream offsetDesc;
		offsetDesc << itsvFundIndex[eventIdx] - itsNbPastFixings << endl;

		rowDescVec[ColAlias::Offset] = offsetDesc.str();
		rowTypeVec[ColAlias::Offset] = ARM_DOUBLE;

		size_t width = itsvFundIndex[eventIdx+1] -itsvFundIndex[eventIdx];
		CC_Ostringstream widthDesc; 
		widthDesc << CC_NS(std, fixed) << width;

		rowDescVec[ColAlias::Width] = widthDesc.str();
		rowTypeVec[ColAlias::Width] = ARM_DOUBLE;
	}

	if( !(isArrearsCpn && isFirstEvent) )
	{
		CC_Ostringstream	startDateDesc;
		double	vFlowStartDate = (*(vDateStrip->GetFlowStartDates()))[eventIdx];
		startDateDesc << CC_NS(std, fixed) << vFlowStartDate;
		rowDescVec[ColAlias::StartDate] = startDateDesc.str();
		rowTypeVec[ColAlias::StartDate] = ARM_DATE_TYPE;

		CC_Ostringstream	endDateDesc;
		endDateDesc << CC_NS(std, fixed) << (*(vDateStrip->GetFlowEndDates()))[eventIdx];
		rowDescVec[ColAlias::EndDate] = endDateDesc.str();
		rowTypeVec[ColAlias::EndDate] = ARM_DATE_TYPE;

		CC_Ostringstream	ITDesc;
		ITDesc << CC_NS(std, fixed) << (*(vDateStrip->GetInterestTerms()))[eventIdx];
		rowDescVec[ColAlias::IT] = ITDesc.str();
		rowTypeVec[ColAlias::IT] = ARM_DOUBLE;

		CC_Ostringstream	redemptionResetDateDesc;
		if(isLastEvent)
		{
			redemptionResetDateDesc << CC_NS(std, fixed) << itsRedemptionResetDate.GetJulian();
			rowDescVec[ColAlias::RedemptionResetDate] = redemptionResetDateDesc.str();
			rowTypeVec[ColAlias::RedemptionResetDate] = ARM_DATE_TYPE;
		}

		CC_Ostringstream	minCpnDesc;
		minCpnDesc << CC_NS(std, fixed) << itsMinCpnCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::MinCpn] = minCpnDesc.str();
		rowTypeVec[ColAlias::MinCpn] = ARM_DOUBLE;

		CC_Ostringstream	maxCpnDesc;
		maxCpnDesc << CC_NS(std, fixed) << itsMaxCpnCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::MaxCpn] = maxCpnDesc.str();
		rowTypeVec[ColAlias::MaxCpn] = ARM_DOUBLE;

		CC_Ostringstream	notionalDesc;
		notionalDesc << CC_NS(std, fixed) << itsCpnNominalCv.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::Notional] = notionalDesc.str();
		rowTypeVec[ColAlias::Notional] = ARM_DOUBLE;

		for (i=0; i<itsNbFX; i++)
		{
			CC_Ostringstream	fxSpotDesc;
			fxSpotDesc	<< const_cast<ARM_HybridIRFXCalculator*>(this)->FXValue(vDateStrip, eventIdx, i);
			rowDescVec[ColAlias::FxSpot1 + i] = fxSpotDesc.str();
			rowTypeVec[ColAlias::FxSpot1 + i] = ARM_STRING;
		}
		
		const_cast<ARM_HybridIRFXCalculator*>(this)->DoPastReset(vDateStrip, eventIdx, rowDescVec, rowTypeVec);

		CC_Ostringstream	feesDesc;
		feesDesc << CC_NS(std, fixed) << itsFees.Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::Fees] = feesDesc.str();
		rowTypeVec[ColAlias::Fees] = ARM_DOUBLE;

		CC_Ostringstream	DFDesc;
		DFDesc	<< "DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
				<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")";

		rowDescVec[ColAlias::DF] = DFDesc.str();
		rowTypeVec[ColAlias::DF] = ARM_STRING;

		CC_Ostringstream	initialFXDesc;
		initialFXDesc << CC_NS(std, fixed) << itsInitialFXCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::InitFX] = initialFXDesc.str();
		rowTypeVec[ColAlias::InitFX] = ARM_DOUBLE;

		CC_Ostringstream	domesticCpnDesc;
		domesticCpnDesc << CC_NS(std, fixed) << itsDomesticCpnCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::DomCpn] = domesticCpnDesc.str();
		rowTypeVec[ColAlias::DomCpn] = ARM_DOUBLE;

		CC_Ostringstream	foreignCpnDesc;
		foreignCpnDesc << CC_NS(std, fixed) << itsForeignCpnCv[0].Interpolate(vFlowStartDate-asOfDate);
		rowDescVec[ColAlias::FgnCpn] = foreignCpnDesc.str();
		rowTypeVec[ColAlias::FgnCpn] = ARM_DOUBLE;

		CC_Ostringstream	couponDesc;
		if ((asOfDate == payDate)
			&&
			(const_cast<ARM_HybridIRFXCalculator*>(this)->GetDiscPricingMode() == ARM_DISC_ACCOUNTING_METH))
		{
			couponDesc << "0";
		}
		else
		{
			couponDesc	<< "Min(Max(" ;
			if ((itsNbFX == 1) || (itsFXChoice == "First") || itsPayOffName == ARM_TARNFXPayoffType::TARNFX)
			{
				couponDesc	<< ColNamesTable[ColAlias::FgnCpn] << "[i]"
							<< "*" << ColNamesTable[ColAlias::FxSpot1] << "[i]"
							<< "/" << ColNamesTable[ColAlias::InitFX] << "[i]" 
							<< "-" << ColNamesTable[ColAlias::DomCpn] << "[i]";
			}
			else if (itsFXChoice == "Second")
			{
				couponDesc	<< ColNamesTable[ColAlias::FgnCpn2] << "[i]" 
							<< "*" << ColNamesTable[ColAlias::FxSpot2] << "[i]"
							<< "/" << ColNamesTable[ColAlias::InitFX2] << "[i]" 
							<< "-" << ColNamesTable[ColAlias::DomCpn2] << "[i]";
			}
			else
			{
				string compare = (itsFXChoice == "Worst" ?	"Min(" : "Max(" );
				couponDesc << compare;

				couponDesc	<< ColNamesTable[ColAlias::FgnCpn] << "[i]" 
							<< "*" << ColNamesTable[ColAlias::FxSpot1] << "[i]"
							<< "/" << ColNamesTable[ColAlias::InitFX] << "[i]" 
							<< "-" << ColNamesTable[ColAlias::DomCpn] << "[i],";

				for (i=1; i<itsNbFX; i++)
				{
					if (i != itsNbFX-1)
						couponDesc << compare;

					couponDesc	<< ColNamesTable[ColAlias::FgnCpn2 + 4*(i-1)] << "[i]" 
								<< "*" << ColNamesTable[ColAlias::FxSpot2 + (i-1)] << "[i]"
								<< "/" << ColNamesTable[ColAlias::InitFX2 + 4*(i-1)] << "[i]" 
								<< "-" << ColNamesTable[ColAlias::DomCpn2 + 4*(i-1)] << "[i]";

					couponDesc << ((i == itsNbFX-1) ? ")" : ",");
				}
				for (i=1; i<itsNbFX-1; i++)
					couponDesc << ")";
			}
			couponDesc	<< "," << ColNamesTable[ColAlias::MinCpn] << "[i]" << "),"
						<< ColNamesTable[ColAlias::MaxCpn] << "[i]" << ")*";

			// NB : check if there is an effective date !
			if (isFirstCpnEvent && (itsEffectiveDate != ARM_DEFAULT_DATE))
				couponDesc	<< CountYears(itsCpnDayCount, itsEffectiveDate.GetJulian(), (*(vDateStrip->GetFlowEndDates()))[eventIdx]);
			else
				couponDesc	<< ColNamesTable[ColAlias::IT] << "[i]";
		}

		rowDescVec[ColAlias::Coupon] = couponDesc.str();
		rowTypeVec[ColAlias::Coupon] = ARM_STRING;

		CC_Ostringstream	paidCouponDesc;
		paidCouponDesc	<< ColNamesTable[ColAlias::Coupon] << "[i]"
						<< "*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")*"
						<< ColNamesTable[ColAlias::Notional] << "[i]";

		if(isLastEvent)
		{
			paidCouponDesc << "+" << ColNamesTable[ColAlias::Redemption] << "[i]";
		}

		rowDescVec[ColAlias::PaidCoupon] = paidCouponDesc.str();
		rowTypeVec[ColAlias::PaidCoupon] = ARM_STRING;

		CC_Ostringstream	sumCouponsDesc;
		sumCouponsDesc << ColNamesTable[ColAlias::Coupon] << "[i]";
		if (!isFirstCpnEvent)
		{
			sumCouponsDesc	<< "+" << ColNamesTable[ColAlias::SumCoupons] << "[i-1]";
		}
		rowDescVec[ColAlias::SumCoupons] = sumCouponsDesc.str();
		rowTypeVec[ColAlias::SumCoupons] = ARM_STRING;
		
		bool isPRDKO = itsPayOffName == ARM_TARNFXPayoffType::PRDKO;
		int index  = isPRDKO ? ColAlias::FxSpot1 : ColAlias::SumCoupons;

		CC_Ostringstream	isAliveDesc;
		isAliveDesc	<< "IF(" << ColNamesTable[index] << "[i]+" 
						<< ColNamesTable[ColAlias::Fees] << "[i]<"
						<< ColNamesTable[ColAlias::Target] << "[i],1,0)";
		if(!isFirstCpnEvent)
		{
			isAliveDesc	<< "*" << ColNamesTable[ColAlias::IsAlive] << "[i-1]";
		}

		rowDescVec[ColAlias::IsAlive] = isAliveDesc.str();
		rowTypeVec[ColAlias::IsAlive] = ARM_STRING;

		CC_Ostringstream	realCouponDesc;
		realCouponDesc	<< ColNamesTable[ColAlias::PaidCoupon] << "[i]";
		if(!isFirstCpnEvent)
		{
			realCouponDesc	<< "*" << ColNamesTable[ColAlias::IsAlive] << "[i-1]";
		}

		rowDescVec[ColAlias::RealCoupon] = realCouponDesc.str();
		rowTypeVec[ColAlias::RealCoupon] = ARM_STRING;

		CC_Ostringstream	realFundingDesc;

		if(!isFirstCpnEvent)
		{
			realFundingDesc	<< ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*";
		}

		if (isArrearsCpn)
		{
			if (itsFirstEventIdx < eventIdx)
				realFundingDesc	<< ColNamesTable[ColAlias::DiscountFunding] << "[i-1]";
			else
				realFundingDesc	<< "0";
		}
		else
			realFundingDesc	<< ColNamesTable[ColAlias::DiscountFunding] << "[i]";

		rowDescVec[ColAlias::RealFunding] = realFundingDesc.str();
		rowTypeVec[ColAlias::RealFunding] = ARM_STRING;

		CC_Ostringstream	capDesc;

		if(!isFirstCpnEvent)
		{
			capDesc	<< ColNamesTable[ColAlias::IsAlive] << "[i-1]*";
		}

		capDesc	<< "Max("
				<< ColNamesTable[ColAlias::SumCoupons] << "[i]" << "-"
				<< ColNamesTable[ColAlias::Target] << "[i]"
				<< ",0)*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
				<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")*"
				<< ColNamesTable[ColAlias::Notional] << "[i]";

		rowDescVec[ColAlias::Cap] = capDesc.str();
		rowTypeVec[ColAlias::Cap] = ARM_STRING;

		CC_Ostringstream	floorDesc;		
		if(isLastEvent)
		{
			floorDesc	<< ColNamesTable[ColAlias::IsAlive] << "[i]*"
						<< "Max("
						<< ColNamesTable[ColAlias::Target] << "[i]" << "-"
						<< ColNamesTable[ColAlias::SumCoupons] << "[i]"
						<< ",0)*DF(" << YC_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
						<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")*"
						<< ColNamesTable[ColAlias::Notional] << "[i]";
			rowDescVec[ColAlias::Floor] = floorDesc.str();
			rowTypeVec[ColAlias::Floor] = ARM_STRING;
		}
		
		CC_Ostringstream	redemptionStrikeDesc;
		double	vRedemptionStrike = 0.;
		if(isLastEvent)
			vRedemptionStrike = itsRedemptionStrike.Elt(0);
		redemptionStrikeDesc << CC_NS(std, fixed) << vRedemptionStrike;
		rowDescVec[ColAlias::RedemptionStrike] = redemptionStrikeDesc.str();
		rowTypeVec[ColAlias::RedemptionStrike] = ARM_DOUBLE;

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
										<< ")/" << ColNamesTable[ColAlias::RedemptionStrike] << "[i]" << "-1)";
					}
					else if (itsFXChoice == "Second")
					{
						redemptionDesc	<< "(Spot(" << FOREX_KEY_NAME << itsForeignCcy[1].GetCcyName()
										<< "/" << itsDomesticCcy.GetCcyName()
										<< ")/" << ColNamesTable[ColAlias::RedemptionStrike2] << "[i]" << "-1)";
					}
					else
					{
						redemptionDesc	<< "Min(Spot(" << FOREX_KEY_NAME << itsForeignCcy[0].GetCcyName()
										<< "/" << itsDomesticCcy.GetCcyName()
										<< ")/" << ColNamesTable[ColAlias::RedemptionStrike] << "[i]" << "-1,";
						
						for (i=1; i<itsNbFX; i++)
						{
							if (i != itsNbFX-1)
								redemptionDesc	<< "Min(";

							redemptionDesc	<< "Spot(" << FOREX_KEY_NAME << itsForeignCcy[i].GetCcyName()
											<< "/" << itsDomesticCcy.GetCcyName()
											<< ")/" << ColNamesTable[ColAlias::RedemptionStrike2+4*(i-1)] << "[i]" << "-1";

							if (i == itsNbFX-1)
								redemptionDesc	<< ")";
							else
								redemptionDesc	<< ",";
						}
						for (i=1; i<itsNbFX-1; i++)
							couponDesc << ")";
					}
					redemptionDesc	<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
									<< ColNamesTable[ColAlias::EndDate] << "[i])*"
									<< ColNamesTable[ColAlias::Notional] << "[i]";
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
					redemptionDesc	<< "(FWD(" << FOREX_KEY_NAME << itsForeignCcy[0].GetCcyName()
									<< "/" << itsDomesticCcy.GetCcyName() << ","
									<< ColNamesTable[ColAlias::RedemptionResetDate] << "[i]" << ", ,"
									<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")"
									<< "/" << ColNamesTable[ColAlias::RedemptionStrike] << "[i]" << "-1)"
									<< "*DF(" << YC_BASIS_KEY_NAME << itsDomesticCcy.GetCcyName() << ","
									<< ColNamesTable[ColAlias::EndDate] << "[i]" << ")"
									<< "*" << ColNamesTable[ColAlias::Notional] << "[i]";
				}
				else
				{
					redemptionDesc	<< "0";
				}
			}
			rowDescVec[ColAlias::Redemption] = redemptionDesc.str();
			rowTypeVec[ColAlias::Redemption] = ARM_STRING;
		}

		CC_Ostringstream	realRedemptionDesc;
		if(isLastEvent)
		{
			if ((itsRedemptionType == ARM_PRCSRedemptionType::mandatoryRedemption) || (itsNbFX > 1))
			{
				if (itsNbFX > 1)
					realRedemptionDesc	<< "-";
				
				realRedemptionDesc	<< ColNamesTable[ColAlias::IsAlive] << "[i]*"
									<< ColNamesTable[ColAlias::Redemption] << "[i]";
			}
			else
			{
				realRedemptionDesc	<< "0";
			}

			rowDescVec[ColAlias::RealRedemption] = realRedemptionDesc.str();
			rowTypeVec[ColAlias::RealRedemption] = ARM_STRING;
		}
		

		CC_Ostringstream	durationDesc;
		if(isFirstCpnEvent)
		{
			durationDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IT] << "[i]" << ")";
		}
		else
		{
			durationDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*"
							<< ColNamesTable[ColAlias::IT] << "[i]" << ")";
		}
		rowDescVec[ColAlias::Duration] = durationDesc.str();
		rowTypeVec[ColAlias::Duration] = ARM_STRING;

		CC_Ostringstream	probaDesc;
		if(isFirstCpnEvent)
		{
			probaDesc	<< "UNPAY(1-" << ColNamesTable[ColAlias::IsAlive] << "[i]" << ")";
		}
		else if(isLastEvent)
		{
			probaDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*(1-"
						<< ColNamesTable[ColAlias::IsAlive] << "[i]" << ")"
						<< "+" << ColNamesTable[ColAlias::IsAlive] << "[i])";
		}
		else
		{
			probaDesc	<< "UNPAY(" << ColNamesTable[ColAlias::IsAlive] << "[i-1]" << "*(1-"
							<< ColNamesTable[ColAlias::IsAlive] << "[i]" << "))";
		}
		rowDescVec[ColAlias::Proba] = probaDesc.str();
		rowTypeVec[ColAlias::Proba] = ARM_STRING;

		CC_Ostringstream	tarnDesc;
		if(itsPayRec == K_PAY)
		{
			tarnDesc	<< ColNamesTable[ColAlias::RealFunding] << "[i]-" << ColNamesTable[ColAlias::RealCoupon] << "[i]";
		}
		else
		{
			tarnDesc	<< ColNamesTable[ColAlias::RealCoupon] << "[i]-" << ColNamesTable[ColAlias::RealFunding] << "[i]";
		}
		rowDescVec[ColAlias::TARN] = tarnDesc.str();
		rowTypeVec[ColAlias::TARN] = ARM_STRING;

		int i;
		for (i=1; i<itsNbFX; i++)
		{
			CC_Ostringstream	initialFXDesc;
			initialFXDesc << CC_NS(std, fixed) << itsInitialFXCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[ColAlias::InitFX2 + 4*(i-1)] = initialFXDesc.str();
			rowTypeVec[ColAlias::InitFX2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	domesticCpnDesc;
			domesticCpnDesc << CC_NS(std, fixed) << itsDomesticCpnCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[ColAlias::DomCpn2 + 4*(i-1)] = domesticCpnDesc.str();
			rowTypeVec[ColAlias::DomCpn2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	foreignCpnDesc;
			foreignCpnDesc << CC_NS(std, fixed) << itsForeignCpnCv[i].Interpolate(vFlowStartDate-asOfDate);
			rowDescVec[ColAlias::FgnCpn2 + 4*(i-1)] = foreignCpnDesc.str();
			rowTypeVec[ColAlias::FgnCpn2 + 4*(i-1)] = ARM_DOUBLE;

			CC_Ostringstream	redemptionStrikeDesc;
			double	vRedemptionStrike = 0.;
			if(isLastEvent)
				vRedemptionStrike = itsRedemptionStrike.Elt(i);
			redemptionStrikeDesc << CC_NS(std, fixed) << vRedemptionStrike;
			rowDescVec[ColAlias::RedemptionStrike2 + 4*(i-1)] = redemptionStrikeDesc.str();
			rowTypeVec[ColAlias::RedemptionStrike2 + 4*(i-1)] = ARM_DOUBLE;
		}
	}

	return ARM_RowInfo(rowDescVec, rowTypeVec);
}

/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateModel2IRFX
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

	if ((!itsIsNormal) && (correlMatrix->size() == 9))
	{
		correlMatrix->Elt(0, 2) *= -1.0;
		correlMatrix->Elt(1, 2) *= -1.0;
		correlMatrix->Elt(2, 0) *= -1.0;
		correlMatrix->Elt(2, 1) *= -1.0;
	}

//	if (!itsIsNormal)
//		ForexKey = FOREX_KEY_NAME + itsDomesticCcy.GetCcyName() + "/" + itsForeignCcy[0].GetCcyName();
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateModel1IRFX
///	Returns: ARM_PricingModelPtr
///	Action : create the 1IR+FX model
/////////////////////////////////////////////////////////////////
ARM_1IRFXModel* ARM_HybridIRFXCalculator::CreateModel1IRFX()
{
	/// Create the 1IR+FX model

	/// Multi-assets avec modèles locaux sur le Fx si nécessaire
	int nbModels = ARM_1IRFXModel::DomBasisModel+1;

	string YcDomKey = YC_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcForKey = YC_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string YcBasisDomKey = YC_BASIS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string YcBasisForKey = YC_BASIS_KEY_NAME + itsForeignCcy[0].GetCcyName();
	string ForexKey = FOREX_KEY_NAME + itsForeignCcy[0].GetCcyName() + "/" + itsDomesticCcy.GetCcyName();
	string MrsDomKey = MRS_KEY_NAME + itsDomesticCcy.GetCcyName();
	string MrsForKey = MRS_KEY_NAME + itsForeignCcy[0].GetCcyName();

	ARM_StringVector names(nbModels);
	ARM_StringVectorVector depends(nbModels);
	if (itsIsNormal)
	{
		names[ARM_1IRFXModel::DomModel]        = YcDomKey;
		names[ARM_1IRFXModel::FxModel]         = ForexKey;	
		names[ARM_1IRFXModel::DomBasisModel]   = YcBasisDomKey;
	}
	else
	{
		names[ARM_1IRFXModel::DomModel]        = YcForKey;
//		names[ARM_1IRFXModel::FxModel]         = ForexKey;	
		names[ARM_1IRFXModel::FxModel]         = FOREX_KEY_NAME + itsDomesticCcy.GetCcyName() + "/" + itsForeignCcy[0].GetCcyName();	
		names[ARM_1IRFXModel::DomBasisModel]   = YcBasisForKey;
	}
	depends[ARM_1IRFXModel::FxModel]       = ARM_StringVector(1);
	depends[ARM_1IRFXModel::FxModel][ARM_1IRFXModel::DomModel]	   = names[ARM_1IRFXModel::DomBasisModel];
	depends[ARM_1IRFXModel::DomBasisModel] = ARM_StringVector(1,names[ARM_1IRFXModel::DomModel]);

	vector< ARM_PricingModelPtr > models(nbModels);
	/// Q vol without any curve because will be bootstrapped latter
	ARM_CurveModelParam HWVolParam( ARM_ModelParamType::QVol,HWVOL_LOWER_BOUND,"QVOL");
	ARM_ModelParamVector IRModelParams(3);
	IRModelParams[0]  = &HWVolParam;
	IRModelParams[1] = dynamic_cast<ARM_CurveModelParam*>(GetMktDataManager()->GetData(MrsDomKey));
	IRModelParams[2] = new ARM_CurveModelParam(ARM_ModelParamType::QParameter, 0.0,"Q");

	/// Create the Q1F model for domestic IR market
	ARM_QModel1F* ycModel = NULL;
	if (itsIsNormal)
	{
		ARM_ZeroCurve* ycDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcDomKey));
		ycModel = new ARM_QModel1F( CreateClonedPtr(ycDomCurve),ARM_ModelParamsQ1F(IRModelParams),true);
	}
	else
	{
		ARM_ZeroCurve* ycForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcForKey));
		ycModel = new ARM_QModel1F( CreateClonedPtr(ycForCurve),ARM_ModelParamsQ1F(IRModelParams),true);
	}
	models[ARM_1IRFXModel::DomModel] = ARM_PricingModelPtr( ycModel );

	/// Create the smiledFX model
	ARM_ZeroCurve* ycBasisDomCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisDomKey));
	ARM_ZeroCurve* ycBasisForCurve = dynamic_cast<ARM_ZeroCurve*>(GetMktDataManager()->GetData(YcBasisForKey));
//	if (!itsIsNormal)
//		ForexKey = FOREX_KEY_NAME + itsDomesticCcy.GetCcyName() + "/" + itsForeignCcy[0].GetCcyName();
	ARM_Forex* forex = dynamic_cast<ARM_Forex*>(GetMktDataManager()->GetData(ForexKey));

	double FXSpot = forex->GetMarketPrice();

	ARM_ModelParamVector FXModelParams;
	ARM_CurveModelParam volParam( ARM_ModelParamType::QVol,0.0,"QVol");
	ARM_CurveModelParam humpParam( ARM_ModelParamType::Hump,1.0,"Hump");
	ARM_CurveModelParam correlParam( ARM_ModelParamType::BetaCorrelation,0.0,"Beta");
	
	if ((itsPayOffName == ARM_TARNFXPayoffType::TARNFX)
		||
		(itsPayOffName == ARM_TARNFXPayoffType::CHOOSERFX))
		FXModelParams.push_back(&volParam);

	FXModelParams.push_back(&humpParam);
	FXModelParams.push_back(&correlParam);

	ARM_ModelParamsSmiled_Fx modelParamsObj(
			FXModelParams,
			ARM_ZeroCurvePtr(CreateClonedPtr(ycBasisDomCurve)),
			ARM_ZeroCurvePtr(CreateClonedPtr(ycBasisForCurve)),
			FXSpot,
			(itsCorrelType==ARM_ModelParamsSmiled::Fwd?2:itsFactorsNb),
			(ARM_ModelParamsSmiled::CorrelType) itsCorrelType,
			0.0);

	ARM_PricingModelPtr mod( new ARM_SmiledModel_Fx(CreateClonedPtr( ycBasisDomCurve ),&modelParamsObj,itsTimeStepNb,itsSpaceStepNb,itsStdDevNb,(itsSkipPDE?true:false),(itsRescalling?true:false)));

	models[ARM_1IRFXModel::FxModel] = mod;

	/// Create both domestic & foreign forward margin models
	ARM_ForwardMarginBasis* basisModel = NULL;
	if (itsIsNormal)
	{
		basisModel = new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisDomCurve));
	}
	else
	{
		basisModel = new ARM_ForwardMarginBasis(CreateClonedPtr(ycBasisForCurve));
	}
	models[ARM_1IRFXModel::DomBasisModel] = ARM_PricingModelPtr(basisModel);
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateModelNP1IRNFX
///	Returns: ARM_PricingModelPtr
///	Action : create the 1IR+nFX model
/////////////////////////////////////////////////////////////////
ARM_NP1IRNFXModel* ARM_HybridIRFXCalculator::CreateModelNP1IRNFX()
{
	const double initIRVol = 0.01;
	const double initFXVol = 0.1;
	const double initIRQ = 0.0;
	const double initFXQ = 1.0;
	std::vector<double> lowerBound(1, 0.0001);
	std::vector<double> upperBound(1, 0.05);
	std::vector<double> lowerQ(1, -10.0);
	std::vector<double> upperQ(1, 10.0);

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
	std::vector<double> breakPointTimes(1, 0.0);
	std::vector<double> values(1, initIRVol );

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


/////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateAndSetModel
///	Returns: void
///	Action : create the IR+FX model
/////////////////////////////////////////////////////////////////
void ARM_HybridIRFXCalculator::CreateAndSetModel()
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
	int scheme = ARM_PathScheme::BrownianBridge;
	if ((itsPayOffName == ARM_TARNFXPayoffType::INDIANFX) || (itsPayOffName == ARM_TARNFXPayoffType::SWITCHER))
		scheme = ARM_PathScheme::Incremental;

	ARM_PathSchemePtr pathScheme(ARM_PathSchemeFactory.Instance()->CreatePathScheme(scheme));
	ARM_MCMethod* mcMethod = new ARM_MCMethod(itsNbSimul,randGen,&sampler,itsBucketSize,ARM_ImpSamplerPtr(NULL),pathScheme);

	numModel->SetNumMethod(ARM_NumMethodPtr( mcMethod ) );
	ARM_NumerairePtr numeraire( ARM_NumeraireFactory.Instance()->CreateNumeraire( (itsModelType==ModelNP1IRNFX)||(itsModelType==Model2IRFX)||(itsCorrelType==ARM_ModelParamsSmiled::Fwd)?ARM_Numeraire::RollingCash:ARM_Numeraire::TerminalZc ) );
    numModel->SetNumeraire(numeraire);
	
	SetPricingModel(numModel);
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
    std::vector<double> stdExp(nbExp+1);
	size_t i;
    stdExp[0] = asOfDate.GetJulian();
    for(i=0;i<nbExp;++i)
        stdExp[i+1] = asOfDate.GetJulian() + K_YEAR_LEN * (*(oswBSVol->GetExpiryTerms()))[i];

	// Get the model dates
	int nbFund = itsFundDateStrip->size();

	std::vector<double>& resetDates = itsFundDateStrip->GetResetDates();
	std::vector<double>& startDates = itsFundDateStrip->GetFlowStartDates();

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


    std::vector<double> initTimes(nbOSW);
    std::vector<double> initVols(nbOSW);
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
    std::vector<double> volLowerBound(nbOSW,HWVOL_LOWER_BOUND);
    std::vector<double> volUpperBound(nbOSW,HWVOL_UPPER_BOUND);
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
	
	std::vector<double>& resetDates = itsCpnDateStrip->GetResetDates();
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

    std::vector<double> initTimes(nbFx);
    std::vector<double> initVols(nbFx);

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
    std::vector<double> volLowerBound(nbFx,QVOL_LOWER_BOUND);
    std::vector<double> volUpperBound(nbFx,QVOL_UPPER_BOUND);
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


ARM_CalibMethod* ARM_HybridIRFXCalculator::CreateNumFxCalibMethod()
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

	std::vector<double>& startDates = itsCpnDateStrip->GetFlowStartDates();
	std::vector<double>& endDates = itsCpnDateStrip->GetFlowEndDates();
	std::vector<double>& fwdStartDates = itsCpnDateStrip->GetFwdRateStartDates();
	std::vector<double>& fwdEndDates = itsCpnDateStrip->GetFwdRateEndDates();
	std::vector<double>& resetDates = itsCpnDateStrip->GetResetDates();
	std::vector<double>& payDates = itsCpnDateStrip->GetPaymentDates();
	std::vector<double>& interestDays = itsCpnDateStrip->GetInterestDays();
	std::vector<double>& interestTerms = itsCpnDateStrip->GetInterestTerms();

	double volATMValue, decVolValue, shiftValue, lambdaValue;

	double fxSpot = forex->GetMarketPrice();
	double fwd, maturity;

	ARM_VanillaSecDensityPtrVector vanillaSecDensities;

	std::vector<double> calibStartDates;
	std::vector<double> calibEndDates;
	std::vector<double> calibFwdStartDates;
	std::vector<double> calibFwdEndDates;
	std::vector<double> calibResetDates;
	std::vector<double> calibPayDates;
	std::vector<double> calibInterestDays;
	std::vector<double> calibInterestTerms;

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
				lambdaValue,
				itsIsNormal));
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
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateCalibration
///	Returns: ARM_CalibMethod*
///	Action : create the calibration method
/////////////////////////////////////////////////////////////////
ARM_CalibMethod* ARM_HybridIRFXCalculator::CreateCalibration(ModelType modelType)
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

///////////////////////////////////////////////////////////////////////////
///	Class  : ARM_HybridIRFXCalculator
///	Routine: CreateDensityFunctor
///	Returns: vector<ARM_DensityFunctor*>
///	Action : create a vector of density functors for NP1IRNFX calibration
///////////////////////////////////////////////////////////////////////////
vector<ARM_DensityFunctor*> ARM_HybridIRFXCalculator::CreateDensityFunctor()
{
	double asOfDate = GetMktDataManager()->GetAsOfDate().GetJulian();
	int nbReset = itsCpnDateStrip->size();

	int start_i = (itsCpnTiming == K_ADVANCE) ? 0 : 1;
	start_i += itsNbPastFixings;
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

		std::vector<double>& startDates = itsCpnDateStrip->GetFlowStartDates();
		std::vector<double>& endDates = itsCpnDateStrip->GetFlowEndDates();
		std::vector<double>& fwdStartDates = itsCpnDateStrip->GetFwdRateStartDates();
		std::vector<double>& fwdEndDates = itsCpnDateStrip->GetFwdRateEndDates();
		std::vector<double>& resetDates = itsCpnDateStrip->GetResetDates();
		std::vector<double>& payDates = itsCpnDateStrip->GetPaymentDates();
		std::vector<double>& interestDays = itsCpnDateStrip->GetInterestDays();
		std::vector<double>& interestTerms = itsCpnDateStrip->GetInterestTerms();

		double volATMValue, decVolValue, shiftValue, lambdaValue;

		double fxSpot = forex->GetMarketPrice();
		double fwd, maturity;

		std::vector<double> calibStartDates;
		std::vector<double> calibEndDates;
		std::vector<double> calibFwdStartDates;
		std::vector<double> calibFwdEndDates;
		std::vector<double> calibResetDates;
		std::vector<double> calibPayDates;
		std::vector<double> calibInterestDays;
		std::vector<double> calibInterestTerms;

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

void ARM_HybridIRFXCalculator::Calibrate()
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
		i += itsNbPastFixings;
		std::vector<double> newResetTimes((&*itsCpnDateStrip)->GetResetDates()->begin()+i, (&*itsCpnDateStrip)->GetResetDates()->end());
		ARM_GP_VectorPtr ResetTimes( (std::vector<double>&)newResetTimes.Clone() );
		
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

	if ((itsPayOffName == ARM_TARNFXPayoffType::TARNFX)
		||
		(itsPayOffName == ARM_TARNFXPayoffType::CHOOSERFX)
	   )
		ComputeDomesticBasis();
}


////////////////////////////////////////////////////
///	Class   : ARM_HybridIRFXCalculator
///	Routines: toString
///	Returns :
///	Action  : Object dump
////////////////////////////////////////////////////
double ARM_HybridIRFXCalculator::Price()
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

		if (itsColumnsToPrice.begin()->find("Proba") && itsIntermediatePrices) 
		{
			GetPricingData()["Probas"] = genPricer->GetPricerInfo()->GetContents( ColNamesTable[ ProductToPriceColumns[ProbaValue] ]).GetData("IntermediatePrices").GetVector();
		}
	}
	itsHasBeenComputed = true;

    return firstColumnPrice;
}


void ARM_HybridIRFXCalculator::ComputePricingData() const
{
	if (!itsHasBeenComputed)
		const_cast<ARM_HybridIRFXCalculator*>(this)->PriceAndTimeIt();
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

